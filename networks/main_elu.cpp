#include "LinearOT/linear-ot.h"
#include "utils/emp-tool.h"
#include <iostream>

// #include "FloatingPoint/floating-point.h"
// #include "FloatingPoint/fp-math.h"
#include <random>
#include <limits>
// #include "float_utils.h"

#include "Millionaire/millionaire.h"
#include "Millionaire/millionaire_with_equality.h"
#include "BuildingBlocks/truncation.h"
#include "BuildingBlocks/aux-protocols.h"
#include <chrono>
// #include <matplotlibcpp.h>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include "OT/ot-utils.h"
#include "OT/emp-ot.h"

#define MAX_THREADS 4
using namespace sci;
using namespace std;
int party, port = 32000;
string address = "127.0.0.1";
int bwL = 21;
uint64_t mask_bwL = (bwL == 64 ? -1 : ((1ULL << bwL) - 1));
int bwL_1 = bwL - 1;
uint64_t mask_bwL_1 = (bwL_1 == 64 ? -1 : ((1ULL << bwL_1) - 1));
bool signed_B = true;
bool accumulate = true;
bool precomputed_MSBs = false;
MultMode mode = MultMode::None;
uint64_t lb = 10;
uint64_t la = 10; // la=5 f=5,la=14,f=12
uint64_t f = 12;
uint64_t s = 7;
// int dim = 4096 * 16;
uint64_t dim = 524288;
// uint64_t dim = 1048576;
// int dim = 100;
uint64_t acc = 2;
// uint64_t init_input = 1032192-10;
uint64_t init_input = 2064384;
// uint64_t init_input = 2080768;
// uint64_t init_input = 100;
uint64_t step_size = 1;
uint64_t correct = 1;

uint64_t h = f + 3;
uint64_t d = f + 3;
uint64_t Tk = f - 1;
uint64_t alpha = 3.5 * pow(2, f);
uint64_t mask_lb = (lb == 64 ? -1 : ((1ULL << lb) - 1));
uint64_t mask_l_Tk = (bwL == 64 ? -1 : ((1ULL << (bwL - Tk)) - 1));
uint64_t mask_lah1 = ((la + h + 1) == 64 ? -1 : ((1ULL << (la + h + 1)) - 1));
uint64_t mask_lla = ((la + bwL) == 64 ? -1 : ((1ULL << (la + bwL)) - 1));
// uint64_t s = 7;
uint64_t mask_s = ((s) == 64 ? -1 : ((1ULL << (s)) - 1));
uint64_t mask_h = (h == 64) ? ~0ULL : (1ULL << h) - 1;

LinearOT *prod;
XTProtocol *ext;
AuxProtocols *aux;
MillionaireWithEquality<sci::NetIO> *mill_eq;
Equality<sci::NetIO> *eq;
Truncation *trunc_oracle;

sci::NetIO *io;
sci::OTPack<sci::NetIO> *otpack;

double calculate_GELU(uint64_t value)
{
    const int64_t shift_amount = 64 - bwL;
    int64_t signed_value = static_cast<int64_t>(value << shift_amount) >> shift_amount;
    const double pow_2_f = static_cast<double>(1ULL << f);
    double x = static_cast<double>(signed_value) / pow_2_f;
    return 0.5 * x + 0.5 * x * std::erf(x / 1.414);
}

double calculate_tanh(uint64_t value, uint64_t f_tanh)
{

    const int64_t shift_amount = 64 - bwL;
    int64_t signed_value = static_cast<int64_t>(value << shift_amount) >> shift_amount;
    const double pow_2_f = static_cast<double>(1ULL << f_tanh);
    double x = static_cast<double>(signed_value) / pow_2_f;
    return std::tanh(x);
}

double calculate_sigmoid(uint64_t value, uint64_t f_sigmoid)
{
    const int64_t shift_amount = 64 - bwL;
    int64_t signed_value = static_cast<int64_t>(value << shift_amount) >> shift_amount;
    const double pow_2_f = static_cast<double>(1ULL << f_sigmoid);
    double x = static_cast<double>(signed_value) / pow_2_f;
    return 1.0 / (1.0 + std::exp(-x));
}

double calculate_ELU(uint64_t value, uint64_t f_ELU, double alpha = 1.0)
{

    const int64_t shift_amount = 64 - bwL;

    int64_t signed_value = static_cast<int64_t>(value << shift_amount) >> shift_amount;

    const double pow_2_f = static_cast<double>(1ULL << f_ELU);
    double x = static_cast<double>(signed_value) / pow_2_f;

    if (x > 0.0)
    {
        return x;
    }
    else
    {
        return alpha * (std::exp(x) - 1.0);
    }
}

int64_t decode_ring(uint64_t input, uint64_t bw)
{
    uint64_t mask = (bw == 64) ? ~0ULL : (1ULL << bw) - 1;
    uint64_t half = 1ULL << (bw - 1);

    // std::cout << "input = " << input << std::endl;
    // std::cout << "half = " << half << std::endl;
    if (input == 1048576)
    {
        return 0;
    }
    if (input < half)
    {
        return input;
    }
    else
    {
        return -((1ULL << (bw)) - input);
    }
}

void assign_lower_h_bits(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *input_lower_h, int32_t h)
{
    // Create a mask that has the lowest h bits set to 1
    uint64_t mask = (h == 64) ? ~0ULL : (1ULL << h) - 1;

    // Assign the lower h bits from inA to inA_

    if (party == sci::ALICE) // Assign the lower h bits from inB to inB_
        for (int i = 0; i < dim; i++)
        {
            input_lower_h[i] = inA[i] & mask;
        }
    else
    {
        for (int i = 0; i < dim; i++)
        {
            input_lower_h[i] = inB[i] & mask;
        }
    }
}

void select_share(uint8_t *sel, uint64_t *x, uint64_t *y, uint64_t *output, int32_t dim, int32_t h)
{
    // Create a mask that has the lowest h bits set to 1
    uint64_t mask = (h == 64) ? ~0ULL : (1ULL << h) - 1;
    uint64_t *mid = new uint64_t[dim];
    for (int i = 0; i < dim; i++)
    {
        mid[i] = (x[i] - y[i]) & mask;
    }

    aux->multiplexer(sel, mid, output, dim, h, h);

    for (int i = 0; i < dim; i++)
    {
        output[i] = (output[i] + y[i]) & mask;
    }
}

void DReLU_Eq(uint64_t *inA, uint8_t *b, uint8_t *b_, uint64_t dim, uint64_t bwl)
{
    uint8_t *m = new uint8_t[dim];
    uint64_t *y = new uint64_t[dim];
    uint64_t mask_l_sub1 = ((bwl - 1) == 64) ? ~0ULL : (1ULL << (bwl - 1)) - 1;
    for (int i = 0; i < dim; i++)
    {
        m[i] = inA[i] >> (bwl - 1);
        // std::cout << "m[" << i << "] = " << static_cast<int>(m[i]) << std::endl;
        y[i] = inA[i] & mask_l_sub1;
    }
    std::cout << "mask_l_sub1 = " << mask_l_sub1 << std::endl;
    uint64_t *comp_eq_input = new uint64_t[dim];

    if (party == sci::ALICE)
    {
        for (int i = 0; i < dim; i++)
        {
            comp_eq_input[i] = (mask_l_sub1 - y[i]) & mask_l_sub1;
            // std::cout << "inA[" << i << "] = " << inA[i] << std::endl;
            // std::cout << "y[" << i << "] = " << y[i] << std::endl;
            // std::cout << "comp_eq_input[" << i << "] = " << comp_eq_input[i] << std::endl;
        }
    }
    else
    {
        for (int i = 0; i < dim; i++)
        {
            comp_eq_input[i] = y[i] & mask_l_sub1;
            // std::cout << "inA[" << i << "] = " << inA[i] << std::endl;
            // std::cout << "y[" << i << "] = " << y[i] << std::endl;
            // std::cout << "comp_eq_input[" << i << "] = " << comp_eq_input[i] << std::endl;
        }
    }

    uint8_t *carry = new uint8_t[dim];
    uint8_t *res_eq = new uint8_t[dim];
    mill_eq->compare_with_eq(carry, res_eq, comp_eq_input, dim, bwl - 1, false);
    for (int i = 0; i < dim; i++)
    {
        // std::cout << "carry[" << i << "] = " << static_cast<int>(carry[i]) << std::endl;
        // std::cout << "res_eq[" << i << "] = " << static_cast<int>(res_eq[i]) << std::endl;
    }
    if (party == sci::ALICE)
    {
        for (int i = 0; i < dim; i++)
        {
            b[i] = carry[i] ^ 1 ^ m[i];
            // b[i] = carry[i] ^ m[i];
        }
    }
    else
    {
        for (int i = 0; i < dim; i++)
        {
            b[i] = carry[i] ^ m[i];
        }
    }

    aux->AND(res_eq, m, b_, dim);
}
int second_interval(uint64_t *input_data, uint8_t *res_drelu_cmp, uint8_t *res_drelu_eq)
{
    mill_eq = new MillionaireWithEquality(party, io, otpack);
    trunc_oracle = new Truncation(party, io, otpack);
    aux = new AuxProtocols(party, io, otpack);
    eq = new Equality(party, io, otpack);
    ext = new XTProtocol(party, io, otpack);
    uint64_t *comp_eq_input = new uint64_t[dim];
    uint64_t *outtrunc = new uint64_t[dim];
    uint8_t *res_eq = new uint8_t[dim];
    trunc_oracle->truncate_and_reduce(dim, input_data, outtrunc, h, bwL); // test comm
    DReLU_Eq(outtrunc, res_drelu_cmp, res_drelu_eq, dim, bwL - h);
    return 1;
}

void third_interval(uint64_t *input_data, uint8_t *res_drelu_cmp, uint8_t *res_drelu_eq, uint8_t *res_eq)
{
    // mill_eq = new MillionaireWithEquality(party, io, otpack);
    trunc_oracle = new Truncation(party, io, otpack);
    aux = new AuxProtocols(party, io, otpack);
    // eq = new Equality(party, io, otpack);
    ext = new XTProtocol(party, io, otpack);
    uint64_t *comp_eq_input = new uint64_t[dim];
    uint64_t *outtrunc = new uint64_t[dim];
    // uint8_t *res_drelu_cmp = new uint8_t[dim];
    // uint8_t *res_drelu_eq = new uint8_t[dim];
    // uint8_t *res_eq = new uint8_t[dim];
    uint8_t *res_cmp = new uint8_t[dim];
    // TR
    uint64_t Comm_start = io->counter;
    auto time_start = std::chrono::high_resolution_clock::now();
    // if (party == sci::ALICE)
    // {
    //     for (int i = 0; i < dim; i++)
    //     {
    //         input_data[i] = (input_data[i] - alpha) & mask_bwL;
    //     }
    // }
    uint64_t trun_start = io->counter;
    trunc_oracle->truncate_and_reduce(dim, input_data, outtrunc, d, bwL); // test comm
    DReLU_Eq(outtrunc, res_drelu_cmp, res_drelu_eq, dim, bwL - d);
    uint64_t mask_l_sub1 = ((bwL - d) == 64) ? ~0ULL : (1ULL << (bwL - d)) - 1;
    if (party == sci::ALICE)
    {
        for (int i = 0; i < dim; i++)
        {
            comp_eq_input[i] = (mask_l_sub1 + 1 - outtrunc[i]) & mask_l_sub1;
        }
    }
    else
    {
        for (int i = 0; i < dim; i++)
        {
            comp_eq_input[i] = outtrunc[i] & mask_l_sub1;
        }
    }
    eq->check_equality(res_eq, comp_eq_input, dim, bwL - d);
}

uint64_t computeULPErr(double calc, double actual, int SCALE)
{
    int64_t calc_fixed = (double(calc) * (1ULL << SCALE));
    int64_t actual_fixed = (double(actual) * (1ULL << SCALE));
    uint64_t ulp_err = (calc_fixed - actual_fixed) > 0
                           ? (calc_fixed - actual_fixed)
                           : (actual_fixed - calc_fixed);
    return ulp_err;
}

void cross_term_mul_batch(
    uint64_t batch_size,
    int party,
    const uint64_t *x,
    uint64_t *z_share, // 输出

    int m, int n, int ell,
    sci::OTPack<sci::NetIO> *otpack)
{
    // std::cout << "cross_term_mul_batch" << std::endl;
    sci::NetIO *io = otpack->io;
    uint64_t mask_m = (m == 64 ? -1 : ((1ULL << m) - 1));
    uint64_t mask_n = (n == 64 ? -1 : ((1ULL << n) - 1));
    uint64_t mask_ell = (ell == 64 ? -1 : ((1ULL << ell) - 1));
    for (int i = 0; i < batch_size; i++)
    {
        z_share[i] = 0;
    }
    if (m <= n)
    {
        std::cout << "m<=n" << std::endl;
        for (int i = 0; i < m; ++i)
        {
            uint64_t *data0 = new uint64_t[batch_size];
            uint64_t *corr = new uint64_t[batch_size];
            uint64_t *recv_data = new uint64_t[batch_size];
            uint64_t mask_ell_i = (ell - i == 64 ? -1 : ((1ULL << (ell - i)) - 1));
            uint8_t *choice = new uint8_t[batch_size];
            // std::cout << "o->counter;" << std::endl;
            size_t comm_before = io->counter;
            // std::cout << "comm_before" << std::endl;
            for (int b = 0; b < batch_size; b++)
            {
                data0[b] = 0;
                corr[b] = 0;
                recv_data[b] = 0;
            }
            if (party == sci::ALICE)
            {
                for (int b = 0; b < batch_size; b++)
                {
                    // corr[b] = x[b];
                    choice[b] = (x[b] >> i) & 1;
                }
                otpack->iknp_straight->recv_cot(recv_data, (bool *)choice, batch_size, ell - i);
                for (int b = 0; b < 5; b++)
                {
                    std::cout << "x[" << i << "] = " << x[b] << std::endl;
                    std::cout << "choice[" << i << "] = " << (static_cast<int>(choice[b])) << std::endl;
                    std::cout << "choice[" << i << "] = " << (bool)(choice[b]) << std::endl;
                    std::cout << "recv_data[" << i << "] = " << recv_data[b] << std::endl;
                    std::cout << "z_share[" << i << "] = " << z_share[b] << std::endl;
                }
                for (int b = 0; b < batch_size; b++)
                {
                    recv_data[b] = (-recv_data[b]) & mask_ell_i;
                    z_share[b] = (z_share[b] + (recv_data[b] << i)) & mask_ell;
                }
            }
            else
            {
                // std::cout << "ALICE" << std::endl;
                for (int b = 0; b < batch_size; b++)
                {
                    corr[b] = x[b];
                    // choice[b] = (x[b] >> i) & 1; // 只 BOB 用 x
                }
                otpack->iknp_straight->send_cot(data0, corr, batch_size, ell - i);
                // otpack->iknp_straight->recv_cot(recv_data, choice, batch_size, ell - i);
                for (int b = 0; b < 5; b++)
                {
                    std::cout << "x[" << i << "] = " << x[b] << std::endl;
                    std::cout << "corr[" << i << "] = " << corr[b] << std::endl;
                    std::cout << "data0[" << i << "] = " << data0[b] << std::endl;
                    std::cout << "z_share[" << i << "] = " << z_share[b] << std::endl;
                }
                for (int b = 0; b < batch_size; b++)
                {
                    // data0[b] = (-data0[b]) & mask_ell_i;
                    z_share[b] = (z_share[b] + (data0[b] << i)) & mask_ell;
                }
            }
            io->flush();
        }
    }
    else
    {
        std::cout << "n<m" << std::endl;
        for (int i = 0; i < n; ++i)
        {
            uint64_t *data0 = new uint64_t[batch_size];
            uint64_t *corr = new uint64_t[batch_size];
            uint64_t *recv_data = new uint64_t[batch_size];
            uint64_t mask_ell_i = (ell - i == 64 ? -1 : ((1ULL << (ell - i)) - 1));
            uint8_t *choice = new uint8_t[batch_size];
            // std::cout << "o->counter;" << std::endl;
            size_t comm_before = io->counter;
            // std::cout << "comm_before" << std::endl;
            for (int b = 0; b < batch_size; b++)
            {
                data0[b] = 0;
                corr[b] = 0;
                recv_data[b] = 0;
            }
            if (party != sci::ALICE)
            {
                // std::cout << "ALICE" << std::endl;
                std::cout << "party = " << party << std::endl;
                for (int b = 0; b < batch_size; b++)
                {
                    choice[b] = (x[b] >> i) & 1;
                }
                otpack->iknp_straight->recv_cot(recv_data, (bool *)choice, batch_size, ell - i);
                for (int b = 0; b < 5; b++)
                {
                    std::cout << "x[" << i << "] = " << x[b] << std::endl;
                    std::cout << "choice[" << i << "] = " << (static_cast<int>(choice[b])) << std::endl;
                    std::cout << "recv_data[" << i << "] = " << recv_data[b] << std::endl;
                    std::cout << "z_share[" << i << "] = " << z_share[b] << std::endl;
                }
                for (int b = 0; b < batch_size; b++)
                {
                    // recv_data[b] = (-recv_data[b]) & mask_ell_i;
                    z_share[b] = (z_share[b] + (recv_data[b] << i)) & mask_ell;
                }
            }
            else
            {
                for (int b = 0; b < batch_size; b++)
                {
                    corr[b] = x[b];
                }
                // std::cout << "send_cot" << std::endl;
                otpack->iknp_straight->send_cot(data0, corr, batch_size, ell - i);
                for (int b = 0; b < 5; b++)
                {
                    std::cout << "x[" << i << "] = " << x[b] << std::endl;
                    std::cout << "corr[" << i << "] = " << corr[b] << std::endl;
                    std::cout << "data0[" << i << "] = " << data0[b] << std::endl;
                    std::cout << "z_share[" << i << "] = " << z_share[b] << std::endl;
                }
                std::cout << "send_cot" << std::endl;
                for (int b = 0; b < batch_size; b++)
                {
                    data0[b] = (-data0[b]) & mask_ell_i;
                    z_share[b] = (z_share[b] + (data0[b] << i)) & mask_ell;
                    // std::cout << "z_share[" << i << "] = " << z_share[b] << std::endl;
                }
            }
            io->flush();
        }
    }
    std::cout << "batch_size = " << batch_size << std::endl;
    std::cout << "m = " << m << std::endl;
    std::cout << "n = " << n << std::endl;
    std::cout << "ell = " << ell << std::endl;
    for (int i = 0; i < 5; i++)
    {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
        std::cout << "z_share[" << i << "] = " << z_share[i] << std::endl;
    }
}

void cross_term(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                int32_t bwA, int32_t bwB, int32_t bwC)
{
    uint64_t mask = (1ULL << bwC) - 1;
    if (party == sci::ALICE)
    {
        sci::PRG128 prg;
        for (int i = 0; i < bwB; i++)
        {
            auto *data0 = new uint64_t[dim];
            prg.random_data(data0, dim * sizeof(uint64_t));
            for (int j = 0; j < dim; j++)
            {
                data0[j] = data0[j] & ((1ULL << (bwC - i)) - 1);
            }
            otpack->iknp_straight->send_cot(data0, inA, dim, bwC - i);
            for (int j = 0; j < dim; j++)
            {
                outC[j] += (-data0[j] * (1ULL << i));
                outC[j] &= mask;
            }
            delete[] data0;
        }
    }
    else
    {
        bool choice[bwB][dim];
        for (int i = 0; i < dim; i++)
        {
            bool *temp = new bool[bwB];
            sci::int64_to_bool(temp, inB[i], bwB);
            for (int j = 0; j < bwB; j++)
            {
                choice[j][i] = temp[j];
            }
            delete[] temp;
        }
        for (int i = 0; i < bwB; i++)
        {
            auto *data = new uint64_t[dim];
            bool *c = new bool[dim];
            for (int j = 0; j < dim; j++)
            {
                bool *temp = new bool[bwB];
                sci::int64_to_bool(temp, inB[j], bwB);
                c[j] = temp[i];
            }
            otpack->iknp_straight->recv_cot(data, c, dim, bwC - i);
            for (int j = 0; j < dim; j++)
            {
                outC[j] += (data[j] * (1ULL << i));
                outC[j] &= mask;
            }
            delete[] data;
            delete[] c;
        }
    }
}

void cross_term_reverse(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                        int32_t bwA, int32_t bwB, int32_t bwC)
{
    uint64_t mask = (1ULL << bwC) - 1;
    if (party == sci::ALICE)
    {
        bool choice[bwA][dim];
        for (int i = 0; i < dim; i++)
        {
            bool *temp = new bool[bwA];
            sci::int64_to_bool(temp, inA[i], bwA);
            for (int j = 0; j < bwA; j++)
            {
                choice[j][i] = temp[j];
            }
            delete[] temp;
        }
        for (int i = 0; i < bwA; i++)
        {
            auto *data = new uint64_t[dim];
            bool *c = new bool[dim];
            for (int j = 0; j < dim; j++)
            {
                bool *temp = new bool[bwA];
                sci::int64_to_bool(temp, inA[j], bwA);
                c[j] = temp[i];
            }
            otpack->iknp_reversed->recv_cot(data, c, dim, bwC - i);
            for (int j = 0; j < dim; j++)
            {
                outC[j] += ((data[j] * (1ULL << i)) & mask);
                outC[j] &= mask;
            }
            delete[] data;
            delete[] c;
        }
    }
    else
    {
        sci::PRG128 prg;
        for (int i = 0; i < bwA; i++)
        {
            auto *data0 = new uint64_t[dim];
            prg.random_data(data0, dim * sizeof(uint64_t));
            otpack->iknp_reversed->send_cot(data0, inB, dim, bwC - i);
            for (int j = 0; j < dim; j++)
            {
                outC[j] += ((-data0[j] * (1ULL << i)) & mask);
                outC[j] &= mask;
            }
            delete[] data0;
        }
    }
}

void sirnn_unsigned_mul(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                        int32_t bwA, int32_t bwB, int32_t bwC)
{
    auto *c = new uint64_t[dim]();
    auto *d = new uint64_t[dim];
    for (int i = 0; i < dim; i++)
    {
        c[i] = 0;
        d[i] = 0;
    }
    if (bwA <= bwB)
    {
        cross_term_reverse(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term(dim, inB, inA, d, bwB, bwA, bwC);
    }
    else
    {
        cross_term(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term_reverse(dim, inB, inA, d, bwB, bwA, bwC);
    }
    auto *wx = new uint8_t[dim];
    auto *wy = new uint8_t[dim];
    aux->wrap_computation(inA, wx, dim, bwA);
    aux->wrap_computation(inB, wy, dim, bwB);
    auto *h = new uint64_t[dim];
    auto *g = new uint64_t[dim];
    aux->multiplexer(wx, inB, h, dim, bwB, bwB);
    aux->multiplexer(wy, inA, g, dim, bwA, bwA);
    uint64_t mask = (1ULL << bwC) - 1;
    for (int i = 0; i < dim; i++)
    {
        outC[i] = (inA[i] * inB[i] + c[i] + d[i] - (g[i] * (1ULL << bwB)) - (h[i] * (1ULL << bwA))) & mask;
    }
    delete[] c;
    delete[] d;
}

void sirnn_signed_mul(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
                      int32_t bwA, int32_t bwB, int32_t bwC)
{
    auto *c = new uint64_t[dim]();
    auto *d = new uint64_t[dim];
    for (int i = 0; i < dim; i++)
    {
        c[i] = 0;
        d[i] = 0;
    }
    uint64_t *inx = new uint64_t[dim];
    uint64_t *iny = new uint64_t[dim];
    for (int i = 0; i < dim; i++)
    {
        inx[i] = (inA[i] + (1ULL << (bwA - 2))) & (bwA == 64 ? -1 : ((1ULL << bwA) - 1));
        iny[i] = (inB[i] + (1ULL << (bwB - 2))) & (bwB == 64 ? -1 : ((1ULL << bwB) - 1));
    }
    if (bwA <= bwB)
    {
        cross_term_reverse(dim, inx, iny, c, bwA, bwB, bwC);
        cross_term(dim, iny, inx, d, bwB, bwA, bwC);
    }
    else
    {
        cross_term(dim, inx, iny, c, bwA, bwB, bwC);
        cross_term_reverse(dim, iny, inx, d, bwB, bwA, bwC);
    }
    auto *wx = new uint8_t[dim];
    auto *wy = new uint8_t[dim];
    aux->wrap_computation(inx, wx, dim, bwA);
    aux->wrap_computation(iny, wy, dim, bwB);
    auto *h = new uint64_t[dim];
    auto *g = new uint64_t[dim];
    aux->multiplexer(wx, iny, h, dim, bwB, bwB);
    aux->multiplexer(wy, inx, g, dim, bwA, bwA);
    uint64_t mask = (1ULL << bwC) - 1;
    uint64_t *unsigned_mul = new uint64_t[dim];
    for (int i = 0; i < dim; i++)
    {
        unsigned_mul[i] = (inx[i] * iny[i] + c[i] + d[i] - (g[i] * (1ULL << bwB)) - (h[i] * (1ULL << bwA))) & mask;
    }
    for (int i = 0; i < dim; i++)
    {
        outC[i] = (unsigned_mul[i] - (1ULL << (bwA - 1)) * (iny[i] - (1ULL << (bwB)) * wy[i]) 
        - (1ULL << (bwB - 1)) * (inx[i] - (1ULL << (bwA)) * wx[i]) + (1ULL << (bwA + bwB - 3)) ) & mask;
    }
    delete[] c;
    delete[] d;
}

void Unsigned_Multiplication(uint64_t *x, uint64_t *y, uint64_t *output, uint8_t *wrapx, uint8_t *wrapy, int m, int n, int ell, sci::NetIO *io)
{
    uint64_t *z_share = new uint64_t[dim];
    uint64_t *x0y1 = new uint64_t[dim];
    uint64_t *x1y0 = new uint64_t[dim];
    // uint8_t *wrapx = new uint8_t[dim];
    // uint8_t *wrapy = new uint8_t[dim];
    uint64_t *g = new uint64_t[dim];
    uint64_t *h = new uint64_t[dim];
    uint64_t mask_ell = (ell == 64 ? -1 : ((1ULL << ell) - 1));
    if (party == sci::ALICE) // x0y1
    {
        std::cout << "ALICE m is " << m << " n is " << n << std::endl;
        cross_term_mul_batch(dim, party, x, x0y1, m, n, ell, otpack); // 要求m<=n
    }
    else
    {
        std::cout << "BOB m is " << m << " n is " << n << std::endl;
        cross_term_mul_batch(dim, party, y, x0y1, m, n, ell, otpack);
    }

    if (party == sci::ALICE) // x1y0
    {
        std::cout << "ALICE m is " << m << " n is " << n << std::endl;
        cross_term_mul_batch(dim, party, y, x1y0, n, m, ell, otpack); // 这里的mn和cross_term_mul_batch内相反。。。
    }
    else
    {
        std::cout << "BOB m is " << m << " n is " << n << std::endl;
        cross_term_mul_batch(dim, party, x, x1y0, n, m, ell, otpack);
    }
    aux->wrap_computation(x, wrapx, dim, m);
    aux->wrap_computation(y, wrapy, dim, n);
    aux->multiplexer(wrapx, y, g, dim, m, m);
    aux->multiplexer(wrapy, x, h, dim, n, n);
    for (int i = 0; i < dim; i++)
    {
        output[i] = (x[i] * y[i] + x0y1[i] + x1y0[i] - (1ULL << n) * g[i] - (1ULL << m) * h[i]) & mask_ell;
    }
    for (int i = 0; i < 5; i++)
    {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
        std::cout << "y[" << i << "] = " << y[i] << std::endl;
        std::cout << "x0y1[" << i << "] = " << x0y1[i] << std::endl;
        std::cout << "x1y0[" << i << "] = " << x1y0[i] << std::endl;
        std::cout << "g[" << i << "] = " << g[i] << std::endl;
        std::cout << "h[" << i << "] = " << h[i] << std::endl;
        std::cout << "output[" << i << "] = " << output[i] << std::endl;
    }
}

void signed_multiplication(uint64_t dim, uint64_t *x, uint64_t *y, uint64_t *output, int m, int n, int ell, sci::NetIO *io)
{
    uint64_t *xyoutput = new uint64_t[dim];
    uint8_t *mwrapx = new uint8_t[dim];
    uint8_t *nwrapy = new uint8_t[dim];
    uint64_t pow_m_2 = (1ull << (m - 1));
    uint64_t pow_n_2 = (1ull << (n - 1));
    uint64_t pow_m = 1ull << m;
    uint64_t pow_n = 1ull << n;
    uint64_t mask_ell = (ell == 64 ? -1 : ((1ULL << ell) - 1));

    for (int i = 0; i < dim; i++)
    {
        if (party == sci::ALICE)
        {
            x[i] = (x[i] + pow_m_2) & (m == 64 ? -1 : ((1ULL << m) - 1));
            y[i] = (y[i] + pow_n_2) & (n == 64 ? -1 : ((1ULL << n) - 1));
        }
    }

    Unsigned_Multiplication(x, y, xyoutput, mwrapx, nwrapy, m, n, ell, io);
    for (int i = 0; i < dim; i++)
    {
        if (party == sci::ALICE)
            output[i] = (xyoutput[i] - (pow_m_2 * (y[i]) - pow_n * mwrapx[i]) - (pow_n_2 * (x[i]) - pow_m * nwrapy[i])) & mask_ell;
        else
            output[i] = (xyoutput[i] - (pow_m_2 * (y[i]) - pow_n * mwrapx[i]) - (pow_n_2 * (x[i]) - pow_m * nwrapy[i]) + pow_m_2 * pow_n_2) & mask_ell;
    }
    // for(int i=0;i<5;i++)
    // {
    //     std::cout<<"m = "<<m<<", n = "<<n<<std::endl;
    //     std::cout<<"mwrapx["<<i<<"] = "<<(static_cast<int>(mwrapx[i]))<<std::endl;
    //     std::cout<<"nwrapy["<<i<<"] = "<<(static_cast<int>(nwrapy[i]))<<std::endl;
    //     std::cout<<"x["<<i<<"] = "<<x[i]<<std::endl;
    //     std::cout<<"y["<<i<<"] = "<<y[i]<<std::endl;
    //     std::cout<<"xyoutput["<<i<<"] = "<<xyoutput[i]<<std::endl;
    //     std::cout<<"output["<<i<<"] = "<<output[i]<<std::endl;
    // }
}

int init_test(uint64_t i, uint64_t j, uint64_t k, uint64_t l)
{
    mill_eq = new MillionaireWithEquality<sci::NetIO>(party, io, otpack);
    eq = new Equality<sci::NetIO>(party, io, otpack);
    uint64_t la = i;
    uint64_t lb = j;
    // la=5 f=5,la=14,f=12
    uint64_t s = k;
    uint64_t f = l;
    uint64_t *a_alice = new uint64_t[dim];
    uint64_t *b_alice = new uint64_t[dim];
    for (size_t i = 0; i < dim; i++)
    {
        a_alice[i] = 0;
        b_alice[i] = 0;
    }

    ////////////////////////////////////////////
    std::vector<std::vector<uint64_t>> data;
    if (party == sci::ALICE)
    {
        std::ifstream file("/root/Cheetah/scripts/elu_la_ld_s7.csv");
        if (!file.is_open())
        {
            std::cerr << "fail to open the file!" << std::endl;
            return 1;
        }

        std::string line;
        int target_line = 24 * (la - 2) + 2 * (lb - 1);
        int current_line = 0;

        while (std::getline(file, line))
        {
            current_line++;

            if (current_line == target_line)
            {
                std::size_t start_pos = line.find("{{");
                std::size_t end_pos = line.find("}}");

                if (start_pos != std::string::npos && end_pos != std::string::npos)
                {

                    std::string data_part = line.substr(start_pos + 2, end_pos - start_pos - 2);

                    std::stringstream ss(data_part);
                    std::string pair_str;

                    while (std::getline(ss, pair_str, '}'))
                    {

                        std::size_t open_bracket_pos = pair_str.find('{');
                        if (open_bracket_pos != std::string::npos)
                        {
                            pair_str = pair_str.substr(open_bracket_pos + 1);
                        }
                        std::stringstream pair_stream(pair_str);
                        std::string number_str;
                        std::vector<uint64_t> pair;
                        while (std::getline(pair_stream, number_str, ','))
                        {
                            if (!number_str.empty())
                            {
                                pair.push_back(static_cast<uint64_t>(std::stoull(number_str)));
                            }
                        }
                        if (pair.size() == 2)
                        {
                            data.push_back(pair);
                        }
                    }
                }
            }
        }

        file.close();
    }
    uint64_t h = f + 3;
    uint64_t alpha = 8 * pow(2, f);
    uint64_t mask_lb = (lb == 64 ? -1 : ((1ULL << lb) - 1));
    prod = new LinearOT(party, io, otpack);
    uint64_t *inA = new uint64_t[dim];
    uint64_t *inB = new uint64_t[dim];

    uint64_t *outax = new uint64_t[dim];
    // uint8_t *outb_sharp = new uint8_t[dim];
    std::cout << "\n=========STEP2 use third_interval to learn [[b]]^B===========" << std::endl;

    uint64_t STEP2_comm_start = io->counter;
    for (int i = 0; i < dim; i++)
    {
        inA[i] = (0 + i * 0) & mask_bwL;
        inB[i] = (init_input + i * step_size) & mask_bwL;
    }
    uint64_t comm_start = io->counter;
    auto time_start = std::chrono::high_resolution_clock::now();
    uint8_t *outb = new uint8_t[dim];
    uint8_t *outb_star = new uint8_t[dim];
    if (party == sci::ALICE)
    {
        second_interval(inA, outb, outb_star); // step 5SSS
    }
    else
    {
        second_interval(inB, outb, outb_star);
    }

    uint64_t STEP2_comm_end = io->counter;

    std::cout << "\n=========STEP4 use EMUX to learn [[|x|]]in L ring===========" << std::endl;
    uint64_t STEP4_comm_start = io->counter;
    aux = new AuxProtocols(party, io, otpack);
    uint64_t *EMUX_output_x = new uint64_t[dim];
    uint64_t *neg_inA = new uint64_t[dim];
    uint64_t *neg_inB = new uint64_t[dim];
    // for (int i = 0; i < dim; i++)
    // {
    //     Drelu[i] = outb[i];
    // }
    // if (party == sci::ALICE)
    // {
    //     for (int i = 0; i < dim; i++)
    //     {
    //         neg_inA[i] = ((-inA[i]+1048575) & mask_bwL); //
    //     }
    //     select_share(Drelu, inA, neg_inA, EMUX_output_x, dim, bwL); // step 10
    //     // aux->multiplexerabs(Drelu, inA, EMUX_output_x, dim, bwL, bwL);
    // }
    // else
    // {
    //     for (int i = 0; i < dim; i++)
    //     {
    //         neg_inB[i] = ((-inB[i]+1048576) & mask_bwL); //
    //     }
    //     select_share(Drelu, inB, neg_inB, EMUX_output_x, dim, bwL);
    //     // aux->multiplexerabs(Drelu, inB, EMUX_output_x, dim, bwL, bwL);
    // }
    // uint64_t *EMUX_output_x1 = new uint64_t[dim];
    // for (int i = 0; i < dim; i++)
    // {
    //     EMUX_output_x1[i] = EMUX_output_x[i];
    // }
    // uint64_t STEP4_comm_end = io->counter;
    // std::cout << "\n=========STEP7 extract the lower h bits===========" << std::endl;
    // std::cout << "inB[" << 0 << "] = " << inB[0] << std::endl;

    uint64_t *input_lower_h = new uint64_t[dim];

    uint64_t *inputx = new uint64_t[dim];
    if (party == sci::ALICE)
    {
        for (int i = 0; i < dim; i++)
        {
            inputx[i] = (inA[i] + 16384) & mask_h;
        }
    }
    else
    {
        for (int i = 0; i < dim; i++)
        {
            inputx[i] = (inB[i] + 16384) & mask_h;
        }
    }

    // uint64_t *EMUX_output_x = new uint64_t[dim];
    uint64_t *neg_x = new uint64_t[dim];
    for (int i = 0; i < dim; i++)
    {
        neg_x[i] = ((-inputx[i] + 1048576) & mask_bwL); // 取反
    }
    if (party == sci::ALICE)
    {
        for (int i = 0; i < dim; i++)
        {
            neg_x[i] = (neg_x[i] - 1) & mask_bwL;
        }
    }
    uint64_t *outtrunc = new uint64_t[dim];
    // step6 check
    std::cout << "\n=========STEP7 get mid s bit for LUT===========" << std::endl;

    trunc_oracle = new Truncation(party, io, otpack);

    uint64_t STEP5_comm_start = io->counter;

    trunc_oracle->truncate_and_reduce(dim, inputx, outtrunc, h - s, h);

    uint64_t STEP5_comm_end = io->counter;
    std::cout << "\n=========STEP6 LookUp Table   ===========" << std::endl;

    //////////////////////////////////////////

    uint64_t *outtrunc1 = new uint64_t[dim];
    uint64_t *outtrunc_a = new uint64_t[dim];
    if (party == sci::ALICE)
    {
        io->send_data(outtrunc, dim * sizeof(uint64_t));
    }
    else
    { // party == BOB
        io->recv_data(outtrunc1, dim * sizeof(uint64_t));

        for (int i = 0; i < dim; i++)
        {
            outtrunc_a[i] = (outtrunc[i] + outtrunc1[i]) & ((1ULL << s) - 1);
            // std::cout << "outtrunc_a[" << i << "] = " << outtrunc_a[i] << std::endl;
        }
    }
    uint64_t N = 1ULL << s; // LUT size

    uint64_t **spec_a = new uint64_t *[dim];

    if (party == sci::ALICE)
        for (int i = 0; i < dim; i++)
        {
            spec_a[i] = new uint64_t[N];
            for (int j = 0; j < N; j++)
            {
                spec_a[i][j] = data[j][0];
                // std::cout << "i = " << i << ", j = " << j << ", data = " << data[j][i] << std::endl;
            }
        }

    uint64_t STEP6_comm_start = io->counter;
    std::cout << "lookup table" << std::endl;
    uint64_t *a_bob = new uint64_t[dim];
    if (party == sci::ALICE)
    {
        aux->lookup_table<uint64_t>(spec_a, nullptr, nullptr, dim, s, la); // step 12 lut
    }
    else
    {                                                                        // party == BOB
        aux->lookup_table<uint64_t>(nullptr, outtrunc_a, a_bob, dim, s, la); //
    }
    std::cout << "lookup table done" << std::endl;
    // if (party != sci::ALICE)
    //     for (int i = 0; i < dim; i++)
    //     {
    //         // std::cout << "a_bob[" << i << "] = " << a_bob[i] << std::endl;
    //     }
    /////选择截距
    uint64_t **spec_b = new uint64_t *[dim];
    uint64_t *b_bob = new uint64_t[dim];
    if (party == sci::ALICE)
        for (int i = 0; i < dim; i++)
        {
            spec_b[i] = new uint64_t[N];
            for (int j = 0; j < N; j++)
            {
                spec_b[i][j] = data[j][1];
            }
        }
    if (party == sci::ALICE)
    {
        aux->lookup_table<uint64_t>(spec_b, nullptr, nullptr, dim, s, lb); // step 12 lut
    }
    else
    {                                                                        // party == BOB
        aux->lookup_table<uint64_t>(nullptr, outtrunc_a, b_bob, dim, s, lb); //
    }
    // if (party != sci::ALICE)
    // std::cout << "b_bob[" << 0 << "] = " << b_bob[0] << std::endl;

    uint64_t STEP6_comm_end = io->counter;
    // cout << "LUT Bytes Sent: " << (comm_end_lut - comm_start_lut) << "bytes" << endl;

    ext = new XTProtocol(party, io, otpack);

    std::cout << "\n=========STEP7 multiplication to get a|x| l+la ===========" << std::endl;
    uint64_t STEP7_comm_start = io->counter;

    //   test_matrix_multiplication(inA, inB, outC, false);
    // test_matrix_multiplication(inA, inB, outC, true);
    if (correct == 1)
    {
        if (party == sci::ALICE)
        {
            // std::cout << "inA_h[" << 0 << "] = " << inA_h[0] << std::endl;
            // std::cout << "a_alice[" << 1 << "] = " << a_alice[1] << std::endl;
            // std::cout << "EMUX_output_x[" << 1 << "] = " << EMUX_output_x[1] << std::endl;
            // prod->hadamard_product(dim, a_alice, EMUX_output_x, outax, la, bwL, la + bwL, true, true, mode, msb1, msb2); // step 13 mul
            // signed_multiplication(dim, a_alice, EMUX_output_x, outax, la, bwL, la + bwL, io);

            sirnn_signed_mul(dim, a_alice, inA, outax, la, bwL, la + bwL);

            // io->send_data(a_alice, dim * sizeof(uint64_t));
            // io->send_data(EMUX_output_x, dim * sizeof(uint64_t));
            // for (int i = 0; i < dim; i++)
            // {
            //     outax[i] = 0;
            // }

            // cross_term_mul(dim,party, a_alice, EMUX_output_x, outax, la, bwL,la + bwL, otpack);
            // std::cout << "outax[" << 0 << "] = " << outax[0] << std::endl;
        }
        else
        {
            // std::cout << "inB_h[" << 0 << "] = " << inB_h[0] << std::endl;
            // std::cout << "a_bob[" << 1 << "] = " << a_bob[1] << std::endl;
            // std::cout << "EMUX_output_x[" << 1 << "] = " << EMUX_output_x[1] << std::endl;
            // prod->hadamard_product(dim, a_bob, EMUX_output_x, outax, la, bwL, la + bwL, true, true, mode, msb1, msb2);
            // signed_multiplication(dim, a_bob, EMUX_output_x, outax, la, bwL, la + bwL, io);
            sirnn_signed_mul(dim, a_bob, inB, outax, la, bwL, la + bwL);

            // uint64_t *recv_a_alice = new uint64_t[dim];
            // io->recv_data(recv_a_alice, dim * sizeof(uint64_t));
            // uint64_t *recv_EMUX_output_x = new uint64_t[dim];
            // io->recv_data(recv_EMUX_output_x, dim * sizeof(uint64_t));
            // uint64_t *a = new uint64_t[dim];
            // uint64_t *b = new uint64_t[dim];
            // uint64_t mask_la = (la == 64 ? -1 : ((1ULL << la) - 1));
            // uint64_t mask_labwL = (la + bwL == 64 ? -1 : ((1ULL << (la + bwL)) - 1));
            // for (int i = 0; i < dim; i++)
            // {
            //     a[i] = (recv_a_alice[i] + a_bob[i]) & mask_la;
            //     b[i] = (recv_EMUX_output_x[i] + EMUX_output_x[i]) & mask_bwL;
            //     outax[i] =( a[i] * b[i]) & mask_labwL;
            // }

            // cross_term_mul(dim,party, a_bob, EMUX_output_x, outax, la, bwL,la + bwL, otpack);
            // std::cout << "outax[" << 0 << "] = " << outax[0] << std::endl;
        }
    }
    uint64_t STEP7_comm_end = io->counter;

    std::cout << "\n=========STEP8 ax truncate from l+la to l+1  ===========" << std::endl; //
                                                                                            ////////////////////////////////////////////////////////
    // uint8_t *msb_zero = new uint8_t[dim];
    // for (int i = 0; i < dim; i++)
    // {
    //     msb_zero[i] = 0;
    // }
    uint64_t *mid_ax = new uint64_t[dim];
    uint64_t STEP8_comm_start = io->counter;

    trunc_oracle->truncate_and_reduce(dim, outax, mid_ax, la - 1, bwL + la); // step 10 tr

    uint64_t STEP8_comm_end = io->counter;
    for (int i = 0; i < dim; i++)
    {
        outax[i] = mid_ax[i];
    }
    std::cout << "\n=========STEP11 d SExt with MSB from f+1 to l   ===========" << std::endl;
    uint64_t *b_SExt = new uint64_t[dim];

    uint8_t *msb_b_extend = new uint8_t[dim];

    uint64_t s_extend_comm_start = io->counter;
    if (party == sci::ALICE)
    {
        for (int i = 0; i < dim; i++)
        {
            msb_b_extend[i] = 1;
            b_alice[i] = (b_alice[i] + 10) & mask_lb;
        }
        ext->s_extend_msb(dim, b_alice, b_SExt, lb, bwL, msb_b_extend); // step 13 s_extend
    }
    else
    {
        for (int i = 0; i < dim; i++)
        {
            msb_b_extend[i] = 1;
            // std::cout << "b_bob[" << i << "] = " << b_bob[i] << std::endl;
            b_bob[i] = (b_bob[i] - 10) & mask_lb;
        }
        ext->s_extend_msb(dim, b_bob, b_SExt, lb, bwL, msb_b_extend);
    }
    uint64_t s_extend_comm_end = io->counter;

    // std::cout << "b_SExt[" << 0 << "] = " << b_SExt[0] << std::endl;

    std::cout << "\n=========STEP12 Caculate z=ax+b   ===========" << std::endl;
    uint64_t *z = new uint64_t[dim];

    for (int i = 0; i < dim; i++)
        z[i] = ((outax[i] + b_SExt[i] * static_cast<uint64_t>(std::pow(2, f - lb + 1))) & mask_bwL); // step 18 add

    std::cout << "\n=========STEP15 get x_half ===========" << std::endl;
    uint64_t *y = new uint64_t[dim];
    uint64_t *non_negative1_part = new uint64_t[dim];

    uint64_t *neg1 = new uint64_t[dim];

    for (int i = 0; i < dim; i++)
    {
        neg1[i] = 1046528;
    }

    int64_t STEP15_comm_start = io->counter;
    if (party == sci::ALICE)
    {
        select_share(outb, inA, z, non_negative1_part, dim, bwL); // step 16
    }
    else
    {
        select_share(outb, inB, z, non_negative1_part, dim, bwL);
    }

    uint8_t *choose_negative_part = new uint8_t[dim];

    for (int i = 0; i < dim; i++)
    {
        choose_negative_part[i] = (outb[i] ^ outb_star[i]) & 1;
    }

    int64_t STEP15_comm_end = io->counter;
    select_share(choose_negative_part, non_negative1_part, neg1, y, dim, bwL);
    // std::cout << "xhalf[" << 0 << "] = " << xhalf[0] << std::endl;
    // std::cout << "abs_xhalf[" << 0 << "] = " << abs_xhalf[0] << std::endl;

    std::cout << "\n=========STEP16 get delta = z-x_half ===========" << std::endl;

    std::cout << "\n=========STEP17 |g|=delta_ + x_half ===========" << std::endl;
    // uint64_t *delta_ = new uint64_t[dim];
    // for (int i = 0; i < dim; i++)
    // {
    //     delta_[i] = 0;
    // }
    // aux->multiplexer(Drelu_, delta, delta_, dim, bwL, bwL);
    // std::cout << "MUX_output_u[" << 0 << "] =" << MUX_output_u[0] << std::endl;
    uint64_t *MUX_output_g = new uint64_t[dim];
    // int64_t STEP21_comm_start = io->counter;

    // for (int i = 0; i < dim; i++)
    // {
    //     Drelu_[i] = (outb_star[i] + outb_sharp[i]) & 1;
    //     if (party == sci::ALICE)
    //     {
    //         Drelu_[i] = Drelu_[i] ^ 1;
    //     }
    // }

    // select_share(Drelu_, abs_xhalf, z, MUX_output_g, dim, bwL); // step 20 ss

    // int64_t STEP21_comm_end = io->counter;
    // for (int i = 0; i < dim; i++)
    // {
    //     std::cout << "outb[" << i << "] = " << static_cast<int>(outb[i]) << std::endl;
    //     std::cout << "outb_star[" << i << "] = " << static_cast<int>(outb_star[i]) << std::endl;
    //     std::cout << "outb_sharp[" << i << "] = " << static_cast<int>(outb_sharp[i]) << std::endl;
    //     std::cout << "Drelu_[" << i << "] = " << static_cast<int>(Drelu_[i]) << std::endl;
    //     // std::cout << "delta_[" << i << "] = " << delta_[i] << std::endl;
    //     std::cout << "z[" << i << "] = " << z[i] << std::endl;
    //     // MUX_output_g[i] = (delta_[i] + z[i]) & mask_bwL;abs_xhalf
    //     std::cout << "abs_xhalf[" << i << "] = " << abs_xhalf[i] << std::endl;
    //     std::cout << "neg_abs_xhalf[" << i << "] = " << neg_abs_xhalf[i] << std::endl;
    //     std::cout << "xhalf[" << i << "] = " << xhalf[i] << std::endl;
    //     std::cout << "MUX_output_g[" << i << "] = " << MUX_output_g[i] << std::endl;
    // }
    uint64_t comm_end = io->counter;
    auto time_end = chrono::high_resolution_clock::now();
    ////////////////////////////////////verfication
    // for (int i = 0; i < dim; i++)
    // {
    //     std::cout << "outb[" << i << "] = " << static_cast<int>(outb[i]) << std::endl;
    //     std::cout << "outb_star[" << i << "] = " << static_cast<int>(outb_star[i]) << std::endl;
    //     std::cout << "outb_sharp[" << i << "] = " << static_cast<int>(outb_sharp[i]) << std::endl;
    //     //outb_star[i] ^outb_sharp[i]
    //     std::cout << "total outb[" << i << "] = " << static_cast<int>((outb_star[i] + outb_sharp[i]) & 1) << std::endl;
    //     std::cout << "Drelu_[" << i << "] = " << static_cast<int>(Drelu_[i]) << std::endl;
    //     std::cout << "Drelu[" << i << "] = " << static_cast<int>(Drelu[i]) << std::endl;

    // }

    std::cout << "\n=========END verification ===========" << std::endl;
    if (party == sci::ALICE)
    {
        io->send_data(y, dim * sizeof(uint64_t));
        uint64_t *comm = new uint64_t[1];
        comm[0] = (comm_end - comm_start) / dim * 8;
        io->send_data(comm, 1 * sizeof(uint64_t));
    }
    else
    {
        uint64_t *recv_y = new uint64_t[dim];
        io->recv_data(recv_y, dim * sizeof(uint64_t));
        // double recv_Total_MSBytes_ALICE;
        // iopack->io->recv_data(&recv_Total_MSBytes_ALICE, sizeof(double));
        double *ULPs = new double[dim];
        double f_pow = pow(2, f);
        int s_y = 12;
        for (int i = 0; i < dim; i++)
        {
            // std::cout << "dim [" << i << "]total y = y0 + y1 =  " << ((y[i] + recv_y[i]) & mask_bwL) << ", real num: " << (double)decode_ring((y[i] + recv_y[i]) & mask_bwL, bwL) / f_pow << std::endl;
            // std::cout << "The result " << inA[i] + inB[i] << " should be calculate_ELU = " << calculate_ELU(inA[i] + inB[i], f) << std::endl;
            // ULPs[i] = abs((((double)decode_ring((y[i] + recv_y[i]) & mask_bwL, bwL) / f_pow) - calculate_ELU(inA[i] + inB[i], f)) / 0.000244140625);
            ULPs[i] = computeULPErr(((double)decode_ring((y[i] + recv_y[i]) & mask_bwL, bwL) / f_pow), calculate_ELU(inA[i] + inB[i], 12), s_y);
            // std::cout << "The ULP is = " << ULPs[i] << std::endl;
        }
        double sum = 0.0;
        for (size_t i = 0; i < dim; ++i)
        {
            sum += (ULPs[i]);
            // std::cout << "ULPs[" << i << "] = " << ULPs[i] << std::endl;
        }
        double average = 0.0;
        double max_val = 0.0;
        double min_val = 0.0;
        average = sum / static_cast<double>(dim);
        // std::cout << "sum: " << sum << std::endl;
        // std::cout << "static_cast<double>(dim): " << static_cast<double>(dim) << std::endl;
        max_val = *std::max_element(ULPs, ULPs + dim); //
        min_val = *std::min_element(ULPs, ULPs + dim); //
        std::cout << "average: " << average << std::endl;
        std::cout << "max_val: " << max_val << std::endl;
        std::cout << "min_val: " << min_val << std::endl;
        uint64_t *alice_comm = new uint64_t[1];
        io->recv_data(alice_comm, 1 * sizeof(uint64_t));

        uint64_t bob_comm = (comm_end - comm_start) / dim * 8;
        uint64_t total_comm = bob_comm + alice_comm[0];
        std::ofstream file("/home/ubuntu/EzPC/elu_output_data.csv", std::ios_base::app);

        if (!file.is_open())
        {
            std::cerr << "Error: Could not open file for writing." << std::endl;
            // return 1;
        }

        auto total_time = chrono::duration_cast<chrono::milliseconds>(time_end - time_start).count();
        file << la << "," << lb << ",  " << average << ",   " << max_val << "  , " << total_comm << ",    " << total_time << "\n";
    }

    ///////////

    cout << "STEP2 Third_interval Bytes Sent: " << (STEP2_comm_end - STEP2_comm_start) / dim * 8 << " bits" << endl;
    // cout << "STEP4 Select_share Bytes Sent: " << (STEP4_comm_end - STEP4_comm_start) / dim * 8 << " bits" << endl;
    cout << "STEP5 TR Bytes Sent: " << (STEP5_comm_end - STEP5_comm_start) / dim * 8 << " bits" << endl;
    cout << "STEP6 LUT*2 Bytes Sent: " << (STEP6_comm_end - STEP6_comm_start) / dim * 8 << " bits" << endl;
    cout << "STEP7 hadamard_product Bytes Sent: " << (STEP7_comm_end - STEP7_comm_start) / dim * 8 << " bits" << endl;
    cout << "STEP8 truncate_and_reduce Bytes Sent: " << (STEP8_comm_end - STEP8_comm_start) / dim * 8 << " bits" << endl;
    std::cout << "s_extend_comm: " << (s_extend_comm_end - s_extend_comm_start) / dim * 8 << std::endl;
    // cout << "STEP14 DRELUsec Bytes Sent: " << (STEP14_comm_end - STEP14_comm_start) / dim * 8 << " bits" << endl;
    cout << "STEP15 clear_MSB_to_Wrap_bitMul and one trunc Bytes Sent: " << (STEP15_comm_end - STEP15_comm_start) / dim * 8 << " bits" << endl;
    // cout << "STEP21 select_share Bytes Sent: " << (STEP21_comm_end - STEP21_comm_start) / dim * 8 << " bits" << endl;
    // cout << "STEP3 Bytes Sent: " << (comm_end - comm_start) << "bytes" << endl;
    // cout << "STEP3 Bytes Sent: " << (comm_end - comm_start) << "bytes" << endl;
    // cout << "STEP3 Bytes Sent: " << (comm_end - comm_start) << "bytes" << endl;
    // cout << "STEP3 Bytes Sent: " << (comm_end - comm_start) << "bytes" << endl;
    // cout << "STEP3 Bytes Sent: " << (comm_end - comm_start) << "bytes" << endl;
    // uint64_t comm_end = io->counter;
    cout << "Total Bytes Sent: " << (comm_end - comm_start) / dim * 8 << " bits" << endl;

    cout << "Total time: "
         << chrono::duration_cast<chrono::milliseconds>(time_end - time_start).count()
         << " ms" << endl;
    delete[] a_alice;
    delete[] b_alice;
    delete[] inA;
    delete[] inB;
    delete[] outax;
    delete[] outb;
    delete[] outb_star;
    delete[] EMUX_output_x;
    delete[] neg_inA;
    delete[] neg_inB;
    delete[] outtrunc;
    delete[] outtrunc1;
    delete[] outtrunc_a;
    delete[] a_bob;
    delete[] b_bob;
    delete[] mid_ax;
    delete[] b_SExt;
    delete[] msb_b_extend;
    delete[] z;
    delete[] MUX_output_g;
    delete[] y;

    if (party == sci::ALICE)
    {
        for (int i = 0; i < dim; i++)
        {
            delete[] spec_a[i];
            delete[] spec_b[i];
        }
    }
    delete[] spec_a;
    delete[] spec_b;

    delete prod;
    delete aux;
    delete trunc_oracle;
    delete ext;
    return 0;
}

int main(int argc, char **argv)
{
    std::cout << "start" << std::endl;
    ArgMapping amap;
    std::cout << "start1" << std::endl;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("m", precomputed_MSBs, "MSB_to_Wrap Optimization?");
    amap.arg("a", ::accumulate, "Accumulate?");
    amap.arg("dim", dim, "Dimension parameter for accumulation");
    amap.arg("init_input", init_input, "init_input for accumulation");
    amap.arg("step_size", step_size, "step_size for accumulation");
    amap.arg("acc", acc, "acc=0 low, acc=1 general (default), acc =2 high");
    amap.arg("correct", correct, "correct=1 or communication=2");

    amap.parse(argc, argv);
    std::cout << "Parsed dimension (dim) = " << dim << std::endl;

    io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    otpack = new sci::OTPack<sci::NetIO>(io, party);
    XTProtocol *ext = new XTProtocol(party, io, otpack);
    aux = new AuxProtocols(party, io, otpack);
    Truncation *trunc_oracle = new Truncation(party, io, otpack);
    MillionaireWithEquality<sci::NetIO> *mill_eq = new MillionaireWithEquality<sci::NetIO>(party, io, otpack);
    LinearOT *prod = new LinearOT(party, io, otpack);
    Equality<sci::NetIO> *eq = new Equality<sci::NetIO>(party, io, otpack);

    std::vector<std::pair<uint64_t, uint64_t>> la_lb_pairs = {
        {8,13},{6,13},{6,12},{5,12},{5,11},{5,10},{4,11},{4,9}};
    // {8, 12}, {7, 12}, {6, 12}, {6, 11}, {5, 12}, {5, 10}, {4, 12}, {4, 10}};
    // {4, 10}};

    for (const auto &pair : la_lb_pairs)
    {
        uint64_t la = pair.first;
        uint64_t lb = pair.second;

        for (uint64_t s = 7; s < 8; s++)
        {
            for (uint64_t k = 12; k < 13; k++)
            {
                // if ((la <= k) & (lb <= k))
                // {
                // std::cout << "la = " << la << ", lb = " << lb << ", s = " << s << ", k = " << k << std::endl;
                for (int i = 0; i < 1; i++)
                {
                    init_test(la, lb, s, k);
                }
                std::cout << "la = " << la << ", lb = " << lb << ", s = " << s << ", k = " << k << std::endl;
                // }
            }
        }
    }
    delete prod;
}
