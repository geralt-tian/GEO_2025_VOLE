/*
 * Test file for exp_nagx function - SEAF-VOLE-AF Port
 * This file tests the exp_nagx function which computes exp(-x) where x > 0
 * Ported from EzPC/SCI to use VOLE (Silent OT) instead of IKNP-OT
 * Including ULP error analysis, communication measurement, and timing
 */

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <random>
#include <cstring>

#include "Math/math-functions.h"
#include "BuildingBlocks/geometric_perspective_protocols.h"
#include "utils/emp-tool.h"

using namespace sci;
using namespace std;

#define MAX_THREADS 4

int party, port = 32000;
int num_threads = 1;  // Start with single thread for simplicity
string address = "127.0.0.1";

// Test parameters for exp_nagx (matching EzPC test)
int dim = 1024;      // Default: 1024 elements
int32_t in_bw = 37;  // Input bit width
int32_t in_f = 12;   // Input fractional bits

int32_t localexp_f = 10;      // Local exp fractional bits
int32_t localexp_bw = localexp_f + 13;  // Local exp bit width

int32_t locallut_f = 32;      // Local LUT fractional bits
int32_t locallut_bw = locallut_f + 2;   // Local LUT bit width

uint64_t mask_in = (in_bw == 64 ? -1 : ((1ULL << in_bw) - 1));

sci::NetIO *ioArr[MAX_THREADS];
sci::OTPack<sci::NetIO> *otpackArr[MAX_THREADS];
GeometricPerspectiveProtocols *gpArr[MAX_THREADS];

uint64_t computeULPErr(double calc, double actual, int SCALE) {
    int64_t calc_fixed = (double(calc) * (1ULL << SCALE));
    int64_t actual_fixed = (double(actual) * (1ULL << SCALE));
    uint64_t ulp_err = (calc_fixed - actual_fixed) > 0
                           ? (calc_fixed - actual_fixed)
                           : (actual_fixed - calc_fixed);
    return ulp_err;
}

void exp_thread(int tid, uint64_t *x, uint64_t *y, int num_exp) {
    cout << "[Thread " << tid << "] Calling exp_softmaxx..." << endl;
    cout << "[Thread " << tid << "] Params: dim=" << num_exp << ", in_bw=" << in_bw << ", in_f=" << in_f << endl;

    // Call exp_softmaxx from GeometricPerspectiveProtocols (exp_nagx implementation)
    gpArr[tid]->exp_softmaxx(num_exp, x, y, in_bw, in_f,
                             localexp_bw, localexp_f, locallut_bw, locallut_f);

    cout << "[Thread " << tid << "] exp_softmaxx completed" << endl;
}

int main(int argc, char **argv) {
    /************* Argument Parsing  ************/
    /********************************************/
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("d", dim, "Number of elements to test");
    amap.arg("nt", num_threads, "Number of threads");
    amap.arg("ip", address, "IP Address of server (ALICE)");

    amap.parse(argc, argv);

    assert(num_threads <= MAX_THREADS);

    cout << "=============================================" << endl;
    cout << "exp_nagx Test Suite - SEAF-VOLE-AF Port" << endl;
    cout << "=============================================" << endl;
    cout << "Testing exp(-x) where x > 0" << endl;
    cout << "Protocol: VOLE (Silent OT)" << endl;
    cout << "Party: " << (party == sci::ALICE ? "ALICE" : "BOB") << endl;
    cout << "Dimensions: " << dim << endl;
    cout << "Input format: " << in_bw << "." << in_f << endl;
    cout << "Threads: " << num_threads << endl;
    cout << "=============================================" << endl;

    /********** Setup IO and Base OTs ***********/
    /********************************************/
    for (int i = 0; i < num_threads; i++) {
        ioArr[i] = new sci::NetIO(party == 1 ? nullptr : address.c_str(), port + i);
        if (i & 1) {
            otpackArr[i] = new sci::OTPack<sci::NetIO>(ioArr[i], 3 - party);
            gpArr[i] = new GeometricPerspectiveProtocols(3 - party, ioArr[i], otpackArr[i]);
        } else {
            otpackArr[i] = new sci::OTPack<sci::NetIO>(ioArr[i], party);
            gpArr[i] = new GeometricPerspectiveProtocols(party, ioArr[i], otpackArr[i]);
        }
    }
    cout << "All Base OTs Done" << endl;

    /************ Generate Test Data ************/
    /********************************************/
    PRG128 prg;

    uint64_t *x = new uint64_t[dim];
    uint64_t *y = new uint64_t[dim];

    prg.random_data(x, dim * sizeof(uint64_t));

    cout << "Generating test data for " << dim << " elements..." << endl;

    if (party == sci::ALICE) {
        ioArr[0]->send_data(x, dim * sizeof(uint64_t));
    } else {
        uint64_t *x0 = new uint64_t[dim];
        ioArr[0]->recv_data(x0, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            // x is always negative
            x[i] = ((1ULL << (in_bw - 1)) + (x[i] & (mask_in >> 1))) - x0[i];
        }
        delete[] x0;
    }
    for (int i = 0; i < dim; i++) {
        x[i] &= mask_in;
    }

    /************** Fork Threads ****************/
    /********************************************/
    uint64_t total_comm = 0;
    uint64_t thread_comm[num_threads];
    for (int i = 0; i < num_threads; i++) {
        thread_comm[i] = ioArr[i]->counter;
    }

    cout << "Starting exp_nagx computation..." << endl;
    auto start_time = chrono::high_resolution_clock::now();

    std::thread exp_threads[num_threads];
    int chunk_size = dim / num_threads;
    for (int i = 0; i < num_threads; ++i) {
        int offset = i * chunk_size;
        int lnum_exp;
        if (i == (num_threads - 1)) {
            lnum_exp = dim - offset;
        } else {
            lnum_exp = chunk_size;
        }
        exp_threads[i] = std::thread(exp_thread, i, x + offset, y + offset, lnum_exp);
    }
    for (int i = 0; i < num_threads; ++i) {
        exp_threads[i].join();
    }

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);

    for (int i = 0; i < num_threads; i++) {
        thread_comm[i] = ioArr[i]->counter - thread_comm[i];
        total_comm += thread_comm[i];
    }

    cout << "Computation completed in " << duration.count() << " ms" << endl;

    /************** Verification ****************/
    /********************************************/
    if (party == sci::ALICE) {
        ioArr[0]->send_data(x, dim * sizeof(uint64_t));
        ioArr[0]->send_data(y, dim * sizeof(uint64_t));
    } else { // party == BOB
        uint64_t *x0 = new uint64_t[dim];
        uint64_t *y0 = new uint64_t[dim];
        ioArr[0]->recv_data(x0, dim * sizeof(uint64_t));
        ioArr[0]->recv_data(y0, dim * sizeof(uint64_t));

        uint64_t total_err = 0;
        uint64_t max_ULP_err = 0;
        for (int i = 0; i < dim; i++) {
            double dbl_x = (signed_val(x0[i] + x[i], in_bw)) / double(1LL << in_f);
            double dbl_y = (signed_val(y0[i] + y[i], in_bw)) / double(1ULL << in_f);
            double exp_x = exp(dbl_x);
            uint64_t err = computeULPErr(dbl_y, exp_x, in_f);
            total_err += err;
            max_ULP_err = std::max(max_ULP_err, err);
        }

        cerr << "Average ULP error: " << total_err / dim << endl;
        cerr << "Max ULP error: " << max_ULP_err << endl;
        cerr << "Number of tests: " << dim << endl;

        delete[] x0;
        delete[] y0;
    }
    cout << "Number of Exp/s:\t" << (double(dim) / duration.count()) * 1000.0 << std::endl;
    cout << "Exp Time\t" << duration.count() << " ms" << endl;
    cout << "Exp Bytes Sent\t" << total_comm << " bytes" << endl;

    /******************* Cleanup ****************/
    /********************************************/
    delete[] x;
    delete[] y;

    for (int i = 0; i < num_threads; i++) {
        delete gpArr[i];
        delete otpackArr[i];
        delete ioArr[i];
    }

    return 0;
}
