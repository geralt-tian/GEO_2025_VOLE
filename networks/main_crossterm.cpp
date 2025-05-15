// Algorithm CrossTermMultiplication(x, y):
//     # 输入: P0 持有 x, P1 持有 y
//     # 输出: P0, P1 各自获得 z = x * y 的 share

//     1. P0: 将 x 解析为 m 位二进制串 (x_{m-1}, ..., x_0), 其中 x_i ∈ {0,1}
//     2. 对于每个 i in {0, ..., m-1}:
//         a. P0 和 P1 共同调用 (2, 1)-COT_{l-i} 协议:
//             - P0 作为 sender，输入 x_i
//             - P1 作为 receiver，输入 y
//             - 获得 t_i^{l-i}（长度为 l-i 的 share）
//     3. 对 b ∈ {0,1}（即每个参与方）:
//         - P_b 计算:
//             <z>_b = ∑_{i=0}^{m-1} 2^i * <t_i>_b^{l-i}
//             # 其中 <t_i>_b^{l-i} 是第 i 个 COT 返回给 P_b 的 share，左移 2^i
//     4. 返回 <z>_b

// party: sci::ALICE (1) or sci::BOB (2)
// x: P0输入（仅ALICE用）
// y: P1输入（仅BOB用）
// m: 比特宽度
// n: y的比特宽
// otpack: OT包
// 返回：本方的 share

#include <iostream>
#include <string>
#include <vector>
#include "utils/emp-tool.h"
#include "OT/ot-utils.h"
#include "OT/emp-ot.h"
#include <chrono>

using namespace std;
using namespace sci;

void cross_term_mul_batch(
    uint64_t batch_size,
    int party,
    const uint64_t *x, // ALICE: x[batch_size]，BOB: nullptr
    // const uint64_t *y, // BOB: y[batch_size]，ALICE: nullptr
    uint64_t *z_share, // 输出

    int m, int n, int ell,
    sci::OTPack<sci::NetIO> *otpack)
{
    // std::cout << "cross_term_mul_batch" << std::endl;
    sci::NetIO *io = otpack->io;
    uint64_t mask_m = (m == 64 ? -1 : ((1ULL << m) - 1));
    uint64_t mask_n = (n == 64 ? -1 : ((1ULL << n) - 1));
    
    for (int i = 0; i < m; ++i)
    {
        uint64_t *data0 = new uint64_t[batch_size];
        uint64_t *corr = new uint64_t[batch_size];
        uint64_t *recv_data = new uint64_t[batch_size];
        uint64_t mask_ell_i = (ell - i == 64 ? -1 : ((1ULL << (ell - i)) - 1));
        bool *choice = new bool[batch_size];
        std::cout << "o->counter;" << std::endl;
        size_t comm_before = io->counter;
        std::cout << "comm_before" << std::endl;
        if (party == sci::ALICE)
        {
            std::cout << "ALICE" << std::endl;
            for (int b = 0; b < batch_size; b++)
            {
                corr[b] = x[b]; // 只 ALICE 用 y
            }
            std::cout << "send_cot" << std::endl;
            otpack->iknp_straight->send_cot(data0, corr, batch_size, ell - i);
            std::cout << "send_cot" << std::endl;
            for (int b = 0; b < batch_size; b++)
            {
                data0[b] = (-data0[b]) & mask_ell_i;
                z_share[b] += (data0[b] << i);
                std::cout << "z_share[" << i << "] = " << z_share[b] << std::endl;
            }
        }
        else
        {
            for (int b = 0; b < batch_size; b++)
            {
                choice[b] = (x[b] >> i) & 1; // 只 BOB 用 x
            }
            otpack->iknp_straight->recv_cot(recv_data, choice, batch_size, ell - i);
            for (int b = 0; b < batch_size; b++)
            {
                z_share[b] += (recv_data[b] << i);
            }
        }
        io->flush();
    }
}
// }

// CrossTerm 乘法协议实现（严格2,1-COT版本）
uint64_t cross_term_mul(int party, uint64_t x, uint64_t y, int m, int n, sci::OTPack<sci::NetIO> *otpack)
{
    int ell = m + n;
    uint64_t z_share = 0;
    sci::NetIO *io = otpack->io;
    uint64_t mask_m = (1ULL << m) - 1;
    uint64_t mask_n = (1ULL << n) - 1;
    for (int i = 0; i < m; i++)
    {
        uint64_t data0[1], corr[1], recv_data[1];
        bool choice[1];
        // size_t comm_before = io->counter;
        if (party == sci::ALICE)
        {
            corr[0] = y;
            // std::cout << "[COT] i=" << i << ", before send_cot, comm=" << comm_before << std::endl;
            otpack->iknp_straight->send_cot(data0, corr, 1, ell - i);
            // size_t comm_after = io->counter;
            // std::cout << "[COT] i=" << i << ", after send_cot, comm=" << comm_after
            //           << ", delta=" << (comm_after - comm_before) << " bytes" << std::endl;
            data0[0] = (-data0[0]) ;
            z_share += (data0[0] << (i))% (1ULL << (ell - i));
        }
        else
        {
            choice[0] = (x >> i) & 1;
            // std::cout << "[COT] i=" << i << ", before recv_cot, comm=" << comm_before << std::endl;
            otpack->iknp_straight->recv_cot(recv_data, choice, 1, ell - i);
            // size_t comm_after = io->counter;
            // std::cout << "[COT] i=" << i << ", after recv_cot, comm=" << comm_after
            //           << ", delta=" << (comm_after - comm_before) << " bytes" << std::endl;
            z_share += (recv_data[0] << (i));
        }
        io->flush();
    }
    return z_share;
}

int main(int argc, char **argv)
{
    int party = 0, port = 8000;
    string address = "127.0.0.1";
    int n = 10; // y的比特宽
    int m = 8; // x的比特宽

    // 参数解析
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("m", m, "Bit width of x");
    amap.arg("n", n, "Bit width of y");
    amap.parse(argc, argv);

    // 网络与OT初始化
    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);

    // 统计开始
    auto start = std::chrono::high_resolution_clock::now();
    size_t comm_start = io->counter;

    // 输入生成
    uint64_t x = 122, y = 77;
    if (party == sci::ALICE)
    {
        x = 122 % (1ULL << m); // 你可以自定义输入
        // std::cout << "ALICE input x = " << x << std::endl;
    }
    else
    {
        y = 77 % (1ULL << n); // 你可以自定义输入
        // std::cout << "BOB input y = " << y << std::endl;
    }
    x = 122, y = 77;
    // 执行 CrossTerm 乘法

    // for (int i = 0; i < 65536; i++)
    // {

    if (party == sci::ALICE)
    {
        uint64_t z_share = cross_term_mul(party, 0, y, m, n, otpack);
        std::cout << "z_share = " << z_share << std::endl;
    }
    else
    {
        uint64_t z_share = cross_term_mul(party, x, 0, m, n, otpack);
        std::cout << "z_share = " << z_share << std::endl;
    }
    // }

    // 输出本方 share
    // if (party == sci::ALICE) std::cout << "ALICE share: " << z_share << std::endl;
    // if (party == sci::BOB) std::cout << "BOB share: " << z_share << std::endl;
    const int batch_size = 1024;
    uint64_t *x1 = new uint64_t[batch_size];
    uint64_t *y1 = new uint64_t[batch_size];
    uint64_t *z_share1 = new uint64_t[batch_size];
    uint64_t mask_m = (m == 64 ? -1 : ((1ULL << m) - 1));
    uint64_t mask_n = (n == 64 ? -1 : ((1ULL << n) - 1));
    for (int i = 0; i < batch_size; ++i)
    {
        x1[i] = (122 + i) & mask_m;
        y1[i] = (77 + i) & mask_n;
        z_share1[i] = 0;
    }
    if (party == sci::ALICE)
    {
        cross_term_mul_batch(batch_size, party,  x1, z_share1, m, n, m + n, otpack);
    }
    else
    {
        cross_term_mul_batch(batch_size, party, y1,  z_share1, m, n, m + n, otpack);
    }
    for (int i = 0; i < batch_size; ++i)
    {
        std::cout << "x1[" << i << "] = " << x1[i] << std::endl;
        std::cout << "y1[" << i << "] = " << y1[i] << std::endl;
        std::cout << "z_share1[" << i << "] = " << z_share1[i] << std::endl;
    }
    // 统计结束
    auto end = std::chrono::high_resolution_clock::now();
    size_t comm_end = io->counter;
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double comm_kb = (comm_end - comm_start) / 1024.0;
    std::cout << "Total communication: " << comm_kb << " KB" << std::endl;
    std::cout << "Elapsed time: " << elapsed_ms << " ms" << std::endl;

    delete otpack;
    delete io;
    return 0;
}
