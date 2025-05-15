#include <iostream>
#include <string>
#include "LinearOT/linear-ot.h"
#include "utils/emp-tool.h"
#include "SCI/src/BuildingBlocks/aux-protocols.h"
#include <chrono>

using namespace std;
using namespace sci;

int party = 0, port = 8000;
    string address = "127.0.0.1";
int dim = 1024;
    int bwA = 8;
    int bwB = 8;
    int bwC = 16;

int main(int argc, char **argv) {
    // 解析命令行参数
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Dimension");
    amap.arg("bwA", bwA, "Bitwidth of A");
    amap.arg("bwB", bwB, "Bitwidth of B");
    amap.arg("bwC", bwC, "Bitwidth of C");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    cout << "[Test] Creating NetIO..." << endl;
    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    cout << "[Test] Creating OTPack..." << endl;
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    cout << "[Test] Creating LinearOT..." << endl;
    LinearOT *prod = new LinearOT(party, io, otpack);

    cout << "[Test] Allocating and initializing input arrays..." << endl;
    uint64_t *inA = new uint64_t[dim];
    uint64_t *inB = new uint64_t[dim];
    uint64_t *outC = new uint64_t[dim];

    PRG128 prg;
    if (party == sci::ALICE) {
        prg.random_data(inA, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; ++i) {
            inA[i] &= (1ULL << bwA) - 1;
            inB[i] = 0;
        }
    } else {
        prg.random_data(inB, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; ++i) {
            inB[i] &= (1ULL << bwB) - 1;
            inA[i] = 0;
        }
    }
    for (int i = 0; i < dim; ++i) outC[i] = 0;

    auto start = std::chrono::high_resolution_clock::now();
    size_t comm_start = io->counter;

    cout << "[Test] Calling matmul_cross_terms..." << endl;
    MultMode mode = MultMode::None;
    prod->matmul_cross_terms(1, dim, 1, inA, inB, outC, bwA, bwB, bwC, false, mode);

    cout << "[Test] Output outC:" << endl;
    for (int i = 0; i < dim; ++i) {
        cout << "outC[" << i << "] = " << outC[i] << endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    size_t comm_end = io->counter;
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double comm_kb = (comm_end - comm_start) / 1024.0;
    cout << "Total communication: " << comm_kb << " KB" << endl;
    cout << "Elapsed time: " << elapsed_ms << " ms" << endl;

    delete[] inA;
    delete[] inB;
    delete[] outC;
    delete prod;
    delete otpack;
    delete io;
    cout << "[Test] Done." << endl;
    return 0;
} 