#include <iostream>
#include <string>
#include "utils/emp-tool.h"
#include "SCI/src/BuildingBlocks/aux-protocols.h"
#include "Millionaire/millionaire_with_equality.h"
#include <chrono>

using namespace std;
using namespace sci;

int party = 0, port = 8000;
string address = "127.0.0.1";
int dim = 1024;
int bw = 16;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of elements");
    amap.arg("bw", bw, "Bitwidth of input");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    MillionaireWithEquality<sci::NetIO> *mill_eq = new MillionaireWithEquality<sci::NetIO>(party, io, otpack);

    uint64_t *data = new uint64_t[dim];
    uint8_t *res_cmp = new uint8_t[dim];
    uint8_t *res_eq = new uint8_t[dim];

    PRG128 prg;
    prg.random_data(data, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; ++i) {
        data[i] &= (1ULL << bw) - 1;
    }

    auto start = std::chrono::high_resolution_clock::now();
    size_t comm_start = io->counter;
    mill_eq->compare_with_eq(res_cmp, res_eq, data, dim, bw, true);
    auto end = std::chrono::high_resolution_clock::now();
    size_t comm_end = io->counter;
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double comm_kb = (comm_end - comm_start) / 1024.0;
    cout << "Total communication: " << comm_kb << " KB" << endl;
    cout << "Elapsed time: " << elapsed_ms << " ms" << endl;

    for (int i = 0; i < std::min(dim, 10); ++i) {
        cout << "data[" << i << "]=" << data[i] << ", cmp[" << i << "]=" << int(res_cmp[i]) << ", eq[" << i << "]=" << int(res_eq[i]) << endl;
    }

    delete[] data;
    delete[] res_cmp;
    delete[] res_eq;
    delete mill_eq;
    delete otpack;
    delete io;
    return 0;
} 