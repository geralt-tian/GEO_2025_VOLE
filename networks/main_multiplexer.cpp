#include <iostream>
#include <string>
#include "utils/emp-tool.h"
#include "SCI/src/BuildingBlocks/aux-protocols.h"
#include <chrono>

using namespace std;
using namespace sci;

int party = 0, port = 8000;
string address = "127.0.0.1";
int dim = 1024;
int bw_x = 8;
int bw_y = 8;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of elements");
    amap.arg("bw_x", bw_x, "Bitwidth of x");
    amap.arg("bw_y", bw_y, "Bitwidth of y");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    AuxProtocols *aux = new AuxProtocols(party, io, otpack);

    uint8_t *sel = new uint8_t[dim];
    uint64_t *x = new uint64_t[dim];
    uint64_t *y = new uint64_t[dim];

    PRG128 prg;
    prg.random_data(sel, dim * sizeof(uint8_t));
    prg.random_data(x, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; ++i) {
        sel[i] &= 1;
        x[i] &= (1ULL << bw_x) - 1;
    }

    auto start = std::chrono::high_resolution_clock::now();
    size_t comm_start = io->counter;
    aux->multiplexer(sel, x, y, dim, bw_x, bw_y);

    for (int i = 0; i < std::min(dim, 10); ++i) {
        cout << "sel[" << i << "]=" << int(sel[i]) << ", x[" << i << "]=" << x[i] << ", y[" << i << "]=" << y[i] << endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    size_t comm_end = io->counter;
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double comm_kb = (comm_end - comm_start) / 1024.0;
    cout << "Total communication: " << comm_kb << " KB" << endl;
    cout << "Elapsed time: " << elapsed_ms << " ms" << endl;

    delete[] sel;
    delete[] x;
    delete[] y;
    delete aux;
    delete otpack;
    delete io;
    return 0;
} 