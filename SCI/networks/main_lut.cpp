#include <iostream>
#include <string>
#include <chrono>
#include "LinearOT/linear-ot.h"
#include "utils/emp-tool.h"
#include "SCI/src/BuildingBlocks/aux-protocols.h"

using namespace std;
using namespace sci;

int party = 0, port = 8000;
string address = "127.0.0.1";
int dim = 1024;
int bw_x = 8;   // 必须 <= 8
int bw_y = 16;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of LUTs");
    amap.arg("bw_x", bw_x, "Bitwidth of x");
    amap.arg("bw_y", bw_y, "Bitwidth of y");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    uint64_t mask_x = (1ULL << bw_x) - 1;

    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    AuxProtocols *aux = new AuxProtocols(party, io, otpack);

    auto start = std::chrono::high_resolution_clock::now();
    size_t comm_start = io->counter;

    if (party == sci::ALICE) {
        uint64_t **spec = new uint64_t*[dim];
        for (int i = 0; i < dim; ++i) {
            spec[i] = new uint64_t[1 << bw_x];
            for (uint64_t j = 0; j < (1ULL << bw_x); ++j) {
                spec[i][j] = (i * 17 + j * 31) % (1ULL << bw_y);
            }
        }
        aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, bw_x, bw_y);
        for (int i = 0; i < dim; ++i) delete[] spec[i];
        delete[] spec;
    } else if (party == sci::BOB) {
        uint64_t *x = new uint64_t[dim];
        uint64_t *y = new uint64_t[dim];
        PRG128 prg;
        prg.random_data(x, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; ++i) x[i] &= mask_x; // 保证范围
        aux->lookup_table<uint64_t>(nullptr, x, y, dim, bw_x, bw_y);
        for (int i = 0; i < std::min(dim, 10); ++i) {
            std::cout << "x[" << i << "]=" << x[i] << ", y[" << i << "]=" << y[i] << std::endl;
        }
        delete[] x;
        delete[] y;
    }

    auto end = std::chrono::high_resolution_clock::now();
    size_t comm_end = io->counter;
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    cout << "LUT protocol done." << endl;
    double comm_kb = (comm_end - comm_start) / 1024.0;
    cout << "Total communication: " << comm_kb << " KB" << endl;
    cout << "Elapsed time: " << elapsed_ms << " ms" << endl;

    delete aux;
    delete otpack;
    delete io;
    return 0;
}