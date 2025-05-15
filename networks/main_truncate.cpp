#include <iostream>
#include <string>
#include "utils/emp-tool.h"
#include "SCI/src/BuildingBlocks/aux-protocols.h"
#include "SCI/src/BuildingBlocks/truncation.h"

using namespace std;
using namespace sci;

int party = 0, port = 8000;
string address = "127.0.0.1";
int dim = 1024;
int bw = 16;
int shift = 4;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of elements");
    amap.arg("bw", bw, "Bitwidth of input");
    amap.arg("shift", shift, "Shift amount");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    AuxProtocols *aux = new AuxProtocols(party, io, otpack);
    Truncation *trunc = new Truncation(party, io, otpack);

    uint64_t *A = new uint64_t[dim];
    uint64_t *B = new uint64_t[dim];

    PRG128 prg;
    prg.random_data(A, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; ++i) {
        A[i] &= (1ULL << bw) - 1;
    }

    trunc->truncate_and_reduce(dim, A, B, shift, bw);

    for (int i = 0; i < std::min(dim, 10); ++i) {
        cout << "A[" << i << "]=" << A[i] << ", B[" << i << "]=" << B[i] << endl;
    }

    delete[] A;
    delete[] B;
    delete trunc;
    delete aux;
    delete otpack;
    delete io;
    return 0;
} 