#include <iostream>
#include <string>
#include "utils/emp-tool.h"
#include "SCI/src/BuildingBlocks/aux-protocols.h"

using namespace std;
using namespace sci;

int party = 0, port = 8000;
string address = "127.0.0.1";
int dim = 1024;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of elements");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    AuxProtocols *aux = new AuxProtocols(party, io, otpack);

    uint8_t *x = new uint8_t[dim];
    uint8_t *y = new uint8_t[dim];
    uint8_t *z = new uint8_t[dim];

    PRG128 prg;
    prg.random_data(x, dim * sizeof(uint8_t));
    prg.random_data(y, dim * sizeof(uint8_t));
    for (int i = 0; i < dim; ++i) {
        x[i] &= 1;  // Ensure binary value
        y[i] &= 1;  // Ensure binary value
    }

    aux->AND(x, y, z, dim);

    for (int i = 0; i < std::min(dim, 10); ++i) {
        cout << "x[" << i << "]=" << int(x[i]) << ", y[" << i << "]=" << int(y[i]) << ", z[" << i << "]=" << int(z[i]) << endl;
    }

    delete[] x;
    delete[] y;
    delete[] z;
    delete aux;
    delete otpack;
    delete io;
    return 0;
} 