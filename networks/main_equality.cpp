#include <iostream>
#include <string>
#include "utils/emp-tool.h"
#include "SCI/src/BuildingBlocks/aux-protocols.h"
#include "Millionaire/equality.h"

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
    Equality<sci::NetIO> *eq = new Equality<sci::NetIO>(party, io, otpack);

    uint64_t *data = new uint64_t[dim];
    uint8_t *res_eq = new uint8_t[dim];

    PRG128 prg;
    prg.random_data(data, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; ++i) {
        data[i] &= (1ULL << bw) - 1;
        // 设置一些相等的值用于测试
        if (i % 5 == 0) {
            data[i] = 42; // 固定一些值为42，这样两方都有相同的值
        }
    }

    eq->check_equality(res_eq, data, dim, bw);

    for (int i = 0; i < std::min(dim, 10); ++i) {
        cout << "data[" << i << "]=" << data[i] << ", eq[" << i << "]=" << int(res_eq[i]) << endl;
    }

    delete[] data;
    delete[] res_eq;
    delete eq;
    delete otpack;
    delete io;
    return 0;
} 