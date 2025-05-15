#include <iostream>
#include <string>
#include "utils/emp-tool.h"
#include "SCI/src/BuildingBlocks/value-extension.h"

using namespace std;
using namespace sci;

int party = 0, port = 8000;
string address = "127.0.0.1";
int dim = 1024;
int bwA = 8;
int bwB = 16;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of elements");
    amap.arg("bwA", bwA, "Bitwidth of input");
    amap.arg("bwB", bwB, "Bitwidth of output");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    XTProtocol *xt = new XTProtocol(party, io, otpack);

    uint64_t *inA = new uint64_t[dim];
    uint64_t *outB = new uint64_t[dim];
    uint8_t *msbA = new uint8_t[dim];

    PRG128 prg;
    prg.random_data(inA, dim * sizeof(uint64_t));
    for (int i = 0; i < dim; ++i) {
        inA[i] &= (1ULL << bwA) - 1;
        // 设置一些负数和正数
        if (i % 2 == 0) {
            // 设置MSB=1表示负数
            msbA[i] = 1;
            // 设置最高位为1，使其成为负数的分享
            inA[i] |= (1ULL << (bwA - 1));
        } else {
            // 设置MSB=0表示正数
            msbA[i] = 0;
            // 清除最高位，使其成为正数的分享
            inA[i] &= ~(1ULL << (bwA - 1));
        }
    }

    // 测试带MSB的符号扩展
    xt->s_extend(dim, inA, outB, bwA, bwB, msbA);

    for (int i = 0; i < std::min(dim, 10); ++i) {
        cout << "inA[" << i << "]=" << inA[i] << ", msbA[" << i << "]=" << int(msbA[i]) 
             << ", outB[" << i << "]=" << outB[i] << endl;
    }

    delete[] inA;
    delete[] outB;
    delete[] msbA;
    delete xt;
    delete otpack;
    delete io;
    return 0;
} 