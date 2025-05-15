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
int bw = 16;

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of elements");
    amap.arg("bw", bw, "Bitwidth of values");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);

    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    AuxProtocols *aux = new AuxProtocols(party, io, otpack);

    uint64_t *inA = new uint64_t[dim];
    uint8_t *msb = new uint8_t[dim];
    uint64_t *outC = new uint64_t[dim];

    PRG128 prg;
    prg.random_data(inA, dim * sizeof(uint64_t));
    
    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t half = 1ULL << (bw - 1);
    
    // 初始化数据，创建多种测试场景
    for (int i = 0; i < dim; ++i) {
        inA[i] &= mask_bw;
        
        // 生成多种测试场景:
        // 1. MSB=0, 值<half
        // 2. MSB=0, 值>=half
        // 3. MSB=1, 值<half
        // 4. MSB=1, 值>=half
        
        if (i % 4 == 0) {
            // 场景1: MSB=0, 值<half
            msb[i] = 0;
            inA[i] = (inA[i] % half);
        } 
        else if (i % 4 == 1) {
            // 场景2: MSB=0, 值>=half
            msb[i] = 0;
            inA[i] = half + (inA[i] % half);
        }
        else if (i % 4 == 2) {
            // 场景3: MSB=1, 值<half
            msb[i] = 1;
            inA[i] = (inA[i] % half);
        }
        else {
            // 场景4: MSB=1, 值>=half
            msb[i] = 1;
            inA[i] = half + (inA[i] % half);
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    size_t comm_start = io->counter;

    // 测试clear_MSB_to_Wrap_bitMul函数
    aux->clear_MSB_to_Wrap_bitMul(dim, inA, msb, outC, bw);

    // 打印前20个元素，包含所有场景
    for (int i = 0; i < std::min(dim, 20); ++i) {
        cout << "Case " << (i % 4 + 1) << ": inA[" << i << "]=" << inA[i] 
             << ", msb[" << i << "]=" << int(msb[i])
             << ", half=" << half
             << ", outC[" << i << "]=" << outC[i] << endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    size_t comm_end = io->counter;
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double comm_kb = (comm_end - comm_start) / 1024.0;
    cout << "Total communication: " << comm_kb << " KB" << endl;
    cout << "Elapsed time: " << elapsed_ms << " ms" << endl;

    delete[] inA;
    delete[] msb;
    delete[] outC;
    delete aux;
    delete otpack;
    delete io;
    return 0;
} 