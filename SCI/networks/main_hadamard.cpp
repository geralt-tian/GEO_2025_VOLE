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
    std::cout << "Step 1: Starting hadamard test..." << std::endl;
    
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("N", dim, "Number of elements");
    amap.arg("bwA", bwA, "Bitwidth of A");
    amap.arg("bwB", bwB, "Bitwidth of B");
    amap.arg("bwC", bwC, "Bitwidth of C");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);
    
    std::cout << "Step 2: Creating NetIO..." << std::endl;
    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    
    std::cout << "Step 3: Creating OTPack..." << std::endl;
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    
    std::cout << "Step 4: Creating LinearOT..." << std::endl;
    LinearOT *prod = new LinearOT(party, io, otpack);

    std::cout << "Step 5: Allocating memory..." << std::endl;
    uint64_t *A = new uint64_t[dim];
    uint64_t *B = new uint64_t[dim];
    uint64_t *C = new uint64_t[dim];
    // 为MSB分配内存
    uint8_t *msb_A = new uint8_t[dim];
    uint8_t *msb_B = new uint8_t[dim];

    std::cout << "Step 6: Generating random data..." << std::endl;
    PRG128 prg;
    prg.random_data(A, dim * sizeof(uint64_t));
    prg.random_data(B, dim * sizeof(uint64_t));
    
    // 初始化MSB
    for (int i = 0; i < dim; ++i) {
        A[i] &= (1ULL << bwA) - 1;
        B[i] &= (1ULL << bwB) - 1;
        msb_A[i] = (A[i] >> (bwA - 1)) & 1; // 设置MSB
        msb_B[i] = (B[i] >> (bwB - 1)) & 1; // 设置MSB
    }
    
    // 初始化C为0
    for (int i = 0; i < dim; ++i) {
        C[i] = 0;
    }

    // 显示指针地址和内存状态
    std::cout << "DEBUG: Memory addresses - A=" << (void*)A << ", B=" << (void*)B << ", C=" << (void*)C << std::endl;
    
    // 先尝试更小的维度
    int test_dim = 2;  // 只使用非常小的测试维度
    std::cout << "Step 7: Calling hadamard_product with very small dim=" << test_dim << "..." << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    size_t comm_start = io->counter;
    
    try {
        // 使用正确的参数调用hadamard_product
        MultMode mode = MultMode::None; // 使用默认模式
        prod->hadamard_product(test_dim, A, B, C, bwA, bwB, bwC, true, true, mode, msb_A, msb_B);
        
        std::cout << "Step 8: Small test hadamard_product completed successfully!" << std::endl;
        
        // 显示小测试的结果
        for (int i = 0; i < test_dim; ++i) {
            cout << "Small test - A[" << i << "]=" << A[i] << ", B[" << i << "]=" << B[i] 
                 << ", C[" << i << "]=" << C[i] << endl;
        }
        
        // 如果小测试成功，再尝试原始维度
        std::cout << "Step 9: Trying with full dimension " << dim << "..." << std::endl;
        
        // 重置C
        for (int i = 0; i < dim; ++i) {
            C[i] = 0;
        }
        
        prod->hadamard_product(dim, A, B, C, bwA, bwB, bwC, true, true, mode, msb_A, msb_B);
        std::cout << "Step 10: Full hadamard_product completed successfully!" << std::endl;
    }
    catch (const std::exception& e) {
        std::cout << "Exception during hadamard_product: " << e.what() << std::endl;
    }
    catch (...) {
        std::cout << "Unknown exception during hadamard_product!" << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    size_t comm_end = io->counter;
    double elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double comm_kb = (comm_end - comm_start) / 1024.0;
    cout << "Total communication: " << comm_kb << " KB" << endl;
    cout << "Elapsed time: " << elapsed_ms << " ms" << endl;

    std::cout << "Step 11: Printing results..." << std::endl;
    for (int i = 0; i < std::min(dim, 10); ++i) {
        cout << "A[" << i << "]=" << A[i] << ", B[" << i << "]=" << B[i] << ", C[" << i << "]=" << C[i] << endl;
    }

    std::cout << "Step 12: Cleaning up..." << std::endl;
    delete[] msb_A;
    delete[] msb_B;
    delete[] A;
    delete[] B;
    delete[] C;
    delete prod;
    delete otpack;
    delete io;
    
    std::cout << "Step 13: Test completed!" << std::endl;
    return 0;
} 