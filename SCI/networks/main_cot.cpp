#include <iostream>
#include <string>
#include <vector>
#include "utils/emp-tool.h"
#include "OT/ot-utils.h"
#include "OT/emp-ot.h"

using namespace std;
using namespace sci;

int party = 0, port = 8000;
string address = "127.0.0.1";
int dim = 1024;  // 使用小维度

int main(int argc, char **argv) {
    std::cout << "Step 1: Starting basic OT test..." << std::endl;
    
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.parse(argc, argv);
    
    std::cout << "Step 2: Creating NetIO..." << std::endl;
    sci::NetIO *io = new sci::NetIO(party == sci::ALICE ? nullptr : address.c_str(), port);
    std::cout << "NetIO created successfully" << std::endl;
    
    // 延迟初始化OTPack，并检查iknp_straight是否存在
    std::cout << "Step 3: Creating OTPack..." << std::endl;
    sci::OTPack<sci::NetIO> *otpack = new sci::OTPack<sci::NetIO>(io, party);
    std::cout << "OTPack created successfully" << std::endl;
    
    // 检查iknp_straight是否初始化
    if (otpack->iknp_straight == nullptr) {
        std::cout << "ERROR: otpack->iknp_straight is NULL" << std::endl;
        delete otpack;
        delete io;
        return 1;
    }
    std::cout << "otpack->iknp_straight is valid" << std::endl;
    
    // 基本的COT通信量测试
    try {
        std::cout << "Step 4: Testing basic COT (non-batched) with comm stats..." << std::endl;
        
        const int num_ot = 1;  // 只做几个OT
        const int bitlen = 21;  // 8位长度
        
        uint64_t *data0 = new uint64_t[num_ot];  // ALICE的数据
        uint64_t *recv_data = new uint64_t[num_ot];  // BOB接收的数据
        uint64_t *corr = new uint64_t[num_ot];  // 相关性数据
        bool *choices = new bool[num_ot];  // BOB的选择位
        
        // 初始化数据
        for(int i = 0; i < num_ot; i++) {
            data0[i] = (i+10);  // 简单值
            corr[i] = (i+20);   // 简单值
            choices[i] = (i % 2 == 0);  // 交替的选择位
        }
        
        // 基本的COT测试 (不使用批处理版本)
        if(party == sci::ALICE) {
            std::cout << "ALICE: Sending basic COT..." << std::endl;
            size_t comm_before = io->counter;
            otpack->iknp_straight->send_cot(data0, corr, num_ot, bitlen);
            size_t comm_after = io->counter;
            std::cout << "ALICE: send_cot comm: before=" << comm_before << ", after=" << comm_after
                      << ", delta=" << (comm_after - comm_before) << " bytes" << std::endl;
            std::cout << "ALICE: Basic COT sent successfully" << std::endl;
        } else {
            std::cout << "BOB: Receiving basic COT..." << std::endl;
            size_t comm_before = io->counter;
            otpack->iknp_straight->recv_cot(recv_data, choices, num_ot, bitlen);
            size_t comm_after = io->counter;
            std::cout << "BOB: recv_cot comm: before=" << comm_before << ", after=" << comm_after
                      << ", delta=" << (comm_after - comm_before) << " bytes" << std::endl;
            std::cout << "BOB: Basic COT received successfully" << std::endl;
            // 打印收到的数据
            for(int i = 0; i < num_ot; i++) {
                std::cout << "BOB received: data[" << i << "] = " << recv_data[i] 
                          << ", choice = " << (choices[i] ? "1" : "0") << std::endl;
            }
        }
        
        // 清理内存
        delete[] data0;
        delete[] recv_data;
        delete[] corr;
        delete[] choices;
        
        std::cout << "Basic OT comm test completed successfully" << std::endl;
        
        // batched cot相关代码全部注释掉
        /*
        // 如果基本OT正常，尝试一个小型的batched_cot测试
        std::cout << "Step 5: Testing batched COT with minimal parameters..." << std::endl;
        ...
        */
    } catch(const std::exception& e) {
        std::cout << "Exception during OT test: " << e.what() << std::endl;
    } catch(...) {
        std::cout << "Unknown exception during OT test!" << std::endl;
    }
    
    std::cout << "Step 6: Cleaning up..." << std::endl;
    delete otpack;
    delete io;
    
    std::cout << "Test completed!" << std::endl;
    return 0;
} 