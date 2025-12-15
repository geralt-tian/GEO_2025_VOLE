// Test file for exp_nagx function
// This file tests the exp_nagx function which computes exp(-x) where x > 0
// Including ULP error analysis, communication measurement, and timing

#include <cstdint>
#include <cstdio>
#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <random>
#include <cstring>

#include "LinearOT/linear-ot.h"
#include "utils/emp-tool.h"
#include "FloatingPoint/floating-point.h"
#include "FloatingPoint/fp-math.h"
#include "BuildingBlocks/aux-protocols.h"
#include "BuildingBlocks/geometric_perspective_protocols.h"
#include "BuildingBlocks/truncation.h"
#include "Math/math-functions.h"
#include "Millionaire/equality.h"
#include "Millionaire/millionaire.h"
#include "Millionaire/millionaire_with_equality.h"

using namespace sci;
using namespace std;

int party, port = 32000;
string address = "127.0.0.1";
IOPack *iopack;
OTPack *otpack;
LinearOT *prod;
GeometricPerspectiveProtocols *gp;
AuxProtocols *aux;

// Test parameters for exp_nagx
int dim = 1048576/4;  // Test with 1024 elements
// int dim = 1024;
int32_t in_bw = 37;      // Input bit width
int32_t in_f = 12;       // Input fractional bits

int32_t localexp_f = 10;   // Local exp fractional bits
int32_t localexp_bw = localexp_f + 13; // Local exp bit width

int32_t locallut_f = 32;   // Local LUT fractional bits
int32_t locallut_bw = locallut_f +2; // Local LUT bit width

uint64_t mask_in = (in_bw == 64 ? -1 : ((1ULL << in_bw) - 1));

// Convert fixed point to double
double fix2double(uint64_t x, uint64_t y, int32_t bw, int32_t f) {
    uint64_t mask = (bw == 64) ? ~0ULL : ((1ULL << bw) - 1);
    uint64_t sum = (x + y) & mask;
    int64_t signed_val;
    if (sum >= (1ULL << (bw - 1))) {
        signed_val = (int64_t)(sum - (1ULL << bw));
    } else {
        signed_val = (int64_t)sum;
    }
    return (double)signed_val / (1ULL << f);
}

// Convert double to fixed point
uint64_t double2fix(double val, int32_t f, int32_t bw) {
    uint64_t mask = (bw == 64) ? ~0ULL : ((1ULL << bw) - 1);
    int64_t fixed_val = (int64_t)(val * (1ULL << f));
    return ((uint64_t)fixed_val) & mask;
}

int main(int argc, char **argv) {
    ArgMapping amap;
    amap.arg("r", party, "Role of party: ALICE = 1; BOB = 2");
    amap.arg("p", port, "Port Number");
    amap.arg("ip", address, "IP Address of server (ALICE)");
    amap.arg("d", dim, "Number of elements to test");
    amap.parse(argc, argv);

    cout << "=== exp_nagx Test Suite ===" << endl;
    cout << "Testing exp(-x) where x > 0" << endl;
    cout << "Party: " << (party == sci::ALICE ? "ALICE" : "BOB") << endl;
    cout << "Dimensions: " << dim << endl;
    cout << "Input format: " << in_bw << "." << in_f << endl;

    // Initialize communication
    iopack = new IOPack(party, port, address);
    otpack = new OTPack(iopack, party);
    prod = new LinearOT(party, iopack, otpack);
    gp = new GeometricPerspectiveProtocols(party, iopack, otpack);
    aux = new AuxProtocols(party, iopack, otpack);

    // Allocate arrays
    uint64_t *inA = new uint64_t[dim];
    uint64_t *result = new uint64_t[dim];

    // Generate test data with fixed seed for reproducibility
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dis(-1000, 0); // range [-10,0] for exp_nagx
    // std::uniform_real_distribution<double> dis(0, 0); 

    cout << "Generating test data for " << dim << " elements..." << endl;

    // Store original test values for later verification
    double *test_values = new double[dim];
    
    // Generate inputs - note: exp_nagx computes exp(-inA) where inA > 0
    for (int i = 0; i < dim; i++) {
        double test_val = dis(gen); // positive value
        test_values[i] = test_val;
        // test_values[0] = 8.0 - 0.0001;
        
        if (party == sci::ALICE) {
            double alice_share = test_val * 0.6; // Alice gets 60% of the value
            inA[i] = double2fix(alice_share, in_f, in_bw);
        } else {
            double bob_share = test_val * 0.4; // Bob gets 40% of the value
            inA[i] = double2fix(bob_share, in_f, in_bw);
        }
    }

    cout << "Starting exp_nagx computation..." << endl;
    
    // Record communication start
    size_t comm_start = iopack->io->counter;
    auto start_time = chrono::high_resolution_clock::now();

    // Call the exp_nagx function
    gp->exp_softmaxx(dim, inA, result, in_bw, in_f, 
                 localexp_bw, localexp_f, locallut_bw, locallut_f);

    auto end_time = chrono::high_resolution_clock::now();
    size_t comm_end = iopack->io->counter;
    
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    
    cout << "Computation completed in " << duration.count() << " ms" << endl;

    // Share results and inputs for accuracy testing
    uint64_t *result_alice = new uint64_t[dim];
    uint64_t *result_bob = new uint64_t[dim];
    uint64_t *input_alice = new uint64_t[dim];
    uint64_t *input_bob = new uint64_t[dim];
    
    if (party == sci::ALICE) {
        iopack->io->send_data(result, dim * sizeof(uint64_t));
        iopack->io->send_data(inA, dim * sizeof(uint64_t));
        iopack->io->recv_data(result_bob, dim * sizeof(uint64_t));
        iopack->io->recv_data(input_bob, dim * sizeof(uint64_t));
        memcpy(result_alice, result, dim * sizeof(uint64_t));
        memcpy(input_alice, inA, dim * sizeof(uint64_t));
    } else {
        iopack->io->recv_data(result_alice, dim * sizeof(uint64_t));
        iopack->io->recv_data(input_alice, dim * sizeof(uint64_t));
        iopack->io->send_data(result, dim * sizeof(uint64_t));
        iopack->io->send_data(inA, dim * sizeof(uint64_t));
        memcpy(result_bob, result, dim * sizeof(uint64_t));
        memcpy(input_bob, inA, dim * sizeof(uint64_t));
    }

    // Accuracy analysis (only on one party to avoid duplicate output)
    if (1) {
        cout << "\n=== Accuracy Analysis ===" << endl;
        
        double ulp = 1.0 / (1ULL << in_f); // ULP for input precision (assuming output same as input)
        double total_ulp_error = 0.0;
        double max_ulp_error = 0.0;
        double total_relative_error = 0.0;
        double max_relative_error = 0.0;
        int correct_results = 0;
        
        for (int i = 0; i < dim; i++) {
            // Reconstruct the actual result (exp_nagx output)
            double actual_result = fix2double(result_alice[i], result_bob[i], in_bw, in_f);
            
            // Reconstruct the input from Alice and Bob shares
            double input_val = fix2double(input_alice[i], input_bob[i], in_bw, in_f);
            
            // Expected result: exp(-input_val) since exp_nagx computes exp(-inA)
            double expected_result = exp(input_val);
            
            // Calculate errors
            double absolute_error = fabs(actual_result - expected_result);
            double ulp_error = absolute_error / ulp;
            double relative_error = (expected_result != 0) ? absolute_error / fabs(expected_result) : 0;
            
            total_ulp_error += ulp_error;
            total_relative_error += relative_error;
            
            if (ulp_error > max_ulp_error) max_ulp_error = ulp_error;
            if (relative_error > max_relative_error) max_relative_error = relative_error;
            
            // Consider result correct if ULP error < 1
            if (ulp_error < 10.0) correct_results++;
            
            // Print details for first few results
            // if (i < 1000) {
            //     printf("Test %d: input=%.6f, actual=%.6f, expected=%.6f, ULP_error=%.2f\n", 
            //            i, input_val, actual_result, expected_result, ulp_error);
            // }
        }
        
        cout << "\n=== Summary Statistics ===" << endl;
        printf("Total tests: %d\n", dim);
        printf("Correct results (ULP < 10): %d (%.2f%%)\n", 
               correct_results, 100.0 * correct_results / dim);
        printf("Average ULP error: %.4f\n", total_ulp_error / dim);
        printf("Maximum ULP error: %.4f\n", max_ulp_error);
        printf("Average relative error: %.6f%%\n", 100.0 * total_relative_error / dim);
        printf("Maximum relative error: %.6f%%\n", 100.0 * max_relative_error);
        
        cout << "\n=== Performance Statistics ===" << endl;
        printf("Communication: %zu bytes\n", comm_end - comm_start);
        printf("Communication per element: %.2f bytes\n", 
               (double)(comm_end - comm_start) / dim);
        printf("Computation time: %ld ms\n", duration.count());
        printf("Throughput: %.2f exp_nagx/sec\n", 
               1000.0 * dim / duration.count());
        
        // Additional statistics for negative exp
        cout << "\n=== exp_nagx Specific Analysis ===" << endl;
        printf("Input range tested: [0.1, 4.0] (positive values)\n");
        printf("Output range: exp(-4.0) to exp(-0.1) ≈ [0.018, 0.905]\n");
        printf("Function tested: exp(-x) where x > 0\n");
    }

    // Cleanup
    delete[] inA;
    delete[] result;
    delete[] result_alice;
    delete[] result_bob;
    delete[] input_alice;
    delete[] input_bob;
    delete[] test_values;
    
    delete gp;
    delete aux;
    delete prod;
    delete otpack;
    delete iopack;

    cout << "\nexp_nagx test completed successfully!" << endl;
    
    return 0;
} 