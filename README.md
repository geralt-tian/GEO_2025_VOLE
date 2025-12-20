# Efficient and High-Accuracy Secure Two-Party Protocols for a Class of Functions with Real-number Inputs



## Setup

For setup instructions, please refer to the README file located in the `SCI` folder.

We successfully completed the compilation on Ubuntu 22.04.5 LTS with Intel(R) Xeon(R) Platinum.

### Requirements

* openssl 
* c++ compiler (>= 8.0 for the better performance on AVX512)
* cmake >= 3.13
* git
* make
* OpenMP (optional, only needed by CryptFlow2 for multi-threading)

### Building Dependencies
* Run `bash scripts/build-deps.sh` which will build the following dependencies
	* [emp-tool](https://github.com/emp-toolkit/emp-tool) We follow the implementation in SCI that using emp-tool for network io and pseudo random generator.
	* [emp-ot](https://github.com/emp-toolkit/emp-ot) We use Ferret in emp-ot as our VOLE-style OT.
	* [Eigen](https://github.com/libigl/eigen) We use Eigen for tensor operations.
	* [SEAL](https://github.com/microsoft/SEAL) We use SEAL's implementation for the BFV homomorphic encryption scheme.
	* [zstd](https://github.com/facebook/zstd) We use zstd for compressing the ciphertext in SEAL which can be replaced by any other compression library.
	* [hexl](https://github.com/intel/hexl/tree/1.2.2) We need hexl's AVX512 acceleration for achieving the reported numbers in our paper.

* The generated objects are placed in the `build/deps/` folder.
* Build has passed on the following setting
  * MacOS 11.6 with clang 13.0.0, Intel Core i5, cmake 3.22.1
  * Red Hat 7.2.0 with gcc 7.2.1, Intel(R) Xeon(R), cmake 3.12.0
  * Ubuntu 18.04 with gcc 7.5.0 Intel(R) Xeon(R),  cmake 3.13
  * Ubuntu 20.04 with gcc 9.4.0 Intel(R) Xeon(R),  cmake 3.16.3
  
### Building the Project

* Run `bash scripts/build-deps.sh` to build dependencies first
* Then build the project:
  ```bash
  cd build
  cmake ..
  make -j
  ```
* This will build the `exp_nagx-cheetah` executable in the `build/bin` folder

## Code Structure

The project is organized as follows:

- **/SCI/src/BuildingBlocks**
  Contains core protocol implementations:
  - geometric_perspective_protocols.cpp/h - MW protocol implementation

- **/networks**
  Contains our evaluation programs:
  - main_exp_nagx.cpp - exp(-x) implementation using VOLE-style OT


## Benchmark

To run the exp_nagx benchmark from 2^10 to 2^18 elements:

```bash
bash test_exp_nagx_benchmark.sh
```

## Reference Repository

**Project webpage:** [OpenCheetah](https://github.com/Alibaba-Gemini-Lab/OpenCheetah)