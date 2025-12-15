# Efficient and High-Accuracy Secure Two-Party Protocols for a Class of Functions with Real-number Inputs



## Setup

For setup instructions, please refer to the README file located in the `SCI` folder.

We successfully completed the compilation on Ubuntu 22.04.5 LTS with Intel(R) Xeon(R) Platinum.


## Code Structure

The project is organized as follows:

- **/SCI/tests**
  Contains all our related code, including implementations of activation functions and models.

- **/SCI/src/BuildingBlocks**
  Contains core protocol implementations:
  - geometric_perspective_protocols.cpp/h - MW protocol implementation

- **/networks**
  Contains our evaluation programs:
  - main_exp_nagx.cpp - exp(-x) implementation using VOLE
  - main_exp_nagx_iknp.cpp - exp(-x) implementation using IKNP-OT
  - main_gelu.cpp - GELU activation function


## Benchmark

To run the exp_nagx benchmark from 2^10 to 2^18 elements:

```bash
bash test_exp_nagx_benchmark.sh
```

## Reference Repository

**Project webpage:** [BOLT](https://github.com/Clive2312/EzPC/tree/bert/SCI)

**Project webpage:** [SEAF](https://github.com/geralt-tian/SEAF)
