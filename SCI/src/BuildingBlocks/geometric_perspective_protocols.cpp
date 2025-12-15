// anonymous authors

#include "geometric_perspective_protocols.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <math.h>

GeometricPerspectiveProtocols::GeometricPerspectiveProtocols(
    int party, sci::NetIO *io, sci::OTPack<sci::NetIO> *otpack) {
    this->party = party;
    this->io = io;
    this->otpack = otpack;
    this->aux = new AuxProtocols(party, io, otpack);
    this->mill = this->aux->mill;
    this->mill_eq = new MillionaireWithEquality<sci::NetIO>(party, io, otpack);
    this->eq = new Equality<sci::NetIO>(party, io, otpack);
    this->triple_gen = this->mill->triple_gen;
}

GeometricPerspectiveProtocols::~GeometricPerspectiveProtocols() {
    delete this->aux;
    delete this->mill_eq;
    delete this->eq;
}

void GeometricPerspectiveProtocols::new_truncate(int32_t dim, uint64_t *inA,
                                                 uint64_t *outB, int32_t shift,
                                                 int32_t bw) {
    if (shift == 0) {
        memcpy(outB, inA, sizeof(uint64_t) * dim);
        return;
    }
    assert((bw - shift) > 0 && "Truncation shouldn't truncate the full bitwidth");
    assert(bw - shift - 1 >= 0);
    assert(inA != outB);

    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t mask_shift = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_upper =
            ((bw - shift) == 64 ? -1 : ((1ULL << (bw - shift)) - 1));

    uint64_t *inA_orig = new uint64_t[dim];
    uint64_t *inA_lower = new uint64_t[dim];
    uint64_t *inA_upper = new uint64_t[dim];
    uint8_t *wrap_lower = new uint8_t[dim];
    uint8_t *wrap_upper = new uint8_t[dim];
    uint8_t *eq_upper = new uint8_t[dim];
    uint8_t *and_upper = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        inA_lower[i] = inA[i] & mask_shift;
        inA_upper[i] = (inA[i] >> shift) & mask_upper;
        if (party == sci::BOB) {
            inA_upper[i] = (mask_upper - inA_upper[i]) & mask_upper;
        }
    }
    uint64_t *inA_prime = new uint64_t[dim];
    uint64_t quarter = 1ULL << (bw - 2);
    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = (inA[i] - quarter) & mask_bw;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = inA[i];
        }
    }
    uint64_t *r = new uint64_t[dim];
    uint64_t half = 1ULL << (bw - 1);
    for (int i = 0; i < dim; i++) {
        if (inA_prime[i] >= half) {
            r[i] = 1;
        } else {
            r[i] = 0;
        }
    }
    uint64_t *bit_mul = new uint64_t[dim];
    if (party == sci::ALICE) {
        sci::PRG128 prg;
        uint64_t *data0 = new uint64_t[dim];
        prg.random_data(data0, dim * sizeof(uint64_t));
        otpack->iknp_straight->send_cot(data0, r, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = ((1ULL << shift) - data0[i]) & mask_shift;
        }
        delete[] data0;
    } else { // party == BOB
        bool *choice = new bool[dim];
        for (int i = 0; i < dim; i++) {
            choice[i] = r[i];
        }
        uint64_t *data = new uint64_t[dim];
        otpack->iknp_straight->recv_cot(data, choice, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = data[i];
        }
        delete[] choice;
    }
    uint64_t *arith_wrap_lower = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
      arith_wrap_lower[i] = 0;
    }
    this->aux->wrap_computation(inA_lower, wrap_lower, dim, shift);
    this->aux->B2A(wrap_lower, arith_wrap_lower, dim, bw);
    uint64_t offset = 1ULL << (bw - shift - 2);
    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
      outB[i] =
          (((inA_prime[i] >> shift) & mask_upper) + offset +
           arith_wrap_lower[i] - (1ULL << (bw - shift)) * (bit_mul[i] + 1)) &
                      mask_bw;
        } else {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) + arith_wrap_lower[i] -
                       (1ULL << (bw - shift)) * bit_mul[i]) &
                      mask_bw;
        }
    }
    delete[] inA_orig;
    delete[] inA_lower;
    delete[] inA_upper;
    delete[] wrap_lower;
    delete[] wrap_upper;
    delete[] eq_upper;
    delete[] and_upper;
    delete[] arith_wrap_lower;

    return;
}

void GeometricPerspectiveProtocols::mw(int32_t dim, uint64_t *input,
                                       uint64_t *output, int32_t in_bw,
                                       int32_t out_bw) {
    uint64_t mask_in = (in_bw == 64 ? -1 : ((1ULL << in_bw) - 1));
    uint64_t mask_out = (out_bw == 64 ? -1 : ((1ULL << out_bw) - 1));
    uint64_t *input_prime = new uint64_t[dim];
    uint64_t N = 1ULL << (in_bw);
    uint64_t quarter = 1ULL << (in_bw - 2);

    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
      input_prime[i] = static_cast<uint64_t>(
          floor((N - ((input[i] - quarter) & mask_in)) / (N - 2 * quarter)));
        }
    } else {
        for (int i = 0; i < dim; i++) {
      input_prime[i] =
          static_cast<uint64_t>(floor(input[i] / (N - 2 * quarter)));
        }
    }
    uint64_t *r = new uint64_t[dim];
    uint64_t half = 1ULL << (in_bw - 1);
    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
            if (input_prime[i] == 0) {
                r[i] = 1;
            } else {
                r[i] = 0;
            }
        }
  } else {
        for (int i = 0; i < dim; i++) {
            if (input_prime[i] > 0) {
                r[i] = 1;
            } else {
                r[i] = 0;
            }
        }
    }
    uint64_t *mul = new uint64_t[dim];
    bit_mul(dim, r, mul, out_bw);

    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
            output[i] = mul[i] + 1;
            if (input[i] < quarter) {
                output[i] = (output[i] - 1) & mask_out;
            }
        } else {
            output[i] = mul[i] & mask_out;
        }
    }

    // for (int i = 0; i < 1000; i++){
    //     printf("input[%d]: %llu\n", i, input[i]);
    //     printf("input_prime[%d]: %llu\n", i, input_prime[i]);
    //     printf("mul[%d]: %llu\n", i, mul[i]);
    //     printf("output[%d]: %llu\n", i, output[i]);
    // }
    // printf ("in_bw: %d\n", in_bw);
    // printf ("N- 2*quarter: %llu\n", N- 2*quarter);
}

void GeometricPerspectiveProtocols::mw_conversion(int32_t dim, uint64_t *input,
                                       uint64_t *output, int32_t in_bw, int32_t act_l_bw,
                                       int32_t out_bw) {
  uint64_t mask_lplus1 = (act_l_bw + 1 == 64 ? -1 : ((1ULL << (act_l_bw + 1)) - 1));
  uint64_t *input_lplus1 = new uint64_t[dim];
  uint64_t *lplus1_output = new uint64_t[dim];
  for (int i = 0; i < dim; i++) {
    input_lplus1[i] = input[i] & mask_lplus1;
  }
  mw(dim, input_lplus1, lplus1_output, act_l_bw + 1, out_bw);
  uint64_t *y_heat = new uint64_t[dim];
  
  if (party == sci::ALICE) {
    for (int i = 0; i < dim; i++) {
      y_heat[i] = (input_lplus1[i] + (1ULL << act_l_bw)) & mask_lplus1;
    }
  } else {
    for (int i = 0; i < dim; i++) {
      y_heat[i] = (input_lplus1[i] + (1ULL << act_l_bw)) & mask_lplus1;
    }
  }
  uint64_t *y_star = new uint64_t[dim];
  uint64_t *delta = new uint64_t[dim];
  if (party == sci::ALICE) {
    for (int i = 0; i < dim; i++) {
      y_star[i] = (y_heat[i] - (1ULL << (act_l_bw-1))) & mask_lplus1;
      if (y_heat[i] < (1ULL << (act_l_bw-1))) {
        delta[i] = 1;
      } else {
        delta[i] = 0;
      }
    }
  } else {
    for (int i = 0; i < dim; i++) {
      delta[i] = 0;
    }
  }
  uint64_t *bit_mul_input = new uint64_t[dim];
  if (party == sci::ALICE) {
    for (int i = 0; i < dim; i++) {
      if (y_star[i] >= (1ULL << (act_l_bw))) {
        bit_mul_input[i] = 1;
      } else {
        bit_mul_input[i] = 0;
      }
    }
  } else {
    for (int i = 0; i < dim; i++) {
      if (y_heat[i] >= (1ULL << (act_l_bw))) {
        bit_mul_input[i] = 1;
      } else {
        bit_mul_input[i] = 0;
      }
    }
  }
  uint64_t *bit_mul_output = new uint64_t[dim];
  bit_mul(dim, bit_mul_input, bit_mul_output, out_bw);
  uint64_t *c = new uint64_t[dim];
  uint64_t mask_out = (out_bw == 64 ? -1 : ((1ULL << out_bw) - 1));
  for (int i = 0; i < dim; i++) {
    c[i] = (delta[i] - bit_mul_output[i]) & mask_out;
  }
  for (int i = 0; i < dim; i++) {
    output[i] = (lplus1_output[i] - c[i]) & mask_out;
  }

  // for (int i = 0; i < dim; i++) {
  //   printf("input[%d]: %llu\n", i, input[i]);
  //   printf("input_lplus1[%d]: %llu\n", i, input_lplus1[i]);
  //   printf("y_heat[%d]: %llu\n", i, y_heat[i]);
  //   printf("y_star[%d]: %llu\n", i, y_star[i]);
  //   printf("delta[%d]: %llu\n", i, delta[i]);
  //   printf("bit_mul_input[%d]: %llu\n", i, bit_mul_input[i]);
  //   printf("bit_mul_output[%d]: %llu\n", i, bit_mul_output[i]);
  //   printf("c[%d]: %llu\n", i, c[i]);
  //   printf("lplus1_output[%d]: %llu\n", i, lplus1_output[i]);
  //   printf("output[%d]: %llu\n", i, output[i]);
  //   printf("\n");
  // }
  delete[] input_lplus1;
  delete[] y_heat;
  delete[] y_star;
  delete[] delta;
  delete[] bit_mul_input;
  delete[] bit_mul_output;
  delete[] c;
  delete[] lplus1_output;
}

void GeometricPerspectiveProtocols::mwwithB(int32_t dim, uint64_t B,
                                            uint64_t *input, uint64_t *output,
                                            int32_t in_bw, int32_t out_bw) {
uint64_t mask_in = (in_bw == 64 ? -1 : ((1ULL << in_bw) - 1));
  uint64_t N = pow(2, in_bw);
uint64_t l_star;
uint64_t mask_out = (out_bw == 64 ? -1 : ((1ULL << out_bw) - 1));
uint64_t *compare_input = new uint64_t[dim];
  if (B == N / 2) {
    printf("B == N / 2\n");
    l_star = in_bw;
  } else {
    l_star = ceil(log2(floor(N / (N - 2 * B))));
}

uint8_t *M = new uint8_t[dim];
uint8_t *M_eq = new uint8_t[dim];
uint8_t *delta = new uint8_t[dim];
uint64_t *x0_star = new uint64_t[dim];
uint64_t *M_result = new uint64_t[dim];
  if (party == sci::ALICE) {
    for (int i = 0; i < dim; i++) {
      if (input[i] > B) { 
    delta[i] = 1;
      } else {
    delta[i] = 0;
    }
    x0_star[i] = (input[i] - B) & mask_in;
   }
  } else {
    for (int i = 0; i < dim; i++) {
        delta[i] = 0;
    }
}

  if (B <= (3 * N / 8)) {
    // mw(dim, input, output, in_bw, out_bw);
    return;
  } else { // B > 3N/8
    if (B == (N / 2)) {
        // compare(N - x0,x1)
        // uint64_t *compare_input = new uint64_t[dim];
        printf("B == N / 2\n");
        printf("l_star: %llu\n", l_star);
      if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
                compare_input[i] = N - x0_star[i];
            }
        this->mill_eq->compare_with_eq(M, M_eq, compare_input, dim, l_star,
                                       false);
      } else {
        for (int i = 0; i < dim; i++) {
                compare_input[i] = input[i];
            }
        this->mill_eq->compare_with_eq(M, M_eq, compare_input, dim, l_star,
                                       false);
        }
        
        // delete[] compare_input;
    } else {
      printf("B != N / 2\n");
      if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
          compare_input[i] = std::floor((N - x0_star[i]) / (N - 2 * B));
        }
        this->mill_eq->compare_with_eq(M, M_eq, compare_input, dim, l_star,
                                       false);
      } else {
        for (int i = 0; i < dim; i++) {
          compare_input[i] = std::floor((input[i]) / (N - 2 * B));
        }
        this->mill_eq->compare_with_eq(M, M_eq, compare_input, dim, l_star,
                                       false);
        }
        
        // delete[] compare_input;
    }
    
    for (int i = 0; i < dim; i++) {
        // M[i] = M[i] ^ M_eq[i];
    }

    this->aux->B2A(M, M_result, dim, out_bw);
}
//   for (int i = 0; i < dim; i++) {
//     output[i] = (M_result[i] + delta[i]) & mask_out;
//     if (party == sci::ALICE){
//         printf("input[%d]: %llu\n", i, input[i]);
//         printf("input_star[%d]: %llu\n", i, x0_star[i]);
//         printf("compare_input[%d]: %llu\n", i, compare_input[i]);
//         // printf("input_delta[%d]: %llu\n", i, delta[i]);
//     }
//     else{
//         printf("input[%d]: %llu\n", i, input[i]);
//         printf("compare_input[%d]: %llu\n", i, compare_input[i]);
//         // printf("input_star[%d]: %llu\n", i, x0_star[i]);
//         // printf("input_delta[%d]: %llu\n", i, delta[i]);
//     }
//     printf("M[%d]: %lld\n", i   , M[i]);    
//     printf("M_result[%d]: %llu\n", i, M_result[i]);
//     printf("delta[%d]: %llu\n", i, delta[i]);
//     printf("output[%d]: %llu\n", i, output[i]);
// }
}

void GeometricPerspectiveProtocols::mux_3(int32_t dim, uint64_t *inA,
                                          uint64_t *inC, uint64_t *out,
                                          int32_t bw) {
    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint8_t *c0 = new uint8_t[dim];
    uint8_t *c1 = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        bool *temp = new bool[2];
        sci::int64_to_bool(temp, inC[i], 2);
        c1[i] = temp[0];
        c0[i] = temp[1];
        delete[] temp;
    }
    uint8_t *carry = new uint8_t[dim];
    if (party == sci::ALICE) {
        uint8_t *dummy = new uint8_t[dim];
        for (int i = 0; i < dim; i++) {
            dummy[i] = 0;
        }
        aux->AND(c1, dummy, carry, dim);
    } else {
        uint8_t *dummy = new uint8_t[dim];
        for (int i = 0; i < dim; i++) {
            dummy[i] = 0;
        }
        aux->AND(dummy, c1, carry, dim);
    }
    uint64_t *temp_inA = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        c0[i] = c0[i] ^ carry[i];
        temp_inA[i] = (inA[i] * 2) & mask_bw;
    }
    uint64_t *t1 = new uint64_t[dim];
    uint64_t *t2 = new uint64_t[dim];
    aux->multiplexer(c1, inA, t1, dim, bw, bw);
    aux->multiplexer(c0, temp_inA, t2, dim, bw, bw);
    for (int i = 0; i < dim; i++) {
        out[i] = (t1[i] + t2[i]) & mask_bw;
    }
}

void GeometricPerspectiveProtocols::bit_mul(int32_t dim, uint64_t *input,
                                            uint64_t *output,
                                            int32_t output_bw) {
    uint64_t mask_out = (output_bw == 64 ? -1 : ((1ULL << output_bw) - 1));
    if (party == sci::ALICE) {
        sci::PRG128 prg;
        uint64_t *data0 = new uint64_t[dim];
        prg.random_data(data0, dim * sizeof(uint64_t));
        otpack->iknp_straight->send_cot(data0, input, dim, output_bw);
        for (int i = 0; i < dim; i++) {
            output[i] = ((1ULL << output_bw) - data0[i]) & mask_out;
        }
        delete[] data0;
    } else { // party == BOB
        bool *choice = new bool[dim];
        for (int i = 0; i < dim; i++) {
            choice[i] = input[i];
        }
        uint64_t *data = new uint64_t[dim];
        otpack->iknp_straight->recv_cot(data, choice, dim, output_bw);
        for (int i = 0; i < dim; i++) {
            output[i] = data[i];
        }
        delete[] choice;
    }
}

void GeometricPerspectiveProtocols::cross_term(int32_t dim, uint64_t *inA,
                                               uint64_t *inB, uint64_t *outC,
                                               int32_t bwA, int32_t bwB,
                                               int32_t bwC) {
  for (int i = 0; i < dim; i++) {
        outC[i] = 0;
    }
    uint64_t mask = (1ULL << bwC) - 1;
    if (party == sci::ALICE) {
        sci::PRG128 prg;
        for (int i = 0; i < bwB; i++) {
            auto *data0 = new uint64_t[dim];
            prg.random_data(data0, dim * sizeof(uint64_t));
            for (int j = 0; j < dim; j++) {
                data0[j] = data0[j] & ((1ULL << (bwC - i)) - 1);
            }
            otpack->iknp_straight->send_cot(data0, inA, dim, bwC - i);
            for (int j = 0; j < dim; j++) {
                outC[j] += (-data0[j] * (1ULL << i));
                outC[j] &= mask;
            }
            delete[] data0;
        }
    } else {
        bool choice[bwB][dim];
        for (int i = 0; i < dim; i++) {
            bool *temp = new bool[bwB];
            sci::int64_to_bool(temp, inB[i], bwB);
            for (int j = 0; j < bwB; j++) {
                choice[j][i] = temp[j];
            }
            delete[] temp;
        }
        for (int i = 0; i < bwB; i++) {
            auto *data = new uint64_t[dim];
            bool *c = new bool[dim];
            for (int j = 0; j < dim; j++) {
                bool *temp = new bool[bwB];
                sci::int64_to_bool(temp, inB[j], bwB);
                c[j] = temp[i];
            }
            otpack->iknp_straight->recv_cot(data, c, dim, bwC - i);
            for (int j = 0; j < dim; j++) {
                outC[j] += (data[j] * (1ULL << i));
                outC[j] &= mask;
            }
            delete[] data;
            delete[] c;
        }
    }
}

void GeometricPerspectiveProtocols::cross_term_reverse(
    int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA,
    int32_t bwB, int32_t bwC) {
    uint64_t mask = (1ULL << bwC) - 1;
    if (party == sci::ALICE) {
        bool choice[bwA][dim];
        for (int i = 0; i < dim; i++) {
            bool *temp = new bool[bwA];
            sci::int64_to_bool(temp, inA[i], bwA);
            for (int j = 0; j < bwA; j++) {
                choice[j][i] = temp[j];
            }
            delete[] temp;
        }
        for (int i = 0; i < bwA; i++) {
            auto *data = new uint64_t[dim];
            bool *c = new bool[dim];
            for (int j = 0; j < dim; j++) {
                bool *temp = new bool[bwA];
                sci::int64_to_bool(temp, inA[j], bwA);
                c[j] = temp[i];
            }
            otpack->iknp_reversed->recv_cot(data, c, dim, bwC - i);
            for (int j = 0; j < dim; j++) {
                outC[j] += ((data[j] * (1ULL << i)) & mask);
                outC[j] &= mask;
            }
            delete[] data;
            delete[] c;
        }
    } else {
        sci::PRG128 prg;
        for (int i = 0; i < bwA; i++) {
            auto *data0 = new uint64_t[dim];
            prg.random_data(data0, dim * sizeof(uint64_t));
            otpack->iknp_reversed->send_cot(data0, inB, dim, bwC - i);
            for (int j = 0; j < dim; j++) {
                outC[j] += ((-data0[j] * (1ULL << i)) & mask);
                outC[j] &= mask;
            }
            delete[] data0;
        }
    }
}

void GeometricPerspectiveProtocols::signed_mul(int32_t dim, uint64_t *inA,
                                               uint64_t *inB, uint64_t *outC,
                                               int32_t bwA, int32_t bwB,
                                               int32_t bwC) {
    auto *c = new uint64_t[dim];
    auto *d = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        c[i] = 0;
        d[i] = 0;
    }
    if (bwA <= bwB) {
        cross_term_reverse(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term(dim, inB, inA, d, bwB, bwA, bwC);
    } else {
        cross_term(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term_reverse(dim, inB, inA, d, bwB, bwA, bwC);
    }
    auto *m_x = new uint64_t[dim];
    auto *m_y = new uint64_t[dim];
    mw(dim, inA, m_x, bwA, 2);
    mw(dim, inB, m_y, bwB, 2);
    auto *g = new uint64_t[dim];
    auto *h = new uint64_t[dim];
    mux_3(dim, inA, m_y, g, bwA);
    mux_3(dim, inB, m_x, h, bwB);
    uint64_t mask = (1ULL << bwC) - 1;
    for (int i = 0; i < dim; i++) {
    outC[i] = (inA[i] * inB[i] + c[i] + d[i] - (g[i] * (1ULL << bwB)) -
               (h[i] * (1ULL << bwA))) &
              mask;
    }
    delete[] c;
    delete[] d;
}

void GeometricPerspectiveProtocols::signed_crossterm(int32_t dim, uint64_t *inA,
                                                     uint64_t *inB,
                                                     uint64_t *outC,
                                                     int32_t bwA, int32_t bwB,
                                                     int32_t bwC) {
    uint64_t maskA = (1ULL << bwA) - 1;
    uint64_t maskB = (1ULL << bwB) - 1;
    uint64_t maskC = (1ULL << bwC) - 1;
    // uint64_t *cross_term_output = new uint64_t[dim];
  if (party == sci::ALICE) {
        cross_term(dim, inA, inB, outC, bwA, bwB, bwC);
  } else {
        cross_term(dim, inA, inB, outC, bwA, bwB, bwC);
    }
    // for (int i = 0; i < dim; i++){
    //     if (party == sci::ALICE){
  //         outC[i] = (cross_term_output[i] - (1ULL << bwB - 1)*inA[i] + (1ULL
  //         << bwA +bwB - 2)) & maskC;
    //     }
    //     else{
  //         outC[i] = (cross_term_output[i] - (1ULL << bwA - 1)*inB[i]) &
  //         maskC;
    //     }
    // }

    // for (int i = 0; i < dim; i++){
    //     printf("inA[%d]: %llu\n", i, inA[i]);
    //     printf("inB[%d]: %llu\n", i, inB[i]);
    //     // printf("cross_term_output[%d]: %llu\n", i, cross_term_output[i]);
    //     printf("outC[%d]: %llu\n", i, outC[i]);
    // }
}

void GeometricPerspectiveProtocols::signed_cipher_plainc(
    int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA,
    int32_t bwB, int32_t bwC) {
    uint64_t maskA = (1ULL << bwA) - 1;
    uint64_t maskB = (1ULL << bwB) - 1;
    uint64_t maskC = (1ULL << bwC) - 1;
  uint64_t *inA_p = new uint64_t[dim]; // inA'
    uint8_t *wrap_A_p = new uint8_t[dim];
  uint64_t *inB_p = new uint64_t[dim]; // inB'

    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
            inA_p[i] = inA[i] + (1ULL << bwA - 1) & maskA;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            inA_p[i] = inA[i] & maskA;
        }
    }

    for (int i = 0; i < dim; i++) {
        inB_p[i] = inB[i] + (1ULL << bwB - 1) & maskB;
    }
    // uint64_t inB_p = new uint64_t[dim]; //inB'
    this->aux->wrap_computation(inA_p, wrap_A_p, dim, bwA);

    if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
      outC[i] =
          (inA_p[i] * inB_p[i] -
           (1ULL << bwA) * wrap_A_p[i] * inB_p[i]  // x0*y'  -2^m * w(x) *y'
           - (1ULL << bwA - 1) * inB_p[i]          //-2^(m-1)*y'
           - (1ULL << bwB - 1) * inA_p[i]          //-2^(n-1)*x0
            + (1ULL << bwA + bwB - 1) * wrap_A_p[i] //+2^(m+n-1) * w(x)
           + (1ULL << bwA + bwB - 2)) &
          maskC; // + 2^(m+n-2)
        }
  } else {
        for (int i = 0; i < dim; i++) {
      outC[i] =
          (inA_p[i] * inB_p[i] -
           (1ULL << bwA) * wrap_A_p[i] * inB_p[i] // x1*y'  -2^m * w(x) *y'
           - (1ULL << bwB - 1) * inA_p[i]         //-2^(n-1)*x1
           + (1ULL << bwA + bwB - 1) * wrap_A_p[i]) &
          maskC; //+2^(m+n-1) * w(x)
    }
  }
}

void GeometricPerspectiveProtocols::sign_extension(int32_t dim, uint64_t *in,
                                                   uint64_t *out, int32_t in_bw,
                                                   int32_t out_bw) {
    if (in_bw == out_bw) {
      memcpy(out, in, sizeof(uint64_t) * dim);
      return;
    }
    assert(in_bw < out_bw);
    uint64_t *c = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        c[i] = 0;
    }
    mw(dim, in, c, in_bw, out_bw - in_bw);
    uint64_t M = 1ULL << in_bw;
    uint64_t N = 1ULL << out_bw;
    uint64_t mask = N - 1;
    for (int i = 0; i < dim; i++) {
        out[i] = (in[i] + c[i] * (N - M)) & mask;
    }
}

void GeometricPerspectiveProtocols::truncate_with_one_bit_error(
    int32_t dim, uint64_t *inA, uint64_t *outB, int32_t shift, int32_t bw) {
    if (shift == 0) {
        memcpy(outB, inA, sizeof(uint64_t) * dim);
        return;
    }
    assert((bw - shift) > 0 && "Truncation shouldn't truncate the full bitwidth");
    assert(inA != outB);
    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t mask_shift = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_upper =
            ((bw - shift) == 64 ? -1 : ((1ULL << (bw - shift)) - 1));

    uint64_t *inA_orig = new uint64_t[dim];
    uint64_t *inA_lower = new uint64_t[dim];
    uint64_t *inA_upper = new uint64_t[dim];
    uint8_t *wrap_lower = new uint8_t[dim];
    uint8_t *wrap_upper = new uint8_t[dim];
    uint8_t *eq_upper = new uint8_t[dim];
    uint8_t *and_upper = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        inA_lower[i] = inA[i] & mask_shift;
        inA_upper[i] = (inA[i] >> shift) & mask_upper;
        if (party == sci::BOB) {
            inA_upper[i] = (mask_upper - inA_upper[i]) & mask_upper;
        }
    }
    uint64_t *inA_prime = new uint64_t[dim];
    uint64_t quarter = 1ULL << (bw - 2);
    if (party == sci::ALICE) {
        // x - L/4
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = (inA[i] - quarter) & mask_bw;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = inA[i];
        }
    }
    uint64_t *r = new uint64_t[dim];
    uint64_t half = 1ULL << (bw - 1);
    for (int i = 0; i < dim; i++) {
        if (inA_prime[i] >= half) {
            r[i] = 1;
        } else {
            r[i] = 0;
        }
    }
    uint64_t *bit_mul = new uint64_t[dim];
    if (party == sci::ALICE) {
        uint64_t mask_shift = (1ULL << shift) - 1;
        sci::PRG128 prg;
        uint64_t *data0 = new uint64_t[dim];
        prg.random_data(data0, dim * sizeof(uint64_t));
        otpack->iknp_straight->send_cot(data0, r, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = ((1ULL << shift) - data0[i]) & mask_shift;
        }
        delete[] data0;
    } else { // party == BOB
        bool *choice = new bool[dim];
        for (int i = 0; i < dim; i++) {
            choice[i] = r[i];
        }
        uint64_t *data = new uint64_t[dim];
        otpack->iknp_straight->recv_cot(data, choice, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = data[i];
        }
        delete[] choice;
    }
    uint64_t offset = 1ULL << (bw - shift - 2);
    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) + offset -
                       (1ULL << (bw - shift)) * (bit_mul[i] + 1)) &
                      mask_bw;
        } else {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) -
                       (1ULL << (bw - shift)) * bit_mul[i]) &
                      mask_bw;
        }
    }
    delete[] inA_orig;
    delete[] inA_lower;
    delete[] inA_upper;
    delete[] wrap_lower;
    delete[] wrap_upper;
    delete[] eq_upper;
    delete[] and_upper;
    return;
}

void GeometricPerspectiveProtocols::sirnn_unsigned_mul(
    int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA,
    int32_t bwB, int32_t bwC) {
    auto *c = new uint64_t[dim]();
    auto *d = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        c[i] = 0;
        d[i] = 0;
    }
    if (bwA <= bwB) {
        cross_term_reverse(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term(dim, inB, inA, d, bwB, bwA, bwC);
    } else {
        cross_term(dim, inA, inB, c, bwA, bwB, bwC);
        cross_term_reverse(dim, inB, inA, d, bwB, bwA, bwC);
    }
    auto *wx = new uint8_t[dim];
    auto *wy = new uint8_t[dim];
    aux->wrap_computation(inA, wx, dim, bwA);
    aux->wrap_computation(inB, wy, dim, bwB);
    auto *h = new uint64_t[dim];
    auto *g = new uint64_t[dim];
    aux->multiplexer(wx, inB, h, dim, bwB, bwB);
    aux->multiplexer(wy, inA, g, dim, bwA, bwA);
    uint64_t mask = (1ULL << bwC) - 1;
    for (int i = 0; i < dim; i++) {
    outC[i] = (inA[i] * inB[i] + c[i] + d[i] - (g[i] * (1ULL << bwB)) -
               (h[i] * (1ULL << bwA))) &
              mask;
    }
    delete[] c;
    delete[] d;
}

void GeometricPerspectiveProtocols::msb0_truncation(int32_t dim, uint64_t *inA,
                                                    uint64_t *outB,
                                                    int32_t shift, int32_t bw) {
    if (shift == 0) {
        memcpy(outB, inA, sizeof(uint64_t) * dim);
        return;
    }
    assert((bw - shift) > 0 && "Truncation shouldn't truncate the full bitwidth");
    assert(inA != outB);
    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t mask_shift = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_upper =
            ((bw - shift) == 64 ? -1 : ((1ULL << (bw - shift)) - 1));
    uint64_t *inA_orig = new uint64_t[dim];
    uint64_t *inA_lower = new uint64_t[dim];
    uint64_t *inA_upper = new uint64_t[dim];
    uint8_t *wrap_lower = new uint8_t[dim];
    uint8_t *wrap_upper = new uint8_t[dim];
    uint8_t *eq_upper = new uint8_t[dim];
    uint8_t *and_upper = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        inA_lower[i] = inA[i] & mask_shift;
        inA_upper[i] = (inA[i] >> shift) & mask_upper;
    }
    uint8_t *r = new uint8_t[dim];
    uint64_t half = 1ULL << (bw - 1);
    for (int i = 0; i < dim; i++) {
        if (inA[i] >= half) {
            r[i] = 1;
        } else {
            r[i] = 0;
        }
    }
    uint8_t *and_result = new uint8_t[dim];
    uint8_t *dummy = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        dummy[i] = 0;
    }
    if (party == sci::ALICE) {
        this->aux->AND(r, dummy, and_result, dim);
    } else {
        this->aux->AND(dummy, r, and_result, dim);
    }
    for (int i = 0; i < dim; i++) {
        and_result[i] = and_result[i] ^ r[i];
    }
  uint64_t *mw = new uint64_t[dim];
    this->aux->B2A(and_result, mw, dim, shift);
    uint64_t *arith_wrap_lower = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
      arith_wrap_lower[i] = 0;
    }
//    this->aux->wrap_computation(inA_lower, wrap_lower, dim, shift);
//    this->aux->B2A(wrap_lower, arith_wrap_lower, dim, bw);
    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
      outB[i] = (((inA[i] >> shift) & mask_upper) + arith_wrap_lower[i] -
                 (1ULL << (bw - shift)) * (mw[i])) &
                mask_bw;
        } else {
      outB[i] = (((inA[i] >> shift) & mask_upper) + arith_wrap_lower[i] -
                 (1ULL << (bw - shift)) * mw[i]) &
                mask_bw;
        }
    }
    delete[] inA_orig;
    delete[] inA_lower;
    delete[] inA_upper;
    delete[] wrap_lower;
    delete[] wrap_upper;
    delete[] eq_upper;
    delete[] and_upper;
    delete[] arith_wrap_lower;
    return;
}

void GeometricPerspectiveProtocols::new_msb0_truncation(
    int32_t dim, uint64_t *inA, uint64_t *outB, int32_t shift, int32_t bw) {
    if (shift == 0) {
        memcpy(outB, inA, sizeof(uint64_t) * dim);
        return;
    }
    assert((bw - shift) > 0 && "Truncation shouldn't truncate the full bitwidth");
    assert(inA != outB);
    uint64_t mask_bw = (bw == 64 ? -1 : ((1ULL << bw) - 1));
    uint64_t mask_shift = (shift == 64 ? -1 : ((1ULL << shift) - 1));
    uint64_t mask_upper =
            ((bw - shift) == 64 ? -1 : ((1ULL << (bw - shift)) - 1));

    uint64_t *inA_orig = new uint64_t[dim];
    uint64_t *inA_lower = new uint64_t[dim];
    uint64_t *inA_upper = new uint64_t[dim];
    uint8_t *wrap_lower = new uint8_t[dim];
    uint8_t *wrap_upper = new uint8_t[dim];
    uint8_t *eq_upper = new uint8_t[dim];
    uint8_t *and_upper = new uint8_t[dim];
    for (int i = 0; i < dim; i++) {
        inA_lower[i] = inA[i] & mask_shift;
        inA_upper[i] = (inA[i] >> shift) & mask_upper;
        if (party == sci::BOB) {
            inA_upper[i] = (mask_upper - inA_upper[i]) & mask_upper;
        }
    }
    uint64_t *inA_prime = new uint64_t[dim];
    uint64_t quarter = 1ULL << (bw - 2);
    if (party == sci::ALICE) {
        // x - L/4
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = (inA[i] - quarter) & mask_bw;
        }
    } else {
        for (int i = 0; i < dim; i++) {
            inA_prime[i] = inA[i];
        }
    }
    uint64_t *r = new uint64_t[dim];
    uint64_t half = 1ULL << (bw - 1);
    for (int i = 0; i < dim; i++) {
        if (inA_prime[i] >= half) {
            r[i] = 1;
        } else {
            r[i] = 0;
        }
    }
    uint64_t *bit_mul = new uint64_t[dim];
    if (party == sci::ALICE) {
        uint64_t mask_shift = (1ULL << shift) - 1;
        sci::PRG128 prg;
        uint64_t *data0 = new uint64_t[dim];
        prg.random_data(data0, dim * sizeof(uint64_t));
        otpack->iknp_straight->send_cot(data0, r, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = ((1ULL << shift) - data0[i]) & mask_shift;
        }
        delete[] data0;
    } else { // party == BOB
        bool *choice = new bool[dim];
        for (int i = 0; i < dim; i++) {
            choice[i] = r[i];
        }
        uint64_t *data = new uint64_t[dim];
        otpack->iknp_straight->recv_cot(data, choice, dim, shift);
        for (int i = 0; i < dim; i++) {
            bit_mul[i] = data[i];
        }
        delete[] choice;
    }
    uint64_t offset = 1ULL << (bw - shift - 2);
    for (int i = 0; i < dim; i++) {
        if (party == sci::ALICE) {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) + offset -
                       (1ULL << (bw - shift)) * (bit_mul[i] + 1)) &
                      mask_bw;
        } else {
            outB[i] = (((inA_prime[i] >> shift) & mask_upper) -
                       (1ULL << (bw - shift)) * bit_mul[i]) &
                      mask_bw;
        }
    }
    delete[] inA_orig;
    delete[] inA_lower;
    delete[] inA_upper;
    delete[] wrap_lower;
    delete[] wrap_upper;
    delete[] eq_upper;
    delete[] and_upper;
    return;
}

void GeometricPerspectiveProtocols::ring_exp(const uint64_t *inA, uint64_t *out,
                                             size_t len, uint64_t f,
                                             uint64_t bwL, uint64_t in_f) {
uint64_t mask_bwL = (bwL == 64 ? -1 : ((1ULL << bwL) - 1));
double pow_f_input = std::pow(2.0, in_f);
double pow_f = std::pow(2.0, f);
for (size_t i = 0; i < len; ++i) {
double x0_real = static_cast<double>(inA[i]) / pow_f_input;
double exp_val = std::exp(x0_real);
uint64_t exp_fixed =
static_cast<uint64_t>(std::round(exp_val * pow_f)) & mask_bwL;
out[i] = exp_fixed;
}
}

void GeometricPerspectiveProtocols::exp(int32_t dim, uint64_t *inA,
                                        uint64_t *result, int32_t in_bw,
                                        int32_t in_f, int32_t out_bw,
                                        int32_t out_f, int32_t localexp_bw,
                                        int32_t localexp_f, int32_t locallut_bw,
                                        int32_t locallut_f) {
    uint64_t *MW = new uint64_t[dim];
    uint64_t *exp_inA = new uint64_t[dim];
    uint64_t *MW_lut = new uint64_t[4];
    uint64_t N_input = 1ULL << in_bw;
  uint64_t MW_B = N_input / 2;
  uint64_t mask_locallut_bw =
      (locallut_bw == 64 ? -1 : ((1ULL << locallut_bw) - 1));
  uint64_t mask_exp_bw = (1ULL << (2 * localexp_bw)) - 1;
    uint64_t *outB = new uint64_t[dim];

    ///////////////////////计算MW_lut///////////////////////
    for (int i = 0; i < 3; i++) {
        double pow_f_input = std::pow(2.0, in_f);
        double pow_f = std::pow(2.0, locallut_f);
    double exp_val = std::exp(-i * static_cast<double>(N_input) /
                              pow_f_input); 
    MW_lut[i] =
        static_cast<uint64_t>(std::round(exp_val * pow(2, locallut_f))) &
        mask_locallut_bw;
      }
      MW_lut[3] = 100;

      ///////////////////////计算MW///////////////////////
  mwwithB(dim, MW_B, inA, MW, in_bw, 2); // compute MW

  if (party == sci::ALICE) { // send MW to bob. Bob get whole MW. Alice get 0.
        io->send_data(MW, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
            MW[i] = 0;
        }
      } else {
        uint64_t *recv_MW = new uint64_t[dim];
        io->recv_data(recv_MW, dim * sizeof(uint64_t));
        for (int i = 0; i < dim; i++) {
          MW[i] = (MW[i] + recv_MW[i]) & ((1ULL << 2) - 1);
        }
        delete[] recv_MW;
      }


    ring_exp(inA, exp_inA, dim, localexp_f, localexp_bw, in_f);
    uint64_t *zero = new uint64_t[dim];
    for (int i = 0; i < dim; i++) {
        zero[i] = 0;
    }


    if (party == sci::ALICE) {
    cross_term(dim, exp_inA, zero, outB, localexp_bw, localexp_bw,
               localexp_bw + localexp_bw);
      } else {
    cross_term(dim, zero, exp_inA, outB, localexp_bw, localexp_bw,
               localexp_bw + localexp_bw);
      }
    

      uint8_t *msb_0 = new uint8_t[dim];
      for (int i = 0; i < dim; i++) {
        msb_0[i] = 0;
      }
    
      uint8_t *wrap_outB = new uint8_t[dim];
    

  this->aux->MSB_to_Wrap(outB, msb_0, wrap_outB, dim, 2 * localexp_bw);
      uint8_t *wrap_outB_alice_send = new uint8_t[dim];
      uint8_t *wrap_outB_bob_recv = new uint8_t[dim];
  if (party == sci::ALICE) {
            io->send_data(wrap_outB, dim * sizeof(uint8_t));
      } else {
        io->recv_data(wrap_outB_bob_recv, dim * sizeof(uint8_t));
      }
      if (party != sci::ALICE) {
        for (int i = 0; i < dim; i++) {
          if (wrap_outB_bob_recv[i] == 1 && wrap_outB[i] == 1) {
            wrap_outB_bob_recv[i] = 0;
            wrap_outB[i] = 0;
          }
        }
      }
      if (party != sci::ALICE) {
        io->send_data(wrap_outB_bob_recv, dim * sizeof(uint8_t));
      } else {
        io->recv_data(wrap_outB_alice_send, dim * sizeof(uint8_t));
        for (int i = 0; i < dim; i++) {
            wrap_outB[i] = wrap_outB_alice_send[i];
        }
      }

      uint64_t *buffer = new uint64_t[dim * 4];
      uint64_t random_numberbob = 0;
      uint64_t random_numberalice = 0;
  uint64_t mask_res_exp =
      (1ULL << (2 * localexp_bw - 2 * localexp_f + 2 + out_f)) - 1;

      if (party != sci::ALICE) {
        for (int i = 0; i < dim; ++i) {
          for (uint8_t j = 0; j < 4; ++j) {

            __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outB[i];

        __int128_t diff = (__int128_t)prod - wrap_outB[i] *
                                                 (__int128_t)MW_lut[j] *
                                                 (1ULL << (2 * localexp_bw));

            uint64_t shift = locallut_f + localexp_f * 2 - out_f; 
        uint64_t res =
            (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) &
                       mask_res_exp);
            buffer[i * 4 + j] = (res + random_numberbob) & mask_res_exp;
          }
        }
        io->send_data(buffer, dim * 4 * sizeof(uint64_t));
      } else {
        io->recv_data(buffer, dim * 4 * sizeof(uint64_t));
      }

      uint64_t *y = new uint64_t[dim];
      if (party == sci::ALICE) {
        uint64_t **spec = new uint64_t *[dim];
    
        for (int i = 0; i < dim; ++i) {
          spec[i] = new uint64_t[4];
          for (uint8_t j = 0; j < 4; ++j) {
            __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outB[i];
        __int128_t diff = (__int128_t)prod - wrap_outB[i] *
                                                 (__int128_t)MW_lut[j] *
                                                 (1ULL << (2 * localexp_bw));

        uint64_t shift = locallut_f + localexp_f * 2 -
                         out_f; 
        uint64_t res =
            (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) &
                       mask_res_exp);
        spec[i][j] =
            (res + random_numberalice + buffer[i * 4 + j]) & mask_res_exp;
          }
        }

      this->aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, 2,
                                      2 * localexp_bw - 2 * localexp_f + 2 +
                                          out_f);
        for (int i = 0; i < dim; ++i)
        delete[] spec[i];
        delete[] spec;
        } else if (party == sci::BOB) {
    this->aux->lookup_table<uint64_t>(
        nullptr, MW, y, dim, 2, 2 * localexp_bw - 2 * localexp_f + 2 + out_f);
        }

        uint64_t mask_result = (1ULL << (out_bw)) - 1;
        if (party != sci::ALICE) {
            for (int i = 0; i < dim; i++) {
                result[i] =
                  (y[i] - random_numberbob) & mask_result; 
            }
          } else {
            for (int i = 0; i < dim; i++) {
      result[i] = -random_numberalice & mask_result; 
            }
          } 

          delete[] MW;
          delete[] exp_inA;
          delete[] MW_lut;
          delete[] outB;
          delete[] zero;
          delete[] msb_0;
          delete[] wrap_outB;
          delete[] wrap_outB_alice_send;
          delete[] wrap_outB_bob_recv;
          delete[] buffer;
          delete[] y;
}

void GeometricPerspectiveProtocols::exp_nag4(
    int32_t dim, uint64_t *inA, uint64_t *result, uint64_t *MW, int32_t in_bw,
    int32_t in_f, int32_t out_bw, int32_t out_f, int32_t localexp_bw,
    int32_t localexp_f, int32_t locallut_bw, int32_t locallut_f) {
//    uint64_t *MW = new uint64_t[dim];
   uint64_t *exp_inA = new uint64_t[dim];
   uint64_t *MW_lut = new uint64_t[4];
   uint64_t N_input = 1ULL << in_bw;
  uint64_t MW_B = N_input / 2;
  uint64_t mask_locallut_bw =
      (locallut_bw == 64 ? -1 : ((1ULL << locallut_bw) - 1));
  uint64_t mask_exp_bw = (1ULL << (2 * localexp_bw)) - 1;
   uint64_t *outB = new uint64_t[dim];


   for (int i = 0; i < 3; i++) {
       double pow_f_input = std::pow(2.0, in_f);
       double pow_f = std::pow(2.0, locallut_f);
    double exp_val = std::exp(-i * static_cast<double>(N_input) / pow_f_input) *
                     std::exp(-4); 

    MW_lut[i] =
        static_cast<uint64_t>(std::round(exp_val * pow(2, locallut_f))) &
        mask_locallut_bw;
     }
     MW_lut[3] = 100;


    size_t comm_start = io->counter;

    //  mw(dim, inA, MW, MW_bw, 2); //compute MW B=N/4
     
    //  mwwithB(dim, MW_B, inA, MW, in_bw, 2); //compute MW
     size_t comm_end = io->counter;
    //  printf("MW communication cost: %zu\n", comm_end - comm_start);

  if (party == sci::ALICE) { // send MW to bob. Bob get whole MW. Alice get 0.
       io->send_data(MW, dim * sizeof(uint64_t));
       for (int i = 0; i < dim; i++) {
           MW[i] = 0;
       }
     } else {
       uint64_t *recv_MW = new uint64_t[dim];
       io->recv_data(recv_MW, dim * sizeof(uint64_t));
       for (int i = 0; i < dim; i++) {
         MW[i] = (MW[i] + recv_MW[i]) & ((1ULL << 2) - 1);
       }
       delete[] recv_MW;
     }


   ring_exp(inA, exp_inA, dim, localexp_f, localexp_bw, in_f);
   uint64_t *zero = new uint64_t[dim];
   for (int i = 0; i < dim; i++) {
       zero[i] = 0;
   }


   if (party == sci::ALICE) {
    cross_term(dim, exp_inA, zero, outB, localexp_bw, localexp_bw,
               localexp_bw + localexp_bw);
     } else {
    cross_term(dim, zero, exp_inA, outB, localexp_bw, localexp_bw,
               localexp_bw + localexp_bw);
     }
    //  printf("outB[384] = %llu\n", outB[384]);
    //  printf("MW_lut[0] = %llu\n", MW_lut[0]);
    //  printf("MW_lut[1] = %llu\n", MW_lut[1]);
    //  printf("MW_lut[2] = %llu\n", MW_lut[2]);
   

     uint8_t *msb_0 = new uint8_t[dim];
     for (int i = 0; i < dim; i++) {
       msb_0[i] = 0;
     }
   
     uint8_t *wrap_outB = new uint8_t[dim];
   

  this->aux->MSB_to_Wrap(outB, msb_0, wrap_outB, dim, 2 * localexp_bw);
     uint8_t *wrap_outB_alice_send = new uint8_t[dim];
     uint8_t *wrap_outB_bob_recv = new uint8_t[dim];
  if (party == sci::ALICE) {
           io->send_data(wrap_outB, dim * sizeof(uint8_t));
     } else {
       io->recv_data(wrap_outB_bob_recv, dim * sizeof(uint8_t));
     }
     if (party != sci::ALICE) {
       for (int i = 0; i < dim; i++) {
         if (wrap_outB_bob_recv[i] == 1 && wrap_outB[i] == 1) {
           wrap_outB_bob_recv[i] = 0;
           wrap_outB[i] = 0;
         }
       }
     }
     if (party != sci::ALICE) {
       io->send_data(wrap_outB_bob_recv, dim * sizeof(uint8_t));
     } else {
       io->recv_data(wrap_outB_alice_send, dim * sizeof(uint8_t));
       for (int i = 0; i < dim; i++) {
           wrap_outB[i] = wrap_outB_alice_send[i];
       }
     }

     uint64_t *buffer = new uint64_t[dim * 4];
     uint64_t random_numberbob = 0;
     uint64_t random_numberalice = 0;
  uint64_t mask_res_exp =
      (1ULL << (2 * localexp_bw - 2 * localexp_f + 2 + out_f)) - 1;

     if (party != sci::ALICE) {
       for (int i = 0; i < dim; ++i) {
         for (uint8_t j = 0; j < 4; ++j) {

           __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outB[i];

        __int128_t diff = (__int128_t)prod - wrap_outB[i] *
                                                 (__int128_t)MW_lut[j] *
                                                 (1ULL << (2 * localexp_bw));

           uint64_t shift = locallut_f + localexp_f * 2 - out_f; 
        uint64_t res =
            (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) &
                       mask_res_exp);
           buffer[i * 4 + j] = (res + random_numberbob) & mask_res_exp;
         }
       }
       io->send_data(buffer, dim * 4 * sizeof(uint64_t));
     } else {
       io->recv_data(buffer, dim * 4 * sizeof(uint64_t));
     }

     uint64_t *y = new uint64_t[dim];
     if (party == sci::ALICE) {
       uint64_t **spec = new uint64_t *[dim];
   
       for (int i = 0; i < dim; ++i) {
         spec[i] = new uint64_t[4];
         for (uint8_t j = 0; j < 4; ++j) {
           __uint128_t prod = (__uint128_t)MW_lut[j] * (__uint128_t)outB[i];
        __int128_t diff = (__int128_t)prod - wrap_outB[i] *
                                                 (__int128_t)MW_lut[j] *
                                                 (1ULL << (2 * localexp_bw));

        uint64_t shift = locallut_f + localexp_f * 2 -
                         out_f; 
        uint64_t res =
            (uint64_t)(((diff + (__int128_t)(1ULL << (shift - 1))) >> shift) &
                       mask_res_exp);
        spec[i][j] =
            (res + random_numberalice + buffer[i * 4 + j]) & mask_res_exp;
         }
       }

     this->aux->lookup_table<uint64_t>(spec, nullptr, nullptr, dim, 2,
                                      2 * localexp_bw - 2 * localexp_f + 2 +
                                          out_f);
       for (int i = 0; i < dim; ++i)
       delete[] spec[i];
       delete[] spec;
       } else if (party == sci::BOB) {
    this->aux->lookup_table<uint64_t>(
        nullptr, MW, y, dim, 2, 2 * localexp_bw - 2 * localexp_f + 2 + out_f);
       }
      //  printf("y[384] = %llu\n", y[384]);
      //  printf("MW[384] = %llu\n", MW[384]);
       uint64_t mask_result = (1ULL << (out_bw)) - 1;
       if (party != sci::ALICE) {
           for (int i = 0; i < dim; i++) {
               result[i] =
                 (y[i] - random_numberbob) & mask_result; 
           }
         } else {
           for (int i = 0; i < dim; i++) {
      result[i] = -random_numberalice & mask_result; 
           }
         } 

         delete[] MW;
         delete[] exp_inA;
         delete[] MW_lut;
         delete[] outB;
         delete[] zero;
         delete[] msb_0;
         delete[] wrap_outB;
         delete[] wrap_outB_alice_send;
         delete[] wrap_outB_bob_recv;
         delete[] buffer;
         delete[] y;
}

void GeometricPerspectiveProtocols::exp_softmaxx( 
  int32_t dim, uint64_t *inA, uint64_t *result, int32_t in_bw, int32_t in_f,
  int32_t localexp_bw, int32_t localexp_f, int32_t locallut_bw,
  int32_t locallut_f) // compute exp( inA), where inA >0
  {
    if (party == sci::ALICE) {
      for (int i = 0; i < dim; i++) {
        inA[i] = inA[i] - 1;
      }
    }
      uint8_t *drelu = new uint8_t[dim];
      uint64_t *drelu_input = new uint64_t[dim];
      uint64_t mask_in_bw = (in_bw == 64 ? -1 : ((1ULL << in_bw) - 1));
      uint64_t mask_in_fplus3 = (in_f + 3 == 64 ? -1 : ((1ULL << (in_f + 3)) - 1));
      uint64_t pow_f_input = std::pow(2, in_f);
      uint64_t *inA_expinput = new uint64_t[dim];
      uint64_t *exp_result = new uint64_t[dim];
      uint64_t *MW = new uint64_t[dim];

      if (party == sci::ALICE) {
      for (int i = 0; i < dim; i++) {
    drelu_input[i] = (inA[i] + 4 * pow_f_input) & mask_in_bw;
        }
      } else {
        for (int i = 0; i < dim; i++) {
    // drelu_input[i] = (inA[i] ) & mask_in_bw;
    drelu_input[i] = (inA[i] + 4 * pow_f_input) & mask_in_bw;
    // drelu_input[i] = (inA[i] + 4 * pow_f_input - 1) & mask_in_bw;
        }
      }
      

      aux->MSB(drelu_input, drelu, dim, in_bw);

if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
    drelu[i] = 1 - drelu[i];
        }
        }

        if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
    inA_expinput[i] = (inA[i] + 2 * pow_f_input) & mask_in_fplus3;
        }
      } else {
        for (int i = 0; i < dim; i++) {
    inA_expinput[i] = (inA[i] + 2 * pow_f_input) & mask_in_fplus3;
        }
      }

      size_t comm_start_exp = otpack->io->counter;

      uint64_t *MW_input = new uint64_t[dim];
      if (party == sci::ALICE) {
          for (int i = 0; i < dim; i++) {
    MW_input[i] = (inA[i] + 2 * pow_f_input) & mask_in_bw;
          }
        } else {
          for (int i = 0; i < dim; i++) {
    MW_input[i] = (inA[i] + 2 * pow_f_input) & mask_in_bw;
          }
        }
        uint64_t N_input = 1ULL << (in_f + 3);
       uint64_t MW_B = N_input / 2;
      // mw(dim, MW_input, MW, in_bw, 2);
      // mwwithB(dim, MW_B, inA_expinput, MW, in_bw, 2);
      // mwwithB(dim, MW_B, MW_input, MW, in_f + 3, 2);
      mw_conversion(dim, MW_input, MW, in_bw, in_f + 3, 2);
//         printf("MW[384] = %llu\n", MW[384]);
//       printf("inA[384] = %llu\n", inA[384]);
//       printf("pow_f_input = %llu\n", pow_f_input);
// printf("inA_expinput[384] = %llu\n", inA_expinput[384]);
// printf("MW_input[384] = %llu\n", MW_input[384]);

exp_nag4(dim, inA_expinput, exp_result, MW, in_f + 3, in_f, in_f + 3, in_f,
         localexp_bw, localexp_f, locallut_bw, locallut_f);
         
  // printf("exp_result[384] = %llu\n", exp_result[384]);
        size_t comm_end_exp = otpack->io->counter;
      std::cout << "exp nag4 comm: " << comm_end_exp - comm_start_exp << std::endl;

        this->aux->multiplexer(drelu, exp_result, result, dim, in_bw, in_bw);

      // for (int i = 0; i < 10; i++) {
      //     printf("inA[%d]: %llu\n", i, inA[i]);
      //     printf("drelu_input[%d]: %d\n", i, drelu_input[i]);
      //     printf("drelu[%d]: %d\n", i, drelu[i]);
      //     printf("MW_input[%d]: %llu\n", i, MW_input[i]);
      //     printf("MW[%d]: %llu\n", i, MW[i]);
      //     printf("exp_result[%d]: %llu\n", i, exp_result[i]);
      //     printf("result[%d]: %llu\n", i, result[i]);
      //     printf("\n");
      // }
      printf("mask_in_bw: %llu\n", mask_in_bw);
  }

void GeometricPerspectiveProtocols::exp_nagx(
    int32_t dim, uint64_t *inA, uint64_t *result, int32_t in_bw, int32_t in_f,
    int32_t localexp_bw, int32_t localexp_f, int32_t locallut_bw,
    int32_t locallut_f) // compute exp( - inA), where inA >0
    {
        uint8_t *drelu = new uint8_t[dim];
        uint64_t *drelu_input = new uint64_t[dim];
        uint64_t mask_in_bw = (in_bw == 64 ? -1 : ((1ULL << in_bw) - 1));
        uint64_t mask_in_fplus3 = (in_f + 3 == 64 ? -1 : ((1ULL << (in_f + 3)) - 1));
        uint64_t pow_f_input = std::pow(2, in_f);
        uint64_t *inA_expinput = new uint64_t[dim];
        uint64_t *exp_result = new uint64_t[dim];
        uint64_t *MW = new uint64_t[dim];


        if (party == sci::ALICE) {
        for (int i = 0; i < dim; i++) {
      drelu_input[i] = (-inA[i] + 4 * pow_f_input) & mask_in_bw;
          }
        } else {
          for (int i = 0; i < dim; i++) {
      drelu_input[i] = (-inA[i] + 4 * pow_f_input) & mask_in_bw;
          }
        }

        aux->MSB(drelu_input, drelu, dim, in_bw);

  if (party == sci::ALICE) {
          for (int i = 0; i < dim; i++) {
      drelu[i] = 1 - drelu[i];
          }
          }


          if (party == sci::ALICE) {
          for (int i = 0; i < dim; i++) {
      inA_expinput[i] = (-inA[i] + 2 * pow_f_input) & mask_in_fplus3;
          }
        } else {
          for (int i = 0; i < dim; i++) {
      inA_expinput[i] = (-inA[i] + 2 * pow_f_input) & mask_in_fplus3;
          }
        }


        size_t comm_start_exp = otpack->io->counter;

        uint64_t *MW_input = new uint64_t[dim];
        if (party == sci::ALICE) {
            for (int i = 0; i < dim; i++) {
      MW_input[i] = (-inA[i] + 2 * pow_f_input) & mask_in_bw;
            }
          } else {
            for (int i = 0; i < dim; i++) {
      MW_input[i] = (-inA[i] + 2 * pow_f_input) & mask_in_bw;
            }
          }
  
        mw(dim, MW_input, MW, in_bw, 2);
  exp_nag4(dim, inA_expinput, exp_result, MW, in_f + 3, in_f, in_f + 3, in_f,
           localexp_bw, localexp_f, locallut_bw, locallut_f);
          size_t comm_end_exp = otpack->io->counter;
        std::cout << "exp nag4 comm: " << comm_end_exp - comm_start_exp << std::endl;
        
          this->aux->multiplexer(drelu, exp_result, result, dim, in_bw, in_bw);

        // for (int i = 0; i < 1000; i++) {
        //     printf("inA[%d]: %llu\n", i, inA[i]);
        //     printf("drelu_input[%d]: %d\n", i, drelu_input[i]);
        //     printf("drelu[%d]: %d\n", i, drelu[i]);
        //     printf("exp_result[%d]: %llu\n", i, exp_result[i]);
        //     printf("result[%d]: %llu\n", i, result[i]);
        //     printf("\n");
        // }
        // printf("mask_in_bw: %llu\n", mask_in_bw);
    }

// vector_scalar multiplication bwA>=bwB
void GeometricPerspectiveProtocols::vector_bit_mul(int32_t dim, uint64_t *inA,
                                                   uint8_t ot_choice,
                                                   uint64_t *outC,
                                                   int32_t bwA) {
  uint64_t mask_A = (1ULL << bwA) - 1;
  // uint64_t mask_C = (1ULL << bwC) - 1;

  int32_t values_per_block = 128 / bwA;
  int32_t total_blocks = (dim + values_per_block - 1) / values_per_block;

  sci::block128 *encoded_blocks = nullptr;
  sci::block128 *encoded_blocks_f = nullptr;

  unsigned char key_bytes[16];
  unsigned char key_bytes_f[16];

  if (party == sci::ALICE) {
    unsigned char alice_key_0[16] = {0, 9, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 9, 9};
    unsigned char alice_key_1[16] = {1, 1, 1, 'a', 0, 0, 1, 0,
                                     0, 0, 0, 5,   0, 0, 4, 0};
    memcpy(key_bytes, alice_key_0, 16);
    memcpy(key_bytes_f, alice_key_1, 16);

        uint64_t *A_random = new uint64_t[dim];
        uint64_t *random_number = new uint64_t[dim];

    srand(1);
        for (int i = 0; i < dim; i++) {
            random_number[i] = (rand() % (mask_A + 1));
      outC[i] = (-random_number[i]) & mask_A;
        }
        for (int i = 0; i < dim; i++) {
            A_random[i] = (inA[i] + random_number[i]) & mask_A;
        }
    for (int i = 0; i < 100; i++) {
      printf("random_number[%d]: %llu\n", i, random_number[i]);
      printf("A_random[%d]: %llu\n", i, A_random[i]);
    }

    sci::AESNI_KEY encrypt_key, encrypt_key_f;
    sci::AESNI_set_encrypt_key(&encrypt_key, key_bytes, 16);
    sci::AESNI_set_encrypt_key(&encrypt_key_f, key_bytes_f, 16);

    encoded_blocks = encode_A_random_to_blocks(A_random, dim, bwA);
    encoded_blocks_f = encode_A_random_to_blocks(random_number, dim, bwA);

    uint64_t block_data[2];
    _mm_storeu_si128((__m128i *)block_data, encoded_blocks[0]);
    printf("encoded_blocks[0]: 0x%016llx%016llx (contains %d values)\n",
           block_data[1], block_data[0], values_per_block);

    sci::AESNI_ecb_encrypt_blks(encoded_blocks, total_blocks, &encrypt_key);
    sci::AESNI_ecb_encrypt_blks(encoded_blocks_f, total_blocks, &encrypt_key_f);

    _mm_storeu_si128((__m128i *)block_data, encoded_blocks[0]);
    printf("encrypt encoded_blocks[0]: 0x%016llx%016llx\n", block_data[1],
           block_data[0]);

    delete[] A_random;
    delete[] random_number;
  }

  sci::block128 *received_encrypt_blocks = new sci::block128[total_blocks];
  sci::block128 *received_encrypt_blocks_f = new sci::block128[total_blocks];

  if (party == sci::ALICE) {
    io->send_data(encoded_blocks, total_blocks * sizeof(sci::block128));
    io->send_data(encoded_blocks_f,
                          total_blocks * sizeof(sci::block128));
  } else {
    io->recv_data(received_encrypt_blocks,
                          total_blocks * sizeof(sci::block128));
    io->recv_data(received_encrypt_blocks_f,
                          total_blocks * sizeof(sci::block128));
  }

  const int num_ots = 8; 
  sci::block128 key_block_0, key_block_1;
  sci::block128 received_key_block[num_ots];
  bool ot_choice_bool[num_ots];

  if (party == sci::ALICE) {
    key_block_0 = _mm_loadu_si128((__m128i *)key_bytes_f);
    key_block_1 = _mm_loadu_si128((__m128i *)key_bytes);

    sci::block128 ot_data0[num_ots];
    sci::block128 ot_data1[num_ots];

    for (int i = 0; i < num_ots; i++) {
      ot_data0[i] = key_block_0; 
      ot_data1[i] = key_block_1; 
    }

    otpack->iknp_straight->send(ot_data0, ot_data1, num_ots);
    printf("Alice sent %d OTs with 128-bit keys\n", num_ots);
    for (int i = 0; i < num_ots; i++) {
      uint64_t data0_parts[2], data1_parts[2];
      _mm_storeu_si128((__m128i *)data0_parts, ot_data0[i]);
      _mm_storeu_si128((__m128i *)data1_parts, ot_data1[i]);
      printf("ot_data0[%d]: 0x%016llx%016llx\n", i, data0_parts[1],
             data0_parts[0]);
      printf("ot_data1[%d]: 0x%016llx%016llx\n", i, data1_parts[1],
             data1_parts[0]);
    }
  } else {
    for (int i = 0; i < num_ots; i++) {
      ot_choice_bool[i] = (ot_choice == 1); 
    }

    otpack->iknp_straight->recv(received_key_block, ot_choice_bool, num_ots);
    printf("Bob received %d OTs, using first result\n", num_ots);
    for (int i = 0; i < num_ots; i++) {
      uint64_t received_parts[2];
      _mm_storeu_si128((__m128i *)received_parts, received_key_block[i]);
      printf("received_key_block[%d]: 0x%016llx%016llx\n", i, received_parts[1],
             received_parts[0]);
    }
  }

  if (party == sci::BOB) {
    unsigned char received_key[16];
    _mm_storeu_si128((__m128i *)received_key, received_key_block[0]);

    sci::AESNI_KEY received_decrypt_key;
    sci::AESNI_set_decrypt_key(&received_decrypt_key, received_key, 16);

    if (ot_choice == 1) {
      sci::AESNI_ecb_decrypt_blks(received_encrypt_blocks, total_blocks,
                                  &received_decrypt_key);
      uint64_t *decoded_data =
          decode_blocks_to_A_random(received_encrypt_blocks, dim, bwA);
      memcpy(outC, decoded_data, dim * sizeof(uint64_t));
      delete[] decoded_data;
    } else if (ot_choice == 0) {
      sci::AESNI_ecb_decrypt_blks(received_encrypt_blocks_f, total_blocks,
                                  &received_decrypt_key);
      uint64_t *decoded_data =
          decode_blocks_to_A_random(received_encrypt_blocks_f, dim, bwA);
      memcpy(outC, decoded_data, dim * sizeof(uint64_t));
      delete[] decoded_data;
    }
    for (int i = 0; i < 100; i++) {
      printf("outC[%d]: %llu \n", i, outC[i]);
    }
    // printf("outC[0]: %llu \n",outC[0]);
  }

  if (party == sci::ALICE) {
    delete[] encoded_blocks;
    delete[] encoded_blocks_f;
  }
  delete[] received_encrypt_blocks;
  delete[] received_encrypt_blocks_f;
}

sci::block128 *GeometricPerspectiveProtocols::encode_matrix_to_blocks(
    const uint64_t *matrix, int32_t m, int32_t n, int32_t bwA) {

  int32_t values_per_block = 128 / bwA;
  int32_t total_blocks_per_row = (n + values_per_block - 1) / values_per_block;
  int32_t total_blocks = m * total_blocks_per_row;
        
  sci::block128 *block_array = new sci::block128[total_blocks];

  for (int32_t row = 0; row < m; row++) {
    for (int32_t block_idx = 0; block_idx < total_blocks_per_row; block_idx++) {
      uint64_t block_data_low = 0;
      uint64_t block_data_high = 0;

      for (int32_t i = 0; i < values_per_block; i++) {
        int32_t col = block_idx * values_per_block + i;

        if (col < n) {
          uint64_t value = matrix[row * n + col]; 
          int32_t bit_offset = i * bwA;

          if (bit_offset < 64) {
            if (bit_offset + bwA <= 64) {
              block_data_low |= (value << bit_offset);
            } else {
              int32_t bits_in_low = 64 - bit_offset;
              int32_t bits_in_high = bwA - bits_in_low;
              block_data_low |= (value << bit_offset);
              block_data_high |= (value >> bits_in_low);
            }
          } else {
            block_data_high |= (value << (bit_offset - 64));
          }
        }
      }

      int32_t block_global_idx = row * total_blocks_per_row + block_idx;
      block_array[block_global_idx] =
          sci::makeBlock128(block_data_high, block_data_low);
    }
  }

  return block_array;
}

uint64_t *GeometricPerspectiveProtocols::decode_blocks_to_matrix(
    const sci::block128 *block_array, int32_t m, int32_t n, int32_t bwA) {

  int32_t values_per_block = 128 / bwA;
  int32_t total_blocks_per_row = (n + values_per_block - 1) / values_per_block;

  uint64_t *matrix = new uint64_t[m * n];

  for (int32_t row = 0; row < m; row++) {
    for (int32_t block_idx = 0; block_idx < total_blocks_per_row; block_idx++) {
      int32_t block_global_idx = row * total_blocks_per_row + block_idx;

      uint64_t block_data[2];
      _mm_storeu_si128((__m128i *)block_data, block_array[block_global_idx]);
      uint64_t block_data_low = block_data[0];
      uint64_t block_data_high = block_data[1];

      for (int32_t i = 0; i < values_per_block; i++) {
        int32_t col = block_idx * values_per_block + i;

        if (col < n) {
          int32_t bit_offset = i * bwA;
          uint64_t value = 0;

          if (bit_offset < 64) {
            if (bit_offset + bwA <= 64) {
              uint64_t mask =
                  (bwA == 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << bwA) - 1);
              value = (block_data_low >> bit_offset) & mask;
                        } else {
              int32_t bits_in_low = 64 - bit_offset;
              int32_t bits_in_high = bwA - bits_in_low;
              uint64_t low_part = (block_data_low >> bit_offset);
              uint64_t high_mask = (bits_in_high == 64)
                                       ? 0xFFFFFFFFFFFFFFFFULL
                                       : ((1ULL << bits_in_high) - 1);
              uint64_t high_part = (block_data_high & high_mask) << bits_in_low;
              value = low_part | high_part;
                        }
                    } else {
            uint64_t mask =
                (bwA == 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << bwA) - 1);
            value = (block_data_high >> (bit_offset - 64)) & mask;
          }

          matrix[row * n + col] = value; 
        }
      }
    }
  }

  return matrix;
}

sci::block128 *GeometricPerspectiveProtocols::encode_A_random_to_blocks(
    const uint64_t *A_random, int32_t dim, int32_t bwA) {

  int32_t values_per_block = 128 / bwA; 
  int32_t total_blocks =
      (dim + values_per_block - 1) / values_per_block; 

  sci::block128 *block_array = new sci::block128[total_blocks];

  for (int32_t block_idx = 0; block_idx < total_blocks; block_idx++) {
    uint64_t block_data_low = 0;
    uint64_t block_data_high = 0;

    for (int32_t i = 0; i < values_per_block; i++) {
      int32_t value_idx = block_idx * values_per_block + i;

      if (value_idx < dim) {
        uint64_t value = A_random[value_idx];
        int32_t bit_offset = i * bwA;

        if (bit_offset < 64) {
          if (bit_offset + bwA <= 64) {
            block_data_low |= (value << bit_offset);
          } else {
            int32_t bits_in_low = 64 - bit_offset;
            int32_t bits_in_high = bwA - bits_in_low;
            block_data_low |= (value << bit_offset);
            block_data_high |= (value >> bits_in_low);
          }
        } else {
          block_data_high |= (value << (bit_offset - 64));
        }
      }
    }

    block_array[block_idx] = sci::makeBlock128(block_data_high, block_data_low);
  }
        
        return block_array;
    }
    
uint64_t *GeometricPerspectiveProtocols::decode_blocks_to_A_random(
    const sci::block128 *block_array, int32_t dim, int32_t bwA) {

  int32_t values_per_block = 128 / bwA; 
  int32_t total_blocks =
      (dim + values_per_block - 1) / values_per_block; 
        
  uint64_t *A_random = new uint64_t[dim];

  for (int32_t block_idx = 0; block_idx < total_blocks; block_idx++) {
    uint64_t block_data[2];
    _mm_storeu_si128((__m128i *)block_data, block_array[block_idx]);
    uint64_t block_data_low = block_data[0];
    uint64_t block_data_high = block_data[1];

    for (int32_t i = 0; i < values_per_block; i++) {
      int32_t value_idx = block_idx * values_per_block + i;

      if (value_idx < dim) {
        int32_t bit_offset = i * bwA;
        uint64_t value = 0;

        if (bit_offset < 64) {
          if (bit_offset + bwA <= 64) {
            uint64_t mask =
                (bwA == 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << bwA) - 1);
            value = (block_data_low >> bit_offset) & mask;
          } else {
            int32_t bits_in_low = 64 - bit_offset;
            int32_t bits_in_high = bwA - bits_in_low;
            uint64_t low_part = (block_data_low >> bit_offset);
            uint64_t high_mask = (bits_in_high == 64)
                                     ? 0xFFFFFFFFFFFFFFFFULL
                                     : ((1ULL << bits_in_high) - 1);
            uint64_t high_part = (block_data_high & high_mask) << bits_in_low;
            value = low_part | high_part;
          }
        } else {
          uint64_t mask =
              (bwA == 64) ? 0xFFFFFFFFFFFFFFFFULL : ((1ULL << bwA) - 1);
          value = (block_data_high >> (bit_offset - 64)) & mask;
        }

        A_random[value_idx] = value;
      }
    }
  }

  return A_random;
}

void GeometricPerspectiveProtocols::matrix_vector_crossterm(
    int32_t m, int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC,
    int32_t bwA, int32_t bwB) {
  int32_t bwC = bwA + bwB;
  uint64_t maskA = (1ULL << bwA) - 1;
  uint64_t maskC = (1ULL << bwC) - 1; 

  uint64_t *random = nullptr;
  uint64_t *inA_random = nullptr;
  sci::block128 *random_blocks = nullptr;
  sci::block128 *inA_random_blocks = nullptr;
  sci::AESNI_KEY *key_random_array = nullptr;
  sci::AESNI_KEY *key_inA_random_array = nullptr;
  unsigned char(*key_bytes_random)[16] = nullptr;
  unsigned char(*key_bytes_inA_random)[16] = nullptr;

  for (int i = 0; i < m * n; i++) {
    outC[i] = 0;
  }

  for (int bit_inB = 0; bit_inB < bwB; bit_inB++) {
    int32_t target_bw = bwC - bit_inB;
    uint64_t target_mask =
        (target_bw >= 64) ? -1ULL : ((1ULL << target_bw) - 1);

    int32_t values_per_block = 128 / target_bw; 
    int32_t total_blocks_per_row =
        (n + values_per_block - 1) / values_per_block; 
    int32_t total_blocks = m * total_blocks_per_row; 

    if (party == sci::ALICE) {
      if (bit_inB == 0) {
        random = new uint64_t[m * n];
        inA_random = new uint64_t[m * n];

        key_random_array = new sci::AESNI_KEY[m];
        key_inA_random_array = new sci::AESNI_KEY[m];
        key_bytes_random = new unsigned char[m][16];
        key_bytes_inA_random = new unsigned char[m][16];

        sci::PRG128 prg;

        for (int row = 0; row < m; row++) {
          prg.random_data(key_bytes_random[row], 16);
          prg.random_data(key_bytes_inA_random[row], 16);

          sci::AESNI_set_encrypt_key(&key_random_array[row],
                                     key_bytes_random[row], 16);
          sci::AESNI_set_encrypt_key(&key_inA_random_array[row],
                                     key_bytes_inA_random[row], 16);
        }
      }

      srand(1234 + bit_inB * 100 +
            party * 10); 
      for (int i = 0; i < m * n; i++) {
        random[i] = (rand() % (target_mask + 1));

        uint64_t extended_inA = inA[i] & maskA;
        inA_random[i] = (extended_inA + random[i]) & target_mask;
      }


      if (bit_inB > 0) {
        delete[] random_blocks;
        delete[] inA_random_blocks;
      }
      random_blocks = encode_matrix_to_blocks(random, m, n, target_bw);
      inA_random_blocks = encode_matrix_to_blocks(inA_random, m, n, target_bw);

      for (int row = 0; row < m; row++) {
        int row_start_idx = row * total_blocks_per_row;

        sci::AESNI_ecb_encrypt_blks(&random_blocks[row_start_idx],
                                    total_blocks_per_row,
                                    &key_random_array[row]);

        sci::AESNI_ecb_encrypt_blks(&inA_random_blocks[row_start_idx],
                                    total_blocks_per_row,
                                    &key_inA_random_array[row]);
      }

      printf("Generated %d key pairs and encrypted %d total blocks for "
             "bit_inB=%d\n",
             m, total_blocks, bit_inB);
    }

    sci::block128 *received_enc_random =
        new sci::block128[total_blocks]; 
    sci::block128 *received_enc_inA_random = new sci::block128[total_blocks];
    if (party == sci::ALICE) {
      io->send_data(random_blocks,
                            total_blocks * sizeof(sci::block128));
      io->send_data(inA_random_blocks,
                            total_blocks * sizeof(sci::block128));
                        } else {
      io->recv_data(received_enc_random,
                            total_blocks * sizeof(sci::block128));
      io->recv_data(received_enc_inA_random,
                            total_blocks * sizeof(sci::block128));
    }

    sci::block128 received_key_block[m];
    bool ot_choice_bool[m];
    if (party == sci::ALICE) {
      sci::block128 ot_data0[m];
      sci::block128 ot_data1[m];

      for (int i = 0; i < m; i++) {
        ot_data0[i] =
            _mm_loadu_si128((__m128i *)key_bytes_random[i]); 
        ot_data1[i] =
            _mm_loadu_si128((__m128i *)key_bytes_inA_random[i]); 
      }
      otpack->iknp_straight->send(ot_data0, ot_data1, m);
    } else {
      for (int i = 0; i < m; i++) {
        ot_choice_bool[i] =
            ((inB[i] >> bit_inB) & 1); 
      }
      otpack->iknp_straight->recv(received_key_block, ot_choice_bool, m);
    }

    sci::block128 *decrypted_blocks = new sci::block128[total_blocks];
    sci::AESNI_KEY received_decrypt_key[m];
    uint64_t *decoded_data = nullptr; 

    if (party == sci::BOB) {
      for (int row = 0; row < m; row++) {
        unsigned char key_bytes[16];
        _mm_storeu_si128((__m128i *)key_bytes, received_key_block[row]);
        sci::AESNI_set_decrypt_key(&received_decrypt_key[row], key_bytes, 16);
      }

      for (int row = 0; row < m; row++) {
        int row_start_idx = row * total_blocks_per_row;

        if (ot_choice_bool[row] == 0) {
          for (int j = 0; j < total_blocks_per_row; j++) {
            decrypted_blocks[row_start_idx + j] =
                received_enc_random[row_start_idx + j];
          }
          sci::AESNI_ecb_decrypt_blks(&decrypted_blocks[row_start_idx],
                                      total_blocks_per_row,
                                      &received_decrypt_key[row]);
        } else {
          for (int j = 0; j < total_blocks_per_row; j++) {
            decrypted_blocks[row_start_idx + j] =
                received_enc_inA_random[row_start_idx + j];
          }
          sci::AESNI_ecb_decrypt_blks(&decrypted_blocks[row_start_idx],
                                      total_blocks_per_row,
                                      &received_decrypt_key[row]);
        }
      }

      decoded_data = decode_blocks_to_matrix(decrypted_blocks, m, n, target_bw);
    }
    if (party == sci::ALICE) {
      uint64_t val0 = 0;
      for (int i = 0; i < m * n; i++) {
        uint64_t val = (-random[i]) & target_mask;

        if (i == 0) {
          val0 = val;
        }
        val = (val * (1ULL << bit_inB)) & maskC;
        outC[i] += val;
        outC[i] &= maskC; 
      }

      printf("ALICE //////////////////////////////////////\n");
    } else {
      uint64_t val0 = 0;
      for (int i = 0; i < m * n; i++) {
        uint64_t val = decoded_data[i] & target_mask;

        if (i == 0) {
          val0 = val;
        }
        val = (val * (1ULL << bit_inB)) & maskC;

        outC[i] += val;
        outC[i] &= maskC; 
      }

    //   printf("BOB bit_inB: %d\n", bit_inB);
    //   printf("bob ot_choice_bool[0]: %d\n", ot_choice_bool[0]);
    //   unsigned char received_key_bytes[16];
    //   _mm_storeu_si128((__m128i *)received_key_bytes, received_key_block[0]);
    //   printf("bob received_key: 0x%02x%02x%02x%02x\n", received_key_bytes[0],
    //          received_key_bytes[1], received_key_bytes[2],
    //          received_key_bytes[3]);
    //   
    //   if (decoded_data != nullptr) {
    //     printf("bob decoded_data[0]: %lu\n", decoded_data[0]);
    //   }
    //   printf("bob val0: %lu\n", val0);
    //   // printf("bob (val * (1ULL << bit_inB)): %lu\n", (val * (1ULL <<
    //   // bit_inB)));
    //   printf("bob outC[0]: %lu\n", outC[0]);
      printf("BOB //////////////////////////////////////\n");
    }

    if (party == sci::BOB && decoded_data != nullptr) {
      delete[] decoded_data;
    }
    delete[] decrypted_blocks;

  }

  if (party == sci::ALICE) {
    if (random != nullptr)
      delete[] random;
    if (inA_random != nullptr)
      delete[] inA_random;
    if (random_blocks != nullptr)
      delete[] random_blocks;
    if (inA_random_blocks != nullptr)
      delete[] inA_random_blocks;
    if (key_random_array != nullptr)
      delete[] key_random_array;
    if (key_inA_random_array != nullptr)
      delete[] key_inA_random_array;
    if (key_bytes_random != nullptr)
      delete[] key_bytes_random;
    if (key_bytes_inA_random != nullptr)
      delete[] key_bytes_inA_random;
  }
}

void GeometricPerspectiveProtocols::matrix_vector_mul(
    int32_t m, int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC,
    int32_t bwA, int32_t bwB, int32_t bwC) {}

void GeometricPerspectiveProtocols::matrix_vector_crossterm_reverse(
    int32_t m, int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC,
    int32_t bwA, int32_t bwB) {
  int32_t bwC = bwA + bwB;
  uint64_t maskA = (1ULL << bwA) - 1;
  uint64_t maskC = (1ULL << bwC) - 1;

  uint64_t *random = nullptr;
  uint64_t *inA_random = nullptr;
  sci::block128 *random_blocks = nullptr;
  sci::block128 *inA_random_blocks = nullptr;
  sci::AESNI_KEY *key_random_array = nullptr;
  sci::AESNI_KEY *key_inA_random_array = nullptr;
  unsigned char(*key_bytes_random)[16] = nullptr;
  unsigned char(*key_bytes_inA_random)[16] = nullptr;

  for (int i = 0; i < m * n; i++) {
    outC[i] = 0;
  }

  for (int bit_inB = 0; bit_inB < bwB; bit_inB++) {
    int32_t target_bw = bwC - bit_inB;
    uint64_t target_mask =
        (target_bw >= 64) ? -1ULL : ((1ULL << target_bw) - 1);

    int32_t values_per_block = 128 / target_bw;
    int32_t total_blocks_per_row =
        (n + values_per_block - 1) / values_per_block;
    int32_t total_blocks = m * total_blocks_per_row;

    if (party == sci::BOB) {
      if (bit_inB == 0) {
        random = new uint64_t[m * n];
        inA_random = new uint64_t[m * n];

        key_random_array = new sci::AESNI_KEY[m];
        key_inA_random_array = new sci::AESNI_KEY[m];
        key_bytes_random = new unsigned char[m][16];
        key_bytes_inA_random = new unsigned char[m][16];

        sci::PRG128 prg;

        for (int row = 0; row < m; row++) {
          prg.random_data(key_bytes_random[row], 16);
          prg.random_data(key_bytes_inA_random[row], 16);

          sci::AESNI_set_encrypt_key(&key_random_array[row],
                                     key_bytes_random[row], 16);
          sci::AESNI_set_encrypt_key(&key_inA_random_array[row],
                                     key_bytes_inA_random[row], 16);
        }
      }

      srand(1234 + bit_inB * 100 + party * 10);
      for (int i = 0; i < m * n; i++) {
        random[i] = (rand() % (target_mask + 1));

        uint64_t extended_inA = inA[i] & maskA;
        inA_random[i] = (extended_inA + random[i]) & target_mask;
      }

      if (bit_inB > 0) {
        delete[] random_blocks;
        delete[] inA_random_blocks;
      }
      random_blocks = encode_matrix_to_blocks(random, m, n, target_bw);
      inA_random_blocks = encode_matrix_to_blocks(inA_random, m, n, target_bw);

      for (int row = 0; row < m; row++) {
        int row_start_idx = row * total_blocks_per_row;

        sci::AESNI_ecb_encrypt_blks(&random_blocks[row_start_idx],
                                    total_blocks_per_row,
                                    &key_random_array[row]);

        sci::AESNI_ecb_encrypt_blks(&inA_random_blocks[row_start_idx],
                                    total_blocks_per_row,
                                    &key_inA_random_array[row]);
      }

      printf("BOB Generated %d key pairs and encrypted %d total blocks for "
             "bit_inB=%d\n",
             m, total_blocks, bit_inB);
    }

    sci::block128 *received_enc_random =
        new sci::block128[total_blocks]; 
    sci::block128 *received_enc_inA_random = new sci::block128[total_blocks];
    if (party == sci::BOB) {
      io->send_data(random_blocks,
                            total_blocks * sizeof(sci::block128));
      io->send_data(inA_random_blocks,
                            total_blocks * sizeof(sci::block128));
    } else {
      io->recv_data(received_enc_random,
                            total_blocks * sizeof(sci::block128));
      io->recv_data(received_enc_inA_random,
                            total_blocks * sizeof(sci::block128));
    }

    sci::block128 received_key_block[m];
    bool ot_choice_bool[m];
    if (party == sci::BOB) {
      sci::block128 ot_data0[m];
      sci::block128 ot_data1[m];

      for (int i = 0; i < m; i++) {
        ot_data0[i] =
            _mm_loadu_si128((__m128i *)key_bytes_random[i]); 
        ot_data1[i] =
            _mm_loadu_si128((__m128i *)key_bytes_inA_random[i]); 
      }
      otpack->iknp_reversed->send(ot_data0, ot_data1, m);
    } else {
      for (int i = 0; i < m; i++) {
        ot_choice_bool[i] =
            ((inB[i] >> bit_inB) & 1); 
      }
      otpack->iknp_reversed->recv(received_key_block, ot_choice_bool, m);
    }

    sci::block128 *decrypted_blocks = new sci::block128[total_blocks];
    sci::AESNI_KEY received_decrypt_key[m];
    uint64_t *decoded_data = nullptr; 

    if (party == sci::ALICE) {
      for (int row = 0; row < m; row++) {
        unsigned char key_bytes[16];
        _mm_storeu_si128((__m128i *)key_bytes, received_key_block[row]);
        sci::AESNI_set_decrypt_key(&received_decrypt_key[row], key_bytes, 16);
      }

      for (int row = 0; row < m; row++) {
        int row_start_idx = row * total_blocks_per_row;

        if (ot_choice_bool[row] == 0) {
          for (int j = 0; j < total_blocks_per_row; j++) {
            decrypted_blocks[row_start_idx + j] =
                received_enc_random[row_start_idx + j];
          }
          sci::AESNI_ecb_decrypt_blks(&decrypted_blocks[row_start_idx],
                                      total_blocks_per_row,
                                      &received_decrypt_key[row]);
        } else {
          for (int j = 0; j < total_blocks_per_row; j++) {
            decrypted_blocks[row_start_idx + j] =
                received_enc_inA_random[row_start_idx + j];
          }
          sci::AESNI_ecb_decrypt_blks(&decrypted_blocks[row_start_idx],
                                      total_blocks_per_row,
                                      &received_decrypt_key[row]);
        }
      }

      decoded_data = decode_blocks_to_matrix(decrypted_blocks, m, n, target_bw);
    }

    if (party == sci::BOB) {
      uint64_t val0 = 0;
      for (int i = 0; i < m * n; i++) {
        uint64_t val = (-random[i]) & target_mask;

        if (i == 0) {
          val0 = val;
        }
        val = (val * (1ULL << bit_inB)) & maskC;
        outC[i] += val;
        outC[i] &= maskC; 
      }
    //   printf("BOB_REVERSE bit_inB: %d\n", bit_inB);
    //   printf("bob random[0]: %lu\n", random[0]);
    //   printf("bob inA_random[0]: %lu\n", inA_random[0]);
    //   printf("bob val0: %lu\n", val0);
    //   printf("bob outC[0]: %lu\n", outC[0]);
    //   printf("BOB_REVERSE //////////////////////////////////////\n");
    } else {
      uint64_t val0 = 0;
      for (int i = 0; i < m * n; i++) {
        uint64_t val = decoded_data[i] & target_mask;

        if (i == 0) {
          val0 = val;
        }
        val = (val * (1ULL << bit_inB)) & maskC;

        outC[i] += val;
        outC[i] &= maskC; 
      }

    //   printf("ALICE_REVERSE bit_inB: %d\n", bit_inB);
    //   printf("alice ot_choice_bool[0]: %d\n", ot_choice_bool[0]);
    //   unsigned char received_key_bytes[16];
    //   _mm_storeu_si128((__m128i *)received_key_bytes, received_key_block[0]);
    //   printf("alice received_key: 0x%02x%02x%02x%02x\n", received_key_bytes[0],
    //          received_key_bytes[1], received_key_bytes[2],
    //          received_key_bytes[3]);
    //   if (decoded_data != nullptr) {
    //     printf("alice decoded_data[0]: %lu\n", decoded_data[0]);
    //   }
    //   printf("alice val0: %lu\n", val0);
    //   printf("alice outC[0]: %lu\n", outC[0]);
    //   printf("ALICE_REVERSE //////////////////////////////////////\n");
    }

    if (party == sci::ALICE && decoded_data != nullptr) {
      delete[] decoded_data;
    }
    delete[] decrypted_blocks;

  }

  if (party == sci::BOB) {
    if (random != nullptr)
      delete[] random;
    if (inA_random != nullptr)
      delete[] inA_random;
    if (random_blocks != nullptr)
      delete[] random_blocks;
    if (inA_random_blocks != nullptr)
      delete[] inA_random_blocks;
    if (key_random_array != nullptr)
      delete[] key_random_array;
    if (key_inA_random_array != nullptr)
      delete[] key_inA_random_array;
    if (key_bytes_random != nullptr)
      delete[] key_bytes_random;
    if (key_bytes_inA_random != nullptr)
      delete[] key_bytes_inA_random;
  }
}

void GeometricPerspectiveProtocols::matrix_vector_unsigned_mul(
    int32_t m, int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA,
    int32_t bwB, int32_t bwC) {
  auto *inB_expanded = new uint64_t[m*n];
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      inB_expanded[i * n + j] = inB[i];
    }
  }
  
  auto *c = new uint64_t[m*n]();
  auto *d = new uint64_t[m*n];
  for (int i = 0; i < m*n; i++) {
    c[i] = 0;
    d[i] = 0;
  }
    matrix_vector_crossterm(m, n, inA, inB, d, bwA, bwB);
    matrix_vector_crossterm_reverse(m, n, inA, inB, c, bwA, bwB);
  auto *wx = new uint8_t[m*n];
  auto *wy = new uint8_t[m*n];
  uint8_t *msb0 = new uint8_t[m*n]();
  for (int i = 0; i < m*n; i++) {
    msb0[i] = 0;
  }
  aux->MSB_to_Wrap(inA, msb0, wx, m*n, bwA);
  aux->MSB_to_Wrap(inB_expanded, msb0, wy, m*n, bwB);  
  auto *h = new uint64_t[m*n];
  auto *g = new uint64_t[m*n];
  aux->multiplexer(wx, inB_expanded, h, m*n, bwB, bwB);  
  aux->multiplexer(wy, inA, g, m*n, bwA, bwA);
  uint64_t mask = (1ULL << bwC) - 1;
  for (int i = 0; i < m*n; i++) {
    outC[i] = (inA[i] * inB_expanded[i] + c[i] + d[i] - (g[i] * (1ULL << bwB)) -
               (h[i] * (1ULL << bwA))) &
              mask;
  }
  delete[] inB_expanded;  
  delete[] c;
  delete[] d;
  delete[] wx;
  delete[] wy;
  delete[] msb0;
  delete[] h;
  delete[] g;
}


// void GeometricPerspectiveProtocols::division(uint64_t *input, uint64_t *output, int32_t dim, uint64_t divisor, uint_32_bw) {
//     int l_d = ceil(log2(divisor * 1.0));
//     int mod_mask = (1ull << bw) - 1;
//     int div_mask = (1ull << l_d) - 1;
//     uint64_t *mw = new uint64_t[dim];
    
//     auto relu = new ReLURingProtocol<sci::NetIO, int64_t>(party, 0, this->io, l_d + 1, t, this->otpack);

//     mwwithB(dim, 1ULL << (bw - 1), input, mw, bw, 2);

//     if (party == sci::BOB)
//     {
//         std::unique_ptr<uint64_t*[]> t_d;
// 	for (uint64_t i = 0; i < dims; i++) {
// 	    t_d[i] = new uint64_t[4];
// 	    t_d[0] = (input[i]) / divisior;
// 	    t_d[1] = (input[i] - mod) / divisior;
// 	    t_d[2] = (input[i] - 2 * mod) / divisor;
// 	    t_d[3] = 0;
// 	}
// 	uint_64_t *a = new uint64_t[dim];
// 	aux->lookup_table(t_d.get(), mw, a, dim, bw, bw);

//         std::unique_ptr<uint64_t*[]> t_e;
//         for (uint64_t i = 0; i < dims; i++) {
//             t_e[i] = new uint64_t[4];
//             t_e[0] = (input[i]) & div_mask;
//             t_e[1] = (input[i] - mod) & div_mask;
//             t_e[2] = (input[i] - 2 * mod) & div_mask;
//             t_e[3] = 0;
//         }
	
//         uint_64_t *b = new uint64_t[dim];
//         aux->lookup_table(t_e.get(), mw, b, dim, l_d+1, l_d+1);

// 	int64_t *temp = new int64_t[dim];
// 	for (int i = 0; i < dim; i++) {
// 	    temp[i] = (input[i] % divisor) + b[i] - divisor;
// 	}

// 	int64_t *relu_res = new int64_t[dim];
// 	uint8_t *drelu_res = new uint8_t[dim];
// 	relu->relu(relu_resm wm dim, drelu_res);
// 	uint64_t *arith_res = new uint64_t[dim];
// 	aux->B2A(drelu_res, arith_res, dim, bw);
//     } else {
// 	std::unique_ptr<uint64_t*[]> t_d;
//         for (uint64_t i = 0; i < dims; i++) {
//             t_d[i] = new uint64_t[4];
//             t_d[0] = 0;//(input[i]) / divisior;
//             t_d[1] = 0;//(input[i] - mod) / divisior;
//             t_d[2] = 0;//(input[i] - 2 * mod) / divisor;
//             t_d[3] = 0;
//         }
//         uint_64_t *a = new uint64_t[dim];
//         aux->lookup_table(t_d.get(), mw, a, dim, bw, bw);
        
//         std::unique_ptr<uint64_t*[]> t_e;
//         for (uint64_t i = 0; i < dims; i++) {
//             t_e[i] = new uint64_t[4];
//             t_e[0] = 0;//(input[i]) & div_mask;
//             t_e[1] = 0;//(input[i] - mod) & div_mask;
//             t_e[2] = 0;//(input[i] - 2 * mod) & div_mask;
//             t_e[3] = 0;
//         }
        
//         uint_64_t *b = new uint64_t[dim];
//         aux->lookup_table(t_e.get(), mw, b, dim, l_d+1, l_d+1);

//         int64_t *temp = new int64_t[dim];
//         for (int i = 0; i < dim; i++) {
//             temp[i] = b[i];
//         }
// 	int64_t *relu_res = new int64_t[dim];
//         uint8_t *drelu_res = new uint8_t[dim];
//         relu->relu(relu_resm wm dim, drelu_res);
//         uint64_t *arith_res = new uint64_t[dim];
//         aux->B2A(drelu_res, arith_res, dim, bw);

//     }    


// }
