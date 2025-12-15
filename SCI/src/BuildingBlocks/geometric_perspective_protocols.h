// anonymous authors

#ifndef GEOMTRIC_PERSPECTIVE_PROTOCOLS_H
#define GEOMTRIC_PERSPECTIVE_PROTOCOLS_H

#include "BuildingBlocks/aux-protocols.h"
#include "Millionaire/equality.h"
#include "Millionaire/millionaire.h"
#include "Millionaire/millionaire_with_equality.h"
#include "NonLinear/relu-ring.h"
#include "OT/emp-ot.h"
#include "utils/aes-ni.h"

class GeometricPerspectiveProtocols {

    public:
        sci::NetIO *io;  // Changed from IOPack to NetIO for VOLE compatibility
        sci::OTPack<sci::NetIO> *otpack;
        TripleGenerator<sci::NetIO> *triple_gen;
        MillionaireProtocol<sci::NetIO> *mill;
        MillionaireWithEquality<sci::NetIO> *mill_eq;
        Equality<sci::NetIO> *eq;
        AuxProtocols *aux;
        int party;

        // Constructor - updated for VOLE (no IOPack, use NetIO directly)
        GeometricPerspectiveProtocols(int party, sci::NetIO *io, sci::OTPack<sci::NetIO> *otpack);

        // Destructor
        ~GeometricPerspectiveProtocols();

// void GeometricPerspectiveProtocols::division(uint64_t *input, uint64_t *output, int32_t dim, uint64_t divisor, uint_32_bw);

        // Bit Multiplication Protocol
        void bit_mul(int32_t dim, uint64_t *input, uint64_t *output, int32_t output_bw);

        // MW(x0, x1, L) = Wrap(x0, x1, L) + MSB(x)
        void mw(int32_t dim, uint64_t *input, uint64_t *output, int32_t in_bw, int32_t out_bw);

        void mw_conversion(int32_t dim, uint64_t *input, uint64_t *output, int32_t in_bw, int32_t act_l_bw, int32_t out_bw);

        void mwwithB(int32_t dim,uint64_t B, uint64_t *input, uint64_t *output, int32_t in_bw, int32_t out_bw);

        // Multiplexer protocol with two-bit choice, c is a two-bit number
        void mux_3(int32_t dim, uint64_t *inA, uint64_t *inC, uint64_t *out, int32_t bw);

        // new truncation protocol
        void new_truncate(int32_t dim, uint64_t *inA, uint64_t *outB, int32_t shift, int32_t bw);

        void cross_term(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA, int32_t bwB, int32_t bwC);

        void signed_crossterm(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA, int32_t bwB, int32_t bwC);

        void signed_mul(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA, int32_t bwB, int32_t bwC);

        void sirnn_unsigned_mul(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA, int32_t bwB, int32_t bwC);

        void signed_cipher_plainc(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA, int32_t bwB, int32_t bwC);

        void cross_term_reverse(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC, int32_t bwA, int32_t bwB, int32_t bwC);

        void sign_extension(int32_t dim, uint64_t *in, uint64_t *out, int32_t in_bw, int32_t out_bw);

        void truncate_with_one_bit_error(int32_t dim, uint64_t *inA, uint64_t *outB, int32_t shift, int32_t bw);

        void msb0_truncation(int32_t dim, uint64_t *inA, uint64_t *outB, int32_t shift, int32_t bw);

        void new_msb0_truncation(int32_t dim, uint64_t *inA, uint64_t *outB, int32_t shift, int32_t bw);

        void ring_exp(const uint64_t *inA, uint64_t *out, size_t len, uint64_t f,
            uint64_t bwL, uint64_t in_f);

        void exp(int32_t dim, uint64_t *inA, uint64_t *result, int32_t in_bw,int32_t in_f, int32_t out_bw, int32_t out_f,
            int32_t localexp_bw, int32_t localexp_f, int32_t locallut_bw, int32_t locallut_f);

        void exp_nag4(int32_t dim, uint64_t *inA, uint64_t *result, uint64_t *MW, int32_t in_bw,int32_t in_f, int32_t out_bw, int32_t out_f,
            int32_t localexp_bw, int32_t localexp_f, int32_t locallut_bw, int32_t locallut_f);

        void exp_nagx(int32_t dim, uint64_t *inA, uint64_t *result, int32_t in_bw,int32_t in_f, 
            int32_t localexp_bw, int32_t localexp_f, int32_t locallut_bw, int32_t locallut_f);

        void exp_softmaxx(int32_t dim, uint64_t *inA, uint64_t *result, int32_t in_bw, int32_t in_f,
            int32_t localexp_bw, int32_t localexp_f, int32_t locallut_bw, int32_t locallut_f);

        void vector_bit_mul(int32_t dim, uint64_t *inA, uint8_t ot_choice, uint64_t *outC,
        int32_t bwA);

        void vector_bit_mul_reverse(int32_t dim, uint64_t *inA, uint8_t ot_choice, uint64_t *outC,
        int32_t bwA);

        // 🔥 新增：A_random数组与block128转换函数
        
        // 将A_random数组编码为block128数组（支持任意位宽bwA，末尾补零）
        sci::block128* encode_A_random_to_blocks(
            const uint64_t* A_random, int32_t dim, int32_t bwA);
        
        // 从block128数组解码回A_random数组
        uint64_t* decode_blocks_to_A_random(
            const sci::block128* block_array, int32_t dim, int32_t bwA);

        void vector_scalar_mul(int32_t dim, uint64_t *inA, uint64_t *inB, uint64_t *outC,
        int32_t bwA, int32_t bwB, int32_t bwC);

        sci::block128* encode_matrix_to_blocks(const uint64_t* matrix, int32_t m, int32_t n, int32_t bwA);

        uint64_t* decode_blocks_to_matrix(const sci::block128* block_array, int32_t m, int32_t n, int32_t bwA);


        
        void matrix_vector_crossterm(int32_t m,int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC,
        int32_t bwA, int32_t bwB);

        void matrix_vector_crossterm_reverse(int32_t m,int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC,
        int32_t bwA, int32_t bwB);

        void matrix_vector_mul(int32_t m,int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC,
        int32_t bwA, int32_t bwB, int32_t bwC);
        
        void matrix_vector_unsigned_mul(int32_t m,int32_t n, uint64_t *inA, uint64_t *inB, uint64_t *outC,
        int32_t bwA, int32_t bwB, int32_t bwC);


};
#endif //GEOMTRIC_PERSPECTIVE_PROTOCOLS_H
