/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny parameters and generation of functions for P434
*********************************************************************************************/  

#include "api.h" 
#include "P434_internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 434-bit field element is represented with Ceil(434 / 64) = 7 64-bit digits or Ceil(434 / 32) = 14 32-bit digits.

//
// Curve isogeny system "SIDHp434". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p434^2), where A=6, B=1, C=1 and p434 = 2^216*3^137-1
//
         
const uint64_t p434[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFDC1767AE2FFFFFF, 
                                                     0x7BC65C783158AEA3, 0x6CFC5FD681C52056, 0x0002341F27177344 };
const uint64_t p434p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xFDC1767AE3000000,
                                                     0x7BC65C783158AEA3, 0x6CFC5FD681C52056, 0x0002341F27177344 };
const uint64_t p434x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFB82ECF5C5FFFFFF,
                                                     0xF78CB8F062B15D47, 0xD9F8BFAD038A40AC, 0x0004683E4E2EE688 }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000001000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0x58AEA3FDC1767AE3, 0xC520567BC65C7831, 0x1773446CFC5FD681, 0x0000000002341F27 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p434^2), expressed in unconventinal radix representation
const uint64_t A_gen[12*6]           = { 0x096962979646cb48, 0x0e67ee98f251f25b, 													 0x0b7d224121c1a9aa, 0x0991571076d39726, 
                                                     0x0ceea48ec3be23a8, 0x0a206df75b516cca, 
													 
													 0x00000000000014fc, 0x0000000000000927,
													 0x0000000000000382, 0x0000000000000f64,
													 0x00000000000004e5, 0x00000000000000c9,
													 // XPA0
                                                     0x039baa6869e3ea50, 0x0a45defc3ee851bb, 0x0da0449ed03b3327, 0x064228b3c61617a8, 
                                                     0x0217d7cd13150585, 0x0cc612e004966e27, 
													 
													 0x0000000000001333, 0x00000000000008ff,
													 0x0000000000001374, 0x00000000000010cb,
													 0x0000000000000d3b, 0x000000000000058e,
													 // XPA1
                                                     0x0251f4fc0deb0c6c, 0x0a8b0564ddde3894, 0x069429e9a44b1b9b, 0x04c038b42c12f356, 
                                                     0x02908e16850c6112, 0x0bb4bc3b8833654f,
													 
													 0x00000000000003f2, 0x0000000000000b91,
													 0x00000000000010ce, 0x0000000000000b20,
													 0x0000000000001156, 0x0000000000000294,
													 // XQA0
                                                     0x0a3cb5a2d61c85f5, 0x0180eef9c1afaafc, 0x084d9d1487e31eaa, 0x0daf5ff63046d504,
                                                     0x07fde99eeaea854c, 0x08f15e7f261199ae,
													 
													 0x0000000000000f63, 0x00000000000006f4,
													 0x000000000000031c, 0x0000000000000144,
													 0x00000000000009e5, 0x000000000000007d,
													 // XQA1
                                                     0x07668483c5fd2516, 0x01a1297d946f7f7b, 0x01b113d928608430, 0x0f9ad9907433375a, 
                                                     0x03700428e3d99674, 0x04dc041edc4c8aaf,
													 
													 0x0000000000001566, 0x0000000000000b47,
													 0x0000000000000e25, 0x000000000000131b,
													 0x0000000000000779, 0x0000000000000327,
													 // XRA0
                                                     0x017a88806036c5ac, 0x0abdef2f5b467b30, 0x04152ef8af2872c6, 0x03059bd31b6f8aaa, 
                                                     0x06c96384028796d1, 0x0cab7f7263d95ab6,
													 
													 0x0000000000000b9e, 0x00000000000000ec,
													 0x00000000000010a8, 0x0000000000000bd2,
													 0x00000000000010c1, 0x0000000000000544};  // XRA1
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p434^2), expressed in unconventinal radix representation
const uint64_t B_gen[12*6]           = {0x0064045b0637b7a9, 0x0fb736147d4c4eca, 													 0x0d2374cb30cd6792, 0x06c3f3b558e1d1eb, 
                                                     0x0f4ed8a21d2260fa, 0x09affc6ffae74eac,
													 
													 0x000000000000060d, 0x0000000000001358,
													 0x000000000000118e, 0x0000000000001103,
													 0x0000000000000421, 0x00000000000001bd,
													 // XPB0 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000,
													 
													 0x0000000000000000, 0x0000000000000000,
													 0x0000000000000000, 0x0000000000000000,
													 0x0000000000000000, 0x0000000000000000,
													 // XPB1
                                                     0x0b511d56da3d3846, 0x097b418fca638161, 0x0e95195ff5fd8151, 0x0635cd752e3efd48, 
                                                     0x07b5d058513bf99d, 0x00fabe49422f1702,
													 
													 0x00000000000008f8, 0x0000000000000692,
													 0x0000000000000a87, 0x000000000000151d,
													 0x000000000000113c, 0x00000000000003eb,
													 // XQB0
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000,
													 
													 0x0000000000000000, 0x0000000000000000,
													 0x0000000000000000, 0x0000000000000000,
													 0x0000000000000000, 0x0000000000000000,
													 // XQB1
                                                     0x023776418a52b76c, 0x05b1060f2312a1e3, 0x0439c294bc702ca2, 0x08056c7de32f673b, 
                                                     0x06fcea038f3cfcf7, 0x00f263f076fcba03,
													 
													 0x0000000000001556, 0x00000000000011bd,
													 0x0000000000000873, 0x0000000000000983,
													 0x00000000000002b2, 0x00000000000005f9,
													 // XRB0
                                                     0x0d56f37404736ae4, 0x0dd07f866b369f6f, 0x0506110dfa14de1d, 0x060a47672da434cc, 
                                                     0x03a729930d819d9c, 0x05389ecfa4b34b62,
													 
													 0x0000000000000ce3, 0x0000000000000033,
													 0x00000000000014cb, 0x0000000000000ed7,
													 0x0000000000000f75, 0x000000000000043c};  // XRB1
// Montgomery constant Montgomery_R2 = (2^448)^2 mod p434
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0x28E55B65DCD69B30, 0xACEC7367768798C2, 0xAB27973F8311688D, 0x175CC6AF8D6C7C0B,
                                                     0xABCD92BF2DDE347E, 0x69E16A61C7686D9A, 0x000025A89BCDD12A };                                                   
// Value one in unconventinal radix representation 
const uint64_t uRadix_one[12]    = { 0x0000000000000001, 0x0000000000000000, 													 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 
													 
													 0x0000000000000000, 0x0000000000000000,0x0000000000000000, 0x0000000000000000,0x0000000000000000, 0x0000000000000000};


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
48, 28, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 
1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
66, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 
2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 3, 1, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 
1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };
           
// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy434
#define fpzero                        fpzero434
#define fpadd                         fpadd434
#define fpsub                         fpsub434
#define fpneg                         fpneg434
#define fpdiv2                        fpdiv2_434
#define fpcorrection                  fpcorrection434
#define fpmul_mont                    fpmul434_mont
#define fpsqr_mont                    fpsqr434_mont
#define fpinv_mont                    fpinv434_mont
#define fpinv_chain_mont              fpinv434_chain_mont
#define fpinv_mont_bingcd             fpinv434_mont_bingcd
#define fp2copy                       fp2copy434
#define fp2zero                       fp2zero434
#define fp2add                        fp2add434
#define fp2sub                        fp2sub434
#define fp2neg                        fp2neg434
#define fp2div2                       fp2div2_434
#define fp2correction                 fp2correction434
#define fp2mul_mont                   fp2mul434_mont
#define fp2sqr_mont                   fp2sqr434_mont
#define fp2inv_mont                   fp2inv434_mont
#define fp2inv_mont_bingcd            fp2inv434_mont_bingcd
#define fpequal_non_constant_time     fpequal434_non_constant_time
#define mp_add_asm                    mp_add434_asm
#define mp_subx2_asm                  mp_sub434x2_asm
#define mp_dblsubx2_asm               mp_dblsub434x2_asm

#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"
