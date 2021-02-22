/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny parameters and generation of functions for P610
*********************************************************************************************/  

#include "api.h" 
#include "P610_internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 610-bit field element is represented with Ceil(610 / 64) = 10 64-bit digits or Ceil(610 / 32) = 20 32-bit digits.

//
// Curve isogeny system "SIDHp610". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p610^2), where A=6, B=1, C=1 and p610 = 2^305*3^192-1
//
         
const uint64_t p610[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x6E01FFFFFFFFFFFF, 
                                                     0xB1784DE8AA5AB02E, 0x9AE7BF45048FF9AB, 0xB255B2FA10C4252A, 0x819010C251E7D88C, 0x000000027BF6A768 };
const uint64_t p610p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x6E02000000000000,
                                                     0xB1784DE8AA5AB02E, 0x9AE7BF45048FF9AB, 0xB255B2FA10C4252A, 0x819010C251E7D88C, 0x000000027BF6A768 };
const uint64_t p610x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xDC03FFFFFFFFFFFF,
                                                     0x62F09BD154B5605C, 0x35CF7E8A091FF357, 0x64AB65F421884A55, 0x03202184A3CFB119, 0x00000004F7ED4ED1 }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0002000000000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0x26F4552D58173701, 0xDFA28247FCD5D8BC, 0xD97D086212954D73, 0x086128F3EC46592A, 0x00013DFB53B440C8 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p610^2), expressed in unconventinal radix representation

const uint64_t A_gen[12 * NWORDS64_FIELD]         = { 0x0ca4ca3070a5f26c, 0x0e6b2c8558468b54, 		  											 0x01e999c7300db042, 0x0cc4e2afe399642d, 													 0x03d0fb18f50b6772, 0x0ab4451b7eb1f2e8,//XPA0
                                                     0x000000f0b412421e, 0x000001493d9bdc06, 
													 0x0000024f23d9655b, 0x00000124e22cd520, 
													 0x000000d6dd043dda, 0x0000012075c5cd65,//XPA01
                                                     0x083601f0637c1aa2, 0x069cb9db46198255, 0x099c5ec0dee0cd95, 0x0d26d3af15f94488, 0x0e7a12e6a2a8959a, 0x00f535e118769642,//XPA1
                                                     0x0000031fea2e569a, 0x000002c3027b181f, 0x00000227212a4040, 0x000002924e11f883, 0x0000029cd2b626b1, 0x000000d7b4026ab6,//XPA11
                                                     0x045438b283798b8a, 0x0646f6f5daeff2f3, 0x0679c48546772a05, 0x0d421c113cfa7d0b, 0x0575f0fb4d7a9e1e, 0x0a0764b948791417,//XQA0
                                                     0x000002eb9eafc363, 0x00000217088359fc, 0x0000020a8cdb3190, 0x000002ec1dc6ff3e, 0x000002b547a84a6b, 0x0000001913cef264,//XQA01
                                                     0x0309fb8c98f0a0ed, 0x06ee815903ee7507, 0x091c6b014e1a217d, 0x0d483453efc109ad, 0x0581f2b4b769bf6b, 0x09ed28ee2e40376c,//XQA1
                                                     0x00000175acb6d2ad, 0x00000263dcf4c2a0, 0x000000a644db9d16, 0x000000b7a8a53b5b, 0x000000e65a11e3a8, 0x0000016bc039eeee,//XQA11
                                                     0x0f0af07fda163fb5, 0x04179990f496b49d, 0x040c1e965bdbdd86, 0x02a9ed9a82edd398, 0x0dd13c4a4ed39d85, 0x02f14649d4dfb969,//XRA0
                                                     0x000002bba5006525, 0x0000019387e73330, 0x00000241f1d429ca, 0x0000008520c588ca, 0x000002da615a4926, 0x00000120769c7ca5,//XRA01
                                                     0x092f37747518a992, 0x068c158cf9457694, 0x0e98cd2ac1e303c1, 0x084bfad6b8149cb2, 0x09b71fe2795a7f06, 0x0c0bad4d85cdf168,//XRA1
                                                     0x000002514fab2148, 0x0000023275d0d43e, 0x0000004d63d46547, 0x00000335e4f1fadf, 0x00000119d65c0a7c, 0x00000100e926888c};//XRA11
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p610^2), expressed in unconventinal radix representation



const uint64_t B_gen[12 * NWORDS64_FIELD]         = {0x0e3417a56f7fefeb, 0x0283eb9cb640875c, 	  												 0x00028d86d413ba6e, 0x09a0ad948c5aa614, 													 0x016522515e4bc59b, 0x0f75f8b00e29d778,//XPB0
                                                     0x00000111697e528a, 0x000001b4eb6e0370, 0x0000002861c72679, 0x0000030248bda00b, 0x000000c058261906, 0x000000e4364b87b7,//XPB01
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,//XPB1
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,//XPB11
                                                     0x0fc48de126b3a047, 0x02af9d6d1fec2266, 0x06a2866a62e0b283, 0x0708a59385c8f5b7, 0x032b76114552baae, 0x08b76d041f3c2cd1,//XQB0
                                                     0x00000258dfa0ce53, 0x000000107a16b44b, 0x000001aaa433819e, 0x00000236acc644f2, 0x000002d405a29bd0, 0x000000dd8944bca1,//XQB01
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,//XQB1
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,//XQB11
                                                     0x038e5f52ded74222, 0x0bb148653d8a3f7c, 0x0006770399f79311, 0x0431d74efc15cf60, 0x0ef92d11342db38b, 0x0956ad5f2c212d64,//XRB0
                                                     0x0000002b8badc62e, 0x000000599cf69b0f, 0x00000206d0359310, 0x000000179ffe7e59, 0x000000c025ab2762, 0x0000013afd185733,//XRB01
                                                     0x0a9c6edb4c238b48, 0x0d8a0cc39cb8d6f9, 0x078f9ec4f37e86c9, 0x07b2cca02c6fe64b, 0x008b08fb25ea2c40, 0x0a0096dee343be55,//XRB1
                                                     0x00000134f6d77b69, 0x00000145be5966f5, 0x0000026b8925b121, 0x000002638afdffb7, 0x0000008d19ac7232, 0x0000012007f9da5d};//XRB11

// Montgomery constant Montgomery_R2 = (2^640)^2 mod p610
const uint64_t R2[NWORDS64_FIELD]     = { 0x0e6f5d201a197727, 0x0d7705e906e70e1e, 0x07c8e4e5ede44e84, 0x0cba1effe2e07a86, 0x06b024f6085a3c7c, 0x07628b29b48c3d39,
                                                     0x00000238f029e43d, 0x0000012c135d12d0, 0x000000e66ea784e6, 0x000000f65ca808ca, 0x000000cc922d4664, 0x0000005e2e878fb9 };                                                    
// Value one in unconventinal radix representation 
const uint64_t uRadix_one[NWORDS64_FIELD]    = { 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                                     0x0000000000000000  };


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
67, 37, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 
2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 
5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 
1, 4, 2, 1, 1, 2, 1, 1, 33, 16, 8, 5, 2, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 
1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 1, 1, 
4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
86, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 
1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 38, 
21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 
9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 
1, 1 };

// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy610
#define fpzero                        fpzero610
#define fpadd                         fpadd610
#define fpsub                         fpsub610
#define fpneg                         fpneg610
#define fpdiv2                        fpdiv2_610
#define fpcorrection                  fpcorrection610
#define fpmul_mont                    fpmul610_mont
#define fpsqr_mont                    fpsqr610_mont
#define fpinv_mont                    fpinv610_mont
#define fpinv_chain_mont              fpinv610_chain_mont
#define fpinv_mont_bingcd             fpinv610_mont_bingcd
#define fp2copy                       fp2copy610
#define fp2zero                       fp2zero610
#define fp2add                        fp2add610
#define fp2sub                        fp2sub610
#define fp2neg                        fp2neg610
#define fp2div2                       fp2div2_610
#define fp2correction                 fp2correction610
#define fp2mul_mont                   fp2mul610_mont
#define fp2sqr_mont                   fp2sqr610_mont
#define fp2inv_mont                   fp2inv610_mont
#define fp2inv_mont_bingcd            fp2inv610_mont_bingcd
#define fpequal_non_constant_time     fpequal610_non_constant_time
#define mp_add_asm                    mp_add610_asm
#define mp_subx2_asm                  mp_sub610x2_asm
#define mp_dblsubx2_asm               mp_dblsub610x2_asm

#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"