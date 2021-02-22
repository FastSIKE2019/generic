/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny parameters and generation of functions for P503
*********************************************************************************************/  

#include "api.h" 
#include "P503_internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 

// Curve isogeny system "SIDHp503". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p503^2), where A=6, B=1, C=1 and p503 = 2^250*3^159-1

// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0400000000000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0xC216F6888479E82B, 0xE6FDB21EDF9F6BC4, 0x1171AF769DE93406, 0x1019BD5060478798 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + XQA1*i, XRA0 + XRA1*i} in GF(p503^2), expressed in unconventinal radix representation
const uint64_t A_gen[6*10]                       = { 0x0003df2ded851678 ,0x00024544dfb8c7ef ,0x0001dfa18972f305 ,0x000233bb246cf157 ,0x0001157e7fa244ac ,
                                                     0x00032be5da1065fe ,0x0003ec8b51728cad ,0x0001885676035346 ,0x0002426ba21f4b99 ,0x000013e60ac1c454 , // XPA0
													  
													 0x00005b8a85f915d2 ,0x00023e890ad093fe ,0x00023f12edf076ff ,0x0004e2c925e62dba ,0x000252f8a2ef9a2a ,
                                                     0x0004cfdd29a5a184 ,0x0002c24f60966610 ,0x0000ab7e9d082fd6 ,0x000470d348fb2a27 ,0x0000d20f742a58e9 , // XPA1
					
													 0x0002aa35c7616461 ,0x0002b993431c19d0 ,0x00003282e0874e73 ,0x00014ff0b6e2d203 ,0x0000384bd7dcfac0 ,
                                                     0x0000bacc66306cbd ,0x0002e21fce444ab3 ,0x00023638e30354d9 ,0x0001287044d40fc4 ,0x0001566ff51d05b9 , // XQA0

													 0x00027741f9e57db1 ,0x00051fb871935d0e ,0x00017a6b584f288e ,0x00039f5caece39ce ,0x00028dad81748dc2 ,
                                                     0x00018773e39342a3 ,0x00013355460e6007 ,0x000349f03a3e8858 ,0x00016f79d9113c44 ,0x0001a8d3f4a77554 , // XQA1
 
													 0x000317e75ec8ce9f ,0x000175981bc5a756 ,0x000484d46d070e2e ,0x0001a57613cd21e9 ,0x000318ecead636a5 ,
                                                     0x0003fb25ae2ab4f7 ,0x0002c64b244d0e9f ,0x0002ed51a5dcf7ae ,0x0000de13eca6c308 ,0x00019fbd2735eb4d , // XRA0

													 0x0004f3f1e0029498 ,0x0002cc80e2674fde ,0x000301850f195fda ,0x0001a985bcbc72d5 ,0x00046a933be272e8 ,
                                                     0x0002fe9649c12a2b ,0x00047a245f40ccc0 ,0x0002a2f3190134c4 ,0x00048d324b28d6bc ,0x00002c5f26e59042   // XRA1
													 }; 
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p503^2), expressed in unconventinal radix representation
const uint64_t B_gen[6*10]                       = { 0x00036340d578badf ,0x000261fcfc29f0e6 ,0x0001fa8b960f576c ,0x00046d5045fdc31e ,0x0004d2c3c3b0ce1a ,
                                                     0x000513897666bb97 ,0x0003ee80a8176b2a ,0x0000e0c455443c6e ,0x0001a4a0ef410baa ,0x0001597fd3d869aa , // XPB0
 
													 0x0000000000000000 ,0x0000000000000000 ,0x0000000000000000 ,0x0000000000000000 ,0x0000000000000000 , 
													 0x0000000000000000 ,0x0000000000000000 ,0x0000000000000000 ,0x0000000000000000 ,0x0000000000000000 , // XPB1
                                                     
													 0x0001d54bc2c24d53 ,0x0003b85786e0a030 ,0x0002e28728056d42 ,0x000233c582ff3a9c ,0x0000d28fff64cdba ,
                                                     0x00040a94f3612966 ,0x000413be14c3add1 ,0x0000e3b0d1f960f0 ,0x00045ee6cc8fc769 ,0x000183991a441c4d , // XQB0												 
													 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
													 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,  // XQB1											 
                                             
													 0x000096034f2cb522 ,0x00021f8236247f86 ,0x00010dfd999e3695 ,0x0001a5f69a933a9b ,0x00036e4eb67e6296 ,
                                                     0x00030648296a8263 ,0x00044d6d2dec66d9 ,0x0003406b3c7235dd ,0x0003af0f701a337c ,0x00000522257f771d , // XRB0
												 
													 0x00037e560b981758 ,0x000518798b8570cb ,0x0004ca0942c2b871 ,0x0001d74497b8c869 ,0x0002b567ea67cf42 ,
                                                     0x000363bdd17ed88c ,0x0001269fc910a0d8 ,0x0000b77c59500a68 ,0x00008f66afd0f29a ,0x000113694cb2fdc8   // XRB1								
													 }; 
                                               
// Value one in unconventinal radix representation 
const uint64_t uRadix_one[10]                    = { 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000};


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
61, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 
4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 
1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 29, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 
1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 
1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
71, 38, 21, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 
1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 
5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 
1, 4, 2, 1, 1, 2, 1, 1, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 
2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 
1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };
           
// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy503
#define fpzero                        fpzero503
#define fp2zero                       fp2zero503
#define fp2neg                        fp2neg503
#define fp2inv_mont_bingcd            fp2inv503_mont_bingcd
#define fpequal_non_constant_time     fpequal503_non_constant_time
#define mp_add_asm                    mp_add503_asm
#define mp_subx2_asm                  mp_sub503x2_asm
#define mp_dblsubx2_asm               mp_dblsub503x2_asm

#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"
