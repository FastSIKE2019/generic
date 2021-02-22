/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny parameters and generation of functions for P751
*********************************************************************************************/  

#include "api.h" 
#include "P751_internal.h"

// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 

//
// Curve isogeny system "SIDHp751". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p751^2), where A=6, B=1, C=1 and p751 = 2^372*3^239-1
//
                  
const uint64_t p751[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xEEAFFFFFFFFFFFFF,
                                                     0xE3EC968549F878A8, 0xDA959B1A13F7CC76, 0x084E9867D6EBE876, 0x8562B5045CB25748, 0x0E12909F97BADC66, 0x00006FE5D541F71C };
const uint64_t p751p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xEEB0000000000000,
                                                     0xE3EC968549F878A8, 0xDA959B1A13F7CC76, 0x084E9867D6EBE876, 0x8562B5045CB25748, 0x0E12909F97BADC66, 0x00006FE5D541F71C };
const uint64_t p751x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xDD5FFFFFFFFFFFFF, 
                                                     0xC7D92D0A93F0F151, 0xB52B363427EF98ED, 0x109D30CFADD7D0ED, 0x0AC56A08B964AE90, 0x1C25213F2F75B8CD, 0x0000DFCBAA83EE38 }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0010000000000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0xC968549F878A8EEB, 0x59B1A13F7CC76E3E, 0xE9867D6EBE876DA9, 0x2B5045CB25748084, 0x2909F97BADC66856, 0x06FE5D541F71C0E1 };

// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p751^2), expressed in unconventinal radix representation
const uint64_t A_gen[6 * NWORDS64_FIELD]         = { 0x2f5d062b87f8b6fa, 0x6544a182b221f3a4, 0x4bc20b5d4729e702, 0x1b5707c0c0a1da49, 0x34fa219056b2c7dc, 0x2a2c6ea37339da66, 0x5f317c8685314ee6, 0x333bbf442c4800a0, 0x1b58f51a8aa31ebc, 0x22a80a4637118bf1, 0x0ac4d392791245f6, 0x156267611b91addc, // XPA0
													 0x384d3a260af02041, 0x350effe35b5729c7, 0x2ff76d119da11654, 0x37c44560cb68c0e8, 0x0c33e02d08e5dbf4, 0x167e625e59713733, 0x53b2de48b2317d83, 0x49eefcc9af4bd428, 0x4a5febe0915fe4cf, 0x0c1b792310e78455, 0x2e88095041a4b9d8, 0x06ab18b504404815, // XPA1
													 0x136194726276290e, 0x111823b28d741f6a, 0x05310302eebb3a0a, 0x27617a14b8b94a79, 0x0757263459c316e3, 0x208a03e63c1b31cd, 0x55179019c1a7665d, 0x2e1f826dd01294d2, 0x5834ba229ec9f156, 0x2ef67e640f2e0ff0, 0x13e038ee0416bf35, 0x0729b9f840802f74, // XQA0
													 0x02d2332cdb8c35c7, 0x23d67536103d11a2, 0x028a5d53fa46186c, 0x11fb3f086abef25b, 0x40795b0eff744f00, 0x33e8c13742473147, 0x41b7e74af1b24690, 0x116ac97c47b6eb78, 0x0b99b01f440d641b, 0x65fd7f2f3725004a, 0x604deb08adae52e2, 0x0b94d52b1b649b86, // XQA1
													 0x174ccf79967525bb, 0x209406fecb4be93f, 0x04eaadd94825d63e, 0x257aae548e486065, 0x0d0e933f1104c424, 0x13ebd49b4463fd34, 0x3d608d9ba9e51747, 0x37c77e6277bdddd9, 0x2f3d61dae887ca05, 0x574f56e0986ccfa0, 0x3c107e0b99d1c33f, 0x1dd760178377b8e5, // XRA0
													 0x0fd19abbfada5165, 0x6304aa18a1aaa718, 0x0f740c49e7b08ba6, 0x4e652915c86b9d76, 0x41965f55c93b68d2, 0x158681bbbdfa4ab8, 0x2bbcdac1997f215e, 0x664ff9ac2949a105, 0x60caa828811dcb49, 0x4151d83126f772f2, 0x6302fe978996e4d9, 0x1909e412299fccbe};// XRA1
													 
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p751^2), expressed in unconventinal radix representation
const uint64_t B_gen[6 * NWORDS64_FIELD]         = { 0x44979607453432fe, 0x01e2f4d8ebec602e, 0x5b007d5a412a7e86, 0x10a3e0ac1d7241dc, 0x04f091103f921956, 0x47a869fe3be79063, 
                                                     0x138d57e15a4c7c03, 0x5c4c4c946335f533, 0x43dfb9862929a724, 0x46761dab003c7e7d, 0x08df4fc98cafa2d3, 0x1dd4673e9bc75809,//XPB0 
													 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,//XPB1
													 
													 0x261b8ce38f8cd0f1, 0x1e1fb6753dd28dd5, 0x55024919131879fb, 0x395e74866a5dc1e1, 0x5209adbf484e91b7, 0x5bd0658abf690a92, 
                                                     0x2ac3049608467766, 0x4136fd83a32e6e6f, 0x380f29bbe89e1c85, 0x0f0b1a0080691396, 0x04d886c7021ad74d, 0x1c787bfcf0a4d3c2,//XQB0 
													 
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,//XQB1
                                             		 
													 0x30bdbed906a79fe3, 0x3de7d0d6cf84a4c1, 0x1f3e0070eace9b5a, 0x3e88cd7c1399e730, 0x05219d0a408a17d3, 0x17d7f1ecf5b3bace, 
                                                     0x647d47ead56307cf, 0x279809ec160c4ece, 0x5db6e1f21019053f, 0x06688a45ec544dc1, 0x11170e1d19b1e9e2, 0x1a96be2f0ddd4367,//XRB0 
												 
													 0x0fa7270dba06c821, 0x5dd688eeb3fd6608, 0x0b2e50d4272f7e80, 0x47a42cd75b9b3859, 0x1c978f457d58fecd, 0x58cdce8a89ab615b, 
                                                     0x158e07ca899cdbd4, 0x53315aebf9aa1c9a, 0x2e142a8dd386d7f3, 0x5f6248b7d740bcb1, 0x63ac74ec0947c08d, 0x1c19319fd1d52636};//XRB1
                                                    
// Value one in unconventinal radix representation 
const uint64_t uRadix_one[NWORDS64_FIELD]        = { 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000,
                                                     0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000};

// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
80, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 
1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 
1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 
33, 20, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 
1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 
1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
112, 63, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 
1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 
1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 31, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 
1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 
2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 49, 31, 16, 8, 4, 2, 
1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 
15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 
1, 1, 1, 21, 12, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 3, 2, 1, 1, 1, 1, 
2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy751
#define fpzero                        fpzero751
#define fp2zero                       fp2zero751
#define fp2neg                        fp2neg751
#define fp2inv_mont_bingcd            fp2inv751_mont_bingcd
#define fpequal_non_constant_time     fpequal751_non_constant_time
#define mp_add_asm                    mp_add751_asm
#define mp_subx2_asm                  mp_sub751x2_asm
#define mp_dblsubx2_asm               mp_dblsub751x2_asm

#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"