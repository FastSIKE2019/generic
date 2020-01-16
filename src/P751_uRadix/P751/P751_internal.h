/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: internal header file for P751
*********************************************************************************************/  

#ifndef __P751_INTERNAL_H__
#define __P751_INTERNAL_H__

#include "../../config.h"
 

#if (TARGET == TARGET_AMD64)
    #define NWORDS_FIELD    12              // Number of words of a 751-bit field element
    #define p751_ZERO_WORDS 5               // Number of "0" digits in the least significant part of p751 + 1     
#elif (TARGET == TARGET_x86)
    #define NWORDS_FIELD    24 
    #define p751_ZERO_WORDS 11
#elif (TARGET == TARGET_ARM)
    #define NWORDS_FIELD    24
    #define p751_ZERO_WORDS 11
#elif (TARGET == TARGET_ARM64)
    #define NWORDS_FIELD    12
    #define p751_ZERO_WORDS 5
#endif 
    

// Basic constants

#define NBITS_FIELD             751  
#define MAXBITS_FIELD           768                
#define MAXWORDS_FIELD          ((MAXBITS_FIELD+RADIX-1)/RADIX)     // Max. number of words to represent field elements
#define NWORDS64_FIELD          ((NBITS_FIELD+63)/64)               // Number of 64-bit words of a 751-bit field element 
#define NBITS_ORDER             384
#define NWORDS_ORDER            ((NBITS_ORDER+RADIX-1)/RADIX)       // Number of words of oA and oB, where oA and oB are the subgroup orders of Alice and Bob, resp.
#define NWORDS64_ORDER          ((NBITS_ORDER+63)/64)               // Number of 64-bit words of a 384-bit element 
#define MAXBITS_ORDER           NBITS_ORDER                         
#define MAXWORDS_ORDER          ((MAXBITS_ORDER+RADIX-1)/RADIX)     // Max. number of words to represent elements in [1, oA-1] or [1, oB].
#define ALICE                   0
#define BOB                     1 
#define OALICE_BITS             372  
#define OBOB_BITS               379    
#define OBOB_EXPON              239 
#define MASK_ALICE              0x0F
#define MASK_BOB                0x03  
#define PRIME                   p751  
#define PARAM_A                 6  
#define PARAM_C                 1
// Fixed parameters for isogeny tree computation
#define MAX_INT_POINTS_ALICE    8      
#define MAX_INT_POINTS_BOB      10 
#define MAX_Alice               186
#define MAX_Bob                 239
#define MSG_BYTES               32
#define SECRETKEY_A_BYTES       (OALICE_BITS + 7) / 8
#define SECRETKEY_B_BYTES       (OBOB_BITS - 1 + 7) / 8
#define FP2_ENCODED_BYTES       2*((NBITS_FIELD + 7) / 8)

// SIDH's basic element definitions and point representations

typedef digit_t felm_t[NWORDS_FIELD];                                 // Datatype for representing 751-bit field elements (768-bit max.)
typedef digit_t dfelm_t[2*NWORDS_FIELD];                              // Datatype for representing double-precision 2x751-bit field elements (2x768-bit max.) 
typedef felm_t  f2elm_t[2];                                           // Datatype for representing quadratic extension field elements GF(p751^2)
        
typedef struct { f2elm_t X; f2elm_t Z; } point_proj;                  // Point representation in projective XZ coordinates.
typedef point_proj point_proj_t[1]; 


/**************** Function prototypes ****************/

// Copy wordsize digits, c = a, where lng(a) = nwords
void copy_words(const digit_t* a, digit_t* c, const unsigned int nwords);

// addition
void digit_a_digit_2(digit_t a, digit_t b, digit_t* c);
void digit_a_digit_1(const digit_t* a, const digit_t* b, digit_t* c);

// Digit multiplication
void digit_x_digit(const digit_t a, const digit_t b, digit_t* c);
void digit_x_digit_2(const digit_t a, const digit_t b, digit_t* co);

//
void mp_mul_new4(const digit_t* a, const digit_t* b, digit_t* c_0, digit_t* c_1, digit_t* ca, const unsigned int nwords);
void improved_BR(const digit_t* t, const digit_t s, digit_t* q, digit_t* r);
void improved_IFFM(digit_t* c_0, digit_t* c_1, digit_t* ca, digit_t* C, const unsigned int nwords);

// Multiprecision comba multiply, c = a*b, where lng(a) = lng(b) = nwords.
void mp_mul_old(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords);

/************ Field arithmetic functions *************/

// Copy of a field element, c = a
void fpcopy751(const digit_t* a, digit_t* c);

// Zeroing a field element, a = 0
void fpzero751(digit_t* a);

//Modular Addition based on the new data representation , with the unconditional radix
extern void fp_add_new(digit_t* a, digit_t* b, digit_t* c, const unsigned int nwords);

//Modular Subtraction based on the new data representation , with the unconditional radix
extern void fpsub_new(digit_t* a, digit_t* b, digit_t* c, const unsigned int nwords);

//Modular Negation based on the new data representation , with the unconditional radix  
extern void fpneg_new(digit_t* a,digit_t* b,const unsigned int nwords);

//Modular Division by Two based on the new data representation , with the unconditional radix
extern void fpdiv2_751_new(digit_t* a, digit_t* c);
            
//Modular Multipication based on the new data representation , with the unconditional radix
void fpmul_new(const digit_t* a, digit_t* b, digit_t* c);
   
// Field squaring using uRadix arithmetic
void mp_sqr_new(const digit_t* a, digit_t* c_0, digit_t* c_1, digit_t* ca, const unsigned int nwords);
void fpsqr_new(const digit_t* a, digit_t* c);

// Conversion to uRadix representation
void IBR_uRadix(const digit_t* a, digit_t* q, digit_t* r, const unsigned int unwords, const digit_t* lamda);
void to_uRadix(digit_t* a, digit_t* uc, const unsigned int nwords);
    
// Conversion from uRadix representation to standard representation
void from_uRadix(digit_t* a, digit_t* c, const unsigned int nwords);

// Field inversion, a = a^-1 in GF(p751)
void fpinv_new(digit_t* a, digit_t* c, const unsigned int nwords);

//
void fpinv_chain_new(digit_t* a , const unsigned int nwords);

/************ GF(p^2) arithmetic functions *************/
    
// Copy of a GF(p751^2) element, c = a
void fp2copy751(const f2elm_t a, f2elm_t c);

// Zeroing a GF(p751^2) element, a = 0
void fp2zero751(f2elm_t a);

// GF(p751^2) negation, a = -a in GF(p751^2)
void fp2neg751(f2elm_t a);

// GF(p751^2) addition, c = a+b in GF(p751^2)
extern void fp2add751(f2elm_t a, f2elm_t b, f2elm_t c);           

// GF(p751^2) subtraction, c = a-b in GF(p751^2)
extern void fp2sub751(f2elm_t a, f2elm_t b, f2elm_t c); 

// GF(p751^2) inversion using uRadix arithmetic, a = (a0-i*a1)/(a0^2+a1^2)
void fp2inv_uRadix(f2elm_t a);

/************ Elliptic curve and isogeny functions *************/

// Differential addition.
void xADD(point_proj_t P, const point_proj_t Q, const f2elm_t xPQ);

// Evaluates the isogeny at the point (X:Z) in the domain of the isogeny.
void eval_4_isog(point_proj_t P, f2elm_t* coeff);

// 3-way simultaneous inversion
void inv_3_way(f2elm_t z1, f2elm_t z2, f2elm_t z3);

#endif
