/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: internal header file for P434
*********************************************************************************************/  

#ifndef __P434_INTERNAL_H__
#define __P434_INTERNAL_H__

#include "../../config.h"
 

#if (TARGET == TARGET_AMD64)
    #define NWORDS_FIELD    6               // Number of words of a 434-bit field element
    #define p434_ZERO_WORDS 3               // Number of "0" digits in the least significant part of p434 + 1     
#elif (TARGET == TARGET_x86)
    #define NWORDS_FIELD    14 
    #define p434_ZERO_WORDS 6
#elif (TARGET == TARGET_ARM)
    #define NWORDS_FIELD    14
    #define p434_ZERO_WORDS 6
#elif (TARGET == TARGET_ARM64)
    #define NWORDS_FIELD    7
    #define p434_ZERO_WORDS 3
#endif
    

// Basic constants

#define NBITS_FIELD             434  
#define MAXBITS_FIELD           448                
#define MAXWORDS_FIELD          ((MAXBITS_FIELD+RADIX-1)/RADIX)     // Max. number of words to represent field elements
#define NWORDS64_FIELD          ((NBITS_FIELD+63)/64)               // Number of 64-bit words of a 434-bit field element 
#define NBITS_ORDER             224
#define NWORDS_ORDER            ((NBITS_ORDER+RADIX-1)/RADIX)       // Number of words of oA and oB, where oA and oB are the subgroup orders of Alice and Bob, resp.
#define NWORDS64_ORDER          ((NBITS_ORDER+63)/64)               // Number of 64-bit words of a 224-bit element 
#define MAXBITS_ORDER           NBITS_ORDER                         
#define MAXWORDS_ORDER          ((MAXBITS_ORDER+RADIX-1)/RADIX)     // Max. number of words to represent elements in [1, oA-1] or [1, oB].
#define ALICE                   0
#define BOB                     1 
#define OALICE_BITS             216  
#define OBOB_BITS               218     
#define OBOB_EXPON              137    
#define MASK_ALICE              0xFF 
#define MASK_BOB                0x01 
#define PRIME                   p434 
#define PARAM_A                 6  
#define PARAM_C                 1
// Fixed parameters for isogeny tree computation
#define MAX_INT_POINTS_ALICE    7        
#define MAX_INT_POINTS_BOB      8      
#define MAX_Alice               108
#define MAX_Bob                 137
#define MSG_BYTES               16
#define SECRETKEY_A_BYTES       (OALICE_BITS + 7) / 8
#define SECRETKEY_B_BYTES       (OBOB_BITS - 1 + 7) / 8
#define FP2_ENCODED_BYTES       2*((NBITS_FIELD + 7) / 8)


// SIDH's basic element definitions and point representations

typedef digit_t felm_t[NWORDS_FIELD];                                 // Datatype for representing 434-bit field elements (448-bit max.)
typedef digit_t dfelm_t[2*NWORDS_FIELD];                              // Datatype for representing double-precision 2x434-bit field elements (448-bit max.) 
typedef felm_t  f2elm_t[2];                                           // Datatype for representing quadratic extension field elements GF(p434^2)
        
typedef struct {f2elm_t X; f2elm_t Xi; f2elm_t Z; f2elm_t Zi;} point_proj;                  
typedef point_proj point_proj_t[1]; 

typedef digit_t normal_t[7];
typedef digit_t dnormal_t[2*7];
typedef normal_t normal2_t[2];

/**************** Function prototypes ****************/
/************* Multiprecision functions **************/ 

// Copy wordsize digits, c = a, where lng(a) = nwords
void copy_words(const digit_t* a, digit_t* c, const unsigned int nwords);

// Digit multiplication, digit * digit -> 2-digit result
void digit_x_digit_MUL(const digit_t a, const digit_t b, digit_t* c);

// Multiprecision comba multiply, c = a*b, where lng(a) = lng(b) = nwords.
void mp_mul_old(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords);

/************ Field arithmetic functions *************/

// Copy of a field element, c = a
void fpcopy434(const digit_t* a, digit_t* c);

// Zeroing a field element, a = 0
void fpzero434(digit_t* a);

//Modular Addition based on the new data representation , with the unconditional radix
extern void fp_add_new(digit_t* a1, digit_t* a0, digit_t* b1, digit_t* b0, digit_t* c0, digit_t* c1);

//Modular Subtraction based on the new data representation , with the unconditional radix
extern void fpsub_new(digit_t* a1, digit_t* a0, digit_t* b1, digit_t* b0, digit_t* c0, digit_t* c1);

//Modular Negation based on the new data representation , with the unconditional radix  
extern void fpneg_new(digit_t* a1, digit_t* a0, digit_t* b0, digit_t* b1);  

//Modular Division by Two based on the new data representation , with the unconditional radix
void fpdiv2_434_new(digit_t* a1, digit_t* a0, digit_t* c0, digit_t* c1);

//Modular Multipication based on the new data representation , with the unconditional radix
void fpmul( digit_t* a1,  digit_t* a0,  digit_t* b1,  digit_t* b0, digit_t* d0, digit_t* d1);
   
// Field squaring using uRadix arithmetic
void fpsqr( digit_t* a1,  digit_t* a0, digit_t* d0, digit_t* d1);

// Conversion to uRadix representation
void to_uRadix(digit_t* a, digit_t* c0, digit_t* c1);
    
// Conversion from uRadix representation to standard representation
void from_uRadix(digit_t* a1, digit_t* a0, digit_t* co);

// Field inversion, a = a^-1 in GF(p434)
void fpinv_new_u(digit_t* a1, digit_t* a0, digit_t* c0, digit_t* c1);


void fpinv_chain_new_u(digit_t* a1, digit_t* a0, digit_t* c0, digit_t* c1);

/************ GF(p^2) arithmetic functions *************/
    
// Copy of a GF(p434^2) element, c = a
void fp2copy434(const f2elm_t a, f2elm_t c);

// Zeroing a GF(p434^2) element, a = 0
void fp2zero434(f2elm_t a);

// GF(p434^2) negation, a = -a in GF(p434^2)
void fp2neg_new(f2elm_t a, f2elm_t ai, f2elm_t b, f2elm_t bi);

// GF(p434^2) addition, c = a+b in GF(p434^2)
extern void fp2add_new(f2elm_t a, f2elm_t ai, f2elm_t b, f2elm_t bi, f2elm_t c, f2elm_t ci);           

// GF(p434^2) subtraction, c = a-b in GF(p434^2)
extern void fp2sub_new(f2elm_t a, f2elm_t ai, f2elm_t b, f2elm_t bi, f2elm_t c, f2elm_t ci); 

// GF(p434^2) division by two, c = a/2  in GF(p434^2) 
void fp2div2(const f2elm_t a, const f2elm_t ai, f2elm_t c, f2elm_t ci) ;

            
// GF(p434^2) squaring using arithmetic, c = a^2 in GF(p434^2)
void fp2sqr_uRadix(f2elm_t a, f2elm_t ai, f2elm_t c, f2elm_t ci);
 
// GF(p434^2) multiplication using arithmetic, c = a*b in GF(p434^2)
void fp2mul_uRadix(f2elm_t a, f2elm_t ai, f2elm_t b, f2elm_t bi, f2elm_t c, f2elm_t ci);
    
// Conversion of a GF(p434^2) element to uRadix representation
void to_fp2uRadix(const normal2_t a, f2elm_t mc, f2elm_t mci);

// Conversion of a GF(p434^2) element from uRadix representation to standard representation
void from_fp2uRadix(const f2elm_t a, const f2elm_t ai, normal2_t c);

// GF(p434^2) inversion, a = (a0-i*a1)/(a0^2+a1^2)
void fp2inv_uRadix(f2elm_t a, f2elm_t ai);

/************ Elliptic curve and isogeny functions *************/

// Computes the j-invariant of a Montgomery curve with projective constant.
void j_inv(const f2elm_t A, const f2elm_t Ai, const f2elm_t C, const f2elm_t Ci, f2elm_t jinv, f2elm_t jinvi);

// Simultaneous doubling and differential addition.
void xDBLADD(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t xPQi, const f2elm_t A24, const f2elm_t A24i);

// Doubling of a point in projective coordinates (X:Z).
void xDBL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t A24plusi, const f2elm_t C24, const f2elm_t C24i);

// Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.
void xDBLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t A24plusi, const f2elm_t C24, const f2elm_t C24i, const int e);

// Computes the corresponding 2-isogeny of a projective point (X2:Z2) of order 2.
void get_2_isog(const point_proj_t P, f2elm_t A, f2elm_t Ai, f2elm_t C, f2elm_t Ci);

// Evaluates the isogeny at the point (X:Z) in the domain of the isogeny.
void eval_2_isog(point_proj_t P, point_proj_t Q);

// Computes the corresponding 4-isogeny of a projective point (X4:Z4) of order 4.
void get_4_isog(const point_proj_t P, f2elm_t A24plus, f2elm_t A24plusi, f2elm_t C24, f2elm_t C24i, f2elm_t* coeff);

// Evaluates the isogeny at the point (X:Z) in the domain of the isogeny.
void eval_4_isog(point_proj_t P, f2elm_t* coeff);

// Tripling of a point in projective coordinates (X:Z).
void xTPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24minusi, const f2elm_t A24plus, const f2elm_t A24plusi)  ;

// Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
void xTPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24minusi, const f2elm_t A24plus, const f2elm_t A24plusi, const int e);

// Computes the corresponding 3-isogeny of a projective point (X3:Z3) of order 3.
void get_3_isog(const point_proj_t P, f2elm_t A24minus, f2elm_t A24minusi, f2elm_t A24plus, f2elm_t A24plusi, f2elm_t* coeff);

// Computes the 3-isogeny R=phi(X:Z), given projective point (X3:Z3) of order 3 on a Montgomery curve and a point P with coefficients given in coeff.
void eval_3_isog(point_proj_t Q, const f2elm_t* coeff);

// 3-way simultaneous inversion
void inv_3_way(f2elm_t z1, f2elm_t z1i, f2elm_t z2, f2elm_t z2i, f2elm_t z3, f2elm_t z3i);

// Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
void get_A(const f2elm_t xP, const f2elm_t xPi, const f2elm_t xQ, const f2elm_t xQi, const f2elm_t xR, const f2elm_t xRi, f2elm_t A, f2elm_t Ai);


#endif
