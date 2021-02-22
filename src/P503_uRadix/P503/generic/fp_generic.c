/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: portable modular arithmetic for P503
*********************************************************************************************/

#include "../P503_internal.h"
#include<stdio.h>

void digit_x_digit(const digit_t a, const digit_t b, digit_t* c)
{// Digit multiplication, digit * digit -> 2-digit result   
 // Inputs : 64bits a & b
 // Output : 128bits c = a * b , c is divided into 2 x 64bits parts 
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4), k = sizeof(digit_t) * 4;

    al = a & mask_low;                        // Low part
    ah = a >> k;                              // High part
    bl = b & mask_low;
    bh = b >> k;

    albl = al*bl;
    albh = al*bh;
    ahbl = ah*bl;
    ahbh = ah*bh;
    c[0] = albl & mask_low;                   // C00

    res1 = albl >> k;
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;  
    temp = res1 + res2 + res3;
    carry = temp >> k;
    c[0] ^= temp << k;                        // C01   
	
    res1 = ahbl >> k;
    res2 = albh >> k;
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    c[1] = temp & mask_low;                   // C10 
    carry = temp & mask_high; 
    c[1] ^= (ahbh & mask_high) + carry;       // C11
}

void digit_x_digit_2(const digit_t a, const digit_t b, digit_t* co)
{// Digit multiplication, digit * digit -> 2-digit result   
 // Inputs : 60bits a & b
 // Output : 120bits c = a * b , c is divided into 2 x 60bits parts
 // log(a)=log(b)=51 (21+30->30+30) make the numbers presented with 60 bits to compute the arithmetics
    register digit_t al, ah, bl, bh, temp;  
    digit_t albl, albh_ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> 34, mask_high = (digit_t)(-1) << 30;
	digit_t k = 30; 
	digit_t h = sizeof(digit_t) * 2;

    al = a & mask_low;    // Low part
    ah = a >> k;          // High part
    bl = b & mask_low;
    bh = b >> k;

    albl = al*bl;
	ahbh = ah*bh;
	albh_ahbl = (al+ah)*(bl+bh) - albl - ahbh;
    
    co[0] = albl & mask_low;            // C00

    res1 = albl >> k;
    temp = res1 + albh_ahbl;
    carry  = temp >> k;
    co[0] ^= (temp & mask_low) << 30;   // C01

    co[1] = ahbh + carry;
}

void mp_mul_old(const digit_t* a, const digit_t* b, digit_t* c, const unsigned int nwords)
{// Multiprecision comba multiply, c = a*b, where lng(a) = lng(b) = nwords.   
    unsigned int i, j;
    digit_t t = 0, u = 0, v = 0, UV[2];
    unsigned int carry = 0;
    
    for (i = 0; i < nwords; i++) {
        for (j = 0; j <= i; j++) {
            MUL(a[j], b[i-j], UV+1, UV[0]); 
            ADDC(0, UV[0], v, carry, v); 
            ADDC(carry, UV[1], u, carry, u); 
            t += carry;
        }
        c[i] = v;
        v = u; 
        u = t;
        t = 0;
    }

    for (i = nwords; i < 2*nwords-1; i++) {
        for (j = i-nwords+1; j < nwords; j++) {
            MUL(a[j], b[i-j], UV+1, UV[0]); 
            ADDC(0, UV[0], v, carry, v); 
            ADDC(carry, UV[1], u, carry, u); 
            t += carry;
        }
        c[i] = v;
        v = u; 
        u = t;
        t = 0;
    }

    c[2*nwords-1] = v; 
}

void improved_BR(const digit_t c_0, const digit_t c_1, digit_t* q, digit_t* r)        
{// Improved Barrett reduction (IBR) 
 // Inputs : 0≤c<2^(2N+lamda), c is divided into 2x60bits parts, input c = {c_1,c_0};
 // Output : q = floor(c / p) , r = c mod p , for p = 2^N_1 * R_3
 // N_1 = 25, N_2=log(R_3) = 26, R_3=3^16=0290D741   //log(t) -> 107 -> 82 (input) -> 58
    digit_t res1, res2, res3[2], carry, s;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*5), mask_low_1 = (digit_t)(-1) >> 37, mask_low_2 = (digit_t)(-1) >> (sizeof(digit_t)*2);
	digit_t R_3 = 0x000000000290D741;                // R_3=3^16=0290D741
	digit_t k = sizeof(digit_t) * 3;
	// N + gamma + 3 = 59
	// High part b = 2^108 / P  -->  031E 32483C 518FAE
	digit_t b = 0x031E32483C518FAE;

	s = c_0 & ((digit_t)(-1) >>39);      // low 25 bits

	res1 = c_0 >> 49;                    // N_1+N_2-2=49
    res2 = c_1 << 11;                    // width of c_0 is 60 --> 60 - 49 = 11
	res2 = res2 ^ res1;
	digit_x_digit_2(res2, b, res3);
	
	//r_h = t - q*R_3
	q[0] = (res3[1]<<1) ^ (res3[0]>>59); // N + gamma + 3 = 59 --> the reset of res3 is 1 bit
	res1 = q[0] & mask_low_1;            // low 27 bits -> N2+1
	res1 = res1 * R_3;                   // compute t1
	r[0] = (((c_0>>25) & mask_low_1) - (res1 & mask_low_1)) & mask_low_1;

	res1 = r[0] - R_3;
	res2 = res1 >> 63;
	carry = 1 - res2;
	r[0] = res1 + (R_3 * res2);

	// if(res2==0)
	// 	r[0] = res1;
	// else
	// 	r[0] = r[0];	

	res1 = r[0] << 25;
	r[0] = res1 ^ s;
	q[0] = q[0] + carry;	
}

void improved_IFFM(digit_t* c_0, digit_t* c_1, digit_t* C, const unsigned int nwords)
{// Apply the IBR function & Adjust all of the coefficients into the standard ranges
 // Inputs : the result from mp_mul_new4 or mp_sqr_new, raw coefficients c_0[i], c_1[i] are all presented in 60 bits ,for 0≤i≤n-1;
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1]
	unsigned int i;
    digit_t t[2], s = 0, q[1]={0}, r[1], temp, t1, s1;
	digit_t mask_low = (digit_t)(-1) >> 39;  // low 25 bits
	digit_t mask_low_1 = (digit_t)(-1) >> 4; // low 60 bits
	digit_t R_div_3 = 0x0001B5E4D6000000;    // R/3
	digit_t R = 0x000521AE82000000;          // R = 2^25 * 3^16
	digit_t k = 60;
	digit_t h = nwords - 1;

	improved_BR(c_0[8], c_1[8], q, r);	
    C[8] = r[0];                             //Apply the IBR function to the n coefficients serially , start at num 8

	temp = q[0] + c_0[9];
	c_0[9] = temp & mask_low_1;
	s = temp >> k;
	c_1[9] = c_1[9] + s;		
	improved_BR(c_0[9], c_1[9], q, r);	
    C[9] = r[0];	

	t1 = C[h] - R_div_3;
	s1 = t1>>63;
	// s = 1 - s1;
	s = 1 ^ s1;
	i = s;
	
	C[h] = t1 + R_div_3 * s1; // slow!
	// if(s1==0)
	// 	C[h] = t1;
	// else 
	// 	C[h] = C[h];
	
	t1 = C[h] - R_div_3;
	s1 = t1>>63;
	// s = 1 - s1;
	s = 1 ^ s1;
	
	C[h] = t1 + R_div_3 * s1; // slow!
	// if(s1==0)
	// 	C[h] = t1;
	// else 
	// 	C[h] = C[h];
	
	s += i;                                   //Adjust the coefficients into the standard ranges

	temp = q[0] + (q[0]<<1) + c_0[0] + s;    

	c_0[0] = temp & mask_low_1;
	s = temp >> k;
	c_1[0] = c_1[0] + s;
	improved_BR(c_0[0], c_1[0], q, r);	
	C[0] = r[0];
	
	for (i = 1; i < h-1; i++){
		temp = q[0] + c_0[i];
		c_0[i] = temp & mask_low_1;
		s = temp >> k;
		c_1[i] = c_1[i] + s;
		improved_BR(c_0[i], c_1[i], q, r);	
        C[i] = r[0];		
	}                                        //Apply the IBR function to the n coefficients serially
	
	c_0[8] = q[0] + C[8];
	c_1[8] = 0;
	improved_BR(c_0[8], c_1[8], q, r);	
	C[8] = r[0];
	
	C[h] += q[0];
	
	t1 = C[h] - R_div_3;
	s1 = t1>>63;

	C[h] = t1 + R_div_3 * s1; // slow!

	// if(s1==0)
	// 	C[h] = t1;
	// else 
	// 	C[h] = C[h];

	C[0] += 1 - s1;                         //Adjust the coefficients into the standard ranges
}

void digit_a_digit_1(const digit_t* a, const digit_t* b, digit_t* c)
{// Digit addition, 2digit + 2digit -> 2-digit+1bit  
 // Inputs : 102 bits a & b , a & b are divided into 2 x 60bits parts
 // Output : 103 bits c = a + b, c is divided into 2 x 60bits parts
    register digit_t temp;
	digit_t  albl, res;
	digit_t mask_low_1 = (digit_t)(-1) >> 4, k = 60;
	
	albl = a[0] + b[0];
	c[0] = albl & mask_low_1;
	res  = albl >> k;  
	c[1] = a[1] + b[1] + res;
}

void digit_s_digit_1(const digit_t* a, const digit_t* b, digit_t* c)
{// Digit subtraction, 2digit - 2digit -> 2-digit+1bit  a-b = c   
 // Inputs : 102 bits a & b , a & b are divided into 2 x 60bits parts
 // Output : 103 bits c = a - b, c is divided into 2 x 60bits parts
    digit_t albl, res;
	digit_t mask_low_1 = (digit_t)(-1) >> 4, k = 63;
	
	albl = a[0] - b[0];
	c[0] = albl & mask_low_1;
	res  = albl >> k;  
	c[1] = a[1] - b[1] - res;
}

void mp_mul_new4(const digit_t* a, const digit_t* b, digit_t* c_0, digit_t* c_1, const unsigned int nwords)
{// Compute multiply-accumulate terms and get the raw coefficients  
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], a[i] & b[i] are presented in 60 bits (51bits only)
 // Output : raw coefficients c_0[i], c_1[i] are all presented in 60 bits ,for 0≤i≤n-1;
	unsigned int i, j, k, s;
    digit_t UV_2[2]={0}, UV_3[2]={0}, p = 60;
	digit_t ci[nwords][2];
	unsigned int carryOut;
	digit_t mask_low = (digit_t)(-1) >> 4;    //low 60 bits
	digit_t temp;

	for (i = 0; i < nwords; i++) {
		digit_x_digit_2(a[i], b[i], ci[i]);
	} 

	c_0[0] = ci[0][0];
    c_1[0] = ci[0][1];
	
	for(i = 2; i < nwords; i = i + 2){
		k = i>>1;
		for(j=0;j<k;j++){
			s = i-j;
			digit_x_digit_2(a[j]+a[s], b[j]+b[s], UV_2);
			digit_s_digit_1(UV_2,ci[j],UV_2);
			digit_s_digit_1(UV_2,ci[s],UV_2);
			digit_a_digit_1(UV_2, UV_3, UV_3);          //one-level Karatsuba to reduce the complexity in the modular multiplication
		}  
		digit_a_digit_1(ci[k],UV_3,UV_3);
		c_0[i] = UV_3[0];
		c_1[i] = UV_3[1];
		UV_3[0] = 0;
		UV_3[1] = 0;
	}   //e.g.  c4 = a4*b0 + a3*b1 + a2*b2 + a1*b3 + a0*b4; c is divided into 2 x 60 bits parts

	for(i = 1; i < nwords; i = i + 2){
		for(j=0;j<=((i-1)>>1);j++){
			s = i-j;
			digit_x_digit_2(a[j]+a[s], b[j]+b[s], UV_2);
			digit_s_digit_1(UV_2,ci[j],UV_2);
			digit_s_digit_1(UV_2,ci[s],UV_2);
			digit_a_digit_1(UV_2, UV_3, UV_3);          //one-level Karatsuba to reduce the complexity in the modular multiplication
		}
		c_0[i] = UV_3[0];
		c_1[i] = UV_3[1];
		UV_3[0] = 0;
		UV_3[1] = 0;
	}   //e.g.  c5 = a5*b0 + a4*b1 + a3*b2 + a2*b3 + a1*b4 + a0*b5; c is divided into 2 x 60 bits parts

	for (i = nwords; i < 2*nwords-2; i = i + 2) {
		k = i>>1;
		for (j = i-nwords+1; j < k; j++) {
			s = i-j;
			digit_x_digit_2(a[j]+a[s], b[j]+b[s], UV_2);
			digit_s_digit_1(UV_2,ci[j],UV_2);
			digit_s_digit_1(UV_2,ci[s],UV_2);
			digit_a_digit_1(UV_2, UV_3, UV_3);          //one-level Karatsuba to reduce the complexity in the modular multiplication
        }
		digit_a_digit_1(ci[k],UV_3,UV_3);               //e.g.  c14 = a9*b5 + a8*b6 + a7*b7 + a6*b8 + a5*b9 ;c is divided into 2 x 60 bits parts
		s = i-nwords;
		temp = (UV_3[0]<<1) + UV_3[0] + c_0[s];	
		c_0[s] = temp & mask_low;
		carryOut = temp >> p;
		c_1[s] = (UV_3[1]<<1) + UV_3[1] + c_1[s] + carryOut;
		UV_3[0] = 0;
		UV_3[1] = 0;
	}   // get the raw coefficients.

	for (i = nwords + 1; i < 2*nwords-2; i = i + 2) {
		for (j = i-nwords+1; j <= ((i-1)>>1); j++) {
			s = i - j;
			digit_x_digit_2(a[j]+a[s], b[j]+b[s], UV_2);
			digit_s_digit_1(UV_2,ci[j],UV_2);
			digit_s_digit_1(UV_2,ci[s],UV_2);
			digit_a_digit_1(UV_2, UV_3, UV_3);          //one-level Karatsuba to reduce the complexity in the modular multiplication
        }                                               //e.g.  c15 = a9*b6 + a8*b7 + a7*b8 + a6*b9 ;c is divided into 2 x 60 bits parts
		s = i-nwords;
		temp = (UV_3[0]<<1) + UV_3[0] + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> p;
		c_1[s] = (UV_3[1]<<1) + UV_3[1] + c_1[s] + carryOut;
		UV_3[0] = 0;
		UV_3[1] = 0;
	}   // get the raw coefficients.

	temp = (ci[nwords-1][0]<<1) + ci[nwords-1][0] + c_0[nwords-2];
	c_0[nwords-2] = temp & mask_low;
	carryOut = temp >> p;
	c_1[nwords-2] = (ci[nwords-1][1]<<1) + ci[nwords-1][1] + c_1[nwords-2] + carryOut; // get the raw coefficients.
}

void fp_add_new(digit_t* a, digit_t* b, digit_t* c, const unsigned int nwords)
{// Modular Addition based on the new data representation , with the unconditional radix
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], a[i] & b[i] are presented in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 60 bits , to compute c = a + b
	unsigned int i,j;
	digit_t ca = 0, ca1, t;
	digit_t R  = 0x000521AE82000000;      // R = 2^25 * 3^16
	digit_t R_div_3 = 0x0001B5E4D6000000; // R/3
 
    for(i = 0; i < nwords; i++){
		c[i] = a[i] + b[i];
	}

	for (i = 1; i < nwords; i++){
		t = c[i-1] - R;
		ca1 = t>>63;
		// ca = 1 - ca1;
		ca = 1 ^ ca1;
		c[i-1] = t  + R * ca1;
		c[i] = c[i] + ca;
	} 

	t = c[nwords-1] - R_div_3;
	ca1 = t >> 63;
	// ca = 1 - ca1;
	ca = 1 ^ ca1;
	c[nwords-1] = t + R_div_3 * ca1;
	c[0] += ca;
}

void fpneg_new(digit_t* a,digit_t* b,const unsigned int nwords)
{// Modular Negation based on the new data representation , with the unconditional radix
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], a[i] all in 60 bits
 // Output : b[i] in [0, R-1] for 0≤i<n-1; b[n-1] in [0, R/3 - 1], b[i] all in 60 bits, to compute b = -a
	unsigned int i;
	digit_t R = 0x000521AE82000000;       // R = 2^25 * 3^16
	digit_t R_div_3 = 0x0001B5E4D6000000; // R/3
	digit_t k = nwords - 1;
	
	for( i = 0; i < k; i++){
		b[i] = R - a[i] - 1;
	}

	b[k] =  R_div_3 - a[k] - 1;
}

void fpsub_new(digit_t* a, digit_t* b, digit_t* c, const unsigned int nwords)
{// Modular Subtraction based on the new data representation , with the unconditional radix
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], a[i] & b[i] are presented in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 60 bits , to compute c = a - b
	unsigned int i, j = 0, t;
	digit_t R = 0x000521AE82000000;       // R = 2^25 * 3^16
	digit_t R_div_3 = 0x0001B5E4D6000000; // R/3
	
	for(i=0;i<nwords;i++){
		c[i] = a[i] - b[i];
	}
	// t = (c[0] == 0);
	// // t = is_digit_zero_ct(c[0]);
	// c[1] = c[1] - t;
	// c[0] = R * t + c[0] * (1 ^ t);
	
	// if(c[0] == 0){
	// 	c[0] = R;
	// 	c[1] = c[1] - 1;
	// }
	// else {
	// 	c[0] = c[0];
	// 	c[1] = c[1] - 0;
	// }
	t = (c[0] == 0);
	j = c[0] >> 63;
	t = t | j;
	c[0] = c[0] + t * R;
	c[1] = c[1] - t;
	
	for(i=1;i<nwords-1;i++){
		j = c[i]>>63;
		c[i] = c[i] + j * R;
		c[i+1] = c[i+1] - j;
	}
	
	j = c[nwords-1]>>63;
	c[nwords-1] = c[nwords-1] + j * R_div_3;
	c[0] = c[0] - j;
}

void fpdiv2_503_new(digit_t* a, digit_t* c)
{// Modular Division by Two based on the new data representation , with the unconditional radix
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], a[i] all in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], c[i] all in 60 bits, to compute c = a/2
    unsigned int i;
    digit_t mask;
	digit_t R_div_2 = 0x000290D741000000; //R/2
	digit_t R_div_6 = 0x0000DAF26B000000; //R/3/2
	digit_t b[NWORDS_FIELD] ={0};

	for (i = 1; i < NWORDS_FIELD; i++) {
		mask = (digit_t)(a[i] & 1);		
		a[i] -= mask;
		b[i-1] = R_div_2 & ((digit_t)(0-mask));
		a[i] >>= 1;
	}

	mask = (digit_t)(a[0] & 1);		
	a[0] -= mask;
	b[NWORDS_FIELD-1] = R_div_6 & ((digit_t)(0-mask));
	a[0] >>= 1;

    for (i = 0; i < NWORDS_FIELD; i++) {
		c[i] = a[i] + b[i];
	}
}

void fpmul_new(const digit_t* a, digit_t* b, digit_t* c)
{// Modular Multipication based on the new data representation , with the unconditional radix
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], all in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 60 bits , to compute c = a * b
    digit_t c1_0[10];
	digit_t c1_1[10];

	mp_mul_new4(a, b, c1_0, c1_1, 10);
	improved_IFFM(c1_0, c1_1, c, 10);
}

void digit_sqr(const digit_t a, digit_t* co)
{// Digit multiplication, digit * digit -> 2-digit result     
 // Inputs : 60bits a & b (51bits at most)
 // Output : 120bits c = a * b , c is divided into 2 x 60bits parts
 // log(a)=log(b)=51 (21+30->30+30) make the numbers presented with 60 bits to compute the arithmetics
    register digit_t al, ah, temp;  
    digit_t alal, alah, ahah, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> 34, mask_high = (digit_t)(-1) << 30;
	digit_t k = 30; digit_t h = sizeof(digit_t) * 2;

    al = a & mask_low;    // Low part
    ah = a >> k;          // High part

    alal = al*al;
	ahah = ah*ah;
	alah = al*ah;
    
    co[0] = alal & mask_low;                         // C00

    res1 = alal >> k;
    temp = res1 + (alah<<1);
    carry = temp >> 30;
    co[0] ^= (temp & mask_low) << 30;                // C01

    co[1] = ahah + carry;
}

void mp_sqr_new(const digit_t* a, digit_t* c_0, digit_t* c_1, const unsigned int nwords)
{// Compute multiply-accumulate terms and get the raw coefficients  
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], all in 60 bits (51bits at most)
 // Output : raw coefficients c_0[i], c_1[i] are all in 60 bits ,for 0≤i≤n-1; 
    unsigned int i, j;
    digit_t UV[2], UV_2[2] = {0},temp;
    digit_t carryOut, k, s;
	digit_t h = 60;
    digit_t mask_low = (digit_t)(-1) >> 4; //low 60 bits
	
	digit_sqr(a[0], UV);
    c_0[0] = UV[0];
    c_1[0] = UV[1];
	
	for (i = 2; i < nwords; i = i + 2) {
		k = i>>1;
		for (j = 0; j < k; j++) {
             digit_x_digit_2(a[j], a[i-j], UV);	
             digit_a_digit_1(UV, UV_2, UV_2);
        }
		digit_a_digit_1(UV_2, UV_2, UV);
		digit_sqr(a[k], UV_2);	
		digit_a_digit_1(UV, UV_2, UV_2);
		c_0[i] = UV_2[0];
		c_1[i] = UV_2[1];
		UV_2[0] = 0;
		UV_2[1] = 0;
	}	//e.g. c4 = 2 * a0 * a4 + 2 * a1 * a3 + a2 ^ 2; c is divided into 2 x 60 bits parts
	
	for (i = 1; i < nwords; i = i + 2) {
		 for (j = 0; j <= (i-1)/2; j++) {
             digit_x_digit_2(a[j], a[i-j], UV);	
             digit_a_digit_1(UV, UV_2, UV_2);
        }
		digit_a_digit_1(UV_2, UV_2, UV);
		c_0[i] = UV[0];
		c_1[i] = UV[1];
		UV_2[0] = 0;
		UV_2[1] = 0;
	}	//e.g. c5 = 2 * a0 * a5 + 2 * a1 * a4 + 2 * a2 * a3 ;c is divided into 2 x 60 bits parts
	
	for (i = nwords; i < 2*nwords-2; i = i + 2) {
		k = i>>1;
		for (j = i-nwords+1; j < k; j++) {
			digit_x_digit_2(a[j], a[i-j], UV);			
            digit_a_digit_1(UV, UV_2, UV_2);
        }
		digit_a_digit_1(UV_2, UV_2, UV);
		digit_sqr(a[k], UV_2);	
		digit_a_digit_1(UV, UV_2, UV_2);                     //e.g. c14 = 2* a9 * a5 + 2 * a8 * a6 + a7 ^ 2;c is divided into 2 x 60 bits parts
		s = i-nwords;
		temp = (UV_2[0]<<1) + UV_2[0] + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> h;
		c_1[s] = (UV_2[1]<<1) + UV_2[1] + c_1[s] + carryOut;
		UV_2[0] = 0;
		UV_2[1] = 0;
	}   // get the raw coefficients.
	
	for (i = nwords + 1; i < 2*nwords-2; i = i + 2) {
		for (j = i-nwords+1; j <= (i-1)/2; j++) {
			digit_x_digit_2(a[j], a[i-j], UV);			
            digit_a_digit_1(UV, UV_2, UV_2);
        }
		digit_a_digit_1(UV_2, UV_2, UV);                     //e.g. c15 = 2 * a9 * a6 + 2 * a8 * a7;
		s = i-nwords;
		temp = (UV[0]<<1) + UV[0] + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> h;
		c_1[s] = (UV[1]<<1) + UV[1] + c_1[s] + carryOut;
		UV_2[0] = 0;
		UV_2[1] = 0;
	}   // get the raw coefficients.

	digit_sqr(a[nwords-1], UV_2);
	temp = (UV_2[0]<<1) + UV_2[0] + c_0[nwords-2];
	c_0[nwords-2] = temp & mask_low;
	carryOut = temp >> h;
	c_1[nwords-2] = (UV_2[1]<<1) + UV_2[1] + c_1[nwords-2] + carryOut; // get the raw coefficients.
}

void fpsqr_new(const digit_t* a, digit_t* c)
{// Modular Squaring based on the new data representation , with the unconditional radix
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], all in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 60 bits , to compute c = a^2
	digit_t c1_0[10];
	digit_t c1_1[10];

	mp_sqr_new(a, c1_0, c1_1, 10);   //performed better than mp_mul_new4
	improved_IFFM(c1_0, c1_1, c, 10);
}

void fpinv_chain_new(digit_t* a , const unsigned int nwords)
{// Modular Inversion
	unsigned int i, j;
	digit_t tt[nwords];
	digit_t t0[nwords], t1[nwords], t2[nwords], t3[nwords], t4[nwords], t5[nwords], t6[nwords], t7[nwords], t8[nwords];
	digit_t t9[nwords], t10[nwords],t11[nwords],t12[nwords],t13[nwords],t14[nwords];

    fpsqr_new(a, tt);
    fpmul_new(a, tt, t0);
	fpmul_new(t0, tt, t1);
	fpmul_new(t1, tt, t2);
	fpmul_new(t2, tt, t3);
	fpmul_new(t3, tt, t4);
	fpmul_new(t4, tt, t5);
	fpmul_new(t5, tt, t6);
	fpmul_new(t6, tt, t7);
	fpmul_new(t7, tt, t8);
	fpmul_new(t8, tt, t9);
	fpmul_new(t9, tt, t10);
	fpmul_new(t10, tt, t11);
	fpmul_new(t11, tt, t12);
	fpmul_new(t12, tt, t13);
	fpmul_new(t13, tt, t14);

	for (i = 0; i < nwords; i++)
		tt[i] = a[i];
	
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(a, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t8, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t6, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t9, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t0, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(a, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t6, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t2, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t8, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(a, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t10, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t0, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t10, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t10, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t5, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t2, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t6, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t3, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t5, tt, tt);
    for (i = 0; i < 12; i++) fpsqr_new(tt, tt);
    fpmul_new(t12, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t8, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t6, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t12, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t11, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t6, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t5, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t14, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t14, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t5, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t6, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t8, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(a, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t4, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t6, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t5, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t7, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(a, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t0, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t11, tt, tt);
    for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
    fpmul_new(t13, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t1, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t10, tt, tt);
    for (j = 0; j < 49; j++) {
        for (i = 0; i < 5; i++) fpsqr_new(tt, tt);
        fpmul_new(t14, tt, tt);
    }

	for (i = 0; i < nwords; i++)
		a[i] = tt[i];
}

void fpinv_new(digit_t* a, digit_t* c, const unsigned int nwords)
{// Modular Inversion for p503
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], all in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 60 bits
    digit_t tt[nwords];
	unsigned int i;

    for (i = 0; i < nwords; i++)
		tt[i] = a[i];
	
    fpinv_chain_new(tt, nwords);
    fpsqr_new(tt, tt);
    fpsqr_new(tt, tt);
    fpmul_new(a, tt, c);
}

void digit_a_digit_2(digit_t a, digit_t b, digit_t* c)
{// Digit addition, digit + digit -> 1-digit+1bit     
 // Inputs : 64 bits a & b
 // Output : 65 bits c = a + b , c[0] c[1]
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, res;
	unsigned int carryOut;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), k = sizeof(digit_t) * 4;
	
    al = a & mask_low;    // Low part
    ah = a >> k;          // High part
    bl = b & mask_low;
    bh = b >> k;
    albl = al+bl;
    c[0] = albl & mask_low;                   
    res = albl >> k;      //carryout
    temp = res + ah + bh;
    res = temp >> k;      //carryout
	c[0] ^= temp << k; 
	c[1] = res;
}


void from_uRadix(digit_t* a, digit_t* co, const unsigned int nwords)
{// From Unconventional Radix Back to Normal
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], all in 60 bits based on the new data representation , with the unconditional radix
 // Output : c in F(p) , presented in normal presentation 
	digit_t R = 0x000521AE82000000;            // R = 2^25 * 3^16
	digit_t R_sub_1 = 0x000521AE81FFFFFF;      // R - 1
	digit_t R_div_3_sub1 = 0x0001B5E4D5FFFFFF; // R/3 - 1
	unsigned int i,j,k;
	digit_t c1[10][2] = {0}; digit_t c[10] = {0};
	digit_t c0[2] = {0}, c01[2] = {0};
	digit_t c2[2] = {0}, c3[2] = {0}, c4[2] = {0};
	digit_t ca = 0;
	digit_t c_next[nwords];
	
	j = 1;
	for(i=1;i<nwords-1;i++){
		// if (a[i] == R_sub_1)
		// 	k = 1;
		// else
		// 	k = 0;

		k = (a[i] == R_sub_1);
		
		j = j & k;
	}

	// if(a[nwords-1] == R_div_3_sub1)
	// 	k = 1;
	// else
	// 	k = 0;

	k = (a[nwords-1] == R_div_3_sub1);
	
	j = j & k;
	
	// if ((a[0] == R)|(a[0] == R_sub_1))
	// 	k = 1;
	// else 
	// 	k = 0;
	k = ((a[0] == R)|(a[0] == R_sub_1));
	
	j = j & k;
	
	k = (digit_t)(a[0] & 1); 
	
	// if(j == 1){
	// 	for(i=1;i<nwords;i++){
	// 		a[i] = 0;
	// 	}
		
	// 	if(k == 1)
	// 		a[0] = 0;
	// 	else
	// 		a[0] = 1;
	// }
	// else{
	// 	for(i=1;i<nwords;i++){
	// 		a[i] = a[i];
	// 	}
		
	// 	if(k == 1)
	// 		a[0] = a[0];
	// 	else
	// 		a[0] = a[0];
	// }
	//Adjust the coefficients into the standard ranges

	for(i=1;i<nwords;i++){
			a[i] = a[i] * (1 ^ j);
	}
	a[0] = a[0] * (1 ^ j) + (1 ^ k) * j;

	digit_x_digit(a[nwords - 1], R, c0);
	digit_a_digit_2(c0[0], a[nwords - 2], c01);
	c[0] = c01[0];
	c[1] = c0[1] + c01[1];

	for(i=nwords-2;i>0;i--){
		digit_x_digit(c[0], R, c3);
		digit_a_digit_2(a[i-1], c3[0], c2);
		c_next[0] = c2[0];
		ca = c2[1];

		for(j=0;j<nwords-i;j++){
			digit_x_digit(c[j], R, c1[j]);
		}

		for(j=0;j<nwords-1-i;j++){
			digit_a_digit_2(c1[j][1], c1[j+1][0], c3);
			digit_a_digit_2(c3[0], ca, c4);
			c_next[j+1] = c4[0];
			ca = c4[1] + c3[1];
		}

		digit_a_digit_2(c1[nwords-1-i][1], ca, c3);
		c[nwords-i] = c3[0];

		for(k=0;k<=nwords-1-i;k++){
			c[k] = c_next[k];
		}
	}
	
	for(i=0;i<8;i++){
		co[i] = c[i];
	}
}

void fp2mul_uRadix(f2elm_t a, f2elm_t b, f2elm_t c)                             
{// GF(p^2) multiplication, c = a*b in GF(p^2).
 // Inputs : a = a0 + a1*i, b = b0 + b1*i. a0, a1, b0, b1 all meet the same requirement with fpmul_new & fpsqr_new
 // Output : c = c0 + c1*i. c = a * b, c0 and c1 meet the same requirement with fpmul_new & fpsqr_new
 // a-b+2^5*R*P mod P = (a_9-b_9+2^5*(R/3-1)*R)*R^9 + (a_8-b_8+2^5*(R-1)*R)*R^8 + ... + a_0-b_0+2^5*(R-1)*R
    digit_t PR_9[2] = {0x04DBCA2FC0000000, 0x0000118E2A5AFB51}; // 2^5*(R/3-1)*R) //      0118E2A5AFB51 4DBCA2FC0000000
    digit_t PR[2] =   {0x0FDBCA2FC0000000, 0x000034AA7F10F1F3}; // 2^5*(R-1)*R    //      034AA7F10F1F3 FDBCA2FC0000000
    digit_t t1[10], t2[10];
    digit_t tt1[2][10], tt2[2][10], tt3[2][10]; 
	int i, j;
	digit_t temp, res;
	digit_t mask_low_1 = (digit_t)(-1) >> 4;
	digit_t k = 60;
	
	for (i=0; i<10; i++){
		t1[i] = a[0][i] + a[1][i];
		t2[i] = b[0][i] + b[1][i];
	}
	
	//for c[0]
    mp_mul_new4(a[0], b[0], tt1[0], tt1[1], 10); // tt1 = a0*b0	
	mp_mul_new4(a[1], b[1], tt2[0], tt2[1], 10); // tt2 = a1*b1

	for (i=0; i<9; i++) {
		//a0*b0 + 2^5*R*P - a1*b1
		temp = tt1[0][i] + PR[0] - tt2[0][i];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> k;
		tt3[1][i] = tt1[1][i] + PR[1] - tt2[1][i] + (res&(digit_t)1) - (res&(digit_t)2);
	}
	
	//a0*b0 + 2^5*R*P - a1*b1
	temp = tt1[0][9] + PR_9[0] - tt2[0][9];
	tt3[0][9] = temp & mask_low_1;
	res = temp >> k;
	tt3[1][9] = tt1[1][9] + PR_9[1] - tt2[1][9] + (res&(digit_t)1) - (res&(digit_t)2);
	
	improved_IFFM(tt3[0], tt3[1], c[0], 10);// c[0]
	
	//for c[1]
	
	mp_mul_new4(t1, t2, tt3[0], tt3[1], 10); // tt1 = (a0+a1)*(b0+b1)	
	
	// c[1] = (a0+a1)*(b0+b1) - (a0*b0 + a1*b1) 
	for (i=0; i<10; i++) {
		//a0*b0 + a1*b1
		res = 0;
		temp = tt2[0][i] + tt1[0][i] + res;
		tt2[0][i] = temp & mask_low_1;
		res = temp >> k;
		tt2[1][i] = tt2[1][i] + tt1[1][i] + res;
		//(a0+a1)*(b0+b1) - a0*b0 - a1*b1
		res = 0;
		temp = tt3[0][i] - tt2[0][i] - res;
		tt3[0][i] = temp & mask_low_1;
		res = temp >> (sizeof(digit_t) * 8 - 1);		
		tt3[1][i] = tt3[1][i] - tt2[1][i] - res;
	}
	
    improved_IFFM(tt3[0], tt3[1], c[1], 10); 
}

void IBR_uRadix(const digit_t* a, digit_t* q, digit_t* r, const unsigned int unwords, const digit_t* lamda)
{// IBR 
	digit_t t[unwords], s, t_r, temp[2*unwords], res1, res2, res3, mask_low_1 = (digit_t)(-1) >> (sizeof(digit_t)*4+5);//low 27 bits n2+1
	digit_t R_3 = 0x000000000290D741;           // R_3=3^16=0290D741
	unsigned int i, j;
	
	s = a[0] & ((digit_t)(-1) >> 39);
	t_r = a[0] >> 25;
	
	if(unwords == 9){
	     for (i = 0; i < unwords-2; i++){
		     t[i] = (a[i] >> 49) ^ (a[i+1] << 15);
	     }
	     t[unwords-2] = a[unwords -2] >> 49;
	     mp_mul_old(t, lamda, temp, (unwords-1));
	}
	else if(unwords == 4){
	     for (i = 0; i < unwords-1; i++){
		     t[i] = (a[i] >> 49) ^ (a[i+1] << 15);
	     }
	     t[unwords-1] = a[unwords -1] >> 49;
	     mp_mul_old(t, lamda, temp, unwords);
	}
	else if(unwords > 4){
	     for (i = 0; i < unwords-1; i++){
		     t[i] = (a[i] >> 49) ^ (a[i+1] << 15);
	     }
	     mp_mul_old(t, lamda, temp, unwords-1);
	}
	else {
	     for (i = 0; i < unwords; i++){
		     t[i] = (a[i] >> 49) ^ (a[i+1] << 15);
	     }
	     mp_mul_old(t, lamda, temp, unwords);
	}
	
	if(unwords>4){
		for(i = 0; i < unwords-1; i++){
			res1 = temp[unwords+i-2] >> (7+14*(9-unwords));
			res2 = temp[unwords+i-1] << (57-14*(9-unwords));
			q[i] = res1 ^ res2;
		}
		
		for(i = unwords -1; i<8; i++){
			q[i] = 0;
		}
	}
	else {
		for(i = 0; i < unwords; i++){
			res1 = temp[unwords+i-1] >> (13+14*(4-unwords));
			res2 = temp[unwords+i]  << (51-14*(4-unwords));
			q[i] = res1 ^ res2;
		}
		
		for(i = unwords; i<8; i++){
			q[i] = 0;
		}
	}
	
	res1 = q[0] & mask_low_1;
	res1 = res1 * R_3;
	r[0] = ((t_r & mask_low_1) - (res1 & mask_low_1)) & mask_low_1;
	
	res1 = r[0] - R_3;
	res2 = res1 >> 63;
	res3 = 1 - res2;
	
	r[0] = res1 + R_3 * res2;
	// if(res2==0)
	// 	r[0] = res1;
	// else 
	// 	r[0] = r[0];
	
	res1 = q[0] + res3;
	
	// if (q[0]==(digit_t)(-1)){
	// 	q[0] = 0;
	// 	j = 1;
	// }
	// else{
	// 	q[0] = res1;
	// 	j = 0;
	// }

	j = (q[0]==(digit_t)(-1));
	q[0] = res1 * (1 ^ j);
	
	int h;
	if(unwords>4){
		h=unwords-2;
	}
	else 
		h=unwords-1;
	
	
	for (i = 1; i < h; i++){	
		res1 = q[i] + j;	
		// if (q[i]==(digit_t)(-1)){
		// 	q[i] = 0;
		// 	j = 1;
		// }
		// else{
		// 	q[i] = res1;
		// 	j = 0;
		// }	
		j = (q[i]==(digit_t)(-1));
		q[i] = res1 * (1 ^ j);	
	}
	q[h] = q[h] + j;
	res1 = r[0] << 25;
	r[0] = res1 ^ s;
}

void to_uRadix(digit_t* a, digit_t* uc, const unsigned int nwords)
{// Conversion to uncoventional radix representation,
 // Inputs : a in F(p) normal field
 // Output : uc[i] in [0, R-1] for 0≤i<n-1; uc[n-1] in [0, R/3 - 1], uc[i] all in 60 bits based on the new data representation , with the unconditional radix
 // mc = a_9*R^9+a_8*R^8+...+a_1*R+a_0.
 // a_i= a mod R, a = (a/R). 
	int i, j;
	digit_t q[nwords], temp[2];
	digit_t lamda[9][8]={{0xbc23af523924bbf5,0x4dacb4c5fb3fe41d,0x9f6bf693516e84e1,0x3c900bf59cc18e4c,0x70cce2743634abb5,0xdbd44c71852bd663,0xe32483c518faeae5,0x0000000000000031},
	                   	 {0x2d317ecff9076f08,0xfda4d45ba138536b,0x02fd6730639327da,0x389d0d8d2aed4f24,0x131c614af598dc33,0x20f1463ebab976f5,0x00000000000c78c9},
						 {0x3516e84e14dacb4c,0x59cc18e4c9f6bf69,0x43634abb53c900bf,0x1852bd66370cce27,0x518faeae5dbd44c7,0x000000031e32483c},
						 {0x0639327dafda4d45,0xd2aed4f2402fd673,0xaf598dc33389d0d8,0xebab976f5131c614,0x0000c78c920f1463},						 
						 {0xb53c900bf59cc18e,0x6370cce2743634ab,0xe5dbd44c71852bd6,0x31e32483c518faea},						 
						 {0x33389d0d8d2aed4f,0xf5131c614af598dc,0xc920f1463ebab976,0x0000000000000c78},						 
						 {0xc71852bd66370cce,0x3c518faeae5dbd44,0x00000000031e3248},						 
						 {0x63ebab976f5131c6,0x000000c78c920f14},
						 {0x0031e32483c518fa}};

    for (i=0; i<nwords-1; i++){
		IBR_uRadix(a, q, temp, nwords-i-1, lamda[i]);
		uc[i] = temp[0];
		
		if(nwords-i-1>4){
		for (j=0; j<nwords-i-2; j++){
			a[j] = q[j];
		}
		for(j = nwords - i -2; j < 8; j++){
			a[j] = 0;
		}
		}
		else{
			for (j=0; j<nwords-i-1; j++){
			a[j] = q[j];
		    }
		    for(j = nwords - i -1; j < 8; j++){
			a[j] = 0;
		    }
		}
	}
	uc[nwords-1] = q[0];
}

void fp2sqr_uRadix(f2elm_t a, f2elm_t c)
{// GF(p^2) squaring, c = a^2 in GF(p^2).
 // Inputs : a = a0 + a1*i, a0 and a1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new
 // Output : c = c0 + c1*i, c = a^2; co & c1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new
    felm_t t1, t2, t3;
	int i, j;

	for (i=0; i<NWORDS_FIELD; i++){
		t1[i] = a[0][i] + a[1][i];
		t3[i] = a[0][i] << 1; 
	}
    fpsub_new (a[0], a[1], t2, NWORDS_FIELD);        // t2 = a0-a1
    fpmul_new(t1, t2, c[0]);                         // c0 = (a0+a1)(a0-a1)
    fpmul_new(t3, a[1], c[1]);                       // c1 = 2a0*a1
}

__inline void fp2add(f2elm_t a, f2elm_t b, f2elm_t c)           
{// GF(p^2) addition, c = a+b in GF(p^2).
 // Inputs : a = a0 + a1*i, b = b0 + b1*i, a0, a1, b0, b1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new
 // Output : c = c0 + c1*i, c = a + b; c0 & c1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new

    fp_add_new(a[0], b[0], c[0], NWORDS_FIELD);
	fp_add_new(a[1], b[1], c[1], NWORDS_FIELD);
}


__inline void fp2sub(f2elm_t a, f2elm_t b, f2elm_t c)          
{// GF(p^2) subtraction, c = a-b in GF(p^2).
 // Inputs : a = a0 + a1*i, b = b0 + b1*i, a0, a1, b0, b1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new
 // Output : c = c0 + c1*i, c = a - b; c0 & c1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new

    fpsub_new(a[0], b[0], c[0], NWORDS_FIELD);
	fpsub_new(a[1], b[1], c[1], NWORDS_FIELD);
}


void fp2neg(f2elm_t a)
{// GF(p^2) negation, a = -a in GF(p^2).
 // Inputs : a = a0 + a1*i, a0 and a1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new
 // Output : a = a0 + a1*i, a = -a, a0 & a1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new

    fpneg_new(a[0], a[0], NWORDS_FIELD);
	fpneg_new(a[1], a[1], NWORDS_FIELD);
}