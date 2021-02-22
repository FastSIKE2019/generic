/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: portable modular arithmetic for P751
*********************************************************************************************/

#include "../P751_internal.h"
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
 // Inputs : 64bits a & b
 // Output : 128bits c = a * b , c is divided into 3 x 48bits parts
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4);
	digit_t k = sizeof(digit_t) * 4; digit_t h = sizeof(digit_t) * 2;

    al = a & mask_low;                              // Low part
    ah = a >> k;                                    // High part
    bl = b & mask_low;                              
    bh = b >> k;

    albl = al*bl;
    albh = al*bh;
    ahbl = ah*bl;
    ahbh = ah*bh;
    co[0] = albl & mask_low;                        

    res1 = albl >> k;
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;  
    temp = res1 + res2 + res3;
    carry = temp >> k;
	co[1] = (temp & mask_low) >> h;
    co[0] ^= temp << k; 
	co[0] &= ((digit_t)(-1) >> h);                  //low 48 bits -- C0
	                                               
    res1 = ahbl >> k;
    res2 = albh >> k;
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    co[1] ^= ((temp & mask_low) << h);              // C1 
    carry = temp >> k;
    co[2] = ((ahbh & mask_high) >> k) + carry;      // C2
}

void digit_a_digit_1(const digit_t* a, const digit_t* b, digit_t* c)
{// Digit addition, 2digit + 2digit -> 2-digit+1bit   
 // Inputs : 128 bits a & b , a & b are divided into 3 x 48bits parts
 // Output : 129 bits c = a + b, c is divided into 3 x 48bits parts
    register digit_t temp;
	digit_t  albl , ambm , ahbh , res;
	digit_t mask_low_1 = (digit_t)(-1) >> (sizeof(digit_t)*2), k = sizeof(digit_t) * 6;
	
	albl = a[0] + b[0];
	c[0] = albl & mask_low_1;                       //low 48 bits
	res  = albl >> k;  
	
	ambm = a[1] + b[1] + res;
	c[1] = ambm & mask_low_1;
	res  = ambm >> k; 
	
	ahbh = a[2] + b[2] + res;
	c[2] = ahbh;
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

void mp_mul_new5(const digit_t* a, const digit_t* b, digit_t* c_0, digit_t* c_1, digit_t* ca, const unsigned int nwords)
{// Compute multiply-accumulate terms and get the raw coefficients  
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], a[i] & b[i] are 64 bits
 // Output : raw coefficients c_0[i], c_1[i], ca[i] are all 48 bits ,for 0≤i≤n-1;
    unsigned int i, j;
    digit_t UV[3], UV_2[3] = {0}, temp;
    digit_t carryOut, k;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*2);  //low 48 bits
	
	for (i = 0; i < nwords; i++) {
		 for (j = 0; j <= i; j++) {
            digit_x_digit_2(a[j], b[i-j], UV);
			digit_a_digit_1(UV, UV_2, UV_2);
        }
		c_0[i] = UV_2[0];
		c_1[i] = UV_2[1];
		ca[i]  = UV_2[2];
		UV_2[0] = 0;
		UV_2[1] = 0;
		UV_2[2] = 0;
	}	

	for (i = nwords; i < 2*nwords-1; i++) {
		for (j = i-nwords+1; j < nwords; j++) {
			digit_x_digit_2(a[j], b[i-j], UV);
			digit_a_digit_1(UV, UV_2, UV_2);
        }
		k = i-nwords;
		temp = (UV_2[0]<<1) + UV_2[0] + c_0[k];
		c_0[k] = temp & mask_low;
		carryOut = temp >> (sizeof(digit_t) * 6);
		temp = (UV_2[1]<<1) + UV_2[1] + c_1[k] + carryOut;
		c_1[k] = temp & mask_low;
		carryOut = temp >> (sizeof(digit_t) * 6);
		ca[k] = (UV_2[2]<<1) + UV_2[2] + ca[k] + carryOut;
		UV_2[0] = 0;
		UV_2[1] = 0;
		UV_2[2] = 0;
	}
}

void digit_s_digit_1(const digit_t* a, const digit_t* b, digit_t* c)//128bits - 128bits 
{// Digit subtraction, 2digit - 2digit -> 2-digit+1bit  
 // Inputs : 128 bits a & b , a & b are divided into 3 x 48bits parts
 // Output : 129 bits c = a - b, c is divided into 3 x 48bits parts
    digit_t albl, ahbh, ambm, res;
	digit_t mask_low_1 = (digit_t)(-1) >> (sizeof(digit_t)*2), k = sizeof(digit_t) * 8 - 1;
	
	albl = a[0] - b[0];
	c[0] = albl & mask_low_1;
	res  = albl >> k;             // get the sign bit
	
	ambm = a[1] - b[1] - res;
	c[1] = ambm & mask_low_1;
	res  = ambm >> k; 
	
	ahbh = a[2] - b[2] - res;
	c[2] = ahbh;
}

void mp_mul_new4(const digit_t* a, const digit_t* b, digit_t* c_0, digit_t* c_1, digit_t* ca, const unsigned int nwords)
{// Compute multiply-accumulate terms and get the raw coefficients  
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], a[i] & b[i] are 64 bits
 // Output : raw coefficients c_0[i], c_1[i], ca[i] are all 48 bits ,for 0≤i≤n-1;
	unsigned int i, j, k, s;
    digit_t UV[2]={0}, UV_1[2]={0}, UV_2[3]={0}, UV_3[3]={0}, p = sizeof(digit_t) * 6;
	digit_t ci[nwords][3];
	unsigned int carryOut;
	digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*2); //low 48 bits
	digit_t temp;

	for (i = 0; i < nwords; i++) {
		digit_x_digit_2(a[i], b[i], ci[i]);
	} 

	c_0[0] = ci[0][0];
    c_1[0] = ci[0][1];
    ca[0]  = ci[0][2];
	
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
		ca[i]  = UV_3[2];
		UV_3[0] = 0;
		UV_3[1] = 0;
		UV_3[2] = 0;
	}   //e.g.  c4 = a4*b0 + a3*b1 + a2*b2 + a1*b3 + a0*b4; c is divided into 3 x 48 bits parts

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
		ca[i]  = UV_3[2];
		UV_3[0] = 0;
		UV_3[1] = 0;
		UV_3[2] = 0;
	}   //e.g.  c5 = a5*b0 + a4*b1 + a3*b2 + a2*b3 + a1*b4 + a0*b5; c is divided into 3 x 48 bits parts

	for (i = nwords; i < 2*nwords-2; i = i + 2) {
		k = i>>1;
		for (j = i-nwords+1; j < k; j++) {
			s = i-j;
			digit_x_digit_2(a[j]+a[s], b[j]+b[s], UV_2);
			digit_s_digit_1(UV_2,ci[j],UV_2);
			digit_s_digit_1(UV_2,ci[s],UV_2);
			digit_a_digit_1(UV_2, UV_3, UV_3);           //one-level Karatsuba to reduce the complexity in the modular multiplication
        }
		digit_a_digit_1(ci[k],UV_3,UV_3);                //e.g.  c18 = a11*b7 + a10*b8 + a9*b9 + a8*b10 + a7*b11 ;c is divided into 3 x 48 bits parts
		temp = (UV_3[0]<<1) + UV_3[0] + c_0[i-nwords];   
		s = i-nwords;
		c_0[s] = temp & mask_low;
		carryOut = temp >> p;
		temp = (UV_3[1]<<1) + UV_3[1] + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> p;
		ca[s] = (UV_3[2]<<1) + UV_3[2] + ca[s] + carryOut;
		UV_3[0] = 0;
		UV_3[1] = 0;
		UV_3[2] = 0;
	}                                                    // get the raw coefficients.

	for (i = nwords + 1; i < 2*nwords-2; i = i + 2) {
		for (j = i-nwords+1; j <= ((i-1)>>1); j++) {
			s = i - j;
			digit_x_digit_2(a[j]+a[s], b[j]+b[s], UV_2);
			digit_s_digit_1(UV_2,ci[j],UV_2);
			digit_s_digit_1(UV_2,ci[s],UV_2);
			digit_a_digit_1(UV_2, UV_3, UV_3);           //one-level Karatsuba to reduce the complexity in the modular multiplication
        }                                                //e.g.  c19 = a11*b8 + a10*b9 + a9*b10 + a8*b11 ;c is divided into 3 x 48 bits parts
		s = i-nwords;
		temp = (UV_3[0]<<1) + UV_3[0] + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> p;//
		temp = (UV_3[1]<<1) + UV_3[1] + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> p;//
		ca[s] = (UV_3[2]<<1) + UV_3[2] + ca[s] + carryOut;
		UV_3[0] = 0;
		UV_3[1] = 0;
		UV_3[2] = 0;
	}                                                    // get the raw coefficients.

	temp = (ci[11][0]<<1) + ci[11][0] + c_0[10];
	c_0[10] = temp & mask_low;
	carryOut = temp >> p;
	temp = (ci[11][1]<<1) + ci[11][1] + c_1[10] + carryOut;
	c_1[10] = temp & mask_low;
	carryOut = temp >> p;
	ca[10] = (ci[11][2]<<1) + ci[11][2] + ca[10] + carryOut;   // get the final raw coefficients.
}

void improved_BR(const digit_t* t, const digit_t s, digit_t* q, digit_t* r)
{// Improved Barrett reduction (IBR) 
 // Inputs : t = floor(c / 2^N_1) , s = c mod 2^N_1 for 0≤c<2^(2N+lamda)
 // Output : q = floor(c / p) , r = c mod p , for p = 2^N_1 * R_3
 // N_1 = 31, N_2=log(R_3) = 32, R_3=3^20=CFD41B91   //log(t) -> 132 ->  101(input) -> 71
    register digit_t al, am, ah, bl, bm, bh, temp;
    digit_t albl, ambm, ahbh, albm_ambl, albh_ahbl, ambh_ahbm, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*5), mask_low_1 = (digit_t)(-1) >> (sizeof(digit_t)*4-1);
	digit_t R_3 = 0x00000000CFD41B91;                // R_3=3^20=CFD41B91
	digit_t k = sizeof(digit_t) * 3;
	
	res1 = t[0] >> 30;     // 2.  t >> (N2-2)
    al = res1 & mask_low;  // low 24bits
	res2 = t[1] << 10;
    am = (res2 & mask_low) ^ (res1 >> (sizeof(digit_t) * 3)); 
	ah = res2 >> (sizeof(digit_t) * 3);
	// High part b = 2^133 / P  --> lamda
    bl = 0x0000000000CBB8BA;
	bm = 0x0000000000091E2D;
    bh = 0x00000000004ED58F;
    // a & b all divided into 3 x 24 bits parts
    albl = al*bl;
	ambm = am*bm;
	ahbh = ah*bh;
	albm_ambl = (al+am)*(bm+bl)-albl-ambm;
    albh_ahbl = (al+ah)*(bh+bl)-albl-ahbh;
	ambh_ahbm = (am+ah)*(bh+bm)-ambm-ahbh;	

    res1 = albl >> k;
    res2 = albm_ambl & mask_low;
    temp = res1 + res2;
    carry = temp >> k;

    res1 = albm_ambl >> k;
    res2 = albh_ahbl & mask_low;
	res3 = ambm & mask_low;	
    temp = res1 + res2 + res3 + carry;
	carry = temp >> k;
	
	res1 = albh_ahbl >> k;
    res2 = ambm >> k;
	res3 = ambh_ahbm & mask_low;	
    temp = res1 + res2 + res3 + carry;
	carry = temp >> k;
	q[0] = temp & mask_low;
	
	res1 = ambh_ahbm >> k;
	res2 = ahbh & mask_low;	
    temp = res1 + res2 + carry;
	carry = temp >> k;
	q[0] ^= (temp & mask_low) << k;
	
    res1 = ahbh >> k;
	q[1] = res1 + carry;
	
	//r_h = t - q*R_3
	res1 = q[0] & mask_low_1;    // low 33 bits -> N2+1
	res1 = res1 * R_3;           // compute t1
	r[0] = ((t[0] & mask_low_1) - (res1 & mask_low_1)) & mask_low_1;  //compute r0
	
	res1 = r[0] - R_3;
	res2 = res1 >> 63;
	res3 = 1 - res2;
	// res3 = 1 ^ res2; //slow!

	// r[0] = (res1 & (0-res3)) + (r[0] & (0-res2));
	
	// r[0] = res1 * res3 + r[0] * res2;

	// res2 = (0 - res2);

	r[0] = res1 + (R_3 * res2);

	// if(res2==0)
	// 	r[0] = res1;
	// else
	// 	r[0] = r[0];
	
	res1 = q[0] + res3;
	q[0] = res1 & ((digit_t)(-1)>>16);
	q[1] = q[1] + (res1>>48);
	
	res1 = r[0] << 31;
	r[0] = res1 ^ s;	
}

void improved_IFFM(digit_t* c_0, digit_t* c_1, digit_t* ca, digit_t* C, const unsigned int nwords)
{// Apply the IBR function & Adjust all of the coefficients into the standard ranges
 // Inputs : the result from mp_mul_new4 or mp_mul_new5 or mp_sqr_new, raw coefficients c_0[i], c_1[i], ca[i] are all 48 bits ,for 0≤i≤n-1;
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 64 bits
	unsigned int i;
    digit_t t[2], s = 0, q[2]={0}, temp, t1, s1;
	digit_t mask_low = (digit_t)(-1) >> 33;                    // low 31 bits
	digit_t mask_low_1 = (digit_t)(-1) >> (sizeof(digit_t)*2); // low 48 bits
	digit_t R_div_3 = 0x22A359ED80000000;                      // R/3
	digit_t R = 0x67EA0DC880000000;                            // R = 2^31 * 3^20
	digit_t B[1] = {0};
	digit_t k = sizeof(digit_t) * 6;                           // 48
	digit_t h = nwords - 1;
		
	t[0] = (c_0[10] >> 31) ^ (c_1[10] << 17);
	t[1] = (ca[10] << 1) ^ (c_1[10] >> 47);
	s = c_0[10] & mask_low;
	improved_BR(t, s, q, B);	
    C[10] = B[0];		

	temp = q[0] + c_0[11];
	c_0[11] = temp & mask_low_1;
	s = temp >> k;
	temp = q[1] + c_1[11] + s;
	c_1[11] = temp & mask_low_1;
	s = temp >> k;
	ca[11] = s + ca[11];		
	t[0] = (c_0[11] >> 31) ^ (c_1[11] << 17);
	t[1] = (ca[11] << 1) ^ (c_1[11] >> 47);
	s = c_0[11] & mask_low;
	improved_BR(t, s, q, B);	
    C[11] = B[0];	            //Apply the IBR function to the n coefficients serially , start at num 10
	
	t1 = C[h] - R_div_3;
	s1 = t1>>63;
	// s = 1 - s1; //slow!
	s = 1 ^ s1;
	i = s;
	
	C[h] = t1 + (R_div_3 * s1);
	// if(s1==0)
	// 	C[h] = t1;
	// else 
	// 	C[h] = C[h];
	
	t1 = C[h] - R_div_3;
	s1 = t1>>63;
	// s = 1 - s1; //slow!
	s = 1 ^ s1;
	
	C[h] = t1 + (R_div_3 * s1);
	// if(s1==0)
	// 	C[h] = t1;
	// else 
	// 	C[h] = C[h];
	
	s += i;                     //Adjust the coefficients into the standard ranges

	temp = q[0] + (q[0]<<1) + c_0[0] + s;
	c_0[0] = temp & mask_low_1;
	s = temp >> k;
	temp = q[1] + (q[1]<<1) + c_1[0] + s;
	c_1[0] = temp & mask_low_1;
	s = temp >> k;
	ca[0] = ca[0] + s;		
	t[0] = (c_0[0] >> 31) ^ (c_1[0] << 17);
	t[1] = (ca[0] << 1) ^ (c_1[0] >> 47);
	s = c_0[0] & mask_low;
	improved_BR(t, s, q, B);	
	C[0] = B[0];
	
	for (i = 1; i < h-1; i++){
		temp = q[0] + c_0[i];
		c_0[i] = temp & mask_low_1;
		s = temp >> k;
		temp = q[1] + c_1[i] + s;
		c_1[i] = temp & mask_low_1;
		s = temp >> k;
		ca[i] = s + ca[i];		
		t[0] = (c_0[i] >> 31) ^ (c_1[i] << 17);
		t[1] = (ca[i] << 1) ^ (c_1[i] >> 47);
		s = c_0[i] & mask_low;
		improved_BR(t, s, q, B);	
        C[i] = B[0];		
	}                           //Apply the IBR function to the n coefficients serially
	
	temp = q[0] + C[10];
	c_0[10] = temp & mask_low_1;
	s = temp >> k;
	temp = q[1] + s;
	c_1[10] = temp & mask_low_1;
	s = temp >> k;
	ca[10] = s;		
	t[0] = (c_0[10] >> 31) ^ (c_1[10] << 17);
	t[1] = (ca[10] << 1) ^ (c_1[10] >> 47);
	s = c_0[10] & mask_low;
	improved_BR(t, s, q, B);	
	C[10] = B[0];
	
	C[h] += q[0];               //Apply the IBR function to the n coefficients serially
	
	t1 = C[h] - R_div_3;
	s1 = t1>>63;

	C[h] = t1 + (R_div_3 * s1);
	// if(s1==0)
	// 	C[h] = t1;
	// else 
	// 	C[h] = C[h];
	
	C[0] = C[0] + 1 - s1;       //Adjust the coefficients into the standard ranges
}

void fp_add_new(digit_t* a, digit_t* b, digit_t* c, const unsigned int nwords)
{// Modular Addition based on the new data representation , with the unconditional radix for p751
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], a[i] & b[i] are 64 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 64 bits , to compute c = a + b
	unsigned int i,j;
	digit_t ca = 0, ca1, t;
	digit_t R  = 0x67EA0DC880000000;       // R = 2^31 * 3^20
	digit_t R_div_3 = 0x22A359ED80000000;  // R/3
 
    for(i = 0; i < nwords; i++){
		c[i] = a[i] + b[i];
	}

	for (i = 1; i < nwords; i++){
		t = c[i-1] - R;
		ca1 = t>>63;
		// ca = 1 - ca1; //slow!
		ca = 1 ^ ca1;
		c[i-1] = t  + R * ca1;
		c[i] = c[i] + ca;
	} 

	t = c[nwords-1] - R_div_3;
	ca1 = t >> 63;
	//ca = 1 - ca1;
	ca = 1 ^ ca1;
	c[nwords-1] = t + R_div_3 * ca1;
	c[0] += ca;

	// for (i = 1; i < nwords; i++){
	// 	t = c[i-1] - R;
	// 	ca = t>>63;
	// 	c[i-1] = t  + R * ca;
	// 	ca = 1 ^ ca;
	// 	c[i] = c[i] + ca;
	// } 
	// t = c[nwords-1] - R_div_3;
 	// ca = t >> 63;
 	// c[nwords-1] = t + R_div_3 * ca;
 	// ca = 1 ^ ca;
 	// c[0] += ca;
}

void fpneg_new(digit_t* a,digit_t* b,const unsigned int nwords)
{// Modular Negation based on the new data representation , with the unconditional radix for p751
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], a[i] all in 64 bits
 // Output : b[i] in [0, R-1] for 0≤i<n-1; b[n-1] in [0, R/3 - 1], b[i] all in 64 bits, to compute b = -a
	unsigned int i;
	digit_t R = 0x67EA0DC880000000;       // R = 2^31 * 3^20
	digit_t R_div_3 = 0x22A359ED80000000; // R/3
	digit_t k = nwords - 1;
	
	for( i = 0; i < k; i++){
		b[i] = R - a[i] - 1;
	}

	b[k] =  R_div_3 - a[k] - 1;
}

void fpsub_new(digit_t* a, digit_t* b, digit_t* c, const unsigned int nwords)
{// Modular Subtraction based on the new data representation , with the unconditional radix for p751
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], a[i] & b[i] are 64 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 64 bits , to compute c = a - b
	unsigned int i, j = 0, t;
	digit_t R = 0x67EA0DC880000000;       // R = 2^31 * 3^20
	digit_t R_div_3 = 0x22A359ED80000000;
	
	for(i=0;i<nwords;i++){
		c[i] = a[i] - b[i];
	}
	
	// t = 1 ^ ((c[0] | (0 - c[0])) >> 63); //slow
	// t = (c[0] == 0);
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

void fpdiv2_751_new(digit_t* a, digit_t* c)
{// Modular Division by Two based on the new data representation , with the unconditional radix
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], a[i] all in 64 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], c[i] all in 64 bits, to compute c = a/2
    unsigned int i;
    digit_t mask;
	digit_t R_div_2 = 0x33F506E440000000;  //  R/2
	digit_t R_div_6 = 0x1151ACF6C0000000;  //  R/3/2
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

void fpmul_new(const digit_t* a, digit_t* b, digit_t* c)
{// Modular Multipication based on the new data representation , with the unconditional radix for p751
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], all in 64 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 64 bits , to compute c = a * b
    digit_t c1_0[NWORDS_FIELD];
	digit_t c1_1[NWORDS_FIELD];
	digit_t ca[NWORDS_FIELD];

	mp_mul_new4(a, b, c1_0, c1_1, ca, NWORDS_FIELD);
	improved_IFFM(c1_0, c1_1, ca, c, NWORDS_FIELD);
}

void mp_sqr_new(const digit_t* a, digit_t* c_0, digit_t* c_1, digit_t* ca, const unsigned int nwords)
{// Compute multiply-accumulate terms and get the raw coefficients of a^2
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], all in 64 bits
 // Output : raw coefficients c_0[i], c_1[i], ca[i] are all 48 bits ,for 0≤i≤n-1; 
    unsigned int i, j;
    digit_t UV[3], UV_2[3] = {0},temp;
    digit_t carryOut, k, s;
	digit_t h = sizeof(digit_t) * 6;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*2); //low 48 bits
	
	digit_x_digit_2(a[0], a[0], UV);
    c_0[0] = UV[0];
    c_1[0] = UV[1];
    ca[0]  = UV[2];
	
	for (i = 2; i < nwords; i = i + 2) {
		k = i>>1;
		for (j = 0; j < k; j++) {
             digit_x_digit_2(a[j], a[i-j], UV);	
             digit_a_digit_1(UV, UV_2, UV_2);
        }
		digit_a_digit_1(UV_2, UV_2, UV);
		digit_x_digit_2(a[k], a[k], UV_2);	
		digit_a_digit_1(UV, UV_2, UV_2);
		c_0[i] = UV_2[0];
		c_1[i] = UV_2[1];
		ca[i]  = UV_2[2];
		UV_2[0] = 0;
		UV_2[1] = 0;
		UV_2[2] = 0;
	}	//e.g. c4 = 2 * a0 * a4 + 2 * a1 * a3 + a2 ^ 2; c is divided into 3 x 48 bits parts
	
	for (i = 1; i < nwords; i = i + 2) {
		 for (j = 0; j <= (i-1)/2; j++) {
             digit_x_digit_2(a[j], a[i-j], UV);	
             digit_a_digit_1(UV, UV_2, UV_2);
        }
		digit_a_digit_1(UV_2, UV_2, UV);
		c_0[i] = UV[0];
		c_1[i] = UV[1];
		ca[i]  = UV[2];
		UV_2[0] = 0;
		UV_2[1] = 0;
		UV_2[2] = 0;
	}	//e.g. c5 = 2 * a0 * a5 + 2 * a1 * a4 + 2 * a2 * a3 ;c is divided into 3 x 48 bits parts
	
	for (i = nwords; i < 2*nwords-2; i = i + 2) {
		k = i>>1;
		for (j = i-nwords+1; j < k; j++) {
			digit_x_digit_2(a[j], a[i-j], UV);			
            digit_a_digit_1(UV, UV_2, UV_2);
        }
		digit_a_digit_1(UV_2, UV_2, UV);
		digit_x_digit_2(a[k], a[k], UV_2);	
		digit_a_digit_1(UV, UV_2, UV_2);                 //e.g. c18= 2* a11 * a7 + 2 * a10 * a8 + a9 ^ 2;c is divided into 3 x 48 bits parts
		s = i-nwords;
		temp = (UV_2[0]<<1) + UV_2[0] + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> h;
		temp = (UV_2[1]<<1) + UV_2[1] + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> h;
		ca[s] = (UV_2[2]<<1) + UV_2[2] + ca[s] + carryOut;
		UV_2[0] = 0;
		UV_2[1] = 0;
		UV_2[2] = 0;
	}                                                    // get the raw coefficients.
	
	for (i = nwords + 1; i < 2*nwords-2; i = i + 2) {
		for (j = i-nwords+1; j <= (i-1)/2; j++) {
			digit_x_digit_2(a[j], a[i-j], UV);			
            digit_a_digit_1(UV, UV_2, UV_2);
        }
		digit_a_digit_1(UV_2, UV_2, UV);                 //e.g. c19 = 2 * a11 * a8 + 2 * a10 * a9;
		s = i-nwords;
		temp = (UV[0]<<1) + UV[0] + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> h;
		temp = (UV[1]<<1) + UV[1] + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> h;
		ca[s] = (UV[2]<<1) + UV[2] + ca[s] + carryOut;
		UV_2[0] = 0;
		UV_2[1] = 0;
		UV_2[2] = 0;
	}                                                    // get the raw coefficients.

	digit_x_digit_2(a[11], a[11], UV_2);
	temp = (UV_2[0]<<1) + UV_2[0] + c_0[10];
	c_0[10] = temp & mask_low;
	carryOut = temp >> h;
	temp = (UV_2[1]<<1) + UV_2[1] + c_1[10] + carryOut;
	c_1[10] = temp & mask_low;
	carryOut = temp >> h;
	ca[10] = (UV_2[2]<<1) + UV_2[2] + ca[10] + carryOut; // get the raw coefficients.
}

void fpsqr_new(const digit_t* a, digit_t* c)
{// Modular Squaring based on the new data representation , with the unconditional radix
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], all in 64 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 64 bits , to compute c = a^2
	digit_t c1_0[NWORDS_FIELD];
	digit_t c1_1[NWORDS_FIELD];
	digit_t ca[NWORDS_FIELD];

	mp_sqr_new(a, c1_0, c1_1, ca, NWORDS_FIELD);         //perform better than mp_mul_new4
	improved_IFFM(c1_0, c1_1, ca, c, NWORDS_FIELD);
}

void fpinv_chain_new(digit_t* a , const unsigned int nwords)
{// Modular Inversion
	unsigned int i, j;
	digit_t tt[nwords];
	digit_t t0[nwords], t1[nwords], t2[nwords], t3[nwords], t4[nwords], t5[nwords], t6[nwords], t7[nwords], t8[nwords];
	digit_t t9[nwords], t10[nwords],t11[nwords],t12[nwords],t13[nwords],t14[nwords],t15[nwords],t16[nwords],t17[nwords];
	digit_t t18[nwords],t19[nwords],t20[nwords],t21[nwords],t22[nwords],t23[nwords],t24[nwords],t25[nwords],t26[nwords];

//  Precomputed table
    fpsqr_new(a, tt);
    fpmul_new(a, tt, t0);
    fpmul_new(t0, tt, t1);
    fpmul_new(t1, tt, t2);
    fpmul_new(t2, tt, t3); 
    fpmul_new(t3, tt, t3);
    fpmul_new(t3, tt, t4);
	fpmul_new(t4, tt, t5);
	fpmul_new(t5, tt, t6);
	fpmul_new(t6, tt, t7);
	fpmul_new(t7, tt, t8);
	fpmul_new(t8, tt, t9);
    fpmul_new(t9, tt, t9);
    fpmul_new(t9 , tt, t10);
	fpmul_new(t10, tt, t11);
	fpmul_new(t11, tt, t12);
	fpmul_new(t12, tt, t13);
	fpmul_new(t13, tt, t14);
	fpmul_new(t14, tt, t15);
	fpmul_new(t15, tt, t16);
	fpmul_new(t16, tt, t17);
	fpmul_new(t17, tt, t18);
	fpmul_new(t18, tt, t19);
	fpmul_new(t19, tt, t20);
	fpmul_new(t20, tt, t21);
    fpmul_new(t21, tt, t21); 
    fpmul_new(t21, tt, t22);
	fpmul_new(t22, tt, t23);
	fpmul_new(t23, tt, t24);
	fpmul_new(t24, tt, t25);
    fpmul_new(t25, tt, t25);
    fpmul_new(t25, tt, t26);

    for (i = 0; i < nwords; i++)
		tt[i] = a[i];
 
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t20, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t24, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t11, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t8 , tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t2 , tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t23, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t2 , tt, tt);
    for (i = 0; i < 9; i++) fpsqr_new(tt, tt);
    fpmul_new(t2 , tt, tt);
    for (i = 0; i < 10; i++) fpsqr_new(tt, tt);
    fpmul_new(t15, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t13, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t26, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t20, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t11, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t10, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t14, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t4 , tt, tt);
    for (i = 0; i < 10; i++) fpsqr_new(tt, tt);
    fpmul_new(t18, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t1 , tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t22, tt, tt);
    for (i = 0; i < 10; i++) fpsqr_new(tt, tt);
    fpmul_new(t6 , tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t24, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t9 , tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t18, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t17, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(a  , tt, tt);
    for (i = 0; i < 10; i++) fpsqr_new(tt, tt);
    fpmul_new(t16, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t7 , tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t0 , tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t12, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t19, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t22, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t25, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t2 , tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t10, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t22, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t18, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t4 , tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t14, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t13, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t5 , tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t23, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t21, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t2 , tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t23, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t12, tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t9 , tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t3 , tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t13, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t17, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t26, tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t5 , tt, tt);
    for (i = 0; i < 8; i++) fpsqr_new(tt, tt);
    fpmul_new(t8 , tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t2 , tt, tt);
    for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
    fpmul_new(t11, tt, tt);
    for (i = 0; i < 7; i++) fpsqr_new(tt, tt);
    fpmul_new(t20, tt, tt);

    for (j = 0; j < 61; j++) {
        for (i = 0; i < 6; i++) fpsqr_new(tt, tt);
        fpmul_new(t26, tt, tt);
    }

    for (i = 0; i < nwords; i++)
		a[i] = tt[i];
}

void fpinv_new(digit_t* a, digit_t* c, const unsigned int nwords)
{// Modular Inversion for p751
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], all in 64 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 64 bits
    digit_t tt[nwords];
	unsigned int i;

    for (i = 0; i < nwords; i++)
		tt[i] = a[i];
	
    fpinv_chain_new(tt, nwords);
    fpsqr_new(tt, tt);
    fpsqr_new(tt, tt);
    fpmul_new(a, tt, c);
}

void fp2mul_uRadix(f2elm_t a, f2elm_t b, f2elm_t c)
{// GF(p^2) multiplication, c = a*b in GF(p^2).
 // Inputs : a = a0 + a1*i, b = b0 + b1*i. a0, a1, b0, b1 all meet the same requirement with fpmul_new & fpsqr_new
 // Output : c = c0 + c1*i. c = a * b, c0 and c1 meet the same requirement with fpmul_new & fpsqr_new
 // a-b+2^5*R*P mod P = (a_11-b_11+2^5*(R/3-1)*R)*R^11 + (a_10-b_10+2^5*(R-1)*R)*R^10 + ... + a_0-b_0+2^5*(R-1)*R
 
    digit_t PR_11[3] = {0x000046F000000000, 0x0000C2FFC04B02BE, 0x00000001C1EC8B85}; //2^5*(R/3-1)*R)
    digit_t PR[3]    = {0x000046F000000000, 0x000048FF40FB02BE, 0x0000000545C5A291}; //2^5*(R-1)*R
    felm_t t1, t2;
    felm_t tt1[3], tt2[3], tt3[3]; 
	int i, j;
	digit_t temp, res;
	digit_t mask_low_1 = (digit_t)(-1) >> (sizeof(digit_t)*2);
	digit_t k = sizeof(digit_t) * 6;
	f2elm_t d;

	for (i=0; i<NWORDS_FIELD; i++){
		t1[i] = a[0][i] + a[1][i];
		t2[i] = b[0][i] + b[1][i];
	}
	
	//for c[0]
    mp_mul_new4(a[0], b[0], tt1[0], tt1[1], tt1[2], NWORDS_FIELD); // tt1 = a0*b0	
	mp_mul_new4(a[1], b[1], tt2[0], tt2[1], tt2[2], NWORDS_FIELD); // tt2 = a1*b1

	for (i=0; i<NWORDS_FIELD-1; i++) {
		//a0*b0 + 2^5*R*P - a1*b1
		res = 0;
		for (j=0; j<2; j++){
			temp = tt1[j][i] + PR[j] - tt2[j][i] + (res&(digit_t)1) - (res&(digit_t)2);
			tt3[j][i] = temp & mask_low_1;
			res = temp >> k;
		}
		tt3[2][i] = tt1[2][i] + PR[2] - tt2[2][i] + (res&(digit_t)1) - (res&(digit_t)2);
	}
	
	//a0*b0 + 2^5*R*P - a1*b1
	res = 0;
	for (j=0; j<2; j++){
		temp = tt1[j][11] + PR_11[j] - tt2[j][11] + (res&(digit_t)1) - (res&(digit_t)2);
		tt3[j][11] = temp & mask_low_1;
		res = temp >> k;
	}

	tt3[2][11] = tt1[2][11] + PR_11[2] - tt2[2][11] + (res&(digit_t)1) - (res&(digit_t)2);

	improved_IFFM(tt3[0], tt3[1], tt3[2], c[0], NWORDS_FIELD);// c[0]
	
	//for c[1]
	
	mp_mul_new5(t1, t2, tt3[0], tt3[1], tt3[2], NWORDS_FIELD); // tt1 = (a0+a1)*(b0+b1)	
	
	// c[1] = (a0+a1)*(b0+b1) - (a0*b0 + a1*b1) 
	for (i=0; i<NWORDS_FIELD; i++) {
		//a0*b0 + a1*b1
		res = 0;
		for (j=0; j<2; j++){
			temp = tt2[j][i] + tt1[j][i] + res;
			tt2[j][i] = temp & mask_low_1;
			res = temp >> k;
		}

		tt2[2][i] = tt2[2][i] + tt1[2][i] + res;
		//(a0+a1)*(b0+b1) - a0*b0 - a1*b1
		res = 0;
		for (j=0; j<2; j++){
			temp = tt3[j][i] - tt2[j][i] - res;
			tt3[j][i] = temp & mask_low_1;
			res = temp >> (sizeof(digit_t) * 8 - 1);
		}

		tt3[2][i] = tt3[2][i] - tt2[2][i] - res;
	}
	
    improved_IFFM(tt3[0], tt3[1], tt3[2], c[1], NWORDS_FIELD); 
}

void from_uRadix(digit_t* a, digit_t* c, const unsigned int nwords)
{// From Unconventional Radix Back to Normal
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], all in 64 bits based on the new data representation , with the unconditional radix
 // Output : c in F(p) , presented in normal presentation 
	digit_t R = 0x67EA0DC880000000;
	digit_t R_sub_1 = 0x67EA0DC87FFFFFFF;     // R - 1
	digit_t R_div_3_sub1 = 0x22A359ED7FFFFFFF;// R/3 - 1
	unsigned int i,j,k;
	digit_t c1[12][2] = {0};
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
	//     }
		
	// 	if(k == 1)
	// 		a[0] = a[0];
	// 	else
	// 		a[0] = a[0];
	// }   //Adjust the coefficients into the standard ranges

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
}

void IBR_uRadix(const digit_t* a, digit_t* q, digit_t* r, const unsigned int unwords, const digit_t* lamda)
{// IBR 
	digit_t t[unwords], s, t_r, temp[2*unwords], res1, res2, res3, mask_low_1 = (digit_t)(-1) >> (sizeof(digit_t)*4-1);
	digit_t R_3 = 0x00000000CFD41B91;            // R_3=3^20=CFD41B91
	unsigned int i, j;
	
	s = a[0] & ((digit_t)(-1) >> 33);
	t_r = a[0] >> 31;
	for (i = 0; i < unwords; i++){
		t[i] = (a[i] >> 61) ^ (a[i+1] << 3);
	}

	mp_mul_old(t, lamda, temp, unwords);
	
	for (i = 0; i < unwords; i++){
		res1 = temp[unwords+i-1] >> (63-unwords);
		res2 = temp[unwords+i] << (unwords+1);
		q[i] = res1 ^ res2;
	}

	res1 = q[0] & mask_low_1;
	res1 = res1 * R_3;
	r[0] = (t_r - (res1 & mask_low_1)) & mask_low_1;

	res1 = r[0] - R_3;
	res2 = res1 >> 63;
	res3 = 1 - res2;
	// res3 = 1 ^ res2;
	r[0] = res1 + R_3 * res2;
	
	res1 = q[0] + res3;
	// if (q[0]==(digit_t)(-1)){
	// 		q[0] = 0;
	// 		j = 1;
	// }
	// else{
	// 		q[0] = res1;
	// 		j = 0;
	// }

	j = (q[0]==(digit_t)(-1));
	q[0] = res1 * (1 ^ j);

	for (i = 1; i < unwords-1; i++){	
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
	
	res1 = r[0] << 31;
	r[0] = res1 ^ s;
}

void to_uRadix(digit_t* a, digit_t* uc, const unsigned int nwords)
{// Conversion to uncoventional radix representation,
 // Inputs : a in F(p) normal field
 // Output : uc[i] in [0, R-1] for 0≤i<n-1; uc[n-1] in [0, R/3 - 1], uc[i] all in 64 bits based on the new data representation , with the unconditional radix
 // mc = a_11*R^11+a_10*R^10+...+a_1*R+a_0.
 // a_i= a mod R, a = (a/R). 

	int i, j;
	digit_t q[nwords], temp[2];
 // 752 minus 63  every time to compute gamma
	digit_t lamda[11][11]={{0x61fe8a242ce4800f, 0xdf46b0c3072f7125, 0xc7845ed41934002d, 0x77b30fdc5c0958b8, 0x2ca7484f6053f900, 0xf70020f49cea77a1, 0xe4f7c3e4bcb2e546, 0x9289ffa9174dfc84, 0xe29b343a29542b47, 0xbb8ba3ef0ab5e0e7,0x0004ed58f091e2dc}, 
						   {0xbe8d61860e5ee24a, 0x8f08bda83268005b, 0xef661fb8b812b171, 0x594e909ec0a7f200, 0xee0041e939d4ef42, 0xc9ef87c97965ca8d, 0x2513ff522e9bf909, 0xc536687452a8568f, 0x771747de156bc1cf, 0x0009dab1e123c5b9}, 
						   {0x1e117b5064d000b7, 0xdecc3f71702562e3, 0xb29d213d814fe401, 0xdc0083d273a9de84, 0x93df0f92f2cb951b, 0x4a27fea45d37f213, 0x8a6cd0e8a550ad1e, 0xee2e8fbc2ad7839f, 0x0013b563c2478b72}, 
						   {0xbd987ee2e04ac5c6, 0x653a427b029fc803, 0xb80107a4e753bd09, 0x27be1f25e5972a37, 0x944ffd48ba6fe427, 0x14d9a1d14aa15a3c, 0xdc5d1f7855af073f, 0x00276ac7848f16e5}, 
						   {0xca7484f6053f9007, 0x70020f49cea77a12, 0x4f7c3e4bcb2e546f, 0x289ffa9174dfc84e, 0x29b343a29542b479, 0xb8ba3ef0ab5e0e7e, 0x004ed58f091e2dcb}, 
						   {0xe0041e939d4ef425, 0x9ef87c97965ca8de, 0x513ff522e9bf909c, 0x536687452a8568f2, 0x71747de156bc1cfc, 0x009dab1e123c5b97}, 
						   {0x3df0f92f2cb951bd, 0xa27fea45d37f2139, 0xa6cd0e8a550ad1e4, 0xe2e8fbc2ad7839f8, 0x013b563c2478b72e}, 
						   {0x44ffd48ba6fe4272, 0x4d9a1d14aa15a3c9, 0xc5d1f7855af073f1, 0x0276ac7848f16e5d}, 
						   {0x9b343a29542b4792, 0x8ba3ef0ab5e0e7e2, 0x04ed58f091e2dcbb}, 
						   {0x1747de156bc1cfc5, 0x09dab1e123c5b977}, 
						   {0x13b563c2478b72ee}};

    for (i=0; i<nwords-1; i++){
		IBR_uRadix(a, q, temp, nwords-i-1, lamda[i]);
		uc[i] = temp[0];
		for (j=0; j<nwords-i-1; j++){
			a[j] = q[j];
		}
	}

	uc[nwords-1] = q[0];
}

void fp2sqr_uRadix(f2elm_t a, f2elm_t c)
{// GF(p^2) squaring using, c = a^2 in GF(p^2).
 // Inputs : a = a0 + a1*i, a0 and a1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new
 // Output : c = c0 + c1*i, c = a^2; co & c1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new
    felm_t t1, t2, t3;
    
    fp_add_new(a[0], a[1], t1, NWORDS_FIELD);        // t1 = a0+a1 
    fpsub_new (a[0], a[1], t2, NWORDS_FIELD);        // t2 = a0-a1
    fp_add_new(a[0], a[0], t3, NWORDS_FIELD);        // t3 = 2a0
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

