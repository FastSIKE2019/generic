/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: portable modular arithmetic for P610
*********************************************************************************************/

#include "../P610_internal.h"


void digit_x_digit_1( digit_t a,  digit_t b, digit_t* co)
{// Digit multiplication, digit * digit -> 2-digit result   
 // Inputs : 60bits a & b
 // Output : 120bits c = a * b , c is divided into 2 x 60bits parts
    register digit_t al, ah, bl, bh, temp;  
    digit_t albl, albh_ahbl, ahbh, res1, carry;
    digit_t mask_low = (digit_t)(-1) >> 34;
	digit_t k = 30; 

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

void improved_BR( digit_t c_3,  digit_t c_2,  digit_t c_1,  digit_t c_0, digit_t* q, digit_t* r)
{// Improved Barrett reduction (IBR) 
 // t = floor(c / 2^N_1) , s = c mod 2^N_1 for 0≤c<2^(2N+lamda) (2*102+4=208)
 // Output : q = floor(c / p) , r = c mod p , for p = 2^N_1 * R_3
 // N_1 = 51, N_2=log(R_3) = 51, R_3=3^32=0x6954FE21E3E81   //log(t) -> 208bit ->  157 -> 108
 // 2*102+4+1=209  2^209 / p = 9B8BD843858 6D9D3E741CCA5D76   
 // p = 34AA7F10F1F408000000000000    
 // 102+3+4=109 
 
	digit_t albl[2], ahbh[2], ahbl[2], albh[2];
	digit_t res1[2], res2, res3;
	digit_t carry, tmp, t;
	digit_t mask_low = (digit_t)(-1) >> 12;  //low 52
	digit_t mask_low_1 = (digit_t)(-1) >> 13;//low 51
	digit_t mask_low_2 = (digit_t)(-1) >> 4; //low 60
	
//	digit_t b[2] = {0x0D9D3E741CCA5D76, 0x00009B8BD8438586};

	digit_t b[4] = {0xCCA5D76, 0x9D3E741, 0x438586D, 0x9B8BD8}; //9B8BD8 438586D 9D3E741 CCA5D76
	digit_t m = 0x0006954FE21E3E81;
	
	digit_t t_r[4];
	digit_t mask_low_3 = (digit_t)(-1) >> 36; //low 28
	t_r[0] = ((c_1>>40)^(c_2<<20)) & mask_low_3;
	t_r[1] = (c_2 >> 8) & mask_low_3;
	t_r[2] = ((c_2 >> 36)^(c_3 << 24)) & mask_low_3;
	t_r[3] = c_3 >> 4;
	
	digit_t a0b0, a1b1, a2b2, a3b3, a0b1_a1b0, a0b2_a2b0, a0b3_a3b0, a1b2_a2b1, a1b3_a3b1, a2b3_a3b2;
	a0b0 = t_r[0] * b[0];
	a1b1 = t_r[1] * b[1];
	a2b2 = t_r[2] * b[2];
	a3b3 = t_r[3] * b[3];
	a0b1_a1b0 = (t_r[0] + t_r[1])*(b[0] + b[1]) - a0b0 - a1b1;
	a0b2_a2b0 = (t_r[0] + t_r[2])*(b[0] + b[2]) - a0b0 - a2b2;
	a0b3_a3b0 = (t_r[0] + t_r[3])*(b[0] + b[3]) - a0b0 - a3b3;
	a1b2_a2b1 = (t_r[1] + t_r[2])*(b[1] + b[2]) - a1b1 - a2b2;
	a1b3_a3b1 = (t_r[1] + t_r[3])*(b[1] + b[3]) - a1b1 - a3b3;
	a2b3_a3b2 = (t_r[2] + t_r[3])*(b[2] + b[3]) - a2b2 - a3b3;
	
	res2 = (a0b0 >> 28) + a0b1_a1b0;
	res3 = a1b1 + a0b2_a2b0 + (res2 >> 28); 
	
	res2 = a1b2_a2b1 + a0b3_a3b0 + (res3 >> 28);
	q[0] = (res2 & mask_low_3) >> 25;
	
	res3 = a2b2 + a1b3_a3b1 + (res2 >> 28);
	q[0] ^= ((res3 & mask_low_3) << 3);
	res2 = a2b3_a3b2 + (res3 >> 28);
	q[0] ^= ((res2 & mask_low_3) << 31);
	
	res3 = a3b3 + (res2 >> 28);
	q[0] ^= ((res3 << 59) & mask_low_2);
	q[1] = res3 >> 1;
	
	tmp = q[0] & mask_low;
	digit_x_digit_1(tmp, m, res1);
	
	t = (c_0 >> 51)^(c_1 << 9);
	r[0] = ((t & mask_low) - (res1[0] & mask_low)) & mask_low;
	
	t = r[0] - m;
	tmp = t >> 63;
	carry = 1 - tmp;
	r[0] = t + (m * tmp);
	
	tmp = q[0] + carry;
	q[0] = tmp & mask_low_2;
	carry = tmp >> 60;
	q[1] = q[1] + carry;
	
	t = c_0 & mask_low_1;
	
	r[1] = r[0] >> 9;
	r[0] = ((r[0] << 51) ^ t) & mask_low_2;
}

void last_improved_BR(digit_t c_1,  digit_t c_0, digit_t* q, digit_t* r)
{
	digit_t b[2] = {0x0D9D3E741CCA5D76, 0x00009B8BD8438586};
	digit_t m = 0x0006954FE21E3E81;
	digit_t t_r[1]; digit_t abl[2], abh[2], res1, res2, t, tmp, carry;
	digit_t mask_low = (digit_t)(-1) >> 12;  //low 52
	digit_t mask_low_1 = (digit_t)(-1) >> 13;//low 51
	digit_t mask_low_2 = (digit_t)(-1) >> 4; //low 60
	
	t_r[0] = c_1 >> 40;
	
	res1 = t_r[0] * b[1];
	
	q[0] = res1 >> 49; q[1] = 0;
	
	res1 = q[0] * m;
	t = (c_0 >> 51)^(c_1 << 9);
	r[0] = ((t & mask_low) - (res1 & mask_low)) & mask_low;
	
	t = r[0] - m;
	tmp = t >> 63;
	carry = 1 - tmp;
	r[0] = t + (m * tmp);
	
	q[0] = q[0] + carry;
	
	t = c_0 & mask_low_1;
	
	r[1] = r[0] >> 9;
	r[0] = ((r[0] << 51) ^ t) & mask_low_2;
}

void improved_IFFM(digit_t* c_0, digit_t* c_1, digit_t* c_2, digit_t* c_3, digit_t* C0, digit_t* C1)
{
	digit_t q[2], r[2];
	digit_t tmp1, tmp2, tmp3, tmp, c_r, carry;
	digit_t h, i, j, k;
	unsigned int m;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	digit_t R_div_2[2] = {0xA04000000000000, 0x1A553F8878F};//1A553F8878F A04000000000000
	
	improved_BR(c_3[4], c_2[4], c_1[4], c_0[4], q, r);	
	C0[4] = r[0];
	C1[4] = r[1];
	
	tmp = c_0[5] + q[0];
	tmp1 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_1[5] + q[1] + carry;
	tmp2 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_2[5] + carry;
	tmp3 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_3[5] + carry;
	
	improved_BR(tmp, tmp3, tmp2, tmp1, q, r);
	C0[5] = r[0];
	C1[5] = r[1];
	
	tmp = C0[5] - R_div_2[0];
	j = tmp >> 63;
	c_r = C1[5] - R_div_2[1] - j;
	k = c_r >> 63;
	i = 1 ^ k;
	h = i;
	
	C0[5] = (tmp + R_div_2[0] * k) & mask_low;
	C1[5] = c_r + R_div_2[1] * k + j*k;
	
	tmp = C0[5] - R_div_2[0];
	j = tmp >> 63;
	c_r = C1[5] - R_div_2[1] - j;
	k = c_r >> 63;
	i = 1 ^ k;
	
	C0[5] = (tmp + R_div_2[0] * k) & mask_low;
	C1[5] = c_r + R_div_2[1] * k + j*k;
	
	h += i;//carry
	
	tmp = c_0[0] + (q[0] << 1) + h;
	tmp1 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_1[0] + (q[1] << 1) + carry;
	tmp2 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_2[0] + carry;
	tmp3 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_3[0] + carry;
	
	improved_BR(tmp, tmp3, tmp2, tmp1, q, r);
	C0[0] = r[0];
	C1[0] = r[1];
	
	for (m = 1; m < 4; m++){
		tmp = c_0[m] + q[0];
		tmp1 = tmp & mask_low;
		carry = tmp >> 60;
		tmp = c_1[m] + q[1] + carry;
		tmp2 = tmp & mask_low;
		carry = tmp >> 60;
		tmp = c_2[m] + carry;
		tmp3 = tmp & mask_low;
		carry = tmp >> 60;
		tmp = c_3[m] + carry;
	
		improved_BR(tmp, tmp3, tmp2, tmp1, q, r);
		C0[m] = r[0];
		C1[m] = r[1];
	}
	tmp = C0[4] + q[0];
	tmp1 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = C1[4] + q[1] + carry;
	tmp2 = tmp & mask_low;

	last_improved_BR(tmp2, tmp1, q, r);

	C0[4] = r[0];
	C1[4] = r[1];
	
	tmp = C0[5] + q[0];
	C0[5] = tmp & mask_low;
	carry = tmp >> 60;
	C1[5] = C1[5] + q[1] + carry;
	
	tmp = C0[5] - R_div_2[0];
	j = tmp >> 63;
	c_r = C1[5] - R_div_2[1] - j;
	k = c_r >> 63;
	i = 1 ^ k;
	
	C0[5] = (tmp + R_div_2[0] * k) & mask_low;
	C1[5] = c_r + R_div_2[1] * k + j*k;
	
	tmp = C0[0] + i;
	C0[0] = tmp & mask_low;
	carry = tmp >> 60;
	C1[0] = C1[0] + carry;
}

void digit_x_digit_2( digit_t a1,  digit_t a0,  digit_t b1,  digit_t b0, digit_t* c)
{
	digit_t a0b0, a1b1, a2b2, a3b3, a0b1_a1b0, a0b2_a2b0, a0b3_a3b0, a1b2_a2b1, a1b3_a3b1, a2b3_a3b2, t_r[4], b[4], res2, res3;
	digit_t mask_low = (digit_t)(-1) >> 34; //low 30
	digit_t mask_low_1 = (digit_t)(-1) >> 4; //low 60
	
	t_r[0] = a0 & mask_low; t_r[1] = a0 >> 30; t_r[2] = a1 & mask_low; t_r[3] = a1 >> 30;
	b[0] = b0 & mask_low; b[1] = b0 >> 30; b[2] = b1 & mask_low; b[3] = b1 >> 30;
	
	a0b0 = t_r[0] * b[0];
	a1b1 = t_r[1] * b[1];
	a2b2 = t_r[2] * b[2];
	a3b3 = t_r[3] * b[3];
	a0b1_a1b0 = (t_r[0] + t_r[1])*(b[0] + b[1]) - a0b0 - a1b1;
	a0b2_a2b0 = (t_r[0] + t_r[2])*(b[0] + b[2]) - a0b0 - a2b2;
	a0b3_a3b0 = (t_r[0] + t_r[3])*(b[0] + b[3]) - a0b0 - a3b3;
	a1b2_a2b1 = (t_r[1] + t_r[2])*(b[1] + b[2]) - a1b1 - a2b2;
	a1b3_a3b1 = (t_r[1] + t_r[3])*(b[1] + b[3]) - a1b1 - a3b3;
	a2b3_a3b2 = (t_r[2] + t_r[3])*(b[2] + b[3]) - a2b2 - a3b3;
	
	c[0] = a0b0 & mask_low;
	res2 = (a0b0 >> 30) + a0b1_a1b0;
	c[0] ^= ((res2 & mask_low) << 30);
	res3 = a1b1 + a0b2_a2b0 + (res2 >> 30); 
	c[1] = (res3 & mask_low);
	res2 = a1b2_a2b1 + a0b3_a3b0 + (res3 >> 30);
	c[1] ^= ((res2 & mask_low) << 30);
	res3 = a2b2 + a1b3_a3b1 + (res2 >> 30);
	c[2] = (res3 & mask_low); 
	res2 = a2b3_a3b2 + (res3 >> 30);
	c[2] ^= ((res2 & mask_low) << 30);
	c[3] = a3b3 + (res2 >> 30);
}

void digit_x_digit_2b( digit_t a1,  digit_t a0,  digit_t b1,  digit_t b0, digit_t* c)//60bit
{// 128x128 max 103X103 -> 4x60bit 
 // input 60bits x 2
 // output 60bits x 4
	digit_t albl[2], ahbh[2], ahbl[2], albh[2];
	digit_t tmp, tmp2[2], tmp3[2], carry ;
	digit_t mask_low = (digit_t)(-1) >> 4; // 4x60bit
	
	digit_x_digit_1(a0, b0, albl);
	digit_x_digit_1(a0, b1, albh);
	digit_x_digit_1(a1, b0, ahbl);
	digit_x_digit_1(a1, b1, ahbh);
	
	c[0] = albl[0];
	
	tmp = albl[1] + albh[0] + ahbl[0];
	c[1] = tmp & mask_low;
	carry = tmp >> 60;
	
	tmp = albh[1] + ahbl[1] + ahbh[0] + carry;
	c[2] = tmp & mask_low;
	carry = tmp >> 60;
	
	c[3] = ahbh[1] + carry;
}

void digit_a_digit_2( digit_t a1,  digit_t a0,  digit_t b1,  digit_t b0, digit_t* c)//60bit
{// 102+102  ->  2x60bits
 // input 60bits x 2
 // output 60bits x 2
	digit_t tmp, carry;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	tmp = a0 + b0;
	c[0] = tmp & mask_low;
	carry = tmp >> 60;
	c[1] = a1 + b1 + carry;
}

void digit_s_digit_4( digit_t* a,  digit_t* b, digit_t* c)
{// input 60x4
 // output 60x4
	digit_t tmp1, tmp2, tmp3, tmp4, borrow;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	tmp1 = a[0] - b[0];
	borrow = tmp1 >> 63;
	tmp2 = a[1] - b[1] - borrow;
	borrow = tmp2 >> 63;
	tmp3 = a[2] - b[2] - borrow;
	borrow = tmp3 >> 63;
	tmp4 = a[3] - b[3] - borrow;
	c[0] = tmp1 & mask_low;
	c[1] = tmp2 & mask_low;
	c[2] = tmp3 & mask_low;
	c[3] = tmp4;// & mask_low;
}

void digit_a_digit_4( digit_t* a, digit_t* b, digit_t* c)
{// input 60x4
 // output 60x4
	digit_t tmp1, tmp2, tmp3, tmp4, carry;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	tmp1 = a[0] + b[0];
	c[0] = tmp1 & mask_low;
	carry = tmp1 >> 60;
	tmp2 = a[1] + b[1] + carry;
	c[1] = tmp2 & mask_low;
	carry = tmp2 >> 60;
	tmp3 = a[2] + b[2] + carry;
	c[2] = tmp3 & mask_low;
	carry = tmp3 >> 60;
	tmp4 = a[3] + b[3] + carry;
	c[3] = tmp4 & mask_low;
}

void mp_mul_new( digit_t* a1,  digit_t* a0,  digit_t* b1,  digit_t* b0, digit_t* c_0, digit_t* c_1, digit_t* c_2, digit_t* c_3)
{// input 64bits  output 60bits x 4
	digit_t ci[6][4];
	digit_t d1[2], d2[2], e[4], f[4] = {0,0,0,0};
	digit_t mask_low = (digit_t)(-1) >> 4; 
	unsigned int carryOut, i, k, j, s; digit_t temp;
	

	for (i = 0; i < 6; i++) {
		digit_x_digit_2(a1[i], a0[i], b1[i], b0[i], ci[i]);//102x102  204bits -> 4x60bits
	} 
	c_0[0] = ci[0][0];
    c_1[0] = ci[0][1];
	c_2[0] = ci[0][2];
	c_3[0] = ci[0][3];
	
	for(i = 2; i < 6; i = i + 2){
		k = i>>1;
		for(j=0;j<k;j++){
			s = i-j;
			digit_a_digit_2(a1[j], a0[j], a1[s], a0[s], d1);//102+102 or 103+103  -> 64bitx2
			digit_a_digit_2(b1[j], b0[j], b1[s], b0[s], d2);
			digit_x_digit_2(d1[1], d1[0], d2[1], d2[0], e); // -> 4x60bits
			digit_s_digit_4(e, ci[j], e);
			digit_s_digit_4(e, ci[s], e);
			digit_a_digit_4(e, f, f);                       //one-level Karatsuba to reduce the complexity in the modular multiplication
		}
		digit_a_digit_4(ci[k],f,f);
		c_0[i] = f[0];
		c_1[i] = f[1];
		c_2[i] = f[2];
		c_3[i] = f[3];
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
		f[3] = 0;
	}   //e.g.  c4 = a4*b0 + a3*b1 + a2*b2 + a1*b3 + a0*b4;

	for(i = 1; i < 6; i = i + 2){
		for(j=0;j<=((i-1)>>1);j++){
			s = i-j;
			digit_a_digit_2(a1[j], a0[j], a1[s], a0[s], d1);//102+102 or 103+103  -> 64bitx2
			digit_a_digit_2(b1[j], b0[j], b1[s], b0[s], d2);
			digit_x_digit_2(d1[1], d1[0], d2[1], d2[0], e); // -> 4x60bits
			digit_s_digit_4(e, ci[j], e);
			digit_s_digit_4(e, ci[s], e);
			digit_a_digit_4(e, f, f);                       //one-level Karatsuba to reduce the complexity in the modular multiplication
		}
		c_0[i] = f[0];
		c_1[i] = f[1];
		c_2[i] = f[2];
		c_3[i] = f[3];
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
		f[3] = 0;
	}  

	for (i = 6; i < 10; i = i + 2) {
		k = i>>1;
		for (j = i-5; j < k; j++) {
			s = i-j;
			digit_a_digit_2(a1[j], a0[j], a1[s], a0[s], d1);//102+102 or 103+103  -> 64bitx2
			digit_a_digit_2(b1[j], b0[j], b1[s], b0[s], d2);
			digit_x_digit_2(d1[1], d1[0], d2[1], d2[0], e); // -> 4x60bits
			digit_s_digit_4(e, ci[j], e);
			digit_s_digit_4(e, ci[s], e);
			digit_a_digit_4(e, f, f);                       //one-level Karatsuba to reduce the complexity in the modular multiplication
        }
		digit_a_digit_4(ci[k],f,f);                        
		temp = (f[0]<<1) + c_0[i-6];   
		s = i-6;
		c_0[s] = temp & mask_low;
		carryOut = temp >> 60;
		temp = (f[1]<<1) + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> 60;
		temp = (f[2]<<1) + c_2[s] + carryOut;
		c_2[s] = temp & mask_low;
		carryOut = temp >> 60;
		c_3[s] = (f[3]<<1) + c_3[s] + carryOut;
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
		f[3] = 0;
	}                                                       // get the raw coefficients.

	for (i = 7; i < 10; i = i + 2) {
		for (j = i-5; j <= ((i-1)>>1); j++) {
			s = i - j;
			digit_a_digit_2(a1[j], a0[j], a1[s], a0[s], d1);//102+102 or 103+103  -> 64bitx2
			digit_a_digit_2(b1[j], b0[j], b1[s], b0[s], d2);
			digit_x_digit_2(d1[1], d1[0], d2[1], d2[0], e); // -> 4x60bits
			digit_s_digit_4(e, ci[j], e);
			digit_s_digit_4(e, ci[s], e);
			digit_a_digit_4(e, f, f);                       //one-level Karatsuba to reduce the complexity in the modular multiplication
        }                                              
		s = i-6;
		temp = (f[0]<<1) + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> 60;//
		temp = (f[1]<<1) + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> 60;//
		temp = (f[2]<<1) + c_2[s] + carryOut;
		c_2[s] = temp & mask_low;
		carryOut = temp >> 60;//
		c_3[s] = (f[3]<<1) + c_3[s] + carryOut;
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
		f[3] = 0;
	}                                                    // get the raw coefficients.

	temp = (ci[5][0]<<1) + c_0[4];
	c_0[4] = temp & mask_low;
	carryOut = temp >> 60;
	temp = (ci[5][1]<<1) + c_1[4] + carryOut;
	c_1[4] = temp & mask_low;
	carryOut = temp >> 60;
	temp = (ci[5][2]<<1) + c_2[4] + carryOut;   
	c_2[4] = temp & mask_low;                       
    carryOut = temp >> 60;            
	c_3[4] = (ci[5][3]<<1) + c_3[4] + carryOut;          // get the final raw coefficients.
}

void digit_sqr( digit_t a, digit_t* co)
{// Digit multiplication, digit * digit -> 2-digit result     
 // Inputs : 60bits a & b 
 // Output : 120bits c = a * b , c is divided into 2 x 60bits parts
    register digit_t al, ah, temp;  
    digit_t alal, alah, ahah, res1, carry;
    digit_t mask_low = (digit_t)(-1) >> 34;
	digit_t mask_low_1 = (digit_t)(-1) >> 35;
	digit_t mask_low_2 = (digit_t)(-1) >> 4;
	digit_t k = 30;

    al = a & mask_low;    // Low part
    ah = a >> k;          // High part

    alal = al*al;
	ahah = ah*ah;
	alah = al*ah;
    
    //co[0] = alal & mask_low;                         // C00

   // res1 = alal >> k;
    temp = alal + ((alah & mask_low_1)<<31);
    carry = temp >> 60;
    co[0] = temp & mask_low_2;                // C01

    co[1] = ahah + carry + (alah>>29);
}

void digit_sqr_2( digit_t a1,  digit_t a0, digit_t* c)
{
 // input 60bits x 2
 // output 60bits x 4
	digit_t alal[2], ahah[2], alah[2];
	digit_t tmp, carry ;
	digit_t mask_low = (digit_t)(-1) >> 4; // 4x60bit
	
	digit_sqr(a0, alal);
	digit_x_digit_1(a0, a1, alah);
	digit_sqr(a1, ahah);
	
	c[0] = alal[0];
	
	tmp = alal[1] + (alah[0] << 1);
	c[1] = tmp & mask_low;
	carry = tmp >> 60;
	
	tmp = (alah[1] << 1) + ahah[0] + carry;
	c[2] = tmp & mask_low;
	carry = tmp >> 60;
	
	c[3] = ahah[1] + carry;
}

void mp_sqr_new( digit_t* a1,  digit_t* a0, digit_t* c_0, digit_t* c_1, digit_t* c_2, digit_t* c_3)
{// input 64bits  output 60bits x 4
	digit_t ci[6][4];
	digit_t f[4] = {0,0,0,0}, e[4] = {0,0,0,0};
	digit_t mask_low = (digit_t)(-1) >> 4; 
	unsigned int carryOut, i, k, j, s; digit_t temp;

	for (i = 0; i < 6; i++) {
		digit_sqr_2(a1[i], a0[i], ci[i]);//102x102  204bits -> 4x60bits
	} 
	c_0[0] = ci[0][0];
    c_1[0] = ci[0][1];
	c_2[0] = ci[0][2];
	c_3[0] = ci[0][3];
	
	for(i = 2; i < 6; i = i + 2){
		k = i>>1;
		for(j=0;j<k;j++){
			s = i-j;
			digit_x_digit_2(a1[j], a0[j], a1[s], a0[s], e);
			digit_a_digit_4(e, f, f);
		}
		digit_a_digit_4(f, f, e);
		digit_a_digit_4(ci[k],e,f);
		c_0[i] = f[0];
		c_1[i] = f[1];
		c_2[i] = f[2];
		c_3[i] = f[3];
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
		f[3] = 0;
	}   //e.g.  c4 = a4*b0 + a3*b1 + a2*b2 + a1*b3 + a0*b4;

	for(i = 1; i < 6; i = i + 2){
		for(j=0;j<=((i-1)>>1);j++){
			s = i-j;
			digit_x_digit_2(a1[j], a0[j], a1[s], a0[s], e);
			digit_a_digit_4(e, f, f);
		}
		digit_a_digit_4(f, f, e);
		c_0[i] = e[0];
		c_1[i] = e[1];
		c_2[i] = e[2];
		c_3[i] = e[3];
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
		f[3] = 0;
	}  

	for (i = 6; i < 10; i = i + 2) {
		k = i>>1;
		for (j = i-5; j < k; j++) {
			s = i-j;
			digit_x_digit_2(a1[j], a0[j], a1[s], a0[s], e);
			digit_a_digit_4(e, f, f);
        }
		digit_a_digit_4(f, f, e);
		digit_a_digit_4(ci[k],e,f);                        
		temp = (f[0]<<1) + c_0[i-6];   
		s = i-6;
		c_0[s] = temp & mask_low;
		carryOut = temp >> 60;
		temp = (f[1]<<1) + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> 60;
		temp = (f[2]<<1) + c_2[s] + carryOut;
		c_2[s] = temp & mask_low;
		carryOut = temp >> 60;
		c_3[s] = (f[3]<<1) + c_3[s] + carryOut;
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
		f[3] = 0;
	}                                                       // get the raw coefficients.

	for (i = 7; i < 10; i = i + 2) {
		for (j = i-5; j <= ((i-1)>>1); j++) {
			s = i - j;
			digit_x_digit_2(a1[j], a0[j], a1[s], a0[s], e);
			digit_a_digit_4(e, f, f);
        }                                                 
		digit_a_digit_4(f, f, e);
		s = i-6;
		temp = (e[0]<<1) + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> 60;//
		temp = (e[1]<<1) + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> 60;//
		temp = (e[2]<<1) + c_2[s] + carryOut;
		c_2[s] = temp & mask_low;
		carryOut = temp >> 60;//
		c_3[s] = (e[3]<<1) + c_3[s] + carryOut;
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
		f[3] = 0;
	}                                                    // get the raw coefficients.
	temp = (ci[5][0]<<1) + c_0[4];
	c_0[4] = temp & mask_low;
	carryOut = temp >> 60;
	temp = (ci[5][1]<<1) + c_1[4] + carryOut;
	c_1[4] = temp & mask_low;
	carryOut = temp >> 60;
	temp = (ci[5][2]<<1) + c_2[4] + carryOut;   
	c_2[4] = temp & mask_low;                       
    carryOut = temp >> 60;            
	c_3[4] = (ci[5][3]<<1) + c_3[4] + carryOut;          // get the final raw coefficients.
}

void fpmul( digit_t* a1,  digit_t* a0,  digit_t* b1,  digit_t* b0, digit_t* d0, digit_t* d1)
{
	digit_t c0[6], c1[6], c2[6], c3[6];
	mp_mul_new(a1, a0, b1, b0, c0, c1, c2, c3);
	improved_IFFM(c0, c1, c2, c3, d0, d1);
}

void fpsqr( digit_t* a1,  digit_t* a0, digit_t* d0, digit_t* d1)
{
	digit_t c0[6], c1[6], c2[6], c3[6];
	mp_sqr_new(a1,a0,c0, c1, c2, c3);
	improved_IFFM(c0, c1, c2, c3, d0, d1);
}

void fpinv_chain_new_u(digit_t* a1, digit_t* a0, digit_t* c0, digit_t* c1)
{
	digit_t t[31][2][6], tt[2][6];
	int i,j;

    // Precomputed table
    fpsqr(a1, a0, tt[0], tt[1]);
    fpmul(a1, a0, tt[1], tt[0], t[0][0], t[0][1]);
    for (i = 0; i <= 29; i++) fpmul(t[i][1], t[i][0], tt[1], tt[0], t[i+1][0], t[i+1][1]);

//    fpcopy(a, tt);

	for (i=0;i<6;i++){
		tt[0][i] = a0[i];
		tt[1][i] = a1[i];
	}

    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[6][1], t[6][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[30][1], t[30][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[25][1], t[25][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[28][1], t[28][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[7][1], t[7][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 11; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[11][1], t[11][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(a1, a0, tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[0][1], t[0][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[3][1], t[3][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[16][1], t[16][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[24][1], t[24][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[28][1], t[28][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[16][1], t[16][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[4][1], t[4][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[3][1], t[3][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[20][1], t[20][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[11][1], t[11][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[14][1], t[14][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[15][1], t[15][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[0][1], t[0][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[15][1], t[15][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[19][1], t[19][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[9][1], t[9][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[5][1], t[5][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[27][1], t[27][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[28][1], t[28][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[29][1], t[29][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[1][1], t[1][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[3][1], t[3][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[2][1], t[2][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[30][1], t[30][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[25][1], t[25][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[28][1], t[28][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[22][1], t[22][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[3][1], t[3][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[22][1], t[22][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[7][1], t[7][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[9][1], t[9][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[4][1], t[4][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[20][1], t[20][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 11; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[10][1], t[10][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[26][1], t[26][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 11; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[2][1], t[2][0], tt[1], tt[0], tt[0], tt[1]);
    for (j = 0; j < 50; j++) {
        for (i = 0; i < 6; i++)fpsqr(tt[1], tt[0], tt[0], tt[1]);
        fpmul(t[30][1], t[30][0], tt[1], tt[0], tt[0], tt[1]);
    }
//    fpcopy(tt, a);    
	for (i=0;i<6;i++){
		c0[i] = tt[0][i];
		c1[i] = tt[1][i];
	}
}

void fpinv_new_u(digit_t* a1, digit_t* a0, digit_t* c0, digit_t* c1)
{
	digit_t tt[2][6], d0[6], d1[6], e0[6], e1[6];
	int i;
	
	for (i = 0; i < 6; i++){
		tt[0][i] = a0[i];
		tt[1][i] = a1[i];
	}
	
	fpinv_chain_new_u(tt[1], tt[0], d0, d1);
	
	fpsqr(d1, d0, e0, e1);
	fpsqr(e1, e0, d0, d1);
	fpmul(a1, a0, d1, d0, c0, c1);
}

void fp_add_new(digit_t* a1, digit_t* a0, digit_t* b1, digit_t* b0, digit_t* c0, digit_t* c1)
{// Modular Addition based on the new data representation , with the unconditional radix
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/2 - 1], a[i] & b[i] are presented in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/2 - 1], all in 60 bits , to compute c = a + b
	unsigned int i,j;
	digit_t ca = 0, ca1, t[2], d[6][2];
	digit_t R[2]  = {0x408000000000000, 0x34AA7F10F1F};//34AA7F10F1F 408000000000000
	digit_t R_div_2[2] = {0xA04000000000000, 0x1A553F8878F};//1A553F8878F A04000000000000
 
	digit_t mask_low = (digit_t)(-1) >> 4;
	
    for(i = 0; i < 6; i++){
		digit_a_digit_2(a1[i], a0[i], b1[i], b0[i], d[i]);
	}

	for (i = 1; i < 6; i++){
		t[0] = d[i-1][0] - R[0];
		ca1 = t[0] >> 63;
		t[1] = d[i-1][1] - R[1] - ca1;
		j = t[1] >> 63;
		ca = 1 ^ j;
		
		d[i-1][0] = (t[0] + R[0]*j) & mask_low;
		d[i-1][1] = t[1] + R[1] * j + ca1*j;
		
		ca1 = d[i][0] + ca;
		d[i][0] = ca1 & mask_low;
		ca = ca1 >> 60;
		d[i][1] = d[i][1] + ca;
		
		c0[i-1] = d[i-1][0];
		c1[i-1] = d[i-1][1];
		c0[i] = d[i][0];
		c1[i] = d[i][1];
	} 
	t[0] = d[5][0] - R_div_2[0];
	ca1 = t[0] >> 63;
	t[1] = d[5][1] - R_div_2[1] - ca1;
	j = t[1] >> 63;
	ca = 1 ^ j;
	
	c0[5] = (t[0] + R_div_2[0]*j) & mask_low;
	c1[5] = t[1] + R_div_2[1] * j + ca1*j;
	ca1 = c0[0] + ca;
	c0[0] = ca1 & mask_low;
	ca = ca1 >> 60;
	c1[0] = c1[0] + ca;
}

void fpneg_new(digit_t* a1, digit_t* a0, digit_t* b0, digit_t* b1)
{// Modular Negation based on the new data representation , with the unconditional radix
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/2 - 1], a[i] all in 60 bits
 // Output : b[i] in [0, R-1] for 0≤i<n-1; b[n-1] in [0, R/2 - 1], b[i] all in 60 bits, to compute b = -a
	unsigned int i;
	digit_t R[2]  = {0x408000000000000, 0x34AA7F10F1F};//34AA7F10F1F 408000000000000
	digit_t R_div_2[2] = {0xA04000000000000, 0x1A553F8878F};//1A553F8878F A04000000000000
	digit_t k = 5, tmp, carry;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	for( i = 0; i < k; i++){
		tmp = R[0] - a0[i] - 1;
		b0[i] = tmp & mask_low;
		carry = tmp >> 63;
		b1[i] = R[1] - a1[i] - carry;
	}

	tmp = R_div_2[0] - a0[k] - 1;
	b0[k] = tmp & mask_low;
	carry = tmp >> 63;
	b1[k] = R_div_2[1] - a1[k] - carry;
}

void fpsub_new(digit_t* a1, digit_t* a0, digit_t* b1, digit_t* b0, digit_t* c0, digit_t* c1)
{// Modular Subtraction based on the new data representation , with the unconditional radix
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/2 - 1], a[i] & b[i] are presented in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/2 - 1], all in 60 bits , to compute c = a - b
	unsigned int i, j = 0, t;
	digit_t R[2]  = {0x408000000000000, 0x34AA7F10F1F};//34AA7F10F1F 408000000000000
	digit_t R_div_2[2] = {0xA04000000000000, 0x1A553F8878F};//1A553F8878F A04000000000000
	digit_t mask_low = (digit_t)(-1) >> 4, tmp, borrow;
	
	for(i=0;i<6;i++){
		tmp = a0[i] - b0[i];
		c0[i] = tmp & mask_low;
		borrow = tmp >> 63;
		c1[i] = a1[i] - b1[i] - borrow;
	}
	t = (c0[0] == 0)&(c1[0] == 0);
	j = c1[0] >> 63;
	t = t | j;
	tmp = c0[0] + t * R[0];
	c0[0] = tmp & mask_low;
	borrow = tmp >> 60;
	c1[0] = c1[0] + t * R[1] + borrow;

	tmp = c0[1] - t;
	c0[1] = tmp & mask_low;
	borrow = tmp >> 63;
	c1[1] = c1[1] - borrow;
	
	for(i=1;i<5;i++){
		j = c1[i]>>63;
		tmp = c0[i] + j * R[0] ;
		c0[i] = tmp & mask_low;
		borrow = tmp >> 60;
		c1[i] = c1[i] + j * R[1] + borrow;
		
		tmp = c0[i+1] - j;
		c0[i+1] = tmp & mask_low;
		borrow = tmp >> 63;
		c1[i+1] = c1[i+1] - borrow;
	}
	
	j = c1[5]>>63;
	tmp = c0[5] + j * R_div_2[0] ;
	c0[5] = tmp & mask_low;
	borrow = tmp >> 60;
	c1[5] = c1[5] + j * R_div_2[1] + borrow;

	c0[0] = c0[0] - j;
}

void fpdiv2_610_new(digit_t* a1, digit_t* a0, digit_t* c0, digit_t* c1)
{// Modular Division by Two based on the new data representation , with the unconditional radix
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/2 - 1], a[i] all in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/2 - 1], c[i] all in 60 bits, to compute c = a/2
    unsigned int i;
    digit_t mask;
	digit_t R_div_2[2] = {0xA04000000000000, 0x1A553F8878F};//1A553F8878F A04000000000000
	digit_t R_div_4[2] =  {0xD02000000000000, 0x0D2A9FC43C7};
	digit_t b[6][2] ={0}, tmp, carry;
	digit_t mask_low = (digit_t)(-1) >> 4;

	for (i = 1; i < 6; i++) {
		mask = (digit_t)(a0[i] & 1);		
//		a0[i] -= mask;
		b[i-1][0] = R_div_2[0] & ((digit_t)(0-mask));
		b[i-1][1] = R_div_2[1] & ((digit_t)(0-mask));
		tmp = (a0[i] >> 1)^(a1[i] << 59);
		a0[i] = tmp & mask_low;
		a1[i] = a1[i] >> 1;
	}

	mask = (digit_t)(a0[0] & 1);		
//	a0[0] -= mask;
	b[5][0] = R_div_4[0] & ((digit_t)(0-mask));
	b[5][1] = R_div_4[1] & ((digit_t)(0-mask));
	tmp = (a0[0] >> 1)^(a1[0] << 59);
	a0[0] = tmp & mask_low;
	a1[0] = a1[0] >> 1;

    for (i = 0; i < 6; i++) {
		tmp = a0[i] + b[i][0];
		c0[i] = tmp & mask_low;
		carry = tmp >> 60;
		c1[i] = a1[i] + b[i][1] + carry;
	}
}

void digit_a_digit_1( digit_t a1,  digit_t a0,  digit_t b1,  digit_t b0, digit_t* c)
{// input 60bits x 2
 // output 60bits x 3
	digit_t tmp, carry;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	tmp = a0 + b0;
	c[0] = tmp & mask_low;
	carry = tmp >> 60;
	tmp = a1 + b1 + carry;
	c[1] = tmp & mask_low;
	c[2] = tmp >> 60;
}

void from_uRadix(digit_t* a1, digit_t* a0, digit_t* co)
{
	digit_t R[2]  = {0x408000000000000, 0x34AA7F10F1F};//34AA7F10F1F 408000000000000
	digit_t R_sub_1[2] = {0x407FFFFFFFFFFFF, 0x34AA7F10F1F};//34AA7F10F1F 407FFFFFFFFFFFF
	digit_t R_div_2_sub1[2] = {0xA03FFFFFFFFFFFF, 0x1A553F8878F};//1A553F8878F A03FFFFFFFFFFFF
	
	digit_t b[4] = {0}, d[4], e[3], f[3], g[3], h[3], zero[1] = {0}, c[11], carry[1];
	digit_t mask_low_1 = (digit_t)(-1) >> 4; int i, j, k;
	
	j = 1;
	for(i=1;i<5;i++){
		k = (a1[i] == R_sub_1[1])&(a0[i]==R_sub_1[0]);
		j = j & k;
	}
	
	k = (a1[5] == R_div_2_sub1[1])&(a0[5]==R_div_2_sub1[0]);
	j = j & k;
	
	k = ((a1[0] == R[1])&(a0[0]==R[0]))|((a1[0] == R_sub_1[1])&(a0[0]==R_sub_1[0]));
	j = j & k;
	
	k = (digit_t)(a0[0] & 1); 
	
	for(i=1;i<6;i++){
			a1[i] = a1[i] * (1 ^ j);
			a0[i] = a0[i] * (1 ^ j);
	}
	a1[0] = a1[0] * (1 ^ j) + (1 ^ k) * j;
	a0[0] = a0[0] * (1 ^ j) + (1 ^ k) * j;
	
	digit_x_digit_2(a1[5], a0[5], R[1], R[0], b);
	digit_a_digit_1(b[1], b[0], a1[4], a0[4], f);
	digit_a_digit_1(b[3], b[2], zero[0], f[2], e);
	c[0] = f[0]; c[1] = f[1];
	c[2] = e[0]; c[3] = e[1];
	
	digit_x_digit_2(c[1], c[0], R[1], R[0], b);
	digit_x_digit_2(c[3], c[2], R[1], R[0], d);
	digit_a_digit_1(b[3], b[2], d[1], d[0], e);
	digit_a_digit_1(d[3], d[2], zero[0], e[2], g);
	//c[0] = b[0]; c[1] = b[1]; c[2] = e[0]; c[3] = e[1]; c[4] = g[0]; c[5] = g[1];
	digit_a_digit_1(b[1], b[0], a1[3], a0[3], f);
	c[0] = f[0]; c[1] = f[1]; 
	digit_a_digit_1(e[1], e[0], zero[0], f[2], h);
	c[2] = h[0]; c[3] = h[1];
	digit_a_digit_1(g[1], g[0], zero[0], h[2], e);
	c[4] = e[0]; c[5] = e[1];
	
	digit_x_digit_2(c[1], c[0], R[1], R[0], b);
	digit_x_digit_2(c[3], c[2], R[1], R[0], d);
	c[0] = b[0]; c[1] = b[1];
	digit_a_digit_1(b[3], b[2], d[1], d[0], e);
	c[2] = e[0]; c[3] = e[1];
	digit_x_digit_2(c[5], c[4], R[1], R[0], b);
	digit_a_digit_1(d[3], d[2], b[1], b[0], f);
	digit_a_digit_1(f[1], f[0], zero[0], e[2], g);
	c[4] = g[0]; c[5] = g[1]; c[6] = b[2] + f[2] + g[2];
	digit_a_digit_1(c[1], c[0], a1[2], a0[2], f);
	c[0] = f[0]; c[1] = f[1]; 
	digit_a_digit_1(c[3], c[2], zero[0], f[2], h);
	c[2] = h[0]; c[3] = h[1];
	digit_a_digit_1(c[5], c[4], zero[0], h[2], e);
	c[4] = e[0]; c[5] = e[1]; c[6] = c[6] + e[2];
	
	digit_x_digit_2(c[1], c[0], R[1], R[0], b);
	digit_x_digit_2(c[3], c[2], R[1], R[0], d);
	c[0] = b[0]; c[1] = b[1];
	digit_a_digit_1(b[3], b[2], d[1], d[0], e);
	c[2] = e[0]; c[3] = e[1];
	digit_x_digit_2(c[5], c[4], R[1], R[0], b);
	digit_a_digit_1(d[3], d[2], b[1], b[0], f);
	digit_a_digit_1(f[1], f[0], zero[0], e[2], g);
	c[4] = g[0]; c[5] = g[1];
	digit_x_digit_2(zero[0], c[6], R[1], R[0], d);
	digit_a_digit_1(b[3], b[2], d[1], d[0], h);
	carry[0] = f[2] + g[2];
	digit_a_digit_1(h[1], h[0], zero[0], carry[0], g);
	c[6] = g[0]; c[7] = g[1]; c[8] = d[2] + g[2] + h[2];
	digit_a_digit_1(c[1], c[0], a1[1], a0[1], f);
	c[0] = f[0]; c[1] = f[1]; 
	digit_a_digit_1(c[3], c[2], zero[0], f[2], h);
	c[2] = h[0]; c[3] = h[1];
	digit_a_digit_1(c[5], c[4], zero[0], h[2], e);
	c[4] = e[0]; c[5] = e[1]; 
	digit_a_digit_1(c[7], c[6], zero[0], e[2], f);
	c[6] = f[0]; c[7] = f[1]; c[8] = c[8] + f[2];
	
	digit_x_digit_2(c[1], c[0], R[1], R[0], b);
	digit_x_digit_2(c[3], c[2], R[1], R[0], d);
	c[0] = b[0]; c[1] = b[1];
	digit_a_digit_1(b[3], b[2], d[1], d[0], e);
	c[2] = e[0]; c[3] = e[1];
	digit_x_digit_2(c[5], c[4], R[1], R[0], b);
	digit_a_digit_1(d[3], d[2], b[1], b[0], f);
	digit_a_digit_1(f[1], f[0], zero[0], e[2], g);
	c[4] = g[0]; c[5] = g[1];
	digit_x_digit_2(c[7], c[6], R[1], R[0], d);
	digit_a_digit_1(b[3], b[2], d[1], d[0], h);
	carry[0] = f[2] + g[2];
	digit_a_digit_1(h[1], h[0], zero[0], carry[0], g);
	c[6] = g[0]; c[7] = g[1];
	digit_x_digit_2(zero[0], c[8], R[1], R[0], b);
	digit_a_digit_1(d[3], d[2], b[1], b[0], f);
	carry[0] = h[2] + g[2];
	digit_a_digit_1(f[1], f[0], zero[0], carry[0], e);
	c[8] = e[0]; c[9] = e[1]; c[10] = b[2] + f[2] + e[2];
	digit_a_digit_1(c[1], c[0], a1[0], a0[0], f);
	c[0] = f[0]; c[1] = f[1]; 
	digit_a_digit_1(c[3], c[2], zero[0], f[2], h);
	c[2] = h[0]; c[3] = h[1];
	digit_a_digit_1(c[5], c[4], zero[0], h[2], e);
	c[4] = e[0]; c[5] = e[1]; 
	digit_a_digit_1(c[7], c[6], zero[0], e[2], f);
	c[6] = f[0]; c[7] = f[1];
	digit_a_digit_1(c[9], c[8], zero[0], f[2], h);
	c[8] = h[0]; c[9] = h[1]; c[10] = c[10] + h[2];
	
	for (i=0;i<10;i++){
		co[i] = (c[i] >> (4*i))^(c[i+1] << (60-4*i));
	}
}

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

void IBR_uRadix(const digit_t* a, digit_t* q, digit_t* r, const unsigned int unwords, const digit_t* lamda)
{
	digit_t t[8] = {0}, temp[16], res1, res2, res3, tmp[2], s, t_r;
	unsigned int i, j;
	digit_t mask_low = (digit_t)(-1) >> 12;  //low 52
	digit_t m = 0x0006954FE21E3E81;
	digit_t mask_low_2 = (digit_t)(-1) >> 4; //low 60
	
	s = a[0] & ((digit_t)(-1) >> 13);
	t_r = (a[0] >> 51)^(a[1] << 13);
	
	if(unwords == 5){
		for(i=0;i < 8;i++){
			t[i] = (a[i+1] >> 36)^(a[i+2] << 28);//t8 -> lamda8
		}
		mp_mul_old(t, lamda, temp, 8);
	}
	else if (unwords == 4){
		for (i=0;i < 6;i++){
			t[i] = (a[i+1] >> 36)^(a[i+2] << 28);//t7 -> lamda7
		}
		t[6] = a[7] >> 36;
		mp_mul_old(t, lamda, temp, 7);
	}
	else if (unwords == 3){
		for (i=0;i < 5;i++){
			t[i] = (a[i+1] >> 36)^(a[i+2] << 28);//t5 -> lamda5
		}
		mp_mul_old(t, lamda, temp, 5);
	}
	else if (unwords == 2){
		for (i=0;i < 3;i++){
			t[i] = (a[i+1] >> 36)^(a[i+2] << 28);//t4 -> lamda4
		}
		t[3] = a[4] >> 36;
		mp_mul_old(t, lamda, temp, 4);
	}
	else /*if (unwords == 1)*/{
		for (i=0;i < 2;i++){
			t[i] = (a[i+1] >> 36)^(a[i+2] << 28);//t2 -> lamda2
		}
		mp_mul_old(t, lamda, temp, 2);
	}
//	else{
//		t[0] = a[1] >> 36;                       //t1 -> lamda1
//		mp_mul_old(t, lamda, temp, 1);
//	}
	
	if(unwords == 5){
		for(i=0;i<8;i++){
			res1 = temp[i+7]>>63;//
			res2 = temp[i+8]<<1;//
			q[i] = res1 ^ res2;
		}//8
		for(i=8;i<10;i++){
			q[i] = 0;
		}
	}
	else if (unwords == 4){
		for(i=0;i<7;i++){
			res1 = temp[i+6]>>26;
			res2 = temp[i+7]<<38;
			q[i] = res1 ^ res2;
		}//7
		for(i=7;i<10;i++){
			q[i] = 0;
		}
	}
	else if (unwords == 3){
		for(i=0;i<5;i++){
			res1 = temp[i+4]>>53;
			res2 = temp[i+5]<<11;
			q[i] = res1 ^ res2;
		}//5
		for(i=5;i<10;i++){
			q[i] = 0;
		}
	}
	else if (unwords == 2){
		for(i=0;i<4;i++){
			res1 = temp[i+3]>>16;
			res2 = temp[i+4]<<48;
			q[i] = res1 ^ res2;
		}//4
		for(i=4;i<10;i++){
			q[i] = 0;
		}
	}
	else /*if (unwords == 1)*/{
		for(i=0;i<2;i++){
			res1 = temp[i+1]>>43;
			res2 = temp[i+2]<<21;
			q[i] = res1 ^ res2;
		}//2
		for(i=2;i<10;i++){
			q[i] = 0;
		}
	}
//	else {
//		res1 = temp[0]>>6;
//		res2 = temp[1]<<58;
//		q[0] = res1 ^ res2;
//		for(i=1;i<8;i++){
//			q[i] = 0;
//		}//1
//	}
	
	res1 = q[0] & mask_low;
	digit_x_digit_1(res1, m, tmp);
	r[0] = ((t_r & mask_low) - (tmp[0] & mask_low)) & mask_low;
	
	res1 = r[0] - m;
	res2 = res1 >> 63;
	res3 = 1 - res2;
	r[0] = res1 + (m * res2);
	//q+res3
	res1 = q[0] + res3;
	j = (q[0]==(digit_t)(-1));
	q[0] = res1 * (1 ^ j);
	
	int h;
	if (unwords == 5){
		h = 8;
	}
	else if (unwords == 4){
		h = 7;
	}
	else if (unwords == 3){
		h = 5;
	}
	else if (unwords == 2){
		h = 4;
	}
	else /*if (unwords == 1)*/{
		h = 2;
	}
//	else{
//		h = 1;
//	}
	
	
	for (i = 1; i < h; i++){	
		res1 = q[i] + j;	
		j = (q[i]==(digit_t)(-1));
		q[i] = res1 * (1 ^ j);	
	}
	q[h] = q[h] + j;
	
	res1 = r[0] << 51;
	r[1] = r[0] >> 9;
	r[0] = (res1 ^ s) & mask_low_2;
}

void to_uRadix(digit_t* a, digit_t* c0, digit_t* c1){
	int i, j;
	digit_t q[10], temp[2], tmp;
	digit_t mask_low = (digit_t)(-1) >> 4; //low 60
	digit_t lamda[6][8]={{0xaad989d5c112ce6c, 0x28f7a5113900cbd5, 0xbe19312a7b3d59be, 0xa10328bcee4a71f4, 0x5a6afbd1aa397d1d, 0x1b65816415aa11f4, 0xf9d0732975d88635, 0x26e2f610e161b674},{0x53d9eacdf147bd28, 0xe772538fa5f0c989, 0x8d51cbe8ed081945, 0x20ad508fa2d357de, 0x4baec431a8db2c0b, 0x870b0db3a7ce8399, 0x13717b0},{0x476840ca2f3b929c, 0x7d169abef46a8e5f, 0x8d46d96059056a84, 0x9d3e741cca5d7621, 0x9b8bd8438586d},{0x02c82b5423e8b4d5, 0xe652ebb10c6a36cb, 0xec21c2c36ce9f3a0, 0x4dc5},{0x1b674f9d0732975d, 0x26e2f610e16},{0x13}};
	
    for (i=0; i<5; i++){
		IBR_uRadix(a, q, temp, 5-i, lamda[i]);
		c0[i] = temp[0];
		c1[i] = temp[1];
		
		if (i==0){
			for(j=0;j<8;j++){
				a[j] = q[j];
			}
			for(j=8;j<10;j++){
				a[j] = 0;
			}
		}
		else if (i==1){
			for(j=0;j<7;j++){
				a[j] = q[j];
			}
			for(j=7;j<10;j++){
				a[j] = 0;
			}
		}
		else if (i==2){
			for(j=0;j<5;j++){
				a[j] = q[j];
			}
			for(j=5;j<10;j++){
				a[j] = 0;
			}
		}
		else if (i==3){
			for(j=0;j<4;j++){
				a[j] = q[j];
			}
			for(j=4;j<10;j++){
				a[j] = 0;
			}
		}
		else{
			for(j=0;j<2;j++){
				a[j] = q[j];
			}
			for(j=2;j<10;j++){
				a[j] = 0;
			}
		}
	}
	tmp = q[0] >> 60;
	c1[5] = (q[1] << 4)^tmp;
	c0[5] = q[0] & mask_low;
}

void fp2sqr_uRadix(f2elm_t a, f2elm_t ai, f2elm_t c, f2elm_t ci){
	f2elm_t t0;
	felm_t tt1, tt2, tt3, tt4;
	fp_add_new(a[1], a[0], ai[1], ai[0], t0[0], t0[1]);//a+ai
	fpsub_new(a[1], a[0], ai[1], ai[0], tt1, tt2);     //a-ai
	fp_add_new(a[1], a[0], a[1], a[0], tt3, tt4);      //2a
	fpmul(t0[1], t0[0], tt2, tt1, c[0], c[1]);         
	fpmul(tt4, tt3, ai[1], ai[0], ci[0], ci[1]);
}

void fp2add_new(f2elm_t a, f2elm_t ai, f2elm_t b, f2elm_t bi, f2elm_t c, f2elm_t ci)
{
	fp_add_new(a[1], a[0], b[1], b[0], c[0], c[1]);
	fp_add_new(ai[1], ai[0], bi[1], bi[0], ci[0], ci[1]);
}

void fp2sub_new(f2elm_t a, f2elm_t ai, f2elm_t b, f2elm_t bi, f2elm_t c, f2elm_t ci)
{
	fpsub_new(a[1], a[0], b[1], b[0], c[0], c[1]);
	fpsub_new(ai[1], ai[0], bi[1], bi[0], ci[0], ci[1]);
}

void fp2neg_new(f2elm_t a, f2elm_t ai, f2elm_t b, f2elm_t bi)
{
	fpneg_new(a[1], a[0], b[0], b[1]);
	fpneg_new(ai[1], ai[0], bi[0], bi[1]);
}

void fp2mul_uRadix(f2elm_t a, f2elm_t ai, f2elm_t b, f2elm_t bi, f2elm_t c, f2elm_t ci)
{// GF(p^2) multiplication, c = a*b in GF(p^2).
 // Inputs : a = a + ai*i, b = b + bi*i. a, ai, b, bi all meet the same requirement with fpmul_new & fpsqr_new
 // Output : c = c + ci*i. c = a * b, c and ci meet the same requirement with fpmul_new & fpsqr_new
 // a-b+2^4*R*P mod P = (a_5-b_5+2^4*(R/2-1)*R)*R^5 + (a_4-b_4+2^4*(R-1)*R)*R^4 + ... + a_0-b_0+2^4*(R-1)*R
    digit_t PR_5[4] = {0x0F80000000000000, 0x079FEB5580EF0E0B, 0x0F1EF2644F18F2FD, 0x56ADA95}; // 2^4*(R/2-1)*R) 56ADA95 F1EF2644F18F2FD 79FEB5580EF0E0B F80000000000000
    digit_t PR[4] =   {0x0F80000000000000, 0x0F400B5580EF0E0B, 0x0E3DE4C89E31E5FA, 0xAD5B52B}; // 2^4*(R-1)*R AD5B52B E3DE4C89E31E5FA F400B5580EF0E0B F80000000000000
    digit_t tt1[4][6], tt2[4][6], tt3[4][6]; digit_t t0[6], t1[6], h0[6], h1[6];
	int i, j;
	digit_t temp, res;
	digit_t mask_low_1 = (digit_t)(-1) >> 4;
	digit_t k = 60;
	
	for (i=0; i<6; i++){
		temp = a[0][i] + ai[0][i];
		t0[i] = temp & mask_low_1;
		res = temp >> k;
		t1[i] = a[1][i] + ai[1][i] + res; //a+ai
		
		temp = b[0][i] + bi[0][i];
		h0[i] = temp & mask_low_1;
		res = temp >> k;
		h1[i] = b[1][i] + bi[1][i] + res; //b+bi
	}
	
	//for c[0]
    mp_mul_new(a[1], a[0], b[1], b[0], tt1[0], tt1[1], tt1[2], tt1[3]); // tt1 = a*b
	mp_mul_new(ai[1], ai[0], bi[1], bi[0], tt2[0], tt2[1], tt2[2], tt2[3]); // tt2 = ai*bi

	for (i=0; i<5; i++){
		//a*b + 2^4*R*P - ai*bi
		temp = tt1[0][i] + PR[0];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> k;
		temp = tt1[1][i] + PR[1] + res;
		tt3[1][i] = temp & mask_low_1;
		res = temp >> k;
		temp = tt1[2][i] + PR[2] + res;
		tt3[2][i] = temp & mask_low_1;
		res = temp >> k;
		tt3[3][i] = tt1[3][i] + PR[3] + res;
		
		temp = tt3[0][i] - tt2[0][i];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> 63;
		temp = tt3[1][i] - tt2[1][i] - res;
		tt3[1][i] = temp & mask_low_1;
		res = temp >> 63;
		temp = tt3[2][i] - tt2[2][i] - res;
		tt3[2][i] = temp & mask_low_1;
		res = temp >> 63;
		tt3[3][i] = tt3[3][i] - tt2[3][i] - res;
	}
	
	//a*b + 2^4*R*P - ai*bi
	temp = tt1[0][5] + PR_5[0];
	tt3[0][5] = temp & mask_low_1;
	res = temp >> k;
	temp = tt1[1][5] + PR_5[1] + res;
	tt3[1][5] = temp & mask_low_1;
	res = temp >> k;
	temp = tt1[2][5] + PR_5[2] + res;
	tt3[2][5] = temp & mask_low_1;
	res = temp >> k;
	tt3[3][5] = tt1[3][5] + PR_5[3] + res;
		
	temp = tt3[0][5] - tt2[0][5];
	tt3[0][5] = temp & mask_low_1;
	res = temp >> 63;
	temp = tt3[1][5] - tt2[1][5] - res;
	tt3[1][5] = temp & mask_low_1;
	res = temp >> 63;
	temp = tt3[2][5] - tt2[2][5] - res;
	tt3[2][5] = temp & mask_low_1;
	res = temp >> 63;
	tt3[3][5] = tt3[3][5] - tt2[3][5] - res;
	
	improved_IFFM(tt3[0], tt3[1], tt3[2], tt3[3], c[0], c[1]);// c[0]
	
	//for c[1]
	
	mp_mul_new(t1, t0, h1, h0, tt3[0], tt3[1], tt3[2], tt3[3]); //   (a+ai)*(b+bi)	
	
	// c[1] = (a+ai)*(b+bi) - (a*b + ai*bi) 
	for (i=0; i<6; i++) {
		//a*b + ai*bi
		temp = tt2[0][i] + tt1[0][i];
		tt2[0][i] = temp & mask_low_1;
		res = temp >> k;
		temp = tt2[1][i] + tt1[1][i] + res;
		tt2[1][i] = temp & mask_low_1;
		res = temp >> k;
		temp = tt2[2][i] + tt1[2][i] + res;
		tt2[2][i] = temp & mask_low_1;
		res = temp >> k;
		tt2[3][i] = tt2[3][i] + tt1[3][i] + res;
		//(a+ai)*(b+bi) - a*b - ai*bi
		temp = tt3[0][i] - tt2[0][i];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> 63;		
		temp = tt3[1][i] - tt2[1][i] - res;
		tt3[1][i] = temp & mask_low_1;
		res = temp >> 63;
		temp = tt3[2][i] - tt2[2][i] - res;
		tt3[2][i] = temp & mask_low_1;
		res = temp >> 63;
		tt3[3][i] = tt3[3][i] - tt2[3][i] - res;
	}
	
    improved_IFFM(tt3[0], tt3[1], tt3[2], tt3[3], ci[0], ci[1]); 
}