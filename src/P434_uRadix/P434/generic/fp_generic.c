/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: portable modular arithmetic for P434
*********************************************************************************************/

#include "../P434_internal.h"

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

void digit_a_digit_2( digit_t a1,  digit_t a0,  digit_t b1,  digit_t b0, digit_t* c)//60bit
{// 73+73  ->  2x60bits
 // input 60bits x 2
 // output 60bits x 2
	digit_t tmp, carry;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	tmp = a0 + b0;
	c[0] = tmp & mask_low;
	carry = tmp >> 60;
	c[1] = a1 + b1 + carry;
}

void fp_add_new(digit_t* a1, digit_t* a0, digit_t* b1, digit_t* b0, digit_t* c0, digit_t* c1)
{// Modular Addition based on the new data representation , with the unconditional radix
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], a[i] & b[i] are presented in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 60 bits , to compute c = a + b
	unsigned int i,j;
	digit_t ca = 0, ca1, t[2], d[6][2];
	digit_t R[2]  = {0x5EE84B000000000, 0x15EB};// 15EB 5EE84B000000000
	digit_t R_div_3[2] = {0x74F819000000000, 0x74E};//74E 74F819000000000
 
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
	t[0] = d[5][0] - R_div_3[0];
	ca1 = t[0] >> 63;
	t[1] = d[5][1] - R_div_3[1] - ca1;
	j = t[1] >> 63;
	ca = 1 ^ j;
	
	c0[5] = (t[0] + R_div_3[0]*j) & mask_low;
	c1[5] = t[1] + R_div_3[1] * j + ca1*j;
	ca1 = c0[0] + ca;
	c0[0] = ca1 & mask_low;
	ca = ca1 >> 60;
	c1[0] = c1[0] + ca;
}


void fpneg_new(digit_t* a1, digit_t* a0, digit_t* b0, digit_t* b1)
{// Modular Negation based on the new data representation , with the unconditional radix
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], a[i] all in 60 bits
 // Output : b[i] in [0, R-1] for 0≤i<n-1; b[n-1] in [0, R/3 - 1], b[i] all in 60 bits, to compute b = -a
	unsigned int i;
	digit_t R[2]  = {0x5EE84B000000000, 0x15EB};// 15EB 5EE84B000000000
	digit_t R_div_3[2] = {0x74F819000000000, 0x74E};//74E 74F819000000000
	digit_t k = 5, tmp, carry;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	for( i = 0; i < k; i++){
		tmp = R[0] - a0[i] - 1;
		b0[i] = tmp & mask_low;
		carry = tmp >> 63;
		b1[i] = R[1] - a1[i] - carry;
	}

	tmp = R_div_3[0] - a0[k] - 1;
	b0[k] = tmp & mask_low;
	carry = tmp >> 63;
	b1[k] = R_div_3[1] - a1[k] - carry;
}


void fpsub_new(digit_t* a1, digit_t* a0, digit_t* b1, digit_t* b0, digit_t* c0, digit_t* c1)
{// Modular Subtraction based on the new data representation , with the unconditional radix
 // Inputs : a[i], b[i] in [0, R-1] for 0≤i<n-1; a[n-1], b[n-1] in [0, R/3 - 1], a[i] & b[i] are presented in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], all in 60 bits , to compute c = a - b
	unsigned int i, j = 0, t;
	digit_t R[2]  = {0x5EE84B000000000, 0x15EB};// 15EB 5EE84B000000000
	digit_t R_div_3[2] = {0x74F819000000000, 0x74E};//74E 74F819000000000
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
	tmp = c0[5] + j * R_div_3[0] ;
	c0[5] = tmp & mask_low;
	borrow = tmp >> 60;
	c1[5] = c1[5] + j * R_div_3[1] + borrow;

	c0[0] = c0[0] - j;
}


void fpdiv2_434_new(digit_t* a1, digit_t* a0, digit_t* c0, digit_t* c1)
{// Modular Division by Two based on the new data representation , with the unconditional radix
 // Inputs : a[i] in [0, R-1] for 0≤i<n-1; a[n-1] in [0, R/3 - 1], a[i] all in 60 bits
 // Output : c[i] in [0, R-1] for 0≤i<n-1; c[n-1] in [0, R/3 - 1], c[i] all in 60 bits, to compute c = a/2
    unsigned int i;
    digit_t mask;
	digit_t R_div_2[2] = {0xAF7425800000000, 0xAF5};//AF5 AF7425800000000
	digit_t R_div_6[2] = {0x3A7C0C800000000, 0x3A7};
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
	b[5][0] = R_div_6[0] & ((digit_t)(0-mask));
	b[5][1] = R_div_6[1] & ((digit_t)(0-mask));
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

void improved_BR(digit_t c_2,  digit_t c_1,  digit_t c_0, digit_t* q, digit_t* r)
{
	digit_t b[3] = {0x7C1F1B, 0x9F06DE, 0x5D6EE269};//5D6EE269 9F06DE 7C1F1B
	digit_t m[3] = {0x5EE84B, 0x15EB, 0x5EFE36};//15EB 5EE84B
	digit_t n = 0x15EB5EE84B;
	
	digit_t mask_low = (digit_t)(-1) >> 40;  //low 24
	digit_t mask_low_1 = (digit_t)(-1) >> 26;  //low 38
	digit_t mask_low_2 = (digit_t)(-1) >> 28;  //low 36
	digit_t mask_low_3 = (digit_t)(-1) >> 4;   //low 60
	digit_t t_r[3], albl, ambm, ahbh, albm_ambl, albh_ahbl, ambh_ahbm, res1, res2, res3, carry, temp, t;
	
	t_r[0] = (c_1>>11) & mask_low;//24
	t_r[1] = (c_1>>35) & mask_low;//24  
	t_r[2] = (c_1>>59)^(c_2<<1);
	
	albl = t_r[0] * b[0];
	ambm = t_r[1] * b[1];
	ahbh = t_r[2] * b[2];
	albm_ambl = (t_r[0]+t_r[1])*(b[1]+b[0])-albl-ambm;
    albh_ahbl = (t_r[0]+t_r[2])*(b[2]+b[0])-albl-ahbh;
	ambh_ahbm = (t_r[1]+t_r[2])*(b[2]+b[1])-ambm-ahbh;	
	
    res2 = albh_ahbl + ambm;
	temp = (albl >> 24) + albm_ambl + ((res2 & mask_low) << 24);
	carry = temp >> 48;
   
	temp = (res2 >> 24) + ambh_ahbm + carry;
	temp = (temp >> 8) + ((ahbh & ((digit_t)(-1) >> 20))<<16);
	q[0] = temp & mask_low_3;

	q[1] = (temp>>60) + (ahbh >> 44);

	temp = q[0] & mask_low;
	carry = (q[0] & mask_low_1)>>24;
	
	albl = temp * m[0];
	ahbh = carry * m[1];
	albh_ahbl = (temp+carry)*m[2]-albl-ahbh;
	
	res1 = albl + (albh_ahbl<<24);
	
	t = (c_0 >> 36)^(c_1 << 24);
	r[0] = ((t & mask_low_1) - (res1 & mask_low_1)) & mask_low_1;
	
	t = r[0] - n;
	temp = t >> 63;
	carry = 1 - temp;
	r[0] = t + (n * temp);
	
	temp = q[0] + carry;
	q[0] = temp & mask_low_3;
	carry = temp >> 60;
	q[1] = q[1] + carry;
	
	t = c_0 & mask_low_2;
	
	r[1] = r[0] >> 24;
	r[0] = ((r[0] << 36) ^ t) & mask_low_3;
}

void last_improved_BR(digit_t c_1,  digit_t c_0, digit_t* q, digit_t* r)
{
	digit_t b[2] = {0x02699F06DE7C1F1B, 0x00005D6EE};
	digit_t m = 0x15EB5EE84B;
	digit_t t_r[1]; digit_t abl[2], abh[2], res1, res2, t, tmp, carry;
	digit_t mask_low = (digit_t)(-1) >> 26;  //low 38
	digit_t mask_low_1 = (digit_t)(-1) >> 28;//low 36
	digit_t mask_low_2 = (digit_t)(-1) >> 4; //low 60
	
	t_r[0] = c_1 >> 11;
	
	res1 = t_r[0] * b[1];
	
	q[0] = res1 >> 20; q[1] = 0;
	
	res1 = q[0] * m;
	t = (c_0 >> 36)^(c_1 << 24);
	r[0] = ((t & mask_low) - (res1 & mask_low)) & mask_low;
	
	t = r[0] - m;
	tmp = t >> 63;
	carry = 1 - tmp;
	r[0] = t + (m * tmp);
	
	q[0] = q[0] + carry;
	
	t = c_0 & mask_low_1;
	
	r[1] = r[0] >> 24;
	r[0] = ((r[0] << 36) ^ t) & mask_low_2;
}


void improved_IFFM(digit_t* c_0, digit_t* c_1, digit_t* c_2, digit_t* C0, digit_t* C1)
{
	digit_t q[2], r[2];
	digit_t tmp1, tmp2, tmp, c_r, carry; //tmp3;
	digit_t h, i, j, k;
	unsigned int m;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	digit_t R_div_3[2] = {0x74F819000000000, 0x74E};//74E 74F819000000000
	
	improved_BR(c_2[4], c_1[4], c_0[4], q, r);	
	C0[4] = r[0];
	C1[4] = r[1];
	
	tmp = c_0[5] + q[0];
	tmp1 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_1[5] + q[1] + carry;
	tmp2 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_2[5] + carry;
	
	improved_BR(tmp, tmp2, tmp1, q, r);
	C0[5] = r[0];
	C1[5] = r[1];
	
	tmp = C0[5] - R_div_3[0];
	j = tmp >> 63;
	c_r = C1[5] - R_div_3[1] - j;
	k = c_r >> 63;
	i = 1 ^ k;
	h = i;
	
	C0[5] = (tmp + R_div_3[0] * k) & mask_low;
	C1[5] = c_r + R_div_3[1] * k + j*k;
	
	tmp = C0[5] - R_div_3[0];
	j = tmp >> 63;
	c_r = C1[5] - R_div_3[1] - j;
	k = c_r >> 63;
	i = 1 ^ k;
	
	C0[5] = (tmp + R_div_3[0] * k) & mask_low;
	C1[5] = c_r + R_div_3[1] * k + j*k;
	
	h += i;//carry
	
	tmp = c_0[0] + (q[0] << 1) + q[0] + h;
	tmp1 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_1[0] + (q[1] << 1) + q[1] + carry;
	tmp2 = tmp & mask_low;
	carry = tmp >> 60;
	tmp = c_2[0] + carry;
	
	improved_BR(tmp, tmp2, tmp1, q, r);
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
	
		improved_BR(tmp, tmp2, tmp1, q, r);
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
	
	tmp = C0[5] - R_div_3[0];
	j = tmp >> 63;
	c_r = C1[5] - R_div_3[1] - j;
	k = c_r >> 63;
	i = 1 ^ k;
	
	C0[5] = (tmp + R_div_3[0] * k) & mask_low;
	C1[5] = c_r + R_div_3[1] * k + j*k;
	
	tmp = C0[0] + i;
	C0[0] = tmp & mask_low;
	carry = tmp >> 60;
	C1[0] = C1[0] + carry;
}


void digit_x_digit_2( digit_t a1,  digit_t a0,  digit_t b1,  digit_t b0, digit_t* c)//60bit
{// 128x128 max 73X73 -> 4x60bit 
 // input 60bits x 2
 // output 60bits x 3
	digit_t t_r[3], t1[3], albl, ambm, ahbh, albm_ambl, albh_ahbl, ambh_ahbm, temp;
	digit_t tmp, tmp2[2], tmp3[2], carry ;
	digit_t mask_low = (digit_t)(-1) >> 34; // 30bit
	digit_t mask_low_1 = (digit_t)(-1) >> 4; //60
	
	//for a
	t_r[0] = a0 & mask_low;//30
	t_r[1] = a0>>30;//30  
	t_r[2] = a1;
	//for b
	t1[0] = b0 & mask_low;//30
	t1[1] = b0>>30;//30  
	t1[2] = b1;
	
	albl = t_r[0] * t1[0];
	ambm = t_r[1] * t1[1];
	ahbh = a1 * b1;
	albm_ambl = (t_r[0]+t_r[1])*(t1[1]+t1[0])-albl-ambm;
    albh_ahbl = (t_r[0]+t_r[2])*(t1[2]+t1[0])-albl-ahbh;
	ambh_ahbm = (t_r[1]+t_r[2])*(t1[2]+t1[1])-ambm-ahbh;	
	
	temp = albl + ((albm_ambl & mask_low) <<30);
	c[0] = temp & mask_low_1;
	carry = (temp >> 60);
	
    temp = (albm_ambl>>30) + ambm + albh_ahbl + carry + ((ambh_ahbm & mask_low) <<30);
	c[1] = temp & mask_low_1;
	carry = (temp >> 60);

	c[2] = carry + (ambh_ahbm >> 30) + ahbh;
}

void digit_x_digit_2b( digit_t a1,  digit_t a0,  digit_t b1,  digit_t b0, digit_t* c)//60bit
{// 128x128 max 73X73 -> 4x60bit 
 // input 60bits x 2
 // output 60bits x 3
	digit_t albl[2], ahbh[2], ahbl[2], albh[2]; //ahbl_albh[2], d0, d1;
	digit_t tmp, tmp2[2], tmp3[2], carry ;
	digit_t mask_low = (digit_t)(-1) >> 4; // 4x60bit
	
	digit_x_digit_1(a0, b0, albl);
	digit_x_digit_1(a0, b1, albh);
	digit_x_digit_1(a1, b0, ahbl);
	digit_x_digit_1(a1, b1, ahbh);
	
	c[0] = albl[0];
	
	tmp = albl[1] + albh[0] + ahbl[0];//d0;
	c[1] = tmp & mask_low;
	carry = tmp >> 60;
	
	c[2] = ahbh[0] + carry + albh[1] + ahbl[1]; //d1; 
}

void digit_s_digit_4( digit_t* a,  digit_t* b, digit_t* c)
{// input 60x3
 // output 60x3
	digit_t tmp1, tmp2, tmp3, borrow;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	tmp1 = a[0] - b[0];
	borrow = tmp1 >> 63;
	tmp2 = a[1] - b[1] - borrow;
	borrow = tmp2 >> 63;
	tmp3 = a[2] - b[2] - borrow;
	c[0] = tmp1 & mask_low;
	c[1] = tmp2 & mask_low;
	c[2] = tmp3;// & mask_low;
}

void digit_a_digit_4( digit_t* a, digit_t* b, digit_t* c)
{// input 60x3
 // output 60x3
	digit_t tmp1, tmp2, tmp3, carry;
	digit_t mask_low = (digit_t)(-1) >> 4;
	
	tmp1 = a[0] + b[0];
	c[0] = tmp1 & mask_low;
	carry = tmp1 >> 60;
	tmp2 = a[1] + b[1] + carry;
	c[1] = tmp2 & mask_low;
	carry = tmp2 >> 60;
	c[2] = a[2] + b[2] + carry;
}

void mp_mul_new( digit_t* a1,  digit_t* a0,  digit_t* b1,  digit_t* b0, digit_t* c_0, digit_t* c_1, digit_t* c_2)
{// input 64bits  output 60bits x 3
	digit_t ci[6][3];
	digit_t d1[2], d2[2], e[3], f[3] = {0,0,0};
	digit_t mask_low = (digit_t)(-1) >> 4; 
	unsigned int carryOut, i, k, j, s; digit_t temp;

	for (i = 0; i < 6; i++) {
		digit_x_digit_2(a1[i], a0[i], b1[i], b0[i], ci[i]);//73x73  146bits -> 3x60bits
	} 
	c_0[0] = ci[0][0];
    c_1[0] = ci[0][1];
	c_2[0] = ci[0][2];
	
	for(i = 2; i < 6; i = i + 2){
		k = i>>1;
		for(j=0;j<k;j++){
			s = i-j;
			digit_a_digit_2(a1[j], a0[j], a1[s], a0[s], d1);//73+73 or 74+74  -> 60bitx2
			digit_a_digit_2(b1[j], b0[j], b1[s], b0[s], d2);
			digit_x_digit_2(d1[1], d1[0], d2[1], d2[0], e); // -> 3x60bits
			digit_s_digit_4(e, ci[j], e);
			digit_s_digit_4(e, ci[s], e);
			digit_a_digit_4(e, f, f);                       //one-level Karatsuba to reduce the complexity in the modular multiplication
		}
		digit_a_digit_4(ci[k],f,f);
		c_0[i] = f[0];
		c_1[i] = f[1];
		c_2[i] = f[2];
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
	}   //e.g.  c4 = a4*b0 + a3*b1 + a2*b2 + a1*b3 + a0*b4;

	for(i = 1; i < 6; i = i + 2){
		for(j=0;j<=((i-1)>>1);j++){
			s = i-j;
			digit_a_digit_2(a1[j], a0[j], a1[s], a0[s], d1);//73+73 or 74+74  -> 60bitx2
			digit_a_digit_2(b1[j], b0[j], b1[s], b0[s], d2);
			digit_x_digit_2(d1[1], d1[0], d2[1], d2[0], e); // -> 3x60bits
			digit_s_digit_4(e, ci[j], e);
			digit_s_digit_4(e, ci[s], e);
			digit_a_digit_4(e, f, f);                       //one-level Karatsuba to reduce the complexity in the modular multiplication
		}
		c_0[i] = f[0];
		c_1[i] = f[1];
		c_2[i] = f[2];
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
	}  

	for (i = 6; i < 10; i = i + 2) {
		k = i>>1;
		for (j = i-5; j < k; j++) {
			s = i-j;
			digit_a_digit_2(a1[j], a0[j], a1[s], a0[s], d1);//73+73 or 74+74  -> 60bitx2
			digit_a_digit_2(b1[j], b0[j], b1[s], b0[s], d2);
			digit_x_digit_2(d1[1], d1[0], d2[1], d2[0], e); // -> 3x60bits
			digit_s_digit_4(e, ci[j], e);
			digit_s_digit_4(e, ci[s], e);
			digit_a_digit_4(e, f, f);                       //one-level Karatsuba to reduce the complexity in the modular multiplication
        }
		digit_a_digit_4(ci[k],f,f);   
		s = i-6;		
		temp = (f[0]<<1) + f[0] + c_0[s];   
		c_0[s] = temp & mask_low;
		carryOut = temp >> 60;
		temp = (f[1]<<1) + f[1] + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> 60;
		c_2[s] = (f[2]<<1) + f[2] + c_2[s] + carryOut;
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
	}                                                       // get the raw coefficients.

	for (i = 7; i < 10; i = i + 2) {
		for (j = i-5; j <= ((i-1)>>1); j++) {
			s = i - j;
			digit_a_digit_2(a1[j], a0[j], a1[s], a0[s], d1);//73+73 or 74+74  -> 60bitx2
			digit_a_digit_2(b1[j], b0[j], b1[s], b0[s], d2);
			digit_x_digit_2(d1[1], d1[0], d2[1], d2[0], e); // -> 3x60bits
			digit_s_digit_4(e, ci[j], e);
			digit_s_digit_4(e, ci[s], e);
			digit_a_digit_4(e, f, f);                       //one-level Karatsuba to reduce the complexity in the modular multiplication
        }                                              
		s = i-6;
		temp = (f[0]<<1) + f[0] + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> 60;//
		temp = (f[1]<<1) + f[1] + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> 60;//
		c_2[s] = (f[2]<<1) + f[2] + c_2[s] + carryOut;
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
	}                                                    // get the raw coefficients.

	temp = (ci[5][0]<<1) + ci[5][0] + c_0[4];
	c_0[4] = temp & mask_low;
	carryOut = temp >> 60;
	temp = (ci[5][1]<<1) + ci[5][1] + c_1[4] + carryOut;
	c_1[4] = temp & mask_low;
	carryOut = temp >> 60;
	c_2[4] = (ci[5][2]<<1) + ci[5][2] + c_2[4] + carryOut;  // get the final raw coefficients.
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

void digit_sqr_2b( digit_t a1,  digit_t a0, digit_t* c)
{
 // input 60bits x 2
 // output 60bits x 3
	digit_t alal[2], ahah, alah[2];
	digit_t tmp, carry ;
	digit_t mask_low = (digit_t)(-1) >> 4; // 4x60bit
	
	digit_sqr(a0, alal);
	digit_x_digit_1(a0, a1, alah);
//	digit_sqr(a1, ahah);
	ahah = a1 * a1;
	
	c[0] = alal[0];
	
	tmp = alal[1] + (alah[0] << 1);
	c[1] = tmp & mask_low;
	carry = tmp >> 60;
	
	c[2] = (alah[1] << 1) + ahah + carry;
}

void digit_sqr_2( digit_t a1,  digit_t a0, digit_t* c)
{
 // input 60bits x 2
 // output 60bits x 3
	digit_t t_r[2], albl, ambm, ahbh, albm_ambl, albh_ahbl, ambh_ahbm, temp;
	digit_t tmp, tmp2[2], tmp3[2], carry ;
	digit_t mask_low = (digit_t)(-1) >> 34; // 30bit
	digit_t mask_low_1 = (digit_t)(-1) >> 4; //60
	
	//for a
	t_r[0] = a0 & mask_low;//30
	t_r[1] = a0>>30;//30  
//	t_r[2] = a1;
	
	albl = t_r[0] * t_r[0];
	ambm = t_r[1] * t_r[1];
	ahbh = a1 * a1;
	albm_ambl = (t_r[0]+t_r[1])*(t_r[1]+t_r[0])-albl-ambm;
    albh_ahbl = (t_r[0]+a1)*(a1+t_r[0])-albl-ahbh;
	ambh_ahbm = (t_r[1]+a1)*(a1+t_r[1])-ambm-ahbh;	
	
	temp = albl + ((albm_ambl & mask_low) <<30);
	c[0] = temp & mask_low_1;
	carry = (temp >> 60);
	
    temp = (albm_ambl>>30) + ambm + albh_ahbl + carry + ((ambh_ahbm & mask_low) <<30);
	c[1] = temp & mask_low_1;
	carry = (temp >> 60);

	c[2] = carry + (ambh_ahbm >> 30) + ahbh;
}

void mp_sqr_new(digit_t* a1,  digit_t* a0, digit_t* c_0, digit_t* c_1, digit_t* c_2)
{// input 60bits  output 60bits x 3
	digit_t ci[6][3];
	digit_t f[3] = {0,0,0}, e[3] = {0,0,0};
	digit_t mask_low = (digit_t)(-1) >> 4; 
	unsigned int carryOut, i, k, j, s; digit_t temp;

	for (i = 0; i < 6; i++) {
		digit_sqr_2(a1[i], a0[i], ci[i]);//73x73  146bits -> 3x60bits
	} 
	c_0[0] = ci[0][0];
    c_1[0] = ci[0][1];
	c_2[0] = ci[0][2];
	
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
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
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
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
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
		s = i-6;		
		temp = (f[0]<<1) + f[0] + c_0[s];   
		c_0[s] = temp & mask_low;
		carryOut = temp >> 60;
		temp = (f[1]<<1) + f[1] + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> 60;
		c_2[s] = (f[2]<<1) + f[2] + c_2[s] + carryOut;
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
	}                                                       // get the raw coefficients.

	for (i = 7; i < 10; i = i + 2) {
		for (j = i-5; j <= ((i-1)>>1); j++) {
			s = i - j;
			digit_x_digit_2(a1[j], a0[j], a1[s], a0[s], e);
			digit_a_digit_4(e, f, f);
        }                                                 
		digit_a_digit_4(f, f, e);
		s = i-6;
		temp = (e[0]<<1) + e[0] + c_0[s];
		c_0[s] = temp & mask_low;
		carryOut = temp >> 60;//
		temp = (e[1]<<1) + e[1] + c_1[s] + carryOut;
		c_1[s] = temp & mask_low;
		carryOut = temp >> 60;//
		c_2[s] = (e[2]<<1) + e[2] + c_2[s] + carryOut;
		f[0] = 0;
		f[1] = 0;
		f[2] = 0;
	}                                                    // get the raw coefficients.
	temp = (ci[5][0]<<1) + ci[5][0] + c_0[4];
	c_0[4] = temp & mask_low;
	carryOut = temp >> 60;
	temp = (ci[5][1]<<1) + ci[5][1] + c_1[4] + carryOut;
	c_1[4] = temp & mask_low;
	carryOut = temp >> 60;
	c_2[4] = (ci[5][2]<<1) + ci[5][2] + c_2[4] + carryOut;     // get the final raw coefficients.
}

void fpmul( digit_t* a1,  digit_t* a0,  digit_t* b1,  digit_t* b0, digit_t* d0, digit_t* d1)
{
	digit_t c0[6], c1[6], c2[6];
	mp_mul_new(a1, a0, b1, b0, c0, c1, c2);
	improved_IFFM(c0, c1, c2, d0, d1);
}

void fpsqr( digit_t* a1,  digit_t* a0, digit_t* d0, digit_t* d1)
{
	digit_t c0[6], c1[6], c2[6];
	mp_sqr_new(a1,a0,c0, c1, c2);
	improved_IFFM(c0, c1, c2, d0, d1);
}


void fpinv_chain_new_u(digit_t* a1, digit_t* a0, digit_t* c0, digit_t* c1)
{ // Chain to compute a^(p-3)/4 using Montgomery arithmetic.
    unsigned int i, j;
	
    digit_t t[31][2][6], tt[2][6];

    // Precomputed table
    fpsqr(a1, a0, tt[0], tt[1]);
    fpmul(a1, a0, tt[1], tt[0], t[0][0], t[0][1]);
    for (i = 0; i <= 29; i++) fpmul(t[i][1], t[i][0], tt[1], tt[0], t[i+1][0], t[i+1][1]);

//    fpcopy(a, tt);

	for (i=0;i<6;i++){
		tt[0][i] = a0[i];
		tt[1][i] = a1[i];
	}

    for (i = 0; i < 7; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[5][1], t[5][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 10; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[14][1], t[14][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[3][1], t[3][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[23][1], t[23][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[13][1], t[13][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[24][1], t[24][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[7][1], t[7][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[12][1], t[12][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[30][1], t[30][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[1][1], t[1][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[30][1], t[30][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[21][1], t[21][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[2][1], t[2][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[19][1], t[19][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[1][1], t[1][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[24][1], t[24][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[26][1], t[26][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[16][1], t[16][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[10][1], t[10][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[6][1], t[6][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[0][1], t[0][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[20][1], t[20][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 8; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[9][1], t[9][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[25][1], t[25][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[30][1], t[30][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[26][1], t[26][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(a1, a0, tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 7; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[28][1], t[28][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[6][1], t[6][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[10][1], t[10][0], tt[1], tt[0], tt[0], tt[1]);
    for (i = 0; i < 9; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
    fpmul(t[22][1], t[22][0], tt[1], tt[0], tt[0], tt[1]);
    for (j = 0; j < 35; j++) {
        for (i = 0; i < 6; i++) fpsqr(tt[1], tt[0], tt[0], tt[1]);
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

void fp2mul_uRadix(f2elm_t a, f2elm_t ai, f2elm_t b, f2elm_t bi, f2elm_t c, f2elm_t ci)
{// GF(p^2) multiplication, c = a*b in GF(p^2).
 // Inputs : a = a + ai*i, b = b + bi*i. a, ai, b, bi all meet the same requirement with fpmul_new & fpsqr_new
 // Output : c = c + ci*i. c = a * b, c and ci meet the same requirement with fpmul_new & fpsqr_new
 // a-b+2^4*R*P mod P = (a_5-b_5+2^4*(R/3-1)*R)*R^5 + (a_4-b_4+2^4*(R-1)*R)*R^4 + ... + a_0-b_0+2^4*(R-1)*R
    digit_t PR_5[3] = {0x0117B50000000000, 0x09FD09495751A14A, 0x0A027532}; // 2^4*(R/3-1)*R) 0A027532 9FD09495751A14A 117B50000000000
    digit_t PR[3] =   {0x0117B50000000000, 0x0DF71BDC05F7A14A, 0x1E075F97}; // 2^4*(R-1)*R 1E075F97 DF71BDC05F7A14A 117B50000000000
    digit_t tt1[3][6], tt2[3][6], tt3[3][6]; digit_t t0[6], t1[6], h0[6], h1[6];
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
    mp_mul_new(a[1], a[0], b[1], b[0], tt1[0], tt1[1], tt1[2]); // tt1 = a*b
	mp_mul_new(ai[1], ai[0], bi[1], bi[0], tt2[0], tt2[1], tt2[2]); // tt2 = ai*bi

	for (i=0; i<5; i++){
		//a*b + 2^4*R*P - ai*bi
		temp = tt1[0][i] + PR[0];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> k;
		temp = tt1[1][i] + PR[1] + res;
		tt3[1][i] = temp & mask_low_1;
		res = temp >> k;
		tt3[2][i] = tt1[2][i] + PR[2] + res;
		
		temp = tt3[0][i] - tt2[0][i];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> 63;
		temp = tt3[1][i] - tt2[1][i] - res;
		tt3[1][i] = temp & mask_low_1;
		res = temp >> 63;
		tt3[2][i] = tt3[2][i] - tt2[2][i] - res;
	}
	
	//a*b + 2^4*R*P - ai*bi
	temp = tt1[0][5] + PR_5[0];
	tt3[0][5] = temp & mask_low_1;
	res = temp >> k;
	temp = tt1[1][5] + PR_5[1] + res;
	tt3[1][5] = temp & mask_low_1;
	res = temp >> k;
	tt3[2][5] = tt1[2][5] + PR_5[2] + res;
		
	temp = tt3[0][5] - tt2[0][5];
	tt3[0][5] = temp & mask_low_1;
	res = temp >> 63;
	temp = tt3[1][5] - tt2[1][5] - res;
	tt3[1][5] = temp & mask_low_1;
	res = temp >> 63;
	tt3[2][5] = tt3[2][5] - tt2[2][5] - res;
	
	improved_IFFM(tt3[0], tt3[1], tt3[2], c[0], c[1]);// c[0]
	
	//for c[1]
	
	mp_mul_new(t1, t0, h1, h0, tt3[0], tt3[1], tt3[2]); //   (a+ai)*(b+bi)	
	
	// c[1] = (a+ai)*(b+bi) - (a*b + ai*bi) 
	for (i=0; i<6; i++) {
		//a*b + ai*bi
		temp = tt2[0][i] + tt1[0][i];
		tt2[0][i] = temp & mask_low_1;
		res = temp >> k;
		temp = tt2[1][i] + tt1[1][i] + res;
		tt2[1][i] = temp & mask_low_1;
		res = temp >> k;
		tt2[2][i] = tt2[2][i] + tt1[2][i] + res;

		//(a+ai)*(b+bi) - a*b - ai*bi
		temp = tt3[0][i] - tt2[0][i];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> 63;		
		temp = tt3[1][i] - tt2[1][i] - res;
		tt3[1][i] = temp & mask_low_1;
		res = temp >> 63;
		tt3[2][i] = tt3[2][i] - tt2[2][i] - res;
	}
	
    improved_IFFM(tt3[0], tt3[1], tt3[2], ci[0], ci[1]); 
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

void digit_x_digit_4( digit_t a1,  digit_t a0,  digit_t b1,  digit_t b0, digit_t* c)//60bit
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

void from_uRadix(digit_t* a1, digit_t* a0, digit_t* co)
{
	digit_t R[2]  = {0x5EE84B000000000, 0x015EB};//0x015EB 5EE84B000000000
	digit_t R_sub_1[2] = {0x5EE84AFFFFFFFFF, 0x015EB};//015EB 5EE84AFFFFFFFFF
	digit_t R_div_3_sub1[2] = {0x74F818FFFFFFFFF, 0x74E};//74E 74F818FFFFFFFFF
	
	digit_t b[3] = {0}, d[4], e[3], f[3], g[3], h[3], zero[1] = {0}, c[11], carry[1], tmp, m[4], n[4];
	digit_t mask_low_1 = (digit_t)(-1) >> 4; int i, j, k;
	
	j = 1;
	for(i=1;i<5;i++){
		k = (a1[i] == R_sub_1[1])&(a0[i]==R_sub_1[0]);
		j = j & k;
	}
	
	k = (a1[5] == R_div_3_sub1[1])&(a0[5]==R_div_3_sub1[0]);
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
	c[0] = f[0]; c[1] = f[1]; c[2] = b[2] + f[2];
	
	digit_x_digit_4(c[1], c[0], R[1], R[0], d);
	digit_x_digit_2(zero[0], c[2], R[1], R[0], g);
	digit_a_digit_1(d[1], d[0], a1[3], a0[3], e);
	c[0] = e[0]; c[1] = e[1];
	tmp = d[2] + g[0] + e[2]; c[2] = tmp & mask_low_1;
	c[3] = d[3] + g[1] + (tmp >> 60);
	
	digit_x_digit_4(c[1], c[0], R[1], R[0], d);
	digit_x_digit_2b(c[3], c[2], R[1], R[0], h); //
	digit_a_digit_1(d[1], d[0], a1[2], a0[2], e);
	c[0] = e[0]; c[1] = e[1];
	tmp = d[2] + h[0] + e[2]; c[2] = tmp & mask_low_1;
	tmp = d[3] + h[1] + (tmp >> 60); c[3] = tmp & mask_low_1;
	c[4] = h[2] + (tmp >> 60);
	
	digit_x_digit_4(c[1], c[0], R[1], R[0], d);
	digit_x_digit_4(c[3], c[2], R[1], R[0], m);
	digit_x_digit_2(zero[0], c[4], R[1], R[0], g);
	digit_a_digit_1(d[1], d[0], a1[1], a0[1], e);
	c[0] = e[0]; c[1] = e[1];
	tmp = d[2] + m[0] + e[2]; c[2] = tmp & mask_low_1;
	tmp = d[3] + m[1] + (tmp >> 60); c[3] = tmp & mask_low_1;
	tmp = m[2] + g[0] + (tmp >> 60); c[4] = tmp & mask_low_1;
	tmp = m[3] + g[1] + (tmp >> 60); c[5] = tmp & mask_low_1;
	c[6] = g[2] + (tmp >> 60);
	
	digit_x_digit_4(c[1], c[0], R[1], R[0], d);
	digit_x_digit_4(c[3], c[2], R[1], R[0], m);
	digit_x_digit_4(c[5], c[4], R[1], R[0], n);
	digit_x_digit_2(zero[0], c[6], R[1], R[0], g);
	digit_a_digit_1(d[1], d[0], a1[0], a0[0], e);
	c[0] = e[0]; c[1] = e[1];
	tmp = d[2] + m[0] + e[2]; c[2] = tmp & mask_low_1;
	tmp = d[3] + m[1] + (tmp >> 60); c[3] = tmp & mask_low_1;
	tmp = m[2] + n[0] + (tmp >> 60); c[4] = tmp & mask_low_1;
	tmp = m[3] + n[1] + (tmp >> 60); c[5] = tmp & mask_low_1;
	tmp = n[2] + g[0] + (tmp >> 60); c[6] = tmp & mask_low_1;
	c[7] = n[3] + g[1] + (tmp >> 60);
	
	for (i=0;i<7;i++){
		co[i] = (c[i] >> (4*i))^(c[i+1] << (60-4*i));
	}
}

void digit_x_digit(const digit_t a, const digit_t b, digit_t* c)
{ // Digit multiplication, digit * digit -> 2-digit result    
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4);

    al = a & mask_low;                        // Low part
    ah = a >> (sizeof(digit_t) * 4);          // High part
    bl = b & mask_low;
    bh = b >> (sizeof(digit_t) * 4);

    albl = al*bl;
    albh = al*bh;
    ahbl = ah*bl;
    ahbh = ah*bh;
    c[0] = albl & mask_low;                   // C00

    res1 = albl >> (sizeof(digit_t) * 4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;  
    temp = res1 + res2 + res3;
    carry = temp >> (sizeof(digit_t) * 4);
    c[0] ^= temp << (sizeof(digit_t) * 4);    // C01   

    res1 = ahbl >> (sizeof(digit_t) * 4);
    res2 = albh >> (sizeof(digit_t) * 4);
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
	digit_t m = 0x15EB5EE84B;
	digit_t t[6] = {0}, temp[12], res1, res2, tmp[2], res3, s, t_r;
	digit_t mask_low = (digit_t)(-1) >> 26;  //low 38
	digit_t mask_low_2 = (digit_t)(-1) >> 4; //low 60
	unsigned int i, j;
	
	s = a[0] & ((digit_t)(-1) >> 28);
	t_r = (a[0] >> 36)^(a[1] << 28);
	
	for (i=0;i<unwords;i++){
		t[i] = (a[i+1] >> 7)^(a[i+2] << 57);
	}
	t[unwords] = a[unwords+1] >> 7;
	
	mp_mul_old(t, lamda, temp, unwords+1);
	
	for (i=0;i<unwords+1;i++){
		res1 = temp[unwords+i] >> (4+unwords*8);
		res2 = temp[unwords+i+1] << (60-unwords*8);
		q[i] = res1 ^ res2;
	}
	for (i=unwords+1;i<7;i++){
		q[i] = 0;
	}
	
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
	
	for (i = 1; i < unwords+1; i++){	
		res1 = q[i] + j;	
		j = (q[i]==(digit_t)(-1));
		q[i] = res1 * (1 ^ j);	
	}
	q[unwords+1] = q[unwords+1] + j;
	
	res1 = r[0] << 36;
	r[1] = r[0] >> 24;
	r[0] = (res1 ^ s) & mask_low_2;
}

void to_uRadix(digit_t* a, digit_t* c0, digit_t* c1){
	int i, j;
	digit_t lamda[6][6]={{0xb2ebda3a05779aa0, 0x810f5e895af29b22, 0xd9913b41f7cfc4b9, 0xe36b629cb6b5b916, 0x6de7c1f1b303365f, 0x5d6ee2699f0},{0xb9810f5e895af29b, 0x16d9913b41f7cfc4, 0x5fe36b629cb6b5b9, 0xf06de7c1f1b30336, 0x5d6ee2699},{0xb916d9913b41f7cf, 0x365fe36b629cb6b5, 0x99f06de7c1f1b303, 0x5d6ee26},{0x03365fe36b629cb6, 0x2699f06de7c1f1b3, 0x5d6ee},{0xee2699f06de7c1f1, 0x5d6},{0x5}};
	
	digit_t q[7], temp[2], tmp;
	digit_t mask_low = (digit_t)(-1) >> 4; //low 60
	
	for (i=0; i<5; i++){
		IBR_uRadix(a, q, temp, 5-i, lamda[i]);
		c0[i] = temp[0];
		c1[i] = temp[1];
		for(j=0;j<6-i;j++){
			a[j] = q[j];
		}
		for(j=6-i;j<7;j++){
			a[j] = 0;
		}
	}
	tmp = q[0] >> 60;
	c1[5] = (q[1] << 4)^tmp;
	c0[5] = q[0] & mask_low;
}

void fp2sqr_uRadix(f2elm_t a, f2elm_t ai, f2elm_t c, f2elm_t ci)
{// GF(p^2) multiplication, c = a*b in GF(p^2).
 // Inputs : a = a + ai*i, b = b + bi*i. a, ai, b, bi all meet the same requirement with fpmul_new & fpsqr_new
 // Output : c = c + ci*i. c = a * b, c and ci meet the same requirement with fpmul_new & fpsqr_new
 // a-b+2^4*R*P mod P = (a_5-b_5+2^4*(R/3-1)*R)*R^5 + (a_4-b_4+2^4*(R-1)*R)*R^4 + ... + a_0-b_0+2^4*(R-1)*R
    digit_t PR_5[3] = {0x0117B50000000000, 0x09FD09495751A14A, 0x0A027532}; // 2^4*(R/3-1)*R) 0A027532 9FD09495751A14A 117B50000000000
    digit_t PR[3] =   {0x0117B50000000000, 0x0DF71BDC05F7A14A, 0x1E075F97}; // 2^4*(R-1)*R 1E075F97 DF71BDC05F7A14A 117B50000000000
    digit_t tt1[3][6], tt2[3][6], tt3[3][6]; digit_t t0[6], t1[6], h0[6], h1[6];
	int i, j;
	digit_t temp, res;
	digit_t mask_low_1 = (digit_t)(-1) >> 4;
	digit_t k = 60;
	
	for (i=0; i<6; i++){
		temp = a[0][i] + ai[0][i];
		t0[i] = temp & mask_low_1;
		res = temp >> k;
		t1[i] = a[1][i] + ai[1][i] + res; //a+ai
	}
	
	//for c[0]
    mp_sqr_new(a[1], a[0], tt1[0], tt1[1], tt1[2]); // tt1 = a*b
	mp_sqr_new(ai[1], ai[0], tt2[0], tt2[1], tt2[2]); // tt2 = ai*bi

	for (i=0; i<5; i++){
		//a*b + 2^4*R*P - ai*bi
		temp = tt1[0][i] + PR[0];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> k;
		temp = tt1[1][i] + PR[1] + res;
		tt3[1][i] = temp & mask_low_1;
		res = temp >> k;
		tt3[2][i] = tt1[2][i] + PR[2] + res;
		
		temp = tt3[0][i] - tt2[0][i];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> 63;
		temp = tt3[1][i] - tt2[1][i] - res;
		tt3[1][i] = temp & mask_low_1;
		res = temp >> 63;
		tt3[2][i] = tt3[2][i] - tt2[2][i] - res;
	}
	
	//a*b + 2^4*R*P - ai*bi
	temp = tt1[0][5] + PR_5[0];
	tt3[0][5] = temp & mask_low_1;
	res = temp >> k;
	temp = tt1[1][5] + PR_5[1] + res;
	tt3[1][5] = temp & mask_low_1;
	res = temp >> k;
	tt3[2][5] = tt1[2][5] + PR_5[2] + res;
		
	temp = tt3[0][5] - tt2[0][5];
	tt3[0][5] = temp & mask_low_1;
	res = temp >> 63;
	temp = tt3[1][5] - tt2[1][5] - res;
	tt3[1][5] = temp & mask_low_1;
	res = temp >> 63;
	tt3[2][5] = tt3[2][5] - tt2[2][5] - res;
	
	improved_IFFM(tt3[0], tt3[1], tt3[2], c[0], c[1]);// c[0]
	
	//for c[1]
	
	mp_sqr_new(t1, t0, tt3[0], tt3[1], tt3[2]); //   (a+ai)*(b+bi)	
	
	// c[1] = (a+ai)*(b+bi) - (a*b + ai*bi) 
	for (i=0; i<6; i++) {
		//a*b + ai*bi
		temp = tt2[0][i] + tt1[0][i];
		tt2[0][i] = temp & mask_low_1;
		res = temp >> k;
		temp = tt2[1][i] + tt1[1][i] + res;
		tt2[1][i] = temp & mask_low_1;
		res = temp >> k;
		tt2[2][i] = tt2[2][i] + tt1[2][i] + res;

		//(a+ai)*(b+bi) - a*b - ai*bi
		temp = tt3[0][i] - tt2[0][i];
		tt3[0][i] = temp & mask_low_1;
		res = temp >> 63;		
		temp = tt3[1][i] - tt2[1][i] - res;
		tt3[1][i] = temp & mask_low_1;
		res = temp >> 63;
		tt3[2][i] = tt3[2][i] - tt2[2][i] - res;
	}
	
    improved_IFFM(tt3[0], tt3[1], tt3[2], ci[0], ci[1]); 
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