/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: ephemeral supersingular isogeny Diffie-Hellman key exchange (SIDH)
*********************************************************************************************/ 

#include "../random/random.h"


static void clear_words(void* mem, digit_t nwords)
{ // Clear digits from memory. "nwords" indicates the number of digits to be zeroed.
  // This function uses the volatile type qualifier to inform the compiler not to optimize out the memory clearing.
    unsigned int i;
    volatile digit_t *v = mem; 

    for (i = 0; i < nwords; i++) {
        v[i] = 0;
    }
}


static void init_basis(digit_t *gen, f2elm_t XP0, f2elm_t XP1, f2elm_t XQ0, f2elm_t XQ1, f2elm_t XR0, f2elm_t XR1)
{ // Initialization of basis points
    
    fpcopy(gen,                  XP0[0]);
    fpcopy(gen +   NWORDS_FIELD, XP0[1]);
    fpcopy(gen + 2*NWORDS_FIELD, XP1[0]);
    fpcopy(gen + 3*NWORDS_FIELD, XP1[1]);
    fpcopy(gen + 4*NWORDS_FIELD, XQ0[0]);
    fpcopy(gen + 5*NWORDS_FIELD, XQ0[1]);
	fpcopy(gen + 6*NWORDS_FIELD, XQ1[0]);
    fpcopy(gen + 7*NWORDS_FIELD, XQ1[1]);
	fpcopy(gen + 8*NWORDS_FIELD, XR0[0]);
    fpcopy(gen + 9*NWORDS_FIELD, XR0[1]);
	fpcopy(gen +10*NWORDS_FIELD, XR1[0]);
    fpcopy(gen +11*NWORDS_FIELD, XR1[1]);
}


static void fp2_encode(const f2elm_t x, const f2elm_t xi, unsigned char *enc)
{ // Conversion of GF(p^2) element from Montgomery to standard representation, and encoding by removing leading 0 bytes
    unsigned int i;
    normal2_t t;

    from_fp2uRadix(x, xi, t);
    for (i = 0; i < FP2_ENCODED_BYTES / 2; i++) {
        enc[i] = ((unsigned char*)t)[i];
        enc[i + FP2_ENCODED_BYTES / 2] = ((unsigned char*)t)[i + MAXBITS_FIELD / 8];
    }
}


static void fp2_decode(const unsigned char *enc, f2elm_t y, f2elm_t yi)
{ // Parse byte sequence back into GF(p^2) element, and conversion to Montgomery representation
    unsigned int i; normal2_t x;

    for (i = 0; i < 2*(MAXBITS_FIELD / 8); i++) ((unsigned char *)x)[i] = 0;
    for (i = 0; i < FP2_ENCODED_BYTES / 2; i++) {
        ((unsigned char*)x)[i] = enc[i];
        ((unsigned char*)x)[i + MAXBITS_FIELD / 8] = enc[i + FP2_ENCODED_BYTES / 2];
    }
    to_fp2uRadix(x, y, yi);
}


void random_mod_order_A(unsigned char* random_digits)
{  // Generation of Alice's secret key  
   // Outputs random value in [0, 2^eA - 1]
    unsigned long long nbytes = NBITS_TO_NBYTES(OALICE_BITS);

    clear_words((void*)random_digits, MAXWORDS_ORDER);
    randombytes(random_digits, nbytes);
    random_digits[nbytes-1] &= MASK_ALICE;    // Masking last byte 
}


void random_mod_order_B(unsigned char* random_digits)
{  // Generation of Bob's secret key  
   // Outputs random value in [0, 2^Floor(Log(2, oB)) - 1]
    unsigned long long nbytes = NBITS_TO_NBYTES(OBOB_BITS-1);

    clear_words((void*)random_digits, MAXWORDS_ORDER);
    randombytes(random_digits, nbytes);
    random_digits[nbytes-1] &= MASK_BOB;     // Masking last byte 
}


int EphemeralKeyGeneration_A(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA)
{ // Alice's ephemeral public key generation
  // Input:  a private key PrivateKeyA in the range [0, 2^eA - 1]. 
  // Output: the public key PublicKeyA consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_ALICE];
    f2elm_t XPA0, XPA1, XQA0, XQA1, XRA0, XRA1, coeff[6], A24plus = {0}, C24 = {0}, A = {0}, A24plusi = {0}, C24i = {0}, Ai = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;

    // Initialize basis points
    init_basis((digit_t*)A_gen, XPA0, XPA1, XQA0, XQA1, XRA0, XRA1);
    init_basis((digit_t*)B_gen, phiP->X, phiP->Xi, phiQ->X, phiQ->Xi, phiR->X, phiR->Xi);
    fpcopy((digit_t*)&uRadix_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&uRadix_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&uRadix_one, (phiR->Z)[0]);

    // Initialize constants: A24plus = A+2C, C24 = 4C, where A=6, C=1
    fpcopy((digit_t*)&uRadix_one, A24plus[0]);
    fp2add_new(A24plus, A24plusi, A24plus, A24plusi, A24plus, A24plusi);
    fp2add_new(A24plus, A24plusi, A24plus, A24plusi, C24, C24i);
    fp2add_new(A24plus, A24plusi, C24, C24i, A, Ai);
    fp2add_new(C24, C24i, C24, C24i, A24plus, A24plusi);

    // Retrieve kernel point
    LADDER3PT(XPA0, XPA1, XQA0, XQA1, XRA0, XRA1, (digit_t*)PrivateKeyA, ALICE, R, A, Ai);  

#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

    xDBLe(R, S, A24plus, A24plusi, C24, C24i, (int)(OALICE_BITS-1));
    get_2_isog(S, A24plus, A24plusi, C24, C24i); 
    eval_2_isog(phiP, S); 
    eval_2_isog(phiQ, S); 
    eval_2_isog(phiR, S);
    eval_2_isog(R, S);
#endif

    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
			fp2copy(R->Xi, pts[npts]->Xi);
            fp2copy(R->Zi, pts[npts]->Zi);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(R, R, A24plus, A24plusi, C24, C24i, (int)(2*m));
            index += m;
        }
        get_4_isog(R, A24plus, A24plusi, C24, C24i, coeff);       

        for (i = 0; i < npts; i++) {
            eval_4_isog(pts[i], coeff);
        }
        eval_4_isog(phiP, coeff);
        eval_4_isog(phiQ, coeff);
        eval_4_isog(phiR, coeff);

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
		fp2copy(pts[npts-1]->Xi, R->Xi); 
        fp2copy(pts[npts-1]->Zi, R->Zi);
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(R, A24plus, A24plusi, C24, C24i, coeff); 
    eval_4_isog(phiP, coeff);
    eval_4_isog(phiQ, coeff);
    eval_4_isog(phiR, coeff);

    inv_3_way(phiP->Z, phiP->Zi, phiQ->Z, phiQ->Zi, phiR->Z, phiR->Zi);
    fp2mul_uRadix(phiP->X, phiP->Xi, phiP->Z, phiP->Zi, phiP->X, phiP->Xi);
    fp2mul_uRadix(phiQ->X, phiQ->Xi, phiQ->Z, phiQ->Zi, phiQ->X, phiQ->Xi);
    fp2mul_uRadix(phiR->X, phiR->Xi, phiR->Z, phiR->Zi, phiR->X, phiR->Xi);
                
    // Format public key                   
    fp2_encode(phiP->X, phiP->Xi, PublicKeyA);
    fp2_encode(phiQ->X, phiQ->Xi, PublicKeyA + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, phiR->Xi, PublicKeyA + 2*FP2_ENCODED_BYTES);

    return 0;
}


int EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB)
{ // Bob's ephemeral public key generation
  // Input:  a private key PrivateKeyB in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  // Output: the public key PublicKeyB consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, phiP = {0}, phiQ = {0}, phiR = {0}, pts[MAX_INT_POINTS_BOB];
    f2elm_t XPB0, XPB1, XQB0, XQB1, XRB0, XRB1, coeff[6], A24plus = {0}, A24plusi = {0}, A24minus = {0}, A24minusi = {0}, A = {0}, Ai = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;

    // Initialize basis points
    init_basis((digit_t*)B_gen, XPB0, XPB1, XQB0, XQB1, XRB0, XRB1);
    init_basis((digit_t*)A_gen, phiP->X, phiP->Xi, phiQ->X, phiQ->Xi, phiR->X, phiR->Xi);
    fpcopy((digit_t*)&uRadix_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&uRadix_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&uRadix_one, (phiR->Z)[0]);

    // Initialize constants: A24minus = A-2C, A24plus = A+2C, where A=6, C=1
    fpcopy((digit_t*)&uRadix_one, A24plus[0]);
    fp2add_new(A24plus, A24plusi, A24plus, A24plusi, A24plus, A24plusi);
    fp2add_new(A24plus, A24plusi, A24plus, A24plusi, A24minus, A24minusi);
    fp2add_new(A24plus, A24plusi, A24minus, A24minusi, A, Ai);
    fp2add_new(A24minus, A24minusi, A24minus, A24minusi, A24plus, A24plusi);
	
    // Retrieve kernel point
    LADDER3PT(XPB0, XPB1, XQB0, XQB1, XRB0, XRB1, (digit_t*)PrivateKeyB, BOB, R, A, Ai);
    
//	normal2_t c;
//	from_fp2uRadix(R->X, R->Xi, c);
//	printf("\n");
//	for(i=0;i<10;i++){
//		printf("R->X0 c[%d] = %#018llx \n;",i,c[0][i]);
//		printf("R->X1 c[%d] = %#018llx \n;",i,c[1][i]);
//	};
//	printf("\n");
//	from_fp2uRadix(R->Z, R->Zi, c);
//	printf("\n");
//	for(i=0;i<10;i++){
//		printf("R->Z0 c[%d] = %#018llx \n;",i,c[0][i]);
//		printf("R->Z1 c[%d] = %#018llx \n;",i,c[1][i]);
//	};
//	printf("\n");
	
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
			fp2copy(R->Xi, pts[npts]->Xi);
            fp2copy(R->Zi, pts[npts]->Zi);
            pts_index[npts++] = index;
            m = strat_Bob[ii++];
            xTPLe(R, R, A24minus, A24minusi, A24plus, A24plusi, (int)m);
            index += m;
        } 
        get_3_isog(R, A24minus, A24minusi, A24plus, A24plusi, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog(pts[i], coeff);
        }     
        eval_3_isog(phiP, coeff);
        eval_3_isog(phiQ, coeff);
        eval_3_isog(phiR, coeff);

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
		fp2copy(pts[npts-1]->Xi, R->Xi); 
        fp2copy(pts[npts-1]->Zi, R->Zi);
        index = pts_index[npts-1];
        npts -= 1;
    }
    
    get_3_isog(R, A24minus, A24minusi, A24plus, A24plusi, coeff);
    eval_3_isog(phiP, coeff);
    eval_3_isog(phiQ, coeff);
    eval_3_isog(phiR, coeff);

    inv_3_way(phiP->Z, phiP->Zi, phiQ->Z, phiQ->Zi, phiR->Z, phiR->Zi);
    fp2mul_uRadix(phiP->X, phiP->Xi, phiP->Z, phiP->Zi, phiP->X, phiP->Xi);
    fp2mul_uRadix(phiQ->X, phiQ->Xi, phiQ->Z, phiQ->Zi, phiQ->X, phiQ->Xi);
    fp2mul_uRadix(phiR->X, phiR->Xi, phiR->Z, phiR->Zi, phiR->X, phiR->Xi);

    // Format public key
    fp2_encode(phiP->X, phiP->Xi, PublicKeyB);
    fp2_encode(phiQ->X, phiQ->Xi, PublicKeyB + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, phiR->Xi, PublicKeyB + 2*FP2_ENCODED_BYTES);

    return 0;
}

int EphemeralSecretAgreement_A(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyB, unsigned char* SharedSecretA)
{ // Alice's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretA using her secret key PrivateKeyA and Bob's public key PublicKeyB
  // Inputs: Alice's PrivateKeyA is an integer in the range [0, oA-1]. 
  //         Bob's PublicKeyB consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretA that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeff[6], PKB[6], jinv, jinvi;
    f2elm_t A24plus = {0}, A24plusi = {0}, C24 = {0}, C24i = {0}, A = {0}, Ai = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
	felm_t zero[1] = {0};
    
    // Initialize images of Bob's basis
    fp2_decode(PublicKeyB, PKB[0], PKB[1]);
    fp2_decode(PublicKeyB + FP2_ENCODED_BYTES, PKB[2], PKB[3]);
    fp2_decode(PublicKeyB + 2*FP2_ENCODED_BYTES, PKB[4], PKB[5]);
	
	normal2_t c;
//	from_fp2uRadix(PKB[4], PKB[5], c);
//	printf("\n");
//	for(i=0;i<10;i++){
//		printf("PKB0 c[%d] = %#018llx \n;",i,c[0][i]);
//		printf("PKB1 c[%d] = %#018llx \n;",i,c[1][i]);
//	};
//	printf("\n");
	

    // Initialize constants: A24plus = A+2C, C24 = 4C, where C=1
    get_A(PKB[0], PKB[1], PKB[2], PKB[3], PKB[4], PKB[5], A, Ai);
	

//	from_fp2uRadix(A, Ai, c);
//	printf("\n");
//	for(i=0;i<10;i++){
//		printf("A0 c[%d] = %#018llx \n;",i,c[0][i]);
//		printf("A1 c[%d] = %#018llx \n;",i,c[1][i]);
//	};
//	printf("\n");
	
    fp_add_new(zero[0], (digit_t*)&uRadix_one, zero[0], (digit_t*)&uRadix_one, C24[0], C24[1]);
    fp2add_new(A, Ai, C24, C24i, A24plus, A24plusi);
    fp_add_new(C24[1], C24[0], C24[1], C24[0], C24[0], C24[1]);

    // Retrieve kernel point
    LADDER3PT(PKB[0], PKB[1], PKB[2], PKB[3], PKB[4], PKB[5], (digit_t*)PrivateKeyA, ALICE, R, A, Ai);   

#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

    xDBLe(R, S, A24plus, A24plusi, C24, C24i, (int)(OALICE_BITS-1));
    get_2_isog(S, A24plus, A24plusi, C24, C24i);
    eval_2_isog(R, S);
#endif

    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
			fp2copy(R->Xi, pts[npts]->Xi);
            fp2copy(R->Zi, pts[npts]->Zi);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(R, R, A24plus, A24plusi, C24, C24i, (int)(2*m));
            index += m;
        }
        get_4_isog(R, A24plus, A24plusi, C24, C24i, coeff);        

        for (i = 0; i < npts; i++) {
            eval_4_isog(pts[i], coeff);
        }

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
		fp2copy(pts[npts-1]->Xi, R->Xi); 
        fp2copy(pts[npts-1]->Zi, R->Zi);
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(R, A24plus, A24plusi, C24, C24i, coeff); 
    fp2add_new(A24plus, A24plusi, A24plus, A24plusi, A24plus, A24plusi);                                    
    fp2sub_new(A24plus, A24plusi, C24, C24i, A24plus, A24plusi); 
    fp2add_new(A24plus, A24plusi, A24plus, A24plusi, A24plus, A24plusi);                    
    j_inv(A24plus, A24plusi, C24, C24i, jinv, jinvi);
    fp2_encode(jinv, jinvi, SharedSecretA);    // Format shared secret

    return 0;
}


int EphemeralSecretAgreement_B(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedSecretB)
{ // Bob's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
  // Inputs: Bob's PrivateKeyB is an integer in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  //         Alice's PublicKeyA consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretB that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, pts[MAX_INT_POINTS_BOB];
    f2elm_t coeff[6], PKB[6], jinv, jinvi;
    f2elm_t A24plus = {0}, A24plusi = {0}, A24minus = {0}, A24minusi = {0}, A = {0}, Ai = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
	felm_t zero[1] = {0};
      
    // Initialize images of Alice's basis
    fp2_decode(PublicKeyA, PKB[0], PKB[1]);
    fp2_decode(PublicKeyA + FP2_ENCODED_BYTES, PKB[2], PKB[3]);
    fp2_decode(PublicKeyA + 2*FP2_ENCODED_BYTES, PKB[4], PKB[5]);

    // Initialize constants: A24plus = A+2C, A24minus = A-2C, where C=1
    get_A(PKB[0], PKB[1], PKB[2], PKB[3], PKB[4], PKB[5], A, Ai);
    fp_add_new(zero[0], (digit_t*)&uRadix_one, zero[0], (digit_t*)&uRadix_one, A24minus[0], A24minus[1]);
    fp2add_new(A, Ai, A24minus, A24minusi, A24plus, A24plusi);
    fp2sub_new(A, Ai, A24minus, A24minusi, A24minus, A24minusi);

    // Retrieve kernel point
    LADDER3PT(PKB[0], PKB[1], PKB[2], PKB[3], PKB[4], PKB[5], (digit_t*)PrivateKeyB, BOB, R, A, Ai);
    
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
			fp2copy(R->Xi, pts[npts]->Xi);
            fp2copy(R->Zi, pts[npts]->Zi);
            pts_index[npts++] = index;
            m = strat_Bob[ii++];
            xTPLe(R, R, A24minus, A24minusi, A24plus, A24plusi, (int)m);
            index += m;
        }
        get_3_isog(R, A24minus, A24minusi, A24plus, A24plusi, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog(pts[i], coeff);
        } 

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
		fp2copy(pts[npts-1]->Xi, R->Xi); 
        fp2copy(pts[npts-1]->Zi, R->Zi);
        index = pts_index[npts-1];
        npts -= 1;
    }
     
    get_3_isog(R, A24minus, A24minusi, A24plus, A24plusi, coeff);    
    fp2add_new(A24plus, A24plusi, A24minus, A24minusi, A, Ai);                 
    fp2add_new(A, Ai, A, Ai, A, Ai);
    fp2sub_new(A24plus, A24plusi, A24minus, A24minusi, A24plus, A24plusi);                   
    j_inv(A, Ai, A24plus, A24plusi, jinv, jinvi);
    fp2_encode(jinv, jinvi, SharedSecretB);    // Format shared secret

    return 0;
}