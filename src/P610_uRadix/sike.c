/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: supersingular isogeny key encapsulation (SIKE) protocol
*********************************************************************************************/ 

#include <string.h>
#include "../sha3/fips202.h"


int crypto_kem_keypair(unsigned char *pk, unsigned char *sk)
{ // SIKE's key generation
  // Outputs: secret key sk (CRYPTO_SECRETKEYBYTES = MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes)
  //          public key pk (CRYPTO_PUBLICKEYBYTES bytes) 
int i;
    // Generate lower portion of secret key sk <- s||SK
    randombytes(sk, MSG_BYTES);
    random_mod_order_B(sk + MSG_BYTES);
	
//	sk[61] = 0x05 ;sk[60] = 0xb7 ;sk[59] = 0x9c ;sk[58] = 0x88 ;sk[57] = 0x2b ;sk[56] = 0x4f ;sk[55] = 0xe4 ;sk[54] = 0x55 ;sk[53] = 0xc6 ;sk[52] = 0xe5 ;sk[51] = 0xa5 ;sk[50] = 0x3f ;sk[49] = 0xb9 ;sk[48] = 0x0f ;sk[47] = 0x31 ;sk[46] = 0x66 ;sk[45] = 0xc2 ;sk[44] = 0xa4 ;sk[43] = 0x91 ;sk[42] = 0x84 ;sk[41] = 0xc8 ;sk[40] = 0xb9 ;sk[39] = 0x0e ;sk[38] = 0x39 ;sk[37] = 0xc3 ;sk[36] = 0xb0 ;sk[35] = 0x0d ;sk[34] = 0x40 ;sk[33] = 0x66 ;sk[32] = 0x09 ;sk[31] = 0x99 ;sk[30] = 0xe9 ;sk[29] = 0x4f ;sk[28] = 0xef ;sk[27] = 0xb5 ;sk[26] = 0x74 ;sk[25] = 0x3c ;sk[24] = 0x64 ;sk[23] = 0x4a ;sk[22] = 0x5a ;sk[21] = 0xa6 ;sk[20] = 0xd7 ;sk[19] = 0x38 ;sk[18] = 0x8c ;sk[17] = 0x43 ;sk[16] = 0x16 ;sk[15] = 0xb7 ;sk[14] = 0x79 ;sk[13] = 0xc9 ;sk[12] = 0x69 ;sk[11] = 0x10 ;sk[10] = 0xe8 ;sk[9] = 0xfa ;sk[8] = 0x97 ;sk[7] = 0x54 ;sk[6] = 0x90 ;sk[5] = 0xfd ;sk[4] = 0xdb ;sk[3] = 0x80 ;sk[2] = 0x16 ;sk[1] = 0xa0 ;sk[0] = 0xb9 ;

    // Generate public key pk
    EphemeralKeyGeneration_B(sk + MSG_BYTES, pk);
	
//	printf("\n");
//	printf("pk = ");
//	for(i=461;i>=0;i--){
//		printf("%02x",pk[i]);
//	};
//	printf("\n");

    // Append public key pk to secret key sk
    memcpy(&sk[MSG_BYTES + SECRETKEY_B_BYTES], pk, CRYPTO_PUBLICKEYBYTES);

    return 0;
}


int crypto_kem_enc(unsigned char *ct, unsigned char *ss, const unsigned char *pk)
{ // SIKE's encapsulation
  // Input:   public key pk         (CRYPTO_PUBLICKEYBYTES bytes)
  // Outputs: shared secret ss      (CRYPTO_BYTES bytes)
  //          ciphertext message ct (CRYPTO_CIPHERTEXTBYTES = CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes)
    unsigned char ephemeralsk[SECRETKEY_A_BYTES];
    unsigned char jinvariant[FP2_ENCODED_BYTES];
    unsigned char h[MSG_BYTES];
    unsigned char temp[CRYPTO_CIPHERTEXTBYTES+MSG_BYTES];
	int i;

    // Generate ephemeralsk <- G(m||pk) mod oA 
    randombytes(temp, MSG_BYTES);
	
//	temp[23] = 0xe2;temp[22] = 0xc9;temp[21] = 0xa0;temp[20] = 0x62;temp[19] = 0x65;temp[18] = 0x45;temp[17] = 0x5d;temp[16] = 0x12;temp[15] = 0x14;temp[14] = 0xf9;temp[13] = 0xf9;temp[12] = 0xf4;temp[11] = 0xaa;temp[10] = 0x40;temp[9] = 0x8b;temp[8] = 0x09;temp[7] = 0xf4;temp[6] = 0x01;temp[5] = 0x75;temp[4] = 0xc5;temp[3] = 0x68;temp[2] = 0xbc;temp[1] = 0xa3;temp[0] = 0xfb;
	
    memcpy(&temp[MSG_BYTES], pk, CRYPTO_PUBLICKEYBYTES);
    shake256(ephemeralsk, SECRETKEY_A_BYTES, temp, CRYPTO_PUBLICKEYBYTES+MSG_BYTES);
    ephemeralsk[SECRETKEY_A_BYTES - 1] &= MASK_ALICE;

    // Encrypt
    EphemeralKeyGeneration_A(ephemeralsk, ct);
	
//	printf("\n");
//	printf("ct = ");
//	for(i=485;i>=0;i--){
//		printf("%02x",ct[i]);
//	};
//	printf("\n");

    EphemeralSecretAgreement_A(ephemeralsk, pk, jinvariant);
	
//	printf("\n");
//	printf("jinvariant = ");
//	for(i=153;i>=0;i--){
//		printf("%02x",jinvariant[i]);
//	};
//	printf("\n");

    shake256(h, MSG_BYTES, jinvariant, FP2_ENCODED_BYTES);
    for (i = 0; i < MSG_BYTES; i++) ct[i + CRYPTO_PUBLICKEYBYTES] = temp[i] ^ h[i];

    // Generate shared secret ss <- H(m||ct)
    memcpy(&temp[MSG_BYTES], ct, CRYPTO_CIPHERTEXTBYTES);
    shake256(ss, CRYPTO_BYTES, temp, CRYPTO_CIPHERTEXTBYTES+MSG_BYTES);

    return 0;
}


int crypto_kem_dec(unsigned char *ss, const unsigned char *ct, const unsigned char *sk)
{ // SIKE's decapsulation
  // Input:   secret key sk         (CRYPTO_SECRETKEYBYTES = MSG_BYTES + SECRETKEY_B_BYTES + CRYPTO_PUBLICKEYBYTES bytes)
  //          ciphertext message ct (CRYPTO_CIPHERTEXTBYTES = CRYPTO_PUBLICKEYBYTES + MSG_BYTES bytes) 
  // Outputs: shared secret ss      (CRYPTO_BYTES bytes)
    unsigned char ephemeralsk_[SECRETKEY_A_BYTES];
    unsigned char jinvariant_[FP2_ENCODED_BYTES];
    unsigned char h_[MSG_BYTES];
    unsigned char c0_[CRYPTO_PUBLICKEYBYTES];
    unsigned char temp[CRYPTO_CIPHERTEXTBYTES+MSG_BYTES];
	int i;

    // Decrypt
    EphemeralSecretAgreement_B(sk + MSG_BYTES, ct, jinvariant_);
    shake256(h_, MSG_BYTES, jinvariant_, FP2_ENCODED_BYTES);
    for (i = 0; i < MSG_BYTES; i++) temp[i] = ct[i + CRYPTO_PUBLICKEYBYTES] ^ h_[i];

    // Generate ephemeralsk_ <- G(m||pk) mod oA
    memcpy(&temp[MSG_BYTES], &sk[MSG_BYTES + SECRETKEY_B_BYTES], CRYPTO_PUBLICKEYBYTES);
    shake256(ephemeralsk_, SECRETKEY_A_BYTES, temp, CRYPTO_PUBLICKEYBYTES+MSG_BYTES);
    ephemeralsk_[SECRETKEY_A_BYTES - 1] &= MASK_ALICE;
    
    // Generate shared secret ss <- H(m||ct) or output ss <- H(s||ct)
    EphemeralKeyGeneration_A(ephemeralsk_, c0_);
    if (memcmp(c0_, ct, CRYPTO_PUBLICKEYBYTES) != 0) {
        memcpy(temp, sk, MSG_BYTES);
    }
    memcpy(&temp[MSG_BYTES], ct, CRYPTO_CIPHERTEXTBYTES);
    shake256(ss, CRYPTO_BYTES, temp, CRYPTO_CIPHERTEXTBYTES+MSG_BYTES);

    return 0;
}