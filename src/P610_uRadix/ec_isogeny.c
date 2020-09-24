/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: elliptic curve and isogeny functions
*********************************************************************************************/


void xDBL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t A24plusi, const f2elm_t C24, const f2elm_t C24i)
{ // Doubling of a Montgomery point in projective coordinates (X:Z).
  // Input: projective Montgomery x-coordinates P = (X1:Z1), where x1=X1/Z1 and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2).
    f2elm_t t0, t1, t2, t3;
    
    fp2sub_new(P->X, P->Xi, P->Z, P->Zi, t0, t1);        // t0  t1->i
    fp2add_new(P->X, P->Xi, P->Z, P->Zi, t2, t3);        // t2  t3->i
    fp2sqr_uRadix(t0, t1, t0, t1);                            // t0 = (X1-Z1)^2 
    fp2sqr_uRadix(t2, t3, t2, t3);                            // t1 = (X1+Z1)^2 
    fp2mul_uRadix(C24, C24i, t0, t1, Q->Z, Q->Zi);   // Z2 = C24*(X1-Z1)^2   
    fp2mul_uRadix(t2, t3, Q->Z, Q->Zi, Q->X, Q->Xi); // X2 = C24*(X1-Z1)^2*(X1+Z1)^2
    fp2sub_new(t2, t3, t0, t1, t2, t3);              // t1 = (X1+Z1)^2-(X1-Z1)^2 
    fp2mul_uRadix(A24plus, A24plusi, t2, t3, t0, t1);// t0 = A24plus*[(X1+Z1)^2-(X1-Z1)^2]
    fp2add_new(Q->Z, Q->Zi, t0, t1, Q->Z, Q->Zi);    // Z2 = A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2
    fp2mul_uRadix(Q->Z, Q->Zi, t2, t3, Q->Z, Q->Zi); // Z2 = [A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2]*[(X1+Z1)^2-(X1-Z1)^2]
}


void xDBLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t A24plusi, const f2elm_t C24, const f2elm_t C24i, const int e)
{ // Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q <- (2^e)*P.
    int i;
    
    copy_words((digit_t*)P, (digit_t*)Q, 2*2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xDBL(Q, Q, A24plus, A24plusi, C24, C24i);
    }
}

#if (OALICE_BITS % 2 == 1)

void get_2_isog(const point_proj_t P, f2elm_t A, f2elm_t Ai, f2elm_t C, f2elm_t Ci)
{ // Computes the corresponding 2-isogeny of a projective Montgomery point (X2:Z2) of order 2.
  // Input:  projective point of order two P = (X2:Z2).
  // Output: the 2-isogenous Montgomery curve with projective coefficients A/C.
    
    fp2sqr_uRadix(P->X, P->Xi, A, Ai);                           // A = X2^2
    fp2sqr_uRadix(P->Z, P->Zi, C, Ci);                           // C = Z2^2
    fp2sub_new(C, Ci, A, Ai, A, Ai);                             // A = Z2^2 - X2^2
}


void eval_2_isog(point_proj_t P, point_proj_t Q)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 2-isogeny phi.
  // Inputs: the projective point P = (X:Z) and the 2-isogeny kernel projetive point Q = (X2:Z2).
  // Output: the projective point P = phi(P) = (X:Z) in the codomain. 
    f2elm_t t0, t1, t2, t3, t4, t5, t6, t7;
    
    fp2add_new(Q->X, Q->Xi, Q->Z, Q->Zi, t0, t1);      // t0 t1->i = X2+Z2
    fp2sub_new(Q->X, Q->Xi, Q->Z, Q->Zi, t2, t3);      // t2 t3->i = X2-Z2
    fp2add_new(P->X, P->Xi, P->Z, P->Zi, t4, t5);      // t4 t5->i = X+Z
    fp2sub_new(P->X, P->Xi, P->Z, P->Zi, t6, t7);      // t6 t7->i = X-Z
    fp2mul_uRadix(t0, t1, t6, t7, t0, t1);             // t0 = (X2+Z2)*(X-Z)
    fp2mul_uRadix(t2, t3, t4, t5, t2, t3);             // t1 = (X2-Z2)*(X+Z)
    fp2add_new(t0, t1, t2, t3, t4, t5);                // t2 = (X2+Z2)*(X-Z) + (X2-Z2)*(X+Z)
    fp2sub_new(t0, t1, t2, t3, t6, t7);                // t3 = (X2+Z2)*(X-Z) - (X2-Z2)*(X+Z)
    fp2mul_uRadix(P->X, P->Xi, t4, t5, P->X, P->Xi);   // Xfinal
    fp2mul_uRadix(P->Z, P->Zi, t6, t7, P->Z, P->Zi);   // Zfinal
}

#endif

void get_4_isog(const point_proj_t P, f2elm_t A24plus, f2elm_t A24plusi, f2elm_t C24, f2elm_t C24i, f2elm_t* coeff)
{ // Computes the corresponding 4-isogeny of a projective Montgomery point (X4:Z4) of order 4.
  // Input:  projective point of order four P = (X4:Z4).
  // Output: the 4-isogenous Montgomery curve with projective coefficients A+2C/4C and the 3 coefficients 
  //         that are used to evaluate the isogeny at a point in eval_4_isog().
    
    fp2sub_new(P->X, P->Xi, P->Z, P->Zi, coeff[2], coeff[3]);         // coeff[1] = X4-Z4
    fp2add_new(P->X, P->Xi, P->Z, P->Zi, coeff[4], coeff[5]);         // coeff[2] = X4+Z4
    fp2sqr_uRadix(P->Z, P->Zi, coeff[0], coeff[1]);                   // coeff[0] = Z4^2
    fp2add_new(coeff[0], coeff[1], coeff[0], coeff[1], coeff[0], coeff[1]);   // coeff[0] = 2*Z4^2
    fp2sqr_uRadix(coeff[0], coeff[1], C24, C24i);                     // C24 = 4*Z4^4
    fp2add_new(coeff[0], coeff[1], coeff[0], coeff[1], coeff[0], coeff[1]);   // coeff[0] = 4*Z4^2
    fp2sqr_uRadix(P->X, P->Xi, A24plus, A24plusi);                    // A24plus = X4^2
    fp2add_new(A24plus, A24plusi, A24plus, A24plusi, A24plus, A24plusi);      // A24plus = 2*X4^2
    fp2sqr_uRadix(A24plus, A24plusi, A24plus, A24plusi);              // A24plus = 4*X4^4
}


void eval_4_isog(point_proj_t P, f2elm_t* coeff)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 4-isogeny phi defined 
  // by the 3 coefficients in coeff (computed in the function get_4_isog()).
  // Inputs: the coefficients defining the isogeny, and the projective point P = (X:Z).
  // Output: the projective point P = phi(P) = (X:Z) in the codomain. 
    f2elm_t t0, t1, t2, t3;
    
    fp2add_new(P->X, P->Xi, P->Z, P->Zi, t0, t1);              // t0 = X+Z
    fp2sub_new(P->X, P->Xi, P->Z, P->Zi, t2, t3);              // t1 = X-Z
    fp2mul_uRadix(t0, t1, coeff[2], coeff[3], P->X, P->Xi);    // X = (X+Z)*coeff[1]
    fp2mul_uRadix(t2, t3, coeff[4], coeff[5], P->Z, P->Zi);    // Z = (X-Z)*coeff[2]
    fp2mul_uRadix(t0, t1, t2, t3, t0, t1);                     // t0 = (X+Z)*(X-Z)
    fp2mul_uRadix(t0, t1, coeff[0], coeff[1], t0, t1);         // t0 = coeff[0]*(X+Z)*(X-Z)
    fp2add_new(P->X, P->Xi, P->Z, P->Zi, t2, t3);              // t1 = (X-Z)*coeff[2] + (X+Z)*coeff[1]
    fp2sub_new(P->X, P->Xi, P->Z, P->Zi, P->Z, P->Zi);         // Z = (X-Z)*coeff[2] - (X+Z)*coeff[1]
    fp2sqr_uRadix(t2, t3, t2, t3);                             // t1 = [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    fp2sqr_uRadix(P->Z, P->Zi, P->Z, P->Zi);                   // Z = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2
    fp2add_new(t2, t3, t0, t1, P->X, P->Xi);                   // X = coeff[0]*(X+Z)*(X-Z) + [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    fp2sub_new(P->Z, P->Zi, t0, t1, t0, t1);                   // t0 = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2 - coeff[0]*(X+Z)*(X-Z)
    fp2mul_uRadix(P->X, P->Xi, t2, t3, P->X, P->Xi);                    // Xfinal
    fp2mul_uRadix(P->Z, P->Zi, t0, t1, P->Z, P->Zi);                    // Zfinal
}


void xTPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24minusi, const f2elm_t A24plus, const f2elm_t A24plusi)              
{ // Tripling of a Montgomery point in projective coordinates (X:Z).
  // Input: projective Montgomery x-coordinates P = (X:Z), where x=X/Z and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q = 3*P = (X3:Z3).
    f2elm_t t0, t1, t2, t3, t4, t5, t6;
	f2elm_t t7, t8, t9, t10, t11, t12, t13;
                                    
    fp2sub_new(P->X, P->Xi, P->Z, P->Zi, t0, t1);               // t0 = X-Z 
    fp2sqr_uRadix(t0, t1, t4, t5);                              // t2 = (X-Z)^2           
    fp2add_new(P->X, P->Xi, P->Z, P->Zi, t2, t3);               // t1 = X+Z 
    fp2sqr_uRadix(t2, t3, t6, t7);                              // t3 = (X+Z)^2
    fp2add_new(t0, t1, t2, t3, t8, t9);                         // t4 = 2*X
    fp2sub_new(t2, t3, t0, t1, t0, t1);                         // t0 = 2*Z 
    fp2sqr_uRadix(t8, t9, t2, t3);                              // t1 = 4*X^2
    fp2sub_new(t2, t3, t6, t7, t2, t3);                         // t1 = 4*X^2 - (X+Z)^2 
    fp2sub_new(t2, t3, t4, t5, t2, t3);                         // t1 = 4*X^2 - (X+Z)^2 - (X-Z)^2
    fp2mul_uRadix(t6, t7, A24plus, A24plusi, t10, t11);         // t5 = A24plus*(X+Z)^2 
    fp2mul_uRadix(t6, t7, t10, t11, t6, t7);                    // t3 = A24plus*(X+Z)^3
    fp2mul_uRadix(A24minus, A24minusi, t4, t5, t12, t13);       // t6 = A24minus*(X-Z)^2
    fp2mul_uRadix(t4, t5, t12, t13, t4, t5);                    // t2 = A24minus*(X-Z)^3
    fp2sub_new(t4, t5, t6, t7, t6, t7);                         // t3 = A24minus*(X-Z)^3 - coeff*(X+Z)^3
    fp2sub_new(t10, t11, t12, t13, t4, t5);                     // t2 = A24plus*(X+Z)^2 - A24minus*(X-Z)^2
    fp2mul_uRadix(t2, t3, t4, t5, t2, t3);                      // t1 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    fp2add_new(t6, t7, t2, t3, t4, t5);                         // t2 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2] + A24minus*(X-Z)^3 - coeff*(X+Z)^3
    fp2sqr_uRadix(t4, t5, t4, t5);                              // t2 = t2^2
    fp2mul_uRadix(t8, t9, t4, t5, Q->X, Q->Xi);                 // X3 = 2*X*t2
    fp2sub_new(t6, t7, t2, t3, t2, t3);                         // t1 = A24minus*(X-Z)^3 - A24plus*(X+Z)^3 - [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    fp2sqr_uRadix(t2, t3, t2, t3);                              // t1 = t1^2
    fp2mul_uRadix(t0, t1, t2, t3, Q->Z, Q->Zi);                 // Z3 = 2*Z*t1
}


void xTPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24minusi, const f2elm_t A24plus, const f2elm_t A24plusi, const int e)
{ // Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q <- (3^e)*P.
    int i;
        
    copy_words((digit_t*)P, (digit_t*)Q, 2*2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xTPL(Q, Q, A24minus, A24minusi, A24plus, A24plusi);
    }
}


void get_3_isog(const point_proj_t P, f2elm_t A24minus, f2elm_t A24minusi, f2elm_t A24plus, f2elm_t A24plusi, f2elm_t* coeff)
{ // Computes the corresponding 3-isogeny of a projective Montgomery point (X3:Z3) of order 3.
  // Input:  projective point of order three P = (X3:Z3).
  // Output: the 3-isogenous Montgomery curve with projective coefficient A/C. 
    f2elm_t t0, t1, t2, t3, t4;
	f2elm_t t5, t6, t7, t8, t9;
    
    fp2sub_new(P->X, P->Xi, P->Z, P->Zi, coeff[0], coeff[1]);         // coeff0 = X-Z
    fp2sqr_uRadix(coeff[0], coeff[1], t0, t1);                        // t0 = (X-Z)^2
    fp2add_new(P->X, P->Xi, P->Z, P->Zi, coeff[2], coeff[3]);         // coeff1 = X+Z
    fp2sqr_uRadix(coeff[2], coeff[3], t2, t3);                        // t1 = (X+Z)^2
    fp2add_new(t0, t1, t2, t3, t4, t5);                               // t2 = (X+Z)^2 + (X-Z)^2
    fp2add_new(coeff[0], coeff[1], coeff[2], coeff[3], t6, t7);       // t3 = 2*X
    fp2sqr_uRadix(t6, t7, t6, t7);                                    // t3 = 4*X^2
    fp2sub_new(t6, t7, t4, t5, t6, t7);                               // t3 = 4*X^2 - (X+Z)^2 - (X-Z)^2 
    fp2add_new(t2, t3, t6, t7, t4, t5);                               // t2 = 4*X^2 - (X-Z)^2 
    fp2add_new(t6, t7, t0, t1, t6, t7);                               // t3 = 4*X^2 - (X+Z)^2
    fp2add_new(t0, t1, t6, t7, t8, t9);                               // t4 = 4*X^2 - (X+Z)^2 + (X-Z)^2 
    fp2add_new(t8, t9, t8, t9, t8, t9);                               // t4 = 2(4*X^2 - (X+Z)^2 + (X-Z)^2) 
    fp2add_new(t2, t3, t8, t9, t8, t9);                               // t4 = 8*X^2 - (X+Z)^2 + 2*(X-Z)^2
    fp2mul_uRadix(t4, t5, t8, t9, A24minus, A24minusi);               // A24minus = [4*X^2 - (X-Z)^2]*[8*X^2 - (X+Z)^2 + 2*(X-Z)^2]
    fp2add_new(t2, t3, t4, t5, t8, t9);                               // t4 = 4*X^2 + (X+Z)^2 - (X-Z)^2
    fp2add_new(t8, t9, t8, t9, t8, t9);                               // t4 = 2(4*X^2 + (X+Z)^2 - (X-Z)^2) 
    fp2add_new(t0, t1, t8, t9, t8, t9);                               // t4 = 8*X^2 + 2*(X+Z)^2 - (X-Z)^2
    fp2mul_uRadix(t6, t7, t8, t9, A24plus, A24plusi);                 // A24plus = [4*X^2 - (X+Z)^2]*[8*X^2 + 2*(X+Z)^2 - (X-Z)^2]
}


void eval_3_isog(point_proj_t Q, const f2elm_t* coeff)
{ // Computes the 3-isogeny R=phi(X:Z), given projective point (X3:Z3) of order 3 on a Montgomery curve and 
  // a point P with 2 coefficients in coeff (computed in the function get_3_isog()).
  // Inputs: projective points P = (X3:Z3) and Q = (X:Z).
  // Output: the projective point Q <- phi(Q) = (X3:Z3). 
    f2elm_t t0, t1, t2;
	f2elm_t t3, t4, t5;

    fp2add_new(Q->X, Q->Xi, Q->Z, Q->Zi, t0, t1);                     // t0 = X+Z
    fp2sub_new(Q->X, Q->Xi, Q->Z, Q->Zi, t2, t3);                     // t1 = X-Z
    fp2mul_uRadix(t0, t1, coeff[0], coeff[1], t0, t1);                // t0 = coeff0*(X+Z)
    fp2mul_uRadix(t2, t3, coeff[2], coeff[3], t2, t3);                // t1 = coeff1*(X-Z)
    fp2add_new(t0, t1, t2, t3, t4, t5);                           // t2 = coeff0*(X+Z) + coeff1*(X-Z)
    fp2sub_new(t2, t3, t0, t1, t0, t1);                           // t0 = coeff1*(X-Z) - coeff0*(X+Z)
    fp2sqr_uRadix(t4, t5, t4, t5);                          // t2 = [coeff0*(X+Z) + coeff1*(X-Z)]^2
    fp2sqr_uRadix(t0, t1, t0, t1);                          // t0 = [coeff1*(X-Z) - coeff0*(X+Z)]^2
    fp2mul_uRadix(Q->X, Q->Xi, t4, t5, Q->X, Q->Xi);        // X3final = X*[coeff0*(X+Z) + coeff1*(X-Z)]^2        
    fp2mul_uRadix(Q->Z, Q->Zi, t0, t1, Q->Z, Q->Zi);        // Z3final = Z*[coeff1*(X-Z) - coeff0*(X+Z)]^2
}


void inv_3_way(f2elm_t z1, f2elm_t z1i, f2elm_t z2, f2elm_t z2i, f2elm_t z3, f2elm_t z3i)
{ // 3-way simultaneous inversion
  // Input:  z1,z2,z3
  // Output: 1/z1,1/z2,1/z3 (override inputs).
    f2elm_t t0, t1, t2, t3;
	f2elm_t t4, t5, t6, t7;

    fp2mul_uRadix(z1, z1i, z2, z2i, t0, t1);                      // t0 = z1*z2
    fp2mul_uRadix(z3, z3i, t0, t1, t2, t3);                       // t1 = z1*z2*z3
    fp2inv_uRadix(t2, t3);                              // t1 = 1/(z1*z2*z3)
    fp2mul_uRadix(z3, z3i, t2, t3, t4, t5);                       // t2 = 1/(z1*z2) 
    fp2mul_uRadix(t4, t5, z2, z2i, t6, t7);                       // t3 = 1/z1
    fp2mul_uRadix(t4, t5, z1, z1i, z2, z2i);                      // z2 = 1/z2
    fp2mul_uRadix(t0, t1, t2, t3, z3, z3i);                       // z3 = 1/z3
    fp2copy(t6, z1);                              // z1 = 1/z1
	fp2copy(t7, z1i); 
}


void get_A(const f2elm_t xP, const f2elm_t xPi, const f2elm_t xQ, const f2elm_t xQi, const f2elm_t xR, const f2elm_t xRi, f2elm_t A, f2elm_t Ai)
{ // Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
  // Input:  the x-coordinates xP, xQ, and xR of the points P, Q and R.
  // Output: the coefficient A corresponding to the curve E_A: y^2=x^3+A*x^2+x.
    f2elm_t t0, t1, one = {0}, t2, t3, zero = {0};
    
    fpcopy((digit_t*)&uRadix_one, one[0]);
    fp2add_new(xP, xPi, xQ, xQi, t2, t3);                           // t1 = xP+xQ
    fp2mul_uRadix(xP, xPi, xQ, xQi, t0, t1);                        // t0 = xP*xQ
    fp2mul_uRadix(xR, xRi, t2, t3, A, Ai);                          // A = xR*t1
    fp2add_new(t0, t1, A, Ai, A, Ai);                               // A = A+t0
    fp2mul_uRadix(t0, t1, xR, xRi, t0, t1);                         // t0 = t0*xR
    fp2sub_new(A, Ai, one, zero, A, Ai);                       // A = A-1
    fp2add_new(t0, t1, t0, t1, t0, t1);                             // t0 = t0+t0
    fp2add_new(t2, t3, xR, xRi, t2, t3);                            // t1 = t1+xR
    fp2add_new(t0, t1, t0, t1, t0, t1);                             // t0 = t0+t0
    fp2sqr_uRadix(A, Ai, A, Ai);                                    // A = A^2
    fp2inv_uRadix(t0, t1);                                          // t0 = 1/t0
    fp2mul_uRadix(A, Ai, t0, t1, A, Ai);                            // A = A*t0
    fp2sub_new(A, Ai, t2, t3, A, Ai);                               // Afinal = A-t1
}


void j_inv(const f2elm_t A, const f2elm_t Ai, const f2elm_t C, const f2elm_t Ci, f2elm_t jinv, f2elm_t jinvi)
{ // Computes the j-invariant of a Montgomery curve with projective constant.
  // Input: A,C in GF(p^2).
  // Output: j=256*(A^2-3*C^2)^3/(C^4*(A^2-4*C^2)), which is the j-invariant of the Montgomery curve B*y^2=x^3+(A/C)*x^2+x or (equivalently) j-invariant of B'*y^2=C*x^3+A*x^2+C*x.
    f2elm_t t0, t1, t2, t3;
    
    fp2sqr_uRadix(A, Ai, jinv, jinvi);                           // jinv = A^2        
    fp2sqr_uRadix(C, Ci, t2, t3);                                // t1 = C^2
    fp2add_new(t2, t3, t2, t3, t0, t1);                          // t0 = t1+t1
    fp2sub_new(jinv, jinvi, t0, t1, t0, t1);                     // t0 = jinv-t0
    fp2sub_new(t0, t1, t2, t3, t0, t1);                          // t0 = t0-t1
    fp2sub_new(t0, t1, t2, t3, jinv, jinvi);                     // jinv = t0-t1
    fp2sqr_uRadix(t2, t3, t2, t3);                               // t1 = t1^2
    fp2mul_uRadix(jinv, jinvi, t2, t3, jinv, jinvi);             // jinv = jinv*t1
    fp2add_new(t0, t1, t0, t1, t0, t1);                          // t0 = t0+t0
    fp2add_new(t0, t1, t0, t1, t0, t1);                          // t0 = t0+t0
    fp2sqr_uRadix(t0, t1, t2, t3);                               // t1 = t0^2
    fp2mul_uRadix(t0, t1, t2, t3, t0, t1);                       // t0 = t0*t1
    fp2add_new(t0, t1, t0, t1, t0, t1);                          // t0 = t0+t0
    fp2add_new(t0, t1, t0, t1, t0, t1);                          // t0 = t0+t0
    fp2inv_uRadix(jinv, jinvi);                                  // jinv = 1/jinv 
    fp2mul_uRadix(jinv, jinvi, t0, t1, jinv, jinvi);             // jinv = t0*jinv
}


void xDBLADD(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t xPQi, const f2elm_t A24, const f2elm_t A24i)
{ // Simultaneous doubling and differential addition.
  // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP. 
    f2elm_t t0, t1, t2, t3, t4, t5;

    fp2add_new(P->X, P->Xi, P->Z, P->Zi, t0, t1);                         // t0 = XP+ZP
    fp2sub_new(P->X, P->Xi, P->Z, P->Zi, t2, t3);                         // t1 = XP-ZP
    fp2sqr_uRadix(t0, t1, P->X, P->Xi);                                   // XP = (XP+ZP)^2
    fp2sub_new(Q->X, Q->Xi, Q->Z, Q->Zi, t4, t5);                         // t2 = XQ-ZQ
    fp2add_new(Q->X, Q->Xi, Q->Z, Q->Zi, Q->X, Q->Xi);                    // XQ = XQ+ZQ
    fp2mul_uRadix(t0, t1, t4, t5, t0, t1);                                // t0 = (XP+ZP)*(XQ-ZQ)
    fp2sqr_uRadix(t2, t3, P->Z, P->Zi);                                   // ZP = (XP-ZP)^2
    fp2mul_uRadix(t2, t3, Q->X, Q->Xi, t2, t3);                           // t1 = (XP-ZP)*(XQ+ZQ)
    fp2sub_new(P->X, P->Xi, P->Z, P->Zi, t4, t5);                         // t2 = (XP+ZP)^2-(XP-ZP)^2
    fp2mul_uRadix(P->X, P->Xi, P->Z, P->Zi, P->X, P->Xi);                 // XP = (XP+ZP)^2*(XP-ZP)^2
    fp2mul_uRadix(t4, t5, A24, A24i, Q->X, Q->Xi);                        // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sub_new(t0, t1, t2, t3, Q->Z, Q->Zi);                      // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    fp2add_new(Q->X, Q->Xi, P->Z, P->Zi, P->Z, P->Zi);            // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    fp2add_new(t0, t1, t2, t3, Q->X, Q->Xi);                      // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    fp2mul_uRadix(P->Z, P->Zi, t4, t5, P->Z, P->Zi);              // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sqr_uRadix(Q->Z, Q->Zi, Q->Z, Q->Zi);                      // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
    fp2sqr_uRadix(Q->X, Q->Xi, Q->X, Q->Xi);                      // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
    fp2mul_uRadix(Q->Z, Q->Zi, xPQ, xPQi, Q->Z, Q->Zi);           // ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
}


static void swap_points(point_proj_t P, point_proj_t Q, const digit_t option)
{ // Swap points.
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    digit_t temp;
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++) {
        temp = option & (P->X[0][i] ^ Q->X[0][i]);
        P->X[0][i] = temp ^ P->X[0][i]; 
        Q->X[0][i] = temp ^ Q->X[0][i]; 
        temp = option & (P->Z[0][i] ^ Q->Z[0][i]);
        P->Z[0][i] = temp ^ P->Z[0][i]; 
        Q->Z[0][i] = temp ^ Q->Z[0][i]; 
        temp = option & (P->X[1][i] ^ Q->X[1][i]);
        P->X[1][i] = temp ^ P->X[1][i]; 
        Q->X[1][i] = temp ^ Q->X[1][i]; 
        temp = option & (P->Z[1][i] ^ Q->Z[1][i]);
        P->Z[1][i] = temp ^ P->Z[1][i]; 
        Q->Z[1][i] = temp ^ Q->Z[1][i]; 
		
		temp = option & (P->Xi[0][i] ^ Q->Xi[0][i]);
        P->Xi[0][i] = temp ^ P->Xi[0][i]; 
        Q->Xi[0][i] = temp ^ Q->Xi[0][i]; 
        temp = option & (P->Zi[0][i] ^ Q->Zi[0][i]);
        P->Zi[0][i] = temp ^ P->Zi[0][i]; 
        Q->Zi[0][i] = temp ^ Q->Zi[0][i]; 
        temp = option & (P->Xi[1][i] ^ Q->Xi[1][i]);
        P->Xi[1][i] = temp ^ P->Xi[1][i]; 
        Q->Xi[1][i] = temp ^ Q->Xi[1][i]; 
        temp = option & (P->Zi[1][i] ^ Q->Zi[1][i]);
        P->Zi[1][i] = temp ^ P->Zi[1][i]; 
        Q->Zi[1][i] = temp ^ Q->Zi[1][i]; 
    }
}


static void LADDER3PT(const f2elm_t xP, const f2elm_t xPi, const f2elm_t xQ, const f2elm_t xQi, const f2elm_t xPQ, const f2elm_t xPQi, const digit_t* m, const unsigned int AliceOrBob, point_proj_t R, const f2elm_t A, const f2elm_t Ai)
{
    point_proj_t R0 = {0}, R2 = {0};
    f2elm_t A24 = {0};f2elm_t A24i = {0};
    digit_t mask;
    int i, nbits, bit, swap, prevbit = 0;

    if (AliceOrBob == ALICE) {
        nbits = OALICE_BITS;
    } else {
        nbits = OBOB_BITS - 1;
    }

    // Initializing constant
    fpcopy((digit_t*)&uRadix_one, A24[0]);
    fp2add_new(A24, A24i, A24, A24i, A24, A24i);
    fp2add_new(A, Ai, A24, A24i, A24, A24i);
    fp2div2(A24, A24i, A24, A24i);  
    fp2div2(A24, A24i, A24, A24i);  // A24 = (A+2)/4

    // Initializing points
    fp2copy(xQ, R0->X);
	fp2copy(xQi, R0->Xi);
    fpcopy((digit_t*)&uRadix_one, (digit_t*)R0->Z);
    fp2copy(xPQ, R2->X);
	fp2copy(xPQi, R2->Xi);
    fpcopy((digit_t*)&uRadix_one, (digit_t*)R2->Z);
    fp2copy(xP, R->X);
	fp2copy(xPi, R->Xi);
    fpcopy((digit_t*)&uRadix_one, (digit_t*)R->Z);
	fpzero((digit_t*)(R->Z)[1]);
	fpzero((digit_t*)(R->Zi)[0]);
	fpzero((digit_t*)(R->Zi)[1]);

    // Main loop
    for (i = 0; i < nbits; i++) {//i < nbits
        bit = (m[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(R, R2, mask);
		
		//	normal2_t c;
	    //    from_fp2uRadix(R->Z, R->Zi, c);
	    //    printf("\n");
	    //    for(i=0;i<10;i++){
	    //    	printf("R->Z0 c[%d] = %#018llx \n;",i,c[0][i]);
	    //    	printf("R->Z1 c[%d] = %#018llx \n;",i,c[1][i]);
	    //    };
	    //    printf("\n");
		//	
		//	from_fp2uRadix(R2->Z, R2->Zi, c);
	    //    printf("\n");
	    //    for(i=0;i<10;i++){
	    //    	printf("R2->Z0 c[%d] = %#018llx \n;",i,c[0][i]);
	    //    	printf("R2->Z1 c[%d] = %#018llx \n;",i,c[1][i]);
	    //    };
	    //    printf("\n");
		
        xDBLADD(R0, R2, R->X, R->Xi, A24, A24i);
        fp2mul_uRadix(R2->X, R2->Xi, R->Z, R->Zi, R2->X, R2->Xi);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(R, R2, mask);
}
