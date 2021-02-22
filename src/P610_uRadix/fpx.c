/********************************************************************************************
* Supersingular Isogeny Key Encapsulation Library
*
* Abstract: core functions over GF(p) and GF(p^2)
*********************************************************************************************/


__inline void fpcopy(const felm_t a, felm_t c)
{ // Copy a field element, c = a.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        c[i] = a[i];
}


__inline void fpzero(felm_t a)
{ // Zero a field element, a = 0.
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++)
        a[i] = 0;
}


void copy_words(const digit_t* a, digit_t* c, const unsigned int nwords)
{ // Copy wordsize digits, c = a, where lng(a) = nwords.
    unsigned int i;
        
    for (i = 0; i < nwords; i++) {                      
        c[i] = a[i];
    }
}


void fp2copy(const f2elm_t a, f2elm_t c)
{ // Copy a GF(p^2) element, c = a.
    fpcopy(a[0], c[0]);
    fpcopy(a[1], c[1]);
}


void fp2zero(f2elm_t a)
{ // Zero a GF(p^2) element, a = 0.
    fpzero(a[0]);
    fpzero(a[1]);
}



void fp2div2(const f2elm_t a, const f2elm_t ai, f2elm_t c, f2elm_t ci)          
{ // GF(p^2) division by two, c = a/2  in GF(p^2).
    fpdiv2_610_new(a[1], a[0], c[0], c[1]);
    fpdiv2_610_new(ai[1], ai[0], ci[0], ci[1]);
}


void fp2inv_uRadix(f2elm_t a, f2elm_t ai)
{// GF(p^2) inversion in GF(p^2), a = (a0-i*a1)/(a0^2+a1^2).
 // Inputs : a = a0 + a1*i, a0 and a1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new
 // Output : a = a0 + a1*i, a0 & a1 in uncoventional radix representation, meet the same requirement with fpmul_new & fpsqr_new
    f2elm_t t1, t2;

    fpsqr(a[1], a[0], t1[0], t1[1]);                      // t10 = a0^2
    fpsqr(ai[1], ai[0], t2[0], t2[1]);                    // t11 = a1^2
    fp_add_new(t1[1], t1[0], t2[1], t2[0], t1[0], t1[1]); // t10 = a0^2+a1^2
    fpinv_new_u(t1[1], t1[0], t1[0], t1[1]);              // t10 = (a0^2+a1^2)^-1
    fpneg_new(ai[1], ai[0], ai[0], ai[1]);                // a = a0-i*a1
    fpmul(a[1], a[0], t1[1], t1[0], a[0], a[1]);
    fpmul(ai[1], ai[0], t1[1], t1[0], ai[0], ai[1]);      // a = (a0-i*a1)*(a0^2+a1^2)^-1
}


void to_fp2uRadix(const normal2_t a, f2elm_t mc, f2elm_t mci)
{// Conversion of a GF(p^2) element from the normal representation to the new data representation based on the unconventional radix.

    to_uRadix(a[0], mc[0], mc[1]);
    to_uRadix(a[1], mci[0], mci[1]);
}


void from_fp2uRadix(const f2elm_t a, const f2elm_t ai, normal2_t c)
{ // Conversion of a GF(p^2) element from the new data representation based on the unconventional radix to normal representation.

    from_uRadix(a[1], a[0], c[0]);
    from_uRadix(ai[1], ai[0], c[1]);
}