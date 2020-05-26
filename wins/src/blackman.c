//There is an 'exact' blackman with:
//a0 = 7938/18608 ≈ 0.42659, a1 = 9240/18608 ≈ 0.49656, and a2 = 1430/18608 ≈ 0.076849
//However:
//a0 = 0.42, a1 = 0.5, a2 = 0.08
//are more-often used (they have improved 18 db/oct fall-off but do not null the sidelobes as well [Wikepedia]).
//However, this later window has negative values at the edges, so I include the exact version here as an option.

//To test compile:
//gcc -c blackman.c -O2 -std=c99 -Wall -Wextra
//clang -c blackman.c -O2 -std=c99 -Weverything
//g++ -c blackman.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c blackman.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int blackman_s (float *X, const int L, const char normalize, const char exact)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    if (L<2) { fprintf(stderr,"error in blackman_s: L must be > 1 \n"); return 1; }
    if (exact) { for (int l=0; l<L; l++) { X[l] = 7938.0f/18608.0f - (9240.0f/18608.0f)*cosf(p*l) + (1430.0f/18608.0f)*cosf(2.0f*p*l); } }
    else { for (int l=0; l<L; l++) { X[l] = 0.42f - 0.5f*cosf(p*l) + 0.08f*cosf(2.0f*p*l); } }
    if (normalize)
    {
        const float d = 1.0f;
        float sm = cblas_sdot(L,&X[0],1,&d,0);
        cblas_sscal(L,1.0f/sm,&X[0],1);
        sm = cblas_sdot(L,&X[0],1,&d,0);
        X[L/2] += 1.0f - sm;
    }
    return 0;
}


int blackman_d (double *X, const int L, const char normalize, const char exact)
{
    const double p = 2.0*M_PI/(L-1.0);
    if (L<2) { fprintf(stderr,"error in blackman_d: L must be > 1 \n"); return 1; }
    if (exact) { for (int l=0; l<L; l++) { X[l] = 7938.0/18608.0 - (9240.0/18608.0)*cos(p*l) + (1430.0/18608.0)*cos(2.0*p*l); } }
    else { for (int l=0; l<L; l++) { X[l] = 0.42 - 0.5*cos(p*l) + 0.08*cos(2.0*p*l); } }
    if (normalize)
    {
        const double d = 1.0;
        double sm = cblas_ddot(L,&X[0],1,&d,0);
        cblas_dscal(L,1.0/sm,&X[0],1);
        sm = cblas_ddot(L,&X[0],1,&d,0);
        X[L/2] += 1.0 - sm;
    }
    return 0;
}


int blackman_c (float *X, const int L, const char normalize, const char exact)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    if (L<2) { fprintf(stderr,"error in blackman_c: L must be > 1 \n"); return 1; }
    if (exact) { for (int l=0; l<L; l++) { X[2*l] = 7938.0f/18608.0f - (9240.0f/18608.0f)*cosf(p*l) + (1430.0f/18608.0f)*cosf(2.0f*p*l); } }
    else { for (int l=0; l<L; l++) { X[2*l] = 0.42f - 0.5f*cosf(p*l) + 0.08f*cosf(2.0f*p*l); } }
    if (normalize)
    {
        const float d = 1.0f;
        float sm = cblas_sdot(L,&X[0],2,&d,0);
        cblas_sscal(L,1.0f/sm,&X[0],2);
        sm = cblas_sdot(L,&X[0],2,&d,0);
        X[2*(L/2)] += 1.0f - sm;
    }
    return 0;
}


int blackman_z (double *X, const int L, const char normalize, const char exact)
{
    const double p = 2.0*M_PI/(L-1.0);
    if (L<2) { fprintf(stderr,"error in blackman_z: L must be > 1 \n"); return 1; }
    if (exact) { for (int l=0; l<L; l++) { X[2*l] = 7938.0/18608.0 - (9240.0/18608.0)*cos(p*l) + (1430.0/18608.0)*cos(2.0*p*l); } }
    else { for (int l=0; l<L; l++) { X[2*l] = 0.42 - 0.5*cos(p*l) + 0.08*cos(2.0*p*l); } }
    if (normalize)
    {
        const double d = 1.0;
        double sm = cblas_ddot(L,&X[0],2,&d,0);
        cblas_dscal(L,1.0/sm,&X[0],2);
        sm = cblas_ddot(L,&X[0],2,&d,0);
        X[2*(L/2)] += 1.0 - sm;
    }
    return 0;
}


#ifdef __cplusplus
}
#endif

