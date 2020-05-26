//Note that this is a continuation of a series of generalized symmetric cosine windows,
//with [Wikipedia]: a0 = 0.35875, a1 = 0.48829, a2 = 0.14128, a3 = 0.01168;

//To test compile:
//gcc -c blackmanharris.c -O2 -std=c99 -Wall -Wextra
//clang -c blackmanharris.c -O2 -std=c99 -Weverything
//g++ -c blackmanharris.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c blackmanharris.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int blackmanharris_s (float *X, const int L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    if (L<2) { fprintf(stderr,"error in blackmanharris_s: L must be > 1 \n"); return 1; }
    for (int l=0; l<L; l++) { X[l] = 0.35875f - 0.48829f*cosf(p*l) + 0.14128f*cosf(2.0f*p*l) - 0.01168f*cosf(3.0f*p*l); }
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


int blackmanharris_d (double *X, const int L, const char normalize)
{
    const double p = 2.0*M_PI/(L-1.0);
    if (L<2) { fprintf(stderr,"error in blackmanharris_d: L must be > 1 \n"); return 1; }
    for (int l=0; l<L; l++) { X[l] = 0.35875 - 0.48829*cos(p*l) + 0.14128*cos(2.0*p*l) - 0.01168*cos(3.0*p*l); }
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


int blackmanharris_c (float *X, const int L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    if (L<2) { fprintf(stderr,"error in blackmanharris_c: L must be > 1 \n"); return 1; }
    for (int l=0; l<L; l++) { X[2*l] = 0.35875f - 0.48829f*cosf(p*l) + 0.14128f*cosf(2.0f*p*l) - 0.01168f*cosf(3.0f*p*l); }
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


int blackmanharris_z (double *X, const int L, const char normalize)
{
    const double p = 2.0*M_PI/(L-1.0);
    if (L<2) { fprintf(stderr,"error in blackmanharris_z: L must be > 1 \n"); return 1; }
    for (int l=0; l<L; l++) { X[2*l] = 0.35875 - 0.48829*cos(p*l) + 0.14128*cos(2.0*p*l) - 0.01168*cos(3.0*p*l); }
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

