//Hamming window of length L

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifndef M_PIf
   #define M_PIf 3.141592653589793238462643383279502884f
#endif

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int hamming_s (float *X, const int L, const char normalize);
int hamming_d (double *X, const int L, const char normalize);
int hamming_c (float *X, const int L, const char normalize);
int hamming_z (double *X, const int L, const char normalize);


int hamming_s (float *X, const int L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in hamming_s: L must be > 1 \n"); return 1; }

    while (l<L) { X[l] = 0.54f - 0.46f*cosf(p*l); l++; }

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


int hamming_d (double *X, const int L, const char normalize)
{
    const double p = 2.0*M_PI/(L-1.0);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in hamming_d: L must be > 1 \n"); return 1; }

    while (l<L) { X[l] = 0.54 - 0.46*cos(p*l); l++; }

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


int hamming_c (float *X, const int L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in hamming_c: L must be > 1 \n"); return 1; }

    while (l<L) { X[2*l] = 0.54f - 0.46f*cosf(p*l); X[2*l+1] = 0.0f; l++; }

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


int hamming_z (double *X, const int L, const char normalize)
{
    const double p = 2.0*M_PI/(L-1.0);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in hamming_z: L must be > 1 \n"); return 1; }

    while (l<L) { X[2*l] = 0.54 - 0.46*cos(p*l); X[2*l+1] = 0.0; l++; }

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
}
#endif

