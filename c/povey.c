//From Kaldi: "like hamming but goes to zero at edges".

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

int povey_s (float *X, const int L, const char normalize);
int povey_d (double *X, const int L, const char normalize);
int povey_c (float *X, const int L, const char normalize);
int povey_z (double *X, const int L, const char normalize);


int povey_s (float *X, const int L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in povey_s: L must be > 1 \n"); return 1; }

    while (l<L) { X[l] = powf(0.5f-0.5f*cosf(p*l),0.85f); l++; }

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


int povey_d (double *X, const int L, const char normalize)
{
    const double p = 2.0*M_PI/(L-1.0);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in povey_d: L must be > 1 \n"); return 1; }

    while (l<L) { X[l] = pow(0.5-0.5*cos(p*l),0.85); l++; }

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


int povey_c (float *X, const int L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in povey_c: L must be > 1 \n"); return 1; }

    while (l<L) { X[2*l] = powf(0.5f-0.5f*cosf(p*l),0.85f); X[2*l+1] = 0.0f; l++; }

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


int povey_z (double *X, const int L, const char normalize)
{
    const double p = 2.0*M_PI/(L-1.0);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in povey_z: L must be > 1 \n"); return 1; }

    while (l<L) { X[2*l] = pow(0.5-0.5*cos(p*l),0.85); X[2*l+1] = 0.0; l++; }

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

