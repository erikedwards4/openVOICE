//Note that this is a continuation of a series of generalized symmetric cosine windows,
//with [Wikipedia]: a0 = 0.21557895, a1 = 0.41663158, a2 = 0.277263158, a3 = 0.083578947, a4 = 0.006947368.
//This window has negative values at the edges.

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

int flattop_s (float *X, const int L, const char normalize);
int flattop_d (double *X, const int L, const char normalize);
int flattop_c (float *X, const int L, const char normalize);
int flattop_z (double *X, const int L, const char normalize);


int flattop_s (float *X, const int L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in flattop_s: L must be > 1 \n"); return 1; }

    while (l<L/2) { X[l] = 0.21557895f - 0.41663158f*cosf(p*l) + 0.277263158f*cosf(2.0f*p*l) - 0.083578947f*cosf(3.0f*p*l) + 0.006947368f*cosf(4.0f*p*l); l++; }
    if (L%2) { X[l] = 1.0f; l++; }
    while (l<L) { X[l] = X[L-l-1]; l++; }

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


int flattop_d (double *X, const int L, const char normalize)
{
    const double p = 2.0*M_PI/(L-1.0);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in flattop_d: L must be > 1 \n"); return 1; }

    while (l<L/2) { X[l] = 0.21557895 - 0.41663158*cos(p*l) + 0.277263158*cos(2.0*p*l) - 0.083578947*cos(3.0*p*l) + 0.006947368*cos(4.0*p*l); l++; }
    if (L%2) { X[l] = 1.0; l++; }
    while (l<L) { X[l] = X[L-l-1]; l++; }

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


int flattop_c (float *X, const int L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in flattop_c: L must be > 1 \n"); return 1; }

    while (l<L/2) { X[2*l] = X[2*l+1] = 0.21557895f - 0.41663158f*cosf(p*l) + 0.277263158f*cosf(2.0f*p*l) - 0.083578947f*cosf(3.0f*p*l) + 0.006947368f*cosf(4.0f*p*l); l++; }
    if (L%2) { X[2*l] = X[2*l+1] = 1.0f; l++; }
    while (l<L) { X[2*l] = X[2*l+1] = X[2*(L-l-1)]; l++; }

    if (normalize)
    {
        const float d = 1.0f;
        float sm = cblas_sdot(L,&X[0],2,&d,0);
        cblas_sscal(2*L,1.0f/sm,&X[0],1);
        sm = cblas_sdot(L,&X[0],2,&d,0);
        X[2*(L/2)] += 1.0f - sm; X[2*(L/2)+1] += 1.0f - sm;
    }

    return 0;
}


int flattop_z (double *X, const int L, const char normalize)
{
    const double p = 2.0*M_PI/(L-1.0);
    int l = 0;

    if (L<2) { fprintf(stderr,"error in flattop_z: L must be > 1 \n"); return 1; }

    while (l<L/2) { X[2*l] = X[2*l+1] = 0.21557895 - 0.41663158*cos(p*l) + 0.277263158*cos(2.0*p*l) - 0.083578947*cos(3.0*p*l) + 0.006947368*cos(4.0*p*l); l++; }
    if (L%2) { X[2*l] = X[2*l+1] = 1.0; l++; }
    while (l<L) { X[2*l] = X[2*l+1] = X[2*(L-l-1)]; l++; }

    if (normalize)
    {
        const double d = 1.0;
        double sm = cblas_ddot(L,&X[0],2,&d,0);
        cblas_dscal(2*L,1.0/sm,&X[0],1);
        sm = cblas_ddot(L,&X[0],2,&d,0);
        X[2*(L/2)] += 1.0 - sm; X[2*(L/2)+1] += 1.0 - sm;
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif

