//Tukey (tapered cosine) window with length L and cosine fraction r

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

int planck_s (float *X, const int L, const float epsilon, const char normalize);
int planck_d (double *X, const int L, const double epsilon, const char normalize);
int planck_c (float *X, const int L, const float epsilon, const char normalize);
int planck_z (double *X, const int L, const double epsilon, const char normalize);


int planck_s (float *X, const int L, const float epsilon, const char normalize)
{
    const float p = epsilon*L;
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in planck_s: L must be > 1 \n"); return 1; }
    if (epsilon<0.0f || epsilon>0.5f) { fprintf(stderr,"error in planck_s: epsilon must be in [0.0 0.5] \n"); return 1; }

    X[l++] = 0.0f;
    while (l<p && l<L/2) { X[l] = 1.0f/(1.0f+expf(p/l-p/(p-l))); l++; }
    while (l<L/2) { X[l] = 1.0f; l++; }
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


int planck_d (double *X, const int L, const double epsilon, const char normalize)
{
    const double p = epsilon*L;
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in planck_d: L must be > 1 \n"); return 1; }
    if (epsilon<0.0 || epsilon>0.5) { fprintf(stderr,"error in planck_d: epsilon must be in [0.0 0.5] \n"); return 1; }

    X[l++] = 0.0;
    while (l<p && l<L/2) { X[l] = 1.0/(1.0+exp(p/l-p/(p-l))); l++; }
    while (l<L/2) { X[l] = 1.0; l++; }
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


int planck_c (float *X, const int L, const float epsilon, const char normalize)
{
    const float p = epsilon*L;
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in planck_c: L must be > 1 \n"); return 1; }
    if (epsilon<0.0f || epsilon>0.5f) { fprintf(stderr,"error in planck_c: epsilon must be in [0.0 0.5] \n"); return 1; }

    X[l++] = 0.0f; X[1] = 0.0f;
    while (l<p && l<L/2) { X[2*l] = X[2*l+1] = 1.0f/(1.0f+expf(p/l-p/(p-l))); l++; }
    while (l<L/2) { X[2*l] = X[2*l+1] = 1.0f; l++; }
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


int planck_z (double *X, const int L, const double epsilon, const char normalize)
{
    const double p = epsilon*L;
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in planck_z: L must be > 1 \n"); return 1; }
    if (epsilon<0.0 || epsilon>0.5) { fprintf(stderr,"error in planck_z: epsilon must be in [0.0 0.5] \n"); return 1; }

    X[l++] = 0.0; X[1] = 0.0;
    while (l<p && l<L/2) { X[2*l] = X[2*l+1] = 1.0/(1.0+exp(p/l-p/(p-l))); l++; }
    while (l<L/2) { X[2*l] = X[2*l+1] = 1.0; l++; }
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

