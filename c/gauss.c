//Gauss window with length L and stdev param a = (L-1)/(2*stdev)

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

int gauss_s (float *X, const int L, const float r, const char normalize);
int gauss_d (double *X, const int L, const double r, const char normalize);
int gauss_c (float *X, const int L, const float r, const char normalize);
int gauss_z (double *X, const int L, const double r, const char normalize);


int gauss_s (float *X, const int L, const float a, const char normalize)
{
    const float p = -2.0f*(a*a/((L-1)*(L-1)));
    const float m = 0.5f*(L-1);
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in gauss_s: L must be > 1 \n"); return 1; }
    if (a<=0.0f) { fprintf(stderr,"error in gauss_s: a must be positive \n"); return 1; }

    while (l<L/2) { X[l] = expf((l-m)*(l-m)*p); l++; }
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


int gauss_d (double *X, const int L, const double a, const char normalize)
{
    const double p = -2.0*(a*a/((L-1)*(L-1)));
    const double m = 0.5*(L-1);
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in gauss_d: L must be > 1 \n"); return 1; }
    if (a<=0.0) { fprintf(stderr,"error in gauss_d: a must be positive \n"); return 1; }

    while (l<L/2) { X[l] = exp((l-m)*(l-m)*p); l++; }
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


int gauss_c (float *X, const int L, const float a, const char normalize)
{
    const float p = -2.0f*(a*a/((L-1)*(L-1)));
    const float m = 0.5f*(L-1);
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in gauss_c: L must be > 1 \n"); return 1; }
    if (a<=0.0f) { fprintf(stderr,"error in gauss_c: a must be positive \n"); return 1; }

    while (l<L/2) { X[2*l] = X[2*l+1] = expf((l-m)*(l-m)*p); l++; }
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


int gauss_z (double *X, const int L, const double a, const char normalize)
{
    const double p = -2.0*(a*a/((L-1)*(L-1)));
    const double m = 0.5*(L-1);
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in gauss_z: L must be > 1 \n"); return 1; }
    if (a<=0.0) { fprintf(stderr,"error in gauss_z: a must be psoitive \n"); return 1; }

    while (l<L/2) { X[2*l] = X[2*l+1] = exp((l-m)*(l-m)*p); l++; }
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

