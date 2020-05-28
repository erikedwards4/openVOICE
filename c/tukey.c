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

int tukey_s (float *X, const int L, const float r, const char normalize);
int tukey_d (double *X, const int L, const double r, const char normalize);
int tukey_c (float *X, const int L, const float r, const char normalize);
int tukey_z (double *X, const int L, const double r, const char normalize);


int tukey_s (float *X, const int L, const float r, const char normalize)
{
    const float p = 2.0f*M_PIf/((L-1)*r);
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in tukey_s: L must be > 1 \n"); return 1; }
    if (r<0.0f || r>1.0f) { fprintf(stderr,"error in tukey_s: r must be in [0.0 1.0] \n"); return 1; }

    while (l<0.5f*r*L) { X[l] = 0.5f - 0.5f*cosf(p*l); l++; }
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


int tukey_d (double *X, const int L, const double r, const char normalize)
{
    const double p = 2.0*M_PI/((L-1)*r);
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in tukey_d: L must be > 1 \n"); return 1; }
    if (r<0.0 || r>1.0) { fprintf(stderr,"error in tukey_d: r must be in [0.0 1.0] \n"); return 1; }

    while (l<0.5*r*L) { X[l] = 0.5 - 0.5*cos(p*l); l++; }
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


int tukey_c (float *X, const int L, const float r, const char normalize)
{
    const float p = 2.0f*M_PIf/((L-1)*r);
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in tukey_s: L must be > 1 \n"); return 1; }
    if (r<0.0f || r>1.0f) { fprintf(stderr,"error in tukey_s: r must be in [0.0 1.0] \n"); return 1; }

    while (l<0.5f*r*L) { X[2*l] = 0.5f - 0.5f*cosf(p*l); X[2*l+1] = 0.0f; l++; }
    while (l<L/2) { X[2*l] = 1.0f; X[2*l+1] = 0.0f; l++; }
    if (L%2) { X[2*l] = 1.0f; X[2*l+1] = 0.0f; l++; }
    while (l<L) { X[2*l] = X[2*(L-l-1)]; X[2*l+1] = 0.0f; l++; }

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


int tukey_z (double *X, const int L, const double r, const char normalize)
{
    const double p = 2.0*M_PI/((L-1)*r);
    int l = 0;

    //Checks
    if (L<2) { fprintf(stderr,"error in tukey_d: L must be > 1 \n"); return 1; }
    if (r<0.0 || r>1.0) { fprintf(stderr,"error in tukey_d: r must be in [0.0 1.0] \n"); return 1; }

    while (l<0.5*r*L) { X[2*l] = 0.5 - 0.5*cos(p*l); X[2*l+1] = 0.0; l++; }
    while (l<L/2) { X[2*l] = 1.0; X[2*l+1] = 0.0; l++; }
    if (L%2) { X[2*l] = 1.0; X[2*l+1] = 0.0; l++; }
    while (l<L) { X[2*l] = X[2*(L-l-1)]; X[2*l+1] = 0.0; l++; }

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

