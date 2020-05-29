//For Barlett, I use Octave/Matlab convention,
//for which they cite Openheim & Shafer: Discrete-Time Signal Processing.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int bartlett_s (float *X, const int L, const char normalize);
int bartlett_d (double *X, const int L, const char normalize);
int bartlett_c (float *X, const int L, const char normalize);
int bartlett_z (double *X, const int L, const char normalize);


int bartlett_s (float *X, const int L, const char normalize)
{
    const float p = 1.0f/(L-(L%2)*(L/2)-1);
    int l = 0;

    if (L<1) { fprintf(stderr,"error in bartlett_s: L must be > 0 \n"); return 1; }

    while (l<L/2) { X[l] = (l*(2-L%2))*p; l++; }
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


int bartlett_d (double *X, const int L, const char normalize)
{
    const double p = 1.0/(L-(L%2)*(L/2)-1);
    int l = 0;

    if (L<1) { fprintf(stderr,"error in bartlett_d: L must be > 0 \n"); return 1; }

    while (l<L/2) { X[l] = (l*(2-L%2))*p; l++; }
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


int bartlett_c (float *X, const int L, const char normalize)
{
    const float p = 1.0f/(L-(L%2)*(L/2)-1);
    int l = 0;

    if (L<1) { fprintf(stderr,"error in bartlett_c: L must be > 0 \n"); return 1; }

    while (l<L/2) { X[2*l] = X[2*l+1] = (l*(2-L%2))*p; l++; }
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


int bartlett_z (double *X, const int L, const char normalize)
{
    const double p = 1.0/(L-(L%2)*(L/2)-1);
    int l = 0;

    if (L<1) { fprintf(stderr,"error in bartlett_z: L must be > 0 \n"); return 1; }

    while (l<L/2) { X[2*l] = X[2*l+1] = (l*(2-L%2))*p; l++; }
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

