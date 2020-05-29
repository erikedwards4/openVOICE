//For triangular, I used Wikipedia definition, but this leads a vanishing 1st sample.
//So, instead I just follow the Octave/Matlab convention.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int triangular_s (float *X, const int L, const char normalize);
int triangular_d (double *X, const int L, const char normalize);
int triangular_c (float *X, const int L, const char normalize);
int triangular_z (double *X, const int L, const char normalize);


int triangular_s (float *X, const int L, const char normalize)
{
    const float p = 1.0f/(L-(L%2)*(L/2));
    //const float den = (L+1.0f)/2.0f;
    int l = 0;

    if (L<1) { fprintf(stderr,"error in triangular_s: L must be > 0 \n"); return 1; }
    
    while (l<L/2) { X[l] = (l*(2-L%2)+1)*p; l++; }
    //while (l<L/2) { X[l] = 1.0f - fabsf((l+1-0.5f*(L+1))/den); l++; }
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


int triangular_d (double *X, const int L, const char normalize)
{
    const double p = 1.0/(L-(L%2)*(L/2));
    int l = 0;

    if (L<1) { fprintf(stderr,"error in triangular_d: L must be > 0 \n"); return 1; }
    
    while (l<L/2) { X[l] = (l*(2-L%2)+1)*p; l++; }
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


int triangular_c (float *X, const int L, const char normalize)
{
    const float p = 1.0f/(L-(L%2)*(L/2));
    int l = 0;

    if (L<1) { fprintf(stderr,"error in triangular_c: L must be > 0 \n"); return 1; }
    
    while (l<L/2) { X[2*l] = X[2*l+1] = (l*(2-L%2)+1)*p; l++; }
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


int triangular_z (double *X, const int L, const char normalize)
{
    const double p = 1.0/(L-(L%2)*(L/2));
    int l = 0;

    if (L<1) { fprintf(stderr,"error in triangular_z: L must be > 0 \n"); return 1; }
    
    while (l<L/2) { X[2*l] = X[2*l+1] = (l*(2-L%2)+1)*p; l++; }
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

