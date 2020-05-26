//For triangular, I use Wikipedia definition, except this leads a vanishing 1st sample,
//so instead I increase L to L+1 and disclude the 1st sample, hence l+1 and (L+1) below.

//To test compile:
//gcc -c triangular.c -O2 -std=c99 -Wall -Wextra
//clang -c triangular.c -O2 -std=c99 -Weverything
//g++ -c triangular.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c triangular.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int triangular_s (float *X, const int L, const char normalize)
{
    const float den = (L+1.0f)/2.0f;
    if (L<1) { fprintf(stderr,"error in triangular_s: L must be > 0 \n"); return 1; }
    for (int l=0; l<L; l++) { X[l] = 1.0f - fabsf((l+1.0f-(L+1.0f)/2.0f)/den); }
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
    const double den = (L+1.0)/2.0;
    if (L<1) { fprintf(stderr,"error in triangular_d: L must be > 0 \n"); return 1; }
    for (int l=0; l<L; l++) { X[l] = 1.0 - fabs((l+1.0-(L+1.0)/2.0)/den); }
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
    const float den = (L+1.0f)/2.0f;
    if (L<1) { fprintf(stderr,"error in triangular_c: L must be > 0 \n"); return 1; }
    for (int l=0; l<L; l++) { X[2*l] = 1.0f - fabsf((l+1.0f-(L+1.0f)/2.0f)/den); }
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


int triangular_z (double *X, const int L, const char normalize)
{
    const double den = (L+1.0)/2.0;
    if (L<1) { fprintf(stderr,"error in triangular_z: L must be > 0 \n"); return 1; }
    for (int l=0; l<L; l++) { X[2*l] = 1.0 - fabs((l+1.0-(L+1.0)/2.0)/den); }
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

