//To test compile:
//gcc -c rectangular.c -O2 -std=c99 -Wall -Wextra
//clang -c rectangular.c -O2 -std=c99 -Weverything
//g++ -c rectangular.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c rectangular.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int rectangular_s (float *X, const int L, const char normalize)
{
    if (L<1) { fprintf(stderr,"error in rectangular_s: L must be > 0 \n"); return 1; }
    for (int l=0; l<L; l++) { X[l] = 1.0f; }
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


int rectangular_d (double *X, const int L, const char normalize)
{
    if (L<1) { fprintf(stderr,"error in rectangular_d: L must be > 0 \n"); return 1; }
    for (int l=0; l<L; l++) { X[l] = 1.0; }
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


int rectangular_c (float *X, const int L, const char normalize)
{
    if (L<1) { fprintf(stderr,"error in rectangular_c: L must be > 0 \n"); return 1; }
    for (int l=0; l<L; l++) { X[2*l] = 1.0f; X[2*l+1] = 0.0f; }
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


int rectangular_z (double *X, const int L, const char normalize)
{
    if (L<1) { fprintf(stderr,"error in rectangular_z: L must be > 0 \n"); return 1; }
    for (int l=0; l<L; l++) { X[2*l] = 1.0; X[2*l+1] = 0.0; }
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

