//To test compile:
//gcc -c hamming.c -O2 -std=c99 -Wall -Wextra
//clang -c hamming.c -O2 -std=c99 -Weverything
//g++ -c hamming.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c hamming.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int hamming_s (float *X, const int L, const char normalize)
{
    const float p = 2.0f*M_PIf/(L-1.0f);
    if (L<2) { fprintf(stderr,"error in hamming_s: L must be > 1 \n"); return 1; }
    for (int l=0; l<L; l++) { X[l] = 0.54f - 0.46f*cosf(p*l); }
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
    if (L<2) { fprintf(stderr,"error in hamming_d: L must be > 1 \n"); return 1; }
    for (int l=0; l<L; l++) { X[l] = 0.54 - 0.46*cos(p*l); }
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
    if (L<2) { fprintf(stderr,"error in hamming_c: L must be > 1 \n"); return 1; }
    for (int l=0; l<L; l++) { X[2*l] = 0.54f - 0.46f*cosf(p*l); }
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
    if (L<2) { fprintf(stderr,"error in hamming_z: L must be > 1 \n"); return 1; }
    for (int l=0; l<L; l++) { X[2*l] = 0.54 - 0.46*cos(p*l); }
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

