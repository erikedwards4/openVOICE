//This just gets absolute value of input X element-wise.
//This does not use in-place, since not good for complex -> real.

//To test compile:
//gcc -c abs.c -O2 -std=c99 -Wall -Wextra
//clang -c abs.c -O2 -std=c99 -Weverything
//g++ -c abs.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c abs.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int abs_s (float *Y, const float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in abs_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = fabsf(X[n]); }

    return 0;
}


int abs_d (double *Y, const double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in abs_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = fabs(X[n]); }
    
    return 0;
}


int abs_c (float *Y, const float *X, const int N)
{
    int n, n2;

    //Checks
    if (N<0) { fprintf(stderr,"error in abs_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0, n2=0; n<2*N; n++, n2+=2) { Y[n] = sqrtf(X[n2]*X[n2] + X[n2+1]*X[n2+1]); }
    
    return 0;
}


int abs_z (double *Y, const double *X, const int N)
{
    int n, n2;

    //Checks
    if (N<0) { fprintf(stderr,"error in abs_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0, n2=0; n<N; n++, n2+=2) { Y[n] = sqrt(X[n2]*X[n2] + X[n2+1]*X[n2+1]); }
    
    return 0;
}


#ifdef __cplusplus
}
#endif

