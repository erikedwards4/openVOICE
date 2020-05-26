//Gets reflection coefficients (RCs) from log area ratios (LARs) along cols or rows of X.

//The complex-valued versions are defined just in case (although I don't know of their use in practice).
//They make it impossible to compile this directly under C++.

//To test compile:
//gcc -c lar2rc.c -O2 -std=c99 -Wall -Wextra
//clang -c lar2rc.c -O2 -std=c99 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif


int lar2rc_s (float *Y, const float *X, const int N)
{
    int n;
    float ex;

    //for (n=0; n<N; n++) { Y[n] = tanhf(-0.5f*X[n]); }
    for (n=0; n<N; n++) { ex = expf(-X[n]); Y[n] = (ex-1.0f)/(ex+1.0f); }

    return 0;
}


int lar2rc_d (double *Y, const double *X, const int N)
{
    int n;
    double ex;

    //for (n=0; n<N; n++) { Y[n] = tanh(-0.5*X[n]); }
    for (n=0; n<N; n++) { ex = exp(-X[n]); Y[n] = (ex-1.0)/(ex+1.0); }

    return 0;
}


int lar2rc_c (float *Y, const float *X, const int N)
{
    int n;
    _Complex float z;

    for (n=0; n<N; n++)
    {
        z = ctanhf(-0.5f*(X[2*n]+_Complex_I*X[2*n+1]));
        Y[2*n] = crealf(z); Y[2*n+1] = cimagf(z);
    }

    return 0;
}


int lar2rc_z (double *Y, const double *X, const int N)
{
    int n;
    _Complex double z;

    for (n=0; n<N; n++)
    {
        z = ctanh(-0.5*(X[2*n]+_Complex_I*X[2*n+1]));
        Y[2*n] = creal(z); Y[2*n+1] = cimag(z);
    }

    return 0;
}


#ifdef __cplusplus
}
#endif

