//Gets log area ratios (LARs) from the reflection coefficients (RCs) along cols or rows of X.

//The complex-valued versions are defined just in case (although I don't know of their use in practice).
//They make it impossible to compile this directly under C++.

//To test compile:
//gcc -c rc2lar.c -O2 -std=c99 -Wall -Wextra
//clang -c rc2lar.c -O2 -std=c99 -Weverything
//g++ -c rc2lar.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c rc2lar.c -O2 -std=c++11 -Weverything -Wno-old-style-cast

#include "/home/erik/codee/openvoice/openvoice.h"
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif


int rc2lar_s (float *Y, const float *X, const int N)
{
    int n;

    //for (n=0; n<N; n++) { Y[n] = -2.0f * atanhf(X[n]); }
    for (n=0; n<N; n++) { Y[n] = logf((1.0f-X[n])/(1.0f+X[n])); }

    return 0;
}


int rc2lar_d (double *Y, const double *X, const int N)
{
    int n;

    //for (n=0; n<N; n++) { Y[n] = -2.0 * atanh(X[n]); }
    for (n=0; n<N; n++) { Y[n] = log((1.0-X[n])/(1.0+X[n])); }

    return 0;
}


int rc2lar_c (float *Y, const float *X, const int N)
{
    int n;
    _Complex float z;

    for (n=0; n<N; n++)
    {
        z = -2.0f * catanhf(X[2*n]+_Complex_I*X[2*n+1]);
        Y[2*n] = crealf(z); Y[2*n+1] = cimagf(z);
    }

    return 0;
}


int rc2lar_z (double *Y, const double *X, const int N)
{
    int n;
    _Complex double z;

    for (n=0; n<N; n++)
    {
        z = -2.0 * catanh(X[2*n]+_Complex_I*X[2*n+1]);
        Y[2*n] = creal(z); Y[2*n+1] = cimag(z);
    }

    return 0;
}


#ifdef __cplusplus
}
#endif

