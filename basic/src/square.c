//This just squares input X element-wise.
//For complex data, this is |X|.^2, i.e. Xr*Xr + Xi*Xi.

//This does not use in-place, since not good for complex -> real.

//To test compile:
//gcc -c square.c -O2 -std=c99 -Wall -Wextra
//clang -c square.c -O2 -std=c99 -Weverything
//g++ -c square.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c square.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int square_s (float *Y, const float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in square_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = X[n]*X[n]; }

    return 0;
}


int square_d (double *Y, const double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in square_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = X[n]*X[n]; }
    
    return 0;
}


int square_c (float *Y, const float *X, const int N)
{
    int n, n2;

    //Checks
    if (N<0) { fprintf(stderr,"error in square_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0, n2=0; n<N; n++, n2+=2) { Y[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; }
    
    return 0;
}


int square_z (double *Y, const double *X, const int N)
{
    int n, n2;

    //Checks
    if (N<0) { fprintf(stderr,"error in square_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0, n2=0; n<N; n++, n2+=2) { Y[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; }
    
    return 0;
}


#ifdef __cplusplus
}
#endif

