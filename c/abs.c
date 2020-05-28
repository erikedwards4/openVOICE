//This just gets absolute value of input X element-wise.
//This has in-place and not-in-place versions, since in-place not good for complex -> real.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int abs_s (float *Y, const float *X, const int N);
int abs_d (double *Y, const double *X, const int N);
int abs_c (float *Y, const float *X, const int N);
int abs_z (double *Y, const double *X, const int N);

int abs_inplace_s (float *X, const int N);
int abs_inplace_d (double *X, const int N);
int abs_inplace_c (float *X, const int N);
int abs_inplace_z (double *X, const int N);


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


int abs_inplace_s (float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in abs_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n] = fabsf(X[n]); }

    return 0;
}


int abs_inplace_d (double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in abs_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n] = fabs(X[n]); }
    
    return 0;
}


int abs_inplace_c (float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in abs_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<2*N; n+=2) { X[n] = sqrtf(X[n]*X[n] + X[n+1]*X[n+1]); X[n+1] = 0.0f; }
    
    return 0;
}


int abs_inplace_z (double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in abs_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<2*N; n+=2) { X[n] = sqrt(X[n]*X[n] + X[n+1]*X[n+1]); X[n+1] = 0.0; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif

