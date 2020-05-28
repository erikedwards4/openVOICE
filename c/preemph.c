//1st-order pre-emphasis of each row or col of X according to dim.

#include <stdio.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int preemph_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const float p);
int preemph_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const double p);


int preemph_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const float p)
{
    const int N = R*C;
    int ns = N, n = N;

    //Checks
    if (R<1) { fprintf(stderr,"error in preemph_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in preemph_s: C (ncols X) must be positive\n"); return 1; }
    if (p<=0.0f) { fprintf(stderr,"error in preemph_s: p (preemph coeff) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            while (ns>0)
            {
                ns -= R;
                while (--n>ns) { X[n] -= p*X[n-1]; }
                X[n] -= p*X[n];
            }
        }
        else
        {
            while (--n>=C) { X[n] -= p*X[n-C]; }
            while (--n>=0) { X[n] -= p*X[n]; }  //Kaldi
            //while (--n>=0) { X[n] = -X[n+C]; }  //my own
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            while (--n>=R) { X[n] -= p*X[n-R]; }
            while (--n>=0) { X[n] -= p*X[n]; }  //Kaldi
            //while (--n>=0) { X[n] = -X[n+R]; }  //my own
        }
        else
        {
            while (ns>0)
            {
                ns -= C;
                while (--n>ns) { X[n] -= p*X[n-1]; }
                X[n] -= p*X[n];
            }
        }
    }
    else
    {
        fprintf(stderr,"error in preemph_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int preemph_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const double p)
{
    const int N = R*C;
    int ns = N, n = N;

    //Checks
    if (R<1) { fprintf(stderr,"error in preemph_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in preemph_d: C (ncols X) must be positive\n"); return 1; }
    if (p<=0.0) { fprintf(stderr,"error in preemph_d: p (preemph coeff) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            while (ns>0)
            {
                ns -= R;
                while (--n>ns) { X[n] -= p*X[n-1]; }
                X[n] -= p*X[n];
            }
        }
        else
        {
            while (--n>=C) { X[n] -= p*X[n-C]; }
            while (--n>=0) { X[n] -= p*X[n]; }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            while (--n>=R) { X[n] -= p*X[n-R]; }
            while (--n>=0) { X[n] -= p*X[n]; }
        }
        else
        {
            while (ns>0)
            {
                ns -= C;
                while (--n>ns) { X[n] -= p*X[n-1]; }
                X[n] -= p*X[n];
            }
        }
    }
    else
    {
        fprintf(stderr,"error in preemph_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

