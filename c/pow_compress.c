//This takes the output of apply_spectrogram_T_mat,
//and compresses the power in-place.

//This also works for the output of stft if desired,
//or for any non-negative data (N is total num samples in X).

//The power compression method is determined from p in [0.0 1.0].
//p==0   -> log
//p==1/3 -> cbrt
//p==1/2 -> sqrt (so changes power to amplitude)
//p==1   -> do nothing (leave as power)

//A small regularization value, preg, is added.
//This is highly recommended for log; set preg to 0.0 to have no effect.

#include <stdio.h>
#include <float.h>
#include <math.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int pow_compress_s (float *X, const int N, const float p, const float preg);
int pow_compress_d (double *X, const int N, const double p, const double preg);


int pow_compress_s (float *X, const int N, const float p, const float preg)
{
    int n;
    
    //Checks
    if (p!=p || p<0.0f || p>1.0f) { fprintf(stderr,"error in pow_compress_s: p must be in [0.0 1.0]\n"); return 1; }
    if (preg!=preg || preg<0.0f) { fprintf(stderr,"error in pow_compress_s: preg must be nonnegative\n"); return 1; }

    if (1.0f-p>FLT_EPSILON)
    {
        if (p>FLT_EPSILON)
        {
            if (fabsf(p-0.5f)<=FLT_EPSILON)
            {
                for (n=0; n<N; n++) { X[n] = sqrtf(X[n]+preg); }
            }
            else if (fabsf(p-1.0f/3.0f)<=FLT_EPSILON)
            {
                for (n=0; n<N; n++) { X[n] = cbrtf(X[n]+preg); }
            }
            else
            {
                for (n=0; n<N; n++) { X[n] = powf(X[n]+preg,p); }
            }
        }
        else
        {
            for (n=0; n<N; n++) { X[n] = logf(X[n]+preg); }
        }
    }
    else
    {
        for (n=0; n<N; n++) { X[n] += preg; }
    }
    
    return 0;
}


int pow_compress_d (double *X, const int N, const double p, const double preg)
{
    int n;
    
    //Checks
    if (p!=p || p<0.0 || p>1.0) { fprintf(stderr,"error in pow_compress_d: p must be in [0.0 1.0]\n"); return 1; }
    if (preg!=preg || preg<0.0) { fprintf(stderr,"error in pow_compress_d: preg must be nonnegative\n"); return 1; }
    
    if (1.0-p>DBL_EPSILON)
    {
        if (p>DBL_EPSILON)
        {
            if (fabs(p-0.5)<=DBL_EPSILON)
            {
                for (n=0; n<N; n++) { X[n] = sqrt(X[n]+preg); }
            }
            else if (fabs(p-1.0/3.0)<=DBL_EPSILON)
            {
                for (n=0; n<N; n++) { X[n] = cbrt(X[n]+preg); }
            }
            else
            {
                for (n=0; n<N; n++) { X[n] = pow(X[n]+preg,p); }
            }
        }
        else
        {
            for (n=0; n<N; n++) { X[n] = log(X[n]+preg); }
        }
    }
    else
    {
        for (n=0; n<N; n++) { X[n] += preg; }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif

