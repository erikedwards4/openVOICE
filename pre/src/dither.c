//Based on Kaldi, dithering for floating-point numbers is
//just adding Gaussian random noise: x[n] += d*randn,
//where d is a dither weight, for which they suggest 0.1 or 1.0.

//However, I don't currently know of a good random number generator for C.
//So, here, I use drand48 and an idea from Knuth, which only requires
//include <time.h>, which I have included anyway.
//Turns out, this is identical to Kaldi (except more robust 1.0-drand48 within log).

//To test compile:
//(Note that gcc must be used with no -std flag, and clang doesn't work.)
//gcc -c dither.c -O2 -Wall -Wextra
//g++ -c dither.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c dither.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int dither_s (float *X, const int N, const float d)
{
    const float pi2 = 2.0f*M_PIf;
    int n;
    struct timespec ts;

    //Checks
    if (N<1) { fprintf(stderr,"error in dither_s: N (num elements X) must be positive\n"); return 1; }
    if (d<0.0f) { fprintf(stderr,"error in dither_s: d (dither weight) must be positive\n"); return 1; }

    if (d>FLT_EPSILON)
    {
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr,"error in dither_s: timespec_get failed. "); perror("timespec_get"); return 1; }
	    srand48(ts.tv_nsec^ts.tv_sec); //seed random number generator
    
        for (n=0; n<N; n++)
        {
            X[n] += d * sqrtf(-2.0f*logf(1.0f-(float)drand48())) * cosf(pi2*(float)drand48());
        }
    }

    return 0;
}


int dither_d (double *X, const int N, const double d)
{
    const double pi2 = 2.0*M_PI;
    int n;
    struct timespec ts;

    //Checks
    if (N<1) { fprintf(stderr,"error in dither_d: N (num elements X) must be positive\n"); return 1; }
    if (d<0.0) { fprintf(stderr,"error in dither_d: d (dither weight) must be positive\n"); return 1; }

    if (d>DBL_EPSILON)
    {
	    if (timespec_get(&ts,TIME_UTC)==0) { fprintf(stderr,"error in dither_d: timespec_get failed. "); perror("timespec_get"); return 1; }
	    srand48(ts.tv_nsec^ts.tv_sec); //seed random number generator
    
        for (n=0; n<N; n++)
        {
            X[n] += d * sqrt(-2.0*log(1.0-drand48())) * cos(pi2*drand48());
        }
    }

    return 0;
}


#ifdef __cplusplus
}
#endif

