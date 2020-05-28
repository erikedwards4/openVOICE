//This applies a lifter (cepstral domain window) to one dimension of X.
//The lifter has a single parameter, Q, which is default 22.0 in HTK, Kaldi and Ellis.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifndef M_PIf
   #define M_PIf 3.141592653589793238462643383279502884f
#endif

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int lifter_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const float Q);
int lifter_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const double Q);


int lifter_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const float Q)
{
    const float Q2 = 0.5f*Q, PQ = M_PIf/Q;
    //const int K = (dim==0) ? R : C;
    int r, c; //k;
    //float *lift;
    //struct timespec tic, toc;

    //Checks
    if (R<1) { fprintf(stderr,"error in lifter_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lifter_s: C (ncols X) must be positive\n"); return 1; }
    if (Q<=0.0f) { fprintf(stderr,"error in lifter_s: Q must be positive\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);

    //Initialize lifter
    //if (!(lift=(float *)malloc((size_t)ndct*sizeof(float)))) { fprintf(stderr,"error in get_ccs_s: problem with malloc for lifter. "); perror("malloc"); return 1; }
    //for (n=0; n<ndct; n++) { lift[n] = fmaf(Q2,sinf(n*PQ),1.0f); }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++) { cblas_sscal(C,fmaf(Q2,sinf(PQ*r),1.0f),&X[r],R); }
            //for (k=0; k<K; k++) { cblas_sscal(C,lift[k],&X[k],K); }
        }
        else
        {
            for (r=0; r<R; r++) { cblas_sscal(C,fmaf(Q2,sinf(PQ*r),1.0f),&X[r*C],1); }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++) { cblas_sscal(R,fmaf(Q2,sinf(PQ*c),1.0f),&X[c*R],1); }
        }
        else
        {
            for (c=0; c<C; c++) { cblas_sscal(R,fmaf(Q2,sinf(PQ*c),1.0f),&X[c],C); }
        }
    }
    else
    {
        fprintf(stderr,"error in lifter_s: dim must be 0 or 1.\n"); return 1;
    }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int lifter_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const double Q)
{
    const double Q2 = 0.5*Q, PQ = M_PI/Q;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in lifter_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lifter_d: C (ncols X) must be positive\n"); return 1; }
    if (Q<=0.0) { fprintf(stderr,"error in lifter_d: Q must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++) { cblas_dscal(C,fma(Q2,sin(PQ*r),1.0),&X[r],R); }
        }
        else
        {
            for (r=0; r<R; r++) { cblas_dscal(C,fma(Q2,sin(PQ*r),1.0),&X[r*C],1); }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++) { cblas_dscal(R,fma(Q2,sin(PQ*c),1.0),&X[c*R],1); }
        }
        else
        {
            for (c=0; c<C; c++) { cblas_dscal(R,fma(Q2,sin(PQ*c),1.0),&X[c],C); }
        }
    }
    else
    {
        fprintf(stderr,"error in lifter_d: dim must be 0 or 1.\n"); return 1;
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif

