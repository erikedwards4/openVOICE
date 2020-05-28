//Gets cepstral coefficients (CCs). This is the same as dct + lifter.
//This computes "the" DCT (discrete cosine transformation) along one dimension of matrix X.
//This uses fftw3 to compute the DCT-II.
//Then it applies a lifter (cepstral domain window) to the same dimension of X.
//The lifter has a single parameter, Q, which is default 22.0 in HTK, Kaldi and Ellis.
//Only the first K CCs are output.

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>

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

int get_ccs_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct, const float Q, const int K);
int get_ccs_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct, const double Q, const int K);


int get_ccs_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct, const float Q, const int K)
{
    const float z = 0.0f;
    const float Q2 = 0.5f*Q, PQ = M_PIf/Q;
    int r, c, k;
    float *X1, *Y1; //*lift;
    fftwf_plan plan;
    //struct timespec tic, toc;

    //Checks
    if (R<1) { fprintf(stderr,"error in get_ccs_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in get_ccs_s: C (ncols X) must be positive\n"); return 1; }
    if (ndct<R && dim==0) { fprintf(stderr,"error in get_ccs_s: ndct must be >= R for dim==0\n"); return 1; }
    if (ndct<C && dim==1) { fprintf(stderr,"error in get_ccs_s: ndct must be >= C for dim==1\n"); return 1; }
    if (Q<0.0f) { fprintf(stderr,"error in get_ccs_s: Q must be positive (or 0 to skip lifter)\n"); return 1; }
    if (K<1) { fprintf(stderr,"error in get_ccs_s: K must be positive\n"); return 1; }
    if (K>ndct) { fprintf(stderr,"error in get_ccs_s: K must be <= ndct>\n"); return 1; }

    //Initialize fftwf
    X1 = fftwf_alloc_real((size_t)ndct);
    Y1 = fftwf_alloc_real((size_t)ndct);
    plan = fftwf_plan_r2r_1d(ndct,X1,Y1,FFTW_REDFT10,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in get_ccs_s: problem creating fftw plan\n"); return 1; }
    cblas_scopy(ndct,&z,0,&X1[0],1); //zero-pad

    //clock_gettime(CLOCK_REALTIME,&tic);

    //Initialize lifter
    //if (!(lift=(float *)malloc((size_t)K*sizeof(float)))) { fprintf(stderr,"error in get_ccs_s: problem with malloc for lifter. "); perror("malloc"); return 1; }
    //for (k=0; k<K; k++) { lift[k] = fmaf(Q2,sinf(k*PQ),1.0f); }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,&X1[0],1);
                fftwf_execute(plan);
                //cblas_ssbmv(CblasColMajor,CblasUpper,K,0,1.0f,&lift[0],1,&Y1[0],1,0.0f,&Y[c*K],1);
                cblas_scopy(K,&Y1[0],1,&Y[c*K],1);
            }
            if (Q>FLT_EPSILON)
            {
                for (k=0; k<K; k++) { cblas_sscal(C,fmaf(Q2,sinf(PQ*k),1.0f),&Y[k],K); }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c],C,&X1[0],1);
                fftwf_execute(plan);
                //cblas_ssbmv(CblasRowMajor,CblasUpper,K,0,1.0f,&lift[0],1,&Y1[0],1,0.0f,&Y[c],C);
                cblas_scopy(K,&Y1[0],1,&Y[c],C);
            }
            if (Q>FLT_EPSILON)
            {
                for (k=0; k<K; k++) { cblas_sscal(C,fmaf(Q2,sinf(PQ*k),1.0f),&Y[k*C],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r],R,&X1[0],1);
                fftwf_execute(plan);
                //cblas_ssbmv(CblasColMajor,CblasUpper,K,0,1.0f,&lift[0],1,&Y1[0],1,0.0f,&Y[r],R);
                cblas_scopy(K,&Y1[0],1,&Y[r],R);
            }
            if (Q>FLT_EPSILON)
            {
                for (k=0; k<K; k++) { cblas_sscal(R,fmaf(Q2,sinf(PQ*k),1.0f),&Y[k*R],1); }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,&X1[0],1);
                fftwf_execute(plan);
                //cblas_ssbmv(CblasRowMajor,CblasUpper,K,0,1.0f,&lift[0],1,&Y1[0],1,0.0f,&Y[r*K],1);
                cblas_scopy(K,&Y1[0],1,&Y[r*K],1);
            }
            if (Q>FLT_EPSILON)
            {
                for (k=0; k<K; k++) { cblas_sscal(R,fmaf(Q2,sinf(PQ*k),1.0f),&Y[k],K); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in get_ccs_s: dim must be 0 or 1.\n"); return 1;
    }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int get_ccs_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct, const double Q, const int K)
{
    const double z = 0.0;
    const double Q2 = 0.5*Q, PQ = M_PI/Q;
    int r, c, k;
    double *X1, *Y1;
    fftw_plan plan;

    //Checks
    if (R<1) { fprintf(stderr,"error in get_ccs_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in get_ccs_d: C (ncols X) must be positive\n"); return 1; }
    if (ndct<R && dim==0) { fprintf(stderr,"error in get_ccs_d: ndct must be >= R for dim==0\n"); return 1; }
    if (ndct<C && dim==1) { fprintf(stderr,"error in get_ccs_d: ndct must be >= C for dim==1\n"); return 1; }
    if (Q<0.0) { fprintf(stderr,"error in get_ccs_d: Q must be positive (or 0 to skip lifter)\n"); return 1; }
    if (K<1) { fprintf(stderr,"error in get_ccs_d: K must be positive\n"); return 1; }
    if (K>ndct) { fprintf(stderr,"error in get_ccs_d: K must be <= ndct>\n"); return 1; }

    //Initialize fftw
    X1 = fftw_alloc_real((size_t)ndct);
    Y1 = fftw_alloc_real((size_t)ndct);
    plan = fftw_plan_r2r_1d(ndct,X1,Y1,FFTW_REDFT10,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in get_ccs_s: problem creating fftw plan\n"); return 1; }
    cblas_dcopy(ndct,&z,0,&X1[0],1); //zero-pad

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,&X1[0],1);
                fftw_execute(plan);
                cblas_dcopy(K,&Y1[0],1,&Y[c*K],1);
            }
            if (Q>DBL_EPSILON)
            {
                for (k=0; k<K; k++) { cblas_dscal(C,fma(Q2,sin(PQ*k),1.0),&Y[k],K); }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c],C,&X1[0],1);
                fftw_execute(plan);
                cblas_dcopy(K,&Y1[0],1,&Y[c],C);
            }
            if (Q>DBL_EPSILON)
            {
                for (k=0; k<K; k++) { cblas_dscal(C,fma(Q2,sin(PQ*k),1.0),&Y[k*C],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r],R,&X1[0],1);
                fftw_execute(plan);
                cblas_dcopy(K,&Y1[0],1,&Y[r],R);
            }
            if (Q>DBL_EPSILON)
            {
                for (k=0; k<K; k++) { cblas_dscal(R,fma(Q2,sin(PQ*k),1.0),&Y[k*R],1); }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,&X1[0],1);
                fftw_execute(plan);
                cblas_dcopy(K,&Y1[0],1,&Y[r*K],1);
            }
            if (Q>DBL_EPSILON)
            {
                for (k=0; k<K; k++) { cblas_dscal(R,fma(Q2,sin(PQ*k),1.0),&Y[k],K); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in get_ccs_d: dim must be 0 or 1.\n"); return 1;
    }

    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

