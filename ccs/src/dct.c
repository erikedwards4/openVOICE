//This computes "the" DCT (discrete cosine transformation) along one dimension of matrix X.
//This uses fftw3 to compute the DCT-II.

//To test compile:
//gcc -c dct.c -O2 -std=c99 -Wall -Wextra
//clang -c dct.c -O2 -std=c99 -Weverything
//g++ -c dct.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c dct.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int dct_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct)
{
    const float z = 0.0f;
    int r, c;
    float *X1, *Y1;
    fftwf_plan plan;
    //struct timespec tic, toc;

    //Checks
    if (R<1) { fprintf(stderr,"error in dct_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in dct_s: C (ncols X) must be positive\n"); return 1; }
    if (ndct<R && dim==0) { fprintf(stderr,"error in dct_s: ndct must be >= R for dim==0\n"); return 1; }
    if (ndct<C && dim==1) { fprintf(stderr,"error in dct_s: ndct must be >= C for dim==1\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);

    //Initialize fftwf
    X1 = fftwf_alloc_real((size_t)ndct);
    Y1 = fftwf_alloc_real((size_t)ndct);
    plan = fftwf_plan_r2r_1d(ndct,X1,Y1,FFTW_REDFT10,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in dct_s: problem creating fftw plan\n"); return 1; }
    cblas_scopy(ndct,&z,0,&X1[0],1); //zero-pad

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,&X1[0],1);
                fftwf_execute(plan);
                cblas_scopy(ndct,&Y1[0],1,&Y[c*ndct],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c],C,&X1[0],1);
                fftwf_execute(plan);
                cblas_scopy(ndct,&Y1[0],1,&Y[c],C);
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
                cblas_scopy(ndct,&Y1[0],1,&Y[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,&X1[0],1);
                fftwf_execute(plan);
                cblas_scopy(ndct,&Y1[0],1,&Y[r*ndct],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in dct_s: dim must be 0 or 1.\n"); return 1;
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    return 0;
}


int dct_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct)
{
    const double z = 0.0;
    int r, c;
    double *X1, *Y1;
    fftw_plan plan;

    //Checks
    if (R<1) { fprintf(stderr,"error in dct_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in dct_d: C (ncols X) must be positive\n"); return 1; }
    if (ndct<R && dim==0) { fprintf(stderr,"error in dct_d: ndct must be >= R for dim==0\n"); return 1; }
    if (ndct<C && dim==1) { fprintf(stderr,"error in dct_d: ndct must be >= C for dim==1\n"); return 1; }

    //Initialize fftw
    X1 = fftw_alloc_real((size_t)ndct);
    Y1 = fftw_alloc_real((size_t)ndct);
    plan = fftw_plan_r2r_1d(ndct,X1,Y1,FFTW_REDFT10,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in dct_s: problem creating fftw plan\n"); return 1; }
    cblas_dcopy(ndct,&z,0,&X1[0],1); //zero-pad

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,&X1[0],1);
                fftw_execute(plan);
                cblas_dcopy(ndct,&Y1[0],1,&Y[c*ndct],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c],C,&X1[0],1);
                fftw_execute(plan);
                cblas_dcopy(ndct,&Y1[0],1,&Y[c],C);
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
                cblas_dcopy(ndct,&Y1[0],1,&Y[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,&X1[0],1);
                fftw_execute(plan);
                cblas_dcopy(ndct,&Y1[0],1,&Y[r*ndct],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in dct_d: dim must be 0 or 1.\n"); return 1;
    }

    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
#endif

