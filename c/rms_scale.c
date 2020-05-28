//Normalizes RMS of each row or col of X according to dim.
//See C++ command-line code for more info.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include "cmp_ascend.c"

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int cmp_ascend_s (const void *a, const void *b);
int cmp_ascend_d (const void *a, const void *b);

int rms_scale_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const float fs, const float tau, const float target_dB_SPL, const float max_dB_SPL);
int rms_scale_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const double fs, const double tau, const double target_dB_SPL, const double max_dB_SPL);


int rms_scale_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const float fs, const float tau, const float target_dB_SPL, const float max_dB_SPL)
{
    const float rms_prctile = 0.9f;
	const float target_RMS = 20e-6f * powf(10.0f, target_dB_SPL/20.0f);
	const float max_RMS = 20e-6f * powf(10.0f, max_dB_SPL/20.0f);
	const float a = expf(-1.0f/(fs*tau));
	const float b = 1.0f - a;
    const int tranz = 5; //num samps to zero in case transient
    float p90_rms, scale, scalemx;
    int r, c, n;
    float *Y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rms_scale_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rms_scale_s: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (!(Y=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in rms_scale_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                Y[0] = b*X[c*R]*X[c*R]; r = 1;
                for (n=c*R+1; n<c*R+R; n++, r++) { Y[r] = b*X[n]*X[n] + a*Y[r-1]; }
                cblas_sscal(tranz,0.0f,&Y[0],1); //in case transient
                qsort(Y,(size_t)(R),sizeof(float),cmp_ascend_s); //gsl_sort_float(Y,1,R);
                p90_rms = sqrtf(Y[(int)(floorf(rms_prctile*R))]);
                scale = target_RMS/p90_rms; scalemx = max_RMS/sqrtf(Y[R-1]);
                scale = fminf(scale, scalemx);
                cblas_sscal(R,scale,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                Y[0] = b*X[c]*X[c]; r = 1;
                for (n=c+C; n<c+C*(R-1); n+=C, r++) { Y[r] = b*X[n]*X[n] + a*Y[r-1]; }
                cblas_sscal(tranz,0.0f,&Y[0],1); //in case transient
                qsort(Y,(size_t)(R),sizeof(float),cmp_ascend_s); //gsl_sort_float(Y,1,R);
                p90_rms = sqrtf(Y[(int)(floorf(rms_prctile*R))]);
                scale = target_RMS/p90_rms; scalemx = max_RMS/sqrtf(Y[R-1]);
                scale = fminf(scale, scalemx);
                cblas_sscal(R,scale,&X[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (!(Y=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in rms_scale_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                Y[0] = b*X[r]*X[r]; c = 1;
                for (n=r+R; n<r+R*(C-1); n+=R, c++) { Y[c] = b*X[n]*X[n] + a*Y[c-1]; }
                cblas_sscal(tranz,0.0f,&Y[0],1); //in case transient
                qsort(Y,(size_t)(C),sizeof(float),cmp_ascend_s); //gsl_sort_float(Y,1,C);
                p90_rms = sqrtf(Y[(int)(floorf(rms_prctile*C))]);
                scale = target_RMS/p90_rms; scalemx = max_RMS/sqrtf(Y[C-1]);
                scale = fminf(scale, scalemx);
                cblas_sscal(C,scale,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                Y[0] = b*X[r*C]*X[r*C]; c = 1;
                for (n=r*C+1; n<r*C+C; n++, c++) { Y[c] = b*X[n]*X[n] + a*Y[c-1]; }
                cblas_sscal(tranz,0.0f,&Y[0],1); //in case transient
                qsort(Y,(size_t)(C),sizeof(float),cmp_ascend_s); //gsl_sort_float(Y,1,C);
                p90_rms = sqrtf(Y[(int)(floorf(rms_prctile*C))]);
                scale = target_RMS/p90_rms; scalemx = max_RMS/sqrtf(Y[C-1]);
                scale = fminf(scale, scalemx);
                cblas_sscal(C,scale,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rms_scale_s: dim must be 0 or 1.\n"); return 1;
    }

    free(Y);
    return 0;
}


int rms_scale_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const double fs, const double tau, const double target_dB_SPL, const double max_dB_SPL)
{
    const double rms_prctile = 0.9;
	const double target_RMS = 20e-6 * pow(10.0, target_dB_SPL/20.0);
	const double max_RMS = 20e-6 * pow(10.0, max_dB_SPL/20.0);
	const double a = exp(-1.0/(fs*tau));
	const double b = 1.0 - a;
    const int tranz = 5; //num samps to zero in case transient
    double p90_rms, scale, scalemx;
    int r, c, n;
    double *Y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rms_scale_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rms_scale_d: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (!(Y=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in rms_scale_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                Y[0] = b*X[c*R]*X[c*R]; r = 1;
                for (n=c*R+1; n<c*R+R; n++, r++) { Y[r] = b*X[n]*X[n] + a*Y[r-1]; }
                cblas_dscal(tranz,0.0,&Y[0],1); //in case transient
                qsort(Y,(size_t)(R),sizeof(double),cmp_ascend_d);
                p90_rms = sqrt(Y[(int)(floor(rms_prctile*R))]);
                scale = target_RMS/p90_rms; scalemx = max_RMS/sqrt(Y[R-1]);
                scale = fmin(scale, scalemx);
                cblas_dscal(R,scale,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                Y[0] = b*X[c]*X[c]; r = 1;
                for (n=c+C; n<c+C*(R-1); n+=C, r++) { Y[r] = b*X[n]*X[n] + a*Y[r-1]; }
                cblas_dscal(tranz,0.0,&Y[0],1); //in case transient
                qsort(Y,(size_t)(R),sizeof(double),cmp_ascend_d);
                p90_rms = sqrt(Y[(int)(floor(rms_prctile*R))]);
                scale = target_RMS/p90_rms; scalemx = max_RMS/sqrt(Y[R-1]);
                scale = fmin(scale, scalemx);
                cblas_dscal(R,scale,&X[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (!(Y=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in rms_scale_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                Y[0] = b*X[r]*X[r]; c = 1;
                for (n=r+R; n<r+R*(C-1); n+=R, c++) { Y[c] = b*X[n]*X[n] + a*Y[c-1]; }
                cblas_dscal(tranz,0.0,&Y[0],1); //in case transient
                qsort(Y,(size_t)(C),sizeof(double),cmp_ascend_d);
                p90_rms = sqrt(Y[(int)(floor(rms_prctile*C))]);
                scale = target_RMS/p90_rms; scalemx = max_RMS/sqrt(Y[C-1]);
                scale = fmin(scale, scalemx);
                cblas_dscal(C,scale,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                Y[0] = b*X[r*C]*X[r*C]; c = 1;
                for (n=r*C+1; n<r*C+C; n++, c++) { Y[c] = b*X[n]*X[n] + a*Y[c-1]; }
                cblas_dscal(tranz,0.0,&Y[0],1); //in case transient
                qsort(Y,(size_t)(C),sizeof(double),cmp_ascend_d);
                p90_rms = sqrt(Y[(int)(floor(rms_prctile*C))]);
                scale = target_RMS/p90_rms; scalemx = max_RMS/sqrt(Y[C-1]);
                scale = fmin(scale, scalemx);
                cblas_dscal(C,scale,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rms_scale_d: dim must be 0 or 1.\n"); return 1;
    }

    free(Y);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

