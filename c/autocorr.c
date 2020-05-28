//Autocorrelation of rows or cols of X according to dim.
//This does not subtract the row or col means (see commented code below).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int autocorr_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L);
int autocorr_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L);
int autocorr_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L);
int autocorr_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L);


int autocorr_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L)
{
    const float o = 1.0f;
    float s;
    int r, c, l;
    float *X1;  //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in autocorr_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in autocorr_s: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in autocorr_s: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in autocorr_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            cblas_scopy(C,&o,0,&Y[0],L); //lag 0
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,&X1[0],1);
                //m = cblas_sdot(R,&X1[0],1,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X1[0],1);
                s = 1.0f / cblas_sdot(R,&X1[0],1,&X1[0],1);
                for (l=1; l<L; l++) { Y[l+c*L] = s * cblas_sdot(R-l,&X1[0],1,&X1[l],1); }
            }
        }
        else
        {
            cblas_scopy(C,&o,0,&Y[0],1); //lag 0
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c],C,&X1[0],1);
                //m = cblas_sdot(R,&X1[0],1,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X1[0],1);
                s = 1.0f / cblas_sdot(R,&X1[0],1,&X1[0],1);
                for (l=1; l<L; l++) { Y[c+l*C] = s * cblas_sdot(R-l,&X1[0],1,&X1[l],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in autocorr_s: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in autocorr_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            cblas_scopy(R,&o,0,&Y[0],1); //lag 0
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r],R,&X1[0],1);
                //m = cblas_sdot(C,&X1[0],1,&o,0) / C;
                //cblas_saxpy(C,-m,&o,0,&X1[0],1);
                s = 1.0f / cblas_sdot(C,&X1[0],1,&X1[0],1);
                for (l=1; l<L; l++) { Y[r+l*R] = s * cblas_sdot(C-l,&X1[0],1,&X1[l],1); }
            }
        }
        else
        {
            cblas_scopy(R,&o,0,&Y[0],L); //lag 0
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,&X1[0],1);
                //m = cblas_sdot(C,&X1[0],1,&o,0) / C;
                //cblas_saxpy(C,-m,&o,0,&X1[0],1);
                s = 1.0f / cblas_sdot(C,&X1[0],1,&X1[0],1);
                for (l=1; l<L; l++) { Y[l+r*L] = s * cblas_sdot(C-l,&X1[0],1,&X1[l],1); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in autocorr_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int autocorr_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L)
{
    const double o = 1.0;
    double s;
    int r, c, l;
    double *X1;  //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in autocorr_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in autocorr_d: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in autocorr_d: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in autocorr_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            cblas_dcopy(C,&o,0,&Y[0],L); //lag 0
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,&X1[0],1);
                //m = cblas_ddot(R,&X1[0],1,&o,0) / R;
                //cblas_daxpy(R,-m,&o,0,&X1[0],1);
                s = 1.0 / cblas_ddot(R,&X1[0],1,&X1[0],1);
                for (l=1; l<L; l++) { Y[l+c*L] = s * cblas_ddot(R-l,&X1[0],1,&X1[l],1); }
            }
        }
        else
        {
            cblas_dcopy(C,&o,0,&Y[0],1); //lag 0
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c],C,&X1[0],1);
                //m = cblas_ddot(R,&X1[0],1,&o,0) / R;
                //cblas_daxpy(R,-m,&o,0,&X1[0],1);
                s = 1.0 / cblas_ddot(R,&X1[0],1,&X1[0],1);
                for (l=1; l<L; l++) { Y[c+l*C] = s * cblas_ddot(R-l,&X1[0],1,&X1[l],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in autocorr_d: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in autocorr_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            cblas_dcopy(R,&o,0,&Y[0],1); //lag 0
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r],R,&X1[0],1);
                //m = cblas_ddot(C,&X1[0],1,&o,0) / C;
                //cblas_daxpy(C,-m,&o,0,&X1[0],1);
                s = 1.0 / cblas_ddot(C,&X1[0],1,&X1[0],1);
                for (l=1; l<L; l++) { Y[r+l*R] = s * cblas_ddot(C-l,&X1[0],1,&X1[l],1); }
            }
        }
        else
        {
            cblas_dcopy(R,&o,0,&Y[0],L); //lag 0
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,&X1[0],1);
                //m = cblas_ddot(C,&X1[0],1,&o,0) / C;
                //cblas_daxpy(C,-m,&o,0,&X1[0],1);
                s = 1.0 / cblas_ddot(C,&X1[0],1,&X1[0],1);
                for (l=1; l<L; l++) { Y[l+r*L] = s * cblas_ddot(C-l,&X1[0],1,&X1[l],1); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in autocorr_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int autocorr_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L)
{
    const float o[2] =  {1.0f,0.0f};
    _Complex float dot;
    float s;
    int r, c, l;
    float *X1;  //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in autocorr_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in autocorr_c: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in autocorr_c: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)(2*R)*sizeof(float)))) { fprintf(stderr,"error in autocorr_c: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            cblas_ccopy(C,&o[0],0,&Y[0],L); //lag 0
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[2*c*R],1,&X1[0],1);
                //m = -cblas_cdotu(R,&X1[0],1,&o[0],0) / R;
                //cblas_caxpy(R,(float *)&m,&o[0],0,&X1[0],1);
                dot = cblas_cdotu(R,&X1[0],1,&X1[0],1);
                memcpy(&s,(float *)&dot,sizeof(float));
                s = 1.0f / s;
                for (l=1; l<L; l++)
                {
                    dot = cblas_cdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_csscal(1,s,(float *)&dot,1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(l+c*L)],1);
                }
            }
        }
        else
        {
            cblas_ccopy(C,&o[0],0,&Y[0],1); //lag 0
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[2*c],2*C,&X1[0],1);
                //m = -cblas_cdotu(R,&X1[0],1,&o[0],0) / R;
                //cblas_caxpy(R,(float *)&m,&o[0],0,&X1[0],1);
                dot = cblas_cdotu(R,&X1[0],1,&X1[0],1);
                memcpy(&s,(float *)&dot,sizeof(float));
                s = 1.0f / s;
                for (l=1; l<L; l++)
                {
                    dot = cblas_cdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_csscal(1,s,(float *)&dot,1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(c+l*C)],1);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in autocorr_c: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)(2*C)*sizeof(float)))) { fprintf(stderr,"error in autocorr_c: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            cblas_ccopy(R,&o[0],0,&Y[0],1); //lag 0
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[2*r],2*R,&X1[0],1);
                //m = -cblas_cdotu(C,&X1[0],1,&o[0],0) / C;
                //cblas_caxpy(C,(float *)&m,&o[0],0,&X1[0],1);
                dot = cblas_cdotu(C,&X1[0],1,&X1[0],1);
                memcpy(&s,(float *)&dot,sizeof(float));
                s = 1.0f / s;
                for (l=1; l<L; l++)
                {
                    dot = cblas_cdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_csscal(1,s,(float *)&dot,1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(r+l*R)],1);
                }
            }
        }
        else
        {
            cblas_ccopy(R,&o[0],0,&Y[0],L); //lag 0
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[2*r*C],1,&X1[0],1);
                //m = -cblas_cdotu(C,&X1[0],1,&o[0],0) / C;
                //cblas_caxpy(C,(float *)&m,&o[0],0,&X1[0],1);
                dot = cblas_cdotu(C,&X1[0],1,&X1[0],1);
                memcpy(&s,(float *)&dot,sizeof(float));
                s = 1.0f / s;
                for (l=1; l<L; l++)
                {
                    dot = cblas_cdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_csscal(1,s,(float *)&dot,1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(l+r*L)],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in autocorr_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int autocorr_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L)
{
    const double o[2] =  {1.0,0.0};
    _Complex double dot;
    double s;
    double *X1;  //1 row or col of X
    int r, c, l;

    //Checks
    if (R<1) { fprintf(stderr,"error in autocorr_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in autocorr_z: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in autocorr_z: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)(2*R)*sizeof(double)))) { fprintf(stderr,"error in autocorr_z: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            cblas_zcopy(C,&o[0],0,&Y[0],L); //lag 0
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[2*c*R],1,&X1[0],1);
                //m = -cblas_zdotu(R,&X1[0],1,&o[0],0) / R;
                //cblas_zaxpy(R,(double *)&m,&o[0],0,&X1[0],1);
                dot = cblas_zdotu(R,&X1[0],1,&X1[0],1);
                memcpy(&s,(double *)&dot,sizeof(double));
                s = 1.0 / s;
                for (l=1; l<L; l++)
                {
                    dot = cblas_zdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_zdscal(1,s,(double *)&dot,1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(l+c*L)],1);
                }
            }
        }
        else
        {
            cblas_zcopy(C,&o[0],0,&Y[0],1); //lag 0
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[2*c],2*C,&X1[0],1);
                //m = -cblas_zdotu(R,&X1[0],1,&o[0],0) / R;
                //cblas_zaxpy(R,(double *)&m,&o[0],0,&X1[0],1);
                dot = cblas_zdotu(R,&X1[0],1,&X1[0],1);
                memcpy(&s,(double *)&dot,sizeof(double));
                s = 1.0 / s;
                for (l=1; l<L; l++)
                {
                    dot = cblas_zdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_zdscal(1,s,(double *)&dot,1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(c+l*C)],1);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in autocorr_z: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)(2*C)*sizeof(double)))) { fprintf(stderr,"error in autocorr_z: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            cblas_zcopy(R,&o[0],0,&Y[0],1); //lag 0
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[2*r],2*R,&X1[0],1);
                //m = -cblas_zdotu(C,&X1[0],1,&o[0],0) / C;
                //cblas_zaxpy(C,(double *)&m,&o[0],0,&X1[0],1);
                dot = cblas_zdotu(C,&X1[0],1,&X1[0],1);
                memcpy(&s,(double *)&dot,sizeof(double));
                s = 1.0 / s;
                for (l=1; l<L; l++)
                {
                    dot = cblas_zdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_zdscal(1,s,(double *)&dot,1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(r+l*R)],1);
                }
            }
        }
        else
        {
            cblas_zcopy(R,&o[0],0,&Y[0],L); //lag 0
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[2*r*C],1,&X1[0],1);
                //m = -cblas_zdotu(C,&X1[0],1,&o[0],0) / C;
                //cblas_zaxpy(C,(double *)&m,&o[0],0,&X1[0],1);
                dot = cblas_zdotu(C,&X1[0],1,&X1[0],1);
                memcpy(&s,(double *)&dot,sizeof(double));
                s = 1.0 / s;
                for (l=1; l<L; l++)
                {
                    dot = cblas_zdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_zdscal(1,s,(double *)&dot,1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(l+r*L)],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in autocorr_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

