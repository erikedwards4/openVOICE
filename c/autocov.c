//The "unbiased" version uses N-l in the denominator instead of N.
//It is actually just "less biased", but is slower,
//has larger mean-squared error, and doesn't match FFT estimate.

//To test compile:
//gcc -c autocov.c -O2 -std=c99 -Wall -Wextra
//clang -c autocov.c -O2 -std=c99 -Weverything
//g++ -c autocov.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c autocov.c -O2 -std=c++11 -Weverything

#include "../include/openvoice.h"

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif


int autocov_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased)
{
    //const float o = 1.0f;
    //float m;
    int r, c, l;
    float *X1;  //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in autocov_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in autocov_s: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in autocov_s: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in autocov_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,&X1[0],1);
                //m = cblas_sdot(R,&X1[0],1,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[l+c*L] = cblas_sdot(R-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_sscal(C,(float)R/(R-l),&Y[l],L); } }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c],C,&X1[0],1);
                //m = cblas_sdot(R,&X1[0],1,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[c+l*C] = cblas_sdot(R-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_sscal(C,(float)R/(R-l),&Y[l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in autocov_s: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in autocov_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r],R,&X1[0],1);
                //m = cblas_sdot(C,&X1[0],1,&o,0) / C;
                //cblas_saxpy(C,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[r+l*R] = cblas_sdot(C-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_sscal(R,(float)C/(C-l),&Y[l*R],1); } }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,&X1[0],1);
                //m = cblas_sdot(C,&X1[0],1,&o,0) / C;
                //cblas_saxpy(C,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[l+r*L] = cblas_sdot(C-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_sscal(R,(float)C/(C-l),&Y[l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in autocov_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int autocov_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased)
{
    int r, c, l;
    double *X1;  //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in autocov_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in autocov_d: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in autocov_d: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in autocov_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,&X1[0],1);
                //m = cblas_ddot(R,&X1[0],1,&o,0) / R;
                //cblas_daxpy(R,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[l+c*L] = cblas_ddot(R-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_dscal(C,(double)R/(R-l),&Y[l],L); } }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c],C,&X1[0],1);
                //m = cblas_ddot(R,&X1[0],1,&o,0) / R;
                //cblas_daxpy(R,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[c+l*C] = cblas_ddot(R-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_dscal(C,(double)R/(R-l),&Y[l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in autocov_d: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in autocov_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r],R,&X1[0],1);
                //m = cblas_ddot(C,&X1[0],1,&o,0) / C;
                //cblas_daxpy(C,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[r+l*R] = cblas_ddot(C-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_dscal(R,(double)C/(C-l),&Y[l*R],1); } }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,&X1[0],1);
                //m = cblas_ddot(C,&X1[0],1,&o,0) / C;
                //cblas_daxpy(C,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[l+r*L] = cblas_ddot(C-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_dscal(R,(double)C/(C-l),&Y[l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in autocov_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int autocov_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased)
{
    _Complex float dot;
    int r, c, l;
    float *X1;  //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in autocov_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in autocov_c: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in autocov_c: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)(2*R)*sizeof(float)))) { fprintf(stderr,"error in autocov_c: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[2*c*R],1,&X1[0],1);
                //m = -cblas_cdotu(R,&X1[0],1,&o[0],0) / R;
                //cblas_caxpy(R,(float *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_cdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(l+c*L)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_csscal(C,(double)R/(R-l),&Y[2*l],L); } }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[2*c],2*C,&X1[0],1);
                //m = -cblas_cdotu(R,&X1[0],1,&o[0],0) / R;
                //cblas_caxpy(R,(float *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_cdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(c+l*C)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_csscal(C,(double)R/(R-l),&Y[2*l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in autocov_c: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)(2*C)*sizeof(float)))) { fprintf(stderr,"error in autocov_c: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[2*r],2*R,&X1[0],1);
                //m = -cblas_cdotu(C,&X1[0],1,&o[0],0) / C;
                //cblas_caxpy(C,(float *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_cdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(r+l*R)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_csscal(R,(float)C/(C-l),&Y[2*l*R],1); } }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[2*r*C],1,&X1[0],1);
                //m = -cblas_cdotu(C,&X1[0],1,&o[0],0) / C;
                //cblas_caxpy(C,(float *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_cdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(l+r*L)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_csscal(R,(float)C/(C-l),&Y[2*l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in autocov_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int autocov_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased)
{
    _Complex double dot;
    int r, c, l;
    double *X1;  //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in autocov_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in autocov_z: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in autocov_z: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)(2*R)*sizeof(double)))) { fprintf(stderr,"error in autocov_z: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[2*c*R],1,&X1[0],1);
                //m = -cblas_zdotu(R,&X1[0],1,&o[0],0) / R;
                //cblas_zaxpy(R,(double *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_zdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(l+c*L)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_zdscal(C,(double)R/(R-l),&Y[2*l],L); } }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[2*c],2*C,&X1[0],1);
                //m = -cblas_zdotu(R,&X1[0],1,&o[0],0) / R;
                //cblas_zaxpy(R,(double *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_zdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(c+l*C)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_zdscal(C,(double)R/(R-l),&Y[2*l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in autocov_z: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)(2*C)*sizeof(double)))) { fprintf(stderr,"error in autocov_z: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[2*r],2*R,&X1[0],1);
                //m = -cblas_zdotu(C,&X1[0],1,&o[0],0) / C;
                //cblas_zaxpy(C,(double *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_zdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(r+l*R)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_zdscal(R,(double)C/(C-l),&Y[2*l*R],1); } }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[2*r*C],1,&X1[0],1);
                //m = -cblas_zdotu(C,&X1[0],1,&o[0],0) / C;
                //cblas_zaxpy(C,(double *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_zdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(l+r*L)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_zdscal(R,(double)C/(C-l),&Y[2*l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in autocov_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

