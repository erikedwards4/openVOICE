//Gets median of each row or col of X according to dim.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include "cmp_ascend.c"

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int median_s (float *X, const char iscolmajor, const int R, const int C, const int dim, float *Y);
int median_d (double *X, const char iscolmajor, const int R, const int C, const int dim, double *Y);
int median_c (float *X, const char iscolmajor, const int R, const int C, const int dim, float *Y);
int median_z (double *X, const char iscolmajor, const int R, const int C, const int dim, double *Y);


int median_s (float *X, const char iscolmajor, const int R, const int C, const int dim, float *Y)
{
    int r, c;
    float *X1; //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in median_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in median_s: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (!(X1=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in median_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,X,1,X1,1);
                qsort(X1,(size_t)(R),sizeof(float),cmp_ascend_s);
                Y[c] = (R%2) ? X1[R/2] : 0.5f*(X1[R/2]+X1[R/2-1]);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,X,C,X1,1);
                qsort(X1,(size_t)(R),sizeof(float),cmp_ascend_s);
                Y[c] = (R%2) ? X1[R/2] : 0.5f*(X1[R/2]+X1[R/2-1]);
            }
        }
    }
    else if (dim==1)
    {
        if (!(X1=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in median_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,X,R,X1,1);
                qsort(X1,(size_t)(C),sizeof(float),cmp_ascend_s);
                Y[r] = (C%2) ? X1[C/2] : 0.5f*(X1[C/2]+X1[C/2-1]);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,X,1,X1,1);
                qsort(X1,(size_t)(C),sizeof(float),cmp_ascend_s);
                Y[r] = (C%2) ? X1[C/2] : 0.5f*(X1[C/2]+X1[C/2-1]);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in median_s: dim must be 0 or 1.\n"); return 1;
    }

    free(X1);
    return 0;
}


int median_d (double *X, const char iscolmajor, const int R, const int C, const int dim, double *Y)
{
    int r, c;
    double *X1; //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in median_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in median_d: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (!(X1=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in median_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,X,1,X1,1);
                qsort(X1,(size_t)(R),sizeof(double),cmp_ascend_s);
                Y[c] = (R%2) ? X1[R/2] : 0.5*(X1[R/2]+X1[R/2-1]);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,X,C,X1,1);
                qsort(X1,(size_t)(R),sizeof(double),cmp_ascend_s);
                Y[c] = (R%2) ? X1[R/2] : 0.5*(X1[R/2]+X1[R/2-1]);
            }
        }
    }
    else if (dim==1)
    {
        if (!(X1=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in median_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,X,R,X1,1);
                qsort(X1,(size_t)(C),sizeof(double),cmp_ascend_s);
                Y[r] = (C%2) ? X1[C/2] : 0.5*(X1[C/2]+X1[C/2-1]);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,X,1,X1,1);
                qsort(X1,(size_t)(C),sizeof(double),cmp_ascend_s);
                Y[r] = (C%2) ? X1[C/2] : 0.5*(X1[C/2]+X1[C/2-1]);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in median_d: dim must be 0 or 1.\n"); return 1;
    }

    free(X1);
    return 0;
}


int median_c (float *X, const char iscolmajor, const int R, const int C, const int dim, float *Y)
{
    int r, c;
    float *X1; //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in median_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in median_c: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (!(X1=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in median_c: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[0],2,X1,1);
                qsort(X1,(size_t)(R),sizeof(float),cmp_ascend_s);
                Y[2*c] = (R%2) ? X1[R/2] : 0.5f*(X1[R/2]+X1[R/2-1]);
                cblas_ccopy(R,&X[1],2,X1,1);
                qsort(X1,(size_t)(R),sizeof(float),cmp_ascend_s);
                Y[2*c+1] = (R%2) ? X1[R/2] : 0.5f*(X1[R/2]+X1[R/2-1]);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[0],2*C,X1,1);
                qsort(X1,(size_t)(R),sizeof(float),cmp_ascend_s);
                Y[2*c] = (R%2) ? X1[R/2] : 0.5f*(X1[R/2]+X1[R/2-1]);
                cblas_ccopy(R,&X[1],2*C,X1,1);
                qsort(X1,(size_t)(R),sizeof(float),cmp_ascend_s);
                Y[2*c+1] = (R%2) ? X1[R/2] : 0.5f*(X1[R/2]+X1[R/2-1]);
            }
        }
    }
    else if (dim==1)
    {
        if (!(X1=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in median_c: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[0],2*R,X1,1);
                qsort(X1,(size_t)(C),sizeof(float),cmp_ascend_s);
                Y[2*r] = (C%2) ? X1[C/2] : 0.5f*(X1[C/2]+X1[C/2-1]);
                cblas_ccopy(C,&X[1],2*R,X1,1);
                qsort(X1,(size_t)(C),sizeof(float),cmp_ascend_s);
                Y[2*r+1] = (C%2) ? X1[C/2] : 0.5f*(X1[C/2]+X1[C/2-1]);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[0],1,X1,1);
                qsort(X1,(size_t)(C),sizeof(float),cmp_ascend_s);
                Y[2*r] = (C%2) ? X1[C/2] : 0.5f*(X1[C/2]+X1[C/2-1]);
                cblas_ccopy(C,&X[1],1,X1,1);
                qsort(X1,(size_t)(C),sizeof(float),cmp_ascend_s);
                Y[2*r+1] = (C%2) ? X1[C/2] : 0.5f*(X1[C/2]+X1[C/2-1]);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in median_c: dim must be 0 or 1.\n"); return 1;
    }

    free(X1);
    return 0;
}


int median_z (double *X, const char iscolmajor, const int R, const int C, const int dim, double *Y)
{
    int r, c;
    double *X1; //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in median_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in median_z: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (!(X1=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in median_z: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[0],2,X1,1);
                qsort(X1,(size_t)(R),sizeof(double),cmp_ascend_s);
                Y[2*c] = (R%2) ? X1[R/2] : 0.5*(X1[R/2]+X1[R/2-1]);
                cblas_zcopy(R,&X[1],2,X1,1);
                qsort(X1,(size_t)(R),sizeof(double),cmp_ascend_s);
                Y[2*c+1] = (R%2) ? X1[R/2] : 0.5*(X1[R/2]+X1[R/2-1]);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[0],2*C,X1,1);
                qsort(X1,(size_t)(R),sizeof(double),cmp_ascend_s);
                Y[2*c] = (R%2) ? X1[R/2] : 0.5*(X1[R/2]+X1[R/2-1]);
                cblas_zcopy(R,&X[1],2*C,X1,1);
                qsort(X1,(size_t)(R),sizeof(double),cmp_ascend_s);
                Y[2*c+1] = (R%2) ? X1[R/2] : 0.5*(X1[R/2]+X1[R/2-1]);
            }
        }
    }
    else if (dim==1)
    {
        if (!(X1=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in median_z: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[0],2*R,X1,1);
                qsort(X1,(size_t)(C),sizeof(double),cmp_ascend_s);
                Y[2*r] = (C%2) ? X1[C/2] : 0.5*(X1[C/2]+X1[C/2-1]);
                cblas_zcopy(C,&X[1],2*R,X1,1);
                qsort(X1,(size_t)(C),sizeof(double),cmp_ascend_s);
                Y[2*r+1] = (C%2) ? X1[C/2] : 0.5*(X1[C/2]+X1[C/2-1]);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[0],1,X1,1);
                qsort(X1,(size_t)(C),sizeof(double),cmp_ascend_s);
                Y[2*r] = (C%2) ? X1[C/2] : 0.5*(X1[C/2]+X1[C/2-1]);
                cblas_zcopy(C,&X[1],1,X1,1);
                qsort(X1,(size_t)(C),sizeof(double),cmp_ascend_s);
                Y[2*r+1] = (C%2) ? X1[C/2] : 0.5*(X1[C/2]+X1[C/2-1]);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in median_z: dim must be 0 or 1.\n"); return 1;
    }

    free(X1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

