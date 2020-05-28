//Gets reflection coefficients (RCs) from autoregressive (AR) parameters along rows or cols of X.
//Input (X) and output (Y) have the same size.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int ar2rc_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2rc_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2rc_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2rc_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);


int ar2rc_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    float sc;
    float *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2rc_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2rc_s: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(float *)malloc((size_t)(P)*sizeof(float)))) { fprintf(stderr,"error in ar2rc_s: problem with malloc. "); perror("malloc"); return 1; }
    cblas_scopy(R*C,X,1,Y,1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_scopy(p+1,&Y[c*R],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c*R+q] += sc * y[p-1-q]; }
                    cblas_sscal(p,1.0f/fmaf(sc,-sc,1.0f),&Y[c*R],1);
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_scopy(p+1,&Y[c],C,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c+q*C] += sc * y[p-1-q]; }
                    cblas_sscal(p,1.0f/fmaf(sc,-sc,1.0f),&Y[c],C);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_scopy(p+1,&Y[r],R,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r+q*R] += sc * y[p-1-q]; }
                    cblas_sscal(p,1.0f/fmaf(sc,-sc,1.0f),&Y[r],R);
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_scopy(p+1,&Y[r*C],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r*C+q] += sc * y[p-1-q]; }
                    cblas_sscal(p,1.0f/fmaf(sc,-sc,1.0f),&Y[r*C],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2rc_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int ar2rc_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    double sc;
    double *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2rc_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2rc_d: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(double *)malloc((size_t)(P)*sizeof(double)))) { fprintf(stderr,"error in ar2rc_d: problem with malloc. "); perror("malloc"); return 1; }
    cblas_dcopy(R*C,X,1,Y,1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_dcopy(p+1,&Y[c*R],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c*R+q] += sc * y[p-1-q]; }
                    cblas_dscal(p,1.0/fma(sc,-sc,1.0),&Y[c*R],1);
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_dcopy(p+1,&Y[c],C,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c+q*C] += sc * y[p-1-q]; }
                    cblas_dscal(p,1.0/fma(sc,-sc,1.0),&Y[c],C);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_dcopy(p+1,&Y[r],R,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r+q*R] += sc * y[p-1-q]; }
                    cblas_dscal(p,1.0/fma(sc,-sc,1.0),&Y[r],R);
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_dcopy(p+1,&Y[r*C],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r*C+q] += sc * y[p-1-q]; }
                    cblas_dscal(p,1.0/fma(sc,-sc,1.0),&Y[r*C],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2rc_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int ar2rc_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    float sc[2] = {0.0f,0.0f};
    float *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2rc_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2rc_c: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(float *)malloc((size_t)(2*P)*sizeof(float)))) { fprintf(stderr,"error in ar2rc_c: problem with malloc. "); perror("malloc"); return 1; }
    cblas_ccopy(R*C,X,1,Y,1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_ccopy(p+1,&Y[2*c*R],1,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c*R+q)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c*R+q)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_csscal(p,1.0f/(1.0f-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*c*R],1);
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_ccopy(p+1,&Y[2*c],C,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c+q*C)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c+q*C)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_csscal(p,1.0f/(1.0f-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*c],C);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_ccopy(p+1,&Y[2*r],R,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r+q*R)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r+q*R)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_csscal(p,1.0f/(1.0f-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*r],R);
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_ccopy(p+1,&Y[2*r*C],1,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r*C+q)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r*C+q)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_csscal(p,1.0f/(1.0f-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*r*C],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2rc_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int ar2rc_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    double sc[2] = {0.0,0.0};
    double *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2rc_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2rc_z: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(double *)malloc((size_t)(2*P)*sizeof(double)))) { fprintf(stderr,"error in ar2rc_z: problem with malloc. "); perror("malloc"); return 1; }
    cblas_zcopy(R*C,X,1,Y,1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_zcopy(p+1,&Y[2*c*R],1,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c*R+q)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c*R+q)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_zdscal(p,1.0/(1.0-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*c*R],1);
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_zcopy(p+1,&Y[2*c],C,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c+q*C)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c+q*C)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_zdscal(p,1.0/(1.0-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*c],C);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_zcopy(p+1,&Y[2*r],R,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r+q*R)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r+q*R)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_zdscal(p,1.0/(1.0-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*r],R);
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (p=P-1; p>0; p--)
                {
                    cblas_zcopy(p+1,&Y[2*r*C],1,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r*C+q)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r*C+q)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_zdscal(p,1.0/(1.0-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*r*C],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2rc_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


#ifdef __cplusplus
}
}
#endif

