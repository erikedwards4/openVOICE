//Gets power spectral densities (PSDs) from autoregressive (AR) parameters along rows or cols of X.
//The 2nd input is the vector V of variances (prediction errors) for each row or col of X.
//The 3rd input is a vector W of F freqs (in radians) at which to get the PSD.

//Following convention of Octave signal package ar_psd.m, I double the power for real-valued X.
//See Eq. (2.38) of Kay and Marple [1981].

//To test compile:
//gcc -c ar2psd.c -O2 -std=c99 -Wall -Wextra
//clang -c ar2psd.c -O2 -std=c99 -Weverything
//g++ -c ar2psd.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c ar2psd.c -O2 -std=c++11 -Weverything -Wno-old-style-cast

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int ar2psd_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *V, const float *W, const int F, const int dim)
{
    const int P = (dim==0) ? R : C;
    int r, c, p, f, n;
    float *Er, *Ei, *yr, *yi;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2psd_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2psd_s: ncols X must be positive\n"); return 1; }
    if (F<1) { fprintf(stderr,"error in ar2psd_s: F (length W) must be positive\n"); return 1; }

    //Allocate
    if (!(yr=(float *)malloc((size_t)(F)*sizeof(float)))) { fprintf(stderr,"error in ar2psd_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(yi=(float *)malloc((size_t)(F)*sizeof(float)))) { fprintf(stderr,"error in ar2psd_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Er=(float *)malloc((size_t)(F*P)*sizeof(float)))) { fprintf(stderr,"error in ar2psd_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Ei=(float *)malloc((size_t)(F*P)*sizeof(float)))) { fprintf(stderr,"error in ar2psd_s: problem with malloc. "); perror("malloc"); return 1; }

    //Make complex-valued E matrix
    for (n=0; n<F*P; n++)
    {
        if (iscolmajor) { f = n%F; p = n/F; } else { f = n/P; p = n%P; }
        Er[n] = -cosf(W[f]*(p+1)); Ei[n] = sinf(W[f]*(p+1));
    }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_sgemv(CblasColMajor,CblasNoTrans,F,P,1.0f,Er,F,&X[c*R],1,0.0f,yr,1);
                cblas_sgemv(CblasColMajor,CblasNoTrans,F,P,1.0f,Ei,F,&X[c*R],1,0.0f,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0f; Y[c*F+f] = 2.0f*V[c]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_sgemv(CblasRowMajor,CblasNoTrans,F,P,1.0f,Er,P,&X[c],C,0.0f,yr,1);
                cblas_sgemv(CblasRowMajor,CblasNoTrans,F,P,1.0f,Ei,P,&X[c],C,0.0f,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0f; Y[c+f*C] = 2.0f*V[c]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_sgemv(CblasColMajor,CblasNoTrans,F,P,1.0f,Er,F,&X[r],R,0.0f,yr,1);
                cblas_sgemv(CblasColMajor,CblasNoTrans,F,P,1.0f,Ei,F,&X[r],R,0.0f,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0f; Y[r+f*R] = 2.0f*V[r]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_sgemv(CblasRowMajor,CblasNoTrans,F,P,1.0f,Er,P,&X[r*C],1,0.0f,yr,1);
                cblas_sgemv(CblasRowMajor,CblasNoTrans,F,P,1.0f,Ei,P,&X[r*C],1,0.0f,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0f; Y[r*F+f] = 2.0f*V[r]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2psd_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(Er); free(Ei); free(yr); free(yi);
    return 0;
}


int ar2psd_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *V, const double *W, const int F, const int dim)
{
    const int P = (dim==0) ? R : C;
    int r, c, p, f, n;
    double *Er, *Ei, *yr, *yi;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2psd_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2psd_d: ncols X must be positive\n"); return 1; }
    if (F<1) { fprintf(stderr,"error in ar2psd_d: F (length W) must be positive\n"); return 1; }

    //Allocate
    if (!(yr=(double *)malloc((size_t)(F)*sizeof(double)))) { fprintf(stderr,"error in ar2psd_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(yi=(double *)malloc((size_t)(F)*sizeof(double)))) { fprintf(stderr,"error in ar2psd_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Er=(double *)malloc((size_t)(F*P)*sizeof(double)))) { fprintf(stderr,"error in ar2psd_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Ei=(double *)malloc((size_t)(F*P)*sizeof(double)))) { fprintf(stderr,"error in ar2psd_d: problem with malloc. "); perror("malloc"); return 1; }

    //Make complex-valued E matrix
    for (n=0; n<F*P; n++)
    {
        if (iscolmajor) { f = n%F; p = n/F; } else { f = n/P; p = n%P; }
        Er[n] = -cos(W[f]*(p+1)); Ei[n] = sin(W[f]*(p+1));
    }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dgemv(CblasColMajor,CblasNoTrans,F,P,1.0,Er,F,&X[c*R],1,0.0,yr,1);
                cblas_dgemv(CblasColMajor,CblasNoTrans,F,P,1.0,Ei,F,&X[c*R],1,0.0,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0; Y[c*F+f] = 2.0*V[c]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dgemv(CblasRowMajor,CblasNoTrans,F,P,1.0,Er,P,&X[c],C,0.0,yr,1);
                cblas_dgemv(CblasRowMajor,CblasNoTrans,F,P,1.0,Ei,P,&X[c],C,0.0,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0; Y[c+f*C] = 2.0*V[c]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_dgemv(CblasColMajor,CblasNoTrans,F,P,1.0,Er,F,&X[r],R,0.0,yr,1);
                cblas_dgemv(CblasColMajor,CblasNoTrans,F,P,1.0,Ei,F,&X[r],R,0.0,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0; Y[r+f*R] = 2.0*V[r]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dgemv(CblasRowMajor,CblasNoTrans,F,P,1.0,Er,P,&X[r*C],1,0.0,yr,1);
                cblas_dgemv(CblasRowMajor,CblasNoTrans,F,P,1.0,Ei,P,&X[r*C],1,0.0,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0; Y[r*F+f] = 2.0*V[r]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2psd_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(Er); free(Ei); free(yr); free(yi);
    return 0;
}


int ar2psd_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *V, const float *W, const int F, const int dim)
{
    const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
    const int P = (dim==0) ? R : C;
    int r, c, p, f, n;
    float *E, *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2psd_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2psd_c: ncols X must be positive\n"); return 1; }
    if (F<1) { fprintf(stderr,"error in ar2psd_c: F (length W) must be positive\n"); return 1; }

    //Allocate
    if (!(y=(float *)malloc((size_t)(2*F)*sizeof(float)))) { fprintf(stderr,"error in ar2psd_c: problem with malloc. "); perror("malloc"); return 1; }
    if (!(E=(float *)malloc((size_t)(2*F*P)*sizeof(float)))) { fprintf(stderr,"error in ar2psd_c: problem with malloc. "); perror("malloc"); return 1; }

    //Make complex-valued E matrix
    for (n=0; n<F*P; n++)
    {
        if (iscolmajor) { f = n%F; p = n/F; } else { f = n/P; p = n%P; }
        E[2*n] = -cosf(W[f]*(p+1)); E[2*n+1] = sinf(W[f]*(p+1));
    }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_cgemv(CblasColMajor,CblasNoTrans,F,P,o,E,F,&X[2*c*R],1,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0f; Y[c*F+f] = V[c] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_cgemv(CblasRowMajor,CblasNoTrans,F,P,o,E,P,&X[2*c],C,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0f; Y[c+f*C] = V[c] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_cgemv(CblasColMajor,CblasNoTrans,F,P,o,E,F,&X[2*r],R,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0f; Y[r+f*R] = V[r] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_cgemv(CblasRowMajor,CblasNoTrans,F,P,o,E,P,&X[2*r*C],1,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0f; Y[r*F+f] = V[r] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2psd_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(E); free(y);
    return 0;
}


int ar2psd_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *V, const double *W, const int F, const int dim)
{
    const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
    const int P = (dim==0) ? R : C;
    int r, c, p, f, n;
    double *E, *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2psd_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2psd_z: ncols X must be positive\n"); return 1; }
    if (F<1) { fprintf(stderr,"error in ar2psd_z: F (length W) must be positive\n"); return 1; }

    //Allocate
    if (!(y=(double *)malloc((size_t)(2*F)*sizeof(double)))) { fprintf(stderr,"error in ar2psd_z: problem with malloc. "); perror("malloc"); return 1; }
    if (!(E=(double *)malloc((size_t)(2*F*P)*sizeof(double)))) { fprintf(stderr,"error in ar2psd_z: problem with malloc. "); perror("malloc"); return 1; }

    //Make complex-valued E matrix
    for (n=0; n<F*P; n++)
    {
        if (iscolmajor) { f = n%F; p = n/F; } else { f = n/P; p = n%P; }
        E[2*n] = -cos(W[f]*(p+1)); E[2*n+1] = sin(W[f]*(p+1));
    }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_zgemv(CblasColMajor,CblasNoTrans,F,P,o,E,F,&X[2*c*R],1,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0; Y[c*F+f] = V[c] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_zgemv(CblasRowMajor,CblasNoTrans,F,P,o,E,P,&X[2*c],C,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0; Y[c+f*C] = V[c] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_zgemv(CblasColMajor,CblasNoTrans,F,P,o,E,F,&X[2*r],R,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0; Y[r+f*R] = V[r] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_zgemv(CblasRowMajor,CblasNoTrans,F,P,o,E,P,&X[2*r*C],1,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0; Y[r*F+f] = V[r] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2psd_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(E); free(y);
    return 0;
}


#ifdef __cplusplus
}
#endif

