//Gets line spectral frequencies (LSFs).
//from the autoregressive (AR) coefficients along cols or rows of X.

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include "cmp_ascend.c"

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int ar2lsf_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2lsf_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2lsf_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2lsf_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);


int ar2lsf_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const float z = 0.0f, o = 1.0f;
    const int P = (dim==0) ? R : C;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P+1, n = P+1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    float *poly, *prev, *omegas, *compan, *wr, *wi, zz[1];
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2lsf_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2lsf_s: ncols X must be positive\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in ar2lsf_s: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    //Allocate
    if (!(poly=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(prev=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(omegas=(float *)malloc((size_t)(2*n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wr=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wi=(float *)malloc((size_t)(n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(float *)malloc((size_t)(n*n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_s: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        for (c=0; c<C; c++)
        {
            if (iscolmajor) { cblas_scopy(n,&X[c*R],1,poly,1); }
            else { cblas_scopy(n,&X[c],C,poly,1); }
            poly[P] = z; prev[P] = -o;
            for (r=0; r<P; r++) { prev[r] = poly[P-r-1]; }
            for (r=0; r<n; r++) { poly[r] -= prev[r]; }
            cblas_scopy(n*n,&z,0,compan,1); cblas_scopy(P,&o,0,&compan[1],n+1); cblas_scopy(n,poly,1,compan,n);
            info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in ar2lsf_s: lapacke decomposition failed\n"); return 1; }
            for (r=0; r<n; r++) { omegas[r] = (wi[r]>10.0f*FLT_EPSILON) ? atan2f(wi[r],wr[r]) : 0.0f; poly[r] += prev[r] + prev[r]; }
            cblas_scopy(n*n,&z,0,compan,1); cblas_scopy(P,&o,0,&compan[1],n+1); cblas_scopy(n,poly,1,compan,n);
            info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in ar2lsf_s: lapacke decomposition failed\n"); return 1; }
            for (r=0; r<n; r++) { omegas[n+r] = (wi[r]>10.0f*FLT_EPSILON) ? atan2f(wi[r],wr[r]) : 0.0f; }
            qsort(omegas,(size_t)(2*n),sizeof(float),cmp_ascend_s);
            if (iscolmajor) { cblas_scopy(P,&omegas[n+1],1,&Y[c*P],1); }
            else { cblas_scopy(P,&omegas[n+1],1,&Y[c],C); }
        }
    }
    else if (dim==1)
    {
        for (r=0; r<R; r++)
        {
            if (iscolmajor) { cblas_scopy(n,&X[r],R,poly,1); }
            else { cblas_scopy(n,&X[r*C],1,poly,1); }
            poly[P] = z; prev[P] = -o;
            for (c=0; c<P; c++) { prev[c] = poly[P-c-1]; }
            for (c=0; c<n; c++) { poly[c] -= prev[c]; }
            cblas_scopy(n*n,&z,0,compan,1); cblas_scopy(P,&o,0,&compan[1],n+1); cblas_scopy(n,poly,1,compan,n);
            info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in ar2lsf_s: lapacke decomposition failed\n"); return 1; }
            for (c=0; c<n; c++) { omegas[c] = (wi[c]>10.0f*FLT_EPSILON) ? atan2f(wi[c],wr[c]) : 0.0f; poly[c] += prev[c] + prev[c]; }
            cblas_scopy(n*n,&z,0,compan,1); cblas_scopy(P,&o,0,&compan[1],n+1); cblas_scopy(n,poly,1,compan,n);
            info = LAPACKE_shseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in ar2lsf_s: lapacke decomposition failed\n"); return 1; }
            for (c=0; c<n; c++) { omegas[n+c] = (wi[c]>10.0f*FLT_EPSILON) ? atan2f(wi[c],wr[c]) : 0.0f; }
            qsort(omegas,(size_t)(2*n),sizeof(float),cmp_ascend_s);
            if (iscolmajor) { cblas_scopy(P,&omegas[n+1],1,&Y[r],R); }
            else { cblas_scopy(P,&omegas[n+1],1,&Y[r*P],1); }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2lsf_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(compan); free(wr); free(wi); free(poly); free(prev); free(omegas);
    return 0;
}


int ar2lsf_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const double z = 0.0, o = 1.0;
    const int P = (dim==0) ? R : C;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P+1, n = P+1, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    double *poly, *prev, *omegas, *compan, *wr, *wi, zz[1];
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2lsf_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2lsf_d: ncols X must be positive\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in ar2lsf_d: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    //Allocate
    if (!(poly=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(prev=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(omegas=(double *)malloc((size_t)(2*n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wr=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wi=(double *)malloc((size_t)(n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(double *)malloc((size_t)(n*n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_d: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        for (c=0; c<C; c++)
        {
            if (iscolmajor) { cblas_dcopy(n,&X[c*R],1,poly,1); }
            else { cblas_dcopy(n,&X[c],C,poly,1); }
            poly[P] = z; prev[P] = -o;
            for (r=0; r<P; r++) { prev[r] = poly[n-r]; }
            for (r=0; r<n; r++) { poly[r] -= prev[r]; }
            cblas_dcopy(n*n,&z,0,compan,1); cblas_dcopy(P,&o,0,&compan[1],n+1); cblas_dcopy(n,poly,1,compan,n);
            info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in ar2lsf_d: lapacke decomposition failed\n"); return 1; }
            for (r=0; r<n; r++) { omegas[r] = (wi[r]>10.0*DBL_EPSILON) ? atan2(wi[r],wr[r]) : 0.0; poly[r] += prev[r] + prev[r]; }
            cblas_dcopy(n*n,&z,0,compan,1); cblas_dcopy(P,&o,0,&compan[1],n+1); cblas_dcopy(n,poly,1,compan,n);
            info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in ar2lsf_d: lapacke decomposition failed\n"); return 1; }
            for (r=0; r<n; r++) { omegas[n+r] = (wi[r]>10.0*DBL_EPSILON) ? atan2(wi[r],wr[r]) : 0.0; }
            qsort(omegas,(size_t)(2*n),sizeof(double),cmp_ascend_d);
            if (iscolmajor) { cblas_dcopy(P,&omegas[n+1],1,&Y[c*P],1); }
            else { cblas_dcopy(P,&omegas[n+1],1,&Y[c],C); }
        }
    }
    else if (dim==1)
    {
        for (r=0; r<R; r++)
        {
            if (iscolmajor) { cblas_dcopy(n,&X[r],R,poly,1); }
            else { cblas_dcopy(n,&X[r*C],1,poly,1); }
            poly[P] = z; prev[P] = -o;
            for (c=0; c<P; c++) { prev[c] = poly[n-c]; }
            for (c=0; c<n; c++) { poly[c] -= prev[c]; }
            cblas_dcopy(n*n,&z,0,compan,1); cblas_dcopy(P,&o,0,&compan[1],n+1); cblas_dcopy(n,poly,1,compan,n);
            info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in ar2lsf_d: lapacke decomposition failed\n"); return 1; }
            for (c=0; c<n; c++) { omegas[c] = (wi[c]>10.0*DBL_EPSILON) ? atan2(wi[c],wr[c]) : 0.0; poly[c] += prev[c] + prev[c]; }
            cblas_dcopy(n*n,&z,0,compan,1); cblas_dcopy(P,&o,0,&compan[1],n+1); cblas_dcopy(n,poly,1,compan,n);
            info = LAPACKE_dhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,compan,ldh,wr,wi,zz,ldz);  //eig
            if (info) { fprintf(stderr,"error in ar2lsf_d: lapacke decomposition failed\n"); return 1; }
            for (c=0; c<n; c++) { omegas[n+c] = (wi[c]>10.0*DBL_EPSILON) ? atan2(wi[c],wr[c]) : 0.0; }
            qsort(omegas,(size_t)(2*n),sizeof(double),cmp_ascend_d);
            if (iscolmajor) { cblas_dcopy(P,&omegas[n+1],1,&Y[r],R); }
            else { cblas_dcopy(P,&omegas[n+1],1,&Y[r*P],1); }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2lsf_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(compan); free(wr); free(wi); free(poly); free(prev); free(omegas);
    return 0;
}


int ar2lsf_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const float z[2] =  {0.0f,0.0f}, o[2] =  {1.0f,0.0f};
    const int P = (dim==0) ? R : C;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P, n = P, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    float *poly, *prev, *omegas, *compan, *wc, zz[2];
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2lsf_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2lsf_c: ncols X must be positive\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in ar2lsf_c: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    //Allocate
    if (!(poly=(float *)malloc((size_t)(2*n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_c: problem with malloc. "); perror("malloc"); return 1; }
    if (!(prev=(float *)malloc((size_t)(2*n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_c: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wc=(float *)malloc((size_t)(2*n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_c: problem with malloc. "); perror("malloc"); return 1; }
    if (!(omegas=(float *)malloc((size_t)(2*n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_c: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(float *)malloc((size_t)(4*n*n)*sizeof(float)))) { fprintf(stderr,"error in ar2lsf_c: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        for (c=0; c<C; c++)
        {
            if (iscolmajor) { cblas_ccopy(n,&X[2*c*R],1,poly,1); }
            else { cblas_ccopy(n,&X[2*c],C,poly,1); }
            poly[2*P] = poly[2*P+1] = prev[2*P+1] = 0.0f; prev[2*P] = -1.0f;
            for (r=0; r<P; r++) { prev[2*r] = poly[2*(P-r-1)]; prev[2*r+1] = poly[2*(P-r-1)+1]; }
            for (r=0; r<n; r++) { poly[2*r] -= prev[2*r]; poly[2*r+1] -= prev[2*r+1]; }
            cblas_ccopy(n*n,z,0,compan,1); cblas_ccopy(P,o,0,&compan[2],n+1); cblas_ccopy(n,poly,1,compan,n);
            info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)wc,(lapack_complex_float *)zz,ldz);
            if (info) { fprintf(stderr,"error in ar2lsf_c: lapacke decomposition failed\n"); return 1; }
            for (r=0; r<n; r++) { omegas[r] = atan2f(wc[2*r+1],wc[2*r]); poly[2*r] += prev[2*r] + prev[2*r]; poly[2*r+1] += prev[2*r+1] + prev[2*r+1]; }
            cblas_ccopy(n*n,z,0,compan,1); cblas_ccopy(P,o,0,&compan[2],n+1); cblas_ccopy(n,poly,1,compan,n);
            info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)wc,(lapack_complex_float *)zz,ldz);
            if (info) { fprintf(stderr,"error in ar2lsf_c: lapacke decomposition failed\n"); return 1; }
            for (r=0; r<n; r++) { omegas[n+r] = atan2f(wc[2*r+1],wc[2*r]); }
            qsort(omegas,(size_t)(2*n),sizeof(float),cmp_ascend_s);
            if (iscolmajor) { cblas_scopy(2*n,omegas,1,&Y[2*c*n],1); }
            else { cblas_scopy(2*n,omegas,1,&Y[c],C); }
        }
    }
    else if (dim==1)
    {
        for (r=0; r<R; r++)
        {
            if (iscolmajor) { cblas_ccopy(n,&X[2*r],R,poly,1); }
            else { cblas_ccopy(n,&X[2*r*C],1,poly,1); }
            poly[2*P] = poly[2*P+1] = prev[2*P+1] = 0.0f; prev[2*P] = -1.0f;
            for (c=0; c<P; c++) { prev[2*c] = poly[2*(P-c-1)]; prev[2*c+1] = poly[2*(P-c-1)+1]; }
            for (c=0; c<n; c++) { poly[2*c] -= prev[2*c]; poly[2*c+1] -= prev[2*c+1]; }
            cblas_ccopy(n*n,z,0,compan,1); cblas_ccopy(P,o,0,&compan[2],n+1); cblas_ccopy(n,poly,1,compan,n);
            info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)wc,(lapack_complex_float *)zz,ldz);
            if (info) { fprintf(stderr,"error in ar2lsf_c: lapacke decomposition failed\n"); return 1; }
            for (c=0; c<n; c++) { omegas[c] = atan2f(wc[2*c+1],wc[2*c]); poly[2*c] += prev[2*c] + prev[2*c]; poly[2*c+1] += prev[2*c+1] + prev[2*c+1]; }
            cblas_ccopy(n*n,z,0,compan,1); cblas_ccopy(P,o,0,&compan[2],n+1); cblas_ccopy(n,poly,1,compan,n);
            info = LAPACKE_chseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_float *)compan,ldh,(lapack_complex_float *)wc,(lapack_complex_float *)zz,ldz);
            if (info) { fprintf(stderr,"error in ar2lsf_c: lapacke decomposition failed\n"); return 1; }
            for (c=0; c<n; c++) { omegas[n+c] = atan2f(wc[2*c+1],wc[2*c]); }
            qsort(omegas,(size_t)(2*n),sizeof(float),cmp_ascend_s);
            if (iscolmajor) { cblas_scopy(2*n,omegas,1,&Y[r],R); }
            else { cblas_scopy(2*n,omegas,1,&Y[2*r*n],1); }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2lsf_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(compan); free(wc); free(poly); free(prev); free(omegas);
    return 0;
}


int ar2lsf_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const double z[2] =  {0.0,0.0}, o[2] =  {1.0,0.0};
    const int P = (dim==0) ? R : C;
    const char job = 'E', compz = 'N';  //eigenvalues only
    const lapack_int ldh = P, n = P, ldz = 1;
    const lapack_int ilo = 1, ihi = n;  //avoids balancing
    lapack_int info;
    double *poly, *prev, *omegas, *compan, *wc, zz[2];
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2lsf_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2lsf_z: ncols X must be positive\n"); return 1; }
    if (P<1) { fprintf(stderr,"error in ar2lsf_z: P (length of polynomial coeffs including a0=1) must be positive\n"); return 1; }

    //Allocate
    if (!(poly=(double *)malloc((size_t)(2*n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_z: problem with malloc. "); perror("malloc"); return 1; }
    if (!(prev=(double *)malloc((size_t)(2*n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_z: problem with malloc. "); perror("malloc"); return 1; }
    if (!(wc=(double *)malloc((size_t)(2*n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_z: problem with malloc. "); perror("malloc"); return 1; }
    if (!(omegas=(double *)malloc((size_t)(2*n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_z: problem with malloc. "); perror("malloc"); return 1; }
    if (!(compan=(double *)malloc((size_t)(4*n*n)*sizeof(double)))) { fprintf(stderr,"error in ar2lsf_z: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        for (c=0; c<C; c++)
        {
            if (iscolmajor) { cblas_zcopy(n,&X[2*c*R],1,poly,1); }
            else { cblas_zcopy(n,&X[2*c],C,poly,1); }
            poly[2*P] = poly[2*P+1] = prev[2*P+1] = 0.0; prev[2*P] = -1.0;
            for (r=0; r<P; r++) { prev[2*r] = poly[2*(P-r-1)]; prev[2*r+1] = poly[2*(P-r-1)+1]; }
            for (r=0; r<n; r++) { poly[2*r] -= prev[2*r]; poly[2*r+1] -= prev[2*r+1]; }
            cblas_zcopy(n*n,z,0,compan,1); cblas_zcopy(P,o,0,&compan[2],n+1); cblas_zcopy(n,poly,1,compan,n);
            info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)wc,(lapack_complex_double *)zz,ldz);
            if (info) { fprintf(stderr,"error in ar2lsf_z: lapacke decomposition failed\n"); return 1; }
            for (r=0; r<n; r++) { omegas[r] = atan2(wc[2*r+1],wc[2*r]); poly[2*r] += prev[2*r] + prev[2*r]; poly[2*r+1] += prev[2*r+1] + prev[2*r+1]; }
            cblas_zcopy(n*n,z,0,compan,1); cblas_zcopy(P,o,0,&compan[2],n+1); cblas_zcopy(n,poly,1,compan,n);
            info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)wc,(lapack_complex_double *)zz,ldz);
            if (info) { fprintf(stderr,"error in ar2lsf_z: lapacke decomposition failed\n"); return 1; }
            for (r=0; r<n; r++) { omegas[n+r] = atan2(wc[2*r+1],wc[2*r]); }
            qsort(omegas,(size_t)(2*n),sizeof(double),cmp_ascend_d);
            if (iscolmajor) { cblas_dcopy(2*n,omegas,1,&Y[2*c*n],1); }
            else { cblas_dcopy(2*n,omegas,1,&Y[c],C); }
        }
    }
    else if (dim==1)
    {
        for (r=0; r<R; r++)
        {
            if (iscolmajor) { cblas_zcopy(n,&X[2*r],R,poly,1); }
            else { cblas_zcopy(n,&X[2*r*C],1,poly,1); }
            poly[2*P] = poly[2*P+1] = prev[2*P+1] = 0.0; prev[2*P] = -1.0;
            for (c=0; c<P; c++) { prev[2*c] = poly[2*(P-c-1)]; prev[2*c+1] = poly[2*(P-c-1)+1]; }
            for (c=0; c<n; c++) { poly[2*c] -= prev[2*c]; poly[2*c+1] -= prev[2*c+1]; }
            cblas_zcopy(n*n,z,0,compan,1); cblas_zcopy(P,o,0,&compan[2],n+1); cblas_zcopy(n,poly,1,compan,n);
            info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)wc,(lapack_complex_double *)zz,ldz);
            if (info) { fprintf(stderr,"error in ar2lsf_z: lapacke decomposition failed\n"); return 1; }
            for (c=0; c<n; c++) { omegas[c] = atan2(wc[2*c+1],wc[2*c]); poly[2*c] += prev[2*c] + prev[2*c]; poly[2*c+1] += prev[2*c+1] + prev[2*c+1]; }
            cblas_zcopy(n*n,z,0,compan,1); cblas_zcopy(P,o,0,&compan[2],n+1); cblas_zcopy(n,poly,1,compan,n);
            info = LAPACKE_zhseqr(LAPACK_COL_MAJOR,job,compz,n,ilo,ihi,(lapack_complex_double *)compan,ldh,(lapack_complex_double *)wc,(lapack_complex_double *)zz,ldz);
            if (info) { fprintf(stderr,"error in ar2lsf_z: lapacke decomposition failed\n"); return 1; }
            for (c=0; c<n; c++) { omegas[n+c] = atan2(wc[2*c+1],wc[2*c]); }
            qsort(omegas,(size_t)(2*n),sizeof(double),cmp_ascend_d);
            if (iscolmajor) { cblas_dcopy(2*n,omegas,1,&Y[r],R); }
            else { cblas_dcopy(2*n,omegas,1,&Y[2*r*n],1); }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2lsf_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(compan); free(wc); free(poly); free(prev); free(omegas);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

