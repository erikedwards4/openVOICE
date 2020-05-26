#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int ar2lsf_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2lsf_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2lsf_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2lsf_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int lsf2ar_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int lsf2ar_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int poly2lsf_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2lsf_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2lsf_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2lsf_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int lsf2poly_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int lsf2poly_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int rc2lar_s (float *Y, const float *X, const int N);
int rc2lar_d (double *Y, const double *X, const int N);
int rc2lar_c (float *Y, const float *X, const int N);
int rc2lar_z (double *Y, const double *X, const int N);

int lar2rc_s (float *Y, const float *X, const int N);
int lar2rc_d (double *Y, const double *X, const int N);
int lar2rc_c (float *Y, const float *X, const int N);
int lar2rc_z (double *Y, const double *X, const int N);

#ifdef __cplusplus
}
#endif
