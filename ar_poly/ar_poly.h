#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int roots_s (float *roots, const float *AS, const int P);
int roots_d (double *roots, const double *AS, const int P);
int roots_c (float *roots, const float *AS, const int P);
int roots_z (double *roots, const double *AS, const int P);

int poly_s (float *AS, const float *roots, const int R);
int poly_d (double *AS, const double *roots, const int R);
int poly_c (float *AS, const float *roots, const int R);
int poly_z (double *AS, const double *roots, const int R);

int poly2roots_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2roots_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2roots_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2roots_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int roots2poly_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int roots2poly_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int roots2poly_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int roots2poly_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int poly2ar_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2ar_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2ar_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2ar_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int ar2poly_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2poly_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2poly_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2poly_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int ar2rc_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2rc_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2rc_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ar2rc_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int rc2ar_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int rc2ar_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int rc2ar_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int rc2ar_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int poly2rc_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2rc_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2rc_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2rc_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int rc2poly_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int rc2poly_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int rc2poly_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int rc2poly_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int ar2psd_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *V, const float *W, const int F, const int dim);
int ar2psd_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *V, const double *W, const int F, const int dim);
int ar2psd_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *V, const float *W, const int F, const int dim);
int ar2psd_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *V, const double *W, const int F, const int dim);

int poly2psd_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *V, const float *W, const int F, const int dim);
int poly2psd_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *V, const double *W, const int F, const int dim);
int poly2psd_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *V, const float *W, const int F, const int dim);
int poly2psd_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *V, const double *W, const int F, const int dim);

#ifdef __cplusplus
}
#endif
