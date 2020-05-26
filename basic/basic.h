#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int mean0_s (float *X, const char iscolmajor, const int R, const int C, const int dim);
int mean0_d (double *X, const char iscolmajor, const int R, const int C, const int dim);
int mean0_c (float *X, const char iscolmajor, const int R, const int C, const int dim);
int mean0_z (double *X, const char iscolmajor, const int R, const int C, const int dim);

int stdev1_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const char biased);
int stdev1_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const char biased);
int stdev1_c (float *X, const char iscolmajor, const int R, const int C, const int dim, const char biased);
int stdev1_z (double *X, const char iscolmajor, const int R, const int C, const int dim, const char biased);

int zscore_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const char biased);
int zscore_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const char biased);
int zscore_c (float *X, const char iscolmajor, const int R, const int C, const int dim, const char biased);
int zscore_z (double *X, const char iscolmajor, const int R, const int C, const int dim, const char biased);

int square_s (float *Y, const float *X, const int N);
int square_d (double *Y, const double *X, const int N);
int square_c (float *Y, const float *X, const int N);
int square_z (double *Y, const double *X, const int N);

int abs_s (float *Y, const float *X, const int N);
int abs_d (double *Y, const double *X, const int N);
int abs_c (float *Y, const float *X, const int N);
int abs_z (double *Y, const double *X, const int N);

int cmp_ascend_s (const void *a, const void *b);
int cmp_ascend_d (const void *a, const void *b);
int cmp_ascend_c (const void *a, const void *b);
int cmp_ascend_z (const void *a, const void *b);

int cmp_descend_s (const void *a, const void *b);
int cmp_descend_d (const void *a, const void *b);
int cmp_descend_c (const void *a, const void *b);
int cmp_descend_z (const void *a, const void *b);

int median_s (float *X, const char iscolmajor, const int R, const int C, const int dim);
int median_d (double *X, const char iscolmajor, const int R, const int C, const int dim);

#ifdef __cplusplus
}
#endif
