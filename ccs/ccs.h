#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int dct_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct);
int dct_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct);

int dct_inplace_s (float *X, const char iscolmajor, const int R, const int C, const int dim);
int dct_inplace_d (double *X, const char iscolmajor, const int R, const int C, const int dim);

int lifter_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const float Q);
int lifter_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const double Q);

int get_ccs_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct, const float Q, const int K);
int get_ccs_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int ndct, const double Q, const int K);

int mfccs_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const float *H, const int B, const int dim, const int c0, const float stp, const char mn0, const int nfft, const float preg, const int ndct, const float Q, const int K);
int mfccs_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const double *H, const int B, const int dim, const int c0, const double stp, const char mn0, const int nfft, const double preg, const int ndct, const double Q, const int K);

#ifdef __cplusplus
}
#endif
