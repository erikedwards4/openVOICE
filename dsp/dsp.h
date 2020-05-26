#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int fft_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);

int fft_matmul_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_matmul_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_matmul_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_matmul_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);

int fir_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *B, const int L, const int dim, const int stride);
int fir_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *B, const int L, const int dim, const int stride);
int fir_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *B, const int L, const int dim, const int stride);
int fir_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *B, const int L, const int dim, const int stride);

int iir_s (float *X, const char iscolmajor, const int R, const int C, const float *A, const int N, const int dim);
int iir_d (double *X, const char iscolmajor, const int R, const int C, const double *A, const int N, const int dim);
int iir_c (float *X, const char iscolmajor, const int R, const int C, const float *A, const int N, const int dim);
int iir_z (double *X, const char iscolmajor, const int R, const int C, const double *A, const int N, const int dim);

#ifdef __cplusplus
}
#endif
