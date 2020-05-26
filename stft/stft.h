#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int fft_hc_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_hc_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);

int hc_square_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int nfft, const int dim);
int hc_square_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int nfft, const int dim);

int fft_squared_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_squared_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);

int stft_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const int dim, const int c0, const float stp, const char mn0, const int nfft);
int stft_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const int dim, const int c0, const double stp, const char mn0, const int nfft);

int stft2_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const int dim, const int c0, const float stp, const char mn0, const int nfft);
int stft2_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const int dim, const int c0, const double stp, const char mn0, const int nfft);

#ifdef __cplusplus
}
#endif
