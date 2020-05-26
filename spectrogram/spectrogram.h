#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int get_cns_s (int *cns, const float *freqs, const int F, const float *cfs, const int B);
int get_cns_d (int *cns, const double *freqs, const int F, const double *cfs, const int B);

int get_spectrogram_T_mat_s (float *H, const char iscolmajor, const float *freqs, const int F, const float *cfs, const int B, const char normalize);
int get_spectrogram_T_mat_d (double *H, const char iscolmajor, const double *freqs, const int F, const double *cfs, const int B, const char normalize);

int apply_spectrogram_T_mat_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *T, const int B, const int F);
int apply_spectrogram_T_mat_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *T, const int B, const int F);

int pow_compress_s (float *X, const int N, const float p, const float preg);
int pow_compress_d (double *X, const int N, const double p, const double preg);

int spectrogram_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const float *H, const int dim, const int c0, const float stp, const char mn0, const int nfft, const float p, const float preg);
int spectrogram_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const double *H, const int dim, const int c0, const double stp, const char mn0, const int nfft, const double p, const double preg);

#ifdef __cplusplus
}
#endif
