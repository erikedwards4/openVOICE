#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int autocorr_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L);
int autocorr_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L);
int autocorr_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L);
int autocorr_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L);

int autocorr_fft_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L);
int autocorr_fft_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L);
int autocorr_fft_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L);
int autocorr_fft_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L);

int sig2ac_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased);
int sig2ac_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased);
int sig2ac_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased);
int sig2ac_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased);

int sig2ac_fft_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased);
int sig2ac_fft_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased);
int sig2ac_fft_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased);
int sig2ac_fft_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased);

int ac2ar_levdurb_s (float *Y, float *V, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int P);
int ac2ar_levdurb_d (double *Y, double *V, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int P);

int ac2poly_levdurb_s (float *Y, float *V, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int P);
int ac2poly_levdurb_d (double *Y, double *V, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int P);

int sig2ar_levdurb_s (float *Y, float *V, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0);
int sig2ar_levdurb_d (double *Y, double *V, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0);

int sig2poly_levdurb_s (float *Y, float *V, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0);
int sig2poly_levdurb_d (double *Y, double *V, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0);

int sig2ar_burg_s (float *Y, float *V, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0);
int sig2ar_burg_d (double *Y, double *V, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0);

int sig2poly_burg_s (float *Y, float *V, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0);
int sig2poly_burg_d (double *Y, double *V, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0);

int ac2rc_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int ac2rc_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);

int ac2cc_s (float *CC, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int K, const float preg);
int ac2cc_d (double *CC, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int K, const double preg);

int ac2mvdr_s (float *MVDR, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int F, const float preg);
int ac2mvdr_d (double *MVDR, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int F, const double preg);

#ifdef __cplusplus
}
#endif
