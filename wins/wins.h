#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int bartlett_s (float *X, const int L, const char normalize);
int bartlett_d (double *X, const int L, const char normalize);
int bartlett_c (float *X, const int L, const char normalize);
int bartlett_z (double *X, const int L, const char normalize);

int blackman_s (float *X, const int L, const char normalize, const char exact);
int blackman_d (double *X, const int L, const char normalize, const char exact);
int blackman_c (float *X, const int L, const char normalize, const char exact);
int blackman_z (double *X, const int L, const char normalize, const char exact);

int blackmanharris_s (float *X, const int L, const char normalize);
int blackmanharris_d (double *X, const int L, const char normalize);
int blackmanharris_c (float *X, const int L, const char normalize);
int blackmanharris_z (double *X, const int L, const char normalize);

int flattop_s (float *X, const int L, const char normalize);
int flattop_d (double *X, const int L, const char normalize);
int flattop_c (float *X, const int L, const char normalize);
int flattop_z (double *X, const int L, const char normalize);

int hamming_s (float *X, const int L, const char normalize);
int hamming_d (double *X, const int L, const char normalize);
int hamming_c (float *X, const int L, const char normalize);
int hamming_z (double *X, const int L, const char normalize);

int hann_s (float *X, const int L, const char normalize);
int hann_d (double *X, const int L, const char normalize);
int hann_c (float *X, const int L, const char normalize);
int hann_z (double *X, const int L, const char normalize);

int rectangular_s (float *X, const int L, const char normalize);
int rectangular_d (double *X, const int L, const char normalize);
int rectangular_c (float *X, const int L, const char normalize);
int rectangular_z (double *X, const int L, const char normalize);

int triangular_s (float *X, const int L, const char normalize);
int triangular_d (double *X, const int L, const char normalize);
int triangular_c (float *X, const int L, const char normalize);
int triangular_z (double *X, const int L, const char normalize);

int povey_s (float *X, const int L, const char normalize);
int povey_d (double *X, const int L, const char normalize);
int povey_c (float *X, const int L, const char normalize);
int povey_z (double *X, const int L, const char normalize);

#ifdef __cplusplus
}
#endif
