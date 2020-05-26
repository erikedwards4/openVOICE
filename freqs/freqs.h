#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "ihcs.h"
#include "sgcs.h"

int interp1q_s (float *Xo, const int No, const float *Xi, const float *Yi, const int Ni);
int interp1q_d (double *Xo, const int No, const double *Xi, const double *Yi, const int Ni);

int convert_freqs_s (float *frqs, const int F, const char in_scale[], const char out_scale[]);
int convert_freqs_d (double *frqs, const int F, const char in_scale[], const char out_scale[]);

int get_cfs_s (float *cfs, const int B, const float lofreq, const float hifreq, const char freq_scale[], const char band_edges);
int get_cfs_d (double *cfs, const int B, const double lofreq, const double hifreq, const char freq_scale[], const char band_edges);

int get_cfs_T_s (float *cfs, const int B, const float fs, const char freq_scale[]);
int get_cfs_T_d (double *cfs, const int B, const double fs, const char freq_scale[]);

int get_stft_freqs_s (float *freqs, const int F, const float fs, const int nfft);
int get_stft_freqs_d (double *freqs, const int F, const double fs, const int nfft);

#ifdef __cplusplus
}
#endif
