#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int rms_scale_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const float fs, const float tau, const float target_dB_SPL, const float max_dB_SPL);
int rms_scale_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const double fs, const double tau, const double target_dB_SPL, const double max_dB_SPL);

int preemph_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const float p);
int preemph_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const double p);

int dither_s (float *X, const int N, const float d);
int dither_d (double *X, const int N, const double d);

#ifdef __cplusplus
}
#endif
