#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int get_deltas_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int N);
int get_deltas_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int N);

int get_delta_deltas_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int N);
int get_delta_deltas_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int N);

int add_deltas_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const int N);
int add_deltas_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const int N);

int add_delta_deltas_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const int N);
int add_delta_deltas_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const int N);

#ifdef __cplusplus
}
#endif
