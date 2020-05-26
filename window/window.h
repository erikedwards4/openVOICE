#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int frame_univar_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const int dim, const int c0, const float stp);
int frame_univar_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const int dim, const int c0, const double stp);
int frame_univar_c (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const int dim, const int c0, const float stp);
int frame_univar_z (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const int dim, const int c0, const double stp);

int apply_win_s (float *X, const char iscolmajor, const int R, const int C, const float *W, const int dim);
int apply_win_d (double *X, const char iscolmajor, const int R, const int C, const double *W, const int dim);
int apply_win_c (float *X, const char iscolmajor, const int R, const int C, const float *W, const int dim);
int apply_win_z (double *X, const char iscolmajor, const int R, const int C, const double *W, const int dim);

int window_univar_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int dim, const int c0, const float stp, const char mn0);
int window_univar_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int dim, const int c0, const double stp, const char mn0);
int window_univar_c (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int dim, const int c0, const float stp, const char mn0);
int window_univar_z (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int dim, const int c0, const double stp, const char mn0);

#ifdef __cplusplus
}
#endif
