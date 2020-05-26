#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int zcs_s (char *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int going);
int zcs_d (char *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int going);
int zcs_c (char *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int going);
int zcs_z (char *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int going);

int lcs_s (char *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int going, const float level);
int lcs_d (char *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int going, const double level);

int mcs_s (char *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int going);
int mcs_d (char *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int going);

int zcr_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const float stp, const int going);
int zcr_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const double stp, const int going);

int lcr_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const float stp, const int going, const float level);
int lcr_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const double stp, const int going, const double level);

int mcr_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const float stp, const int going);
int mcr_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int L, const int dim, const int c0, const double stp, const int going);

int zcr_windowed_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *W, const int L, const int dim, const int c0, const float stp, const int going);
int zcr_windowed_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *W, const int L, const int dim, const int c0, const double stp, const int going);

int lcr_windowed_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *W, const int L, const int dim, const int c0, const float stp, const int going, const float level);
int lcr_windowed_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *W, const int L, const int dim, const int c0, const double stp, const int going, const double level);

int mcr_windowed_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *W, const int L, const int dim, const int c0, const float stp, const int going);
int mcr_windowed_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *W, const int L, const int dim, const int c0, const double stp, const int going);

#ifdef __cplusplus
}
#endif
