#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int moments_s (float *X, const char iscolmajor, const int R, const int C, const int dim);
int moments_d (double *X, const char iscolmajor, const int R, const int C, const int dim);

#ifdef __cplusplus
}
#endif
