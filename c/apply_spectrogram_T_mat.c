//This takes the output of STFT, X, and a transform matrix, T,
//and outputs a T*X or X*T', depending on orientation of X.
//T must be BxF, as from get_spectrogram_T_mat.

#include <stdio.h>
#include <float.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int apply_spectrogram_T_mat_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *T, const int B, const int F);
int apply_spectrogram_T_mat_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *T, const int B, const int F);


int apply_spectrogram_T_mat_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *T, const int B, const int F)
{
    const float z = 0.0f;
    int b, f, n = 0;
    
    //Checks
    if (R<1) { fprintf(stderr,"error in apply_spectrogram_T_mat_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in apply_spectrogram_T_mat_s: C (ncols X) must be positive\n"); return 1; }
    if (B<1) { fprintf(stderr,"error in apply_spectrogram_T_mat_s: B (nrows T) must be positive\n"); return 1; }
    if (F<1) { fprintf(stderr,"error in apply_spectrogram_T_mat_s: F (ncols T) must be positive\n"); return 1; }

    if (R==F) //dim==0
    {
        cblas_scopy(B*C,&z,0,&Y[0],1); //make sure Y is initialized to 0
        if (iscolmajor)
        {
            for (f=0; f<F; f++)
            {
                for (b=0; b<B; b++, n++)
                {
                    if (T[n]>FLT_EPSILON || -T[n]>FLT_EPSILON) { cblas_saxpy(C,T[n],&X[f],F,&Y[b],B); }
                }
            }
        }
        else
        {
            for (b=0; b<B; b++)
            {
                for (f=0; f<F; f++, n++)
                {
                    if (T[n]>FLT_EPSILON || -T[n]>FLT_EPSILON) { cblas_saxpy(C,T[n],&X[f*C],1,&Y[b*C],1); }
                }
            }
        }
    }
    else if (C==F) //dim==1
    {
        cblas_scopy(B*R,&z,0,&Y[0],1); //make sure Y is initialized to 0
        if (iscolmajor)
        {
            for (f=0; f<F; f++)
            {
                for (b=0; b<B; b++, n++)
                {
                    if (T[n]>FLT_EPSILON || -T[n]>FLT_EPSILON) { cblas_saxpy(R,T[n],&X[f*R],1,&Y[b*R],1); }
                }
            }
            
        }
        else
        {
            for (b=0; b<B; b++)
            {
                for (f=0; f<F; f++, n++)
                {
                    if (T[n]>FLT_EPSILON || -T[n]>FLT_EPSILON) { cblas_saxpy(R,T[n],&X[f],F,&Y[b],B); }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in apply_spectrogram_T_mat_s: T mat not size compatible with X\n"); return 1;
    }
    
    return 0;
}


int apply_spectrogram_T_mat_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *T, const int B, const int F)
{
    const double z = 0.0;
    int b, f, n = 0;
    
    //Checks
    if (R<1) { fprintf(stderr,"error in apply_spectrogram_T_mat_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in apply_spectrogram_T_mat_d: C (ncols X) must be positive\n"); return 1; }
    if (B<1) { fprintf(stderr,"error in apply_spectrogram_T_mat_d: B (nrows T) must be positive\n"); return 1; }
    if (F<1) { fprintf(stderr,"error in apply_spectrogram_T_mat_d: F (ncols T) must be positive\n"); return 1; }

    if (R==F) //dim==0
    {
        cblas_dcopy(B*C,&z,0,&Y[0],1); //make sure Y is initialized to 0
        if (iscolmajor)
        {
            for (f=0; f<F; f++)
            {
                for (b=0; b<B; b++, n++)
                {
                    if (T[n]>DBL_EPSILON || -T[n]>DBL_EPSILON) { cblas_daxpy(C,T[n],&X[f],F,&Y[b],B); }
                }
            }
        }
        else
        {
            for (b=0; b<B; b++)
            {
                for (f=0; f<F; f++, n++)
                {
                    if (T[n]>DBL_EPSILON || -T[n]>DBL_EPSILON) { cblas_daxpy(C,T[n],&X[f*C],1,&Y[b*C],1); }
                }
            }
        }
    }
    else if (C==F) //dim==1
    {
        cblas_dcopy(B*R,&z,0,&Y[0],1); //make sure Y is initialized to 0
        if (iscolmajor)
        {
            for (f=0; f<F; f++)
            {
                for (b=0; b<B; b++, n++)
                {
                    if (T[n]>DBL_EPSILON || -T[n]>DBL_EPSILON) { cblas_daxpy(R,T[n],&X[f*R],1,&Y[b*R],1); }
                }
            }
            
        }
        else
        {
            for (b=0; b<B; b++)
            {
                for (f=0; f<F; f++, n++)
                {
                    if (T[n]>DBL_EPSILON || -T[n]>DBL_EPSILON) { cblas_daxpy(R,T[n],&X[f],F,&Y[b],B); }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in apply_spectrogram_T_mat_d: T mat not size compatible with X\n"); return 1;
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif

