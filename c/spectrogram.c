//This takes univariate X, windows X with W, applies the 1D FFT to each frame,
//takes power, transforms with T mat to B cfs, and then compresses power.
//The output Y is power (real-valued), i.e. the FFT squared.
//The size of Y is BxC if dim==0, and RxB if dim==1.

//This is equivalent to stft -> apply_spectrogram_T_mat -> pow_compress.

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int spectrogram_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const float *H, const int dim, const int c0, const float stp, const char mn0, const int nfft, const float p, const float preg);
int spectrogram_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const double *H, const int dim, const int c0, const double stp, const char mn0, const int nfft, const double p, const double preg);


int spectrogram_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const float *H, const int dim, const int c0, const float stp, const char mn0, const int nfft, const float p, const float preg)
{
    const float z = 0.0f, o = 1.0f;
    const int F = nfft/2 + 1;  //num non-negative FFT freqs
    const int Lpre = L/2;      //nsamps before center samp
    int ss = c0 - Lpre;        //start samp of current frame
    int r, c, f, n = 0;
    float *X1, *Y1;
    fftwf_plan plan;
    //struct timespec tic, toc;

    //Checks
    if (N<1) { fprintf(stderr,"error in spectrogram_s: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in spectrogram_s: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in spectrogram_s: C (ncols Y) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in spectrogram_s: L (winlength) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in spectrogram_s: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in spectrogram_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in spectrogram_s: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0f) { fprintf(stderr,"error in spectrogram_s: stp (step size) must be positive\n"); return 1; }
    if (nfft<L) { fprintf(stderr,"error in spectrogram_s: nfft must be >= L (winlength)\n"); return 1; }
    if (p!=p || p<0.0f || p>1.0f) { fprintf(stderr,"error in spectrogram_s: p must be in [0.0 1.0]\n"); return 1; }
    if (preg!=preg || preg<0.0f) { fprintf(stderr,"error in spectrogram_s: preg must be nonnegative\n"); return 1; }

    //Initialize fftwf
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_real((size_t)nfft);
    plan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in spectrogram_s: problem creating fftw plan\n"); return 1; }
    cblas_scopy(nfft,&z,0,&X1[0],1); //zero-pad

    //clock_gettime(CLOCK_REALTIME,&tic);

    //Initialize Y to preg
    cblas_scopy(R*C,&preg,0,&Y[0],1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                //Window
                ss = (int)(roundf(c*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_scopy(-ss,&z,0,&X1[0],1);
                    cblas_ssbmv(CblasColMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                    cblas_scopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_scopy(L,&z,0,&X1[0],1);
                }

                //Subtract mean
                if (mn0) { cblas_saxpy(L,-cblas_sdot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }

                //FFT
                fftwf_execute(plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_sgemv(CblasColMajor,CblasNoTrans,R,F,1.0f,&H[0],R,&Y1[0],1,1.0f,&Y[c*R],1);
                //for (b=0; b<R; b++) { Y[b+c*R] += cblas_sdot(f2s[b],&H[b+f1s[b]*R],R,&Y1[f1s[b]],1); }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                //Window
                ss = (int)(roundf(c*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_scopy(-ss,&z,0,&X1[0],1);
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                    cblas_scopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_scopy(L,&z,0,&X1[0],1);
                }

                //Subtract mean
                if (mn0) { cblas_saxpy(L,-cblas_sdot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }

                //FFT
                fftwf_execute(plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_sgemv(CblasRowMajor,CblasNoTrans,R,F,1.0f,&H[0],F,&Y1[0],1,1.0f,&Y[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                //Window
                ss = (int)(roundf(r*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_scopy(-ss,&z,0,&X1[0],1);
                    cblas_ssbmv(CblasColMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                    cblas_scopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_scopy(L,&z,0,&X1[0],1);
                }

                //Subtract mean
                if (mn0) { cblas_saxpy(L,-cblas_sdot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }

                //FFT
                fftwf_execute(plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_sgemv(CblasColMajor,CblasNoTrans,C,F,1.0f,&H[0],C,&Y1[0],1,1.0f,&Y[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                //Window
                ss = (int)(roundf(r*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_scopy(-ss,&z,0,&X1[0],1);
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                    cblas_scopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_scopy(L,&z,0,&X1[0],1);
                }

                //Subtract mean
                if (mn0) { cblas_saxpy(L,-cblas_sdot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }

                //FFT
                fftwf_execute(plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_sgemv(CblasRowMajor,CblasNoTrans,C,F,1.0f,&H[0],F,&Y1[0],1,1.0f,&Y[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in spectrogram_s: dim must be 0 or 1.\n"); return 1;
    }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    //Power compress
    if (1.0f-p>FLT_EPSILON)
    {
        if (p>FLT_EPSILON)
        {
            if (fabsf(p-0.5f)<=FLT_EPSILON) { for (n=0; n<R*C; n++) { Y[n] = sqrtf(Y[n]+preg); } }
            else if (fabsf(p-1.0f/3.0f)<=FLT_EPSILON) { for (n=0; n<R*C; n++) { Y[n] = cbrtf(Y[n]+preg); } }
            else { for (n=0; n<R*C; n++) { Y[n] = powf(Y[n]+preg,p); } }
        }
        else { for (n=0; n<R*C; n++) { Y[n] = logf(Y[n]+preg); } }
    }
    else { for (n=0; n<R*C; n++) { Y[n] += preg; } }
    
    //Exit
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int spectrogram_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const double *H, const int dim, const int c0, const double stp, const char mn0, const int nfft, const double p, const double preg)
{
    const double z = 0.0, o = 1.0;
    const int F = nfft/2 + 1;  //num non-negative FFT freqs
    const int Lpre = L/2;      //nsamps before center samp
    int ss = c0 - Lpre;        //start samp of current frame
    int r, c, f, n = 0;
    double *X1, *Y1;
    fftw_plan plan;

    //Checks
    if (N<1) { fprintf(stderr,"error in spectrogram_d: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in spectrogram_d: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in spectrogram_d: C (ncols Y) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in spectrogram_d: L (winlength) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in spectrogram_d: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in spectrogram_d: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in spectrogram_d: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0) { fprintf(stderr,"error in spectrogram_d: stp (step size) must be positive\n"); return 1; }
    if (nfft<L) { fprintf(stderr,"error in spectrogram_d: nfft must be >= L (winlength)\n"); return 1; }
    if (p!=p || p<0.0 || p>1.0) { fprintf(stderr,"error in spectrogram_d: p must be in [0.0 1.0]\n"); return 1; }
    if (preg!=preg || preg<0.0) { fprintf(stderr,"error in spectrogram_d: preg must be nonnegative\n"); return 1; }

    //Initialize fftw
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_real((size_t)nfft);
    plan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in spectrogram_d: problem creating fftw plan\n"); return 1; }
    cblas_dcopy(nfft,&z,0,&X1[0],1); //zero-pad

    //Initialize Y to preg
    cblas_dcopy(R*C,&preg,0,&Y[0],1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                //Window
                ss = (int)(round(c*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&X1[0],1);
                    cblas_dsbmv(CblasColMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                    cblas_dcopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_dcopy(L,&z,0,&X1[0],1);
                }

                //Subtract mean
                if (mn0) { cblas_daxpy(L,-cblas_ddot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }

                //FFT
                fftw_execute(plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_dgemv(CblasColMajor,CblasNoTrans,R,F,1.0,&H[0],R,&Y1[0],1,1.0,&Y[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                //Window
                ss = (int)(round(c*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&X1[0],1);
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                    cblas_dcopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_dcopy(L,&z,0,&X1[0],1);
                }

                //Subtract mean
                if (mn0) { cblas_daxpy(L,-cblas_ddot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }

                //FFT
                fftw_execute(plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_dgemv(CblasRowMajor,CblasNoTrans,R,F,1.0,&H[0],F,&Y1[0],1,1.0,&Y[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                //Window
                ss = (int)(round(r*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&X1[0],1);
                    cblas_dsbmv(CblasColMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                    cblas_dcopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_dcopy(L,&z,0,&X1[0],1);
                }

                //Subtract mean
                if (mn0) { cblas_daxpy(L,-cblas_ddot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }

                //FFT
                fftw_execute(plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_dgemv(CblasColMajor,CblasNoTrans,C,F,1.0,&H[0],C,&Y1[0],1,1.0,&Y[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                //Window
                ss = (int)(round(r*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&X1[0],1);
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                    cblas_dcopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_dcopy(L,&z,0,&X1[0],1);
                }

                //Subtract mean
                if (mn0) { cblas_daxpy(L,-cblas_ddot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }

                //FFT
                fftw_execute(plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_dgemv(CblasRowMajor,CblasNoTrans,C,F,1.0,&H[0],F,&Y1[0],1,1.0,&Y[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in spectrogram_d: dim must be 0 or 1.\n"); return 1;
    }

    //Power compress
    if (1.0-p>DBL_EPSILON)
    {
        if (p>DBL_EPSILON)
        {
            if (fabs(p-0.5)<=DBL_EPSILON) { for (n=0; n<R*C; n++) { Y[n] = sqrt(Y[n]+preg); } }
            else if (fabs(p-1.0/3.0)<=DBL_EPSILON) { for (n=0; n<R*C; n++) { Y[n] = cbrt(Y[n]+preg); } }
            else { for (n=0; n<R*C; n++) { Y[n] = pow(Y[n]+preg,p); } }
        }
        else { for (n=0; n<R*C; n++) { Y[n] = log(Y[n]+preg); } }
    }
    else { for (n=0; n<R*C; n++) { Y[n] += preg; } }
    
    //Exit
    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif


// //Get Hb (banded H) //this worked to get banded matrix storage (but was slower!):
// if (!(Hb=(float *)calloc((size_t)((kl+ku+1)*B),sizeof(float)))) { fprintf(stderr,"error in spectrogram_s: problem with calloc for Hb. "); perror("calloc"); return 1; }
// if (iscolmajor)
// {
//     for (r=0; r<kl+1+ku; r++)
//     {
//         if (r<ku) { b = 0; f = ku - r; }
//         else { b = r - ku; f = 0; }
//         while (b<B && f<F) { Hb[r+f*(kl+1+ku)] = H[b+f*B]; b++; f++; }
//     }
// }
// else
// {
//     for (c=0; c<kl+1+ku; c++)
//     {
//         if (c<kl) { b = kl - c; f = 0; }
//         else { b = 0; f = c - kl; }
//         while (b<B && f<F) { Hb[c+b*(kl+1+ku)] = H[f+b*F]; b++; f++; }
//     }
// }
// ...
// cblas_sgbmv(CblasRowMajor,CblasNoTrans,B,F,kl,ku,1.0f,&Hb[0],kl+ku+1,&Y1[0],1,1.0f,&Y[c],C);

