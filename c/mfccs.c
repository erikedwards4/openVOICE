//This takes univariate X, windows X with W, applies the 1D FFT to each frame,
//takes power, transforms with T mat to B cfs, and then compresses power with log.
//Then the cepstral coeffs (CCs) are obtained by the 1D DCT-II of each frame,
//followed by a lifter (weighting in the cepstral domain).

//This is equivalent to :                                                                                 spectrogram -> get_ccs.
//which is equivalent to:                                                                                 spectrogram -> dct -> lifter.
//which is equivalent to:                                             stft -> apply_spectrogram_T_mat -> pow_compress -> dct -> lifter.
//which is equivalent to:                     window_univar -> fft_squared -> apply_spectrogram_T_mat -> pow_compress -> dct -> lifter.
//which is equivalent to:             window_univar -> fft_hc -> hc_square -> apply_spectrogram_T_mat -> pow_compress -> dct -> lifter.
//which is equivalent to: frame_univar -> apply_win -> fft_hc -> hc_square -> apply_spectrogram_T_mat -> pow_compress -> dct -> lifter.

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>

#ifndef M_PIf
   #define M_PIf 3.141592653589793238462643383279502884f
#endif

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int mfccs_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const float *H, const int B, const int dim, const int c0, const float stp, const char mn0, const int nfft, const float preg, const int ndct, const float Q, const int K);
int mfccs_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const double *H, const int B, const int dim, const int c0, const double stp, const char mn0, const int nfft, const double preg, const int ndct, const double Q, const int K);


//int mfccs_s (float *Y, sigi Yi, const float *X, const sigi Xi, const float *W, const sigi Wi, const float *H, const sigi H, const int dim, const int c0, const float stp, const float preg, const int ndct, const float Q)
int mfccs_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const float *H, const int B, const int dim, const int c0, const float stp, const char mn0, const int nfft, const float preg, const int ndct, const float Q, const int K)
{
    const float z = 0.0f, o = 1.0f;
    const float Q2 = 0.5f*Q, PQ = M_PIf/Q;
    //const char normalize_H = 0;
    const int F = nfft/2 + 1;   //num non-negative FFT freqs
    const int F2 = F-1+nfft%2;  //for square loop below
    const int Lpre = L/2;       //nsamps before center samp
    int ss = c0 - Lpre;         //start samp of current frame
    int r, c, f, b, k;
    float *X1, *Y1, *X1d, *Y1d; //, *freqs, *cfs, *H;
    fftwf_plan fft_plan, dct_plan;
    //struct timespec tic, toc;

    //Checks
    if (N<1) { fprintf(stderr,"error in mfccs_s: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in mfccs_s: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mfccs_s: C (ncols Y) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in mfccs_s: L (winlength) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in mfccs_s: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in mfccs_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in mfccs_s: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0f) { fprintf(stderr,"error in mfccs_s: stp (step size) must be positive\n"); return 1; }
    if (nfft<L) { fprintf(stderr,"error in mfccs_s: nfft must be >= L (winlength)\n"); return 1; }
    if (preg!=preg || preg<0.0f) { fprintf(stderr,"error in mfccs_s: preg must be nonnegative\n"); return 1; }
    if (Q<0.0f) { fprintf(stderr,"error in mfccs_s: Q must be positive (or 0 to skip lifter)\n"); return 1; }
    if (K<1) { fprintf(stderr,"error in mfccs_s: K must be positive\n"); return 1; }
    if (K>ndct) { fprintf(stderr,"error in mfccs_s: K must be <= ndct>\n"); return 1; }

    //Initialize fftwf for FFT
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_real((size_t)nfft);
    fft_plan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fft_plan) { fprintf(stderr,"error in mfccs_s: problem creating fftw plan for fft\n"); return 1; }
    cblas_scopy(nfft,&z,0,&X1[0],1); //zero-pad

    //Initialize fftwf for DCT
    X1d = fftwf_alloc_real((size_t)ndct);
    Y1d = fftwf_alloc_real((size_t)ndct);
    dct_plan = fftwf_plan_r2r_1d(ndct,X1d,Y1d,FFTW_REDFT10,FFTW_ESTIMATE);
    if (!dct_plan) { fprintf(stderr,"error in mfccs_s: problem creating fftw plan for dct\n"); return 1; }
    cblas_scopy(ndct,&z,0,&X1d[0],1); //zero-pad

    //Get frequencies and H
    //if (get_stft_freqs_s(&freqs[0],F,fs,nfft)) { fprintf(stderr,"error in mfccs_s: problem getting STFT freqs\n"); return 1; }
    //if (get_cfs_T_s(&cfs[0],B,fs,'mel')) { fprintf(stderr,"error in mfccs_s: problem getting Mel cfs\n"); return 1; }
    //if (get_spectrogram_T_mat_s(&H[0],iscolmajor,&freqs[0],F,&cfs[0],B,normalize_H))
    //{ fprintf(stderr,"error in mfccs_s: problem getting STFT freq -> Mel cfs transform matrix\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);

    //Initialize lifter
    //float *lift; int n;
    //if (!(lift=(float *)malloc((size_t)K*sizeof(float)))) { fprintf(stderr,"error in get_ccs_s: problem with malloc for lifter. "); perror("malloc"); return 1; }
    //for (k=0; k<K; k++) { lift[k] = fmaf(Q2,sinf(k*PQ),1.0f); }

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
                fftwf_execute(fft_plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_sgemv(CblasColMajor,CblasNoTrans,B,F,1.0f,&H[0],B,&Y1[0],1,0.0f,&X1d[0],1);

                //Power compress
                for (b=0; b<B; b++) { X1d[b] = logf(X1d[b]+preg); }

                //DCT
                fftwf_execute(dct_plan);
                cblas_scopy(K,&Y1d[0],1,&Y[c*K],1);
                //cblas_ssbmv(CblasColMajor,CblasUpper,K,0,1.0f,&lift[0],1,&Y1d[0],1,0.0f,&Y[c*K],1);
                //for (k=0; k<K; k++) { Y[k+c*K] = lift[k]*Y1d[k]; }
            }
            //Lifter
            if (Q>FLT_EPSILON) { for (k=0; k<K; k++) { cblas_sscal(C,fmaf(Q2,sinf(PQ*k),1.0f),&Y[k],K); } }
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
                fftwf_execute(fft_plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_sgemv(CblasRowMajor,CblasNoTrans,B,F,1.0f,&H[0],F,&Y1[0],1,0.0f,&X1d[0],1);

                //Power compress
                for (b=0; b<B; b++) { X1d[b] = logf(X1d[b]+preg); }

                //DCT
                fftwf_execute(dct_plan);
                cblas_scopy(K,&Y1d[0],1,&Y[c],C);
            }
            //Lifter
            if (Q>FLT_EPSILON) { for (k=0; k<K; k++) { cblas_sscal(C,fmaf(Q2,sinf(PQ*k),1.0f),&Y[k*C],1); } }
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
                fftwf_execute(fft_plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_sgemv(CblasColMajor,CblasNoTrans,B,F,1.0f,&H[0],B,&Y1[0],1,0.0f,&X1d[0],1);

                //Power compress
                for (b=0; b<B; b++) { X1d[b] = logf(X1d[b]+preg); }

                //DCT
                fftwf_execute(dct_plan);
                cblas_scopy(K,&Y1d[0],1,&Y[r],R);
            }
            //Lifter
            if (Q>FLT_EPSILON) { for (k=0; k<K; k++) { cblas_sscal(R,fmaf(Q2,sinf(PQ*k),1.0f),&Y[k*R],1); } }
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
                fftwf_execute(fft_plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_sgemv(CblasRowMajor,CblasNoTrans,B,F,1.0f,&H[0],F,&Y1[0],1,0.0f,&X1d[0],1);

                //Power compress
                for (b=0; b<B; b++) { X1d[b] = logf(X1d[b]+preg); }

                //DCT
                fftwf_execute(dct_plan);
                cblas_scopy(K,&Y1d[0],1,&Y[r*K],1);
            }
            //Lifter
            if (Q>FLT_EPSILON) { for (k=0; k<K; k++) { cblas_sscal(R,fmaf(Q2,sinf(PQ*k),1.0f),&Y[k],K); } }
        }
    }
    else
    {
        fprintf(stderr,"error in mfccs_s: dim must be 0 or 1.\n"); return 1;
    }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    //Exit
    fftwf_destroy_plan(fft_plan); fftwf_free(X1); fftwf_free(Y1);
    fftwf_destroy_plan(dct_plan); fftwf_free(X1d); fftwf_free(Y1d);
    return 0;
}


int mfccs_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const double *H, const int B, const int dim, const int c0, const double stp, const char mn0, const int nfft, const double preg, const int ndct, const double Q, const int K)
{
    const double z = 0.0, o = 1.0;
    const double Q2 = 0.5*Q, PQ = M_PI/Q;
    const int F = nfft/2 + 1;   //num non-negative FFT freqs
    const int F2 = F-1+nfft%2;  //for square loop below
    const int Lpre = L/2;       //nsamps before center samp
    int ss = c0 - Lpre;         //start samp of current frame
    int r, c, f, b, k;
    double *X1, *Y1, *X1d, *Y1d;
    fftw_plan fft_plan, dct_plan;

    //Checks
    if (N<1) { fprintf(stderr,"error in mfccs_d: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in mfccs_d: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mfccs_d: C (ncols Y) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in mfccs_d: L (winlength) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in mfccs_d: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in mfccs_d: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in mfccs_d: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0) { fprintf(stderr,"error in mfccs_d: stp (step size) must be positive\n"); return 1; }
    if (nfft<L) { fprintf(stderr,"error in mfccs_d: nfft must be >= L (winlength)\n"); return 1; }
    if (preg!=preg || preg<0.0) { fprintf(stderr,"error in mfccs_d: preg must be nonnegative\n"); return 1; }
    if (Q<0.0) { fprintf(stderr,"error in mfccs_d: Q must be positive (or 0 to skip lifter)\n"); return 1; }
    if (K<1) { fprintf(stderr,"error in mfccs_d: K must be positive\n"); return 1; }
    if (K>ndct) { fprintf(stderr,"error in mfccs_d: K must be <= ndct>\n"); return 1; }

    //Initialize fftw for FFT
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_real((size_t)nfft);
    fft_plan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fft_plan) { fprintf(stderr,"error in mfccs_d: problem creating fftw plan for fft\n"); return 1; }
    cblas_dcopy(nfft,&z,0,&X1[0],1); //zero-pad

    //Initialize fftw for DCT
    X1d = fftw_alloc_real((size_t)ndct);
    Y1d = fftw_alloc_real((size_t)ndct);
    dct_plan = fftw_plan_r2r_1d(ndct,X1d,Y1d,FFTW_REDFT10,FFTW_ESTIMATE);
    if (!dct_plan) { fprintf(stderr,"error in mfccs_d: problem creating fftw plan for dct\n"); return 1; }
    cblas_dcopy(ndct,&z,0,&X1d[0],1); //zero-pad

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
                fftw_execute(fft_plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_dgemv(CblasColMajor,CblasNoTrans,B,F,1.0,&H[0],B,&Y1[0],1,0.0,&X1d[0],1);

                //Power compress
                for (b=0; b<B; b++) { X1d[b] = log(X1d[b]+preg); }

                //DCT
                fftw_execute(dct_plan);
                cblas_dcopy(K,&Y1d[0],1,&Y[c*K],1);
            }
            //Lifter
            if (Q>DBL_EPSILON) { for (k=0; k<K; k++) { cblas_dscal(C,fma(Q2,sin(PQ*k),1.0),&Y[k],K); } }
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
                fftw_execute(fft_plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_dgemv(CblasRowMajor,CblasNoTrans,B,F,1.0,&H[0],F,&Y1[0],1,0.0,&X1d[0],1);

                //Power compress
                for (b=0; b<B; b++) { X1d[b] = log(X1d[b]+preg); }

                //DCT
                fftw_execute(dct_plan);
                cblas_dcopy(K,&Y1d[0],1,&Y[c],C);
            }
            //Lifter
            if (Q>DBL_EPSILON) { for (k=0; k<K; k++) { cblas_dscal(C,fma(Q2,sin(PQ*k),1.0),&Y[k*C],1); } }
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
                fftw_execute(fft_plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_dgemv(CblasColMajor,CblasNoTrans,B,F,1.0,&H[0],B,&Y1[0],1,0.0,&X1d[0],1);

                //Power compress
                for (b=0; b<B; b++) { X1d[b] = log(X1d[b]+preg); }

                //DCT
                fftw_execute(dct_plan);
                cblas_dcopy(K,&Y1d[0],1,&Y[r],R);
            }
            //Lifter
            if (Q>DBL_EPSILON) { for (k=0; k<K; k++) { cblas_dscal(R,fma(Q2,sin(PQ*k),1.0),&Y[k*R],1); } }
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
                fftw_execute(fft_plan);

                //Power
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F2; f++) { Y1[f] += Y1[nfft-f]; }
                
                //Transform freq scale
                cblas_dgemv(CblasRowMajor,CblasNoTrans,B,F,1.0,&H[0],F,&Y1[0],1,0.0,&X1d[0],1);

                //Power compress
                for (b=0; b<B; b++) { X1d[b] = log(X1d[b]+preg); }

                //DCT
                fftw_execute(dct_plan);
                cblas_dcopy(K,&Y1d[0],1,&Y[r*K],1);
            }
            //Lifter
            if (Q>DBL_EPSILON) { for (k=0; k<K; k++) { cblas_dscal(R,fma(Q2,sin(PQ*k),1.0),&Y[k],K); } }
        }
    }
    else
    {
        fprintf(stderr,"error in mfccs_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    fftw_destroy_plan(fft_plan); fftw_free(X1); fftw_free(Y1);
    fftw_destroy_plan(dct_plan); fftw_free(X1d); fftw_free(Y1d);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

