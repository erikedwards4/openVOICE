//This takes the output of window_univar, and applies the 1D FFT to each row or col.
//The input dim is which dim to take the 1D FFT along (0->along cols, 1->along rows).
//The output Y is power (real-valued), i.e. the FFT squared, Y = X*conj(X) element-wise.
//The size of Y must be FxC if dim==0, and RxF if dim==1, where F = nfft/2 + 1.

//This v2 (originally v1) uses window_univar and fft_squared quite literally in order,
//but requires a large intermdiate matrix XW to hold the windowed data, and cblas_?copy operations.
//Nonetheless, it can be slightly faster than v1 (which is frame-by-frame).
//However, it is just as often slightly slower than v1 (depends on dim and iscolmajor).
//But, it requires much more code and will be less adaptable to other methods, so I use v1 as default.

//To test compile:
//gcc -c stft2.c -O2 -std=c99 -Wall -Wextra
//clang -c stft2.c -O2 -std=c99 -Weverything
//g++ -c stft2.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c stft2.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int stft2_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const int dim, const int c0, const float stp, const char mn0, const int nfft)
{
    const float z = 0.0f, o = 1.0f;
    const int T = (dim==0) ? C : R;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    const int F = nfft/2 + 1;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;
    float *XW;                 //intermediate output (as of window_univar)
    int r, c, f, n = 0;
    float *X1, *Y1;
    fftwf_plan plan;

    //Checks
    if (N<1) { fprintf(stderr,"error in stft2_s: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in stft2_s: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in stft2_s: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in stft2_s: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in stft2_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in stft2_s: L (winlength) must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in stft2_s: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0f) { fprintf(stderr,"error in stft2_s: stp (step size) must be positive\n"); return 1; }
    if (nfft<L) { fprintf(stderr,"error in stft2_s: nfft must be >= L (winlength)\n"); return 1; }
    if (dim==0 && F!=R) { fprintf(stderr,"error in stft2_s: F must equal nrows in X for dim==0\n"); return 1; }
    if (dim==1 && F!=C) { fprintf(stderr,"error in stft2_s: F must equal ncols in X for dim==1\n"); return 1; }

    //Initialize fftwf
    XW = fftwf_alloc_real((size_t)(T*L));
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_real((size_t)nfft);
    while (n<nfft) { X1[n++] = 0.0f; }
    plan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in stft2_s: problem creating fftw plan\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&XW[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&XW[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&XW[l+L*Tpre],L); cblas_sscal(Tmid,W[l],&XW[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&XW[t*L],1);
                    cblas_ssbmv(CblasColMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&XW[t*L-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&XW[t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&XW[t*L],1);
                    cblas_scopy(L-N+ss,&z,0,&XW[t*L+N-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_saxpy(L,-cblas_sdot(L,&XW[t*L],1,&o,0)/L,&o,0,&XW[t*L],1); } }
            for (c=0; c<C; c++)
            {
                cblas_scopy(L,&XW[c*L],1,&X1[0],1);
                fftwf_execute(plan);
                n = c*F; Y[n++] = Y1[0]*Y1[0];
                for (f=1; f<F-1; f++, n++) { Y[n] = Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f]; }
                Y[n] = (nfft%2) ? Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f] : Y1[f]*Y1[f];
            }
        }
        else
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&XW[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&XW[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&XW[l*T+Tpre],1); cblas_sscal(Tmid,W[l],&XW[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&XW[t],T);
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&XW[-T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&XW[t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&XW[t],T);
                    cblas_scopy(L-N+ss,&z,0,&XW[t+T*(N-ss)],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_saxpy(L,-cblas_sdot(L,&XW[t],T,&o,0)/L,&o,0,&XW[t],T); } }
            for (c=0; c<C; c++)
            {
                cblas_scopy(L,&XW[c],C,&X1[0],1);
                fftwf_execute(plan);
                n = c; Y[n] = Y1[0]*Y1[0]; n += C;
                for (f=1; f<F-1; f++, n+=C) { Y[n] = Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f]; }
                Y[n] = (nfft%2) ? Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f] : Y1[f]*Y1[f];
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&XW[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&XW[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&XW[l*T+Tpre],1); cblas_sscal(Tmid,W[l],&XW[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&XW[t],T);
                    cblas_ssbmv(CblasColMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&XW[-T*ss],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&XW[t],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&XW[t],T);
                    cblas_scopy(L-N+ss,&z,0,&XW[t+T*(N-ss)],T);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_saxpy(L,-cblas_sdot(L,&XW[t],T,&o,0)/L,&o,0,&XW[t],T); } }
            for (r=0; r<R; r++)
            {
                cblas_scopy(L,&XW[r],R,&X1[0],1);
                fftwf_execute(plan);
                n = r; Y[n] = Y1[0]*Y1[0]; n += R;
                for (f=1; f<F-1; f++, n+=R) { Y[n] = Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f]; }
                Y[n] = (nfft%2) ? Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f] : Y1[f]*Y1[f];
            }
        }
        else
        {
            if (stp==floorf(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_scopy(Tpre,&z,0,&XW[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_scopy(Tpost,&z,0,&XW[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_scopy(Tmid,&X[ssl+Tpre*istp],istp,&XW[l+L*Tpre],L); cblas_sscal(Tmid,W[l],&XW[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_scopy(-ss,&z,0,&XW[t*L],1);
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&XW[t*L-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&XW[t*L],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&XW[t*L],1);
                    cblas_scopy(L-N+ss,&z,0,&XW[t*L+N-ss],1);
                    t++; ss = (int)(roundf(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_saxpy(L,-cblas_sdot(L,&XW[t*L],1,&o,0)/L,&o,0,&XW[t*L],1); } }
            for (r=0; r<R; r++)
            {
                cblas_scopy(L,&XW[r*L],1,&X1[0],1);
                fftwf_execute(plan);
                n = r*F; Y[n++] = Y1[0]*Y1[0];
                for (f=1; f<F-1; f++, n++) { Y[n] = Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f]; }
                Y[n] = (nfft%2) ? Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f] : Y1[f]*Y1[f];
            }
        }
    }
    else
    {
        fprintf(stderr,"error in stft2_s: dim must be 0 or 1.\n"); return 1;
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(XW); fftwf_free(Y1);
    return 0;
}


int stft2_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const int dim, const int c0, const double stp, const char mn0, const int nfft)
{
    const double z = 0.0, o = 1.0;
    const int T = (dim==0) ? C : R;
    const int Lpre = L/2;      //nsamps before center samp
    const int istp = (int)stp;
    const int F = nfft/2 + 1;
    int ss = c0 - Lpre;        //start samp of current frame
    int l, t = 0;
    int Tpre, Tmid, Tpost;
    int ssl = ss, esl = ss + (T-1)*istp;
    double *XW;                 //intermediate output (as of window_univar)
    int r, c, f, n = 0;
    double *X1, *Y1;
    fftw_plan plan;

    //Checks
    if (N<1) { fprintf(stderr,"error in stft2_d: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in stft2_d: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in stft2_d: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in stft2_d: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in stft2_d: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in stft2_d: L (winlength) must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in stft2_d: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0) { fprintf(stderr,"error in stft2_d: stp (step size) must be positive\n"); return 1; }
    if (nfft<L) { fprintf(stderr,"error in stft2_d: nfft must be >= L (winlength)\n"); return 1; }
    if (dim==0 && F!=R) { fprintf(stderr,"error in stft2_d: F must equal nrows in X for dim==0\n"); return 1; }
    if (dim==1 && F!=C) { fprintf(stderr,"error in stft2_d: F must equal ncols in X for dim==1\n"); return 1; }

    //Initialize fftw
    XW = fftw_alloc_real((size_t)(T*L));
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_real((size_t)nfft);
    while (n<nfft) { X1[n++] = 0.0; }
    plan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in stft2_d: problem creating fftw plan\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&XW[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&XW[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&XW[l+L*Tpre],L); cblas_dscal(Tmid,W[l],&XW[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&XW[t*L],1);
                    cblas_dsbmv(CblasColMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&XW[t*L-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&XW[t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&XW[t*L],1);
                    cblas_dcopy(L-N+ss,&z,0,&XW[t*L+N-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_daxpy(L,-cblas_ddot(L,&XW[t*L],1,&o,0)/L,&o,0,&XW[t*L],1); } }
            for (c=0; c<C; c++)
            {
                cblas_dcopy(L,&XW[c*L],1,&X1[0],1);
                fftw_execute(plan);
                n = c*F; Y[n++] = Y1[0]*Y1[0];
                for (f=1; f<F-1; f++, n++) { Y[n] = Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f]; }
                Y[n] = (nfft%2) ? Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f] : Y1[f]*Y1[f];
            }
        }
        else
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&XW[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&XW[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&XW[l*T+Tpre],1); cblas_dscal(Tmid,W[l],&XW[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&XW[t],T);
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&XW[-T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&XW[t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&XW[t],T);
                    cblas_dcopy(L-N+ss,&z,0,&XW[t+T*(N-ss)],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_daxpy(L,-cblas_ddot(L,&XW[t],T,&o,0)/L,&o,0,&XW[t],T); } }
            for (c=0; c<C; c++)
            {
                cblas_dcopy(L,&XW[c],C,&X1[0],1);
                fftw_execute(plan);
                n = c; Y[n] = Y1[0]*Y1[0]; n += C;
                for (f=1; f<F-1; f++, n+=C) { Y[n] = Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f]; }
                Y[n] = (nfft%2) ? Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f] : Y1[f]*Y1[f];
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&XW[l*T],1); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&XW[l*T+Tpre+Tmid],1); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&XW[l*T+Tpre],1); cblas_dscal(Tmid,W[l],&XW[l*T+Tpre],1); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&XW[t],T);
                    cblas_dsbmv(CblasColMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&XW[-T*ss],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&XW[t],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&XW[t],T);
                    cblas_dcopy(L-N+ss,&z,0,&XW[t+T*(N-ss)],T);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_daxpy(L,-cblas_ddot(L,&XW[t],T,&o,0)/L,&o,0,&XW[t],T); } }
            for (r=0; r<R; r++)
            {
                cblas_dcopy(L,&XW[r],R,&X1[0],1);
                fftw_execute(plan);
                n = r; Y[n] = Y1[0]*Y1[0]; n += R;
                for (f=1; f<F-1; f++, n+=R) { Y[n] = Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f]; }
                Y[n] = (nfft%2) ? Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f] : Y1[f]*Y1[f];
            }
        }
        else
        {
            if (stp==floor(stp))
            {
                for (l=0; l<L; l++)
                {
                    Tmid = T;
                    if (ssl<0) { Tpre = 1 - ssl/istp; Tmid -= Tpre; cblas_dcopy(Tpre,&z,0,&XW[l],L); } else { Tpre = 0; }
                    if (esl>=N) { Tpost = 1 + (esl-N)/istp; Tmid -= Tpost; cblas_dcopy(Tpost,&z,0,&XW[l+L*(Tpre+Tmid)],L); }
                    if (Tmid>0) { cblas_dcopy(Tmid,&X[ssl+Tpre*istp],istp,&XW[l+L*Tpre],L); cblas_dscal(Tmid,W[l],&XW[l+L*Tpre],L); }
                    ssl++; esl++;
                }
            }
            else
            {
                while (t<T && ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&XW[t*L],1);
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&XW[t*L-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss+L<=N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&XW[t*L],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
                while (t<T && ss<N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&XW[t*L],1);
                    cblas_dcopy(L-N+ss,&z,0,&XW[t*L+N-ss],1);
                    t++; ss = (int)(round(t*stp)) + c0 - Lpre;
                }
            }
            if (mn0) { for (t=0; t<T; t++) { cblas_daxpy(L,-cblas_ddot(L,&XW[t*L],1,&o,0)/L,&o,0,&XW[t*L],1); } }
            for (r=0; r<R; r++)
            {
                cblas_dcopy(L,&XW[r*L],1,&X1[0],1);
                fftw_execute(plan);
                n = r*F; Y[n++] = Y1[0]*Y1[0];
                for (f=1; f<F-1; f++, n++) { Y[n] = Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f]; }
                Y[n] = (nfft%2) ? Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f] : Y1[f]*Y1[f];
            }
        }
    }
    else
    {
        fprintf(stderr,"error in stft2_d: dim must be 0 or 1.\n"); return 1;
    }
    
    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(XW); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
#endif

