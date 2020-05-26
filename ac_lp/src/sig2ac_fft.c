//The "unbiased" version uses N-l in the denominator instead of N.
//It is actually just "less biased", but is slower,
//has larger mean-squared error, and doesn't match FFT estimate.

//Recall that FFT is order N*log(N), instead of N*N for sig2ac.
//Since 2 FFTs (and set up time), this means use this if 2*log(N) < N.
//However, this is the case for N>2, so only set up time would prevent using this.
//On quick local-system test, this is worth using for N>?.

//To test compile:
//gcc -c sig2ac_fft.c -O2 -std=c99 -Wall -Wextra
//clang -c sig2ac_fft.c -O2 -std=c99 -Weverything
//g++ -c sig2ac_fft.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c sig2ac_fft.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int sig2ac_fft_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased)
{
    const float z = 0.0f; //, o = 1.0f;
    //float m;
    int r, c, l, f, nfft, F;
    float *X1, *Y1;
    fftwf_plan fplan, iplan;

    //Checks
    if (R<1) { fprintf(stderr,"error in sig2ac_fft_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sig2ac_fft_s: ncols X must be positive\n"); return 1; }
    if (dim==0 && L>R) { fprintf(stderr,"error in sig2ac_fft_s: nlags must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && L>C) { fprintf(stderr,"error in sig2ac_fft_s: nlags must be < ncols X for dim==1\n"); return 1; }

    //Get nfft
    nfft = (dim==0) ? R+L : C+L;
    if (nfft>16384) { nfft += nfft%2; }
    else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }
    F = nfft/2 + 1;

    //Initialize fftwf
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_real((size_t)nfft);
    fplan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    iplan = fftwf_plan_r2r_1d(nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_s: problem creating fftw plan"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_scopy(R,&X[c*R],1,&X1[0],1);
                //m = cblas_sdot(R,&X1[0],1,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X1[0],1);
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                cblas_scopy(L,&X1[0],1,&Y[c*L],1);
            }
            cblas_sscal(L*C,1.0f/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_sscal(C,(float)R/(float)(R-l),&Y[l],L); } }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_scopy(R,&X[c],C,&X1[0],1);
                //m = cblas_sdot(R,&X1[0],1,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X1[0],1);
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                cblas_scopy(L,&X1[0],1,&Y[c],1);
            }
            cblas_sscal(L*C,1.0f/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_sscal(C,(float)R/(float)(R-l),&Y[l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_scopy(C,&X[r],R,&X1[0],1);
                //m = cblas_sdot(C,&X1[0],1,&o,0) / C;
                //cblas_saxpy(C,-m,&o,0,&X1[0],1);
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                cblas_scopy(L,&X1[0],1,&Y[r],R);
            }
            cblas_sscal(L*R,1.0f/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_sscal(R,(float)C/(float)(C-l),&Y[l*R],1); } }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_scopy(C,&X[r*C],1,&X1[0],1);
                //m = cblas_sdot(C,&X1[0],1,&o,0) / C;
                //cblas_saxpy(C,-m,&o,0,&X1[0],1);
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                cblas_scopy(L,&X1[0],1,&Y[r*L],1);
            }
            cblas_sscal(L*R,1.0f/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_sscal(R,(float)C/(float)(C-l),&Y[l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2ac_fft_s: dim must be 0 or 1.\n"); return 1;
    }

    fftwf_destroy_plan(fplan); fftwf_destroy_plan(iplan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int sig2ac_fft_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased)
{
    const double z = 0.0;
    int r, c, l, f, nfft, F;
    double *X1, *Y1;
    fftw_plan fplan, iplan;

    //Checks
    if (R<1) { fprintf(stderr,"error in sig2ac_fft_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sig2ac_fft_d: ncols X must be positive\n"); return 1; }
    if (dim==0 && L>R) { fprintf(stderr,"error in sig2ac_fft_d: nlags must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && L>C) { fprintf(stderr,"error in sig2ac_fft_d: nlags must be < ncols X for dim==1\n"); return 1; }

    //Get nfft
    nfft = (dim==0) ? R+L : C+L;
    if (nfft>16384) { nfft += nfft%2; }
    else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }
    F = nfft/2 + 1;

    //Initialize fftw
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_real((size_t)nfft);
    fplan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    iplan = fftw_plan_r2r_1d(nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_d: problem creating fftw plan"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_dcopy(R,&X[c*R],1,&X1[0],1);
                //m = cblas_ddot(R,&X1[0],1,&o,0) / R;
                //cblas_daxpy(R,-m,&o,0,&X1[0],1);
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                cblas_dcopy(L,&X1[0],1,&Y[c*L],1);
            }
            cblas_dscal(L*C,1.0/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_dscal(C,(double)R/(double)(R-l),&Y[l],L); } }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_dcopy(R,&X[c],C,&X1[0],1);
                //m = cblas_ddot(R,&X1[0],1,&o,0) / R;
                //cblas_daxpy(R,-m,&o,0,&X1[0],1);
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                cblas_dcopy(L,&X1[0],1,&Y[c],1);
            }
            cblas_dscal(L*C,1.0/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_dscal(C,(double)R/(double)(R-l),&Y[l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_dcopy(C,&X[r],R,&X1[0],1);
                //m = cblas_ddot(C,&X1[0],1,&o,0) / C;
                //cblas_daxpy(C,-m,&o,0,&X1[0],1);
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                cblas_dcopy(L,&X1[0],1,&Y[r],R);
            }
            cblas_dscal(L*R,1.0/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_dscal(R,(double)C/(double)(C-l),&Y[l*R],1); } }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_dcopy(C,&X[r*C],1,&X1[0],1);
                //m = cblas_ddot(C,&X1[0],1,&o,0) / C;
                //cblas_daxpy(C,-m,&o,0,&X1[0],1);
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                cblas_dcopy(L,&X1[0],1,&Y[r*L],1);
            }
            cblas_dscal(L*R,1.0/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_dscal(R,(double)C/(double)(C-l),&Y[l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2ac_fft_d: dim must be 0 or 1.\n"); return 1;
    }

    fftw_destroy_plan(fplan); fftw_destroy_plan(iplan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


int sig2ac_fft_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased)
{
    const float z[2] = {0.0f,0.0f}; //, o[2] = {1.0f,0.0f};
    //_Complex float m;
    int r, c, l, f, nfft;
    fftwf_complex *X1, *Y1;
    fftwf_plan fplan, iplan;

    //Checks
    if (R<1) { fprintf(stderr,"error in sig2ac_fft_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sig2ac_fft_c: ncols X must be positive\n"); return 1; }
    if (dim==0 && L>R) { fprintf(stderr,"error in sig2ac_fft_c: nlags must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && L>C) { fprintf(stderr,"error in sig2ac_fft_c: nlags must be < ncols X for dim==1\n"); return 1; }

    //Get nfft
    nfft = (dim==0) ? R+L : C+L;
    if (nfft>16384) { nfft += nfft%2; }
    else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }

    //Initialize fftw
    X1 = fftwf_alloc_complex((size_t)nfft);
    Y1 = fftwf_alloc_complex((size_t)nfft);
    fplan = fftwf_plan_dft_1d(nfft,X1,Y1,FFTW_FORWARD,FFTW_ESTIMATE);
    iplan = fftwf_plan_dft_1d(nfft,Y1,X1,FFTW_BACKWARD,FFTW_ESTIMATE);
    if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_c: problem creating fftw plan"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(nfft-R,&z[0],0,(float *)&X1[R],1); //zero-pad
                cblas_ccopy(R,&X[2*c*R],1,(float *)&X1[0],1);
                //m = -cblas_cdotu(R,(float *)&X1[0],1,&o[0],0) / R;
                //cblas_caxpy(R,(float *)&m,&o[0],0,(float *)&X1[0],1);
                fftwf_execute(fplan);
                #ifdef __cplusplus
                    for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0f; }
                #else
                    for (f=0; f<nfft; f++) { Y1[f] *= conjf(Y1[f]); }
                #endif
                fftwf_execute(iplan);
                cblas_ccopy(L,(float *)&X1[0],1,&Y[2*c*L],1);
            }
            cblas_csscal(L*C,1.0f/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_csscal(C,(float)R/(float)(R-l),&Y[2*l],L); } }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(nfft-R,&z[0],0,(float *)&X1[R],1); //zero-pad
                cblas_ccopy(R,&X[2*c],2*C,(float *)&X1[0],1);
                //m = -cblas_cdotu(R,(float *)&X1[0],1,&o[0],0) / R;
                //cblas_caxpy(R,(float *)&m,&o[0],0,(float *)&X1[0],1);
                fftwf_execute(fplan);
                #ifdef __cplusplus
                    for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0f; }
                #else
                    for (f=0; f<nfft; f++) { Y1[f] *= conjf(Y1[f]); }
                #endif
                fftwf_execute(iplan);
                cblas_ccopy(L,(float *)&X1[0],1,&Y[2*c],C);
            }
            cblas_csscal(L*C,1.0f/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_csscal(C,(float)R/(float)(R-l),&Y[2*l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(nfft-C,&z[0],0,(float *)&X1[C],1); //zero-pad
                cblas_ccopy(C,&X[2*r],2*R,(float *)&X1[0],1);
                //m = -cblas_cdotu(C,(float *)&X1[0],1,&o[0],0) / C;
                //cblas_caxpy(C,(float *)&m,&o[0],0,(float *)&X1[0],1);
                fftwf_execute(fplan);
                #ifdef __cplusplus
                    for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0f; }
                #else
                    for (f=0; f<nfft; f++) { Y1[f] *= conjf(Y1[f]); }
                #endif
                fftwf_execute(iplan);
                cblas_ccopy(L,(float *)&X1[0],1,&Y[2*r],R);
            }
            cblas_csscal(L*R,1.0f/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_csscal(R,(float)C/(float)(C-l),&Y[2*l*R],1); } }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(nfft-C,&z[0],0,(float *)&X1[C],1); //zero-pad
                cblas_ccopy(C,&X[2*r*C],1,(float *)&X1[0],1);
                //m = -cblas_cdotu(C,(float *)&X1[0],1,&o[0],0) / C;
                //cblas_caxpy(C,(float *)&m,&o[0],0,(float *)&X1[0],1);
                fftwf_execute(fplan);
                #ifdef __cplusplus
                    for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0f; }
                #else
                    for (f=0; f<nfft; f++) { Y1[f] *= conjf(Y1[f]); }
                #endif
                fftwf_execute(iplan);
                cblas_ccopy(L,(float *)&X1[0],1,&Y[2*r*L],1);
            }
            cblas_csscal(L*R,1.0f/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_csscal(R,(float)C/(float)(C-l),&Y[2*l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2ac_fft_c: dim must be 0 or 1.\n"); return 1;
    }

    fftwf_destroy_plan(fplan); fftwf_destroy_plan(iplan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int sig2ac_fft_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int L, const char unbiased)
{
    const double z[2] = {0.0,0.0};
    int r, c, l, f, nfft;
    fftw_complex *X1, *Y1;
    fftw_plan fplan, iplan;

    //Checks
    if (R<1) { fprintf(stderr,"error in sig2ac_fft_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sig2ac_fft_z: ncols X must be positive\n"); return 1; }
    if (dim==0 && L>R) { fprintf(stderr,"error in sig2ac_fft_z: nlags must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && L>C) { fprintf(stderr,"error in sig2ac_fft_z: nlags must be < ncols X for dim==1\n"); return 1; }

    //Get nfft
    nfft = (dim==0) ? R+L : C+L;
    if (nfft>16384) { nfft += nfft%2; }
    else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }

    //Initialize fftw
    X1 = fftw_alloc_complex((size_t)nfft);
    Y1 = fftw_alloc_complex((size_t)nfft);
    fplan = fftw_plan_dft_1d(nfft,X1,Y1,FFTW_FORWARD,FFTW_ESTIMATE);
    iplan = fftw_plan_dft_1d(nfft,Y1,X1,FFTW_BACKWARD,FFTW_ESTIMATE);
    if (!fplan || !iplan) { fprintf(stderr,"error in sig2ac_fft_z: problem creating fftw plan"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_zcopy(nfft-R,&z[0],0,(double *)&X1[R],1); //zero-pad
                cblas_zcopy(R,&X[2*c*R],1,(double *)&X1[0],1);
                //m = -cblas_zdotu(R,(double *)&X1[0],1,&o[0],0) / R;
                //cblas_zaxpy(R,(double *)&m,&o[0],0,(double *)&X1[0],1);
                fftw_execute(fplan);
                #ifdef __cplusplus
                    for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0; }
                #else
                    for (f=0; f<nfft; f++) { Y1[f] *= conj(Y1[f]); }
                #endif
                fftw_execute(iplan);
                cblas_zcopy(L,(double *)&X1[0],1,&Y[2*c*L],1);
            }
            cblas_zdscal(L*C,1.0/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_zdscal(C,(double)R/(double)(R-l),&Y[2*l],L); } }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_zcopy(nfft-R,&z[0],0,(double *)&X1[R],1); //zero-pad
                cblas_zcopy(R,&X[2*c],2*C,(double *)&X1[0],1);
                //m = -cblas_zdotu(R,(double *)&X1[0],1,&o[0],0) / R;
                //cblas_zaxpy(R,(double *)&m,&o[0],0,(double *)&X1[0],1);
                fftw_execute(fplan);
                #ifdef __cplusplus
                    for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0; }
                #else
                    for (f=0; f<nfft; f++) { Y1[f] *= conj(Y1[f]); }
                #endif
                fftw_execute(iplan);
                cblas_zcopy(L,(double *)&X1[0],1,&Y[2*c],C);
            }
            cblas_zdscal(L*C,1.0/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_zdscal(C,(double)R/(double)(R-l),&Y[2*l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_zcopy(nfft-C,&z[0],0,(double *)&X1[C],1); //zero-pad
                cblas_zcopy(C,&X[2*r],2*R,(double *)&X1[0],1);
                //m = -cblas_zdotu(C,(double *)&X1[0],1,&o[0],0) / C;
                //cblas_zaxpy(C,(double *)&m,&o[0],0,(double *)&X1[0],1);
                fftw_execute(fplan);
                #ifdef __cplusplus
                    for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0; }
                #else
                    for (f=0; f<nfft; f++) { Y1[f] *= conj(Y1[f]); }
                #endif
                fftw_execute(iplan);
                cblas_zcopy(L,(double *)&X1[0],1,&Y[2*r],R);
            }
            cblas_zdscal(L*R,1.0/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_zdscal(R,(double)C/(double)(C-l),&Y[2*l*R],1); } }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_zcopy(nfft-C,&z[0],0,(double *)&X1[C],1); //zero-pad
                cblas_zcopy(C,&X[2*r*C],1,(double *)&X1[0],1);
                //m = -cblas_zdotu(C,(double *)&X1[0],1,&o[0],0) / C;
                //cblas_zaxpy(C,(double *)&m,&o[0],0,(double *)&X1[0],1);
                fftw_execute(fplan);
                #ifdef __cplusplus
                    for (f=0; f<nfft; f++) { Y1[f][0] = Y1[f][0]*Y1[f][0] + Y1[f][1]*Y1[f][1]; Y1[f][1] = 0.0; }
                #else
                    for (f=0; f<nfft; f++) { Y1[f] *= conj(Y1[f]); }
                #endif
                fftw_execute(iplan);
                cblas_zcopy(L,(double *)&X1[0],1,&Y[2*r*L],1);
            }
            cblas_zdscal(L*R,1.0/nfft,&Y[0],1);
            if (unbiased) { for (l=1; l<L; l++) { cblas_zdscal(R,(double)C/(double)(C-l),&Y[2*l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2ac_fft_z: dim must be 0 or 1.\n"); return 1;
    }

    fftw_destroy_plan(fplan); fftw_destroy_plan(iplan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
#endif

