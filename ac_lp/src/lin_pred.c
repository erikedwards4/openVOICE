//This does lin_pred followed by lev_durb for each row/col of X.
//The autocov is the biased version (e.g. default of Octave, etc.).
//The inputs/outputs are the same as lev_durb,
//except works directly from X rather than from AC of X.

//An opt mean0 is added to zero the mean of each row/col of X first.
//In this case, this is mean0 -> autocov_fft -> lev_durb.

//To test compile:
//gcc -c lin_pred.c -O2 -std=c99 -Wall -Wextra
//clang -c lin_pred.c -O2 -std=c99 -Weverything
//g++ -c lin_pred.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c lin_pred.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int lin_pred_s (float *AS, float *V, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0)
{
    const float z = 0.0f, o = 1.0f;
    float m, g, sc;
    int r, c, f, nfft, F, p, q;
    float *X1, *Y1, *AStmp;
    fftwf_plan fplan, iplan;

    //Checks
    if (R<1) { fprintf(stderr,"error in lin_pred_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lin_pred_s: ncols X must be positive\n"); return 1; }
    if (dim==0 && P>R) { fprintf(stderr,"error in lin_pred_s: P must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && P>C) { fprintf(stderr,"error in lin_pred_s: P must be < ncols X for dim==1\n"); return 1; }

    //Get nfft
    nfft = (dim==0) ? R+P : C+P;
    if (nfft>16384) { nfft += nfft%2; }
    else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }
    F = nfft/2 + 1;
    sc = 1.0f/nfft;

    //Initialize
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_real((size_t)nfft);
    fplan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    iplan = fftwf_plan_r2r_1d(nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fplan || !iplan) { fprintf(stderr,"error in lin_pred_s: problem creating fftw plan"); return 1; }
    if (!(AStmp=(float *)malloc((size_t)(P-1)*sizeof(float)))) { fprintf(stderr,"error in lin_pred_s: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_scopy(C,&o,0,&AS[0],P);
            for (c=0; c<C; c++)
            {
                cblas_scopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_scopy(R,&X[c*R],1,&X1[0],1);
                if (mean0)
                {
                    m = cblas_sdot(R,&X1[0],1,&o,0) / R;
                    cblas_saxpy(R,-m,&o,0,&X1[0],1);
                }
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                AS[c*P+1] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_sscal(P,sc,&X1[0],1);
                V[c] = fmaf(X1[1],g,X1[0]);
                for (p=2; p<P; p++)
                {
                    g = X1[p];
                    for (q=1; q<p; q++) { g = fmaf(X1[q],AS[p-q+c*P],g); }
                    AS[p+c*P] = g = -g/V[c];
                    for (q=1; q<p; q++) { AS[q+c*P] = fmaf(g,AStmp[p-q-1],AS[q+c*P]); }
                    cblas_scopy(p,&AS[1+c*P],1,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
            }
        }
        else
        {
            cblas_scopy(C,&o,0,&AS[0],1);
            for (c=0; c<C; c++)
            {
                cblas_scopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_scopy(R,&X[c],C,&X1[0],1);
                if (mean0)
                {
                    m = cblas_sdot(R,&X1[0],1,&o,0) / R;
                    cblas_saxpy(R,-m,&o,0,&X1[0],1);
                }
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                AS[c+C] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_sscal(P,sc,&X1[0],1);
                V[c] = fmaf(X1[1],g,X1[0]);
                for (p=2; p<P; p++)
                {
                    g = X1[p];
                    for (q=1; q<p; q++) { g = fmaf(X1[q],AS[c+(p-q)*C],g); }
                    AS[c+p*C] = g = -g/V[c];
                    for (q=1; q<p; q++) { AS[c+q*C] = fmaf(g,AStmp[p-q-1],AS[c+q*C]); }
                    cblas_scopy(p,&AS[c+C],C,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy(R,&o,0,&AS[0],1);
            for (r=0; r<R; r++)
            {
                cblas_scopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_scopy(C,&X[r],R,&X1[0],1);
                if (mean0)
                {
                    m = cblas_sdot(C,&X1[0],1,&o,0) / C;
                    cblas_saxpy(C,-m,&o,0,&X1[0],1);
                }
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                AS[r+R] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_sscal(P,sc,&X1[0],1);
                V[r] = fmaf(X1[1],g,X1[0]);
                for (p=2; p<P; p++)
                {
                    g = X1[p];
                    for (q=1; q<p; q++) { g = fmaf(X1[q],AS[r+(p-q)*R],g); }
                    AS[r+p*R] = g = -g/V[r];
                    for (q=1; q<p; q++) { AS[r+q*R] = fmaf(g,AStmp[p-q-1],AS[r+q*R]); }
                    cblas_scopy(p,&AS[r+R],R,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
            }
        }
        else
        {
            cblas_scopy(R,&o,0,&AS[0],P);
            for (r=0; r<R; r++)
            {
                cblas_scopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_scopy(C,&X[r*C],1,&X1[0],1);
                if (mean0)
                {
                    m = cblas_sdot(C,&X1[0],1,&o,0) / C;
                    cblas_saxpy(C,-m,&o,0,&X1[0],1);
                }
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                AS[r*P+1] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_sscal(P,sc,&X1[0],1);
                V[r] = fmaf(X1[1],g,X1[0]);
                for (p=2; p<P; p++)
                {
                    g = X1[p];
                    for (q=1; q<p; q++) { g = fmaf(X1[q],AS[p-q+r*P],g); }
                    AS[p+r*P] = g = -g/V[r];
                    for (q=1; q<p; q++) { AS[q+r*P] = fmaf(g,AStmp[p-q-1],AS[q+r*P]); }
                    cblas_scopy(p,&AS[1+r*P],1,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in lin_pred_s: dim must be 0 or 1.\n"); return 1;
    }

    fftwf_destroy_plan(fplan); fftwf_destroy_plan(iplan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int lin_pred_d (double *AS, double *V, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0)
{
    const double z = 0.0, o = 1.0;
    double m, g, sc;
    int r, c, f, nfft, F, p, q;
    double *X1, *Y1, *AStmp;
    fftw_plan fplan, iplan;

    //Checks
    if (R<1) { fprintf(stderr,"error in lin_pred_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lin_pred_d: ncols X must be positive\n"); return 1; }
    if (dim==0 && P>R) { fprintf(stderr,"error in lin_pred_d: P must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && P>C) { fprintf(stderr,"error in lin_pred_d: P must be < ncols X for dim==1\n"); return 1; }

    //Get nfft
    nfft = (dim==0) ? R+P : C+P;
    if (nfft>16384) { nfft += nfft%2; }
    else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }
    F = nfft/2 + 1;
    sc = 1.0/nfft;

    //Initialize
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_real((size_t)nfft);
    fplan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    iplan = fftw_plan_r2r_1d(nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fplan || !iplan) { fprintf(stderr,"error in lin_pred_d: problem creating fftw plan"); return 1; }
    if (!(AStmp=(double *)malloc((size_t)(P-1)*sizeof(double)))) { fprintf(stderr,"error in lin_pred_d: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_dcopy(C,&o,0,&AS[0],P);
            for (c=0; c<C; c++)
            {
                cblas_dcopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_dcopy(R,&X[c*R],1,&X1[0],1);
                if (mean0)
                {
                    m = cblas_ddot(R,&X1[0],1,&o,0) / R;
                    cblas_daxpy(R,-m,&o,0,&X1[0],1);
                }
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                AS[c*P+1] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_dscal(P,sc,&X1[0],1);
                V[c] = fma(X1[1],g,X1[0]);
                for (p=2; p<P; p++)
                {
                    g = X1[p];
                    for (q=1; q<p; q++) { g = fma(X1[q],AS[p-q+c*P],g); }
                    AS[p+c*P] = g = -g/V[c];
                    for (q=1; q<p; q++) { AS[q+c*P] = fma(g,AStmp[p-q-1],AS[q+c*P]); }
                    cblas_dcopy(p,&AS[1+c*P],1,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
            }
        }
        else
        {
            cblas_dcopy(C,&o,0,&AS[0],1);
            for (c=0; c<C; c++)
            {
                cblas_dcopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_dcopy(R,&X[c],C,&X1[0],1);
                if (mean0)
                {
                    m = cblas_ddot(R,&X1[0],1,&o,0) / R;
                    cblas_daxpy(R,-m,&o,0,&X1[0],1);
                }
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                AS[c+C] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_dscal(P,sc,&X1[0],1);
                V[c] = fma(X1[1],g,X1[0]);
                for (p=2; p<P; p++)
                {
                    g = X1[p];
                    for (q=1; q<p; q++) { g = fma(X1[q],AS[c+(p-q)*C],g); }
                    AS[c+p*C] = g = -g/V[c];
                    for (q=1; q<p; q++) { AS[c+q*C] = fma(g,AStmp[p-q-1],AS[c+q*C]); }
                    cblas_dcopy(p,&AS[c+C],C,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy(R,&o,0,&AS[0],1);
            for (r=0; r<R; r++)
            {
                cblas_dcopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_dcopy(C,&X[r],R,&X1[0],1);
                if (mean0)
                {
                    m = cblas_ddot(C,&X1[0],1,&o,0) / C;
                    cblas_daxpy(C,-m,&o,0,&X1[0],1);
                }
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                AS[r+R] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_dscal(P,sc,&X1[0],1);
                V[r] = fma(X1[1],g,X1[0]);
                for (p=2; p<P; p++)
                {
                    g = X1[p];
                    for (q=1; q<p; q++) { g = fma(X1[q],AS[r+(p-q)*R],g); }
                    AS[r+p*R] = g = -g/V[r];
                    for (q=1; q<p; q++) { AS[r+q*R] = fma(g,AStmp[p-q-1],AS[r+q*R]); }
                    cblas_dcopy(p,&AS[r+R],R,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
            }
        }
        else
        {
            cblas_dcopy(R,&o,0,&AS[0],P);
            for (r=0; r<R; r++)
            {
                cblas_dcopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_dcopy(C,&X[r*C],1,&X1[0],1);
                if (mean0)
                {
                    m = cblas_ddot(C,&X1[0],1,&o,0) / C;
                    cblas_daxpy(C,-m,&o,0,&X1[0],1);
                }
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                AS[r*P+1] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_dscal(P,sc,&X1[0],1);
                V[r] = fma(X1[1],g,X1[0]);
                for (p=2; p<P; p++)
                {
                    g = X1[p];
                    for (q=1; q<p; q++) { g = fma(X1[q],AS[p-q+r*P],g); }
                    AS[p+r*P] = g = -g/V[r];
                    for (q=1; q<p; q++) { AS[q+r*P] = fma(g,AStmp[p-q-1],AS[q+r*P]); }
                    cblas_dcopy(p,&AS[1+r*P],1,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in lin_pred_d: dim must be 0 or 1.\n"); return 1;
    }

    fftw_destroy_plan(fplan); fftw_destroy_plan(iplan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
#endif

