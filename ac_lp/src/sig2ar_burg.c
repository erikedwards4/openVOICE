//This gets polynomial coeffs for each row/col of X using the Burg method.

//An opt mean0 is added to zero the mean of each row/col of X first.
//In this case, this is mean0 -> autocov_fft -> sig2ar_burg.

//To test compile:
//gcc -c sig2ar_burg.c -O2 -std=c99 -Wall -Wextra
//clang -c sig2ar_burg.c -O2 -std=c99 -Weverything
//g++ -c sig2ar_burg.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c sig2ar_burg.c -O2 -std=c++11 -Weverything -Wno-old-style-cast

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int sig2ar_burg_s (float *Y, float *V, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0)
{
    const float o = 1.0f;
    const int N = (dim==0) ? R : C;
    int r, c, p, q;
    float m, g, b2, f2;
    float *b, *f, *oldf, *AStmp;

    //Checks
    if (R<1) { fprintf(stderr,"error in sig2ar_burg_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sig2ar_burg_s: ncols X must be positive\n"); return 1; }
    if (dim==0 && P>=R) { fprintf(stderr,"error in sig2ar_burg_s: P must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && P>=C) { fprintf(stderr,"error in sig2ar_burg_s: P must be < ncols X for dim==1\n"); return 1; }

    //Allocate
    if (!(b=(float *)malloc((size_t)(N)*sizeof(float)))) { fprintf(stderr,"error in sig2ar_burg_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(f=(float *)malloc((size_t)(N)*sizeof(float)))) { fprintf(stderr,"error in sig2ar_burg_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(oldf=(float *)malloc((size_t)(N)*sizeof(float)))) { fprintf(stderr,"error in sig2ar_burg_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AStmp=(float *)malloc((size_t)(P)*sizeof(float)))) { fprintf(stderr,"error in sig2ar_burg_s: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,b,1);
                if (mean0)
                {
                    m = cblas_sdot(R,b,1,&o,0) / R;
                    cblas_saxpy(R,-m,&o,0,b,1);
                }
                cblas_scopy(R,b,1,f,1);
                V[c] = cblas_sdot(R,b,1,f,1) / R;
                for (p=1; p<=P; p++)
                {
                    b2 = cblas_sdot(R-p,b,1,b,1);
                    f2 = cblas_sdot(R-p,&f[p],1,&f[p],1);
                    Y[c*P+p-1] = g = (-2.0f/(b2+f2)) * cblas_sdot(R-p,b,1,&f[p],1);
                    V[c] *= fmaf(g,-g,1.0f);
                    if (p>1) { for (q=0; q<p; q++) { Y[c*P+q] += g * AStmp[p-q-2]; } }
                    if (p<P)
                    {
                        cblas_scopy(p,&Y[c*P],1,AStmp,1);
                        cblas_scopy(N-p,&f[p],1,&oldf[p],1);
                        cblas_saxpy(N-p,g,b,1,&f[p],1);
                        cblas_saxpy(N-p,g,&oldf[p],1,b,1);
                    }
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c],C,b,1);
                if (mean0)
                {
                    m = cblas_sdot(R,b,1,&o,0) / R;
                    cblas_saxpy(R,-m,&o,0,b,1);
                }
                cblas_scopy(R,b,1,f,1);
                V[c] = cblas_sdot(R,b,1,f,1) / R;
                for (p=1; p<=P; p++)
                {
                    b2 = cblas_sdot(R-p,b,1,b,1);
                    f2 = cblas_sdot(R-p,&f[p],1,&f[p],1);
                    Y[c+(p-1)*C] = g = (-2.0f/(b2+f2)) * cblas_sdot(R-p,b,1,&f[p],1);
                    V[c] *= fmaf(g,-g,1.0f);
                    if (p>1) { for (q=0; q<p; q++) { Y[c+q*C] += g * AStmp[p-q-2]; } }
                    if (p<P)
                    {
                        cblas_scopy(p,&Y[c],C,AStmp,1);
                        cblas_scopy(N-p,&f[p],1,&oldf[p],1);
                        cblas_saxpy(N-p,g,b,1,&f[p],1);
                        cblas_saxpy(N-p,g,&oldf[p],1,b,1);
                    }
                }
            }
        }
        cblas_sscal(P*C,-1.0f,Y,1);
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r],R,b,1);
                if (mean0)
                {
                    m = cblas_sdot(C,b,1,&o,0) / C;
                    cblas_saxpy(C,-m,&o,0,b,1);
                }
                cblas_scopy(C,b,1,f,1);
                V[r] = cblas_sdot(C,b,1,f,1) / C;
                for (p=1; p<=P; p++)
                {
                    b2 = cblas_sdot(C-p,b,1,b,1);
                    f2 = cblas_sdot(C-p,&f[p],1,&f[p],1);
                    Y[r+(p-1)*R] = g = (-2.0f/(b2+f2)) * cblas_sdot(C-p,b,1,&f[p],1);
                    V[r] *= fmaf(g,-g,1.0f);
                    if (p>1) { for (q=0; q<p; q++) { Y[r+q*R] += g * AStmp[p-q-2]; } }
                    if (p<P)
                    {
                        cblas_scopy(p,&Y[r],R,AStmp,1);
                        cblas_scopy(N-p,&f[p],1,&oldf[p],1);
                        cblas_saxpy(N-p,g,b,1,&f[p],1);
                        cblas_saxpy(N-p,g,&oldf[p],1,b,1);
                    }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,b,1);
                if (mean0)
                {
                    m = cblas_sdot(C,b,1,&o,0) / C;
                    cblas_saxpy(C,-m,&o,0,b,1);
                }
                cblas_scopy(C,b,1,f,1);
                V[r] = cblas_sdot(C,b,1,f,1) / C;
                for (p=1; p<=P; p++)
                {
                    b2 = cblas_sdot(C-p,b,1,b,1);
                    f2 = cblas_sdot(C-p,&f[p],1,&f[p],1);
                    Y[r*P+p-1] = g = (-2.0f/(b2+f2)) * cblas_sdot(C-p,b,1,&f[p],1);
                    V[r] *= fmaf(g,-g,1.0f);
                    if (p>1) { for (q=0; q<p; q++) { Y[r*P+q] += g * AStmp[p-q-2]; } }
                    if (p<P)
                    {
                        cblas_scopy(p,&Y[r*P],1,AStmp,1);
                        cblas_scopy(N-p,&f[p],1,&oldf[p],1);
                        cblas_saxpy(N-p,g,b,1,&f[p],1);
                        cblas_saxpy(N-p,g,&oldf[p],1,b,1);
                    }
                }
            }
        }
        cblas_sscal(R*P,-1.0f,Y,1);
    }
    else
    {
        fprintf(stderr,"error in sig2ar_burg_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(b); free(f); free(oldf); free(AStmp);
    return 0;
}


int sig2ar_burg_d (double *Y, double *V, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int P, const char mean0)
{
    const double o = 1.0;
    const int N = (dim==0) ? R : C;
    int r, c, p, q;
    double m, g, b2, f2;
    double *b, *f, *oldf, *AStmp;

    //Checks
    if (R<1) { fprintf(stderr,"error in sig2ar_burg_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sig2ar_burg_d: ncols X must be positive\n"); return 1; }
    if (dim==0 && P>=R) { fprintf(stderr,"error in sig2ar_burg_d: P must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && P>=C) { fprintf(stderr,"error in sig2ar_burg_d: P must be < ncols X for dim==1\n"); return 1; }

    //Allocate
    if (!(b=(double *)malloc((size_t)(N)*sizeof(double)))) { fprintf(stderr,"error in sig2ar_burg_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(f=(double *)malloc((size_t)(N)*sizeof(double)))) { fprintf(stderr,"error in sig2ar_burg_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(oldf=(double *)malloc((size_t)(N)*sizeof(double)))) { fprintf(stderr,"error in sig2ar_burg_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AStmp=(double *)malloc((size_t)(P)*sizeof(double)))) { fprintf(stderr,"error in sig2ar_burg_d: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,b,1);
                if (mean0)
                {
                    m = cblas_ddot(R,b,1,&o,0) / R;
                    cblas_daxpy(R,-m,&o,0,b,1);
                }
                cblas_dcopy(R,b,1,f,1);
                V[c] = cblas_ddot(R,b,1,f,1) / R;
                for (p=1; p<=P; p++)
                {
                    b2 = cblas_ddot(R-p,b,1,b,1);
                    f2 = cblas_ddot(R-p,&f[p],1,&f[p],1);
                    Y[c*P+p-1] = g = (-2.0/(b2+f2)) * cblas_ddot(R-p,b,1,&f[p],1);
                    V[c] *= fma(g,-g,1.0);
                    if (p>1) { for (q=0; q<p; q++) { Y[c*P+q] += g * AStmp[p-q-2]; } }
                    if (p<P)
                    {
                        cblas_dcopy(p,&Y[c*P],1,AStmp,1);
                        cblas_dcopy(N-p,&f[p],1,&oldf[p],1);
                        cblas_daxpy(N-p,g,b,1,&f[p],1);
                        cblas_daxpy(N-p,g,&oldf[p],1,b,1);
                    }
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c],C,b,1);
                if (mean0)
                {
                    m = cblas_ddot(R,b,1,&o,0) / R;
                    cblas_daxpy(R,-m,&o,0,b,1);
                }
                cblas_dcopy(R,b,1,f,1);
                V[c] = cblas_ddot(R,b,1,f,1) / R;
                for (p=1; p<=P; p++)
                {
                    b2 = cblas_ddot(R-p,b,1,b,1);
                    f2 = cblas_ddot(R-p,&f[p],1,&f[p],1);
                    Y[c+(p-1)*C] = g = (-2.0/(b2+f2)) * cblas_ddot(R-p,b,1,&f[p],1);
                    V[c] *= fma(g,-g,1.0);
                    if (p>1) { for (q=0; q<p; q++) { Y[c+q*C] += g * AStmp[p-q-2]; } }
                    if (p<P)
                    {
                        cblas_dcopy(p,&Y[c],C,AStmp,1);
                        cblas_dcopy(N-p,&f[p],1,&oldf[p],1);
                        cblas_daxpy(N-p,g,b,1,&f[p],1);
                        cblas_daxpy(N-p,g,&oldf[p],1,b,1);
                    }
                }
            }
        }
        cblas_dscal(P*C,-1.0,Y,1);
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r],R,b,1);
                if (mean0)
                {
                    m = cblas_ddot(C,b,1,&o,0) / C;
                    cblas_daxpy(C,-m,&o,0,b,1);
                }
                cblas_dcopy(C,b,1,f,1);
                V[r] = cblas_ddot(C,b,1,f,1) / C;
                for (p=1; p<=P; p++)
                {
                    b2 = cblas_ddot(C-p,b,1,b,1);
                    f2 = cblas_ddot(C-p,&f[p],1,&f[p],1);
                    Y[r+(p-1)*R] = g = (-2.0/(b2+f2)) * cblas_ddot(C-p,b,1,&f[p],1);
                    V[r] *= fma(g,-g,1.0);
                    if (p>1) { for (q=0; q<p; q++) { Y[r+q*R] += g * AStmp[p-q-2]; } }
                    if (p<P)
                    {
                        cblas_dcopy(p,&Y[r],R,AStmp,1);
                        cblas_dcopy(N-p,&f[p],1,&oldf[p],1);
                        cblas_daxpy(N-p,g,b,1,&f[p],1);
                        cblas_daxpy(N-p,g,&oldf[p],1,b,1);
                    }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,b,1);
                if (mean0)
                {
                    m = cblas_ddot(C,b,1,&o,0) / C;
                    cblas_daxpy(C,-m,&o,0,b,1);
                }
                cblas_dcopy(C,b,1,f,1);
                V[r] = cblas_ddot(C,b,1,f,1) / C;
                for (p=1; p<=P; p++)
                {
                    b2 = cblas_ddot(C-p,b,1,b,1);
                    f2 = cblas_ddot(C-p,&f[p],1,&f[p],1);
                    Y[r*P+p-1] = g = (-2.0/(b2+f2)) * cblas_ddot(C-p,b,1,&f[p],1);
                    V[r] *= fma(g,-g,1.0);
                    if (p>1) { for (q=0; q<p; q++) { Y[r*P+q] += g * AStmp[p-q-2]; } }
                    if (p<P)
                    {
                        cblas_dcopy(p,&Y[r*P],1,AStmp,1);
                        cblas_dcopy(N-p,&f[p],1,&oldf[p],1);
                        cblas_daxpy(N-p,g,b,1,&f[p],1);
                        cblas_daxpy(N-p,g,&oldf[p],1,b,1);
                    }
                }
            }
        }
        cblas_dscal(R*P,-1.0,Y,1);
    }
    else
    {
        fprintf(stderr,"error in sig2ar_burg_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(b); free(f); free(oldf); free(AStmp);
    return 0;
}


#ifdef __cplusplus
}
#endif

