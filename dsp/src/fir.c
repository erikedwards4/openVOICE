//for (n=L-1; n<N; n++) { Y[n] = cblas_sdot(L,&B[0],1,&X[n-L+1],1); }  //slower even for 1-chan case

//To test compile:
//gcc -c fir.c -O2 -std=c99 -Wall -Wextra
//clang -c fir.c -O2 -std=c99 -Weverything
//g++ -c fir.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c fir.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int fir_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *B, const int L, const int dim, const int stride)
{
    const float z = 0.0f;
    const int Ro = R/stride, Co = C/stride;
    int r, c, l;
    int *ss, *Rs, *Cs;

    //Checks
    if (R<1) { fprintf(stderr,"error in fir_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fir_s: C (ncols X) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in fir_s: L (filter IR length) must be positive\n"); return 1; }
    if (stride<1) { fprintf(stderr,"error in fir_s: stride must be positive\n"); return 1; }

    //Initialize the start-samp offsets and lengths
    if (!(ss=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Rs=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Cs=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_s: problem with malloc. "); perror("malloc"); return 1; }
    for (l=0; l<L; l++) { ss[l] = (stride-l%stride)%stride; Rs[l] = (R-l)/stride; Cs[l] = (C-l)/stride; }
    
    if (dim==0)
    {
        cblas_scopy(Ro*C,&z,0,&Y[0],1);
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (l=0; l<L; l++) { cblas_saxpy(Rs[l],B[l],&X[c*R+ss[l]],stride,&Y[c*Ro+Ro-Rs[l]],1); }
            }
        }
        else
        {
            if (stride==1)
            {
                for (l=0; l<L; l++) { cblas_saxpy(C*(R-l),B[l],&X[0],1,&Y[l*C],1); }
            }
            else
            {
                for (c=0; c<C; c++)
                {
                    for (l=0; l<L; l++) { cblas_saxpy(Rs[l],B[l],&X[c+C*ss[l]],C*stride,&Y[c+C*(Ro-Rs[l])],C); }
                }
            }
        }
    }
    else if (dim==1)
    {
        cblas_scopy(R*Co,&z,0,&Y[0],1);
        if (iscolmajor)
        {
            if (stride==1)
            {
                for (l=0; l<L; l++) { cblas_saxpy(R*(C-l),B[l],&X[0],1,&Y[l*R],1); }
            }
            else
            {
                for (r=0; r<R; r++)
                {
                    for (l=0; l<L; l++) { cblas_saxpy(Cs[l],B[l],&X[r+R*ss[l]],R*stride,&Y[r+R*(Co-Cs[l])],R); }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (l=0; l<L; l++) { cblas_saxpy(Cs[l],B[l],&X[r*C+ss[l]],stride,&Y[r*Co+Co-Cs[l]],1); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fir_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(ss); free(Rs); free(Cs);
    return 0;
}


int fir_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *B, const int L, const int dim, const int stride)
{
    const double z = 0.0;
    const int Ro = R/stride, Co = C/stride;
    int r, c, l;
    int *ss, *Rs, *Cs;

    //Checks
    if (R<1) { fprintf(stderr,"error in fir_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fir_d: C (ncols X) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in fir_d: L (filter IR length) must be positive\n"); return 1; }
    if (stride<1) { fprintf(stderr,"error in fir_d: stride must be positive\n"); return 1; }

    //Initialize the start-samp offsets and lengths
    if (!(ss=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Rs=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Cs=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_d: problem with malloc. "); perror("malloc"); return 1; }
    for (l=0; l<L; l++) { ss[l] = (stride-l%stride)%stride; Rs[l] = (R-l)/stride; Cs[l] = (C-l)/stride; }
    
    if (dim==0)
    {
        cblas_dcopy(Ro*C,&z,0,&Y[0],1);
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (l=0; l<L; l++) { cblas_daxpy(Rs[l],B[l],&X[c*R+ss[l]],stride,&Y[c*Ro+Ro-Rs[l]],1); }
            }
        }
        else
        {
            if (stride==1)
            {
                for (l=0; l<L; l++) { cblas_daxpy(C*(R-l),B[l],&X[0],1,&Y[l*C],1); }
            }
            else
            {
                for (c=0; c<C; c++)
                {
                    for (l=0; l<L; l++) { cblas_daxpy(Rs[l],B[l],&X[c+C*ss[l]],C*stride,&Y[c+C*(Ro-Rs[l])],C); }
                }
            }
        }
    }
    else if (dim==1)
    {
        cblas_dcopy(R*Co,&z,0,&Y[0],1);
        if (iscolmajor)
        {
            if (stride==1)
            {
                for (l=0; l<L; l++) { cblas_daxpy(R*(C-l),B[l],&X[0],1,&Y[l*R],1); }
            }
            else
            {
                for (r=0; r<R; r++)
                {
                    for (l=0; l<L; l++) { cblas_daxpy(Cs[l],B[l],&X[r+R*ss[l]],R*stride,&Y[r+R*(Co-Cs[l])],R); }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (l=0; l<L; l++) { cblas_daxpy(Cs[l],B[l],&X[r*C+ss[l]],stride,&Y[r*Co+Co-Cs[l]],1); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fir_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(ss); free(Rs); free(Cs);
    return 0;
}


int fir_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const float *B, const int L, const int dim, const int stride)
{
    const int Ro = R/stride, Co = C/stride;
    int r, c, l;
    int *ss, *Rs, *Cs;

    //Checks
    if (R<1) { fprintf(stderr,"error in fir_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fir_c: C (ncols X) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in fir_c: L (filter IR length) must be positive\n"); return 1; }
    if (stride<1) { fprintf(stderr,"error in fir_c: stride must be positive\n"); return 1; }

    //Initialize the start-samp offsets and lengths
    if (!(ss=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_c: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Rs=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_c: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Cs=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_c: problem with malloc. "); perror("malloc"); return 1; }
    for (l=0; l<L; l++) { ss[l] = (stride-l%stride)%stride; Rs[l] = (R-l)/stride; Cs[l] = (C-l)/stride; }
    
    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (l=0; l<L; l++) { cblas_caxpy(Rs[l],&B[2*l],&X[2*(c*R+ss[l])],stride,&Y[2*(c*Ro+Ro-Rs[l])],1); }
            }
        }
        else
        {
            if (stride==1)
            {
                for (l=0; l<L; l++) { cblas_caxpy(C*(R-l),&B[2*l],&X[0],1,&Y[2*l*C],1); }
            }
            else
            {
                for (c=0; c<C; c++)
                {
                    for (l=0; l<L; l++) { cblas_caxpy(Rs[l],&B[2*l],&X[2*(c+C*ss[l])],C*stride,&Y[2*(c+C*(Ro-Rs[l]))],C); }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            if (stride==1)
            {
                for (l=0; l<L; l++) { cblas_caxpy(R*(C-l),&B[2*l],&X[0],1,&Y[2*l*R],1); }
            }
            else
            {
                for (r=0; r<R; r++)
                {
                    for (l=0; l<L; l++) { cblas_caxpy(Cs[l],&B[2*l],&X[2*(r+R*ss[l])],R*stride,&Y[2*(r+R*(Co-Cs[l]))],R); }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (l=0; l<L; l++) { cblas_caxpy(Cs[l],&B[2*l],&X[2*(r*C+ss[l])],stride,&Y[2*(r*Co+Co-Cs[l])],1); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fir_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(ss); free(Rs); free(Cs);
    return 0;
}


int fir_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const double *B, const int L, const int dim, const int stride)
{
    const int Ro = R/stride, Co = C/stride;
    int r, c, l;
    int *ss, *Rs, *Cs;

    //Checks
    if (R<1) { fprintf(stderr,"error in fir_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fir_z: C (ncols X) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in fir_z: L (filter IR length) must be positive\n"); return 1; }
    if (stride<1) { fprintf(stderr,"error in fir_z: stride must be positive\n"); return 1; }

    //Initialize the start-samp offsets and lengths
    if (!(ss=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_z: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Rs=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_z: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Cs=(int *)malloc((size_t)L*sizeof(int)))) { fprintf(stderr,"error in fir_z: problem with malloc. "); perror("malloc"); return 1; }
    for (l=0; l<L; l++) { ss[l] = (stride-l%stride)%stride; Rs[l] = (R-l)/stride; Cs[l] = (C-l)/stride; }
    
    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (l=0; l<L; l++) { cblas_zaxpy(Rs[l],&B[2*l],&X[2*(c*R+ss[l])],stride,&Y[2*(c*Ro+Ro-Rs[l])],1); }
            }
        }
        else
        {
            if (stride==1)
            {
                for (l=0; l<L; l++) { cblas_zaxpy(C*(R-l),&B[2*l],&X[0],1,&Y[2*l*C],1); }
            }
            else
            {
                for (c=0; c<C; c++)
                {
                    for (l=0; l<L; l++) { cblas_zaxpy(Rs[l],&B[2*l],&X[2*(c+C*ss[l])],C*stride,&Y[2*(c+C*(Ro-Rs[l]))],C); }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            if (stride==1)
            {
                for (l=0; l<L; l++) { cblas_zaxpy(R*(C-l),&B[2*l],&X[0],1,&Y[2*l*R],1); }
            }
            else
            {
                for (r=0; r<R; r++)
                {
                    for (l=0; l<L; l++) { cblas_zaxpy(Cs[l],&B[2*l],&X[2*(r+R*ss[l])],R*stride,&Y[2*(r+R*(Co-Cs[l]))],R); }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (l=0; l<L; l++) { cblas_zaxpy(Cs[l],&B[2*l],&X[2*(r*C+ss[l])],stride,&Y[2*(r*Co+Co-Cs[l])],1); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fir_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(ss); free(Rs); free(Cs);
    return 0;
}


#ifdef __cplusplus
}
#endif

