//Gets deltas (first order differences) of X.
//Only the deltas are output in Y.

//I implement this just like FIR for speed and shorter code,
//except non-causal and mid-sample of B is 0,
//so I don't explicitly make B (e.g., B[n] just equals sc*n).

//It seems that Kaldi and others use the same method for edge samples,
//which is numpy.edge, so that is used here.
//Comment out those lines to set out-of-range samps to 0.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int get_deltas_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int N);
int get_deltas_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int N);


int get_deltas_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int N)
{
    const float z = 0.0f;
    int r, c, n;
    float sc = 1.0f;

    //Checks
    if (R<1) { fprintf(stderr,"error in get_deltas_s: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in get_deltas_s: C (ncols Y) must be positive\n"); return 1; }
    if (N<1) { fprintf(stderr,"error in get_deltas_s: N (delta winlength) must be positive\n"); return 1; }

    //Get sc (normalizer)
    for (n=2; n<=N; n++) { sc += n*n; }
    sc = 0.5f/sc;

    //Initialize Y
    cblas_scopy(R*C,&z,0,&Y[0],1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[c*R],0,&Y[c*R],1);        //beg edge samps
                    cblas_saxpy(R-n,-sc*n,&X[c*R],1,&Y[c*R+n],1);    //past samps
                    cblas_saxpy(R-n,sc*n,&X[c*R+n],1,&Y[c*R],1);     //future samps
                    cblas_saxpy(n,sc*n,&X[c*R+R-1],0,&Y[c*R+R-n],1); //end edge samps
                }
            }
        }
        else
        {
            for (n=1; n<=N; n++)
            {
                cblas_saxpy(C*(R-n),-sc*n,&X[0],1,&Y[n*C],1); //past samps
                cblas_saxpy(C*(R-n),sc*n,&X[n*C],1,&Y[0],1);  //future samps
                for (c=0; c<C; c++)
                {
                    cblas_saxpy(n,-sc*n,&X[c],0,&Y[c],C);                //beg edge samps
                    cblas_saxpy(n,sc*n,&X[c+C*(R-1)],0,&Y[c+C*(R-n)],C); //end edge samps
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (n=1; n<=N; n++)
            {
                cblas_saxpy(R*(C-n),-sc*n,&X[0],1,&Y[n*R],1); //past samps
                cblas_saxpy(R*(C-n),sc*n,&X[n*R],1,&Y[0],1);  //future samps
                for (r=0; r<R; r++)
                {
                    cblas_saxpy(n,-sc*n,&X[r],0,&Y[r],R);                //beg edge samps
                    cblas_saxpy(n,sc*n,&X[r+R*(C-1)],0,&Y[r+R*(C-n)],R); //end edge samps
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[r*C],0,&Y[r*C],1);        //beg edge samps
                    cblas_saxpy(C-n,-sc*n,&X[r*C],1,&Y[r*C+n],1);    //past samps
                    cblas_saxpy(C-n,sc*n,&X[r*C+n],1,&Y[r*C],1);     //future samps
                    cblas_saxpy(n,sc*n,&X[r*C+C-1],0,&Y[r*C+C-n],1); //end edge samps
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in get_deltas_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int get_deltas_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int N)
{
    const double z = 0.0;
    int r, c, n;
    double sc = 1.0;

    //Checks
    if (R<1) { fprintf(stderr,"error in get_deltas_d: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in get_deltas_d: C (ncols Y) must be positive\n"); return 1; }
    if (N<1) { fprintf(stderr,"error in get_deltas_d: N (delta winlength) must be positive\n"); return 1; }

    //Get sc (normalizer)
    for (n=2; n<=N; n++) { sc += n*n; }
    sc = 0.5/sc;

    //Initialize Y
    cblas_dcopy(R*C,&z,0,&Y[0],1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[c*R],0,&Y[c*R],1);        //beg edge samps
                    cblas_daxpy(R-n,-sc*n,&X[c*R],1,&Y[c*R+n],1);    //past samps
                    cblas_daxpy(R-n,sc*n,&X[c*R+n],1,&Y[c*R],1);     //future samps
                    cblas_daxpy(n,sc*n,&X[c*R+R-1],0,&Y[c*R+R-n],1); //end edge samps
                }
            }
        }
        else
        {
            for (n=1; n<=N; n++)
            {
                cblas_daxpy(C*(R-n),-sc*n,&X[0],1,&Y[n*C],1); //past samps
                cblas_daxpy(C*(R-n),sc*n,&X[n*C],1,&Y[0],1);  //future samps
                for (c=0; c<C; c++)
                {
                    cblas_daxpy(n,-sc*n,&X[c],0,&Y[c],C);                //beg edge samps
                    cblas_daxpy(n,sc*n,&X[c+C*(R-1)],0,&Y[c+C*(R-n)],C); //end edge samps
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (n=1; n<=N; n++)
            {
                cblas_daxpy(R*(C-n),-sc*n,&X[0],1,&Y[n*R],1); //past samps
                cblas_daxpy(R*(C-n),sc*n,&X[n*R],1,&Y[0],1);  //future samps
                for (r=0; r<R; r++)
                {
                    cblas_daxpy(n,-sc*n,&X[r],0,&Y[r],R);                //beg edge samps
                    cblas_daxpy(n,sc*n,&X[r+R*(C-1)],0,&Y[r+R*(C-n)],R); //end edge samps
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[r*C],0,&Y[r*C],1);        //beg edge samps
                    cblas_daxpy(C-n,-sc*n,&X[r*C],1,&Y[r*C+n],1);    //past samps
                    cblas_daxpy(C-n,sc*n,&X[r*C+n],1,&Y[r*C],1);     //future samps
                    cblas_daxpy(n,sc*n,&X[r*C+C-1],0,&Y[r*C+C-n],1); //end edge samps
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in get_deltas_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

