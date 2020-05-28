//Gets deltas (1st order) and delta-deltas (2nd order) differences of X.
//Only the deltas and delta-deltas are output in Y.
//Thus, Y must be pre-allocated to have twice the size of X.

//I implement this just like FIR for speed and shorter code,
//except non-causal and mid-sample of B is 0,
//so I don't explicitly make B (e.g., B[n] just equals sc*n).

//Note that this may treat edge samples differently than other code.
//But this could be changed with a few lines of additional code here,
//without changing the super-efficient FIR implementation.
//That is, since out-of-range samps were assumed to be 0, nothing was
//added to Y for them, so can just add something later.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int get_delta_deltas_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int N);
int get_delta_deltas_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int N);


int get_delta_deltas_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int N)
{
    const float z = 0.0f;
    int r, c, n;
    float sc = 1.0f;

    //Checks
    if (R<1) { fprintf(stderr,"error in get_delta_deltas_s: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in get_delta_deltas_s: C (ncols Y) must be positive\n"); return 1; }
    if (N<1) { fprintf(stderr,"error in get_delta_deltas_s: N (delta winlength) must be positive\n"); return 1; }

    //Get sc (normalizer)
    for (n=2; n<=N; n++) { sc += n*n; }
    sc = 0.5f/sc;

    //Initialize Y
    cblas_scopy(2*R*C,&z,0,&Y[0],1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[c*R],0,&Y[c*R],1);            //beg edge samps
                    cblas_saxpy(R-n,-sc*n,&X[c*R],1,&Y[c*R+n],1);        //past samps
                    cblas_saxpy(R-n,sc*n,&X[c*R+n],1,&Y[c*R],1);         //future samps
                    cblas_saxpy(n,sc*n,&X[c*R+R-1],0,&Y[c*R+R-n],1);     //end edge samps
                }
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&Y[c*R],0,&Y[R*C+c*R],1);        //beg edge samps
                    cblas_saxpy(R-n,-sc*n,&Y[c*R],1,&Y[R*C+c*R+n],1);    //past samps
                    cblas_saxpy(R-n,sc*n,&Y[c*R+n],1,&Y[R*C+c*R],1);     //future samps
                    cblas_saxpy(n,sc*n,&Y[c*R+R-1],0,&Y[R*C+c*R+R-n],1); //end edge samps
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[c],0,&Y[c],2*C);                      //beg edge samps
                    cblas_saxpy(R-n,-sc*n,&X[c],C,&Y[c+n*2*C],2*C);              //past samps
                    cblas_saxpy(R-n,sc*n,&X[c+n*C],C,&Y[c],2*C);                 //future samps
                    cblas_saxpy(n,sc*n,&X[c+C*(R-1)],0,&Y[c+2*C*(R-n)],2*C);     //end edge samps
                }
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&Y[c],0,&Y[C+c],2*C);                    //beg edge samps
                    cblas_saxpy(R-n,-sc*n,&Y[c],2*C,&Y[C+c+n*2*C],2*C);          //past samps
                    cblas_saxpy(R-n,sc*n,&Y[c+n*2*C],2*C,&Y[C+c],2*C);           //future samps
                    cblas_saxpy(n,sc*n,&Y[c+2*C*(R-1)],0,&Y[C+c+2*C*(R-n)],2*C); //end edge samps
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[r],0,&Y[r],2*R);                      //beg edge samps
                    cblas_saxpy(C-n,-sc*n,&X[r],R,&Y[r+n*2*R],2*R);              //past samps
                    cblas_saxpy(C-n,sc*n,&X[r+n*R],R,&Y[r],2*R);                 //future samps
                    cblas_saxpy(n,sc*n,&X[r+R*(C-1)],0,&Y[r+2*R*(C-n)],2*R);     //end edge samps
                }
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&Y[r],0,&Y[R+r],2*R);                    //beg edge samps
                    cblas_saxpy(C-n,-sc*n,&Y[r],2*R,&Y[R+r+n*2*R],2*R);          //past samps
                    cblas_saxpy(C-n,sc*n,&Y[r+n*2*R],2*R,&Y[R+r],2*R);           //future samps
                    cblas_saxpy(n,sc*n,&Y[r+2*R*(C-1)],0,&Y[R+r+2*R*(C-n)],2*R); //end edge samps
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[r*C],0,&Y[r*C],1);            //beg edge samps
                    cblas_saxpy(C-n,-sc*n,&X[r*C],1,&Y[r*C+n],1);        //past samps
                    cblas_saxpy(C-n,sc*n,&X[r*C+n],1,&Y[r*C],1);         //future samps
                    cblas_saxpy(n,sc*n,&X[r*C+C-1],0,&Y[r*C+C-n],1);     //end edge samps
                }
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&Y[r*C],0,&Y[R*C+r*C],1);        //beg edge samps
                    cblas_saxpy(C-n,-sc*n,&Y[r*C],1,&Y[R*C+r*C+n],1);    //past samps
                    cblas_saxpy(C-n,sc*n,&Y[r*C+n],1,&Y[R*C+r*C],1);     //future samps
                    cblas_saxpy(n,sc*n,&Y[r*C+C-1],0,&Y[R*C+r*C+C-n],1); //end edge samps
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in get_delta_deltas_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int get_delta_deltas_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int N)
{
    const double z = 0.0;
    int r, c, n;
    double sc = 1.0;

    //Checks
    if (R<1) { fprintf(stderr,"error in get_delta_deltas_d: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in get_delta_deltas_d: C (ncols Y) must be positive\n"); return 1; }
    if (N<1) { fprintf(stderr,"error in get_delta_deltas_d: N (delta winlength) must be positive\n"); return 1; }

    //Get sc (normalizer)
    for (n=2; n<=N; n++) { sc += n*n; }
    sc = 0.5/sc;

    //Initialize Y
    cblas_dcopy(2*R*C,&z,0,&Y[0],1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[c*R],0,&Y[c*R],1);            //beg edge samps
                    cblas_daxpy(R-n,-sc*n,&X[c*R],1,&Y[c*R+n],1);        //past samps
                    cblas_daxpy(R-n,sc*n,&X[c*R+n],1,&Y[c*R],1);         //future samps
                    cblas_daxpy(n,sc*n,&X[c*R+R-1],0,&Y[c*R+R-n],1);     //end edge samps
                }
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&Y[c*R],0,&Y[R*C+c*R],1);        //beg edge samps
                    cblas_daxpy(R-n,-sc*n,&Y[c*R],1,&Y[R*C+c*R+n],1);    //past samps
                    cblas_daxpy(R-n,sc*n,&Y[c*R+n],1,&Y[R*C+c*R],1);     //future samps
                    cblas_daxpy(n,sc*n,&Y[c*R+R-1],0,&Y[R*C+c*R+R-n],1); //end edge samps
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[c],0,&Y[c],2*C);                      //beg edge samps
                    cblas_daxpy(R-n,-sc*n,&X[c],C,&Y[c+n*2*C],2*C);              //past samps
                    cblas_daxpy(R-n,sc*n,&X[c+n*C],C,&Y[c],2*C);                 //future samps
                    cblas_daxpy(n,sc*n,&X[c+C*(R-1)],0,&Y[c+2*C*(R-n)],2*C);     //end edge samps
                }
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&Y[c],0,&Y[C+c],2*C);                    //beg edge samps
                    cblas_daxpy(R-n,-sc*n,&Y[c],2*C,&Y[C+c+n*2*C],2*C);          //past samps
                    cblas_daxpy(R-n,sc*n,&Y[c+n*2*C],2*C,&Y[C+c],2*C);           //future samps
                    cblas_daxpy(n,sc*n,&Y[c+2*C*(R-1)],0,&Y[C+c+2*C*(R-n)],2*C); //end edge samps
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[r],0,&Y[r],2*R);                      //beg edge samps
                    cblas_daxpy(C-n,-sc*n,&X[r],R,&Y[r+n*2*R],2*R);              //past samps
                    cblas_daxpy(C-n,sc*n,&X[r+n*R],R,&Y[r],2*R);                 //future samps
                    cblas_daxpy(n,sc*n,&X[r+R*(C-1)],0,&Y[r+2*R*(C-n)],2*R);     //end edge samps
                }
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&Y[r],0,&Y[R+r],2*R);                    //beg edge samps
                    cblas_daxpy(C-n,-sc*n,&Y[r],2*R,&Y[R+r+n*2*R],2*R);          //past samps
                    cblas_daxpy(C-n,sc*n,&Y[r+n*2*R],2*R,&Y[R+r],2*R);           //future samps
                    cblas_daxpy(n,sc*n,&Y[r+2*R*(C-1)],0,&Y[R+r+2*R*(C-n)],2*R); //end edge samps
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[r*C],0,&Y[r*C],1);            //beg edge samps
                    cblas_daxpy(C-n,-sc*n,&X[r*C],1,&Y[r*C+n],1);        //past samps
                    cblas_daxpy(C-n,sc*n,&X[r*C+n],1,&Y[r*C],1);         //future samps
                    cblas_daxpy(n,sc*n,&X[r*C+C-1],0,&Y[r*C+C-n],1);     //end edge samps
                }
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&Y[r*C],0,&Y[R*C+r*C],1);        //beg edge samps
                    cblas_daxpy(C-n,-sc*n,&Y[r*C],1,&Y[R*C+r*C+n],1);    //past samps
                    cblas_daxpy(C-n,sc*n,&Y[r*C+n],1,&Y[R*C+r*C],1);     //future samps
                    cblas_daxpy(n,sc*n,&Y[r*C+C-1],0,&Y[R*C+r*C+C-n],1); //end edge samps
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in get_delta_deltas_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

