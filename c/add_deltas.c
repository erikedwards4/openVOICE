//Gets deltas (first order differences) of X,
//which is assumed to be pre-allocated of double size
//so that the appended rows/cols will hold the deltas.
//Thus, for dim==0, C must be even, with C/2 the original ncols of X.
//And,  for dim==1, R must be even, with R/2 the original nrows of X.

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

int add_deltas_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const int N);
int add_deltas_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const int N);


int add_deltas_s (float *X, const char iscolmajor, const int R, const int C, const int dim, const int N)
{
    const float z = 0.0f;
    const int No = R*C/2;
    int r, c, n;
    float sc = 1.0f;

    //Checks
    if (R<1) { fprintf(stderr,"error in add_deltas_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in add_deltas_s: C (ncols X) must be positive\n"); return 1; }
    if (N<1) { fprintf(stderr,"error in add_deltas_s: N (delta winlength) must be positive\n"); return 1; }
    if (dim==0 && C%2!=0) { fprintf(stderr,"error in add_deltas_s: C (ncols X) must be even for dim==0\n"); return 1; }
    if (dim==1 && R%2!=0) { fprintf(stderr,"error in add_deltas_s: R (nrows X) must be even for dim==1\n"); return 1; }

    //Get sc (normalizer)
    for (n=2; n<=N; n++) { sc += n*n; }
    sc = 0.5f/sc;

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_scopy(No,&z,0,&X[No],1);
            for (c=0; c<C/2; c++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[c*R],0,&X[No+c*R],1);        //beg edge samps
                    cblas_saxpy(R-n,-sc*n,&X[c*R],1,&X[No+c*R+n],1);    //past samps
                    cblas_saxpy(R-n,sc*n,&X[c*R+n],1,&X[No+c*R],1);     //future samps
                    cblas_saxpy(n,sc*n,&X[c*R+R-1],0,&X[No+c*R+R-n],1); //end edge samps
                }
            }
        }
        else
        {
            for (c=0; c<C/2; c++)
            {
                cblas_scopy(R,&z,0,&X[C/2+c],C);
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[c],0,&X[C/2+c],C);                //beg edge samps
                    cblas_saxpy(R-n,-sc*n,&X[c],C,&X[C/2+c+n*C],C);          //past samps
                    cblas_saxpy(R-n,sc*n,&X[c+n*C],C,&X[C/2+c],C);           //future samps
                    cblas_saxpy(n,sc*n,&X[c+C*(R-1)],0,&X[C/2+c+C*(R-n)],C); //end edge samps
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R/2; r++)
            {
                cblas_scopy(C,&z,0,&X[R/2+r],R);
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[r],0,&X[R/2+r],R);                //beg edge samps
                    cblas_saxpy(C-n,-sc*n,&X[r],R,&X[R/2+r+n*R],R);          //past samps
                    cblas_saxpy(C-n,sc*n,&X[r+n*R],R,&X[R/2+r],R);           //future samps
                    cblas_saxpy(n,sc*n,&X[r+R*(C-1)],0,&X[R/2+r+R*(C-n)],R); //end edge samps
                }
            }
        }
        else
        {
            cblas_scopy(No,&z,0,&X[No],1);
            for (r=0; r<R/2; r++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_saxpy(n,-sc*n,&X[r*C],0,&X[No+r*C],1);        //beg edge samps
                    cblas_saxpy(C-n,-sc*n,&X[r*C],1,&X[No+r*C+n],1);    //past samps
                    cblas_saxpy(C-n,sc*n,&X[r*C+n],1,&X[No+r*C],1);     //future samps
                    cblas_saxpy(n,sc*n,&X[r*C+C-1],0,&X[No+r*C+C-n],1); //end edge samps
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in add_deltas_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int add_deltas_d (double *X, const char iscolmajor, const int R, const int C, const int dim, const int N)
{
    const double z = 0.0;
    const int No = R*C/2;
    int r, c, n;
    double sc = 1.0;

    //Checks
    if (R<1) { fprintf(stderr,"error in add_deltas_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in add_deltas_d: C (ncols X) must be positive\n"); return 1; }
    if (N<1) { fprintf(stderr,"error in add_deltas_d: N (delta winlength) must be positive\n"); return 1; }
    if (dim==0 && C%2!=0) { fprintf(stderr,"error in add_deltas_d: C (ncols X) must be even for dim==0\n"); return 1; }
    if (dim==1 && R%2!=0) { fprintf(stderr,"error in add_deltas_d: R (nrows X) must be even for dim==1\n"); return 1; }

    //Get sc (normalizer)
    for (n=2; n<=N; n++) { sc += n*n; }
    sc = 0.5/sc;

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_dcopy(No,&z,0,&X[No],1);
            for (c=0; c<C/2; c++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[c*R],0,&X[No+c*R],1);        //beg edge samps
                    cblas_daxpy(R-n,-sc*n,&X[c*R],1,&X[No+c*R+n],1);    //past samps
                    cblas_daxpy(R-n,sc*n,&X[c*R+n],1,&X[No+c*R],1);     //future samps
                    cblas_daxpy(n,sc*n,&X[c*R+R-1],0,&X[No+c*R+R-n],1); //end edge samps
                }
            }
        }
        else
        {
            for (c=0; c<C/2; c++)
            {
                cblas_dcopy(R,&z,0,&X[C/2+c],C);
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[c],0,&X[C/2+c],C);                //beg edge samps
                    cblas_daxpy(R-n,-sc*n,&X[c],C,&X[C/2+c+n*C],C);          //past samps
                    cblas_daxpy(R-n,sc*n,&X[c+n*C],C,&X[C/2+c],C);           //future samps
                    cblas_daxpy(n,sc*n,&X[c+C*(R-1)],0,&X[C/2+c+C*(R-n)],C); //end edge samps
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R/2; r++)
            {
                cblas_dcopy(C,&z,0,&X[R/2+r],R);
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[r],0,&X[R/2+r],R);                //beg edge samps
                    cblas_daxpy(C-n,-sc*n,&X[r],R,&X[R/2+r+n*R],R);          //past samps
                    cblas_daxpy(C-n,sc*n,&X[r+n*R],R,&X[R/2+r],R);           //future samps
                    cblas_daxpy(n,sc*n,&X[r+R*(C-1)],0,&X[R/2+r+R*(C-n)],R); //end edge samps
                }
            }
        }
        else
        {
            cblas_dcopy(No,&z,0,&X[No],1);
            for (r=0; r<R/2; r++)
            {
                for (n=1; n<=N; n++)
                {
                    cblas_daxpy(n,-sc*n,&X[r*C],0,&X[No+r*C],1);        //beg edge samps
                    cblas_daxpy(C-n,-sc*n,&X[r*C],1,&X[No+r*C+n],1);    //past samps
                    cblas_daxpy(C-n,sc*n,&X[r*C+n],1,&X[No+r*C],1);     //future samps
                    cblas_daxpy(n,sc*n,&X[r*C+C-1],0,&X[No+r*C+C-n],1); //end edge samps
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in add_deltas_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

