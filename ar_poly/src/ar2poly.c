//Gets polynomials from autoregressive (AR) parameters along rows or cols of X.
//If the polynomial is a0 a1 a2..., then the AR coeffs are -a1/a0 -a2/a0...
//Since a0 cannot be recovered from AR coeffs, a0 is always set to 1 here.

//To test compile:
//gcc -c ar2poly.c -O2 -std=c99 -Wall -Wextra
//clang -c ar2poly.c -O2 -std=c99 -Weverything
//g++ -c ar2poly.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c ar2poly.c -O2 -std=c++11 -Weverything -Wno-old-style-cast

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int ar2poly_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const float o = 1.0f;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2poly_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2poly_s: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_scopy(C,&o,0,Y,R+1);
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,&Y[c*(R+1)+1],1);
                cblas_sscal(R,-1.0f,&Y[c*(R+1)+1],1);
            }
        }
        else
        {
            cblas_scopy(C,&o,0,Y,1);
            cblas_scopy(R*C,X,1,&Y[C],1);
            cblas_sscal(R*C,-1.0f,&Y[C],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy(R,&o,0,Y,1);
            cblas_scopy(R*C,X,1,&Y[R],1);
            cblas_sscal(R*C,-1.0f,&Y[R],1);
        }
        else
        {
            cblas_scopy(R,&o,0,Y,C+1);
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,&Y[r*(C+1)+1],1);
                cblas_sscal(C,-1.0f,&Y[r*(C+1)+1],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2poly_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int ar2poly_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const double o = 1.0;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2poly_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2poly_d: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_dcopy(C,&o,0,Y,R+1);
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,&Y[c*(R+1)+1],1);
                cblas_dscal(R,-1.0,&Y[c*(R+1)+1],1);
            }
        }
        else
        {
            cblas_dcopy(C,&o,0,Y,1);
            cblas_dcopy(R*C,X,1,&Y[C],1);
            cblas_dscal(R*C,-1.0,&Y[C],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy(R,&o,0,Y,1);
            cblas_dcopy(R*C,X,1,&Y[R],1);
            cblas_dscal(R*C,-1.0,&Y[R],1);
        }
        else
        {
            cblas_dcopy(R,&o,0,Y,C+1);
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,&Y[r*(C+1)+1],1);
                cblas_dscal(C,-1.0,&Y[r*(C+1)+1],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2poly_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int ar2poly_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const float o[2] =  {1.0f,0.0f};
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2poly_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2poly_c: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_ccopy(C,&o[0],0,Y,R+1);
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[2*c*R],1,&Y[2*c*(R+1)+2],1);
                cblas_csscal(R,-1.0,&Y[2*c*(R+1)+2],1);
            }
        }
        else
        {
            cblas_ccopy(C,&o[0],0,Y,1);
            cblas_ccopy(R*C,X,1,&Y[2*C],1);
            cblas_csscal(R*C,-1.0,&Y[2*C],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_ccopy(R,&o[0],0,Y,1);
            cblas_ccopy(R*C,X,1,&Y[2*R],1);
            cblas_csscal(R*C,-1.0,&Y[2*R],1);
        }
        else
        {
            cblas_ccopy(R,&o[0],0,Y,C+1);
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[2*r*C],1,&Y[2*r*(C+1)+2],1);
                cblas_csscal(C,-1.0,&Y[2*r*(C+1)+2],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2poly_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int ar2poly_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const double o[2] =  {1.0,0.0};
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in ar2poly_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ar2poly_z: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_zcopy(C,&o[0],0,Y,R+1);
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[2*c*R],1,&Y[2*c*(R+1)+2],1);
                cblas_zdscal(R,-1.0,&Y[2*c*(R+1)+2],1);
            }
        }
        else
        {
            cblas_zcopy(C,&o[0],0,Y,1);
            cblas_zcopy(R*C,X,1,&Y[2*C],1);
            cblas_zdscal(R*C,-1.0,&Y[2*C],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_zcopy(R,&o[0],0,Y,1);
            cblas_zcopy(R*C,X,1,&Y[2*R],1);
            cblas_zdscal(R*C,-1.0,&Y[2*R],1);
        }
        else
        {
            cblas_zcopy(R,&o[0],0,Y,C+1);
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[2*r*C],1,&Y[2*r*(C+1)+2],1);
                cblas_zdscal(C,-1.0,&Y[2*r*(C+1)+2],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ar2poly_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


#ifdef __cplusplus
}
#endif

