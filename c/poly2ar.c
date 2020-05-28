//Gets autoregressive (AR) parameters from polynomials along rows or cols of X.
//If the polynomial is a0 a1 a2..., then the AR coeffs are -a1/a0 -a2/a0...
//This is invertible with ar2poly only if a0==1.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int poly2ar_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2ar_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2ar_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int poly2ar_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);


int poly2ar_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim)
{
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in poly2ar_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in poly2ar_s: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R-1,&X[c*R+1],1,&Y[c*(R-1)],1);
                cblas_sscal(R-1,-1.0f/X[c*R],&Y[c*(R-1)],1);
            }
        }
        else
        {
            cblas_scopy((R-1)*C,&X[C],1,&Y[0],1);
            for (c=0; c<C; c++)
            {
                cblas_sscal(R-1,-1.0f/X[c],&Y[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy(R*(C-1),&X[R],1,&Y[0],1);
            for (r=0; r<R; r++)
            {
                cblas_sscal(C-1,-1.0f/X[r],&Y[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C-1,&X[r*C+1],1,&Y[r*(C-1)],1);
                cblas_sscal(C-1,-1.0f/X[r*C],&Y[r*(C-1)],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2ar_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int poly2ar_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim)
{
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in poly2ar_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in poly2ar_d: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R-1,&X[c*R+1],1,&Y[c*(R-1)],1);
                cblas_dscal(R-1,-1.0/X[c*R],&Y[c*(R-1)],1);
            }
        }
        else
        {
            cblas_dcopy((R-1)*C,&X[C],1,&Y[0],1);
            for (c=0; c<C; c++)
            {
                cblas_dscal(R-1,-1.0/X[c],&Y[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy(R*(C-1),&X[R],1,&Y[0],1);
            for (r=0; r<R; r++)
            {
                cblas_dscal(C-1,-1.0/X[r],&Y[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C-1,&X[r*C+1],1,&Y[r*(C-1)],1);
                cblas_dscal(C-1,-1.0/X[r*C],&Y[r*(C-1)],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2ar_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int poly2ar_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim)
{
    int r, c;
    float sc[2] = {0.0f,0.0f};

    //Checks
    if (R<1) { fprintf(stderr,"error in poly2ar_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in poly2ar_c: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R-1,&X[2*c*R+2],1,&Y[2*c*(R-1)],1);
                sc[0] = 1.0f / (X[2*c*R]*X[2*c*R]+X[2*c*R+1]*X[2*c*R+1]);
                sc[1] = sc[0] * X[2*c*R+1]; sc[0] *= -X[2*c*R];
                cblas_cscal(R-1,&sc[0],&Y[2*c*(R-1)],1);
            }
        }
        else
        {
            cblas_ccopy((R-1)*C,&X[2*C],1,&Y[0],1);
            for (c=0; c<C; c++)
            {
                sc[0] = 1.0f / (X[2*c]*X[2*c]+X[2*c+1]*X[2*c+1]);
                sc[1] = sc[0] * X[2*c+1]; sc[0] *= -X[2*c];
                cblas_cscal(R-1,&sc[0],&Y[2*c],2*C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_ccopy(R*(C-1),&X[2*R],1,&Y[0],1);
            for (r=0; r<R; r++)
            {
                sc[0] = 1.0f / (X[2*r]*X[2*r]+X[2*r+1]*X[2*r+1]);
                sc[1] = sc[0] * X[2*r+1]; sc[0] *= -X[2*r];
                cblas_cscal(C-1,&sc[0],&Y[2*r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C-1,&X[2*r*C+2],1,&Y[2*r*(C-1)],1);
                sc[0] = 1.0f / (X[2*r*C]*X[2*r*C]+X[2*r*C+1]*X[2*r*C+1]);
                sc[1] = sc[0] * X[2*r*C+1]; sc[0] *= -X[2*r*C];
                cblas_cscal(C-1,&sc[0],&Y[2*r*(C-1)],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2ar_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int poly2ar_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim)
{
    int r, c;
    double sc[2] = {0.0,0.0};

    //Checks
    if (R<1) { fprintf(stderr,"error in poly2ar_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in poly2ar_z: ncols X must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R-1,&X[2*c*R+2],1,&Y[2*c*(R-1)],1);
                sc[0] = 1.0 / (X[2*c*R]*X[2*c*R]+X[2*c*R+1]*X[2*c*R+1]);
                sc[1] = sc[0] * X[2*c*R+1]; sc[0] *= -X[2*c*R];
                cblas_zscal(R-1,&sc[0],&Y[2*c*(R-1)],1);
            }
        }
        else
        {
            cblas_zcopy((R-1)*C,&X[2*C],1,&Y[0],1);
            for (c=0; c<C; c++)
            {
                sc[0] = 1.0 / (X[2*c]*X[2*c]+X[2*c+1]*X[2*c+1]);
                sc[1] = sc[0] * X[2*c+1]; sc[0] *= -X[2*c];
                cblas_zscal(R-1,&sc[0],&Y[2*c],2*C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_zcopy(R*(C-1),&X[2*R],1,&Y[0],1);
            for (r=0; r<R; r++)
            {
                sc[0] = 1.0 / (X[2*r]*X[2*r]+X[2*r+1]*X[2*r+1]);
                sc[1] = sc[0] * X[2*r+1]; sc[0] *= -X[2*r];
                cblas_zscal(C-1,&sc[0],&Y[2*r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C-1,&X[2*r*C+2],1,&Y[2*r*(C-1)],1);
                sc[0] = 1.0 / (X[2*r*C]*X[2*r*C]+X[2*r*C+1]*X[2*r*C+1]);
                sc[1] = sc[0] * X[2*r*C+1]; sc[0] *= -X[2*r*C];
                cblas_zscal(C-1,&sc[0],&Y[2*r*(C-1)],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2ar_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


#ifdef __cplusplus
}
}
#endif

