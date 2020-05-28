//Converts LSFs (line spectral frequencies) to polynomial representation along cols or rows of X.
//Works by getting p_sym and p_asym and then roots and then lsf2poly.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int lsf2poly_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim);
int lsf2poly_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim);


int lsf2poly_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const float z = 0.0f, o = 1.0f;
    const int F = (dim==0) ? R : C;  //num LSFs
    const int nr = 2*F + 2;          //length roots
    const int P = F + 2;             //length psym and pasym
    const int Po = P - 1;            //num output poly coeffs
    int r, c, p, f;
    float *roots, *psym, *pasym;

    //Checks
    if (R<1) { fprintf(stderr,"error in lsf2poly_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lsf2poly_s: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(psym=(float *)malloc((size_t)(2*P)*sizeof(float)))) { fprintf(stderr,"error in lsf2poly_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(pasym=(float *)malloc((size_t)(2*P)*sizeof(float)))) { fprintf(stderr,"error in lsf2poly_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(roots=(float *)malloc((size_t)(2*nr)*sizeof(float)))) { fprintf(stderr,"error in lsf2poly_s: problem with malloc. "); perror("malloc"); return 1; }
    roots[0] = o; roots[nr] = -o; roots[1] = roots[nr+1] = z;

    if (dim==0)
    {
        if (iscolmajor)
        {
            //cblas_scopy(C,&o,0,Y,Po);
            for (c=0; c<C; c++)
            {
                for (f=0; f<F; f++)
                {
                    roots[2*f+2] = roots[4*F-2*f+2] = cosf(X[f+c*R]);
                    roots[2*f+3] = sinf(X[f+c*R]); roots[4*F-2*f+3] = -roots[2*f+3];
                }
                cblas_scopy(2*P,&z,0,psym,1); cblas_scopy(2*P,&z,0,pasym,1);
                psym[0] = pasym[0] = o;
                for (r=0; r<P-1; r++)
                {
                    for (p=r+1; p>0; p--)
                    {
                        pasym[2*p] = fmaf(-roots[4*r],pasym[2*p-2],fmaf(roots[4*r+1],pasym[2*p-1],pasym[2*p]));
                        pasym[2*p+1] = fmaf(-roots[4*r],pasym[2*p-1],fmaf(-roots[4*r+1],pasym[2*p-2],pasym[2*p+1]));
                        psym[2*p] = fmaf(-roots[4*r+2],psym[2*p-2],fmaf(roots[4*r+3],psym[2*p-1],psym[2*p]));
                        psym[2*p+1] = fmaf(-roots[4*r+2],psym[2*p-1],fmaf(-roots[4*r+3],psym[2*p-2],psym[2*p+1]));
                    }
                }
                for (r=0; r<Po; r++) { Y[r+c*Po] = 0.5f*(psym[2*(Po-r)]-pasym[2*(Po-r)]); }
            }
        }
        else
        {
            //cblas_scopy(C,&o,0,Y,1);
            for (c=0; c<C; c++)
            {
                for (f=0; f<F; f++)
                {
                    roots[2*f+2] = roots[4*F-2*f+2] = cosf(X[c+f*C]);
                    roots[2*f+3] = sinf(X[c+f*C]); roots[4*F-2*f+3] = -roots[2*f+3];
                }
                cblas_scopy(2*P,&z,0,psym,1); cblas_scopy(2*P,&z,0,pasym,1);
                psym[0] = pasym[0] = o;
                for (r=0; r<P-1; r++)
                {
                    for (p=r+1; p>0; p--)
                    {
                        pasym[2*p] = fmaf(-roots[4*r],pasym[2*p-2],fmaf(roots[4*r+1],pasym[2*p-1],pasym[2*p]));
                        pasym[2*p+1] = fmaf(-roots[4*r],pasym[2*p-1],fmaf(-roots[4*r+1],pasym[2*p-2],pasym[2*p+1]));
                        psym[2*p] = fmaf(-roots[4*r+2],psym[2*p-2],fmaf(roots[4*r+3],psym[2*p-1],psym[2*p]));
                        psym[2*p+1] = fmaf(-roots[4*r+2],psym[2*p-1],fmaf(-roots[4*r+3],psym[2*p-2],psym[2*p+1]));
                    }
                }
                for (r=0; r<Po; r++) { Y[c+r*C] = 0.5f*(psym[2*(Po-r)]-pasym[2*(Po-r)]); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            //cblas_scopy(R,&o,0,Y,1);
            for (r=0; r<R; r++)
            {
                for (f=0; f<F; f++)
                {
                    roots[2*f+2] = roots[4*F-2*f+2] = cosf(X[r+f*R]);
                    roots[2*f+3] = sinf(X[r+f*R]); roots[4*F-2*f+3] = -roots[2*f+3];
                }
                cblas_scopy(2*P,&z,0,psym,1); cblas_scopy(2*P,&z,0,pasym,1);
                psym[0] = pasym[0] = o;
                for (c=0; c<P-1; c++)
                {
                    for (p=c+1; p>0; p--)
                    {
                        pasym[2*p] = fmaf(-roots[4*c],pasym[2*p-2],fmaf(roots[4*c+1],pasym[2*p-1],pasym[2*p]));
                        pasym[2*p+1] = fmaf(-roots[4*c],pasym[2*p-1],fmaf(-roots[4*c+1],pasym[2*p-2],pasym[2*p+1]));
                        psym[2*p] = fmaf(-roots[4*c+2],psym[2*p-2],fmaf(roots[4*c+3],psym[2*p-1],psym[2*p]));
                        psym[2*p+1] = fmaf(-roots[4*c+2],psym[2*p-1],fmaf(-roots[4*c+3],psym[2*p-2],psym[2*p+1]));
                    }
                }
                for (c=0; c<Po; c++) { Y[r+c*R] = 0.5f*(psym[2*(Po-c)]-pasym[2*(Po-c)]); }
            }
        }
        else
        {
            //cblas_scopy(R,&o,0,Y,Po);
            for (r=0; r<R; r++)
            {
                for (f=0; f<F; f++)
                {
                    roots[2*f+2] = roots[4*F-2*f+2] = cosf(X[f+r*C]);
                    roots[2*f+3] = sinf(X[f+r*C]); roots[4*F-2*f+3] = -roots[2*f+3];
                }
                cblas_scopy(2*P,&z,0,psym,1); cblas_scopy(2*P,&z,0,pasym,1);
                psym[0] = pasym[0] = o;
                for (c=0; c<P-1; c++)
                {
                    for (p=c+1; p>0; p--)
                    {
                        pasym[2*p] = fmaf(-roots[4*c],pasym[2*p-2],fmaf(roots[4*c+1],pasym[2*p-1],pasym[2*p]));
                        pasym[2*p+1] = fmaf(-roots[4*c],pasym[2*p-1],fmaf(-roots[4*c+1],pasym[2*p-2],pasym[2*p+1]));
                        psym[2*p] = fmaf(-roots[4*c+2],psym[2*p-2],fmaf(roots[4*c+3],psym[2*p-1],psym[2*p]));
                        psym[2*p+1] = fmaf(-roots[4*c+2],psym[2*p-1],fmaf(-roots[4*c+3],psym[2*p-2],psym[2*p+1]));
                    }
                }
                for (c=0; c<Po; c++) { Y[c+r*Po] = 0.5f*(psym[2*(Po-c)]-pasym[2*(Po-c)]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in lsf2poly_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(roots); free(psym); free(pasym);
    return 0;
}


int lsf2poly_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim)
{
    const double z = 0.0, o = 1.0;
    const int F = (dim==0) ? R : C;  //num LSFs
    const int nr = 2*F + 2;          //length roots
    const int P = F + 2;             //length psym and pasym
    const int Po = P - 1;            //num output poly coeffs
    int r, c, p, f;
    double *roots, *psym, *pasym;

    //Checks
    if (R<1) { fprintf(stderr,"error in lsf2poly_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lsf2poly_d: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(psym=(double *)malloc((size_t)(2*P)*sizeof(double)))) { fprintf(stderr,"error in lsf2poly_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(pasym=(double *)malloc((size_t)(2*P)*sizeof(double)))) { fprintf(stderr,"error in lsf2poly_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(roots=(double *)malloc((size_t)(2*nr)*sizeof(double)))) { fprintf(stderr,"error in lsf2poly_d: problem with malloc. "); perror("malloc"); return 1; }
    roots[0] = o; roots[nr] = -o; roots[1] = roots[nr+1] = z;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (f=0; f<F; f++)
                {
                    roots[2*f+2] = roots[4*F-2*f+2] = cos(X[f+c*R]);
                    roots[2*f+3] = sin(X[f+c*R]); roots[4*F-2*f+3] = -roots[2*f+3];
                }
                cblas_dcopy(2*P,&z,0,psym,1); cblas_dcopy(2*P,&z,0,pasym,1);
                psym[0] = pasym[0] = o;
                for (r=0; r<P-1; r++)
                {
                    for (p=r+1; p>0; p--)
                    {
                        pasym[2*p] = fma(-roots[4*r],pasym[2*p-2],fma(roots[4*r+1],pasym[2*p-1],pasym[2*p]));
                        pasym[2*p+1] = fma(-roots[4*r],pasym[2*p-1],fma(-roots[4*r+1],pasym[2*p-2],pasym[2*p+1]));
                        psym[2*p] = fma(-roots[4*r+2],psym[2*p-2],fma(roots[4*r+3],psym[2*p-1],psym[2*p]));
                        psym[2*p+1] = fma(-roots[4*r+2],psym[2*p-1],fma(-roots[4*r+3],psym[2*p-2],psym[2*p+1]));
                    }
                }
                for (r=0; r<Po; r++) { Y[r+c*Po] = 0.5*(psym[2*(Po-r)]-pasym[2*(Po-r)]); }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (f=0; f<F; f++)
                {
                    roots[2*f+2] = roots[4*F-2*f+2] = cos(X[c+f*C]);
                    roots[2*f+3] = sin(X[c+f*C]); roots[4*F-2*f+3] = -roots[2*f+3];
                }
                cblas_dcopy(2*P,&z,0,psym,1); cblas_dcopy(2*P,&z,0,pasym,1);
                psym[0] = pasym[0] = o;
                for (r=0; r<P-1; r++)
                {
                    for (p=r+1; p>0; p--)
                    {
                        pasym[2*p] = fma(-roots[4*r],pasym[2*p-2],fma(roots[4*r+1],pasym[2*p-1],pasym[2*p]));
                        pasym[2*p+1] = fma(-roots[4*r],pasym[2*p-1],fma(-roots[4*r+1],pasym[2*p-2],pasym[2*p+1]));
                        psym[2*p] = fma(-roots[4*r+2],psym[2*p-2],fma(roots[4*r+3],psym[2*p-1],psym[2*p]));
                        psym[2*p+1] = fma(-roots[4*r+2],psym[2*p-1],fma(-roots[4*r+3],psym[2*p-2],psym[2*p+1]));
                    }
                }
                for (r=0; r<Po; r++) { Y[c+r*C] = 0.5*(psym[2*(Po-r)]-pasym[2*(Po-r)]); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (f=0; f<F; f++)
                {
                    roots[2*f+2] = roots[4*F-2*f+2] = cos(X[r+f*R]);
                    roots[2*f+3] = sin(X[r+f*R]); roots[4*F-2*f+3] = -roots[2*f+3];
                }
                cblas_dcopy(2*P,&z,0,psym,1); cblas_dcopy(2*P,&z,0,pasym,1);
                psym[0] = pasym[0] = o;
                for (c=0; c<P-1; c++)
                {
                    for (p=c+1; p>0; p--)
                    {
                        pasym[2*p] = fma(-roots[4*c],pasym[2*p-2],fma(roots[4*c+1],pasym[2*p-1],pasym[2*p]));
                        pasym[2*p+1] = fma(-roots[4*c],pasym[2*p-1],fma(-roots[4*c+1],pasym[2*p-2],pasym[2*p+1]));
                        psym[2*p] = fma(-roots[4*c+2],psym[2*p-2],fma(roots[4*c+3],psym[2*p-1],psym[2*p]));
                        psym[2*p+1] = fma(-roots[4*c+2],psym[2*p-1],fma(-roots[4*c+3],psym[2*p-2],psym[2*p+1]));
                    }
                }
                for (c=0; c<Po; c++) { Y[r+c*R] = 0.5*(psym[2*(Po-c)]-pasym[2*(Po-c)]); }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (f=0; f<F; f++)
                {
                    roots[2*f+2] = roots[4*F-2*f+2] = cos(X[f+r*C]);
                    roots[2*f+3] = sin(X[f+r*C]); roots[4*F-2*f+3] = -roots[2*f+3];
                }
                cblas_dcopy(2*P,&z,0,psym,1); cblas_dcopy(2*P,&z,0,pasym,1);
                psym[0] = pasym[0] = o;
                for (c=0; c<P-1; c++)
                {
                    for (p=c+1; p>0; p--)
                    {
                        pasym[2*p] = fma(-roots[4*c],pasym[2*p-2],fma(roots[4*c+1],pasym[2*p-1],pasym[2*p]));
                        pasym[2*p+1] = fma(-roots[4*c],pasym[2*p-1],fma(-roots[4*c+1],pasym[2*p-2],pasym[2*p+1]));
                        psym[2*p] = fma(-roots[4*c+2],psym[2*p-2],fma(roots[4*c+3],psym[2*p-1],psym[2*p]));
                        psym[2*p+1] = fma(-roots[4*c+2],psym[2*p-1],fma(-roots[4*c+3],psym[2*p-2],psym[2*p+1]));
                    }
                }
                for (c=0; c<Po; c++) { Y[c+r*Po] = 0.5*(psym[2*(Po-c)]-pasym[2*(Po-c)]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in lsf2poly_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(roots); free(psym); free(pasym);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

