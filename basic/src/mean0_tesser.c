//To test compile:
//gcc -c mean0.c -O2 -std=c99 -Wall -Wextra
//clang -c mean0.c -O2 -std=c99 -Weverything
//g++ -c mean0.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c mean0.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int mean0_s (Tesser_s X, const int dim)
{
    const float o = 1.0f;
    float m;
    int r, c;

    //X.Checks
    if (X.R<1) { fprintf(stderr,"error in mean0_s: X.R (nrows X) must be positive\n"); return 1; }
    if (X.C<1) { fprintf(stderr,"error in mean0_s: X.C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (X.iscolmajor)
        {
            for (c=0; c<X.C; c++)
            {
                m = cblas_sdot(X.R,&X.dat[c*X.R],1,&o,0) / X.R;
                cblas_saxpy(X.R,-m,&o,0,&X.dat[c*X.R],1);
            }
        }
        else
        {
            for (c=0; c<X.C; c++)
            {
                m = cblas_sdot(X.R,&X.dat[c],X.C,&o,0) / X.R;
                cblas_saxpy(X.R,-m,&o,0,&X.dat[c],X.C);
            }
        }
    }
    else if (dim==1)
    {
        if (X.iscolmajor)
        {
            for (r=0; r<X.R; r++)
            {
                m = cblas_sdot(X.C,&X.dat[r],X.R,&o,0) / X.C;
                cblas_saxpy(X.C,-m,&o,0,&X.dat[r],X.R);
            }
        }
        else
        {
            for (r=0; r<X.R; r++)
            {
                m = cblas_sdot(X.C,&X.dat[r*X.C],1,&o,0) / X.C;
                cblas_saxpy(X.C,-m,&o,0,&X.dat[r*X.C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in mean0_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int mean0_d (Tesser_d X, const int dim)
{
    const double o = 1.0;
    double m;
    int r, c;

    //X.Checks
    if (X.R<1) { fprintf(stderr,"error in mean0_d: X.R (nrows X) must be positive\n"); return 1; }
    if (X.C<1) { fprintf(stderr,"error in mean0_d: X.C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (X.iscolmajor)
        {
            for (c=0; c<X.C; c++)
            {
                m = cblas_ddot(X.R,&X.dat[c*X.R],1,&o,0) / X.R;
                cblas_daxpy(X.R,-m,&o,0,&X.dat[c*X.R],1);
            }
        }
        else
        {
            for (c=0; c<X.C; c++)
            {
                m = cblas_ddot(X.R,&X.dat[c],X.C,&o,0) / X.R;
                cblas_daxpy(X.R,-m,&o,0,&X.dat[c],X.C);
            }
        }
    }
    else if (dim==1)
    {
        if (X.iscolmajor)
        {
            for (r=0; r<X.R; r++)
            {
                m = cblas_ddot(X.C,&X.dat[r],X.R,&o,0) / X.C;
                cblas_daxpy(X.C,-m,&o,0,&X.dat[r],X.R);
            }
        }
        else
        {
            for (r=0; r<X.R; r++)
            {
                m = cblas_ddot(X.C,&X.dat[r*X.C],1,&o,0) / X.C;
                cblas_daxpy(X.C,-m,&o,0,&X.dat[r*X.C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in mean0_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int mean0_c (Tesser_s X, const int dim)
{
    //const float o = 1.0f; float m;
    const float o[2] =  {1.0f,0.0f};
    _Complex float m;
    int r, c;

    //X.Checks
    if (X.R<1) { fprintf(stderr,"error in mean0_c: X.R (nrows X) must be positive\n"); return 1; }
    if (X.C<1) { fprintf(stderr,"error in mean0_c: X.C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (X.iscolmajor)
        {
            for (c=0; c<2*X.C; c+=2)
            {
                //m = cblas_sdot(X.R,&X.dat[c*X.R],2,&o,0) / X.R;
                //cblas_saxpy(X.R,-m,&o,0,&X.dat[c*X.R],2);
                //m = cblas_sdot(X.R,&X.dat[c*X.R+1],2,&o,0) / X.R;
                //cblas_saxpy(X.R,-m,&o,0,&X.dat[c*X.R+1],2);

                //cblas_saxpy(X.R,-cblas_sdotu(X.R,&X.dat[c*X.R],2,&o,0)/X.R,&o,0,&X.dat[c*X.R],2);
                //cblas_saxpy(X.R,-cblas_sdotu(X.R,&X.dat[c*X.R+1],2,&o,0)/X.R,&o,0,&X.dat[c*X.R+1],2);

                m = -cblas_cdotu(X.R,&X.dat[c*X.R],1,&o[0],0) / X.R;
                cblas_caxpy(X.R,(float *)&m,&o[0],0,&X.dat[c*X.R],1);
            }
        }
        else
        {
            for (c=0; c<2*X.C; c+=2)
            {
                //m = cblas_sdot(X.R,&X.dat[2*c],2*X.C,&o,0) / X.R;
                //cblas_saxpy(X.R,-m,&o,0,&X.dat[2*c],2*X.C);
                //m = cblas_sdot(X.R,&X.dat[2*c+1],2*X.C,&o,0) / X.R;
                //cblas_saxpy(X.R,-m,&o,0,&X.dat[2*c+1],2*X.C);
                
                m = -cblas_cdotu(X.R,&X.dat[c],X.C,&o[0],0) / X.R;
                cblas_caxpy(X.R,(float *)&m,&o[0],0,&X.dat[c],X.C);
            }
        }
    }
    else if (dim==1)
    {
        if (X.iscolmajor)
        {
            for (r=0; r<2*X.R; r+=2)
            {
                m = -cblas_cdotu(X.C,&X.dat[r],X.R,&o[0],0) / X.C;
                cblas_caxpy(X.C,(float *)&m,&o[0],0,&X.dat[r],X.R);
            }
        }
        else
        {
            for (r=0; r<2*X.R; r+=2)
            {
                m = -cblas_cdotu(X.C,&X.dat[r*X.C],1,&o[0],0) / X.C;
                cblas_caxpy(X.C,(float *)&m,&o[0],0,&X.dat[r*X.C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in mean0_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int mean0_z (Tesser_d X, const int dim)
{
    const double o[2] =  {1.0,0.0};
    _Complex double m;
    int r, c;

    //X.Checks
    if (X.R<1) { fprintf(stderr,"error in mean0_z: X.R (nrows X) must be positive\n"); return 1; }
    if (X.C<1) { fprintf(stderr,"error in mean0_z: X.C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (X.iscolmajor)
        {
            for (c=0; c<2*X.C; c+=2)
            {
                m = -cblas_zdotu(X.R,&X.dat[c*X.R],1,&o[0],0) / X.R;
                cblas_zaxpy(X.R,(double *)&m,&o[0],0,&X.dat[c*X.R],1);
            }
        }
        else
        {
            for (c=0; c<2*X.C; c+=2)
            {
                m = -cblas_zdotu(X.R,&X.dat[c],X.C,&o[0],0) / X.R;
                cblas_zaxpy(X.R,(double *)&m,&o[0],0,&X.dat[c],X.C);
            }
        }
    }
    else if (dim==1)
    {
        if (X.iscolmajor)
        {
            for (r=0; r<2*X.R; r+=2)
            {
                m = -cblas_zdotu(X.C,&X.dat[r],X.R,&o[0],0) / X.C;
                cblas_zaxpy(X.C,(double *)&m,&o[0],0,&X.dat[r],X.R);
            }
        }
        else
        {
            for (r=0; r<2*X.R; r+=2)
            {
                m = -cblas_zdotu(X.C,&X.dat[r*X.C],1,&o[0],0) / X.C;
                cblas_zaxpy(X.C,(double *)&m,&o[0],0,&X.dat[r*X.C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in mean0_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
#endif

