//Does Levinson-Durbin recursion from autocorrelation (AC) values.
//This was based originally on Matlab/Octave levinson.m, (except no R\ac option).
//The output array "AS" should be initialized before calling AS double as[P+1].
//The AC must be at least of length P+1 for AS of length P+1.
//Reflection coefficients (RCs) can also be obtained from this.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int ac2poly_levdurb_s (float *Y, float *V, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int P);
int ac2poly_levdurb_d (double *Y, double *V, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int P);


int ac2poly_levdurb_s (float *Y, float *V, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int P)
{
    const float o = 1.0f;
    const int P1 = P + 1; //size of Y along dim
    int r, c, p, q;
	float g;
    float *AStmp; //, *rcs

    //Checks
    if (dim==0 && P>=R) { fprintf(stderr,"error in ac2poly_levdurb_s: P must be < nrows AC for dim==0\n"); return 1; }
    if (dim==1 && P>=C) { fprintf(stderr,"error in ac2poly_levdurb_s: P must be < ncols AC for dim==1\n"); return 1; }
	
    //Initialize
    if (!(AStmp=(float *)malloc((size_t)(P)*sizeof(float)))) { fprintf(stderr,"error in ac2poly_levdurb_s: problem with malloc. "); perror("malloc"); return 1; }
	
    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_scopy(C,&o,0,Y,P1);
            for (c=0; c<C; c++)
            {
                Y[c*P1+1] = AStmp[0] = g = -AC[1+c*R]/AC[c*R];
                V[c] = fmaf(AC[1+c*R],g,AC[c*R]);
                for (p=2; p<=P; p++)
                {
                    g = AC[p+c*R];
                    for (q=1; q<p; q++) { g = fmaf(AC[q+c*R],Y[p-q+c*P1],g); }
                    Y[p+c*P1] = g = -g/V[c]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { Y[q+c*P1] = fmaf(g,AStmp[p-q-1],Y[q+c*P1]); }
                    cblas_scopy(p,&Y[1+c*P1],1,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
            }
        }
        else
        {
            cblas_scopy(C,&o,0,Y,1);
            for (c=0; c<C; c++)
            {
                Y[c+C] = AStmp[0] = g = -AC[c+C]/AC[c];
                V[c] = fmaf(AC[c+C],g,AC[c]);
                for (p=2; p<=P; p++)
                {
                    g = AC[c+p*C];
                    for (q=1; q<p; q++) { g = fmaf(AC[c+q*C],Y[c+(p-q)*C],g); }
                    Y[c+p*C] = g = -g/V[c]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { Y[c+q*C] = fmaf(g,AStmp[p-q-1],Y[c+q*C]); }
                    cblas_scopy(p,&Y[c+C],C,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy(R,&o,0,Y,1);
            for (r=0; r<R; r++)
            {
                Y[r+R] = AStmp[0] = g = -AC[r+R]/AC[r];
                V[r] = fmaf(AC[r+R],g,AC[r]);
                for (p=2; p<=P; p++)
                {
                    g = AC[r+p*R];
                    for (q=1; q<p; q++) { g = fmaf(AC[r+q*R],Y[r+(p-q)*R],g); }
                    Y[r+p*R] = g = -g/V[r]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { Y[r+q*R] = fmaf(g,AStmp[p-q-1],Y[r+q*R]); }
                    cblas_scopy(p,&Y[r+R],R,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
            }
        }
        else
        {
            cblas_scopy(R,&o,0,Y,P1);
            for (r=0; r<R; r++)
            {
                Y[r*P1+1] = AStmp[0] = g = -AC[1+r*C]/AC[r*C];
                V[r] = fmaf(AC[1+r*C],g,AC[r*C]);
                for (p=2; p<=P; p++)
                {
                    g = AC[p+r*C];
                    for (q=1; q<p; q++) { g = fmaf(AC[q+r*C],Y[p-q+r*P1],g); }
                    Y[p+r*P1] = g = -g/V[r]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { Y[q+r*P1] = fmaf(g,AStmp[p-q-1],Y[q+r*P1]); }
                    cblas_scopy(p,&Y[1+r*P1],1,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ac2poly_levdurb_s: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


int ac2poly_levdurb_d (double *Y, double *V, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int P)
{
	const double o = 1.0;
    const int P1 = P + 1; //size of Y along dim
    int r, c, p, q;
	double g;
    double *AStmp; //, *rcs

    //Checks
    if (dim==0 && P>=R) { fprintf(stderr,"error in ac2poly_levdurb_d: P must be < nrows AC for dim==0\n"); return 1; }
    if (dim==1 && P>=C) { fprintf(stderr,"error in ac2poly_levdurb_d: P must be < ncols AC for dim==1\n"); return 1; }
	
    //Initialize
    if (!(AStmp=(double *)malloc((size_t)(P)*sizeof(double)))) { fprintf(stderr,"error in ac2poly_levdurb_d: problem with malloc. "); perror("malloc"); return 1; }
    //if (!(rcs=(double *)malloc((size_t)(P)*sizeof(double)))) { fprintf(stderr,"error in ac2poly_levdurb_d: problem with malloc. "); perror("malloc"); return 1; }
	
    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_dcopy(C,&o,0,Y,P1);
            for (c=0; c<C; c++)
            {
                Y[c*P1+1] = AStmp[0] = g = -AC[1+c*R]/AC[c*R];
                V[c] = fma(AC[1+c*R],g,AC[c*R]);
                for (p=2; p<=P; p++)
                {
                    g = AC[p+c*R];
                    for (q=1; q<p; q++) { g = fma(AC[q+c*R],Y[p-q+c*P1],g); }
                    Y[p+c*P1] = g = -g/V[c]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { Y[q+c*P1] = fma(g,AStmp[p-q-1],Y[q+c*P1]); }
                    cblas_dcopy(p,&Y[1+c*P1],1,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
            }
        }
        else
        {
            cblas_dcopy(C,&o,0,Y,1);
            for (c=0; c<C; c++)
            {
                Y[c+C] = AStmp[0] = g = -AC[c+C]/AC[c];
                V[c] = fma(AC[c+C],g,AC[c]);
                for (p=2; p<=P; p++)
                {
                    g = AC[c+p*C];
                    for (q=1; q<p; q++) { g = fma(AC[c+q*C],Y[c+(p-q)*C],g); }
                    Y[c+p*C] = g = -g/V[c]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { Y[c+q*C] = fma(g,AStmp[p-q-1],Y[c+q*C]); }
                    cblas_dcopy(p,&Y[c+C],C,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy(R,&o,0,Y,1);
            for (r=0; r<R; r++)
            {
                Y[r+R] = AStmp[0] = g = -AC[r+R]/AC[r];
                V[r] = fma(AC[r+R],g,AC[r]);
                for (p=2; p<=P; p++)
                {
                    g = AC[r+p*R];
                    for (q=1; q<p; q++) { g = fma(AC[r+q*R],Y[r+(p-q)*R],g); }
                    Y[r+p*R] = g = -g/V[r]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { Y[r+q*R] = fma(g,AStmp[p-q-1],Y[r+q*R]); }
                    cblas_dcopy(p,&Y[r+R],R,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
            }
        }
        else
        {
            cblas_dcopy(R,&o,0,Y,P1);
            for (r=0; r<R; r++)
            {
                Y[r*P1+1] = AStmp[0] = g = -AC[1+r*C]/AC[r*C];
                V[r] = fma(AC[1+r*C],g,AC[r*C]);
                for (p=2; p<=P; p++)
                {
                    g = AC[p+r*C];
                    for (q=1; q<p; q++) { g = fma(AC[q+r*C],Y[p-q+r*P1],g); }
                    Y[p+r*P1] = g = -g/V[r]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { Y[q+r*P1] = fma(g,AStmp[p-q-1],Y[q+r*P1]); }
                    cblas_dcopy(p,&Y[1+r*P1],1,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in ac2poly_levdurb_d: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


#ifdef __cplusplus
}
}
#endif

