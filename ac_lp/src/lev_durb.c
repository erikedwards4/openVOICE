//Does Levinson-Durbin recursion from autocorrelation (AC) values.
//This was based originally on Matlab/Octave levinson.m, (except no R\ac option).
//The output array "as" should be initialized before calling as double as[P+1].
//The AC must be at least of length P+1 for as of length P+1.
//Reflection coefficients (rcs) can also be obtained from this.

//To test compile:
//gcc -c lev_durb.c -O2 -std=c99 -Wall -Wextra
//clang -c lev_durb.c -O2 -std=c99 -Weverything
//g++ -c lev_durb.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c lev_durb.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int lev_durb_s (float *AS, float *V, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int P)
{
    const float o = 1.0f;
    int r, c, p, q;
	float g;
    float *AStmp; //, *rcs

    //Checks
    if (dim==0 && P>R) { fprintf(stderr,"error in lev_durb_s: P must be < nrows AC for dim==0\n"); return 1; }
    if (dim==1 && P>C) { fprintf(stderr,"error in lev_durb_s: P must be < ncols AC for dim==1\n"); return 1; }
	
    //Initialize
    if (!(AStmp=(float *)malloc((size_t)(P-1)*sizeof(float)))) { fprintf(stderr,"error in lev_durb_s: problem with malloc. "); perror("malloc"); return 1; }
    //if (!(rcs=(float *)malloc((size_t)(P-1)*sizeof(float)))) { fprintf(stderr,"error in lev_durb_s: problem with malloc. "); perror("malloc"); return 1; }
	
    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_scopy(C,&o,0,&AS[0],P);
            for (c=0; c<C; c++)
            {
                AS[c*P+1] = AStmp[0] = g = -AC[1+c*R]/AC[c*R];
                V[c] = fmaf(AC[1+c*R],g,AC[c*R]);
                for (p=2; p<P; p++)
                {
                    g = AC[p+c*R];
                    for (q=1; q<p; q++) { g = fmaf(AC[q+c*R],AS[p-q+c*P],g); }
                    AS[p+c*P] = g = -g/V[c]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { AS[q+c*P] = fmaf(g,AStmp[p-q-1],AS[q+c*P]); }
                    cblas_scopy(p,&AS[1+c*P],1,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
            }
        }
        else
        {
            cblas_scopy(C,&o,0,&AS[0],1);
            for (c=0; c<C; c++)
            {
                AS[c+C] = AStmp[0] = g = -AC[c+C]/AC[c];
                V[c] = fmaf(AC[c+C],g,AC[c]);
                for (p=2; p<P; p++)
                {
                    g = AC[c+p*C];
                    for (q=1; q<p; q++) { g = fmaf(AC[c+q*C],AS[c+(p-q)*C],g); }
                    AS[c+p*C] = g = -g/V[c]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { AS[c+q*C] = fmaf(g,AStmp[p-q-1],AS[c+q*C]); }
                    cblas_scopy(p,&AS[c+C],C,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy(R,&o,0,&AS[0],1);
            for (r=0; r<R; r++)
            {
                AS[r+R] = AStmp[0] = g = -AC[r+R]/AC[r];
                V[r] = fmaf(AC[r+R],g,AC[r]);
                for (p=2; p<P; p++)
                {
                    g = AC[r+p*R];
                    for (q=1; q<p; q++) { g = fmaf(AC[r+q*R],AS[r+(p-q)*R],g); }
                    AS[r+p*R] = g = -g/V[r]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { AS[r+q*R] = fmaf(g,AStmp[p-q-1],AS[r+q*R]); }
                    cblas_scopy(p,&AS[r+R],R,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
            }
        }
        else
        {
            cblas_scopy(R,&o,0,&AS[0],P);
            for (r=0; r<R; r++)
            {
                AS[r*P+1] = AStmp[0] = g = -AC[1+r*C]/AC[r*C];
                V[r] = fmaf(AC[1+r*C],g,AC[r*C]);
                for (p=2; p<P; p++)
                {
                    g = AC[p+r*C];
                    for (q=1; q<p; q++) { g = fmaf(AC[q+r*C],AS[p-q+r*P],g); }
                    AS[p+r*P] = g = -g/V[r]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { AS[q+r*P] = fmaf(g,AStmp[p-q-1],AS[q+r*P]); }
                    cblas_scopy(p,&AS[1+r*P],1,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in lev_durb_s: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


int lev_durb_d (double *AS, double *V, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int P)
{
	const double o = 1.0;
    int r, c, p, q;
	double g;
    double *AStmp; //, *rcs

    //Checks
    if (dim==0 && P>R) { fprintf(stderr,"error in lev_durb_d: P must be < nrows AC for dim==0\n"); return 1; }
    if (dim==1 && P>C) { fprintf(stderr,"error in lev_durb_d: P must be < ncols AC for dim==1\n"); return 1; }
	
    //Initialize
    if (!(AStmp=(double *)malloc((size_t)(P-1)*sizeof(double)))) { fprintf(stderr,"error in lev_durb_d: problem with malloc. "); perror("malloc"); return 1; }
    //if (!(rcs=(double *)malloc((size_t)(P-1)*sizeof(double)))) { fprintf(stderr,"error in lev_durb_d: problem with malloc. "); perror("malloc"); return 1; }
	
    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_dcopy(C,&o,0,&AS[0],P);
            for (c=0; c<C; c++)
            {
                AS[c*P+1] = AStmp[0] = g = -AC[1+c*R]/AC[c*R];
                V[c] = fma(AC[1+c*R],g,AC[c*R]);
                for (p=2; p<P; p++)
                {
                    g = AC[p+c*R];
                    for (q=1; q<p; q++) { g = fma(AC[q+c*R],AS[p-q+c*P],g); }
                    AS[p+c*P] = g = -g/V[c]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { AS[q+c*P] = fma(g,AStmp[p-q-1],AS[q+c*P]); }
                    cblas_dcopy(p,&AS[1+c*P],1,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
            }
        }
        else
        {
            cblas_dcopy(C,&o,0,&AS[0],1);
            for (c=0; c<C; c++)
            {
                AS[c+C] = AStmp[0] = g = -AC[c+C]/AC[c];
                V[c] = fma(AC[c+C],g,AC[c]);
                for (p=2; p<P; p++)
                {
                    g = AC[c+p*C];
                    for (q=1; q<p; q++) { g = fma(AC[c+q*C],AS[c+(p-q)*C],g); }
                    AS[c+p*C] = g = -g/V[c]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { AS[c+q*C] = fma(g,AStmp[p-q-1],AS[c+q*C]); }
                    cblas_dcopy(p,&AS[c+C],C,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy(R,&o,0,&AS[0],1);
            for (r=0; r<R; r++)
            {
                AS[r+R] = AStmp[0] = g = -AC[r+R]/AC[r];
                V[r] = fma(AC[r+R],g,AC[r]);
                for (p=2; p<P; p++)
                {
                    g = AC[r+p*R];
                    for (q=1; q<p; q++) { g = fma(AC[r+q*R],AS[r+(p-q)*R],g); }
                    AS[r+p*R] = g = -g/V[r]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { AS[r+q*R] = fma(g,AStmp[p-q-1],AS[r+q*R]); }
                    cblas_dcopy(p,&AS[r+R],R,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
            }
        }
        else
        {
            cblas_dcopy(R,&o,0,&AS[0],P);
            for (r=0; r<R; r++)
            {
                AS[r*P+1] = AStmp[0] = g = -AC[1+r*C]/AC[r*C];
                V[r] = fma(AC[1+r*C],g,AC[r*C]);
                for (p=2; p<P; p++)
                {
                    g = AC[p+r*C];
                    for (q=1; q<p; q++) { g = fma(AC[q+r*C],AS[p-q+r*P],g); }
                    AS[p+r*P] = g = -g/V[r]; //rcs[p-1] = g;
                    for (q=1; q<p; q++) { AS[q+r*P] = fma(g,AStmp[p-q-1],AS[q+r*P]); }
                    cblas_dcopy(p,&AS[1+r*P],1,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in lev_durb_d: dim must be 0 or 1.\n"); return 1;
    }
	
	return 0;
}


#ifdef __cplusplus
}
#endif

