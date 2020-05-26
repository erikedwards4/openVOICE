//Does Levinson-Durbin recursion from autocorrelation (AC) values.
//The output array "AS" should be initialized before calling as double as[P].
//The AC must be at least of length P for AS of length P.

//To test compile:
//gcc -c ac2ar_levdurb.c -O2 -std=c99 -Wall -Wextra
//clang -c ac2ar_levdurb.c -O2 -std=c99 -Weverything
//g++ -c ac2ar_levdurb.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c ac2ar_levdurb.c -O2 -std=c++11 -Weverything -Wno-old-style-cast

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int ac2ar_levdurb_s (float *Y, float *V, const float *AC, const char iscolmajor, const int R, const int C, const int dim, const int P)
{
    int r, c, p, q;
	float g;
    float *AS, *AStmp;

    //Checks
    if (dim==0 && P>=R) { fprintf(stderr,"error in ac2ar_levdurb_s: P must be < nrows AC for dim==0\n"); return 1; }
    if (dim==1 && P>=C) { fprintf(stderr,"error in ac2ar_levdurb_s: P must be < ncols AC for dim==1\n"); return 1; }
	
    //Initialize
    if (!(AStmp=(float *)malloc((size_t)(P)*sizeof(float)))) { fprintf(stderr,"error in ac2ar_levdurb_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AS=(float *)malloc((size_t)(P+1)*sizeof(float)))) { fprintf(stderr,"error in ac2ar_levdurb_s: problem with malloc. "); perror("malloc"); return 1; }
	AS[0] = 1.0f;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                AS[1] = AStmp[0] = g = -AC[1+c*R]/AC[c*R];
                V[c] = fmaf(AC[1+c*R],g,AC[c*R]);
                for (p=2; p<=P; p++)
                {
                    g = AC[p+c*R];
                    for (q=1; q<p; q++) { g = fmaf(AC[q+c*R],AS[p-q],g); }
                    AS[p] = g = -g/V[c];
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],1,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
                cblas_scopy(P,&AS[1],1,&Y[c*P],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                AS[1] = AStmp[0] = g = -AC[c+C]/AC[c];
                V[c] = fmaf(AC[c+C],g,AC[c]);
                for (p=2; p<=P; p++)
                {
                    g = AC[c+p*C];
                    for (q=1; q<p; q++) { g = fmaf(AC[c+q*C],AS[p-q],g); }
                    AS[p] = g = -g/V[c];
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],1,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
                cblas_scopy(P,&AS[1],1,&Y[c],C);
            }
        }
        cblas_sscal(P*C,-1.0f,Y,1);
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                AS[1] = AStmp[0] = g = -AC[r+R]/AC[r];
                V[r] = fmaf(AC[r+R],g,AC[r]);
                for (p=2; p<=P; p++)
                {
                    g = AC[r+p*R];
                    for (q=1; q<p; q++) { g = fmaf(AC[r+q*R],AS[p-q],g); }
                    AS[p] = g = -g/V[r];
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],1,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
                cblas_scopy(P,&AS[1],1,&Y[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                AS[1] = AStmp[0] = g = -AC[1+r*C]/AC[r*C];
                V[r] = fmaf(AC[1+r*C],g,AC[r*C]);
                for (p=2; p<=P; p++)
                {
                    g = AC[p+r*C];
                    for (q=1; q<p; q++) { g = fmaf(AC[q+r*C],AS[p-q],g); }
                    AS[p] = g = -g/V[r];
                    for (q=1; q<p; q++) { AS[q] = fmaf(g,AStmp[p-q-1],AS[q]); }
                    cblas_scopy(p,&AS[1],1,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
                cblas_scopy(P,&AS[1],1,&Y[r*P],1);
            }
        }
        cblas_sscal(R*P,-1.0f,Y,1);
    }
    else
    {
        fprintf(stderr,"error in ac2ar_levdurb_s: dim must be 0 or 1.\n"); return 1;
    }

    //Finish
	return 0;
}


int ac2ar_levdurb_d (double *Y, double *V, const double *AC, const char iscolmajor, const int R, const int C, const int dim, const int P)
{
    int r, c, p, q;
	double g;
    double *AS, *AStmp;

    //Checks
    if (dim==0 && P>=R) { fprintf(stderr,"error in ac2ar_levdurb_d: P must be < nrows AC for dim==0\n"); return 1; }
    if (dim==1 && P>=C) { fprintf(stderr,"error in ac2ar_levdurb_d: P must be < ncols AC for dim==1\n"); return 1; }
	
    //Initialize
    if (!(AStmp=(double *)malloc((size_t)(P)*sizeof(double)))) { fprintf(stderr,"error in ac2ar_levdurb_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(AS=(double *)malloc((size_t)(P+1)*sizeof(double)))) { fprintf(stderr,"error in ac2ar_levdurb_d: problem with malloc. "); perror("malloc"); return 1; }
	AS[0] = 1.0;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                AS[1] = AStmp[0] = g = -AC[1+c*R]/AC[c*R];
                V[c] = fma(AC[1+c*R],g,AC[c*R]);
                for (p=2; p<=P; p++)
                {
                    g = AC[p+c*R];
                    for (q=1; q<p; q++) { g = fma(AC[q+c*R],AS[p-q],g); }
                    AS[p] = g = -g/V[c];
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],1,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
                cblas_dcopy(P,&AS[1],1,&Y[c*P],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                AS[1] = AStmp[0] = g = -AC[c+C]/AC[c];
                V[c] = fma(AC[c+C],g,AC[c]);
                for (p=2; p<=P; p++)
                {
                    g = AC[c+p*C];
                    for (q=1; q<p; q++) { g = fma(AC[c+q*C],AS[p-q],g); }
                    AS[p] = g = -g/V[c];
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],1,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
                cblas_dcopy(P,&AS[1],1,&Y[c],C);
            }
        }
        cblas_dscal(P*C,-1.0,Y,1);
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                AS[1] = AStmp[0] = g = -AC[r+R]/AC[r];
                V[r] = fma(AC[r+R],g,AC[r]);
                for (p=2; p<=P; p++)
                {
                    g = AC[r+p*R];
                    for (q=1; q<p; q++) { g = fma(AC[r+q*R],AS[p-q],g); }
                    AS[p] = g = -g/V[r];
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],1,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
                cblas_dcopy(P,&AS[1],1,&Y[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                AS[1] = AStmp[0] = g = -AC[1+r*C]/AC[r*C];
                V[r] = fma(AC[1+r*C],g,AC[r*C]);
                for (p=2; p<=P; p++)
                {
                    g = AC[p+r*C];
                    for (q=1; q<p; q++) { g = fma(AC[q+r*C],AS[p-q],g); }
                    AS[p] = g = -g/V[r];
                    for (q=1; q<p; q++) { AS[q] = fma(g,AStmp[p-q-1],AS[q]); }
                    cblas_dcopy(p,&AS[1],1,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
                cblas_dcopy(P,&AS[1],1,&Y[r*P],1);
            }
        }
        cblas_dscal(R*P,-1.0,Y,1);
    }
    else
    {
        fprintf(stderr,"error in ac2ar_levdurb_d: dim must be 0 or 1.\n"); return 1;
    }

    //Finish
	return 0;
}


#ifdef __cplusplus
}
#endif

