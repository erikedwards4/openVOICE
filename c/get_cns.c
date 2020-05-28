//Get B center-numbers (cns), which are the closest-fit integer indices into the STFT freqs.
//Before using this, cns should be initialized as an array with B members.
//This assumes that the STFT freqs are sorted ascending.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int get_cns_s (int *cns, const float *freqs, const int F, const float *cfs, const int B);
int get_cns_d (int *cns, const double *freqs, const int F, const double *cfs, const int B);


int get_cns_s (int *cns, const float *freqs, const int F, const float *cfs, const int B)
{
    int b;

    //Checks
    if (B<2) { fprintf(stderr,"error in get_cns_s: B (num cfs) must be > 1\n"); return 1; }
    if (F<2) { fprintf(stderr,"error in get_cns_s: F (num STFT freqs) must > 1\n"); return 1; }
	if (cfs[0]<0.0f) { fprintf(stderr,"error in get_cns_s: cfs must be nonnegative\n"); return 1; }
    if (freqs[0]<0.0f) { fprintf(stderr,"error in get_cns_s: STFT freqs must be nonnegative\n"); return 1; }
    if (freqs[1]<=freqs[0] || freqs[F-1]<=freqs[F-2]) { fprintf(stderr,"error in get_cns_s: STFT freqs must be sorted ascending\n"); return 1; }
	
    for (b=0; b<B; b++)
    {
        cns[b] = 0;
        while (cns[b]<F-1 && fabsf(cfs[b]-freqs[cns[b]])>fabsf(cfs[b]-freqs[cns[b]+1]))
		{
			cns[b]++;
		}
        //cfs[b] = freqs[cns[b]]; //this would reassign cf to nearest FFT freq if desired
    }
	
    return 0;
}


int get_cns_d (int *cns, const double *freqs, const int F, const double *cfs, const int B)
{
    int b;

    //Checks
    if (B<2) { fprintf(stderr,"error in get_cns_d: B (num cfs) must be > 1\n"); return 1; }
    if (F<2) { fprintf(stderr,"error in get_cns_d: F (num STFT freqs) must > 1\n"); return 1; }
	if (cfs[0]<0.0) { fprintf(stderr,"error in get_cns_d: cfs must be nonnegative\n"); return 1; }
    if (freqs[0]<0.0) { fprintf(stderr,"error in get_cns_d: STFT freqs must be nonnegative\n"); return 1; }
    if (freqs[1]<=freqs[0] || freqs[F-1]<=freqs[F-2]) { fprintf(stderr,"error in get_cns_d: STFT freqs must be sorted ascending\n"); return 1; }
	
    for (b=0; b<B; b++)
    {
        cns[b] = 0;
        while (cns[b]<F-1 && fabs(cfs[b]-freqs[cns[b]])>fabs(cfs[b]-freqs[cns[b]+1]))
		{
			cns[b]++;
		}
        //cfs[b] = freqs[cns[b]]; //this would reassign cf to nearest FFT freq if desired
    }
	
    return 0;
}


#ifdef __cplusplus
}
}
#endif

