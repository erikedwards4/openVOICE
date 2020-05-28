//Get B center-frequencies (cfs) for a given frequency scale (freq_scale).
//Before using this, cfs should be initialized as an array with B members.

//This is for filter bank purposes, for example, where a set of cfs from lofreq to hifreq
//is desired that are linearly spaced on the specified frequency scale (freq_scale).
//However, the lofreq and hifreq are specified in Hz, and the output cfs are also in Hz.

//For the frequency range, I always use Nyquist as hifreq, and 0 Hz as lofreq (except log scales use 1 Hz).
//Note that lofreq/hifreq are band *edges*, where the first/last triangle filters start/end, and not the first/last cfs.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include "convert_freqs.c"

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int get_cfs_s (float *cfs, const int B, const float lofreq, const float hifreq, const char freq_scale[], const char band_edges);
int get_cfs_d (double *cfs, const int B, const double lofreq, const double hifreq, const char freq_scale[], const char band_edges);


int get_cfs_s (float *cfs, const int B, const float lofreq, const float hifreq, const char freq_scale[], const char band_edges)
{
    int b;
    float zrng[2], zinc;

	//Checks
    if (B<2) { fprintf(stderr,"error in get_cfs_s: B must be > 1\n"); return 1; }
	if (lofreq<0.0f) { fprintf(stderr,"error in get_cfs_s: lofreq must be nonnegative\n"); return 1; }
	if (hifreq<=lofreq) { fprintf(stderr,"error in get_cfs_s: hifreq must be > lofreq\n"); return 1; }
    if (hifreq>=100000.0f) { fprintf(stderr,"error in get_cfs_s: hifreq must be < 100 kHz\n"); return 1; }
	
	//Get zrng, which is 2-member array with loz and hiz
    zrng[0] = lofreq; zrng[1] = hifreq;
	if (convert_freqs_s(zrng,2,"hz",freq_scale)) { return 1; }
    
    //Get zinc, the increment in units of freq_scale
    if (band_edges) { zinc = (zrng[1]-zrng[0])/(1.0f+B); }
    else { zinc = (zrng[1]-zrng[0])/(B-1.0f); }
	
	//Initially, set cfs as the "czs", i.e. the linear-spaced center position in freq_scale units
    if (band_edges) { cfs[0] = zrng[0] + zinc; cfs[B-1] = zrng[1] - zinc; }
    else { cfs[0] = zrng[0]; cfs[B-1] = zrng[1]; }
	for (b=1; b<B-1; b++) { cfs[b] = fmaf(b,zinc,cfs[0]); }
	
	//Now convert these from freq_scale back to Hz
	if (convert_freqs_s(cfs,B,freq_scale,"hz")) { return 1; }
    if (!band_edges) { cfs[0] = lofreq; cfs[B-1] = hifreq; }
	
    return 0;
}


int get_cfs_d (double *cfs, const int B, const double lofreq, const double hifreq, const char freq_scale[], const char band_edges)
{
    int b;
    double zrng[2], zinc;

	//Checks
    if (B<2) { fprintf(stderr,"error in get_cfs_d: B must be > 1"); return 1; }
	if (lofreq<0.0) { fprintf(stderr,"error in get_cfs_d: lofreq must be nonnegative\n"); return 1; }
	if (hifreq<=lofreq) { fprintf(stderr,"error in get_cfs_d: hifreq must be > lofreq\n"); return 1; }
    if (hifreq>=100000.0) { fprintf(stderr,"error in get_cfs_d: hifreq must be < 100 kHz\n"); return 1; }
	
	//Get zrng, which is 2-member array with loz and hiz
    zrng[0] = lofreq; zrng[1] = hifreq;
	if (convert_freqs_d(zrng,2,"hz",freq_scale)) { return 1; }
    
    //Get zinc, the increment in units of freq_scale
    if (band_edges) { zinc = (zrng[1]-zrng[0])/(1.0+B); }
    else { zinc = (zrng[1]-zrng[0])/(B-1.0); }
	
	//Initially, set cfs as the "czs", i.e. the linear-spaced center position in freq_scale units
    if (band_edges) { cfs[0] = zrng[0] + zinc; cfs[B-1] = zrng[1] - zinc; }
    else { cfs[0] = zrng[0]; cfs[B-1] = zrng[1]; }
	for (b=1; b<B-1; b++) { cfs[b] = fma(b,zinc,cfs[0]); }
	
	//Now convert these from freq_scale back to Hz
	if (convert_freqs_d(cfs,B,freq_scale,"hz")) { return 1; }
    if (!band_edges) { cfs[0] = lofreq; cfs[B-1] = hifreq; }
	
    return 0;
}


#ifdef __cplusplus
}
}
#endif

