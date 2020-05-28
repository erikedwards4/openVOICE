//Get B center-frequencies (cfs) for a given frequency scale (freq_scale).
//Before using this, cfs should be initialized as an array with B members.

//This is for purposes of the spectrogram transform matrix T, i.e. for a triangular filter re-weighting.
//It is therefore simpler to use with fewer options than get_cfs (which is more general purpose).

//For the frequency range, I always use Nyquist as hifreq, and 0 Hz as lofreq (except log scales use 1 Hz).
//Note that lofreq/hifreq are band *edges*, where the first/last triangle filters start/end, and not the first/last cfs.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include "convert_freqs.c"

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int get_cfs_T_s (float *cfs, const int B, const float fs, const char freq_scale[]);
int get_cfs_T_d (double *cfs, const int B, const double fs, const char freq_scale[]);


int get_cfs_T_s (float *cfs, const int B, const float fs, const char freq_scale[])
{
    int b;
    float zrng[2], zinc;

    //Checks
    if (B<1) { fprintf(stderr,"error in get_cfs_T_s: B must be positive\n"); return 1; }
    if (fs<=0.0f) { fprintf(stderr,"error in get_cfs_T_s: fs (sampling frequency) must be positive\n"); return 1; }

	//Put lofreq at 1.0 for logarithmic scales
	if (strncmp(freq_scale,"log",3)==0 || strncmp(freq_scale,"midi",4)==0 || strncmp(freq_scale,"octave",6)==0 || strncmp(freq_scale,"piano",5)==0)
	{
		zrng[0] = 1.0f; zrng[1] = fs/2.0f;
	}
    else
    {
        zrng[0] = 0.0f; zrng[1] = fs/2.0f;
    }
	
	//Get zrng, which is 2-member array with loz and hiz
	if (convert_freqs_s(zrng,2,"hz",freq_scale)) { return 1; }
    zinc = (zrng[1]-zrng[0])/(1.0f+B);
	
	//Initially, set cfs as the "czs", i.e. the linear-spaced center position in freq_scale units
    //cfs[0] = zrng[0] + zinc; for (b=1; b<B; b++) { cfs[b] = cfs[b-1] + zinc; }  //fastest
	for (b=0; b<B; b++) { cfs[b] = fmaf(1.0f+b,zinc,zrng[0]); }                   //most accurate
	
	//Now convert these from freq_scale back to Hz
	if (convert_freqs_s(cfs,B,freq_scale,"hz")) { return 1; }
	
    return 0;
}


int get_cfs_T_d (double *cfs, const int B, const double fs, const char freq_scale[])
{
    int b;
    double zrng[2], zinc;

    //Checks
    if (B<1) { fprintf(stderr,"error in get_cfs_T_d: B must be positive\n"); return 1; }
    if (fs<=0.0) { fprintf(stderr,"error in get_cfs_T_d: fs (sampling frequency) must be positive\n"); return 1; }

	//Put lofreq at 1.0 for logarithmic scales
	if (strncmp(freq_scale,"log",3)==0 || strncmp(freq_scale,"midi",4)==0 || strncmp(freq_scale,"octave",6)==0 || strncmp(freq_scale,"piano",5)==0)
	{
		zrng[0] = 1.0; zrng[1] = fs/2.0;
	}
    else
    {
        zrng[0] = 0.0; zrng[1] = fs/2.0;
    }
	
	//Get zrng, which is 2-member array with loz and hiz
	if (convert_freqs_d(zrng,2,"hz",freq_scale)) { return 1; }
    zinc = (zrng[1]-zrng[0])/(1.0+B);
	
	//Initially, set cfs as the "czs", i.e. the linear-spaced center position in freq_scale units
    //cfs[0] = zrng[0] + zinc; for (b=1; b<B; b++) { cfs[b] = cfs[b-1] + zinc; }  //fastest
	for (b=0; b<B; b++) { cfs[b] = fma(1.0+b,zinc,zrng[0]); }                     //most accurate
	
	//Now convert these from freq_scale back to Hz
	if (convert_freqs_d(cfs,B,freq_scale,"hz")) { return 1; }
	
    return 0;
}


#ifdef __cplusplus
}
}
#endif

