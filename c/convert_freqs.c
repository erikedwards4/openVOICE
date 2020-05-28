//This converts between several known requency scales:
//bark, cochlea, erb, hz, mel, midi, octave, piano.
//
//I also include a couple of bio-frequency scales of my own,
//based on anatomical and physiological studies of the auditory periphery:
//ihcs (inner hair cells), sgcs (spiral ganglion cells).
//
//The later require a quick linear interpolation,
//so I made my own quick functions, interp1q_s and interp1q_d.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include "ihcs.h"
#include "sgcs.h"
#include "interp1q.c"

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int convert_freqs_s (float *frqs, const int F, const char in_scale[], const char out_scale[]);
int convert_freqs_d (double *frqs, const int F, const char in_scale[], const char out_scale[]);


int convert_freqs_s (float *frqs, const int F, const char in_scale[], const char out_scale[])
{
	//Declaration
	int f, cf;

	//Checks
    if (F<1) { fprintf(stderr,"error in convert_freqs_s: F (num freqs to convert) must be positive\n"); return 1; }
	if (strlen(in_scale)<2) { fprintf(stderr,"error in convert_freqs_s: input freq scale must be string with length > 1\n"); return 1; }
	if (strlen(out_scale)<2) { fprintf(stderr,"error in convert_freqs_s: output freq scale must be string with length > 1\n"); return 1; }
	
	//Convert values in frqs
	if (strncmp(in_scale,"bark",4)==0) //Bark
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10f((600.0f*sinhf(frqs[f]/6.0f))/165.4f+0.88f) / 2.1f; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680f * logf(1.0f+27639.2274f*sinhf(frqs[f]/6.0f))/(600.0f*sinhf(frqs[f]/6.0f)+14678.494617f); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 600.0f * sinhf(frqs[f]/6.0f); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0f * log10f(1.0f + (6.0f/7.0f)*sinhf(frqs[f]/6.0f)); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f((600.0f/440.0f)*sinhf(frqs[f]/6.0f)) + 69.0f; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2f(600.0f*sinhf(frqs[f]/6.0f)); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f((600.0f/440.0f)*sinhf(frqs[f]/6.0f)) + 49.0f; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_s(frqs,F,"bark","hz"); convert_freqs_s(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_s(frqs,F,"bark","hz"); convert_freqs_s(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"cochlea",7)==0) //Cochlea
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0f * asinhf((165.4f*(powf(10.0f,2.1f*frqs[f])-0.88f))/600.0f); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680f * logf(1.0f + 7619.213687f/(165.4f+14678.494617f/(powf(10.0f,2.1f*frqs[f])-0.88f))); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 165.4f * (powf(10.0f,2.1f*frqs[f]) - 0.88f); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0f * log10f(1.0f + (165.4f/700.0f)*(powf(10.0f,2.1f*frqs[f])-0.88f)); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f((165.4f/440.0f)*(powf(10.0f,2.1f*frqs[f])-0.88f)) + 69.0f; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2f(165.4f*(powf(10.0f,2.1f*frqs[f])-0.88f)); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f((165.4f/440.0f)*(powf(10.0f,2.1f*frqs[f])-0.88f)) + 49.0f; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_s(frqs,F,"cochlea","hz"); convert_freqs_s(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_s(frqs,F,"cochlea","hz"); convert_freqs_s(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"erb",3)==0) //ERB
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0f * asinhf((676170.419311f/(47.065379f-expf(0.089504f*frqs[f]))-14678.494617f)/600.0f); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10f((676170.419311f/(47.065379f-expf(0.089504f*frqs[f]))-14678.494617f)/165.4f+0.88f) / 2.1f; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 676170.419311f/(47.065379f-expf(0.089504f*frqs[f])) - 14678.494617f; }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0f * log10f(1.0f+(676170.419311f/(47.065379f-expf(0.089504f*frqs[f]))-14678.494617f)/700.0f); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f((676170.419311f/(47.065379f-expf(0.089504f*frqs[f]))-14678.494617f)/440.0f) + 69.0f; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2f((676170.419311f/(47.065379f-expf(0.089504f*frqs[f]))-14678.494617f)); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f((676170.419311f/(47.065379f-expf(0.089504f*frqs[f]))-14678.494617f)/440.0f) + 49.0f; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_s(frqs,F,"erb","hz"); convert_freqs_s(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_s(frqs,F,"erb","hz"); convert_freqs_s(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"hz",2)==0) //Hz
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0f * asinhf(frqs[f]/600.0f); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10f(frqs[f]/165.4f+0.88f) / 2.1f; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680f * logf(1.0f+46.065379f*frqs[f]/(frqs[f]+14678.494617f)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0f * log10f(1.0f+frqs[f]/700.0f); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f(frqs[f]/440.0f) + 69.0f; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2f(frqs[f]); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f(frqs[f]/440.0f) + 49.0f; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			float yi[6803];//, yo[F];
            for (cf=0; cf<6803; cf++) { yi[cf] = 0.5f*cf; }
			interp1q_s(&frqs[0],F,&ihc_cfs_s[0],&yi[0],6803);
			//for (f=0; f<F; f++) { frqs[f] = yo[f]; }
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			float yi[64003];
			for (cf=0; cf<64003; cf++) { yi[cf] = 0.5f*cf; }
			interp1q_s(&frqs[0],F,&sgc_cfs_s[0],&yi[0],64003);
			//for (f=0; f<F; f++) { frqs[f] = yo[f]; }
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"mel",3)==0) //Mel
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0f * asinhf((7.0f/6.0f)*(powf(10.0f,frqs[f]/2595.0f)-1.0f)); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10f(700.0f*(powf(10.0f,frqs[f]/2595.0f)-1.0f)/165.4f + 0.88f) / 2.1f; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680f * logf(1.0f+46.065379f*(700.0f*powf(10.0f,frqs[f]/2595.0f)-700.0f)/((700.0f*powf(10.0f,frqs[f]/2595.0f)-700.0f)+14678.494617f)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 700.0f * (powf(10.0f,frqs[f]/2595.0f) - 1.0f); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f((700.0f*powf(10.0f,frqs[f]/2595.0f)-700.0f)/440.0f) + 69.0f; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2f(700.0f*powf(10.0f,frqs[f]/2595.0f)-700.0f); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*log2f((700.0f*powf(10.0f,frqs[f]/2595.0f)-700.0f)/440.0f) + 49.0f; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_s(frqs,F,"mel","hz"); convert_freqs_s(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_s(frqs,F,"mel","hz"); convert_freqs_s(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"midi",4)==0) //MIDI
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0f * asinhf(440.0f*powf(2.0f,(frqs[f]-69.0f)/12.0f)/600.0f); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10f(440.0f*powf(2.0f,(frqs[f]-69.0f)/12.0f)/165.4f+0.88f) / 2.1f; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680f * logf(1.0f+46.065379f*440.0f*powf(2.0f,(frqs[f]-69.0f)/12.0f)/(440.0f*powf(2.0f,(frqs[f]-69.0f)/12.0f)+14678.494617f)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 440.0f * powf(2.0f,(frqs[f]-69.0f)/12.0f); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0f * log10f(1.0f+440.0f*powf(2.0f,(frqs[f]-69.0f)/12.0f)/700.0f); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2f(440.0f) + (frqs[f]-69.0f)/12.0f; }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = frqs[f] - 20.0f; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_s(frqs,F,"midi","hz"); convert_freqs_s(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_s(frqs,F,"midi","hz"); convert_freqs_s(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"octave",6)==0) //Octaves
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0f * asinhf(powf(2.0f,frqs[f])/600.0f); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10f(powf(2.0f,frqs[f])/165.4f+0.88f) / 2.1f; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680f * logf(1.0f+46.065379f*powf(2.0f,frqs[f])/(powf(2.0f,frqs[f])+14678.494617f)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = powf(2.0f,frqs[f]); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0f * log10f(1.0f+powf(2.0f,frqs[f])/700.0f); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*(frqs[f]-log2f(440.0f)) + 69.0f; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0f*(frqs[f]-log2f(440.0f)) + 49.0f; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_s(frqs,F,"octave","hz"); convert_freqs_s(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_s(frqs,F,"octave","hz"); convert_freqs_s(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"piano",5)==0) //Piano
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0f * asinhf(440.0f*powf(2.0f,(frqs[f]-49.0f)/12.0f)/600.0f); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10f(440.0f*powf(2.0f,(frqs[f]-49.0f)/12.0f)/165.4f+0.88f) / 2.1f; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680f * logf(1.0f+46.065379f*440.0f*powf(2.0f,(frqs[f]-49.0f)/12.0f)/(440.0f*powf(2.0f,(frqs[f]-49.0f)/12.0f)+14678.494617f)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 440.0f * powf(2.0f,(frqs[f]-49.0f)/12.0f); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0f * log10f(1.0f+440.0f*powf(2.0f,(frqs[f]-49.0f)/12.0f)/700.0f); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = frqs[f] + 20.0f; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = (frqs[f]-49.0f)/12.0f + log2f(440.0f); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			return 0;
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_s(frqs,F,"piano","hz"); convert_freqs_s(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_s(frqs,F,"piano","hz"); convert_freqs_s(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
    else if (strncmp(in_scale,"ihcs",4)==0) //IHCs
	{
		if (strncmp(out_scale,"cochlea",7)==0)
		{
			convert_freqs_s(frqs,F,"ihcs","hz"); convert_freqs_s(frqs,F,"hz","cochlea");
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			float yi[6803];
			for (cf=0; cf<6803; cf++) { yi[cf] = 0.5f*cf; }
			interp1q_s(&frqs[0],F,&yi[0],&ihc_cfs_s[0],6803);
		}
		else if (strncmp(out_scale,"ihcs",4)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_s(frqs,F,"ihcs","hz"); convert_freqs_s(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"sgcs",4)==0) //SGCs
	{
		if (strncmp(out_scale,"cochlea",7)==0)
		{
			convert_freqs_s(frqs,F,"sgcs","hz"); convert_freqs_s(frqs,F,"hz","cochlea");
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			float yi[64003];
			for (cf=0; cf<64003; cf++) { yi[cf] = 0.5f*cf; }
            interp1q_s(&frqs[0],F,&yi[0],&sgc_cfs_s[0],64003);
		}
		else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_s(frqs,F,"sgcs","hz"); convert_freqs_s(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			return 0;
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else
	{
		fprintf(stderr,"error in convert_freqs: input frequency scale not recognized"); return 1;
	}
	
	return 0;
}


int convert_freqs_d (double *frqs, const int F, const char in_scale[], const char out_scale[])
{
	//Declaration
	int f, cf;

	//Checks
    if (F<1) { fprintf(stderr,"error in convert_freqs_d: F (num freqs to convert) must be positive\n"); return 1; }
	if (strlen(in_scale)<2) { fprintf(stderr,"error in convert_freqs_s: input freq scale must be string with length > 1\n"); return 1; }
	if (strlen(out_scale)<2) { fprintf(stderr,"error in convert_freqs_s: output freq scale must be string with length > 1\n"); return 1; }
	
	//Convert values in frqs
	if (strncmp(in_scale,"bark",4)==0) //Bark
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10((600.0*sinh(frqs[f]/6.0))/165.4+0.88) / 2.1; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680 * log(1.0+27639.2274*sinh(frqs[f]/6.0))/(600.0*sinh(frqs[f]/6.0)+14678.494617); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 600.0 * sinh(frqs[f]/6.0); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0 * log10(1.0 + (6.0/7.0)*sinh(frqs[f]/6.0)); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2((600.0/440.0)*sinh(frqs[f]/6.0)) + 69.0; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2(600.0*sinh(frqs[f]/6.0)); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2((600.0/440.0)*sinh(frqs[f]/6.0)) + 49.0; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_d(frqs,F,"bark","hz"); convert_freqs_d(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_d(frqs,F,"bark","hz"); convert_freqs_d(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"cochlea",7)==0) //Cochlea
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0 * asinh((165.4*(pow(10.0,2.1*frqs[f])-0.88))/600.0); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680 * log(1.0 + 7619.213687/(165.4+14678.494617/(pow(10.0,2.1*frqs[f])-0.88))); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 165.4 * (pow(10.0,2.1*frqs[f]) - 0.88); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0 * log10(1.0 + (165.4/700.0)*(pow(10.0,2.1*frqs[f])-0.88)); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2((165.4/440.0)*(pow(10.0,2.1*frqs[f])-0.88)) + 69.0; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2(165.4*(pow(10.0,2.1*frqs[f])-0.88)); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2((165.4/440.0)*(pow(10.0,2.1*frqs[f])-0.88)) + 49.0; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_d(frqs,F,"cochlea","hz"); convert_freqs_d(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_d(frqs,F,"cochlea","hz"); convert_freqs_d(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"erb",3)==0) //ERB
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0 * asinh((676170.419311/(47.065379-exp(0.089504*frqs[f]))-14678.494617)/600.0); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10((676170.419311/(47.065379-exp(0.089504*frqs[f]))-14678.494617)/165.4+0.88) / 2.1; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 676170.419311/(47.065379-exp(0.089504*frqs[f])) - 14678.494617; }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0 * log10(1.0+(676170.419311/(47.065379-exp(0.089504*frqs[f]))-14678.494617)/700.0); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2((676170.419311/(47.065379-exp(0.089504*frqs[f]))-14678.494617)/440.0) + 69.0; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2((676170.419311/(47.065379-exp(0.089504*frqs[f]))-14678.494617)); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2((676170.419311/(47.065379-exp(0.089504*frqs[f]))-14678.494617)/440.0) + 49.0; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_d(frqs,F,"erb","hz"); convert_freqs_d(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_d(frqs,F,"erb","hz"); convert_freqs_d(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"hz",2)==0) //Hz
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0 * asinh(frqs[f]/600.0); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10(frqs[f]/165.4+0.88) / 2.1; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680 * log(1.0+46.065379*frqs[f]/(frqs[f]+14678.494617)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0 * log10(1.0+frqs[f]/700.0); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2(frqs[f]/440.0) + 69.0; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2(frqs[f]); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2(frqs[f]/440.0) + 49.0; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			double yi[6803];
            for (cf=0; cf<6803; cf++) { yi[cf] = 0.5*cf; }
			interp1q_d(&frqs[0],F,&ihc_cfs_d[0],&yi[0],6803);
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			double yi[64003];
			for (cf=0; cf<64003; cf++) { yi[cf] = 0.5*cf; }
            interp1q_d(&frqs[0],F,&sgc_cfs_d[0],&yi[0],64003);
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"mel",3)==0) //Mel
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0 * asinh((7.0/6.0)*(pow(10.0,frqs[f]/2595.0)-1.0)); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10(700.0*(pow(10.0,frqs[f]/2595.0)-1.0)/165.4 + 0.88) / 2.1; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680 * log(1.0+46.065379*(700.0*pow(10.0,frqs[f]/2595.0)-700.0)/((700.0*pow(10.0,frqs[f]/2595.0)-700.0)+14678.494617)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 700.0 * (pow(10.0,frqs[f]/2595.0) - 1.0); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2((700.0*pow(10.0,frqs[f]/2595.0)-700.0)/440.0) + 69.0; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2(700.0*pow(10.0,frqs[f]/2595.0)-700.0); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*log2((700.0*pow(10.0,frqs[f]/2595.0)-700.0)/440.0) + 49.0; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_d(frqs,F,"mel","hz"); convert_freqs_d(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_d(frqs,F,"mel","hz"); convert_freqs_d(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"midi",4)==0) //MIDI
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0 * asinh(440.0*pow(2.0,(frqs[f]-69.0)/12.0)/600.0); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10(440.0*pow(2.0,(frqs[f]-69.0)/12.0)/165.4+0.88) / 2.1; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680 * log(1.0+46.065379*440.0*pow(2.0,(frqs[f]-69.0)/12.0)/(440.0*pow(2.0,(frqs[f]-69.0)/12.0)+14678.494617)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 440.0 * pow(2.0,(frqs[f]-69.0)/12.0); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0 * log10(1.0+440.0*pow(2.0,(frqs[f]-69.0)/12.0)/700.0); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log2(440.0) + (frqs[f]-69.0)/12.0; }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = frqs[f] - 20.0; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_d(frqs,F,"midi","hz"); convert_freqs_d(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_d(frqs,F,"midi","hz"); convert_freqs_d(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"octave",6)==0) //Octaves
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0 * asinh(pow(2.0,frqs[f])/600.0); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10(pow(2.0,frqs[f])/165.4+0.88) / 2.1; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680 * log(1.0+46.065379*pow(2.0,frqs[f])/(pow(2.0,frqs[f])+14678.494617)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = pow(2.0,frqs[f]); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0 * log10(1.0+pow(2.0,frqs[f])/700.0); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*(frqs[f]-log2(440.0)) + 69.0; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 12.0*(frqs[f]-log2(440.0)) + 49.0; }
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_d(frqs,F,"octave","hz"); convert_freqs_d(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_d(frqs,F,"octave","hz"); convert_freqs_d(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"piano",5)==0) //Piano
	{
		if (strncmp(out_scale,"bark",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 6.0 * asinh(440.0*pow(2.0,(frqs[f]-49.0)/12.0)/600.0); }
		}
		else if (strncmp(out_scale,"cochlea",7)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = log10(440.0*pow(2.0,(frqs[f]-49.0)/12.0)/165.4+0.88) / 2.1; }
		}
		else if (strncmp(out_scale,"erb",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 11.172680 * log(1.0+46.065379*440.0*pow(2.0,(frqs[f]-49.0)/12.0)/(440.0*pow(2.0,(frqs[f]-49.0)/12.0)+14678.494617)); }
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 440.0 * pow(2.0,(frqs[f]-49.0)/12.0); }
		}
		else if (strncmp(out_scale,"mel",3)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = 2595.0 * log10(1.0+440.0*pow(2.0,(frqs[f]-49.0)/12.0)/700.0); }
		}
		else if (strncmp(out_scale,"midi",4)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = frqs[f] + 20.0; }
		}
		else if (strncmp(out_scale,"octave",6)==0)
		{
			for (f=0; f<F; f++) { frqs[f] = (frqs[f]-49.0)/12.0 + log2(440.0); }
		}
		else if (strncmp(out_scale,"piano",5)==0)
		{
			return 0;
		}
        else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_d(frqs,F,"piano","hz"); convert_freqs_d(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_d(frqs,F,"piano","hz"); convert_freqs_d(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
    else if (strncmp(in_scale,"ihcs",4)==0) //IHCs
	{
		if (strncmp(out_scale,"cochlea",7)==0)
		{
			convert_freqs_d(frqs,F,"ihcs","hz"); convert_freqs_d(frqs,F,"hz","cochlea");
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			double yi[6803];
			for (cf=0; cf<6803; cf++) { yi[cf] = 0.5*cf; }
			interp1q_d(&frqs[0],F,&yi[0],&ihc_cfs_d[0],6803);
		}
		else if (strncmp(out_scale,"ihcs",4)==0)
		{
			return 0;
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			convert_freqs_d(frqs,F,"ihcs","hz"); convert_freqs_d(frqs,F,"hz","sgcs");
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else if (strncmp(in_scale,"sgcs",4)==0) //SGCs
	{
		if (strncmp(out_scale,"cochlea",7)==0)
		{
			convert_freqs_d(frqs,F,"sgcs","hz"); convert_freqs_d(frqs,F,"hz","cochlea");
		}
		else if (strncmp(out_scale,"hz",2)==0)
		{
			double yi[64003];
			for (cf=0; cf<64003; cf++) { yi[cf] = 0.5*cf; }
			interp1q_d(&frqs[0],F,&yi[0],&sgc_cfs_d[0],64003);
		}
		else if (strncmp(out_scale,"ihcs",4)==0)
		{
			convert_freqs_d(frqs,F,"sgcs","hz"); convert_freqs_d(frqs,F,"hz","ihcs");
		}
		else if (strncmp(out_scale,"sgcs",4)==0)
		{
			return 0;
		}
		else
		{
			fprintf(stderr,"error in convert_freqs: output frequency scale not recognized"); return 1;
		}
	}
	else
	{
		fprintf(stderr,"error in convert_freqs: input frequency scale not recognized"); return 1;
	}
	
	return 0;
}


#ifdef __cplusplus
}
}
#endif

