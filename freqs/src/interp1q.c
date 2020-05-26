//These do a quick linear interpolation (required by convert_freqs for biofreq scales).
//These assume no extrapolation, and assume that Xi is sorted ascending.
//These are also in-place (overwrite Xo with Yo).

//To test compile:
//gcc -c interp1q.c -O2 -std=c99 -Wall -Wextra
//clang -c interp1q.c -O2 -std=c99 -Weverything
//g++ -c interp1q.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c interp1q.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int interp1q_s (float *Xo, const int No, const float *Xi, const float *Yi, const int Ni)
{
	int ni = 0, no = 0;

	//Checks
    if (Ni<2) { fprintf(stderr,"error in interp1q_s: Ni (input length) must be > 1\n"); return 1; }
    if (No<1) { fprintf(stderr,"error in interp1q_s: No (output length) must be positive\n"); return 1; }

    while (no<No)
    {
        ni = 0; //disclude this for slightly faster, but require Xo to be ascending
        while (ni<Ni && Xi[ni]<Xo[no]) { ni++; }
        if (fabsf(Xi[ni]-Xo[no])<FLT_EPSILON) { Xo[no] = Yi[ni]; }
        else { Xo[no] = Yi[ni-1] + (Yi[ni]-Yi[ni-1])*(Xo[no]-Xi[ni-1])/(Xi[ni]-Xi[ni-1]); }
        no++;
    }
	
	return 0;
}


int interp1q_d (double *Xo, const int No, const double *Xi, const double *Yi, const int Ni)
{
	int ni = 0, no = 0;

	//Checks
    if (Ni<2) { fprintf(stderr,"error in interp1q_s: Ni (input length) must be > 1\n"); return 1; }
    if (No<1) { fprintf(stderr,"error in interp1q_s: No (output length) must be positive\n"); return 1; }

    while (no<No)
    {
        ni = 0; //disclude this for slightly faster, but require Xo to be ascending
        while (ni<Ni && Xi[ni]<Xo[no]) { ni++; }
        if (fabs(Xi[ni]-Xo[no])<DBL_EPSILON) { Xo[no] = Yi[ni]; }
        else { Xo[no] = Yi[ni-1] + (Yi[ni]-Yi[ni-1])*(Xo[no]-Xi[ni-1])/(Xi[ni]-Xi[ni-1]); }
        no++;
    }
	
	return 0;
}


#ifdef __cplusplus
}
#endif

