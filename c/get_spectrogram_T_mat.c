//This makes the transformation matrix T for converting from a single frame of raw STFT power
//on a linear scale (from freqs) to a single frame of power on the new frequency scale.

//H is of size B x F, where B is the number of output frequency bands,
//and F is the number of non-negative STFT freqs (nfft/2+1).

//Each row of H has one triangle-shaped weighting function to make the weighted average.
//The width of the triangle is determined by B, the number of output bands;
//here I assume the usual default behavior where the lower edge of the first triangle ends on 0 Hz,
//and the upper edge of the Bth triangle ends on the Nyquist frq.
//Each row can be normalized to sum to 1.0, so that it represents a weighted average, by setting normalize to 1.

//The output freq_scale is previously determined by calling get_cfs_T,
//and then the cns (ints of length B) are determined here.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include "get_cns.c"

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int get_spectrogram_T_mat_s (float *H, const char iscolmajor, const float *freqs, const int F, const float *cfs, const int B, const char normalize);
int get_spectrogram_T_mat_d (double *H, const char iscolmajor, const double *freqs, const int F, const double *cfs, const int B, const char normalize);


int get_spectrogram_T_mat_s (float *H, const char iscolmajor, const float *freqs, const int F, const float *cfs, const int B, const char normalize)
{
	const float z = 0.0f;
    float Hsm;
    int b, f, f1, f2;
    int *cns;

    //Checks
    if (B<2) { fprintf(stderr,"error in get_spectrogram_T_mat_s: B (num cfs) must be > 1\n"); return 1; }
    if (F<2) { fprintf(stderr,"error in get_spectrogram_T_mat_s: F (num STFT freqs) must > 1\n"); return 1; }
    if (freqs[0]<0.0f) { fprintf(stderr,"error in get_spectrogram_T_mat_s: STFT freqs must be nonnegative\n"); return 1; }
    if (freqs[1]<=freqs[0] || freqs[F-1]<=freqs[F-2]) { fprintf(stderr,"error in get_spectrogram_T_mat_s: STFT freqs must be sorted ascending\n"); return 1; }
	
    //Get cns (center numbers, i.e. indices to freqs)
    if (!(cns=(int *)malloc((size_t)B*sizeof(int)))) { fprintf(stderr,"error in get_spectrogram_T_mat_s: problem with malloc for cns. "); perror("malloc"); return 1; }
	if (get_cns_s(cns,freqs,F,cfs,B)) { fprintf(stderr,"error in get_spectrogram_T_mat_s: problem getting cns\n"); return 1; }

    //Make sure H is initialized to 0
    cblas_scopy(B*F,&z,0,&H[0],1);

	//Make each row of H
    for (b=0; b<B; b++)
    {
        if (cns[b]<0 || cns[b]>=F) { fprintf(stderr,"error in get_spectrogram_T_mat_s: cns must be ints in [0 F-1]\n"); return 1; }
        if (b==0) { f1 = (2*cns[0]<1+cns[1]) ? 1 : 2*cns[0]-cns[1]; }
        else
        {
            if (cns[b]<=cns[b-1]) { fprintf(stderr,"error in get_spectrogram_T_mat_s: cns must be sorted ascending\n"); return 1; }
            f1 = cns[b-1];
        }
        if (b==B-1) { f2 = F - 1; } else { f2 = cns[b+1]; }

        if (iscolmajor)
        {
            for (f=f1+1; f<cns[b]; f++) { H[b+f*B] = (f-f1)*(1.0f/(cns[b]-f1)); }
            for (f=cns[b]; f<f2; f++) { H[b+f*B] = (f2-f)*(1.0f/(f2-cns[b])); }
            if (normalize)
            {
                Hsm = cblas_sasum(f2-f1-1,&H[b+(f1+1)*B],B);
                cblas_sscal(f2-f1-1,1.0f/Hsm,&H[b+(f1+1)*B],B);
            }
        }
        else
        {
            for (f=f1+1; f<cns[b]; f++) { H[f+b*F] = (f-f1)*(1.0f/(cns[b]-f1)); }
            for (f=cns[b]; f<f2; f++) { H[f+b*F] = (f2-f)*(1.0f/(f2-cns[b])); }
            if (normalize)
            {
                Hsm = cblas_sasum(f2-f1-1,&H[f1+1+b*F],1);
                cblas_sscal(f2-f1-1,1.0f/Hsm,&H[f1+1+b*F],1);
            }
        }
    }

    //Exit
	return 0;
}


int get_spectrogram_T_mat_d (double *H, const char iscolmajor, const double *freqs, const int F, const double *cfs, const int B, const char normalize)
{
	const double z = 0.0;
    double Hsm;
    int b, f, f1, f2;
    int *cns;

    //Checks
    if (B<2) { fprintf(stderr,"error in get_spectrogram_T_mat_d: B (num cfs) must be > 1\n"); return 1; }
    if (F<2) { fprintf(stderr,"error in get_spectrogram_T_mat_d: F (num STFT freqs) must > 1\n"); return 1; }
    if (freqs[0]<0.0) { fprintf(stderr,"error in get_spectrogram_T_mat_d: STFT freqs must be nonnegative\n"); return 1; }
    if (freqs[1]<=freqs[0] || freqs[F-1]<=freqs[F-2]) { fprintf(stderr,"error in get_spectrogram_T_mat_d: STFT freqs must be sorted ascending\n"); return 1; }

    //Get cns (center numbers, i.e. indices to freqs)
    if (!(cns=(int *)malloc((size_t)B*sizeof(int)))) { fprintf(stderr,"error in get_spectrogram_T_mat_d: problem with malloc for cns. "); perror("malloc"); return 1; }
	if (get_cns_d(cns,freqs,F,cfs,B)) { fprintf(stderr,"error in get_spectrogram_T_mat_d: problem getting cns\n"); return 1; }
	
    //Make sure H is initialized to 0
    cblas_dcopy(B*F,&z,0,&H[0],1);

	//Make each row of H
    for (b=0; b<B; b++)
    {
        if (cns[b]<0 || cns[b]>=F) { fprintf(stderr,"error in get_spectrogram_T_mat_s: cns must be ints in [0 F-1]\n"); return 1; }
        if (b==0) { f1 = (2*cns[0]<1+cns[1]) ? 1 : 2*cns[0]-cns[1]; }
        else
        {
            if (cns[b]<=cns[b-1]) { fprintf(stderr,"error in get_spectrogram_T_mat_s: cns must be sorted ascending\n"); return 1; }
            f1 = cns[b-1];
        }
        if (b==B-1) { f2 = F - 1; } else { f2 = cns[b+1]; }

        if (iscolmajor)
        {
            for (f=f1+1; f<cns[b]; f++) { H[b+f*B] = (f-f1)*(1.0/(cns[b]-f1)); }
            for (f=cns[b]; f<f2; f++) { H[b+f*B] = (f2-f)*(1.0/(f2-cns[b])); }
            if (normalize)
            {
                Hsm = cblas_dasum(f2-f1-1,&H[b+(f1+1)*B],B);
                cblas_dscal(f2-f1-1,1.0/Hsm,&H[b+(f1+1)*B],B);
            }
        }
        else
        {
            for (f=f1+1; f<cns[b]; f++) { H[f+b*F] = (f-f1)*(1.0/(cns[b]-f1)); }
            for (f=cns[b]; f<f2; f++) { H[f+b*F] = (f2-f)*(1.0/(f2-cns[b])); }
            if (normalize)
            {
                Hsm = cblas_dasum(f2-f1-1,&H[f1+1+b*F],1);
                cblas_dscal(f2-f1-1,1.0/Hsm,&H[f1+1+b*F],1);
            }
        }
    }
	
    //Exit
	return 0;
}


#ifdef __cplusplus
}
}
#endif

