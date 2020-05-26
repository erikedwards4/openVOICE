//Gets F frequencies in Hz for the STFT of length nfft with sample rate fs.
//The frequencies are just linearly spaced from 0 to Nyquist (fs/2).
//If F>nfft/2+1, then the negative frequencies are included.
//Only up to F frequencies are output, and freqs must be of length F.

//To test compile:
//gcc -c get_stft_freqs.c -O2 -std=c99 -Wall -Wextra
//clang -c get_stft_freqs.c -O2 -std=c99 -Weverything
//g++ -c get_stft_freqs.c -O2 -std=c++11 -Wall -Wextra
//clang++ -c get_stft_freqs.c -O2 -std=c++11 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int get_stft_freqs_s (float *freqs, const int F, const float fs, const int nfft)
{
    const int Nyq = nfft/2;
    const float finc = 0.5f*fs/Nyq;
    int f;

    //Checks
    if (F<1) { fprintf(stderr,"error in get_stft_freqs_s: F must be positive\n"); return 1; }
    if (nfft<1) { fprintf(stderr,"error in get_stft_freqs_s: nfft must be positive\n"); return 1; }
    if (F>nfft) { fprintf(stderr,"error in get_stft_freqs_s: F must be <= nfft (STFT length)\n"); return 1; }
    if (fs<=0.0f) { fprintf(stderr,"error in get_stft_freqs_s: fs must be positive\n"); return 1; }
	
    freqs[0] = 0.0f;
    for (f=1; f<F && f<Nyq; f++) { freqs[f] = f*finc; }
    if (F>Nyq) { freqs[Nyq] = 0.5f*fs; }
    for (f=Nyq+1; f<F; f++) { freqs[f] = -finc*(nfft-f); }
	
    return 0;
}


int get_stft_freqs_d (double *freqs, const int F, const double fs, const int nfft)
{
    const int Nyq = nfft/2;
    const double finc = 0.5*fs/Nyq;
    int f;

    //Checks
    if (F<1) { fprintf(stderr,"error in get_stft_freqs_d: F must be positive\n"); return 1; }
    if (nfft<1) { fprintf(stderr,"error in get_stft_freqs_d: nfft must be positive\n"); return 1; }
    if (F>nfft) { fprintf(stderr,"error in get_stft_freqs_d: F must be <= nfft (STFT length)\n"); return 1; }
    if (fs<=0.0) { fprintf(stderr,"error in get_stft_freqs_d: fs must be positive\n"); return 1; }
	
    freqs[0] = 0.0;
    for (f=1; f<F && f<Nyq; f++) { freqs[f] = f*finc; }
    if (F>Nyq) { freqs[Nyq] = 0.5*fs; }
    for (f=Nyq+1; f<F; f++) { freqs[f] = -finc*(nfft-f); }
	
    return 0;
}


#ifdef __cplusplus
}
#endif

