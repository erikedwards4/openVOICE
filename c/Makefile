#@author Erik Edwards
#@date 2019-2020

#openVOICE is my own library of C functions for voice processing in C.
#This is the Makefile used during development of the C functions.
#Since these are not stand-alone programs (no main functions),
#the only outputs are object files to obj.

#This document is also important since it defines the project structure,
#that is the groupings of the C functions.

SHELL=/bin/bash

CC=g++

ifeq ($(CC),$(filter $(CC),clang gcc))
	STD=-std=c99
else
	STD=-std=c++11
endif

ifeq ($(CC),clang)
	WFLAGS=-Weverything -Wno-old-style-cast
else ifeq ($(CC),clang++)
	WFLAGS=-Weverything -Wno-old-style-cast -Wno-deprecated
else
	WFLAGS=-Wall -Wextra
endif

CFLAGS=$(WFLAGS) -O3 $(STD) -march=native -fopenmp -fPIC


all: Basic DSP Pre Wins Freqs Window STFT Spectrogram CCs Deltas AR_Poly AC_LP LSF_LAR ZCs_LCs #Functionals

Basic: mean0 stdev1 zscore abs square cmp_ascend cmp_descend median
mean0: mean0.c; $(CC) -c mean0.c -o obj/mean0.o $(CFLAGS)
stdev1: stdev1.c; $(CC) -c stdev1.c -o obj/stdev1.o $(CFLAGS)
zscore: zscore.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
abs: abs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
square: square.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
cmp_ascend: cmp_ascend.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
cmp_descend: cmp_descend.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
median: median.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

DSP: fft fir iir #medfilt fft_matmul
fft: fft.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
fft_matmul: fft_matmul.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
fir: fir.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
iir: iir.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
medfilt: medfilt.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

Pre: rms_scale preemph #dither
rms_scale: rms_scale.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
preemph: preemph.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
#dither: dither.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

Wins: rectangular triangular bartlett hann hamming blackman blackmanharris flattop povey #tukey
rectangular: rectangular.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
triangular: triangular.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
bartlett: bartlett.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
hann: hann.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
hamming: hamming.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
blackman: blackman.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
blackmanharris: blackmanharris.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
flattop: flattop.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
povey: povey.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
tukey: tukey.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

Freqs: interp1q convert_freqs get_cfs get_stft_freqs get_cfs_T #get_cns
interp1q: interp1q.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
convert_freqs: convert_freqs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
get_cfs: get_cfs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
get_stft_freqs: get_stft_freqs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
get_cfs_T: get_cfs_T.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

Window: frame_univar apply_win window_univar
frame_univar: frame_univar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
apply_win: apply_win.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
window_univar: window_univar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

STFT: fft_hc hc_square fft_squared stft #stft2
fft_hc: fft_hc.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
hc_square: hc_square.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
fft_squared: fft_squared.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
stft: stft.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
stft2: stft2.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

Spectrogram: get_cns get_spectrogram_T_mat apply_spectrogram_T_mat pow_compress spectrogram
get_cns: get_cns.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
abs2: abs2.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
get_spectrogram_T_mat: get_spectrogram_T_mat.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
apply_spectrogram_T_mat: apply_spectrogram_T_mat.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
pow_compress: pow_compress.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
spectrogram: spectrogram.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

CCs: dct dct_inplace lifter get_ccs mfccs
dct: dct.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
dct_inplace: dct_inplace.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
lifter: lifter.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
get_ccs: get_ccs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
mfccs: mfccs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

Deltas: get_deltas add_deltas get_delta_deltas add_delta_deltas
get_deltas: get_deltas.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
add_deltas: add_deltas.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
get_delta_deltas: get_delta_deltas.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
add_delta_deltas: add_delta_deltas.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

AR_Poly: roots poly poly2roots roots2poly ar2poly poly2ar ar2rc rc2ar poly2rc rc2poly ar2psd poly2psd
roots: roots.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
poly: poly.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
poly2roots: poly2roots.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
roots2poly: roots2poly.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
ar2poly: ar2poly.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
poly2ar: poly2ar.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
ar2rc: ar2rc.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
rc2ar: rc2ar.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
poly2rc: poly2rc.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
rc2poly: rc2poly.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
ar2psd: ar2psd.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
poly2psd: poly2psd.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

AC_LP: autocorr autocorr_fft sig2ac sig2ac_fft ac2ar_levdurb ac2poly_levdurb sig2ar_levdurb sig2poly_levdurb sig2ar_burg sig2poly_burg ac2rc ac2cc ac2mvdr
autocorr: autocorr.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
autocorr_fft: autocorr_fft.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
sig2ac: sig2ac.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
sig2ac_fft: sig2ac_fft.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
ac2ar_levdurb: ac2ar_levdurb.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
ac2poly_levdurb: ac2poly_levdurb.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
sig2ar_levdurb: sig2ar_levdurb.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
sig2poly_levdurb: sig2poly_levdurb.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
sig2ar_burg: sig2ar_burg.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
sig2poly_burg: sig2poly_burg.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
ac2rc: ac2rc.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
ac2cc: ac2cc.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
ac2mvdr: ac2mvdr.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

LSF_LAR: ar2lsf lsf2ar poly2lsf lsf2poly rc2lar lar2rc
ar2lsf: ar2lsf.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
lsf2ar: lsf2ar.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
poly2lsf: poly2lsf.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
lsf2poly: lsf2poly.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
rc2lar: rc2lar.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
lar2rc: lar2rc.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

ZCs_LCs: zcs lcs mcs zcr lcr mcr zcr_windowed lcr_windowed mcr_windowed
zcs: zcs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
lcs: lcs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
mcs: mcs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
zcr: zcr.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
lcr: lcr.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
mcr: mcr.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
zcr_windowed: zcr_windowed.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
lcr_windowed: lcr_windowed.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
mcr_windowed: mcr_windowed.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)

Functionals: moments #cumulants prctiles lcrs
moments: moments.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
cumulants: cumulants.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
prctiles: prctiles.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)
lcrs: lcrs.c; $(CC) -c $^ -o obj/$@.o $(CFLAGS)


clean:
	find ./obj -type f -name *.o | xargs rm -f

