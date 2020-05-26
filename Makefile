#@author Erik Edwards
#@date 2019-2020

#openVOICE is my own library of C functions for voice processing in C.
#This is the Makefile to make static (.a) and dynamic (.so) libraries.

#To use the .o object files, or the .a or .so libraries:
#Do this once: sudo ln -s $(OV)/libopenvoice.a /usr/lib/libopenvoice.a
#Or this once: sudo ln -s $(OV)/libopenvoice.so /usr/lib/libopenvoice.so

#For libopenvoice.a:
#ar r creates new archive or adds file(s) to archive
#ar c suppresses stderr message when creating archive
#ar s forces regeneration of symbol table
#ar v is for verbose

SHELL=/bin/bash

OV=/home/erik/codee/openvoice

CC=gcc

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

#Use these (via make a or make so) to make a final library
AR=ar crs libopenvoice.a $(OV)/{basic,dsp,pre,wins,freqs,window,stft,spectrogram,ccs,deltas,ar_poly,ac_lp,lsf_lar,zcs_lcs}/obj/*.o
SO=$(CC) -shared -o libopenvoice.so $(OV)/{basic,dsp,pre,wins,freqs,window,stft,spectrogram,ccs,deltas,ar_poly,ac_lp,lsf_lar,zcs_lcs}/obj/*.o

openVOICE: Basic DSP Pre Wins Freqs Window STFT Spectrogram CCs Deltas AR_Poly AC_LP LSF_LAR ZCs_LCs #Functionals 
	$(AR)
	$(SO)

Basic: mean0 stdev1 zscore abs square cmp_ascend cmp_descend #median
mean0: basic/src/mean0.c; $(CC) -c $^ -o basic/obj/$@.o $(CFLAGS)
stdev1: basic/src/stdev1.c; $(CC) -c $^ -o basic/obj/$@.o $(CFLAGS)
zscore: basic/src/zscore.c; $(CC) -c $^ -o basic/obj/$@.o $(CFLAGS)
abs: basic/src/abs.c; $(CC) -c $^ -o basic/obj/$@.o $(CFLAGS)
square: basic/src/square.c; $(CC) -c $^ -o basic/obj/$@.o $(CFLAGS)
cmp_ascend: basic/src/cmp_ascend.c; $(CC) -c $^ -o basic/obj/$@.o $(CFLAGS)
cmp_descend: basic/src/cmp_descend.c; $(CC) -c $^ -o basic/obj/$@.o $(CFLAGS)
median: basic/src/median.c; $(CC) -c $^ -o basic/obj/$@.o $(CFLAGS)

DSP: fft fir iir #medfilt fft_matmul
fft: dsp/src/fft.c; $(CC) -c $^ -o dsp/obj/$@.o $(CFLAGS)
fft_matmul: dsp/src/fft_matmul.c; $(CC) -c $^ -o dsp/obj/$@.o $(CFLAGS)
fir: dsp/src/fir.c; $(CC) -c $^ -o dsp/obj/$@.o $(CFLAGS)
iir: dsp/src/iir.c; $(CC) -c $^ -o dsp/obj/$@.o $(CFLAGS)
medfilt: dsp/src/medfilt.c; $(CC) -c $^ -o dsp/obj/$@.o $(CFLAGS)

Pre: rms_scale preemph #dither
rms_scale: pre/src/rms_scale.c; $(CC) -c $^ -o pre/obj/$@.o $(CFLAGS)
preemph: pre/src/preemph.c; $(CC) -c $^ -o pre/obj/$@.o $(CFLAGS)
#dither: pre/src/dither.c; $(CC) -c $^ -o pre/obj/$@.o $(CFLAGS)

Wins: rectangular triangular bartlett hann hamming blackman blackmanharris flattop povey #tukey
rectangular: wins/src/rectangular.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)
triangular: wins/src/triangular.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)
bartlett: wins/src/bartlett.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)
hann: wins/src/hann.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)
hamming: wins/src/hamming.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)
blackman: wins/src/blackman.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)
blackmanharris: wins/src/blackmanharris.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)
flattop: wins/src/flattop.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)
povey: wins/src/povey.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)
tukey: wins/src/tukey.c; $(CC) -c $^ -o wins/obj/$@.o $(CFLAGS)

Freqs: interp1q convert_freqs get_cfs get_stft_freqs get_cfs_T #get_cns
interp1q: freqs/src/interp1q.c; $(CC) -c $^ -o freqs/obj/$@.o $(CFLAGS)
convert_freqs: freqs/src/convert_freqs.c; $(CC) -c $^ -o freqs/obj/$@.o $(CFLAGS)
get_cfs: freqs/src/get_cfs.c; $(CC) -c $^ -o freqs/obj/$@.o $(CFLAGS)
get_stft_freqs: freqs/src/get_stft_freqs.c; $(CC) -c $^ -o freqs/obj/$@.o $(CFLAGS)
get_cfs_T: freqs/src/get_cfs_T.c; $(CC) -c $^ -o freqs/obj/$@.o $(CFLAGS)

Window: frame_univar apply_win window_univar
frame_univar: window/src/frame_univar.c; $(CC) -c window/src/$@.c -o window/obj/$@.o $(CFLAGS)
apply_win: window/src/apply_win.c; $(CC) -c window/src/$@.c -o window/obj/$@.o $(CFLAGS)
window_univar: window/src/window_univar.c; $(CC) -c window/src/$@.c -o window/obj/$@.o $(CFLAGS)

STFT: fft_hc hc_square fft_squared stft #stft2
fft_hc: stft/src/fft_hc.c; $(CC) -c $^ -o stft/obj/$@.o $(CFLAGS)
hc_square: stft/src/hc_square.c; $(CC) -c $^ -o stft/obj/$@.o $(CFLAGS)
fft_squared: stft/src/fft_squared.c; $(CC) -c $^ -o stft/obj/$@.o $(CFLAGS)
stft: stft/src/stft.c; $(CC) -c $^ -o stft/obj/$@.o $(CFLAGS)
stft2: stft/src/stft2.c; $(CC) -c $^ -o stft/obj/$@.o $(CFLAGS)

Spectrogram: get_cns get_spectrogram_T_mat apply_spectrogram_T_mat pow_compress spectrogram
get_cns: spectrogram/src/get_cns.c; $(CC) -c $^ -o spectrogram/obj/$@.o $(CFLAGS)
abs2: spectrogram/src/abs2.c; $(CC) -c $^ -o spectrogram/obj/$@.o $(CFLAGS)
get_spectrogram_T_mat: spectrogram/src/get_spectrogram_T_mat.c; $(CC) -c $^ -o spectrogram/obj/$@.o $(CFLAGS)
apply_spectrogram_T_mat: spectrogram/src/apply_spectrogram_T_mat.c; $(CC) -c $^ -o spectrogram/obj/$@.o $(CFLAGS)
pow_compress: spectrogram/src/pow_compress.c; $(CC) -c $^ -o spectrogram/obj/$@.o $(CFLAGS)
spectrogram: spectrogram/src/spectrogram.c; $(CC) -c $^ -o spectrogram/obj/$@.o $(CFLAGS)

CCs: dct dct_inplace lifter get_ccs mfccs
dct: ccs/src/dct.c; $(CC) -c $^ -o ccs/obj/$@.o $(CFLAGS)
dct_inplace: ccs/src/dct_inplace.c; $(CC) -c $^ -o ccs/obj/$@.o $(CFLAGS)
lifter: ccs/src/lifter.c; $(CC) -c $^ -o ccs/obj/$@.o $(CFLAGS)
get_ccs: ccs/src/get_ccs.c; $(CC) -c $^ -o ccs/obj/$@.o $(CFLAGS)
mfccs: ccs/src/mfccs.c; $(CC) -c $^ -o ccs/obj/$@.o $(CFLAGS)

Deltas: get_deltas add_deltas get_delta_deltas add_delta_deltas
get_deltas: deltas/src/get_deltas.c; $(CC) -c $^ -o deltas/obj/$@.o $(CFLAGS)
add_deltas: deltas/src/add_deltas.c; $(CC) -c $^ -o deltas/obj/$@.o $(CFLAGS)
get_delta_deltas: deltas/src/get_delta_deltas.c; $(CC) -c $^ -o deltas/obj/$@.o $(CFLAGS)
add_delta_deltas: deltas/src/add_delta_deltas.c; $(CC) -c $^ -o deltas/obj/$@.o $(CFLAGS)

AR_Poly: roots poly poly2roots roots2poly ar2poly poly2ar ar2rc rc2ar poly2rc rc2poly ar2psd poly2psd
roots: ar_poly/src/roots.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
poly: ar_poly/src/poly.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
poly2roots: ar_poly/src/poly2roots.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
roots2poly: ar_poly/src/roots2poly.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
ar2poly: ar_poly/src/ar2poly.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
poly2ar: ar_poly/src/poly2ar.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
ar2rc: ar_poly/src/ar2rc.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
rc2ar: ar_poly/src/rc2ar.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
poly2rc: ar_poly/src/poly2rc.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
rc2poly: ar_poly/src/rc2poly.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
ar2psd: ar_poly/src/ar2psd.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)
poly2psd: ar_poly/src/poly2psd.c; $(CC) -c $^ -o ar_poly/obj/$@.o $(CFLAGS)

AC_LP: autocorr autocorr_fft sig2ac sig2ac_fft ac2ar_levdurb ac2poly_levdurb sig2ar_levdurb sig2poly_levdurb sig2ar_burg sig2poly_burg ac2rc ac2cc ac2mvdr
autocorr: ac_lp/src/autocorr.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
autocorr_fft: ac_lp/src/autocorr_fft.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
sig2ac: ac_lp/src/sig2ac.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
sig2ac_fft: ac_lp/src/sig2ac_fft.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
ac2ar_levdurb: ac_lp/src/ac2ar_levdurb.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
ac2poly_levdurb: ac_lp/src/ac2poly_levdurb.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
sig2ar_levdurb: ac_lp/src/sig2ar_levdurb.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
sig2poly_levdurb: ac_lp/src/sig2poly_levdurb.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
sig2ar_burg: ac_lp/src/sig2ar_burg.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
sig2poly_burg: ac_lp/src/sig2poly_burg.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
ac2rc: ac_lp/src/ac2rc.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
ac2cc: ac_lp/src/ac2cc.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)
ac2mvdr: ac_lp/src/ac2mvdr.c; $(CC) -c $^ -o ac_lp/obj/$@.o $(CFLAGS)

LSF_LAR: ar2lsf lsf2ar poly2lsf lsf2poly rc2lar lar2rc
ar2lsf: lsf_lar/src/ar2lsf.c; $(CC) -c $^ -o lsf_lar/obj/$@.o $(CFLAGS)
lsf2ar: lsf_lar/src/lsf2ar.c; $(CC) -c $^ -o lsf_lar/obj/$@.o $(CFLAGS)
poly2lsf: lsf_lar/src/poly2lsf.c; $(CC) -c $^ -o lsf_lar/obj/$@.o $(CFLAGS)
lsf2poly: lsf_lar/src/lsf2poly.c; $(CC) -c $^ -o lsf_lar/obj/$@.o $(CFLAGS)
rc2lar: lsf_lar/src/rc2lar.c; $(CC) -c $^ -o lsf_lar/obj/$@.o $(CFLAGS)
lar2rc: lsf_lar/src/lar2rc.c; $(CC) -c $^ -o lsf_lar/obj/$@.o $(CFLAGS)

ZCs_LCs: zcs lcs mcs zcr lcr mcr zcr_windowed lcr_windowed mcr_windowed
zcs: zcs_lcs/src/zcs.c; $(CC) -c $^ -o zcs_lcs/obj/$@.o $(CFLAGS)
lcs: zcs_lcs/src/lcs.c; $(CC) -c $^ -o zcs_lcs/obj/$@.o $(CFLAGS)
mcs: zcs_lcs/src/mcs.c; $(CC) -c $^ -o zcs_lcs/obj/$@.o $(CFLAGS)
zcr: zcs_lcs/src/zcr.c; $(CC) -c $^ -o zcs_lcs/obj/$@.o $(CFLAGS)
lcr: zcs_lcs/src/lcr.c; $(CC) -c $^ -o zcs_lcs/obj/$@.o $(CFLAGS)
mcr: zcs_lcs/src/mcr.c; $(CC) -c $^ -o zcs_lcs/obj/$@.o $(CFLAGS)
zcr_windowed: zcs_lcs/src/zcr_windowed.c; $(CC) -c $^ -o zcs_lcs/obj/$@.o $(CFLAGS)
lcr_windowed: zcs_lcs/src/lcr_windowed.c; $(CC) -c $^ -o zcs_lcs/obj/$@.o $(CFLAGS)
mcr_windowed: zcs_lcs/src/mcr_windowed.c; $(CC) -c $^ -o zcs_lcs/obj/$@.o $(CFLAGS)

Functionals: moments #cumulants prctiles lcrs
moments: functionals/src/moments.c; $(CC) -c $^ -o functionals/obj/$@.o $(CFLAGS)
cumulants: functionals/src/cumulants.c; $(CC) -c $^ -o functionals/obj/$@.o $(CFLAGS)
prctiles: functionals/src/prctiles.c; $(CC) -c $^ -o functionals/obj/$@.o $(CFLAGS)
lcrs: functionals/src/lcrs.c; $(CC) -c $^ -o functionals/obj/$@.o $(CFLAGS)

#Use these (i.e. make a or make so) to make a final library
a: libopenvoice.a; $(AR)
so: libopenvoice.so; $(SO)

clean:
	find $(OV) -type f -name *.o | xargs rm -f

