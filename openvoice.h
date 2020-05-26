#pragma once


#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>
#include <time.h>


#ifndef M_PIf
   #define M_PIf 3.141592653589793238462643383279502884f
#endif

#ifndef M_PI
   #define M_PI 3.141592653589793238462643383279502884
#endif


//Define struct for minimum input info (3 ints)
// struct iii
// {
//    int iscolmajor;
// }

//Define struct for 4D tensors with critical info
//This is an intentional play on words ("Tesseract", "Tensor", "Time-series/signal + extra info")
// struct Tesser_s
// {
//    float *dat;
//    char iscolmajor;
//    int R;
//    int C;
//    int S;
//    int H;
// };

// struct Tesser_d
// {
//    double *dat;
//    char iscolmajor;
//    int R;
//    int C;
//    int S;
//    int H;
// };


#ifdef __cplusplus
namespace ov {
#endif

#include "/home/erik/codee/openvoice/basic/basic.h"
#include "/home/erik/codee/openvoice/dsp/dsp.h"
#include "/home/erik/codee/openvoice/pre/pre.h"
#include "/home/erik/codee/openvoice/wins/wins.h"
#include "/home/erik/codee/openvoice/freqs/freqs.h"
#include "/home/erik/codee/openvoice/window/window.h"
#include "/home/erik/codee/openvoice/stft/stft.h"
#include "/home/erik/codee/openvoice/spectrogram/spectrogram.h"
#include "/home/erik/codee/openvoice/ccs/ccs.h"
#include "/home/erik/codee/openvoice/deltas/deltas.h"
#include "/home/erik/codee/openvoice/ar_poly/ar_poly.h"
#include "/home/erik/codee/openvoice/ac_lp/ac_lp.h"
#include "/home/erik/codee/openvoice/lsf_lar/lsf_lar.h"
#include "/home/erik/codee/openvoice/zcs_lcs/zcs_lcs.h"
#include "/home/erik/codee/openvoice/functionals/functionals.h"

#ifdef __cplusplus
}
#endif

