//Only compiles under C, since complex.

//To test compile:
//gcc -c cmp_ascend.c -O2 -std=c99 -Wall -Wextra
//clang -c cmp_ascend.c -O2 -std=c99 -Weverything

#include "/home/erik/codee/openvoice/openvoice.h"

#ifdef __cplusplus
extern "C" {
#endif


int cmp_ascend_s (const void *a, const void *b)
{
	float x1 = *(const float*)a, x2 = *(const float*)b;
	if (x1!=x1) { return 1; }
    else if (x2!=x2) { return -1; }
    else if (x1>x2) { return 1; }
    else if (x2>x1) { return -1; }
    else { return 0; }
}


int cmp_ascend_d (const void *a, const void *b)
{
	double x1 = *(const double*)a, x2 = *(const double*)b;
	if (x1!=x1) { return 1; }
    else if (x2!=x2) { return -1; }
    else if (x1>x2) { return 1; }
    else if (x2>x1) { return -1; }
    else { return 0; }
}


int cmp_ascend_c (const void *a, const void *b)
{
	float complex x1 = *(const float complex*)a, x2 = *(const float complex*)b;
	if (cabsf(x1)!=cabsf(x1)) { return 1; }
    else if (cabsf(x2)!=cabsf(x2)) { return -1; }
    else if (cabsf(x1)>cabsf(x2)) { return 1; }
    else if (cabsf(x2)>cabsf(x1)) { return -1; }
	else if (cargf(x1)>cargf(x2)) { return 1; }
    else if (cargf(x2)>cargf(x1)) { return -1; }
    else { return 0; }
}


int cmp_ascend_z (const void *a, const void *b)
{
	double complex x1 = *(const double complex*)a, x2 = *(const double complex*)b;
	if (cabs(x1)!=cabs(x1)) { return 1; }
    else if (cabs(x2)!=cabs(x2)) { return -1; }
    else if (cabs(x1)>cabs(x2)) { return 1; }
    else if (cabs(x2)>cabs(x1)) { return -1; }
	else if (carg(x1)>carg(x2)) { return 1; }
    else if (carg(x2)>carg(x1)) { return -1; }
    else { return 0; }
}


#ifdef __cplusplus
}
#endif

