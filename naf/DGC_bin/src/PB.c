#include "PB.h"
#include "FFT_solve.h"
#include "gets.h"

void PB_STEP(void) {
	get_EPS_PROF();

double tic = omp_get_wtime();
	if (HALF_BOOL) PB_COS(phA, PHA, PHA_T, fftw_inc, fftw_outc_r);
	else           PB_FFT(phA, PHA, PHA_T, fftw_in,  fftw_out_c);
double toc = omp_get_wtime();
//printf("\t[PB.c/PB_COS]: %f\n", toc-tic);

	get_ELFIELD();
}
