#include "common.h"
#include </share/apps/software/fftw-3.3.8/include/fftw3.h> // Pitzer cluster fftw3.3.8 location

void FFTW_SETUP(void);
	fftw_plan p, ip, sfp, sip, cfp, cip;
	double *fftw_in, *fftw_inc, *fftw_ins;
	fftw_complex *fftw_out_c;
	double *fftw_outs_r, *fftw_outc_r;
void FFT_solve(void);
void MDE_FFT(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, fftw_complex *fftw_out);
void MDE_COS(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out);
void MDE_SIN(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out);
void PB_FFT(double *phA, double **PHA, double **PHA_T, double *fftw_in, fftw_complex *fftw_out);
void PB_COS(double *phA, double **PHA, double **PHA_T, double *fftw_in, double *fftw_out);
void FFTW_CLEAN(void);
