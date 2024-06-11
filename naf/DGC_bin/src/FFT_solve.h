#include "common.h"
#include <fftw3.h>

void FFTW_SETUP(void);
void DGC_COS(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out, double *PROPGAUSS, double *DG_in, double *DG_out);
void FFT_solve(void);
void MDE_FFT(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, fftw_complex *fftw_out);
void MDE_COS(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out);
void MDE_COSxy(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out);
void MDE_SIN(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out);
void PB_FFT(double *phA, double **PHA, double **PHA_T, double *fftw_in, fftw_complex *fftw_out);
void PB_COS(double *phA, double **PHA, double **PHA_T, double *fftw_in, double *fftw_out);
void FFTW_CLEAN(void);
