#include "common.h"
double get_freeE(void);
void get_3D_GAUSSPROP(long int ijk);
void get_2D_GAUSSPROP(long int ij, double *mnGAUSSPROP, double *mGAUSSDIST, int N_m, double *nGAUSSDIST, int N_n);
void get_1D_GAUSSPROP(int i, double *mGAUSSPROP, double *mGAUSSDIST, int N);
void get_PHI(void);
void get_EPS_PROF(void);
void int_PHA_Reimm(int k);
void int_PHA_Quad(int k);

void get_int_Q1_Quen(int X);
void get_ELFIELD(void);

double And_mix(double **WA, double *wB); // No other place to put

void solve(int ndim, double **a, double *b, double *x);
void mult(int ndim, double **A, double *B, double *x);
void inv(int ndim, double **A, double **invA);
