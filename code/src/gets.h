#include "common.h"

void get_freeE(void);
	double freeEnergy, freeDiff, inCompMax;
	int MAXMAX;

void get_PHI(void);
void int_PHA_Reimm(int k);
void int_PHA_Quad(int k);

void get_ELFIELD(void);

double And_mix(double **WA, double *wB); // No other place to put
	double and_err;
	double **WAdiff, **WAnew, **DAnew;
	double **u, **u_temp, *v;
	double *wBdiff, *wBnew, *DBnew, *Cs;

void solve(int ndim, double **a, double *b, double *x);
void mult(int ndim, double **A, double *B, double *x);
void inv(int ndim, double **A, double **invA);
