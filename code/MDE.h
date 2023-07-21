#include "common.h"

#define _IJS(i,j,s) (j*(Nm_i[k][j]+1) + s)*Nx + i

void MDE_STEP(void);

void MDE_FFT_FORWARD(int k);
void MDE_FFT_BACKWARD(int k);
void MDE_STEP_INIT(int k);
	double *qInt, *wA, **qA, **qcA, **QA, **QcA, *Q1;
	double Q2;
	int *Ns;
void MDE_STEP_CLEAN(void);

void get_PHI(void);
void int_PHA_Reimm(int k);
void int_PHA_Quad(int k);

void qInt_Gauss(void);
void qInt_Disc(void);
void qInt_Free(void);

void MDE_SETUP(void);
void MDE_CLEAN(void);
