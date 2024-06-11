#include "common.h"
void DGC_STEP(void);

void DGC_FFT_FORWARD(int X);
void DGC_FFT_BACKWARD(int X);
void DGC_STEP_INIT(int x);
void DGC_STEP_CLEAN(void);

void int_PHA_Reimm(int k);
void int_PHA_Quad(int k);

void qInt_Graft(int x, double **qstar);
void qInt_Free(void);

void DGC_SETUP(void);
void DGC_CLEAN(void);

