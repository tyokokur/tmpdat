#include "common.h"
void MDE_STEP(void);

void MDE_FFT_FORWARD(int k);
void MDE_FFT_BACKWARD(int k);
void MDE_STEP_INIT(int k);
void MDE_STEP_CLEAN(void);

void get_PHI(void);
void int_PHA_Reimm(int k);
void int_PHA_Quad(int k);

void qInt_Graft(int x, double **qstar);
void qInt_Free(void);


void MDE_SETUP(void);
void MDE_CLEAN(void);

