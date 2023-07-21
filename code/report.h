#include "common.h"

void REPORT_SETUP(void);
	char pattern[30];
	char phname[30];
	char phname_elec[30];
	char Wname[30];
	char itname[30];
void write_ph(double *phA, double **PHA_T, double **PHA, double *phB);
void write_elec(double *pot_elec,double *rho_elec_plus,double *rho_elec_minus,double *rho_elec_polym);
void write_W(double **WA,double *wB,double *eta); 
void write_qA(char *name, int N, int Ns, double **qA);
