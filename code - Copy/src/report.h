#include "common.h"

void REPORT_SETUP(void);
	char pattern[30];
	char phname[30];
	char phname_elec[30];
	char Wname[30];
	char itname[30];
void write_ph(void);
void write_elec(void);
void write_W(void);
void write_qA(char *name, int N, int Ns, double **qA);
void write_it(void);
