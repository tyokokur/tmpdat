#include "common.h"
void report_vec(int i, double *vec);
void stdout_redir_setup(void);
int  stdout_redir_start(void);
void stdout_redir_refresh(void);
void stdout_redir_end(void);
void printout(void);
void REPORT_SETUP(void);
void write_ph(void);
void write_elec(void);
void write_W(void);
void write_qA(char *name, int N, int Ns, double **qA);
void write_qA_s(int s, char *name, int N, int Ns, double **qA);
void write_it(void);
void write_finish(void);
