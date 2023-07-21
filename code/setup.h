#include "common.h"

void SETUP(void);
	char pattern[30], phname[30], phname_elec[30], Wname[30], itname[30];
void ASSIGN_CONSTS(void);
	double fugac_plus, fugac_minus, L_Deby_S, L_Deby_P;
	double eps_P, eps_S, zs, freeEnergy_bulk;
	int Nx, Nx2, Nxc, Ns0;
	double integ_cons;
void CALLOC_FIELDS(void);
	double *rx, *rx_sq;
	double *pot_elec, *eps_prof, *field_elec, *field_elec_sq;
	double *rho_elec_polym, *rho_elec_plus, *rho_elec_minus;
	double *phA, **PHA, **PHA_T, **WA;
	double *phB, *wB, *eta;
void LOAD_AA(char *seq_file);
	int NF_N, K, *K_i, Nm, **Nm_i, **Ns_i;
	double **Alpha_i, **Chi_i, *sigma_i;
void READ_AA(int *i, FILE *fp);
void INIT_FIELDS(int init_opt);
	double xCmax, xCltot, xCneck;
void ANALY_PB(void);
void get_ELFIELD(void);
void CLEAN(void);
