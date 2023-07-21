/*Stand Libs*/
#ifndef _STANDARDS
	#define _STANDARDS
	#include <stdlib.h>
	#include <stdio.h>
	#include <string.h>
	#define _USE_MATH_DEFINES
	#include <math.h>
#endif

/*Memory Macros*/
#ifndef TJY_MEM
	#define TJY_MEM
	#define INIT(size) (double *)calloc(size, sizeof(double))
	#define INIT_2D(NAME, size1, size2) (double **)calloc(size1, sizeof(double)); for (i=0; i<size1; i++) NAME[i] = INIT(size2)
#endif

/*Print Macros*/
#ifndef TJY_PRINT
	#define TJY_PRINT
	#define CHECK_VAL(VAR) printf("%.6e ", VAR)
	#define CHECK_VAL_N(VAR) printf(#VAR": %.6e\n", VAR)
#endif

// Physical
extern double eps_vacuum, kB, Length, Charge, avo_num;
// System 
extern double freeEnergy;
extern int Temp;
extern double lx;
extern int Nx, Nx2, Nxc;
extern double *rx, *rx_sq;
// Solution
extern int Z_plus, Z_minus;
extern double eps_r_P, eps_r_S, Q2, freeEnergy_bulk;
extern double c0_plus, c0_minus, fugac_plus, fugac_minus; 
extern double eps_P, eps_S, L_Deby_S, L_Deby_P, mu_s, zs;
extern double *eps_prof, *pot_elec, *rho_elec_minus, *rho_elec_plus, *rho_elec_polym, *field_elec, *field_elec_sq;
extern double *phB, *wB;
// Block Polymer 
extern char seq_file[30];
extern int NF_N, K, *K_i, Nm, **Nm_i;
extern double b0, v0;
extern double *qInt, *wA, **qA, **qcA, **QA, **QcA, *Q1;
extern double **Alpha_i, **Chi_i, *sigma_i;
extern double **PHA, **PHA_T, *phA, **WA;
//Numerical 
extern double dx, ds0, integ_cons;
extern int init_opt, HALF_BOOL;
extern char Win[30];
extern double xCmax, xCltot, xCneck;
extern int Ns0, MaxIT_PB;
extern double wopt_PB, Sm_PB;
extern int iter, iter_PB;
extern int **Ns_i, *Ns;
