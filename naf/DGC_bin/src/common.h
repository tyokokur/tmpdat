#include <fftw3.h>
/* OMP */
#ifndef _OMP
	#define _OMP
	#include <omp.h>
#endif

/*Stand Libs*/
#ifndef _STANDARDS
	#define _STANDARDS
	#include <stdlib.h>
	#include <stdio.h>
	#include <string.h>
	#define _USE_MATH_DEFINES
	#include <math.h>
	#include <time.h>
	#include <unistd.h>
	#include <fcntl.h>
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

/*IJS Macro*/
#ifndef _IJSMACRO
	#define _IJSMACRO
	#define _IJS(ijk,Y,s) (Y*(Nm_i[X][Y]+1) + s)*NxNyNz + ijk
	#define _Yijk Y*NxNyNz + ijk
#endif

// Physical
extern double eps_vacuum, kB, Length, Charge, avo_num;
// System 
extern double freeE, freeDiff;
extern int Temp;
extern double lx, ly, lz, lz_cutoff;
extern int Nx, Nx2, Nxc;
extern int Ny, Nz, Nz_cutoff;
extern int NxNyNz, NxNy, NxNyNz_cutoff;
extern int Nzg_2;
extern double *rx, *rx_sq, *ry, *ry_sq, *rz, *rz_sq;
// Solution
extern int Z_plus, Z_minus;
extern double eps_r_P, eps_r_S, Q2, freeEnergy_bulk;
extern double c0_plus_n, c0_plus, c0_minus, fugac_plus, fugac_minus; 
extern double eps_P, eps_S, L_Deby_S, L_Deby_P, mu_s, zs;
extern double *eps_prof, *pot_elec, *rho_elec_minus, *rho_elec_plus, *rho_elec_polym, *field_elec, *field_elec_sq;
extern double *phB, *wB, *eta;
// Block Polymer 
extern char seq_file[30];
extern int NF_N, K, *K_i, Nm, **Nm_i;
extern double b0, v0;
extern double *qInt, *wA, **qA, **qcA, **QA, **QcA, *int_Q1_Quen;
extern double **Alpha_i, **Chi_i, *sigma_i, *pSURF;
extern double **PHA, **PHA_T, *phA, **WA;
extern double *GAUSSPROP, *xGAUSSDIST, *yGAUSSDIST, *zGAUSSDIST;
//Numerical 
extern int and_NrMax, MAXMAX;
extern double dx, dy, dz, ds0, integ_cons;
extern int init_opt, HALF_BOOL;
extern int win_box, win_nx, win_ny, win_nz;
extern char Win[30];
extern double xCmax, xCltot, xCneck;
extern int Ns0, MaxIT_PB;
extern double wopt, wcmp, wopt_PB, Sm_PB, wand;
extern double and_err, inCompMax, A_r;
extern int iter, iter_PB, and_it;
extern int **Ns_i, *Ns;
extern double *ksq, *deltaFun;
extern double ***DAs, ***WAs, **WBs, **DBs;
extern int **IJK_PERP;
extern int N_PERP;
//Logs
extern char outlog_name[30], errlog_name[30];
extern int  STDOUT_REDIR;
//Files
extern char pattern[30], phname[30], phname_elec[30], Wname[30], itname[30];
//FFTW
extern double *fftw_in, *fftw_inc, *fftw_inc_xy, *fftw_ins, *fftw_inc_cutoff;
extern double *fftw_outs_r, *fftw_outc_r, *fftw_outc_r_xy, *fftw_outc_r_cutoff;
extern fftw_complex *fftw_out_c;
extern double *ksq, *ksq_xy, *ksq_cutoff;
extern double *GAUSSPROP, *zGAUSSPROP, *FG_in, *FG_out;
extern int OMP_NUM_THREADS;
