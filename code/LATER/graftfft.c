/*UPDATED 2023618*/
//Discrete sine transform (MDE)
//phA(s) quadrature
//fftw (cos/full) MDE, PBE

/**********************************************************************/
/************************* length unit b0 *****************************/

#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>
#include <string.h>
//#include "mpi.h"
#include</share/apps/software/fftw-3.3.8/include/fftw3.h> // Pitzer cluster fftw3.3.8 location

#define MaxIT 2e4           //Maximum iteration steps
int Nx = 150;	//grid size
int Nx2, Nxc; //used in fftw
//#define _USE_BULK

#define Pi 3.141592653589
#define Pi_2 Pi*2.0
#define Pi_4 Pi*4.0

void PhA_quad(int k, double delta, double **PHA, double **PHA_T, double **QA, double **QcA);
void MDE_FFT(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, fftw_complex *fftw_out);
void PB_FFT(double *phA, double **PHA, double **PHA_T, double *fftw_in, fftw_complex *fftw_out);
void MDE_SIN(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out);
void MDE_COS(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out);
void PB_COS(double *phA, double **PHA, double **PHA_T, double *fftw_in, double *fftw_out);

void initW(double *wA1,double *wA2,double *wA3,double *wA4,double *wA5,double *wB);
double freeE(double *phA, double **PHA, double **PHA_T, double *phB, double **WA, double *wB, double *eta);
double getConc(double *phA, double **PHA, double **PHA_T, double *phB, double **WA, double *wB);
void sovDif_CR(double *in,double **g, int K, double *W, int *Ns, int sign);
void thomas(int n, double **a, double *b, double *x);
void solve_eq_tridiag(int n,double *in,double *out,double **aij);
void solve_PB(double *phA, double **PHA, double **PHA_T);
void write_ph(double *phA, double **PHA_T, double **PHA, double *phB);
void write_sim(double *phA, double **PHA_T, double **PHA, double *phB);
void write_W(double **WA,double *wB,double *eta);
void write_elec(double *pot_elec,double *rho_elec_plus,double *rho_elec_minus,double *rho_elec_polym);
void write_qA(char *name, int N, int Ns, double **qA); 
double And_mix(double **WA, double *wB); 
void inv(int ndim, double **A, double **invA);
void solve(int ndim, double **a, double *b, double *x);
void mult(int ndim, double **A, double *B, double *x);

double *sigma_i; //grafting density of each L/M/H
int NF_N; //Number of chain types (= 3 for L/M/H)
int KMAX; //Max number of blocks in a chain type
int *K_i; //Number of blocks in a chain type
double chi;  // Flory-Huiggins chi //
double lx; // size along x-axis //
double c1; // averaged volume fraction of polymer //
double Vbox;
double Nm; //Number of Kuhn segments for longest chain type
double b0,v0,p0;  // Kuhn segment lengh, segment volume, stiffness parameter //
double Nm_02,Nm_03; 

double mu_s;  // chemical potential od solvent //
double zs;  // fugacity of solvent //

int Ns0;  // number of discretization of polymer chain //
int **Nm_i, **Ns_i;
double **Alpha_i, **Chi_i;
int KMAX; 
int block_N;
double ds0;  // difference step of discretization of polymer chain //

double R_p,V_p; // estimated radius of polymeric globule //
double *rx,*rx_sq; // position //
double dx;

double *Q1,Q2; // partition function of polymer and solvent, respectively //

double q0,q0_sq;  // q0=1 //
double A_r;  // random amplitude //
double wand,wopt,wcmp;  // coefficent of iteration //
int and_it;
double Sm1,Sm2; // error tolerance //


double wa0,wb0;
double freeEnergy_bulk;

double rsqt;
double integ_cons;

int iter;
/// electric ///
int Z_minus,Z_plus;
double **PHA,*phA, *phB, **WA, *wB, *eta;
double *pot_elec,*rho_elec;
double *rho_elec_plus,*rho_elec_minus,*rho_elec_polym;
double *field_elec,*field_elec_sq;
double *eps_prof;
double fugac_plus,fugac_minus;
double c0_plus,c0_minus;
double eps_0,eps_r_P,eps_r_S,eps_P,eps_S;

double delt_PS_v0;

double L_Deby_P,L_Deby_S,L_Bjer_S;
double multip_Deby,multip_Deby_y;

double Length,Temp,Charge,kbT;  // SI unit //


double rho_fix;
double a11,a12,a21,a22;
double b1,b2;
double det,det1,det2;
double A1,A2;

////// PB iteration //////
double wopt_PB,Sm_PB,err_PB;
int MaxIT_PB, iter_PB;
//////////////////////////


double tt_mask;  // interfacial width of tanh() //

//////////////////////
	double freeEnergy,freeOld;
	double freeW,freeU,freeS,free_elec,freeDiff;
	double free_elec_polym,free_elec_laplace,free_elec_ion;
	//double free_factorial;

/// cycle of chain number ///
int i_np;
int np_m, np_M;
/////////////////////////////

//AndersonMixing Histories//
int and_NrMax=1; // History record
double ***DAs, ***WAs;
double **DBs, **WBs;
double *Cs;
double **WAdiff, *wBdiff;

// FFTW //
int FFTW_HALF = 1;
fftw_plan p, ip, cfp, cip, sfp, sip;
double *fftw_in;
double *fftw_inc, *fftw_ins;
fftw_complex *fftw_out_c;
double *fftw_outs_r, *fftw_outc_r;

/*Memory Allocation Macros*/
#define INIT(size) (double *)calloc(size, sizeof(double))
#define INIT_2D(NAME, size) (double **)calloc(size, sizeof(double)); for (i=0; i<size; i++) NAME[i] = INIT(size)
#define CHECK_VAL(VAR) printf("%.6e ", VAR)
#define CHECK_VAL_N(VAR) printf(#VAR": %.6e\n", VAR)

/*Matrix Collapse Reference Macros*/
#define _IJS(i,j,s) (j*(Nm_i[k][j]+1) + s)*Nx + i

char phname[30],phname_elec[30],Win[30],Wname[30],itname[30];

int main(int argc, char **argv)
{

	const double *TEST;
	double START[10] = {1,2,3,4,5,6,7,8,9,10};
	memcpy(TEST, START, sizeof(START));
	printf("%f\n", TEST[0]);
	exit(1);	

	double **PHA_T;
	double *wA1,*wA2,*wA3,*wA4,*wA5;
	double *phA1,*phA2,*phA3,*phA4,*phA5;
	double ves_dx, ves_r, vcust;
	double e1,e2,e3,e4,e5,e6,e7,e8;
	double xCneck, xCltot, xCmax;
	int i,j,k,in,iseed=-3;
	
	int i1,i2;
	FILE *fp, *fp1;
	char seq_file[50], seq[30], pattern[20];
	int vopt;

	time_t ts;
	iseed=time(&ts); srand(iseed);

	fp1=fopen("graft1.txt","r");
	fscanf(fp1,"%d",&in);		//in=1: inputing configuration is given;
	fscanf(fp1, "%s", seq_file);

	fp = fopen(seq_file, "r");
	fscanf(fp, "%s\n", seq);
	NF_N = strlen(seq);
	sigma_i = INIT(NF_N);
	Q1 = INIT(NF_N);
	K_i = (int *)calloc(NF_N, sizeof(double));
	Nm_i = (int **)calloc(NF_N, sizeof(int));
	Ns_i = (int **)calloc(NF_N, sizeof(int));
	Alpha_i = (double **)calloc(NF_N, sizeof(double));
	Chi_i = (double **)calloc(NF_N, sizeof(double));

	for (i = 0; i < NF_N-1; i++) fscanf(fp1, "%lf,", &sigma_i[i]); 
	fscanf(fp1, "%lf\n", &sigma_i[NF_N-1]); 
	
	fscanf(fp1,"%d,%d",&np_m,&np_M);
	fscanf(fp1,"%lf,%lf,%lf",&chi,&lx); //scaff PROTEIN_i_CONFIG.txt
	fscanf(fp1,"%lf,%lf",&Nm,&p0); //Nm not read, taken from seq_file
	fscanf(fp1,"%d",&Ns0); 
	fscanf(fp1,"%lf,%lf",&eps_r_P,&eps_r_S); //scaff PROTEIN_i_CONFIG.txt
	fscanf(fp1,"%lf,%lf",&c0_plus,&c0_minus); 
	fscanf(fp1,"%d,%d",&Z_plus,&Z_minus); 
	fscanf(fp1,"%lf,%lf",&Length,&Temp);
	fscanf(fp1,"%s",Win);		//input file name for W field;
	fscanf(fp1,"%s",Wname);		//output file name for W field;
	fscanf(fp1,"%lf",&A_r);
	fscanf(fp1,"%lf,%lf",&wopt,&wcmp);
	fscanf(fp1,"%lf,%d,%d",&wand, &and_it, &and_NrMax);
	fscanf(fp1,"%lf,%lf",&Sm1,&Sm2);
	fscanf(fp1,"%lf,%lf,%d",&wopt_PB,&Sm_PB,&MaxIT_PB);
	fscanf(fp1,"%lf",&tt_mask);
	fscanf(fp1, "%lf, %lf, %lf", &xCmax, &xCltot, &xCneck);
	fscanf(fp1, "%lf, %d, %lf", &b0, &vopt, &vcust);
	fclose(fp1);

	if (vopt == -1) v0 = vcust;
	if (vopt == 0) v0 = 4.0/3.0*M_PI*pow(b0, 3);
	if (vopt == 1) v0 = 1.00*pow(b0, 3);
	if (vopt == 2) v0 = 0.50*pow(b0, 3);
	if (vopt == 3) v0 = 0.33*pow(b0, 3);
	if (vopt == 4) v0 = 0.15*pow(b0, 3);
	printf("b0: %.2f\n", b0);
	printf("vopt: %d; v0 = %.2f\n", vopt, v0);

	V_p=np_M*Nm*v0;
	R_p=pow(0.75/Pi*V_p,1.0/3);
	/////////////////////////////////////////
	int orig_Nm; 	

	k = 0;
	Nm = 0; //Initialize max chain length

for (i=0;i<NF_N;i++) printf("%.2f ", sigma_i[i]);
	if (strstr(seq, "L") != NULL){
		fscanf(fp, "%ld\n", &K_i[k]);
		Nm_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Ns_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Alpha_i[k] = INIT(K_i[k]);
		Chi_i[k] = INIT(K_i[k]);
		for (j=0; j<K_i[k]; j++){
			fscanf(fp, "[%*d, %d], %lf, %lf\n", &Nm_i[k][j], &Alpha_i[k][j], &Chi_i[k][j]);
			orig_Nm = Nm_i[k][j];
			Nm_i[k][j] = round(Nm_i[k][j]/b0);
		}
		if (Nm_i[k][K_i[k]-1] > Nm) Nm = Nm_i[k][K_i[k]-1]; //Check last at last value of i = K_i[j]-1

		printf("Sigma_L: %.2f\n", sigma_i[k]);
		printf("L: %d Kuhn N =  %d AA  / %.2f [nm/b0]\n \tsigma_L = %.2f\n", Nm_i[k][K_i[k]-1], orig_Nm, b0,  sigma_i[k]);
		printf("\tAlpha[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Alpha_i[k][j]);
		printf("%.4f", Alpha_i[k][j]);
		printf("]\n");
		printf("\tChi[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Chi_i[k][j]);
		printf("%.4f", Chi_i[k][j]);
		printf("]\n");
		k += 1;
	}
	
	
	if (strstr(seq, "M") != NULL){
		fscanf(fp, "%ld\n", &K_i[k]);
		Nm_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Ns_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Alpha_i[k] = INIT(K_i[k]);
		Chi_i[k] = INIT(K_i[k]);
		for (j=0; j<K_i[k]; j++){
			fscanf(fp, "[%*d, %d], %lf, %lf\n", &Nm_i[k][j], &Alpha_i[k][j], &Chi_i[k][j]);
			orig_Nm = Nm_i[k][j];
			Nm_i[k][j] = round(Nm_i[k][j]/b0);
		}
		if (Nm_i[k][K_i[k]-1] > Nm) Nm = Nm_i[k][K_i[k]-1]; //Check last at last value of i = K_i[j]-1

		printf("Sigma_M: %.2f\n", sigma_i[k]);
		printf("M: %d Kuhn N =  %d AA / %.2f [nm/b0]\n \tsigma_M = %.2f\n", Nm_i[k][K_i[k]-1], orig_Nm, b0,  sigma_i[k]);
		printf("\tAlpha[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Alpha_i[k][j]);
		printf("%.4f", Alpha_i[k][j]);
		printf("]\n");
		printf("\tChi[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Chi_i[k][j]);
		printf("%.4f", Chi_i[k][j]);
		printf("]\n");
		k += 1;
	}

	if (strstr(seq, "H") != NULL){
		fscanf(fp, "%ld\n", &K_i[k]);
		Nm_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Ns_i[k] = (int *)calloc(K_i[k], sizeof(int));
		Alpha_i[k] = INIT(K_i[k]);
		Chi_i[k] = INIT(K_i[k]);
		for (j=0; j<K_i[k]; j++){
			fscanf(fp, "[%*d, %d], %lf, %lf\n", &Nm_i[k][j], &Alpha_i[k][j], &Chi_i[k][j]);
			orig_Nm = Nm_i[k][j];
			Nm_i[k][j] = round(Nm_i[k][j]/b0);
		}
		if (Nm_i[k][K_i[k]-1] > Nm) Nm = Nm_i[k][K_i[k]-1]; //Check last at last value of i = K_i[j]-1

		printf("Sigma_H: %.2f\n", sigma_i[k]);
		printf("H: %d Kuhn N =  %d AA  / %.2f [nm/b0]\n \tsigma_H = %.2f\n", Nm_i[k][K_i[k]-1], orig_Nm, b0,  sigma_i[k]);

		printf("\tAlpha[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Alpha_i[k][j]);
		printf("%.4f", Alpha_i[k][j]);
		printf("]\n");
		printf("\tChi[%d]: [", K_i[k]);
		for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Chi_i[k][j]);
		printf("%.4f", Chi_i[k][j]);
		printf("]\n");
		k += 1;
	}
	fclose(fp);

	snprintf(pattern, sizeof(char)*20, "_a%03.0fc%03.0f_x%03.0fs%03.0f", Alpha_i[0][0]*100, c0_plus*1000, Chi_i[0][0]*100, sigma_i[0]*100);
	snprintf(phname, sizeof(char)*30, "ph%s.dat", pattern);
	snprintf(phname_elec, sizeof(char)*30, "el%s.dat", pattern);
	snprintf(Wname, sizeof(char)*30, "W%s.dat", pattern);
	snprintf(itname, sizeof(char)*30, "it%s.dat", pattern);

	printf("%s\n", pattern);

	Charge=1.6e-19;
	kbT=1.38e-23*Temp;

	c0_plus*=6.022e23*1e3*pow(Length,3.0);  // dimensionless cation concentration //
	c0_minus*=6.022e23*1e3*pow(Length,3.0);  // dimensionless anion concentration //
	eps_0=8.854e-12*kbT*Length/(Charge*Charge);  // dimensionless permittivity in vaccum //
	eps_P=eps_0*eps_r_P;  // dimensionless dielectric constant of polymer //
	eps_S=eps_0*eps_r_S;  // dimensionless dielectric constant of solvent //
	L_Bjer_S=1.0/(Pi_4*eps_0*eps_r_S);  // dimensionless Bjerrum length in solvent //
	L_Deby_P=sqrt(eps_0*eps_r_P/(Z_plus*Z_plus*c0_plus+Z_minus*Z_minus*c0_minus));  // dimensionless Debye length in polymer //
	L_Deby_S=sqrt(eps_0*eps_r_S/(Z_plus*Z_plus*c0_plus+Z_minus*Z_minus*c0_minus));  // dimensionless Debye length in solvent //


	//ds0=1.0*Nm*0.360/b0/Ns0; // s=1,2,3,...,Ns0 //
	ds0 = 1.0*Nm/Ns0;
	CHECK_VAL_N(ds0);
	
	for (k=0;k<NF_N;k++) {
		printf("%d : [", k);
		for (j=0;j<K_i[k]-1;j++) {
			Ns_i[k][j] = roundl(Nm_i[k][j] / ds0); 
			printf("%d ",Ns_i[k][j]);
		}
		Ns_i[k][j] = roundl(Nm_i[k][j] / ds0); 
		printf("%d]\n",Ns_i[k][j]);
	}
	
	///////////////////////
	
	///////////////// length unit Nm^(1/3)b /////////////////
	/*
	b0/=Nm_03;
	ds0/=Nm;
	c0_plus*=Nm;
	c0_minus*=Nm;
	eps_0*=Nm_03;
	L_Bjer_S/=Nm_03;
	L_Deby_P/=Nm_03;
	L_Deby_S/=Nm_03;
	*/
	/////////////////////////////////////////////////////////

	mu_s=-1.0;
	zs=exp(mu_s);   //  exp(beta*mu_s)  //

	delt_PS_v0=0.5*(eps_P-eps_S)*v0;

	//q0_sq=q0*q0;

	fugac_plus=c0_plus; //*exp(u2_plus);
	fugac_minus=c0_minus; //*exp(u2_minus);


	//////////////////////// size of subvolume ////////////////////////////
	/*
	if(c0_plus!=0)
	{
		lx=R_p+multip_Deby*L_Deby_S;
	}
	//printf("lx=%10.4e\n",lx);  exit(0);

	if(c0_plus==0)
	{
		lx=3.0*R_p;
	}
	*/


	//Vbox=Pi*pow(lx,2.0)*2.0*ly;
	Vbox=4.0/3*Pi*pow(lx,3.0);
	c1=V_p/Vbox;

	dx=17.0/150; //Arbitrary dx that gives reasonable results
	Nx = roundl(lx/dx);
	printf("Nx: %d\n", Nx);

	//integ_cons=2.0*Pi_2*dx*dy/v0;  /// double the semisphere ///
	integ_cons=dx/v0; 
	
	PHA = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) PHA[i] = calloc(Nx*K_i[i], sizeof(double));
	PHA_T = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) PHA_T[i] = calloc(Nx, sizeof(double));
	WA = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WA[i] = calloc(Nx*K_i[i], sizeof(double));

	wB = INIT(Nx); phA = INIT(Nx); phB = INIT(Nx); eta = INIT(Nx);
	rx = INIT(Nx); rx_sq = INIT(Nx);
	rho_elec = INIT(Nx); rho_elec_minus = INIT(Nx); rho_elec_plus = INIT(Nx); rho_elec_polym = INIT(Nx); pot_elec = INIT(Nx);
	field_elec = INIT(Nx); field_elec_sq = INIT(Nx); eps_prof = INIT(Nx);

	for(i=0;i<Nx;i++)
	{
		rx[i]=i*dx;
		rx_sq[i]=rx[i]*rx[i];
	}


	/*
	freeEnergy_bulk=0.0;
	for(i=0;i<Nx;i++)
	{
		freeEnergy_bulk+=integ_cons[i];
	}
	freeEnergy_bulk*=-(mu_s+1.0)*Pi_4*dx/v0;
	*/

	freeEnergy_bulk=0.0;  // bulk free energy of ions //
	//if(fc_P_0!=0)
	//{
		for(i=0;i<Nx;i++)  
		{
			freeEnergy_bulk-=c0_plus+c0_minus;
		}
		freeEnergy_bulk*=integ_cons*v0;
	//}



	////////// for one polymer //////
	V_p=Nm*v0;
	R_p=pow(0.75/Pi*V_p,1.0/3);

	R_p=pow(0.75/Pi*1*Nm*v0,1.0/3);  // actual radius for each np //
	////// find log(np!) //////
	//free_factorial = (np+0.5)*log(np)-np+0.5*log(2*M_PI); //Stirling's Approx. Takashi

	//**************************definition of surface field and confinement***********************

	/***************initialize the electric field ***************/
	rho_fix=Alpha_i[0][0]/v0;

	a11=2.0*eps_P*(sinh(R_p/L_Deby_P)-R_p/L_Deby_P*cosh(R_p/L_Deby_P));
	a12=eps_S*(R_p/L_Deby_S+1.0)*exp(-R_p/L_Deby_S);
	a21=2.0*sinh(R_p/L_Deby_P);
	a22=exp(-R_p/L_Deby_S);

	b1=0.0;
	b2=0.5*rho_fix*R_p/c0_plus;

	det=a11*a22-a12*a21;
	det1=b1*a22-b2*a12;
	det2=a11*b2-a21*b1;

	A1=det1/det;
	A2=det2/det;

	for(i=1;i<Nx;i++)
	{
		if(rx[i]<R_p)
		{
			pot_elec[i]=A1/rx[i]*(exp(-rx[i]/L_Deby_P)-exp(rx[i]/L_Deby_P))+0.5*rho_fix/c0_plus; //scaff spherical analytical
		}
		else
		{
			pot_elec[i]=A2/rx[i]*exp(-rx[i]/L_Deby_S);
		}
	}
	pot_elec[0]=pot_elec[1];

	printf("Analytical over\n");


	for(i=0;i<Nx;i++)
	{
		pot_elec[i]=-pot_elec[i]; // negative charge on polymer //
	}

	/***************Initialize wA, wB******************/
	if (in==-2){
		fp = fopen(Win, "r");
		for (i=0; i<Nx; i++) {
			fscanf(fp, "%*lf %lf ", &phA[i]); // lx, phA
			for (k=0; k<NF_N; k++) {
				fscanf(fp, "%*lf [ "); // PHA_T
				for (j=0; j<K_i[k]-1; j++) {
					fscanf(fp, "%lf ", &PHA[k][j*Nx + i]);
				}
				fscanf(fp, "%lf ] ", &PHA[k][j*Nx + i]); 
				fscanf(fp, "%*lf\n"); // phB
			}
		}
		fclose(fp);

		for (i=0; i<Nx; i++){
			pot_elec[i]=0.0;
			phB[i] = 1.0 - phA[i];
			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++){	
				WA[k][j*Nx + i] = Chi_i[k][j]*phB[i] - Alpha_i[k][j]*pot_elec[i];
			}
			wB[i] = 0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) wB[i] += Chi_i[k][j]*PHA[k][j*Nx + i];
			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++)	WA[k][j*Nx + i] *= (1.0+A_r*(rand()/RAND_MAX-0.5));
			wB[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
		
			field_elec[i]=0.0;
			field_elec_sq[i]=0.0;
		}
	}

	else if(in==-1) { //READ W_INIT
		fp=fopen(Win,"r"); 
		for(i=0;i<Nx;i++) {
			for (k=0; k<NF_N; k++) {
				fscanf(fp, "[ ");
				for (j=0; j<K_i[k]-1; j++) {fscanf(fp, "%lf ", &WA[k][j*Nx + i]); WA[k][j*Nx + i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));}
				fscanf(fp, "%lf ] ", &WA[k][j*Nx + i]); 
				WA[k][j*Nx + i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
			}
			fscanf(fp, "%lf %lf %lf\n", &wB[i], &eta[i], &pot_elec[i]);
			wB[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
			eta[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
			pot_elec[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
		}
		fclose(fp);
		///// get electric field /////
		i=0;i2=i+1;
		field_elec[i]=field_elec[i2];

		for(i=1;i<Nx-1;i++)
		{
			i1=i-1;
			i2=i+1;
			field_elec[i]=0.5*(pot_elec[i1]-pot_elec[i2])/dx;
		}

		i=Nx-1;
		field_elec[i]=field_elec[i2];


		for(i=0;i<Nx;i++)
		{
			field_elec_sq[i]=field_elec[i]*field_elec[i];
		}

	}
	else if (in == 0) { //Linear using xCmax
		if (xCmax == 0)  xCmax = 0.80; //default
		for (i=1; i<Nx; i++){ //Respect no pen at x = 0
			for (k=0; k<NF_N; k++) {
				PHA_T[k][i] = 0;
				PHA[k][0*Nx + i] = xCmax * Nm_i[k][0]/Nm;
				for (j=1; j<K_i[k]; j++){
					PHA[k][j*Nx + i] = xCmax*((Nm_i[k][j]-Nm_i[k][j-1])/Nm);
					PHA_T[k][i] += PHA[k][j*Nx + i];
				}
				PHA_T[k][i] += PHA[k][0*Nx + i];
			}
		}
	
		for (i=0; i<Nx; i++){
			phA[i] = 0;
			for (k=0; k<NF_N; k++) phA[i] += PHA_T[k][i]; 

    		phB[i]=1.0-phA[i]; 
			eps_prof[i]=phA[i]*eps_P+(1.0-phA[i])*eps_S;

			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++){	
				WA[k][j*Nx + i] = Chi_i[k][j]*phB[i] - Alpha_i[k][j]*pot_elec[i];
			}
			
			wB[i] = 0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) wB[i] += Chi_i[k][j]*PHA[k][j*Nx + i];			

			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++)	WA[k][j*Nx + i] *= (1.0+A_r*(rand()/RAND_MAX-0.5));
			wB[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
		
			pot_elec[i]=0.0;
			field_elec[i]=0.0;
			field_elec_sq[i]=0.0;
		}
		//write_ph(phA, PHA_T, PHA, phB); write_W(WA, wB, eta); write_elec(pot_elec, rho_elec, rho_elec_minus, rho_elec_polym); exit(0);
	}
	else { //xC init, in==no. peaks; not super useful for xC = 1
		if (xCmax == 0)  xCmax = 0.80; //defaults
		if (xCltot==0) xCltot = lx/2.0;
		if (xCneck==0) xCneck = xCltot/(2.0*in);
		double w;
		int c;

		w = (xCltot - (in-1)*xCneck)/in; 
		for (i=0; i<Nx; i++){
			for (k=0; k<NF_N; k++) {
				PHA_T[k][i] = 0;
				PHA[k][0*Nx + i] = 0; for (c=1; c<=in; c++) PHA[k][0*Nx + i] += xCmax*(Nm_i[k][0]/Nm) * exp(-1.0/2* pow((i*dx-c*1.25* (w+xCneck)/2), 2) / sqrt(w/xCltot));
				for (j=1; j<K_i[k]; j++){
					PHA[k][j*Nx + i] = 0;
					for (c=1; c<=in; c++) PHA[k][j*Nx + i] += xCmax*((Nm_i[k][j]-Nm_i[k][j-1])/Nm) * exp(-1.0/2* pow((i*dx-c*1.25* (w+xCneck)/2), 2) / sqrt(w/xCltot));
					PHA_T[k][i] += PHA[k][j*Nx + i];
				}
				PHA_T[k][i] += PHA[k][0*Nx + i];
			}
			phA[i] = 0;
			for (k=0; k<NF_N; k++) phA[i] += PHA_T[k][i]; 

			phB[i]=1.0-phA[i]; 
			eps_prof[i]=phA[i]*eps_P+(1.0-phA[i])*eps_S;

			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++){	
				WA[k][j*Nx + i] = Chi_i[k][j]*phB[i] - Alpha_i[k][j]*pot_elec[i];
			}
			
			wB[i] = 0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) wB[i] += Chi_i[k][j]*PHA[k][j*Nx + i];			

			for (k=0; k<NF_N; k++) for (j=0; j<K_i[k]; j++)	WA[k][j*Nx + i] *= (1.0+A_r*(rand()/RAND_MAX-0.5));
			wB[i]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
		
			pot_elec[i]=0.0;
			field_elec[i]=0.0;
			field_elec_sq[i]=0.0;
		}
		//printf("for out\n");
		//write_ph(phA, PHA_T, PHA, phB); write_W(WA, wB, eta); write_elec(pot_elec, rho_elec, rho_elec_minus, rho_elec_polym); exit(0);
	}
	printf("Init over\n");
	freeE(phA, PHA, PHA_T, phB, WA, wB, eta);

}// End main()

//********************Output configuration******************************

void write_ph(double *phA, double **PHA_T, double **PHA, double *phB)
{
	int i, j, k;

	FILE *fp=fopen(phname,"w");

	for(i=0;i<Nx;i++) 
	{
		fprintf(fp,"%10.3e ",rx[i]);
		fprintf(fp, "%10.4e ", phA[i]);
		for (k=0; k<NF_N; k++) {
			fprintf(fp, "%10.4e ", PHA_T[k][i]);	
			fprintf(fp, "[ ");
			for (j=0; j<K_i[k]-1; j++) fprintf(fp, "%10.4e ", PHA[k][j*Nx + i]);
			fprintf(fp, "%10.4e ] ", PHA[k][j*Nx + i]);
		}
		fprintf(fp, "%10.4e\n", phB[i]);
	}

	fclose(fp);
}
void write_sim(double *phA, double **PHA_T, double **PHA, double *phB)
{
	int i, j, k;
	char phnew[30] = "SIM.dat";
	
	FILE *fp=fopen(phnew,"a");
	for(i=0;i<Nx;i++) 
	{
		fprintf(fp, "%10.4e\n", phA[i]);
	}
	fclose(fp);
}

void write_elec(double *pot_elec,double *rho_elec_plus,double *rho_elec_minus,double *rho_elec_polym)
{
	int i;

	FILE *fp=fopen(phname_elec,"w");

	for(i=0;i<Nx;i++)
	{	
		fprintf(fp,"%10.3e %10.4e %10.4e %10.4e %10.4e\n",rx[i],-pot_elec[i],rho_elec_plus[i],rho_elec_minus[i],rho_elec_polym[i]);
	}

	fclose(fp);
}


void write_W(double **WA,double *wB,double *eta) 
{
	int i,j,k;

	FILE *fp=fopen(Wname,"w");

	for(i=0;i<Nx;i++) 
	{
		for (k=0; k<NF_N; k++) {
			fprintf(fp, "[ ");
			for (j=0; j<K_i[k]-1; j++) fprintf(fp, "%10.4e ", WA[k][j*Nx + i]);
			fprintf(fp, "%10.4e ] ", WA[k][j*Nx + i]);
		}
		fprintf(fp, "%10.5e %10.5e %10.5e\n", wB[i], eta[i], pot_elec[i]);
	}

	fclose(fp);
}

void write_qA(char *name, int N, int Ns, double **qA) {
	int i,j,s;
	FILE *fp=fopen(name,"w");
	for (s=0; s<=Ns; s++) for(i=0;i<N;i++) {
		fprintf(fp, "%10.5e\n", qA[i][s]);
	}
	fclose(fp);
	printf("wrote qA\n");
}

void initW(double *wA1,double *wA2,double *wA3,double *wA4,double *wA5,double *wB)
{
	int i;

	wA1[0] = 0.0; wA2[0] = 0.0; wA3[0] = 0.0; wA4[0] = 0.0; wA5[0] = 0.0; wB[0] = 1.0; 
	for(i=1;i<Nx;i++)
	{
		wA1[i]=chi*(1-c1)*(1.0+A_r*(rand()/RAND_MAX-0.5));
		wA2[i]=chi*(1-c1)*(1.0+A_r*(rand()/RAND_MAX-0.5));
		wA3[i]=chi*(1-c1)*(1.0+A_r*(rand()/RAND_MAX-0.5));
		wA4[i]=chi*(1-c1)*(1.0+A_r*(rand()/RAND_MAX-0.5));
		wA5[i]=chi*(1-c1)*(1.0+A_r*(rand()/RAND_MAX-0.5));
		wB[i]=chi*c1*(1.0+A_r*(rand()/RAND_MAX-0.5));
	}

}


//*************************************main loop****************************************

double freeE(double *phA, double **PHA, double **PHA_T, double *phB, double **WA, double *wB, double *eta)
{
	int i,j,k,n,maxIter;
	int is, stab = 0;
	int MAXMAX; //scaff
	int chi_step = 0;
	double chi_kj;
	double and_err;
//
	double *pot_old;
	pot_old = INIT(Nx);
//

	double inCompMax;
	double beta,psum,fpsum,wpref, free_elec_polym_step,eta1,eta2,eta3,eta4,eta5,freeU_step;

	double **WAnew;
	WAdiff = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WAdiff[i] = calloc(Nx*K_i[i], sizeof(double));
	WAnew = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WAnew[i] = calloc(Nx*K_i[i], sizeof(double));

	double *wBnew;

	wBdiff=(double *)malloc(sizeof(double)*Nx);
	wBnew=(double *)malloc(sizeof(double)*Nx);

	DAs = (double ***)calloc(and_NrMax+1, sizeof(double)); 
	for (n=0; n<and_NrMax+1; n++){ 
		DAs[n] = (double **)calloc(NF_N, sizeof(double));
		for (k=0; k<NF_N; k++)	DAs[n][k] = (double *)calloc(K_i[k]*Nx, sizeof(double));	
	}
	WAs = (double ***)calloc(and_NrMax+1, sizeof(double)); 
	for (n=0; n<and_NrMax+1; n++){
		WAs[n] = (double **)calloc(NF_N, sizeof(double));
		for (k=0; k<NF_N; k++) WAs[n][k] = (double *)calloc(K_i[k]*Nx, sizeof(double));
	}
	DBs = (double **)calloc(and_NrMax+1, sizeof(double));
	for (n=0; n<and_NrMax+1; n++) DBs[n] = (double *)calloc(Nx, sizeof(double));
	WBs = (double **)calloc(and_NrMax+1, sizeof(double));
	for (n=0; n<and_NrMax+1; n++) WBs[n] = (double *)calloc(Nx, sizeof(double));

	// FFTW //
	//global bool FFTW_HALF
	Nx2 = Nx * 2;
	Nxc = Nx + 1; //  = (2 * Nx / 2) + 1 ; // fft space (2 Nx), / 2 + 1 in fftw r2c alg
	if (FFTW_HALF){
		printf("--- FFTW Halfspace ---\n");
		fftw_inc    = INIT(Nx);
		fftw_outc_r = INIT(Nx);
		fftw_ins    = INIT(Nx-1); // 
		fftw_outs_r = INIT(Nx-1);
		cfp = fftw_plan_r2r_1d(Nx, fftw_inc, fftw_outc_r, FFTW_REDFT00, FFTW_MEASURE); 
		cip = fftw_plan_r2r_1d(Nx, fftw_outc_r, fftw_inc, FFTW_REDFT00, FFTW_MEASURE);
		sfp = fftw_plan_r2r_1d(Nx-2, fftw_ins, fftw_outs_r, FFTW_RODFT00, FFTW_MEASURE); 
		sip = fftw_plan_r2r_1d(Nx-2, fftw_outs_r, fftw_ins, FFTW_RODFT00, FFTW_MEASURE);
	}
	else {
		printf("--- FFTW Fullspace ---\n");
		fftw_in =  INIT(Nx2);
		fftw_out_c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nxc);
		p  = fftw_plan_dft_r2c_1d(Nx2, fftw_in, fftw_out_c, FFTW_ESTIMATE); 
		ip = fftw_plan_dft_c2r_1d(Nx2, fftw_out_c, fftw_in, FFTW_ESTIMATE);
	}
	//////////
	
	maxIter=MaxIT;
	beta=1.0;

    iter=0;		

	freeEnergy=0.0;

	printf("freeE in\n");
	do
	{
        iter=iter+1;
        getConc(phA, PHA, PHA_T, phB, WA, wB);
//	if (iter==1) write_sim(phA, PHA_T, PHA,phB);
        for(i=0;i<Nx;i++) eps_prof[i]=phA[i]*eps_P+(1-phA[i])*eps_S;

	solve_PB(phA, PHA, PHA_T); 
	//if (FFTW_COS) PB_COS(phA, PHA, PHA_T, fftw_in, fftw_out_r);
	//else PB_FFT(phA, PHA, PHA_T, fftw_in, fftw_out_c);

		freeW=0.0;
		freeU=0.0;
		freeS=0.0;
		free_elec_polym=0.0;
		free_elec_laplace=0.0;
		free_elec_ion=0.0;
		inCompMax=0.0;

		for(i=0;i<Nx;i++)
		{
			psum=1.0-phA[i]-phB[i];

			fpsum=fabs(psum);
			if(fpsum>inCompMax){inCompMax=fpsum; MAXMAX = i;}

				eta1 = 0; eta2 = 0; 
				for (k=0;k<NF_N;k++){	
					eta1 += K_i[k];
					for (j=0;j<K_i[k];j++){
						eta2 += (1-phA[i])*Chi_i[k][j] - pot_elec[i]*Alpha_i[k][j] + Chi_i[k][j]*PHA[k][j*Nx + i] - WA[k][j*Nx + i];
					}	
				}
				eta[i] = 1/(eta1 + 1) * (eta2 - wB[i]);
				
				freeU_step = 0;
				for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) freeU_step += Chi_i[k][j]*PHA[k][j*Nx+i];
				freeU += freeU_step * phB[i];
			

			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) freeW -= WA[k][j*Nx+i] * PHA[k][j*Nx+i];
			freeW -= wB[i]*phB[i];

			free_elec_polym_step = 0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) free_elec_polym_step += Alpha_i[k][j]*PHA[k][j*Nx+i];
			free_elec_polym -= free_elec_polym_step*pot_elec[i]/v0;

			free_elec_laplace-=0.5*eps_prof[i]*field_elec_sq[i];
			free_elec_ion-=(rho_elec_plus[i]+rho_elec_minus[i]);
		}
		
		//CHECK_VAL_N(eta4);
		//printf("ETA: "); CHECK_VAL(eta[Nx-1]); printf("\n");

		//CHECK_VAL_N(integ_cons);

		freeU*=integ_cons; //CHECK_VAL_N(freeU);
		freeW*=integ_cons; //CHECK_VAL_N(freeW);

		free_elec_polym*=integ_cons*v0; //CHECK_VAL_N(free_elec_polym); //scaff check this
		free_elec_laplace*=integ_cons*v0; //CHECK_VAL_N(free_elec_laplace); //scaff check this
		free_elec_ion*=integ_cons*v0; //CHECK_VAL_N(free_elec_ion);
		free_elec_ion-=freeEnergy_bulk; //CHECK_VAL_N(free_elec_ion);
		free_elec=free_elec_polym+free_elec_laplace+free_elec_ion; //CHECK_VAL_N(free_elec); 

		freeS = -zs*Q2;
		for (k=0;k<NF_N;k++) freeS -= sigma_i[k]*log(Q1[k]);  //CHECK_VAL_N(freeS);

		freeOld=freeEnergy;
		freeEnergy=freeU+freeW+freeS+free_elec+freeEnergy_bulk; //+free_factorial;
		freeDiff=fabs((freeEnergy-freeOld)/freeEnergy);

		and_err = And_mix(WA, wB);

		if(iter%10==0)printf(" %5d : %.8e, %.8e, %.8e, %.8e, %d\n", //scaff
			iter, freeEnergy, freeDiff, inCompMax, and_err, MAXMAX);

		if(iter%100==0||iter>=maxIter)
		{

		//get_shape_parameter(phA);

		FILE *fp=fopen("printout.dat","w");
		fprintf(fp,"%d %10.3e %10.5e %10.5e %10.5e %10.5e\n",
			iter,lx,freeEnergy,freeDiff,inCompMax,and_err);
		fclose(fp);

		fp=fopen(itname,"a");
		fprintf(fp,"%d %10.5e %10.5e %10.5e %10.5e\n",
			iter,freeEnergy,freeDiff,inCompMax, and_err);
		fclose(fp);
		}
		if(iter%100==0)
		//if(iter%20==0)
		{
			write_ph(phA, PHA_T, PHA,phB);
//			write_sim(phA, PHA_T, PHA,phB);
			write_W(WA,wB,eta);
			write_elec(pot_elec,rho_elec_plus,rho_elec_minus,rho_elec_polym);
			
		}

/////////////////

		//write_ph(phA, PHA_T, PHA, phB);
		//write_W(WA,wB,eta);
		//write_elec(pot_elec,rho_elec_plus,rho_elec_minus,rho_elec_polym);
		//if (iter == 110) exit(0);
//////////////////

		/*if ((iter_PB != MaxIT_PB || err_PB == 0.0) && inCompMax <= 1e-03) {
			stab += 1; //stability of prev. chi
			if (stab > 500 && chi_step != 10) {
				chi_step += 1; 
				printf("\t New chi_step: %d\n", chi_step);
				stab = 0;
			} //If stable for more than 100 iters
		}
		else {
			stab = 0;
		}*/
		}while(iter<maxIter&&(and_err>1e-3||inCompMax>Sm1||freeDiff>Sm2||iter<123));

		if (FFTW_HALF){
			fftw_destroy_plan(cfp); fftw_destroy_plan(cip);
			fftw_destroy_plan(sfp); fftw_destroy_plan(sip);
		}
		else{
			fftw_destroy_plan(p);
			fftw_destroy_plan(ip);
		}
		free(fftw_in); 
		if (FFTW_HALF == 1){
			free(fftw_inc);
			free(fftw_ins);
			free(fftw_outs_r);
			free(fftw_outc_r);
		}
		else fftw_free(fftw_out_c);

		write_ph(phA, PHA_T, PHA, phB);
		write_W(WA,wB,eta);
		write_elec(pot_elec,rho_elec_plus,rho_elec_minus,rho_elec_polym);

		//get_shape_parameter(phA);
		
		FILE *fp=fopen("printout.dat","w");
		fprintf(fp,"%d %10.3e %10.5e %10.5e %10.5e %10.5e\n",
			iter,lx,freeEnergy,freeDiff,inCompMax,and_err);
		fclose(fp);

		#if !defined(_USE_BULK)
		fp = fopen("(0)_ProgramOver.dat", "w");
		fclose(fp);
		#endif

		fp=fopen(itname,"a");
		fprintf(fp,"%d %10.5e %10.5e %10.5e %10.5e\n",
			iter,freeEnergy,freeDiff,inCompMax,and_err);
		fclose(fp);
	}



double getConc(double *phA, double **PHA, double **PHA_T, double *phB, double **WA, double *wB)
{
	int i,j,k,iz,s,s1;
	double fflA,fflB;

	double *qInt;
	double **QA, **QcA;
	double **qA, **qcA, *wA;
	int s0;


	//////////////////////// Polymer //////////////////////////
	QA = (double **)calloc(NF_N, sizeof(double));
	for (k=0; k<NF_N; k++) QA[k] = INIT(K_i[k] * (Ns_i[k][K_i[k]-1]+1) * Nx); // 
	QcA = (double **)calloc(NF_N, sizeof(double));
	for (k=0; k<NF_N; k++) QcA[k] = INIT(K_i[k] * (Ns_i[k][K_i[k]-1]+1) * Nx);

	qInt = INIT(Nx);
	qA = (double **)calloc(Nx, sizeof(double));
	qcA = (double **)calloc(Nx, sizeof(double));
	int *Ns;
	double sigma = 1.0e-03; // 1.0e=02; From Chao PEB.c
	
	for (k=0; k<NF_N; k++){ 
		wA = INIT(Nx*K_i[k]);
		Ns = calloc(Ns_i[k][K_i[k]-1], sizeof(int));
		for (j=0; j<K_i[k]; j++) {
			for (i=0; i<Nx; i++) wA[j*Nx + i] = WA[k][j*Nx+i];
			Ns[j] = Ns_i[k][j];
		}
		//Begin forwards prop
		for (i=0;i<Nx;i++){
			qInt[i] = 0.0; // DBC wall
			//qInt[i] = 1.0/sqrt(2*M_PI * sigma) * exp(-0.5 * (i*dx)*(i*dx) / sigma);  //Gauss approx for \delta(z-eps), NBC wall
			qA[i] = INIT(Ns[K_i[k]-1] + 1); 
		} 
		qInt[1] = 1.0/dx; // Discretized \delta(z-eps), DBC wall

		// sovDif_CR(qInt, qA, K_i[k], wA, Ns, 1); //Use previous CR algorithm		
		// write_qA("CRqA.dat", Nx, Ns[K_i[k]-1], qA);

		if (FFTW_HALF) MDE_SIN(qInt, qA, K_i[k], wA, Ns, 1, fftw_ins, fftw_outs_r); 	
		else MDE_FFT(qInt, qA, K_i[k], wA, Ns, 1, fftw_in, fftw_out_c); 	

		for (j=0; j<K_i[k]; j++){
			if (j==0) s0 = 0;
			else s0 = Ns_i[k][j-1];

			for (i=0; i<Nx; i++) for (s=s0; s<=Ns_i[k][j]; s++){
				QA[k][_IJS(i,j,s)] = qA[i][s];
			}
		}

		//Begin backwards prop
		for (i=0;i<Nx;i++){
			qInt[i] = 1.0; //Graft at epsilon = 1*dx, discretized Dirac delta distribution
			qcA[i] = INIT(Ns[K_i[k]-1] + 1); 
		} 
		qInt[0] = 0.0; // DBC wall

		// sovDif_CR(qInt, qcA, K_i[k], wA, Ns, -1); //Use previous CR algorithm		
		// write_qA("CRqcA.dat", Nx, Ns[K_i[k]-1], qcA);
		if (FFTW_HALF) MDE_SIN(qInt, qcA, K_i[k], wA, Ns, -1, fftw_ins, fftw_outs_r);
		else MDE_FFT(qInt, qcA, K_i[k], wA, Ns, -1, fftw_in, fftw_out_c);

		for (j=0; j<K_i[k]; j++){
			if (j==0) s0 = 0;
			else s0 = Ns_i[k][j-1];

			for (i=0; i<Nx; i++) for (s=s0; s<=Ns_i[k][j]; s++){
				QcA[k][_IJS(i,j,s)] = qcA[i][s];
			}
		}

		for (i=0;i<Nx;i++) free(qcA[i]);
		for (i=0;i<Nx;i++) free(qA[i]);
		free(wA); free(Ns);
	}

	for (k=0;k<NF_N;k++){

		/*Q1[k] = QcA[k][_IJS(0,0,0)]/v0; 
		fflA = sigma_i[k] * ds0 / Q1[k];
		PhA_quad(k, fflA, PHA, PHA_T, QA, QcA);*/

		Q1[k]=0.0;
		for(i=0;i<Nx;i++)
		{
			Q1[k]+=QA[k][_IJS(i, (K_i[k]-1) , (Ns_i[k][K_i[k]-1]) )]; //integral of q(r, N)
		}
		Q1[k]*=integ_cons; 
		
		fflA=sigma_i[k]*ds0/Q1[k];
		for(i=0;i<Nx;i++){
			for (j=0;j<K_i[k];j++){
				PHA[k][j*Nx + i] = 0;

				if (j==0) {
					PHA[k][0*Nx+i] += 0.50*QA[k][_IJS(i,0, 0 )]*QcA[k][_IJS(i,0, 0)];
					for (s=1;s<Ns_i[k][0];s++) PHA[k][0*Nx+i] += QA[k][_IJS(i,0,s)]*QcA[k][_IJS(i,0,s)];
					PHA[k][0*Nx+i] += 0.50*QA[k][_IJS( i, 0, Ns_i[k][0] )]*QcA[k][_IJS( i, 0, Ns_i[k][0] )];
				}
				else{
					PHA[k][j*Nx+i] += 0.50*QA[k][_IJS(i,j, Ns_i[k][j-1] )]*QcA[k][_IJS(i,j, Ns_i[k][j-1] )];
					for (s=Ns_i[k][j-1]+1;s<Ns_i[k][j];s++) PHA[k][j*Nx+i] += QA[k][_IJS(i,j,s)]*QcA[k][_IJS(i,j,s)];
					PHA[k][j*Nx+i] += 0.50*QA[k][_IJS( i, j, Ns_i[k][j] )]*QcA[k][_IJS( i, j, Ns_i[k][j] )];
				}
			}
		}
		for (i=0; i<Nx;i++){
			PHA_T[k][i] = 0;
			for (j=0;j<K_i[k];j++){
				PHA[k][j*Nx+i] *= fflA;
				PHA_T[k][i] += PHA[k][j*Nx+i];
			}
		}
	}
	for (i=0;i<Nx;i++) {
		phA[i] = 0;
		for(k=0;k<NF_N;k++) phA[i] += PHA_T[k][i];	
	}

//for (i=0; i<10; i++) CHECK_VAL(phA[i]); printf("\n");
//exit(1);

	//if (iter==10) {write_ph(phA, PHA_T, PHA, phB); exit(0);}
  	///////////////////////////////////////////////////////////


  	//////////////////////// Solvent //////////////////////////
  	Q2=0.0;
  	for(i=0;i<Nx;i++)
  	{
  		phB[i]=zs*exp(-wB[i]);
  		Q2+=exp(-wB[i]);
  	}
  	Q2*=integ_cons;

  	///////////////////////////////////////////////////////////

	free(qInt); 
	for(k=0;k<NF_N;k++){free(QA[k]); free(QcA[k]);}
	free(QA); free(QcA); free(qA); free(qcA);
}

void sovDif_CR(double *in,double **g, int K, double *W, int *Ns, int sign){
	int i,j,k,is;
	double *Aq,*Bq;
	double *func_in,*func_out,**Mij;

	Aq=(double *)malloc(sizeof(double)*Nx);
	Bq=(double *)malloc(sizeof(double)*Nx);
	func_in=(double *)malloc(sizeof(double)*Nx);
	func_out=(double *)malloc(sizeof(double)*Nx);

	Mij=(double **)malloc(sizeof(double)*Nx);
  	for(i=0;i<Nx;i++)
  	{
    	Mij[i]=(double *)malloc(sizeof(double)*3);
  	}


	for(i=0;i<Nx;i++)
	{
		Aq[i] = -ds0*b0*b0; //Cartesian 1D
	}

    if(sign==1)  // forward //
    {

    	is=0;
		for(i=0;i<Nx;i++)
    	{
        	g[i][is]=in[i];  // initial value //
    	}
		
    	for(is=1;is<=Ns[K-1];is++)
    	{
			for (j=0; j<K; j++) if (is <= Ns[j]){ //First less than is belonging block
				for (i=0; i<Nx; i++) Bq[i] = 2*ds0*b0*b0 + 6*ds0*dx*dx*W[j*Nx+i]; 
				break;
			}
        	for(i=0;i<Nx;i++)
        	{
        		if(i==0)
            	{
                	Mij[i][0]= 0.0;
	                Mij[i][1]= 1.0; //Dirichlet
    	            Mij[i][2]= 0.0;
	            }
    	        else if(i==Nx-1)
        	    {
            	    Mij[i][0]= 0.0;
                	Mij[i][1]= 1.0; //Dirichlet
                	Mij[i][2]= 0.0;
	            }
    	        else
        	    {
            	    Mij[i][0] = Aq[i];
                	Mij[i][1] = 12*dx*dx + Bq[i];
	                Mij[i][2] = Aq[i];
        	    }
        	}

	        for(i=0;i<Nx;i++)
    	    {
        	    if(i==0)
            	{
            		func_in[i]=0.0; //No penetration (via Neumann)              									
	            }
    	        else if(i==Nx-1)
        	    {
            	    func_in[i]=0.0; //Infinite 
            	}
	            else
    	        {
					func_in[i] = -Aq[i]*g[i-1][is-1]  //Cartesian 1D
									+ (12*dx*dx - Bq[i]) *g[i][is-1]
									-Aq[i]*g[i+1][is-1];			
            	}
        	}

	    	thomas(Nx,Mij,func_in,func_out);
			

    		for(i=0;i<Nx;i++)
    		{
    			g[i][is]=func_out[i];
	    	}

    	}
    }
    if(sign==-1)  // backward //
    {

    	is=Ns[K-1]; //scaffffffff
		for(i=0;i<Nx;i++)
    	{ 
        	g[i][is]=in[i];  // initial value //
    	}
		
    	for(is=Ns[K-1]-1;is>=0;is--)
    	{
			for (j=0; j<K; j++) if (is <= Ns[j]){ //First less than is belonging block
				for (i=0; i<Nx; i++) Bq[i] = 2*ds0*b0*b0 + 6*ds0*dx*dx*W[j*Nx+i]; 
				break;
			}
        	for(i=0;i<Nx;i++)
        	{
        		if(i==0)
            	{
                	Mij[i][0]= 0.0;
	                Mij[i][1]= 1.0; //Dirichlet
    	            Mij[i][2]= 0.0;
	            }
    	        else if(i==Nx-1)
        	    {
            	    Mij[i][0]= 0.0;
                	Mij[i][1]= 1.0; //Dirichlet
                	Mij[i][2]= 0.0;
	            }
    	        else
        	    {
            	    Mij[i][0] = Aq[i];
                	Mij[i][1] = 12*dx*dx + Bq[i];
	                Mij[i][2] = Aq[i];
        	    }
        	}

	        for(i=0;i<Nx;i++)
    	    {
        	    if(i==0)
            	{
            		func_in[i]=0.0; //No penetration                									
	            }
    	        else if(i==Nx-1)
        	    {
            	    func_in[i]=0.0; //Infinite 
            	}
	            else
    	        {
					func_in[i] = -Aq[i]*g[i-1][is+1]  //Cartesian 1D
									+ (12*dx*dx - Bq[i]) *g[i][is+1]
									-Aq[i]*g[i+1][is+1];			
            	}
        	}

	    	thomas(Nx,Mij,func_in,func_out);

    		for(i=0;i<Nx;i++)
    		{
    			g[i][is]=func_out[i];
	    	}

	    }


    }

	free(Aq);
	free(Bq);

	free(func_in);
	free(func_out);

	for(i=0;i<Nx;i++)
	{
		free(Mij[i]);
	}
	free(Mij);

}

void solve_PB(double *phA, double **PHA, double **PHA_T)
{
    int i, j, k; 
    int i1,i2;

    double err_max,err_check;

    double *pot_elec_old;

    double *V_in,*V_out,**Aij;

	pot_elec_old = INIT(Nx);
	V_in = INIT(Nx); V_out = INIT(Nx);

    Aij = (double **)malloc(sizeof(double)*Nx);
    for(i=0;i<Nx;i++)
    {
    	Aij[i] = INIT(3);
    }


    iter_PB=0;

	do  // iteration for solving PB eq. //
	{
        iter_PB=iter_PB+1;

        for(i=0;i<Nx;i++)
        {
        	pot_elec_old[i]=pot_elec[i];
        }

		for(i=0;i<Nx;i++)
		{
			rho_elec_plus[i] = Z_plus*fugac_plus*exp(-Z_plus*pot_elec_old[i]); //SCAFF NEW
			rho_elec_minus[i] = Z_minus*fugac_minus*exp(Z_minus*pot_elec_old[i]); //SCAFF NEW
			rho_elec_polym[i] = 0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) rho_elec_polym[i] += Alpha_i[k][j]*PHA[k][j*Nx+i];
			rho_elec_polym[i] *= 1/v0;

		}

		Aij[0][1] = 1.0; Aij[0][2] = -1.0; //Neumann
    	for(i=1;i<Nx-1;i++)
    	{
    		//Aij[i][0] = 1 - dx/2; //SCAFF bad approx
			//Aij[i][1] = -2.0; //SCAFF bad approx
			//Aij[i][2] = 1 + dx/2; //SCAFF bad approx
			
			Aij[i][0] = 0.25*(eps_prof[i+1]-eps_prof[i-1]) - eps_prof[i]; //SCAFF USE FOR VALID W_INIT
   			Aij[i][1] = 2.0*eps_prof[i] + dx*dx* (Z_minus*rho_elec_minus[i] + Z_plus*rho_elec_plus[i]);
    		Aij[i][2] = -0.25*(eps_prof[i+1]-eps_prof[i-1]) - eps_prof[i];
    	}
   		Aij[Nx-1][1] = 1.0; //Dirichlet

   		V_in[0] = 0.0; //Symmetry
    	for(i=1;i<Nx-1;i++)
    	{
    		V_in[i] = dx*dx * (rho_elec_plus[i]*(1+Z_plus*pot_elec_old[i]) - rho_elec_minus[i]*(1-Z_minus*pot_elec_old[i])+rho_elec_polym[i]); //SCAFF USE FOR VALID w_INIT
    	}
    	V_in[Nx-1] = 0.0; //Infinite

		thomas(Nx,Aij,V_in,V_out);

		for(i=0;i<Nx;i++)
		{
    		pot_elec[i]=wopt_PB*V_out[i]+(1.0-wopt_PB)*pot_elec_old[i];  // simople linear superposit //
   		}
		
		err_max=0.0;
		for(i=0;i<Nx;i++)
		{
			err_check= (pot_elec[i]-pot_elec_old[i])*(pot_elec[i]-pot_elec_old[i]);
			if(err_check>err_max)err_max=err_check;
		}
		err_max=sqrt(err_max);

	}while(err_max>Sm_PB&&iter_PB<MaxIT_PB);
	
	//printf("iter=%d,err_max_PB=%10.5e,iter_PB=%d\n",iter,err_max,iter_PB);  // check //

	err_PB=err_max;

	///// get electric field /////
	for(i=0;i<Nx;i++)
	{
		if((i>0)&&(i<Nx-1))
		{
			i1=i-1;
			i2=i+1;
			field_elec[i]=0.5*(pot_elec[i1]-pot_elec[i2])/dx;
		}
	}

	i=0;i2=i+1;
	field_elec[i]=field_elec[i2];

	i=Nx-1;i2=i-1;
	field_elec[i]=field_elec[i2];


	for(i=0;i<Nx;i++)
	{
		field_elec_sq[i]=field_elec[i]*field_elec[i];
	}


	free(pot_elec_old);


	free(V_out);
    free(V_in);
    for(i=0;i<Nx;i++)
	{
		free(Aij[i]);
	}
	free(Aij);

}


void thomas(int n, double **a, double *b, double *x) ///////Takashi
{
    /*Thomas Algorithm for tridiagonal systems -- Hoffman 1.8.3*/
    //Input: 
        //n: number of rows in a
        //a: reduced nx3 matrix (changed in-place)
        //b: homogeneity vector from ax = b
    //Output:
        //a: nx3 simplified through Gaussain elim (changed in-place)
        //x: nx1 solution vector to ax = b (changed in-place)
    int i, f = n-1;
    double em; 

    //Forw1rd Elimination
    for (i = 1; i < n; i++)
    {
        em = a[i][0]/a[i-1][1];
        a[i][0] = em;
        a[i][1] = a[i][1] - em*a[i-1][2];
        b[i] = b[i] - a[i][0]*b[i-1];
    }

    //Back Substitution
    x[f] = b[f]/a[f][1];
    for (i = f-1; i >= 0; i--)
    {
        x[i] = (b[i] - a[i][2]*x[i+1])/a[i][1];
    }
    return;
}

void solve_eq_tridiag(int n,double *in,double *out,double **aij)   //  not considering pivoting //
{
    int i, j, k;
    double multip;


    multip=1.0/aij[0][1];
    aij[0][2]*=multip;
    in[0]*=multip;

    //aij[0][1]=1.0;


   	for(i=1;i<n-1;i++)
   	//for(i=1;i<n;i++)
    {
        multip=1.0/aij[i][0];
        for(j=1;j<3;j++)
        {
            aij[i][j]*=multip;
        }
        in[i]*=multip;

        //aij[i][0]=1.0;
    }

    for(i=1;i<n-1;i++)
   	//for(i=1;i<n;i++)
    {
        multip=1.0/(aij[i][1]-aij[i-1][2]);

        aij[i][2]*=multip;
        in[i]=(in[i]-in[i-1])*multip;

        //aij[i][1]=1.0;
    }

    for(i=n-2;i>=0;i--)
    {
        in[i]-=aij[i][2]*in[i+1];
    }


    for(i=0;i<n;i++)
    {
        out[i]=in[i];
    }

}

double And_mix(double **WA, double *wB){ 
//Anderson Mixing
	//Input globals: PHA, Chi_i, Alpha_i, pot_elec, eta 
	//Output updated WA and wB
	//Updated globals DAs WAs
	double **WAdiff, **WAnew, *wBdiff, *wBnew, **DAnew, *DBnew;
	double **u, **u_temp, *v;
	double psum, lambda;
	double amax, bmax;
	double and_err;

	int i,j,k, n, m;
	int and_Nr, end;
	WAdiff = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WAdiff[i] = calloc(Nx*K_i[i], sizeof(double));
	WAnew = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) WAnew[i] = calloc(Nx*K_i[i], sizeof(double));
	wBdiff=(double *)malloc(sizeof(double)*Nx);
	wBnew=(double *)malloc(sizeof(double)*Nx);
	DAnew = (double **)calloc(NF_N, sizeof(double)); for (i=0; i<NF_N; i++) DAnew[i] = calloc(Nx*K_i[i], sizeof(double));
	DBnew = (double *)calloc(Nx, sizeof(double)); 

	for (n=1; n<=and_NrMax; n++){ 
		for (k=0; k<NF_N; k++){ 
			for (j=0;j<K_i[k];j++){
				for (i=0;i<Nx;i++){ 
					DAs[n-1][k][j*Nx+i] = DAs[n][k][j*Nx+i];
					DBs[n-1][i] = DBs[n][i];
					WAs[n-1][k][j*Nx+i] = WAs[n][k][j*Nx+i];
					WBs[n-1][i] = WBs[n][i];
				}
			}
		}	
	}

	//Calculate and Update
	end = and_NrMax;

	amax = 0.0; 
	bmax = 0.0;
	for (i=0;i<Nx;i++){
		for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) {
			WAnew[k][j*Nx+i]    = Chi_i[k][j]*phB[i] - Alpha_i[k][j]*pot_elec[i] - eta[i];
			WAdiff[k][j*Nx+i]   = WAnew[k][j*Nx+i] - WA[k][j*Nx+i];
			DAs[end][k][j*Nx+i] = WAdiff[k][j*Nx+i]; //Update last val
			WAs[end][k][j*Nx+i] = WA[k][j*Nx+i]; //Update last val
			amax += WAdiff[k][j*Nx+i]*WAdiff[k][j*Nx+i] / (WA[k][j*Nx+i]*WA[k][j*Nx+i]);
			}
		wBnew[i] = 0;
		for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) wBnew[i] += Chi_i[k][j]*PHA[k][j*Nx + i];
		wBnew[i] -= eta[i];
		wBdiff[i] = wBnew[i]-wB[i];		
		DBs[end][i] = wBdiff[i]; //Update last val 
		WBs[end][i] = wB[i]; //Update last val
		bmax += wBdiff[i]*wBdiff[i] / (wB[i]*wB[i]);
	}	
	and_err = amax+bmax;
	and_err = pow(and_err, 0.50);
	
	if (iter%and_it!=0 || and_NrMax == 0){// || and_err>1e-02 || bmax>1e-03 || amax>1e-3){ 
		for (i=0; i<Nx; i++){
			psum = 1.0-phA[i]-phB[i];
			for (k=0;k<NF_N;k++) for (j=0; j<K_i[k]; j++) WA[k][j*Nx+i] += wopt*(WAdiff[k][j*Nx+i]- wcmp*psum);
			wB[i] += wopt*(wBdiff[i]-wcmp*psum);
		}
		for (k=0;k<NF_N;k++){free(WAdiff[k]); free(WAnew[k]); free(DAnew[k]);}
		free(WAdiff); free(WAnew); free(DAnew);
		free(wBdiff); free(wBnew); free(DBnew);
		return and_err;
	}	
	//if (iter%10==0) printf("err: %.5e, amax: %.3e, bmax: %.3e\n", and_err, amax, bmax);
	and_Nr = fmin(iter-1, and_NrMax); 
	Cs = (double *)calloc(and_Nr, sizeof(double)); 
	//Store histories
	u = INIT_2D(u, and_Nr);
	u_temp = INIT_2D(u_temp, and_Nr);
	v = INIT(and_Nr);
	for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) for (i=0;i<Nx;i++){
		for (m=1;m<=and_Nr;m++) {
			for (n=1;n<=and_Nr;n++){
				u[m-1][n-1] += (DAs[end][k][j*Nx+i] - DAs[end-m][k][j*Nx+i]) * (DAs[end][k][j*Nx+i] - DAs[end-n][k][j*Nx+i]);
			}	
			v[m-1] += (DAs[end][k][j*Nx+i] - DAs[end-m][k][j*Nx+i])*DAs[end][k][j*Nx+i];
		}
	}	
	for (i=0;i<Nx;i++) for (m=1;m<=and_Nr;m++){
 		for (n=1;n<=and_Nr;n++) u[m-1][n-1] += (DBs[end][i]-DBs[end-m][i])*(DBs[end][i]-DBs[end-n][i]);
		v[m-1] += (DBs[end][i]-DBs[end-m][i])*DBs[end][i];
	}

	//printf("Umn    ");
	//for (m=0;m<and_Nr;m++) for (n=0;n<and_Nr;n++) printf("%.6e ", u[m][n]);
	//printf("\n");
	//printf("Vm    ");
	//for (m=0;m<and_Nr;m++) printf("%.6e ", v[m]);
	//printf("\n");
	// Compute Coeffs	
	inv(and_Nr, u, u_temp); //u ^-1	
	for (m=0;m<and_Nr;m++) for (n=0;n<and_Nr;n++) u[n][m] = u_temp[m][n]; //uT^-1
	mult(and_Nr, u, v, Cs); //CN = uT^-1 * v
printf("CAS: ");
for (n=0;n<and_Nr;n++) printf("%.6e ", Cs[n]);
printf("\n");
	//Half step
	for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) for (i=0;i<Nx;i++){
		WAnew[k][j*Nx+i] = WAs[end][k][j*Nx+i]; //k+1/2
		DAnew[k][j*Nx+i] = DAs[end][k][j*Nx+i];
		for (n=1;n<=and_Nr;n++) {
			WAnew[k][j*Nx+i] += wand*Cs[end-n]*(WAs[end-n][k][j*Nx+i]-WAs[end][k][j*Nx+i]);
			DAnew[k][j*Nx+i] += wand*Cs[end-n]*(DAs[end-n][k][j*Nx+i]-DAs[end][k][j*Nx+i]);
		}
	}
	for (i=0;i<Nx;i++){
		wBnew[i] = WBs[end][i]; //k+1/2
		DBnew[i] = DBs[end][i];
		for (n=1; n<=and_Nr; n++){
			wBnew[i] += wand*Cs[end-n]*(WBs[end-n][i] - WBs[end][i]);
			DBnew[i] += wand*Cs[end-n]*(DBs[end-n][i] - DBs[end][i]);
		}
	}	
	//Update fields
	lambda =1.0; //1.0-pow(0.9, iter/200);
	for (i=0;i<Nx;i++){
		for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++){
			WA[k][j*Nx+i] = WAnew[k][j*Nx+i] + lambda * DAnew[k][j*Nx+i];
		}
		wB[i] = wBnew[i] + lambda * DBnew[i];
	}

	//Free
	for (k=0;k<NF_N;k++){free(WAdiff[k]); free(WAnew[k]); free(DAnew[k]);}
	free(WAdiff); free(WAnew); free(DAnew);
	free(wBdiff); free(wBnew); free(DBnew);
	for (n=0;n<and_Nr;n++) {free(u[n]); free(u_temp[n]);}
	free(v);
	return and_err;
}


void inv(int ndim, double **A, double **invA){
//Inverse using Doolittle LU Factorization
//Adapted from Hoffman
//Output inverse of A into invA
	int i,j,k;
	double em;
	double **l, **u;
	double *b, *x, *col; 
	b = INIT(ndim); x = INIT(ndim); col = INIT(ndim);
	u = INIT_2D(u, ndim);
	l = INIT_2D(l, ndim);
	//LU Decomp
	for (k=0;k<ndim-1;k++){
		for (i=k+1;i<ndim;i++){
			em = A[i][k] / A[k][k]; 	
			A[i][k] = em;
			for (j=k+1;j<ndim;j++){
				A[i][j] = A[i][j] - em*A[k][j];
			}
		}	
	}	
	//A-1 = L-1 * U-1
	for (i=0;i<ndim;i++){
		l[i][i] = 1.0;
		for (j=i-1;j>=0;j--) l[i][j] = A[i][j];
		for (j=i;j<ndim;j++) u[i][j] = A[i][j];
	}
	for (i=0;i<ndim;i++){
		for (j=0;j<ndim;j++){ 
			if (j==i) b[j] = 1.0;
			else      b[j] = 0.0;
		}
		solve(ndim, l, b, x);	
		solve(ndim, u, x, col);
		for (j=0;j<ndim;j++) invA[i][j] = col[j];
	}
	return;
}
void solve(int ndim, double **a, double *b, double *x){
	double *bp; bp = INIT(ndim);
	int i,j,k;
	int n = ndim-1;
	bp[0]= b[0];
	for (i=1;i<ndim;i++){
		bp[i] = b[i];
		for (j=0;j<i;j++){
			bp[i] = bp[i] - a[i][j]*bp[j];
		}
	}	
	x[n] = bp[n]/a[n][n];
	for (i=n-1;i>=0;i--){
		x[i] = bp[i];
		for (j=n;j>i;j--){
			x[i] = x[i]-a[i][j]*x[j];
		}
		x[i] = x[i]/a[i][i];
	}
}
void mult(int ndim, double **A, double *B, double *x){
//x = Ab: Square matrix by vector
	int i, j;
	for (i=0;i<ndim;i++){
		x[i] = 0.0;
		for (j=0;j<ndim;j++) x[i] += A[i][j]*B[j];
	}
	return;
}

void MDE_FFT(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, fftw_complex *fftw_out){
	double *rpref, *ksq;
	double ds, _mirror;
	int s, is, i, j, Nsf;
	int start;

	// FFTW requires PBC; Calculate with z mirror and return only [0, Nz)
	ksq = INIT(Nx2);
	for (i=0; i<Nx2; i++)   ksq[i] = (i * 2.0*M_PI/(Nx2*dx)) * (i * 2.0*M_PI/(Nx2*dx)); // Fourier transform change of vars
	rpref = INIT(Nx);

	Nsf = Ns[K-1]; // Last monomer
	ds = 1.0; // Discrete Gaussian chain

	if (sign == 1)  start = 0;   // Forwards
	if (sign == -1) start = Nsf; // Backwards

	for (i=0; i<Nx; i++) g[i][start] = init[i];

	for (is = 1; is <= Nsf; is++){
		s = start + sign*is; // Contour length s = [1, Nsf]
		for (j=0; j<K; j++) if (s <= Ns[j]){ // Piecewise W field based on is
			for (i=0; i<Nx; i++) rpref[i] = exp(-0.5 * ds * W[j*Nx + i]); 
			break; // Only do this once
		}
		for (i=0; i<Nx; i++) {
			_mirror = rpref[i] * g[i][s - sign];
			fftw_in[i]       = _mirror; // Real step 1
			fftw_in[Nx2-1-i] = _mirror;
		}
		fftw_execute(p); // FFT forward (in --> out)
		for (i=0; i<Nxc; i++) for(j=0;j<2;j++) fftw_out[i][j] *= exp( -ds * b0*b0/6.0 * ksq[i]); // Im step 2
		fftw_execute(ip); // FFT back   (out --> in)
		// Return only first halfspace
		for (i=0; i<Nx; i++) g[i][s] = fftw_in[i] * rpref[i] / Nx2; // Real step 3
	}
	free(rpref); free(ksq);
}

void MDE_SIN(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out){
	double *rpref, *ksq;
	double ds, _mirror;
	int s, is, i, j, Nsf;
	int start;

	// FFTW requires PBC; Calculate with z mirror and return only [0, Nz)
	ksq = INIT(Nx);
	for (i=0; i<Nx; i++)   ksq[i] = (i * 2.0*M_PI/(Nx2*dx)) * (i * 2.0*M_PI/(Nx2*dx)); // Fourier transform change of vars
	rpref = INIT(Nx);

	Nsf = Ns[K-1]; // Last monomer
	ds = 1.0; // Discrete Gaussian chain

	if (sign == 1)  start = 0;   // Forwards
	if (sign == -1) start = Nsf; // Backwards

	for (i=0; i<Nx; i++) g[i][start] = init[i];

	for (is = 1; is <= Nsf; is++){
		s = start + sign*is; // Contour length s = [1, Nsf]
		for (j=0; j<K; j++) if (s <= Ns[j]){ // Piecewise W field based on is
			for (i=0; i<Nx; i++) rpref[i] = exp(-0.5 * ds * W[j*Nx + i]); 
			break; // Only do this once
		}
		for (i=1; i<Nx; i++) fftw_in[i-1] = rpref[i] * g[i][s - sign];
		fftw_execute(sfp); // FFT forward (in --> out)
		for (i=1; i<Nx; i++) fftw_out[i-1] *= exp( -ds * b0*b0/6.0 * ksq[i]); // Im step 2
		fftw_execute(sip); // FFT back   (out --> in)
		// Return only first halfspace
		g[0][s] = 0.0;
		for (i=1; i<Nx; i++) g[i][s] = fftw_in[i-1] * rpref[i] / (Nx2); // Real step 3
}
	free(rpref); free(ksq);
}

void MDE_COS(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out){
	double *rpref, *ksq;
	double ds, _mirror;
	int s, is, i, j, Nsf;
	int start;

	// FFTW requires PBC; Calculate with z mirror and return only [0, Nz)
	ksq = INIT(Nx);
	for (i=0; i<Nx; i++)   ksq[i] = (i * 2.0*M_PI/(Nx2*dx)) * (i * 2.0*M_PI/(Nx2*dx)); // Fourier transform change of vars
	rpref = INIT(Nx);

	Nsf = Ns[K-1]; // Last monomer
	ds = 1.0; // Discrete Gaussian chain

	if (sign == 1)  start = 0;   // Forwards
	if (sign == -1) start = Nsf; // Backwards

	for (i=0; i<Nx; i++) g[i][start] = init[i];

	for (is = 1; is <= Nsf; is++){
		s = start + sign*is; // Contour length s = [1, Nsf]
		for (j=0; j<K; j++) if (s <= Ns[j]){ // Piecewise W field based on is
			for (i=0; i<Nx; i++) rpref[i] = exp(-0.5 * ds * W[j*Nx + i]); 
			break; // Only do this once
		}
		for (i=0; i<Nx; i++) fftw_in[i] = rpref[i] * g[i][s - sign];
		fftw_execute(cfp); // FFT forward (in --> out)
		for (i=0; i<Nx; i++) fftw_out[i] *= exp( -ds * b0*b0/6.0 * ksq[i]); // Im step 2
		fftw_execute(cip); // FFT back   (out --> in)
		// Return only first halfspace
		for (i=0; i<Nx; i++) g[i][s] = fftw_in[i] * rpref[i] / (Nx2-2); // Real step 3
	}
	free(rpref); free(ksq);
}

void PB_FFT(double *phA, double **PHA, double **PHA_T, double *fftw_in, fftw_complex *fftw_out){
	double *ksq, *pot_elec_old;
	double _mirror, c, err_check, err_max;
	int k, i1, i2, i, j;

	pot_elec_old = INIT(Nx);
	c = 1.0 / (L_Deby_S * L_Deby_S);
	// FFTW requires PBC; Calculate with z mirror and return only [0, Nz)
	ksq = INIT(Nx2);
	for (i=0; i<Nx; i++)     ksq[i] = (i * 2.0*M_PI/(Nx2*dx)) * (i * 2.0*M_PI/(Nx2*dx)); // Fourier transform change of vars
	for (i=Nx; i<Nx2; i++) ksq[i] = ((Nx2-i) * 2.0*M_PI/(Nx2*dx)) * ((Nx2-i) * 2.0*M_PI/(Nx2*dx)); // Fourier transform change of vars

	// Simple mixing it solver
	iter_PB = 0;
	do { 
		iter_PB += 1;

		for(i=0;i<Nx;i++){
			pot_elec_old[i]   = pot_elec[i];
			rho_elec_plus[i]  = Z_plus *fugac_plus *exp(-Z_plus*pot_elec[i]); 
			rho_elec_minus[i] = Z_minus*fugac_minus*exp( Z_minus*pot_elec[i]);
			rho_elec_polym[i] = 0.0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) rho_elec_polym[i] += Alpha_i[k][j]*PHA[k][j*Nx+i];
			rho_elec_polym[i] *= 1.0/v0;
		}
//	CHECK_VAL_N(rho_elec_plus[0]);
//	CHECK_VAL_N(rho_elec_minus[0]);
//	CHECK_VAL_N(rho_elec_polym[0]);

		for (i=0; i<Nx; i++){
			_mirror = (rho_elec_plus[i] - rho_elec_minus[i]  + rho_elec_polym[i]) / eps_prof[i] + c * pot_elec[i];
			fftw_in[i]       = _mirror;
			fftw_in[Nx2-1-i] = _mirror;
		}

		fftw_execute(p); // FFT forward (in --> out)

	//CHECK_VAL_N(fftw_in[0]);
		for (i=0; i<Nxc; i++) for(j=0;j<2;j++) fftw_out[i][j] /= (ksq[i] + c); 
//	for (i=0;i<Nx;i++) CHECK_VAL(fftw_out[i][1]);
//	printf("\n");
		fftw_execute(ip); // FFT back   (out --> in)
//	for (i=0;i<Nx;i++) CHECK_VAL(fftw_in[i]);
//	printf("\n");
//exit(1);
	//CHECK_VAL_N(fftw_in[0]);
		// Return only first halfspace
		for (i=0; i<Nx; i++) pot_elec[i] = fftw_in[i] / Nx2; 
//	CHECK_VAL_N(pot_elec[0]);

		for(i=0;i<Nx;i++) pot_elec[i]=wopt_PB*pot_elec[i]+(1.0-wopt_PB)*pot_elec_old[i];  // simople linear superposit // 
		
		err_max=0.0;
		for(i=0;i<Nx;i++)
		{
			err_check= (pot_elec[i]-pot_elec_old[i])*(pot_elec[i]-pot_elec_old[i]);
			if(err_check>err_max)err_max=err_check;
		}
		err_max=sqrt(err_max);

	} while(err_max>Sm_PB && iter_PB<MaxIT_PB);
	if (iter_PB == MaxIT_PB){ printf("Max_PB it %d. ", iter_PB); CHECK_VAL_N(err_max);}


	/*if (iter == 20){
		if (sign == 1)  write_qA("FFTqA.dat", Nx, Nsf, g);
		if (sign == -1) write_qA("FFTqcA.dat",Nx, Nsf, g);
		for (i=0; i<Nx; i++) { g2[i][start] = init[i]; g2[Nx2-1-i][start] = init[i]; }
		if (sign == 1)  write_qA("FFT2qA.dat", Nx2, Nsf, g2);
		if (sign == -1) write_qA("FFT2qcA.dat", Nx2, Nsf, g2);
	}*/

	//for(i=0; i<Nx; i++) free(g2[i]); free(g2);
	
	///// get electric field /////
	for(i=0;i<Nx;i++)
	{
		if((i>0)&&(i<Nx-1))
		{
			i1=i-1;
			i2=i+1;
			field_elec[i]=0.5*(pot_elec[i1]-pot_elec[i2])/dx;
		}
	}

	i=0;i2=i+1;
	field_elec[i]=field_elec[i2];

	i=Nx-1;i2=i-1;
	field_elec[i]=field_elec[i2];


	for(i=0;i<Nx;i++) field_elec_sq[i]=field_elec[i]*field_elec[i];
	free(ksq);
	free(pot_elec_old);
}

void PB_COS(double *phA, double **PHA, double **PHA_T, double *fftw_in, double *fftw_out){
	double *ksq, *pot_elec_old;
	double _mirror, c, err_check, err_max;
	int k, i1, i2, i, j;

	pot_elec_old = INIT(Nx);
	c = 1.0 / (L_Deby_S * L_Deby_S);
	// FFTW requires PBC; Calculate with z mirror and return only [0, Nz)
	ksq = INIT(Nx);
	for (i=0; i<Nx; i++)   ksq[i] = (i * 2.0*M_PI/(Nx2*dx)) * (i * 2.0*M_PI/(Nx2*dx)); // Fourier transform change of vars

	// Simple mixing it solver
	iter_PB = 0;
	do { 
		iter_PB += 1;

		for(i=0;i<Nx;i++){
			pot_elec_old[i]   = pot_elec[i];
			rho_elec_plus[i]  = Z_plus *fugac_plus *exp(-Z_plus*pot_elec[i]); 
			rho_elec_minus[i] = Z_minus*fugac_minus*exp( Z_minus*pot_elec[i]);
			rho_elec_polym[i] = 0.0;
			for (k=0;k<NF_N;k++) for (j=0;j<K_i[k];j++) rho_elec_polym[i] += Alpha_i[k][j]*PHA[k][j*Nx+i];
			rho_elec_polym[i] *= 1.0/v0;
		}
	//CHECK_VAL_N(rho_elec_plus[0]);
	//CHECK_VAL_N(rho_elec_minus[0]);
	//CHECK_VAL_N(rho_elec_polym[0]);

		for (i=0; i<Nx; i++) fftw_in[i] = (rho_elec_plus[i] - rho_elec_minus[i]  + rho_elec_polym[i]) / eps_prof[i] + c * pot_elec[i];
		fftw_execute(cfp); // FFT forward (in --> out)

	//CHECK_VAL_N(fftw_in[0]);
		for (i=0; i<Nx; i++) fftw_out[i] /= (ksq[i] + c); 
//	for (i=0;i<Nx;i++) CHECK_VAL(fftw_out[i][1]);
//	printf("\n");
		fftw_execute(cip); // FFT back   (out --> in)
//	for (i=0;i<Nx;i++) CHECK_VAL(fftw_in[i]);
//	printf("\n");
//exit(1);
	//CHECK_VAL_N(fftw_in[0]);
		// Return only first halfspace
		for (i=0; i<Nx; i++) pot_elec[i] = fftw_in[i] / (Nx2-2); 
	//CHECK_VAL_N(pot_elec[0]);

		for(i=0;i<Nx;i++) pot_elec[i]=wopt_PB*pot_elec[i]+(1.0-wopt_PB)*pot_elec_old[i];  // simople linear superposit // 
		
		err_max=0.0;
		for(i=0;i<Nx;i++)
		{
			err_check= (pot_elec[i]-pot_elec_old[i])*(pot_elec[i]-pot_elec_old[i]);
			if(err_check>err_max)err_max=err_check;
		}
		err_max=sqrt(err_max);

	} while(err_max>Sm_PB && iter_PB<MaxIT_PB);
	if (iter_PB == MaxIT_PB){ printf("Max_PB it %d. ", iter_PB); CHECK_VAL_N(err_max);}


	/*if (iter == 20){
		if (sign == 1)  write_qA("FFTqA.dat", Nx, Nsf, g);
		if (sign == -1) write_qA("FFTqcA.dat",Nx, Nsf, g);
		for (i=0; i<Nx; i++) { g2[i][start] = init[i]; g2[Nx-1-i][start] = init[i]; }
		if (sign == 1)  write_qA("FFT2qA.dat", Nx, Nsf, g2);
		if (sign == -1) write_qA("FFT2qcA.dat", Nx, Nsf, g2);
	}*/

	//for(i=0; i<Nx; i++) free(g2[i]); free(g2);
	
	///// get electric field /////
	for(i=0;i<Nx;i++)
	{
		if((i>0)&&(i<Nx-1))
		{
			i1=i-1;
			i2=i+1;
			field_elec[i]=0.5*(pot_elec[i1]-pot_elec[i2])/dx;
		}
	}

	i=0;i2=i+1;
	field_elec[i]=field_elec[i2];

	i=Nx-1;i2=i-1;
	field_elec[i]=field_elec[i2];


	for(i=0;i<Nx;i++) field_elec_sq[i]=field_elec[i]*field_elec[i];
	free(ksq);
	free(pot_elec_old);
}

void PhA_quad(int k, double delta, double **PHA, double **PHA_T, double **QA, double **QcA){
	int i,j,s;
	int ns;	
	
	for (i=0; i<Nx; i++) for (j=0; j<K_i[k]; j++){
		// Fourth order Gauss Quad without endpoints (Chantawansri 2011 Fredrickson)
		ns = Ns_i[k][j]; // End monomer
		PHA[k][j*Nx+i]  = 55.0/24 * QA[k][_IJS(i,j,1)]*QcA[k][_IJS(i,j,1)];
		PHA[k][j*Nx+i] += -1.0/6  * QA[k][_IJS(i,j,2)]*QcA[k][_IJS(i,j,2)];
		PHA[k][j*Nx+i] += 11.0/8  * QA[k][_IJS(i,j,3)]*QcA[k][_IJS(i,j,3)];
		for (s=4; s<=ns-4; s++){ 
			PHA[k][j*Nx+i]   += QA[k][_IJS(i,j,s)]*QcA[k][_IJS(i,j,s)];
		}
		PHA[k][j*Nx+i] += 11.0/8  * QA[k][_IJS(i,j,ns-3)]*QcA[k][_IJS(i,j,ns-3)];
		PHA[k][j*Nx+i] += -1.0/6  * QA[k][_IJS(i,j,ns-2)]*QcA[k][_IJS(i,j,ns-2)];
		PHA[k][j*Nx+i] += 55.0/24 * QA[k][_IJS(i,j,ns-1)]*QcA[k][_IJS(i,j,ns-1)];
		PHA[k][j*Nx+i] *= delta;
	}
//CHECK_VAL(PHA[k][0]);
	
	for (i=0; i<Nx;i++){
		PHA_T[k][i] = 0;
		for (j=0;j<K_i[k];j++){
			PHA_T[k][i] += PHA[k][j*Nx+i];
		}
	}
} 
