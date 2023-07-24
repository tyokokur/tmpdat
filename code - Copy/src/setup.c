#include "setup.h"
#include "gets.h"
#include "report.h"

void SETUP(void){
	ASSIGN_CONSTS();
	LOAD_AA(seq_file);
	REPORT_SETUP();
	CALLOC_FIELDS();
	INIT_FIELDS(init_opt);
}

void ASSIGN_CONSTS(void){
	int i;

	// Phys Props
	c0_plus *= 1e-03 * avo_num * pow(Length, 3); 
	c0_minus*= 1e-03 * avo_num * pow(Length, 3); 
	fugac_plus  = c0_plus;
	fugac_minus = c0_minus;
	L_Deby_S=sqrt(eps_vacuum*eps_r_S/(Z_plus*Z_plus*c0_plus+Z_minus*Z_minus*c0_minus));  
	L_Deby_P=sqrt(eps_vacuum*eps_r_P/(Z_plus*Z_plus*c0_plus+Z_minus*Z_minus*c0_minus));  
	eps_P = eps_vacuum * eps_r_P;
	eps_S = eps_vacuum * eps_r_S;

	// Fields
	zs = exp(mu_s);

	freeEnergy_bulk = 0.0;
	for (i=0; i<Nx; i++) freeEnergy_bulk -= c0_plus + c0_minus;
	freeEnergy_bulk *= integ_cons*v0;

	// Numerical 
	Nx = roundl(lx/dx);
	Nx2 = Nx * 2;	
	Nxc = Nx + 1;
	Ns0 = roundl(1.0*Nm*ds0);
	integ_cons = dx/v0;
}

void CALLOC_FIELDS(void){
	int i, k, n;
	// System
	rx = INIT(Nx); rx_sq = INIT(Nx);

	// Electric
	pot_elec = INIT(Nx); eps_prof = INIT(Nx); field_elec = INIT(Nx); field_elec_sq = INIT(Nx);
	rho_elec_polym = INIT(Nx); rho_elec_plus = INIT(Nx); rho_elec_minus = INIT(Nx);

	// Polymer
	phA = INIT(Nx); PHA_T = INIT_2D(PHA_T, NF_N, Nx); PHA = INIT_2D(PHA, NF_N, K_i[i]*Nx); 
	WA = INIT_2D(WA, NF_N, K_i[i]*Nx);

	// Solvent
	phB = INIT(Nx); wB = INIT(Nx); eta = INIT(Nx);

	// Anderson Mixing Histories
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
}

void LOAD_AA(char *seq_file){
	char seq[30];
	int i;
	FILE *fp = fopen(seq_file, "r");

	fscanf(fp, "%s\n", seq);
	NF_N = strlen(seq);
	sigma_i = INIT(NF_N);
	K_i = (int *)calloc(NF_N, sizeof(double));
	Nm_i = (int **)calloc(NF_N, sizeof(int));
	Ns_i = (int **)calloc(NF_N, sizeof(int));
	Alpha_i = (double **)calloc(NF_N, sizeof(double));
	Chi_i = (double **)calloc(NF_N, sizeof(double));

	i = 0;
	if (strstr(seq, "L") != NULL){ printf("NFL:\n"); READ_AA(&i, fp); }
	if (strstr(seq, "M") != NULL){ printf("NFM:\n"); READ_AA(&i, fp); }
	if (strstr(seq, "H") != NULL){ printf("NFH:\n"); READ_AA(&i, fp); }

	fclose(fp);
}

void READ_AA(int *i, FILE *fp){
	int k = *i, j, orig_Nm;
	fscanf(fp, "%lf\n", &sigma_i[k]);
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

	printf("\t%d Kuhn N =  %d AA  / %.2f [nm/b0]\n", Nm_i[k][K_i[k]-1], orig_Nm, b0);

	printf("\tAlpha[%d]: [", K_i[k]);
	for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Alpha_i[k][j]);
	printf("%.4f", Alpha_i[k][j]); printf("]\n");

	printf("\tChi[%d]: [", K_i[k]);
	for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Chi_i[k][j]);
	printf("%.4f", Chi_i[k][j]); printf("]\n");
	
	*i += 1;
}

void INIT_FIELDS(int init_opt){
	int i, j, k, i1, i2;
	double A_r = 2.3e-16; // Field perturbation constant
	FILE *fp;
	
	for (i=0; i<Nx; i++) {
		rx[i] = i * dx;
		rx_sq[i] = rx[i] * rx[i];
	}
	
	printf("INIT_OPT: %d\n", init_opt);
	switch(init_opt) {
		case -2: // Init from pha
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
			break;

		case -1: // Init from wA
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

			get_ELFIELD();

			break;

		case 0 : // Init as const
			if (xCmax == 0)  xCmax = 0.80; //default
			for (i=1; i<Nx; i++){ //Respect no pen at x = 0
				for (k=0; k<NF_N; k++) {
					PHA_T[k][i] = 0;
					PHA[k][0*Nx + i] = xCmax * Nm_i[k][0]/Nm;
					PHA_T[k][i] += PHA[k][0*Nx + i];
					for (j=1; j<K_i[k]; j++){
						PHA[k][j*Nx + i] = xCmax*(1.0*(Nm_i[k][j]-Nm_i[k][j-1])/Nm);
						PHA_T[k][i] += PHA[k][j*Nx + i];
					}
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
			}
			ANALY_PB();
			//write_ph(phA, PHA_T, PHA, phB); write_W(WA, wB, eta); write_elec(pot_elec, rho_elec_plus, rho_elec_minus, rho_elec_polym); exit(0);
			break;

		default : 
			if (xCmax == 0) xCmax = 0.80; //defaults
			if (xCltot==0)  xCltot = lx/2.0;
			if (xCneck==0)  xCneck = xCltot/(2.0*init_opt);
			double w;
			int c;

			w = (xCltot - (init_opt-1)*xCneck)/init_opt; 
			for (i=0; i<Nx; i++){
				for (k=0; k<NF_N; k++) {
					PHA_T[k][i] = 0;
					PHA[k][0*Nx + i] = 0; for (c=1; c<=init_opt; c++) PHA[k][0*Nx + i] += xCmax*(Nm_i[k][0]/Nm) * exp(-1.0/2* pow((i*dx-c*1.25* (w+xCneck)/2), 2) / sqrt(w/xCltot));
					for (j=1; j<K_i[k]; j++){
						PHA[k][j*Nx + i] = 0;
						for (c=1; c<=init_opt; c++) PHA[k][j*Nx + i] += xCmax*((Nm_i[k][j]-Nm_i[k][j-1])/Nm) * exp(-1.0/2* pow((i*dx-c*1.25* (w+xCneck)/2), 2) / sqrt(w/xCltot));
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
			}
			ANALY_PB();
			//write_ph(phA, PHA_T, PHA, phB); write_W(WA, wB, eta); write_elec(pot_elec, rho_elec, rho_elec_minus, rho_elec_polym); exit(0);
	}

}

void ANALY_PB(void){
	// Analytical PB solution (spherical)
	double rho_fix = Alpha_i[0][0]/v0;
	double R_p=pow(0.75/M_PI*Nm*v0, 1.0/3);
	double a11=2.0*eps_P*(sinh(R_p/L_Deby_P)-R_p/L_Deby_P*cosh(R_p/L_Deby_P));
	double a12=eps_S*(R_p/L_Deby_S+1.0)*exp(-R_p/L_Deby_S);
	double a21=2.0*sinh(R_p/L_Deby_P);
	double a22=exp(-R_p/L_Deby_S);
	double b1=0.0;
	double b2=0.5*rho_fix*R_p/c0_plus;
	double det=a11*a22-a12*a21;
	double det1=b1*a22-b2*a12;
	double det2=a11*b2-a21*b1;
	double A1=det1/det;
	double A2=det2/det;

	int i;
	for(i=1;i<Nx;i++) {
		if(rx[i]<R_p)  pot_elec[i]=A1/rx[i]*(exp(-rx[i]/L_Deby_P)-exp(rx[i]/L_Deby_P))+0.5*rho_fix/c0_plus; 
		else  pot_elec[i]=A2/rx[i]*exp(-rx[i]/L_Deby_S);
	}
	pot_elec[0]=pot_elec[1]; // Neumann

	get_ELFIELD();
}

void CLEAN(void){
	int k;

	// System
	free(rx); free(rx_sq);

	// Electric
	free(pot_elec); free(eps_prof); free(field_elec); free(field_elec_sq);
	free(rho_elec_polym); free(rho_elec_plus); free(rho_elec_minus);

	// Polymer
	free(phA); 
	for (k=0; k<NF_N; k++) { free(PHA_T[k]); free(PHA[k]); free(WA[k]); }

	// Solvent
	free(phB); free(wB); free(eta);
}
