#include "setup.h"
#include "gets.h"
#include "report.h"

/* GLOBALS */
//SETUP()
	char pattern[30], phname[30], phname_elec[30], Wname[30], itname[30];
//ASSIGN_CONSTS()
	double fugac_plus, fugac_minus, L_Deby_S, L_Deby_P;
	double eps_P, eps_S, zs, freeEnergy_bulk;
	int Nx, Nx_2, Nx2, Nxc, Ns0;
	int Ny, Ny_2, Nz, Nz_2, Nz_cutoff;
	int NxNyNz, NxNyNz_8, NxNy, NxNy_4, NxNyNz_cutoff;
	double integ_cons;
	int OMP_NUM_THREADS;
//CALLOC_FIELDS()
	double *rx, *rx_sq, *ry, *ry_sq, *rz, *rz_sq;
	double *pot_elec, *eps_prof, *field_elec, *field_elec_sq;
	double *rho_elec_polym, *rho_elec_plus, *rho_elec_minus;
	double *phA, **PHA, **PHA_T, **WA;
	double *phB, *wB, *eta;
	double ***DAs, ***WAs, **DBs, **WBs;
//LOAD_AA()
	int NF_N, K, *K_i, Nm, **Nm_i, **Ns_i;
	double **Alpha_i, **Chi_i, *sigma_i;
//xy_NBC()
	int **IJK_PERP; int N_PERP;
	double *pSURF;

void SETUP(void){
	double tic = omp_get_wtime();
	stdout_redir_setup();
	stdout_redir_start();
	LOAD_AA(seq_file);
	ASSIGN_CONSTS();
	REPORT_SETUP();
	CALLOC_FIELDS();
	INIT_FIELDS(init_opt);
	double toc = omp_get_wtime();
	//printf("[setup.c/SETUP]: %f\n", toc-tic);	
}

void ASSIGN_CONSTS(void){
	long int ijk;
	double eps_0;

	// OMP
	if (getenv("OMP_NUM_THREADS")!=NULL) OMP_NUM_THREADS = atoi(getenv("OMP_NUM_THREADS")); //From environment
	else { 
		printf("*** OMP_NUM_THREADS NOT FOUND; DEFAULT TO ONE ***\n");
		OMP_NUM_THREADS = 1;
	}
	omp_set_num_threads(OMP_NUM_THREADS);

	// Names for report
	c0_plus_n = c0_plus;

	// Phys Props
	c0_plus *= 1e3 * avo_num * pow(Length, 3.0); 
	c0_minus*= 1e3 * avo_num * pow(Length, 3.0); 
	fugac_plus  = c0_plus;
	fugac_minus = c0_minus;
	eps_0 = eps_vacuum * kB * Temp * Length / Charge / Charge;
	eps_P = eps_0 * eps_r_P;
	eps_S = eps_0 * eps_r_S;
	L_Deby_S=sqrt(eps_S/(Z_plus*Z_plus*c0_plus+Z_minus*Z_minus*c0_minus));  
	L_Deby_P=sqrt(eps_P/(Z_plus*Z_plus*c0_plus+Z_minus*Z_minus*c0_minus));  

	// Numerical 
	Nz = roundl(lz/dz);
	//Nz_cutoff = roundl(lz_cutoff/dz);

	/*if (HALF_BOOL) xy_NBC_HALF(); // Set Nx, Ny, pSURF based on sigma, dx
	else           xy_NBC(); */
	
	/* WHEN NOT USING xy_NBC() SPECIFYING GRAFTING POINTS */
	Nx = roundl(lx/dx);
	Ny = roundl(ly/dy);
	NxNyNz = Nx*Ny*Nz;
	NxNy = Nx*Ny; 

	NxNyNz_cutoff = Nx*Ny*Nz;

	Nx2 = Nx * 2; Nxc = Nx + 1;
	integ_cons = dx*dy*dz/v0;
	
printf("Nx, Ny, Nz: %d %d %d\n", Nx, Ny, Nz);
printf("NxNy: %d, NxNyNz: %d\n", NxNy, NxNyNz);

	Ns0 = roundl(1.0*Nm*ds0);

	// Fields
	zs = exp(mu_s);

	freeEnergy_bulk = 0.0;
	#pragma omp parallel for reduction(+:freeEnergy_bulk)
	for (ijk=0; ijk<NxNyNz; ijk++) freeEnergy_bulk -= c0_plus + c0_minus;
	freeEnergy_bulk *= integ_cons*v0;
}

void CALLOC_FIELDS(void){
	int i, k, n, X;
	// System
	rx = INIT(Nx); rx_sq = INIT(Nx);
	ry = INIT(Ny); ry_sq = INIT(Ny);
	rz = INIT(Nz); rz_sq = INIT(Nz);

	// Electric
	pot_elec = INIT(NxNyNz); eps_prof = INIT(NxNyNz); field_elec = INIT(NxNyNz); field_elec_sq = INIT(NxNyNz);
	rho_elec_polym = INIT(NxNyNz); rho_elec_plus = INIT(NxNyNz); rho_elec_minus = INIT(NxNyNz);

	// Polymer
	phA = INIT(NxNyNz); PHA_T = INIT_2D(PHA_T, NF_N, NxNyNz); PHA = INIT_2D(PHA, NF_N, K_i[i]*NxNyNz); 
	WA = INIT_2D(WA, NF_N, K_i[i]*NxNyNz);

	// Solvent
	phB = INIT(NxNyNz); wB = INIT(NxNyNz); eta = INIT(NxNyNz);

	// Anderson Mixing Histories
	DAs = (double ***)calloc(and_NrMax+1, sizeof(double)); 
	for (n=0; n<and_NrMax+1; n++){ 
		DAs[n] = (double **)calloc(NF_N, sizeof(double));
		for (X=0; X<NF_N; X++)	DAs[n][X] = (double *)calloc(K_i[X]*NxNyNz, sizeof(double));	
	}

	WAs = (double ***)calloc(and_NrMax+1, sizeof(double)); 
	for (n=0; n<and_NrMax+1; n++){
		WAs[n] = (double **)calloc(NF_N, sizeof(double));
		for (X=0; X<NF_N; X++) WAs[n][X] = (double *)calloc(K_i[X]*NxNyNz, sizeof(double));
	}
	DBs = (double **)calloc(and_NrMax+1, sizeof(double));
	for (n=0; n<and_NrMax+1; n++) DBs[n] = (double *)calloc(NxNyNz, sizeof(double));
	WBs = (double **)calloc(and_NrMax+1, sizeof(double));
	for (n=0; n<and_NrMax+1; n++) WBs[n] = (double *)calloc(NxNyNz, sizeof(double));
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
	fscanf(fp, "%d\n", &K_i[k]);
	Nm_i[k] = (int *)calloc(K_i[k], sizeof(int));
	Ns_i[k] = (int *)calloc(K_i[k], sizeof(int));
	Alpha_i[k] = INIT(K_i[k]);
	Chi_i[k] = INIT(K_i[k]);
	for (j=0; j<K_i[k]; j++){
		fscanf(fp, "[%*d, %d], %lf, %lf\n", &Nm_i[k][j], &Alpha_i[k][j], &Chi_i[k][j]);
	        orig_Nm = Nm_i[k][j];
		Nm_i[k][j] = round(Nm_i[k][j]);// TO MAINTAIN END-END: /b0);
	}
	if (Nm_i[k][K_i[k]-1] > Nm) Nm = Nm_i[k][K_i[k]-1]; //Check last at last value of i = K_i[j]-1

	printf("\t%d Kuhn N =  %d [nm]\n", Nm_i[k][K_i[k]-1], orig_Nm);

	printf("\tAlpha[%d]: [", K_i[k]);
	for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Alpha_i[k][j]);
	printf("%.4f", Alpha_i[k][j]); printf("]\n");

	printf("\tChi[%d]: [", K_i[k]);
	for (j=0;j<K_i[k]-1;j++) printf("%.4f, ", Chi_i[k][j]);
	printf("%.4f", Chi_i[k][j]); printf("]\n");
	
	*i += 1;
}

void INIT_FIELDS(int init_opt){
	int i, j, k, Y, X;
	long int ij, ijk;
	
	for (i=0; i<Nx; i++) {
		rx[i] = i * dx;
		rx_sq[i] = rx[i] * rx[i];
	}
	
	for (j=0; j<Ny; j++) {
		ry[j] = j * dy;
		ry_sq[j] = ry[j] * ry[j];
	}
	for (k=0; k<Nz; k++) {
		rz[k] = k * dz;
		rz_sq[k] = rz[k] * rz[k];
	}

	printf("INIT_OPT: %d\n", init_opt);
	switch(init_opt) {
		case -2: 
			printf("  PHAz\n");
			ANALY_PB();
			INIT_PHA_PHAz(Win);			
			break;

		case -1:
			printf("  WAz\n");
			INIT_PHA_WAz(Win, win_box);
			get_ELFIELD();
			break;

		case 0 : 
			printf("  CONST\n");
			ANALY_PB();
			INIT_PHA_CONST(xCmax);
			break;

		case 1:
			printf("  WA spread\n");
			INIT_PHA_WA_SPREAD(Win, win_box);
			get_ELFIELD();
			break;

		default : 
			printf("  WA empty\n");
			INIT_PHA_WA(Win, win_box);
			get_ELFIELD();
	}

}

void INIT_PHA_WAz(char *Win, int win_box) {
	FILE *fp;
	long int ijk, ij;
	long int ijk0, Yijk0;
	int k, X, Y;
	double eta1, eta2;

	if (win_box == 0) {   // Win boxsize same as current
		win_nx = Nx; 
		win_ny = Ny; 
		win_nz = Nz;
	}

	fp=fopen(Win,"r"); 
	for (k=0; k<Nz; k++){
		ijk0 = (0*Ny + 0)*Nz + k; // x = 0, y = 0, z = k*dz
		ijk = ijk0; // Copy to ijk for use in _Yijk macros

		if (k >= win_nz) { // If empty
			wB[ijk]  = mu_s; 
			eta[ijk] = -mu_s;
			for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++) WA[X][_Yijk] = Chi_i[X][Y] + mu_s;
		}
		else {
			// READ WAz
			for (X=0; X<NF_N; X++) {
				fscanf(fp, "[ ");
				for (Y=0; Y<K_i[X]-1; Y++) {fscanf(fp, "%lf ", &WA[X][_Yijk]); WA[X][_Yijk]*=(1.0+A_r*(rand()/RAND_MAX-0.5));}
				fscanf(fp, "%lf ] ", &WA[X][_Yijk]); 
				//WA[X][_Yijk]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
			}
			fscanf(fp, "%lf %lf %lf\n", &wB[ijk], &eta[ijk], &pot_elec[ijk]);
			//wB[ijk] *= 1.0+A_r*(rand()/RAND_MAX-0.5);
			//eta[ijk]*= 1.0+A_r*(rand()/RAND_MAX-0.5);
			//pot_elec[ijk]*= 1.0+A_r*(rand()/RAND_MAX-0.5);
		}
		
		// COPY TO XY
		for (ij=0; ij<NxNy; ij++){
			ijk = ij*Nz + k;
			wB[ijk] = wB[ijk0];
			eta[ijk]= eta[ijk0];
			pot_elec[ijk] = pot_elec[ijk0];
		}
		for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++) {
			Yijk0 = Y*NxNyNz + ijk0;
			for (ij=0; ij<NxNy; ij++) {
				ijk = ij*Nz + k;
				WA[X][_Yijk] = WA[X][Yijk0];	
			}
		}
	}

	INIT_ELECS();

	fclose(fp);
}

void INIT_PHA_WA_SPREAD(char *Win, int win_box) {
	FILE *fp;
	long int ijk, ijk0;
	int X, Y, i, j, k, i0, j0;

	if (win_box == 0) {   // Win boxsize same as current
		win_nx = Nx; 
		win_ny = Ny; 
		win_nz = Nz;
	}

	fp=fopen(Win,"r"); 
	for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) for (k=0; k<Nz; k++) {
		ijk = (i*Ny + j)*Nz + k;

		if (k >= win_nz){
			wB[ijk]  = mu_s; 
			eta[ijk] = -mu_s;
			for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++) WA[X][_Yijk] = Chi_i[X][Y] + mu_s;
		}

		else if (i >= win_nx || j >= win_ny) { // If not known 
			ijk0 = (i0*Ny + j0)*Nz + k; // Last known value
			wB[ijk]  = wB[ijk0];
			eta[ijk] = eta[ijk0];
			pot_elec[ijk] = pot_elec[ijk0];
			for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++) WA[X][_Yijk] = WA[X][Y*NxNyNz+ijk0];
		}
		
		else {
			i0 = i; j0 = j; // Store i,j of last known k value
			for (X=0; X<NF_N; X++) {
				fscanf(fp, "[ ");
				for (Y=0; Y<K_i[X]-1; Y++) {fscanf(fp, "%lf ", &WA[X][_Yijk]); WA[X][_Yijk]*=(1.0+A_r*(rand()/RAND_MAX-0.5));}
				fscanf(fp, "%lf ] ", &WA[X][_Yijk]); 
				WA[X][_Yijk]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
			}
			fscanf(fp, "%lf %lf %lf\n", &wB[ijk], &eta[ijk], &pot_elec[ijk]);
			wB[ijk] *= 1.0+A_r*(rand()/RAND_MAX-0.5);
			eta[ijk]*= 1.0+A_r*(rand()/RAND_MAX-0.5);
			pot_elec[ijk]*= 1.0+A_r*(rand()/RAND_MAX-0.5);
		}
	}

	INIT_ELECS();

	fclose(fp);
}


void INIT_PHA_WA(char *Win, int win_box) {
	FILE *fp;
	long int ijk;
	int X, Y, i, j, k;

	if (win_box == 0) {   // Win boxsize same as current
		win_nx = Nx; 
		win_ny = Ny; 
		win_nz = Nz;
	}

	fp=fopen(Win,"r"); 
	for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) for (k=0; k<Nz; k++) {
		ijk = (i*Ny + j)*Nz + k;

		if (i >= win_nx || j >= win_ny || k >= win_nz) { // If empty
			wB[ijk]  = mu_s; 
			eta[ijk] = -mu_s;
			for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++) WA[X][_Yijk] = Chi_i[X][Y] + mu_s;
		}
		else {
			for (X=0; X<NF_N; X++) {
				fscanf(fp, "[ ");
				for (Y=0; Y<K_i[X]-1; Y++) {fscanf(fp, "%lf ", &WA[X][_Yijk]); WA[X][_Yijk]*=(1.0+A_r*(rand()/RAND_MAX-0.5));}
				fscanf(fp, "%lf ] ", &WA[X][_Yijk]); 
				WA[X][_Yijk]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
			}
			fscanf(fp, "%lf %lf %lf\n", &wB[ijk], &eta[ijk], &pot_elec[ijk]);
			wB[ijk] *= 1.0+A_r*(rand()/RAND_MAX-0.5);
			eta[ijk]*= 1.0+A_r*(rand()/RAND_MAX-0.5);
			pot_elec[ijk]*= 1.0+A_r*(rand()/RAND_MAX-0.5);
		}
	}

	INIT_ELECS();

	fclose(fp);
}

void INIT_PHA_CONST(double xCmax) {
	long int ijk, ij;
	int i, j, k, X, Y;

	if (xCmax == 0)  xCmax = 0.80; //default
	//for (k=1;k<Nz-1;k++) { //no pen edges
	for (k=0;k<Nz;k++) {
		for (ij=0; ij<NxNy; ij++) {
			ijk = ij*Nz + k;
			for (X=0; X<NF_N; X++) {
				PHA_T[X][ijk] = 0;
				PHA[X][0*NxNyNz + ijk] = xCmax * Nm_i[X][0]/Nm;
				PHA_T[X][ijk] += PHA[X][0*NxNyNz + ijk];
				for (Y=1; Y<K_i[X]; Y++){
					PHA[X][_Yijk] = xCmax*(1.0*(Nm_i[X][Y]-Nm_i[X][Y-1])/Nm);
					PHA_T[X][ijk] += PHA[X][_Yijk];
				}
			}
		}
	}
for (ijk=0;ijk<NxNyNz;ijk++) pot_elec[ijk]=0;
//write_elec();exit(1);
/*
	for (ij=0; ij<NxNy; ij++) {
		for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++) {
			ijk = ij*Nz + 0; // No pen wall
			PHA[X][_Yijk] = 0.0;
			PHA_T[X][ijk] = 0.0;
			ijk = ij*Nz + Nz-1; // No pen wall
			PHA[X][_Yijk] = 0.0;
			PHA_T[X][ijk] = 0.0;
		}
	}
*/
	for (ijk=0; ijk<NxNyNz; ijk++){
		phA[ijk] = 0;
		for (X=0; X<NF_N; X++) phA[ijk] += PHA_T[X][ijk]; 

		phB[ijk]=1.0-phA[ijk]; 
		eps_prof[ijk]=phA[ijk]*eps_P+(1.0-phA[ijk])*eps_S;

		for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++){	
			WA[X][_Yijk] = Chi_i[X][Y]*phB[ijk] - Alpha_i[X][Y]*pot_elec[ijk];
		}
		
		wB[ijk] = 0;
		for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) wB[ijk] += Chi_i[X][Y]*PHA[X][_Yijk];			

		//for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++)	WA[X][_Yijk] *= (1.0+A_r*(rand()/RAND_MAX-0.5));
		//wB[ijk]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
	}

}

void INIT_PHA_xC(double xCmax, double xCltot, double xCneck) {
printf("NOT TRANSFERRED TO 3D\n"); exit(1);
	int i, j, k;

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
}

void INIT_PHA_PHAz(char *Win) {
	FILE *fp;
	long int ijk, ij;
	long int ijk0, Yijk0;
	int k, X, Y;
	double eta1, eta2;

	fp = fopen(Win, "r");
	for (k=0; k<Nz; k++) {
		ijk0 = (0*Ny + 0)*Nz + k; // x = 0, y = 0, z = k*dz
		ijk = ijk0; // Copy to ijk for use in _Yijk macros
		// READ PHA
		fscanf(fp, "%*f %lf ", &phA[ijk0]); 
		for (X=0; X<NF_N; X++) {
			fscanf(fp, "%*lf [ "); // PHA_T
			for (Y=0; Y<K_i[X]-1; Y++) {
				fscanf(fp, "%lf ", &PHA[X][_Yijk]);
			}
			fscanf(fp, "%lf ] ", &PHA[X][_Yijk]); 
			fscanf(fp, "%*f\n"); // phB
		}
		// COPY TO XY
		for (ij=0; ij<NxNy; ij++) {
			ijk = ij*Nz + k;
			phA[ijk] = phA[ijk0];
		}
		for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++) {
			Yijk0 = Y*NxNyNz + ijk0;
			for (ij=0; ij<NxNy; ij++) {
				ijk = ij*Nz + k;	
				PHA[X][_Yijk] = PHA[X][Yijk0];
				PHA_T[X][ijk] += PHA[X][_Yijk];
			}
		}
	}
	fclose(fp);

	for (ijk=0; ijk<NxNyNz; ijk++){
		phB[ijk] = 1.0 - phA[ijk];
		wB[ijk] = -log(phB[ijk] / zs);
//
		eta1 = 0; eta2 = 0; 
		for (X=0;X<NF_N;X++){	
			eta1 += K_i[X];
			for (Y=0;Y<K_i[X];Y++){
				eta2 += (1-phA[ijk])*Chi_i[X][Y] - pot_elec[ijk]*Alpha_i[X][Y] + Chi_i[X][Y]*PHA[X][_Yijk] - WA[X][_Yijk];
			}	
		}
		eta[ijk] = 1/(eta1 + 1) * (eta2 - wB[ijk]);

//
		for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++){	
			WA[X][_Yijk] = Chi_i[X][Y]*phB[ijk] - Alpha_i[X][Y]*pot_elec[ijk]-eta[ijk];
		}
		//for (X=0; X<NF_N; X++) for (Y=0; Y<K_i[X]; Y++)	WA[X][_Yijk] *= (1.0+A_r*(rand()/RAND_MAX-0.5));
		//wB[ijk]*=(1.0+A_r*(rand()/RAND_MAX-0.5));
	}
}

void INIT_PHA_PHA(char *Win) {
printf("NOT TRANSFERRED TO 3D\n"); exit(1);
	FILE *fp;
	int i,j,k;
	fp = fopen(Win, "r");
	for (i=0; i<Nx; i++) {
		fscanf(fp, "%*f %lf ", &phA[i]); // lx, phA
		for (k=0; k<NF_N; k++) {
			fscanf(fp, "%*f [ "); // PHA_T
			for (j=0; j<K_i[k]-1; j++) {
				fscanf(fp, "%lf ", &PHA[k][j*Nx + i]);
			}
			fscanf(fp, "%lf ] ", &PHA[k][j*Nx + i]); 
			fscanf(fp, "%*f\n"); // phB
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

	int i,j,k;
	long int ijk, ij0, ij1;
	for(k=1;k<Nz;k++) {
		for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) {
			ijk = (i*Ny + j)*Nz + k;
			if(rz[k]<R_p)  pot_elec[ijk]=A1/rz[k]*(exp(-rz[k]/L_Deby_P)-exp(rz[k]/L_Deby_P))+0.5*rho_fix/c0_plus; 
			else  pot_elec[ijk]=A2/rz[k]*exp(-rz[k]/L_Deby_S);
		}
		
	}

	for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) {
		ij0 = (i*Ny + j)*Nz + 0;  // k = 0
		ij1 = (i*Ny + j)*Nz + 1;  // k = 1 
		pot_elec[ij0] = pot_elec[ij1]; // Neumann
	}

	get_ELFIELD();
}

void INIT_ELECS(void) {
	int ijk, X, Y;
	for(ijk=0;ijk<NxNyNz;ijk++){
		rho_elec_plus[ijk]  = Z_plus *fugac_plus *exp(-Z_plus*pot_elec[ijk]); 
		rho_elec_minus[ijk] = Z_minus*fugac_minus*exp( Z_minus*pot_elec[ijk]);
	}
	ELEC_POLYM();
}

void ELEC_POLYM(void) {
	int ijk, X, Y;
	for (ijk=0; ijk<NxNyNz; ijk++){
		rho_elec_polym[ijk] = 0.0;
		for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) {
			#pragma omp atomic
			rho_elec_polym[ijk] += Alpha_i[X][Y]*PHA[X][_Yijk];
		}
		rho_elec_polym[ijk] *= 1.0/v0;
	}
}


void xy_NBC(void) {
	// ONLY IMPLEMENTED FOR PURE (NF_N = 1) BRUSH
	// ONLY IMPLEMENTED FOR DX = DY
	// Return binary array of input to signify locations of grafted chains
	int x1, y1, x2, y2;
	double R = floor( sqrt(2.0/(3.0*sqrt(3)*sigma_i[0])) / dx) * dx;
	double s3R = floor( sqrt(3)*R / dx) * dx;
	int Deltax = round(3.0*R / dx);
	int Deltay = round(s3R / dy);

	int nx = 2 * floor(lx / (3.0*R));
	if (nx % 2 == 0) nx += 1; // nx must be odd
	lx = Deltax*dx*(nx-1)/2 + dx;

	int ny = floor(ly / s3R);
	if (ny % 2 == 1) ny += 1; // ny must be even
	ly = Deltay*dy*(ny-1) + dy;

	Nx = floor(Deltax * (nx-1)/2.0) + 1;
	Ny = floor(Deltay * (ny-1)) + 1;

	double s_calc = 2.0/(3.0*sqrt(3)*R*R);
	printf("xy_NBC CHANGES:\n");
	printf("Sigma_target: %.4f, Sigma_calc: %.4f (%.2f %), R: %.2f nm\n", sigma_i[0], s_calc, 100*(s_calc/sigma_i[0]-1.0), R);
	printf("nx (odd) : %d, Lx: %.2f\n", nx, lx);
	printf("ny (even): %d, Ly: %.2f\n", ny, ly);
	printf("Nx: %d, Ny: %d\n", Nx, Ny);

	sigma_i[0] = s_calc;

	pSURF = INIT(NxNy);
	double seed1x = 0, seed1y = ly-dy;
	x1 = round(seed1x/dx); 
	while (x1*dx <= lx) {
		y1 = round(seed1y/dy);
		pSURF[x1*Ny+y1] = 1.0;
		while (y1*dy >= 0.0) {
			pSURF[x1*Ny+y1] = 1.0;
			y1 -= Deltay; 
		}
		x1 += Deltax;
	}
	
	double seed2x = floor(3.0*R/2.0/dx)*dx, seed2y = floor(s3R/2.0/dy)*dy;
	x2 = round(seed2x/dx);
	while (x2*dx <= lx) { 
		y2 = round(seed2y/dy);
		pSURF[x2*Ny+y2] = 1.0;
		while (y2*dy <= ly) {
			pSURF[x2*Ny+y2] = 1.0;
			y2 += Deltay;
		}
		x2 += Deltax;
	}
}

void xy_NBC_HALF(void) {
	// ONLY IMPLEMENTED FOR PURE (NF_N = 1) BRUSH
	// ONLY IMPLEMENTED FOR DX = DY
	// Return input binary array to signify locations of grafted chains
	// Even at x, y boundaries, PBC at x0 = 2*Lx, y0 = 2*Ny
	int n;
	int x1, y1, x2, y2;
	double R = floor( sqrt(2.0/(3.0*sqrt(3)*sigma_i[0])) / dx) * dx;
	double s3R = floor( sqrt(3)*R / dx) * dx;
	int Deltax = round(3.0*R / dx), Deltay = round(s3R / dy);
//printf("    "); CHECK_VAL_N(s3R);
//printf("Deltax: %d, Deltay: %d\n", Deltax, Deltay);

//printf("lx %f, ly %f\n", lx, ly);
	
	int nx = 2 * floor(lx / (3.0*R));
	if (nx % 2 == 0) nx += 1; // nx must be odd
	int nx_2 = (nx-1) / 2 + 1;
	lx = Deltax*dx*(nx_2-1)/2 + dx;
	int ny = floor(ly / s3R);
	if (ny % 2 == 1) ny += 1; // ny must be even

	int ny_2 = ny/2;
	ly = 1.5*Deltay*dy*(ny_2-1) + dy;

	Nx = floor(Deltax * (nx_2-1)/2.0) + 1;
	Ny = floor(1.5*Deltay * (ny_2-1)) + 1;

	if (floor(ly/s3R) == 2) {
		ny = floor(2*ly/s3R) / 2;
		ly = Deltay*dy*(ny-1)+dy;
		Ny = floor(Deltay*(ny-1))+1;
	}

	NxNy = Nx*Ny;
	NxNyNz = Nx*Ny*Nz;

	IJK_PERP = (int **)calloc(1, sizeof(int));
	IJK_PERP[0] = (int *)calloc(nx_2*ny_2, sizeof(int)); 

	double s_calc = 2.0/(3.0*sqrt(3)*R*R);
	printf("xy_NBC_HALF CHANGES:\n");
	printf("\tSigma_target: %.4f, Sigma_calc: %.4f (%.2f %), R: %.2f nm\n", sigma_i[0], s_calc, 100*(s_calc/sigma_i[0]-1.0), R);
	printf("\tnx (odd, %d) --> nx_2: %d, Lx (/2): %.2f\n", nx, nx_2, lx);
	printf("\tny (even, %d) --> ny_2: %d, Ly (/2): %.2f\n", ny, ny_2, ly);
	printf("\tNx: %d, Ny: %d\n", Nx, Ny);

	sigma_i[0] = s_calc;

	pSURF = INIT(NxNy);
	double seed1x = 0, seed1y = 0;
	x1 = round(seed1x/dx);   (x1*Ny+y1) * Nz + 1;
	n = -1;
	while (x1*dx <= lx) {
		y1 = round(seed1y/dy);
		while (y1*dy <= ly) {
			pSURF[x1*Ny+y1] = 1.0;
			printf("x1, y1: %d, %d\n", x1, y1);
			n += 1;
			IJK_PERP[0][n] = (x1*Ny+y1)*Nz + 1;
			y1 += Deltay; 
		}
		x1 += Deltax;
	}
	
	double seed2x = floor(3.0*R/2.0/dx)*dx, seed2y = floor(s3R/2.0/dy)*dy;
	x2 = round(seed2x/dx);
	while (x2*dx <= lx) { 
		y2 = round(seed2y/dy);
		while (y2*dy <= ly) {
			pSURF[x2*Ny+y2] = 1.0;
			printf("x2, y2: %d, %d\n", x2, y2);
			if (IJK_PERP[0][n] != IJK_PERP[0][n-1]) {
				n += 1;
				IJK_PERP[0][n] = (x2*Ny+y2)*Nz + 1;
			}
			y2 += Deltay;
		}
		x2 += Deltax;
	}

	N_PERP = n+1; //Add to transfrom index to count
	if (N_PERP != nx_2*ny_2) {printf("ERROR: N_PERP (%d) != nx_2 * ny_2 (%d), Exiting...\n", N_PERP, nx_2*ny_2); exit(1);}

	for (n=0;n<N_PERP;n++) printf("%d ", IJK_PERP[0][n]);
	printf("\n");
//exit(1);
}

void CLEAN(void){
	int k;
	
	stdout_redir_end();

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
