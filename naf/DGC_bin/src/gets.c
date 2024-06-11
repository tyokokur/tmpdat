#include "gets.h"
/* GLOBALS */
//get_freeE()
	double freeE, freeDiff, inCompMax;
	int MAXMAX;
//get_PHI()
	double Q2;
//And_mix()
	double and_err;
	double **WAdiff, **WAnew, **DAnew;
	double **u, **u_temp, *v;
	double *wBdiff, *wBnew, *DBnew, *Cs;

void get_3D_GAUSSPROP(long int ijk){
	/* TODO: PARALLELIZE EACH DIMENSION */
	// Given some inputted ijk
	// return the Gaussian chain distribution centered around ijk
	int i, j, k, ip, ig;
	double *xGAUSSPROP, *yGAUSSPROP, *zGAUSSPROP;
	xGAUSSPROP = INIT(Nx); yGAUSSPROP = INIT(Ny); zGAUSSPROP = INIT(Nz);
	
	k = ijk % Nz;
	j = (ijk - k)/Nz % Ny;
	i = ((ijk - k)/Nz - j)/Ny;
//printf("for ijk: %ld, k: %d\n", ijk, k);

//printf("Input i,j,k = %d,%d,%d\n", i,j,k);
	// X
	get_1D_GAUSSPROP(i, xGAUSSPROP, xGAUSSDIST, Nx);
//printf("xGAUSSPROP:\n");
//for (i=0;i<Nx;i++) printf("%f ", xGAUSSPROP[i]);
//printf("\n");

	// Y
	get_1D_GAUSSPROP(j, yGAUSSPROP, yGAUSSDIST, Ny);
//printf("yGAUSSPROP:\n");
//for (j=0;j<Ny;j++) printf("%f ", yGAUSSPROP[j]);
//printf("\n");
	
	// Z
	get_1D_GAUSSPROP(k, zGAUSSPROP, zGAUSSDIST, Nz);
//printf("zGAUSSPROP:\n");
//for (k=0;k<Nz;k++) printf("%f ", zGAUSSPROP[k]);
//printf("\n");

	for (i=0;i<Nx;i++) for (j=0;j<Ny;j++) for (k=0;k<Nz;k++){
		ijk = (i*Ny+j)*Nz+k;
		GAUSSPROP[ijk] = xGAUSSPROP[i]*yGAUSSPROP[j]*zGAUSSPROP[k];
	}

	free(xGAUSSPROP); free(yGAUSSPROP); free(zGAUSSPROP);
}

void get_2D_GAUSSPROP(long int ij, double *mnGAUSSPROP, double *mGAUSSDIST, int N_m, double *nGAUSSDIST, int N_n){
	int i, j;
	j = ij % N_n;
	i = (ij - j) / N_n; 
	double *mGAUSSPROP, *nGAUSSPROP;
	mGAUSSPROP = INIT(N_m); nGAUSSPROP = INIT(N_n);

	get_1D_GAUSSPROP(i, mGAUSSPROP, mGAUSSDIST, N_m);
	get_1D_GAUSSPROP(j, nGAUSSPROP, nGAUSSDIST, N_n);
	
	for (i=0;i<N_m;i++) for (j=0;j<N_n;j++) mnGAUSSPROP[i*N_n+j] = mGAUSSPROP[i]*nGAUSSPROP[j];

	free(mGAUSSPROP); free(nGAUSSPROP);
}

void get_1D_GAUSSPROP(int i, double *mGAUSSPROP, double *mGAUSSDIST, int N){
	// Return 1D filter centered around index i
	// Output mGAUSSPROP, input mGAUSSDIST
	int ip, ig;

	mGAUSSPROP[i] = mGAUSSDIST[0];
	for (ip=i+1, ig=1; ip<N ; ip++,ig++) mGAUSSPROP[ip] = mGAUSSDIST[ig]; 
	for (ip=i-1, ig=1; ip>=0; ip--,ig++) mGAUSSPROP[ip] = mGAUSSDIST[ig]; 
}

void get_PHI(void) {
	long int ijk;
	int X, Y;

	// Polymer
	for (X=0; X<NF_N; X++) {
		int_PHA_Quad(X);
		//int_PHA_Reimm(X);
		get_int_Q1_Quen(X);
	}
	
	for (ijk=0; ijk<NxNyNz; ijk++) {
		phA[ijk] = 0.0;
		for (X=0; X<NF_N; X++) phA[ijk] += PHA_T[X][ijk];
	}

	// Solvent
	Q2 = 0.0;
	#pragma parallel for reduction(+:Q2) \
		shared(phB, NxNyNz, wB, zs) \
		private(ijk)
  	for(ijk=0;ijk<NxNyNz;ijk++) {
  		Q2 += exp(-wB[ijk]);
  		phB[ijk] = zs*exp(-wB[ijk]);
  	}
  	Q2 *= integ_cons;
}

double get_freeE(void) {
	double freeEnergy, freeW, freeU, freeS, free_elec; 
	double free_elec_polym, free_elec_laplace, free_elec_ion;
	double freeU_step, free_elec_polym_step;
	double psum, fpsum, eta1, eta2;
	int X, Y;
	long int ijk;

int nthreads;
double tic = omp_get_wtime();
	freeW = 0.0; freeU = 0.0; freeS = 0.0; free_elec_polym = 0.0; free_elec_ion = 0.0; free_elec_laplace = 0.0; inCompMax = 0.0;
	#pragma omp parallel for \
		reduction(+:freeW, freeU, freeS, free_elec_polym, free_elec_ion, free_elec_laplace) \
		shared(inCompMax, MAXMAX)\
		private(eta1, eta2, fpsum, free_elec_polym_step, freeU_step, ijk, psum, X, Y) 
	for(ijk=0;ijk<NxNyNz;ijk++)
	{
		psum=1.0-phA[ijk]-phB[ijk];

		fpsum=fabs(psum);
		#pragma omp critical
		{
		if(fpsum>inCompMax){inCompMax=fpsum; MAXMAX = ijk;}
		}

		eta1 = 0; eta2 = 0; 
		for (X=0;X<NF_N;X++){	
			eta1 += K_i[X];
			for (Y=0;Y<K_i[X];Y++){
				eta2 += (1-phA[ijk])*Chi_i[X][Y] - pot_elec[ijk]*Alpha_i[X][Y] + Chi_i[X][Y]*PHA[X][_Yijk] - WA[X][_Yijk];
			}	
		}

		eta[ijk] = 1/(eta1 + 1) * (eta2 - wB[ijk]);

/*
if (ijk==(10*Ny+10)*Nz+20) { 
	X = 0; Y = 0;
	CHECK_VAL_N(phA[ijk]);
	CHECK_VAL_N(Chi_i[X][Y]);
	CHECK_VAL_N(pot_elec[ijk]);
	CHECK_VAL_N(Alpha_i[X][Y]);
	CHECK_VAL_N(PHA[X][_Yijk]);
	CHECK_VAL_N(WA[X][_Yijk]);
	CHECK_VAL_N(eta1);
	CHECK_VAL_N(eta2);
	CHECK_VAL_N(eta[ijk]);
	exit(1);
}
*/
		
		freeU_step = 0;
		for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) freeU_step += Chi_i[X][Y]*PHA[X][_Yijk];
		freeU += freeU_step * phB[ijk];
	
		for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) freeW -= WA[X][_Yijk] * PHA[X][_Yijk];
		freeW -= wB[ijk]*phB[ijk];

		free_elec_polym_step = 0;
		for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) free_elec_polym_step += Alpha_i[X][Y]*PHA[X][_Yijk];
		free_elec_polym   -= free_elec_polym_step * pot_elec[ijk] / v0;
		free_elec_laplace -= 0.5 * eps_prof[ijk] * field_elec_sq[ijk];
		free_elec_ion     -= (rho_elec_plus[ijk] + rho_elec_minus[ijk]);
		nthreads = omp_get_num_threads();
	}
double toc = omp_get_wtime();
//printf("\t[gets.c/get_freeEnergy] (%d threads): %f\n", nthreads, toc-tic);

	freeU *= integ_cons; 
	freeW *= integ_cons; 

	free_elec_polym   *= integ_cons*v0; 
	free_elec_laplace *= integ_cons*v0; 
	free_elec_ion     *= integ_cons*v0; 
	free_elec_ion     -= freeEnergy_bulk; 
	free_elec          = free_elec_polym + free_elec_laplace + free_elec_ion; 

	freeS = -zs*Q2; 
	for (X=0;X<NF_N;X++) freeS -= int_Q1_Quen[X]; 

	freeEnergy=freeU+freeW+freeS+free_elec+freeEnergy_bulk; //+free_factorial; 

/*
double A = Nx*dx*Ny*dy;
CHECK_VAL_N(freeU/A);
CHECK_VAL_N(freeW/A);
CHECK_VAL_N(Q2);
CHECK_VAL_N(int_Q1_Quen[0]);
CHECK_VAL_N(freeS/A);
CHECK_VAL_N(free_elec_polym/A);
CHECK_VAL_N(free_elec_laplace/A);
CHECK_VAL_N(free_elec_ion/A);
CHECK_VAL_N(free_elec/A);

write_ph(); write_W();
exit(1);
*/

/*if (iter%10==0){
printf("\t"); CHECK_VAL_N(int_Q1_Quen[0]);
printf("\t"); CHECK_VAL_N(freeU);
printf("\t"); CHECK_VAL_N(freeW);
printf("\t"); CHECK_VAL_N(freeS);
}*/
//CHECK_VAL_N(free_elec_polym); 
//CHECK_VAL_N(free_elec_laplace);
//CHECK_VAL_N(free_elec_ion);
//CHECK_VAL_N(free_elec_ion);
//CHECK_VAL_N(free_elec); 
//CHECK_VAL_N(freeEnergy);

	return freeEnergy;
}

void get_EPS_PROF(void) {
	long int ijk;
	for (ijk=0; ijk<NxNyNz; ijk++) eps_prof[ijk] = phA[ijk]*eps_P + phB[ijk]*eps_S;
}

void int_PHA_Reimm(int X){
	long int ijk;
	int Y, s;
	for (ijk=0, PHA[X][ijk]=0.0;ijk<NxNyNz;ijk++) {
		for (s=0; s<=Ns_i[X][0]; s++) for (Y=0;Y<K_i[X];Y++) PHA[X][ijk] += QA[X][_IJS(ijk, Y, s)] * QcA[X][_IJS(ijk, Y, s)] * exp(WA[X][_Yijk]);
	}
//for (ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(PHA[X][ijk]);
//write_qA_s(0, "qA.dat",  NxNyNz, Ns_i[X][0], qA);
//write_qA_s(0, "qcA.dat", NxNyNz, Ns_i[X][0], qcA);
//exit(1);

	for (ijk=0; ijk<NxNyNz;ijk++){
		PHA_T[X][ijk] = 0;
		for (Y=0;Y<K_i[X];Y++){
			PHA_T[X][ijk] += PHA[X][_Yijk];
		}
	}	
}
/*
void int_PHA_Reimm(int X) {
	long int ijk;
	int Y, s;

	for(ijk=0;ijk<NxNyNz;ijk++){ 
		
		for (Y=0;Y<K_i[X];Y++){
			PHA[X][_Yijk] = 0;

			if (Y==0) {
				PHA[X][0*NxNyNz+ijk] += 0.50*QA[X][_IJS(ijk,0, 0 )]*QcA[X][_IJS(ijk,0, 0)];
				for (s=1;s<Ns_i[X][0];s++) PHA[X][0*NxNyNz+ijk] += QA[X][_IJS(ijk,0,s)]*QcA[X][_IJS(ijk,0,s)];
				PHA[X][0*NxNyNz+ijk] += 0.50*QA[X][_IJS( ijk, 0, Ns_i[X][0] )]*QcA[X][_IJS( ijk, 0, Ns_i[X][0] )];
			}
			else{
				PHA[X][_Yijk] += 0.50*QA[X][_IJS(ijk,Y, Ns_i[X][Y-1] )]*QcA[X][_IJS(ijk,Y, Ns_i[X][Y-1] )];
				for (s=Ns_i[X][Y-1]+1;s<Ns_i[X][Y];s++) PHA[X][_Yijk] += QA[X][_IJS(ijk,Y,s)]*QcA[X][_IJS(ijk,Y,s)];
				PHA[X][_Yijk] += 0.50*QA[X][_IJS(ijk,Y, Ns_i[X][Y] )]*QcA[X][_IJS(ijk, Y, Ns_i[X][Y] )];
			}
			
			PHA[X][_Yijk] *= ds0;
		}
	}
	for (ijk=0; ijk<NxNyNz;ijk++){
		PHA_T[X][ijk] = 0;
		for (Y=0;Y<K_i[X];Y++){
			PHA_T[X][ijk] += PHA[X][_Yijk];
		}
	}	
}*/

void int_PHA_Quad(int X) {
	// Fourth order Gauss Quad without endpoints (Chantawansri 2011 Fredrickson)
	long int ijk;
	int Y;
	int s, ni, nf;	
	double sum1, sum2;

int nthreads;
double tic = omp_get_wtime();
	for (Y=0; Y<K_i[X]; Y++) {
		if (Y!=0) ni = Ns_i[X][Y-1];
		else      ni = 0;
		nf = Ns_i[X][Y];

		#pragma omp parallel \
			shared(nf, ni, Y)\
			private(ijk, s, sum1, sum2)
		{ 
		#pragma omp for 
		for (ijk=0; ijk<NxNyNz; ijk++) {
			sum2 = 0.0;
			for (s=4; s<=nf-ni-4; s++){ 
				sum2 += QA[X][_IJS(ijk,Y,ni+s)]*QcA[X][_IJS(ijk,Y,ni+s)];
			}
			sum1  =      55.0/24 * QA[X][_IJS(ijk,Y,ni+1)]*QcA[X][_IJS(ijk,Y,ni+1)]
				   -  1.0/6  * QA[X][_IJS(ijk,Y,ni+2)]*QcA[X][_IJS(ijk,Y,ni+2)]
				   + 11.0/8  * QA[X][_IJS(ijk,Y,ni+3)]*QcA[X][_IJS(ijk,Y,ni+3)]
				   + 11.0/8  * QA[X][_IJS(ijk,Y,nf-3)]*QcA[X][_IJS(ijk,Y,nf-3)]
				   -  1.0/6  * QA[X][_IJS(ijk,Y,nf-2)]*QcA[X][_IJS(ijk,Y,nf-2)]
				   + 55.0/24 * QA[X][_IJS(ijk,Y,nf-1)]*QcA[X][_IJS(ijk,Y,nf-1)];

			PHA[X][_Yijk] = (sum1+sum2)*exp(WA[X][_Yijk]); // *ds0;
		}
		if (omp_get_thread_num() == 0) nthreads = omp_get_num_threads();
		}
	}


	#pragma omp parallel for \
		shared(PHA_T, PHA) \
		private(ijk, Y)
	for (ijk=0; ijk<NxNyNz;ijk++){
		PHA_T[X][ijk] = 0;
		for (Y=0;Y<K_i[X];Y++){
			PHA_T[X][ijk] += PHA[X][_Yijk];
		}
	}

double toc = omp_get_wtime();
//printf("\t[gets.c/int_PHA_Quad] (%d threads): %f\n", nthreads, toc-tic);
}

void get_int_Q1_Quen(int X) {

	int ij;
	int_Q1_Quen[X] = 0.0;
	int Y = 0; //TODO DGC TMP
	for (ij=0;ij<NxNy;ij++) int_Q1_Quen[X] += log( QcA[X][_IJS(ij*Nz+1, Y, 0)] ); 
	int_Q1_Quen[X] *= sigma_i[X]*dx*dy; 
	//for (ij=0;ij<NxNy;ij++) int_Q1_Quen[X] += log( QcA[X][_IJS(ij*Nz+1, 0, 0)]); // OLD 03/12
	//int_Q1_Quen[X] *= sigma_i[X]*dx*dy/v0; // OLD 03/12
			
}


// ORIGINAL
/*void get_int_Q1_Quen(int X) {
	long int ij;
	int k, Y, s;
	double Q1;

	int_Q1_Quen[X] = 0.0;

double tic = omp_get_wtime();
	for (ij=0; ij<NxNy; ij++) { // Could parallelize this but already fast enough
		Y = 0; s = 0;
		Q1 = 0.0;
		#pragma omp parallel for reduction(+:Q1)
		for (k=0; k<Nz; k++) Q1 += QcA[X][_IJS(ij*Nz+k, Y, s)] * deltaFun[k];
		#pragma omp atomic
		int_Q1_Quen[X] += log( Q1*dz/v0 );
	}
	int_Q1_Quen[X] *= dx*dy;
double toc = omp_get_wtime();
//printf("\t[gets.c/get_int_Q1_Quen]: %f\n", toc-tic);
}*/

void get_ELFIELD(void){
	int i,j,k;
	long int ijk, p, m;
	double dpdx, dpdy, dpdz;

	for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) for (k=0; k<Nz; k++) {
		ijk = (i*Ny + j)*Nz + k;

		if (i==0 || i==Nx-1 || j==0 || j==Ny-1 || k==0 || k==Nz-1) {
			field_elec[ijk] = 0.0; // Neumann 
		}
		else {
			p = ( (i+1) *Ny + j)*Nz + k; m = ( (i-1) *Ny + j)*Nz + k;
			dpdx = 0.5*(pot_elec[p]-pot_elec[m])/dx;
			p = (i*Ny + (j+1) )*Nz + k; m = (i*Ny + (j-1) )*Nz + k;
			dpdy = 0.5*(pot_elec[p]-pot_elec[m])/dy;
			p = (i*Ny + j)*Nz + k+1; m = (i*Ny + j)*Nz + k-1;
			dpdz = 0.5*(pot_elec[p]-pot_elec[m])/dz;

			field_elec[ijk] = sqrt( dpdx*dpdx + dpdy*dpdy + dpdz*dpdz );
		}

		field_elec_sq[ijk] = field_elec[ijk] * field_elec[ijk];
	}
}

double And_mix(double **WA, double *wB){ 
//Anderson Mixing
	//Input globals: PHA, Chi_i, Alpha_i, pot_elec, eta 
	//Output updated WA and wB
	//Updated globals DAs WAs DBs WBs
	double **WAdiff, **WAnew, *wBdiff, *wBnew, **DAnew, *DBnew;
	double **u, **u_temp, *v;
	double psum, lambda;
	double amax, bmax;
	double and_err;

	long int ijk;
	int X, Y, n, m;
	int and_Nr, end;

	int i;
	WAdiff = INIT_2D(WAdiff, NF_N, NxNyNz*K_i[i]);
	WAnew  = INIT_2D(WAnew , NF_N, NxNyNz*K_i[i]);
	wBdiff = INIT(NxNyNz); wBnew = INIT(NxNyNz);
	DAnew  = INIT_2D(DAnew, NF_N, NxNyNz*K_i[i]);
	DBnew  = INIT(NxNyNz);

	for (n=1; n<=and_NrMax; n++){ 
		for (X=0; X<NF_N; X++){ 
			for (Y=0;Y<K_i[X];Y++){
				#pragma omp parallel for \
					shared(n, X, Y) \
					private(ijk)
				for (ijk=0;ijk<NxNyNz;ijk++){ 
					DAs[n-1][X][_Yijk] = DAs[n][X][_Yijk];
					DBs[n-1][ijk] = DBs[n][ijk];
					WAs[n-1][X][_Yijk] = WAs[n][X][_Yijk];
					WBs[n-1][ijk] = WBs[n][ijk];
				}
			}
		}	
	}

	//Calculate and Update
	end = and_NrMax;

	amax = 0.0; 
	bmax = 0.0;
	#pragma omp parallel for \
		shared(amax, bmax)\
		private(ijk, X, Y)
	for (ijk=0;ijk<NxNyNz;ijk++){
		for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) {
			WAnew[X][_Yijk]    = Chi_i[X][Y]*phB[ijk] - Alpha_i[X][Y]*pot_elec[ijk] - eta[ijk];
			WAdiff[X][_Yijk]   = WAnew[X][_Yijk] - WA[X][_Yijk];
			DAs[end][X][_Yijk] = WAdiff[X][_Yijk]; //Update last val
			WAs[end][X][_Yijk] = WA[X][_Yijk]; //Update last val
			#pragma omp atomic
			amax += WAdiff[X][_Yijk]*WAdiff[X][_Yijk] / (WA[X][_Yijk]*WA[X][_Yijk]);
		}
		wBnew[ijk] = 0;
		for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) wBnew[ijk] += Chi_i[X][Y]*PHA[X][_Yijk];
		wBnew[ijk] -= eta[ijk];
		wBdiff[ijk] = wBnew[ijk]-wB[ijk];		
		DBs[end][ijk] = wBdiff[ijk]; //Update last val 
		WBs[end][ijk] = wB[ijk]; //Update last val
		#pragma omp atomic
		bmax += wBdiff[ijk]*wBdiff[ijk] / (wB[ijk]*wB[ijk]);
	}	
	and_err = amax+bmax;
	and_err = pow(and_err, 0.50);
	
	if (iter%and_it!=0 || and_NrMax == 0){// || and_err>1e-02 || bmax>1e-03 || amax>1e-3){ 
		#pragma omp parallel for \
			private(ijk, psum, X, Y)
		for (ijk=0; ijk<NxNyNz; ijk++){
			psum = 1.0-phA[ijk]-phB[ijk];
			for (X=0;X<NF_N;X++) for (Y=0; Y<K_i[X]; Y++) WA[X][_Yijk] += wopt*(WAdiff[X][_Yijk]- wcmp*psum);
			wB[ijk] += wopt*(wBdiff[ijk]-wcmp*psum);
		}
		for (X=0;X<NF_N;X++){free(WAdiff[X]); free(WAnew[X]); free(DAnew[X]);}
		free(WAdiff); free(WAnew); free(DAnew);
		free(wBdiff); free(wBnew); free(DBnew);
		return and_err;
	}	
	//if (iter%10==0) printf("err: %.5e, amax: %.3e, bmax: %.3e\n", and_err, amax, bmax);
	and_Nr = fmin(iter-1, and_NrMax); 
	Cs = (double *)calloc(and_Nr, sizeof(double)); 

	//Store histories
	u = INIT_2D(u, and_Nr, and_Nr);
	u_temp = INIT_2D(u_temp, and_Nr, and_Nr);
	v = INIT(and_Nr);
	for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) for (ijk=0;ijk<NxNyNz;ijk++){
		for (m=1;m<=and_Nr;m++) {
			for (n=1;n<=and_Nr;n++){
				u[m-1][n-1] += (DAs[end][X][_Yijk] - DAs[end-m][X][_Yijk]) * (DAs[end][X][_Yijk] - DAs[end-n][X][_Yijk]);
			}	
			v[m-1] += (DAs[end][X][_Yijk] - DAs[end-m][X][_Yijk])*DAs[end][X][_Yijk];
		}
	}	
	for (ijk=0;ijk<NxNyNz;ijk++) for (m=1;m<=and_Nr;m++){
 		for (n=1;n<=and_Nr;n++) u[m-1][n-1] += (DBs[end][ijk]-DBs[end-m][ijk])*(DBs[end][ijk]-DBs[end-n][ijk]);
		v[m-1] += (DBs[end][ijk]-DBs[end-m][ijk])*DBs[end][ijk];
	}

	// Compute Coeffs	
	inv(and_Nr, u, u_temp); //u ^-1	
	for (m=0;m<and_Nr;m++) for (n=0;n<and_Nr;n++) u[n][m] = u_temp[m][n]; //uT^-1
	mult(and_Nr, u, v, Cs); //CN = uT^-1 * v
printf("CAS: ");
for (n=0;n<and_Nr;n++) printf("%.6e ", Cs[n]);
printf("\n");

	//Half step
	for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) for (ijk=0;ijk<NxNyNz;ijk++){
		WAnew[X][_Yijk] = WAs[end][X][_Yijk]; //k+1/2
		DAnew[X][_Yijk] = DAs[end][X][_Yijk];
		for (n=1;n<=and_Nr;n++) {
			WAnew[X][_Yijk] += wand*Cs[end-n]*(WAs[end-n][X][_Yijk]-WAs[end][X][_Yijk]);
			DAnew[X][_Yijk] += wand*Cs[end-n]*(DAs[end-n][X][_Yijk]-DAs[end][X][_Yijk]);
		}
	}
	for (ijk=0;ijk<NxNyNz;ijk++){
		wBnew[ijk] = WBs[end][ijk]; //k+1/2
		DBnew[ijk] = DBs[end][ijk];
		for (n=1; n<=and_Nr; n++){
			wBnew[ijk] += wand*Cs[end-n]*(WBs[end-n][ijk] - WBs[end][ijk]);
			DBnew[ijk] += wand*Cs[end-n]*(DBs[end-n][ijk] - DBs[end][ijk]);
		}
	}	
	//Update fields
	lambda = 1.0; //1.0-pow(0.9, iter/200);
	for (ijk=0;ijk<NxNyNz;ijk++){
		for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++){
			WA[X][_Yijk] = WAnew[X][_Yijk] + lambda * DAnew[X][_Yijk];
		}
		wB[ijk] = wBnew[ijk] + lambda * DBnew[ijk];
	}

	//Free
	for (X=0;X<NF_N;X++){free(WAdiff[X]); free(WAnew[X]); free(DAnew[X]);}
	free(WAdiff); free(WAnew); free(DAnew);
	free(wBdiff); free(wBnew); free(DBnew);
	free(Cs);
	for (n=0;n<and_Nr;n++) {free(u[n]); free(u_temp[n]);}
	free(u); free(u_temp);
	free(v);
	return and_err;
}

void solve(int ndim, double **a, double *b, double *x){
	double *bp; bp = INIT(ndim);
	int i,j;
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

void inv(int ndim, double **A, double **invA){
//Inverse using Doolittle LU Factorization
//Adapted from Hoffman
//Output inverse of A into invA
	int i,j,k;
	double em;
	double **l, **u;
	double *b, *x, *col; 
	b = INIT(ndim); x = INIT(ndim); col = INIT(ndim);
	u = INIT_2D(u, ndim, ndim);
	l = INIT_2D(l, ndim, ndim);
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
