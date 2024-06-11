#include "MDE.h"
#include "FFT_solve.h"
#include "gets.h"

#ifndef _IJSMACRO
	#define _IJSMACRO
	#define _IJS(ijk,Y,s) (Y*(Nm_i[k][Y]+1) + s)*NxNyNz_8 + ijk
#endif

/* GLOBALS */
//MDE_STEP_INIT();
	double *wA;
	double Q2;
	int *Ns;
//MDE_SETUP();
	double *qInt, **qA, **qcA, **QA, **QcA; 
	double *int_Q1_Quen, *deltaFun;

void MDE_STEP(void) {
	int X;
	for (X=0; X<NF_N; X++) { 
		MDE_STEP_INIT(X);

		qInt_Free();

int nthreads;
#pragma omp parallel
{ nthreads = omp_get_num_threads(); }
double tic = omp_get_wtime();
		MDE_FFT_BACKWARD(X);	
double toc = omp_get_wtime();
//printf("\t[MDE.c/MDE_FFT_BACKWARD] (%d threads) %f\n", nthreads, toc-tic);

		qInt_Graft(X, qcA);
tic = omp_get_wtime();
		MDE_FFT_FORWARD(X);
toc = omp_get_wtime();
//printf("\t[MDE.c/MDE_FFT_FORWARD] (%d threads) %f\n", nthreads, toc-tic);

		get_PHI();

		MDE_STEP_CLEAN();
	}
}

void MDE_STEP_INIT(int x) {
	long int ijk;
	int Y; 

	wA = INIT(NxNyNz*K_i[x]);
	Ns = (int *)calloc(Ns_i[x][K_i[x]-1], sizeof(int));
	for (Y=0; Y<K_i[x]; Y++) {
		for (ijk=0; ijk<NxNyNz; ijk++) wA[_Yijk] = WA[x][_Yijk];
		Ns[Y] = Ns_i[x][Y];
	}
	for (ijk=0; ijk<NxNyNz; ijk++) {
		qA[ijk]  = INIT(Ns[K_i[x]-1]+1);
		qcA[ijk] = INIT(Ns[K_i[x]-1]+1);
	}
}

void MDE_FFT_FORWARD(int X) {
	long int ijk;
	int Y;
	int s, s0;
	// Solve MDE
	//if (HALF_BOOL) MDE_COS(qInt, qA, K_i[X], wA, Ns, 1, fftw_inc, fftw_outc_r); 
	if (HALF_BOOL) MDE_COS_CUTOFF(qInt, qA, K_i[X], wA, Ns, 1, fftw_inc_cutoff, fftw_outc_r_cutoff); 
	else           MDE_FFT(qInt, qA, K_i[X], wA, Ns, 1, fftw_in,  fftw_out_c);

	// Get QA
	for (Y=0; Y<K_i[X]; Y++) {
		if (Y==0) s0 = 0;
		else s0 = Ns_i[X][Y-1];

		for (ijk=0; ijk<NxNyNz; ijk++) for (s=s0; s<=Ns_i[X][Y]; s++) {
			QA[X][_IJS(ijk,Y,s)] = qA[ijk][s];
		}
	}
}

void MDE_FFT_BACKWARD(int X) {
	long int ijk;
	int s, s0;
	int Y;
	// Solve MDE
	if (HALF_BOOL) MDE_COS_CUTOFF(qInt, qcA, K_i[X], wA, Ns, -1, fftw_inc_cutoff, fftw_outc_r_cutoff); 
	//if (HALF_BOOL) MDE_COS(qInt, qcA, K_i[X], wA, Ns, -1, fftw_inc, fftw_outc_r); 
	else           MDE_FFT(qInt, qcA, K_i[X], wA, Ns, -1, fftw_in, fftw_out_c);
	// Get QcA
	for (Y=0; Y<K_i[X]; Y++) {
		if (Y==0) s0 = 0;
		else s0 = Ns_i[X][Y-1];

		for (ijk=0; ijk<NxNyNz; ijk++) for (s=s0; s<=Ns_i[X][Y]; s++) {
			QcA[X][_IJS(ijk,Y,s)] = qcA[ijk][s];
		}
	}
}

void MDE_STEP_CLEAN(void) {
	long int ijk;
	for (ijk=0; ijk<NxNyNz; ijk++) { free(qA[ijk]); free(qcA[ijk]); }
	free(wA); free(Ns);
}

void MDE_SETUP(void) {
	qInt = INIT(NxNyNz);
	qA   = (double **)calloc(NxNyNz, sizeof(double));
	qcA  = (double **)calloc(NxNyNz, sizeof(double));
	int_Q1_Quen = INIT(NF_N);
	
	int i, k, X, Y;
	printf("Ns0: \n");
	for (X=0;X<NF_N;X++) {
		printf("\t%d : [", X);
		for (Y=0;Y<K_i[X]-1;Y++) {
			Ns_i[X][Y] = roundl(Nm_i[X][Y] / ds0); 
			printf("%d ",Ns_i[X][Y]);
		}
		Ns_i[X][Y] = roundl(Nm_i[X][Y] / ds0); 
		printf("%d]\n",Ns_i[X][Y]);
	}

	// Delta function
	deltaFun = INIT(Nz);

	//// Gaussian (Chantawansri Fredrickson 2011)
	double gamma = 0.20; // variance, recommend > 0.10 based on accuracy on of phi_p quadrature
	double sum=0.0;
	for (k=0; k<Nz; k++) {
		deltaFun[k] = 1.0/sqrt(M_PI*gamma/2.0) * exp(-0.5 * (k*dz-dz)*(k*dz-dz) / gamma);
		//if (k!=0) sum += deltaFun[k]; // int from z* to inf
		sum += deltaFun[k]; 
		if (k>= Nz-2) deltaFun[k] = 0.0; // Force neumann at opposite boundary
	}
	sum *= dz;
	for (k=0;k<Nz;k++) deltaFun[k] /= sum; // Ensure int = 1
	
	/*
	double csum=0.0;
	for (k=0;k<Nz;k++) csum += deltaFun[k];
	CHECK_VAL_N(csum*dz); exit(1);
	*/
	

	//// Disc
	//for (k=0; k<Nz; k++) deltaFun[k] = 0.0;
	//deltaFun[1] = 1.0/dz;
	//printf("TEMP deltaFun[1]!!!!\n");
	//deltaFun[1] = 1e-35; 

	QA = INIT_2D(QA, NF_N, K_i[i] * (Ns_i[i][K_i[i]-1]+1) * NxNyNz);
	QcA= INIT_2D(QcA,NF_N, K_i[i] * (Ns_i[i][K_i[i]-1]+1) * NxNyNz);
}

void MDE_CLEAN(void){
	free(qInt); free(qA); free(qcA);
	free(int_Q1_Quen); free(QA); free(QcA);
	free(deltaFun);
}

void qInt_Graft(int X, double **qstar){
	long int ij, ijk, ij1;
	int k;

	for (ij=0; ij<NxNy; ij++) {
		//qInt[ij*Nz+0] = 0.0; //DBC
		//for (k=1; k<Nz-1; k++) { //DBC
		for (k=0; k<Nz; k++) {
			ijk = ij*Nz + k;
			ij1 = ij*Nz + 1;
			qInt[ijk] = sigma_i[X]*v0 / qstar[ij1][0] * deltaFun[k];
			//if ((ij==0) && (k<10)) {CHECK_VAL(deltaFun[k]); CHECK_VAL_N(qInt[ijk]);}
		}
		//qInt[ij*Nz+Nz-1] = 0.0; //DBC
	}
//exit(1);
}

void qInt_Free(void){
	// NBC
	long int ij, ijk; 
	for (ijk=0; ijk<NxNyNz; ijk++) qInt[ijk] = 1.0;	
//	for (ij=0; ij<NxNy; ij++) qInt[ij*Nz+0]  = 0.0;
}


