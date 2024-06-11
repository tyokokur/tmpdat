#include "DGC.h"
#include "FFT_solve.h"
#include "gets.h"

/* GLOBALS */
// DGC_STEP_INIT()
	double *wA;
	int *Ns;

// DGC_SETUP()
	double *qInt, **qA, **qcA, **QA, **QcA; 
	double *xGAUSSDIST, *yGAUSSDIST, *zGAUSSDIST;
	double *int_Q1_Quen, *deltaFun;
	int Nzg_2;
	

void DGC_STEP(void){
	int X;
	for (X=0; X<NF_N; X++) {
		
		DGC_STEP_INIT(X);
		qInt_Free();
		DGC_FFT_BACKWARD(X);
		qInt_Graft(X, qcA); 
		DGC_FFT_FORWARD(X);
		get_PHI();

//write_qA_s(1, "qAs0.dat", NxNyNz, Ns_i[0], qA);
//write_qA_s(1, "qcAs0.dat",NxNyNz, Ns_i[0], qcA);
//exit(1);
		DGC_STEP_CLEAN();
	}
}

void DGC_STEP_INIT(int x){
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

void DGC_FFT_FORWARD(int X) {
	long int ijk;
	int Y;
	int s, s0;

	if (HALF_BOOL) DGC_COS(qInt, qA, K_i[X], wA, Ns, 1, fftw_inc_xy, fftw_outc_r_xy, GAUSSPROP, FG_in, FG_out);
	else           {printf("ERROR: ONLY CODED FOR HALF_BOOL\n\n"); exit(1);}

	// Get QA
	for (Y=0; Y<K_i[X]; Y++) {
		if (Y==0) s0 = 0;
		else s0 = Ns_i[X][Y-1];

		for (ijk=0; ijk<NxNyNz; ijk++) for (s=s0; s<=Ns_i[X][Y]; s++) {
			QA[X][_IJS(ijk,Y,s)] = qA[ijk][s];
		}
	}
}

void DGC_FFT_BACKWARD(int X) {
	long int ijk;
	int s, s0;
	int Y;

	if (HALF_BOOL) DGC_COS(qInt, qcA, K_i[X], wA, Ns, -1, fftw_inc_xy, fftw_outc_r_xy, GAUSSPROP, FG_in, FG_out); 
	else           {printf("ERROR: ONLY CODED FOR HALF_BOOL\n\n"); exit(1);}

	// Get QcA
	for (Y=0; Y<K_i[X]; Y++) {
		if (Y==0) s0 = 0;
		else s0 = Ns_i[X][Y-1];

		for (ijk=0; ijk<NxNyNz; ijk++) for (s=s0; s<=Ns_i[X][Y]; s++) {
			QcA[X][_IJS(ijk,Y,s)] = qcA[ijk][s];
		}
	}
}

void DGC_STEP_CLEAN(void) {
	long int ijk;
	for (ijk=0; ijk<NxNyNz; ijk++) { free(qA[ijk]); free(qcA[ijk]); }
	free(wA); free(Ns);
}

void DGC_SETUP(void) {
	qInt = INIT(NxNyNz);
	qA   = (double **)calloc(NxNyNz, sizeof(double));
	qcA  = (double **)calloc(NxNyNz, sizeof(double));
	int_Q1_Quen = INIT(NF_N);
	
	int i, j, k, X, Y;
	printf("Ns0 (= Nm since DGC): \n");
	for (X=0;X<NF_N;X++) {
		printf("\t%d : [", X);
		for (Y=0;Y<K_i[X]-1;Y++) {
			Ns_i[X][Y] = Nm_i[X][Y];
			printf("%d ",Ns_i[X][Y]);
		}
		Ns_i[X][Y] = Nm_i[X][Y];
		printf("%d]\n",Ns_i[X][Y]);
	}

	// Delta function
	deltaFun = INIT(Nz);

	// Gaussian (Chantawansri Fredrickson 2011)
	double gamma = 0.10; // variance, recommend > 0.10 based on accuracy on of phi_p quadrature
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
	

	// Disc
	//for (k=0; k<Nz; k++) deltaFun[k] = 0.0;
	//deltaFun[1] = 1.0/dz;

	// GAUSSDIST
	//double norm = pow(2.0*M_PI*b0*b0/3.0, -0.5)*0.5;
	double norm = 1.0;
double dsum=0.0;
	xGAUSSDIST = INIT(Nx); yGAUSSDIST = INIT(Ny); zGAUSSDIST = INIT(Nz);
	for (i=0;i<Nx;i++) xGAUSSDIST[i] = exp( -3.0/(2*b0*b0) * i*dx*i*dx)*norm;
	for (j=0;j<Ny;j++) yGAUSSDIST[j] = exp( -3.0/(2*b0*b0) * j*dy*j*dy)*norm;
	for (k=0;k<Nz;k++){zGAUSSDIST[k] = exp( -3.0/(2*b0*b0) * k*dz*k*dz)*norm; dsum+=zGAUSSDIST[k];}
//printf("norm: %.5f, DSUM: %.6g\n", norm, dsum);

	Nzg_2 = roundl( 3.5*b0/sqrt(3.0) / dz ); //Number of grid points for 4*std of gaussian chain (for total of 8\sigma)
	printf("Nzg_2 = %d\n", Nzg_2); 
	
	// QAs
	QA = INIT_2D(QA, NF_N, K_i[i] * (Ns_i[i][K_i[i]-1]+1) * NxNyNz);
	QcA= INIT_2D(QcA,NF_N, K_i[i] * (Ns_i[i][K_i[i]-1]+1) * NxNyNz);
}

void DGC_CLEAN(void){
	free(qInt); free(qA); free(qcA);
	free(int_Q1_Quen); free(QA); free(QcA);
	free(deltaFun);
	free(GAUSSPROP); free(xGAUSSDIST); free(yGAUSSDIST); free(zGAUSSDIST);
	free(zGAUSSPROP);
}

void qInt_Graft(int X, double **qstar){
	long int ij, ijk, ij1;
	int k;

	for (ij=0; ij<NxNy; ij++) {
		//qInt[ij*Nz+0] = 0.0; //DBC
		//for (k=1; k<Nz-1; k++) { //DBC
		for (k=0; k<Nz; k++) {
			ijk = ij*Nz + k;
			//ij1 = ij*Nz + 1;
			//qInt[ijk] = sigma_i[X]*v0 / qstar[ij1][0] * deltaFun[k];
			qInt[ijk] = sigma_i[X]*v0 / qstar[ijk][0] * deltaFun[k] * exp(-wA[ijk]);
			//if ((ij==0) && (k<10)) {CHECK_VAL(deltaFun[k]); CHECK_VAL_N(qInt[ijk]);}
		}
		//qInt[ij*Nz+Nz-1] = 0.0; //DBC
	}
//exit(1);
}

void qInt_Free(void){
	long int ij, ijk; 
	for (ijk=0; ijk<NxNyNz; ijk++) qInt[ijk] = exp(-wA[ijk]); // DGC
	//for (ijk=0; ijk<NxNyNz; ijk++) qInt[ijk] = 1.0;	// NBC
	//for (ij=0; ij<NxNy; ij++) qInt[ij*Nz+0]  = 0.0; // DBC
}
