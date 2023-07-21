#include "MDE.h"
#include "FFT_solve.h"

void MDE_STEP(void) {
	int k;
	for (k=0; k<NF_N; k++) { 
		MDE_STEP_INIT(k);

		qInt_Gauss();
		MDE_FFT_FORWARD(k);

		qInt_Free();
		MDE_FFT_BACKWARD(k);	

		get_PHI();

		MDE_STEP_CLEAN();
	}
}

void MDE_STEP_INIT(int k) {
	int i, j; 
	wA = INIT(Nx*K_i[k]);
	Ns = (int *)calloc(Ns_i[k][K_i[k]-1], sizeof(int));
	for (j=0; j<K_i[k]; j++) {
		for (i=0; i<Nx; i++) wA[j*Nx + i] = WA[k][j*Nx+i];
		Ns[j] = Ns_i[k][j];
	}
	for (i=0; i<Nx; i++) {
		qA[i]  = INIT(Ns[K_i[k]-1]+1);
		qcA[i] = INIT(Ns[K_i[k]-1]+1);
	}
}

void MDE_FFT_FORWARD(int k) {
	int i, j, s, s0;
	double _temp;
	// Solve MDE
	if (HALF_BOOL) MDE_COS(qInt, qA, K_i[k], wA, Ns, 1, fftw_inc, fftw_outc_r); 
	else           MDE_FFT(qInt, qA, K_i[k], wA, Ns, 1, fftw_in, fftw_out_c);

	// Get QA
	for (j=0; j<K_i[k]; j++) {
		if (j==0) s0 = 0;
		else s0 = Ns_i[k][j-1];

		for (i=0; i<Nx; i++) for (s=s0; s<=Ns_i[k][j]; s++) {
			QA[k][_IJS(i,j,s)] = qA[i][s];
		}
	}
}

void MDE_FFT_BACKWARD(int k) {
	int i, j, s, s0;
	// Solve MDE
	if (HALF_BOOL) MDE_COS(qInt, qcA, K_i[k], wA, Ns, -1, fftw_inc, fftw_outc_r); 
	else           MDE_FFT(qInt, qcA, K_i[k], wA, Ns, -1, fftw_in, fftw_out_c);
	// Get QcA
	for (j=0; j<K_i[k]; j++) {
		if (j==0) s0 = 0;
		else s0 = Ns_i[k][j-1];

		for (i=0; i<Nx; i++) for (s=s0; s<=Ns_i[k][j]; s++) {
			QcA[k][_IJS(i,j,s)] = qcA[i][s];
		}
	}
}

void MDE_STEP_CLEAN(void) {
	int i;
	for (i=0; i<Nx; i++) { free(qA[i]); free(qcA[i]); }
	free(wA); free(Ns);
}

void get_PHI(void) {
	int i,j,k;

	// Polymer
	for (k=0; k<NF_N; k++) {
		Q1[k] = 0.0;
		for (i=0; i<Nx; i++) Q1[k] += QA[k][_IJS(i, (K_i[k]-1) , (Ns_i[k][K_i[k]-1]) )]; //integral of q(r, N)
		Q1[k] *= integ_cons;
		int_PHA_Reimm(k);
	}
	
	for (i=0; i<Nx; i++) {
		phA[i] = 0.0;
		for (k=0; k<NF_N; k++) phA[i] += PHA_T[k][i];
	}
	
	// Solvent
	Q2 = 0.0;
  	for(i=0;i<Nx;i++) {
  		phB[i] = zs*exp(-wB[i]);
  		Q2 += exp(-wB[i]);
  	}
  	Q2 *= integ_cons;
}

void int_PHA_Reimm(int k) {
	int i, j, s;

	double fflA=sigma_i[k]*ds0/Q1[k];
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

void int_PHA_Quad(int k) {
	// Fourth order Gauss Quad without endpoints (Chantawansri 2011 Fredrickson)
	int i,j,s;
	int ns;	
	
	double fflA=sigma_i[k]*ds0/Q1[k];
	for (i=0; i<Nx; i++) for (j=0; j<K_i[k]; j++){
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
		PHA[k][j*Nx+i] *= fflA;
	}
	
	for (i=0; i<Nx;i++){
		PHA_T[k][i] = 0;
		for (j=0;j<K_i[k];j++){
			PHA_T[k][i] += PHA[k][j*Nx+i];
		}
	}
}

void MDE_SETUP(void) {
	qInt = INIT(Nx);
	qA   = (double **)calloc(Nx, sizeof(double));
	qcA  = (double **)calloc(Nx, sizeof(double));
	Q1   = INIT(NF_N);
	
	int i, j, k;
	printf("Ns0: \n");
	for (k=0;k<NF_N;k++) {
		printf("\t%d : [", k);
		for (j=0;j<K_i[k]-1;j++) {
			Ns_i[k][j] = roundl(Nm_i[k][j] / ds0); 
			printf("%d ",Ns_i[k][j]);
		}
		Ns_i[k][j] = roundl(Nm_i[k][j] / ds0); 
		printf("%d]\n",Ns_i[k][j]);
	}

	QA = INIT_2D(QA, NF_N, K_i[i] * (Ns_i[i][K_i[i]-1]+1) * Nx);
	QcA= INIT_2D(QcA,NF_N, K_i[i] * (Ns_i[i][K_i[i]-1]+1) * Nx);
}

void MDE_CLEAN(void){
	free(qInt); free(qA); free(qcA);
	free(Q1); free(QA); free(QcA);
}

void qInt_Gauss(void){
	// Gauss approx for \delta(z-eps), NBC wall
	int i;
	double sigma = 1.0e-03; 
	for (i=0; i<Nx; i++) qInt[i] = 1.0/sqrt(2*M_PI * sigma) * exp(-0.5 * (i*dx)*(i*dx) / sigma);
}

void qInt_Disc(void){
	// Inpen wall, DBC
	int i;
	for (i=0; i<Nx; i++) qInt[i] = 0.0;
	qInt[1] = 1.0/dx;
}

void qInt_Free(void){
	// Inpen wall, DBC
	int i; 
	for (i=0; i<Nx; i++) qInt[i] = 1.0;	
	qInt[0] = 0.0;
}

