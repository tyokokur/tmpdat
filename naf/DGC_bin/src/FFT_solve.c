#include "FFT_solve.h"
#include "report.h"
#include "gets.h"

/* GLOBALS */
//FFTW_SETUP()
	fftw_plan p, ip, sfp, sip, cfp, cip, cfp_xy, cip_xy, cfp_cutoff, cip_cutoff, cfp_FG, cip_FG;
	double *fftw_in, *fftw_inc, *fftw_inc_xy, *fftw_ins, *fftw_inc_cutoff;
	fftw_complex *fftw_out_c;
	double *fftw_outs_r, *fftw_outc_r, *fftw_outc_r_xy, *fftw_outc_r_cutoff;
	double *ksq, *ksq_xy, *ksq_cutoff;
	double *GAUSSPROP, *zGAUSSPROP, *FG_in, *FG_out;

void FFTW_SETUP(void){
	if (!fftw_init_threads()){
		printf("FFTW THREADS ERROR. Exiting...\n");
		exit(1);
	}
	printf("OMP_NUM_THREADS: %d\n", OMP_NUM_THREADS);
	fftw_plan_with_nthreads(OMP_NUM_THREADS);
	
	if (HALF_BOOL){
		printf("--- FFTW Halfspace ---\n");
		printf("Nx = Nx_real / 2, Ny = Ny_real / 2, Nz = Nz_real\n");
double tic = omp_get_wtime();
		// FFT xy only
		fftw_inc_xy    = INIT(NxNy);
		fftw_outc_r_xy = INIT(NxNy);
		cfp_xy = fftw_plan_r2r_2d(Nx, Ny, fftw_inc_xy, fftw_outc_r_xy, FFTW_REDFT10, FFTW_REDFT10, FFTW_PATIENT); 
		cip_xy = fftw_plan_r2r_2d(Nx, Ny, fftw_outc_r_xy, fftw_inc_xy, FFTW_REDFT01, FFTW_REDFT01, FFTW_PATIENT);

		// FFT full, Nz_cutoff
		/*fftw_inc_cutoff    = INIT(NxNyNz_cutoff);
		fftw_outc_r_cutoff = INIT(NxNyNz_cutoff);
		cfp_cutoff = fftw_plan_r2r_3d(Nx, Ny, Nz_cutoff, fftw_inc_cutoff, fftw_outc_r_cutoff, FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_PATIENT);
		cip_cutoff = fftw_plan_r2r_3d(Nx, Ny, Nz_cutoff, fftw_outc_r_cutoff, fftw_inc_cutoff, FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_PATIENT);*/

		// FFT full
		fftw_inc    = INIT(NxNyNz);
		fftw_outc_r = INIT(NxNyNz);
		cfp = fftw_plan_r2r_3d(Nx, Ny, Nz, fftw_inc, fftw_outc_r, FFTW_REDFT10, FFTW_REDFT10, FFTW_REDFT10, FFTW_PATIENT);
		cip = fftw_plan_r2r_3d(Nx, Ny, Nz, fftw_outc_r, fftw_inc, FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_PATIENT);

		// DGC 
		GAUSSPROP = INIT(NxNy); //EDIT CONVUL IN 2D
		zGAUSSPROP = INIT(Nz); //EDIT CONVUL IN 2D
		FG_in     = INIT(NxNy);
		FG_out    = INIT(NxNy);
		cfp_FG = fftw_plan_r2r_2d(Nx, Ny, GAUSSPROP, FG_out, FFTW_REDFT10, FFTW_REDFT10, FFTW_PATIENT);
		cip_FG = fftw_plan_r2r_2d(Nx, Ny, FG_out, FG_in, FFTW_REDFT01, FFTW_REDFT01, FFTW_PATIENT);

double toc = omp_get_wtime();
printf("FFTW PATIENT: %.3g sec\n", toc - tic);

		ksq = INIT(NxNyNz);
		ksq_xy = INIT(NxNy);
		ksq_cutoff = INIT(NxNyNz_cutoff);
		long int ijk, ij;
		int i, j, k;
		
		for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) for (k=0; k<Nz; k++) {  
			ijk = (i*Ny + j)*Nz + k;
			ksq[ijk] =  (i * 2.0*M_PI/(2.0*Nx*dx)) * (i * 2.0*M_PI/(2.0*Nx*dx)); // Fourier transform change of vars
			ksq[ijk] += (j * 2.0*M_PI/(2.0*Ny*dy)) * (j * 2.0*M_PI/(2.0*Ny*dy)); 
			ksq[ijk] += (k * 2.0*M_PI/(2.0*Nz*dz)) * (k * 2.0*M_PI/(2.0*Nz*dz)); 
		}
		/*
		for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) for (k=0; k<Nz_cutoff; k++) {  
			ijk = (i*Ny + j)*Nz + k;
			ksq_cutoff[ijk] =  (i * 2.0*M_PI/(2.0*Nx*dx)) * (i * 2.0*M_PI/(2.0*Nx*dx)); // Fourier transform change of vars
			ksq_cutoff[ijk] += (j * 2.0*M_PI/(2.0*Ny*dy)) * (j * 2.0*M_PI/(2.0*Ny*dy)); 
			ksq_cutoff[ijk] += (k * 2.0*M_PI/(2.0*Nz_cutoff*dz)) * (k * 2.0*M_PI/(2.0*Nz_cutoff*dz)); 
		}

		for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) {
			ij = i*Ny + j;
			ksq_xy[ij] =  (i * 2.0*M_PI/(2.0*Nx*dx)) * (i * 2.0*M_PI/(2.0*Nx*dx)); // Fourier transform change of vars
			ksq_xy[ij] += (j * 2.0*M_PI/(2.0*Ny*dy)) * (j * 2.0*M_PI/(2.0*Ny*dy)); 
		}*/
	}
	else {
		printf("--- FFTW Fullspace ---\n");
		printf("UNEDITED\n"); exit(1);
		fftw_in =  INIT(NxNy);
		fftw_out_c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NxNy);
		p  = fftw_plan_dft_r2c_1d(NxNy, fftw_in, fftw_out_c, FFTW_ESTIMATE); 
		ip = fftw_plan_dft_c2r_1d(NxNy, fftw_out_c, fftw_in, FFTW_ESTIMATE);
	}
}

/*void DGC_COS(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out, double *GAUSSPROP, double *FG_in, double *FG_out){
	int s, is, start, Nsf;
	long int ijk, ijkp;
	double sum;

	Nsf = Ns[K-1]; // Last monomer

	if (sign == 1)  start = 0;   // Forwards
	if (sign == -1) start = Nsf; // Backwards

	for (ijk=0; ijk<NxNyNz; ijk++) g[ijk][start] = init[ijk];

	for (is=1; is<=Nsf; is++){
		s = start + sign*is; // Contour length s = [1, Nsf]
		
		for (ijk=0; ijk<NxNyNz; ijk++){
			get_GAUSSPROP(ijk);
			sum = 0.0;
			for (ijkp=0; ijkp<NxNyNz; ijkp++) sum += GAUSSPROP[ijkp]*g[ijkp][s-sign];
			g[ijk][s] = sum*dx*dy*dz * exp(-W[ijk])/v0;
		}
	}
}*/
void DGC_COS(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out, double *GAUSSPROP, double *FG_in, double *FG_out){
	int s, is, start, Nsf, k, kp, fsum, bsum;
	long int ijk, ij;
	double *FG_out_perm;
	FG_out_perm = INIT(NxNy);

	Nsf = Ns[K-1]; // Last monomer

	if (sign == 1)  start = 0;   // Forwards
	if (sign == -1) start = Nsf; // Backwards

	for (ijk=0; ijk<NxNyNz; ijk++) g[ijk][start] = init[ijk];
//printf("DGC_COS init done\n");

	get_2D_GAUSSPROP(0, GAUSSPROP, xGAUSSDIST, Nx, yGAUSSDIST, Ny); //TODO: REDFT01 point of reflection? 

//printf("real gauss: ");
//for(ij=0;ij<NxNy;ij++) CHECK_VAL(GAUSSPROP[ij]);
//printf("\n\n");

	fftw_execute(cfp_FG); // GAUSSPROP -> FG_out
	for (ij=0;ij<NxNy;ij++) FG_out_perm[ij] = FG_out[ij];
//printf("FG recip (inv Gauss): ");
//for(ij=0;ij<NxNy;ij++) CHECK_VAL(FG_out[ij]);
//printf("\n\n");

//exit(1);
//
	for (is = 1; is <= Nsf; is++){
		s = start + sign*is; // Contour length s = [1, Nsf]

		//for (ij=0;ij<NxNy;ij++) g[ij*Nz+0][s] = 0.0; //TODO
		for (k=0; k<Nz; k++){
			get_1D_GAUSSPROP(k, zGAUSSPROP, zGAUSSDIST, Nz); 
/*if ((k==1) && (s==1) && (sign=1)){
ij=0;
printf("zGAUSSPROP[%d]: %.5e ", k, zGAUSSPROP[k]); 
for (kp=k+1, isum=0; kp<Nz && isum<Nzg_2; kp++, isum++){ CHECK_VAL(zGAUSSPROP[kp]); }
printf(" zGAUSSPROP f done\n");
for (kp=k-1, isum=0; kp>=0 && isum<Nzg_2; kp--, isum++){ CHECK_VAL(zGAUSSPROP[kp]); }
printf(" zGAUSSPROP b done\n");

printf("COMPARE:\n");
CHECK_VAL(exp(-3.0/2/b0/b0 * pow(k*dz - k*dz, 2.0))); 
for (kp=k+1, isum=0; kp<Nz && isum<Nzg_2; kp++, isum++){ CHECK_VAL(exp(-3.0/2/b0/b0 * pow(k*dz - kp*dz, 2.0))); }
printf(" f done\n");
for (kp=k-1, isum=0; kp>=0 && isum<Nzg_2; kp--, isum++){ CHECK_VAL(exp(-3.0/2/b0/b0 * pow(k*dz - kp*dz, 2.0))); }
printf(" b done\n");
exit(1);
}*/
			//for (ij=0; ij<NxNy; ij++) for (kp=0, fftw_in[ij]=0.0; kp<Nz; kp++) {
				//fftw_in[ij] += zGAUSSPROP[kp]*g[ij*Nz+kp][s-sign]*dz;
			//}
			
			for (ij=0; ij<NxNy; ij++){ 
				fftw_in[ij] = zGAUSSPROP[k]*g[ij*Nz+k][s-sign];
				for (kp=k+1, fsum=0; kp<Nz && fsum<Nzg_2; kp++, fsum++){
					fftw_in[ij] += zGAUSSPROP[kp]*g[ij*Nz+kp][s-sign];
				}
				for (kp=k-1, bsum=0; kp>=0 && bsum<Nzg_2; kp--, bsum++){
					fftw_in[ij] += zGAUSSPROP[kp]*g[ij*Nz+kp][s-sign];
				}
				fftw_in[ij] *= dz;
			} 

			fftw_execute(cfp_xy); // fftw_in --> fftw_out
			
			for (ij=0; ij<NxNy; ij++) fftw_out[ij] = fftw_out[ij]*FG_out_perm[ij];

			fftw_execute(cip_xy); // fftw_out --> fftw_in

			for (ij=0;ij<NxNy;ij++) g[ij*Nz+k][s] = exp(-W[ij*Nz+k])*fftw_in[ij]/(4.0*NxNy)/v0;
			//for (ij=0;ij<NxNy;ij++) g[ij*Nz+k][s] = exp(-W[ij*Nz+k])*fftw_in[ij]/v0;
		}

//printf("real q(s-1): ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(fftw_in[ijk]);
//printf("\n\n");
//printf("inv q: ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(fftw_out[ijk]);
//printf("\n\n");

//printf("FG real (after pointwise): ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(FG_in[ijk]);
//printf("\n\n");

//printf("g[ijk][s]: ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(g[ijk][s]);
//printf("\n\n");
//exit(1);
	}
	free(FG_out_perm);
}

/* DGC 3D FULL
void DGC_COS(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out, double *GAUSSPROP, double *FG_in, double *FG_out){
	int s, is, start, Nsf;
	long int ijk, ijkp;
	double *FG_out_perm;
	FG_out_perm = INIT(NxNyNz);

	Nsf = Ns[K-1]; // Last monomer

	if (sign == 1)  start = 0;   // Forwards
	if (sign == -1) start = Nsf; // Backwards

	for (ijk=0; ijk<NxNyNz; ijk++) g[ijk][start] = init[ijk];

	get_GAUSSPROP(0); //TODO: REDFT01 point of reflection? 
//printf("real gauss: ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(GAUSSPROP[ijk]);
//printf("\n\n");

	fftw_execute(cfp_FG);
	for (ijk=0;ijk<NxNyNz;ijk++) FG_out_perm[ijk] = FG_out[ijk];
//printf("FG recip (inv Gauss): ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(FG_out[ijk]);
//printf("\n\n");

//exit(1);
//
	//for (is = 1; is <= Nsf; is++){
	for (is = 1; is <= Nsf; is++){
		s = start + sign*is; // Contour length s = [1, Nsf]
		for (ijk=0; ijk<NxNyNz; ijk++) fftw_in[ijk] = g[ijk][s-sign];
//printf("real q(s-1): ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(fftw_in[ijk]);
//printf("\n\n");
		fftw_execute(cfp);
//printf("inv q: ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(fftw_out[ijk]);
//printf("\n\n");

		for (ijk=0;ijk<NxNyNz;ijk++) FG_out[ijk] = FG_out_perm[ijk]*fftw_out[ijk];
		fftw_execute(cip_FG);
//printf("FG real (after pointwise): ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(FG_in[ijk]);
//printf("\n\n");

		for (ijk=0;ijk<NxNyNz;ijk++) g[ijk][s] = exp(-W[ijk])*FG_in[ijk]/(8.0*NxNyNz)/v0;
//printf("g[ijk][s]: ");
//for(ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(g[ijk][s]);
//printf("\n\n");
//exit(1);
	}
	free(FG_out_perm);
}*/

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

void MDE_COS_CUTOFF(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out){
	double *rpref;
	double ds, gaussp;
	long int ijk, ij, ijkc;
	int s, is, i, j, k, Y, Nsf;
	int start, Yold=-1;


	rpref = INIT(NxNyNz);
	Nsf = Ns[K-1]; // Last monomer
//	ds = 1.0; // Discrete Gaussian chain
	ds = ds0; 

	if (sign == 1)  start = 0;   // Forwards
	if (sign == -1) start = Nsf; // Backwards

	for (ijk=0; ijk<NxNyNz; ijk++) g[ijk][start] = init[ijk];

//if(sign==1) for (ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(init[ijk]);

	for (is = 1; is <= Nsf; is++){
		s = start + sign*is; // Contour length s = [1, Nsf]

		// Piecewise W field based on is
		for (Y=0; Y<K; Y++) if (s <= Ns[Y]){ 
			if (Y != Yold) {
				#pragma omp parallel for shared(rpref, W, Y) private(ijk)
				for (ijk=0; ijk<NxNyNz; ijk++) rpref[ijk] = exp(-0.5 * ds * W[_Yijk]); 
				Yold = Y;
			}
			break; // Only do this once
		}

		for (ij=0; ij<NxNy; ij++) for (k=0; k < Nz_cutoff; k++) {
			ijk  = ij*Nz + k; // Non-cutoff index
			ijkc = ij*Nz_cutoff + k;
			fftw_in[ijkc] = g[ijk][s-sign] * rpref[ijk];
		}

		fftw_execute(cfp_cutoff);

		for (ijkc=0; ijkc<NxNyNz_cutoff; ijkc++) fftw_out[ijkc] *= exp( -1.0 * ds * b0*b0/6.0 * ksq_cutoff[ijkc]);

		fftw_execute(cip_cutoff);

		for (ij=0; ij<NxNy; ij++) for (k=0; k<Nz; k++) {
			if (k < Nz_cutoff) {
				ijk  = ij*Nz + k; // Non-cutoff index
				ijkc = ij*Nz_cutoff + k;
				g[ijk][s] = fftw_in[ijkc] * rpref[ijk] / (8.0 * NxNyNz_cutoff);
			}
			else    g[ijk][s] = 0.0;
		}

/*if (s==roundl(1) && (sign == 1)){
write_qA_s(s, "ph_a_QA.dat", NxNyNz, s, g);
exit(1);}*/

	}
	free(rpref);
}

void MDE_COS(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out){
	double *rpref;
	double ds, gaussp;
	long int ijk, ij;
	int s, is, i, j, k, Y, Nsf;
	int start, Yold=-1;

	rpref = INIT(NxNyNz);
	Nsf = Ns[K-1]; // Last monomer
//	ds = 1.0; // Discrete Gaussian chain
	ds = ds0; 

	if (sign == 1)  start = 0;   // Forwards
	if (sign == -1) start = Nsf; // Backwards

	for (ijk=0; ijk<NxNyNz; ijk++) g[ijk][start] = init[ijk];

//if(sign==1) for (ijk=0;ijk<NxNyNz;ijk++) CHECK_VAL(init[ijk]);

	for (is = 1; is <= Nsf; is++){
		s = start + sign*is; // Contour length s = [1, Nsf]

		// Piecewise W field based on is
		for (Y=0; Y<K; Y++) if (s <= Ns[Y]){ 
			if (Y != Yold) {
				#pragma omp parallel for shared(rpref, W, Y) private(ijk)
				for (ijk=0; ijk<NxNyNz; ijk++) rpref[ijk] = exp(-0.5 * ds * W[_Yijk]); 
				Yold = Y;
			}
			break; // Only do this once
		}

		for (ijk=0; ijk<NxNyNz; ijk++) fftw_in[ijk] = g[ijk][s-sign] * rpref[ijk];
		fftw_execute(cfp);
		for (ijk=0; ijk<NxNyNz; ijk++) fftw_out[ijk] *= exp( -1.0 * ds * b0*b0/6.0 * ksq[ijk] );
		fftw_execute(cip);
		for (ijk=0; ijk<NxNyNz; ijk++) g[ijk][s] = fftw_in[ijk] * rpref[ijk] / (8.0 * NxNyNz);

/*if (s==roundl(1) && (sign == 1)){
write_qA_s(s, "ph_a_QA.dat", NxNyNz, s, g);
exit(1);}*/

	}
	free(rpref);
}

void MDE_COSxy(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out){
	// xy FFT and z CFD
	double *rpref, *fhat_pref, **fhat;
	double ds, gaussp;
	long int ijk, ij;
	int s, is, i, j, k, Y, Nsf;
	int start, Yold = -1;

	rpref = INIT(NxNyNz);
	fhat_pref = INIT(NxNy);
	fhat  = INIT_2D(fhat, Nz, NxNy);
	Nsf = Ns[K-1]; // Last monomer
//	ds = 1.0; // Discrete Gaussian chain
	ds = ds0; // ASK CHAO WHY HIS IS 1.0

	gaussp = ds * b0*b0 / 6.0;
	#pragma omp parallel for shared(fhat_pref, ksq_xy) private(ij)
	for (ij=0; ij<NxNy; ij++) fhat_pref[ij] = exp(-gaussp*ksq_xy[ij]) - gaussp*2.0/dz/dz;

	if (sign == 1)  start = 0;   // Forwards
	if (sign == -1) start = Nsf; // Backwards

	#pragma omp parallel for shared(g, init, start) private(ijk)
	for (ijk=0; ijk<NxNyNz; ijk++) g[ijk][start] = init[ijk];

	for (is = 1; is <= Nsf; is++){
		s = start + sign*is; // Contour length s = [1, Nsf]

		// Piecewise W field based on is
		for (Y=0; Y<K; Y++) if (s <= Ns[Y]){ 
			if (Y != Yold) {
				#pragma omp parallel for shared(rpref, W, Y) private(ijk)
				for (ijk=0; ijk<NxNyNz; ijk++) rpref[ijk] = exp(-0.5 * ds * W[_Yijk]); 
				Yold = Y;
			}
			break; // Only do this once
		}
		// Init fhat
		for (k=0; k<Nz; k++) {
			#pragma omp parallel for shared(fftw_in, g, k, rpref, s, sign) private(ij, ijk)
			for (ij=0; ij<NxNy; ij++) {
				ijk = ij*Nz + k;
				fftw_in[ij] = rpref[ijk] * g[ijk][s - sign];
			}
double tic, toc;
if (is == 10 && k==10) tic = omp_get_wtime();
			fftw_execute(cfp_xy); // fftw_in[ij] --> fftw_out[ij]
if (is == 10 && k==10) {
	   toc = omp_get_wtime();
	   printf("fftw est time: %.3g\n", (toc - tic) * 2 * Nz * Nsf);
	}
			#pragma omp parallel for shared(fhat, fftw_out, k) private(ij)
			for (ij=0; ij<NxNy; ij++) fhat[k][ij] = fftw_out[ij]; // Store FFT results
		}
		
		#pragma omp parallel for shared(g, s) private(ij)
		for (ij=0; ij<NxNy; ij++){
			g[ij*Nz + 0][s]    = 0.0; // No pen wall @ z = 0 
			g[ij*Nz + Nz-1][s] = 0.0; // Far away DBC @ z = Lz
		}

		// Inner grid
		for (k=1; k<Nz-1; k++){ 
			#pragma omp parallel for shared(fftw_out, fhat, fhat_pref, k) private(ij, ijk)
			for (ij=0; ij<NxNy; ij++) {
				ijk = ij*Nz + k;
				fftw_out[ij] = fhat_pref[ij] * fhat[k][ij] + gaussp/dz/dz * (fhat[k+1][ij] - 2.0*fhat[k][ij] + fhat[k-1][ij]);
			}
			fftw_execute(cip_xy); // fftw_out[ij] --> fftw_in[ij]
			#pragma omp parallel for shared(fftw_in, g, k, rpref, s) private(ij, ijk)
			for (ij=0; ij<NxNy; ij++) {
				ijk = ij*Nz + k;
				//g[ijk][s] = rpref[ijk] * fftw_in[ij] / ( 2.0*(Nx-1) * 2.0*(Ny-1) ); // REDFT00 normalization
				g[ijk][s] = rpref[ijk] * fftw_in[ij] / ( 2.0*Nx * 2.0*Ny ); // REDFT10 normalization
			}
		}
	}

	for (i=0; i<Nz; i++) free(fhat[i]);
	free(fhat); free(fhat_pref);
	free(rpref);
}

void MDE_SIN(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out){
// NOT VALIDATED
	double *rpref, *ksq;
	double ds;
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


void PB_FFT(double *phA, double **PHA, double **PHA_T, double *fftw_in, fftw_complex *fftw_out){
	double *ksq, *pot_elec_old;
	double _mirror, c, err_check, err_max;
	int k, i, j;

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

		for (i=0; i<Nx; i++){
			_mirror = (rho_elec_plus[i] - rho_elec_minus[i]  - rho_elec_polym[i]) / eps_prof[i] + c * pot_elec[i];
			fftw_in[i]       = _mirror;
			fftw_in[Nx2-1-i] = _mirror;
		}

		fftw_execute(p); // FFT forward (in --> out)

		for (i=0; i<Nxc; i++) for(j=0;j<2;j++) fftw_out[i][j] /= (ksq[i] + c); 

		fftw_execute(ip); // FFT back   (out --> in)

		// Return only first halfspace
		for (i=0; i<Nx; i++) pot_elec[i] = fftw_in[i] / Nx2; 

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

	get_ELFIELD();

	free(ksq);
	free(pot_elec_old);
}

void PB_COS(double *phA, double **PHA, double **PHA_T, double *fftw_in, double *fftw_out){
	double *pot_elec_old;
	double c, err_check, err_max;
	int i1, i2, X, Y;
	long int ijk;

	pot_elec_old = INIT(NxNyNz);
	c = 1.0 / (L_Deby_S * L_Deby_S);

	// FFTW requires PBC; Calculate with z mirror and return only [0, Nz)
//printf("iter_PB: ");	
	// Simple mixing it solver
	iter_PB = 0;
	do { 
		iter_PB += 1;
//printf("%d ", iter_PB);

		#pragma omp parallel for \
			shared(rho_elec_plus, rho_elec_minus, rho_elec_polym, pot_elec, pot_elec_old) \
			private(ijk, X, Y)
		for(ijk=0;ijk<NxNyNz;ijk++){
			pot_elec_old[ijk]   = pot_elec[ijk];
			rho_elec_plus[ijk]  = Z_plus *fugac_plus *exp(-Z_plus*pot_elec[ijk]); 
			rho_elec_minus[ijk] = Z_minus*fugac_minus*exp( Z_minus*pot_elec[ijk]);
			rho_elec_polym[ijk] = 0.0;
			for (X=0;X<NF_N;X++) for (Y=0;Y<K_i[X];Y++) {
				#pragma omp atomic
				rho_elec_polym[ijk] += Alpha_i[X][Y]*PHA[X][_Yijk];
			}
			rho_elec_polym[ijk] *= 1.0/v0;
		}

		#pragma omp parallel for \
			shared(eps_prof, fftw_in, rho_elec_plus, rho_elec_minus, rho_elec_polym, pot_elec) \
			private(ijk)
		for (ijk=0; ijk<NxNyNz; ijk++) fftw_in[ijk] = (rho_elec_plus[ijk] - rho_elec_minus[ijk] - rho_elec_polym[ijk]) / eps_prof[ijk] + c * pot_elec[ijk];

		fftw_execute(cfp); // fftw_in[ijk] --> fftw_out[ijk]

		#pragma omp parallel for shared(fftw_out, ksq, c) private(ijk)
		for (ijk=0; ijk<NxNyNz; ijk++) fftw_out[ijk] /= (ksq[ijk] + c); 

		fftw_execute(cip); // fftw_out[ijk] --> fftw_in[ijk]

		// Return only first halfspace
		#pragma omp parallel for shared(fftw_in, pot_elec) private(ijk)
		for (ijk=0; ijk<NxNyNz; ijk++) pot_elec[ijk] = fftw_in[ijk] / (8.0*NxNyNz); 

		#pragma omp parallel for shared(pot_elec, pot_elec_old) private(ijk)
		for (ijk=0; ijk<NxNyNz; ijk++) pot_elec[ijk]=wopt_PB*pot_elec[ijk]+(1.0-wopt_PB)*pot_elec_old[ijk];  // simple linear superposition // 
		

		err_max=0.0;
		for(ijk=0;ijk<NxNyNz;ijk++)
		{
			err_check= (pot_elec[ijk]-pot_elec_old[ijk])*(pot_elec[ijk]-pot_elec_old[ijk]);
			if(err_check>err_max) err_max=err_check;
		}
		err_max=sqrt(err_max);

	} while(err_max>Sm_PB && iter_PB<MaxIT_PB);

	if (iter_PB == MaxIT_PB){ printf("Max_PB it %d. ", iter_PB); CHECK_VAL_N(err_max);}
	
	free(pot_elec_old);
}

void FFTW_CLEAN(void) {
	if (HALF_BOOL){
		fftw_destroy_plan(cfp); fftw_destroy_plan(cip);
		//fftw_destroy_plan(sfp); fftw_destroy_plan(sip);
	}
	else{
		fftw_destroy_plan(p);
		fftw_destroy_plan(ip);
	}
	free(fftw_in); 

	if (HALF_BOOL){
		free(fftw_inc);
		//free(fftw_ins);
		//free(fftw_outs_r);
		free(fftw_outc_r);
	}
	else fftw_free(fftw_out_c);
}

