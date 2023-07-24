#include "FFT_solve.h"
#include "report.h"

void FFTW_SETUP(void){
	if (HALF_BOOL){
		printf("--- FFTW Halfspace ---\n");
		fftw_inc    = INIT(Nx);
		fftw_outc_r = INIT(Nx);
//		fftw_ins    = INIT(Nx-1);  
//		fftw_outs_r = INIT(Nx-1);
		cfp = fftw_plan_r2r_1d(Nx, fftw_inc, fftw_outc_r, FFTW_REDFT10, FFTW_MEASURE); 
		cip = fftw_plan_r2r_1d(Nx, fftw_outc_r, fftw_inc, FFTW_REDFT01, FFTW_MEASURE);
//		sfp = fftw_plan_r2r_1d(Nx-2, fftw_ins, fftw_outs_r, FFTW_RODFT00, FFTW_MEASURE); 
//		sip = fftw_plan_r2r_1d(Nx-2, fftw_outs_r, fftw_ins, FFTW_RODFT00, FFTW_MEASURE);
	}
	else {
		printf("--- FFTW Fullspace ---\n");
		fftw_in =  INIT(Nx2);
		fftw_out_c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nxc);
		p  = fftw_plan_dft_r2c_1d(Nx2, fftw_in, fftw_out_c, FFTW_ESTIMATE); 
		ip = fftw_plan_dft_c2r_1d(Nx2, fftw_out_c, fftw_in, FFTW_ESTIMATE);
	}
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
		for (i=0; i<Nx; i++) fftw_inc[i] = rpref[i] * g[i][s - sign];
		fftw_execute(cfp); // FFT forward (in --> out)
		for (i=0; i<Nx; i++) fftw_out[i] *= exp( -ds * b0*b0/6.0 * ksq[i]); // Im step 2
		fftw_execute(cip); // FFT back   (out --> in)
		// Return only first halfspace
		for (i=0; i<Nx; i++) g[i][s] = fftw_in[i] * rpref[i] / (Nx2-2); // Real step 3
	}
	free(rpref); free(ksq);
}

void MDE_SIN(double *init, double **g, int K, double *W, int *Ns, int sign, double *fftw_in, double *fftw_out){
// NOT VALIDATED
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

		for (i=0; i<Nx; i++){
			_mirror = (rho_elec_plus[i] - rho_elec_minus[i]  + rho_elec_polym[i]) / eps_prof[i] + c * pot_elec[i];
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

		for (i=0; i<Nx; i++) fftw_in[i] = (rho_elec_plus[i] - rho_elec_minus[i]  + rho_elec_polym[i]) / eps_prof[i] + c * pot_elec[i];
		fftw_execute(cfp); // FFT forward (in --> out)

		for (i=0; i<Nx; i++) fftw_out[i] /= (ksq[i] + c); 

		fftw_execute(cip); // FFT back   (out --> in)

		// Return only first halfspace
		for (i=0; i<Nx; i++) pot_elec[i] = fftw_in[i] / (Nx2-2); 

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

