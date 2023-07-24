#include "gets.h"

void get_PHI(void) {
	int i,j,k;

	// Polymer
	for (k=0; k<NF_N; k++) {
		Q1[k] = 0.0;
		for (i=0; i<Nx; i++) Q1[k] += QA[k][_IJS(i, (K_i[k]-1) , (Ns_i[k][K_i[k]-1]) )]; //integral of q(r, N)
		Q1[k] *= integ_cons;
		int_PHA_Quad(k);
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

void get_freeE(void) {
	double freeW, freeU, freeS, free_elec; 
	double free_elec_polym, free_elec_laplace, free_elec_ion;
	double freeU_step, free_elec_polym_step, freeOld;
	double psum, fpsum, eta1, eta2;
	int i,j,k;

	freeW=0.0; freeU=0.0; freeS=0.0; free_elec_polym=0.0; free_elec_laplace=0.0; free_elec_ion=0.0; inCompMax=0.0;

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
	
	freeU*=integ_cons; //CHECK_VAL_N(freeU);
	freeW*=integ_cons; //CHECK_VAL_N(freeW);

	free_elec_polym*=integ_cons*v0; //CHECK_VAL_N(free_elec_polym); //scaff check this
	free_elec_laplace*=integ_cons*v0; //CHECK_VAL_N(free_elec_laplace); //scaff check this
	free_elec_ion*=integ_cons*v0; //CHECK_VAL_N(free_elec_ion);
	free_elec_ion-=freeEnergy_bulk; //CHECK_VAL_N(free_elec_ion);
	free_elec=free_elec_polym+free_elec_laplace+free_elec_ion; //CHECK_VAL_N(free_elec); 

	freeS = -zs*Q2;
	for (k=0;k<NF_N;k++) freeS -= sigma_i[k]*log(Q1[k]);  //CHECK_VAL_N(freeS);

	freeOld = freeEnergy;
	freeEnergy=freeU+freeW+freeS+free_elec+freeEnergy_bulk; //CHECK_VAL_N(freeEnergy); //+free_factorial; 
	freeDiff = fabs((freeEnergy-freeOld)/freeEnergy);
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

void get_ELFIELD(void){
	int i, i1, i2;

	i=0;i2=i+1;
	field_elec[i]=field_elec[i2];

	for(i=1;i<Nx-1;i++) {
		i1=i-1; i2=i+1;
		field_elec[i]=0.5*(pot_elec[i1]-pot_elec[i2])/dx;
	}

	i=Nx-1;
	field_elec[i]=field_elec[i2];

	for(i=0;i<Nx;i++) field_elec_sq[i]=field_elec[i]*field_elec[i];
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

	int i,j,k, n, m;
	int and_Nr, end;

	WAdiff = INIT_2D(WAdiff, NF_N, Nx*K_i[i]);
	WAnew  = INIT_2D(WAnew , NF_N, Nx*K_i[i]);
	wBdiff = INIT(Nx); wBnew = INIT(Nx);
	DAnew  = INIT_2D(DAnew, NF_N, Nx*K_i[i]);
	DBnew  = INIT(Nx);

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
	u = INIT_2D(u, and_Nr, and_Nr);
	u_temp = INIT_2D(u_temp, and_Nr, and_Nr);
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
	lambda = 1.0; //1.0-pow(0.9, iter/200);
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
	free(Cs);
	for (n=0;n<and_Nr;n++) {free(u[n]); free(u_temp[n]);}
	free(u); free(u_temp);
	free(v);
	return and_err;
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
