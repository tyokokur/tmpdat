// Import
#include "common.h"
#include "consts.h"

// Input
#include "input.h"

// Functions
#include "main.h"
#include "setup.h"
#include "MDE.h"


int main(void) {
	SETUP();
	FFTW_SETUP();
	MDE_SETUP();

	do { 
		iter += 1;
printf("ITER: %d\n", iter);
		MDE_STEP();
printf("Post MDE_STEP: ");
CHECK_VAL_N(phA[1]);

double freeOld, freeDiff;

		freeOld=freeEnergy;
		get_freeE();
		freeDiff=fabs((freeEnergy-freeOld)/freeEnergy);

exit(1);

	} while(iter < 100);

	MDE_CLEAN();
	FFTW_CLEAN();
	CLEAN();

	return 0;
}

void get_freeE(void) {
	double freeW, freeU, freeS, free_elec; 
	double free_elec_polym, free_elec_laplace, free_elec_ion;
	double freeU_step, free_elec_polym_step;
	double psum, fpsum, eta1, eta2, inCompMax;
	int MAXMAX, i,j,k;

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
	
	freeU*=integ_cons; CHECK_VAL_N(freeU);
	freeW*=integ_cons; CHECK_VAL_N(freeW);

	free_elec_polym*=integ_cons*v0; CHECK_VAL_N(free_elec_polym); //scaff check this
	free_elec_laplace*=integ_cons*v0; CHECK_VAL_N(free_elec_laplace); //scaff check this
	free_elec_ion*=integ_cons*v0; //CHECK_VAL_N(free_elec_ion);
	free_elec_ion-=freeEnergy_bulk; CHECK_VAL_N(free_elec_ion);
	free_elec=free_elec_polym+free_elec_laplace+free_elec_ion; CHECK_VAL_N(free_elec); 

	freeS = -zs*Q2;
	for (k=0;k<NF_N;k++) freeS -= sigma_i[k]*log(Q1[k]);  CHECK_VAL_N(freeS);

	freeEnergy=freeU+freeW+freeS+free_elec+freeEnergy_bulk; CHECK_VAL_N(freeEnergy); //+free_factorial; 
}
