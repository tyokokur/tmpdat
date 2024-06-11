// Import
#include "consts.h"

// Input
#include "../../input.h"

// Functions
#include "setup.h"
#include "report.h"
#include "FFT_solve.h"
#include "DGC.h"
#include "PB.h"
#include "gets.h"

int main(void) {
	int and_hold; double freeE_old;

	SETUP();
	FFTW_SETUP();
	DGC_SETUP();
	stdout_redir_refresh();

	do { 
		iter += 1;
		DGC_STEP();
		stdout_redir_refresh();
//if(iter==1){
//write_ph();
//exit(1);}
		freeE_old = freeE;

		if ((iter % 10 == 0)){// && (inCompMax < 0.1) ) {
			PB_STEP();		
			stdout_redir_refresh();

			freeE = get_freeE();
			stdout_redir_refresh();

			and_hold = and_NrMax; and_NrMax = 0; // Force superposition
			and_err = And_mix(WA, wB);
			and_NrMax = and_hold;
		}	
		else {
			ELEC_POLYM(); // Update rho_elec_polym if no PB_STEP

			freeE = get_freeE();
			stdout_redir_refresh();

			and_hold = and_NrMax; and_NrMax = 0; // Force superposition
			and_err = And_mix(WA, wB);
			and_NrMax = and_hold;
		}
		freeDiff = fabs((freeE-freeE_old)/freeE);

		if (iter % 1  == 0) {
			printout();
		}
		if (iter % 10 == 0) {
			write_it();
			write_W();
			write_ph();
			write_elec();
		}
		
		stdout_redir_refresh();

	} while(iter<maxIter&&(and_err>1e-3||inCompMax>Sm1||freeDiff>Sm2||iter<123));

	printout();
	write_it();
	write_ph();
	write_W();
	write_elec();
	write_finish();

/*
FILE *fp;
long int ijk; double sum=0;
fp = fopen("VOL.dat", "w");
for (ijk=0;ijk<NxNyNz;ijk++) sum+=PHA[0][ijk]*dz*dx*dy/lx/ly;
fprintf(fp, "PHA SUM                 : %10.5g\n", sum);
fprintf(fp, "EXPECTED sigma * Nm * v0: %10.5g\n", sigma_i[0]*Nm*v0);
fprintf(fp, "ERROR                   : %10.10f", sum/(sigma_i[0]*Nm*v0));
fclose(fp);
*/

	DGC_CLEAN();
	FFTW_CLEAN();
	CLEAN();

	return 0;
}
