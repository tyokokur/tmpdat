// Import
#include "common.h"
#include "consts.h"

// Input
#include "input.h"

// Functions
#include "setup.h"
#include "report.h"
#include "FFT_solve.h"
#include "MDE.h"
#include "gets.h"

int main(void) {
	SETUP();
	FFTW_SETUP();
	MDE_SETUP();

	do { 
		iter += 1;

		MDE_STEP();
		get_freeE();
		and_err = And_mix(WA, wB);

		if (iter % 10  == 0) {
			printout();
			write_it();
		}
		if (iter % 100 == 0) {
			write_ph();
			write_W();
			write_elec();
		}

	} while(iter<maxIter&&(and_err>1e-3||inCompMax>Sm1||freeDiff>Sm2||iter<123));

	printout();
	write_it();
	write_ph();
	write_W();
	write_elec();

	MDE_CLEAN();
	FFTW_CLEAN();
	CLEAN();

	return 0;
}
