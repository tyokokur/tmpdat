/* input.h */
// Parameters that change from code to code //

// System
double lx = 10.0;

// Block Polymer 
double b0 = 1.00;
double v0 = 4.19; // Add decimal points

// Solution
double c0_plus  = 10e-03;
double c0_minus = 10e-03;

// Numerical 
int init_opt   = 0;
int HALF_BOOL  = 1;
double xCmax   = 0.80;
double wopt_PB = 1e-03;
double Sm_PB   = 1e-10;

// Files
char seq_file[30] = "H300.txt";
char Win[30] = "W_in.dat";
