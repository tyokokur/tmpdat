/* input.h */
// Parameters that change from code to code //

// Setup
int init_opt   = -1;
int HALF_BOOL  = 1;
int maxIter    = 10000;

// System
double lx = 12.0;

// Block Polymer 
double b0 = 1.00;
double v0 = 4.19; // Add decimal points

// Solution
double c0_plus  = 10e-03;
double c0_minus = 10e-03;

// Initialization
double xCmax   = 0.80;
double xCneck  = 10;

// Convergence
double wopt = 1e-03;
double wcmp = 30;

double Sm1 = 2e-06;
double Sm2 = 1e-09;

double wopt_PB = 1e-03;
double Sm_PB   = 1e-10;

int and_it    = 1000;
double wand   = 1e-09;
int and_NrMax = 10;

// Files
char seq_file[30] = "H300.txt";
char Win[30] = "W_in.dat";
