/* input.h */
// Parameters that change from code to code //

// Convergence
double wopt = 1e-02;
double wcmp = 1.5;

int and_it    = 1000;
double wand   = 1e-09;
int and_NrMax = 10;

double wopt_PB = 1e-02;
double Sm_PB   = 1e-10;
int MaxIT_PB   = 2000;

double Sm1 = 2e-06;
double Sm2 = 1e-09;

// Setup
int init_opt   =  0;
int win_box    =  0;
int win_nx     = 1;
int win_ny     = 1;
int win_nz     = 400;

int HALF_BOOL  = 1;
int maxIter    = 5000;

// Solution
double c0_plus  = 10e-03;
double c0_minus = 10e-03;

// System
double lx = 6.50;
double ly = 3.75;
double lz = 50.0;

// Block Polymer 
double b0 = 1.00;
double v0 = 4.19; // Add decimal points

// Initialization
double xCmax   = 0.60;
double xCneck  = 10;
double xCltot  = 0; 

// Files
char seq_file[30] = "ac.txt";
char Win[30] = "W_in.dat";
