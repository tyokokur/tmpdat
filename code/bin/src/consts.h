// Physical
double eps_vacuum = 8.854e-12; 
double kB         = 1.38e-23;
double Length     = 1e-09;
double Charge     = 1.6e-19;
double avo_num    = 6.022e23;

// System 
int Temp = 273;

// Solution
int Z_plus     = 1; 
int Z_minus    = 1;
double eps_r_P = 80.0; 
double eps_r_S = 80.0;
double mu_s    = -1.0;

//Numerical 
double dx    = 0.20;
double ds0   = 0.30;
int iter     = 0;
int iter_PB  = 0;
int MaxIT_PB = 0;
