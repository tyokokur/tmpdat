#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>
#include <string.h>

#define INIT(size) (double *)calloc(size, sizeof(double))
#define INIT_2D(NAME, size) (double **)calloc(size, sizeof(double)); for (i=0; i<size; i++) NAME[i] = INIT(size)
#define CHECK_VAL_N(VAR) printf(#VAR": %.6e\n", VAR)
#define CHECK_VAL(VAR) printf("%.6e ", VAR)

void heat_NN(double *x, double *f, double s1);
void thomas(int n, double **a, double *b, double *x);
void write(double *x);
void writep(double *x);
void write_cat(double *x);
void past(void);

double lz = 100.0;
size_t N = 1500;
double T = 5.0;
size_t Nt= 10000;

size_t i,j,step;
size_t MaxIT=1000;
double h, d0, d1, s0, s1, h2, k, wopt=0.25;
double sum, diff, err, err_NL, err_max = 1e-09, err_max_NL = 1e-09;
double *x, *b; 

int main(void){
	printf("N: %d, Nt: %d\n", N, Nt);
	x = INIT(N); b = INIT(N);
	h = lz / (N-1); h2=h*h;
	k = T / (Nt-1);

	int Z_plus = 1;
	double eps_P = 1.0;
	double sigma_cat = 10.0;
	s1 = -sigma_cat / eps_P;

	//for (i=0;i<N;i++) x[i] = sigma_cat/lz; //Initial condition
	for (i=0;i<N;i++) {
		sum  = 2.0*log(1.0+i*h/(2.0/sigma_cat));
		x[N-1-i] = sum; //Initial condition
	}
	sum = x[0];
	for(i=0;i<N;i++) x[i] -= sum;
	write(x);

	size_t it, t;
	double *nG0; nG0 = INIT(N);
	double *rho_elec_polym; rho_elec_polym = INIT(N);
	double **A; A = (double **)calloc(N, sizeof(double)); for(i=0;i<N;i++) A[i] = (double *)calloc(3, sizeof(double));
	double *b; b = INIT(N); 
	double *del; del = INIT(N);
	double *x_old; x_old = INIT(N);
	double *xp; xp = INIT(N);
	double *xp_old; xp_old = INIT(N);
	double M0, M1, N1, f, sum;

	double kh = k/h, kh2 = k/h2, kh2_2 = kh2/2.0;
	for (t=1; t<=Nt; t++){
		for (i=0;i<N;i++) { x_old[i] = x[i]; xp_old[i] = xp[i]; } // N, N+1

		for (i=1;i<N-1;i++) nG0[i] = x_old[i] + eps_P*kh2_2*(x_old[i+1] - 2.0*x_old[i] + x_old[i-1]) + k*rho_elec_polym[i];
		nG0[0] = x_old[0] + eps_P*kh2*(x_old[1] - x_old[0]) + k*rho_elec_polym[0];
		nG0[N-1] = x_old[N-1] + eps_P*kh2*(x_old[N-2] - x_old[N-1] + 2.0*h*s1) + k*rho_elec_polym[N-1];
		for (M0=0.0,i=0;i<N;i++) M0 += exp(-Z_plus*x_old[i]);

		it = 0;
		do {
			it += 1;
			for (M1=0.0,i=0;i<N;i++) M1 += exp(-Z_plus*x[i]);
			f = Z_plus*sigma_cat*kh/(M1 + M0);

			// Interior
			for (i=1;i<N;i++) { 
				N1 = exp(-Z_plus*x[i]) + exp(-Z_plus*x_old[i]);

				A[i][0] = -eps_P*kh2_2;
				A[i][1] = 1.0 + eps_P*kh2 + Z_plus*Z_plus*sigma_cat*kh * N1/M1 * (1 - exp(-Z_plus*x[i])/M1);
				A[i][2] = -eps_P*kh2_2;

				b[i] = -x[i] + eps_P*kh2_2*(x[i+1]-2.0*x[i]+x[i-1]) + f*N1 + nG0[i];
			}

			// Left boundary
			N1 = exp(-Z_plus*x[0]) + exp(-Z_plus*x_old[0]);
			A[0][1] = 1.0 + eps_P*kh2 + Z_plus*Z_plus*sigma_cat*kh * N1/M1 * (1 - exp(-Z_plus*x[0])/M1);
			A[0][2] = -eps_P*kh2_2;
			b[0] = -x[0] + eps_P*kh2*(x[1]-x[0]) + f*N1 + nG0[0];

			// Right boundary
			A[N-1][0] = -eps_P*kh2_2;
			N1 = exp(-Z_plus*x[N-1]) + exp(-Z_plus*x_old[N-1]);
			A[N-1][1] = 1.0 + eps_P*kh2 + Z_plus*Z_plus*sigma_cat*kh * N1/M1 * (1 - exp(-Z_plus*x[1])/M1);
			b[N-1] = -x[N-1] + eps_P*kh2*(x[N-2]-x[N-1]) + f*N1 + nG0[N-1];

			thomas(N, A, b, del);

			for (sum=0.0,i=0;i<N;i++) { x[i] += del[i]; sum += del[i] * del[i]; }
			err_NL = sqrt(sum);

			if (it%100==0) printf("\tit_NL: %d, err_NL: %10.5e\n", it, err_NL);
		} while(err_NL > err_max_NL && it < MaxIT);

		if (t%100==0) {
			printf("t: %d (%10.5e), err: %10.5e\n", t, t*k, err);
			//write_cat(x);
			writep(xp);
		}

		xp[0] = (-3.0*x[0] + 4.0*x[1] - x[2]) / 2.0/h;
		for (i=1;i<N-1;i++) xp[i] = (x[i+1] - x[i-1]) /2.0/h;
		xp[N-1] = (x[N-3] - 4.0*x[N-2] + 3.0*x[N-1])/ 2.0/h;

		for (sum=0.0,i=0;i<N;i++) sum += (xp[i] - xp_old[i]) * (xp[i] - xp_old[i]);
		err = sqrt(sum);

		if (err < err_max) break;
	}

	printf("DONE. t: %d (%10.5e), err: %10.5e\n", t, t*k, err);
	f = x[0];
	for (i=0;i<N;i++) x[i] -= f; //shift x[0] = 0
	write_cat(x);
	writep(xp);

	return 0;
}

void nonlin_heat(double *x, double *f, double s1){
	//Didn't properly check if worked
	// s1 = -1.0 for compatibility
	double *x_old; x_old = INIT(N);
	double *del; del = INIT(N);
	double *G0; G0 = INIT(N);
	double kh = k/h, kh2 = k/h2, kh2_2 = kh2/2.0; 
	double M;
	double **A; A = (double **)calloc(N, sizeof(double)); for(i=0;i<N;i++) A[i] = (double *)calloc(3, sizeof(double));

	size_t it, t;
	
	for (t=1;t<=Nt;t++){
		for(i=0;i<N;i++) x_old[i] = x[i]; // x_old = x(N), x = x(N+1)

		for(i=1;i<N-1;i++) G0[i] = -x_old[i] - kh2_2*(x_old[i+1]-2.0*x_old[i]+x_old[i-1]); //Interior
		G0[0]   = -x_old[0] - kh2*(x_old[1]-x_old[0]);
		G0[N-1] = -x_old[N-1] - kh2*(x_old[N-2]-x_old[N-1]+2.0*h*s1);

		it = 0;
		do { 
			it += 1;

			for(M=0,i=0;i<N;i++) M += x[i]*x[i] + x_old[i]*x_old[i]; // u^2

			for(i=1;i<N-1;i++){
				A[i][0] = -kh2_2;
				A[i][1] = 1.0 + kh2 - 2.0*kh*x[i]/M + 2.0*kh*(x[i]*x[i]+x_old[i]*x_old[i])*x[i]/M/M;
				A[i][2] = -kh2_2;

				b[i] = -x[i] + kh2_2*(x[i+1]-2.0*x[i]+x[i-1]) + kh*(x[i]*x[i]+x_old[i]*x_old[i])/M - G0[i];
			}
			A[0][1] = 1.0 + kh2 - 2.0*kh*x[0]/M + 2.0*kh*(x[0]*x[0]+x_old[0]*x_old[0])*x[0]/M/M;
			A[0][2] = -kh2_2;
			b[0] = -x[0] + kh2*(x[1]-x[0]) + kh*(x[0]*x[0]+x_old[0]*x_old[0])/M - G0[0];

			A[N-1][0] = -kh2_2;
			A[N-1][1] = 1.0 + kh2 - 2.0*kh*x[N-1]/M + 2.0*kh*(x[N-1]*x[N-1]+x_old[N-1]*x_old[N-1])*x[N-1]/M/M;
			b[N-1] = -x[N-1] + kh2*(x[N-2]-x[N-1]) + kh*(x[N-1]*x[N-1]+x_old[N-1]*x_old[N-1])/M - G0[N-1];

			thomas(N, A, b, del);

			for (sum=0.0,i=0;i<N;i++) { x[i] += del[i]; sum += del[i]*del[i]; }
			err_NL = sqrt(sum);

		} while(it<MaxIT && err_NL>err_max_NL);


		for (sum=0.0,i=0;i<N;i++) sum += (x[i] - x_old[i])*(x[i] - x_old[i]);
		err = sqrt(sum);
		if (t%100==0) {
			printf("t: %d, err: %10.5e, it_NL: %d, err_NL: %10.5e\n", t, err, it, err_NL);
			write_cat(x);
			writep(x);
		}

		if (err < err_max) break;
	}
	printf("DONE: t = %d, err %10.5e\n", t, err);
	write_cat(x);

	for(i=0;i<N;i++) free(A[i]); free(A);
	free(x_old); free(del); free(G0);
}

void heat_NN(double *x, double *f, double s1){
	size_t i,t;
	double *x_old; x_old = INIT(N);
	double **A; A = (double **)calloc(N, sizeof(double)); for(i=0;i<N;i++) A[i] = (double *)calloc(3, sizeof(double));

	double kh2 = k/h2, kh2_2 = kh2/2.0;
	double sum;
	for (t=1;t<=Nt;t++){
		//Interior
		for (i=1;i<N-1;i++){
			A[i][0] = -kh2_2;
			A[i][1] = 1+kh2;
			A[i][2] = -kh2_2;

			b[i] = x[i] + kh2_2 * (x[i+1] - 2.0*x[i] + x[i-1]) + k*f[i];
		}
		//Left boundary
		A[0][1] = 1+kh2;
		A[0][2] = -kh2;
		b[0] = x[0] + kh2*(x[1]-x[0]) + k*f[0];
		//Right boundary
		A[N-1][0] = -kh2;
		A[N-1][1] = 1+kh2;
		b[N-1] = x[N-1] + kh2*(x[N-2]-x[N-1]) + k*f[N-1] + 2.0*k*s1/h;

		for (i=0;i<N;i++) x_old[i] = x[i];
		thomas(N, A, b, x);
		write_cat(x);

		for (sum=0,i=0;i<N;i++) sum += (x[i] - x_old[i])*(x[i] - x_old[i]);
		err = sqrt(sum);
		if (err < err_max) break;
	}
	printf("k = %d (t = %6.3f): %10.5e\n", t, t*k, err);

	for(i=0;i<N;i++) free(A[i]);free(A);
	free(x_old);
}


void write_cat(double *x){
	FILE *fp;
	fp = fopen("ph_a.dat", "a");
	for(i=0;i<N;i++) fprintf(fp,"%10.5e %10.5e\n", h*i, x[i]);
	fclose(fp);
}
void write(double *x){
	FILE *fp;
	fp = fopen("ph_a.dat", "w");
	for(i=0;i<N;i++) fprintf(fp,"%10.5e %10.5e\n", h*i, x[i]);
	fclose(fp);
}
void writep(double *xp){
	FILE *fp;
	fp = fopen("ph_p.dat", "w");
	for(i=0;i<N;i++) fprintf(fp,"%10.3e %12.5e\n", h*i, xp[i]);
}

void thomas(int n, double **a, double *b, double *x)
{
    /*Thomas Algorithm for tridiagonal systems -- Hoffman 1.8.3*/
    //Input: 
        //n: number of rows in a
        //a: reduced nx3 matrix (changed in-place)
        //b: homogeneity vector from ax = b
    //Output:
        //a: nx3 simplified through Gaussain elim (changed in-place)
        //x: nx1 solution vector to ax = b (changed in-place)
    int i, f = n-1;
    double em; 

    //Forward Elimination
    for (i = 1; i < n; i++)
    {
        em = a[i][0]/a[i-1][1];
        a[i][0] = em;
        a[i][1] = a[i][1] - em*a[i-1][2];
        b[i] = b[i] - a[i][0]*b[i-1];
    }

    //Back Substitution
    x[f] = b[f]/a[f][1];
    for (i = f-1; i >= 0; i--)
    {
        x[i] = (b[i] - a[i][2]*x[i+1])/a[i][1];
    }
    return;
}

void past(void){
	double *f; f = INIT(N);
	double *in; in = INIT(N);
	double Q;

	// 1A: f = 0, in = const, homogeneous neumann (no pen)
	s1 = 0.0;
	for (i=0;i<N;i++) in[i] = 1.0; //Initial condition
	for (i=0;i<N;i++) f[i]  = 0.0; //Source
	for (i=0;i<N;i++) x[i] = in[i]; write(x);
	heat_NN(x, f, s1);

	// 1B: f = 0, in = lin, homogeneous neumann (no pen)
	s1 = 0.0;
	for (i=0;i<N;i++) in[i] = 0.5*i*h; //Initial condition
	for (i=0;i<N;i++) f[i]  = 0.0; //Source
	for (i=0;i<N;i++) x[i] = in[i]; write(x);
	heat_NN(x, f, s1);

	// 1C: f = 0, in = quad, homogeneous neumann (no pen)
	s1 = 0.0;
	for (i=0;i<N;i++) in[i] = 1.5*(i*h)*(i*h); //Initial condition
	for (i=0;i<N;i++) f[i]  = 0.0; //Source
	for (i=0;i<N;i++) x[i] = in[i]; write(x);
	heat_NN(x, f, s1);

	// 2A: in = const, compat cond: int(f) = -s1 = 1, f = const
	s1 = -1.0;
	for (i=0;i<N;i++) in[i] = 1.0; //Initial condition
	for (i=0;i<N;i++) f[i]  = 1.0/(N-1)/h; //Source
	for (i=0;i<N;i++) x[i] = in[i]; write(x);
	heat_NN(x, f, s1);

	// 2B: in = quad, compat cond: int(f) = -s1 = 1, f = const
	s1 = -1.0;
	for (i=0;i<N;i++) in[i] = 1.0*(i*h)*(i*h); //Initial condition
	for (i=0;i<N;i++) f[i]  = 1.0/(N-1)/h; //Source
	for (i=0;i<N;i++) x[i] = in[i]; write(x);
	heat_NN(x, f, s1);

	// 2C: in = exp, compat cod: int(f) = -s1 = 1, f = const
	s1 = -1.0;
	for (Q=0,i=0;i<N;i++) Q += h*exp(-2.0*i*h);
	for (i=0;i<N;i++) in[i] = exp(-2.0*i*h)/Q; //Initial condition
	for (i=0;i<N;i++) f[i]  = 1.0/(N-1)/h; //Source
	for (i=0;i<N;i++) x[i] = in[i]; write(x);
	heat_NN(x, f, s1);

	// 2D: in = gauss, compat cod: int(f) = -s1 = 1, f = const
	s1 = -1.0;
	for (Q=0,i=0;i<N;i++) Q += h*exp(-1.0*i*h*i*h/0.2);
	for (i=0;i<N;i++) in[i] = exp(-1.0*i*h*i*h/0.2)/Q; //Initial condition
	for (i=0;i<N;i++) f[i]  = 1.0/(N-1)/h; //Source
	for (i=0;i<N;i++) x[i] = in[i]; write(x);
	heat_NN(x, f, s1);

	free(f); free(in);
}
