#include "report.h"

void REPORT_SETUP(void){
	snprintf(pattern, sizeof(char)*20, "_a%03.0fc%03.0f_x%03.0fs%03.0f", Alpha_i[0][0]*100, c0_plus*1000, Chi_i[0][0]*100, sigma_i[0]*100);
	snprintf(phname, sizeof(char)*30, "ph%s.dat", pattern);
	snprintf(phname_elec, sizeof(char)*30, "el%s.dat", pattern);
	snprintf(Wname, sizeof(char)*30, "W%s.dat", pattern);
	snprintf(itname, sizeof(char)*30, "it%s.dat", pattern);
}

void write_ph(double *phA, double **PHA_T, double **PHA, double *phB) {
	int i, j, k;
	FILE *fp=fopen(phname,"w");

	for(i=0;i<Nx;i++) 
	{
		fprintf(fp,"%10.3e ",rx[i]);
		fprintf(fp, "%10.4e ", phA[i]);
		for (k=0; k<NF_N; k++) {
			fprintf(fp, "%10.4e ", PHA_T[k][i]);	
			fprintf(fp, "[ ");
			for (j=0; j<K_i[k]-1; j++) fprintf(fp, "%10.4e ", PHA[k][j*Nx + i]);
			fprintf(fp, "%10.4e ] ", PHA[k][j*Nx + i]);
		}
		fprintf(fp, "%10.4e\n", phB[i]);
	}

	fclose(fp);
}

void write_elec(double *pot_elec,double *rho_elec_plus,double *rho_elec_minus,double *rho_elec_polym) {
	int i;

	FILE *fp=fopen(phname_elec,"w");

	for(i=0;i<Nx;i++)
	{	
		fprintf(fp,"%10.3e %10.4e %10.4e %10.4e %10.4e\n",rx[i],-pot_elec[i],rho_elec_plus[i],rho_elec_minus[i],rho_elec_polym[i]);
	}

	fclose(fp);
}

void write_W(double **WA,double *wB,double *eta) {
	int i,j,k;

	FILE *fp=fopen(Wname,"w");

	for(i=0;i<Nx;i++) 
	{
		for (k=0; k<NF_N; k++) {
			fprintf(fp, "[ ");
			for (j=0; j<K_i[k]-1; j++) fprintf(fp, "%10.4e ", WA[k][j*Nx + i]);
			fprintf(fp, "%10.4e ] ", WA[k][j*Nx + i]);
		}
		fprintf(fp, "%10.5e %10.5e %10.5e\n", wB[i], eta[i], pot_elec[i]);
	}

	fclose(fp);
}

void write_qA(char *name, int N, int Ns, double **qA) {
	int i,j,s;
	FILE *fp=fopen(name,"w");
	for (s=0; s<=Ns; s++) for(i=0;i<N;i++) {
		fprintf(fp, "%10.5e\n", qA[i][s]);
	}
	fclose(fp);
	printf("wrote qA\n");
}
