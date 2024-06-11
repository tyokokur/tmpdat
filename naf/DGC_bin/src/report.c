#include "report.h"

/* GLOBALS */
int STDOUT_REDIR;
int errlog, outlog;
int errlog_save, outlog_save;
double c0_plus_n;

void report_vec(int N, double *vec){
	int i;
	FILE *fp = fopen("ph_a_EXPORT.dat", "w");
	for (i=0; i<N; i++) fprintf(fp, "%.5e\n", vec[i]);
	fclose(fp);
}

void stdout_redir_setup(void){
	if (getenv("STDOUT_REDIR")!=NULL) STDOUT_REDIR=atoi(getenv("STDOUT_REDIR"));
	else {
		printf("*** STDOUT_REDIR NOT FOUND; DEFAULT TO ZERO ***\n");
		STDOUT_REDIR = 0;
	}
}

int stdout_redir_start(void){
// source: https://www.unix.com/programming/268879-c-unix-how-redirect-stdout-file-c-code.html
	if (STDOUT_REDIR) {
		outlog = open(outlog_name, O_RDWR|O_CREAT|O_APPEND, 0600);
		if (-1 == outlog) { perror("opening outlog"); return 255; }

		errlog = open(errlog_name, O_RDWR|O_CREAT|O_APPEND, 0600);
		if (-1 == errlog) { perror("opening errlog"); return 255; }

		outlog_save = dup(fileno(stdout));
		errlog_save = dup(fileno(stderr));

		if (-1 == dup2(outlog, fileno(stdout))) { perror("cannot redirect stdout"); return 255; }
		if (-1 == dup2(errlog, fileno(stderr))) { perror("cannot redirect stderr"); return 255; }

		//printf("\n ... stdout start ...\n\n");
	}
}

void stdout_redir_refresh(void){
	//if (STDOUT_REDIR) printf("\n ... stdout refresh ...\n\n");
	stdout_redir_end();
	stdout_redir_start();
}

void stdout_redir_end(void){
// source: https://www.unix.com/programming/268879-c-unix-how-redirect-stdout-file-c-code.html
	if (STDOUT_REDIR) {
		//printf("\n ... stdout end ...\n\n");

		fflush(stdout); close(outlog);
		fflush(stderr); close(errlog);

		dup2(outlog_save, fileno(stdout));
		dup2(errlog_save, fileno(stderr));

		close(outlog_save);
		close(errlog_save);
	}
}

void printout(void) {
	printf(" %5d : %.8e, %.8e, %.8e, %.8e, %d\n",
		iter, freeE, freeDiff, inCompMax, and_err, MAXMAX);

	FILE *fp = fopen("printout.dat", "w");
	fprintf(fp,"%d %10.3e %10.5e %10.5e %10.5e %10.5e\n",
			iter,lx,freeE,freeDiff,inCompMax,and_err);
	fclose(fp);
}

void REPORT_SETUP(void){
	snprintf(pattern, sizeof(char)*20, "_a%03.0fc%03.0f_x%03.0fs%03.0f", Alpha_i[0][0]*100, c0_plus_n*1000, Chi_i[0][0]*100, sigma_i[0]*100);
	snprintf(phname, sizeof(char)*30, "ph%s.dat", pattern);
	snprintf(phname_elec, sizeof(char)*30, "el%s.dat", pattern);
	snprintf(Wname, sizeof(char)*30, "W%s.dat", pattern);
	snprintf(itname, sizeof(char)*30, "it%s.dat", pattern);
}

void write_it(void) {
	FILE *fp=fopen(itname,"a");
	fprintf(fp,"%d %10.5e %10.5e %10.5e %10.5e\n",
		iter,freeE,freeDiff,inCompMax,and_err);
	fclose(fp);
}

void write_ph(void) {
	int i, j, k, X, Y;
	long int ijk;
	FILE *fp=fopen(phname,"w");

	for (i=0; i<Nx; i++) for (j=0; j<Ny; j++) for (k=0; k<Nz; k++) {	
		ijk = (i*Ny + j)*Nz + k;
		fprintf(fp,"%10.3e %10.3e %10.3e ",rx[i], ry[j], rz[k]);
		fprintf(fp, "%10.4e ", phA[ijk]);
		for (X=0; X<NF_N; X++) {
			fprintf(fp, "%10.4e ", PHA_T[X][ijk]);	
			fprintf(fp, "[ ");
			for (Y=0; Y<K_i[X]-1; Y++) fprintf(fp, "%10.4e ", PHA[X][_Yijk]);
			fprintf(fp, "%10.4e ] ", PHA[X][_Yijk]);
		}
		fprintf(fp, "%10.4e\n", phB[ijk]);
	}

	fclose(fp);
}

void write_elec(void) {
	int i, j, k;
	long int ijk;

	FILE *fp=fopen(phname_elec,"w");

	for(i=0;i<Nx;i++) for(j=0;j<Ny;j++) for(k=0;k<Nz;k++) {
		ijk = (i*Ny + j)*Nz + k;
		fprintf(fp,"%10.3e %10.3e %10.3e %10.4e %10.4e %10.4e %10.4e\n",rx[i],ry[j],rz[k],-pot_elec[ijk],rho_elec_plus[ijk],rho_elec_minus[ijk],rho_elec_polym[ijk]);
	}

	fclose(fp);
}

void write_W(void) {
	int X, Y;
	long int ijk;

	FILE *fp=fopen(Wname,"w");

	for (ijk=0; ijk<NxNyNz; ijk++){
		for (X=0; X<NF_N; X++) {
			fprintf(fp, "[ ");
			for (Y=0; Y<K_i[X]-1; Y++) fprintf(fp, "%10.4e ", WA[X][_Yijk]);
			fprintf(fp, "%10.4e ] ", WA[X][_Yijk]);
		}
		fprintf(fp, "%10.5e %10.5e %10.5e\n", wB[ijk], eta[ijk], pot_elec[ijk]);
	}

	fclose(fp);
}

void write_qA(char *name, int N, int Ns, double **qA) {
	int i,s;
	FILE *fp=fopen(name,"w");
	for (s=0; s<=Ns; s++) for(i=0;i<N;i++) {
		fprintf(fp, "%10.5e %10.5e\n", i*dz, qA[i][s]);
	}
	fclose(fp);
	printf("wrote qA\n");
}

void write_qA_s(int s, char *name, int N, int Ns, double **qA) {
	int i;
	FILE *fp=fopen(name,"w");
	for(i=0;i<N;i++) {
		fprintf(fp, "%10.5e %10.5e\n", i*dz, qA[i][s]);
	}
	fclose(fp);
	printf("wrote qA\n");
}

void write_finish(void) {
	FILE *fp = fopen("(0) Program_Over.dat", "w");
	fclose(fp);
}
	
