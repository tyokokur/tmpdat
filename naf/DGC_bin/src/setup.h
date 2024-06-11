#include "common.h"
void SETUP(void);
void ASSIGN_CONSTS(void);
void CALLOC_FIELDS(void);
void LOAD_AA(char *seq_file);
void READ_AA(int *i, FILE *fp);
void INIT_FIELDS(int init_opt);
void INIT_PHA_WA(char *Win, int win_box);
void INIT_PHA_WA_SPREAD(char *Win, int win_box);
void INIT_PHA_WAz(char *Win, int win_box);
void INIT_PHA_CONST(double xCmax);
void INIT_PHA_xC(double xCmax, double xCltot, double xCneck);
void INIT_PHA_PHAz(char *Win);
void INIT_PHA_PHA(char *Win);
void ANALY_PB(void);
void INIT_ELECS(void);
void ELEC_POLYM(void);
void get_ELFIELD(void);

void xy_NBC(void); 
void xy_NBC_HALF(void); 

void CLEAN(void);
