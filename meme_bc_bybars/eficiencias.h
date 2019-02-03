#ifndef _eficiencias_h_
#define _eficiencias_h_

#include "efi_RF_o.h"
#include "efi_RF_c.h"
#include "efi_RC.h"
#include "efi_CR.h"


//double SF, SM;

void Eficiencias(int nelem, int ndime, int ntipo, int* matnu, int* ntips,
				 double** props, double* fuerc, double* valor, double* xlong,
				 double ***arear, double ***arerf, double ***arerc,
				 double ***arecr, double ****dtcon, int** indfu,
				 double** fuerb, double** fuefl, double* wcarg, int iwrit,
				 int *ubiCat, double *efi_max, int nevab );
void Cambia_Extremos_Barra(int ndime,int ntipo,int nelem, double **fuefl);
#endif
