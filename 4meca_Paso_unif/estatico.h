#ifndef _estatico_h_
#define _estatico_h_

//*****************************
// Rutina Estatico
//*****************************


/******************************************************************************
 Solucion a problema Estatico Lineal
 ******************************************************************************/
int Estatico(int nelem, int npnod, int nevab, int ndime, int ntipo, int ncaso, int npres, int iwrit,
			 int indso, int isale, int ngdln, int nnode, int ntotv, int nincl, int nreso, int neqns,
			 int nwktl, int ishot, int** lnods, int* matnu, int** inpre,  int* iffix,int* ntips,
			 int* maxad, int* nodea, int* nodpr, int** leqns, int* lincl, int* lreso,double* react,
			 double** coord, double** presc, double* fixed, double* despl, double* aslod,double* stiff,
			 double*** srmat, double** props, double* xlong, double** vectr, double** rigid,
			 double** carpp, double** fuepp, double** carga, double** fuemp, double** girom,
			 double** girtm, double** tempr,double* fuerc, double* valor, double ***arear,
			 double ***arerf, double ***arerc, double ***arecr, double ****dtcon, int** indfu,
			 double** fuerb, double** fuefl, double* wcarg, double* xincl, double** resor,
			 double** astif, double* vectp, double* vectu, double* vect1, double* deslo, int *ubiCat,
			 double *efi_max, bool generaMatRigidez, bool Optimiza);

#endif
