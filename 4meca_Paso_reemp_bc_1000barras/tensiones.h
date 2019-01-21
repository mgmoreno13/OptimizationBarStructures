#ifndef _tensiones_h_
#define _tensiones_h_

//*****************************
// Rutina Tensiones
//*****************************


/******************************************************************************
 Calcula las tensiones y esfuerzos en los puntos de Gauss.
 ******************************************************************************/
void Tensiones(int nelem, int nnode, int ngdln, int ndime, int ntipo, int nevab,
			   int ntotv, int iwrit, int isale, int nincl, int* matnu, int** lnods, int* lincl, 
			   double** props, double* despl, double*** srmat, double* fuerc, double** fuepp, 
			   double** fuefl, double** vectr, double* xlong, double** girom, double** girtm, 
			   double** tempr, double** rigid, double* xincl, double* valor);


void Reacci(int nelem, int ndime, int npnod, int ngdln, int nincl, int nevab, int ntotv, int iwrit, int isale, 
			int** lnods, int* nodea, int* iffix, int* lincl,  double* despl, double* fixed, double*** srmat, 
			double** carpp, double* react, double** girom, double** girtm, double** tempr, double** rigid, 
			double* xincl);

#endif
