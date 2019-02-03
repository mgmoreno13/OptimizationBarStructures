#ifndef _rigidez_h_
#define _rigidez_h_

//*****************************
// Rutinas para Calculo d ematrices de Rigidez y Masa
//*****************************


/******************************************************************************
 Calcula matriz de Rigidez de cada elemento.
 ******************************************************************************/
void Rigimat(int nelem, int nevab, int ndime, int ntipo, int isale, int* matnu, int* ntips, double*** srmat, 
			 double** props, double* xlong, double** vectr, double** rigid, double** girom, 
			 double** girtm, double** tempr);


void ApoyoInclinado(int ielem, int ndime, int ngdln, int nevab, int isale, int nincl, int** lnods, int* lincl,double** girom,  
					double** girtm, double** tempr, double** rigid, double* xincl);


/******************************************************************************
 Calcula matriz de Masa Consistente de cada elemento.
 ******************************************************************************/
void Masaconsis(int nelem, int nevab, int ndime, int ntipo, int isale, int* matnu, int* ntips,  
				double*** smmat,double** props, double* xlong, double** vectr, double** xmasa,  
				double** girom,double** girtm, double** tempr);


/******************************************************************************
 Calcula matriz de Masa Concentrada de cada elemento.
 ******************************************************************************/
void Masaconsen(int nelem, int nevab, int ndime, int ntipo, int isale, int* matnu, int* ntips,  
				double*** smmat,double** props, double* xlong, double** vectr, double** xmasa,  
				double** girom,double** girtm, double** tempr);

#endif
