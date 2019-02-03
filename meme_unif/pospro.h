#ifndef _pospro_h_
#define _pospro_h_

//*****************************
// Rutina para posproceso
//*****************************


/******************************************************************************
 Posproceso.
 ******************************************************************************/
void Pospro(int npnod, int nelem, int npres, int ntipo, int nnode, int ngdln, int nmats, 
			int ngaus, int ndime, int ntens, int nevab, int ntotv, int** lnods, int* matnu,
			int* iffix, double** coord, double* aslod, double* despl, double** carga);

#endif
