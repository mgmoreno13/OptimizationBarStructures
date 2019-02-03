#ifndef _posprogid_h_
#define _posprogid_h_

/******************************************************************************
 Posproceso.
 ******************************************************************************/
void PosproGid(int npnod, int nelem, int ndime, int nnode, int ngdln, int ntotv, int nmats,int ncaso,
			   int paso, int** lnods, int* matnu, int* matva, double** coord, double **props,
			   double* despl, double* fuerc, double* aslod, double** carpp, char* titulo,double** mvector);
#endif
