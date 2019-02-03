#ifndef _fuerzas_h_
#define _fuerzas_h_


/******************************************************************************
 Calculo de fuerzas nodales equivalentes.
 ******************************************************************************/
void Fuerzas(int nelem, int npnod, int ndime, int ntipo, int nevab, int nnode, int ngdln,
	         int iwrit, int icarg, int** lnods, int* matnu, int** indfu, int* ntips,
	         double** coord, double* xlong, double** props, double** fuerb, double* wcarg,
	         double** carpp, double** fuepp, double** vectr, double** carga, double** fuemp,
	         bool Optimiza, int &npesp );


/******************************************************************************
 Ensambla fuerzas puntuales
 ******************************************************* **********************/
void Fuerzas_Puntuales(int nnodf, int nelem, int nnode, int ngdln, int iwrit, int** lnods, double** carga);


/******************************************************************************
 Calcula fuerzas equivalentes para movimientos prescritos.
 ******************************************************************************/
void Fuefix(int nelem, int ndime, int ngdln, int nevab, int isale, int nincl, int** lnods,int* lincl,
			double** girom, double** girtm, double** tempr, double** rigid, double* xincl,int* nodea,
			int* iffix, double* fixed, double*** srmat, double** carga);

// Copia el vector de fuerza (de un caso de carga dado) que se empleara para resolver el sistema.
// Se hace una copia porque despues lo empleara el optimizador.
void copia_vector_fuerza( int caso_i, int nelem, int nevab, double **carga, double ***VF );


// Escoge el vector de fuerzas del caso de carga 
void escoge_vector_fuerza( int caso_i, int nelem, int nevab, double **carga, double ***VF );

// Determina el vector de fueras de empotramiento perfecto
void vector_fuer_emp_perf( int nelem, int nevab, double **fuemp, double **carga );

// Incluye las fuerzas de peso propio al vector de fuerzas. Se usa cuando con el optmizador
void incluye_peso_propio( int npesp, int nelem, int ndime, int ntipo, int **lnods, int *matnu,
                          int *ntips, double **coord, double *xlong, double **props, double **vectr,
                          double **carga );

#endif
