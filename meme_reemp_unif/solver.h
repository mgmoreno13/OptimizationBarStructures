#ifndef _solver_h_
#define _solver_h_

/******************************************************************************
 Ensambla y resuelve el sistema de ecuaciones matriciales.
 Solucion del sistema de ecuaciones utilizando eliminacion gaussiana directa.
 ******************************************************************************/
int Solucion(int npnod, int nelem, int ncaso, int ngdln, int nnode, int npres, int ntipo, int ndime, int nincl, 
			 int nreso, int nevab, int ntotv, int iwrit, int isale, int** lnods, int* nodpr, int** inpre, 
			 int* iffix, int* lreso, int* lincl, int* ntips, double** presc, double* despl, 
			 double* fixed, double** astif, double* aslod, double** carpp, double*** srmat, double** resor,
			 double** girom, double** girtm, double** tempr, double** rigid, double* xincl, double* react);


/******************************************************************************
 Ensambla las matrices de rigidez y los vectores de fuerzas.
 ******************************************************************************/
void Ensambla(int* nsvab, int ncaso, int nelem, int nnode, int ngdln, int npres, int ndime, int nincl, 
			  int nreso, int nevab, int ntotv, int iwrit, int isale,  
			  int** lnods, int* nodpr, int** inpre, int* iffix, int* lreso, int* lincl, double** presc, double* despl, 
			  double* fixed, double** astif, double* aslod, double** carpp, double*** srmat, double** resor, 
			  double** girom, double** girtm, double** tempr, double** rigid, double* xincl);

			  
/******************************************************************************
 Reduce el sistema de ecuaciones globales por eliminacion gaussiana directa.
 ******************************************************************************/ 
int Reduce(int neons, int* iffix, double* fixed, double** astif, double* aslod);


/******************************************************************************
 Realiza la sustitucion hacia atras.
 ******************************************************************************/
void Sustituye(int neons, int nelem, int ntipo, int ndime, int ngdln, int* iffix, int* ntips, int** lnods,   
			   double* fixed, double* react, double** astif, double* aslod, double* despl);


/******************************************************************************
 Links with profile solver.
 ******************************************************************************/
void Linkin(int nelem, int nnode, int nevab, int ntotv, int npres, int ngdln, int& neqns, int& nwktl, int& mkoun,  
			int** lnods,int* nodea, int* nodpr, int** inpre, int* iffix, int** leqns, int* maxad, int* mhigh, 
			double** presc, double  *fixed, double*& stiff);

			
/******************************************************************************
 Evaluates the column heights of sitffness matrix.
 ******************************************************************************/
void Colmht(int mevab,int ielem, int** leqns, int* mhigh);


/******************************************************************************
 Evaluates addresses of diagonal elements.
 ******************************************************************************/
void Addres(int neqns, int* nwktl, int* mkoun, int* maxad, int* mhigh);



/******************************************************************************
 Assembly of total stiffness vector.
 ******************************************************************************/
void Addban(int nelem, int ngdln,int ndime, int nnode, int nevab, int isale, int nincl, int nwktl, int nreso, 
			int** lnods,  int** leqns, int* maxad, int* nodea,int* lincl, int* lreso, int* iffix, double** girom, 
			double** girtm, double** tempr, double** rigid,double* stiff, double*** srmat, double* xincl, 
			double** resor);


/******************************************************************************
 Factorises (l)*(d)*(l) transpose of siffness matrix.
 ******************************************************************************/
int Decomp(int neqns, int ishot, int* maxad, double* stiff);


/******************************************************************************
 Calculate the incremental displacements.
 ******************************************************************************/
void Profile(int nelem, int ndime, int npnod, int ngdln, int ncaso, int neqns, int iwrit,int ntotv, 
			 int nincl, int nevab, int isale,    
			 int** lnods, int* nodea, int* maxad, int* iffix,int* lincl, double** carpp, double* despl, 
			 double* aslod, double* fixed, double* stiff, double*** srmat, double* react, 
			 double** girom, double** girtm, double** tempr, double** rigid, double* xincl);

			 
/******************************************************************************
 For reduce and bak-substitute iteration vectors.
 ******************************************************************************/
void Redbak(int neqns, int* maxad, double* stiff, double* aslod);


/******************************************************************************
 Ensambla el vector de fuerzas obtener la para solucion de un sistema de
 ecuaciones con gradiente conjugado.
 ******************************************************************************/
int Gradc(int nelem, int npnod, int npres, int ngdln, int ntotv, int nnode, int ncaso, int iwrit,int ndime, 
		  int nevab,  int isale, int nincl, int nreso, int** lnods, int** inpre, int* nodpr, int* iffix,
		  int* nodea, int* lincl, int* lreso, double* fixed, double** presc, double* despl, double* aslod, 
		  double** carpp, double* vectp, double* vectu, double** girom, double** girtm, double** tempr, 
		  double** rigid, double* xincl, double** resor, double*** srmat, double* react, double* vect1, 
		  double* deslo);

		  

void EvaluaMatrizVector(int nelem, int ngdln,int ndime, int nevab, int isale, int nincl, int nreso, 
                        int ntotv, int** lnods,  int* nodea, int* iffix, int* lincl,int* lreso, 
						double** girom, double** girtm, double** tempr, double** rigid, double* xincl, 
						double** resor, double* vectu, double* vectp, double* vect1, double* deslo, 
						double*** srmat);

#endif
