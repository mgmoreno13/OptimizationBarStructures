#ifndef _ValoresPropios_h_
#define _ValoresPropios_h_

//*****************************
// Rutina Estatico
//*****************************


/******************************************************************************
 Solucion a problema Estatico Lineal
 ******************************************************************************/
int ValoresPropios(int nelem, int npnod, int nevab, int ndime, int ntipo, int ncaso, int npres, int iwrit,
                   int indso, int isale, int ngdln, int nnode, int ntotv, int nincl, int nreso, int neqns,
                   int nwktl, int ishot, int** lnods, int* matnu, int** inpre,  int* iffix,int* ntips,
                   int* maxad, int* nodea, int* nodpr, int** leqns, int* lincl, int* lreso,double* react,
                   double** coord, double** presc, double* fixed, double* despl, double* aslod,double* stiff,
                   double*** srmat, double** props, double* xlong, double** vectr, double** rigid,
                   double** carpp, double** fuepp, double** carga, double** fuemp, double** girom,
                   double** girtm, double** tempr,int mkoun, int *mhigh, double* fuerc, double* valor,
                   double** fuerb, double** fuefl, double* wcarg, double* xincl, double** resor,
                   double** astif, double* vectp, double* vectu, double* vect1, double* deslo, double *tenlo,
                   double *aux, double* aux1, double*** smmat, double** xmasa , double* vecta, int nvpro,
                   double* mvalor, double** mvector, double** vecto1, int threads, double** deslm,
                   double** tenlm, double** vecth,bool generaMatRigidez);

int VP_Bajos(int nelem,int ngdln,int ndime,int nnode,int npnod,int ncaso,int ntipo,int iwrit,int nevab,int isale,
             int nincl,int nwktl,int nreso,int indso,int ntotv,int npres,int neqns,int mkoun,int ishot,int nvpro,
             int *nodpr, int** inpre, int *mhigh, int *ntips,int** lnods,int** leqns,int* maxad,
             int* nodea,int* lincl,int* lreso,int* iffix,double* fixed,double* stiff,double* xincl,
             double** resor,double* mvalor, double* despl,double *deslo, double *tenlo, double** astif,
             double* aslod,double *react, double* vecta, double *vectp, double* vectu, double* vect1, double* aux,
             double* aux1,double **carga, double** girom,double** girtm,double** tempr,double** rigid,
             double** presc,double**mvector, double** vecto1, double*** srmat,double*** smmat,int threads,
             double** deslm, double** tenlm, double** vecth);
int VP_Alto(int nelem,int ngdln,int ndime,int nnode,int npnod,int ncaso,int ntipo,int iwrit,int nevab,int isale,
             int nincl,int nwktl,int nreso,int indso,int ntotv,int npres,int neqns,int mkoun,int ishot,int nvpro,
             int *nodpr, int** inpre, int *mhigh, int *ntips,int** lnods,int** leqns,int* maxad,
             int* nodea,int* lincl,int* lreso,int* iffix,double* fixed,double* stiff,double* xincl,
             double** resor,double* mvalor, double* despl,double *deslo, double *tenlo, double** astif,
             double* aslod,double *react, double* vecta, double *vectp, double* vectu, double* vect1, double* aux,
             double* aux1,double **carga, double** girom,double** girtm,double** tempr,double** rigid,
             double** presc,double**mvector, double** vecto1, double*** srmat,double*** smmat, int threads,
             double** deslm,double** tenlm, double** vecth);
void Inicia(int ntotv,int npres, int ngdln, int indso, int **inpre,int *nodpr, int *iffix, double **presc,
			double *despl, double *vecta, double *fixed);
void Derecho(int nelem, int ntotv, int ndime, int ngdln, int nnode, int nevab, int isale, int nincl, int** lnods,
			 int* iffix, int* lincl,double** girom, double** girtm, double** tempr, double** rigid, double* xincl,
			 double* deslo, double* tenlo, double* vect1, double* aslod, double* despl, double *aux, double* aux1,
             double*** smmat, double* mvalor, double** mvector,int ivpro,int threads,double** deslm,double** tenlm,double** vecth);
void Almacena(double vectn,double* vecta, double* aslod, double* mvalor,double** mvector, double** vecto1, int ivpro, int ntotv);
void MasaVect(int nelem, int nnode, int ngdln, int nevab, int ntotv,int* iffix, int **lnods, double* deslo,
              double* tenlo, double* aux, double* aslod, double*** smmat,int threads,double** deslm,double** tenlm,double** vecth);
#endif
