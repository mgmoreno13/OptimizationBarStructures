#ifndef _memoria_h_
#define _memoria_h_

int PideMemoria(int npnod, int nelem, int npres, int ncarg, int ngdln,int nmats,
                int ndime, int nreso, int nincl, int mnode, int nprop, int nevab, int indso,int threads,
                int**& lnods, int**& inpre, int*& nareas, int*& narerf, int*& narerc, int*& narecon,
                int*& lreso, int*& lincl,int*& ntips, int*& nodpr, int*& matnu, int*& matva, int*& matvm,
                int*& iffix, int**& indfu, int*& nodea, int*& maxad, int*& mhigh, int**& leqns,
                double**& coord, double**& carga, double**& carpp, double**& presc, double**& props,
                double**& fuemp, double**& fuepp, double**& rigid, double**& xmasa, double**& vmatr,
                double**& girom, double**& girtm, double**& tempr, double**& fuefl,double**& fuerb,
                double***& srmat,double***& smmat, double**& resor, double*& xincl,
                double*& fuerc, double*& aslod, double*& despl, double*& fixed, double*& react,
                double*& vecta, double*& vect1, double*& deslo, double*& xlong, double*& valor, double*& valom,
                double*& angulo, double**& angulg, double*& wcarg, double**& vectr, double**& astif, double*& vectp,
                double*& vectu, char**& listaCatalogos, int*& ubiCat, double ***&vectores_fuerzas, int *&npesp,
                double *&efi_max, int &ncaso, int &nvpro, double* &tenlo, double* &mvalor, double* &aux, double* &aux1,
                double** &mvector,double** &vecto1, double** &deslm, double** &tenlm ,double** &vecth);

void LiberaMemoria( int npnod, int nelem, int ncarg, int ngdln, int nmats, int ndime, int nevab, int indso,
				    int**& lnods, int**& inpre, int*& nareas, int*& narerf, int*& narerc, int*& narecon,
				    int*& lreso, int*& lincl,int*& ntips, int*& nodpr, int*& matnu, int*& matva, int*& matvm,
				    int*& iffix, int**& indfu, int*& nodea, int*& maxad, int*& mhigh, int**& leqns,
				    double**& coord, double**& carga, double**& carpp, double**& presc, double**& props,
				    double**& fuemp, double**& fuepp, double**& rigid, double**& xmasa, double**& vmatr,
				    double**& girom, double**& girtm, double**& tempr, double**& fuefl,double**& fuerb,
				    double***& srmat,double***& smmat, double**& resor, double*& xincl,
				    double*& fuerc, double*& aslod, double*& despl, double*& fixed, double*& react,
				    double*& vecta, double*& vect1, double*& deslo, double*& xlong, double*& valor, double*& valom,
				    double*& angulo, double**& angulg, double*& wcarg, double**& vectr, double**& astif, double*& vectp,
				    double*& vectu, char**& listaCatalogos, int*& ubiCat, double ***&vectores_fuerzas,
				    char**& nombreA, int*& tamCatA, double***& arear, char**& nombreRF, int*& tamCatRF, double***& arerf,
					char**& nombreRC, int*& tamCatRC, double***& arerc, char**& nombreC, int*& tamCatC, double***& arecr,
					double****&dtcon, int nCatArm, int nCatAcRF, int nCatAcRC, int nCatConcr, double*& stiff, int *&npesp,
					double *&efi_max );

void mem_catalogos( char**& nombreA, int*& tamCatA, double***& arear,
					char**& nombreRF, int*& tamCatRF, double***& arerf,
					char**& nombreRC, int*& tamCatRC, double***& arerc,
					char**& nombreC, int*& tamCatC, double***& arecr,
					double****&dtcon, int nCatArm, int nCatAcRF,
					int nCatAcRC, int nCatConcr );

void lib_mem_catalogos( char**& nombreA, int*& tamCatA, double***& arear,
					    char**& nombreRF, int*& tamCatRF, double***& arerf,
					    char**& nombreRC, int*& tamCatRC, double***& arerc,
					    char**& nombreC, int*& tamCatC, double***& arecr,
					    double****&dtcon, int nCatArm, int nCatAcRF,
					    int nCatAcRC, int nCatConcr );

double*** matrixx(int i0, int i1, int j0, int j1,  int k0, int k1, char* nombre);

double*** mem_triple_ptr( int i0, int i1, char *nombre );

double**** mem_quadruple_ptr( int i0, int i1, char *nombre );

double** matrix(int i0, int i1, int j0, int j1,char* nombre);

int** imatrix(int i0, int i1, int j0, int j1,char* nombre);

double* vector(int i0,int i1,char* nombre);

char** cmatrix(int i0, int i1, int j0, int j1,char* nombre);

int* ivector(int i0,int i1,char* nombre);

char* cvector( int i0, int i1, char* nombre );

void Error(char* nombre);

void nrerror(char* msg);

void free_vector(double* v, int ini);

void free_ivector(int* v, int ini);

void free_cvector(char* v, int ini);

void free_matrixx( double ***A, int ini_i, int l, int ini_j, int m, int ini_k );

void free_matrix(double** m,int nrl,int nrh,int ncl);

void free_imatrix(int** m,int nrl,int nrh,int ncl);

void free_cmatrix( char** m, int nrl, int nrh, int ncl );

#endif
