#ifndef _datos_h_
#define _datos_h_

typedef struct{
	int nCat;
	char **nombre;
	int *tamCat;
	double ***arear;
}cat_armadura;

typedef struct{
	int nCat;
	char **nombre;
	int *tamCat;
	double ***arerf;
}cat_acero_RF;

typedef struct{
	int nCat;
	char **nombre;
	int *tamCat;
	double ***arerc;
}cat_acero_RC;

typedef struct{
	int nCat;
	char **nombre;
	int *tamCat;
	double ***arecr;
	double ****dtcon;
}cat_concreto;

int Datos(int& npnod, int& nelem, int& npres, int& ncarg, int& ntipo, int& ngdln, int& nmats,
          int& ndime, int& iwrit, int& isale, int& nreso, int& nincl, int& nfami, int& ncaso,
          int& nnode, int& mnode, int& nprop, int& ngaus, int& ntens, int& nevab, int& indso,
          int& ntotv,int &threads,int &nvpro,int**& lnods, int**& inpre, int*& nareas,
          int*& narerf, int*& narerc, int*& narecon,
          int*& lreso, int*& lincl,int*& ntips, int*& nodpr, int*& matnu, int*& matva, int*& matvm,
          int*& iffix, int**& indfu, int*& nodea, int*& maxad, int*& mhigh, int**& leqns,
          double**& coord, double**& carga, double**& carpp, double**& presc, double**& props,
          double**& fuemp, double**& fuepp, double**& rigid, double**& xmasa, double**& vmatr,
          double**& girom, double**& girtm, double**& tempr, double**& fuefl,double**& fuerb,
          double***& srmat,double***& smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
          cat_acero_RC &cAceroRC, cat_concreto &cConcr, double**& resor, double*& xincl,
          double*& fuerc, double*& aslod, double*& despl, double*& fixed, double*& react,
          double*& vecta, double*& vect1, double*& deslo, double*& xlong, double*& valor,
          double*& valom,double*& angulo, double**& angulg, double*& wcarg, double**& vectr,
          double**& astif, double*& vectp, double*& vectu, char **&listaCatalogos, int*& ubiCat,
          double ***&vectores_fuerzas, int *&npesp, double *&efi_max,
          double* &tenlo, double* &mvalor, double* &aux, double* &aux1, double** &mvector,
          double** &vecto1, double** &deslm, double** &tenlm, double** &vecth);

void CargaPropCatalogo( int nmats, int ndime, double** props, double*** arear,
						double*** arerf, double*** arerc, double*** arecr,
						int *ubiCat );

void CargaPropCatOpt( int* matva, int ndime, int nmats, double** props,
					  double*** arear, double*** arerf, double*** arerc,
					  double*** arecr, int *ubiCat );

#endif
