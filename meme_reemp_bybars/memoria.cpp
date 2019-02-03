//*****************************
// Rutinas para manejo de memoria
//*****************************

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "memoria.h"
#include "raros.h"


//************************************************
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
                double** &mvector,double** &vecto1, double** &deslm, double** &tenlm ,double** &vecth)

{
	// Regresa 0 si no hubo ningun error.
	// Regresa -1 si no se encontro memoria.

	int mreso,mincl;

	mreso=nreso ;
	mincl=nincl ;
	if(mreso <= 0) mreso=1;
	if(mincl <= 0) mincl=1;
	ntips = ivector(1, nelem,(char*)"ntips");
	nodpr = ivector(1, npnod,(char*)"nodpr");
	matnu = ivector(1, nelem,(char*)"matnu");
	matva = ivector(1, nmats,(char*)"matva");
	matvm = ivector(1, nmats,(char*)"matvm");
	iffix = ivector(1, npnod*ngdln,(char*)"iffix");
	lreso = ivector(1, mreso,(char*)"lreso");
	lincl = ivector(1, mincl,(char*)"nincl");
	npesp = ivector(1, ncarg,(char*)"npesp");
	fuerc =  vector(1, nelem,(char*)"fuerc");
	aslod =  vector(1, npnod*ngdln,(char*)"aslod");
	despl =  vector(1, npnod*ngdln,(char*)"despl");
	fixed =  vector(1, npnod*ngdln,(char*)"fixed");
	react =  vector(1, npnod*ngdln,(char*)"react");
	xincl =  vector(1, mincl,(char*)"xincl");
	vecta =  vector(1, npnod*ngdln,(char*)"vectp");
	vect1 =  vector(1, nevab,(char*)"vect1");
	deslo =  vector(1, nevab,(char*)"deslo");
    tenlo =  vector(1, nevab,(char*)"tenlo");
	xlong =  vector(1, nelem,(char*)"xlong");
	valor =  vector(1, 100   ,(char*)"valor");
	valom =  vector(1, 100   ,(char*)"valom");
	efi_max = vector(1, nelem,(char*)"efi_max");
	inpre = imatrix(1, npnod,1,ngdln,(char*)"inpre");
	lnods = imatrix(1, nelem,1,mnode,(char*)"lnods");
	indfu = imatrix(1,     4,1,nelem,(char*)"indfu");
	coord =  matrix(1, npnod,1,ndime,(char*)"coord");
	carga =  matrix(1, nelem,1,nevab,(char*)"carga");
	carpp =  matrix(1, nelem,1,nevab,(char*)"carpp");
	props =  matrix(1, nelem,1,nprop,(char*)"props");
	presc =  matrix(1, npnod,1,ngdln,(char*)"presc");
	fuemp =  matrix(1, nelem,1,nevab,(char*)"fuemp");
	fuepp =  matrix(1, nelem,1,nevab,(char*)"fuepp");
	rigid =  matrix(1, nevab,1,nevab,(char*)"rigid");
	xmasa =  matrix(1, nevab,1,nevab,(char*)"xmasa");
	vmatr =  matrix(1, ndime,1,ndime,(char*)"vmatr");
	resor =  matrix(1, ngdln,1,mreso,(char*)"resor");
	girom =  matrix(1, nevab,1,nevab,(char*)"girom");
	girtm =  matrix(1, nevab,1,nevab,(char*)"girtm");
	tempr =  matrix(1, nevab,1,nevab,(char*)"tempr");
	fuefl =  matrix(1, nelem,1,   12,(char*)"fuefl");
	fuerb =  matrix(1,     6,1,nelem,(char*)"fuerb");
	srmat =  matrixx(1,nelem,1,nevab,1,nevab,(char*)"srmat");
	smmat =  matrixx(1,nelem,1,nevab,1,nevab,(char*)"smmat");
	vectores_fuerzas = matrixx(1,ncarg,1,nelem,1,nevab,(char*)"vecfr");
	if(ndime ==3) {
		angulo = vector(1, nelem,(char*)"angulo");
		angulg = matrix(1,     3,1,nelem,(char*)"angulg");
		wcarg =  vector(1,2*nelem,(char*)"wcarg");
		vectr =  matrix(1, nelem,1,    9,(char*)"vectr");
	}
	else{
		wcarg =  vector(1, nelem,(char*)"wcarg");
		vectr =  matrix(1, nelem,1,    4,(char*)"vectr");
	}
//*****catalogos***
	nareas  = ivector(1, 233,(char*)"nareas");
	narerc  = ivector(1,  77,(char*)"narerc");
	narerf  = ivector(1,  22,(char*)"narerf");
	narecon = ivector(1,  59,(char*)"narecon");
	listaCatalogos = cmatrix( 1, nmats, 1, 300, (char*)"Lista catalogos" );
	ubiCat = ivector( 1, nmats, (char*)"Ubicacion catalogos" );

	if(ntips<=0 ||	nodpr<=0 || matnu<=0 || iffix<=0 || lreso<=0 || lincl<=0 ||
	   ntips<=0 || nodpr<=0 || matnu<=0 || iffix<=0 || lreso<=0 || lincl<=0 ||
	   fuerc<=0 || aslod<=0 || despl<=0 || fixed<=0 || react<=0 || xincl<=0 ||
	   props<=0 || presc<=0 || fuemp<=0 || rigid<=0 || resor<=0 || girom<=0 ||
	   girtm<=0 || tempr<=0 || inpre<=0 || lnods<=0 || coord<=0 || carga<=0)
		return -1;
    

    
    if(ncaso ==2){
        mvalor  =  vector(1, nvpro,(char*)"mvalo");
        aux     =  vector(1,npnod*ngdln,(char*)"aux ");
        aux1    =  vector(1,npnod*ngdln,(char*)"aux1");
        mvector =  matrix(1, nvpro, 1, npnod*ngdln,(char*)"mvect");
        vecto1  =  matrix(1, nvpro, 1, npnod*ngdln,(char*)"vecto");
        deslm   =  matrix(1, nelem, 1, mnode*ngdln,(char*)"deslm");
        tenlm   =  matrix(1, nelem, 1, mnode*ngdln,(char*)"tenlm");
        vecth   =  matrix(1,threads,1, npnod*ngdln,(char*)"vecth");
        
    }

	switch(indso) {
		case 1:
			astif = matrix(1, npnod*ngdln, 1, npnod*ngdln,(char*)"astif");
			if(astif<=0) return -1;
			break;
		case 2:
			nodea = ivector(1, nelem+1,(char*)"nodea");
			maxad = ivector(1, npnod*ngdln+1,(char*)"maxad");
			mhigh = ivector(1, npnod*ngdln+1,(char*)"mhigh");
			leqns = imatrix(1, nevab, 1, nelem,(char*)"leqns");

			if(nodea<=0 || maxad<=0 || mhigh<=0 || leqns<=0) return -1;
			break;
		case 3:
			nodea = ivector(1, nelem,(char*)"nodea");
			vectp =  vector(1, npnod*ngdln,(char*)"vectp");
			vectu =  vector(1, npnod*ngdln,(char*)"vectu");

			if(nodea<=0 || vectp<=0 || vectu<=0 || vect1<=0 || deslo<=0)
				return -1;
			break;
	}

	return 0;
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
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
					double *&efi_max )
{
    free_ivector(ntips, 1);
    free_ivector(nodpr, 1);
    free_ivector(matnu, 1);
    free_ivector(matva, 1);
    free_ivector(matvm, 1);
    free_ivector(iffix, 1);
    free_ivector(lreso, 1);
    free_ivector(lincl, 1);
    free_ivector(npesp, 1);
    free_vector(fuerc, 1);
    free_vector(aslod, 1);
    free_vector(despl, 1);
    free_vector(fixed, 1);
    free_vector(react, 1);
    free_vector(xincl, 1);
    free_vector(vecta, 1);
    free_vector(vect1, 1);
    free_vector(deslo, 1);
    free_vector(xlong, 1);
    free_vector(valor, 1);
    free_vector(valom, 1);
    free_vector(efi_max, 1);
    free_imatrix(inpre, 1, npnod, 1);
    free_imatrix(lnods, 1, nelem, 1);
    free_imatrix(indfu, 1, 4, 1);
    free_matrix(coord, 1, npnod, 1);
    free_matrix(carga, 1, nelem, 1);
    free_matrix(carpp, 1, nelem, 1);
    free_matrix(props, 1, nelem, 1);
    free_matrix(presc, 1, npnod, 1);
    free_matrix(fuemp, 1, nelem, 1);
    free_matrix(fuepp, 1, nelem, 1);
    free_matrix(rigid, 1, nevab, 1);
    free_matrix(xmasa, 1, nevab, 1);
    free_matrix(vmatr, 1, ndime, 1);
    free_matrix(resor, 1, ngdln, 1);
    free_matrix(girom, 1, nevab, 1);
    free_matrix(girtm, 1, nevab, 1);
    free_matrix(tempr, 1, nevab, 1);
    free_matrix(fuefl, 1, nelem, 1);
    free_matrix(fuerb, 1, 6, 1);
    free_matrixx(srmat, 1, nelem, 1, nevab, 1 );
    free_matrixx(smmat, 1, nelem, 1, nevab, 1 );
    free_matrixx(vectores_fuerzas, 1, ncarg, 1, nelem, 1 );

    free_vector(stiff, 1);

    if( ndime == 3 ){
        free_vector(angulo, 1);
        free_matrix(angulg, 1, 3, 1 );
        free_vector(wcarg, 1);
        free_matrix(vectr, 1, nelem, 1);
    }
    else{
        free_vector(wcarg, 1);
        free_matrix(vectr, 1, nelem, 1);
    }

    free_ivector(nareas, 1);
    free_ivector(narerc, 1);
    free_ivector(narerf, 1);
    free_ivector(narecon, 1);
    free_cmatrix(listaCatalogos, 1, nmats, 1);
    free_ivector(ubiCat, 1);

    switch( indso ){
        case 1:
            free_matrix(astif, 1, npnod*ngdln, 1);
            break;
        case 2:
            free_ivector(nodea, 1);
            free_ivector(maxad, 1);
            free_ivector(mhigh, 1);
            free_imatrix(leqns, 1, nevab, 1);
            break;
        case 3:
            free_ivector(nodea, 1);
            free_vector(vectp, 1);
            free_vector(vectu, 1);
            break;
    }

    lib_mem_catalogos( nombreA, tamCatA, arear, nombreRF, tamCatRF, arerf,
                       nombreRC, tamCatRC, arerc, nombreC, tamCatC, arecr,
                       dtcon, nCatArm, nCatAcRF, nCatAcRC, nCatConcr );

}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void mem_catalogos( char**& nombreA, int*& tamCatA, double***& arear,
					char**& nombreRF, int*& tamCatRF, double***& arerf,
					char**& nombreRC, int*& tamCatRC, double***& arerc,
					char**& nombreC, int*& tamCatC, double***& arecr,
					double****&dtcon, int nCatArm, int nCatAcRF,
					int nCatAcRC, int nCatConcr ){
	if( nCatArm > 0 ){
		nombreA = cmatrix( 1, nCatArm, 1, 300, (char*)"Cat. Armaduras" );
		tamCatA = ivector( 1, nCatArm, (char*)"Cat. Armaduras" );
		arear   = mem_triple_ptr( 1, nCatArm, (char*)"Cat. Armaduras" );
	}
	if( nCatAcRF > 0 ){
		nombreRF = cmatrix( 1, nCatAcRF, 1, 300, (char*)"Cat. Acero RF" );
		tamCatRF = ivector( 1, nCatAcRF, (char*)"Cat. Acero RF" );
		arerf    = mem_triple_ptr( 1, nCatAcRF, (char*)"Cat. Acero RF" );
	}
	if( nCatAcRC > 0 ){
		nombreRC = cmatrix( 1, nCatAcRC, 1, 300, (char*)"Cat. Acero RC" );
		tamCatRC = ivector( 1, nCatAcRC, (char*)"Cat. Acero RC" );
		arerc    = mem_triple_ptr( 1, nCatAcRC, (char*)"Cat. Acero RC" );
	}
	if( nCatConcr > 0 ){
		nombreC = cmatrix( 1, nCatConcr, 1, 300, (char*)"Cat. Concreto" );
		tamCatC = ivector( 1, nCatConcr, (char*)"Cat. Concreto" );
		arecr   = mem_triple_ptr( 1, nCatConcr, (char*)"Cat. Concreto" );
		dtcon   = mem_quadruple_ptr( 1, nCatConcr, (char*)"Cat. Concreto" );
	}
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑

void lib_mem_catalogos( char**& nombreA, int*& tamCatA, double***& arear,
					    char**& nombreRF, int*& tamCatRF, double***& arerf,
					    char**& nombreRC, int*& tamCatRC, double***& arerc,
					    char**& nombreC, int*& tamCatC, double***& arecr,
					    double****&dtcon, int nCatArm, int nCatAcRF,
					    int nCatAcRC, int nCatConcr )
{
	if( nCatArm > 0 ){
		free_cmatrix( nombreA, 1, nCatArm, 1 );
		for( int i = 1; i <= nCatArm; i++ ) free_matrix( arear[i], 1, tamCatA[i], 1 );
		free( arear + 1 );
		free_ivector( tamCatA, 1 );
	}
	if( nCatAcRF > 0 ){
		free_cmatrix( nombreRF, 1, nCatAcRF, 1 );
		for( int i = 1; i <= nCatAcRF; i++ ) free_matrix( arerf[i], 1, tamCatRF[i], 1 );
		free( arerf + 1 );
		free_ivector( tamCatRF, 1 );
	}
	if( nCatAcRC > 0 ){
		free_cmatrix( nombreRC, 1, nCatAcRC, 1 );
		for( int i = 1; i <= nCatAcRC; i++ ) free_matrix( arerc[i], 1, tamCatRC[i], 1 );
		free( arerc + 1 );
		free_ivector( tamCatRC, 1 );
	}
	if( nCatConcr > 0 ){
		free_cmatrix( nombreC, 1, nCatConcr, 1 );
		for( int i = 1; i <= nCatConcr; i++ ){
			free_matrix( arecr[i], 1, tamCatC[i], 1 );
			for( int j = 1; j <= tamCatC[i]; j++ ) free_matrix( dtcon[i][j], 1, 3, 1 );
			free( dtcon[i] + 1 );
		}
		free( arecr + 1 );
		free( dtcon + 1 );
		free_ivector( tamCatC, 1 );
	}
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
// Reserva mememoria para un cuadruple apuntador
double**** mem_quadruple_ptr( int i0, int i1, char *nombre ){
	double ****a;
	if( ( a = (double****) calloc( i1-i0+1, sizeof(double***) ) ) == NULL )
		Error( nombre );
	a -= i0;
	return a;
}

// Reserva mememoria para un triple apuntador
double*** mem_triple_ptr( int i0, int i1, char *nombre ){
	double ***a;
	if( ( a = (double***) calloc( i1-i0+1, sizeof(double**) ) ) == NULL )
		Error( nombre );
	a -= i0;
	return a;
}

double*** matrixx(int i0, int i1, int j0, int j1, int k0, int k1, char* nombre)
// Reserva memoria para un arreglo tridimensional
{
	register int i,j;
	double*** a;

	if((a = (double***)calloc(i1-i0+1, sizeof(double**))) == NULL) Error(nombre);
	a-=i0;
	for(i=i0; i<=i1; i++) {
		if((a[i] = (double**)calloc(j1-j0+1,sizeof(double*))) == NULL) Error(nombre);
		a[i]-=j0;
	}
	for(i=i0; i<=i1; i++)
		for(j=j0; j<=j1; j++) {
			if((a[i][j] = (double*)calloc(k1-k0+1,sizeof(double))) == NULL) Error(nombre);
			a[i][j]-=k0;
		}

	return a;
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
double** matrix(int i0, int i1, int j0, int j1,char* nombre)
// Allocates memory for a matrix.
{
	register int i;
	double** a;

	if((a = (double**)calloc(i1-i0+1, sizeof(double*))) == NULL)
		Error(nombre);

	a-=i0;

	for(i=i0; i<=i1; i++) {
		if((a[i] = (double*)calloc(j1-j0+1, sizeof(double))) == NULL)
			Error(nombre) ;

		a[i]-=j0;
	}

	return a ;
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
int** imatrix(int i0, int i1, int j0, int j1,char* nombre)
// Allocates memory for a matrix.
{
	register int i;
	int** a;

	if((a = (int**)calloc(i1-i0+1, sizeof(int* ))) == NULL) Error(nombre);

	a-=i0;

	for(i=i0; i<=i1; i++) {
		if((a[i] = (int* )calloc(j1-j0+1,sizeof(int))) == NULL) Error(nombre);
		a[i]-=j0;
	}

	return a;
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
double* vector(int i0,int i1,char* nombre)
// Allocates memory for vector with indices from i0 to i1.
{
	double* v;

	if((v = (double*)calloc(i1-i0+1,sizeof (double))) == NULL) Error(nombre);

	return v-i0;
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
char** cmatrix(int i0, int i1, int j0, int j1,char* nombre)
// Allocates memory for a matrix.
{
	char** a;
	int i;

	if((a = (char**)calloc(i1-i0+1,sizeof(char*))) == NULL) Error(nombre);

	a -= i0;

	for(i=i0; i<=i1; i++) {
		if((a[i] = (char*)calloc(j1-j0+1,sizeof(char))) == NULL) Error(nombre);
		a[i] -= j0;
	}

	return a;
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
int* ivector(int i0,int i1,char* nombre)
// Allocates memory for vector with indices from i0 to i1.
{
	int* v;

	if((v = (int* )calloc(i1-i0+1,sizeof (int))) == NULL) Error(nombre);

	return v-i0;
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
// Allocates memory for vector with indices from i0 to i1.
char* cvector( int i0, int i1, char* nombre ){
	char* v;
	if( ( v = (char*) calloc( i1-i0+1, sizeof (char) ) ) == NULL ) Error(nombre);
	return v-i0;
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Error(char* nombre)
{
	char mensaje[80];

	sprintf(mensaje, "No hay memoria, parado en: %s\n\nIntente resolverlo con \
			otro metodo de solucion -Ver Ayuda-", nombre);
	fprintf(fp16, "%s", mensaje);
	//MessageYesNo("Error de memoria");
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void nrerror(char* msg)
// Prints error message for num. rec.
{
	//MessageYesNo(msg);
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void free_vector(double* v, int ini)
// Frees memory for vector v.
{
	free(v+ini);
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void free_ivector(int* v, int ini)
// Frees memory for vector v.
{
	free(v+ini);
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
// Frees memory for vector v.
void free_cvector(char* v, int ini)
{
	free(v+1);
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
// Libera la memoria de un arreglo triple de tamano l x m x n
void free_matrixx( double ***A, int ini_i, int l, int ini_j, int m, int ini_k ){
	register int i, j;
	for( i = ini_i; i <= l; i++ ){
        for( j = ini_j; j <= m; j++ ){
            free( A[i][j] + ini_k );
        }
        free( A[i] + ini_j );
	}
	free( A + ini_i );
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void free_matrix(double** m,int nrl,int nrh,int ncl)
{
	register int i;

	for(i=nrh; i>=nrl; i--) free((char*)(m[i]+ncl));

	free((char*)(m+nrl));
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void free_imatrix(int** m,int nrl,int nrh,int ncl)
{
	register int i;

	for(i=nrh; i>=nrl; i--) free((char*)(m[i]+ncl));

	free((char*) (m+nrl));
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑

void free_cmatrix( char** m, int nrl, int nrh, int ncl ){
	register int i;
	for( i = nrh; i >= nrl; i-- ) free((char*)(m[i]+ncl));
	free((char*) (m+nrl));
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
