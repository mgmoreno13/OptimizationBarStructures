#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "datos.h"
#include "memoria.h"
#include "raros.h"

// --------------------------------------------------------------------------------
// Rutina Datos
//
// Lee datos de la tipologia de la malla, condiciones de contorno y
// propiedades del material.
// --------------------------------------------------------------------------------
int Datos(int& npnod, int& nelem, int& npres, int& ncarg, int& ntipo, int& ngdln, int& nmats,
		  int& ndime, int& iwrit, int& isale, int& nreso, int& nincl, int& nfami, int& ncaso,
		  int& nnode, int& mnode, int& nprop, int& ngaus, int& ntens,int& nevab, int& indso, int& ntotv,
		  int**& lnods, int**& inpre, int*& nareas, int*& narerf, int*& narerc, int*& narecon,
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
		  double ***&vectores_fuerzas, int *&npesp, double *&efi_max )
{
	/*	Tipo de problema:
	 1  armaduras.
	 2  estructuras reticulares.
	 3  combinacion de armaduras, estructuras reticulares y articulaciones.
	 */
	// Regresa 0 si no hugo ning√∫n error.
	// Regresa -1 si se cuenta con una versi√≥n de evaluaci√≥n.
	// Regresa -2 si no hubo memoria suficiente.

	int ielem, inode, idime, ipres, igdln, imats, iprop, ipnod, ireso, iincl;
	int ndato, nbar, nCatalogos, aux, nCatArm, nCatAcRF, nCatAcRC, nCatConcr;
	int iArm, iRF, iRC, iCon, compa;

	int CargaCatalogo =1;
	double factor;
	char barra[200];
    int i;
    char letrero[200];	//Contenido no importante, solo para lectura de letreros en archivos de datos.
	FILE *fp20;



	// Lee parametros de control:
    fscanf(fp5, "%s %d\n", letrero, &npnod);
    fscanf(fp5, "%s %d", letrero, &nelem);
    fscanf(fp5, "%s %d", letrero, &npres);
    fscanf(fp5, "%s %d", letrero, &ncarg);
    fscanf(fp5, "%s %d", letrero, &ntipo);
    fscanf(fp5, "%s %d", letrero, &ngdln);
    fscanf(fp5, "%s %d", letrero, &nmats);
    fscanf(fp5, "%s %d", letrero, &iwrit);
    fscanf(fp5, "%s %d", letrero, &isale);
    fscanf(fp5, "%s %d", letrero, &ndime);
    fscanf(fp5, "%s %d", letrero, &nreso);
    fscanf(fp5, "%s %d", letrero, &nincl);
    fscanf(fp5, "%s %d", letrero, &nfami);
    fscanf(fp5, "%s %d", letrero, &ncaso);


    //Lectura de letreros para definir barras
    for (i=1; i<=8;i++) fscanf(fp5, "%s", letrero);





	nnode = 2;
	mnode = nnode;
	if(ndime == 3 && ntipo >= 2) mnode=3;
	nprop = (ndime == 2) ? 6 : 9;
	ngaus = 2;
	ntens = ngdln;
	nevab = ngdln*nnode;
	factor= 3.141592654/180.0;
	fprintf(fp16,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d \n",
			npnod, nelem, npres, ncarg, ntipo, ngdln, nmats, iwrit, isale, ndime,
			nreso, nincl, nfami,ncaso);

	if(PideMemoria(npnod,nelem,npres,ncarg,ngdln,nmats,ndime,nreso,nincl,
				   mnode,nprop,nevab,indso, lnods,inpre,nareas,narerf,narerc,
				   narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,iffix,indfu,
				   nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
				   fuepp,rigid,xmasa,vmatr,girom,girtm,tempr,fuefl,fuerb,srmat,
				   smmat,resor, xincl,fuerc,aslod,despl,fixed,react,vecta,
				   vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
				   vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,
				   efi_max)) return -2;
	// Lee y escribe conexiones nodales y numero de propiedad de material
    // Tipo de problema Articulado y Empotrado
	if(ntipo != 3) {
		for(ielem=1; ielem<=nelem; ielem++) {
			fscanf(fp5, "%d", &matnu[ielem]);
			for(inode=1; inode<=mnode; inode++)
				fscanf(fp5,"%d",&lnods[ielem][inode]);
            if(ntipo==1) {
                lnods[ielem][3]==0;
                angulo[ielem]=angulg[1][ielem]=angulg[2][ielem]=angulg[3][ielem]=0.0;
            }
			if (ndime ==3 && ntipo == 2 && lnods[ielem][3]==0) fscanf(fp5,"%lf",&angulo[ielem]);
			if (ndime ==3 && ntipo == 2 && lnods[ielem][3]==-1)fscanf(fp5,"%lf %lf %lf",&angulg[1][ielem],&angulg[2][ielem],&angulg[3][ielem]);
		}

		if(iwrit != 0) {
			for(ielem=1; ielem<=nelem; ielem++) {
				fprintf(fp16, "%d\t",matnu[ielem]);
				for(inode=1; inode<=mnode; inode++)
					fprintf(fp16, "%d\t", lnods[ielem][inode]);
				if (ndime ==3 && ntipo == 2 && lnods[ielem][3]==0) fprintf(fp16,"%lf\t",angulo[ielem]);
				if (ndime ==3 && ntipo == 2 && lnods[ielem][3]==-1)fprintf(fp16,"%lf\t %lf\t %lf\t",angulg[1][ielem],angulg[2][ielem],angulg[3][ielem]);
				fprintf(fp16, "\n");
			}
		}
		for(ielem=1; ielem<=nelem; ielem++)
			ntips[ielem]=ntipo;
	}
    // Tipo de problema combinado (articulados y empotrados)
	else {
		for(ielem=1; ielem<=nelem; ielem++) {
			fscanf(fp5, "%d", &matnu[ielem]);
           	for(inode=1; inode<=mnode; inode++)
				fscanf(fp5, "%d", &lnods[ielem][inode]);
            //ntips guarda el tipo de barra.
            //1=AA; 2=EE; 3=EA; 4=EA, donde A-->Articualdo, E--> Empotrado
			fscanf(fp5, "%d ", &ntips[ielem]);
			if (ndime ==3 && lnods[ielem][3]==0) fscanf(fp5,"%lf",&angulo[ielem]);
			if (ndime ==3 && lnods[ielem][3]==-1)fscanf(fp5,"%lf %lf %lf",&angulg[1][ielem],&angulg[2][ielem],&angulg[3][ielem]);
		}

		if(iwrit != 0) {
			for(ielem=1; ielem<=nelem; ielem++) {
				fprintf(fp16, "%d\t", matnu[ielem]);
				for(inode=1; inode<=mnode; inode++)
					fprintf(fp16, "%d\t", lnods[ielem][inode]);
				fprintf(fp16, "%d\t", ntips[ielem]);
				if (ndime ==3 && lnods[ielem][3]== 0) fprintf(fp16,"%lf\t",angulo[ielem]);
				if (ndime ==3 && lnods[ielem][3]==-1)fprintf(fp16,"%lf\t %lf\t %lf\t",angulg[1][ielem],angulg[2][ielem],angulg[3][ielem]);
				fprintf(fp16, "\n");
			}
		}
	}

	for(ipnod=1; ipnod<=npnod; ipnod++)
		for(idime=1; idime<=ndime; idime++) coord[ipnod][idime] = 0.0;


    //Lectura de letreros de coordenadas de nodos
    for ( i=1; i<=4; i++) fscanf(fp5, "%s", letrero);

	// Lee y escribe coordenadas nodales:
	for(ipnod=1; ipnod<=npnod; ipnod ++){
		for(idime=1; idime<=ndime; idime++)
			fscanf(fp5, "%lf", &coord[ipnod][idime]);
	}

	if(iwrit != 0) {
		for(ipnod=1; ipnod<=npnod; ipnod++) {
			for(idime=1; idime<=ndime; idime++)
				fprintf(fp16, "%lf\t", coord[ipnod][idime]);
			fprintf(fp16, "\n");
		}
	}
    // Lectura de letreros sobre restricciones:
    for (i=1; i<=14;i++) fscanf(fp5, "%s", letrero);

	//Lee y escribe movimientos prescritos:
	for(ipres=1; ipres<=npres; ipres++) {
		fscanf(fp5, "%d", &nodpr[ipres]);
		for(igdln=1; igdln<=ngdln; igdln++)
			fscanf(fp5, "%d", &inpre[ipres][igdln]);
		for(igdln=1; igdln<=ngdln; igdln++)
			fscanf(fp5,"%lf",&presc[ipres][igdln]);
	}

	if(iwrit != 0) {
		for(ipres=1; ipres<=npres; ipres++) {
			fprintf(fp16, "%d\t", nodpr[ipres]);
			for(igdln=1; igdln<=ngdln; igdln++)
				fprintf(fp16, "%d\t", inpre[ipres][igdln]);
			for(igdln=1; igdln<=ngdln; igdln++)
				fprintf(fp16, "%lf\t", presc[ipres][igdln]);
			fprintf(fp16, "\n");
		}
	}
    
    
    // Lectura de letreros sobre propiedades de los materiales:
    for (i=1; i<=11;i++) {
        fscanf(fp5, "%s", letrero);
     //   printf("%s \t", letrero);
    }
	// Lee y escribe propiedades de los materiales:
	for(imats=1; imats<=nmats; imats++) {
		// Lee las propiedades antes del tipo de material
		for(iprop=1; iprop<=nprop-1; iprop++){
			fscanf(fp5, "%lf", &props[imats][iprop]);
		}
        //printf("%s  props=%lf \n",listaCatalogos[imats]+2,props[imats][8]);
		// Coloca en una lista los nombres de los catalogos
   //     if(imats ==1)listaCatalogos[imats]+2=(char*)"seccione.dat";
   //     if(imats ==2)listaCatalogos[imats]+2 =(char*)"seccione.dat";
   //     if(imats ==3)listaCatalogos[imats]+2 =(char*)"CatalogoCONC.dat";
   //     listaCatalogos[imats]+2=arch;
   //     letero= (char*)props[imats][8];
   //     printf("letrero = %s",letrero);
		listaCatalogos[imats][1] = props[imats][8];
        fscanf( fp5, "%s", listaCatalogos[imats]+2 );
        printf("%s  props=%lf \n",listaCatalogos[imats]+2,props[imats][8]);
        printf("%s aqui \n",listaCatalogos[imats]+2);
    /*    compa=1000;
        compa=strcmp(listaCatalogos[imats]+2,"seccione.dat");
        if(compa == 0)props[imats][8]=0;
        printf("compa = %d \n",compa);
        compa=1000;
        compa=strcmp(listaCatalogos[imats]+2,"CatalogoRF.dat");
        if(compa == 0)props[imats][8]=1;
        printf("compa = %d \n",compa);
        compa=1000;
        compa=strcmp(listaCatalogos[imats]+2,"CatalogoRC.dat");
        if(compa == 0)props[imats][8]=2;
        printf("compa = %d \n",compa);
        compa=1000;
        compa=strcmp(listaCatalogos[imats]+2,"CatalogoCONCR.dat");
        if(compa == 0)props[imats][8]=3;
        printf("compa = %d \n",compa);
        printf("%s  props=%lf \n",listaCatalogos[imats]+2,props[imats][8]);
		// Lee la ultima propiedad, que es el numero de seccion
		fscanf( fp5, "%lf", &props[imats][nprop]);
	}

	// Carga los catalogos de la lista
		// De la lista quita los repetidos. nCatalogos es el numero de catalogos
		// distintos
		// nCatalogos = quitaRepetidos_ListaCatalogos( listaCatalogos );
		nCatalogos = nmats;
		// Cuenta cuantos catalogos se necesican para cada material
		nCatArm  = 0;
		nCatAcRF = 0;
		nCatAcRC = 0;
		nCatConcr= 0;
		for( i = 1; i <= nCatalogos; i++ ){
			switch( (int)listaCatalogos[i][1] ){
				case 0: nCatArm++;   break;
				case 1: nCatAcRF++;  break;
				case 2: nCatAcRC++;  break;
				case 3: nCatConcr++; break;
			}
		}
		mem_catalogos( cArmaduras.nombre, cArmaduras.tamCat, cArmaduras.arear,
					   cAceroRF.nombre, cAceroRF.tamCat, cAceroRF.arerf,
					   cAceroRC.nombre, cAceroRC.tamCat, cAceroRC.arerc,
					   cConcr.nombre, cConcr.tamCat, cConcr.arecr, cConcr.dtcon,
					   nCatArm, nCatAcRF, nCatAcRC, nCatConcr );
		iArm = 1;
		iRF  = 1;
		iRC  = 1;
		iCon = 1;
	cArmaduras.nCat = nCatArm;
	cAceroRF.nCat   = nCatAcRF;
	cAceroRC.nCat   = nCatAcRC;
	cConcr.nCat     = nCatConcr;

	for( i = 1; i <= nCatalogos; i++ ){
		switch( (int)listaCatalogos[i][1] ){
			// Armaduras
			case 0:
				//area-inerciax-inerciay-???-???-????
				fp20 = fopen( listaCatalogos[i]+2, "rt" );
				fscanf( fp20, "%s %d", letrero, &ndato );

				cArmaduras.tamCat[iArm] = ndato;
				strcpy( cArmaduras.nombre[iArm]+1, listaCatalogos[i]+2 );
				cArmaduras.arear[iArm] = matrix( 1, ndato, 1, 6, (char*)"Cat. Armaduras" );
				ubiCat[i] = iArm;

				for( imats = 1; imats <=ndato; imats++ ){
					fscanf(fp20, "%d %s", &aux, barra);
					for(iprop=1; iprop<=6; iprop++)
						fscanf(fp20, "%lf", &cArmaduras.arear[iArm][imats][iprop]);
				}
				fclose(fp20);
				iArm++;
                props[imats][8]=0.0;
                printf("Armaduras \n");
				break;
			// Acero rolado en frio
			case 1:
				fp20 = fopen( listaCatalogos[i]+2, "rt" );
				//Á-Ix-Iy-Sx-Sy-Qx-Qy-Tx-Ty-H-B-D-R
				fscanf( fp20, "%s %d", letrero, &ndato);

				cAceroRF.tamCat[iRF] = ndato;
				strcpy( cAceroRF.nombre[iRF]+1, listaCatalogos[i]+2 );
				cAceroRF.arerf[iRF] = matrix( 1, ndato, 1, 13, (char*)"Cat. Ac. RF" );
				ubiCat[i] = iRF;

				for(imats=1; imats<=ndato; imats++) {
					fscanf(fp20, "%d %s", &aux,barra);
					for(iprop=1; iprop<=13; iprop++)
						fscanf(fp20, "%lf", &cAceroRF.arerf[iRF][imats][iprop]);
				}
				fclose(fp20);
                props[imats][8]=1.0;
				iRF++;
                printf("Rolado en Frio \n");
				break;
			// Acero rolado en caliente
			case 2:
				fp20 = fopen( listaCatalogos[i]+2, "rt" );
				// E-Fy-A-Ix-Iy-J-Sx-Sy-Qx-Qy-tx-ty-nu
				fscanf( fp20, "%s %d",letrero, &ndato);

				cAceroRC.tamCat[iRC] = ndato;
				strcpy( cAceroRC.nombre[iRC]+1, listaCatalogos[i]+2 );
				cAceroRC.arerc[iRC] = matrix( 1, ndato, 1, 13, (char*)"Cat. Ac. RC" );
				ubiCat[i] = iRC;

				for(imats=1; imats<=ndato; imats++) {
					fscanf(fp20, "%d %s", &aux,barra);
					for(iprop=1; iprop<=13; iprop++){
						fscanf(fp20, "%lf", &cAceroRC.arerc[iRC][imats][iprop]);
                        printf ("%lf \t", cAceroRC.arerc[iRC][imats][iprop]);
					}
                    printf("\n");
				}
				fclose(fp20);
                props[imats][8]=2.0;
                printf("Rolado en Caliente \n");
				iRC++;
				break;
			// Concreto
			case 3:
				fp20 = fopen( listaCatalogos[i]+2, "rt" );
				//b-h-f'c,fy,numero d ebarras de acero y coordenadas
				fscanf( fp20, "%s %d",letrero, &ndato);

				cConcr.tamCat[iCon] = ndato;
				strcpy( cConcr.nombre[iCon]+1, listaCatalogos[i]+2 );
				cConcr.arecr[iCon] = matrix( 1, ndato, 1, 5, (char*)"Cat. Concreto" );
				cConcr.dtcon[iCon] = mem_triple_ptr( 1, ndato, (char*)"Cat. Concreto" );
				ubiCat[i] = iCon;

				for(imats=1; imats<=ndato; imats++) {
					fscanf(fp20, "%d %s", &aux, barra);
					for(iprop=1; iprop<=5; iprop++)
						fscanf(fp20, "%lf", &cConcr.arecr[iCon][imats][iprop]);
					nbar=(int)(cConcr.arecr[iCon][imats][5]+0.5);

					cConcr.dtcon[iCon][imats] = matrix( 1, 3, 1, nbar, (char*)"dtcon");

					for(iprop=1; iprop<=nbar; iprop++)
						fscanf(fp20, "%lf %lf %lf", &cConcr.dtcon[iCon][imats][1][iprop],
							   &cConcr.dtcon[iCon][imats][2][iprop],
							   &cConcr.dtcon[iCon][imats][3][iprop]);
				   }
				fclose(fp20);
                props[imats][8]=3.0;
                printf("Concreto \n");
				iCon++;
				break;
			default:
				break;
		}
	}

	// Numero de catalogo de cada material
	if(CargaCatalogo ==1) CargaPropCatalogo(nmats,ndime,props,cArmaduras.arear,
											cAceroRF.arerf,cAceroRC.arerc,
											cConcr.arecr,ubiCat);
	if(iwrit != 0) {
		for(imats=1; imats<=nmats; imats++) {
			for(iprop=1; iprop<=nprop; iprop++)
				fprintf(fp16, "%e\t", props[imats][iprop]);
			fprintf(fp16, "\n");
		}
	}

//    printf("termino de leer materiales \n");
    
    // Lectura de letreros sobre apoyos elasticos:
    for (i=1; i<=8;i++) fscanf(fp5, "%s", letrero);
    //Lee y escribe apoyos elasticos
	if(nreso > 0){
		for(ireso=1; ireso<=nreso; ireso++) {
			fscanf(fp5,"%d",&lreso[ireso]);
			for(igdln=1; igdln<=ngdln; igdln++)
				fscanf(fp5, "%lf", &resor[igdln][ireso]);
		}

		if(iwrit != 0) {
			for(ireso=1; ireso<=nreso; ireso++) {
				fprintf(fp16, "%d\t", lreso[ireso]);
				for(igdln=1; igdln<=ngdln; igdln++)
					fprintf(fp16, "%lf\t", resor[igdln][ireso]);
				fprintf(fp16, "\n");
			}
		}
	}
 //   printf("termino apoyos elasticos \n");

    // Lectura de letreros sobre apoyos inclinados:
    for (i=1; i<=3;i++) fscanf(fp5, "%s", letrero);
    //Lee y escribe apoyos inclinados
	if(nincl > 0){
		for(iincl=1; iincl<=nincl; iincl++)
			fscanf(fp5,"%d %lf",&lincl[iincl],&xincl[iincl]);

		if(iwrit != 0) {
			for(iincl=1; iincl<=nincl; iincl++)
				fprintf(fp16,"%d %lf \n",lincl[iincl],xincl[iincl]);
			for(iincl=1; iincl<=nincl; iincl++)
				xincl[iincl]=xincl[iincl]*factor;
		}
	}
 //   printf("termino apoyos inclinados \n");

    ntotv = npnod*ngdln;
	return 0;
}

void CargaPropCatalogo( int nmats, int ndime, double** props, double*** arear,
						double*** arerf, double*** arerc, double*** arecr,
						int *ubiCat ){

	int imats, rolado, jj, iArm, iRF, iRC, iCon;


	for( imats = 1; imats <= nmats; imats++ ){
		if( ndime == 2 ) rolado = (int)props[imats][5]+0.1;
		if( ndime == 3 ) rolado = (int)props[imats][8]+0.1;
        printf("imats=%d   rolado=%d \n",imats, rolado);
		//Armaduras
		//area-inerciax-inerciay-???-???-????
		if( rolado == 0 ){
			iArm = ubiCat[imats];
			if( ndime == 2 ){
				jj = (int)props[imats][6]+0.5;
				props[imats][2] = arear[iArm][jj][1];
				props[imats][3] = arear[iArm][jj][2];
				props[imats][6] = jj;
			}
			if( ndime == 3 ){
				jj = (int)props[imats][9] + 0.5;
				props[imats][2] = arear[iArm][jj][1];
				props[imats][3] = arear[iArm][jj][2];
				props[imats][4] = arear[iArm][jj][3];
				props[imats][5] = arear[iArm][jj][2]+arear[iArm][jj][3];
				props[imats][6] = props[imats][1]/(2.6);
				props[imats][9] = jj;
			}
		}
		//Rolado Frio
	    //Á-Ix-Iy-Sx-Sy-Qx-Qy-Tx-Ty-H-B-D-R
		if( rolado == 1 ){
			iRF = ubiCat[imats];
			if( ndime == 2 ){
				jj = (int)props[imats][6]+0.5;
				props[imats][2] = arerf[iRF][jj][1];
				props[imats][3] = arerf[iRF][jj][2];
				props[imats][6] = jj;
			}
			if( ndime == 3 ){
				 jj = (int)props[imats][9] + 0.5;
				 props[imats][2] = arerf[iRF][jj][1];
				 props[imats][3] = arerf[iRF][jj][2];
				 props[imats][4] = arerf[iRF][jj][3];
				 props[imats][5] = arerf[iRF][jj][2] + arerf[iRF][jj][3];
				 props[imats][6] = props[imats][1]/(2.4);
				 props[imats][9] = jj;
			}
		}
		//Rolado Caliente
	    // E-Fy-A-Ix-Iy-J-Sx-Sy-Qx-Qy-tx-ty-nu
		if( rolado == 2 ){
			iRC = ubiCat[imats];
			if( ndime == 2 ){
				jj = (int)props[imats][6] + 0.5;
				props[imats][1] = arerc[iRC][jj][1];
				props[imats][2] = arerc[iRC][jj][3];
				props[imats][3] = arerc[iRC][jj][4];
				props[imats][6] = jj;
			}
			if(ndime ==3){
				jj= (int)props[imats][9]+0.5;
				props[imats][1] = arerc[iRC][jj][1];
				props[imats][2] = arerc[iRC][jj][3];
				props[imats][3] = arerc[iRC][jj][4];
				props[imats][4] = arerc[iRC][jj][5];
				props[imats][5] = arerc[iRC][jj][6];
				props[imats][6] = props[imats][1]/(2.0+2.0*arerc[iRC][imats][13]); // PENDIENTE
				props[imats][9] = jj;
			}
		}
		//Concreto
		//b-h-f'c,fy,numero de barras de acero y coordenadas
		if( rolado == 3 ){
			iCon = ubiCat[imats];
			if(ndime ==2){
				jj= (int)props[imats][6]+0.5;
				props[imats][2] = arecr[iCon][jj][1]*arecr[iCon][jj][2];
				props[imats][3] = arecr[iCon][jj][1]*arecr[iCon][jj][2]*
								  arecr[iCon][jj][2]*arecr[iCon][jj][2]/12.0;
				props[imats][6] = jj;
			}
			if(ndime ==3){
				jj= (int)props[imats][9]+0.5;
				props[imats][1] = sqrt(arecr[iCon][jj][3])*15000.0;
				props[imats][2] = arecr[iCon][jj][1]*arecr[iCon][jj][2];
				props[imats][3] = arecr[iCon][jj][1]*arecr[iCon][jj][2]*
								  arecr[iCon][jj][2]*arecr[iCon][jj][2]/12.0;
				props[imats][4] = arecr[iCon][jj][2]*arecr[iCon][jj][1]*
								  arecr[iCon][jj][1]*arecr[iCon][jj][1]/12.0;
				props[imats][5] = props[imats][3] + props[imats][4];
				props[imats][6] = props[imats][1]/(2.6);
				props[imats][9] = jj;
			}
		}
	}
}

void CargaPropCatOpt( int* matva, int ndime, int nmats, double** props,
					  double*** arear, double*** arerf, double*** arerc,
					  double*** arecr, int *ubiCat ){

	int imats, jj, rolado, iArm, iRF, iRC, iCon;
	for( imats = 1; imats <= nmats; imats++ ){
		jj = matva[imats];
		if( ndime == 2 ) rolado = (int)props[imats][5] + 0.1;
		if( ndime == 3 ) rolado = (int)props[imats][8] + 0.1;
		//Armaduras
		//area-inerciax-inerciay-???-???-????
		if( rolado == 0 ){
			iArm = ubiCat[imats];
			if( ndime == 2 ){
				props[imats][2] = arear[iArm][jj][1];
				props[imats][3] = arear[iArm][jj][2];
				props[imats][6] = jj;
			}
			if( ndime == 3 ){
				props[imats][2] = arear[iArm][jj][1];
				props[imats][3] = arear[iArm][jj][2];
				props[imats][4] = arear[iArm][jj][3];
				props[imats][5] = arear[iArm][jj][2] + arear[iArm][jj][3];
				props[imats][6] = props[imats][1]/(2.6);
				props[imats][9] = jj;
			}
		}
		//Rolado Frio
	    //Á-Ix-Iy-Sx-Sy-Qx-Qy-Tx-Ty-H-B-D-R
		if( rolado == 1 ){
			iRF = ubiCat[imats];
			if( ndime == 2 ){
				props[imats][2] = arerf[iRF][jj][1];
				props[imats][3] = arerf[iRF][jj][2];
				props[imats][6] = jj;
			}
			if( ndime == 3 ){
				 props[imats][2] = arerf[iRF][jj][1];
				 props[imats][3] = arerf[iRF][jj][2];
				 props[imats][4] = arerf[iRF][jj][3];
				 props[imats][5] = arerf[iRF][jj][2] + arerf[iRF][jj][3];
				 props[imats][6] = props[imats][1]/(2.4);
				 props[imats][9] = jj;
			}
		}
		//Rolado Caliente
	    // E-Fy-A-Ix-Iy-J-Sx-Sy-Qx-Qy-tx-ty-nu
		if( rolado == 2 ){
			iRC = ubiCat[imats];
			if( ndime == 2 ){
				props[imats][1] = arerc[iRC][jj][1];
				props[imats][2] = arerc[iRC][jj][3];
				props[imats][3] = arerc[iRC][jj][4];
				props[imats][6] = jj;
			}
			if(ndime ==3){
				props[imats][1] = arerc[iRC][jj][1];
				props[imats][2] = arerc[iRC][jj][3];
				props[imats][3] = arerc[iRC][jj][4];
				props[imats][4] = arerc[iRC][jj][5];
				props[imats][5] = arerc[iRC][jj][6];
				props[imats][6] = props[imats][1]/(2.0+2.0*arerc[iRC][imats][13]); // PENDIENTE
				props[imats][9] = jj;
			}
		}
		//Concreto
		//b-h-f'c,fy,numero d ebarras de acero y coordenadas
		if( rolado == 3 ){
			iCon = ubiCat[imats];
			if(ndime == 2 ){
				props[imats][2] = arecr[iCon][jj][1]*arecr[iCon][jj][2];
				props[imats][3] = arecr[iCon][jj][1]*arecr[iCon][jj][2]*
								  arecr[iCon][jj][2]*arecr[iCon][jj][2]/12.0;
				props[imats][6] = jj;
			}
			if(ndime ==3){
				props[imats][1] = sqrt(arecr[iCon][jj][3])*15000.0;
				props[imats][2] = arecr[iCon][jj][1]*arecr[iCon][jj][2];
				props[imats][3] = arecr[iCon][jj][1]*arecr[iCon][jj][2]*
								  arecr[iCon][jj][2]*arecr[iCon][jj][2]/12.0;
				props[imats][4] = arecr[iCon][jj][2]*arecr[iCon][jj][1]*
								  arecr[iCon][jj][1]*arecr[iCon][jj][1]/12.0;
				props[imats][5] = props[imats][3] + props[imats][4];
				props[imats][6] = props[imats][1]/(2.6);
				props[imats][9] = jj;
			}
		}
	}
}

int comp_cadenas( char *c1, char *c2 ){
	int igual, i, continua;
	if( c1 != NULL && c2 != NULL ){
		i = 1;
		igual    = 1;
		continua = 1;
		while( continua ){
			if( c1[i] == '\0' && c2[i] == '\0' ){
				continua = 0;
			}
			else{
				if( c1[i] != c2[i] ){
					igual    = 0;
					continua = 0;
				}
			}
			i++;
		}
	}
	return igual;
}

int quitaRepetidos_ListaCatalogos( int nmats, char** listaCatalogos ){
	int i, j, conta, contenida;
	conta = 1;
	for( i = 2; i <= nmats; i++ ){
		contenida = 0;
		for( j = 1; j <= conta; j++ ){
			if( comp_cadenas( listaCatalogos[j], listaCatalogos[i] ) ){
				contenida = 1;
				break;
			}
		}
		if( contenida == 0 ){
			conta++;
			j = 1;
			while( listaCatalogos[i][j] != '\0' ){
				listaCatalogos[conta][j] = listaCatalogos[i][j];
                j++;
			}
			listaCatalogos[conta][j] = '\0';
		}
	}
	return conta;
}

