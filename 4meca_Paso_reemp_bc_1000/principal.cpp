//*****************************
// Rutina Principal
//*****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "datos.h"
#include "estatico.h"
#include "ValoresPropios.h"
#include "fuerzas.h"
#include "main.h"
#include "posprogid.h"
#include "prelim.h"
#include "principal.h"
#include "raros.h"
#include "solver.h"
#include "memoria.h"
#include "Optimizador.h"

//******************************************************************************
//							Principal
//******************************************************************************
int Principal(int indso, int ishot, char* input1, char* output1, bool Optimiza,int prueba)
{
	// Finaliza 0 si no hubo ningun error.
	// Finaliza -1 si se hubo algun error de calculo.
	// Finaliza -2 si falto algun dato es erroneo.
	// Finaliza -3 si no se cuenta con la version completa (o la copia de evaluacion).
	// Finaliza -4 si el archivo no es de datos.
	// Finaliza -5 si no se pudo crear el archivo de resultados.

	int *nareas, *narerc, *narerf, *narecon, *matva, *matvm, *ubiCat, casoCarga, *npesp;
	cat_armadura cArmaduras;
	cat_acero_RF cAceroRF;
	cat_acero_RC cAceroRC;
	cat_concreto cConcr;
	double ***dtcon, *efi_max;
	double *angulo, **angulg;
	bool iwritModificado;
	//jacob.modif aumento de nom[4][256] a nom[5][256]
	char titulo[256], uLongitud[256], uFuerza[256], ver[10], **listaCatalogos;
	int iprob, ncarg, ndime, nelem, nevab, ngaus, ngdln, nmats, nfami, ncaso, nvpro, threads;
	int nnode, npnod, npres, nprob, nprop, ntens, ntipo, isale,ntotv, neqns, mkoun;
	int iwrit, nwktl, nreso,nincl, mnode, i, ival, iv, paso, *ntips, **inpre,  *nodpr;

    
	int **lnods,  *matnu,  *iffix, *maxad, *mhigh, **leqns,  *nodea,  *lreso,  *lincl, **indfu;

	double 	 *xlong, *wcarg, **coord, **carga,  **carpp, *despl, **props, **presc, **fuemp, **fuepp, **astif,
	*aslod,  *fixed,  *react,  *fuerc,  *stiff,  *vectp,  *vectu, *valor, *valom,
	*vecta, *vect1,  *deslo, **resor,  *xincl, **rigid, **girom, **fuerb,
	**girtm, **tempr, **vmatr, **fuefl, **xmasa, **vectr, ***srmat, ***smmat;
    double *tenlo, *mvalor, *aux, *aux1, ** mvector, ** vecto1, ** deslm, ** tenlm, ** vecth; // Se agrego para valores propios
    double ***vectores_fuerzas;           // Se agrego para almacenar los vectores de fuerza para cada caso de carga
	double beta,alfa;                     // Se agrego para optimizador
	double deltaf,dd,de,fact,penaliza;    // Se agrego para optimizador
	int jj,imats,Num_iter;                // Se agrego para optimizador

    char letrero[100]; //Contenido no importante, solo para lectura de letreros en archivos de datos.

    bool generaMatRigidez;
	int auxa=0;
    Optimiza = false;
	paso = 1; // este sirve para posrogid

	// Apertura de archivo de datos
	//int nchar = ;
	fp5  = fopen(input1, "rt");
	 for(i= (int) strlen(input1);i>=0;i--)
	   if(input1[i]=='.') {
			input1[i]='\0';
			break;
        }
	    else {
			input1[i]='\0';
		}
	strcpy(nomFlaviaMesh, input1);
	strcpy(nomFlaviaRes, input1);
//	strcpy(nomPos, input1);
//	strcpy(nomGir, input1);
	strcpy(nomOpt, input1);
	strcpy(nomLog, input1);

	strcat(nomFlaviaMesh, ".flavia.msh");         // para gid
	strcat(nomFlaviaRes, ".flavia.res");          // para gid
	strcat(nomPos, ".pos");
	strcat(nomGir, ".gir");


	// Crea los archivos de resultados
	fp16 = fopen(output1, "wt");
//  fp79 = fopen(nomGir,"wt");
//	fp9  = fopen(nomPos, "wt");

	// Inicia lectura del archivo de datos
	fscanf(fp5, "%s", titulo);
	fscanf(fp5, "%s", ver);
	fprintf(fp16,"%s\n", "ProMECA 1.0");

	ival  = 0;
	nprob = 1;
	for(iprob=1; iprob<=nprob; iprob++){

		// Lee datos del problema
        fscanf(fp5, "%s %s",letrero, titulo);
		fprintf(fp16, "%s\n", titulo);
		fscanf(fp5, "%s %s",letrero, uLongitud);
		fscanf(fp5, "%s %s", letrero, uFuerza);
		fprintf(fp16, "%s\t", uLongitud);
		fprintf(fp16, "%s\n", uFuerza);

		// Lee datos geometricos y del material:
		iv=Datos(npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,nfami,ncaso,
				 nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,threads,nvpro,lnods,inpre,nareas,narerf,narerc,
				 narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,iffix,indfu,nodea,maxad,mhigh,leqns,
				 coord,carga,carpp,presc,props,fuemp,fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,
				 fuerb,srmat,smmat,cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,
				 fixed,react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,vectp,
				 vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,tenlo,mvalor,
                 aux,aux1,mvector,vecto1,deslm,tenlm,vecth);
        
    
    
       // printf("termino de leer los datos");
   

		if( iwrit == 1 && Optimiza ) fprintf(fpLog,"Lectura de datos del problema = %d \n\n",Prob());


        // Calculo de cosenos dierectores de las barras //mg vectr es una matriz que contiene los cosenos directores.
		Prelim(nelem,ndime,ntipo,lnods,coord,xlong,vmatr,vectr,angulo);

		if(indso == 2) Linkin(nelem,nnode,nevab,ntotv,npres,ngdln,neqns,nwktl,mkoun,
							  lnods,nodea,nodpr,inpre,iffix,leqns,maxad,mhigh,
							  presc,fixed,stiff);

        // Lee el letrero de los casos de carga
        for( int i = 1; i <= 10; i++ ) {
            fscanf(fp5, "%s", letrero);
            //printf("%s \n", letrero);

        }

		generaMatRigidez = true;
		for( casoCarga = 1; casoCarga <= ncarg; casoCarga++ ){


            if(ncaso ==1) {
                Fuerzas(nelem,npnod,ndime,ntipo,nevab,nnode,ngdln,iwrit,casoCarga,lnods,matnu,indfu,ntips,coord,
                        xlong,props,fuerb,wcarg,carpp,fuepp,vectr,carga,fuemp,Optimiza,auxa);

                ival = Estatico( nelem,npnod,nevab,ndime,ntipo,ncaso,npres,iwrit,indso,isale,ngdln,nnode,
				                           ntotv,nincl,nreso,neqns,nwktl,ishot,lnods,matnu,inpre,iffix,ntips,maxad,
				                           nodea,nodpr,leqns,lincl,lreso,react,coord,presc,fixed,despl,aslod,stiff,
				                           srmat,props,xlong,vectr,rigid,carpp,fuepp,carga,fuemp,girom,girtm,tempr,
				                           fuerc,valor,cArmaduras.arear,cAceroRF.arerf,cAceroRC.arerc,cConcr.arecr,
				                           cConcr.dtcon,indfu,fuerb,fuefl,wcarg,xincl,resor,astif,vectp,vectu,vect1,
				                           deslo,ubiCat,efi_max,generaMatRigidez,Optimiza);
                nvpro=paso;
                generaMatRigidez = false;
            }
            if(ncaso ==2){
                     printf("Inicia Calculo de Valores y Vectores Propios= %d \n\n",Prob());
                    ival =ValoresPropios( nelem,npnod,nevab,ndime,ntipo,ncaso,npres,iwrit,indso,isale,ngdln,nnode,
                                          ntotv,nincl,nreso,neqns,nwktl,ishot,lnods,matnu,inpre,iffix,ntips,maxad,
                                          nodea,nodpr,leqns,lincl,lreso,react,coord,presc,fixed,despl,aslod,stiff,
                                          srmat,props,xlong,vectr,rigid,carpp,fuepp,carga,fuemp,girom,girtm,tempr,
                                          mkoun,mhigh,fuerc,valor,fuerb,fuefl,wcarg,xincl,resor,astif,vectp,vectu,vect1,
                                          deslo,tenlo,aux,aux1,smmat,xmasa,vecta,nvpro,mvalor,mvector,vecto1,threads,
                                          deslm,tenlm,vecth,generaMatRigidez);
                generaMatRigidez = false;
            }
            if(ncaso ==3){
                strcat(nomOpt, "Opt.txt");  fpOpt =fopen(nomOpt,"wt");
                strcat(nomLog, "Log.txt");  fpLog =fopen(nomLog,"wt");
                ival = Optimizador(ishot, npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,
                                                  nfami,ncaso, nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,lnods,inpre,
                                                  nareas,narerf,narerc,narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,
                                                  iffix,indfu,nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
                                                  fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,fuerb,srmat,smmat,
                                                  cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,fixed,
                                                  react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
                                                  vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,stiff,
                                                  neqns, nwktl ,prueba);
              }

		// Post-proceso con GID
		PosproGid(npnod,nelem,ndime,nnode,ngdln,ntotv,nmats,ncaso,nvpro,lnods,matnu,matva,coord,
				  props,despl,fuerc,aslod,carga, titulo,mvector);
        
        }

		// Libera memoria
		LiberaMemoria( npnod, nelem, ncarg, ngdln, nmats, ndime, nevab, indso, lnods,
		               inpre, nareas, narerf, narerc, narecon, lreso, lincl, ntips,
		               nodpr, matnu, matva, matvm, iffix, indfu, nodea, maxad, mhigh,
		               leqns, coord, carga, carpp, presc, props, fuemp, fuepp, rigid,
		               xmasa, vmatr, girom, girtm, tempr, fuefl, fuerb, srmat, smmat,
		               resor, xincl, fuerc, aslod, despl, fixed, react, vecta, vect1,
		               deslo, xlong, valor, valom, angulo, angulg, wcarg, vectr, astif,
		               vectp, vectu, listaCatalogos, ubiCat, vectores_fuerzas, cArmaduras.nombre,
		               cArmaduras.tamCat, cArmaduras.arear, cAceroRF.nombre, cAceroRF.tamCat,
		               cAceroRF.arerf, cAceroRC.nombre, cAceroRC.tamCat, cAceroRC.arerc,
		               cConcr.nombre, cConcr.tamCat, cConcr.arecr, cConcr.dtcon,
		               cArmaduras.nCat, cAceroRF.nCat, cAceroRC.nCat, cConcr.nCat,
		               stiff, npesp, efi_max );

     }
//printf("Ejecucion del programa = %d \n\n",Prob());
	return ival;

}

//******************************************************************************
void Actualiza(int nmats, int* matva, int* matvm, double* valor, double* valom)
{
	int imats;
	//printf("aceptado \n");
	for(imats=1; imats<=nmats; imats++)
 	    matvm[imats]=matva[imats];
	for(imats=1; imats<=10; imats++)
		valom[imats]=valor[imats];
}
//******************************************************************************
void Recupera(int nmats, int* matva, int* matvm, double* valor, double* valom)
{
	int imats;
	for(imats=1; imats<=nmats; imats++)
 	    matva[imats]=matvm[imats];
	for(imats=1; imats<=10; imats++)
		valor[imats]=valom[imats];
}
//******************************************************************************
void Inicializa_Materiales( int nmats, int ndime, double** props, int* matva,
							int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
							int *Con_tamCat, int *ubiCat ){
	// Inicia los materiales con la seccion ultima, la de mayor area transversal,
	// de cada catalogo. En el catalogo deben estar ordenadas de menor a mayor
	// en cuanto a area.
	int imats, rolado, iCat;
	for( imats = 1; imats <= nmats; imats++ ){
		//if( ndime ==2 ) rolado = (int)props[imats][5] + 0.5;
		//if( ndime ==3 )
        rolado = (int)props[imats][8] + 0.5;
		iCat = ubiCat[imats];
		switch( rolado ){
			case 0:
				matva[imats] = Arm_tamCat[iCat];
				break;
			case 1:
				matva[imats] = RF_tamCat[iCat];
				break;
			case 2:
				matva[imats] = RC_tamCat[iCat];
				break;
			case 3:
				matva[imats] = Con_tamCat[iCat];
				break;
			default:
				break;
		}
	}
}
//******************************************************************************
void inicializa_array1D( int n, double *v, double cte ){
	for( int i = 1; i <= n; i++ ) v[i] = cte;
}
//******************************************************************************
void escribe_efi_max( int nelem, double *efi_max ){
	for( int ielem = 1; ielem <= nelem; ielem++ )
		fprintf(fp16,"elemento = %d eficiencia_max = %f\n", ielem, efi_max[ielem] );
}
//******************************************************************************
//void Recosido_Simulado(int Num_iter)
  /*
	int jj,imats, ival;
	double dd,de,beta,deltaf,fact,penaliza;

	penaliza=1.e07;
	for(jj=1; jj<=Num_iter; jj++){
    for(imats=1; imats<=nmats;imats++)
    printf("matva[%d]= %d \t",imats,matva[imats]);
    printf("\n");
    Genera_Nueva_Poblacion(nmats,ndime,matva,valmax,props);
    CargaPropCatOpt(matva);
    ival =Estatico();
    dd = valor[2]-10.0;
    de = valor[3]-1.0;
    valor[4]=valor[1];
    if (dd >0)valor[4]+=dd*penaliza;
    if (de >0)valor[4]+=de*penaliza;
    if (beta < 1.e50) beta *= alfa ;
    //	printf("caso de estudio = %d peso =%lf  desp. Max=%lf  eficiencia=%lf   funcion=%lf \n",jj,valor[1],valor[2],valor[3],valor[4]);
    //	printf("Materiales 1=%d    2=%d    3=%d    4=%d \n",matva[1],matva[2],matva[3],matva[4]);
    deltaf = valor[4] - valom[4];
    if (deltaf < 0) Actualiza() ;
    else{
    //	printf(" %d deltaf =%lf valor=%lf\n",jj,deltaf,valor[4]);
    deltaf *= beta ;
    deltaf = exp(-deltaf) ;
    fact=(float)rand()/(float)RAND_MAX;
    if (fact < deltaf) {
    Actualiza() ;
    //printf("entro \n");
    //printf("deltaf =%lf \n",deltaf);
    //printf("fact=%lf deltaf=%lf beta=%lf  itera =%d\n",fact,deltaf,beta,jj);
    }
    else Recupera();
    }
	}    */



