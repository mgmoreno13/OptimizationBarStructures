//*****************************
// Rutina Estatico
//*****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "eficiencias.h"
#include "estatico.h"
#include "main.h"
#include "raros.h"
#include "rigidez.h"
#include "solver.h"
#include "tensiones.h"
#include "fuerzas.h"


//************************************************
int Estatico(int nelem, int npnod, int nevab, int ndime, int ntipo, int ncaso, int npres, int iwrit,
			 int indso, int isale, int ngdln, int nnode, int ntotv, int nincl, int nreso, int neqns,
			 int nwktl, int ishot, int** lnods, int* matnu, int** inpre,  int* iffix,int* ntips,
			 int* maxad, int* nodea, int* nodpr, int** leqns, int* lincl, int* lreso,double* react,
			 double** coord, double** presc, double* fixed, double* despl, double* aslod,double* stiff,
			 double*** srmat, double** props, double* xlong, double** vectr, double** rigid,
			 double** carpp, double** fuepp, double** carga, double** fuemp, double** girom,
			 double** girtm, double** tempr,double* fuerc, double* valor, double ***arear,
			 double ***arerf, double ***arerc, double ***arecr, double ****dtcon, int** indfu,
			 double** fuerb, double** fuefl, double* wcarg, double* xincl, double** resor,
			 double** astif, double* vectp, double* vectu, double* vect1, double* deslo, int *ubiCat,
			 double *efi_max, bool generaMatRigidez, bool Optimiza )
{

    int ncasa=1;

	/******************************************************************************
	 Solucion a problema Estatico Lineal
	 ******************************************************************************/
	if(iwrit ==1 && Optimiza == false) {
		printf("\n\n");
		printf("    *******************************************\n");
		printf("    *  Calculo de un problema Estatico Lineal *\n");
		printf("    *******************************************\n");
		printf("\n\n");
	}


	if(generaMatRigidez == true){
		// Calcula la matriz de rigidez de los elementos:
		Rigimat(nelem,nevab,ndime,ntipo,isale,matnu,ntips,srmat,props,xlong,vectr,rigid,
			    girom,girtm,tempr);

		if(iwrit ==1) printf("Calculo de la matriz de rigidez = %d \n\n",Prob()) ;

		// Ensambla y factoriza la matriz de rigidez
		if(indso == 2) {
			Addban(nelem,ngdln,ndime,nnode,nevab,isale,nincl,nwktl,nreso, //mg Montaje del vector de rigidez total.
				   lnods,leqns,maxad,nodea,lincl,lreso,iffix,girom,
				   girtm,tempr,rigid,stiff,srmat,xincl,resor);
			if(iwrit ==1) printf("Ensamblando matriz de rigidez = %d \n\n",Prob()) ;

			if(Decomp(neqns,ishot,maxad,stiff)) {//Factoriza (l) * (d) * (l) transposición de la matriz de rigidez.  (stiff es rigidez)
				fclose(fp1);
				CierraLimpia();
				return -1;
			}
			if(iwrit ==1) printf("Factorizando matriz de rigidez = %d \n\n",Prob()) ;
		}
	}

	switch(indso) {
		case 1:      /* Ensambla y resuelve las ecuaciones de rigidez
	 por eliminacion gaussiana: */
			if(Solucion(npnod,nelem,ncasa,ngdln,nnode,npres,ntipo,ndime,nincl,
						nreso,nevab,ntotv,iwrit,isale,lnods,nodpr,inpre,
						iffix,lreso,lincl,ntips,presc,despl,
						fixed,astif,aslod,carga,srmat,resor,
						girom,girtm,tempr,rigid,xincl,react))return -1;
			break;
		case 2:
			Profile(nelem,ndime,npnod,ngdln,ncasa,neqns,iwrit,ntotv,nincl,nevab,isale,// mg Calcula los desplazamientos incrementales... como indso=2 (establecido en main.cpp)
					lnods,nodea,maxad,iffix,lincl,carga,despl,
					aslod,fixed,stiff,srmat,react,
					girom,girtm,tempr,rigid,xincl);
			break;
		case 3:
			Gradc(nelem,npnod,npres,ngdln,ntotv,nnode,ncasa,iwrit,ndime,
				  nevab,isale,nincl,nreso,lnods,inpre,nodpr,iffix,
				  nodea,lincl,lreso,fixed,presc,despl,aslod,
				  carga,vectp,vectu,girom,girtm,tempr,
				  rigid,xincl,resor,srmat,react,vect1,deslo);
			break;
	}
	if(iwrit ==1) printf("Resolviendo sistema de ecuaciones = %d \n\n",Prob()) ;

	// Calcula los esfuerzos en los elementos:
	Tensiones(nelem,nnode,ngdln,ndime,ntipo,nevab,ntotv,iwrit,isale,nincl,matnu,lnods,lincl,//mg  Calcula las tensiones y esfuerzos en los puntos de Gauss.
			  props,despl,srmat,fuerc,fuemp,fuefl,vectr,xlong,girom,girtm,tempr,rigid,xincl,valor);//mg <---aqui se llena el vector "valor"
   // printf("%lf   %lf  %lf %lf \n",valor[1],valor[2],valor[3],valor[4]);

    if(iwrit ==1) printf("Calculo de tensiones = %d \n\n",Prob()) ;

    
    
	// Calcula las Eficiencias en los elementos:
    //Muttio.Agregado de nevab
	Eficiencias(nelem,ndime,ntipo,matnu,ntips,props,fuerc,valor,xlong,arear,arerf,
				arerc,arecr,dtcon,indfu,fuerb,fuefl,wcarg,iwrit,ubiCat,efi_max,nevab);
   // printf("valor1= %lf   valor2 =%lf  valor3=%lf valor4= %lf \n",valor[1],valor[2],valor[3],valor[4]);
	if(iwrit ==1) printf("Calculo de eficiencias = %d \n\n",Prob());

	//LiberaMemoria();
	return 0;

}

