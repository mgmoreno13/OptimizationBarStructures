//*****************************
// Rutina ValoresPropios
//*****************************
#pragma hdrstop
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include "ValoresPropios.h"
#include "estatico.h"
#include "main.h"
#include "raros.h"
#include "rigidez.h"
#include "solver.h"
#include "tensiones.h"
//************************************************

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
                   double** tenlm, double** vecth,bool generaMatRigidez)
{
	int imasa,ivpti,indic=2,ielem,ievab;
	/******************************************************************************
	 Solucion a problema de Modos de vibracion
	 ******************************************************************************/
	printf("\n\n");
	printf("    *******************************************\n");
	printf("    *  Calculo de Modos de vibracion          *\n");
	printf("    *******************************************\n");
	printf("\n\n");
    imasa=1;
    ivpti=2;

	// Calcula la matriz de rigidezy de masa de los elementos:
	
    printf("valores propios=%d \n",nvpro);
    
    if(generaMatRigidez == true){
	Rigimat(nelem,nevab,ndime,ntipo,isale,matnu,ntips,srmat,props,xlong,vectr,rigid,
		    girom,girtm,tempr);
	if(imasa == 1) Masaconsis(nelem,nevab,ndime,ntipo,isale,matnu,ntips,smmat,props,xlong,vectr,xmasa,
							 girom,girtm,tempr);
	if(imasa == 2) Masaconsen(nelem,nevab,ndime,ntipo,isale,matnu,ntips,smmat,props,xlong,vectr,xmasa,  
						     girom,girtm,tempr);
	printf("Calculo de la matriz de rigidez y de masa= %d \n\n",Prob()) ;
    }
    
    if(ivpti == 1) VP_Bajos(nelem,ngdln,ndime,nnode,npnod,ncaso,ntipo,iwrit,nevab,isale,nincl,nwktl,nreso,indso,ntotv,npres,
                           neqns,mkoun,ishot,nvpro,nodpr,inpre,mhigh,ntips,lnods,leqns,maxad,nodea,lincl,lreso,iffix,
                           fixed,stiff,xincl,resor,mvalor,despl,deslo,tenlo,astif,aslod,react,vecta,vectp,vectu,vect1,aux,
                           aux1,carga,girom,girtm,tempr,rigid,presc,mvector,vecto1,srmat,smmat,threads,deslm,tenlm,vecth);
    else if (ivpti ==2) VP_Alto(nelem,ngdln,ndime,nnode,npnod,ncaso,ntipo,iwrit,nevab,isale,nincl,nwktl,nreso,indso,ntotv,npres,
                                neqns,mkoun,ishot,nvpro,nodpr,inpre,mhigh,ntips,lnods,leqns,maxad,nodea,lincl,lreso,iffix,
                                fixed,stiff,xincl,resor,mvalor,despl,deslo,tenlo,astif,aslod,react,vecta,vectp,vectu,vect1,aux,
                                aux1,carga,girom,girtm,tempr,rigid,presc,mvector,vecto1,srmat,smmat,threads,deslm,tenlm,vecth);
    
    
	return indic;
	
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
int VP_Bajos(int nelem,int ngdln,int ndime,int nnode,int npnod,int ncaso,int ntipo,int iwrit,int nevab,int isale,
             int nincl,int nwktl,int nreso,int indso,int ntotv,int npres,int neqns,int mkoun,int ishot,int nvpro,
             int *nodpr, int** inpre, int *mhigh, int *ntips,int** lnods,int** leqns,int* maxad,
             int* nodea,int* lincl,int* lreso,int* iffix,double* fixed,double* stiff,double* xincl,
             double** resor,double* mvalor, double* despl,double *deslo, double *tenlo, double** astif,
             double* aslod,double *react, double* vecta, double *vectp, double* vectu, double* vect1, double* aux,
             double* aux1,double **carga, double** girom,double** girtm,double** tempr,double** rigid,
             double** presc,double** mvector, double**vecto1, double*** srmat,double*** smmat,int threads, double** deslm,
             double** tenlm, double** vecth)
{
    int ivpro, itotv, itera,tmaxc,indic;
    double toler, vectn;
    double vecto,numer, denom, diffe;
	toler=1e-20;
	tmaxc=1000000;
 
if(indso == 2) {
    Linkin(nelem,nnode,nevab,ntotv,npres,ngdln,neqns,nwktl,mkoun,
                          lnods,nodea,nodpr,inpre,iffix,leqns,maxad,mhigh,
                          presc,fixed,stiff);
    Addban(nelem,ngdln,ndime,nnode,nevab,isale,nincl,nwktl,nreso,
           lnods,leqns,maxad,nodea,lincl,lreso,iffix,girom,
           girtm,tempr,rigid,stiff,srmat,xincl,resor);
    printf("Ensamblando matriz de rigidez = %d \n\n",Prob()) ;
    if(Decomp(neqns,ishot,maxad,stiff)) {
        CierraLimpia();
        return -1;
    }
    printf("Factorizando matriz de rigidez = %d \n\n",Prob()) ;
}
 
 
for(ivpro=1; ivpro<=nvpro; ivpro++){
    mvalor[ivpro]=0;
 #pragma omp parallel for default(shared)
    for(itotv=1; itotv<=ntotv; itotv++){
        mvector[ivpro][itotv]=0;
         vecto1[ivpro][itotv]=0;
       }
    }
for(ivpro=1; ivpro<=nvpro; ivpro++){
   // printf("Calculando Vector Propio Numero %d   en %d  seg.\n",ivpro,Prob());
    Inicia(ntotv,npres,ngdln,indso,inpre,nodpr,iffix,presc,despl,vecta,fixed);

    vecto=1e30;
    itera=0;
    do {
        itera++;
        Derecho(nelem,ntotv,ndime,ngdln,nnode,nevab,isale,nincl,lnods,iffix,lincl,girom,girtm,tempr,rigid,xincl,
                deslo,tenlo,vecta,aslod,despl,aux,aux1,smmat,mvalor,mvector,ivpro,threads,deslm,tenlm,vecth);
		/* Ensambla y resuelve las ecuaciones */
        switch(indso) {
            case 1:
                Solucion(npnod,nelem,ncaso,ngdln,nnode,npres,ntipo,ndime,nincl,
                         nreso,nevab,ntotv,iwrit,isale,lnods,nodpr,inpre,
                         iffix,lreso,lincl,ntips,presc,despl,
                         fixed,astif,aslod,carga,srmat,resor,
                         girom,girtm,tempr,rigid,xincl,react);
                break;
            case 2:
                Profile(nelem,ndime,npnod,ngdln,ncaso,neqns,iwrit,ntotv,nincl,nevab,isale,
                        lnods,nodea,maxad,iffix,lincl,carga,despl,
                        aslod,fixed,stiff,srmat,react,
                        girom,girtm,tempr,rigid,xincl);
                break;
            case 3:
                Gradc(nelem,npnod,npres,ngdln,ntotv,nnode,ncaso,iwrit,ndime,
                      nevab,isale,nincl,nreso,lnods,inpre,nodpr,iffix,
                      nodea,lincl,lreso,fixed,presc,despl,aslod,
                      carga,vectp,vectu,girom,girtm,tempr,
                      rigid,xincl,resor,srmat,react,vect1,deslo);
                break;
        }
        numer=1.0;
        denom=0;
#pragma omp parallel for default(shared)reduction(+:denom)
        for(itotv=1; itotv<=ntotv; itotv++)
            denom+=despl[itotv]*despl[itotv];
        vectn=denom;
        diffe=fabs(vectn-vecto);
       // printf("itera =%d vectn =%e   vecto=%e numer=%e   diffe=%e toler =%e\n",itera,vectn,vecto,numer,diffe,toler);
        denom=sqrt(denom);
#pragma omp parallel for default(shared)
        for(itotv=1; itotv<=ntotv; itotv++)
            vecta[itotv]=despl[itotv]/denom;
        if(diffe<toler) break;
        vecto=vectn;
    }while(itera <= tmaxc);
    MasaVect(nelem,nnode,ngdln,nevab,ntotv,iffix,lnods,deslo,tenlo,vecta,aslod,smmat,threads,deslm,tenlm,vecth);
    denom=0;
#pragma omp parallel for default(shared)reduction(+:denom)
    for(itotv=1; itotv<=ntotv; itotv++)
        denom+=vecta[itotv]*aslod[itotv];
    denom=sqrt(denom);
#pragma omp parallel for default(shared)
    for(itotv=1; itotv<=ntotv; itotv++){
        vecta[itotv]=vecta[itotv]/denom;
        aslod[itotv]=aslod[itotv]/denom;
    }
    MasaVect(nelem,nnode,ngdln,nevab,ntotv,iffix,lnods,deslo,tenlo,vecta,despl,srmat,threads,deslm,tenlm,vecth);
    vectn=0;
#pragma omp parallel for default(shared)reduction(+:vectn)
    for(itotv=1; itotv<=ntotv; itotv++){
        vectn+=vecta[itotv]*despl[itotv];
    }
    Almacena(vectn,vecta,aslod,mvalor,mvector,vecto1,ivpro,ntotv);
    printf("itera=%d \t valorP[%d]=%lf \t diffe=%e  \t  en %d \t seg.\n",itera,ivpro,vectn,diffe,Prob());
    if(itera < tmaxc) indic=0;
    else indic =5;

}
printf("resultados finales \n");
printf("valores propios \n");
for (ivpro=1; ivpro<=nvpro; ivpro++)
printf("valor[%d]=%e \n",ivpro,mvalor[ivpro]);
printf("\n");

/*
    //propuctos punto de vectores
    for (ivpro=1; ivpro<=nvpro; ivpro++){
        for (int jvpro=1; jvpro<=nvpro; jvpro++){
           denom=0;
           for(itotv=1; itotv<=ntotv; itotv++)
               denom+=mvector[jvpro][itotv]*vecto1[ivpro][itotv];
            printf("ivpro=%d  jvpro=%d prod=%e \n",ivpro,jvpro,denom);
        }
    }
    printf("\n\n\n");
    for (ivpro=1; ivpro<=nvpro; ivpro++){
        MasaVect(nelem,nnode,ngdln,nevab,ntotv,iffix,lnods,deslo,tenlo,vecto1[ivpro],aslod,srmat);
        for (int jvpro=1; jvpro<=nvpro; jvpro++){
            denom=0;
            for(itotv=1; itotv<=ntotv; itotv++)
                denom+=aslod[itotv]*vecto1[jvpro][itotv];
            printf("ivpro=%d  jvpro=%d prod=%e \n",ivpro,jvpro,denom);
        }
    }
 */
            return indic;
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
int VP_Alto(int nelem,int ngdln,int ndime,int nnode,int npnod,int ncaso,int ntipo,int iwrit,int nevab,int isale,
             int nincl,int nwktl,int nreso,int indso,int ntotv,int npres,int neqns,int mkoun,int ishot,int nvpro,
             int *nodpr, int** inpre, int *mhigh, int *ntips,int** lnods,int** leqns,int* maxad,
             int* nodea,int* lincl,int* lreso,int* iffix,double* fixed,double* stiff,double* xincl,
             double** resor,double* mvalor, double* despl,double *deslo, double *tenlo, double** astif,
             double* aslod,double *react, double* vecta, double *vectp, double* vectu, double* vect1, double* aux,
             double* aux1,double **carga, double** girom,double** girtm,double** tempr,double** rigid,
             double** presc,double** mvector, double** vecto1, double*** srmat,double*** smmat, int threads,
             double** deslm,double** tenlm, double** vecth)
{
    int ivpro, itotv, itera,tmaxc,indic;
    double vecto, numer, valrr, toler, vectn,diffe;
    toler=1e-20;
    tmaxc=1000000;
    if(indso == 2) {
        Linkin(nelem,nnode,nevab,ntotv,npres,ngdln,neqns,nwktl,mkoun,
               lnods,nodea,nodpr,inpre,iffix,leqns,maxad,mhigh,
               presc,fixed,stiff);
        Addban(nelem,ngdln,ndime,nnode,nevab,isale,nincl,nwktl,nreso,
               lnods,leqns,maxad,nodea,lincl,lreso,iffix,girom,
               girtm,tempr,rigid,stiff,smmat,xincl,resor);
        printf("Ensamblando matriz de masa = %d \n\n",Prob()) ;
        
        if(Decomp(neqns,ishot,maxad,stiff)) {
            CierraLimpia();
            return -1;
        }
        printf("Factorizando matriz de masa = %d \n\n",Prob()) ;
    }
    for(ivpro=1; ivpro<=nvpro; ivpro++){
        mvalor[ivpro]=0;
#pragma omp parallel for default(shared)
        for(itotv=1; itotv<=ntotv; itotv++){
            mvector[ivpro][itotv]=0;
             vecto1[ivpro][itotv]=0;
        }
    }
    for(ivpro=1; ivpro<=nvpro; ivpro++){
        Inicia(ntotv,npres,ngdln,indso,inpre,nodpr,iffix,presc,despl,vecta,fixed);
        vecto=1e-30;
        itera=0;
        do {
            itera++;
            /* calcula lado derecho M-fi */
            Derecho(nelem,ntotv,ndime,ngdln,nnode,nevab,isale,nincl,lnods,iffix,lincl,girom,girtm,tempr,rigid,xincl,
                    deslo,tenlo,vecta,aslod,despl,aux,aux1,srmat,mvalor,mvector,ivpro,threads,deslm,tenlm,vecth);
            switch(indso) {
                case 1:
                    Solucion(npnod,nelem,ncaso,ngdln,nnode,npres,ntipo,ndime,nincl,
                             nreso,nevab,ntotv,iwrit,isale,lnods,nodpr,inpre,
                             iffix,lreso,lincl,ntips,presc,despl,
                             fixed,astif,aslod,carga,smmat,resor,
                             girom,girtm,tempr,rigid,xincl,react);
                    break;
                case 2:
                    Profile(nelem,ndime,npnod,ngdln,ncaso,neqns,iwrit,ntotv,nincl,nevab,isale,
                            lnods,nodea,maxad,iffix,lincl,carga,despl,
                            aslod,fixed,stiff,smmat,react,
                            girom,girtm,tempr,rigid,xincl);
                    break;
                case 3:
                    Gradc(nelem,npnod,npres,ngdln,ntotv,nnode,ncaso,iwrit,ndime,
                          nevab,isale,nincl,nreso,lnods,inpre,nodpr,iffix,
                          nodea,lincl,lreso,fixed,presc,despl,aslod,
                          carga,vectp,vectu,girom,girtm,tempr,
                          rigid,xincl,resor,smmat,react,vect1,deslo);
                    break;
            }
            numer=0;
#pragma omp parallel for default(shared)reduction(+:numer)
            for(itotv=1; itotv<=ntotv; itotv++)
                numer+=despl[itotv]*despl[itotv];
            vectn=fabs(numer);
            diffe=fabs(vectn-vecto);
            //printf("itera =%d vectn =%e   vecto=%e numer=%e   diffe=%e toler =%e\n",itera,vectn,vecto,numer,diffe,toler);
            valrr=sqrt(numer);
#pragma omp parallel for default(shared)
            for(itotv=1; itotv<=ntotv; itotv++)
                vecta[itotv]=despl[itotv]/valrr;
            if(diffe<toler) break;
            vecto=vectn;
        }while(itera <= tmaxc);
        MasaVect(nelem,nnode,ngdln,nevab,ntotv,iffix,lnods,deslo,tenlo,vecta,aslod,srmat,threads,deslm,tenlm,vecth);
        valrr=0;
#pragma omp parallel for default(shared)reduction(+:valrr)
        for(itotv=1; itotv<=ntotv; itotv++)
            valrr+=vecta[itotv]*aslod[itotv];
        valrr=sqrt(valrr);
#pragma omp parallel for default(shared)
        for(itotv=1; itotv<=ntotv; itotv++){
            vecta[itotv]=vecta[itotv]/valrr;
            aslod[itotv]=aslod[itotv]/valrr;
        }
        MasaVect(nelem,nnode,ngdln,nevab,ntotv,iffix,lnods,deslo,tenlo,vecta,despl,smmat,threads,deslm,tenlm,vecth);
        vectn=0;
#pragma omp parallel for default(shared)reduction(+:vectn)
        for(itotv=1; itotv<=ntotv; itotv++)
            vectn+=vecta[itotv]*despl[itotv];
        vectn=1.0/vectn;
        Almacena(vectn,vecta,aslod,mvalor,mvector,vecto1,ivpro,ntotv);
        printf("itera=%d \t valorP[%d]=%lf \t diffe=%e  \t  en %d \t seg.\n",itera,ivpro,vectn,diffe,Prob());
        if(itera < tmaxc) indic=0;
        else indic =5;
    }
    printf("resultados finales \n");
    printf("valores propios \n");
    for (ivpro=1; ivpro<=nvpro; ivpro++)
        printf("valor[%d]=%e \n",ivpro,mvalor[ivpro]);
    printf("\n");
/*
    //propuctos punto de vectores
    for (ivpro=1; ivpro<=nvpro; ivpro++){
        for (int jvpro=1; jvpro<=nvpro; jvpro++){
            denom=0;
            for(itotv=1; itotv<=ntotv; itotv++)
                denom+=mvector[jvpro][itotv]*vecto1[ivpro][itotv];
            printf("ivpro=%d  jvpro=%d prod=%e \n",ivpro,jvpro,denom);
        }
    }
    printf("\n\n\n");
    for (ivpro=1; ivpro<=nvpro; ivpro++){
        MasaVect(nelem,nnode,ngdln,nevab,ntotv,iffix,lnods,deslo,tenlo,vecto1[ivpro],aslod,smmat);
        for (int jvpro=1; jvpro<=nvpro; jvpro++){
            denom=0;
            for(itotv=1; itotv<=ntotv; itotv++)
                denom+=aslod[itotv]*vecto1[jvpro][itotv];
            if(jvpro==ivpro)denom=1/denom;
            printf("ivpro=%d  jvpro=%d prod=%e \n",ivpro,jvpro,denom);
        }
    }
 */
    return indic;
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Inicia(int ntotv,int npres, int ngdln, int indso, int **inpre,int *nodpr, int *iffix, double **presc,
			double *despl, double *vecta, double *fixed)
/******************************************************************************
 Calcula matriz de rigidez de cada elemento.
 ******************************************************************************/
{
	int itotv,ipres,igdln,nloca,mtotv;
    double fact;
	if (indso ==3){
#pragma omp parallel for default(shared)
		for(itotv=1; itotv<=ntotv; itotv++) {
			fixed[itotv]=0.0;
            iffix[itotv]=0;
	    }
		
		for(ipres=1; ipres<=npres; ipres++) {
			for(igdln=1; igdln<=ngdln; igdln++) {
				if(inpre[ipres][igdln] == 1) {
					nloca = (nodpr[ipres]-1)*ngdln+igdln;
					iffix[nloca] = 1;
					fixed[nloca] = presc[ipres][igdln];
				}
				
			}
		}
	}

	mtotv=0;
	for(itotv=1; itotv<=ntotv; itotv++)
        if(iffix[itotv] ==0) mtotv++;
   // printf("mtotv= %d\n",mtotv);
    fact=1.0; ///sqrt((double)mtotv);
    for(itotv=1; itotv<=ntotv; itotv++){
		if(iffix[itotv] ==0)despl[itotv] =vecta[itotv] = fact;
        else despl[itotv] =vecta[itotv] = 0;
   // printf("itotv=%d ,iffix[itotv]= %d,despl[itotv]=%lf,vecta[itotv]=%lf \n",itotv,iffix[itotv],despl[itotv],vecta[itotv]);
    }
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Derecho(int nelem, int ntotv, int ndime, int ngdln, int nnode, int nevab, int isale, int nincl, int** lnods,  
			 int* iffix, int* lincl,double** girom, double** girtm, double** tempr, double** rigid, double* xincl,
             double* deslo, double *tenlo,double* vecta, double* aslod, double* despl, double* aux, double *aux1,
             double*** smmat, double* mvalor,double** mvector,int ivpro,int threads,double** deslm,double** tenlm,
             double** vecth)
/******************************************************************************
 Calcula matriz de masa por vector iterativo
 ******************************************************************************/
{
	int    itotv,jvpro;
    double xnorm;
// inicializa producto matriz-vector
 #pragma omp parallel for default(shared)
    for(itotv=1; itotv<=ntotv; itotv++)
        aux[itotv]=0.0;
    MasaVect(nelem,nnode,ngdln,nevab,ntotv,iffix,lnods,deslo,tenlo,vecta,aslod,smmat,threads,deslm,tenlm,vecth);
        if(ivpro >1){
        // calcula contribuciones de cada vector ya conocido
        for(jvpro=1; jvpro<=(ivpro-1); jvpro++){
            xnorm=0;
            for(itotv=1; itotv<=ntotv; itotv++)
                xnorm+=vecta[itotv]*mvector[jvpro][itotv];
            for(itotv=1; itotv<=ntotv; itotv++)
                aux[itotv]+=xnorm*mvector[jvpro][itotv];
        }
        }
//    printf("paso 1 \n");
    xnorm=0;
#pragma omp parallel for default(shared)reduction(+:xnorm)
    for(itotv=1; itotv<=ntotv; itotv++){
        aslod[itotv]=aslod[itotv]-aux[itotv];
        xnorm+= aslod[itotv]* aslod[itotv];
    }
//    printf("xnorm= %d \n",xnorm);
    xnorm=sqrt(xnorm);
#pragma omp parallel for default(shared)
    for(itotv=1; itotv<=ntotv; itotv++){
        aslod[itotv]=aslod[itotv]/xnorm;
//    printf("aslod[%d]=%lf \n",itotv,aslod[itotv]);
    }
        }
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Almacena(double vectn,double* vecta, double* aslod, double* mvalor,double** mvector, double** vecto1, int ivpro, int ntotv)
{
    int itotv;
    mvalor[ivpro]=vectn;
#pragma omp parallel for default(shared)
    for(itotv=1; itotv<=ntotv; itotv++){
        mvector[ivpro][itotv]=aslod[itotv];
         vecto1[ivpro][itotv]=vecta[itotv];
    }
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void MasaVect(int nelem, int nnode, int ngdln,int nevab, int ntotv,int* iffix, int **lnods, double* deslo, double* tenlo, double* aux,
              double* aslod, double*** smmat,int threads,double** deslm,double** tenlm,double** vecth)
{
    int ielem,inode,nodei,igdln,nrows,nrowe,ievab,jevab,itotv,tid,jt;
    
    #pragma omp parallel for default(shared)
    for(itotv=1; itotv<=ntotv; itotv++){
        aslod[itotv] = 0.0;
        for (jt=1; jt<=threads;jt ++)
            vecth[jt][itotv]=0.0;
    }

#pragma omp parallel for default(shared)private(ielem,inode,nodei,igdln,nrows,nrowe)
    for(ielem=1; ielem<=nelem; ielem++){
        for(inode=1; inode<=nnode; inode++) {
            nodei = lnods[ielem][inode];
  //          printf("nodei = %d \n",nodei);
            for(igdln=1; igdln<=ngdln; igdln++) {
                nrows = ((nodei-1)*ngdln)+igdln;
                nrowe = ((inode-1)*ngdln)+igdln;
                if(iffix[(nodei-1)*ngdln+igdln]==0)
                    deslm[ielem][nrowe] = aux[nrows];
                else deslm[ielem][nrowe] = 0.0;
            }
        }
    }
 //   printf("vamos aqui 1\n");
#pragma omp parallel for default(shared)private(ielem,ievab,jevab)
    for(ielem=1; ielem<=nelem; ielem++){
        for(ievab=1; ievab<=nevab; ievab++) {
            tenlm[ielem][ievab] = 0.0;
            for(jevab=1; jevab<=nevab; jevab++)
                tenlm[ielem][ievab] =tenlm[ielem][ievab]+smmat[ielem][ievab][jevab]*deslm[ielem][jevab];
        }
    }
//    printf("vamos aqui \n");
#pragma omp parallel for default(shared)private(ielem,tid,inode,nodei,igdln,nrows,nrowe)
    for(ielem=1; ielem<=nelem; ielem++) {
        //tid = omp_get_thread_num()+1;
        tid=1;
   //     printf("tid=%d\n",tid);
        for(inode=1; inode<=nnode; inode++) {
            nodei = lnods[ielem][inode];
            for(igdln=1; igdln<=ngdln; igdln++) {
                nrows = ((nodei-1)*ngdln)+igdln;
                nrowe = ((inode-1)*ngdln)+igdln;
                vecth[tid][nrows]+= tenlm[ielem][nrowe];
            }
        }
    }
#pragma omp parallel for default(shared)
	for(itotv=1; itotv<=ntotv; itotv++){
        if(iffix[itotv] !=0) aslod[itotv] = 0.0;
        else for (jt=1; jt<=threads;jt ++)
            aslod[itotv]+=vecth[jt][itotv];
    }
    
}
              
