//*****************************
// Rutina Para Optimizar
//*****************************
#pragma hdrstop
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
//#include <vector>
#include <algorithm> 
#include "datos.h"
#include "estatico.h"
#include "ValoresPropios.h"
#include "fuerzas.h"
#include "main.h"
#include "posprogid.h"
#include "principal.h"
#include "raros.h"
#include "solver.h"
#include "memoria.h"
#include "Optimizador.h"

std::vector<Individuo> poblacion;
std::vector<Individuo> selec;
std::vector<Individuo> hijos;
std::vector<Individuo> newpob;
std::vector<Individuo> NP;
std::vector<Individuo> Elegibles;
std::vector<Individuo> Penalizados;
Individuo vecino;
int dimension;
int tam_pob=50;
std::vector<int> max_secciones;
int evaluaciones;

int secciones_max;
std::vector<std::vector<int> > Grafo;
std::vector<std::vector<int> > Grafobarras;
std::vector<std::vector<int> > MatrixAdy;
std::vector<bool> visitados;
std::vector<int> bars;
std::vector<int> barras_tent;

int cont_bfs=0;
int num1;

//************************************************

int Optimizador( int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl,int prueba )
{
/**********************
Optimizador
**********************/
srand(time(NULL)); //semilla 

Grafo.resize(npnod);
Grafobarras.resize(nelem);


MatrixAdy.assign(npnod,std::vector<int> (npnod,-1));

for(int ielem=1; ielem<=nelem; ielem++) {//elementos=cantidad de barras
    
    MatrixAdy[lnods[ielem][1]-1][lnods[ielem][2]-1]=ielem-1;//Llena matriz de adyacencias
    MatrixAdy[lnods[ielem][2]-1][lnods[ielem][1]-1]=ielem-1;
    Grafo[lnods[ielem][1]-1].push_back(lnods[ielem][2]-1);
    Grafo[lnods[ielem][2]-1].push_back(lnods[ielem][1]-1);
   

}


//Llenado grafo de barras

for(int ielem=1; ielem<=nelem; ielem++) {//elementos=cantidad de barras
    

    for (int i = 0; i < MatrixAdy[lnods[ielem][1]-1].size(); ++i)
    {
       
        if (MatrixAdy[lnods[ielem][1]-1][i]!=-1 && MatrixAdy[lnods[ielem][1]-1][i]!=ielem-1)
        {
            Grafobarras[ielem-1].push_back( MatrixAdy[lnods[ielem][1]-1][i]);
        }
    }

    for (int i = 0; i < MatrixAdy[lnods[ielem][2]-1].size(); ++i)
    {
        
        if (MatrixAdy[lnods[ielem][2]-1][i]!=-1 && MatrixAdy[lnods[ielem][2]-1][i]!=ielem-1)
        {
            Grafobarras[ielem-1].push_back( MatrixAdy[lnods[ielem][2]-1][i]);
        }

    }

 }


char buf[150];
sprintf(buf, "%d_%d_reemp_bc_1000.txt",nelem, prueba);
//sprintf(buf, "/home/maria.moreno/meca_Paso10/%d_reemp.txt", prueba);


 FILE * pFile;
pFile = fopen (buf,"w");

   

evaluaciones=0;
bool Optimiza = true;
bool generaMatRigidez;
int ival, casoCarga, jj, Num_iter, imats;
double beta, alfa, penaliza, dd, de, deltaf, fact;
//srand(time(NULL)); //semilla 
// Vector de fuerza de cada caso de carga

//cada caso de carga?
    for( casoCarga = 1; casoCarga <= ncarg; casoCarga++ ){
        npesp[casoCarga]=0;
    }

    
    //ncarg=1;
    
for( casoCarga = 1; casoCarga <= ncarg; casoCarga++ ){

// Lee caso de carga y al mismo tiempo ensambla el vector de fuerzas
    
    Fuerzas(nelem,npnod,ndime,ntipo,nevab,nnode,ngdln,iwrit,casoCarga,lnods,matnu,indfu,ntips,
            coord,xlong,props,fuerb,wcarg,carpp,fuepp,vectr,carga,fuemp,Optimiza,
            npesp[casoCarga]);

// Copia el vector de fuerzas del caso de carga, en caso de que se vaya a optimizar
    copia_vector_fuerza( casoCarga, nelem, nevab, carga, vectores_fuerzas );
    //mg vectores fuerzas es una matriz

}
    
    if( iwrit ==1 ) fprintf(fpLog,"Calculo de Fuerzas = %d \n\n",Prob());

// Inicializa valores de vector de materiales

    
Inicializa_Materiales( nmats, ndime, props, matva, cArmaduras.tamCat, cAceroRF.tamCat,// ----------> mg funcion que esta en principal.cpp
                       cAceroRC.tamCat, cConcr.tamCat, ubiCat );

	secciones_max=matva[1];//maxima cantidad de secciones

// Valor inicial de la funcion de optimizacion
//valor[4] = 1.e10;

// Primer iteracion
Actualiza(nmats,matva,matvm,valor,valom);

CargaPropCatOpt(matva,ndime,nmats,props,cArmaduras.arear,cAceroRF.arerf,cAceroRC.arerc,
                cConcr.arecr,ubiCat);
// Analiza los casos de carga
valor[2] = 0.0; // Para guardar el desplazamiento maximo
valor[3] = 0.0; // Para guardar la eficiencia maxima
generaMatRigidez = true;
    
for( casoCarga = 1; casoCarga <= ncarg; casoCarga++ ){
// Escoge el vector de fuerzas del caso de carga
    escoge_vector_fuerza(casoCarga,nelem,nevab,carga,vectores_fuerzas);
// Incluye las fuerzas de peso propio al vector de fuerzas
    incluye_peso_propio(npesp[casoCarga],nelem,ndime,ntipo,lnods,matnu,ntips,coord,
                        xlong,props,vectr,carga);
// Determina el vector de fuerzas de empotramiento perfecto
    vector_fuer_emp_perf(nelem,nevab,fuemp,carga);
    // Incluye las fuerzas de peso propio al vector de fuerzas
    incluye_peso_propio(npesp[casoCarga],nelem,ndime,ntipo,lnods,matnu,ntips,coord,
                        xlong,props,vectr,carga);
   // printf("indso = %d \n",indso);//mg indice del solucionador que se utilizara
    ival = Estatico( nelem,npnod,nevab,ndime,ntipo,ncaso,npres,false,indso,isale,ngdln,nnode,
                     ntotv,nincl,nreso,neqns,nwktl,ishot,lnods,matnu,inpre,iffix,ntips,maxad,
                     nodea,nodpr,leqns,lincl,lreso,react,coord,presc,fixed,despl,aslod,stiff,
                     srmat,props,xlong,vectr,rigid,carpp,fuepp,carga,fuemp,girom,girtm,tempr,
                     fuerc,valor,cArmaduras.arear,cAceroRF.arerf,cAceroRC.arerc,cConcr.arecr,
                     cConcr.dtcon,indfu,fuerb,fuefl,wcarg,xincl,resor,astif,vectp,vectu,vect1,
                     deslo,ubiCat,efi_max,generaMatRigidez,Optimiza);
// Ya se genero una vez la matriz de rigidez
generaMatRigidez = false;
}
//    printf("cero \n");
  
// Parametros para el recosido simulado

    beta = 1.0;
    alfa = 1.001;
    //penaliza = 1.e7;
   // Num_iter = 1000;
    dimension=nmats;


    //@mg 
    max_secciones.resize(dimension);

    Genera_Poblacion_Inicial(ishot, npnod,nelem, npres,  ncarg, ntipo, ngdln, nmats,
                 ndime, iwrit, isale, nreso, nincl, nfami, ncaso,
                  nnode, mnode, nprop, ngaus, ntens, nevab, indso, ntotv,
                 lnods, inpre, nareas, narerf, narerc, narecon,
                 lreso, lincl, ntips, nodpr, matnu, matva, matvm,
                 iffix, indfu, nodea, maxad,mhigh, leqns,
                 coord, carga,carpp, presc, props,
                 fuemp, fuepp, rigid, xmasa, vmatr,
                 girom, girtm, tempr, fuefl, fuerb,
                 srmat, smmat, cArmaduras, cAceroRF,
                 cAceroRC,cConcr, resor, xincl,
                 fuerc, aslod, despl, fixed, react,
                 vecta, vect1, deslo, xlong, valor,
                 valom, angulo, angulg, wcarg, vectr,
                 astif, vectp, vectu, listaCatalogos, ubiCat,
                 vectores_fuerzas, npesp, efi_max, stiff, neqns,
                 nwktl);

    


    int inicio=0,gen=0;
    sort(poblacion.begin(),poblacion.end());


        fprintf (pFile,"\nPoblacion Inicial\n");
        fflush(pFile);
       fprintf (pFile,"\n");
       fflush(pFile);


            for (int i = 0; i < poblacion.size(); ++i)
            {
                for (int j = 0; j < poblacion[i].secciones.size(); ++j)
                {
                   fprintf (pFile,"%d ", poblacion[i].secciones[j]);
                   fflush(pFile);
                }
               // std::cout<<poblacion[i].fit<<std::endl;
                fprintf (pFile," fit: %lf desplazamiento: %lf eficiencia: %lf \n", poblacion[i].fit, poblacion[i].desp, poblacion[i].efi);
                fflush(pFile);
               /* for (int k = 0; k < poblacion[i].valores.size(); ++k)
                {
                    printf(" %lf ",poblacion[i].valores[k]);
                    printf("\n");
                }*/
                
            }

        fprintf (pFile,"\n");
        fflush(pFile);
        fprintf (pFile,"\n");
        fflush(pFile);

       // exit(0);
        double Dini,D,Dgeneracional;

do{
    

    parent_selection();
    Crossover();
    mutacion(ishot, npnod,nelem, npres,  ncarg, ntipo, ngdln, nmats,
                 ndime, iwrit, isale, nreso, nincl, nfami, ncaso,
                  nnode, mnode, nprop, ngaus, ntens, nevab, indso, ntotv,
                 lnods, inpre, nareas, narerf, narerc, narecon,
                 lreso, lincl, ntips, nodpr, matnu, matva, matvm,
                 iffix, indfu, nodea, maxad,mhigh, leqns,
                 coord, carga,carpp, presc, props,
                 fuemp, fuepp, rigid, xmasa, vmatr,
                 girom, girtm, tempr, fuefl, fuerb,
                 srmat, smmat, cArmaduras, cAceroRF,
                 cAceroRC,cConcr, resor, xincl,
                 fuerc, aslod, despl, fixed, react,
                 vecta, vect1, deslo, xlong, valor,
                 valom, angulo, angulg, wcarg, vectr,
                 astif, vectp, vectu, listaCatalogos, ubiCat,
                 vectores_fuerzas, npesp, efi_max, stiff, neqns,
                 nwktl);
    new_pob();   
    	//lena elegibles con papas
    	for (int i = 0; i < poblacion.size(); ++i)
         {
         	Elegibles.push_back(poblacion[i]);//se guardan padres en elegibles
         } 

    actualizacion(matvm,valom);//guarda al mejor


    local_search(ishot, npnod,nelem, npres,  ncarg, ntipo, ngdln, nmats,
                 ndime, iwrit, isale, nreso, nincl, nfami, ncaso,
                  nnode, mnode, nprop, ngaus, ntens, nevab, indso, ntotv,
                 lnods, inpre, nareas, narerf, narerc, narecon,
                 lreso, lincl, ntips, nodpr, matnu, matva, matvm,
                 iffix, indfu, nodea, maxad,mhigh, leqns,
                 coord, carga,carpp, presc, props,
                 fuemp, fuepp, rigid, xmasa, vmatr,
                 girom, girtm, tempr, fuefl, fuerb,
                 srmat, smmat, cArmaduras, cAceroRF,
                 cAceroRC,cConcr, resor, xincl,
                 fuerc, aslod, despl, fixed, react,
                 vecta, vect1, deslo, xlong, valor,
                 valom, angulo, angulg, wcarg, vectr,
                 astif, vectp, vectu, listaCatalogos, ubiCat,
                 vectores_fuerzas, npesp, efi_max, stiff, neqns,
                 nwktl);



        for (int i = 0; i < poblacion.size(); ++i)//se guardan hijos en elegibles
         {
         	Elegibles.push_back(poblacion[i]);
         } 
         if (gen==0)
         {
         	Dini=distancia_inicial(Elegibles);
         	
         }
         
         D=Dini*(1-((double)gen/900.0));
        

         reemplazamiento_npb(D);//ddentro se ordena a Elegibles y se hace el reemplazamiento
         
         poblacion.clear();

         for (int i = 0; i < selec.size(); ++i)
         {
         	poblacion.push_back(selec[i]);
         }


         Elegibles.clear();
         Penalizados.clear();
         selec.clear();

    	sort(poblacion.begin(),poblacion.end());//ordena poblacion después de haber pasado por reemplazamiento


    if (gen%10==0)
    {
        Dgeneracional=distancia_inicial(poblacion);
        printf("%d %lf %lf\n",gen,Dgeneracional,poblacion[0].fit);
        fprintf (pFile,"%d %lf %lf\n",gen,Dgeneracional,poblacion[0].fit);
       // fprintf (pFile,"\n");
        fflush(pFile);

    }

    gen++;
    
 }while(gen <=1000);



        /*fprintf (pFile,"Población final: \n");
        fprintf (pFile,"\n");
        fprintf (pFile,"\n");


            for (int i = 0; i < poblacion.size(); ++i)
            {
                for (int j = 0; j < poblacion[i].secciones.size(); ++j)
                {
                   fprintf (pFile,"%d ", poblacion[i].secciones[j]);
                }
                fprintf (pFile," fit: %lf desplazamiento: %lf eficiencia: %lf \n", poblacion[i].fit, poblacion[i].desp, poblacion[i].efi);
                
            }

        fprintf (pFile,"\n");
        fprintf (pFile,"\n");
        fprintf (pFile,"\n");


        fprintf (pFile,"Mejor individuo encontrado: \n");
 for (int j = 0; j < poblacion[0].secciones.size(); ++j)
                {
                   fprintf (pFile,"%d ", poblacion[0].secciones[j]);
                }
                
fprintf (pFile,"\n");
fprintf (pFile,"%lf ", poblacion[0].fit);*/
fclose (pFile);

return ival;

}
//***********************************************************************************************************

void reemplazamiento_npb(double D){

	int dist;

	do{

			sort(Elegibles.begin(),Elegibles.end());
			//Elegibles[0].d=-155000;
			selec.push_back(Elegibles[0]);
			Elegibles.erase(Elegibles.begin()+0);//lo quito poruqe es el de mejor fitness


			for (int i = 0; i < Elegibles.size(); ++i)
			{
				//distancia entre el nuevo de la poblacion y los elegibles para actualizar dicha distancia
				Elegibles[i].d=distancia_individuos(Elegibles[i],selec[selec.size()-1]);
				
			}

			for (int i = 0; i < Penalizados.size(); ++i)
			{
				//distancia entre el nuevo de la poblacion y los penalizados para actualizar dicha distancia
				Penalizados[i].d=distancia_individuos(Penalizados[i],selec[selec.size()-1]);
				
			}
			

			for (int i = 0; i < Elegibles.size(); ++i)
			{
				if (Elegibles[i].d < D)
				{
					Penalizados.push_back(Elegibles[i]);
					Elegibles.erase(Elegibles.begin()+i);
					i--;//para que vuelva a evaluar la posicion que quito ya que ahi debe estar la posicion siguiente a evaluar
				}		
				
			}

	}while(((int)Elegibles.size() != 0) && (selec.size() < poblacion.size()));//hasta que se acaben los elegibles

	int dmayor,indmayor;
	while(selec.size() < poblacion.size())//mientras no se haya llenado nuevamente toda la poblacion empezará a utilizar penalizados
	{
		dmayor=0;
		for (int i = 0; i < Penalizados.size(); ++i)
		{
			if (dmayor < Penalizados[i].d)
			{
				dmayor=Penalizados[i].d;
				indmayor=i;
			}
			
		}

		selec.push_back(Penalizados[indmayor]);
		Penalizados.erase(Penalizados.begin()+indmayor);


	}


}


double distancia_inicial(std::vector<Individuo> Pob){

	int cont=0,acum=0;

	for(int i=0; i<Pob.size(); i++){//ind
		for(int j=0; j<Pob.size(); j++){//itera con los demas ind
			if (i==j)continue;
			acum+=distancia_individuos(Pob[i],Pob[j]);
			//printf("%d\n", acum);
			cont++;
			//printf("%d\n", cont);
		}

	}

	return ((double)acum/(double)cont)*0.5;

}

double distancia_individuos(Individuo a, Individuo b){


	int sum=0;
	for (int i = 0; i < a.secciones.size(); ++i)
	{
		sum+=abs(a.secciones[i]-b.secciones[i]);
	}

	return sum;

}

void local_search(int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl){


    double num,n;
    int secc,cont,j,jj,b=0,modif,mejor=0;
    //std::set<int> myset1;

    for (int i = 0; i < tam_pob; ++i)//recorre  poblacion
    { 
    		b=0;
    		//CREA AL VECINO-ANTES DE MODIFICAR
    		for (int k= 0; k < poblacion[i].secciones.size(); ++k)//recorre secciones 
	        {

	            vecino.secciones.push_back(poblacion[i].secciones[k]);
	        }


	       for (int k= 0; k < poblacion[i].valores.size(); ++k)//recorre secciones 
	        {

	            vecino.valores.push_back(poblacion[i].valores[k]);
	        }


    	do{  

								//std::cout<<"vecino size "<<vecino.secciones.size()<<std::endl;
					    	modif=rand()%(poblacion[i].secciones.size());//cual posicion voy a modificar


					    		//std::cout<<"modif "<<modif<<" secciones size "<<poblacion[i].secciones.size()<<" barras size "<<poblacion[i].barras.size()<<std::endl;
						        

						     //if (modif <= vecino.seccio

							for (int k= 0; k < poblacion[i].secciones.size(); ++k)//regreso al individuo original
					        {


					           matva[k+1]=poblacion[i].secciones[k];
                               vecino.secciones[k]=poblacion[i].secciones[k];
					                  // 	std::cout<<" "<<matva[k+1]<<" ";
					        }

					 //std::cout<<" matnu "<<std::endl;
					       /* for (int k= 0; k < poblacion[i].barras.size(); ++k)//recorre secciones 
					        {
					        	matnu[k+1]=vecino.barras[k];
					        	        //	std::cout<<" "<<matnu[k+1]<<" ";
					        }*/
                            do{
                                secc=1+rand()%(secciones_max+1-1);

                            }while(secc==vecino.secciones[modif]);


					        vecino.secciones[modif]=secc;
					        matva[modif+1]=secc;
					          	
					        
					            
					                
					       
						//printf("secciones %d  %d \n",poblacion[i].secciones[0],poblacion[i].secciones[1]);


							vecino.fit=Funcion_fitness(ishot, npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,
					                                                  nfami,ncaso, nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,lnods,inpre,
					                                                  nareas,narerf,narerc,narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,
					                                                  iffix,indfu,nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
					                                                  fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,fuerb,srmat,smmat,
					                                                  cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,fixed,
					                                                  react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
					                                                  vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,stiff,
					                                                  neqns, nwktl);//@mg
					        
					        /*poblacion[i].fit=Funcion_fitness(ishot, npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,
					                                                  nfami,ncaso, nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,lnods,inpre,
					                                                  nareas,narerf,narerc,narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,
					                                                  iffix,indfu,nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
					                                                  fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,fuerb,srmat,smmat,
					                                                  cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,fixed,
					                                                  react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
					                                                  vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,stiff,
					                                                  neqns, nwktl);//@mg*/

					       // std::cout<<" valores tamñano "<<vecino.valores.size()<<std::endl;
                            vecino.desp=valor[2];
                            vecino.efi=valor[3];

					        vecino.valores.clear();

					       for (int j = 1; j < 5; ++j)//4 porque son 4 valores, peso propio,desplazamiento maximo,eficiencia maxima
					        {
					           vecino.valores.push_back(valor[j]);
					           // printf("valores %lf ",valor[j] );
					        }

					        //printf("\n");
					     


					        if (vecino.fit < poblacion[i].fit)
					        {
					        	//printf("ENCONTRO MEJOR \n");
                                //mejor=i;
					        	actualizar_vecino_poblacion(vecino,i);
					        	//break;
					        }


					        b++;

					       // std::cout<<"buscando becino para "<<i<<std::endl;
					        // std::cout<<std::endl;

    }while(b < 75/*((59*poblacion[0].secciones.size())+(poblacion[0].secciones.size()*poblacion[0].barras.size()))*/);

                              /* printf("Individuo despues de la busqueda local \n");
    for (int i = 0; i < poblacion[0].secciones.size(); ++i)
    {
        printf("%d  ", poblacion[0].secciones[i]);
    }
    printf("%lf  \n", poblacion[0].fit);

    exit(0);*/

                           // exit(0);
        vecino.secciones.clear();
       // vecino.barras.clear();
        vecino.valores.clear();
    
        
        
    }



}

void actualizar_vecino_poblacion(Individuo vecino,int i){

	poblacion.erase(poblacion.begin()+i);

	poblacion.insert(poblacion.begin()+i,vecino);



}


void parent_selection_bnp(double D){

	selec.push_back(poblacion[0]);//pasa el best de primero


	for (int i = 1; i < poblacion.size(); ++i)
	{
		Elegibles.push_back(poblacion[i]);
		Elegibles[i-1].d=99999;
	}

	int dist,mayor_dist,mayor_valor;

	


	do{

				mayor_valor=0;
				for (int i = 0; i < Elegibles.size(); ++i)
				{


					dist=distancia_individuos(selec[selec.size()-1],Elegibles[i]);//distancia entre el ultimo de los individuos selecccionados y todos los elegibles
				
					if (dist < Elegibles[i].d)
					{
						Elegibles[i].d=dist;
					}

					if (Elegibles[i].d < D)
					{
						Penalizados.push_back(Elegibles[i]);
						//borrar.push_back(i);
						
						Elegibles.erase(Elegibles.begin()+i);//pasa a individuo a penalizados 
						i=i-1;

					}

					if (Elegibles[i].d > mayor_valor)
					{
						mayor_valor=Elegibles[i].d;
						mayor_dist=i;//indice con el q tiene mayor distancia
					}
				}


				selec.push_back(Elegibles[mayor_dist]);//se agregrego a selec un nuevo individuo
				Elegibles.erase(Elegibles.begin()+mayor_dist);//pasa a individuo a penalizados 


	}while(Elegibles.size()!=0);


	
	
	while(selec.size()<poblacion.size())
	{
		mayor_valor=0;
		for (int i = 0; i < Penalizados.size(); ++i)
		{	
			if (Penalizados[i].d > mayor_valor)
			{
				mayor_dist=i;
				mayor_valor= Penalizados[i].d;
			}
		}

		selec.push_back(Penalizados[mayor_dist]);
		Penalizados.erase(Penalizados.begin()+mayor_dist);

	}


	Elegibles.clear();
	Penalizados.clear();

}





void parent_selection(){

	int a,b;
	std::vector<int> pseudopob;

	for (int i = 0; i < tam_pob; ++i)
	{
		pseudopob.push_back(i);
	}

	random_shuffle(pseudopob.begin(),pseudopob.end());//vector de posiciones de los individuos randomizados!!
	Individuo Indi;

	for (int i = 0; i < tam_pob; ++i)
	{
		selec.push_back(Indi);
            

        if(poblacion[i].fit > poblacion[pseudopob[i]].fit) {

			for (int j = 0; j < dimension; ++j)
			{
				selec[i].secciones.push_back(poblacion[pseudopob[i]].secciones[j]);

			}
			
		}else{

			for (int j = 0; j < dimension; ++j)
			{
				selec[i].secciones.push_back(poblacion[i].secciones[j]);

			}
		}

		//fitpob.push_back((this->*f)(seleccion[i]));


	}

}

void BC_barras(int v){//para recorrer grafo GRAFO
    visitados.assign(Grafobarras.size(),false);
    //printf("%d\n",v);
   
    bars.push_back(v);
    cont_bfs++;
    for (int i = 0; i < Grafobarras[v].size(); ++i)
    {
        //printf("tentativos\n");
        barras_tent.push_back(Grafobarras[v][i]);
    }
    
    visitados[v]=true;
    double numrand;

    while(cont_bfs<num1){

        numrand=rand()%(barras_tent.size());
        if (!visitados[barras_tent[numrand]]){

            bars.push_back(barras_tent[numrand]);
            visitados[barras_tent[numrand]]=true;
            //printf("%d\n",barras_tent[numrand]);

            for (int i = 0; i < Grafobarras[barras_tent[numrand]].size(); ++i)
            {
                barras_tent.push_back(Grafobarras[barras_tent[numrand]][i]);
            }
            cont_bfs++;

        }




    }

}


void Crossover(){
    double point,num;
    Individuo Indi;

        int aux,unif,n;
        std::vector<int> hijo1;
        std::vector<int> hijo2;


        for (int k = 0; k < tam_pob;k=k+2)//OJO
        {

            num=rand()/(double)RAND_MAX;

            if(num<0.8){//PROBABILIDAD DE CROSSOVER

                hijo1=selec[k].secciones;
                hijo2=selec[k+1].secciones;

                n=rand()%(Grafobarras.size());//en que barra empieza
                num1=rand()%(Grafobarras.size());//maximo numero de barras 
                bars.clear();//barras conexas
                cont_bfs=0;//cuenta los nodos que se van abriendo
                visitados.clear();//reinicia a los visitados.
                barras_tent.clear();

                BC_barras(n);

                

                for (int i = 0; i < bars.size(); ++i)
                {
                    hijo1[bars[i]]=selec[k+1].secciones[bars[i]];
                    hijo2[bars[i]]=selec[k].secciones[bars[i]];
                }



                hijos.push_back(Indi);
                hijos.push_back(Indi);
                for (int l = 0; l < dimension; ++l)
                {
                    hijos[k].secciones.push_back(hijo1[l]);
                    hijos[k+1].secciones.push_back(hijo2[l]);
                }
                //hijos.push_back(hijo1);
                //hijos.push_back(hijo2);
                hijo1.clear();
                hijo2.clear();
            


            }else{


                hijos.push_back(Indi);
                hijos.push_back(Indi);
                for (int l = 0; l < dimension; ++l)
                {
                    hijos[k].secciones.push_back(selec[k].secciones[l]);
                    hijos[k+1].secciones.push_back(selec[k+1].secciones[l]);
                }

            }

        }   


}

void mutacion(int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl){

    double num,n,secc;

    for (int i = 0; i < tam_pob; ++i)//recorre  
    {   
        
        
        for (int j = 0; j < dimension; ++j)//recorre secciones segun sea viga o columna
        {
             num=rand()/(double)RAND_MAX;
            //cout<<num<<endl;
            if (num<0.055)//antes en 0.055-(1/(double)nelem)
            {
                secc=1+rand()%(max_secciones[j]+1-1);//necesito saber  cual es el numero maximo de secciones en el catalogo 1-59
                hijos[i].secciones[j]=(int)secc;
                    
            }


            matva[j+1]=hijos[i].secciones[j];
                
        }
	//printf("secciones %d  %d \n",hijos[i].secciones[0],hijos[i].secciones[1]);


	hijos[i].fit=Funcion_fitness(ishot, npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,//QUEDE AQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                                                  nfami,ncaso, nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,lnods,inpre,
                                                  nareas,narerf,narerc,narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,
                                                  iffix,indfu,nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
                                                  fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,fuerb,srmat,smmat,
                                                  cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,fixed,
                                                  react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
                                                  vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,stiff,
                                                  neqns, nwktl);//@mg
        
        /*hijos[i].fit=Funcion_fitness(ishot, npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,
                                                  nfami,ncaso, nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,lnods,inpre,
                                                  nareas,narerf,narerc,narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,
                                                  iffix,indfu,nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
                                                  fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,fuerb,srmat,smmat,
                                                  cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,fixed,
                                                  react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
                                                  vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,stiff,
                                                  neqns, nwktl);//@mg*/
        hijos[i].desp=valor[2];
        hijos[i].efi=valor[3];


       for (int j = 1; j < 5; ++j)//4 porque son 4 valores, peso propio,desplazamiento maximo,eficiencia maxima
        {
            hijos[i].valores.push_back(valor[j]);
           // printf("valores %lf ",valor[j] );
        }

        //printf("\n");
     

        
        
    }




}


int new_pob(){


    sort(poblacion.begin(),poblacion.end());
    sort(hijos.begin(),hijos.end());


    int cont2=0,bandera=0;

    for (int i = 0; i < hijos.size(); ++i)
    {   
        if(i==0 && (poblacion[i].fit < hijos[i].fit)) //solo si el mejor de los padres es mejor de los hijos este sera conservado para la siguiente generacion 
        {
            newpob.push_back(poblacion[0]);

        }else{

            newpob.push_back(hijos[cont2]);
            cont2++;

        }
        
    }

    

    
    return bandera;


}

void actualizacion(int* matvm,double* valom){

    for (int i = 0; i < tam_pob; ++i)
    {
        for (int j = 0; j< dimension; ++j)
        {
            poblacion[i].secciones[j]=newpob[i].secciones[j];
        }
        for (int j = 0; j< 4; ++j)
        {
            poblacion[i].valores[j]=newpob[i].valores[j];
        }

        poblacion[i].fit=newpob[i].fit;
        poblacion[i].desp=newpob[i].desp;
        poblacion[i].efi=newpob[i].efi;
    }

    newpob.clear();
    hijos.clear();
    selec.clear();

    for (int j = 0; j< dimension; ++j)
    {
        matvm[j+1]=poblacion[0].secciones[j];//copiar el mejor
    }

    for (int j = 0; j< 4; ++j)
    {
       valom[j+1]=poblacion[0].valores[j];
    }



    
}



void Genera_Poblacion_Inicial(int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                 int nwktl){//@mg 

	//int tam_pob=10;
  	poblacion.resize(tam_pob);  

  	
	//poblacion.resize(1);
//
	for (int i = 0; i < tam_pob; ++i)
	//for (int i = 0; i < 1; ++i)
	{
       // printf("\nINDIVIDUO %d\n", i);
		Genera_Nuevo_Individuo(nmats, ndime, matva, props, cArmaduras.tamCat,
                           cAceroRF.tamCat, cAceroRC.tamCat, cConcr.tamCat, ubiCat);

		
		//matva[1]=8;
  		//matva[2]=2;


		poblacion[i].fit=Funcion_fitness(ishot, npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,//QUEDE AQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                                                  nfami,ncaso, nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,lnods,inpre,
                                                  nareas,narerf,narerc,narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,
                                                  iffix,indfu,nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
                                                  fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,fuerb,srmat,smmat,
                                                  cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,fixed,
                                                  react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
                                                  vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,stiff,
                                                  neqns, nwktl);//@mg'

        poblacion[i].desp=valor[2];
        poblacion[i].efi=valor[3];

		//printf(" ALOOOOOO %lf \n",poblacion[i].fit );

    	
		for (int j = 0; j < nmats; ++j)
		{
			poblacion[i].secciones.push_back(matva[j+1]);//haciendo individuo

    		//poblacion[i].secciones.push_back(valor[j+3]);
		}
       // printf("pob secciones %d  %d\n", poblacion[i].secciones[0],poblacion[i].secciones[1]);
		for (int j = 1; j < 5; ++j)//4 porque son 4 valores, peso propio,desplazamiento maximo,eficiencia maxima
		{
			poblacion[i].valores.push_back(valor[j]);//
    		//poblacion[i].secciones.push_back(valor[j+3]);
		}
	}

}



double Funcion_fitness( int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl){//@mg 
	
	double dd,de,fitness;
	double penaliza = 1.e7;
   // bool Optimiza = true;
    bool generaMatRigidez,Optimiza=true;
    int ival, casoCarga;

	CargaPropCatOpt(matva,ndime,nmats,props,cArmaduras.arear,cAceroRF.arerf,cAceroRC.arerc,cConcr.arecr,ubiCat);
		valor[2] = 0.0; // Para guardar el desplazamiento maximo
	    valor[3] = 0.0; // Para guardar la eficiencia maxima
	    generaMatRigidez = true;

	    for( casoCarga = 1; casoCarga <= ncarg; casoCarga++ ){
			// Escoge el vector de fuerzas del caso de carga
	      escoge_vector_fuerza(casoCarga,nelem,nevab,carga,vectores_fuerzas);
			// Incluye las fuerzas de peso propio al vector de fuerzas
	      incluye_peso_propio(npesp[casoCarga],nelem,ndime,ntipo,lnods,matnu,ntips,coord,
	                         xlong,props,vectr,carga);
			// Determina el vector de fueras de empotramiento perfecto
			//      vector_fuer_emp_perf(nelem,nevab,fuemp,carga);
  	      ival = Estatico( nelem,npnod,nevab,ndime,ntipo,ncaso,npres,false,indso,isale,ngdln,nnode,
                     ntotv,nincl,nreso,neqns,nwktl,ishot,lnods,matnu,inpre,iffix,ntips,maxad,
                     nodea,nodpr,leqns,lincl,lreso,react,coord,presc,fixed,despl,aslod,stiff,
                     srmat,props,xlong,vectr,rigid,carpp,fuepp,carga,fuemp,girom,girtm,tempr,
                     fuerc,valor,cArmaduras.arear,cAceroRF.arerf,cAceroRC.arerc,cConcr.arecr,
                     cConcr.dtcon,indfu,fuerb,fuefl,wcarg,xincl,resor,astif,vectp,vectu,vect1,
                     deslo,ubiCat,efi_max,generaMatRigidez,Optimiza);
			// Ya se genero una vez la matriz de rigidez
	      generaMatRigidez = false;
    	}


	dd = valor[2]-10.0;  // Castiga desplazamientos mayores a 10cm. Aparentemente aqui podemos poner claro/norma.
	de = valor[3]-1.0;   // Castiga eficiencia maxima. Ahorita la trata de llevar al 100%.
	valor[4]=valor[1];   // Peso de la estructura total.
         	if (dd >0) valor[4]+=dd*penaliza;// mg se agrega al valor de la funcion
            if (de >0) valor[4]+=de*penaliza;// mg se agrega al valor de la funcion 
    fitness= valor[4];
    evaluaciones++;
    
   /* printf(" En Funcion_fitness \n ",fitness);

    for(int i=1;i<=nelem;i++){
        printf(" %d ",matva[i]);
    }

     printf(" %lf \n ",fitness);*/

	return fitness;


}


void Genera_Nuevo_Individuo( int nmats, int ndime, int* matva, double** props,//@mg 
                             int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
                             int *Con_tamCat, int *ubiCat){


	int rolado, iCat;
	// Escoge un material al azar de la lista de materiales. Entre 1 y nmats.
	//mg otra opcion seria variable = limite_inferior + rand() % (limite_superior +1 - limite_inferior) ;
	//imats=(int)(((float) rand() / (float) RAND_MAX)*(float)nmats-.000000000001)+1;//lo que elige aqui es la columna o viga 
	/*printf("nMATS: %d\n",nmats ); //mg lo comento porque lo que se estaba generando era un vecino, ahora quiero es una poblacion
    
    printf("tam max_secciones %d\n", max_secciones.size() );*/
	//printf("ndime: %d\n",ndime );
	for (int imats = 1; imats <= nmats; ++imats)
	{
       // printf("iMATS: %d\n",imats );
        //printf("1 %lf\n",props[imats][5] );
        //printf("2 %lf\n",props[imats][8] );
		//if (ndime ==2 ) rolado= (int)(props[imats][5]+0.5);
	rolado= (int)(props[imats][8]+0.5);
        // printf("rolado: %d\n",rolado );

	// Del material antes escogido, se le cambia al azar la seccion de las
	// disponibles en el catalogo.
		iCat = ubiCat[imats];

		switch( rolado ){
		case 0://tamaño del catalogo para ese material "imats" en especifico.
		matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Arm_tamCat[iCat]-.000000000001)+1;
        	max_secciones[imats-1]=Arm_tamCat[iCat];
		//printf("------------- matva %d \n",matva[imats]);
        	//printf(" Arm %d\n",Arm_tamCat[iCat]);
		break;
		case 1:
		matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RF_tamCat[iCat]-.000000000001)+1;
       		 max_secciones[imats-1]=RF_tamCat[iCat];
		//printf("------------- matva %d \n",matva[imats]);
        	//printf("RF %d\n",RF_tamCat[iCat]);
		break;
		case 2:
		matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RC_tamCat[iCat]-.000000000001)+1;
        	max_secciones[imats-1]=RC_tamCat[iCat];
        	//printf(" RC %d\n",RC_tamCat[iCat]);
		//printf("------------- matva %d \n",matva[imats]);
		break;
		case 3:
		matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Con_tamCat[iCat]-.000000000001)+1;
		//printf("------------- matva %d \n",matva[imats]);
        	max_secciones[imats-1]=Con_tamCat[iCat];
        	//printf("cat_concreto %d\n",Con_tamCat[iCat]);
		break;
		default:
		break;
		}

	}

	
	/*printf("matva1: %d matva2: %d \n",matva[1],matva[2]);
    printf("max_secciones: %d max_secciones: %d \n",max_secciones[0],max_secciones[1]);*/


}



/*void Genera_Nueva_Poblacion( int nmats, int ndime, int* matva, double** props,//mg LO QUE EN REALIDAD GENERA ES UN VECINO, escoje un numero al azar de cual tipo de estructura va a modificar
                             int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
                             int *Con_tamCat, int *ubiCat )
{
int imats, rolado, iCat;
// Escoge un material al azar de la lista de materiales. Entre 1 y nmats.

//mg otra opcion seria variable = limite_inferior + rand() % (limite_superior +1 - limite_inferior) ;
imats=(int)(((float) rand() / (float) RAND_MAX)*(float)nmats-.000000000001)+1;//lo que elige aqui es cambiar la seccion de la columna o de la viga , imats es viga o columna
printf("IMATS: %d\n",imats );

if (ndime ==2 ) rolado= (int)props[imats][5]+0.5;
if (ndime ==3 ) rolado= (int)props[imats][8]+0.5;

// Del material antes escogido, se le cambia al azar la seccion de las
// disponibles en el catalogo.
iCat = ubiCat[imats];
switch( rolado ){
case 0://tamaño del catalogo para ese material "imats" en especifico.
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Arm_tamCat[iCat]-.000000000001)+1;
break;
case 1:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RF_tamCat[iCat]-.000000000001)+1;
break;
case 2:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RC_tamCat[iCat]-.000000000001)+1;
break;
case 3:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Con_tamCat[iCat]-.000000000001)+1;
break;
default:
break;
}

printf("matva1: %d matva2: %d \n",matva[1],matva[2]);

}*/



/*void Genera_Nuevo_Individuo( int nmats, int ndime, int* matva, double** props,
                             int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
                             int *Con_tamCat, int *ubiCat )
{
int imats, rolado, iCat;
// Escoge un material al azar de la lista de materiales. Entre 1 y nmats.

//mg otra opcion seria variable = limite_inferior + rand() % (limite_superior +1 - limite_inferior) ;
imats=(int)(((float) rand() / (float) RAND_MAX)*(float)nmats-.000000000001)+1;//lo que elige aqui 
printf("IMATS: %d\n",imats );

if (ndime ==2 ) rolado= (int)props[imats][5]+0.5;
if (ndime ==3 ) rolado= (int)props[imats][8]+0.5;

// Del material antes escogido, se le cambia al azar la seccion de las
// disponibles en el catalogo.
iCat = ubiCat[imats];//mg tama;o del catalogo?
switch( rolado ){
case 0://tamaño del catalogo para ese material "imats" en especifico.
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Arm_tamCat[iCat]-.000000000001)+1;
break;
case 1:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RF_tamCat[iCat]-.000000000001)+1;
break;
case 2:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RC_tamCat[iCat]-.000000000001)+1;
break;
case 3:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Con_tamCat[iCat]-.000000000001)+1;
break;
default:
break;
}

printf("matva1: %d matva2: %d \n",matva[1],matva[2]);

}*/


//******************************************************************************
void Recosido_Simulado(int Num_iter)
{  /*
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
}


