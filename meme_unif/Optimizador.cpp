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
#include <string>

std::vector<Individuo> poblacion;
std::vector<Individuo> selec;
std::vector<Individuo> hijos;
std::vector<Individuo> newpob;
Individuo vecino;
int dimension;
int tam_pob=50;
std::vector<int> max_secciones;
int evaluaciones;

int secciones_max;
std::vector<std::vector<int> > Grafo;
std::vector<std::vector<int> > MatrixAdy;
std::vector<bool> visitados;
std::vector<int> rama;
std::vector<int> bars;
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


MatrixAdy.assign(npnod,std::vector<int> (npnod,-1));

for(int ielem=1; ielem<=nelem; ielem++) {
    
    MatrixAdy[lnods[ielem][1]-1][lnods[ielem][2]-1]=ielem-1;//Llena matriz de adyacencias
    MatrixAdy[lnods[ielem][2]-1][lnods[ielem][1]-1]=ielem-1;
    Grafo[lnods[ielem][1]-1].push_back(lnods[ielem][2]-1);
    Grafo[lnods[ielem][2]-1].push_back(lnods[ielem][1]-1);

}



char buf[150];
sprintf(buf, "%d_%d_meme_unif.txt", nelem,prueba);
 FILE * pFile;
pFile = fopen (buf,"w");
 // ofs << valom[4]<<"\n";


  //exit(0);

evaluaciones=0;
bool Optimiza = true;
bool generaMatRigidez;
int ival, casoCarga, jj, Num_iter, imats;
double beta, alfa, penaliza, dd, de, deltaf, fact;

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
   // printf("nmats %d \n",nmats );


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

       

    
   gen=1;
   double Dgeneracional;

do{
    gen++;

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
    sort(poblacion.begin(),poblacion.end());

    if (gen%10==0)
    {
        
        Dgeneracional=distancia_inicial(poblacion);
        printf("%d %lf %lf\n",gen,Dgeneracional,poblacion[0].fit);
        fprintf (pFile,"%d %lf %lf\n",gen,Dgeneracional,poblacion[0].fit);
        fflush(pFile);

    }


    
 }while(gen <=1000);


	fclose (pFile);


return ival;

}
//***********************************************************************************************************

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

void DFS_DT(int v, int padre){//para recorrer grafo GRAFO

    int j;
   // rama.push_back(v);
    if (padre != -1)
    {
        visitados[padre]=true;
        bars.push_back(MatrixAdy[v][padre]);
        cont_bfs++;
        
    }
    
        
    for (int j= 0; j < Grafo[v].size(); ++j)
    {
    	if (cont_bfs>num1)
    	{
    		break;
    	}
    	 if (visitados[ Grafo[v][j]] != true)
        {
            DFS_DT( Grafo[v][j],v);
            //break;
            
        }

            	
            

    }

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

void Crossover(){
	double point,num;
	Individuo Indi;

		int aux,unif;
		std::vector<int> hijo1;
		std::vector<int> hijo2;


		for (int k = 0; k < tam_pob;k=k+2)//OJO
		{

			num=rand()/(double)RAND_MAX;

			if(num<0.8){//PROBABILIDAD DE CROSSOVER

				

				for (int i = 0; i < dimension; i++ )//recorre variables
				{
                    unif= rand() / (double)RAND_MAX;//0 + rand()%(dimension-1+1-1);//de 0 a la dimension-1


		            if (unif <= 0.5)
                    {
                        hijo1.push_back(selec[k].secciones[i]);
                        hijo2.push_back(selec[k+1].secciones[i]);

                    }else{

                        hijo1.push_back(selec[k+1].secciones[i]);
                        hijo2.push_back(selec[k].secciones[i]);

                    }

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


