//*****************************
// Rutina Tensiones
//*****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "raros.h"
#include "rigidez.h"
#include "tensiones.h"


//************************************************
void Tensiones(int nelem, int nnode, int ngdln, int ndime, int ntipo, int nevab,
			   int ntotv, int iwrit, int isale, int nincl, int* matnu, int** lnods, int* lincl,
			   double** props, double* despl, double*** srmat, double* fuerc, double** fuepp,
			   double** fuefl, double** vectr, double* xlong, double** girom, double** girtm,
			   double** tempr, double** rigid, double* xincl, double* valor)
/******************************************************************************
 Calcula las tensiones y esfuerzos en los puntos de Gauss.
 ******************************************************************************/
{
	int   ielem, mposn, inode, lnode, nposn, igdln, ievab, jevab,itotv,imats,rolado;
	double tens1, angux, factr,    pi, toler, pesot, despt,pesob, coset, senot, cosef, senof,
	       desel[13], tensg[13];

	pi = 3.14159265;
	toler = 1.0e-5;
 /*   for(itotv=1; itotv<=ntotv; itotv++){
        printf("despl[%d]=%lf \n",itotv,despl[itotv]);
    }
*/
	factr = 180.0/pi;
	// Ciclo sobre cada elemento:
	for(ielem=1; ielem<=nelem; ielem++) {
		imats=matnu[ielem];
		//if(ndime == 2) rolado=(int)props[imats][5]+0.5;
		//if(ndime == 3)
            rolado=(int)props[imats][8]+0.5;

		// Lee la matriz de rigidez elemental:
		if(nincl >0)ApoyoInclinado(ielem,ndime,ngdln,nevab,isale, nincl,lnods,
								   lincl,girom,girtm,tempr,rigid,xincl);

		// Idenifica los desplazamientos y los puntos nodales del elemento:
		mposn = 0;

		for(inode=1; inode<=nnode; inode++) {
			lnode = lnods[ielem][inode];
			nposn = (lnode-1)*ngdln;

			for(igdln=1; igdln<=ngdln; igdln++)
				desel[++mposn] = despl[++nposn];

		}
        mposn=0;
		for(inode=1; inode<=nnode; inode++)
			for(igdln=1; igdln<=ngdln; igdln++)
				// Calculo de las fuerzas en cada nodo:
				for(ievab=1; ievab<=nevab; ievab++) {
					tensg[ievab] = 0.0;

					for(jevab=1; jevab<=nevab; jevab++){
						tensg[ievab] =tensg[ievab]+srmat[ielem][ievab][jevab]*desel[jevab];
					}
					tensg[ievab] =tensg[ievab]+fuepp[ielem][ievab];
				}

		if(ntipo == 1) {
			if(ndime == 2) {
				tens1 = sqrt(tensg[1]*tensg[1]+tensg[2]*tensg[2]);

				if(tens1 != 0.0) {
					angux = tensg[1]/tens1;
		   			if(fabs(angux) <= toler) angux = 0.0;
		   			angux = acos(angux)*factr;
		   			if(angux >= (180.0-toler)) angux = 0.0;
		   		}
				if(iwrit ==1) {
					for(ievab=1; ievab<=nevab; ievab++)
						fprintf(fp16, "%lf\t", tensg[ievab]);
					fprintf(fp16, "\n");
				}
		   		coset =vectr[ielem][1];
		   		senot =vectr[ielem][2];
		   		fuerc[ielem] =-tensg[1]*coset-tensg[2]*senot;
				if(rolado >0) {
				   fuefl[ielem][1] =-fuerc[ielem]+fuepp[ielem][1]*coset-fuepp[ielem][2]*senot;
				   fuefl[ielem][2] =              fuepp[ielem][1]*senot+fuepp[ielem][2]*coset;
				   fuefl[ielem][3] = fuerc[ielem]+fuepp[ielem][3]*coset-fuepp[ielem][4]*senot;
				   fuefl[ielem][4] =              fuepp[ielem][3]*senot+fuepp[ielem][4]*coset;
			         }
		   	}
		   	else if (ndime==3){
				if(iwrit ==1) {
					for(ievab=1; ievab<=nevab; ievab++)
						fprintf(fp16, "%d %lf\t", ielem, tensg[ievab]);
					fprintf(fp16, "\n");
				}
				fuerc[ielem]=(-tensg[1]*vectr[ielem][1]-tensg[2]*vectr[ielem][2])*vectr[ielem][4]-tensg[3]*vectr[ielem][3];
				if(rolado >0) {
					coset=vectr[ielem][1];
					senot=vectr[ielem][2];
					cosef=vectr[ielem][3];
					senof=vectr[ielem][4];
					fuefl[ielem][1] =-fuerc[ielem]+fuepp[ielem][1]*coset*senof-fuepp[ielem][2]*senot*senof-fuepp[ielem][3]*cosef;
					fuefl[ielem][2] =              fuepp[ielem][1]*senot*senof+fuepp[ielem][2]*coset*senof-fuepp[ielem][3]*cosef;
					fuefl[ielem][3] =              fuepp[ielem][1]*cosef      +fuepp[ielem][2]*cosef      +fuepp[ielem][3]*senof;
					fuefl[ielem][4] =-fuerc[ielem]+fuepp[ielem][4]*coset*senof-fuepp[ielem][5]*senot*senof-fuepp[ielem][6]*cosef;
					fuefl[ielem][5] =              fuepp[ielem][4]*senot*senof+fuepp[ielem][5]*coset*senof-fuepp[ielem][6]*cosef;
					fuefl[ielem][6] =              fuepp[ielem][4]*cosef      +fuepp[ielem][5]*cosef      +fuepp[ielem][6]*senof;
				}
			}
		}

		if(ntipo >= 2) {
			if(ndime ==2){
				/*	fprintf(fp16, "globales en elemento %d\t",ielem);
				 for(ievab=1; ievab<=nevab; ievab++)
				 fprintf(fp16, "%lf\t", tensg[ievab]);
				 fprintf(fp16, "\n"); */
				coset = vectr[ielem][1];
				senot = vectr[ielem][2];
				fuefl[ielem][1] = tensg[1]*coset+tensg[2]*senot;
				fuefl[ielem][2] =-tensg[1]*senot+tensg[2]*coset;
				fuefl[ielem][3] = tensg[3];
				fuefl[ielem][4] = ( tensg[4]*coset+tensg[5]*senot);
				fuefl[ielem][5] = (-tensg[4]*senot+tensg[5]*coset);
				fuefl[ielem][6] = tensg[6];
				if(iwrit ==1) {
                   fprintf(fp16, "ejes locales en elemento %d\t",ielem);
			       for(ievab=1; ievab<=nevab; ievab++)
				       fprintf(fp16, "%lf\t", fuefl[ielem][ievab]);
			       fprintf(fp16, "\n");
				}
			}
			else if(ndime == 3){

				/*fprintf(fp16, "globales en elemento %d\t",ielem);
				 for(ievab=1; ievab<=nevab; ievab++)
				 fprintf(fp16, "%lf\t", tensg[ievab]);
				 fprintf(fp16, "\n"); //*/

				fuefl[ielem][1] = vectr[ielem][1]*tensg[1]+ vectr[ielem][2]*tensg[2]+ vectr[ielem][3]*tensg[3];
				fuefl[ielem][2] = vectr[ielem][4]*tensg[1]+ vectr[ielem][5]*tensg[2]+ vectr[ielem][6]*tensg[3];
				fuefl[ielem][3] = vectr[ielem][7]*tensg[1]+ vectr[ielem][8]*tensg[2]+ vectr[ielem][9]*tensg[3];
				fuefl[ielem][4] = vectr[ielem][1]*tensg[4] +vectr[ielem][2]*tensg[5] +vectr[ielem][3]*tensg[6];
				fuefl[ielem][5] = vectr[ielem][4]*tensg[4] +vectr[ielem][5]*tensg[5] +vectr[ielem][6]*tensg[6];
				fuefl[ielem][6] = vectr[ielem][7]*tensg[4] +vectr[ielem][8]*tensg[5] +vectr[ielem][9]*tensg[6];
				fuefl[ielem][7] = vectr[ielem][1]*tensg[7] +vectr[ielem][2]*tensg[8] +vectr[ielem][3]*tensg[9];
				fuefl[ielem][8] = vectr[ielem][4]*tensg[7] +vectr[ielem][5]*tensg[8] +vectr[ielem][6]*tensg[9];
				fuefl[ielem][9] = vectr[ielem][7]*tensg[7] +vectr[ielem][8]*tensg[8] +vectr[ielem][9]*tensg[9];
				fuefl[ielem][10]= vectr[ielem][1]*tensg[10]+vectr[ielem][2]*tensg[11]+vectr[ielem][3]*tensg[12];
				fuefl[ielem][11]= vectr[ielem][4]*tensg[10]+vectr[ielem][5]*tensg[11]+vectr[ielem][6]*tensg[12];
				fuefl[ielem][12]= vectr[ielem][7]*tensg[10]+vectr[ielem][8]*tensg[11]+vectr[ielem][9]*tensg[12];
				if(iwrit ==1) {
					fprintf(fp16, "Fue.Loc %d\t",ielem);
					for(ievab=1; ievab<=nevab; ievab++)
						fprintf(fp16, "%lf\t", fuefl[ielem][ievab]);
					fprintf(fp16, "\n");
				}
			}
		}
	}
	pesot=despt=0;
	for(ielem=1; ielem<=nelem; ielem++){
	/*	if(ndime ==2) pesob = props[matnu[ielem]][2]*props[matnu[ielem]][4]*xlong[ielem];
		if(ndime ==3) */
        pesob = props[matnu[ielem]][2]*props[matnu[ielem]][7]*xlong[ielem];
       // printf("ielem %d, matnu=%d, pesob=%lf \n", ielem, matnu[ielem],pesob);
		pesot +=pesob;
	}
    for(itotv=1; itotv<=ntotv; itotv++){
		if(despt < fabs(despl[itotv]))despt =fabs(despl[itotv]);
      //  printf("despt=%lf despl[%d]=%lf \n",despt,itotv,despl[itotv]);
    }
	valor[1]=pesot;
	// Verifica si el desplazamiento maximo de este caso de carga supera al anterior
	if( despt > valor[2] ) valor[2] = despt;
	//printf("\n\nPeso total: %lf       despl. max. = %lf \n\n", pesot,despt);
	/*
		if(ntipo == 1){
		fprintf(fp16,"\n\n\n barra      fuerza         esf. act.   eficiencia          area             longitud        peso          pesotot          volumen \n");
		pesot=0.;
		v=0.;
		for(ielem=1; ielem<=nelem; ielem++){
			if(ndime ==2)pesob = props[matnu[ielem]][2]*props[matnu[ielem]][4]*xlong[ielem];
			if(ndime ==3)pesob = props[matnu[ielem]][2]*props[matnu[ielem]][7]*xlong[ielem];
			pesot +=pesob;
			v+=props[matnu[ielem]][2]*xlong[ielem];
			Fi(fuerc[ielem], xlong[ielem], matnu[ielem],ielem);

			fprintf(fp16, "%d\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", ielem, fuerc[ielem],
					fuerc[ielem]/props[matnu[ielem]][2],esfac/esfre,
					props[matnu[ielem]][2],xlong[ielem],pesob,pesot,v);
		}
	}
	 */
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Reacci(int nelem, int ndime, int npnod, int ngdln, int nincl, int nevab, int ntotv, int iwrit, int isale,
			int** lnods, int* nodea, int* iffix, int* lincl,  double* despl, double* fixed, double*** srmat,
			double** carpp, double* react, double** girom, double** girtm, double** tempr, double** rigid,
			double* xincl)
{
	int    ipnod,nloca,ngush,igdln,mcont,ncont,nn1,nposn,mposn,ievab,jevab,
	ielem,ne1,inode,itotv;
	double treac[13],desel[13];

	for(itotv=1; itotv<=ntotv; itotv++) react[itotv]=0.0;

	for(ielem=1; ielem<=nelem; ielem++) {
		nn1 = nodea[ielem];
		ne1 = nn1*ngdln ;
		mposn=0;

		for(inode=1; inode<=nn1; inode++) {
			nposn = (lnods[ielem][inode]-1)*ngdln;

			for(igdln=1; igdln<=ngdln; igdln++){
				nposn++;
				desel[++mposn] = despl[nposn]-fixed[nposn];
			}
		}


		if(nincl >0)ApoyoInclinado(ielem,ndime,ngdln,nevab,isale, nincl,lnods,
								   lincl,girom,girtm,tempr,rigid,xincl);

		for(ievab=1; ievab<=ne1; ievab++){
			treac[ievab]=0.0;

			for(jevab=1; jevab<=ne1; jevab++)
				treac[ievab] += srmat[ielem][ievab][jevab]*desel[jevab];
		}

		mposn=0;

		for(inode=1; inode<=nn1; inode++) {
			nposn = (lnods[ielem][inode]-1)*ngdln ;

			for(igdln=1; igdln<=ngdln; igdln++){
				nposn++;
				mposn++;
				if(iffix[nposn]==1)
					react[nposn]+=treac[mposn]-carpp[ielem][mposn];
			}
		}
	}


	for(ipnod=1; ipnod<=npnod; ipnod ++) {
		mcont = 0;
		nloca = (ipnod-1)*ngdln;

		for(igdln=1; igdln<=ngdln; igdln ++) {
			ngush = nloca+igdln;

			if(iffix[ngush] > 0) mcont = mcont+1;
		}

		if(mcont > 0) {
			for(igdln=1; igdln<=ngdln; igdln ++) {
				ncont = nloca+igdln;
				treac[igdln] = react[ncont];
			}

            if(iwrit ==1){
				fprintf(fp16,"Reac.Nodo\t%d\t",ipnod);
				for(igdln=1; igdln<=ngdln; igdln ++)
					fprintf(fp16,"%le\t", treac[igdln]);
				fprintf(fp16, "\n");
			}
		}
	}
}
