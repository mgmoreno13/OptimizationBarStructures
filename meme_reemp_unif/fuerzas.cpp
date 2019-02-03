//*****************************
// Rutina para caluculo de fuerzas
//*****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "fuerzas.h"
#include "raros.h"
#include "rigidez.h"


//************************************************
void Fuerzas(int nelem, int npnod, int ndime, int ntipo, int nevab, int nnode, int ngdln,
	         int iwrit, int icarg, int** lnods, int* matnu, int** indfu, int* ntips,
	         double** coord, double* xlong, double** props, double** fuerb, double* wcarg,
	         double** carpp, double** fuepp, double** vectr, double** carga, double** fuemp,
	         bool Optimiza, int &npesp )
/******************************************************************************
 Cálculo de fuerzas nodales equivalentes.
 ******************************************************************************/
{
	int   ielem, ievab, nnodf, npunt, ndist, ntria, npara, npunty, ndisty,
	ipunt, jelem, lnod1, lnod2, idist, itria, ipara;
    int i;
    char letrero[100];			//Contenido no importante, solo para lectura de letreros en archivos de datos.
	double fuerx, fuerz, posi , xlon1, ylon1, xlon2, ylon2,
	coset, senot, xmom1, xmom2, reac2, reac1, reax1, reax2, posib,    f1,    f2,
	f3,    f4,    f5,    f6,    f7,     f8,   f9,   f10,   f11,
	f12, cargw, wtotl, fueaz, fuebz, fuecz, ylong, factr, uvalu,
	vvalu, wvalu, reai1, reai2, xloxy;
	char titulo[250];

	// Inicializa vector de fuerzas nodales equivalentes:
	for(ielem=1; ielem<=nelem; ielem++) {
		for(ievab=1; ievab<=4; ievab++)
			indfu[ievab][ielem]=0;
		for(ievab=1; ievab<=6; ievab++)
			fuerb[ievab][ielem]=0.0;
		for(ievab=1; ievab<=nevab; ievab++) {
			carga[ielem][ievab] = 0.0;
			fuemp[ielem][ievab] = 0.0;
		}
	}
    // Lee el titulo del caso de carga
	fscanf(fp5, "%s", titulo);
	if(iwrit ==1) fprintf(fp16, "CASO_DE_CARGA: %d\t%s\n", icarg, titulo);

	// Cargas en elementos:
	if(ndime == 2) {
		if(ntipo == 1) {
			fscanf(fp5, "%d %d %d %d %d %d ", &nnodf, &npunt, &ndist, &ntria,
				   &npara, &npesp);
			if(iwrit ==1) fprintf(fp16, "%d\t%d\t%d\t%d\t%d\t%d\n", nnodf, npunt, ndist,
					ntria, npara, npesp);

			if(nnodf >0) Fuerzas_Puntuales(nnodf,nelem,nnode,ngdln,iwrit,lnods, carga);

			if(npunt > 0) {
				for(ipunt=1; ipunt<=npunt; ipunt++) {
					fscanf(fp5, "%d %lf %lf", &jelem, &fuerz, &posi);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\t%lf\n", jelem, fuerz, posi);
					indfu[2][jelem]=1;
					fuerb[2][jelem]=fuerz;
					fuerb[4][jelem]=posi;
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					posib = xlong[jelem]-posi;
					reac1 =-fuerz*posib/xlong[jelem];
					reac2 =-fuerz*posi/xlong[jelem];
					f1 =-reac1*senot;
					f2 = reac1*coset;
					f3 =-reac2*senot;
					f4 = reac2*coset;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
				}
			}

			if(ndist > 0) {
				for(idist=1; idist<=ndist; idist++) {
					fscanf(fp5, "%d %lf", &jelem, &fuerz);
					indfu[4][jelem]=1;
					fuerb[6][jelem]=fuerz;
					fuerx=0;
					if(iwrit ==1) fprintf(fp16, "%d \t %lf \t %lf \n", jelem, fuerz,fuerx);
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					reac1 =-(fuerz*xlong[jelem]/2.0);
					reac2 = reac1;
					reax1 =-(fuerx*xlong[jelem]/2.0);
					reax2 = reax1;
					f1 =-reac1*senot+reax1*coset;
					f2 = reac1*coset+reax1*senot;
					f3 =-reac2*senot+reax2*coset;
					f4 = reac2*coset+reax2*senot;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
				}
			}

			if(ntria > 0) {
				for(itria=1; itria<=ntria; itria++) {
					fscanf(fp5, "%d %lf %lf", &jelem, &fueaz, &fuebz);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\t%lf\n", jelem, fueaz, fuebz);
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					reac1 =-(    fueaz+2.0*fuebz)*xlong[jelem]/6.0;
					reac2 =-(2.0*fueaz+    fuebz)*xlong[jelem]/6.0;
					f1 =-reac1*senot;
					f2 = reac1*coset;
					f3 =-reac2*senot;
					f4 = reac2*coset;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
				}
			}

			if(npara > 0) {
				for(ipara=1; ipara<=npara; ipara++) {
					fscanf(fp5, "%d %lf %lf %lf %lf", &jelem, &fueaz,
						   &fuebz, &fuecz, &ylong);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\t%lf\t%lf\t%lf\n", jelem, fueaz,
                         	fuebz, fuecz, ylong);
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					factr = fuecz*ylong*ylong+fueaz*(xlong[jelem]*xlong[jelem]-ylong*ylong)
					- fuebz*xlong[jelem]*xlong[jelem] ;
					uvalu = (fuebz-fueaz)/ylong/ylong-factr/(ylong*ylong*
															 xlong[jelem]*(ylong-xlong[jelem]));
					vvalu = factr/(ylong*xlong[jelem]*(ylong-xlong[jelem]));
					wvalu = fueaz  ;
					reac1 =-xlong[jelem]*(    uvalu*xlong[jelem]*xlong[jelem]+3.0*vvalu*xlong[jelem]+
								   6.0*wvalu)/12.0;
					reac2 =-xlong[jelem]*(3.0*uvalu*xlong[jelem]*xlong[jelem]+4.0*vvalu*xlong[jelem]+
								   6.0*wvalu)/12.0;
					f1 =-reac1*senot;
					f2 = reac1*coset;
					f3 =-reac2*senot;
					f4 = reac2*coset;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
				}
			}

			// Carga por peso propio en elementos tipo 1. Si la bandera Optimiza es true
			// no se carga el peso propio de los elementos porque en el optimizador este va variando.
			if(npesp == 1 && Optimiza == false) {
				for(ielem=1; ielem<=nelem; ielem++) {
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					coset = vectr[ielem][1];
					cargw= 1.0;
               		if (coset < 0.0) cargw= -1;
					indfu[4][ielem]=1;
					fuerb[6][ielem]+=-cargw*wtotl;
					reac1 = wtotl*xlong[ielem]/2.0;
					reac2 = reac1;
					f2 = reac1;
					f4 = reac2;
					fuemp[ielem][2] = fuemp[ielem][2]+f2;
					fuemp[ielem][4] = fuemp[ielem][4]+f4;
					carga[ielem][2] = carga[ielem][2]-f2;
					carga[ielem][4] = carga[ielem][4]-f4;
				}
			}
		}

		if(ntipo == 2) {
			fscanf(fp5, "%d %d %d %d %d %d", &nnodf, &npunt, &ndist, &ntria,
				   &npara, &npesp);
			if(iwrit ==1) fprintf(fp16, "%d\t%d\t%d\t%d\t%d\t%d\n", nnodf, npunt, ndist,
					ntria, npara, npesp);

			if(nnodf >0) Fuerzas_Puntuales(nnodf,nelem,nnode,ngdln,iwrit,lnods, carga);

			if(npunt > 0) {
				for(ipunt=1; ipunt<=npunt; ipunt++) {
					fscanf(fp5, " %d %lf %lf", &jelem, &fuerz, &posi);
					if(iwrit ==1) fprintf(fp16, " %d\t%lf\t%lf\n", jelem, fuerz, posi);
					indfu[2][jelem]=1;
					fuerb[2][jelem]=fuerz;
					fuerb[4][jelem]=posi;
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					posib = xlong[jelem]-posi;
					xmom1 =-fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
					xmom2 = fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
					reac1 = (+xmom1+xmom2-fuerz*posib)/xlong[jelem];
					reac2 =-( xmom2+xmom1+fuerz*posi)/xlong[jelem];
 					f1 =-reac1*senot;
					f2 = reac1*coset;
					f3 = xmom1;
					f4 =-reac2*senot;
					f5 = reac2*coset;
					f6 = xmom2;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					fuemp[jelem][5] = fuemp[jelem][5]+f5;
					fuemp[jelem][6] = fuemp[jelem][6]+f6;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
				}
			}

			if(ndist > 0) {
				for(idist=1; idist<=ndist; idist++) {
					fscanf(fp5, "%d %lf", &jelem, &fuerz);
					fuerx=0;
					if(iwrit ==1) fprintf(fp16, "%d \t% lf \t% lf \n", jelem, fuerz,fuerx);
					indfu[4][jelem]=1;
					fuerb[6][jelem]=fuerz;
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					xmom1 =-fuerz*xlong[jelem]*xlong[jelem]/12.0;
					xmom2 =-xmom1;
					reac1 =-(fuerz*xlong[jelem]/2.0);
					reac2 = reac1;
     				reax1 =-(fuerx*xlong[jelem]/2.0);
					reax2 = reax1;
					f1 =-reac1*senot+reax1*coset;
					f2 = reac1*coset+reax1*senot;
					f3 = xmom1;
					f4 =-reac2*senot+reax2*coset;
					f5 = reac2*coset+reax2*senot;
					f6 = xmom2;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					fuemp[jelem][5] = fuemp[jelem][5]+f5;
					fuemp[jelem][6] = fuemp[jelem][6]+f6;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
				}
			}

			if(ntria > 0) {
				for(itria=1; itria<=ntria; itria++) {
					fscanf(fp5, "%d %lf %lf", &jelem, &fueaz, &fuebz);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\t%lf\n", jelem, fueaz, fuebz);
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					xmom1 =-(9.0*fueaz+6.0*fuebz)*xlong[jelem]*xlong[jelem]/180.0;
					xmom2 = (6.0*fueaz+9.0*fuebz)*xlong[jelem]*xlong[jelem]/180.0;
					reac1 =-(63.0*fueaz+27.0*fuebz)*xlong[jelem]/180.0;
					reac2 =-(27.0*fueaz+63.0*fuebz)*xlong[jelem]/180.0;
					f1 =-reac1*senot;
					f2 = reac1*coset;
					f3 = xmom1;
					f4 =-reac2*senot;
					f5 = reac2*coset;
					f6 = xmom2;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					fuemp[jelem][5] = fuemp[jelem][5]+f5;
					fuemp[jelem][6] = fuemp[jelem][6]+f6;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
				}
			}

			if(npara > 0) {
				for(ipara=1; ipara<=npara; ipara++) {
					fscanf(fp5, "%d %lf %lf %lf %lf ", &jelem, &fueaz,
						   &fuebz, &fuecz, &ylong);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\t%lf\t%lf\t%lf\n", jelem, fueaz,
                         	fuebz, fuecz, ylong);
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					factr = fuecz*ylong*ylong+fueaz*(xlong[jelem]*xlong[jelem]-ylong*ylong)
					- fuebz*xlong[jelem]*xlong[jelem] ;
					uvalu = (fuebz-fueaz)/ylong/ylong-factr/(ylong*ylong*
															 xlong[jelem]*(ylong-xlong[jelem]));
					vvalu = factr/(ylong*xlong[jelem]*(ylong-xlong[jelem]));
					wvalu = fueaz  ;
					xmom1 =-( 4.0*uvalu*xlong[jelem]*xlong[jelem]*xlong[jelem]*xlong[jelem]+12.0*vvalu*
							 xlong[jelem]*xlong[jelem]*xlong[jelem]+30.0*wvalu*xlong[jelem]*xlong[jelem])/360.0;
					reai1 =-(uvalu*xlong[jelem]*xlong[jelem]+2.0*vvalu*xlong[jelem]+6.0*wvalu)*
					xlong[jelem]/12.0;
					xmom2 =-(18.0*uvalu*xlong[jelem]*xlong[jelem]*xlong[jelem]*xlong[jelem]+42.0*vvalu*
							 xlong[jelem]*xlong[jelem]*xlong[jelem]+150.0*wvalu*xlong[jelem]*xlong[jelem]+360.0*
							 reai1*xlong[jelem])/360.0;
					reai2 =-(3.0*uvalu*xlong[jelem]*xlong[jelem]+4.0*vvalu*xlong[jelem]+6.0*wvalu)
					*xlong[jelem]/12.0;
					reac1 = (xmom1+xmom2)/xlong[jelem]+reai1;
					reac2 =-(xmom1+xmom2)/xlong[jelem]+reai2;
					f1 =-reac1*senot;
					f2 = reac1*coset;
					f3 = xmom1;
					f4 =-reac2*senot;
					f5 = reac2*coset;
					f6 = xmom2;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					fuemp[jelem][5] = fuemp[jelem][5]+f5;
					fuemp[jelem][6] = fuemp[jelem][6]+f6;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
				}
			}

			// Carga por peso propio en elementos tipo 2. Si la bandera Optimiza es true
			// no se carga el peso propio de los elementos porque en el optimizador este va variando.
			if(npesp == 1 && Optimiza == false) {
				for(ielem=1; ielem<=nelem; ielem++) {
					coset = vectr[ielem][1];
					cargw= 1.0;
               		if (coset < 0.0) cargw= -1;
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					indfu[4][ielem]=1;
					fuerb[6][ielem]+=-cargw*wtotl*coset;
					xmom1 = wtotl*xlong[ielem]*xlong[ielem]*fabs(coset)*cargw/12.0;
					xmom2 =-xmom1;
					reac1 = wtotl*xlong[ielem]/2.0;
					reac2 = reac1;
					f2 = reac1;
					f3 = xmom1;
					f5 = reac2;
					f6 = xmom2;
					fuemp[ielem][2] = fuemp[ielem][2]+f2;
					fuemp[ielem][3] = fuemp[ielem][3]+f3;
					fuemp[ielem][5] = fuemp[ielem][5]+f5;
					fuemp[ielem][6] = fuemp[ielem][6]+f6;
					carga[ielem][2] = carga[ielem][2]-f2;
					carga[ielem][3] = carga[ielem][3]-f3;
					carga[ielem][5] = carga[ielem][5]-f5;
					carga[ielem][6] = carga[ielem][6]-f6;
				}
			}
		}

		if(ntipo == 3) {
			fscanf(fp5, "%d %d %d %d %d %d", &nnodf, &npunt, &ndist, &ntria,
				   &npara, &npesp);
			if(iwrit ==1) fprintf(fp16, "%d\t%d\t%d\t%d\t%d\t%d\n", nnodf, npunt, ndist,
					ntria, npara, npesp);

			if(nnodf >0) Fuerzas_Puntuales(nnodf,nelem,nnode,ngdln,iwrit,lnods, carga);

			if(npunt > 0) {
				for(ipunt=1; ipunt<=npunt; ipunt++) {
					fscanf(fp5, "%d %lf %lf", &jelem, &fuerz, &posi);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\t%lf\n", jelem, fuerz, posi);
					indfu[1][jelem]=1;
					fuerb[2][jelem]=fuerz;
					fuerb[4][jelem]=posi;
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					posib = xlong[jelem]-posi;

					if(ntips[jelem] == 1) {
						reac1 =-fuerz*posib/xlong[jelem];
						reac2 =-fuerz*posi/xlong[jelem];
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
					}

					if(ntips[jelem] == 2) {
						xmom1 =-fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
						xmom2 = fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
						reac1 = (+xmom1+xmom2-fuerz*posib)/xlong[jelem];
						reac2 =-( xmom2+xmom1+fuerz*posi)/xlong[jelem];
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 = xmom1;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						f6 = xmom2;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][3] = carga[jelem][3]-f3;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
						carga[jelem][6] = carga[jelem][6]-f6;
					}

					if(ntips[jelem] == 3) {
						xmom1 =-fuerz*posi*posib*(posib+posi/2)/(xlong[jelem]*
																 xlong[jelem]);
						reac1= (xmom1-fuerz*posib)/xlong[jelem];
						reac2=-(fuerz+reac1);
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 = xmom1;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][3] = carga[jelem][3]-f3;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
					}
					if(ntips[jelem] == 4) {
						xmom2 =fuerz*posi*posib*(posi+posib/2)/(xlong[jelem]*
																 xlong[jelem]);
						reac1= (xmom2-fuerz*posi)/xlong[jelem];
						reac2= -(fuerz+reac1);
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						f6 = xmom2;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
						carga[jelem][6] = carga[jelem][6]-f6;
					}
                }
			}

			if(ndist > 0) {
				for(idist=1; idist<=ndist; idist++) {
					fscanf(fp5, "%d %lf", &jelem, &fuerz);
					fuerx=0;
					if(iwrit ==1) fprintf(fp16, "%d \t %lf  \n", jelem, fuerz);
					indfu[4][jelem]=1;
					fuerb[6][jelem]=fuerz;
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];

					if(ntips[jelem] == 1) {
						reac1 =-(fuerz*xlong[jelem]/2.0);
						reac2 = reac1;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
					}

					if(ntips[jelem] == 2) {
						xmom1 =-fuerz*xlong[jelem]*xlong[jelem]/12.0;
						xmom2 =-xmom1;
						reac1 =-(fuerz*xlong[jelem]/2.0);
						reac2 = reac1;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 = xmom1;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						f6 = xmom2;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][3] = carga[jelem][3]-f3;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
						carga[jelem][6] = carga[jelem][6]-f6;
					}

					if(ntips[jelem] == 3) {
						xmom1 =-fuerz*xlong[jelem]*xlong[jelem]/8.0;
						reac1 =-(5.0*fuerz*xlong[jelem]/8.0);
						reac2 = reac1*.6;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 = xmom1;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][3] = carga[jelem][3]-f3;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
					}
					if(ntips[jelem] == 4) {
						xmom2 =fuerz*xlong[jelem]*xlong[jelem]/8.0;
						reac2 =-(5.0*fuerz*xlong[jelem]/8.0);
						reac1 = reac2*.6;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						f6 = xmom2;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
						carga[jelem][6] = carga[jelem][6]-f6;
					}
				}
			}

			if(ntria > 0){
				for(itria=1; itria<=ntria; itria++) {
					fscanf(fp5, "%d %lf %lf", &jelem, &fueaz, &fuebz);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\t%lf\n", jelem, fueaz, fuebz);
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					if(ntips[jelem] == 1) {
						reac1 =-(    fueaz+2.0*fuebz)*xlong[jelem]/6.0;
						reac2 =-(2.0*fueaz+    fuebz)*xlong[jelem]/6.0;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 =-reac2*senot;
						f4 = reac2*coset;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][4] = fuemp[jelem][4]+f3;
						fuemp[jelem][5] = fuemp[jelem][5]+f4;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][4] = carga[jelem][4]-f3;
						carga[jelem][5] = carga[jelem][5]-f4;
					}
					if(ntips[jelem] == 2) {
						xmom1 =-(9.0*fueaz+6.0*fuebz)*xlong[jelem]*xlong[jelem]/180.0;
						xmom2 = (6.0*fueaz+9.0*fuebz)*xlong[jelem]*xlong[jelem]/180.0;
						reac1 =-(63.0*fueaz+27.0*fuebz)*xlong[jelem]/180.0;
						reac2 =-(27.0*fueaz+63.0*fuebz)*xlong[jelem]/180.0;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 = xmom1;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						f6 = xmom2;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][3] = carga[jelem][3]-f3;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
						carga[jelem][6] = carga[jelem][6]-f6;
					}
					if(ntips[jelem] == 3) {
						xmom1 =-( 8.0*fueaz+ 7.0*fuebz)*xlong[jelem]*xlong[jelem]/120.0;
						reac1 =-(16.0*fueaz+ 9.0*fuebz)*xlong[jelem]/40.0;
						reac2 =-( 4.0*fueaz+11.0*fuebz)*xlong[jelem]/40.0;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 = xmom1;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][3] = carga[jelem][3]-f3;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
					}
				}
			}

			if(npara > 0){
				for(ipara=1; ipara<=npara; ipara++) {
					fscanf(fp5, "%d %lf %lf %lf %lf ", &jelem, &fueaz,
						   &fuebz, &fuecz, &ylong);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\t%lf\t%lf\t%lf\n", jelem, fueaz,
                         	fuebz, fuecz, ylong);
					coset = vectr[jelem][1];
					senot = vectr[jelem][2];
					factr = fuecz*ylong*ylong+fueaz*(xlong[jelem]*xlong[jelem]-ylong*ylong)
					- fuebz*xlong[jelem]*xlong[jelem] ;
					uvalu = (fuebz-fueaz)/ylong/ylong-factr/(ylong*ylong*
															 xlong[jelem]*(ylong-xlong[jelem]));
					vvalu = factr/(ylong*xlong[jelem]*(ylong-xlong[jelem]));
					wvalu = fueaz  ;

					if(ntips[jelem] == 1) {
						reac1 =-xlong[jelem]*(    uvalu*xlong[jelem]*xlong[jelem]+3.0*vvalu*
									   xlong[jelem]+6.0*wvalu)/12.0;
						reac2 =-xlong[jelem]*(3.0*uvalu*xlong[jelem]*xlong[jelem]+4.0*vvalu*
									   xlong[jelem]+6.0*wvalu)/12.0;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 =-reac2*senot;
						f4 = reac2*coset;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][4] = fuemp[jelem][4]+f3;
						fuemp[jelem][5] = fuemp[jelem][5]+f4;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][4] = carga[jelem][4]-f3;
						carga[jelem][5] = carga[jelem][5]-f4;
					}

					if(ntips[jelem] == 2) {
						xmom1 =-(4.0*uvalu*xlong[jelem]*xlong[jelem]*xlong[jelem]*xlong[jelem]+12.0*
								 vvalu*xlong[jelem]*xlong[jelem]*xlong[jelem]+30.0*wvalu*xlong[jelem]*
								 xlong[jelem])/360.0;
						reai1 =-(uvalu*xlong[jelem]*xlong[jelem]+2.0*vvalu*xlong[jelem]+6.0*
								 wvalu)*xlong[jelem]/12.0;
						xmom2 =-(18.0*uvalu*xlong[jelem]*xlong[jelem]*xlong[jelem]*xlong[jelem]+42.0*
								 vvalu*xlong[jelem]*xlong[jelem]*xlong[jelem]+150.0*wvalu*xlong[jelem]*
								 xlong[jelem]+360.0*reai1*xlong[jelem])/360.0;
						reai2 =-(3.0*uvalu*xlong[jelem]*xlong[jelem]+4.0*vvalu*xlong[jelem]+6.0
								 *wvalu)*xlong[jelem]/12.0;
						reac1 = (xmom1+xmom2)/xlong[jelem]+reai1;
						reac2 =-(xmom1+xmom2)/xlong[jelem]+reai2;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 = xmom1;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						f6 = xmom2;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][3] = carga[jelem][3]-f3;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
						carga[jelem][6] = carga[jelem][6]-f6;
					}

					if(ntips[jelem] == 3) {
						reac1 =-(uvalu*xlong[jelem]*xlong[jelem]+2.0*vvalu*xlong[jelem]+6.0*
								 wvalu)*xlong[jelem]/12.0;
						xmom1 = (uvalu*xlong[jelem]*xlong[jelem]*xlong[jelem]*xlong[jelem]+3.0*vvalu*
								 xlong[jelem]*xlong[jelem]*xlong[jelem]+15.0*wvalu*xlong[jelem]*xlong[jelem])/
						120.0+reac1*xlong[jelem]/2.0;
						reac2 = -(uvalu*xlong[jelem]*xlong[jelem]+3.0*vvalu*xlong[jelem]+15.0*
								  wvalu)*xlong[jelem]/120.0-reac1/2.0-(3.0*uvalu*xlong[jelem]
																*xlong[jelem]+4.0*vvalu*xlong[jelem]+6.0*wvalu)*xlong[jelem]/12.0;
						reac1 = reac1*3.0/2.0 +(uvalu*xlong[jelem]*xlong[jelem]*xlong[jelem]+
												3.0*vvalu*xlong[jelem]*xlong[jelem]+15.0*wvalu*xlong[jelem])/120.0;
						f1 =-reac1*senot;
						f2 = reac1*coset;
						f3 = xmom1;
						f4 =-reac2*senot;
						f5 = reac2*coset;
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						carga[jelem][1] = carga[jelem][1]-f1;
						carga[jelem][2] = carga[jelem][2]-f2;
						carga[jelem][3] = carga[jelem][3]-f3;
						carga[jelem][4] = carga[jelem][4]-f4;
						carga[jelem][5] = carga[jelem][5]-f5;
					}
				}
			}

			// Carga por peso propio en elementos subtipos 1, 2 y 3. Si la bandera Optimiza es true
			// no se carga el peso propio de los elementos porque en el optimizador este va variando.
			if(npesp == 1 && Optimiza == false) {
				for(ielem=1; ielem<=nelem; ielem++) {
					coset = vectr[ielem][1];
               		cargw= 1.0;
               		if(coset<0.0) cargw= -1;
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					indfu[4][ielem]=1;
					fuerb[6][ielem]+=-cargw*wtotl*coset;
					if(ntips[ielem] == 1) {
						reac1 = (wtotl*xlong[ielem]/2.0);
						reac2 = reac1;
						f2 = reac1;
						f5 = reac2;
						fuemp[ielem][2] = fuemp[ielem][2]+f2;
						fuemp[ielem][5] = fuemp[ielem][5]+f5;
						carga[ielem][2] = carga[ielem][2]-f2;
						carga[ielem][5] = carga[ielem][5]-f5;
					}

					if(ntips[ielem] == 2) {
						xmom1 = wtotl*xlong[ielem]*xlong[ielem]*fabs(coset)*cargw/12.0;
						xmom2 =-xmom1;
						reac1 = wtotl*xlong[ielem]/2.0;
						reac2 = reac1;
						f2 = reac1;
						f3 = xmom1;
						f5 = reac2;
						f6 = xmom2;
						fuemp[ielem][2] = fuemp[ielem][2]+f2;
						fuemp[ielem][3] = fuemp[ielem][3]+f3;
						fuemp[ielem][5] = fuemp[ielem][5]+f5;
						fuemp[ielem][6] = fuemp[ielem][6]+f6;
						carga[ielem][2] = carga[ielem][2]-f2;
						carga[ielem][3] = carga[ielem][3]-f3;
						carga[ielem][5] = carga[ielem][5]-f5;
						carga[ielem][6] = carga[ielem][6]-f6;
					}

					if(ntips[ielem] == 3) {
						xmom1 = wtotl*xlong[ielem]*xlong[ielem]*fabs(coset)*cargw/8.0;
						reac1 = 5.0*wtotl*xlong[ielem]/8.0;
						reac2 = 3.0*wtotl*xlong[ielem]/8.0;
						f2 = reac1;
						f3 = xmom1;
						f5 = reac2;
						fuemp[ielem][2] = fuemp[ielem][2]+f2;
						fuemp[ielem][3] = fuemp[ielem][3]+f3;
						fuemp[ielem][5] = fuemp[ielem][5]+f5;
						carga[ielem][2] = carga[ielem][2]-f2;
						carga[ielem][3] = carga[ielem][3]-f3;
						carga[ielem][5] = carga[ielem][5]-f5;
					}
					if(ntips[ielem] == 4) {
						xmom2 = -wtotl*xlong[ielem]*xlong[ielem]*fabs(coset)*cargw/8.0;
						reac1 = 3.0*wtotl*xlong[ielem]/8.0;
						reac2 = 5.0*wtotl*xlong[ielem]/8.0;
						f2 = reac1;
						f4 = reac2;
						f5 = xmom1;
						fuemp[ielem][2] = fuemp[ielem][2]+f2;
						fuemp[ielem][5] = fuemp[ielem][5]+f5;
						fuemp[ielem][6] = fuemp[ielem][6]+f6;
						carga[ielem][2] = carga[ielem][2]-f2;
						carga[ielem][3] = carga[ielem][3]-f4;
						carga[ielem][5] = carga[ielem][5]-f5;
					}
				}
			}
		}
	}

	if(ndime == 3) {
		if(ntipo == 1) {
			fscanf(fp5, "%d %d %d %d %d %d ", &nnodf, &npunt, &ndist, &ntria,
				   &npara, &npesp);
			if(iwrit ==1) fprintf(fp16, "%d\t%d\t%d\t%d\t%d\t%d\n", nnodf, npunt, ndist,
					ntria, npara, npesp);

			if(nnodf >0) Fuerzas_Puntuales(nnodf,nelem,nnode,ngdln,iwrit,lnods, carga);

            // Si la bandera Optimiza es true no se carga el peso propio de los elementos
            // porque en el optimizador este va variando.
			if(npesp == 1 && Optimiza == false) {
				for(ielem=1; ielem<=nelem; ielem++) {
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					reac1 = (wtotl*xlong[ielem]/2.0);
					reac2 = reac1;
					indfu[4][ielem]=1;
					fuerb[6][ielem]=-wtotl*vectr[ielem][3];
					indfu[3][ielem]=1;
					fuerb[5][ielem]=0;
					f3 = reac1;
					f6 = reac2;
					fuemp[ielem][3] = fuemp[ielem][3]+f3;
					fuemp[ielem][6] = fuemp[ielem][6]+f6;
					carga[ielem][3] = carga[ielem][3]-f3;
					carga[ielem][6] = carga[ielem][6]-f6;
				}
			}
		}

		if(ntipo == 2) {
			fscanf(fp5, "%d %d %d %d %d %d", &nnodf, &npunt, &npunty,
				   &ndist, &ndisty, &npesp);
			if(iwrit ==1) fprintf(fp16, "%d\t%d\t%d\t%d\t%d\t%d\n", nnodf, npunt, npunty,
					ndist, ndisty, npesp);
			if(nnodf >0) Fuerzas_Puntuales(nnodf,nelem,nnode,ngdln,iwrit,lnods, carga);
			if(npunt > 0) {
				for(ipunt=1; ipunt<=npunt; ipunt++) {
					fscanf(fp5, " %d %lf %lf", &jelem, &fuerz, &posi);
					if(iwrit ==1) fprintf(fp16, " %d\t%lf\t%lf\n", jelem, fuerz, posi);
					indfu[1][jelem]=1;
					fuerb[1][jelem]=fuerz;
					fuerb[3][jelem]=posi;
					posib = xlong[jelem]-posi;
					xmom1 =fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
					xmom2 =-fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
					reac1 = (+xmom1+xmom2-fuerz*posib)/xlong[jelem];
					reac2 =-( xmom2+xmom1+fuerz*posi)/xlong[jelem];
					f1 =vectr[jelem][7]*reac1;
					f2 =vectr[jelem][8]*reac1;
					f3 =vectr[jelem][9]*reac1;
 					f4 =vectr[jelem][4]*xmom1;
					f5 =vectr[jelem][5]*xmom1;
					f6 =vectr[jelem][6]*xmom1;
					f7 =vectr[jelem][7]*reac2;
					f8 =vectr[jelem][8]*reac2;
					f9 =vectr[jelem][9]*reac2;
 					f10=vectr[jelem][4]*xmom2;
					f11=vectr[jelem][5]*xmom2;
					f12=vectr[jelem][6]*xmom2;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					fuemp[jelem][5] = fuemp[jelem][5]+f5;
					fuemp[jelem][6] = fuemp[jelem][6]+f6;
					fuemp[jelem][7] = fuemp[jelem][7]+f7;
					fuemp[jelem][8] = fuemp[jelem][8]+f8;
					fuemp[jelem][9] = fuemp[jelem][9]+f9;
					fuemp[jelem][10]= fuemp[jelem][10]+f10;
					fuemp[jelem][11]= fuemp[jelem][11]+f11;
					fuemp[jelem][12]= fuemp[jelem][12]+f12;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
					carga[jelem][7] = carga[jelem][7]-f7;
					carga[jelem][8] = carga[jelem][8]-f8;
					carga[jelem][9] = carga[jelem][9]-f9;
					carga[jelem][10]= carga[jelem][10]-f10;
					carga[jelem][11]= carga[jelem][11]-f11;
					carga[jelem][12]= carga[jelem][12]-f12;
				}
			}
			if(npunty > 0) {
				for(ipunt=1; ipunt<=npunty; ipunt++) {
					fscanf(fp5, " %d %lf %lf", &jelem, &fuerz, &posi);
					if(iwrit ==1) fprintf(fp16, " %d\t%lf\t%lf\n", jelem, fuerz, posi);
					indfu[1][jelem]=2;
					fuerb[2][jelem]=fuerz;
					fuerb[4][jelem]=posi;
					posib = xlong[jelem]-posi;
					xmom1 =-fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
					xmom2 = fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
					reac1 = (+xmom1+xmom2-fuerz*posib)/xlong[jelem];
					reac2 =-( xmom2+xmom1+fuerz*posi)/xlong[jelem];

					f1 =vectr[jelem][4]*reac1;
					f2 =vectr[jelem][5]*reac1;
					f3 =vectr[jelem][6]*reac1;
 					f4 =vectr[jelem][7]*xmom1;
					f5 =vectr[jelem][8]*xmom1;
					f6 =vectr[jelem][9]*xmom1;
					f7 =vectr[jelem][4]*reac2;
					f8 =vectr[jelem][5]*reac2;
					f9 =vectr[jelem][6]*reac2;
 					f10=vectr[jelem][7]*xmom2;
					f11=vectr[jelem][8]*xmom2;
					f12=vectr[jelem][9]*xmom2;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					fuemp[jelem][5] = fuemp[jelem][5]+f5;
					fuemp[jelem][6] = fuemp[jelem][6]+f6;
					fuemp[jelem][7] = fuemp[jelem][7]+f7;
					fuemp[jelem][8] = fuemp[jelem][8]+f8;
					fuemp[jelem][9] = fuemp[jelem][9]+f9;
					fuemp[jelem][10]= fuemp[jelem][10]+f10;
					fuemp[jelem][11]= fuemp[jelem][11]+f11;
					fuemp[jelem][12]= fuemp[jelem][12]+f12;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
					carga[jelem][7] = carga[jelem][7]-f7;
					carga[jelem][8] = carga[jelem][8]-f8;
					carga[jelem][9] = carga[jelem][9]-f9;
					carga[jelem][10]= carga[jelem][10]-f10;
					carga[jelem][11]= carga[jelem][11]-f11;
					carga[jelem][12]= carga[jelem][12]-f12;
				}
			}

			if(ndist > 0) {
				for(idist=1; idist<=ndist; idist++) {
					fscanf(fp5, "%d %lf", &jelem, &fuerz);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\n", jelem, fuerz);
					indfu[4][jelem]=1;
					fuerb[6][jelem]=fuerz;
					xmom1 = fuerz*xlong[jelem]*xlong[jelem]/12.0;
					xmom2 =-xmom1;
					reac1 =-(fuerz*xlong[jelem]/2.0);
					reac2 = reac1;
					f1 =vectr[jelem][7]*reac1;
					f2 =vectr[jelem][8]*reac1;
					f3 =vectr[jelem][9]*reac1;
 					f4 =vectr[jelem][4]*xmom1;
					f5 =vectr[jelem][5]*xmom1;
					f6 =vectr[jelem][6]*xmom1;
					f7 =vectr[jelem][7]*reac2;
					f8 =vectr[jelem][8]*reac2;
					f9 =vectr[jelem][9]*reac2;
 					f10=vectr[jelem][4]*xmom2;
					f11=vectr[jelem][5]*xmom2;
					f12=vectr[jelem][6]*xmom2;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					fuemp[jelem][5] = fuemp[jelem][5]+f5;
					fuemp[jelem][6] = fuemp[jelem][6]+f6;
					fuemp[jelem][7] = fuemp[jelem][7]+f7;
					fuemp[jelem][8] = fuemp[jelem][8]+f8;
					fuemp[jelem][9] = fuemp[jelem][9]+f9;
					fuemp[jelem][10]= fuemp[jelem][10]+f10;
					fuemp[jelem][11]= fuemp[jelem][11]+f11;
					fuemp[jelem][12]= fuemp[jelem][12]+f12;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
					carga[jelem][7] = carga[jelem][7]-f7;
					carga[jelem][8] = carga[jelem][8]-f8;
					carga[jelem][9] = carga[jelem][9]-f9;
					carga[jelem][10]= carga[jelem][10]-f10;
					carga[jelem][11]= carga[jelem][11]-f11;
					carga[jelem][12]= carga[jelem][12]-f12;
				}
			}
			if(ndisty > 0) {
				for(idist=1; idist<=ndisty; idist++) {
					fscanf(fp5, "%d %lf", &jelem, &fuerz);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\n", jelem, fuerz);
					indfu[3][jelem]=1;
					fuerb[5][jelem]=fuerz;
					xmom1 =-fuerz*xlong[jelem]*xlong[jelem]/12.0;
					xmom2 =-xmom1;
					reac1 =-(fuerz*xlong[jelem]/2.0);
					reac2 = reac1;
					f1 =vectr[jelem][4]*reac1;
					f2 =vectr[jelem][5]*reac1;
					f3 =vectr[jelem][6]*reac1;
 					f4 =vectr[jelem][7]*xmom1;
					f5 =vectr[jelem][8]*xmom1;
					f6 =vectr[jelem][9]*xmom1;
					f7 =vectr[jelem][4]*reac2;
					f8 =vectr[jelem][5]*reac2;
					f9 =vectr[jelem][6]*reac2;
 					f10=vectr[jelem][7]*xmom2;
					f11=vectr[jelem][8]*xmom2;
					f12=vectr[jelem][9]*xmom2;
					fuemp[jelem][1] = fuemp[jelem][1]+f1;
					fuemp[jelem][2] = fuemp[jelem][2]+f2;
					fuemp[jelem][3] = fuemp[jelem][3]+f3;
					fuemp[jelem][4] = fuemp[jelem][4]+f4;
					fuemp[jelem][5] = fuemp[jelem][5]+f5;
					fuemp[jelem][6] = fuemp[jelem][6]+f6;
					fuemp[jelem][7] = fuemp[jelem][7]+f7;
					fuemp[jelem][8] = fuemp[jelem][8]+f8;
					fuemp[jelem][9] = fuemp[jelem][9]+f9;
					fuemp[jelem][10]= fuemp[jelem][10]+f10;
					fuemp[jelem][11]= fuemp[jelem][11]+f11;
					fuemp[jelem][12]= fuemp[jelem][12]+f12;
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
					carga[jelem][7] = carga[jelem][7]-f7;
					carga[jelem][8] = carga[jelem][8]-f8;
					carga[jelem][9] = carga[jelem][9]-f9;
					carga[jelem][10]= carga[jelem][10]-f10;
					carga[jelem][11]= carga[jelem][11]-f11;
					carga[jelem][12]= carga[jelem][12]-f12;
				}
			}
			// Si la bandera Optimiza es true no se carga el peso propio de
			// los elementos porque en el optimizador este va variando.
			if(npesp == 1 && Optimiza == false) {
				for(ielem=1; ielem<=nelem; ielem++) {
					lnod1 = lnods[ielem][1];
					lnod2 = lnods[ielem][2];
					xlon1 = coord[lnod1][1];
					ylon1 = coord[lnod1][2];
					xlon2 = coord[lnod2][1];
					ylon2 = coord[lnod2][2];
					xloxy= sqrt((xlon2-xlon1)*(xlon2-xlon1)+(ylon2-ylon1)*
								(ylon2-ylon1));
					if(xloxy >0){
						coset=(xlon2-xlon1)/xloxy;
						senot=(ylon2-ylon1)/xloxy;
					}
					else {
						coset=senot=0.0;
					}
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					indfu[4][ielem]=1;
					fuerb[6][ielem]=-wtotl*vectr[ielem][9];
					indfu[3][ielem]=1;
					fuerb[5][ielem]=-wtotl*vectr[ielem][8];
					xmom1 =-wtotl*xlong[ielem]*xloxy/12.0;
					xmom2 =-xmom1;
					reac1 = wtotl*xlong[ielem]/2.0;
					reac2 = reac1;
					f1 =0.0;
					f2 =0.0;
					f3 =reac1;
 					f4 =-senot*xmom1;
					f5 = coset*xmom1;
					f6 =0.0;
					f7 =0.0;
					f8 =0.0;
					f9 =reac2;
 					f10=-senot*xmom2;
					f11= coset*xmom2;
					f12=0.0;
					fuemp[ielem][1] = fuemp[ielem][1]+f1;
					fuemp[ielem][2] = fuemp[ielem][2]+f2;
					fuemp[ielem][3] = fuemp[ielem][3]+f3;
					fuemp[ielem][4] = fuemp[ielem][4]+f4;
					fuemp[ielem][5] = fuemp[ielem][5]+f5;
					fuemp[ielem][6] = fuemp[ielem][6]+f6;
					fuemp[ielem][7] = fuemp[ielem][7]+f7;
					fuemp[ielem][8] = fuemp[ielem][8]+f8;
					fuemp[ielem][9] = fuemp[ielem][9]+f9;
					fuemp[ielem][10]= fuemp[ielem][10]+f10;
					fuemp[ielem][11]= fuemp[ielem][11]+f11;
					fuemp[ielem][12]= fuemp[ielem][12]+f12;
					carga[ielem][1] = carga[ielem][1]-f1;
					carga[ielem][2] = carga[ielem][2]-f2;
					carga[ielem][3] = carga[ielem][3]-f3;
					carga[ielem][4] = carga[ielem][4]-f4;
					carga[ielem][5] = carga[ielem][5]-f5;
					carga[ielem][6] = carga[ielem][6]-f6;
					carga[ielem][7] = carga[ielem][7]-f7;
					carga[ielem][8] = carga[ielem][8]-f8;
					carga[ielem][9] = carga[ielem][9]-f9;
					carga[ielem][10]= carga[ielem][10]-f10;
					carga[ielem][11]= carga[ielem][11]-f11;
					carga[ielem][12]= carga[ielem][12]-f12;
				}
			}
		}
		if(ntipo == 3) {
			fscanf(fp5, "%d %d %d %d %d %d", &nnodf, &npunt, &npunty,
				   &ndist, &ndisty, &npesp);
			if(iwrit ==1) fprintf(fp16, "%d\t%d\t%d\t%d\t%d\t%d\n", nnodf, npunt, npunty,
					ndist, ndisty, npesp);
			if(nnodf >0) Fuerzas_Puntuales(nnodf,nelem,nnode,ngdln,iwrit,lnods, carga);
			if(npunt > 0) {
				for(ipunt=1; ipunt<=npunt; ipunt++) {
					fscanf(fp5, " %d %lf %lf", &jelem, &fuerz, &posi);
					if(iwrit ==1) fprintf(fp16, " %d\t%lf\t%lf\n", jelem, fuerz, posi);
					indfu[1][jelem]=1;
					fuerb[1][jelem]=fuerz;
					fuerb[3][jelem]=posi;
					posib = xlong[jelem]-posi;
                    if(ntips[jelem] ==1){xmom1 = xmom2 = 0;}
                    if(ntips[jelem] ==2){
						xmom1 =fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
						xmom2 =-fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
					}
                    if(ntips[jelem] ==3){
						xmom1 =fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
						xmom2 =-fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
						xmom1 =xmom1-xmom2/2.0;
						xmom2 =0;
					}
                    if(ntips[jelem] ==4){
						xmom1 =fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
						xmom2 =-fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
						xmom2 =-xmom1/2.0+xmom2;
						xmom1 = 0;
					}
					reac1 =-(+xmom1+xmom2+fuerz*posib)/xlong[jelem];
					reac2 =+( xmom2+xmom1-fuerz*posi)/xlong[jelem];
					f1 =vectr[jelem][7]*reac1;
					f2 =vectr[jelem][8]*reac1;
					f3 =vectr[jelem][9]*reac1;
 					f4 =vectr[jelem][4]*xmom1;
					f5 =vectr[jelem][5]*xmom1;
					f6 =vectr[jelem][6]*xmom1;
					f7 =vectr[jelem][7]*reac2;
					f8 =vectr[jelem][8]*reac2;
					f9 =vectr[jelem][9]*reac2;
 					f10=vectr[jelem][4]*xmom2;
					f11=vectr[jelem][5]*xmom2;
					f12=vectr[jelem][6]*xmom2;
					// PENDIENTE: verificar que este cambio es correcto
					//if(ntips[ielem]!=1){
					if(ntips[jelem]!=1){
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						fuemp[jelem][7] = fuemp[jelem][7]+f7;
						fuemp[jelem][8] = fuemp[jelem][8]+f8;
						fuemp[jelem][9] = fuemp[jelem][9]+f9;
						fuemp[jelem][10]= fuemp[jelem][10]+f10;
						fuemp[jelem][11]= fuemp[jelem][11]+f11;
						fuemp[jelem][12]= fuemp[jelem][12]+f12;
					}
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
					carga[jelem][7] = carga[jelem][7]-f7;
					carga[jelem][8] = carga[jelem][8]-f8;
					carga[jelem][9] = carga[jelem][9]-f9;
					carga[jelem][10]= carga[jelem][10]-f10;
					carga[jelem][11]= carga[jelem][11]-f11;
					carga[jelem][12]= carga[jelem][12]-f12;
				}
			}
			if(npunty > 0) {
				for(ipunt=1; ipunt<=npunty; ipunt++) {
					fscanf(fp5, " %d %lf %lf", &jelem, &fuerz, &posi);
					if(iwrit ==1) fprintf(fp16, " %d\t%lf\t%lf\n", jelem, fuerz, posi);
					indfu[2][jelem]=1;
					fuerb[2][jelem]=fuerz;
					fuerb[4][jelem]=posi;
					lnod1 = lnods[jelem][1];
					posib = xlong[jelem]-posi;
                    if(ntips[jelem] ==1){xmom1 = xmom2 = 0;}
                    if(ntips[jelem] ==2){
						xmom1 =-fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
						xmom2 = fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
					}
                    if(ntips[jelem] ==3){
						xmom1 =-fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
						xmom2 = fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
						xmom1 = xmom1-xmom2/2;
						xmom2 = 0;
					}
                    if(ntips[jelem] ==4){
						xmom1 =fuerz*posi*posib*posib/(xlong[jelem]*xlong[jelem]);
						xmom2 =-fuerz*posi*posi*posib/(xlong[jelem]*xlong[jelem]);
						xmom2 =-xmom1/2.0+xmom2;
						xmom1 = 0;
					}
					reac1 = ( xmom1+xmom2-fuerz*posib)/xlong[jelem];
					reac2 =-( xmom2+xmom1+fuerz*posi)/xlong[jelem];
					f1 =vectr[jelem][4]*reac1;
					f2 =vectr[jelem][5]*reac1;
					f3 =vectr[jelem][6]*reac1;
 					f4 =vectr[jelem][7]*xmom1;
					f5 =vectr[jelem][8]*xmom1;
					f6 =vectr[jelem][9]*xmom1;
					f7 =vectr[jelem][4]*reac2;
					f8 =vectr[jelem][5]*reac2;
					f9 =vectr[jelem][6]*reac2;
 					f10=vectr[jelem][7]*xmom2;
					f11=vectr[jelem][8]*xmom2;
					f12=vectr[jelem][9]*xmom2;
					// PENDIENTE: verificar que este cambio es correcto
					//if(ntips[ielem]!=1){
					if(ntips[jelem]!=1){
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						fuemp[jelem][7] = fuemp[jelem][7]+f7;
						fuemp[jelem][8] = fuemp[jelem][8]+f8;
						fuemp[jelem][9] = fuemp[jelem][9]+f9;
						fuemp[jelem][10]= fuemp[jelem][10]+f10;
						fuemp[jelem][11]= fuemp[jelem][11]+f11;
						fuemp[jelem][12]= fuemp[jelem][12]+f12;
					}
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
					carga[jelem][7] = carga[jelem][7]-f7;
					carga[jelem][8] = carga[jelem][8]-f8;
					carga[jelem][9] = carga[jelem][9]-f9;
					carga[jelem][10]= carga[jelem][10]-f10;
					carga[jelem][11]= carga[jelem][11]-f11;
					carga[jelem][12]= carga[jelem][12]-f12;
				}
			}

			if(ndist > 0) {
				for(idist=1; idist<=ndist; idist++) {
					fscanf(fp5, "%d %lf", &jelem, &fuerz);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\n", jelem, fuerz);
					indfu[4][jelem]=1;
					fuerb[6][jelem]=fuerz;
                    if(ntips[jelem] ==1){xmom1 =xmom2 = 0;}
                    if(ntips[jelem] ==2){
					    xmom1 = fuerz*xlong[jelem]*xlong[jelem]/12.0;
					    xmom2 =-xmom1;
					}
                    if(ntips[jelem] ==3){
  	    				xmom1 = fuerz*xlong[jelem]*xlong[jelem]/12.0;
		    			xmom2 =-xmom1;
                        xmom1 = xmom1-xmom2/2.0;
                        xmom2 = 0;
					}
                    if(ntips[jelem] ==4){
  	    				xmom1 = fuerz*xlong[jelem]*xlong[jelem]/12.0;
		    			xmom2 =-xmom1;
                        xmom2 = -xmom1/2.0+xmom2;
                        xmom1 = 0;
					}
					reac1 = ( xmom1+xmom2-fuerz*xlong[jelem]*xlong[jelem]/2.0)/xlong[jelem];
					reac2 =-( xmom2+xmom1+fuerz*xlong[jelem]*xlong[jelem]/2.0)/xlong[jelem];
					f1 =vectr[jelem][7]*reac1;
					f2 =vectr[jelem][8]*reac1;
					f3 =vectr[jelem][9]*reac1;
 					f4 =vectr[jelem][4]*xmom1;
					f5 =vectr[jelem][5]*xmom1;
					f6 =vectr[jelem][6]*xmom1;
					f7 =vectr[jelem][7]*reac2;
					f8 =vectr[jelem][8]*reac2;
					f9 =vectr[jelem][9]*reac2;
 					f10=vectr[jelem][4]*xmom2;
					f11=vectr[jelem][5]*xmom2;
					f12=vectr[jelem][6]*xmom2;
					// PENDIENTE: verificar que este cambio es correcto
					//if(ntips[ielem]!=1){
					if(ntips[jelem]!=1){
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						fuemp[jelem][7] = fuemp[jelem][7]+f7;
						fuemp[jelem][8] = fuemp[jelem][8]+f8;
						fuemp[jelem][9] = fuemp[jelem][9]+f9;
						fuemp[jelem][10]= fuemp[jelem][10]+f10;
						fuemp[jelem][11]= fuemp[jelem][11]+f11;
						fuemp[jelem][12]= fuemp[jelem][12]+f12;
					}
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
					carga[jelem][7] = carga[jelem][7]-f7;
					carga[jelem][8] = carga[jelem][8]-f8;
					carga[jelem][9] = carga[jelem][9]-f9;
					carga[jelem][10]= carga[jelem][10]-f10;
					carga[jelem][11]= carga[jelem][11]-f11;
					carga[jelem][12]= carga[jelem][12]-f12;
				}
			}
			if(ndisty > 0) {
				for(idist=1; idist<=ndisty; idist++) {
					fscanf(fp5, "%d %lf", &jelem, &fuerz);
					if(iwrit ==1) fprintf(fp16, "%d\t%lf\n", jelem, fuerz);
					indfu[3][jelem]=1;
					fuerb[5][jelem]=fuerz;
                    if(ntips[jelem] ==1){xmom1 =xmom2 = 0;}
				    if(ntips[jelem] ==2){
					    xmom1 =-fuerz*xlong[jelem]*xlong[jelem]/12.0;
					    xmom2 =-xmom1;
					}
                    if(ntips[jelem] ==3){
  	    				xmom1 = fuerz*xlong[jelem]*xlong[jelem]/12.0;
		    			xmom2 =-xmom1;
                        xmom1 = xmom1-xmom2/2;
                        xmom2 = 0;
                    }
                    if(ntips[jelem] ==4){
  	    				xmom1 = fuerz*xlong[jelem]*xlong[jelem]/12.0;
		    			xmom2 =-xmom1;
                        xmom2 = -xmom1/2.0+xmom2;
                        xmom1 = 0;
					}
					reac1 =-( xmom1+xmom2-fuerz*xlong[jelem]*xlong[jelem]/2.0)/xlong[jelem];
					reac2 = ( xmom2+xmom1+fuerz*xlong[jelem]*xlong[jelem]/2.0)/xlong[jelem];

					f1 =vectr[jelem][4]*reac1;
					f2 =vectr[jelem][5]*reac1;
					f3 =vectr[jelem][6]*reac1;
 					f4 =vectr[jelem][7]*xmom1;
					f5 =vectr[jelem][8]*xmom1;
					f6 =vectr[jelem][9]*xmom1;
					f7 =vectr[jelem][4]*reac2;
					f8 =vectr[jelem][5]*reac2;
					f9 =vectr[jelem][6]*reac2;
 					f10=vectr[jelem][7]*xmom2;
					f11=vectr[jelem][8]*xmom2;
					f12=vectr[jelem][9]*xmom2;
					// PENDIENTE: verificar que este cambio es correcto
					//if(ntips[ielem]!=1){
					if(ntips[jelem]!=1){
						fuemp[jelem][1] = fuemp[jelem][1]+f1;
						fuemp[jelem][2] = fuemp[jelem][2]+f2;
						fuemp[jelem][3] = fuemp[jelem][3]+f3;
						fuemp[jelem][4] = fuemp[jelem][4]+f4;
						fuemp[jelem][5] = fuemp[jelem][5]+f5;
						fuemp[jelem][6] = fuemp[jelem][6]+f6;
						fuemp[jelem][7] = fuemp[jelem][7]+f7;
						fuemp[jelem][8] = fuemp[jelem][8]+f8;
						fuemp[jelem][9] = fuemp[jelem][9]+f9;
						fuemp[jelem][10]= fuemp[jelem][10]+f10;
						fuemp[jelem][11]= fuemp[jelem][11]+f11;
						fuemp[jelem][12]= fuemp[jelem][12]+f12;
					}
					carga[jelem][1] = carga[jelem][1]-f1;
					carga[jelem][2] = carga[jelem][2]-f2;
					carga[jelem][3] = carga[jelem][3]-f3;
					carga[jelem][4] = carga[jelem][4]-f4;
					carga[jelem][5] = carga[jelem][5]-f5;
					carga[jelem][6] = carga[jelem][6]-f6;
					carga[jelem][7] = carga[jelem][7]-f7;
					carga[jelem][8] = carga[jelem][8]-f8;
					carga[jelem][9] = carga[jelem][9]-f9;
					carga[jelem][10]= carga[jelem][10]-f10;
					carga[jelem][11]= carga[jelem][11]-f11;
					carga[jelem][12]= carga[jelem][12]-f12;
				}
			}
			// Si la bandera Optimiza es true no se carga el peso propio de
			// los elementos porque en el optimizador este va variando.
			if(npesp == 1 && Optimiza == false) {
				for(ielem=1; ielem<=nelem; ielem++) {
					lnod1 = lnods[ielem][1];
					lnod2 = lnods[ielem][2];
					xlon1 = coord[lnod1][1];
					ylon1 = coord[lnod1][2];
					xlon2 = coord[lnod2][1];
					ylon2 = coord[lnod2][2];
					xloxy= sqrt((xlon2-xlon1)*(xlon2-xlon1)+(ylon2-ylon1)*
								(ylon2-ylon1));
					if(xloxy >0){
						coset=(xlon2-xlon1)/xloxy;
						senot=(ylon2-ylon1)/xloxy;
					}
					else {
						coset=senot=0.0;
					}
					wtotl =-props[matnu[ielem]][2]*props[matnu[ielem]][7];
					indfu[4][ielem]=1;
					fuerb[6][ielem]=-wtotl*vectr[ielem][9];
					indfu[3][ielem]=1;
					fuerb[5][ielem]=-wtotl*vectr[ielem][8];
                    if(ntips[ielem] ==1) {xmom1 = xmom2 = 0;}
                    if(ntips[ielem] ==2){
						xmom1 =-wtotl*xlong[ielem]*xloxy/12.0;
						xmom2 =-xmom1;
					}
					if(ntips[ielem] ==3){
						xmom1 =-wtotl*xlong[ielem]*xloxy/12.0;
						xmom2 =-xmom1;
						xmom1 = xmom1-xmom2/2;
						xmom2 = 0;
                    }
					if(ntips[ielem] ==4){
						xmom1 =-wtotl*xlong[ielem]*xloxy/12.0;
						xmom2 =-xmom1;
						xmom2 = -xmom1/2.0+xmom2;
						xmom1 = 0;
                    }					if(xloxy >0){
						reac1 = ( xmom1+xmom2-wtotl*xlong[ielem]*xloxy/2.0)/xloxy;
						reac2 =-( xmom2+xmom1+wtotl*xlong[ielem]*xloxy/2.0)/xloxy;
					}
					else {reac1=reac2=-wtotl*xlong[ielem]/2.0;}

					f1 =0.0;
					f2 =0.0;
					f3 =reac1;
 					f4 =-senot*xmom1;
					f5 = coset*xmom1;
					f6 =0.0;
					f7 =0.0;
					f8 =0.0;
					f9 =reac2;
 					f10=-senot*xmom2;
					f11= coset*xmom2;
					f12=0.0;
                    if(ntips[ielem]!=1){
						fuemp[ielem][1] = fuemp[ielem][1]+f1;
						fuemp[ielem][2] = fuemp[ielem][2]+f2;
						fuemp[ielem][3] = fuemp[ielem][3]+f3;
						fuemp[ielem][4] = fuemp[ielem][4]+f4;
						fuemp[ielem][5] = fuemp[ielem][5]+f5;
						fuemp[ielem][6] = fuemp[ielem][6]+f6;
						fuemp[ielem][7] = fuemp[ielem][7]+f7;
						fuemp[ielem][8] = fuemp[ielem][8]+f8;
						fuemp[ielem][9] = fuemp[ielem][9]+f9;
						fuemp[ielem][10]= fuemp[ielem][10]+f10;
						fuemp[ielem][11]= fuemp[ielem][11]+f11;
						fuemp[ielem][12]= fuemp[ielem][12]+f12;
					}
					carga[ielem][1] = carga[ielem][1]-f1;
					carga[ielem][2] = carga[ielem][2]-f2;
					carga[ielem][3] = carga[ielem][3]-f3;
					carga[ielem][4] = carga[ielem][4]-f4;
					carga[ielem][5] = carga[ielem][5]-f5;
					carga[ielem][6] = carga[ielem][6]-f6;
					carga[ielem][7] = carga[ielem][7]-f7;
					carga[ielem][8] = carga[ielem][8]-f8;
					carga[ielem][9] = carga[ielem][9]-f9;
					carga[ielem][10]= carga[ielem][10]-f10;
					carga[ielem][11]= carga[ielem][11]-f11;
					carga[ielem][12]= carga[ielem][12]-f12;
				}
			}
		}
	}


/*
		printf("CARGAS \n");
	 for(ielem=1; ielem<=nelem; ielem++) {
	 for(ievab=1; ievab<=nevab; ievab++) {
	 printf("elemento=%d  GDL=%d FEP=%e Carga=%e\n",ielem,ievab,fuemp[ielem][ievab],carga[ielem][ievab]);
	 }
	 } */
    //printf("termino Fuerzas \n");
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Fuerzas_Puntuales(int nnodf, int nelem, int nnode, int ngdln, int iwrit, int** lnods, double** carga)
/******************************************************************************
 Ensambla fuerzas puntuales
 ******************************************************* **********************/
{
	int ipnod, ielem, inode, nloca, lodpt, igdln, ncont;
 	double pnodt[16];

	for(ipnod=1; ipnod<=nnodf; ipnod++) {
		fscanf(fp5, "%d", &lodpt);
		for(igdln=1; igdln<=ngdln; igdln++) fscanf(fp5, "%lf", &pnodt[igdln]);
		if(iwrit ==1) {
		   fprintf(fp16, "%d\t", lodpt);
		   for(igdln=1; igdln<=ngdln; igdln++) fprintf(fp16, "%lf\t", pnodt[igdln]);
		   fprintf(fp16, "\n");
		}
		// Asocia las cargas puntuales nodales con un elemento:
		for(ielem=1; ielem<=nelem; ielem++) {
			for(inode=1; inode<=nnode; inode++) {
				nloca = lnods[ielem][inode];
				if(lodpt == nloca) break;
			}

			if(lodpt == nloca) break;
		}

		for(igdln=1; igdln<=ngdln; igdln++) {
			ncont = (inode-1)*ngdln+igdln;
			carga[ielem][ncont] = carga[ielem][ncont]+pnodt[igdln];
		}
	}
}


//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Fuefix(int nelem, int ndime, int ngdln, int nevab, int isale, int nincl, int** lnods,int* lincl,
			double** girom, double** girtm, double** tempr, double** rigid, double* xincl,int* nodea,
			int* iffix, double* fixed, double*** srmat, double** carga)
/******************************************************************************
 Calcula fuerzas equivalentes para movimientos prescritos.
 ******************************************************************************/
{
	int   nn1,ne1,igdln,ielem,inode,nloca,ncont,ievab,jevab;
	double deslo[120],tenlo[120],sum;


	for(ielem=1; ielem <= nelem; ielem++) {
		nn1=nodea[ielem] ;
		ne1=nn1*ngdln;
		ncont=0;


		if(nincl >0)ApoyoInclinado(ielem,ndime,ngdln,nevab,isale,nincl,lnods,lincl,girom,girtm,tempr,
								   rigid,xincl);


		for(inode=1; inode <= nn1; inode++) {
			nloca=(lnods[ielem][inode]-1)*ngdln;
			for(igdln=1; igdln <= ngdln; igdln++) {
				if(iffix[nloca+igdln] == 1) ncont++;
			}
		}

		if(ncont > 0) {
			ncont=0;

			for(inode=1; inode <= nn1; inode++) {
				nloca=(lnods[ielem][inode]-1)*ngdln;
				for(igdln=1; igdln <= ngdln; igdln++)
					deslo[++ncont]=fixed[++nloca];
			}
			for (ievab = 1 ; ievab <= ne1 ; ievab++) {
				sum=0.0;

				for (jevab = 1 ; jevab <= ne1 ; jevab++)
					sum+=srmat[ielem][ievab][jevab]*deslo[jevab];
				tenlo[ievab]=sum;
			}

			ncont=0;

			for(inode=1; inode <= nn1; inode++) {
				for(igdln=1; igdln <= ngdln; igdln++){
					ncont++;
					carga[ielem][ncont]-=tenlo[ncont];
				}
			}
		}
	}
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑

// Copia el vector de fuerza (de un caso de carga dado) que se empleara para resolver el sistema.
// Se hace una copia porque despues lo empleara el optimizador.
void copia_vector_fuerza( int caso_i, int nelem, int nevab, double **carga, double ***VF ){
	for( int i = 1; i <= nelem; i++ ){
		for( int j = 1; j <= nevab; j++ ){
			VF[caso_i][i][j] = carga[i][j];
		}
	}
}

// Escoge el vector de fuerzas del caso de carga 
void escoge_vector_fuerza( int caso_i, int nelem, int nevab, double **carga, double ***VF ){
	for( int i = 1; i <= nelem; i++ ){
		for( int j = 1; j <= nevab; j++ ){
			carga[i][j] = VF[caso_i][i][j];
            //printf("carga[%d][%d]= %lf\n",i,j,carga[i][j]);
		}
	}	
}

// Determina el vector de fueras de empotramiento perfecto
void vector_fuer_emp_perf( int nelem, int nevab, double **fuemp, double **carga ){
	for( int i = 1; i <= nelem; i++ ){
		for( int j = 1; j <= nevab; j++ ){
			fuemp[i][j] = -carga[i][j];
		}
	}
}

// Incluye las fuerzas de peso propio al vector de fuerzas. Se usa cuando con el optmizador
void incluye_peso_propio( int npesp, int nelem, int ndime, int ntipo, int **lnods, int *matnu,
                          int *ntips, double **coord, double *xlong, double **props, double **vectr,
                          double **carga ){

	int ielem, lnod1, lnod2;
	double xlon1, ylon1, xlon2, ylon2, coset, senot, xmom1, xmom2, reac2, reac1;
	double f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, cargw, wtotl, xloxy;

	// Cargas en elementos:
	if(ndime == 2) {
		if(ntipo == 1) {
			if(npesp == 1){
				for(ielem=1; ielem<=nelem; ielem++) {
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					//coset = vectr[ielem][1];
					//cargw= 1.0;
               		//if (coset < 0.0) cargw= -1;
					//indfu[4][ielem]=1;
					//fuerb[6][ielem]+=-cargw*wtotl;
					reac1 = wtotl*xlong[ielem]/2.0;
					reac2 = reac1;
					f2 = reac1;
					f4 = reac2;
					//fuemp[ielem][2] = fuemp[ielem][2]+f2;
					//fuemp[ielem][4] = fuemp[ielem][4]+f4;
					carga[ielem][2] = carga[ielem][2]-f2;
					carga[ielem][4] = carga[ielem][4]-f4;
				}
			}
		}

		if(ntipo == 2) {
			if(npesp == 1) {
				for(ielem=1; ielem<=nelem; ielem++) {
					coset = vectr[ielem][1];
					cargw= 1.0;
               		if (coset < 0.0) cargw= -1;
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					//indfu[4][ielem]=1;
					//fuerb[6][ielem]+=-cargw*wtotl*coset;
					xmom1 = wtotl*xlong[ielem]*xlong[ielem]*fabs(coset)*cargw/12.0;
					xmom2 =-xmom1;
					reac1 = wtotl*xlong[ielem]/2.0;
					reac2 = reac1;
					f2 = reac1;
					f3 = xmom1;
					f5 = reac2;
					f6 = xmom2;
					//fuemp[ielem][2] = fuemp[ielem][2]+f2;
					//fuemp[ielem][3] = fuemp[ielem][3]+f3;
					//fuemp[ielem][5] = fuemp[ielem][5]+f5;
					//fuemp[ielem][6] = fuemp[ielem][6]+f6;
					carga[ielem][2] = carga[ielem][2]-f2;
					carga[ielem][3] = carga[ielem][3]-f3;
					carga[ielem][5] = carga[ielem][5]-f5;
					carga[ielem][6] = carga[ielem][6]-f6;
				}
			}
		}

		if(ntipo == 3) {
			if(npesp == 1) {
				for(ielem=1; ielem<=nelem; ielem++) {
					coset = vectr[ielem][1];
               		cargw= 1.0;
               		if(coset<0.0) cargw= -1;
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					//indfu[4][ielem]=1;
					//fuerb[6][ielem]+=-cargw*wtotl*coset;
					if(ntips[ielem] == 1) {
						reac1 = (wtotl*xlong[ielem]/2.0);
						reac2 = reac1;
						f2 = reac1;
						f5 = reac2;
						//fuemp[ielem][2] = fuemp[ielem][2]+f2;
						//fuemp[ielem][5] = fuemp[ielem][5]+f5;
						carga[ielem][2] = carga[ielem][2]-f2;
						carga[ielem][5] = carga[ielem][5]-f5;
					}

					if(ntips[ielem] == 2) {
						xmom1 = wtotl*xlong[ielem]*xlong[ielem]*fabs(coset)*cargw/12.0;
						xmom2 =-xmom1;
						reac1 = wtotl*xlong[ielem]/2.0;
						reac2 = reac1;
						f2 = reac1;
						f3 = xmom1;
						f5 = reac2;
						f6 = xmom2;
						//fuemp[ielem][2] = fuemp[ielem][2]+f2;
						//fuemp[ielem][3] = fuemp[ielem][3]+f3;
						//fuemp[ielem][5] = fuemp[ielem][5]+f5;
						//fuemp[ielem][6] = fuemp[ielem][6]+f6;
						carga[ielem][2] = carga[ielem][2]-f2;
						carga[ielem][3] = carga[ielem][3]-f3;
						carga[ielem][5] = carga[ielem][5]-f5;
						carga[ielem][6] = carga[ielem][6]-f6;
					}

					if(ntips[ielem] == 3) {
						xmom1 = wtotl*xlong[ielem]*xlong[ielem]*fabs(coset)*cargw/8.0;
						reac1 = 5.0*wtotl*xlong[ielem]/8.0;
						reac2 = 3.0*wtotl*xlong[ielem]/8.0;
						f2 = reac1;
						f3 = xmom1;
						f5 = reac2;
						//fuemp[ielem][2] = fuemp[ielem][2]+f2;
						//fuemp[ielem][3] = fuemp[ielem][3]+f3;
						//fuemp[ielem][5] = fuemp[ielem][5]+f5;
						carga[ielem][2] = carga[ielem][2]-f2;
						carga[ielem][3] = carga[ielem][3]-f3;
						carga[ielem][5] = carga[ielem][5]-f5;
					}
					if(ntips[ielem] == 4) {
						xmom2 = -wtotl*xlong[ielem]*xlong[ielem]*fabs(coset)*cargw/8.0;
						reac1 = 3.0*wtotl*xlong[ielem]/8.0;
						reac2 = 5.0*wtotl*xlong[ielem]/8.0;
						f2 = reac1;
						f4 = reac2;
						f5 = xmom1;
						//fuemp[ielem][2] = fuemp[ielem][2]+f2;
						//fuemp[ielem][5] = fuemp[ielem][5]+f5;
						//fuemp[ielem][6] = fuemp[ielem][6]+f6;
						carga[ielem][2] = carga[ielem][2]-f2;
						carga[ielem][3] = carga[ielem][3]-f4;
						carga[ielem][5] = carga[ielem][5]-f5;
					}
				}
			}
		}
	}

	if(ndime == 3) {
		if(ntipo == 1) {
			if(npesp == 1) {
				for(ielem=1; ielem<=nelem; ielem++) {
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					reac1 = (wtotl*xlong[ielem]/2.0);
					reac2 = reac1;
					//indfu[4][ielem]=1;
					//fuerb[6][ielem]=-wtotl*vectr[ielem][3];
					//indfu[3][ielem]=1;
					//fuerb[5][ielem]=0;
					f3 = reac1;
					f6 = reac2;
					//fuemp[ielem][3] = fuemp[ielem][3]+f3;
					//fuemp[ielem][6] = fuemp[ielem][6]+f6;
					carga[ielem][3] = carga[ielem][3]-f3;
					carga[ielem][6] = carga[ielem][6]-f6;
				}
			}
		}

		if(ntipo == 2) {
			if(npesp == 1) {
				for(ielem=1; ielem<=nelem; ielem++) {
					lnod1 = lnods[ielem][1];
					lnod2 = lnods[ielem][2];
					xlon1 = coord[lnod1][1];
					ylon1 = coord[lnod1][2];
					xlon2 = coord[lnod2][1];
					ylon2 = coord[lnod2][2];
					xloxy= sqrt((xlon2-xlon1)*(xlon2-xlon1)+(ylon2-ylon1)*
								(ylon2-ylon1));
					if(xloxy >0){
						coset=(xlon2-xlon1)/xloxy;
						senot=(ylon2-ylon1)/xloxy;
					}
					else {
						coset=senot=0.0;
					}
					wtotl = props[matnu[ielem]][2]*props[matnu[ielem]][7];
					//indfu[4][ielem]=1;
					//fuerb[6][ielem]=-wtotl*vectr[ielem][9];
					//indfu[3][ielem]=1;
					//fuerb[5][ielem]=-wtotl*vectr[ielem][8];
					xmom1 =-wtotl*xlong[ielem]*xloxy/12.0;
					xmom2 =-xmom1;
					reac1 = wtotl*xlong[ielem]/2.0;
					reac2 = reac1;
					f1 =0.0;
					f2 =0.0;
					f3 =reac1;
 					f4 =-senot*xmom1;
					f5 = coset*xmom1;
					f6 =0.0;
					f7 =0.0;
					f8 =0.0;
					f9 =reac2;
 					f10=-senot*xmom2;
					f11= coset*xmom2;
					f12=0.0;
					/*
					fuemp[ielem][1] = fuemp[ielem][1]+f1;
					fuemp[ielem][2] = fuemp[ielem][2]+f2;
					fuemp[ielem][3] = fuemp[ielem][3]+f3;
					fuemp[ielem][4] = fuemp[ielem][4]+f4;
					fuemp[ielem][5] = fuemp[ielem][5]+f5;
					fuemp[ielem][6] = fuemp[ielem][6]+f6;
					fuemp[ielem][7] = fuemp[ielem][7]+f7;
					fuemp[ielem][8] = fuemp[ielem][8]+f8;
					fuemp[ielem][9] = fuemp[ielem][9]+f9;
					fuemp[ielem][10]= fuemp[ielem][10]+f10;
					fuemp[ielem][11]= fuemp[ielem][11]+f11;
					fuemp[ielem][12]= fuemp[ielem][12]+f12;
					*/
					carga[ielem][1] = carga[ielem][1]-f1;
					carga[ielem][2] = carga[ielem][2]-f2;
					carga[ielem][3] = carga[ielem][3]-f3;
					carga[ielem][4] = carga[ielem][4]-f4;
					carga[ielem][5] = carga[ielem][5]-f5;
					carga[ielem][6] = carga[ielem][6]-f6;
					carga[ielem][7] = carga[ielem][7]-f7;
					carga[ielem][8] = carga[ielem][8]-f8;
					carga[ielem][9] = carga[ielem][9]-f9;
					carga[ielem][10]= carga[ielem][10]-f10;
					carga[ielem][11]= carga[ielem][11]-f11;
					carga[ielem][12]= carga[ielem][12]-f12;
				}
			}
		}

		if(ntipo == 3) {
			if(npesp == 1) {
				for(ielem=1; ielem<=nelem; ielem++) {
					lnod1 = lnods[ielem][1];
					lnod2 = lnods[ielem][2];
					xlon1 = coord[lnod1][1];
					ylon1 = coord[lnod1][2];
					xlon2 = coord[lnod2][1];
					ylon2 = coord[lnod2][2];
					xloxy= sqrt((xlon2-xlon1)*(xlon2-xlon1)+(ylon2-ylon1)*
								(ylon2-ylon1));
					if(xloxy >0){
						coset=(xlon2-xlon1)/xloxy;
						senot=(ylon2-ylon1)/xloxy;
					}
					else {
						coset=senot=0.0;
					}
					wtotl =-props[matnu[ielem]][2]*props[matnu[ielem]][7];
					//indfu[4][ielem]=1;
					//fuerb[6][ielem]=-wtotl*vectr[ielem][9];
					//indfu[3][ielem]=1;
					//fuerb[5][ielem]=-wtotl*vectr[ielem][8];
                    if(ntips[ielem] ==1) {xmom1 = xmom2 = 0;}
                    if(ntips[ielem] ==2){
						xmom1 =-wtotl*xlong[ielem]*xloxy/12.0;
						xmom2 =-xmom1;
					}
					if(ntips[ielem] ==3){
						xmom1 =-wtotl*xlong[ielem]*xloxy/12.0;
						xmom2 =-xmom1;
						xmom1 = xmom1-xmom2/2;
						xmom2 = 0;
                    }
					if(ntips[ielem] ==4){
						xmom1 =-wtotl*xlong[ielem]*xloxy/12.0;
						xmom2 =-xmom1;
						xmom2 = -xmom1/2.0+xmom2;
						xmom1 = 0;
                    }					if(xloxy >0){
						reac1 = ( xmom1+xmom2-wtotl*xlong[ielem]*xloxy/2.0)/xloxy;
						reac2 =-( xmom2+xmom1+wtotl*xlong[ielem]*xloxy/2.0)/xloxy;
					}
					else {reac1=reac2=-wtotl*xlong[ielem]/2.0;}

					f1 =0.0;
					f2 =0.0;
					f3 =reac1;
 					f4 =-senot*xmom1;
					f5 = coset*xmom1;
					f6 =0.0;
					f7 =0.0;
					f8 =0.0;
					f9 =reac2;
 					f10=-senot*xmom2;
					f11= coset*xmom2;
					f12=0.0;
					/*
                    if(ntips[ielem]!=1){
						fuemp[ielem][1] = fuemp[ielem][1]+f1;
						fuemp[ielem][2] = fuemp[ielem][2]+f2;
						fuemp[ielem][3] = fuemp[ielem][3]+f3;
						fuemp[ielem][4] = fuemp[ielem][4]+f4;
						fuemp[ielem][5] = fuemp[ielem][5]+f5;
						fuemp[ielem][6] = fuemp[ielem][6]+f6;
						fuemp[ielem][7] = fuemp[ielem][7]+f7;
						fuemp[ielem][8] = fuemp[ielem][8]+f8;
						fuemp[ielem][9] = fuemp[ielem][9]+f9;
						fuemp[ielem][10]= fuemp[ielem][10]+f10;
						fuemp[ielem][11]= fuemp[ielem][11]+f11;
						fuemp[ielem][12]= fuemp[ielem][12]+f12;
					}
					*/
					carga[ielem][1] = carga[ielem][1]-f1;
					carga[ielem][2] = carga[ielem][2]-f2;
					carga[ielem][3] = carga[ielem][3]-f3;
					carga[ielem][4] = carga[ielem][4]-f4;
					carga[ielem][5] = carga[ielem][5]-f5;
					carga[ielem][6] = carga[ielem][6]-f6;
					carga[ielem][7] = carga[ielem][7]-f7;
					carga[ielem][8] = carga[ielem][8]-f8;
					carga[ielem][9] = carga[ielem][9]-f9;
					carga[ielem][10]= carga[ielem][10]-f10;
					carga[ielem][11]= carga[ielem][11]-f11;
					carga[ielem][12]= carga[ielem][12]-f12;
				}
			}
		}
	}
 /*   for( int i = 1; i <= nelem; i++ ){
        for( int j = 1; j <= 12; j++ ){
            printf("carga[%d][%d]= %lf\n",i,j,carga[i][j]);
        }
    }*/
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
