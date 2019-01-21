//*****************************
// Rutina para posproceso
//*****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "pospro.h"
#include "raros.h"


//************************************************
void Pospro(int npnod, int nelem, int npres, int ntipo, int nnode, int ngdln, int nmats,
			int ngaus, int ndime, int ntens, int nevab, int ntotv, int** lnods, int* matnu,
			int* iffix, double** coord, double* aslod, double* despl, double** carga)
/******************************************************************************
 Posproceso.
 ******************************************************************************/
{

	int ipnod, ielem, itotv, ktotv, jelem, igdln, inode, ievab, ncont, ngish,
	igash;

	fprintf(fp9, "%d\t", 1);
	fprintf(fp9, "\n") ;
	fprintf(fp9, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", npnod,
			nelem, npres, 1, ntipo, nnode, ngdln, nmats, ngaus, ndime, ntens,
			nevab);
	jelem = 0;

	for(ielem=1; ielem<=nelem; ielem++) {
		jelem++;
		fprintf(fp9, "%d\t%d\t%d\t%d\t%d\t%d\t%d\n", jelem,	lnods[ielem][1],
				lnods[ielem][2], 0, 0, 0, matnu[ielem]);
	}

	if(ndime == 2)
		for(ipnod=1; ipnod<=npnod; ipnod++)
			fprintf(fp9, "%d\t%lf\t%lf\n", ipnod, coord[ipnod][1], coord[ipnod][2]);

	if(ndime == 3)
		for(ipnod=1; ipnod<=npnod; ipnod++)
			fprintf(fp9, "%d\t%lf\t%lf\n%lf\n",
					ipnod, coord[ipnod][1], coord[ipnod][2], coord[ipnod][3]);

	fprintf(fp9, "%d\n", ntotv);

	for(ipnod=1; ipnod<=npnod; ipnod++) {
		ktotv = (ipnod-1)*ngdln+1;

		for(igdln=1; igdln<=ngdln; igdln++)
			if(iffix[ktotv++] == 0) fprintf(fp9, "%d\n", 1);
			else fprintf(fp9, "%d\t", 0);
	}

	fprintf(fp9, "\n");
	fprintf(fp9, "%d\n", ntotv);

	for(itotv=1; itotv<=ntotv; itotv++)
		aslod[itotv] = 0.0;

	for(ielem=1; ielem<=nelem; ielem++) {
		for(inode=1; inode<=nnode; inode++) {
			ktotv = (lnods[ielem][inode]-1)*ngdln;
			ievab = (inode-1)*ngdln;

			for(igdln=1; igdln<=ngdln; igdln++)
				aslod[ktotv+igdln] = aslod[ktotv+igdln]+
				carga[ielem][ievab+igdln];
		}
	}

	for(ipnod=1; ipnod<=npnod; ipnod++) {
		ktotv = (ipnod-1)*ngdln+1;

		for(igdln=1; igdln<=ngdln; igdln++)
			fprintf(fp9, "%lf\t", aslod[ktotv++]);
	}

	fprintf(fp9, "\n");

	// C√°lculo de desplazamientos:
	for(ipnod=1; ipnod<=npnod; ipnod++) {
		ncont = ipnod*ngdln;
		ngish = ncont-ngdln+1;
		fprintf(fp9, "%d\t", ipnod);

		for(igash=ngish; igash<=ncont; igash++)
			fprintf(fp9, "%lf\t", despl[igash]);

     	fprintf(fp9, "\n");
	}

	// Calculo de tensiones en los nodos:
	fprintf(fp9, "%d\n", 0);
	fclose(fp9);
}
