//*****************************
// Rutina para posproceso
//*****************************
#pragma hdrstop
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "posprogid.h"
#include "raros.h"


//************************************************
void PosproGid(int npnod, int nelem, int ndime, int nnode, int ngdln, int ntotv, int nmats,int ncaso,
			   int paso, int** lnods, int* matnu, int* matva, double** coord, double **props,
			   double* despl, double* fuerc, double* aslod, double** carpp, char* titulo, double** mvector)
/******************************************************************************
 Posproceso.
 ******************************************************************************/
{

	int ipnod, ielem, itotv, ktotv, igdln, inode, ievab, ncont, ngish,
	    igash, idime, imats, ngaus=1, ivpro;
 

    if(ncaso ==1){
    if(paso ==1){
    //Archivo de postproceso para GiD malla
	fpresults = fopen(nomFlaviaMesh, "wt");

	//archivo de posproceso para la escritura de la geometria de la Malla
	fprintf(fpresults, "#  MEFI 2.0\n#\tElaborado por Salvador Botello\n");
	fprintf(fpresults, "#\t Centro de Investigacion en Matematicas\n");
	fprintf(fpresults, "#\t\tGuanajuato, Mexico\t\t   20102n\n");
	fprintf(fpresults, "MESH \"%s\"  ", titulo);
	fprintf(fpresults, "dimension  %d  ", ndime);
	fprintf(fpresults, "ElemType Linear  Nnode    2 \n");

			
	fprintf(fpresults, "Coordinates\n");
	fprintf(fpresults, "#node number\tcoord-X\tcoord-Y");
	if(ndime==3)
		fprintf(fpresults, "\tcoord-Z\n");
	else
		fprintf(fpresults, "\n");
	
	for(ipnod=1; ipnod<=npnod; ipnod++) {
		fprintf(fpresults, "\t%d\t", ipnod);
		for(idime=1; idime<=ndime; idime++)
			fprintf(fpresults, "%e\t", coord[ipnod][idime]);
		fprintf(fpresults, "\n");
	}
	fprintf(fpresults, "end coordinates\n");

	fprintf(fpresults, "Elements\n");
	fprintf(fpresults, "#element  ");
	for(inode=1; inode<=2; inode++)
		fprintf(fpresults, "node_%d  ", inode);
	fprintf(fpresults, "material\n");
	
	for(ielem=1; ielem<=nelem; ielem++) {
		fprintf(fpresults, "\t%d\t",ielem);
		for(inode=1; inode<=2; inode++)
			fprintf(fpresults, "%d\t", lnods[ielem][inode]);
		fprintf(fpresults, "%d\t",matnu[ielem]);
		fprintf(fpresults, "\n");
	}
	fprintf(fpresults, "end elements\n");
	
	fclose(fpresults);

	}
    //Archivo de postproceso para GiD resultados

    if(paso ==1){
	fpresults = fopen(nomFlaviaRes, "wt");
	fprintf(fpresults, "GiD  Post  Results File  1.0\n\n");
    fprintf(fpresults, "GaussPoints  \"Board Internal\"  ElemType Linear \n ");	
	fprintf(fpresults, "Number Of Gauss Points: %d\n", ngaus);
	fprintf(fpresults, "Natural Coordinates: internal\nend gausspoints\n\n");
	}
	fprintf(fpresults, "Result \"Displacements\" \"Load Analysis\" %d Vector OnNodes\n", paso);
	fprintf(fpresults, "ComponentNames \"X-DISP\", \"Y-DISP\", \"Z-DISP\"\n");
	fprintf(fpresults, "Values\n") ;
	for(ipnod=1; ipnod<=npnod; ipnod++) {
		ncont = ipnod*ngdln;
		ngish = ncont-ngdln+1;
		if(ndime==2) ncont=ngish+1;
		if(ndime==3) ncont=ngish+2;
		fprintf(fpresults, "    %d \t", ipnod);
		for(igash=ngish; igash<= ncont; igash++)
			fprintf(fpresults, "%8e\t", despl[igash]);
		if(ndime==2) fprintf(fpresults, "0.000000+e00");
		fprintf(fpresults, "\n");
	}
	fprintf(fpresults,"End Values\n\n");

	for(imats=1; imats<=nmats; imats++) {
		if (ndime ==2) matva[imats]= (int)props[imats][6];
		if (ndime ==3) matva[imats]= (int)props[imats][9];	
	}		
	fprintf(fpresults, "Result \"Eficiencias\" \"Load Analysis\" %d Scalar OnGaussPoints   \"Board Internal\" \n", paso);
	fprintf(fpresults, "ComponentNames \"Efic\"\n");
	fprintf(fpresults, "Values\n") ;
	for(ielem=1; ielem<=nelem; ielem++) 
			fprintf(fpresults, "%d  \t  %8e\n", ielem, fabs(fuerc[ielem]));
	fprintf(fpresults,"End Values\n\n");
	
	fprintf(fpresults, "Result \"Material\" \"Load Analysis\" %d Scalar OnGaussPoints   \"Board Internal\" \n", paso);
	fprintf(fpresults, "ComponentNames \"Mater\"\n");
	fprintf(fpresults, "Values\n") ;
	for(ielem=1; ielem<=nelem; ielem++){ 
		imats=matnu[ielem];
		fprintf(fpresults, "%d  \t  %8e\n", ielem, (float)matva[imats]);
	}
	fprintf(fpresults,"End Values\n\n");
	
	for(itotv=1; itotv<=ntotv; itotv++)
		aslod[itotv] = 0.0;
	
	for(ielem=1; ielem<=nelem; ielem++) {
		for(inode=1; inode<=nnode; inode++) {
			ktotv = (lnods[ielem][inode]-1)*ngdln;
			ievab = (inode-1)*ngdln;			
			for(igdln=1; igdln<=ngdln; igdln++)
				aslod[ktotv+igdln] = aslod[ktotv+igdln]+
				carpp[ielem][ievab+igdln];
		}
	}
	fprintf(fpresults, "Result \"Loads\" \"Load Analysis\"  %d Vector OnNodes\n", paso);
	fprintf(fpresults, "ComponentNames \"X-Load\", \"Y-Load\", \"Z-Load\"\n");
	fprintf(fpresults, "Values\n") ;
	for(ipnod=1; ipnod<=npnod; ipnod++) {
		ncont = ipnod*ngdln;
		ngish = ncont-ngdln+1;
		if(ndime==2) ncont=ngish+1;
		if(ndime==3) ncont=ngish+2;
		fprintf(fpresults, "    %d \t", ipnod);
		for(igash=ngish; igash<= ncont; igash++)
			fprintf(fpresults, "%8e\t", aslod[igash]);
		if(ndime==2) fprintf(fpresults, "0.000000+e00");
		fprintf(fpresults, "\n");
	}
	fprintf(fpresults,"End Values\n\n");
	
	fprintf(fpresults, "Result \"Moments\" \"Load Analysis\"  %d Vector OnNodes\n", paso);
	fprintf(fpresults, "ComponentNames \"Mx-Load\", \"My-Load\", \"Mz-Load\"\n");
	fprintf(fpresults, "Values\n") ;
	for(ipnod=1; ipnod<=npnod; ipnod++) {
		ncont = ipnod*ngdln;
		ngish = ncont-ngdln+4;
		if(ndime==2) ngish=ncont;
		if(ndime==3) ncont=ngish+2;
		fprintf(fpresults, "    %d \t", ipnod);
		for(igash=ngish; igash<= ncont; igash++)
			fprintf(fpresults, "%8e\t", aslod[igash]);
		if(ndime==2) fprintf(fpresults, "0.000000+e00\t 0.000000+e00");
		fprintf(fpresults, "\n");
	}
	fprintf(fpresults,"End Values\n\n");
	for(itotv=1; itotv<=ntotv; itotv++)
		aslod[itotv] = despl[itotv]= 0.0;
    }
    if(ncaso ==2){

        //Archivo de postproceso para GiD malla
        fpresults = fopen(nomFlaviaMesh, "wt");
        
        //archivo de posproceso para la escritura de la geometria de la Malla
        fprintf(fpresults, "#  MEFI 2.0\n#\tElaborado por Salvador Botello\n");
        fprintf(fpresults, "#\t Centro de Investigacion en Matematicas\n");
        fprintf(fpresults, "#\t\tGuanajuato, Mexico\t\t   20102n\n");
        fprintf(fpresults, "MESH \"%s\"  ", titulo);
        fprintf(fpresults, "dimension  %d  ", ndime);
        fprintf(fpresults, "ElemType Linear  Nnode    2 \n");
        
        
        fprintf(fpresults, "Coordinates\n");
        fprintf(fpresults, "#node number\tcoord-X\tcoord-Y");
        if(ndime==3)
            fprintf(fpresults, "\tcoord-Z\n");
        else
            fprintf(fpresults, "\n");
        
        for(ipnod=1; ipnod<=npnod; ipnod++) {
            fprintf(fpresults, "\t%d\t", ipnod);
            for(idime=1; idime<=ndime; idime++)
                fprintf(fpresults, "%e\t", coord[ipnod][idime]);
            fprintf(fpresults, "\n");
        }
        fprintf(fpresults, "end coordinates\n");
        
        fprintf(fpresults, "Elements\n");
        fprintf(fpresults, "#element  ");
        for(inode=1; inode<=2; inode++)
            fprintf(fpresults, "node_%d  ", inode);
        fprintf(fpresults, "material\n");
        
        for(ielem=1; ielem<=nelem; ielem++) {
            fprintf(fpresults, "\t%d\t",ielem);
            for(inode=1; inode<=2; inode++)
                fprintf(fpresults, "%d\t", lnods[ielem][inode]);
            fprintf(fpresults, "%d\t",matnu[ielem]);
            fprintf(fpresults, "\n");
        }
        fprintf(fpresults, "end elements\n");
        
        fclose(fpresults);

        fpresults = fopen(nomFlaviaRes, "wt");
        fprintf(fpresults, "GiD  Post  Results File  1.0\n\n");
        fprintf(fpresults, "GaussPoints  \"Board Internal\"  ElemType Linear \n ");
        fprintf(fpresults, "Number Of Gauss Points: %d\n", ngaus);
        fprintf(fpresults, "Natural Coordinates: internal\nend gausspoints\n\n");
        for(ivpro=1; ivpro<=paso; ivpro++){
            fprintf(fpresults, "Result \"Displacements\" \"Load Analysis\" %d Vector OnNodes\n", ivpro);
        fprintf(fpresults, "ComponentNames \"X-DISP\", \"Y-DISP\", \"Z-DISP\"\n");
        fprintf(fpresults, "Values\n") ;
        for(ipnod=1; ipnod<=npnod; ipnod++) {
            ncont = ipnod*ngdln;
            ngish = ncont-ngdln+1;
            if(ndime==2) ncont=ngish+1;
            if(ndime==3) ncont=ngish+2;
            fprintf(fpresults, "    %d \t", ipnod);
            for(igash=ngish; igash<= ncont; igash++)
                fprintf(fpresults, "%8e\t", mvector[ivpro][igash]);
            if(ndime==2) fprintf(fpresults, "0.000000+e00");
            fprintf(fpresults, "\n");
        }
        fprintf(fpresults,"End Values\n\n");
        }
    }
	return;
}
