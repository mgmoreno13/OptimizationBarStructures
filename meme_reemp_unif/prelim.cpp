//*****************************
// Rutinas para Calculos preliminares
//*****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "prelim.h"
#include "raros.h"


//************************************************void prelim
void Prelim(int nelem, int ndime, int ntipo, int** lnods, double** coord, double* xlong,
			double** vmatr, double** vectr, double* angulo)
{
	int    ielem,lnod1,lnod2,lnod3,indic;
	double xlon1,xlon2,xlon3,ylon1,ylon2,ylon3,zlon1,zlon2,zlon3,xlona,
	coset,senot,cosef,senof,valo1,valo2,valo3,beta,ca,sa,cb,sb,cg,sg;

	indic=1;

    for(ielem=1; ielem<=nelem; ielem++){
		lnod1 = lnods[ielem][1];
		lnod2 = lnods[ielem][2];
		xlon1 = coord[lnod1][1];
	    ylon1 = coord[lnod1][2];
	    xlon2 = coord[lnod2][1];
	    ylon2 = coord[lnod2][2];
	    if(ndime ==2){
			xlong[ielem]=sqrt((xlon2-xlon1)*(xlon2-xlon1)+(ylon2-ylon1)*(ylon2-ylon1));
			vectr[ielem][1] = (xlon2-xlon1)/xlong[ielem]; //coset
			vectr[ielem][2] = (ylon2-ylon1)/xlong[ielem]; //senot
		}
		else if(ndime==3){
			zlon1 = coord[lnod1][3];
			zlon2 = coord[lnod2][3];
			xlong[ielem]=sqrt((xlon2-xlon1)*(xlon2-xlon1)+(ylon2-ylon1)*(ylon2-ylon1)+(zlon2-zlon1)*(zlon2-zlon1));
			if (ntipo ==1) {
				xlona = sqrt((xlon2-xlon1)*(xlon2-xlon1)+(ylon2-ylon1)*
							 (ylon2-ylon1));
				if(xlona > 0) {
					coset = (xlon2-xlon1)/xlona;
					senot = (ylon2-ylon1)/xlona;
	   		    }
				else {
					coset = 0.0;
					senot = 0.0;
   			    }
				cosef = (zlon2-zlon1)/xlong[ielem];
				senof = 1.0-cosef*cosef;
				if(senof > 0) {
					if(senof > 1) senof = 1.0;
	   			    senof = sqrt(senof);
				}
			    else senof = 0.0;
				vectr[ielem][1]=coset;
				vectr[ielem][2]=senot;
				vectr[ielem][3]=cosef;
				vectr[ielem][4]=senof;
			}
			else {
				if(indic ==1 && lnods[ielem][3]!=0){
					lnod3=lnods[ielem][3];
					xlon3=coord[lnod3][1];
					ylon3=coord[lnod3][2];
					zlon3=coord[lnod3][3];
					//V1
					vmatr[1][1]=xlon2-xlon1;
					vmatr[2][1]=ylon2-ylon1;
					vmatr[3][1]=zlon2-zlon1;
					//V2
					vmatr[1][2]=xlon3-xlon1;
					vmatr[2][2]=ylon3-ylon1;
					vmatr[3][2]=zlon3-zlon1;
					// Producto V3=V1xV2
					vmatr[1][3]=vmatr[2][1]*vmatr[3][2]-vmatr[2][2]*vmatr[3][1];
					vmatr[2][3]=vmatr[1][2]*vmatr[3][1]-vmatr[1][1]*vmatr[3][2];
					vmatr[3][3]=vmatr[1][1]*vmatr[2][2]-vmatr[1][2]*vmatr[2][1];
					// Producto V2=V3*V1
					vmatr[1][2]=vmatr[2][3]*vmatr[3][1]-vmatr[2][1]*vmatr[3][3];
					vmatr[2][2]=vmatr[1][1]*vmatr[3][3]-vmatr[1][3]*vmatr[3][1];
					vmatr[3][2]=vmatr[1][3]*vmatr[2][1]-vmatr[1][1]*vmatr[2][3];
					// Norma de V1, V2 y V3
					valo1=sqrt(vmatr[1][1]*vmatr[1][1]+vmatr[2][1]*vmatr[2][1]+vmatr[3][1]
							   *vmatr[3][1]);
					valo2=sqrt(vmatr[1][2]*vmatr[1][2]+vmatr[2][2]*vmatr[2][2]+vmatr[3][2]
							   *vmatr[3][2]);
					valo3=sqrt(vmatr[1][3]*vmatr[1][3]+vmatr[2][3]*vmatr[2][3]+vmatr[3][3]
							   *vmatr[3][3]);
					// Normaliza V1
					vmatr[1][1]=vmatr[1][1]/valo1;
					vmatr[2][1]=vmatr[2][1]/valo1;
					vmatr[3][1]=vmatr[3][1]/valo1;
					// Normaliza V2
					vmatr[1][2]=vmatr[1][2]/valo2;
					vmatr[2][2]=vmatr[2][2]/valo2;
					vmatr[3][2]=vmatr[3][2]/valo2;
					// Normaliza V3
					vmatr[1][3]=vmatr[1][3]/valo3;
					vmatr[2][3]=vmatr[2][3]/valo3;
					vmatr[3][3]=vmatr[3][3]/valo3;
				}
				else if(indic ==1 && lnods[ielem][3]==0){
					beta=angulo[ielem]*0.017453292519943;
					xlona = sqrt((xlon2-xlon1)*(xlon2-xlon1)+(ylon2-ylon1)*(ylon2-ylon1));
					if(xlona > 0) {
						coset = (xlon2-xlon1)/xlona;
						senot = (ylon2-ylon1)/xlona;
					}
					else {
						coset = 1.0;
						senot = 0.0;
					}

                    senof = (zlon2-zlon1)/xlong[ielem];
					cosef = 1.0-senof*senof;

					if(cosef > 0) {
						if(cosef > 1) cosef = 1.0;
						cosef = sqrt(cosef);
					}
					else cosef = 0.0;
					ca=coset;
					sa=senot;
					cg=cosef;
					sg=senof;
					cb=cos(beta);
					sb=sin(beta);
					vmatr[1][1]= cg*ca;
					vmatr[2][1]= cg*sa;
					vmatr[3][1]= sg;
					vmatr[1][2]=-sb*sg*ca-cb*sa;		//
					vmatr[2][2]=-sb*sg*sa+cb*ca;
					vmatr[3][2]= cg*sb;
                    vmatr[1][3]=-cb*sg*ca+sb*sa;		//
					vmatr[2][3]=-cb*sg*sa-sb*ca;		//
					vmatr[3][3]= cb*cg;

				}
				else{
					xlona = sqrt((xlon2-xlon1)*(xlon2-xlon1)+(ylon2-ylon1)*(ylon2-ylon1));
					if(xlona > 0) {
						coset = (xlon2-xlon1)/xlona;
						senot = (ylon2-ylon1)/xlona;
					}
					else {
						coset = 0.0;
						senot = 0.0;
					}
					cosef = (zlon2-zlon1)/xlong[ielem];
					senof = 1.0-cosef*cosef;
					if(senof > 0) {
						if(senof > 1) senof = 1.0;
						senof = sqrt(senof);
					}
					else senof = 0.0;
					vmatr[1][1]= coset*senof;
					vmatr[1][2]=-senot*senof;
					vmatr[1][3]=-cosef;
					vmatr[2][1]= senot*senof;
					vmatr[2][2]= coset*senof;
					vmatr[2][3]=-cosef;
					vmatr[3][1]= cosef;
					vmatr[3][2]= cosef;
					vmatr[3][3]= senof;
				}
				vectr[ielem][1]=vmatr[1][1];
				vectr[ielem][2]=vmatr[2][1];
				vectr[ielem][3]=vmatr[3][1];
				vectr[ielem][4]=vmatr[1][2];
				vectr[ielem][5]=vmatr[2][2];
				vectr[ielem][6]=vmatr[3][2];
				vectr[ielem][7]=vmatr[1][3];
				vectr[ielem][8]=vmatr[2][3];
				vectr[ielem][9]=vmatr[3][3];
			}
		}
	}
}
