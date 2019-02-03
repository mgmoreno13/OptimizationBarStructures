//*****************************
// Rutinas para Calculo d ematrices de Rigidez y Masa
//*****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "main.h"
#include "raros.h"
#include "rigidez.h"



void Rigimat(int nelem, int nevab, int ndime, int ntipo, int isale, int* matnu, int* ntips, double*** srmat,
			 double** props, double* xlong, double** vectr, double** rigid, double** girom,
			 double** girtm, double** tempr)
/******************************************************************************
 Calcula matriz de Rigidez de cada elemento.
 ******************************************************************************/
{
	int ielem, lprop, ievab, jevab;
	double 	young, xarea, xiner, yiner, xjner, gcort, consta, consa, consi, consl, con4y,
			con6y, co12y, con4z, con6z, co12z, conjx, coset, senot, cosef, senof;

	for(ielem=1; ielem<=nelem; ielem++)
		for(ievab=1; ievab<=nevab; ievab++)
			for(jevab=1; jevab<=nevab; jevab++)
				srmat[ielem][ievab][jevab]=0;

    // Ciclo sobre cada elemento:
	for(ielem=1; ielem<=nelem; ielem++) {

		lprop = abs(matnu[ielem]);
		young = props[lprop][1];
		xarea = props[lprop][2];
		xiner = props[lprop][3];
		//if(ndime ==3)
        {
			yiner = props[lprop][4];
			xjner = props[lprop][5];
			gcort = props[lprop][6];
		}
		// Inicializa la matriz de rigidez elemental:
		for(ievab=1; ievab<=nevab; ievab++)
			for(jevab=1; jevab<=nevab; jevab++)
				rigid[ievab][jevab] = 0.0;

		if(ndime == 2) {
			// Para 2D:
			coset = vectr[ielem][1];
			senot = vectr[ielem][2];
			// Calculo de la matriz de rigidez:
			consa = young*xarea/xlong[ielem];

			if(ntipo == 1) {
				rigid[1][1] = consa*coset*coset;
				rigid[1][2] = consa*coset*senot;
				rigid[1][3] =-rigid[1][1];
				rigid[1][4] =-rigid[1][2];
				rigid[2][1] = consa*senot*coset;
				rigid[2][2] = consa*senot*senot;
				rigid[2][3] =-rigid[2][1];
				rigid[2][4] =-rigid[2][2];
				rigid[3][1] =-rigid[1][1];
				rigid[3][2] =-rigid[1][2];
				rigid[3][3] = rigid[1][1];
				rigid[3][4] = rigid[1][2];
				rigid[4][1] =-rigid[2][1];
				rigid[4][2] =-rigid[2][2];
				rigid[4][4] = rigid[2][2];
				rigid[4][3] = rigid[2][1];
			}

			if(ntipo == 2) {
				consi = 4.0*young*xiner/xlong[ielem];
				consta = 1.5*consi/xlong[ielem];
				consl = 2.0*consta/xlong[ielem];
				rigid[1][1] = consa*coset*coset+consl*senot*senot;
				rigid[1][2] = (consa-consl)*senot*coset;
				rigid[1][3] =-consta*senot;
				rigid[1][4] =-rigid[1][1];
				rigid[1][5] =-rigid[1][2];
				rigid[1][6] = rigid[1][3];
				rigid[2][2] = consa*senot*senot+consl*coset*coset;
				rigid[2][3] = consta*coset;
				rigid[2][4] =-rigid[1][2];
				rigid[2][5] =-rigid[2][2];
				rigid[2][6] = rigid[2][3];
				rigid[3][3] = consi;
				rigid[3][4] =-rigid[1][3];
				rigid[3][5] =-rigid[2][3];
				rigid[3][6] = consi/2.0;
				rigid[4][4] = consa*coset*coset+consl*senot*senot;
				rigid[4][5] = rigid[1][2];
				rigid[4][6] =-rigid[1][3];
				rigid[5][5] = consa*senot*senot+consl*coset*coset;
				rigid[5][6] =-rigid[2][3];
				rigid[6][6] = consi;

				for(ievab=1; ievab<=nevab; ievab++) {
					for(jevab=ievab; jevab<=nevab; jevab++)
						rigid[jevab][ievab] = rigid[ievab][jevab];
				}
			}

			if(ntipo == 3) {
				if(ntips[ielem] == 1) {
					rigid[1][1] = consa*coset*coset;
					rigid[1][2] = consa*coset*senot;
					rigid[1][4] =-rigid[1][1];
					rigid[1][5] =-rigid[1][2];
					rigid[2][1] = consa*senot*coset;
					rigid[2][2] = consa*senot*senot;
					rigid[2][4] =-rigid[2][1];
					rigid[2][5] =-rigid[2][2];
					rigid[4][1] =-rigid[1][1];
					rigid[4][2] =-rigid[1][2];
					rigid[4][4] = rigid[1][1];
					rigid[4][5] = rigid[1][2];
					rigid[5][1] =-rigid[2][1];
					rigid[5][2] =-rigid[2][2];
					rigid[5][4] = rigid[2][1];
					rigid[5][5] = rigid[2][2];
				}

				if(ntips[ielem] == 2) {
					consi = 4.0*young*xiner/xlong[ielem];
					consta = 1.5*consi/xlong[ielem];
					consl = 2.0*consta/xlong[ielem];
					rigid[1][1] = consa*coset*coset+consl*senot*senot;
					rigid[1][2] = (consa-consl)*senot*coset;
					rigid[1][3] =-consta*senot;
					rigid[1][4] =-rigid[1][1];
					rigid[1][5] =-rigid[1][2];
					rigid[1][6] = rigid[1][3];
					rigid[2][2] = consa*senot*senot+consl*coset*coset;
					rigid[2][3] = consta*coset;
					rigid[2][4] =-rigid[1][2];
					rigid[2][5] =-rigid[2][2];
					rigid[2][6] = rigid[2][3];
					rigid[3][3] = consi;
					rigid[3][4] =-rigid[1][3];
					rigid[3][5] =-rigid[2][3];
					rigid[3][6] = consi/2.0;
					rigid[4][4] = consa*coset*coset+consl*senot*senot;
					rigid[4][5] = rigid[1][2];
					rigid[4][6] =-rigid[1][3];
					rigid[5][5] = consa*senot*senot+consl*coset*coset;
					rigid[5][6] =-rigid[2][3];
					rigid[6][6] = consi;

					for(ievab=1; ievab<=nevab; ievab++)
						for(jevab=ievab; jevab<=nevab; jevab++)
							rigid[jevab][ievab] = rigid[ievab][jevab];
				}
				if(ntips[ielem] == 3) {
					consi = 3.0*young*xiner/xlong[ielem];
					consta = consi/xlong[ielem];
					consl = consta/xlong[ielem];
					rigid[1][1] = consa*coset*coset+consl*senot*senot;
					rigid[1][2] = (consa-consl)*senot*coset;
					rigid[1][3] =-consta*senot;
					rigid[1][4] =-rigid[1][1];
					rigid[1][5] =-rigid[1][2];
					rigid[2][2] = consa*senot*senot+consl*coset*coset;
					rigid[2][3] = consta*coset;
					rigid[2][4] =-rigid[1][2];
					rigid[2][5] =-rigid[2][2];
					rigid[3][3] = consi;
					rigid[3][4] =-rigid[1][3];
					rigid[3][5] =-rigid[2][3];
					rigid[4][4] = consa*coset*coset+consl*senot*senot;
					rigid[4][5] = rigid[1][2];
					rigid[5][5] = consa*senot*senot+consl*coset*coset;

					for(ievab=1; ievab<=nevab; ievab++)
						for(jevab=ievab; jevab<=nevab; jevab++)
							rigid[jevab][ievab] = rigid[ievab][jevab];
				}
				if(ntips[ielem] == 4) {
					consi = 3.0*young*xiner/xlong[ielem];
					consta = consi/xlong[ielem];
					consl = consta/xlong[ielem];
					rigid[1][1] = consa*coset*coset+consl*senot*senot;
					rigid[1][2] = (consa-consl)*senot*coset;
					rigid[1][4] =-rigid[1][1];
					rigid[1][5] =-rigid[1][2];
					rigid[1][6] =-consta*senot;
					rigid[2][2] = consa*senot*senot+consl*coset*coset;
					rigid[2][4] =-rigid[1][2];
					rigid[2][5] =-rigid[2][2];
					rigid[2][6] =consta*coset;
					rigid[4][4] = rigid[1][1];
					rigid[4][5] = rigid[1][2];
					rigid[4][6] =-rigid[1][6];
					rigid[5][5] = rigid[2][2];
					rigid[5][6] =-rigid[2][6];
					rigid[6][6] = consi;

					for(ievab=1; ievab<=nevab; ievab++)
						for(jevab=ievab; jevab<=nevab; jevab++)
							rigid[jevab][ievab] = rigid[ievab][jevab];
				}

			}
		}
		else {
			// Para 3D:

			consa = young*xarea/xlong[ielem];
			if(ntipo == 1) {
				coset=vectr[ielem][1];
				senot=vectr[ielem][2];
				cosef=vectr[ielem][3];
				senof=vectr[ielem][4];
				rigid[1][1] = consa*coset*coset*senof*senof;
				rigid[1][2] = consa*coset*senot*senof*senof;
				rigid[1][3] = consa*coset*cosef*senof;
				rigid[1][4] =-rigid[1][1];
				rigid[1][5] =-rigid[1][2];
				rigid[1][6] =-rigid[1][3];
				rigid[2][2] = consa*senot*senot*senof*senof;
				rigid[2][3] = consa*senot*cosef*senof;
				rigid[2][4] =-rigid[1][2];
				rigid[2][5] =-rigid[2][2];
				rigid[2][6] =-rigid[2][3];
				rigid[3][3] = consa*cosef*cosef;
				rigid[3][4] =-rigid[1][3];
				rigid[3][5] =-rigid[2][3];
				rigid[3][6] =-rigid[3][3];
				rigid[4][4] = rigid[1][1];
				rigid[4][5] = rigid[1][2];
				rigid[4][6] = rigid[1][3];
				rigid[5][5] = rigid[2][2];
				rigid[5][6] = rigid[2][3];
				rigid[6][6] = rigid[3][3];

				for(ievab=1; ievab<=nevab; ievab++)
					for(jevab=ievab; jevab<=nevab; jevab++)
						rigid[jevab][ievab] = rigid[ievab][jevab];
			}
			if(ntipo == 2) {
				for(ievab=1; ievab<=nevab; ievab++) {
					for(jevab=1; jevab<=nevab; jevab++) {
						rigid[ievab][jevab] = 0.0;
						girom[ievab][jevab] = 0.0;
					}
  				}
				girom[ 1][ 1] = vectr[ielem][1];
				girom[ 1][ 2] = vectr[ielem][2];
				girom[ 1][ 3] = vectr[ielem][3];
				girom[ 2][ 1] = vectr[ielem][4];
				girom[ 2][ 2] = vectr[ielem][5];
				girom[ 2][ 3] = vectr[ielem][6];
				girom[ 3][ 1] = vectr[ielem][7];
				girom[ 3][ 2] = vectr[ielem][8];
				girom[ 3][ 3] = vectr[ielem][9];
				girom[ 4][ 4] = girom[1][1];
				girom[ 4][ 5] = girom[1][2];
				girom[ 4][ 6] = girom[1][3];
				girom[ 5][ 4] = girom[2][1];
				girom[ 5][ 5] = girom[2][2];
				girom[ 5][ 6] = girom[2][3];
				girom[ 6][ 4] = girom[3][1];
				girom[ 6][ 5] = girom[3][2];
				girom[ 6][ 6] = girom[3][3];
				girom[ 7][ 7] = girom[1][1];
				girom[ 7][ 8] = girom[1][2];
				girom[ 7][ 9] = girom[1][3];
				girom[ 8][ 7] = girom[2][1];
				girom[ 8][ 8] = girom[2][2];
				girom[ 8][ 9] = girom[2][3];
				girom[ 9][ 7] = girom[3][1];
				girom[ 9][ 8] = girom[3][2];
				girom[ 9][ 9] = girom[3][3];
				girom[10][10] = girom[1][1];
				girom[10][11] = girom[1][2];
				girom[10][12] = girom[1][3];
				girom[11][10] = girom[2][1];
				girom[11][11] = girom[2][2];
				girom[11][12] = girom[2][3];
				girom[12][10] = girom[3][1];
				girom[12][11] = girom[3][2];
				girom[12][12] = girom[3][3];

				for(ievab=1; ievab<=nevab; ievab++)
					for(jevab=1; jevab<=nevab; jevab++)
						girtm[ievab][jevab] = girom[jevab][ievab];

				con4y = 4.0*young*xiner/xlong[ielem];
				con6y = 1.5*con4y/xlong[ielem];
				co12y = 2.0*con6y/xlong[ielem];
				con4z = 4.0*young*yiner/xlong[ielem];
				con6z = 1.5*con4z/xlong[ielem];
				co12z = 2.0*con6z/xlong[ielem];
				conjx = gcort*xjner/xlong[ielem];
				rigid[ 1][ 1] = consa;
				rigid[ 1][ 7] =-rigid[1][1];
				rigid[ 2][ 2] = co12z;
				rigid[ 2][ 6] = con6z;
				rigid[ 2][ 8] =-rigid[2][2];
				rigid[ 2][12] = rigid[2][6];
				rigid[ 3][ 3] = co12y;
				rigid[ 3][ 5] =-con6y;
				rigid[ 3][ 9] =-rigid[3][3];
				rigid[ 3][11] = rigid[3][5];
				rigid[ 4][ 4] = conjx;
				rigid[ 4][10] =-rigid[4][4];
				rigid[ 5][ 3] =-con6y;
				rigid[ 5][ 5] = con4y;
				rigid[ 5][ 9] =-rigid[3][5];
				rigid[ 5][11] = con4y/2.0;
				rigid[ 6][ 2] = con6z;
				rigid[ 6][ 6] = con4z;
				rigid[ 6][ 8] =-rigid[2][6];
				rigid[ 6][12] = con4z/2.0;
				rigid[ 7][ 1] =-rigid[1][1];
				rigid[ 7][ 7] = consa;
				rigid[ 8][ 2] =-rigid[2][2];
				rigid[ 8][ 6] =-rigid[2][6];
				rigid[ 8][ 8] = co12z;
				rigid[ 8][12] =-rigid[2][6];
				rigid[ 9][ 3] =-rigid[3][3];
				rigid[ 9][ 5] =-rigid[3][5];
				rigid[ 9][ 9] = co12y;
				rigid[ 9][11] =-rigid[3][5];
				rigid[10][ 4] =-rigid[4][4];
				rigid[10][10] = conjx;
				rigid[11][ 3] = rigid[3][5];
				rigid[11][ 5] = con4y/2.0;
				rigid[11][ 9] =-rigid[3][5];
				rigid[11][11] = con4y;
				rigid[12][ 2] = rigid[2][6];
				rigid[12][ 6] = con4z/2.0;
				rigid[12][ 8] =-rigid[2][6];
				rigid[12][12] = con4z;

				Comult(nevab, nevab, nevab, tempr, girtm, rigid);
				Comult(nevab, nevab, nevab, rigid, tempr, girom);
			}
            if(ntipo == 3) {
				for(ievab=1; ievab<=nevab; ievab++) {
					for(jevab=1; jevab<=nevab; jevab++) {
						rigid[ievab][jevab] = 0.0;
						girom[ievab][jevab] = 0.0;
					}
  				}
				girom[ 1][ 1] = vectr[ielem][1];
				girom[ 1][ 2] = vectr[ielem][2];
				girom[ 1][ 3] = vectr[ielem][3];
				girom[ 2][ 1] = vectr[ielem][4];
				girom[ 2][ 2] = vectr[ielem][5];
				girom[ 2][ 3] = vectr[ielem][6];
				girom[ 3][ 1] = vectr[ielem][7];
				girom[ 3][ 2] = vectr[ielem][8];
				girom[ 3][ 3] = vectr[ielem][9];
				girom[ 4][ 4] = girom[1][1];
				girom[ 4][ 5] = girom[1][2];
				girom[ 4][ 6] = girom[1][3];
				girom[ 5][ 4] = girom[2][1];
				girom[ 5][ 5] = girom[2][2];
				girom[ 5][ 6] = girom[2][3];
				girom[ 6][ 4] = girom[3][1];
				girom[ 6][ 5] = girom[3][2];
				girom[ 6][ 6] = girom[3][3];
				girom[ 7][ 7] = girom[1][1];
				girom[ 7][ 8] = girom[1][2];
				girom[ 7][ 9] = girom[1][3];
				girom[ 8][ 7] = girom[2][1];
				girom[ 8][ 8] = girom[2][2];
				girom[ 8][ 9] = girom[2][3];
				girom[ 9][ 7] = girom[3][1];
				girom[ 9][ 8] = girom[3][2];
				girom[ 9][ 9] = girom[3][3];
				girom[10][10] = girom[1][1];
				girom[10][11] = girom[1][2];
				girom[10][12] = girom[1][3];
				girom[11][10] = girom[2][1];
				girom[11][11] = girom[2][2];
				girom[11][12] = girom[2][3];
				girom[12][10] = girom[3][1];
				girom[12][11] = girom[3][2];
				girom[12][12] = girom[3][3];

				for(ievab=1; ievab<=nevab; ievab++)
					for(jevab=1; jevab<=nevab; jevab++)
						girtm[ievab][jevab] = girom[jevab][ievab];
				if (ntips[ielem] == 1) {
					rigid[1][1] = consa;
					rigid[7][7] = rigid[1][1];
					rigid[1][7] =-rigid[1][1];
					rigid[7][1] =-rigid[1][1];
				}
				if (ntips[ielem] == 2){
					con4y = 4.0*young*xiner/xlong[ielem];
					con6y = 1.5*con4y/xlong[ielem];
					co12y = 2.0*con6y/xlong[ielem];
					con4z = 4.0*young*yiner/xlong[ielem];
					con6z = 1.5*con4z/xlong[ielem];
					co12z = 2.0*con6z/xlong[ielem];
					conjx = gcort*xjner/xlong[ielem];
					rigid[ 1][ 1] = consa;
					rigid[ 1][ 7] =-rigid[1][1];
					rigid[ 2][ 2] = co12z;
					rigid[ 2][ 6] = con6z;
					rigid[ 2][ 8] =-rigid[2][2];
					rigid[ 2][12] = rigid[2][6];
					rigid[ 3][ 3] = co12y;
					rigid[ 3][ 5] =-con6y;
					rigid[ 3][ 9] =-rigid[3][3];
					rigid[ 3][11] = rigid[3][5];
					rigid[ 4][ 4] = conjx;
					rigid[ 4][10] =-rigid[4][4];
					rigid[ 5][ 3] =-con6y;
					rigid[ 5][ 5] = con4y;
					rigid[ 5][ 9] =-rigid[3][5];
					rigid[ 5][11] = con4y/2.0;
					rigid[ 6][ 2] = con6z;
					rigid[ 6][ 6] = con4z;
					rigid[ 6][ 8] =-rigid[2][6];
					rigid[ 6][12] = con4z/2.0;
					rigid[ 7][ 1] =-rigid[1][1];
					rigid[ 7][ 7] = consa;
					rigid[ 8][ 2] =-rigid[2][2];
					rigid[ 8][ 6] =-rigid[2][6];
					rigid[ 8][ 8] = co12z;
					rigid[ 8][12] =-rigid[2][6];
					rigid[ 9][ 3] =-rigid[3][3];
					rigid[ 9][ 5] =-rigid[3][5];
					rigid[ 9][ 9] = co12y;
					rigid[ 9][11] =-rigid[3][5];
					rigid[10][ 4] =-rigid[4][4];
					rigid[10][10] = conjx;
					rigid[11][ 3] = rigid[3][5];
					rigid[11][ 5] = con4y/2.0;
					rigid[11][ 9] =-rigid[3][5];
					rigid[11][11] = con4y;
					rigid[12][ 2] = rigid[2][6];
					rigid[12][ 6] = con4z/2.0;
					rigid[12][ 8] =-rigid[2][6];
					rigid[12][12] = con4z;
				}
				if(ntips[ielem] == 3) {
					con4y = 3.0*young*xiner/xlong[ielem];
					con6y = con4y/xlong[ielem];
					co12y = con6y/xlong[ielem];
					con4z = 3.0*young*yiner/xlong[ielem];
					con6z = con4z/xlong[ielem];
					co12z = con6z/xlong[ielem];
                     //jacob. modificacion, valor tenía cero
					conjx = gcort*xjner/xlong[ielem];
					rigid[ 1][ 1] = consa;
					rigid[ 1][ 7] =-rigid[1][1];
					rigid[ 2][ 2] = co12z;
					rigid[ 2][ 6] = con6z;
					rigid[ 2][ 8] =-rigid[2][2];
					rigid[ 3][ 3] = co12y;
					rigid[ 3][ 5] =-con6y;
					rigid[ 3][ 9] =-rigid[3][3];
					rigid[ 4][ 4] = conjx;
					rigid[ 5][ 3] =-con6y;
					rigid[ 5][ 5] = con4y;
					rigid[ 5][ 9] =-rigid[3][5];
					rigid[ 6][ 2] = con6z;
					rigid[ 6][ 6] = con4z;
					rigid[ 6][ 8] =-rigid[2][6];
					rigid[ 7][ 1] =-rigid[1][1];
					rigid[ 7][ 7] = consa;
					rigid[ 8][ 2] =-rigid[2][2];
					rigid[ 8][ 6] =-rigid[2][6];
					rigid[ 8][ 8] = co12z;
					rigid[ 9][ 3] =-rigid[3][3];
					rigid[ 9][ 5] =-rigid[3][5];
					rigid[ 9][ 9] = co12y;
				}
				if (ntips[ielem] == 4){
					con4y = 3.0*young*xiner/xlong[ielem];
					con6y = con4y/xlong[ielem];
					co12y = con6y/xlong[ielem];
					con4z = 3.0*young*yiner/xlong[ielem];
					con6z = con4z/xlong[ielem];
					co12z = con6z/xlong[ielem];
                    //jacob. modificacion, valor tenía cero
					conjx = gcort*xjner/xlong[ielem];
					rigid[ 1][ 1] = consa;
					rigid[ 1][ 7] =-consa;
					rigid[ 2][ 2] = co12z;
					rigid[ 2][ 8] =-co12z;
					rigid[ 2][12] = con6z;
					rigid[ 3][ 3] = co12y;
					rigid[ 3][ 9] =-co12y;
					rigid[ 3][11] =-con6y;
					rigid[ 7][ 1] =-consa;
					rigid[ 7][ 7] = consa;
					rigid[ 8][ 2] =-co12z;
					rigid[ 8][ 8] = co12z;
					rigid[ 8][12] =-con6z;
					rigid[ 9][ 3] =-co12y;
					rigid[ 9][ 9] = co12y;
					rigid[ 9][11] = con6y;
					rigid[10][10] = conjx;
					rigid[11][ 3] =-con6y;
					rigid[11][ 9] = con6y;
					rigid[11][11] = con4y;
					rigid[12][ 2] = con6z;
					rigid[12][ 8] =-con6z;
					rigid[12][12] = con4z;
				}
				Comult(nevab, nevab, nevab, tempr, girtm, rigid);
				Comult(nevab, nevab, nevab, rigid, tempr, girom);
			}
		}

		for(ievab=1; ievab<=nevab; ievab++)
			for(jevab=1; jevab<=nevab; jevab++)
				srmat[ielem][ievab][jevab]=rigid[ievab][jevab];

		if(isale >= 1) {
			fprintf(fp16, "Elemento %d\n", ielem);

			for(ievab=1; ievab<=nevab; ievab++) {
				for(jevab=1; jevab<=nevab; jevab++)
					fprintf(fp16, "%lf\t", rigid[ievab][jevab]);

				fprintf(fp16, "\n");
			}

			fprintf(fp16, "\n");
		}
	}
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void ApoyoInclinado(int ielem, int ndime, int ngdln, int nevab, int isale, int nincl, int** lnods, int* lincl,double** girom,
					double** girtm, double** tempr, double** rigid, double* xincl)
{
	int    lnod1,lnod2,iincl,ievab,jevab;
	double senox,cosen;

	lnod1 = lnods[ielem][1];
	lnod2 = lnods[ielem][2];

	for(ievab=1; ievab <= nevab; ievab++){
		for(jevab=1; jevab <= nevab; jevab++)
			girom[ievab][jevab] = 0.0;

		girom[ievab][ievab] = 1.0;
	}

	if(ndime == 2){
		for(iincl=1; iincl<=nincl; iincl++) {
			if(lnod1 == lincl[iincl]) {
				senox = sin(xincl[iincl]);
				cosen = cos(xincl[iincl]);
				girom[1][1] = cosen;
				girom[1][2] =-senox;
				girom[2][1] = senox;
				girom[2][2] = cosen;
			}

			if(lnod2 == lincl[iincl]){
				senox = sin(xincl[iincl]);
				cosen = cos(xincl[iincl]);
				girom[ngdln+1][ngdln+1] = cosen;
				girom[ngdln+1][ngdln+2] =-senox;
				girom[ngdln+2][ngdln+1] = senox;
				girom[ngdln+2][ngdln+2] = cosen;
			}
		}
	}

	for(ievab=1; ievab <= nevab; ievab++)
		for(jevab=1; jevab <= nevab; jevab++)
			girtm[ievab][jevab] = girom[jevab][ievab];

	Comult(nevab, nevab, nevab, tempr, girtm, rigid);
	Comult(nevab, nevab, nevab, rigid, tempr, girom);

	/*	if(isale >= 1) {
	 fprintf(fp16, "Elemento_%d\n", ielem);
	 for(ievab=1; ievab<=nevab; ievab++) {
	 for(jevab=1; jevab<=nevab; jevab++)
	 fprintf(fp16, "%lf\t", rigid[ievab][jevab]);
	 fprintf(fp16, "\n");
	 }
	 }
	 */
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
/************************************************/
void Masaconsis(int nelem, int nevab, int ndime, int ntipo, int isale, int* matnu, int* ntips,
				double*** smmat,double** props, double* xlong, double** vectr, double** xmasa,
				double** girom,double** girtm, double** tempr)
/******************************************************************************
 Calcula matriz de Masa Consistente de cada elemento.
 ******************************************************************************/
{
	int ielem, lprop, ievab, jevab;
	double pesoe,coset, senot, consa, consb;

	for(ielem=1; ielem<=nelem; ielem++)
		for(ievab=1; ievab<=nevab; ievab++)
			for(jevab=1; jevab<=nevab; jevab++)
				smmat[ielem][ievab][jevab]=0;

	// Ciclo sobre cada elemento:
	for(ielem=1; ielem<=nelem; ielem++) {

		lprop = abs(matnu[ielem]);

		// Inicializa la matriz de rigidez elemental:
		for(ievab=1; ievab<=nevab; ievab++)
			for(jevab=1; jevab<=nevab; jevab++)
				xmasa[ievab][jevab] = 0.0;

		if(ndime == 2) {
            pesoe = props[lprop][7];
			//pesoe = props[lprop][4];
			// Para 2D:
			coset = vectr[ielem][1];
			senot = vectr[ielem][2];
			// Calculo de la matriz de rigidez:
			if(ntipo == 1) {
				consa = pesoe*xlong[ielem]/3.0;
				consb = consa/2.0;
				xmasa[1][1] = consa;
				xmasa[1][3] = consb;
				xmasa[2][2] = consa;
				xmasa[2][4] = consb;
				xmasa[3][1] = consb;
				xmasa[3][3] = consa;
				xmasa[4][2] = consb;
				xmasa[4][4] = consa;
			}

			if(ntipo == 2) {
				consa = pesoe*xlong[ielem]/420.0;
				xmasa[1][1] = consa*(140*coset*coset+156*senot*senot);
				xmasa[1][2] =-consa*16*senot*coset;
				xmasa[1][3] =-consa*22*xlong[ielem]*senot;
				xmasa[1][4] = consa*(70*coset*coset+54*senot*senot);
				xmasa[1][5] =-xmasa[1][2];
				xmasa[1][6] = consa*13*xlong[ielem]*senot;
				xmasa[2][2] = consa*(156*coset*coset+140*senot*senot);
				xmasa[2][3] = consa*22*xlong[ielem]*coset;
				xmasa[2][4] =-xmasa[1][2];
				xmasa[2][5] = consa*(54*coset*coset+70*senot*senot);
				xmasa[2][6] =-consa*13*xlong[ielem]*coset;
				xmasa[3][3] = consa*4*xlong[ielem]*xlong[ielem];
				xmasa[3][4] =-consa*13*senot*xlong[ielem];
				xmasa[3][5] = consa*13*coset*xlong[ielem];
				xmasa[3][6] =-consa*3*xlong[ielem]*xlong[ielem];
				xmasa[4][4] = xmasa[1][1];
				xmasa[4][5] = xmasa[1][2];
				xmasa[4][6] =-xmasa[1][3];
				xmasa[5][5] = xmasa[2][2];
				xmasa[5][6] =-xmasa[2][3];
				xmasa[6][6] = xmasa[3][3];

				for(ievab=1; ievab<=nevab; ievab++)
					for(jevab=ievab; jevab<=nevab; jevab++)
						xmasa[jevab][ievab] = xmasa[ievab][jevab];
			}

			if(ntipo == 3) {
				if(ntips[ielem] == 1) {
					consa = pesoe*xlong[ielem]/3.0;
					consb = consa/2.0;
					xmasa[1][1] = consa;
					xmasa[1][3] = consb;
					xmasa[2][2] = consa;
					xmasa[2][4] = consb;
					xmasa[3][1] = consb;
					xmasa[3][3] = consa;
					xmasa[4][2] = consb;
					xmasa[4][4] = consa;
				}

				if(ntips[ielem] == 2) {
					consa = pesoe*xlong[ielem]/420.0;
					xmasa[1][1] = consa*(140*coset*coset+156*senot*senot);
					xmasa[1][2] =-consa*16*senot*coset;
					xmasa[1][3] =-consa*22*xlong[ielem]*senot;
					xmasa[1][4] = consa*(70*coset*coset+54*senot*senot);
					xmasa[1][5] =-xmasa[1][2];
					xmasa[1][6] = consa*13*xlong[ielem]*senot;
					xmasa[2][2] = consa*(156*coset*coset+140*senot*senot);
					xmasa[2][3] = consa*22*xlong[ielem]*coset;
					xmasa[2][4] =-xmasa[1][2];
					xmasa[2][5] = consa*(54*coset*coset+70*senot*senot);
					xmasa[2][6] =-consa*13*xlong[ielem]*coset;
					xmasa[3][3] = consa*4*xlong[ielem]*xlong[ielem];
					xmasa[3][4] =-consa*13*senot*xlong[ielem];
					xmasa[3][5] = consa*13*coset*xlong[ielem];
					xmasa[3][6] =-consa*3*xlong[ielem]*xlong[ielem];
					xmasa[4][4] = xmasa[1][1];
					xmasa[4][5] = xmasa[1][2];
					xmasa[4][6] =-xmasa[1][3];
					xmasa[5][5] = xmasa[2][2];
					xmasa[5][6] =-xmasa[2][3];
					xmasa[6][6] = xmasa[3][3];

					for(ievab=1; ievab<=nevab; ievab++)
						for(jevab=ievab; jevab<=nevab; jevab++)
							xmasa[jevab][ievab] = xmasa[ievab][jevab];
				}

				if(ntips[ielem] == 3) {
					consa = pesoe*xlong[ielem]/840.0;
					xmasa[1][1] = consa*(280*coset*coset+408*senot*senot);
					xmasa[1][2] =-consa*128*senot*coset;
					xmasa[1][3] =-consa*72*xlong[ielem]*senot;
					xmasa[1][4] = consa*(140*coset*coset+117*senot*senot);
					xmasa[1][5] = consa*23*senot*coset;
					xmasa[2][2] = consa*(408*coset*coset+280*senot*senot);
					xmasa[2][3] = consa*72*xlong[ielem]*coset;
					xmasa[2][4] = consa*23*senot*coset;
					xmasa[2][5] = consa*(117*coset*coset+140*senot*senot);
					xmasa[3][3] = consa*16*xlong[ielem]*xlong[ielem];
					xmasa[3][4] =-consa*33*senot*xlong[ielem];
					xmasa[3][5] = consa*33*coset*xlong[ielem];
					xmasa[4][4] = consa*(280*coset*coset+198*senot*senot);
					xmasa[4][5] = consa*82*senot*coset;
					xmasa[5][5] = consa*(198*coset*coset+280*senot*senot);

					for(ievab=1; ievab<=nevab; ievab++)
						for(jevab=ievab; jevab<=nevab; jevab++)
							xmasa[jevab][ievab] = xmasa[ievab][jevab];
				}
			}
		}
		else {
			// Para 3D:
			pesoe = props[lprop][7];
			if(ntipo == 1) {
				consa = pesoe*xlong[ielem]/3.0;
				consb = consa/2.0;
				xmasa[1][1] = consa;
				xmasa[1][4] = consb;
				xmasa[2][2] = consa;
				xmasa[2][5] = consb;
				xmasa[3][3] = consa;
				xmasa[3][6] = consb;
				xmasa[4][1] = consb;
				xmasa[4][4] = consa;
				xmasa[5][2] = consb;
				xmasa[5][5] = consa;
				xmasa[6][3] = consb;
				xmasa[6][6] = consa;

				for(ievab=1; ievab<=nevab; ievab++)
					for(jevab=ievab; jevab<=nevab; jevab++)
						xmasa[jevab][ievab] = xmasa[ievab][jevab];
			}
			if(ntipo == 2) {
				for(ievab=1; ievab<=nevab; ievab++) {
					for(jevab=1; jevab<=nevab; jevab++) {
						xmasa[ievab][jevab] = 0.0;
						girom[ievab][jevab] = 0.0;
					}
  				}
				girom[ 1][ 1] = vectr[ielem][1];
				girom[ 1][ 2] = vectr[ielem][2];
				girom[ 1][ 3] = vectr[ielem][3];
				girom[ 2][ 1] = vectr[ielem][4];
				girom[ 2][ 2] = vectr[ielem][5];
				girom[ 2][ 3] = vectr[ielem][6];
				girom[ 3][ 1] = vectr[ielem][7];
				girom[ 3][ 2] = vectr[ielem][8];
				girom[ 3][ 3] = vectr[ielem][9];
				girom[ 4][ 4] = girom[1][1];
				girom[ 4][ 5] = girom[1][2];
				girom[ 4][ 6] = girom[1][3];
				girom[ 5][ 4] = girom[2][1];
				girom[ 5][ 5] = girom[2][2];
				girom[ 5][ 6] = girom[2][3];
				girom[ 6][ 4] = girom[3][1];
				girom[ 6][ 5] = girom[3][2];
				girom[ 6][ 6] = girom[3][3];
				girom[ 7][ 7] = girom[1][1];
				girom[ 7][ 8] = girom[1][2];
				girom[ 7][ 9] = girom[1][3];
				girom[ 8][ 7] = girom[2][1];
				girom[ 8][ 8] = girom[2][2];
				girom[ 8][ 9] = girom[2][3];
				girom[ 9][ 7] = girom[3][1];
				girom[ 9][ 8] = girom[3][2];
				girom[ 9][ 9] = girom[3][3];
				girom[10][10] = girom[1][1];
				girom[10][11] = girom[1][2];
				girom[10][12] = girom[1][3];
				girom[11][10] = girom[2][1];
				girom[11][11] = girom[2][2];
				girom[11][12] = girom[2][3];
				girom[12][10] = girom[3][1];
				girom[12][11] = girom[3][2];
				girom[12][12] = girom[3][3];

				for(ievab=1; ievab<=nevab; ievab++)
					for(jevab=1; jevab<=nevab; jevab++)
						girtm[ievab][jevab] = girom[jevab][ievab];

				consa = pesoe*xlong[ielem]/420.0;
				xmasa[ 1][ 1] = consa*140.;
				xmasa[ 1][ 7] = consa*70.;
				xmasa[ 2][ 2] = consa*156.0;
				xmasa[ 2][ 6] = consa*22.0*xlong[ielem];
				xmasa[ 2][ 8] = consa*54.0;
				xmasa[ 2][12] =-consa*13*xlong[ielem];
				xmasa[ 3][ 3] = xmasa[2][2];
				xmasa[ 3][ 5] = xmasa[2][6];
				xmasa[ 3][ 9] = xmasa[2][8];
				xmasa[ 3][11] = xmasa[2][12];
				xmasa[ 4][ 4] = consa*4.0*xlong[ielem]*xlong[ielem];
				xmasa[ 4][10] =-consa*3.0*xlong[ielem]*xlong[ielem];
				xmasa[ 5][ 3] = xmasa[3][5];
				xmasa[ 5][ 5] = xmasa[4][4];
				xmasa[ 5][ 9] = consa*13.*xlong[ielem];
				xmasa[ 5][11] = xmasa[4][10];
				xmasa[ 6][ 2] = xmasa[2][6];
				xmasa[ 6][ 6] = xmasa[5][5];
				xmasa[ 6][ 8] = xmasa[5][9];
				xmasa[ 6][12] = xmasa[5][11];
				xmasa[ 7][ 1] = xmasa[1][7];
				xmasa[ 7][ 7] = xmasa[1][1];
				xmasa[ 8][ 2] = xmasa[2][8];
				xmasa[ 8][ 6] = xmasa[6][8];
				xmasa[ 8][ 8] = xmasa[2][2];
				xmasa[ 8][12] =-xmasa[2][6];
				xmasa[ 9][ 3] = xmasa[3][9];
				xmasa[ 9][ 5] = xmasa[5][9];
				xmasa[ 9][ 9] = xmasa[3][3];
				xmasa[ 9][11] = xmasa[8][12];
				xmasa[10][ 4] = xmasa[4][10];
				xmasa[10][10] = xmasa[4][4];
				xmasa[11][ 3] = xmasa[3][11];
				xmasa[11][ 5] = xmasa[5][11];
				xmasa[11][ 9] = xmasa[9][11];
				xmasa[11][11] = xmasa[5][5];
				xmasa[12][ 2] = xmasa[2][12];
				xmasa[12][ 6] = xmasa[6][12];
				xmasa[12][ 8] = xmasa[8][12];
				xmasa[12][12] = xmasa[6][6];

				Comult(nevab, nevab, nevab, tempr, girtm, xmasa);
				Comult(nevab, nevab, nevab, xmasa, tempr, girom);
			}
            if(ntipo == 3) {
				for(ievab=1; ievab<=nevab; ievab++) {
					for(jevab=1; jevab<=nevab; jevab++) {
						xmasa[ievab][jevab] = 0.0;
						girom[ievab][jevab] = 0.0;
					}
  				}
				girom[ 1][ 1] = vectr[ielem][1];
				girom[ 1][ 2] = vectr[ielem][2];
				girom[ 1][ 3] = vectr[ielem][3];
				girom[ 2][ 1] = vectr[ielem][4];
				girom[ 2][ 2] = vectr[ielem][5];
				girom[ 2][ 3] = vectr[ielem][6];
				girom[ 3][ 1] = vectr[ielem][7];
				girom[ 3][ 2] = vectr[ielem][8];
				girom[ 3][ 3] = vectr[ielem][9];
				girom[ 4][ 4] = girom[1][1];
				girom[ 4][ 5] = girom[1][2];
				girom[ 4][ 6] = girom[1][3];
				girom[ 5][ 4] = girom[2][1];
				girom[ 5][ 5] = girom[2][2];
				girom[ 5][ 6] = girom[2][3];
				girom[ 6][ 4] = girom[3][1];
				girom[ 6][ 5] = girom[3][2];
				girom[ 6][ 6] = girom[3][3];
				girom[ 7][ 7] = girom[1][1];
				girom[ 7][ 8] = girom[1][2];
				girom[ 7][ 9] = girom[1][3];
				girom[ 8][ 7] = girom[2][1];
				girom[ 8][ 8] = girom[2][2];
				girom[ 8][ 9] = girom[2][3];
				girom[ 9][ 7] = girom[3][1];
				girom[ 9][ 8] = girom[3][2];
				girom[ 9][ 9] = girom[3][3];
				girom[10][10] = girom[1][1];
				girom[10][11] = girom[1][2];
				girom[10][12] = girom[1][3];
				girom[11][10] = girom[2][1];
				girom[11][11] = girom[2][2];
				girom[11][12] = girom[2][3];
				girom[12][10] = girom[3][1];
				girom[12][11] = girom[3][2];
				girom[12][12] = girom[3][3];

				for(ievab=1; ievab<=nevab; ievab++)
					for(jevab=1; jevab<=nevab; jevab++)
						girtm[ievab][jevab] = girom[jevab][ievab];
				if (ntips[ielem] == 1) {
					consa = pesoe*xlong[ielem]/3.0;
					consb = consa/2.0;
					xmasa[1][1] = consa;
					xmasa[1][4] = consb;
					xmasa[2][2] = consa;
					xmasa[2][5] = consb;
					xmasa[3][3] = consa;
					xmasa[3][6] = consb;
					xmasa[4][1] = consb;
					xmasa[4][4] = consa;
					xmasa[5][2] = consb;
					xmasa[5][5] = consa;
					xmasa[6][3] = consb;
					xmasa[6][6] = consa;
				}
				if (ntips[ielem] == 2){
					consa = pesoe*xlong[ielem]/420.0;
					xmasa[ 1][ 1] = consa*140.;
					xmasa[ 1][ 7] = consa*70.;
					xmasa[ 2][ 2] = consa*156.0;
					xmasa[ 2][ 6] = consa*22.0*xlong[ielem];
					xmasa[ 2][ 8] = consa*54.0;
					xmasa[ 2][12] =-consa*13*xlong[ielem];
					xmasa[ 3][ 3] = xmasa[2][2];
					xmasa[ 3][ 5] = xmasa[2][6];
					xmasa[ 3][ 9] = xmasa[2][8];
					xmasa[ 3][11] = xmasa[2][12];
					xmasa[ 4][ 4] = consa*4.0*xlong[ielem]*xlong[ielem];
					xmasa[ 4][10] =-consa*3.0*xlong[ielem]*xlong[ielem];
					xmasa[ 5][ 3] = xmasa[3][5];
					xmasa[ 5][ 5] = xmasa[4][4];
					xmasa[ 5][ 9] = consa*13.*xlong[ielem];
					xmasa[ 5][11] = xmasa[4][10];
					xmasa[ 6][ 2] = xmasa[2][6];
					xmasa[ 6][ 6] = xmasa[5][5];
					xmasa[ 6][ 8] = xmasa[5][9];
					xmasa[ 6][12] = xmasa[5][11];
					xmasa[ 7][ 1] = xmasa[1][7];
					xmasa[ 7][ 7] = xmasa[1][1];
					xmasa[ 8][ 2] = xmasa[2][8];
					xmasa[ 8][ 6] = xmasa[6][8];
					xmasa[ 8][ 8] = xmasa[2][2];
					xmasa[ 8][12] =-xmasa[2][6];
					xmasa[ 9][ 3] = xmasa[3][9];
					xmasa[ 9][ 5] = xmasa[5][9];
					xmasa[ 9][ 9] = xmasa[3][3];
					xmasa[ 9][11] = xmasa[8][12];
					xmasa[10][ 4] = xmasa[4][10];
					xmasa[10][10] = xmasa[4][4];
					xmasa[11][ 3] = xmasa[3][11];
					xmasa[11][ 5] = xmasa[5][11];
					xmasa[11][ 9] = xmasa[9][11];
					xmasa[11][11] = xmasa[5][5];
					xmasa[12][ 2] = xmasa[2][12];
					xmasa[12][ 6] = xmasa[6][12];
					xmasa[12][ 8] = xmasa[8][12];
					xmasa[12][12] = xmasa[6][6];              }
				if(ntips[ielem] == 3) {
					consa = pesoe*xlong[ielem]/840.0;
					xmasa[ 1][ 1] = consa*240.0;
					xmasa[ 1][ 7] = consa*140.0;
					xmasa[ 2][ 2] = consa*408.0;
					xmasa[ 2][ 6] = consa*72*xlong[ielem];
					xmasa[ 2][ 8] = consa*117.0;
					xmasa[ 3][ 3] = xmasa[2][2];
					xmasa[ 3][ 5] = xmasa[2][6];
					xmasa[ 3][ 9] = xmasa[2][8];
					xmasa[ 4][ 4] = consa*8.0*xlong[ielem]*xlong[ielem];
					xmasa[ 5][ 3] = xmasa[3][5];
					xmasa[ 5][ 5] = consa*16.0*xlong[ielem]*xlong[ielem];
					xmasa[ 5][ 9] = consa*33.0*xlong[ielem];
					xmasa[ 6][ 2] = xmasa[2][6];
					xmasa[ 6][ 6] = xmasa[5][5];
					xmasa[ 6][ 8] = xmasa[5][9];
					xmasa[ 7][ 1] = xmasa[1][7];
					xmasa[ 7][ 7] = xmasa[1][1];
					xmasa[ 8][ 2] = xmasa[2][8];
					xmasa[ 8][ 6] = xmasa[6][8];
					xmasa[ 8][ 8] = consa*198.0;
					xmasa[ 9][ 3] = xmasa[3][9];
					xmasa[ 9][ 5] = xmasa[5][9];
					xmasa[ 9][ 9] = xmasa[8][8];
				}
  				Comult(nevab, nevab, nevab, tempr, girtm, xmasa);
				Comult(nevab, nevab, nevab, xmasa, tempr, girom);
			}
		}

		for(ievab=1; ievab<=nevab; ievab++)
			for(jevab=1; jevab<=nevab; jevab++)
				smmat[ielem][ievab][jevab]=xmasa[ievab][jevab];

		if(isale >= 1) {
			fprintf(fp16, "Elemento %d\n", ielem);

			for(ievab=1; ievab<=nevab; ievab++) {
				for(jevab=1; jevab<=nevab; jevab++)
					fprintf(fp16, "%lf\t", xmasa[ievab][jevab]);

				fprintf(fp16, "\n");
			}

			fprintf(fp16, "\n");
		}
	}
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
/************************************************/
void Masaconsen(int nelem, int nevab, int ndime, int ntipo, int isale, int* matnu, int* ntips,
				double*** smmat,double** props, double* xlong, double** vectr, double** xmasa,
				double** girom,double** girtm, double** tempr)
/******************************************************************************
 Calcula matriz de Masa Concentrada de cada elemento.
 ******************************************************************************/
{
	int ielem, lprop, ievab, jevab;
	double pesoe, consa;

	for(ielem=1; ielem<=nelem; ielem++)
		for(ievab=1; ievab<=nevab; ievab++)
			for(jevab=1; jevab<=nevab; jevab++)
				smmat[ielem][ievab][jevab]=0;

	// Ciclo sobre cada elemento:
	for(ielem=1; ielem<=nelem; ielem++) {

		lprop = abs(matnu[ielem]);

		// Inicializa la matriz de rigidez elemental:
		for(ievab=1; ievab<=nevab; ievab++)
			for(jevab=1; jevab<=nevab; jevab++)
				xmasa[ievab][jevab] = 0.0;

		if(ndime == 2) {
			// Para 2D:
			pesoe = props[lprop][7];
            //pesoe = props[lprop][4];
			// Calculo de la matriz de masa:
			consa = pesoe*xlong[ielem]/2.0;
			if(ntipo == 1) {
				xmasa[1][1] = consa;
				xmasa[2][2] = consa;
				xmasa[3][3] = consa;
				xmasa[4][4] = consa;
			}

			if(ntipo == 2 || ntipo ==3) {
				xmasa[1][1] = consa;
				xmasa[2][2] = consa;
				xmasa[4][4] = consa;
				xmasa[5][5] = consa;
			}

		}
		else {
			// Para 3D:
			pesoe = props[lprop][7];
			consa = pesoe*xlong[ielem]/2.0;
			if(ntipo == 1) {
				xmasa[1][1] = consa;
				xmasa[2][2] = consa;
				xmasa[3][3] = consa;
				xmasa[4][4] = consa;
				xmasa[5][5] = consa;
				xmasa[6][6] = consa;

			}
			if(ntipo == 2 || ntipo ==3) {
				xmasa[ 1][ 1] = consa;
				xmasa[ 2][ 2] = consa;
				xmasa[ 3][ 3] = consa;
				xmasa[ 7][ 7] = consa;
				xmasa[ 8][ 8] = consa;
				xmasa[ 9][ 9] = consa;
			}
		}

		for(ievab=1; ievab<=nevab; ievab++)
			for(jevab=1; jevab<=nevab; jevab++)
				smmat[ielem][ievab][jevab]=xmasa[ievab][jevab];


		if(isale >= 1) {
			fprintf(fp16, "Elemento %d\n", ielem);

			for(ievab=1; ievab<=nevab; ievab++) {
				for(jevab=1; jevab<=nevab; jevab++)
					fprintf(fp16, "%lf\t", xmasa[ievab][jevab]);

				fprintf(fp16, "\n");
			}

			fprintf(fp16, "\n");
		}
	}
}
