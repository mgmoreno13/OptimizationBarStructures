// *****************************
// Rutinas para resolver sistemas de ecuaciones
// *****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "memoria.h"
#include "fuerzas.h"
#include "raros.h"
#include "rigidez.h"
#include "solver.h"
#include "tensiones.h"


/******************************************************************************
 Ensambla y resuelve el sistema de ecuaciones matriciales.
 Solucion del sistema de ecuaciones utilizando eliminacion gaussiana directa.
 ******************************************************************************/
int Solucion(int npnod, int nelem, int ncaso, int ngdln, int nnode, int npres, int ntipo, int ndime, int nincl,
			 int nreso, int nevab, int ntotv, int iwrit, int isale, int** lnods, int* nodpr, int** inpre,
			 int* iffix, int* lreso, int* lincl, int* ntips, double** presc, double* despl,
			 double* fixed, double** astif, double* aslod, double** carpp, double*** srmat, double** resor,
			 double** girom, double** girtm, double** tempr, double** rigid, double* xincl, double* react)
{
	int ipnod, ncont, ngish, igash, nsvab, mcont, nloca, igdln, ngush;
	double treac[12];


	Ensambla(&nsvab, ncaso, nelem, nnode, ngdln, npres, ndime, nincl, nreso, nevab, ntotv, iwrit, isale,
			 lnods, nodpr, inpre, iffix, lreso, lincl, presc, despl, fixed, astif, aslod, carpp, srmat,
			 resor, girom, girtm, tempr, rigid, xincl);
	if(Reduce(nsvab, iffix, fixed, astif, aslod)) return -1;
	Sustituye(nsvab, nelem, ntipo, ndime, ngdln, iffix, ntips, lnods, fixed, react, astif, aslod, despl);
	for(ipnod=1; ipnod<=npnod; ipnod++) {
		ncont = ipnod*ngdln;
		ngish = ncont-ngdln+1;

        if(iwrit ==1) {
			for(igash=ngish; igash<= ncont; igash++)
				fprintf(fp16, "%e\t", despl[igash]);
			fprintf(fp16, "\n");
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
			if(iwrit ==1) {
				for(igdln=1; igdln<=ngdln; igdln ++)
					fprintf(fp16,"%e\t",treac[igdln]);
				fprintf(fp16, "\n");
			}
		}
	}

	return 0;
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Ensambla(int* nsvab, int ncaso, int nelem, int nnode, int ngdln, int npres, int ndime, int nincl,
			  int nreso, int nevab, int ntotv, int iwrit, int isale,
			  int** lnods, int* nodpr, int** inpre, int* iffix, int* lreso, int* lincl, double** presc, double* despl,
			  double* fixed, double** astif, double* aslod, double** carpp, double*** srmat, double** resor,
			  double** girom, double** girtm, double** tempr, double** rigid, double* xincl)
/******************************************************************************
 Ensambla las matrices de rigidez y los vectores de fuerzas.
 ******************************************************************************/
{
	int  itotv, jtotv, ipres, nloca, igdln, ncont, ielem, inode, nodei,
	nrows, nrowe, jnode, nodej, jgdln, ncols, ncole, ireso;

	for(itotv=1; itotv<=ntotv; itotv++) {
		despl[itotv] = 0.0;
		fixed[itotv] = 0.0;
		iffix[itotv] = 0.0;

		for(jtotv=1; jtotv<=ntotv; jtotv++) astif[itotv][jtotv] = 0.0;
	}
	if(ncaso==1){
		for(itotv=1; itotv<=ntotv; itotv++)
			aslod[itotv] = 0.0;
	}
	// Ensambla vector de fuerzas elementales:
	*nsvab = 1;

	for(ipres=1; ipres<=npres; ipres++) {
		nloca = (nodpr[ipres]-1)*ngdln;

		for(igdln=1; igdln<=ngdln; igdln++) {
			ncont = nloca+igdln;
			iffix[ncont] = inpre[ipres][igdln];
			fixed[ncont] = presc[ipres][igdln];
		}
	}

	for(ielem=1; ielem<=nelem; ielem++) {

		if(nincl >0)ApoyoInclinado(ielem,ndime,ngdln,nevab,isale, nincl,lnods,
								   lincl,girom,girtm,tempr,rigid,xincl);

		for(inode=1; inode<=nnode; inode++) {
			nodei = lnods[ielem][inode];

			for(igdln=1; igdln<=ngdln; igdln++) {
				nrows = ((nodei-1)*ngdln)+igdln;
				nrowe = ((inode-1)*ngdln)+igdln;
				if(ncaso ==1)aslod[nrows] += carpp[ielem][nrowe];
				if(nrows > *nsvab) *nsvab = nrows;

				// Ensambla las matrices de rigideces elementales:
				for(jnode=1; jnode<=nnode; jnode++) {
					nodej = lnods[ielem][jnode];

					for(jgdln=1; jgdln<=ngdln; jgdln++) {
						ncols = (nodej-1)*ngdln+jgdln;
						ncole = (jnode-1)*ngdln+jgdln;
						astif[nrows][ncols] += srmat[ielem][nrowe][ncole];
					}
				}
			}
		}
	}

	if(nreso > 0){
		for(ireso=1; ireso<=nreso; ireso++)
			for(igdln=1; igdln<=ngdln; igdln++) {
				nrows = ((lreso[ireso]-1)*ngdln)+igdln;
				astif[nrows][nrows] += resor[igdln][ireso];
			}
	}

	if(isale == 2) {
		fprintf(fp16,"Matriz de rigidez global:\n");

		for(itotv=1; itotv<=ntotv; itotv++) {
			for(jtotv=1; jtotv<=ntotv; jtotv++)
				fprintf(fp16, "%lf\t", astif[itotv][jtotv]);

			fprintf(fp16, "\n");
		}
		if(iwrit ==1) {
			fprintf(fp16, "Vector de cargas global:\n");
			for(itotv=1; itotv<=ntotv; itotv ++)
				fprintf(fp16, "%lf\n", aslod[itotv]);
		}
	}
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
int Reduce(int neons, int* iffix, double* fixed, double** astif, double* aslod)
/******************************************************************************
 Reduce el sistema de ecuaciones globales por eliminacion gaussiana directa.
 ******************************************************************************/
{
	// Regresa 0 si no ocurrio ningun error.
	// Regresa -1 en caso de error de c√°lculo.

	int ieons, irows, ieon1, icols;
	double pivot, factr;

	for(ieons=1; ieons<=neons; ieons++) {

		if(iffix[ieons] != 1) {
			// Reduce ecuaciones:
			pivot = astif[ieons][ieons];

			if(fabs(pivot) < 1.0e-20) {
				printf("\nPivote incorrecto %lf.\nEcuacion %d.", pivot, ieons);
				fprintf(fp16,"\nPivote incorrecto %lf.\nEcuacion %d.", pivot, ieons);
				return -1;
			}

			if(ieons != neons) {
				ieon1 = ieons+1;

				for(irows=ieon1; irows<=neons; irows++) {
					factr=astif[irows][ieons]/pivot;

					if(factr != 0.0) {
						for(icols=ieons; icols<=neons; icols++)
							astif[irows][icols] = astif[irows][icols]-
							factr*astif[ieons][icols];
						aslod[irows] = aslod[irows]-factr*
						aslod[ieons];
					}
				}
			}
		}
		else {
          	// Ajusta vector de cargas para desplazamientos prescritos:
			for(irows=ieons; irows<=neons; irows ++) {
				aslod[irows] = aslod[irows]-astif[irows][ieons]*
				fixed[ieons];
				astif[irows][ieons] = 0.0;
			}
		}
	}

	return 0;
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Sustituye(int neons, int nelem, int ntipo, int ndime, int ngdln, int* iffix, int* ntips, int** lnods,
			   double* fixed, double* react, double** astif, double* aslod, double* despl)
/******************************************************************************
 Realiza la sustitucion hacia atras.
 ******************************************************************************/
{
	int nloca, ielem, ieons, neon1, nback, nbac1, icols;
	double pivot, resid;

	for(ieons=1; ieons<=neons; ieons++)
		react[ieons] = 0.0;

	neon1 = neons+1;

	for(ieons=1; ieons<=neons; ieons++) {
		nback = neon1-ieons;
		pivot = astif[nback][nback];
		resid = aslod[nback];

		if(nback != neons) {
			nbac1 = nback+1;

			for(icols=nbac1; icols<=neons; icols++)
				resid = resid-astif[nback][icols]*despl[icols];
		}

		if(iffix[nback] == 0) despl[nback] = resid/pivot;
		if(iffix[nback] == 1) despl[nback] = fixed[nback];
		if(iffix[nback] == 1) react[nback] = -resid;
	}

	if(ntipo == 3 && ndime ==2) {
		for(ielem=1; ielem<=nelem; ielem ++) {
			if(ntips[ielem] == 1) {
				nloca = (lnods[ielem][1]-1)*ngdln+3;
				react[nloca] = 0.0;
				nloca = (lnods[ielem][2]-1)*ngdln+3;
				react[nloca] = 0.0;
			}

			if(ntips[ielem] == 3) {
				nloca = (lnods[ielem][2]-1)*ngdln+3;
				react[nloca] = 0.0;
			}
		}
	}
	if(ntipo == 3 && ndime ==3) {
		for(ielem=1; ielem<=nelem; ielem ++) {
			if(ntips[ielem] == 1) {
				nloca = (lnods[ielem][1]-1)*ngdln+4;
				react[nloca] = 0.0;
				react[nloca+1] = 0.0;
				react[nloca+2] = 0.0;
				nloca = (lnods[ielem][2]-1)*ngdln+4;
				react[nloca] = 0.0;
				react[nloca+1] = 0.0;
				react[nloca+2] = 0.0;
			}

			if(ntips[ielem] == 3) {
				nloca = (lnods[ielem][2]-1)*ngdln+4;
				react[nloca] = 0.0;
				react[nloca+1] = 0.0;
				react[nloca+2] = 0.0;
			}
		}
	}
}


//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
/******************************************************************************
 Links with profile solver.
 ******************************************************************************/
void Linkin(int nelem, int nnode, int nevab, int ntotv, int npres, int ngdln, int& neqns, int& nwktl, int& mkoun,
			int** lnods,int* nodea, int* nodpr, int** inpre, int* iffix, int** leqns, int* maxad, int* mhigh,
			double** presc, double  *fixed, double*& stiff)
{
	int ipres, igdln, nloca, itotv, ielem, ievab, nn1, inode, ident, ne1,nposn;

	for(ielem=1; ielem<=nelem; ielem++) nodea[ielem] = nnode;
	for(itotv=1; itotv<=ntotv; itotv++) {fixed[itotv] = 0.0;iffix[itotv]=0;}

	for(ipres=1; ipres<=npres; ipres++) {
		for(igdln=1; igdln<=ngdln; igdln++) {
			if(inpre[ipres][igdln] == 1) {
				nloca = (nodpr[ipres]-1)*ngdln+igdln;
				iffix[nloca] = 1;
				fixed[nloca] = presc[ipres][igdln];
			}
		}
	}

	// Number of unknowns:
	neqns = 0;

	for(itotv=1; itotv<=ntotv; itotv++) {
		if(iffix[itotv] == 0) {
			neqns++;
			iffix[itotv] = neqns;
		}
		else iffix[itotv] = 0;
	}

	// Conectivity array leqns:
	for(ielem=1; ielem<=nelem; ielem++)
		for(ievab=1; ievab<=nevab; ievab++) leqns[ievab][ielem] = 0;

	for(ielem=1; ielem<=nelem; ielem++) {
		nn1 = nodea[ielem];
		ievab = 1;

		for(inode=1; inode<=nn1; inode++) {
			ident = lnods[ielem][inode] ;

			for(igdln=1; igdln<=ngdln; igdln++) {
				itotv = (ident-1)*ngdln+igdln;
				leqns[ievab][ielem] = iffix[itotv];
				ievab++;
			}
		}
	}

	// Loop over all elements:
	for(ielem=1; ielem<=nelem; ielem++) {
		nn1 = nodea[ielem];
		ne1 = nn1*ngdln;
		Colmht(ne1,ielem,leqns,mhigh);

	}

	/* Follows the standard profile technique (see Bathe textbook)
	 addreses of diagonal elements -  maxad array -  */
	Addres(neqns, &nwktl, &mkoun, maxad, mhigh);

	stiff = vector(1, nwktl+1, (char*)"stiff");

	for(itotv=1; itotv<=ntotv; itotv++) {
		nposn = iffix[itotv];

		if(nposn != 0) iffix[itotv] = 0;
		else iffix[itotv] = 1;
	}
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Colmht(int mevab,int ielem, int** leqns, int* mhigh)
/******************************************************************************
 Evaluates the column heights of sitffness matrix.
 ******************************************************************************/
{
	register int maxam, ievab, ieqns, jhigh;

	maxam = 100000000;

	for(ievab=1; ievab<=mevab; ievab++) {
		if(leqns[ievab][ielem] != 0)
			if((leqns[ievab][ielem]-maxam) < 0) maxam = leqns[ievab][ielem];
	}

	for(ievab=1; ievab<=mevab; ievab++) {
		ieqns = leqns[ievab][ielem];

		if(ieqns != 0) {
			jhigh = ieqns-maxam;
			if(jhigh > mhigh[ieqns]) mhigh[ieqns] = jhigh;
		}
	}
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Addres(int neqns, int* nwktl, int* mkoun, int* maxad, int* mhigh)
/******************************************************************************
 Evaluates addresses of diagonal elements.
 ******************************************************************************/
{
	register int neqnn, ieqnn, ieqns,ikoun ;

	neqnn = neqns+1;

	for(ieqnn=1; ieqnn<=neqnn; ieqnn++) maxad[ieqnn] = 0.0;

	maxad[1] = 1;
	maxad[2] = 2;
	ikoun = 0;

	if(neqns != 1) {
		for(ieqns=2; ieqns<=neqns; ieqns++) {
			if(mhigh[ieqns] > ikoun) ikoun = mhigh[ieqns];
			maxad[ieqns+1] = maxad[ieqns]+mhigh[ieqns]+1;
		}
	}

	*mkoun=ikoun+1;
	*nwktl = maxad[neqns+1]-maxad[1];
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Addban(int nelem, int ngdln,int ndime, int nnode, int nevab, int isale, int nincl, int nwktl, int nreso,
			int** lnods,  int** leqns, int* maxad, int* nodea,int* lincl, int* lreso, int* iffix, double** girom,
			double** girtm, double** tempr, double** rigid,double* stiff, double*** srmat, double* xincl,
			double** resor)
/******************************************************************************
 Assembly of total stiffness vector.
 ******************************************************************************/
{
	int 	ielem,   nn1,   ne1, kount, ievab, ieqns, imaxa, kevab, jevab, jeqns,
	ijeqn, isize, iwktl, igdln, ireso, inode, nloca, itotv;

	for(iwktl=1; iwktl<=nwktl+1; iwktl++) stiff[iwktl] = 0.0;

	for(ielem=1; ielem<=nelem; ielem++) {
		nn1 = nodea[ielem];
		ne1 = nn1*ngdln;

		if(nincl >0)ApoyoInclinado(ielem,ndime,ngdln,nevab,isale, nincl,lnods,
								   lincl,girom,girtm,tempr,rigid,xincl);


		kount = 0;

		for(ievab=1; ievab<=ne1; ievab++) {
			ieqns = leqns[ievab][ielem];

			if(ieqns > 0) {
				imaxa = maxad[ieqns];
				kevab = ievab;

				for(jevab=1; jevab<=ne1; jevab++) {
					jeqns = leqns[jevab][ielem];

					if(jeqns > 0) {
						ijeqn = ieqns-jeqns;

						if(ijeqn >=0) {
							isize = imaxa+ijeqn;
							stiff[isize] += srmat[ielem][ievab][jevab];
						}
					}

					kevab = kevab+ne1-jevab;
				}
			}

			kount = kount+ne1-ievab;
		}
	}

 	if(nreso > 0){
		for(ireso=1; ireso<=nreso; ireso++){
			for(ielem=1; ielem<=nelem; ielem++) {
				for(inode=1; inode<=nnode; inode++) {
			     	nloca = lnods[ielem][inode];
					fprintf(fp16,"%d,%d\n", nloca,lreso[ireso]);
					if(lreso[ireso] == nloca) break;
				}
				if(lreso[ireso] == nloca) break;
			}
      		for(igdln=1; igdln<=ngdln; igdln++) {
				ievab = (inode-1)*ngdln+igdln;
				itotv=(lreso[ireso]-1)*ngdln+igdln;
				if(iffix[itotv] == 0) {
					ieqns = leqns[ievab][ielem];
			   		imaxa = maxad[ieqns];
					fprintf(fp16,"%d,%lf\t", imaxa,resor[igdln][ireso]);
					stiff[imaxa] += resor[igdln][ireso];
				}
			}
		}
	}

}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
int Decomp(int neqns, int ishot, int* maxad, double* stiff)
/******************************************************************************
 Factorises (l)*(d)*(l) transpose of siffness matrix.
 ******************************************************************************/
{
	int 	ieqns, imaxa, lower, kuper, khigh, ksize, icoun, juper, jhigh, kmaxa,
	ndiag, ncolm, icolm, jmaxa;

	double count, bsumm, ratio;

	if(neqns == 1) return 0;

	for(ieqns=1; ieqns<=neqns; ieqns++) {
		imaxa = maxad[ieqns];
		lower = imaxa+1;
		kuper = maxad[ieqns+1]-1;
		khigh = kuper-lower;

		if(khigh > 0) {
			ksize = ieqns-khigh;
			icoun = 0;
			juper = kuper;

			for(jhigh=1; jhigh<=khigh; jhigh++) {
				icoun++;
				juper--;
				kmaxa = maxad[ksize];
				ndiag = maxad[ksize+1]-kmaxa-1;

				if(ndiag > 0.0) {
					ncolm = icoun;
					if(ndiag < ncolm) ncolm = ndiag;
					count = 0.0;
					for(icolm=1; icolm<=ncolm; icolm++)
						count = count+stiff[kmaxa+icolm]*
						stiff[juper+icolm];
					stiff[juper] = stiff[juper]-count;
				}

				ksize++;
			}
		}
  //      printf("factorizando ecuacion numero %d \t %d \n",ieqns,neqns);

		if(khigh >= 0) {
			ksize = ieqns;
			bsumm = 0.0;

			for(icolm=lower; icolm<=kuper; icolm++) {
				ksize--;
				jmaxa = maxad[ksize];

				if(stiff[jmaxa] == 0.0) {
					fprintf(fp16,"La matriz de rigidez no es def. positiva; ecuaci√≥n = %d, pivote = %e. \n",	ieqns, stiff[jmaxa]);
					return -1;
				}

				ratio = stiff[icolm]/stiff[jmaxa];
				bsumm+= stiff[icolm]*ratio;
				stiff[icolm] = ratio;
			}

			stiff[imaxa] = stiff[imaxa]-bsumm;
		}

		if(stiff[imaxa] <= 0.0) {
			if(ishot == 0) {
				fprintf(fp16,"La matriz de rigidez no es def. positiva; ecuacion = %d, pivote = %e. \n", ieqns, stiff[imaxa]);
                return -1;
			}

			if(stiff[imaxa] <= 0.0) stiff[imaxa] =-1.e-16 ;
		}
	}

	return 0;
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Profile(int nelem, int ndime, int npnod, int ngdln, int ncaso, int neqns, int iwrit,int ntotv,
			 int nincl, int nevab, int isale,
			 int** lnods, int* nodea, int* maxad, int* iffix,int* lincl, double** carpp, double* despl,
			 double* aslod, double* fixed, double* stiff, double*** srmat, double* react,
			 double** girom, double** girtm, double** tempr, double** rigid, double* xincl)
/******************************************************************************
 Calculate the incremental displacements.
 ******************************************************************************/
{
	int 	itotv, ielem, nn1, iposn, ievab, igdln, npoi, inode, ipoi, ipnod,
	ncont, ngish, igash;

	for(itotv=1; itotv<=ntotv; itotv++) despl[itotv] = 0.0;

	if(ncaso == 1) {
		for(ielem=1; ielem<=nelem; ielem++) {
			nn1 = nodea[ielem];
			for(inode=1; inode<=nn1; inode++) {
				iposn=lnods[ielem][inode];
				for(igdln=1; igdln<=ngdln; igdln++) {
					itotv = (iposn-1)*ngdln+igdln;
					ievab = (inode-1)*ngdln+igdln;
					despl[itotv] += carpp[ielem][ievab];
				}
			}
		}
	}
	else if(ncaso ==2){
		for(itotv=1; itotv <= ntotv; itotv++)
			despl[itotv]=aslod[itotv];
	}

    //Impresion de Cargas Nodales Equivalentes
	if(iwrit ==1) {
		for(inode=1; inode<=npnod; inode++) {
			itotv = (inode-1)*ngdln;
			fprintf(fp16, "%d \t ", inode);
			for(igdln=1; igdln<=ngdln; igdln++)
				fprintf(fp16, "%le \t ", despl[itotv+igdln]);
			fprintf(fp16, "\n");
		}
	}
	npoi = 0;


	for(itotv=1; itotv <= ntotv; itotv++) {
		if(iffix[itotv] != 1) {
			npoi++;
			aslod[npoi] = despl[itotv];
		}
	}

	for(ipoi=1; ipoi<=npoi; ipoi++) despl[ipoi] = aslod[ipoi];

	Redbak(neqns,maxad,stiff,aslod);
	npoi = 0;

	for(itotv=1; itotv<=ntotv; itotv++) {
		despl[itotv] = 0.0;

		if(iffix[itotv] == 1) despl[itotv] = fixed[itotv];
		else {
			npoi++;
			despl[itotv] = aslod[npoi];
		}
	}

	if(iwrit == 1){
		for(ipnod=1; ipnod<=npnod; ipnod++) {
			ncont = ipnod*ngdln;
			ngish = ncont-ngdln+1;
            fprintf(fp16, "Despl.Nodo\t%d\t", ipnod);
			for(igash=ngish; igash<=ncont; igash++)
				fprintf(fp16, "%le\t", despl[igash]);
			fprintf(fp16, " \n");
		}
	}

	Reacci(nelem,ndime,npnod,ngdln,nincl,nevab,ntotv,iwrit,isale,lnods,nodea,iffix,lincl,despl,fixed,srmat,
		   carpp,react,girom,girtm,tempr,rigid, xincl);
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void Redbak(int neqns, int* maxad, double* stiff, double* aslod)
/******************************************************************************
 For reduce and bak-substitute iteration vectors.
 ******************************************************************************/
{
	int ieqns, lower, kuper, jeqns, kmaxa, icolm, keqns;
	double sumcc;

	for(ieqns=1; ieqns<=neqns; ieqns++) {
		lower = maxad[ieqns]+1;
		kuper = maxad[ieqns+1]-1;

		if((kuper-lower) >= 0) {
			jeqns = ieqns;
			sumcc = 0.0 ;

			for(icolm=lower; icolm<=kuper; icolm++) {
				jeqns--;
				sumcc+= stiff[icolm]*aslod[jeqns];
			}

			aslod[ieqns] = aslod[ieqns]-sumcc;
		}
	}

	for(ieqns=1; ieqns<=neqns; ieqns++) {
		kmaxa = maxad[ieqns];
		aslod[ieqns] = aslod[ieqns]/stiff[kmaxa];
	}

	if(neqns == 1) return;

	jeqns = neqns;

	for(ieqns=2; ieqns<=neqns; ieqns++) {
		lower = maxad[jeqns]+1;
		kuper = maxad[jeqns+1]-1;

		if((kuper-lower) >= 0) {
			keqns = jeqns;
			for(icolm=lower; icolm<=kuper; icolm++){
				keqns--;
				aslod[keqns] = aslod[keqns]-stiff[icolm]*aslod[jeqns];
			}
		}

		jeqns--;
	}
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
int Gradc(int nelem, int npnod, int npres, int ngdln, int ntotv, int nnode, int ncaso, int iwrit,int ndime,
		  int nevab,  int isale, int nincl, int nreso, int** lnods, int** inpre, int* nodpr, int* iffix,
		  int* nodea, int* lincl, int* lreso, double* fixed, double** presc, double* despl, double* aslod,
		  double** carpp, double* vectp, double* vectu, double** girom, double** girtm, double** tempr,
		  double** rigid, double* xincl, double** resor, double*** srmat, double* react, double* vect1,
		  double* deslo)
/******************************************************************************
 Ensambla el vector de fuerzas obtener la para solucion de un sistema de
 ecuaciones con gradiente conjugado.
 ******************************************************************************/
{
	// Regresa 0 si no se encontr√≥ ning√∫n error.
	// Regresa -1 si el numero de iteraciones fu√© excedido.

	int 	 ipnod, ipres, igdln, nloca, itotv, ielem, ievab, nn1, kmaxc, inode,
	iposn, kcont, ncont, ngish, igash;
	double wvalu, wvold, wbeta, puval, alfa;
	double toler;

	toler = 1e-10;
	kmaxc = 10000*ntotv;

	for(ielem=1; ielem<=nelem; ielem++) nodea[ielem] = nnode;

	if(ncaso ==1){
		for(itotv=1; itotv<=ntotv; itotv++) {
			fixed[itotv] = 0.0;
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
		Fuefix(nelem,ndime,ngdln,nevab,isale,nincl,lnods,lincl,girom,girtm,tempr,rigid,xincl,nodea,
					iffix,fixed,srmat,carpp);
		for(ielem=1; ielem<=nelem; ielem++) {
			nn1 = nodea[ielem];
			for(inode=1; inode <= nn1; inode++) {
				iposn = lnods[ielem][inode];
				for(igdln=1; igdln<=ngdln; igdln++) {
					itotv = (iposn-1)*ngdln+igdln;
					ievab = (inode-1)*ngdln+igdln;
					aslod[itotv] = aslod[itotv]+carpp[ielem][ievab];
				}
			}
		}
    }

	for(itotv=1; itotv<=ntotv; itotv++){
		despl[itotv] = 0.0;
		if(iffix[itotv] !=0) aslod[itotv] = 0.0;
		vectp[itotv] = aslod[itotv];
	}
    if(iwrit ==1) {
		for(inode=1; inode<=npnod; inode++) {
			itotv = (inode-1)*ngdln;
			fprintf(fp16, "%d \t ", inode);
			for(igdln=1; igdln<=ngdln; igdln++)
				fprintf(fp16, "%le \t ", aslod[itotv+igdln]);
			fprintf(fp16, "\n");
		}
	}
	kcont = 0;
	wvold = 0.0;

	do {
		wvalu = 0.0;

		for(itotv=1; itotv<=ntotv; itotv++)
			wvalu+= aslod[itotv]*aslod[itotv];

		if(kcont != 0) {
			wbeta = wvalu/wvold;

			for(itotv=1; itotv<=ntotv; itotv++)
				vectp[itotv] = aslod[itotv]+wbeta*vectp[itotv];
		}

		EvaluaMatrizVector(nelem,ngdln,ndime,nevab,isale,nincl,nreso,ntotv,lnods,nodea,iffix,lincl,lreso,
						   girom,girtm,tempr,rigid,xincl,resor,vectu,vectp, vect1,deslo,srmat);
		puval = 0.0;

		for(itotv=1; itotv<=ntotv; itotv++)
			puval+= vectu[itotv]*vectp[itotv];

		if(iwrit ==1) printf("kcont=%d wvalu=%e  toler=%e iwrit =%d\n",kcont,wvalu,toler,iwrit);
        if(fabs(puval) < toler)break;
		else {
			alfa = wvalu/puval;
			for(itotv=1; itotv<=ntotv; itotv++) {
				despl[itotv]+= alfa*vectp[itotv];
				aslod[itotv]-= alfa*vectu[itotv];
			}

		}
		if(wvalu < toler) break;
		kcont++;
		wvold = wvalu;
	} while(kcont <= kmaxc);

	if(iwrit == 1) {
		for(ipnod =1; ipnod<=npnod; ipnod++) {
			ncont = ipnod*ngdln;
			ngish = ncont-ngdln+1;
			for(igash=ngish; igash<=ncont; igash++)
				if(iffix[igash] !=0) despl[igash]=fixed[igash];
			for(igash=ngish; igash<=ncont; igash++)
				fprintf(fp16, " %e\t", despl[igash]);
			fprintf(fp16, "\n");
		}
	}

Reacci(nelem,ndime,npnod,ngdln,nincl,nevab,ntotv,iwrit,isale,lnods,nodea,iffix,lincl,despl,fixed,srmat,
	   carpp,react,girom,girtm,tempr,rigid, xincl);
	return 0;
}

//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void EvaluaMatrizVector(int nelem, int ngdln,int ndime, int nevab, int isale, int nincl, int nreso,
                        int ntotv, int** lnods,  int* nodea, int* iffix, int* lincl,int* lreso,
						double** girom, double** girtm, double** tempr, double** rigid, double* xincl,
						double** resor, double* vectu, double* vectp, double* vect1, double* deslo,
						double*** srmat)
{
	int 	ievab, itotv, igdln, inode, nn1, nn2, nodei, nrows, nrowe, jevab,
	ielem, ireso;


	for(itotv=1; itotv<=ntotv; itotv++) vectu[itotv] = 0.0;

	for(ielem=1; ielem<=nelem; ielem++) {
		nn1 = nodea[ielem];
		nn2 = nn1*ngdln;

		if(nincl >0) ApoyoInclinado(ielem,ndime,ngdln,nevab,isale, nincl,lnods,
									lincl,girom,girtm,tempr,rigid,xincl);
		for(inode=1; inode<=nn1;inode++) {
			nodei = lnods[ielem][inode];

			for(igdln=1; igdln<=ngdln; igdln++) {
				nrows = ((nodei-1)*ngdln)+igdln;
				nrowe = ((inode-1)*ngdln)+igdln;
				if(iffix[(nodei-1)*ngdln+igdln]==0)
					deslo[nrowe] = vectp[nrows];
				else deslo[nrowe] = 0.0;
			}
		}

		for(ievab=1; ievab<=nn2; ievab++) {
			vect1[ievab] = 0.0;

			for(jevab=1; jevab<=nn2; jevab++)
				vect1[ievab] =vect1[ievab]+srmat[ielem][ievab][jevab]*deslo[jevab];
		}

		for(inode=1; inode<=nn1; inode++) {
			nodei = lnods[ielem][inode];

			for(igdln=1; igdln<=ngdln; igdln++) {
				nrows = ((nodei-1)*ngdln)+igdln;
				nrowe = ((inode-1)*ngdln)+igdln;
				vectu[nrows]+= vect1[nrowe];
			}
		}
	}

	if(nreso > 0) {
		for(ireso=1; ireso<=nreso; ireso++)
			for(igdln=1; igdln<=ngdln; igdln++) {
				nrows = ((lreso[ireso]-1)*ngdln)+igdln;
				vectu[nrows] += resor[igdln][ireso]*vectp[nrows];
			}
	}


	for(itotv=1; itotv<=ntotv; itotv++)
		if(iffix[itotv] !=0) vectu[itotv] = 0.0;
}
