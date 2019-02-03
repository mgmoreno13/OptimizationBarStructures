//
//  efi_CR.cpp
//  new_mecap
//
//  Created by Jacob Esau Salazar Solano on 5/24/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "efi_CR.h"
#include <stdio.h>

void Eficiencia_MarcosCR(int ielem, int nelem, int ndime, int ntipo, int* matnu,
                         int* ntips, double** props, double* xlong,double* fuerc,
                         double** arecr, int** indfu, double** fuerb,
                         double** fuefl, double* wcarg, double*** dtcon )
{
	int lmats,iprop,nbarr;
	double as[101],x[101],y[101];
	double cero=0;
    ////printf("Entra a eficiencia marcos de CR ielem = % d \n", ielem);

	//if(ndime== 2) lmats=(int)props[matnu[ielem]][6]+0.5;
	//if(ndime ==3) lmats=(int)props[matnu[ielem]][9]+0.5;
    lmats=(int)props[matnu[ielem]][9]+0.5;
	nbarr=(int)arecr[lmats][5]+0.5;
	for(iprop=1; iprop<=100; iprop++){
		as[iprop]=x[iprop]=y[iprop]=0.0;
	}
	for(iprop=1; iprop<=nbarr; iprop++){
		as[iprop]=dtcon[lmats][1][iprop];
		x[iprop]=dtcon[lmats][2][iprop];
		y[iprop]=dtcon[lmats][3][iprop];
	}
/*	if(ndime== 2)
		fuerc[ielem]= Flexcomp_concr(indfu[1][ielem],indfu[2][ielem],indfu[3][ielem],
									 indfu[4][ielem],xlong[ielem],arecr[lmats][1],
									 arecr[lmats][2],arecr[lmats][3],arecr[lmats][4],
									 nbarr,as,x,y,ntips[ielem],cero,fuerb[2][ielem],cero,
									 fuerb[4][ielem],cero,fuerb[6][ielem],
									 fuefl[ielem][1],cero,fuefl[ielem][2],
									 fuefl[ielem][4],cero,fuefl[ielem][5],
									 cero,fuefl[ielem][3],cero,cero,fuefl[ielem][6],cero);
	if(ndime== 3)                 */
		fuerc[ielem]= Flexcomp_concr(indfu[1][ielem],indfu[2][ielem],indfu[3][ielem],
									 indfu[4][ielem],xlong[ielem],arecr[lmats][1],
									 arecr[lmats][2],arecr[lmats][3],arecr[lmats][4],
									 nbarr,as,x,y,ntips[ielem],fuerb[1][ielem],fuerb[2][ielem],
									 fuerb[3][ielem],fuerb[4][ielem],fuerb[5][ielem],
									 fuerb[6][ielem],fuefl[ielem][1],fuefl[ielem][2],
									 fuefl[ielem][3],fuefl[ielem][7],fuefl[ielem][8],
									 fuefl[ielem][9],fuefl[ielem][4],fuefl[ielem][5],
									 fuefl[ielem][6],fuefl[ielem][10],fuefl[ielem][11],
									 fuefl[ielem][12]);
}

// Rutinas de Concreto Reforzado
//------------------------------------------------------------------------------

double Flexcomp_concr(int ifcy,int ifcz,int ifdy,int ifdz,double L, double b, double h, double fc, double fy, int nbarr, double *as, double *x, double *y,
					  int tb, double mfcy, double mfcz,double dfcy,double dfcz, double mfdy, double mfdz, double fxi, double fyi, double fzi, double fxf,
					  double fyf,  double fzf,  double mxi, double myi,  double mzi,  double mxf,  double myf, double mzf)
{
    double A,Iy,Iz,K,fam,rgy,rgz,rgmin,ast,E,EIy,EIz,KLsR,Mmaxy,Mmaxz,Mry,Mrz,P0,efi,Pc,pi=M_PI;
	int    i;

	////printf("ifcy=%d ifcz=%d ifdy=%d ifdz=%d dfcy=%lf dfcz=%lf mfcy=%lf mfcz=%lf mfdy=%lf mfdz=%lf \n",ifcy, ifcz,ifdy,ifdz,dfcy,dfcz,mfcy,mfcz,mfdy,mfdz);
	////printf("L=%lf b=%lf h=%lf fc=%lf fy=%lf nbarr=%d\n",L,b,h,fc,fy,nbarr);
	////printf("fxi =%lf fyi=%lf fzi =%lf fxf =%lf fyf=%lf fzf =%lf mxi=%lf myi=%lf mzi %lf mxf =%lf myf =%lf mzf =%lf\n",fxi,fyi,fzi,fxf,fyf,fzf,mxi,myi,mzi,mxf,myf,mzf);
	Mmaxy=0;
	Mmaxz=0;
	A=b*h;
	Iy=b*h*h*h/12.;
	Iz=h*b*b*b/12.;
	E=15000.*sqrt(fc);
	//jaco.pdr
	EIy=E*Iy;
	EIz=E*Iz;
	////printf("A=%lf Iy=%lf Iz=%lf tb=%d\n",A,Iy,Iz,tb);

    // Factor de Longitud Efectiva.
	switch( tb){
		case 1 : // Art-Art.
			K = 1.00;
			break;
		case 2 : // Emp-Emp.
			K = 0.65;
			break;
		case 3 : // Emp-Art.
			K = 0.80;
			break;
			//jacob.pdr
        	case 4: //  Art-Emp.
            	K= 0.80;
            	break;
	}

	// Obtiene Fuerza Axial Maxima.
	fam = -fxi;
	if( fabs(fam) < fabs( fxf)){
		fam = fxf;
	}

	// Radios de Giro y su Minimo.
	rgy = sqrt( Iy / A );
	rgz = sqrt( Iz / A );
	rgmin = rgy;
	if( rgmin > rgz ){
		rgmin = rgz;
	}


	KLsR=K*L/rgmin;

////printf("KLsR=%lf \n",KLsR);
    if (KLsR>100){printf("ERROR: Requiere analisis de segundo orden. KL/r= = %lf",KLsR); return 0;}
	//  Revision por esbeltez del Momento Maximo en Y.
	//jacob.pdr
	Pc=pi*pi*EIy/(K*L*K*L);
	cfmmcon( &Mmaxy, L, ifcz, ifdz, dfcz, mfcz, mfdz, fzi, fzf, myi, myf,fam,Pc,KLsR);
	//jacob.pdr
	Pc=pi*pi*EIz/(K*L*K*L);
	//  Revision por esbeltez del Momento Maximo en Z.
	cfmmcon( &Mmaxz, L, ifcy, ifdy, dfcy, mfcy, mfdy, fyi, fyf, mzi, mzf,fam,Pc,KLsR);
//printf("fam =%lf Mmaxy=%lf ,Mmaxz = %lf\n",fam,Mmaxy,Mmaxz);

	ast=0.0;
	for(i= 1 ; i<=nbarr; i++) ast+=as[i];

	//jacob.pdr
	//muttio Calculo de PO.
    //Problema: En algunos casos no deberia entrar al calculo de Bressler,
    // Se encontro que en la función de Bressler, PO se calculo distinto ¿?
	//Antes
	//P0 = - 0.08*( 0.85*fc*(b*h-ast) + ast*fy);
	//Despues: Modificacion Neto

	P0 = - 0.8*0.65*( 0.85*fc*(b*h-ast) + ast*fy );

	//jacob.pdr
    //printf("fam=%lf,  P0= %lf \n",fam,P0);
    if(fabs(fam)<0.1*fabs(P0)){
		//	//printf("solo flexion\n");
        ////printf("entro con My \n");
			Mry=flexion_r(b,h,fc,fy,nbarr,as,y,ast);
        ////printf("entro con Mz \n");
	    	Mrz=flexion_r(h,b,fc,fy,nbarr,as,x,ast);
////printf("uno Mry =%lf Mrz=%lf Mmaxy=%lf Mmaxz=%lf \n",Mry,Mrz,Mmaxy,Mmaxz);

		efi=fabs(Mmaxy)/fabs(Mry)+fabs(Mmaxz)/fabs(Mrz);
        if (efi >0.1){
            efi += fabs(fam/P0);
        }
		//jacob.pdr
		//neto.agregado
		//fuerc2[ielem]=0;
		////printf("\nEntra a flexion");
	}
	else{
 //       //printf("dos");
		efi=Flexcomp_Rect_ACI_Bress(fam,Mmaxy,Mmaxz,b,h,fc,fy,nbarr,as,x,y,ast);
		//jacob.pdr
		//neto.agregado
		//fuerc2[ielem]=Flexotension_Rect_ACI_CTHsu(fam,Mmaxy,Mmaxz,b,h,fc,fy,nbarr,as,x,y,ast);
//		//printf("\nEntra a Flexcomp_Rect_ACI_Bress");
	}
   // //printf("efi=%lf \n",efi);
	return efi;
}

//---------------------------------------------------------------------
//jacob.pdr
//neto.agregado

double Flexotension_Rect_ACI_CTHsu(double Pu, double Mux, double Muy, double b, double h, double fc, double fy, double nbarr, double *as, double *x, double *y,double ast)
{
	/*REVISION POR Y FLEXOCOMPRESIÓN FLEXOTENSIÓN DE COLUMAS RECTANGULARES MEDIANTE EL CODIGO ACI CONSIDERANDO INTERACCION MEDIANTE LA FORMULA DE C.T. Hsu*/
	/*SOLO FUNCIONA PARA ARMADOS SIMÉTRICOS*/
	/*FUNCIONA CORRECTAMENTE PARA SIMPLE Y DOBLE FELXIÓN*/
	/*BIBLIOGRAFÍA: Aspectos Fundamentales del Concreto Reforzado, González Cuevas Robles Fernandez*/

	//imprimir
	////printf("\n\nMux = %lf",Mux);
	////printf("\nMuy = %lf\n",Muy);

    double b1,pacer;
	double P0,Pnbx,Pnby,Mnbx,Mnby,Pnb,Mny,Mnx,Pn,d1,d2,M1,M2;
    double SF, SM;

	pacer=ast / (b*h);
    
// Error comentado
//	if(pacer < 0.01) {//printf("ERROR: Refuerzo en columna menor al 1%% = %e", pacer); return 0;}
//	if(pacer > 0.08) {//printf("ERROR: Refuerzo en columna mayor al 8%% = %e", pacer); return 0;}

	b1 = 1.05- fc/1400.0;
	if(b1 < 0.65) b1 = 0.65;
	if(b1 > 0.85) b1 = 0.85;
	Mnx = fabs(Mux);//falta agregar valores mínimos
	Mny = fabs(Muy);
	//
	d1=2.0;
	d2=5.0-2.0*h;
	if(d2 >d1)d1=d2;
	M1=fabs(Pu*d1);
	M2=fabs(Pu*d1);
	if (M1>Mnx) Mnx=M1;
	if (M2>Mny) Mny=M2;
	// también aplica en este caso?
	SumaF2(fc, fy, as, nbarr, b, h, y, b1, &SF, &SM);
	Pnbx = SF;
	Mnbx = fabs(SM);
	SumaF2(fc, fy, as, nbarr, h, b, x, b1, &SF, &SM);
	Pnby = SF;
	Mnby = fabs(SM);

	if ((Pu<=Pnbx)&(Pu<=Pnby)){ //flexocompresión por encima de la falla balanceada de los 2 diagramas de interacción
		P0 = - 0.8*0.65*( 0.85*fc*(b*h-ast) + ast*fy );

	}
	else{
		P0 = 0.9 * ast * fy;//flexotensión y flexocompresión por debajo de la falla balanceada de ambos diagramas

	}

	if (Mny==0){
		Pnb = Pnbx;

	}
	else if (Mnx==0){
		Pnb = Pnby;

	}
	else if (Pnbx>=Pnby){
		Pnb = Pnby;
	}
	else {
		Pnb = Pnbx;

	}

	Pn = Pnb + (P0 - Pnb)*(1 - pow((Mnx/Mnbx),1.5) - pow((Mny/Mnby),1.5));


	return Pu/Pn;
}


void SumaF2(double fc, double fy, double *as, int nbarr, double b, double h, double *y, double b1, double *SF, double *SM )
{
	//Solo funciona para la condición balanceada
	int i;
	double aux1 = 0,aux2 = 0,es,a,d,c,E = 2100000.0,phi=0.65;
	*SF = 0;
	*SM = 0;
	//Se considera que aux1 siempre va a apuntar al acero en tensión de forma que siempre sea negativa
	for (i=1;i<=nbarr;i++){
		if (y[i]<aux1) aux1 = y[i];
		if (y[i]>aux2) aux2 = y[i];
	}
	c = 0.003*(h/2-aux1)/(0.003+fy/E);
	for (i=1;i<=nbarr;i++){
        d = 0.5*h - y[i];
		if (y[i] == aux1){
			es = fy/E;
		}
		else{
			es = d * 0.003/c - 0.003;
		}
		if(fabs(es) < fy/E) {
			*SF = *SF + E*es*as[i];
			*SM = *SM + (d - 0.5*h) * E*es*as[i];
		}
		else {
			if(es > 0){
				*SF = *SF + fy*as[i];
				*SM = *SM + (d - 0.5*h) * fy*as[i];
			}
			else{
				*SF = *SF - fy*as[i];
				*SM = *SM - (d - 0.5*h) * fy*as[i];
			}
		}
	}
	a = b1 * c;
	if(a > h) a = h;
	//El valor de phi para la condición balanceada es de 0.65
	*SF = phi * (*SF - 0.85 * fc * a * b);
	*SM = phi *(*SM + 0.5*( h - a ) * ( 0.85*fc*a*b ));
}



//------------------------------------------------------------------------------


double Flexcomp_Rect_ACI_Bress(double Pu, double Mux, double Muy, double b, double h, double fc, double fy, double nbarr, double* as, double* x, double* y,double ast)
{
	/*
	 REVISION POR FLEXOCOMPRESION DE COLUMAS RECTANGULARES MEDIANTE EL CODIGO ACI CONSIDERANDO INTERACCION MEDIANTE LA FORMULA DE BRESSLER
	 DATOS:
	 Pu  = Carga axial ( Pu<0 Compresion, Pu>0 Tension )			(kg)
	 Mux = Momento ultimo alrededor del eje horizontal			(kg-cm)
	 Muy = Momento ultimo alrededor del eje vertical				(kg-cm)
	 b   = Dimension de la seccion paralela al eje horizontal	(cm)
	 h   = Dimension de la seccion paralela al eje vertical		(cm)
	 fc  = Resistencia a compresion del concreto					(kg/cm2)
	 fy  = Esfuerzo de fluencia del acero de refuerzo			(kg/cm2)
	 {as} = Areas de acero de las varillas de refuerzo			(cm2)
	 {x}  = Coordenadas en el eje horizontal de las varillas de {as} respecto al centroide de la seccion	(cm)
	 {y}  = Coordenadas en el eje vertical de las varillas de {as} respecto al centroide de la seccion	(cm)

	 SALIDA:
	 Relacion Pu / Pr  (Carga ultima / Carga Resistente) sin considerar los efectos de esbeltez
	 */

    double b1, Pnx, Pny, P0, Pn,pacer,M1,M2,d1,d2;
	pacer=ast / (b*h);
	//jacob.pdr revisar Acorde al Dr. Botello estas revisiones seguido dan error. Tal vez sea causa de los catalogos
    // Error comentado
//	if(pacer < 0.01) {//printf("ERROR: Refuerzo en columna menor al 1 porciento = %e",pacer); return 0;}
//	if(pacer > 0.08) {//printf("ERROR: Refuerzo en columna mayor al 8 porciento = %e",pacer); return 0;}

	b1 = 1.05- fc/1400.0;
	if(b1 < 0.65) b1 = 0.65;
	if(b1 > 0.85) b1 = 0.85;
	M1=fabs(Mux);
	d1=2.0;
	d2=5.0-2.0*h;
	//	//printf("d1=%lf d2=%lf\n",d1,d2);
	if(d2 >d1)d1=d2;
	M2=fabs(Pu*d1);
	////printf("M1 =%lf M2=%lf \n",M1,M2);
	if(M2 >M1) M1=M2;
	Pnx = Pnom(Pu, M1, b, h, fc, fy, y, as, nbarr, b1);
	////printf("Pnx=%lf\n",Pnx);
	M1=fabs(Muy);
	d1=2.0;
	d2=5.0-2.0*b;
	//	//printf("d1=%lf d2=%lf\n",d1,d2);
	if(d2 >d1)d1=d2;
	M2=fabs(Pu*d1);
	//	//printf("M1 =%lf M2=%lf \n",M1,M2);
	if(M2 >M1) M1=M2;
	Pny = Pnom(Pu, M1, h, b, fc, fy, x, as, nbarr, b1);
	//		//printf("Pny=%lf\n",Pny);

	//jacob.pdr
	//neto.modificado
	if (Pu<0){
		P0 = - 0.8*0.65*( 0.85*fc*(b*h-ast) + ast*fy );

	}
	else{
		P0 = 0.9 * ast * fy;

	}
	//	//printf("P0=%lf\n",P0);
	//jacob.pdr IMPORTANTE
	//neto.modificado: antes no estaban comentadas y en teoría nunca deberían ser mayores
	//if(Pnx < P0) Pnx = P0;
	//if(Pny < P0) Pny = P0;

	//antes
	//Pn =  0.70/ ( 1.0/Pnx + 1.0/Pny - 1.0/P0 );
	//neto.modificado
    //printf("Pnx=%lf, Pny=%lf, P0= %lf \n",Pnx,Pny,P0);
	Pn =  1/ ( 1.0/Pnx + 1.0/Pny - 1.0/P0 );
	//	    //printf("Pu=%lf Pn0%lf Pu/Pn=%lf\n",Pu,Pn,Pu/Pn);
	return Pu / Pn;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double Pnom(double Pu, double Mu, double b, double h, double fc, double fy, double* y, double* as, int nbarr, double b1)
{
	/*
	 Calcula la carga axial nominal para interaccion carga axial-flexion en un sentido de columnas rectangulares

	 DATOS:
	 Pu = Carga ultima aplicada a la secciÛn		(kg)
	 Mu = Momento alrededor del eje horizontal paralelo a b
	 b  = DimensiÛn horizontal de la secciÛn		(cm)
	 h  = DimensiÛn vertical de la secciÛn		(cm)
	 fc = Resistencia a compresiÛn del concreto	(kg/cm2)
	 fy = Esfuerzo de fluencia del acero			(kg/cm2)
	 {y} = Coordenadas verticales del acero de refuerzo respecto al centroide de la secciÛn		(cm)
	 {as} = area c/u de las varillas o paquetes de acero de refuerzo		(cm2)
	 nbarr = Numero de varillas de acero de refuerzo
	 b1 = Parametro que modifica la profundidad del bloque equivalente de esfuerzos  0.65 <= b1 <= 0.85
	 SALIDA:
	 regresa = Carga axial nominal resistente		(kg)
	 */

	double SF, SM;
	int kcont,kmax=1000;
	double c, c1, c2, e, e1, e2, eu;
	//jacob.pdr
	//Modificado
	double i,paso=0.05;
	int encontrado=0,valmax=1000;
	/////////

	c1 = h;						//Eje neutro supuesto
	c2 = 0.5*h;
	//jacob.pdr
	//Modificado
	eu = Mu / Pu;//se le cambia el signo para que en caso de compresion busque una excentricidad negativa y para tension una positiva

//jacob.pdr revisar dese aqui hasta excentricidades (aprox 15 lineas)
	//Modificado (obtiene valores iniciales de auxc1,auxe1,auxc2 y auxe2 para el ciclo)
	if (Pu>0){
		do{
			for (i=0 ; i<h ; i+=paso){
				SumaF(fc , fy , as , nbarr , b , h , y , i , b1, &SF, &SM );
	e1 = SM / SF;
				if (e1 < (-valmax)){
				c1 = i;
				encontrado = 1;
				break;//duda: esto permite salir de ambos ciclos?
				}
			}
			paso /= 1.2;
		}while (encontrado != 1);
	}
	else{
		do{
			for (i=0 ; i<h ; i+=paso){
				SumaF(fc , fy , as , nbarr , b , h , y , i , b1, &SF, &SM );
				e1 = SM / SF;
				if (e1 > valmax){
				c1 = i;
				encontrado = 1;
				break;
				}
			}
			paso /= 1.2;
		}while (encontrado != 1);
	}
	c2=c1-paso;
	SumaF(fc, fy, as, nbarr, b, h, y, c2, b1, &SF, &SM );
	e2 = SM / SF;
	/*Código anterior
	SumaF(fc, fy, as, nbarr, b, h, y, c1, b1);
	e1 = SM / SF;
	SumaF(fc, fy, as, nbarr, b, h, y, c2, b1);
	e2 = SM / SF;*/

	kcont=0;
	do {
		kcont++;
		//jacob.pdr
		c = c1 - (e1-eu) * (c1 - c2) / (e1 - e2);//TENER CUIDADO CON EL SIGNO PROPIO DE eu (eu - : cuando Pu es +; eu + : cuando Pu es -)
		SumaF(fc, fy, as, nbarr, b, h, y, c, b1, &SF, &SM );
		e = SM / SF;
		////printf("kcont=%d e=%lf eu=%lf \n",kcont,e,eu);
		//jacob.pdr
		if(fabs(fabs(e) - fabs(eu)) < 0.00000001) break;//modificado
		c2 = c1;
		c1 = c;
		e2 = e1;
		e1 = e;
	}while  (kcont <= kmax);
    //printf("Pu=%lf,Mu=%lf,e=%lf,eu=%lf,SM=%lf,SF=%lf  \n",Pu,Mu,e,eu,SM,SF);
	return SF;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void SumaF(double fc, double fy, double* as, int nbarr, double b, double h, double* y, double c, double b1, double *SF, double *SM )
{
	/*
	 DATOS:
	 fc = Ressitencia a la compresion del concreto				(kg/cm2)
	 fy = Esfuerzo de fluencia del acero							(kg/cm2)
	 {as} = Areas de acero de cada varilla o paquete de varillas en la secciÛn		(cm2)
	 nbarr = Numero de varillas
	 b  = DimensiÛn horizontal de la seccion de columna			(cm)
	 h  = Dimension vertical de la columna						(cm)
	 {y} = Coordenadas verticales de las varillas de refuerzo		(cm)
	 c  = Eje neutro propuesto de la seccion						(cm)
	 b1 = Parametro que modifica la profundidad del bloque equivalente de esfuerzos  0.65 <= b1 <= 0.85

	 SALIDA:
	 SF = Suma de fuerzas (Carga nominal Pn), SF < 0 --> CompresiÛn, SF > 0 --> Tension		(kg)
	 SM = Suma de momentos respecto al centroide												(kg-cm)
	 */
	int i;
	double d, es, a,E=2000000.0,phi,est;
	*SF = 0;
	*SM = 0;
	//jacob.pdr
	//neto.agregado
	est=0;
	for (i=1; i<=nbarr; i++){											//De las varillas
		d = 0.5*h - y[i];
		es = d * 0.003/c - 0.003;
		//jacob.pdr
		//neto.agregado
		if (es>est) est=es;
		if(fabs(es) < fy/E) {
			*SF = *SF + E*es*as[i];
			*SM = *SM + (d - 0.5*h) * E*es*as[i];
		}
		else {
			if(es > 0){
				*SF = *SF + fy*as[i];
				*SM = *SM + (d - 0.5*h) * fy*as[i];
			}
			else{
				*SF = *SF - fy*as[i];
				*SM = *SM - (d - 0.5*h) * fy*as[i];
			}
		}
	}
	a = b1 * c;
	if(a > h) a = h;


	//jacob.pdr al parecer el phi viene en archivos del Dr. Alejandro
	//neto.agregado
	phi = ( est - fy/E ) / ( 6 * fy/E ) + 0.65;
	if (phi < 0.65) phi = 0.65;
	if (phi > 0.90) phi = 0.90;

	//jacob.pdr comentario D.B. //OJO valor absoluto.......
	/*
	*SF = fabs(*SF - 0.85 * fc * a * b);					//Compresion en el bloque de concreto
	*SM = fabs(*SM + 0.5*( h - a) * ( 0.85*fc*a*b));
*/
	*SF = phi * (*SF - 0.85 * fc * a * b);
	*SM = phi * (*SM + 0.5*( h - a ) * ( 0.85*fc*a*b ));

}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double flexion_r(double b,double h,double fc,double fy,int nbarr, double* as,double* y,double ast)
{
	double ymin, d,ro,w,asb,a,bb,c,raiz,Mr,Beta1=0.85,fyr[100],beta1=0.85,E=2.1e06,toler=0.00001;
	int    i,ind,tipo=1,kcont,kmax=10;

	ymin=0;
        
    for (i=1; i<=nbarr; i++){
        fyr[i]=0;
		if(y[i]<ymin) ymin=y[i];
    }
	if (tipo ==0) {
		d=0.5*h-ymin;
		asb=0;
		for (i=1; i<=nbarr; i++)
			if(fabs(y[i]-ymin)<toler) asb+=as[i];
		ro=asb/(b*h);
		w=ro*fy/fc;
		Mr=0.9*b*d*d*fc*w*(1-0.59*w);
	} else {
		kcont=0;
		do{
			kcont++;
			a=0.85*fc*b;
			bb=c=0;
			for (i=1; i<=nbarr; i++){
				if(fyr[i]==fy || fyr[i]==-fy)c+=fyr[i]*as[i];
				else {
					bb+=0.003*E*as[i];
					c+=-0.003*E*as[i]*Beta1*(0.5*h-y[i]);
                    //printf("barra =%d, bb=%lf ,  c=%lf \n ",i, bb,c);
			    }
		    }
            raiz =bb*bb-4*a*c;
            if(raiz<1e-10) raiz=0.0;
	        raiz=(-bb+sqrt(raiz))/(2.*a);
//printf("a=%lf bb=%lf c=%lf raiz=%lf kcont=%d \n",a,bb,c,raiz,kcont);
			ind=0;
			for (i=1; i<=nbarr; i++){
                
				if(fyr[i]==fy || fyr[i]==-fy){}
				else {
					bb=(1.0-beta1*(0.5*h-y[i])/raiz)*0.003*E;
//printf("kcont=%d barra=%d bb=%lf fy = %lf y=%lf, raiz =%lf\n",kcont,i,bb,fy, y[i],raiz);
					if(bb> fy){fyr[i]= fy; ind++;}
					if(bb<-fy){fyr[i]=-fy; ind++;}
				}
			}
			if(ind ==0)break;
		}while  (kcont <= kmax);
		c=raiz/0.85;
		Mr=a*raiz*(c-raiz/2);
		for (i=1; i<=nbarr; i++){
			if(fyr[i]==fy || fyr[i]==-fy)Mr+=as[i]*fyr[i]*(c-(0.5*h-y[i]));
			else {
				bb=(1.0-beta1*(0.5*h-y[i])/raiz)*0.003*E;
				Mr+=as[i]*bb*(c-(0.5*h-y[i]));
			}
		}
        //printf("Mr= %lf \n",Mr);
	}
    //printf("va de regreso \n");
	return Mr;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void cfmmcon( double* mmax, double L, int ifc, int ifd, double dfc, double mfc, double mfd, double fi, double ff,
			 double mi, double mf, double Pu,double Pc, double KLsR)
{
	// ----------------------------------------------------------------- //
	// Nombre : cfmmcon Coeficiente de Flexion y Momento Maximo.         //
	// Funcion: Calcula el Coeficiente de Flexion segun la metodologia   //
	//          ACI para vigas de concreto reforzado sometida a          //
	//          fuerza concetrada y uniformemente distribuidad ademas    //
	//          obtiene el Momento Maximo para dise~no.                  //
	// Autor  : Salvador Botello (CIMAT). 2010.                          //
	// Reviso :                                                          //
	// ----------------------------------------------------------------- //

	// --------------------- DATOS DE ENTRADA ---------------------- //
	// mmax  Momento Maximo en la Barra.
	// L     Longitud de la Barra.
	// ifc   Indicador de Fuerza Concentrada.
	// ifd   Indicador de Fuerza Distribuida.
	// dfc   Distancia de Aplicacion de Fuerza Concentrada.
	// mfc   Magnitud de la Fuerza Concentrada.
	// mfd   Magnitud de la Fuerza Distribuida.
	// fi    Fuerza en Nodo Inicial.
	// ff    Fuerza en Nodo Final.
	// mi    Momento en Nodo Inicial.
	// mf    Momento en Nodo Final.
	// Pu    Carga axial maxima actuante
	// Pc    Carga critica de Euler
	// ------------------------------------------------------------- //

	// ------------- DECLARACION DE VARIABLES LOCALES -------------- //
	double m1  = 0.0; // Momento en Extremo de Barra Maximo.
	double m2  = 0.0; // Momento en Extremo de Barra Minimo.
	double m3  = 0.0; // Momento Maximo en interior de barra.
	double xvo = 0.0; // Distancia a la que el Cortante es cero.
	// ------------------------------------------------------------- //
	double coef,Cm,fact,delta;

	// ------------------ COEFICIENTE DE FLEXION ------------------- //
	// Momentos en extremos de barra.
	//jacob.modif antes:
    /*
	 if( mi < mf){
	 m1 = mi;
	 m2 = mf;
	 }
	 else{
	 m1 = mf;
	 m2 = mi;
	 }*/
    //ahora:

    if( fabs(mi) < fabs(mf)){
        // muttio fabs(
		m1 = fabs(mi);
		m2 = fabs(mf);
	}
	else{
        // muttio fabs(
		m1 = fabs(mf);
		m2 = fabs(mi);
	}



	// Coeficiente de Flexion.
	if( ifc == 0){
		if( ifd == 0){
			// No hay ni Fuerzas Concentradas ni Distribuidas.
			// Los momentos maximos estan en extremos de barra.
			// muttio momentos en cero se cambian por momentos maximos en extremo
			//Caso columnas
			//Antes
			//m3 = 0.0;
			//Despues
			 if( fabs(mi) < fabs(mf)){
                m3 = fabs(mf);
			 } else {
                m3 = fabs(mi);
			 }


		}
		else{
			// No hay Fuerzas Concentradas pero Si Distribuidas.
			// Distancia a la que el Cortante es cero.
			//muttio fabs( fi )
			xvo = - fabs(fi) / mfd;
			// Verifica que la distancia sea dentro de la barra.
			if( ( xvo > 0.0) && ( xvo < L)){
				// Momento Maximo dentro de la barra.
				//jacob.modif
                //antes:
                //m3 = mi - 0.5 * fi * fi / mfd;
                //Despues
                //muttio fabs( ) regreso a (-0.5) ahora si hay resultados simetricos
                // Se piensa debe ser negativo porque mfd es negativo, entonces la distancia
                // al centro tiene que ser positiva (Habra que checar los casos con concentradas, etc...)
                m3 =  mi - 0.5 * fabs(fi) * fabs(fi) / mfd;

            }
			else{
            // muttio momentos en cero se cambian por momentos maximos en extremo
			//Caso columnas
			//Antes
			//m3 = 0.0;
			//Despues
			 if( fabs(mi) < fabs(mf)){
                m3 = fabs(mf);
			 } else {
                m3 = fabs(mi);
			 }
			}
		}
	}
	else{
		if( ifd == 0){
			// Si hay Fuerzas Concentradas pero No Distribuidas.
			// Evalua el Momento flexionante en el
			// punto de aplicacion de la carga.
			//jacob.modif
            //antes:
            //m3 = mi + fi * dfc;
            //despues
            m3 = mi - fi * dfc;

        }
		else{
			// Si hay Fuerzas Concentradas y Distribuidas.
			// Evalua el Momento flexionante en el
			// punto de aplicacion de la carga.
            //jacob.modif
            //antes:
			//m3 = 0.5 * mfd * dfc * dfc + mfc * dfc + mi;
            //depues
            //muttio fabs(
            m3 = -0.5 * fabs(mfd) * dfc * dfc - fabs(fi) * dfc + fabs(mi);

        }
	}

	// Coeficiente de Flexion.
	if (ifc >0 || ifd >0 || m2==0.) {Cm =1;}
    else {
        //muttio fabs(
		coef =fabs(m1)/fabs(m2);
		if (coef >1 || coef < -1) coef=1/coef;
		Cm=0.6 -0.4*coef;
		if(Cm<0.4)Cm=0.4;
	}
	fact=34-12*KLsR;
	if(KLsR <40 && KLsR <fact) delta=1.0;
	else delta=Cm/(1.0-Pu/(0.75*Pc));
	if(delta <1.0) delta =1.0;
	////printf("delta =%lf\n",delta);
	// Momento Maximo en la Barra.
	*mmax = fabs( m1);
	if( *mmax < fabs( m2)){
		*mmax = fabs( m2);
	}
	if( *mmax < fabs( m3)){
		*mmax = fabs( m3);
	}
	//muttio.
	*mmax=*mmax;//*delta;
	// ------------------------------------------------------------- //

	return;
}
