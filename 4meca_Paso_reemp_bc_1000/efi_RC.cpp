//
//  efi_RC.cpp
//  new_mecap
//
//  Created by Jacob Esau Salazar Solano on 5/24/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <stdio.h>

#include "efi_RC.h"

void Eficiencia_MarcosRC(int ielem, 	int nelem, 	int ndime, 	int ntipo,  int* matnu, 	int* ntips, double** props,		
                         double* xlong,	double* fuerc,	double** arerc, 
                         int** indfu, 	double** fuerb,	double** fuefl, 
                         double* wcarg	)
{
	int    lmats;
	double cero=0;
	if(ndime== 2){
		lmats=(int)props[matnu[ielem]][6]+0.5;
		if(ntipo ==1)
			fuerc[ielem]=darc(ntips[ielem],cero,
							  indfu[2][ielem],cero,
							  indfu[4][ielem],xlong[ielem],
							  arerc[lmats][1],arerc[lmats][2],arerc[lmats][3],
							  arerc[lmats][4],arerc[lmats][5],arerc[lmats][7],
							  arerc[lmats][8],arerc[lmats][9],arerc[lmats][10],
							  arerc[lmats][11],arerc[lmats][12],cero,
							  fuerb[2][ielem],cero,fuerb[4][ielem],cero,
							  fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],fuefl[ielem][2],cero,
							  fuefl[ielem][3],fuefl[ielem][4],cero,cero,cero,cero,
							  cero,cero,cero);
        else
			fuerc[ielem]=darc(ntips[ielem],cero,
			  			      indfu[2][ielem],cero,
						      indfu[4][ielem],xlong[ielem],
						      arerc[lmats][1],arerc[lmats][2],arerc[lmats][3],
						      arerc[lmats][4],arerc[lmats][5],arerc[lmats][7],
						      arerc[lmats][8],arerc[lmats][9],arerc[lmats][10],
						      arerc[lmats][11],arerc[lmats][12],cero,
						      fuerb[2][ielem],cero,fuerb[4][ielem],cero,
						      fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],cero,
						      fuefl[ielem][2],fuefl[ielem][4],cero,
						      fuefl[ielem][5],cero,fuefl[ielem][3],cero,cero,
						      fuefl[ielem][6],cero);
	}
	else if(ndime==3) {
		lmats=(int)props[matnu[ielem]][9]+0.5;
        if(ntipo ==1)
			fuerc[ielem]= darc(ntips[ielem],indfu[1][ielem],indfu[2][ielem],
							   indfu[3][ielem],indfu[4][ielem],xlong[ielem],
							   arerc[lmats][1],arerc[lmats][2],
							   arerc[lmats][3],arerc[lmats][4],
							   arerc[lmats][5],arerc[lmats][7],
							   arerc[lmats][8],arerc[lmats][9],
							   arerc[lmats][10],arerc[lmats][11],
							   arerc[lmats][12],cero,
							   fuerb[2][ielem],cero,
							   fuerb[4][ielem],fuerb[5][ielem]+wcarg[ielem+nelem],
							   fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],
							   0.0,fuefl[ielem][2],
							   fuefl[ielem][4],0.0,
							   fuefl[ielem][5],cero,fuefl[ielem][3],cero,cero,fuefl[ielem][6],cero);
		else
			fuerc[ielem]= darc(ntips[ielem],indfu[1][ielem],indfu[2][ielem],
                               indfu[3][ielem],indfu[4][ielem],xlong[ielem],
                               arerc[lmats][1],arerc[lmats][2],
                               arerc[lmats][3],arerc[lmats][4],
                               arerc[lmats][5],arerc[lmats][7],
                               arerc[lmats][8],arerc[lmats][9],
                               arerc[lmats][10],arerc[lmats][11],
                               arerc[lmats][12],fuerb[1][ielem],
                               fuerb[2][ielem],fuerb[3][ielem],
                               fuerb[4][ielem],fuerb[5][ielem]+wcarg[ielem+nelem],
                               fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],
                               fuefl[ielem][2],fuefl[ielem][3],
                               fuefl[ielem][7],fuefl[ielem][8],
                               fuefl[ielem][9],fuefl[ielem][4],
                               fuefl[ielem][5],fuefl[ielem][6],
                               fuefl[ielem][10],fuefl[ielem][11],
                               fuefl[ielem][12]);
	}
}



//-------INICIA EFICIENCIAS ARMADURAS-------------------------------------------
//------------------------------------------------------------------------------
void Eficiencia_Armaduras(int ielem, int ndime, int* matnu, double** props, double* xlong, double** arear, double* fuerc)
{
	int    lmats;
    // Seccion del catalogo
	if(ndime == 2)lmats=(int)props[matnu[ielem]][6]+0.5;
	if(ndime == 3)lmats=(int)props[matnu[ielem]][9]+0.5;
	Fi(xlong[ielem],lmats,ielem,props,arear,fuerc);
    
}

//------------------------------------------------------------------------------
void Fi(double l, int lmats, int ielem, double** props, double** arear, double* fuerc)
/*
 calcula restriccion de esfuerzos
 */
{
	double rel, Esbeltez,xmin, rg, Cc,Fy,E,a,pi=M_PI,esfac,esfre ;
	int    icomp=1;
	double sigmamax=3500.;
    
	E = props[1][1];
	a = arear[lmats][1];
	esfac = fuerc[ielem] / a ;
	if (icomp == 0) esfac =fabs(esfac) ;
	if (esfac >= 0.0) esfre=sigmamax;
	else{
		Fy=sigmamax;
		esfac=fabs(esfac) ;
		xmin=arear[lmats][2];
		if(xmin > arear[lmats][3]) xmin = arear[lmats][3];
		rg = sqrt(xmin/a);
		Esbeltez = l/rg;
		Cc = pi*sqrt(2.0*E/Fy);
		rel = Esbeltez/Cc;
		if(Esbeltez <= Cc)
			esfre = Fy*(1.0-0.5*rel*rel)/(5.0/3.0+0.375*rel-0.125*rel*rel*rel);
		else
			esfre = (12.0*pi*pi*E)/(23.0*Esbeltez*Esbeltez);
	}
	fuerc[ielem]=esfac/esfre;
	//fprintf(fp16, "%d\t %lf\t %lf \t %lf \t %lf \n", ielem, esfac, esfre,fuerc[ielem],a);
}
//----FIN EFICIENCIAS ARMADURAS-------------------------------------------------
//------------------------------------------------------------------------------




// Rutinas de rolado en caliente
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

float darc( int tb,		int ifcy,   int ifcz,   int ifdy,  	int ifdz,  	float L,   
           float E,   	float Fy,   float A,    float Iy, 	float Iz,  	float Sy,   
           float Sz,   float Qy,  	float Qz,  	float ty,  	float tz,  	float mfcy,
           float mfcz, float dfcy,	float dfcz, float mfdy, float mfdz, float fxi, 
           float fyi, 	float fzi, 	float fxf, 	float fyf,  float fzf,  float mxi,
           float myi,  float mzi,  float mxf,  float myf, 	float mzf)
{
	// ----------------------------------------------------------------- //
	// Nombre : DARC Dise~no de Acero Rolado en Caliente.                //
	// Funcion: Calcula la eficiencia de una seccion de acero rolada en  //
	//          caliente y sometida a flexion en dos ejes y fuerza axial //
	//          segun la metodologia de la AISC.                         //
	// Autor  : Hector Hernandez (CIMAT). 2010.                          //
	// Reviso :                                                          //
	// ----------------------------------------------------------------- //
	// ----------------- CONSIDERACIONES GENERALES ----------------- //
	// La fuerza de gravedad actua en el eje Z Global.
	// Eje X Local esta alineado con la seccion de forma longitudinal.
	// Eje Y Local coincide con el eje fuerte de la seccion.
	// Eje Z Local coincide con el eje debil de la seccion.
	// ------------------------------------------------------------- //
    
	// --------------------- DATOS DE ENTRADA ---------------------- //
	// tb     Tipo de Barra 1 Art-Art 2 Emp-Emp 3 Emp-Art.
	// ifcy   Indicador de Fuerza Concentrada en Y 1 Si 0 No.
	// ifcz   Indicador de Fuerza Concentrada en Z 1 Si 0 No.
	// ifdy   Indicador de Fuerza Distribuida en Y 1 Si 0 No.
	// ifdz   Indicador de Fuerza Distribuida en Z 1 Si 0 No.
	// L      Longitud de la barra.
	// E      Modulo de Elasticidad.
	// Fy     Modulo de Fluencia.
	// A      Area de la seccion.
	// Iy     Momento de Inercia en Y.
	// Iz     Momento de Inercia en Z.
	// Sy     Modulo de Seccion en Y.
	// Sz     Modulo de Seccion en Z.
	// Qy     Momento de Area en Y.
	// Qz     Momento de Area en Z.
	// ty     Espesor resistente a Cortante en Y.
	// tz     Espesor resistente a Cortante en Z.
	// mfcy   Magnitud de la Fuerza Concentrada en Y.
	// mfcz   Magnitud de la Fuerza Concentrada en Z.
	// dfcy   Distancia de Fuerza Concentrada en Y.
	// dfcz   Distancia de Fuerza Concentrada en Z.
	// mfdy   Magnitud de la Fuerza Distribuida en Y.
	// mfdz   Magnitud de la Fuerza Distribuida en Z.
	// fxi    Fuerza en eje X en nodo Inicial.
	// fyi    Fuerza en eje Y en nodo Inicial.
	// fzi    Fuerza en eje Z en nodo Inicial.
	// fxf    Fuerza en eje X en nodo Final.
	// fyf    Fuerza en eje Y en nodo Final.
	// fzf    Fuerza en eje Z en nodo Final.
	// mxi    Momento en eje X en nodo Inicial.
	// myi    Momento en eje Y en nodo Inicial.
	// mzi    Momento en eje Z en nodo Inicial.
	// mxf    Momento en eje X en nodo Final.
	// myf    Momento en eje Y en nodo Final.
	// mzf    Momento en eje Z en nodo Final.
	// ------------------------------------------------------------- //
    
	// ------------- DECLARACION DE VARIABLES LOCALES -------------- //
	float PICU  = 0.0; // PI al CUadrado.
	float tefl  = 0.0; // TEmporal Flotante.
	float K     = 0.0; // Factor de Longitud Efectiva.
	float fam   = 0.0; // Fuerza Axial Maxima.
	float rgy   = 0.0; // Radio de Giro en Y.
	float rgz   = 0.0; // Radio de Giro en Z.
	float rgmin = 0.0; // Radio de Giro MINimo.
	float faa   = 0.0; // Esfuerzo Axial Actuante.
	float Cby   = 0.0; // Coeficiente de Flexion en Y.
	float Cbz   = 0.0; // Coeficiente de Flexion en Z.
	float Mmaxy = 0.0; // Momento Maximo en Y.
	float Mmaxz = 0.0; // Momento Maximo en Z.
	float Fbpy  = 0.0; // Esfuerzo a Flexion Permisible en Y.
	float Fbpz  = 0.0; // Esfuerzo a Flexion Permisible en Z.
	float fbay  = 0.0; // Esfuerzo a Flexion Actuante en Y.
	float fbaz  = 0.0; // Esfuerzo a Flexion Actuante en Z.
	float Ccc   = 0.0; // Carga Critica a Compresion.
	float Fa    = 0.0; // Esfuerzo Axial Permisible.
	float Fey   = 0.0; // Esfuerzo Critico de Euler en eje Y.
	float Fez   = 0.0; // Esfuerzo Critico de Euler en eje Z.
	float Cmy   = 0.0; // Factor de Reduccion de Momentos en Eje Y.
	float Cmz   = 0.0; // Factor de Reduccion de Momentos en Eje Z.
	float Vy    = 0.0; // Fuerza de Corte en eje Y.
	float Vz    = 0.0; // Fuerza de Corte en eje Z.
	float Vca   = 0.0; // Esfuerzo Cortante Actuante.
	float efi   = 0.0; // EFIciencia de la seccion.
	float efco  = 0.0; // EFiciencia a COrtante.
	// ------------------------------------------------------------- //
    
	/*	 fprintf(fp16,"nuevo elemento\n");
	 fprintf(fp16,"tb = %d\n",tb);
	 fprintf(fp16,"ifcy =%d\n",ifcy);
	 fprintf(fp16,"ifcz =%d\n",ifcz);
	 fprintf(fp16,"ifdy =%d\n",ifdy);
	 fprintf(fp16,"ifdz =%d\n",ifdz);
	 fprintf(fp16,"L =%f\n",L);
	 fprintf(fp16,"E =%f\n",E);
	 fprintf(fp16,"Fy =%f\n",Fy);
	 fprintf(fp16,"A =%f\n",A);
	 fprintf(fp16,"Iy =%f\n",Iy);
	 fprintf(fp16,"Iz =%f\n",Iz);
	 fprintf(fp16,"Sy =%f\n",Sy);
	 fprintf(fp16,"Sz =%f\n",Sz);
	 fprintf(fp16,"Qy =%f\n",Qy);
	 fprintf(fp16,"Qz =%f\n",Qz);
	 fprintf(fp16,"ty =%f\n",ty);
	 fprintf(fp16,"tz =%f\n",tz);
	 fprintf(fp16,"mfcy =%f\n",mfcy);
	 fprintf(fp16,"mfcz =%f\n",mfcz);
	 fprintf(fp16,"dfcy =%f\n",dfcy);
	 fprintf(fp16,"dfcz =%f\n",dfcz);
	 fprintf(fp16,"mfdy =%f\n",mfdy);
	 fprintf(fp16,"mfdz =%f\n",mfdz);
	 fprintf(fp16,"fxi =%f\n",fxi);
	 fprintf(fp16,"fyi =%f\n",fyi);
	 fprintf(fp16,"fzi =%f\n",fzi);
	 fprintf(fp16,"fxf =%f\n",fxf);
	 fprintf(fp16,"fyf =%f\n",fyf);
	 fprintf(fp16,"fzf =%f\n",fzf);
	 fprintf(fp16,"mxi =%f\n",mxi);
	 fprintf(fp16,"myi =%f\n",myi);
	 fprintf(fp16,"mzi =%f\n",mzi);
	 fprintf(fp16,"mxf =%f\n",mxf);
	 fprintf(fp16,"myf =%f\n",myf);
	 fprintf(fp16,"mzf =%f\n",mzf);
	 */
	// Declaracion de Constantes.
	PICU = 9.869604401089358;
    
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
        case 4 : // Art-Emp.
			K = 0.80;
			break;
            
	}
    
	// Obtiene Fuerza Axial Maxima.
	fam = fabs( fxi);
	if( fam < fabs( fxf)){
		fam = fabs( fxf);
	}
    
	// Radios de Giro y su Minimo.
	rgy = sqrt( Iy / A);
	rgz = sqrt( Iz / A);
	rgmin = rgy;
	if( rgmin > rgz){
		rgmin = rgz;
	}
    
	// Esfuerzo Axial Actuante.
	faa = fam / A;
    
	// Coeficiente de Flexion y Momento Maximo en Y.
	cfmm( &Cby, &Mmaxy, L, ifcz, ifdz, dfcz, mfcz, mfdz, fzi, fzf, myi, myf);
    
	// Coeficiente de Flexion y Momento Maximo en Z.
	cfmm( &Cbz, &Mmaxz, L, ifcy, ifdy, dfcy, mfcy, mfdy, fyi, fyf, mzi, mzf);
    
	// Esfuerzo a Flexion Permisible en Y.
	cefp( &Fbpy, rgy, L, Cby, Fy);
    
	// Esfuerzo a Flexion Permisible en Z.
	cefp( &Fbpz, rgz, L, Cbz, Fy);
    
	// Esfuerzo por Flexion actuante en Y.
	fbay = Mmaxy / Sy;
    
	// Esfuerzo por Flexion actuante en Z.
	fbaz = Mmaxz / Sz;
    
	// Verifica si la barra esta sometida a Tension o Compresion.
	if( fxi > 0.0){
		// ---------------------- COMPRESION ----------------------- //
		// Carga Critica a Compresion.
		Ccc = sqrt( 2.0 * PICU * E / Fy);
        
		// Factor de Esbeltez.
		tefl = K * L / rgmin;
        
		// Esfuerzo Axial Permisible.
		if( tefl <= Ccc){
			Fa = ( ( 1.0 - ( ( tefl * tefl) / ( 2.0 * Ccc * Ccc))) / ( 1.666666666666667 + 0.375 * tefl / Ccc - 0.125 * tefl / ( Ccc * Ccc * Ccc))) * Fy;
		}
		else{
			Fa = 0.5217391304347826 * ( ( PICU * E) / ( tefl * tefl));
		}
        
		// Evalua la EFIciencia de la seccion sometida
		// a Compresion y Flexion en los dos ejes.
		if( ( faa / Fa) < 0.15){
			// Evalua la EFIciencia de la seccion.
			efi = ( faa / Fa) + ( fbay / Fbpy) + ( fbaz / Fbpz);
		}
		else{
			// Esfuerzo Critico de Euler en eje Y.
			Fey = PICU * E / ( ( K * L / rgy) * ( K * L / rgy));
            
			// Esfuerzo Critico de Euler en eje Z.
			Fez = PICU * E / ( ( K * L / rgz) * ( K * L / rgz));
            
			// Factor de Reduccion de Momentos en Eje Y.
			frmf( &Cmy, Fey, ifcz, ifdz, myi, myf, tb, faa);
            
			// Factor de Reduccion de Momentos en Eje Z.
			frmf( &Cmz, Fez, ifcy, ifdy, mzi, mzf, tb, faa);
            
			// Evalua la EFIciencia de la seccion.
			efi = ( faa / ( 0.6 * Fy)) + ( fbay / Fbpy) + ( fbaz / Fbpz);
			tefl = ( faa / Fa) + ( Cmy * fbay / ( ( 1.0 - ( faa / Fey)) * Fbpy)) + ( Cmz * fbaz / ( ( 1.0 - ( faa / Fez)) * Fbpz));
			if( efi < tefl){
				efi = tefl;
			}
		}
	}
	else{
		// ----------------------- TENSION ------------------------- //
		// Evalua la EFIciencia de la seccion.
		efi = ( faa / ( 0.6 * Fy)) + ( fbay / Fbpy) + ( fbaz / Fbpz);
	}
	// ------------------------------------------------------------- //
    
	// ------------------- VERIFICA POR CORTANTE ------------------- //
	// Obtiene el maximo cortante en Y.
	Vy = fabs( fyi);
	if( Vy < fabs( fyf)){
		Vy = fabs( fyf);
	}
    
	// Obtiene el maximo cortante en Z.
	Vz = fabs( fzi);
	if( Vz < fabs( fzf)){
		Vz = fabs( fzf);
	}
    
	// Esfuerzo Cortante Actuante.
	Vca = ( ( Vz * Qy) / ( ty * Iy)) + ( ( Vy * Qz) / ( tz * Iz));
    
	// EFiciencia a COrtante.
	efco = Vca / ( 0.4 * Fy);
	// ------------------------------------------------------------- //
    
	// ---------------- VERIFICA MAXIMA EFICIENCIA ----------------- //
	if( efi < efco){
		efi = efco;
	}
    // ------------------------------------------------------------- //
	return efi;
    
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void cefp( float *Fb, float rg, float L, float Cb, float Fy)
{
	// ----------------------------------------------------------------- //
	// Nombre : CEPF Calculo de Esfuerzo Permisible a Flexion.           //
	// Funcion: Calcula el Esfuerzo Permisible a Flexion para una viga   //
	//          de acero rolado en caliente segun la metodologia AISC.   //
	// Autor  : Hector Hernandez (CIMAT). 2010.                          //
	// Reviso :                                                          //
	// ----------------------------------------------------------------- //
	// --------------------- DATOS DE ENTRADA ---------------------- //
	// Fb  Esfuerzo Permisible a Flexion.
	// rg  Radio de Giro.
	// L   Longitud.
	// Cb  Coeficiente de Flexion.
	// Fy  Modulo de Fluencia.
	// ------------------------------------------------------------- //
    
	// ------------- DECLARACION DE VARIABLES LOCALES -------------- //
	float esb  = 0.0; // ESBeltez de la barra.
	float Lim1 = 0.0; // LIMites de esbeltez Inferior.
	float Lim2 = 0.0; // LIMites de esbeltez Superior.
	// ------------------------------------------------------------- //
    
	// ESBeltez de la barra.
	esb = L / rg;
    
	// Limites de esbeltez.
	Lim1 = sqrt( (  7170000.0 * Cb) / Fy);
	Lim2 = sqrt( ( 35900000.0 * Cb) / Fy);
    
	// Esfuerzo permisible a flexion (pag 153 McCormack).
	*Fb = 0.66 * Fy;
	if( ( esb < Lim2) && ( esb > Lim1)){
		*Fb = ( 0.66 - Fy * esb * esb / ( 108000000.0 * Cb)) * Fy;
	}
	if( esb > Lim2){
		*Fb = 12000000.0 * Cb / ( esb * esb);
	}
    
	return;
    
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void cfmm(	float *Cb, float *mmax, float L, int ifc, int ifd, float dfc, 
          float mfc, float mfd, float fi, float ff, float mi, float mf)
{
	// ----------------------------------------------------------------- //
	// Nombre : CFMM Coeficiente de Flexion y Momento Maximo.            //
	// Funcion: Calcula el Coeficiente de Flexion segun la metodologia   //
	//          AISC para vigas de acero rolado en caliente sometida a   //
	//          fuerza concetrada y uniformemente distribuidad ademas    //
	//          obtiene el Momento Maximo para dise~no.                  //
	// Autor  : Hector Hernandez (CIMAT). 2010.                          //
	// Reviso :                                                          //
	// ----------------------------------------------------------------- //
    
	// --------------------- DATOS DE ENTRADA ---------------------- //
	// Cb    Coeficiente de Flexion.
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
	// ------------------------------------------------------------- //
    
	// ------------- DECLARACION DE VARIABLES LOCALES -------------- //
	float m1  = 0.0; // Momento en Extremo de Barra Maximo.
	float m2  = 0.0; // Momento en Extremo de Barra Minimo.
	float m3  = 0.0; // Momento Maximo en interior de barra.
	float xvo = 0.0; // Distancia a la que el Cortante es cero.
	// ------------------------------------------------------------- //
    
	// ------------------ COEFICIENTE DE FLEXION ------------------- //
	// Momentos en extremos de barra.
	if( mi < mf){
		m1 = mi;
		m2 = mf;
	}
	else{
		m1 = mf;
		m2 = mi;
	}
    
	// Coeficiente de Flexion.
	if( ifc == 0){
		if( ifd == 0){
			// No hay ni Fuerzas Concentradas ni Distribuidas.
			// Los momentos maximos estan en extremos de barra.
			m3 = 0.0;
		}
		else{
			// No hay Fuerzas Concentradas pero Si Distribuidas.
			// Distancia a la que el Cortante es cero.
			xvo = - fi / mfd;
			// Verifica que la distancia sea dentro de la barra.
			if( ( xvo > 0.0) && ( xvo < L)){
				// Momento Maximo dentro de la barra.
				//jacob.modif
                //antes:
                //m3 = mi - 0.5 * fi * fi / mfd;
                //despues
                m3 = mi + 0.5 * fi * fi / mfd;
                
            }
			else{
				m3 = 0.0;
			}
		}
	}
	else{
		if( ifd == 0){
			// Si hay Fuerzas Concentradas pero No Distribuidas.
			// Evalua el Momento flexionante en el
			// punto de aplicacion de la carga.
			//jacob.modif
            //antes
            //m3 = mi + fi * dfc;
            //antes:
            m3 = mi - fi * dfc;
            
        }
		else{
			// Si hay Fuerzas Concentradas y Distribuidas.
			// Evalua el Momento flexionante en el
			// punto de aplicacion de la carga.
			//jacob.modif
            //antes
            //m3 = 0.5 * mfd * dfc * dfc + mfc * dfc + mi;
            //despues
            m3 = -0.5 * mfd * dfc * dfc - fi * dfc + mi;
            
        }
	}
    
	// El Momento Flexionante Maximo se encuentra dentro de la Barra.
	if( ( fabs( m3) > fabs( m1)) && ( fabs( m3) > fabs( m2))){
		m1 = 1.0;
		m2 = 1.0;
	}
    
	// Coeficiente de Flexion.
	if( ( fabs( m1) < 0.00000000001) && ( fabs( m2) < 0.00000000001)){
		*Cb = 1.75;
	}
	else{
		*Cb = 1.75 + 1.05 * ( m1 / m2) + 0.3 * ( m1 / m2) * ( m1 / m2);
	}
    
	if( *Cb > 2.3){
		*Cb = 2.3;
	}
    
	// Momento Maximo en la Barra.
	*mmax = fabs( m1);
	if( *mmax < fabs( m2)){
		*mmax = fabs( m2);
	}
	if( *mmax < fabs( m3)){
		*mmax = fabs( m3);
	}
	// ------------------------------------------------------------- //
    
	return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void frmf( float *Cm, float Fe, int ifc, int ifd, float mi, float mf, int tb, float faa)
{
	// ----------------------------------------------------------------- //
	// Nombre : FRMF Factor de Reduccion de Momento Flexionante.         //
	// Funcion: Calcula el Factor de Reduccion de Momentos Flexionantes  //
	//          para una viga de acero rolado en caliente segun la meto- //
	//          dologia AISC.                                            //
	// Autor  : Hector Hernandez (CIMAT). 2010.                          //
	// Reviso :                                                          //
	// ----------------------------------------------------------------- //
    
    
	// --------------------- DATOS DE ENTRADA ---------------------- //
	// Cm   Factor de Reduccion de Momentos.
	// Fe   Esfuerzo de Euler.
	// ifc  Indicador de Fuerza Concentrada.
	// ifd  Indicador de Fuerza Distribuida.
	// mi   Momento en Nodo Inicial.
	// mf   Momento en Nodo Final.
	// tb   Tipo de Barra.
	// faa  Fuerza Axial Actuante.
	// ------------------------------------------------------------- //
    
	float psi = 0.0;
    
	// Factor de Reduccion de Momentos.
	if( ( ifc == 0) && ( ifd == 0)){
		// Momentos Maximos en los extremos de barra
		// no hay cargas transversales y los apoyos
		// pueden desplazarse.
		if( ( fabs( mi) < 1.0e-16) && ( fabs( mf) < 1.0e-16)){
			*Cm = 1.0;
		}
		else{
			if( mi <= mf){
				*Cm = 0.6 + 0.4 * fabs( mi / mf);
			}
			else{
				*Cm = 0.6 + 0.4 * fabs( mf / mi);
			}
		}
	}
	else{
		// Hay cargas transverslaes entre apoyos.
		switch( tb){
			case 1 : // Art-Art.
				if( ifc == 1){
					psi = 0.2;
				}
				if( ifd == 1){
					psi = 0.0;
				}
				break;
			case 2 : // Emp-Emp.
				if( ifc == 1){
					psi = 0.6;
				}
				if( ifd == 1){
					psi = 0.4;
				}
				break;
			case 3 : // Emp-Art.
				if( ifc == 1){
					psi = 0.4;
				}
				if( ifd == 1){
					psi = 0.3;
				}
				break;
            case 4 : // Art-Emp.
				if( ifc == 1){
					psi = 0.4;
				}
				if( ifd == 1){
					psi = 0.3;
				}
				break;
		}
		*Cm = 1.0 - psi * ( faa / Fe);
	}
    
	return;
    
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Fin de Rutinas de rolado en caliente