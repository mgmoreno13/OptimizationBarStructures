//
//  efi_CR.h
//  new_mecap
//
//  Created by Jacob Esau Salazar Solano on 5/24/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//
#include <math.h>

#ifndef new_mecap_efi_CR_h
#define new_mecap_efi_CR_h

void Eficiencia_MarcosCR(int ielem, int nelem, int ndime, int ntipo, int* matnu, int* ntips, double** props,
						 double* xlong,double* fuerc,double** arecr, int** indfu, double** fuerb,
						 double** fuefl, double* wcarg, double*** dtcon);

double Flexcomp_concr(int ifcy,int ifcz,int ifdy,int ifdz,double L, double b, double h, double fc, double fy, int nbarr, double* as, double* x, double* y,
					  int tb, double mfcy, double mfcz,double dfcy,double dfcz, double mfdy, double mfdz, double fxi, double fyi, double fzi, double fxf,
					  double fyf,  double fzf,  double mxi, double myi,  double mzi,  double mxf,  double myf, double mzf);


/*
 REVISION POR FLEXOCOMPRESION DE COLUMAS RECTANGULARES MEDIANTE EL CODIGO ACI CONSIDERANDO INTERACCION MEDIANTE LA FORMULA DE BRESSLER
 DATOS:
 Pu  = Carga axial ( Pu<0 Compresion, Pu>0 Tension)			(kg)
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
double Flexcomp_Rect_ACI_Bress(double Pu, double Mux, double Muy, double b, double h, double fc, double fy, double nbarr, double* as, double* x, double* y,double ast);


/*
 Calcula la carga axial nominal para interaccion carga axial-flexion en un sentido de columnas rectangulares
 
 DATOS:
 Pu = Carga ultima aplicada a la secciÛn		(kg)
 Mu = Momento alrededor del eje horizontal paralelo a b
 b  = DimensiÛn horizontal de la secciÛn		(cm)
 h  = DimensiÛn vertical de la secciÛn		(cm)
 fc = Resistencia a compresiÛn del concreto	(kg/cm≤)
 fy = Esfuerzo de fluencia del acero			(kg/cm≤)
 {y} = Coordenadas verticales del acero de refuerzo respecto al centroide de la secciÛn		(cm)
 {as} = ¡rea c/u de las varillas o paquetes de acero de refuerzo		(cm≤)
 nbarr = N˙mero de varillas de acero de refuerzo
 b1 = Par·metro que modifica la profundidad del bloque equivalente de esfuerzos  0.65 <= b1 <= 0.85
 SALIDA:
 regresa = eficiencia = Pu / Pn
 */
double Pnom(double Pu, double Mu, double b, double h, double fc, double fy, double* y, double* as, int nbarr, double b1);


/*
 DATOS:
 fc = Ressitencia a la compresion del concreto				(kg/cm≤)
 fy = Esfuerzo de fluencia del acero							(kg/cm≤)
 {as} = Areas de acero de cada varilla o paquete de varillas en la secciÛn		(cm≤)
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
void SumaF(double fc, double fy, double* as, int nbarr, double b, double h, double* y, double c, double b1,
           double *SF, double *SM );




// ----------------------------------------------------------------- //
// Nombre : SumaF2											         //
// Funcion: jacob.pdr ERNESTO : falta agregar descripción            //
// Autor  : Ernesto Ortega Trujillo. 2010.                         	 //
// Reviso :                                                          //
// ----------------------------------------------------------------- //
// --------------------- DATOS DE ENTRADA -------------------------- //
// fc = Ressitencia a la compresion del concreto			(kg/cm≤)	//
// fy = Esfuerzo de fluencia del acero						(kg/cm≤)	//
// {as} = Areas de acero de cada varilla o paquete de varillas 			//
//			en la secciÛn		(cm≤)									//
// nbarr = Numero de varillas											//
// b  = DimensiÛn horizontal de la seccion de columna			(cm)	//
// h  = Dimension vertical de la columna						(cm)	//
// {y} = Coordenadas verticales de las varillas de refuerzo		(cm)	//
// c  = Eje neutro propuesto de la seccion						(cm)	//
// b1 = Parametro que modifica la profundidad del bloque equivalente de	//
// 		esfuerzos  0.65 <= b1 <= 0.85									//
// ----------------------------------------------------------------- 	//
// --------------------- DATOS DE SALIDA --------------------------- 	//
// SF = Suma de fuerzas (Carga nominal Pn), 							//
//		SF < 0 --> CompresiÛn, SF > 0 --> Tension		(kg)			//
// SM = Suma de momentos respecto al centroide			(kg-cm)			//
// ----------------------------------------------------------------- //
void SumaF2(double fc, double fy, double *as, int nbarr, double b, double h, double *y, double b1, double *SF, double *SM );





double flexion_r(double b,double h,double fc,double fy,int nbarr, double* as,double* y,double ast);


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
void cfmmcon( double* mmax, double L, int ifc, int ifd, double dfc, double mfc, double mfd, double fi, double ff,
			 double mi, double mf, double Pu,double Pc, double KLsR);


// -----------------------------------------------------------------//
// Nombre : Flexotension_Rect_ACI_CTHsu 							//
// Funcion:	Calcula la carga axial nominal para interaccion carga 	//
//			axial-flexion en un sentido de columnas rectangulares	//
//			Cuando Carga Axial = TENSION.							//
// BIBLIOGRAFÍA: Aspectos Fundamentales del Concreto Reforzado, 	//
// 			González Cuevas Robles Fernandez						//	
// Autor  : Ernesto Ortega Trujillo (CIMAT). 2011.                  //
// Reviso :                                                         //
// -----------------------------------------------------------------//
// --------------------- DATOS DE ENTRADA --------------------------//
// Pu  = Carga axial ( Pu<0 Compresion, Pu>0 Tension)		(kg)	//
// Mux = Momento ultimo alrededor del eje horizontal		(kg-cm)	//
// Muy = Momento ultimo alrededor del eje vertical			(kg-cm)	//
// b   = Dimension de la seccion paralela al eje horizontal	(cm)	//
// h   = Dimension de la seccion paralela al eje vertical	(cm)	//
// fc  = Resistencia a compresion del concreto				(kg/cm2)//
// fy  = Esfuerzo de fluencia del acero de refuerzo			(kg/cm2)//
// {as} = Areas de acero de las varillas de refuerzo		(cm2)	//
// {x}  = Coordenadas en el eje horizontal de las varillas de {as}	// 
// 			respecto al centroide de la seccion	(cm)				//
// {y}  = Coordenadas en el eje vertical de las varillas de {as} 	//
// 			respecto al centroide de la seccion	(cm)				//
// -----------------------------------------------------------------//
// --------------------- DATOS DE SALIDA ---------------------------//
// Relacion Pu / Pr  (Carga ultima / Carga Resistente) sin 			//
//        	considerar los efectos de esbeltez						//
// -----------------------------------------------------------------//


double Flexotension_Rect_ACI_CTHsu(double Pu, double Mux, double Muy, double b, double h, double fc, double fy, double nbarr, double *as, double *x, double *y,double ast);

#endif
