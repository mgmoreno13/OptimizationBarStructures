//
//  efi_RC.h
//  new_mecap
//
//  Created by Jacob Esau Salazar Solano on 5/24/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <math.h>

#ifndef new_mecap_efi_RC_h
#define new_mecap_efi_RC_h


void Eficiencia_MarcosRC(int ielem, int nelem, int ndime, int ntipo, int* matnu, 
                         int* ntips, double** props, double* xlong,double* fuerc,
                         double** arerc, int** indfu, double** fuerb, 
                         double** fuefl, double* wcarg);

void Eficiencia_Armaduras(int ielem, int ndime, int* matnu, double** props, double* xlong, double** arear, double* fuerc);





void Fi(double l, int lmats, int ielem, double** props, double** arear, double* fuerc);





float darc( int tb,   int ifcy,   int ifcz,   int ifdy,  int ifdz,  float L,   float E,   float Fy,   float A,    float Iy,
		   float Iz,   float Sy,   float Sz,   float Qy,  float Qz,  float ty,  float tz,  float mfcy, float mfcz, float dfcy,
		   float dfcz, float mfdy, float mfdz, float fxi, float fyi, float fzi, float fxf, float fyf,  float fzf,  float mxi,
		   float myi,  float mzi,  float mxf,  float myf, float mzf);

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

void cefp( float *Fb, float rg, float L, float Cb, float Fy);




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






void cfmm( float *Cb, float *mmax, float L, int ifc, int ifd, float dfc, float mfc, float mfd, float fi, float ff,
		  float mi, float mf);


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





void frmf( float *Cm, float Fe, int ifc, int ifd, float mi, float mf, int tb, float faa);

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



#endif
