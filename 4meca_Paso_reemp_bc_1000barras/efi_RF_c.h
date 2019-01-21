//
//  efi_RF_c.h
//  new_mecap
//
//  Created by Jacob Esau Salazar Solano on 5/24/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <math.h>

#ifndef new_mecap_efi_RF_c_h
#define new_mecap_efi_RF_c_h

#include "efi_RF_o.h"

void Eficiencia_MarcosRF(int ielem, int nelem, int ndime, int ntipo, int* matnu, int* ntips, double** props,
						 double* xlong,double* fuerc,double** arerf, int** indfu, double** fuerb,
						 double** fuefl, double* wcarg, int tipo_seccion);
double f2esc(double H, double W, double D, double R, double t, double Fy,
             double E, double LongX, double LongY, double LongT, int tb,
             int ifcz, int ifcy, int ifdz, int ifdy, double mfcz, double mfcy,
             double dfcz, double dfcy, double mfdz, double mfdy, double fxi,
             double fyi, double fzi, double fxf,double fyf, double fzf,
             double mxi, double myi, double mzi, double mxf, double myf,
             double mzf, int ielem);


//----------------------------------------------------------------------------//
// Este programa calcula el area efectiva de una seccion C sometida a
// compresion uniforme segun la metodologia de la AISI B2 explicado en el
// capitulo 7 del libro de Hancock
//
// Autor          : Hector Hernandez. 17 de Marzo del 2009.
// Modificaciones : 18 de Marzo del 2009, Inicio de codificacion.
// Reviso         :
//----------------------------------------------------------------------------//

// Datos de Entrada:
//     H  =  Peralte de la seccion                   Deep         cm
//     W  =  Ancho de la seccion                     Flange       cm
//     D  =  Longitud del doblez                     Lip          cm
//     R  =  Radio interno del doblez                Radius       cm
//     t  =  Espesor de la lamina                    tickness     cm
// F_eval =  Esfuerzo al cual se desea calcular el area efectiva  Kg/cm^2
//     E  =  Modulo de elasticidad de acero                       Kg/cm^2
//----------------------------------------------------------------------------//







//-----------------------------------------------------------------------------//
//                                                                             //
// Esta rutuna calcula la eficiencia de una seccion C sometida a Compresion y  //
// Flexion en dos ejes, usando el Sistema Internacional de medicion.           //
//                                                                             //
// Autor  : Hector Hernandez. 05 de Octubre del 2010.                          //
// Reviso :                                                                    //
//                                                                             //
//-----------------------------------------------------------------------------//
//                                                                             //
// Datos de Entrada:                                                           //
//   H     =  Peralte de la seccion.     Deep         cm                       //
//   W     =  Ancho de la seccion.       Flange       cm                       //
//   D     =  Longitud del doblez.       Lip          cm                       //
//   R     =  Radio interno del doblez.               cm                       //
//   t     =  Espesor de la lamina.      tickness     cm                       //
//   Fy    =  Esfuerzo de fluencia del acero.         Kg/cm^2                  //
//   E     =  Modulo de Elasticidad.                  Kg/cm^2                  //
//   LongX =  Longitud sin Arriostrar en X.           cm                       //
//   LongY =  Longitud sin Arriostrar en Y.           cm                       //
//   LongT =  Longitud sin Arriostrar a Torsion.      cm                       //
//   tb    =  Tipo de Barra 1 Art-Art 2 Emp-Emp 3 Emp-Art.                     //
//   ifcy  =  Indicador de Fuerza Concentrada en Y 1 Si 0 No.                  //
//   ifcz  =  Indicador de Fuerza Concentrada en Z 1 Si 0 No.                  //
//   ifdy  =  Indicador de Fuerza Distribuida en Y 1 Si 0 No.                  //
//   ifdz  =  Indicador de Fuerza Distribuida en Z 1 Si 0 No.                  //
//   L     =  Longitud de la barra.                                            //
//   mfcy  =  Magnitud de la Fuerza Concentrada en Y.                          //
//   mfcz  =  Magnitud de la Fuerza Concentrada en Z.                          //
//   dfcy  =  Distancia de Fuerza Concentrada en Y.                            //
//   dfcz  =  Distancia de Fuerza Concentrada en Z.                            //
//   mfdy  =  Magnitud de la Fuerza Distribuida en Y.                          //
//   mfdz  =  Magnitud de la Fuerza Distribuida en Z.                          //
//   fxi   =  Fuerza en eje X en nodo Inicial.                                 //
//   fyi   =  Fuerza en eje Y en nodo Inicial.                                 //
//   fzi   =  Fuerza en eje Z en nodo Inicial.                                 //
//   fxf   =  Fuerza en eje X en nodo Final.                                   //
//   fyf   =  Fuerza en eje Y en nodo Final.                                   //
//   fzf   =  Fuerza en eje Z en nodo Final.                                   //
//   mxi   =  Momento en eje X en nodo Inicial.                                //
//   myi   =  Momento en eje Y en nodo Inicial.                                //
//   mzi   =  Momento en eje Z en nodo Inicial.                                //
//   mxf   =  Momento en eje X en nodo Final.                                  //
//   myf   =  Momento en eje Y en nodo Final.                                  //
//   mzf   =  Momento en eje Z en nodo Final.                                  //
//                                                                             //
//-----------------------------------------------------------------------------//
//                                                                             //
// Datos de Salida:                                                            //
//   EFI  = Eficiencia de la Seccion.                                          //
//                                                                             //
//-----------------------------------------------------------------------------//
//                                                                             //
// Funciones Llamadas:                                                         //
//                                                                             //
//  Area Efectiva a Compresion de una Seccion C.                               //
//   [ A_efe_n ] = aecsc(H, W, D, R, t, Fn, E);                                //
//                                                                             //
//  Momentos Nominales Flexionantes de una Seccion C.                          //
//   [ Mnx , Mny1 , Mny2 ] = mnfsc( H, W, D, R, t, E, Fcx, Fcy1, Fcy2 ,ielem); //
//                                                                             //
// Obtiene Acciones Mecanicas Actuantes.                                       //
//   [ K, P, Vx, Vy, Mx, My, Cb, Ctf ] = oama( tb, ifcy, ifcz, ifdy, ifdz, ... //
//                                                                             //
//-----------------------------------------------------------------------------//


double aecsc( double H, double W, double D, double R, double t, double F_eval, double E);


//-----------------------------------------------------------------------------//
//                                                                             //
// Esta rutuna calcula los momentos nomimales de una seccion C de acero rolado //
// en frio que estara sujeto a flexion en ambos ejes X e Y, usando el Sistema  //
// Internacional de medicion.                                                  //
//                                                                             //
// Autor  : Hector Hernandez. 05 de Marzo del 2009.                            //
// Cambios: 09/03/2009 Inicio de codificacion.                                 //
//          11/03/2009 Se Agrego calculo de propiedades con flexion en eje X.  //
//          16/03/2009 Se Agrego calculo de propiedades con flexion en eje Y.  //
//          26/09/2010 Se Agrego calculo de propiedades con flexion en eje Y   //
//                     con giro de 180 grados (por lo general mas debil).      //
// Reviso :                                                                    //
//                                                                             //
//-----------------------------------------------------------------------------//
//                                                                             //
// Datos de Entrada:                                                           //
//   H    =  Peralte de la seccion.     Deep         cm                        //
//   W    =  Ancho de la seccion.       Flange       cm                        //
//   D    =  Longitud del doblez.       Lip          cm                        //
//   R    =  Radio interno del doblez.               cm                        //
//   t    =  Espesor de la lamina.      tickness     cm                        //
//   E    =  Modulo de Elasticidad.                  Kg/cm^2                   //
//   Fcx  =  Esfuerzo de Evaluacion en X.            Kg/cm^2                   //
//   Fcy1 =  Esfuerzo de Evaluacion en Y1.           Kg/cm^2                   //
//   Fcy2 =  Esfuerzo de Evaluacion en Y2.           Kg/cm^2                   //
//                                                                             //
//-----------------------------------------------------------------------------//
//                                                                             //
// Datos de Salida:                                                            //
//   Mnx  =  Momento Nominal en X.                   Kg*cm                     //
//   Mny1 =  Momento Nominal en Y1.                  Kg*cm                     //
//   Mny2 =  Momento Nominal en Y2.                  Kg*cm                     //
//                                                                             //
//-----------------------------------------------------------------------------//

void mnfsc( double H, double W, double D, double R, double t, double E, double Fcx,
            double Fcy1, double Fcy2, double* Mnx, double* Mny1, double* Mny2, int ielem );


//-----------------------------------------------------------------------------//
//                                                                             //
// Esta rutuna obtiene las acciones mecanicas maximas y sus correspondientes   //
// coeficientes en una barra de acero rolado en frio para su posterior dise~no //
//                                                                             //
// Autor  : Hector Hernandez. 22 de Noviembre del 2010.                        //
// Reviso :                                                                    //
//                                                                             //
//-----------------------------------------------------------------------------//
//                                                                             //
// Datos de Entrada:                                                           //
//    tb   =  Tipo de Barra 1 Art-Art 2 Emp-Emp 3 Emp-Art.                     //
//    ifcy =  Indicador de Fuerza Concentrada en Y 1 Si 0 No.                  //
//    ifcz =  Indicador de Fuerza Concentrada en Z 1 Si 0 No.                  //
//    ifdy =  Indicador de Fuerza Distribuida en Y 1 Si 0 No.                  //
//    ifdz =  Indicador de Fuerza Distribuida en Z 1 Si 0 No.                  //
//    L    =  Longitud de la barra.                                            //
//    mfcy =  Magnitud de la Fuerza Concentrada en Y.                          //
//    mfcz =  Magnitud de la Fuerza Concentrada en Z.                          //
//    dfcy =  Distancia de Fuerza Concentrada en Y.                            //
//    dfcz =  Distancia de Fuerza Concentrada en Z.                            //
//    mfdy =  Magnitud de la Fuerza Distribuida en Y.                          //
//    mfdz =  Magnitud de la Fuerza Distribuida en Z.                          //
//    fxi  =  Fuerza en eje X en nodo Inicial.                                 //
//    fyi  =  Fuerza en eje Y en nodo Inicial.                                 //
//    fzi  =  Fuerza en eje Z en nodo Inicial.                                 //
//    fxf  =  Fuerza en eje X en nodo Final.                                   //
//    fyf  =  Fuerza en eje Y en nodo Final.                                   //
//    fzf  =  Fuerza en eje Z en nodo Final.                                   //
//    mxi  =  Momento en eje X en nodo Inicial.                                //
//    myi  =  Momento en eje Y en nodo Inicial.                                //
//    mzi  =  Momento en eje Z en nodo Inicial.                                //
//    mxf  =  Momento en eje X en nodo Final.                                  //
//    myf  =  Momento en eje Y en nodo Final.                                  //
//    mzf  =  Momento en eje Z en nodo Final.                                  //
//-----------------------------------------------------------------------------//
//                                                                             //
// Datos de Salida:                                                            //
//    K    =  Factor de Longitud Efectiva.                                     //
//    P    =  Fuerza Axial Maxima.                                             //
//    Vx   =  Fuerza Cortante en X.                                            //
//    Vy   =  Fuerza Cortante en Y.                                            //
//    Mx   =  Momento Flexionante en X.                                        //
//    My   =  Momento Flexionante en Y.                                        //
//    Cb   =  Coeficiente de Flexion en X.                                     //
//    Ctf  =  Coeficiente de Flexion en Y.                                      //
//                                                                             //
//-----------------------------------------------------------------------------//
void oama( double* K, double* P, double* Vx, double* Vy, double* Mx, double* My, double* Cb,
          double* Ctf, int tb, int ifcy, int ifcz, int ifdy, int ifdz, double L, double mfcy,
          double mfcz, double dfcy, double dfcz, double mfdy, double mfdz, double fxi, double fyi,
          double fzi, double fxf, double fyf, double fzf, double mxi, double myi, double mzi,
          double mxf, double myf, double mzf, int ielem, double A, double Fy );

double mom_pandeo_distorsional_flex( double E, double Fy, double mu, double Ix,
                                    double t, double ho, double bo, double D );

double esf_pandeo_distorsional_flex( double E, double mu, double t, double ho,
                                     double bo, double D );

double fuer_pandeo_distorsional_comp( double E, double Fy, double mu, double Ag,
                                      double t, double ho, double bo, double D );

double esf_pandeo_distorsional_comp( double E, double mu, double t, double ho,
                                     double bo, double D );


#endif
