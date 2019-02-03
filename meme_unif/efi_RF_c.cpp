//
//  efi_RF_c.cpp
//  new_mecap
//
//  Created by Jacob Esau Salazar Solano on 5/24/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "efi_RF_c.h"
#include <stdio.h>


void Eficiencia_MarcosRF(int ielem, int nelem, int ndime, int ntipo, int* matnu, int* ntips, double** props,
						 double* xlong,double* fuerc,double** arerf, int** indfu, double** fuerb,
						 double** fuefl, double* wcarg, int tipo_seccion) //tipo_seccion=  1-->seccion C , 4-->seccion Omega
{
	int    lmats;
	double cero = 0;
	double Fy = 3515.33, E = 2074044.42, poisson = 0.3, LongX, LongY, LongT;
	LongX=LongY=LongT=xlong[ielem];

    if (tipo_seccion==1) { //Para secciones C
        if(ndime== 2){
            lmats = (int)props[matnu[ielem]][6]+0.5;
            if(ntipo ==1)
                fuerc[ielem]=f2esc(arerf[lmats][10],arerf[lmats][11],
                                   arerf[lmats][12],arerf[lmats][13],
                                   arerf[lmats][8],
                                   Fy,E,LongX,LongY, LongT,ntips[ielem],
                                   cero,indfu[2][ielem],cero, indfu[4][ielem],
                                   cero,fuerb[2][ielem],cero, fuerb[4][ielem],
                                   cero,fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],fuefl[ielem][2],cero,
                                   fuefl[ielem][3],fuefl[ielem][4],cero,cero,cero,cero,
                                   cero,cero,cero,ielem);
            else
                fuerc[ielem]=f2esc(arerf[lmats][10],arerf[lmats][11],
                                   arerf[lmats][12],arerf[lmats][13],
                                   arerf[lmats][8],
                                   Fy,E,LongX,LongY, LongT,ntips[ielem],
                                   cero,indfu[2][ielem],cero, indfu[4][ielem],
                                   cero,fuerb[2][ielem],cero, fuerb[4][ielem],
                                   cero,fuerb[6][ielem]+wcarg[ielem],
                                   fuefl[ielem][1],cero,fuefl[ielem][2],
                                   fuefl[ielem][4],cero,fuefl[ielem][5],
                                   cero, fuefl[ielem][3],cero,
                                   cero, fuefl[ielem][6],cero,ielem);
        }
        else if(ndime ==3){
            lmats=(int)props[matnu[ielem]][9]+0.5;
            if(ntipo ==1)
                fuerc[ielem]=f2esc(arerf[lmats][10],arerf[lmats][11],
                                   arerf[lmats][12],arerf[lmats][13],
                                   arerf[lmats][8],Fy,E,LongX,LongY, LongT,ntips[ielem],
                                   cero,indfu[2][ielem],cero, indfu[4][ielem],
                                   cero,fuerb[2][ielem],cero, fuerb[4][ielem],
                                   fuerb[5][ielem]+wcarg[ielem+nelem],
                                   fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],fuefl[ielem][2],
                                   fuefl[ielem][3],fuefl[ielem][4],fuefl[ielem][5],
                                   fuefl[ielem][6],cero,cero,cero,cero,cero,cero,ielem);
            else
                fuerc[ielem]=f2esc(arerf[lmats][10],arerf[lmats][11],
                                   arerf[lmats][12],arerf[lmats][13],
                                   arerf[lmats][8],Fy,E,LongX,LongY, LongT,
                                   ntips[ielem], indfu[1][ielem],
                                   indfu[2][ielem],indfu[3][ielem],
                                   indfu[4][ielem],fuerb[1][ielem],
                                   fuerb[2][ielem],fuerb[3][ielem],
                                   fuerb[4][ielem],
                                   fuerb[5][ielem]+wcarg[ielem+nelem],
                                   fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],
                                   fuefl[ielem][2],fuefl[ielem][3],
                                   fuefl[ielem][7],fuefl[ielem][8],
                                   fuefl[ielem][9],fuefl[ielem][4],
                                   fuefl[ielem][5],fuefl[ielem][6],
                                   fuefl[ielem][10],fuefl[ielem][11],
                                   fuefl[ielem][12],ielem);
        }
    }
    else {   //Para secciones Omega
        if(ndime== 2){
            lmats = (int)props[matnu[ielem]][6]+0.5;
            if(ntipo ==1)
                fuerc[ielem]=f2eso(arerf[lmats][10],arerf[lmats][11],
                                   arerf[lmats][12],arerf[lmats][13],
                                   arerf[lmats][8],
                                   Fy,E,LongX,LongY, LongT,ntips[ielem],
                                   cero,indfu[2][ielem],cero, indfu[4][ielem],
                                   cero,fuerb[2][ielem],cero, fuerb[4][ielem],
                                   cero,fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],fuefl[ielem][2],cero,
                                   fuefl[ielem][3],fuefl[ielem][4],cero,cero,cero,cero,
                                   cero,cero,cero,ielem);
            else
                fuerc[ielem]=f2eso(arerf[lmats][10],arerf[lmats][11],
                                   arerf[lmats][12],arerf[lmats][13],
                                   arerf[lmats][8],
                                   Fy,E,LongX,LongY, LongT,ntips[ielem],
                                   cero,indfu[2][ielem],cero, indfu[4][ielem],
                                   cero,fuerb[2][ielem],cero, fuerb[4][ielem],
                                   cero,fuerb[6][ielem]+wcarg[ielem],
                                   fuefl[ielem][1],cero,fuefl[ielem][2],
                                   fuefl[ielem][4],cero,fuefl[ielem][5],
                                   cero, fuefl[ielem][3],cero,
                                   cero, fuefl[ielem][6],cero,ielem);
        }
        else if(ndime ==3){
            lmats=(int)props[matnu[ielem]][9]+0.5;
            if(ntipo ==1)
                fuerc[ielem]=f2eso(arerf[lmats][10],arerf[lmats][11],
                                   arerf[lmats][12],arerf[lmats][13],
                                   arerf[lmats][8],Fy,E,LongX,LongY, LongT,ntips[ielem],
                                   cero,indfu[2][ielem],cero, indfu[4][ielem],
                                   cero,fuerb[2][ielem],cero, fuerb[4][ielem],
                                   fuerb[5][ielem]+wcarg[ielem+nelem],
                                   fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],fuefl[ielem][2],
                                   fuefl[ielem][3],fuefl[ielem][4],fuefl[ielem][5],
                                   fuefl[ielem][6],cero,cero,cero,cero,cero,cero,ielem);
            else
                fuerc[ielem]=f2eso(arerf[lmats][10],arerf[lmats][11],
                                   arerf[lmats][12],arerf[lmats][13],
                                   arerf[lmats][8],Fy,E,LongX,LongY, LongT,
                                   ntips[ielem], indfu[1][ielem],
                                   indfu[2][ielem],indfu[3][ielem],
                                   indfu[4][ielem],fuerb[1][ielem],
                                   fuerb[2][ielem],fuerb[3][ielem],
                                   fuerb[4][ielem],
                                   fuerb[5][ielem]+wcarg[ielem+nelem],
                                   fuerb[6][ielem]+wcarg[ielem],fuefl[ielem][1],
                                   fuefl[ielem][2],fuefl[ielem][3],
                                   fuefl[ielem][7],fuefl[ielem][8],
                                   fuefl[ielem][9],fuefl[ielem][4],
                                   fuefl[ielem][5],fuefl[ielem][6],
                                   fuefl[ielem][10],fuefl[ielem][11],
                                   fuefl[ielem][12],ielem);
        }
    }

}
// Rutinas de rolado en frio
//------------------------------------------------------------------------------


double f2esc(double H, double W, double D, double R, double t, double Fy,
             double E, double LongX, double LongY, double LongT, int tb,
             int ifcz, int ifcy, int ifdz, int ifdy, double mfcz, double mfcy,
             double dfcz, double dfcy, double mfdz, double mfdy, double fxi,
             double fyi, double fzi, double fxf,double fyf, double fzf,
             double mxi, double myi, double mzi, double mxf, double myf,
             double mzf, int ielem)
{
    //jacob.modif provisional
    /*se agrego la entrada del numero del elemento (ielem) para revision de rutina.
	 la modificación se realizo en 4 puntos a) esta funcion, b) en declaracion de esta funcion al
	 inicio del archivo, c) en las dos ocasiones en que es llamada
     */

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

	//-------------------------------------------------------------------------//
	//                           VARIABLES LOCALES                             //
	//-------------------------------------------------------------------------//
	double DOVO;       // Un DOceaVO = 1.0 / 12.0.
	double SEXT;       // Un SEXTo   = 1.0 / 06.0.
	double TERC;       // Un TERCio  = 1.0 / 03.0.
	double PICU;       // Pi al cuadrado.
	double tempo;      // Variable TEMPOral.
	double r;          // Radio medio del cuarto de rondana.
	double u;          // Longitud de arco de la rondana.
	double c;          // Distancia del centroidal del cuarto de rondana.
	double h;          // Longitud plana del Web o Peralte.
	double w;          // Longitud plana del Flange o Ancho.
	double d;          // Longitud plana del Lip o Labio.
	double A;          // Area de la seccion completa.
	double Ix;         // Momento de inercia en X.
	double Iy;         // Momento de inercia en Y.
	double J;          // Momento polar de inercia.
	double xbar;       // Distancia del lomo a la coordenada centroidal.
	double xcg;        // Abscisa del centroide de la seccion.
	double xco;        // Abscisa del centro de cortante de la seccion.
	double xo;         // Distancia de los centros de cortante al de gravedad.
	double aa, bb, cc; // Constantes para evaluar el alabeo.
	double Ms, Cw;     // Constantes de alabeo.
	double rx, ry, ro; // Radios de giro.
	double betaw;      // Constante de evaluacion de propiedad j.
	double betaf;      // Constante de evaluacion de propiedad j.
	double beta1;      // Constante de evaluacion de propiedad j.
	double betay;      // Constante de evaluacion de propiedad j.
	double jtofl;      // Propiedad j de torsion.
	double SigmaX;     // Esfuerzo critico de Euler para pandeo en el eje X.
	double SigmaY;     // Esfuerzo critico de Euler para pandeo en el eje Y.
	double SigmaT;     // Esfuerzo critico de Euler para pandeo torsional.
	double FE;         // Esfuerzo critico de pandeo lateral de torsion.
	double Fn;         // Esfuerzo critico a compresion.
	double Fc;         // Esfuerzo inelastico de pandeo lateral de torsion.
	double Fcx;        // Esfuerzo inelastico de pandeo lateral de torsion en X.
	double Fcy1;       // Esfuerzo inelastico de pandeo lateral de torsion en Y1.
	double Fcy2;       // Esfuerzo inelastico de pandeo lateral de torsion en Y2.
	double Lambda;     // Factor de esbeltez.
	double A_efe_n;    // Area Efectiva con el esfuerzo critico Fn.
	double K;          //  Factor de Longitud Efectiva.
	double P;          //  Fuerza Axial Maxima.
	double Vx;         //  Fuerza Cortante en X.
	double Vy;         //  Fuerza Cortante en Y.
	double Mx;         //  Momento Flexionante en X.
	double My;         //  Momento Flexionante en Y.
	double Cb;         //  Coeficiente de Flexion en X.
	double Ctf;        //  Coeficiente de Flexion en Y.
	double Pn;         // Resistencia Axial Nominal.
	double PnD;        // Resistencia Nominal a la Compresion por Pandeo Distorsional
	double Sf;         // Modulo de Seccion.
	double Cs;         // Coeficiente de esfuerzo critico de pandeo.
	double Mnx;        // Momento Nominal en X.
	double Mny;        // Momento Nominal en Y.
    double Mnxt;	   // Momento Nominal de la sección completa = Sftx * Fy
    double Mnyt;	   // Momento Nominal de la sección completa = Sfty * Fy
	double Mny1;       // Momento Nominal en Y1.
	double Mny2;       // Momento Nominal en Y2.
    double MnD;        // Momento Nominal por pandeo distorsional ( solo aplica
                       // cuando la flexion es sobre el eje fuerte)
	double Vny;        // Fuerza Cortante Resistente Nominal sobre el eje Y
	double Vnx;        // Fuerza Cortante Resistente Nominal sobre el eje X
	double omega_b;    // Factor de seguridad a Flexion.
	double omega_c;    // Factor de seguridad a Compresion.
	double omega_v;    // Factor de seguridad a Cortante.
    double omega_t;    // Factor de seguridad a Tensión.
	double ax, ay;     // Factores de esbeltes.
	double Pex, Pey;   // Factores de pandeo.
	double EFIFC;      // EFIciencia de la seccion a Flexion y Compresion.
    double EFIFT;      // EFIciencia de la seccion a Flexion y Tensión.
    double EFIFT2;	   // EFIciencia de la seccion a Flexion y Tensión auxiliar.
	double EFIFV;      // EFIciencia de la seccion a Flexion y Cortante.
	double EFI;        // EFIciencia de la seccion.

	//-------------------------------------------------------------------------//


    /*
	 fprintf(fp16,"nuevo elemento\n");
	 fprintf(fp16,"H = %lf\n",H);
	 fprintf(fp16,"W = %lf\n",W);
	 fprintf(fp16,"D = %lf\n",D);
	 fprintf(fp16,"R = %lf\n",R);
	 fprintf(fp16,"t = %lf\n",t);
	 fprintf(fp16,"Fy =%f\n",Fy);
	 fprintf(fp16,"E =%f\n",E);
	 fprintf(fp16,"LongX =%f\n",LongX);
	 fprintf(fp16,"LongY =%f\n",LongY);
	 fprintf(fp16,"LongT =%f\n",LongT);

	 fprintf(fp16,"tb = %d\n",tb);
	 fprintf(fp16,"ifcy =%d\n",ifcy);
	 fprintf(fp16,"ifcz =%d\n",ifcz);
	 fprintf(fp16,"ifdy =%d\n",ifdy);
	 fprintf(fp16,"ifdz =%d\n",ifdz);
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

	 fprintf(fp16,"L =%f\n",LongT);
     */

    /*
    //Ejemplo 4.2 del Yu. Seccion C sometida a flexion
    H  = 10.0;
    W  = 3.5;
    D  = 0.72;
    R  = 3.0 / 32.0;
    t  = 0.075;
    Fy = 50.0;
    E  = 29500.0;
    LongX = 5 * 12.0;
    LongY = 5 * 12.0;
    LongT = 5 * 12.0;
    tb = 1;
    ifcz = 0;
    ifcy = 0;
    ifdz = 0;
    ifdy = 0;
    fxi  = 3.48;
    fyi  = 0.0;
    fzi  = 0.0;
    mxi  = 0.0;
    myi  = 0.0;
    mzi  = 3.48 * 2.116;
    fxf  = -3.48;
    fyf  = 0.0;
    fzf  = 0.0;
    mxf  = 0.0;
    myf  = 0.0;
    mzf  = -3.48 * 2.116;
    //*/

    /*
    //Ejemplo 6-3 de la pag 399 del Yu (Flexo-compresion)
    H  = 8.0;
    W  = 3.0;
    D  = 0.81;
    R  = 0.1875;
    t  = 0.105;
    Fy = 50.0;
    E  = 29500.0;
    LongX = 5.0 * 12.0;
    LongY = 5.0 * 12.0;
    LongT = 5.0 * 12.0;
    tb = 1;
    ifcz = 0;
    ifcy = 0;
    ifdz = 0;
    ifdy = 0;
    fxi  = 3.48;
    fyi  = 0.0;
    fzi  = 0.0;
    mxi  = 0.0;
    myi  = 0.0;
    mzi  = 3.48 * 2.116;
    fxf  = -3.48;
    fyf  = 0.0;
    fzf  = 0.0;
    mxf  = 0.0;
    myf  = 0.0;
    mzf  = -3.48 * 2.116;
    //*/

    /*
    // Ejemplo de la pag 165 del Yu (Flexion)
    // Seccion omega que se puede simular como una seccion C
    H  = 15.0;
    W  = 10.0;
    D  = 1.34;
    R  = 0.3;
    t  = 0.2;
    Fy = 50.0;
    E  = 29500.0;
    LongX = 100.0;
    LongY = 100.0;
    LongT = 100.0;
    tb = 1;
    ifcz = 0;
    ifcy = 0;
    ifdz = 0;
    ifdy = 0;
    fxi  = 3.48;
    fyi  = 0.0;
    fzi  = 0.0;
    mxi  = 0.0;
    myi  = 0.0;
    mzi  = 3.48 * 2.116;
    fxf  = -3.48;
    fyf  = 0.0;
    fzf  = 0.0;
    mxf  = 0.0;
    myf  = 0.0;
    mzf  = -3.48 * 2.116;
    //*/

    /* Ejemplo 5.6 del Hancock pag. 158
    H  = 8.0;
    W  = 2.75;
    D  = 0.625;
    R  = 0.1875;
    t  = 0.060;
    Fy = 55.0;
    E  = 29500.0;
    LongX = 5.0 * 12.0;
    LongY = 5.0 * 12.0;
    LongT = 5.0 * 12.0;
    tb = 1;
    ifcz = 0;
    ifcy = 0;
    ifdz = 0;
    ifdy = 0;
    fxi  =   2.0;
    fyi  =  -0.1;
    fzi  =   0.2;
    mxi  =   0.0;
    myi  = -12.3;
    mzi  =  -5.7;
    fxf  =  -2.0;
    fyf  =   0.1;
    fzf  =  -0.2;
    mxf  =   0.0;
    myf  =   0.3;
    mzf  =  -0.3;
    tb = 1; // TEMPORAL: SOLO PARA PRUEBAS
    //*/

	//-------------------------------------------------------------------------//
	//                    DECLARACION DE ALGUNAS CONSTANTES                    //
	//-------------------------------------------------------------------------//
	// Un DOceaVO.
	DOVO = 1.0 / 12.0;
	// Un SEXTo.
	SEXT = 1.0 / 6.0;
	// Un TERCio.
	TERC = 1.0 / 3.0;
	// PI al CUadrado.
	PICU = 9.86960440108936;
	// Variable TEMPOral.
	tempo = 0.0;
	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//             PROPIEDADES GEOMETRICAS DE LA SECCION COMPLETA              //
	//                  CALCULADAS POR EL METODO DE LINEAS                     //
	//-------------------------------------------------------------------------//
	// Propiedades de Cuarto de Rondana.
	r = R + 0.5 * t;   // Radio medio.
	u = 1.570 * r;     // Longitud de arco.
	c = 0.637 * r;     // Centroide.



	// Longitudes planas de la seccion.
	h = H - 2.0 * ( R + t);   // Peralte.
	w = W - 2.0 * ( R + t);   // Flange.
	d = D - 1.0 * ( R + t);   // Labio.

	// Area de la seccion completa.
	A = t * ( 2.0 * w + h + 4.0 * u + 2.0 * d);

	// Momento de inercia en X.
	Ix = 2.0 * t * ( 2.0 * u * ( 0.5 * h + c) * ( 0.5 * h + c) +
					w * ( 0.5 * H - 0.5 * t) * ( 0.5 * H - 0.5 * t) + d *
					( 0.5 * H - R - t - 0.5 * d) * ( 0.5 * H - R - t - 0.5 * d)) +
	DOVO * h * h * h * t;

	// Abscisa del centroide de la seccion medido
	// a partir de la parte interna del peralte.
	xcg = ( 0.5 * t * h + 2.0 * u * ( t + R - c) +
		   2.0 * u * ( W - R - t + c) +
		   w * W + 2.0 * d * ( W - 0.5 * t)) * t / A;

	// Momento de inercia en Y.
	Iy = ( 0.25 * t * t * h + 2.0 * u * ( t + R - c) * ( t + R - c) +
		  2.0 * u * ( W - R - t + c) * ( W - R - t + c) + 0.5 * W * W * w +
		  2.0 * d * ( W - 0.5 * t) * ( W - 0.5 * t)) * t +
	SEXT * w * w * w * t - A * xcg * xcg;

	// Momento polar de inercia.
	// jacob.comentario constante torsional de St. Venant.
    // jacob.comentario p 196 Yu
    J = TERC * A * t * t;

	// Centro de cortante y constantes para evaluar el alabeo.
	//jacob.comentario p160 Hancock
    aa = H - t;
	bb = W - t;
	cc = D - 0.5 * t;
	Ms = ( bb * ( 3.0 * aa * aa * bb + cc * ( 6.0 * aa * aa - 8.0 * cc * cc))) /
	( aa * aa * aa + 6.0 * aa * aa * bb +
	 cc * ( 8.0 * cc * cc - 12.0 * aa * cc + 6.0 * aa * aa));

	// Abscisa del centro de cortante de la seccion completa.
	xco = - ( Ms - 0.5 * t) - xcg;

	// Constante de alabeo.
	Cw = DOVO * aa * aa * bb * bb * t * ( ( 2.0 * aa * aa * aa * bb +
										   003.0 * aa * aa * bb * bb + 48.0 * cc * cc * cc * cc +
										   112.0 * bb * cc * cc * cc + 8.0 * aa * cc * cc * cc +
										   48.0 * aa * bb * cc * cc + 12.0 * aa * aa * cc * cc +
										   12.0 * aa * aa * bb * cc + 6.0 * aa * aa * aa * cc) /
										 ( 6.0 * aa * aa * bb + ( aa + 2.0 * cc) * ( aa + 2.0 * cc) *
										  ( aa + 2.0 * cc) - 24.0 * aa * cc * cc));

	// Radios de giro.
	// jacob.comentario p161 Hancock
    rx = sqrt( Ix / A);
	ry = sqrt( Iy / A);
	ro = sqrt( ( Ix + Iy) / A + xco * xco);


	// Distancia del lomo a la coordenada centroidal.
    // jacob.comentario en realidad al parecer es del eje del alma al centro de
    // gravedad, considerando dimensiones como en
    // pag 159 delHancock fig 5.20 (b) ; formula en pag. 645 Appendix C lipped channel.
	xbar = bb * ( bb + 2.0 * cc) / ( aa + 2.0 * bb + 2.0 * cc);

	// Distancia de los centros de cortante al de gravedad.
	xo = - Ms - xbar;

	// Constantes de evaluacion de propiedad j.
    // jacob.comentario Xbar debe ser considerada negativa por
    // eso signos negativos.
    // Appendix C del Yu
	betaw = - DOVO * t * xbar * aa * aa * aa - t * xbar * xbar * xbar * aa;
	tempo = bb - xbar;
	betaf = 0.5 * t * ( tempo * tempo * tempo * tempo -
					   xbar * xbar * xbar * xbar) + 0.25 * aa * aa * t *
	( tempo * tempo - xbar * xbar);
	beta1 = 2.0 * cc * t * tempo * tempo * tempo + (2.0/3.0) * t *
	tempo * ( ( aa * aa * aa / 8.0) - ( ( aa / 2.0) - cc) *
			 ( ( aa / 2.0) - cc) * ( ( aa / 2.0) - cc));
	betay = ( ( betaw + betaf + beta1) / Iy) - 2.0 * xo;

	// Propiedad de Pandeo de Torsion-Flexion.
	// Ec. 6.47 pag 375 del Yu
    jtofl = 0.5 * betay;
	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//      OBTIENE ACCIONES MECANICAS ACTUANTES Y COEFICIENTES DE APOYO       //
	//-------------------------------------------------------------------------//

	oama( &K, &P, &Vx, &Vy, &Mx, &My, &Cb, &Ctf, tb, ifcy, ifcz, ifdy, ifdz, LongT, mfcy, mfcz, dfcy, dfcz, mfdy, mfdz, fxi, fyi, fzi,
		 fxf, fyf, fzf, mxi, myi, mzi, mxf, myf, mzf, ielem, A, Fy );
	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//        ESFUERZOS CRITICOS DE EULER PARA EVALUAR PANDEO Y TORSION        //
	//-------------------------------------------------------------------------//
	// jacob.comentario pag 208 del Yu
    // Esfuerzo critico de Euler para pandeo en el eje X.
	SigmaX = PICU * E / ( ( K * LongX / rx) * ( K * LongX / rx));
	// Esfuerzo critico de Euler para pandeo en el eje Y.
	SigmaY = PICU * E / ( ( K * LongY / ry) * ( K * LongY / ry));
	// Esfuerzo critico de Euler para pandeo torsional.

    //jacob.comentario  formula factorizando E de a cuerdo a ec. 4.69 pag 208 del yu.
    SigmaT = ( E / ( A * ro * ro)) *
	( 0.38461538461538 * J + PICU * Cw / ( K * K * LongT * LongT));

	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//                       CARGA NOMINAL A COMPRESION                        //
	//-------------------------------------------------------------------------//

    // Esfuerzo critico dominante.
    // jacob.comentario pag 336 del yu
	FE = ( (SigmaX + SigmaT) - sqrt( (SigmaX + SigmaT) * (SigmaX + SigmaT) -
									4.0 * ( 1.0 - ( xco / ro) * ( xco / ro)) * SigmaX * SigmaT)) /
	( 2.0 * ( 1.0 - ( xco / ro) * ( xco / ro)));
	if( FE > SigmaY){
		FE = SigmaY;
    }

	// Factor de esbeltez.
	Lambda = sqrt( Fy / FE);

	// Esfuerzo critico.
	if( Lambda <= 1.5){
		Fn = pow(0.658,( Lambda * Lambda)) * Fy;
	}
	else{
		Fn = ( 0.877 / ( Lambda * Lambda)) * Fy;
    }

	// Area Efectiva con el esfuerzo critico Fn.
	A_efe_n = aecsc(H, W, D, R, t, Fn, E);

	// Carga Nominal a Compresion.
	Pn = A_efe_n * Fn;

	// Resistencia Nominal a la Compresion por Pandeo Distorsional
	PnD = fuer_pandeo_distorsional_comp( E, Fy, 0.3, A, t, H, W, D );

	// Escoge la menor resistencia a la compresion
	if( Pn > PnD ) Pn = PnD;

	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//          DETERMINACION DEL ESFUERZO CRITICO DE TORSION LATERAL          //
	//               CUANDO LA FLEXION ACTUA EN EL EJE X (Fcx)                 //
	//-------------------------------------------------------------------------//
	//                                    w                                    //
	//                               <--------->                               //
	//                                                                         //
	//                              /-----------\              Ejes para el    //
	//                             / ----------- \             analisis de la  //
	//                            / /           \ \            seccion. No     //
	//                       ^   | |             | |  ^        tienen nada que //
	//                       |   | |             | |  |        con los ejes    //
	//                       |   | |             | |  | d      locales del elemento //
	//                       |   | |             | |  |             Y (eje debil)   //
	//                       |   | |             | |  v             ^          //
	//                       |   | |                                |          //
	//                       |   | |                                |          //
	//                     h |   | |                                |          //
	//                       |   | |                                -------> X (eje fuerte)//
	//                       |   | |                           centroide       //
	//                       |   | |             | |  ^                        //
	//                       |   | |             | |  |                        //
	//                       |   | |             | |  | d                      //
	//                       |   | |             | |  |                        //
	//                       v   | |             | |  v                        //
	//                            \ \           / /                            //
	//                             \ ----------- /                             //
	//                              \-----------/                              //
	//                                                                         //
	// Inciso a) de articula C3.1.2.1 del AISI. Pagina 207 de Yu.              //
	// Cuando la flexion esta presente en el eje "X" mostado arriba, la        //
	// seccion resiste mejor la flexion que cuando es sobre el otro eje, por   //
	// esto se le conoce como eje fuerte. Cuando la flexion es sobre el eje "Y"//
	// la seccion no es tan resistente, por eso lo de eje debil. Notar que     //
	// cuando la flexion es sobre el eje debil, hay dos casos distintos.       //
	//-------------------------------------------------------------------------//
	// Modulo de seccion elatico de la seccion completa.
	Sf = Ix / ( H / 2.0);

	// Esfuerzo critico de pandeo lateral de torsion.
	FE = ( Cb * ro * A / Sf) * sqrt( SigmaY * SigmaT);
	//jacob.impresion provisional
    //rif(ielem==157||ielem==176){printf("Cb=%f   ro=%f  A=%f   \n",Cb, ro, A);}


	// Esfuerzo inelastico de pandeo lateral de torsion en eje X.
	Fc = 1.111111111111111 * Fy * ( 1.0 - 0.2777777777777778 * ( Fy / FE));

    //jacob.impresion prov
    //if(ielem==157||ielem==176){printf("Fy=%f, FE=%f, Fc=%f\n",Fy,FE,Fc);}

    if( FE >= 2.78 * Fy){
		Fc = Fy;
	}
	if( FE <= 0.56 * Fy){
		Fc = FE;
	}
	Fcx = Fc;


	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//          DETERMINACION DEL ESFUERZO CRITICO DE TORSION LATERAL          //
	//               CUANDO LA FLEXION ACTUA EN EL EJE Y (Fcy1)                //
	//-------------------------------------------------------------------------//
	//                                                                         //
	//                              d                                          //
	//                         <--------->                                     //
	//                                                                         //
	//                        /-----------       -----------\                  //
	//                       / -----------       ----------- \                 //
	//                      / /                             \ \                //
	//                 ^   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//               w |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 v   | |                               | |               //
	//                      \ \                              / /               //
	//                       \ ------------------------------ /                //
	//                        \------------------------------/                 //
	//                                                                         //
	//                         <---------------------------->                  //
	//                                       h                                 //
	//                                                                         //
	// Inciso b) de articula C3.1.2.1 del AISI. Pagina 207 de Yu.              //
	//-------------------------------------------------------------------------//
	// Modulo de seccion elatico de la seccion completa.
	Sf = Iy / ( W - xcg);
	Cs = -1.0;


	// Esfuerzo critico de pandeo lateral de torsion.
	FE = ( ( Cs * A * SigmaX) / ( Ctf * Sf)) * ( jtofl + Cs *
                                                sqrt( jtofl * jtofl + ro * ro * ( SigmaT / SigmaX)));


	//jacob.impresion provisional
    /*printf("Cs=%3f  A=%3f  SigmaX=%3f  Ctf=%3f   Sf=%3f     jtofl=%f\n",Cs,A,SigmaX,Ctf,Sf,jtofl);
	 printf("raiz= %f\n",sqrt( jtofl * jtofl + ro * ro * ( SigmaT / SigmaX)));
	 printf("Fc=%f     Fy=%f     FE=%f  \n",Fc,Fy,FE);
	 */

	// Esfuerzo inelastico de pandeo lateral de torsion en eje X.
	Fc = 1.111111111111111 * Fy * ( 1.0 - 0.2777777777777778 * ( Fy / FE));
	if( FE >= 2.78 * Fy){
		Fc = Fy;
	}
	if( FE <= 0.56 * Fy){
		Fc = FE;
	}
	Fcy1 = Fc;
	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//          DETERMINACION DEL ESFUERZO CRITICO DE TORSION LATERAL          //
	//               CUANDO LA FLEXION ACTUA EN EL EJE Y (Fcy2)                //
	//-------------------------------------------------------------------------//
	//                                                                         //
	//                                       h                                 //
	//                         <---------------------------->                  //
	//                                                                         //
	//                        /-----------------------------\                  //
	//                       / ----------------------------- \                 //
	//                      / /                             \ \                //
	//                 ^   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//               w |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 v   | |                               | |               //
	//                      \ \                              / /               //
	//                       \ -----------        ----------- /                //
	//                        \-----------        -----------/                 //
	//                                                                         //
	//                         <--------->                                     //
	//                              d                                          //
	//                                                                         //
	// Inciso b) de articula C3.1.2.1 del AISI. Pagina 207 de Yu.              //
	//-------------------------------------------------------------------------//
	// Modulo de seccion elatico de la seccion completa.
	Sf = Iy / xcg;
	Cs = 1.0;

	// Esfuerzo critico de pandeo lateral de torsion.
	FE = ( ( Cs * A * SigmaX) / ( Ctf * Sf ) ) * ( jtofl + Cs *
                                       sqrt( jtofl * jtofl + ro * ro * ( SigmaT / SigmaX)));

	// Esfuerzo inelastico de pandeo lateral de torsion en eje X.
	Fc = 1.111111111111111 * Fy * ( 1.0 - 0.2777777777777778 * ( Fy / FE));
	if( FE >= 2.78 * Fy){
		Fc = Fy;
	}
	if( FE <= 0.56 * Fy){
		Fc = FE;
	}
	Fcy2 = Fc;
	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//      OBTIENE MOMENTOS PERMISIBLES EN AMBOS EJES CON LOS ESFUERZOS       //
	//      CRITICOS DE PANDEO LATERAL DE TORSION ACTUANDO EN AMBOS EJES.      //
	//      TAMBIEN SE TOMA EN CUENTA EL PANDEO DISTORSIONAL.                  //
	//-------------------------------------------------------------------------//
	// Funcion Momentos Nominales a Flexion de Seccion C.
	mnfsc( H, W, D, R, t, E, Fcx, Fcy1, Fcy2, &Mnx, &Mny1, &Mny2, ielem );
	MnD = mom_pandeo_distorsional_flex( E, Fy, 0.3, Ix, t, H, W, D );

    // Menor momento permisible sobre eje X, tomando en cuenta el pandeo
    // distorsional (eje fuerte)
    if( Mnx > MnD ) Mnx = MnD;

	// Obtiene el menor de los momentos permisibles en Y (eje debil).
	Mny = Mny1;
	if( Mny > Mny2){
		Mny = Mny2;
	}

	//Mny = Mny2; // TEMPORAL, SOLO PARA PRUEBAS

	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//                   OBTIENE FUERZA CORTANTE PERMISIBLE                    //
	//-------------------------------------------------------------------------//
	// Resistencia Nominal a Cortante cuando hay combinacion con flexion
	// segun capitulo 9 de STRUCTURAL STEEL DESIGNER'S HANDBOOK.
	//pag 170 del Hancock
	// Sobre el alma (web):
    Vny = 0.60 * t * t * sqrt( 5.34 * Fy * E);
	if( (h/t) < sqrt( 5.34 * E / Fy)){
		Vny = 0.60 * Fy * h * t;
	}
	if( (h/t) > ( 1.51 * sqrt( 5.34 * E / Fy))){
		Vny = 0.904 * 5.34 * E * t * t * t / h;
	}
    // Sobre los lomos (flanges)
    Vnx = 2.0 * 0.60 * t * t * sqrt( 5.34 * Fy * E);
	if( (w/t) < sqrt( 5.34 * E / Fy)){
		Vnx = 2.0 * 0.60* Fy * w * t;
	}
	if( (w/t) > ( 1.51 * sqrt( 5.34 * E / Fy))){
		Vnx = 2.0 * 0.904 * 5.34 * E * t * t * t / w;
	}
	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//                        EFICIENCIA DE LA SECCION                         //
	//-------------------------------------------------------------------------//
	// Factores de seguridad ASD.
	omega_c = 1.80;   // Compresion.
	omega_b = 1.67;   // Flexion.
	omega_v = 1.67;   // Corte.
    omega_t = 1.67;	  // Tensión.





    // Pag. 364  del YU.
    // Si P = tensión,
    if ( fxi < 0.0) {
    	//  Ix / ( H / 2.0 )= Sftx 	;
        //  Iy / ( W - xcg )= Sfty1 ; Se utilizará esta por dar el menor Stfy
        //  Iy / xcg 		= Sfty2 ;
        //  Sftx * Fy 		= Mnxt 	;
        //  Sfty * Fy		= Mnyt 	;

        Mnxt = Ix / ( H / 2.0 ) * Fy;
        Mnyt = Iy / ( W - xcg ) * Fy;
    	EFIFT  = omega_b * ( Mx / Mnxt + My / Mnyt) + omega_t * P / (A * Fy);

        EFIFT2 = omega_b * ( Mx / Mnx + My / Mny ) - omega_t * P / (A * Fy);

    	if ( EFIFT2 > EFIFT){
        	EFIFT = EFIFT2;
    	}

        // Solo con el proposito de ahorrar variables. No implica relación real:
        EFIFC = EFIFT;

    }

	else {
        // Flexo Compresión
        // Eficiencia de la seccion sometida a flexion en dos ejes y compresion.
        // Articulo C5.2.1 del AISI con metodologia ASD.
        //jacob.comentario pag 223. Hancock
        if( ( omega_c * P / Pn ) <= 0.15 ){
            Pex = 1.0;
            Pey = 1.0;
            ax  = 1.0;
            ay  = 1.0;
        }
        else{
            Pex = PICU * E * Ix / ( LongX * LongX );
            Pey = PICU * E * Iy / ( LongY * LongY );
            // Siempre Cmx/ax >= y Cmy/ay >= 1            
            ax  = 1.0 - omega_c * P / Pex;
            ay  = 1.0 - omega_c * P / Pey;
            if( ax < 1e-5 ) ax = 1e-5; 
            if( ay < 1e-5 ) ay = 1e-5; 
            if( ax > 1.0 )  ax = 1.0;
            if( ay > 1.0 )  ay = 1.0;
        }

        // Eficiencia de la seccion.
        EFIFC = omega_c * P / Pn + omega_b * ( ( Mx / ( ax * Mnx ) ) + ( My / ( ay * Mny ) ) );

        // Comentado por Ernesto
        /*
        //jacob.modif cambio de eficiencia a 10 cuando ay o ax sea menor a 0
        if( ay < 0.0 || ax < 0.0 ){
          //printf("Elemento %d \n\n\n\n",ielem);
          EFIFC=10;
        }
        //*/

    }

	// Eficiencia de la seccion a Flexion y Cortante.
	//jacob.modif
    //se lesaca raiz a la suma de los cuadrados de las relaciones Actuante/nominal
    // En realidad no importa mucho la reaiz por 0 < efi < 1
    EFIFV = ( omega_b * Mx / Mnx) * ( omega_b * Mx / Mnx)   +
	        ( omega_b * My / Mny) * ( omega_b * My / Mny)   +
	        ( omega_v * Vx / Vnx ) * ( omega_v * Vx / Vnx ) +
	        ( omega_v * Vy / Vny ) * ( omega_v * Vy / Vny );

    EFIFV=sqrt(EFIFV);
	// EFIciencia final de la seccion.
	if( EFIFC < EFIFV){
		EFI = EFIFV;
	}
	else{
		EFI = EFIFC;
	}

	//-------------------------------------------------------------------------//



    /*
     fprintf(fp16,"elemento = %d eficiencia =\t%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t \n",ielem,EFI,EFIFC,EFIFV,
     omega_c*P/Pn       , omega_b*Mx/(ax*Mnx), omega_b*My/(ay*Mny),
     omega_b * Mx / Mnx , omega_b * My / Mny , omega_v * Vx / Vn);
    */
    //jacob.impresion provisional
    //elemento = 1 eficiencia =EFI  EFIFC  EFIFV efiAXIAL efiMX1 efiMY1 efiMX2 efiMY2 efiCORT
    //if(ielem==157||ielem==176)printf("Elemento %d  omega_b=%f  Mx=%f   Mnx=%f  EFi2=%f\n",ielem,omega_b,Mx, Mnx, omega_b * Mx / Mnx);
    //if(ielem==157||ielem==176)printf("Elemento %d  omega_b=%f  My=%f   Mny=%f  EFi2=%f\n",ielem,omega_b,My, Mny, omega_b * My / Mny);

    /*
	if(ielem==1){
        printf("Elemento %d  EFIFC=%f    EFIFV=%f  \n",ielem, EFIFC, EFIFV);
        printf("Elemento %d  omega_b=%f  Mx=%f   Mnx=%f  EFi1=%f\n",ielem,omega_b, Mx, Mnx, omega_b * Mx / Mnx);
        printf("Elemento %d  omega_b=%f  My=%f   Mny=%f  EFi2=%f\n",ielem,omega_b,My, Mny, omega_b * My / Mny);
        printf("Elemento %d  omega_v=%f  Vx=%f   Vn=%f  EFi3=%f\n\n",ielem,omega_v,Vx,Vnx, omega_v * Vx / Vnx);
        printf("Elemento %d  omega_v=%f  Vx=%f   Vn=%f  EFi3=%f\n\n",ielem,omega_v,Vy,Vny, omega_v * Vy / Vny);
	 }
	 //*/

    //printf( "\nEFI_C = %lf", EFI );
	return EFI;

}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
double aecsc( double H, double W, double D, double R, double t, double F_eval, double E)
{
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
	//  Bibliotecas incluidas:
	// #include <iostream.h>
	// #include <fstream.h>
	// #include <stdlib.h>
	// #include <math.h>
	//----------------------------------------------------------------------------//

	//****************************************************************************//
	// Inicio de la ejecucion del programa.
	// float AreaEfectivaCompresion(double H, double W, double D, double R, double t, double F_eval, double E)

	//----------------------------------------------------------------------------//
	//                        Declaracion de variables                            //
	//----------------------------------------------------------------------------//
	double r;     // Radio medio del cuarto de rondana
	double u;     // Longitud de arco de la rondana
	//double c;     // Distancia del centroide al centro del radio del cuarto de rondana
	double h;     // Longitud plana del Web o Peralte
	double w;     // Longitud plana del Flange o Ancho
	double d;     // Longitud plana del Lip o Labio
	double S;     // Coeficiente para evaluar pandeo local
	double K;     // Coeficiente de rigidez del Flange
	double Ka;    // Coeficiente de rigidez del Labio respecto al Flange
	double Ku;    // Coeficiente de pandeo local para el Labio a compresion
	double n;     // Constante para la evaluacion de la rigidez del Flange o Ancho
	double Ia;    // Momento de Inercia requerido del atiesador lateral del Flange (el labio)
	double Is;    // Momento de Inercia presente  del atiesador lateral del Flange respecto a su eje centroidal
	//double C1;    // Coeficiente relacionado con C2
	double C2;    // Razon de momento de inercia del atiesador lateral del Flange presente respecto al requerido
	double Lam_W; // Factor de Esbeltez del Web o Peralte
	double Lam_F; // Factor de Esbeltez del Flange o Ancho
	double Lam_L; // Factor de Esbeltez del Lip o Labio
	double Rho_W; // Factor de Reduccion de longitud efectiva del Web o Peralte
	double Rho_F; // Factor de Reduccion de longitud efectiva del Flange o Ancho
	double Rho_L; // Factor de Reduccion de longitud efectiva del Lip o Labio
	double h_efe; // Longitud efectiva del Web o Peralte
	double w_efe; // Longitud efectiva del Flange o Ancho
	double d_efe; // Longitud efectiva del Lip o Labio
	double A_efe; // Area efectiva de la seccion a compresion
	//----------------------------------------------------------------------------//

	// Propiedades de la seccion completa
	r = R + 0.5 * t;           // Radio medio del cuarto de rondana
	u = 1.570 * r;             // Longitud de arco de la rondana
	//c = 0.637 * r;             // Distancia del centroide al centro del
	// radio del cuarto de rondana
	h = H - 2.0 * ( R + t);   // Longitud plana del peralte
	w = W - 2.0 * ( R + t);   // Longitud plana del Flange
	d = D - 1.0 * ( R + t);   // Longitud plana del Labio

	// Area Efectiva con el esfuerzo al cual se desea calcular el area efectiva

	// Evaluacion de la longitud efectiva del Web o Peralte a compresion uniforme
	// segun la metodologia del AISI B2 (pagina 214 de Hancock).
	// Factor de Esbeltez del Web
	Lam_W = (1.052/sqrt(4.0)) * (h/t) * sqrt(F_eval/E);
	// Factor de Reduccion de longitud efectiva del Web que no puede ser mayor que 1.0
	if(Lam_W <= 0.673){
		Rho_W = 1.0;
	}
	else{
		Rho_W = (1.0 - 0.22/Lam_W) / Lam_W;
	}
	if(Rho_W > 1.0){
		Rho_W = 1.0;
	}
	// Longitud efectiva del Web o Peralte
	h_efe = Rho_W * h;

	// Evaluacion de las longitudes efectivas del Flange y del Labio superiores
	// a compresion segun la metodologia del AISI B4.2 (pagina 128 de Yu).
	S  = 1.28 * sqrt(E/F_eval);  // Coeficiente para evaluar pandeo local
	Ku = 0.43;               // Coeficiente de pandeo local para el Labio a
	// compresion que esta simplemente soportado en
	// tres lados y libre en uno (pagina 94 de Yu)

	// Determinacion del caso a evaluar
	if((w/t) <= (S/3.0)){ //----------------------------------------------- Caso I
		// La longitud efectiva del labio y del ancho de la seccion se
		// evaluara usando el CASO I del articulo B4.2 del AISI.
		// Para este caso la el Flange es totalmente efectivo por lo que no
		// es necesario calcular los factores de reduccion
		n     = 0.0;
		Ia    = 0.0;
		Is    = 0.0;
		C2    = 0.0;
		//C1    = 0.0;
		Ka    = 0.0;
		K     = 0.0;
		Lam_F = 0.0;
		Rho_F = 0.0;
		Lam_L = 0.0;
		Rho_L = 0.0;
		// Longitud efectiva del Flange o Ancho en compresion uniforme
		w_efe =   w;
		// Longitud efectiva del Lip o Labio en compresion uniforme
		d_efe =   d;
	}
	else{
		if((w/t) < S){   //----------------------------------------------- Caso II
			// La longitud efectiva del labio y del ancho de la seccion se
			// evaluara usando el CASO II del articulo B4.2 del AISI.
			// Constante para la evaluacion de la rigidez del Flange
			n  = 1.0 / 2.0;
			// Momento de Inercia requerido del atiesador lateral del Flange
			// (el labio superior)
			Ia = 399.00 * (((w/t)/S) - (sqrt(Ku)/2.0)) * (((w/t)/S) - (sqrt(Ku)/2.0)) * (((w/t)/S) - (sqrt(Ku)/2.0)) * t * t * t * t;
		}
		else{            //----------------------------------------------- Caso III
			// La longitud efectiva del labio y del ancho de la seccion se
			// evaluara usando el CASO III del articulo B4.2 del AISI.
			// Constante para la evaluacion de la rigidez del Flange
			n  = 1.0 / 3.0;
			// Momento de Inercia requerido del atiesador lateral del Flange
			// (el labio superior)
			Ia = ((115.0*(w/t)/S) + 5.0) * t * t * t * t;
		}
		// Momento de Inercia del atiesador lateral del Flange (el labio
		// superior) respecto a su eje centroidal
		Is = d*d*d*t / 12.0;
		// Razon de momento de inercia del atiesador lateral del Flange
		// presente respecto al requerido que no puede ser mayor que 1.0
		C2 = Is / Ia;
		if(C2 > 1.0){
			C2 = 1.0;
		}
		//C1 = 2.0 - C2;
		// Coeficiente de rigidez del Labio respecto al Flange que no puede
		// ser mayor que 4.0
		Ka = 5.25 - 5.0 * (D/w);
		if(Ka > 4.0){
			Ka = 4.0;
		}
		// Coeficiente de rigidez del Flange que no puede ser mayor que 4.0
		K  = pow(C2,n) * (Ka - Ku) + Ku;
		if(K > 4.0){
			K  = 4.0;
		}
		// Determinacion del Factor de Esbeltez del Flange
		Lam_F = (1.052/sqrt(K)) * (w/t) * sqrt(F_eval/E);
		// Factor de Reduccion de longitud efectiva del Flange que no puede
		// ser mayor que 1.0
		if(Lam_F <= 0.673){
			Rho_F = 1.0;
		}
		else{
			Rho_F = (1.0 - 0.22/Lam_F) / Lam_F;
		}
		if(Rho_F > 1.0){
			Rho_F = 1.0;
		}
		// Determinacion del Factor de Esbeltez del Labio
		Lam_L= (1.052/sqrt(Ku)) * (d/t) * sqrt(F_eval/E);
		// Factor de Reduccion de longitud efectiva del Labio que no puede
		// ser mayor que 1.0
		if(Lam_L <= 0.673){
			Rho_L = 1.0;
		}
		else{
			Rho_L = (1.0 - 0.22/Lam_L) / Lam_L;
		}
		if(Rho_L > 1.0){
			Rho_L = 1.0;
		}
		// Longitud efectiva del Flange o Ancho en compresion uniforme
		w_efe = Rho_F * w;
		// Longitud efectiva del Lip o Labio en compresion uniforme
		d_efe = Rho_L * C2 * d;
	}

	// Area efectiva a esfuerzo de evaluacion
	A_efe = t * ( 2.0*w_efe + h_efe + 4.0*u + 2.0*d_efe);

	return A_efe;
}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void mnfsc( double H, double W, double D, double R, double t, double E, double Fcx,
            double Fcy1, double Fcy2, double* Mnx, double* Mny1, double* Mny2, int ielem )
{
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

	// ------------------------- VARIABLES LOCALES --------------------------- //
	// ESCALARES:
	int     i;         // Variables de control de ciclos.
	int     conta;     // CONTAdor de iteraciones.
	double  tempo;     // Variable TEMPOral.
	double  DOVO;      // DOceaVO = 1.0 / 12.0.
	double  ACR;       // Area de Cuarto de Rondana.
	double  CCR;       // Centroide de Cuarto de Rondana.
	double  ICR;       // Momento de Inercia de Cuarto de Rondana.
	double  h;         // Longitud plana del peralte.
	double  w;         // Longitud plana del Flange.
	double  d;         // Longitud plana del Labio.
	double  d_efe;     // Longitud efectiva del Labio.
	double  h_efe;     // Longitud efectiva del Web.
	double  h_efe_1;   // Longitud efectiva del Web encima del centroide.
	double  h_efe_2;   // Longitud efectiva del Web debajo del centroide.
	double  h_inefe;   // Longitud inefectiva del Web.
	double  w_efe;     // Longitud efectiva del Flange.
	double  w_efe_1;   // Longitud efectiva del Flange encima del centroide.
	double  w_efe_2;   // Longitud efectiva del Flange debajo del centroide.
	double  w_inefe;   // Longitud inefectiva del Flange.
	double  sum_efe_a; // Suma de longitudes efectivas Anterior.
	double  sum_efe_s; // Suma de longitudes efectivas Siguiente.
	double  f1;        // Esfuerzo a compresion encima del centroide.
	double  f2;        // Esfuerzo a tension    debajo del centroide.
	double  Psi;       // Razon de esfuerzos a compresion y tension.
	double  n;         // Constante para evaluacion de rigidez de Flange.
	double  K;         // Coeficiente de rigidez del Flange.
	double  Ka;        // Coeficiente de rigidez del Labio respecto a Flange.
	double  Kf;        // Ceoficiente de rigidez del Flange.
	double  Ku;        // Coeficiente de pandeo local para el Labio.
	double  Kw;        // Ceoficiente de rigidez del Web.
	double  S;         // Coeficiente para evaluar pandeo local.
	double  Is;        // Momento de Inercia del atiesador lateral del Flange.
	double  Ia;        // Momento de Inercia requerido de atiesador lateral.
	double  C2;        // Razon de momentos de inercia de atiesador lateral.
	double  Lam_F;     // Factor de Esbeltez del Flange.
	double  Lam_H;     // Factor de Esbeltez del Lomo.
	double  Lam_L;     // Factor de Esbeltez del Labio.
	double  Lam_W;     // Factor de Esbeltez del Web.
	double  Rho_F;     // Factor de Reduccion del Flange.
	double  Rho_H;     // Factor de Reduccion del Lomo.
	double  Rho_L;     // Factor de Reduccion del Labio.
	double  Rho_W;     // Factor de Reduccion del Web.
	double  Area;      // Area total de la seccion.
	double  Mome;      // Momento de area total de la seccion.
	double  CC;        // Coordenada Centroidal.
	double  Cy;        // Distancias maximas de fibras extremas a ejes cent.
	double  Ixx;       // Momento de Inercia en X.
	double  Iyy;       // Momento de Inercia en Y.
	double  Iy2;       // Momento de Inercia en Y con giro de 180.
	double  Qx;        // Momento estatico de Area en X.
	double  Qy;        // Momento estatico de Area en Y.
	double  Qy2;       // Momento estatico de Area en Y con giro de 180.
	double  Sx;        // Modulo de Seccion en X.
	double  Sy;        // Modulo de Seccion en Y.
	double  Sy2;       // Modulo de Seccion en Y con giro de 180.
	double f_real;     // Esfuerzo de compresion en la fibra superior,
	                   // considerando que se alcanza la fluencia en la
					   // parte inferior (Mny2)
	double f_supuesto; // Esfuerzo "f" empleado para calcular la longitud
	                   // efectiva sobre el lado de "h" (alma)

	// ARREGLOS:
	double AREAS[11];  // Areas de las subsecciones.
	double DISTA[11];  // DISTAncias de centroides o ejes locales a globales.
	double MILOC[11];  // Momentos de Inercia LOCales.
	double MELOC[11];  // Momentos Estaticos LOCales.
	// ----------------------------------------------------------------------- //

    //jacob.impresion provisional

    /* if(Fcy1<0){
	 printf("\n\n\n\n\n\n\n\n\n\n\n - - - -     Fcy1=%f     - - - - - - - - -    - - - - - - - - - - - - -  - - - - - - - - - - - - - - -\n",Fcy1,Fcy1);
	 }//*/

    //jacob.impresion provisional
    //printf("Fcy1=%f   Fcy2=%f ",Fcy1, Fcy2);

	// Un DOceaVO.
	DOVO  = 1.0 / 12.0;
	tempo = 0.0;

	// Longitudes planas de la seccion.
	h = H - 2.0 * ( R + t);   // Peralte.
	w = W - 2.0 * ( R + t);   // Flange.
	d = D - 1.0 * ( R + t);   // Labio.

	// Propiedades de Cuarto de Rondana (Area, Centroide y Momento de Inercia).
	ACR = 0.7853981633974483 * ( 2.0*R*t + t*t);
	CCR = 0.4244131815783876 * ( (3.0*R*R + 3*R*t + t*t) / ( 2*R + t));
	ICR = 0.1963495408493621 * ( 4.0*R*R*R*t + 6.0*R*R*t*t + 4.0*R*t*t*t + t*t*t*t);

	//-------------------------------------------------------------------------//
	//        Determinacion de las propiedades de Area de la seccion           //
	//         cuando la Flexion actua en el eje X (Ixx, Sx, Qx, tx)           //
	//-------------------------------------------------------------------------//
	//                                    w                                    //
	//                               <--------->                               //
	//                                                                         //
	//                              /-----------\                              //
	//                             / ----------- \                             //
	//                            / /           \ \                            //
	//                       ^   | |             | |  ^                        //
	//                       |   | |             | |  |                        //
	//                       |   | |             | |  | d                      //
	//                       |   | |             | |  |                        //
	//                       |   | |             | |  v                        //
	//                       |   | |                                           //
	//                       |   | |                                           //
	//                     h |   | |                                           //
	//                       |   | |                                           //
	//                       |   | |                                           //
	//                       |   | |             | |  ^                        //
	//                       |   | |             | |  |                        //
	//                       |   | |             | |  | d                      //
	//                       |   | |             | |  |                        //
	//                       v   | |             | |  v                        //
	//                            \ \           / /                            //
	//                             \ ----------- /                             //
	//                              \-----------/                              //
	//                                                                         //
	// Evaluacion de las longitudes efectivas del Flange y del Labio           //
	// superiores a compresion y del Lomo con gradiente de esfuerzo            //
	// para evaluar las propiedades de la seccion cuando esta sometida         //
	// a flexion en el eje X segun metodologia AISI B4.2 (pag 128 de Yu).      //
	//-------------------------------------------------------------------------//

	// Ceoficientes de rigidez y pandeo del Web.
	S  = 1.28 * sqrt(E/Fcx);
	Ku = 0.43;

	// Determinacion del caso a evaluar segun el articulo B4.2 del AISI.
	if( (w/t) <= (S/3.0)){ // Caso I
		// En este caso el Flange y el Labio son
		// totalmente efectivos por lo que no es
		// necesario calcular factores de reduccion.
		n     = 0.0;
		Ia    = 0.0;
		Is    = 0.0;
		C2    = 0.0;
		Ka    = 0.0;
		K     = 0.0;
		Lam_F = 0.0;
		Rho_F = 0.0;
		Lam_L = 0.0;
		Rho_L = 0.0;

		// Longitud efectiva del Flange
		// a compresion (el superior).
		w_efe = w;

		// Longitud efectiva del Labio
		// a compresion (el superior).
		d_efe = d;
	}
	else{
		if( (w/t) < S){ // Caso II
			// Constante de rigidez del Flange.
			n  = 1.0 / 2.0;
			// Momento de Inercia requerido del atiesador
			// lateral del Flange (el labio superior).
			tempo = ((w/t)/S) - (sqrt(Ku)/2.0);
			tempo = tempo * tempo * tempo;
			Ia = 399.00 * tempo * t * t * t * t;


		}
		else{ // Caso III
			// Constante de rigidez del Flange.
			n  = 1.0 / 3.0;
			// Momento de Inercia requerido del atiesador
			// lateral del Flange (el labio superior).
			Ia = ((115.0*(w/t)/S) + 5.0) * t * t * t * t;
		}

		// Momento de Inercia del atiesador lateral
		// del Flange (el labio superior).
		Is = DOVO * d * d * d * t;

		// Razon de momento de inercia del atiesador
		// lateral del Flange presente respecto al
		// requerido que no puede ser mayor que 1.0.
		C2 = Is / Ia;
		if( C2 > 1.0){
			C2 = 1.0;
		}

		// Coeficiente de rigidez del Labio respecto
		// al Flange que no puede ser mayor que 4.0.
		Ka = 5.25 - 5.0 * (D/w);
		if( Ka > 4.0){
			Ka = 4.0;
		}

		// Coeficiente de rigidez del Flange
		// que no puede ser mayor que 4.0.
		K  = pow(C2,n) * (Ka - Ku) + Ku;
		if( K > 4.0){
			K = 4.0;
		}

		// Factor de Esbeltez del Flange.
		Lam_F = (1.052/sqrt(K)) * (w/t) * sqrt(Fcx/E);

		// Factor de Reduccion de longitud efectiva
		// del Flange que no puede ser mayor que 1.0.
		if( Lam_F <= 0.673){
			Rho_F = 1.0;
		}
		else{
			Rho_F = (1.0 - 0.22/Lam_F) / Lam_F;
		}
		if( Rho_F > 1.0){
			Rho_F = 1.0;
		}

		// Factor de Esbeltez del Labio.
		Lam_L = (1.052/sqrt(Ku)) * (d/t) * sqrt(Fcx/E);

		// Factor de Reduccion de longitud efectiva
		// del Labio que no puede ser mayor que 1.0.
		if( Lam_L <= 0.673){
			Rho_L = 1.0;
		}
		else{
			Rho_L = (1.0 - 0.22/Lam_L) / Lam_L;
		}
		if( Rho_L > 1.0){
			Rho_L = 1.0;
		}

		// Longitud efectiva del Flange
		// a compresion (el superior)
		w_efe = Rho_F * w;

		// Longitud efectiva del Labio
		// a compresion (el superior)
		d_efe = Rho_L * C2 * d;
	}

	// Evaluacion de las longitud efectiva del Web considerandolo como un
	// elemento atiesado sujeto a gradiente de esfuerzo segun la metodo-
	// logia del AISI B2.3 (pagina 113 de Yu).

	// Inicializacion de variables para el proceso iterativo.
	//jacob. modif provisional
    /*
     h_efe_1 = 0.4 * h;
	 h_efe_2 = 0.1*h;
	 h_inefe = 0.2*h;
	 //*/

    h_efe_1 = 0.5 * h;  // Longitud efectiva del Web en la parte superior.
	h_efe_2 = 0.0;      // Longitud efectiva del Web en la parte inferior.
	h_inefe = 0.0;      // Longitud no efectiva del Web que no trabaja.

	// Inicializacion de controladores del ciclo while.
	sum_efe_a = h;                 // Suma de longitudes efectivas del Web.
	sum_efe_s = h_efe_1 + h_efe_2; // Suma de longitudes efectivas del Web.
	conta = 0;

	// Inicio de iteraciones para encontrar la longitud efectiva
	// del Web, esta restringido a 1000 iteraciones.
	while((fabs((sum_efe_s-sum_efe_a)) > 1.0e-6) && (conta < 1000) && (h_inefe >= 0.0)){

        //jacob. agregado para impresion
        //printf("Antes del calculo :\n h_efe_1= %5f     h_efe_2= %5f     h_inefe= %5f  \n",h_efe_1,h_efe_2,h_inefe);



		// Localizacion del eje neutro y calculo del momento de inecia
		// Ixx y Modulo de seccion Sx respecto al eje X tomando como
		// referencia el eje que pasa sobre la fibra superior.

		// Areas de cada sub-seccion.
		AREAS[0]  = d * t;
		AREAS[1]  = ACR;
		AREAS[2]  = w * t;
		AREAS[3]  = AREAS[1];
		AREAS[4]  = (h - h_efe_1 - h_inefe) * t;
		AREAS[5]  = h_efe_1 * t;
		AREAS[6]  = AREAS[3];
		AREAS[7]  = 0.5 * w_efe * t;
		AREAS[8]  = AREAS[7];
		AREAS[9]  = AREAS[6];
		AREAS[10] = d_efe * t;

		// Distancias de centroides de cada
		// sub-seccion a la fibra superior.
		DISTA[0]  = h + R + t - 0.5*d;
		DISTA[1]  = h + R + t + CCR;
		DISTA[2]  = h + 2.0*R + 1.5*t;
		DISTA[3]  = DISTA[1];
		DISTA[4]  = R + t + 0.5*h + 0.5*h_efe_1 + 0.5*h_inefe;
		DISTA[5]  = R + t + 0.5*h_efe_1;
		DISTA[6]  = R + t - CCR;
		DISTA[7]  = 0.5*t;
		DISTA[8]  = DISTA[7];
		DISTA[9]  = DISTA[6];
		DISTA[10] = R + t + 0.5*d_efe;

		// Distancia de fibra superior al eje neutro.
		// Momento de area total de la seccion.
		Mome = 0.0;
		for ( i = 0 ; i < 11 ; i++){
			Mome = Mome + AREAS[i] * DISTA[i];
		}

		// Area total de la seccion.
		Area = 0.0;
		for ( i = 0 ; i < 11 ; i++){
			Area = Area + AREAS[i];
		}

		// Distancia de la fibra superior al eje neutro.
		CC = Mome / Area;

		// Esfuerzos en la parte plana del peralte (Web).
		f1 = +Fcx * ((CC-R-t)/CC);     // Compresion en fibra superior.
		f2 = -Fcx * ((h-(CC-R-t))/CC); // Tension    en fibra inferior.

		// Razon de esfuerzos presentes en el Web.
		Psi = f2 / f1;

		// Ceoficiente de rigidez del Web.
		tempo = 1.0 - Psi;
		Kw = 4.0 + 2.0*tempo*tempo*tempo + 2.0*tempo;

		// Factor de Esbeltez del Web.
		Lam_W = (1.052/sqrt(Kw)) * (h/t) * sqrt(f1/E);

		// Factor de Reduccion de longitud efectiva
		// del Web que no puede ser mayor que 1.0.
		if( Lam_W <= 0.673){
			Rho_W = 1.0;
		}
		else{
			Rho_W = (1.0 - 0.22/Lam_W) / Lam_W;
		}
		if( Rho_W > 1.0){
			Rho_W = 1.0;
		}

		// Longitud efectiva del Web.
		h_efe = Rho_W * h;

		// Longitud efectiva del Web en la parte superior.
		h_efe_1 = h_efe / (3.0-Psi);

		// Longitud efectiva del Web en la parte inferior.
		if( Psi <= -0.236){
			h_efe_2 = h_efe / 2.0;
		}
		else{
			h_efe_2 = h_efe - h_efe_1;
		}

		// Longitud no efectiva del Web.
		h_inefe = (CC-R-t) - (h_efe_1+h_efe_2);

		// Actualizacion de la suma de
		// longitudes efectivas en el Web.
		sum_efe_a = sum_efe_s;
		sum_efe_s = h_efe_1 + h_efe_2;

        //jacob. agregado para impresion
        //printf("despues del calculo :\n h_efe_1= %5f     h_efe_2= %5f     h_inefe= %5f  \n",h_efe_1,h_efe_2,h_inefe);
        //jacob. agregada impresion
        //printf("contador= %d\n",conta);
        //printf("sum_efe_s= %5f       sum_efe_a = %5f    h_infe=%f\n\n",sum_efe_s,sum_efe_a,h_inefe);

		// Actualizacion del contador.
		conta = conta + 1;
	}



	// Correccion de longitud no efectiva del Web
	// o Peralte que no puede se negativa.
	if( h_inefe < 0.0){
		h_inefe = 0.0;
		h_efe_1 = 0.5 * h;
		h_efe_2 = 0.0;
	}

	// Areas de cada sub-seccion.
	AREAS[0]  = d * t;
	AREAS[1]  = ACR;
	AREAS[2]  = w * t;
	AREAS[3]  = AREAS[1];
	AREAS[4]  = (h - h_efe_1 - h_inefe) * t;
	AREAS[5]  = h_efe_1 * t;
	AREAS[6]  = AREAS[3];
	AREAS[7]  = 0.5 * w_efe * t;
	AREAS[8]  = AREAS[7];
	AREAS[9]  = AREAS[6];
	AREAS[10] = d_efe * t;

	// Distancias de centroides de cada
	// sub-seccion a la fibra superior.
	DISTA[0]  = h + R + t - 0.5*d;
	DISTA[1]  = h + R + t + CCR;
	DISTA[2]  = h + 2.0*R + 1.5*t;
	DISTA[3]  = DISTA[1];
	DISTA[4]  = R + t + 0.5*h + 0.5*h_efe_1 + 0.5*h_inefe;
	DISTA[5]  = R + t + 0.5*h_efe_1;
	DISTA[6]  = R + t - CCR;
	DISTA[7]  = 0.5*t;
	DISTA[8]  = DISTA[7];
	DISTA[9]  = DISTA[6];
	DISTA[10] = R + t + 0.5*d_efe;

	// Distancia de fibra superior al eje neutro.
	// Momento de area total de la seccion.
	Mome = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		Mome = Mome + AREAS[i] * DISTA[i];
	}

	// Area total de la seccion.
	Area = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		Area = Area + AREAS[i];
	}

	// Distancia de la fibra superior al eje neutro.
	CC = Mome / Area;

	// Momentos de Inercia en X respecto a
	// ejes locales de cada sub-seccion.
	MILOC[0]  = DOVO * d * d * d * t;
	MILOC[1]  = ICR;
	MILOC[2]  = DOVO * t * t * t * w;
	MILOC[3]  = MILOC[1];
	MILOC[4]  = DOVO * (h - h_efe_1 - h_inefe) * (h - h_efe_1 - h_inefe) * (h - h_efe_1 - h_inefe) * t;
	MILOC[5]  = DOVO * h_efe_1 * h_efe_1 * h_efe_1 * t;
	MILOC[6]  = MILOC[3];
	MILOC[7]  = 0.5 * DOVO * t * t * t * w_efe;
	MILOC[8]  = MILOC[7];
	MILOC[9]  = MILOC[6];
	MILOC[10] = DOVO * d_efe * d_efe * d_efe * t;

	// Distancia paralela al eje Y del
	// centroide de la seccion a los ejes
	// inerciales locales de cada sub-seccion.
	DISTA[0]  = h + R + t - 0.5*d - CC;
	DISTA[1]  = h + R + t - CC;
	DISTA[2]  = h + 2.0*R + 1.5*t - CC;
	DISTA[3]  = DISTA[1];
	DISTA[4]  = h + R + t - 0.5*(h - h_efe_1 - h_inefe) - CC;
	DISTA[5]  = R + t + 0.5*h_efe_1 - CC;
	DISTA[6]  = R + t - CC;
	DISTA[7]  = 0.5*t - CC;
	DISTA[8]  = DISTA[7];
	DISTA[9]  = DISTA[6];
	DISTA[10] = R + t + 0.5*d_efe - CC;

	// Momento de inercia de la seccion respecto al eje X.
	Ixx = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		Ixx = Ixx + MILOC[i] + DISTA[i] * DISTA[i] * AREAS[i];
        //jacob.impresion
        //if(MILOC[i]<0){printf("Error, MILOC[%d] = %f\n",i,MILOC[i]);}
        //if(DISTA[i]<0){printf("Error, DISTA[%d] = %f\n",i,DISTA[i]);}
        //if(AREAS[i]<0){printf("Error, AREAS[%d] = %f\n\n",i,AREAS[i]);}

	}

    //jacob.impresion provisional
    //if(Ixx<0){printf("Error, \t\t\t   Ixx = %f\n\n\n\n",Ixx);}


	// /////////////////////////////////////////////////////////// //
	// Hasta aqui ya esta probado con los ejemplos del libro de Yu //
	// y varia un poco con los resultados del ejemplo 4.2          //
	// /////////////////////////////////////////////////////////// //

	// Distancias maximas de las fibras
	// extremas a los ejes centroidales.
	Cy = CC;
	if( Cy < fabs(H-Cy)){
		Cy = fabs(H-Cy);
	}

	// Modulo de Seccion calculado con la maxima
	// distancia del centroide a fibras extremas.
	Sx = Ixx / Cy;

    //jacob.impresion provisional
    //if(ielem==176||ielem==157){printf("Sx= %f \n",Sx);}
    //printf("Sx=%f\t Ixx=%f\t Cy=%f\n",Sx,Ixx,Cy);

	// Momentos estaticos de area de las porciones de la seccion
	// que quedan arriba del eje neutro y que trabajan a compre-
	// sion respecto al eje neutro (solo los negativos contribui-
	// ran, los positivos trabajan a tension y no cuentan para
	// la evaluacion del esfuerzo cortante en la seccion).
	MELOC[0]  = (h + R + t - 0.5*d - CC) * (d * t);
	MELOC[1]  = (h + R + t + CCR - CC) * ACR;
	MELOC[2]  = (h + 2.0*R + 1.5*t - CC) * (w * t);
	MELOC[3]  = (h + R + t + CCR - CC) * ACR;
	MELOC[4]  = (- 0.5 * h_efe_2) * (h_efe_2 * t);
	MELOC[5]  = (R + t + 0.5*h_efe_1 - CC) * (h_efe_1 * t);
	MELOC[6]  = (R + t - CCR - CC) * ACR;
	MELOC[7]  = (0.5*t - CC) * (0.5 * w_efe * t);
	MELOC[8]  = (0.5*t - CC) * (0.5 * w_efe * t);
	MELOC[9]  = (R + t - CCR - CC) * ACR;
	MELOC[10] = (R + t + 0.5*d_efe - CC) * (d_efe * t);

	// Momento estatico de la seccion
	// por encima del eje neutro.
	Qx = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		if( MELOC[i] < 0.0){
			Qx = Qx - MELOC[i];
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	//        Determinacion de las propiedades de Area de la seccion           //
	//         cuando la Flexion actua en el eje Y (Iyy, Sy, Qy, ty)           //
	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                              d                                          //
	//                         <--------->                                     //
	//                                                                         //
	//                        /-----------       -----------\                  //
	//                       / -----------       ----------- \                 //
	//                      / /                             \ \                //
	//                 ^   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//               w |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 v   | |                               | |               //
	//                      \ \                              / /               //
	//                       \ ------------------------------ /                //
	//                        \------------------------------/                 //
	//                                                                         //
	//                         <---------------------------->                  //
	//                                       h                                 //
	//                                                                         //
	// Evaluacion de las longitudes efectivas de los Labios superiores         //
	// a compresion y del Flange con gradiente de esfuerzo para evaluar        //
	// las propiedades de la seccion cuando este sometida a flexion en         //
	// el eje Y segun la metodologia del AISI B4.2 (pagina 122 de Yu).         //
	//-------------------------------------------------------------------------//

	// Factor de Esbeltez del Labio.
	Lam_L = (1.052/sqrt(Ku)) * (d/t) * sqrt(Fcy1/E);

    //jacob.impresion provisional

    /*
     if(Fcy1<1){
	 printf(" - - - - - - - - - - - - - - - - -");
	 printf("Ku=%f    Fcy1= %f     E=%f  \n",Ku,Fcy1,E);
	 }
	 //*/


	// Factor de Reduccion de longitud efectiva
	// del Labio que no puede ser mayor que 1.0.

    //jacob.impresion provisional
    /*
     printf("Lam_L=%f   \n",Lam_L);
     //*/

    if( Lam_L <= 0.673){
		Rho_L = 1.0;
	}
	else{
		Rho_L = (1.0 - 0.22/Lam_L) / Lam_L;
	}
	if( Rho_L > 1.0){
		Rho_L = 1.0;
	}



	// Longitud efectiva de Labios a compresion.
	d_efe = Rho_L * d;

    //jacob.impresion provisional
    //printf("d_efe=%f    Rho_L=%f     d=%f  \n",d_efe,Rho_L,d);


	// Evaluacion de las longitudes efectivas de los Flanges considerandolos
	// como elementos atiesados sujeto a gradiente de esfuerzo segun la
	// metodologia del AISI B2.3 (pagina 113 de Yu).

	// Inicializacion de variables para el proceso iterativo.
	w_efe_1 = 0.5 * w; // Longitud efectiva del Flange en la parte superior.
	w_efe_2 = 0.0;     // Longitud efectiva del Flange en la parte inferior.
	w_inefe = 0.0;     // Longitud no efectiva del Flange que no trabaja.

	// Inicializacion de controladores del ciclo while.
	sum_efe_a = w;                 // Suma de longitudes efectivas del Flange.
	sum_efe_s = w_efe_1 + w_efe_2; // Suma de longitudes efectivas del Flange.
	conta = 0;

	// Inicio de iteraciones para encontrar la longitud efectiva
	// del Flange, esta restringido a 1000 iteraciones.
	while((fabs((sum_efe_s-sum_efe_a)) > 1.e-6) && (conta < 1000) && (w_inefe >= 0.0)){

		// Localizacion del eje neutro y calculo del momento de inecia
		// Iyy y Modulo de seccion Sy respecto al eje Y tomando como
		// referencia el eje que pasa sobre la fibra superior.

		// Areas de cada sub-seccion.
		AREAS[0]  = d_efe * t;
		AREAS[1]  = ACR;
		AREAS[2]  = w_efe_1 * t;
		AREAS[3]  = (w - w_efe_1 - w_inefe) * t;
		AREAS[4]  = AREAS[1];
		AREAS[5]  = h * t;
		AREAS[6]  = AREAS[4];
		AREAS[7]  = AREAS[3];
		AREAS[8]  = AREAS[2];
		AREAS[9]  = AREAS[6];
		AREAS[10] = AREAS[0];

		// Distancias de centroides de cada
		// sub-seccion a la fibra superior.
		DISTA[0]  = 0.5 * t;
		DISTA[1]  = R + t - CCR;
		DISTA[2]  = R + t + 0.5*w_efe_1;
		DISTA[3]  = R + t + 0.5*w + 0.5*w_efe_1 + 0.5*w_inefe;
		DISTA[4]  = R + t + w + CCR;
		DISTA[5]  = w + 2.0*R + 1.5*t;
		DISTA[6]  = DISTA[4];
		DISTA[7]  = DISTA[3];
		DISTA[8]  = DISTA[2];
		DISTA[9]  = DISTA[1];
		DISTA[10] = DISTA[0];

		// Distancia de fibra superior al eje neutro.
		// Momento de area total de la seccion.
		Mome = 0.0;
		for ( i = 0 ; i < 11 ; i++){
			Mome = Mome + AREAS[i] * DISTA[i];
		}

		// Area total de la seccion.
		Area = 0.0;
		for ( i = 0 ; i < 11 ; i++){
			Area = Area + AREAS[i];
		}

		// Distancia de la fibra superior al eje neutro.
		CC = Mome / Area;

		// Esfuerzos en la parte plana de los Flanges.
		f1 = +Fcy1 * ((CC-R-t)/CC);     // Compresion en fibra superior.
		f2 = -Fcy1 * ((w-(CC-R-t))/CC); // Tension    en fibra inferior.

		// Razon de esfuerzos presentes en el Flange.
		Psi = f2 / f1;

		// Ceoficiente de rigidez del Flange.
		tempo = 1.0 - Psi;
		Kf = 4.0 + 2.0*tempo*tempo*tempo + 2.0*tempo;

		// Factor de Esbeltez del Flange.
		Lam_F = (1.052/sqrt(Kf)) * (w/t) * sqrt(f1/E);

		// Factor de Reduccion de longitud efectiva
		// del Flange que no puede ser mayor que 1.0.
		if( Lam_F <= 0.673){
			Rho_F = 1.0;
		}
		else{
			Rho_F = (1.0 - 0.22/Lam_F) / Lam_F;
		}
		if( Rho_F > 1.0){
			Rho_F = 1.0;
		}

		// Longitud efectiva del Flange.
		w_efe = Rho_F * w;

		// Longitud efectiva del Flange en la parte superior.
		w_efe_1 = w_efe / (3.0-Psi);

		// Longitud efectiva del Flange en la parte inferior.
		if( Psi <= -0.236){
			w_efe_2 = w_efe / 2.0;
		}
		else{
			w_efe_2 = w_efe - w_efe_1;
		}

		// Longitud no efectiva del Flange.
		w_inefe = (CC-R-t) - (w_efe_1+w_efe_2);

		// Actualizacion de la suma de
		// longitudes efectivas en Flange.
		sum_efe_a = sum_efe_s;
		sum_efe_s = w_efe_1 + w_efe_2;

		// Actualizacion del contador.
		conta = conta + 1;
	}

	// Correccion de la longitud no efectiva del
	// Flange o Ancho que no puede se negativa.
	if( w_inefe < 0.0){
		w_inefe = 0.0;
		w_efe_1 = 0.5 * w;
		w_efe_2 = 0.0;
	}

	// Areas de cada sub-seccion.
	AREAS[0]  = d_efe * t;
	AREAS[1]  = ACR;
	AREAS[2]  = w_efe_1 * t;
	AREAS[3]  = (w - w_efe_1 - w_inefe) * t;
	AREAS[4]  = AREAS[1];
	AREAS[5]  = h * t;
	AREAS[6]  = AREAS[4];
	AREAS[7]  = AREAS[3];
	AREAS[8]  = AREAS[2];
	AREAS[9]  = AREAS[6];
	AREAS[10] = AREAS[0];

	// Distancias de centroides de cada
	// sub-seccion a la fibra superior.
	DISTA[0]  = 0.5*t;
	DISTA[1]  = R + t - CCR;
	DISTA[2]  = R + t + 0.5*w_efe_1;
	DISTA[3]  = R + t + 0.5*w + 0.5*w_efe_1 + 0.5*w_inefe;
	DISTA[4]  = R + t + w + CCR;
	DISTA[5]  = w + 2.0*R + 1.5*t;
	DISTA[6]  = DISTA[4];
	DISTA[7]  = DISTA[3];
	DISTA[8]  = DISTA[2];
	DISTA[9]  = DISTA[1];
	DISTA[10] = DISTA[0];

	// Distancia de fibra superior al eje neutro.
	// Momento de area total de la seccion.
	Mome = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		Mome = Mome + AREAS[i] * DISTA[i];
        //jacob.impresion provisional
        //printf("AREAS[%d]= %f   DISTA[%d]=%f    d_efe=%f    t=%f \n",i,AREAS[i],i,DISTA[i],d_efe,t);
	}

	// Area total de la seccion.
	Area = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		Area = Area + AREAS[i];
        //jacob.impresion provional
        //printf("AREAS[%d]= %f \n",i,AREAS[i]);
	}

	// Distancia de la fibra superior al eje neutro.
	CC = Mome / Area;

    //jacob.impresion provisional
    //printf("CC= %f     Mome= %f    Area= %f  \n",CC,Mome,Area);


	// Momentos de Inercia en Y respecto a
	// ejes locales de cada sub-seccion.
	MILOC[0]  = DOVO * t * t * t * d_efe;
	MILOC[1]  = ICR;
	MILOC[2]  = DOVO * w_efe_1 * w_efe_1 * w_efe_1 * t;
	MILOC[3]  = DOVO * (w - w_efe_1 - w_inefe) * (w - w_efe_1 - w_inefe) * (w - w_efe_1 - w_inefe) * t;
	MILOC[4]  = MILOC[1];
	MILOC[5]  = DOVO * t * t * t * h;
	MILOC[6]  = MILOC[4];
	MILOC[7]  = MILOC[3];
	MILOC[8]  = MILOC[2];
	MILOC[9]  = MILOC[6];
	MILOC[10] = MILOC[0];

	// Distancia paralela al eje Y del
	// centroide de la seccion a los ejes
	// inerciales locales de cada sub-seccion.
	DISTA[0]  = 0.5*t - CC;
	DISTA[1]  = R + t - CC;
	DISTA[2]  = R + t + 0.5*w_efe_1 - CC;
	DISTA[3]  = w + R + t - 0.5*(w - w_efe_1 - w_inefe) - CC;
	DISTA[4]  = w + R + t - CC;
	DISTA[5]  = w + 2.0*R + 1.5*t - CC;
	DISTA[6]  = DISTA[4];
	DISTA[7]  = DISTA[3];
	DISTA[8]  = DISTA[2];
	DISTA[9]  = DISTA[1];
	DISTA[10] = DISTA[0];

	// Momento de inercia de la seccion respecto al eje Y.
	Iyy = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		Iyy = Iyy + MILOC[i] + DISTA[i] * DISTA[i] * AREAS[i];
	}

	// Distancias maximas de las fibras
	// extremas a los ejes centroidales.

    //jacob.impresion provional
    //printf("CC= %f    abs(W)-abs(Cy)=%f   \n",CC,fabs(fabs(W)-fabs(Cy)));

    Cy = CC;
	if( Cy < fabs(fabs(W)-fabs(Cy))){
		Cy = fabs(fabs(W)-fabs(Cy));
	}

	// Modulo de Seccion calculado con la maxima
	// distancia del centroide a fibras extremas.
	Sy = Iyy / Cy;


    //jacob.impresion provisional
    //printf("Sy= %f     Iyy = %f  Cy= %f  \n ",Sy, Iyy, Cy);


	// Momentos estaticos de area de las porciones de la seccion
	// que quedan arriba del eje neutro y que trabajan a compre-
	// sion respecto al eje neutro (solo los negativos contribui-
	// ran, los positivos trabajan a tension y no cuentan para
	// la evaluacion del esfuerzo cortante en la seccion).
	MELOC[0]  = (0.5*t - CC) * (d_efe * t);
	MELOC[1]  = (R + t - CCR - CC) * ACR;
	MELOC[2]  = (R + t + 0.5*w_efe_1 - CC) * (w_efe_1 * t);
	MELOC[3]  = (- 0.5 * w_efe_2) * (w_efe_2 * t);
	MELOC[4]  = (w + R + t + CCR - CC) * ACR;
	MELOC[5]  = (w + 2.0*R + 1.5*t - CC) * (h * t);
	MELOC[6]  = MELOC[4];
	MELOC[7]  = MELOC[3];
	MELOC[8]  = MELOC[2];
	MELOC[9]  = MELOC[1];
	MELOC[10] = MELOC[0];

	// Momento estatico de la seccion
	// por encima del eje neutro.
	Qy = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		if( MELOC[i] < 0.0){
			Qy = Qy - MELOC[i];
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	//        Determinacion de las propiedades de Area de la seccion           //
	//         cuando la Flexion actua en el eje Y (Iyy, Sy, Qy, ty)           //
	/////////////////////////////////////////////////////////////////////////////
	//                                                                         //
	//                                       h                                 //
	//                         <---------------------------->                  //
	//                                                                         //
	//                        /-----------------------------\                  //
	//                       / ----------------------------- \                 //
	//                      / /                             \ \                //
	//                 ^   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//               w |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 |   | |                               | |               //
	//                 v   | |                               | |               //
	//                      \ \                              / /               //
	//                       \ -----------        ----------- /                //
	//                        \-----------        -----------/                 //
	//                                                                         //
	//                         <--------->                                     //
	//                              d                                          //
	//                                                                         //
	// Evaluacion de las longitudes efectivas de los Labios superiores         //
	// a compresion y del Flange con gradiente de esfuerzo para evaluar        //
	// las propiedades de la seccion cuando este sometida a flexion en         //
	// el eje Y segun la metodologia del AISI B4.2 (pagina 122 de Yu).         //
	// En este caso se evalua la flexion con un giro de 180 en la seccion.     //
	//-------------------------------------------------------------------------//

	// Dado que el punto de fluencia se alcanza en el parte inferior de la seccion
	// (en el lado de tension), se tiene que iterar sobre el esfuerzo empleado para
	// calcular la longitud efectiva en h. Se considera que el "w" es totalmente
	// efectivo.

		// En este caso los siguientes valores no cambian (a excepcion de AREAS[5]
		// que corresponde al area efetiva sobre h.
		// Areas de cada sub-seccion.
		AREAS[0]  = d * t;
		AREAS[1]  = ACR;
		AREAS[2]  = 0.5 * w * t;
		AREAS[3]  = 0.5 * w * t;
		AREAS[4]  = ACR;
		AREAS[5]  = 0.0;
		AREAS[6]  = ACR;
		AREAS[7]  = 0.5 * w * t;
		AREAS[8]  = 0.5 * w * t;
		AREAS[9]  = ACR;
		AREAS[10] = d * t;

		// Distancias de centroides de cada
		// sub-seccion a la fibra superior.
		DISTA[0]  = W - 0.5 * t;
		DISTA[1]  = W - t - R + CCR;
		DISTA[2]  = W - t - R - 0.25 * w;
		DISTA[3]  = W - t - R - 0.75 * w;
		DISTA[4]  = R + t - CCR;
		DISTA[5]  = 0.5 * t;
		DISTA[6]  = R + t - CCR;
		DISTA[7]  = W - t - R - 0.75 * w;
		DISTA[8]  = W - t - R - 0.25 * w;
		DISTA[9]  = W - t - R + CCR;
		DISTA[10] = W - 0.5 * t;

		Area = 0.0;
		Mome = 0.0;
		for( i = 0; i < 11; i++ ){
			Area += AREAS[i];
			Mome += AREAS[i] * DISTA[i];
		}

		conta = 0;
		f_real     = Fcy2;
		f_supuesto = 0.0;  // Solo para que entre al ciclo
		while( fabs( f_real - f_supuesto ) > 1e-6 && conta < 1000 ){

			f_supuesto = f_real;

			// Factor de Esbeltez del Lomo.
			Lam_H = 0.526 * ( h / t) * sqrt( f_supuesto / E );

			// Factor de Reduccion de longitud efectiva
			// del Lomo que no puede ser mayor que 1.0
			if( Lam_H <= 0.673 ){
				Rho_H = 1.0;
			}
			else{
				Rho_H = ( 1.0 - 0.22 / Lam_H ) / Lam_H;
			}
			if( Rho_H > 1.0){
				Rho_H = 1.0;
			}

			// Longitud efectiva del Lomo a compresion
			h_efe = Rho_H * h;

			AREAS[5] = h_efe * t;

			// Calculo del centroide (con respecto a la fibra superior)
			Area += AREAS[5];
			Mome += AREAS[5] * DISTA[5];
			CC = Mome / Area;

			// Esfuerzo de compresion en la fibra superior
			f_real = Fcy2 * CC / ( W - CC );

			Area -= AREAS[5];
			Mome -= AREAS[5] * DISTA[5];
			conta++;

		}

	// Evaluacion de las longitudes efectivas de los Flanges considerandolos
	// como elementos atiesados sujeto a gradiente de esfuerzo segun la
	// metodologia del AISI B2.3 (pagina 113 de Yu).

	// Inicializacion de variables para el proceso iterativo.
	w_efe_1 = 0.5 * w;  // Longitud efectiva del Flange en la parte superior.
	w_efe_2 = 0.0;      // Longitud efectiva del Flange en la parte inferior.
	w_inefe = 0.0;      // Longitud no efectiva del Flange que no trabaja.

	// Inicializacion de controladores del ciclo while.
	sum_efe_a = w;                 // Suma de longitudes efectivas del Flange.
	sum_efe_s = w_efe_1 + w_efe_2; // Suma de longitudes efectivas del Flange.
	conta = 0;

	// Inicio de iteraciones para encontrar la longitud efectiva
	// del Flange, esta restringido a 1000 iteraciones.
	while((fabs((sum_efe_s-sum_efe_a)) > 1.0e-6) && (conta < 1000) && (w_inefe >= 0.0)){

		// Localizacion del eje neutro y calculo del momento de inecia
		// Iyy y Modulo de seccion Sy respecto al eje Y tomando como
		// referencia el eje que pasa sobre la fibra superior.

		// Areas de cada sub-seccion.
		AREAS[0]  = d * t;
		AREAS[1]  = ACR;
		AREAS[2]  = (w - w_efe_1 - w_inefe) * t;
		AREAS[3]  = w_efe_1 * t;
		AREAS[4]  = AREAS[1];
		AREAS[5]  = h_efe * t;
		AREAS[6]  = AREAS[4];
		AREAS[7]  = AREAS[3];
		AREAS[8]  = AREAS[2];
		AREAS[9]  = AREAS[1];
		AREAS[10] = AREAS[0];

		// Distancias de centroides de cada
		// sub-seccion a la fibra superior.
		DISTA[0]  = W - 0.5*t;
		DISTA[1]  = W - t - R + CCR;
		DISTA[2]  = R + t + 0.5*w + 0.5*w_efe_1 + 0.5*w_inefe;
		DISTA[3]  = R + t + 0.5*w_efe_1;
		DISTA[4]  = R + t - CCR;
		DISTA[5]  = 0.5*t;
		DISTA[6]  = DISTA[4];
		DISTA[7]  = DISTA[3];
		DISTA[8]  = DISTA[2];
		DISTA[9]  = DISTA[1];
		DISTA[10] = DISTA[0];

		// Distancia de fibra superior al eje neutro.
		// Momento de area total de la seccion.
		Mome = 0.0;
		for ( i = 0 ; i < 11 ; i++){
			Mome = Mome + AREAS[i] * DISTA[i];
		}

		// Area total de la seccion.
		Area = 0.0;
		for ( i = 0 ; i < 11 ; i++){
			Area = Area + AREAS[i];
		}

		// Distancia de la fibra superior al eje neutro.
		CC = Mome / Area;

		// Esfuerzos en la parte plana de los Flanges.
		f1 = +Fcy2 * ((CC-R-t)/CC);     // Compresion en fibra superior.
		f2 = -Fcy2 * ((w-(CC-R-t))/CC); // Tension    en fibra inferior.

        f1 =  Fcy2 * ( CC - R - t ) / ( W - CC );
		f2 = -Fcy2 * ( W - CC - R - t ) / ( W - CC );

		// Razon de esfuerzos presentes en el Flange.
		Psi = f2 / f1;

		// Ceoficiente de rigidez del Flange.
		tempo = 1.0 - Psi;
		Kf = 4.0 + 2.0*tempo*tempo*tempo + 2.0*tempo;

		// Factor de Esbeltez del Flange.
		Lam_F = (1.052/sqrt(Kf)) * (w/t) * sqrt(f1/E);

		// Factor de Reduccion de longitud efectiva
		// del Flange que no puede ser mayor que 1.0.
		if( Lam_F <= 0.673){
			Rho_F = 1.0;
		}
		else{
			Rho_F = (1.0 - 0.22/Lam_F) / Lam_F;
		}
		if( Rho_F > 1.0){
			Rho_F = 1.0;
		}

		// Longitud efectiva del Flange.
		w_efe = Rho_F * w;

		// Longitud efectiva del Flange en la parte superior.
		w_efe_1 = w_efe / (3.0-Psi);

		// Longitud efectiva del Flange en la parte inferior.
		if( Psi <= -0.236){
			w_efe_2 = w_efe / 2.0;
		}
		else{
			w_efe_2 = w_efe - w_efe_1;
		}

		// Longitud no efectiva del Flange.
		w_inefe = (CC-R-t) - (w_efe_1+w_efe_2);

		// Actualizacion de la suma de
		// longitudes efectivas en Flange.
		sum_efe_a = sum_efe_s;
		sum_efe_s = w_efe_1+w_efe_2;

		// Actualizacion del contador.
		conta = conta + 1;
	}

	// Correccion de la longitud no efectiva del
	// Flange o Ancho que no puede se negativa.
	if( w_inefe < 0.0){
		w_inefe = 0.0;
		w_efe_1 = 0.5 * w;
		w_efe_2 = 0.0;
	}

	// Areas de cada sub-seccion.
	AREAS[0]  = d * t;
	AREAS[1]  = ACR;
	AREAS[2]  = (w - w_efe_1 - w_inefe) * t;
	AREAS[3]  = w_efe_1 * t;
	AREAS[4]  = AREAS[1];
	AREAS[5]  = h_efe * t;
	AREAS[6]  = AREAS[4];
	AREAS[7]  = AREAS[3];
	AREAS[8]  = AREAS[2];
	AREAS[9]  = AREAS[1];
	AREAS[10] = AREAS[0];

	// Distancias de centroides de cada
	// sub-seccion a la fibra superior.
	DISTA[0]  = W - 0.5*t;
	DISTA[1]  = W - t - R + CCR;
	DISTA[2]  = R + t + 0.5*w + 0.5*w_efe_1 + 0.5*w_inefe;
	DISTA[3]  = R + t + 0.5*w_efe_1;
	DISTA[4]  = R + t - CCR;
	DISTA[5]  = 0.5*t;
	DISTA[6]  = DISTA[4];
	DISTA[7]  = DISTA[3];
	DISTA[8]  = DISTA[2];
	DISTA[9]  = DISTA[1];
	DISTA[10] = DISTA[0];

	// Distancia de fibra superior al eje neutro.
	// Momento de area total de la seccion.
	Mome = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		Mome = Mome + AREAS[i] * DISTA[i];
	}

	// Area total de la seccion.
	Area = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		Area = Area + AREAS[i];
	}

	// Distancia de la fibra superior al eje neutro.
	CC = Mome / Area;

	// Momentos de Inercia en Y respecto a
	// ejes locales de cada sub-seccion.
	MILOC[0]  = DOVO * t * t * t * d;
	MILOC[1]  = ICR;
	MILOC[2]  = DOVO * (w - w_efe_1 - w_inefe) * (w - w_efe_1 - w_inefe) * (w - w_efe_1 - w_inefe) * t;
	MILOC[3]  = DOVO * w_efe_1 * w_efe_1 * w_efe_1 * t;
	MILOC[4]  = MILOC[1];
	MILOC[5]  = DOVO * t * t * t * h_efe;
	MILOC[6]  = MILOC[4];
	MILOC[7]  = MILOC[3];
	MILOC[8]  = MILOC[2];
	MILOC[9]  = MILOC[1];
	MILOC[10] = MILOC[0];

	// Distancia paralela al eje Y del
	// centroide de la seccion a los ejes
	// inerciales locales de cada sub-seccion.
	DISTA[0]  = w + 2.0*R + 1.5*t - CC;
	DISTA[1]  = w + R + t - CC;
	DISTA[2]  = w + R + t - 0.5*(w - w_efe_1 - w_inefe) - CC;
	DISTA[3]  = R + t + 0.5*w_efe_1 - CC;
	DISTA[4]  = R + t - CC;
	DISTA[5]  = 0.5*t - CC;
	DISTA[6]  = DISTA[4];
	DISTA[7]  = DISTA[3];
	DISTA[8]  = DISTA[2];
	DISTA[9]  = DISTA[1];
	DISTA[10] = DISTA[0];

	// Momento de inercia de la seccion respecto al eje Y.
	Iy2 = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		Iy2 = Iy2 + MILOC[i] + DISTA[i] * DISTA[i] * AREAS[i];
	}

	// Distancias maximas de las fibras
	// extremas a los ejes centroidales.
	Cy = CC;
	if( Cy < fabs(fabs(W)-fabs(Cy))){
		Cy = fabs(fabs(W)-fabs(Cy));
	}

	// Modulo de Seccion calculado con la maxima
	// distancia del centroide a fibras extremas.
	Sy2 = Iy2 / Cy;

	// Momentos estaticos de area de las porciones de la seccion
	// que quedan arriba del eje neutro y que trabajan a compre-
	// sion respecto al eje neutro (solo los negativos contribui-
	// ran, los positivos trabajan a tension y no cuentan para
	// la evaluacion del esfuerzo cortante en la seccion).
	MELOC[0]  = (w + 2.0*R + 1.5*t - CC) * (d * t);
	MELOC[1]  = (w + R + t + CCR - CC) * ACR;
	MELOC[2]  = (- 0.5 * w_efe_2) * (w_efe_2 * t);
	MELOC[3]  = (R + t + 0.5*w_efe_1 - CC) * (w_efe_1 * t);
	MELOC[4]  = (R + t - CCR - CC) * ACR;
	MELOC[5]  = (0.5*t - CC) * (h_efe * t);
	MELOC[6]  = MELOC[4];
	MELOC[7]  = MELOC[3];
	MELOC[8]  = MELOC[2];
	MELOC[9]  = MELOC[1];
	MELOC[10] = MELOC[0];

	// Momento estatico de la seccion
	// por encima del eje neutro.
	Qy2 = 0.0;
	for ( i = 0 ; i < 11 ; i++){
		if( MELOC[i] < 0.0){
			Qy2 = Qy2 - MELOC[i];
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	//            Obtine los momentos flexionantes permisibles .               //
	/////////////////////////////////////////////////////////////////////////////



    *Mnx  = Sx  * Fcx;
    *Mny1 = Sy  * Fcy1;
    *Mny2 = Sy2 * Fcy2;

    //jacob. impresion provisional
    //printf("Mnx= %f   Sx=%f      Fcx= %f\n",*Mnx,Sx,Fcx);

    /*
    printf("Mny1= %f   Sy=%f      Fcx= %f\n",*Mny1,Sy,Fcy1);
    printf("Mny2= %f   Sy2=%f      Fcx= %f\n\n",*Mny2,Sy2,Fcy2);
    //*/

    //jacob.impresion provisional
    //if( ielem == 1 ) printf( "\nMnx = %f  Sx=%f  Fcx= %f\n", *Mnx, Sx, Fcx );

    return;

}
//¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑¬∑
void oama( double* K, double* P, double* Vx, double* Vy, double* Mx, double* My, double* Cb,
          double* Ctf, int tb, int ifcy, int ifcz, int ifdy, int ifdz, double L, double mfcy,
          double mfcz, double dfcy, double dfcz, double mfdy, double mfdz, double fxi, double fyi,
          double fzi, double fxf, double fyf, double fzf, double mxi, double myi, double mzi,
          double mxf, double myf, double mzf, int ielem, double A, double Fy )
{
	//-----------------------------------------------------------------------------//
	//                                                                             //
	// Esta rutuna obtiene las acciones mecanicas maximas y sus correspondientes   //
	// coeficientes en una barra de acero rolado en frio para su posterior dise~no //
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
	//    Ctf  =  Coeficiente de Flexion en Y.                                     //
	//                                                                             //
	//-----------------------------------------------------------------------------//


	/** Ejes propios de la sección **/

    // Mismos ejes para perfiles omega y C

	/*		 Y
     ^
     |				  +
     |				  +
     |				  +
     |				  +
     |	              +
     |               +
       + + + + + + + + + +
     +     |
     +		 |
     +		 |
     +		 |
     +		 |
     +		 |
     +		 |
     + - - - - - - - - - - - - - - - - - -> X
     +		 |
     +		 |
     +		 |
     +		 |
     +		 |
     +     |
     + + + + + + + + + +
     |			     +
     |		  		  +
     |   	          +
     |	              +
     |	              +
     |		          +
     */


	// ------------- DECLARACION DE VARIABLES LOCALES -------------- //
	double m1;   // Momento Mayor en extremo de barra.
	double m2;   // Momento Menor en extremo de barra.
	double m3;   // Momento Mayor dentro de barra.
	double Ma;   // Momento a un   cuarto  de barra.
	double Mb;   // Momento a un   medio   de barra.
	double Mc;   // Momento a tres cuartos de barra.
	double xvo;  // Distancia a la que el Cortante es cero.
	double P1;   // Pendiente antes de carga concentrada.
	double P2;   // Pendiente despues de carga concentrada.
	// ------------------------------------------------------------- //


    // Factor de Longitud Efectiva.
	switch( tb){
		case 1 : // Art-Art.
			*K = 1.00;
			break;
		case 2 : // Emp-Emp.
			*K = 0.65;
			break;
		case 3 : // Emp-Art.
			*K = 0.80;
			break;
        case 4 : // Art-Emp.
			*K = 0.80;
			break;
	}
	*K = 3.0; // TEMPORAL, SOLO PARA PRUEBAS

	// Fuerza Axial Maxima.
	*P = fabs( fxi );
	if( *P < fabs( fxf ) ){
		*P = fabs( fxf );
	}


	// Obtiene el maximo cortante en X.
	/*
	*Vx = fabs( fzi);
	if( *Vx < fabs( fyf)){
		*Vx = fabs( fyf);
	}
	*/
	*Vx = fabs( fyi );

	// Obtiene el maximo cortante en Y.
	/*
	*Vy = fabs ( fyi);
	if( *Vy < fabs( fzf ) ){
		*Vy = fabs( fzf );
	}
	*/
	*Vy = fabs( fzi );

	///////////////////////////////////////////////////////////////////
	//            COEFICIENTE DE FLEXION Y MOMENTO                   //
	//            Flexión sobre eje X de la sección                  //
	///////////////////////////////////////////////////////////////////
	// Momentos en extremos de barra.

    //jacob.modif se realiza la comparación con valores absolutos.

    //Antes:
    /*

     if( myi < myf){
     m1 = myi;
     m2 = myf;
     }
     else{
     m1 = myf;
     m2 = myi;
     }

     */

    //ahora:

    if( fabs( myi ) < fabs( myf ) ){
		m1 = myi;
		m2 = myf;
	}
	else{
		m1 = myf;
		m2 = myi;
	}

	// Coeficiente de Flexion.
	if( ifcz == 0){
		if( ifdz == 0){
			// No hay ni Fuerzas Concentradas ni Distribuidas.
			// Los momentos maximos estan en extremos de barra.
			m3 = 0.0;
			// Evalua momentos en los cuartos de la barra.

            //jacob.modif cambio de signos para myf para
            //tener un diagrama de momentos coherente. Esto
            //utilizando la convención de ejes locales para un
            //elemento en 3D.
            Ma = 0.25 * ( 3.0 * myi -       myf);
			Mb = 0.50 * (       myi -       myf);
			Mc = 0.25 * (       myi - 3.0 * myf);

		}
		else{
			// No hay Fuerzas Concentradas pero Si Distribuidas.
			// Distancia a la que el Cortante es cero.
			xvo = - fzi / mfdz;
			// Verifica que la distancia sea dentro de la barra.
			if( ( xvo > 0.0) && ( xvo < L)){
				// Momento Maximo dentro de la barra.
				m3 = myi - 0.5 * fzi * fzi / mfdz;
			}
			else{
				m3 = 0.0;
			}
			// Evalua momentos en los cuartos de la barra.
            //jacob.modif
			Ma = 0.03125 * mfdz * L * L + 0.25 * fzi * L + myi;
			Mb = 0.12500 * mfdz * L * L + 0.50 * fzi * L + myi;
			Mc = 0.28125 * mfdz * L * L + 0.75 * fzi * L + myi;
		}


	}
	else{
		if( ifdz == 0){
			// Si hay Fuerzas Concentradas pero No Distribuidas.
			// Evalua el Momento flexionante en el
			// punto de aplicacion de la carga.
            m3 = myi + fzi * dfcz;

            // Pendientes de las ecuaciones de momento.
			//jacob.modif cambio de signo a myf para ajustar diagramas de momentos
            P1 = ( m3 - myi) / (     dfcz);
			P2 = ( -myf - m3) / ( L - dfcz);
			// Evalua momentos en los cuartos de la barra.
			if( dfcz <= 0.25 * L){
				Ma = P2 * ( 0.25 * L - dfcz) + m3;
				Mb = P2 * ( 0.50 * L - dfcz) + m3;
				Mc = P2 * ( 0.75 * L - dfcz) + m3;
			}
			else if( ( dfcz > 0.25 * L) && ( dfcz <= 0.50 * L)){
				Ma = P1 * ( 0.25 * L) + myi;
				Mb = P2 * ( 0.50 * L - dfcz) + m3;
				Mc = P2 * ( 0.75 * L - dfcz) + m3;
			}
			else if( ( dfcz > 0.50 * L) && ( dfcz <= 0.75 * L)){
				Ma = P1 * ( 0.25 * L) + myi;
				Mb = P1 * ( 0.50 * L) + myi;
				Mc = P2 * ( 0.75 * L - dfcz) + m3;
			}
			else{
				Ma = P1 * ( 0.25 * L) + myi;
				Mb = P1 * ( 0.50 * L) + myi;
				Mc = P1 * ( 0.75 * L) + myi;
			}
		}
		else{
			// Si hay Fuerzas Concentradas y Distribuidas.
			// Evalua el Momento flexionante en el
			// punto de aplicacion de la carga.

			m3 = 0.5 * mfdz * dfcz * dfcz + fzi * dfcz + myi;

            // Evalua momentos en los cuartos de la barra.
			if( dfcz < 0.25 * L){
				Ma = myi + 0.25 * fzi * L + 0.03125 * mfdz * L * L + mfcz * ( 0.25 * L - dfcz );
				Mb = myi + 0.50 * fzi * L + 0.12500 * mfdz * L * L + mfcz * ( 0.50 * L - dfcz );
				Mc = myi + 0.75 * fzi * L + 0.28125 * mfdz * L * L + mfcz * ( 0.75 * L - dfcz );
			}
			else if( ( dfcz >= 0.25 * L) && ( dfcz < 0.50 * L)){
				Ma = myi + 0.25 * fzi * L + 0.03125 * mfdz * L * L;
				Mb = myi + 0.50 * fzi * L + 0.12500 * mfdz * L * L + mfcz * ( 0.50 * L - dfcz );
				Mc = myi + 0.75 * fzi * L + 0.28125 * mfdz * L * L + mfcz * ( 0.75 * L - dfcz );
			}
			else if( ( dfcz >= 0.50 * L) && ( dfcz < 0.75 * L)){
				Ma = myi + 0.25 * fzi * L + 0.03125 * mfdz * L * L;
				Mb = myi + 0.50 * fzi * L + 0.12500 * mfdz * L * L;
				Mc = myi + 0.75 * fzi * L + 0.28125 * mfdz * L * L + mfcz * ( 0.75 * L - dfcz );
			}
			else{
				Ma = myi + 0.25 * fzi * L + 0.03125 * mfdz * L * L;
				Mb = myi + 0.50 * fzi * L + 0.12500 * mfdz * L * L;
				Mc = myi + 0.75 * fzi * L + 0.28125 * mfdz * L * L;
			}
		}
	}
	// Obtiene valores absolutos.
	Ma = fabs( Ma );
	Mb = fabs( Mb );
	Mc = fabs( Mc );
	// Momento Maximo en la Barra.
	*Mx = fabs( m1 );
	if( *Mx < fabs( m2 ) ){
		*Mx = fabs( m2 );
	}
	if( *Mx < fabs( m3 ) ){
		*Mx = fabs( m3 );
	}

	// Coeficiente de Flexion.
	if( fabs( *Mx) < 0.00000000001 ){
		*Cb = 1.0;
	}
	else{
		*Cb = 12.5 * *Mx / ( 2.5 * *Mx + 3.0 * Ma + 4.0 * Mb + 3.0 * Mc );
	}

    // Flexocompresión o flexión
    // Se compara con el valor de aprox. el 10% de la resistencia nominal del perfil
	if ( *P >= 0.05 * A * Fy ) *Cb = 1.0;


    //jacob.impresion prov
    //if(ielem==176||ielem==157){printf("\nMx=%f  Ma=%f  Mb=%f  Mc= %f",*Mx, Ma, Mb, Mc);
	//printf("Distribuidas= %d   Concentradas = %d", ifdy, ifcy);}


	// ------------------------------------------------------------- //

	///////////////////////////////////////////////////////////////////
	//            COEFICIENTE DE FLEXION Y MOMENTO EN Y              //
	///////////////////////////////////////////////////////////////////
	// Momentos en extremos de barra.

    //jacob.modif se realiza la comparación con valores absolutos.

    //Antes:
    /*

	 if( mzi < mzf){
	 m1 = mzi;
	 m2 = mzf;
	 }
	 else{
	 m1 = mzf;
	 m2 = mzi;
	 }
     */

    //ahora:

    if( fabs(mzi) < fabs(mzf)){
		m1 = mzi;
		m2 = mzf;
	}
	else{
		m1 = mzf;
		m2 = mzi;
    }



    //jacob.impresion provisional
    //printf("mzi= %f    mzf=%f   m1=%f    m2=%f  \n",mzi,mzf,m1,m2);

	// Coeficiente de Flexion.
	if( ifcy == 0){
		if( ifdy == 0){
			// No hay ni Fuerzas Concentradas ni Distribuidas.
			// Los momentos maximos estan en extremos de barra.
			m3 = 0.0;
		}
		else{
			// No hay Fuerzas Concentradas pero Si Distribuidas.
			// Distancia a la que el Cortante es cero.
			xvo = - fyi / mfdy;
			// Verifica que la distancia sea dentro de la barra.
			if( ( xvo > 0.0) & ( xvo < L)){
				// Momento Maximo dentro de la barra.
				//jacob.modif
                //antes:
                //m3 = mzi - 0.5 * fyi * fyi / mfdy;
                //después
                m3 = mzi + 0.5 * fyi * fyi / mfdy;

				if(ielem==3){printf(" \n --L=%f  mfdy=%f,  fyi=%f,  fyf=%f,  mzi=%f, mzf=%f m3= %f-- \n",L,mfdy,fyi, fyf,mzi, mzf,m3);}
                if(ielem==3){printf(" \n --Longa=%f-- \n",xvo);}


			}
			else{
				m3 = 0.0;
			}
		}
	}
	else{
		if( ifdy == 0){
			// Si hay Fuerzas Concentradas pero No Distribuidas.
			// Evalua el Momento flexionante en el
			// punto de aplicacion de la carga.
			//jacob.modif
            //antes:
            //m3 = mzi + fyi * dfcy;
            //después
            m3 = mzi - fyi * dfcy;

		}
		else{
			// Si hay Fuerzas Concentradas y Distribuidas.
			// Evalua el Momento flexionante en el
			// punto de aplicacion de la carga.
			//jacob.modif
            //antes
            //m3 = 0.5 * mfdy * dfcy * dfcy + mfcy * dfcy + mzi;
            //después:
            m3 = - 0.5 * mfdy * dfcy * dfcy - fyi * dfcy + mzi;

        }
	}

    // En este caso m3 > m1 y m3 > m2. Hay un momento mayor en medio dentro
    // de la barra
	if( fabs( m3 ) > fabs( m1 ) && fabs( m3 ) > fabs( m2 ) ) *Ctf = 1.0;
	else{
        // En este caso m2 = 0, y por lo tanto m1 = m3 = 0
        if( fabs( m2 ) < 0.000001 ) *Ctf = 0.6;
        // En este caso m2 no es cero
        else                        *Ctf = 0.6 - 0.4 * ( m1 / m2 );
	}

	// Flexocompresión o flexión
    // Se compara con el valor de aprox. el 10% de la resistencia nominal del perfil
	if ( *P >= 0.05 * A * Fy ) *Ctf = 1.0;

    //jacob.impresion
    //printf("  m1=%f   m2=%f   Ctf=%f \n",m1,m2,*Ctf);

	// Momento Maximo en la Barra.
	*My = fabs( m1);
	if( *My < fabs( m2)){
		*My = fabs( m2);
	}
	if( *My < fabs( m3)){
		*My = fabs( m3);
	}
	// ------------------------------------------------------------- //

}

double mom_pandeo_distorsional_flex( double E, double Fy, double mu, double Ix,
                                     double t, double ho, double bo, double D ){

    // Calcula el momento resistente nominal al pandeo distorsional para
    // el caso de FLEXIÓN

    double Fd, Sf, Mcrd, My, lambda_d, Mn;

    // Esfuerzo por pandeo distorsional
    Fd = esf_pandeo_distorsional_flex( E, mu, t, ho, bo, D );

    // Modulo de seccion con respecto a la fibra superior a compresion,
    // de la seccion completa
    Sf = Ix / ( 0.5 * ho );

    // Mcr
    Mcrd = Sf * Fd;

    // En este casoo Sfy = Sf
    My = Sf * Fy;

    // Momento nominal resisntente al pandeo distorsional
    lambda_d = sqrt( My / Mcrd );
    if( lambda_d <= 0.673 ) Mn = My;
    else                    Mn = ( 1.0 - 0.22*sqrt(Mcrd/My) ) * sqrt(Mcrd/My) * My;

    return Mn;

}

double esf_pandeo_distorsional_flex( double E, double mu, double t, double ho,
                                     double bo, double D ){

    // Calcula el esfuerzo por pandeo distorsional ( "distortional
    // buckling stress Fd" ) para el caso de FLEXION

    double G, XIweb;
    double pi, b, d, Af, xo, yo, hx, Jf, Cwf, Ixf, Iyf, Ixyf;
    double beta, Lcr, L, k_phi_fe, k_phi_we, k_phi, k_phi_fg, k_phi_wg, Fd;
    pi = 3.14159265359;
    G  = E / ( 2.0 * ( 1.0 + mu ) );
    // Dado que esto aplica solo cuando la flexion es sobre el eje
    // fuerte (hay simetria)
    XIweb = 2.0;

    // Propiedades (seccion C)
    //h  = ho - t;
    b  = bo - t;
    d  = D - 0.5 * t;
    Af = ( b + d ) * t;
    xo =  b*b / (2.0*(b+d));
    yo = -d*d / (2.0*(b+d));
    hx = -(b*b+2.0*d*b) / (2.0*(b+d));
    Jf = (b*t*t*t+d*t*t*t) / 3.0;
    Cwf  = 0.0;
    Ixf  = t * ( t*t*b*b + 4.0*b*d*d*d + t*t*b*d + d*d*d*d ) / (12.0*(b+d));
    Iyf  = t * ( b*b*b*b + 4.0*b*b*b*d ) / (12.0*(b+d));
    Ixyf = t*b*b*d*d / (4.0*(b+d));

    // Se toma beta = 1, ignorando el gradiente de momentos
    beta = 1.0;

    // Longitud L
        Lcr = pow( 4.0*pi*pi*pi*pi*ho*(1.0-mu*mu)/(t*t*t) * ( Ixf*(xo-hx)*(xo-hx)
                   + Cwf - Ixyf*Ixyf*(xo-hx)*(xo-hx)/Iyf )
                   + pi*pi*pi*pi*ho*ho*ho*ho/720.0, 0.25 );

        // Dado que lo mas probable es que no se coloquen elementos que
        // eviten el pandeo distorsional se hace la siguiente considerecion
        //Lm = Lcr;
    L = Lcr;

    // Rigidez rotacional elastica que proporciona el lomo (flange) a la
    // junta entre el lomo y el alma (web): k_phi_fe
    k_phi_fe = pi*pi*pi*pi/(L*L*L*L) * ( E*Ixf*(xo-hx)*(xo-hx) + E*Cwf
               - E*Ixyf*Ixyf/Iyf*(xo-hx)*(xo-hx) ) + pi*pi/(L*L)*G*Jf;

    // Rigidez rotacional elastica que propociona el alma (web) a la junta
    // entre el lomo (flange) y el alma: k_phi_we
    k_phi_we = E*t*t*t/(12.0*(1.0-mu*mu)) * ( 3.0/ho + pi*pi/(L*L)*19.0*ho/60.0
               + pi*pi*pi*pi/(L*L*L*L)*ho*ho*ho/240.0 );

    // Rigidez rotacional elastica proporcionada por algun elemento
    // que de soporte a la junta entre el lomo y el alma ("brace","panel",
    // "sheating"). En este caso se considera nula.
    k_phi = 0.0;

    // Rigidez rotacional geometrica que se requiere que proporcione el lomo
    // a la junta lomo-web
    k_phi_fg = pi*pi/(L*L) * ( Af * ( (xo-hx)*(xo-hx)*Ixyf*Ixyf/(Iyf*Iyf)
               - 2.0*yo*(xo-hx)*Ixyf/Iyf + hx*hx + yo*yo ) + Ixf + Iyf );

    // Rigidez rotacional geometrica que se requiere que proporcione el alma
    // a la junta lomo-web
    k_phi_wg = ho*t*pi*pi/13440.0 * ( (45360.0*(1.0-XIweb)+62160.0)*L*L/(ho*ho) +
               448.0*pi*pi + ho*ho/(L*L)*(53.0+3.0*(1.0-XIweb))*pi*pi*pi*pi ) /
               ( pi*pi*pi*pi + 28.0*pi*pi*L*L/(ho*ho) + 420.0*L*L*L*L/(ho*ho*ho*ho) );

    Fd = beta * ( k_phi_fe + k_phi_we + k_phi ) / ( k_phi_fg + k_phi_wg );

    return Fd;

}

double fuer_pandeo_distorsional_comp( double E, double Fy, double mu, double Ag,
                                      double t, double ho, double bo, double D ){

    // Calcula la resistencia a la COMPRESION por pandeo distorsinal

    double Fd, Pcrd, Py, lambda_d, Pn;

    // Esfuerzo de pandeo distorsional
    Fd = esf_pandeo_distorsional_comp( E, mu, t, ho, bo, D );

    // Pcrd
    Pcrd = Ag * Fd;

    // Py
    Py = Ag * Fy;

    // Resistencia nominal
    lambda_d = sqrt( Py / Pcrd );
    if( lambda_d <= 0.561 ) Pn = Py;
    else                    Pn = ( 1.0 - 0.25*pow(Pcrd/Py,0.6) ) * pow(Pcrd/Py,0.6) * Py;

    return Pn;

}


double esf_pandeo_distorsional_comp( double E, double mu, double t, double ho,
                                     double bo, double D ){

    // Calcula el esfuerzo por pandeo distorsional ( "distortional
    // buckling stress Fd" ) para el caso de FLEXION

    double G, XIweb;
    double pi, b, d, Af, xo, yo, hx, Jf, Cwf, Ixf, Iyf, Ixyf;
    double Lcr, L, k_phi_fe, k_phi_we, k_phi, k_phi_fg, k_phi_wg, Fd;
    pi = 3.14159265359;
    G  = E / ( 2.0 * ( 1.0 + mu ) );
    // Dado que esto aplica solo cuando la flexion es sobre el eje
    // fuerte (hay simetria)
    XIweb = 2.0;

    // Propiedades (seccion C)
    //h  = ho - t;
    b  = bo - t;
    d  = D - 0.5 * t;
    Af = ( b + d ) * t;
    xo =  b*b / (2.0*(b+d));
    yo = -d*d / (2.0*(b+d));
    hx = -(b*b+2.0*d*b) / (2.0*(b+d));
    Jf = (b*t*t*t+d*t*t*t) / 3.0;
    Cwf  = 0.0;
    Ixf  = t * ( t*t*b*b + 4.0*b*d*d*d + t*t*b*d + d*d*d*d ) / (12.0*(b+d));
    Iyf  = t * ( b*b*b*b + 4.0*b*b*b*d ) / (12.0*(b+d));
    Ixyf = t*b*b*d*d / (4.0*(b+d));

    // Longitud L
        Lcr = pow( 6.0*pi*pi*pi*pi*ho*(1.0-mu*mu)/(t*t*t) * ( Ixf*(xo-hx)*(xo-hx)
                   + Cwf - Ixyf*Ixyf*(xo-hx)*(xo-hx)/Iyf ), 0.25 );

        // Dado que lo mas probable es que no se coloquen elementos que
        // eviten el pandeo distorsional se hace la siguiente considerecion
        //Lm = Lcr;
    L = Lcr;

    // Rigidez rotacional elastica que proporciona el lomo (flange) a la
    // junta entre el lomo y el alma (web): k_phi_fe
    k_phi_fe = pi*pi*pi*pi/(L*L*L*L) * ( E*Ixf*(xo-hx)*(xo-hx) + E*Cwf
               - E*Ixyf*Ixyf/Iyf*(xo-hx)*(xo-hx) ) + pi*pi/(L*L)*G*Jf;

    // Rigidez rotacional elastica que propociona el alma (web) a la junta
    // entre el lomo (flange) y el alma: k_phi_we
    k_phi_we = E*t*t*t  / ( 6.0*ho*(1.0-mu*mu) );

    // Rigidez rotacional elastica proporcionada por algun elemento
    // que de soporte a la junta entre el lomo y el alma ("brace","panel",
    // "sheating"). En este caso se considera nula.
    k_phi = 0.0;

    // Rigidez rotacional geometrica que se requiere que proporcione el lomo
    // a la junta lomo-web
    k_phi_fg = pi*pi/(L*L) * ( Af * ( (xo-hx)*(xo-hx)*Ixyf*Ixyf/(Iyf*Iyf)
               - 2.0*yo*(xo-hx)*Ixyf/Iyf + hx*hx + yo*yo ) + Ixf + Iyf );

    // Rigidez rotacional geometrica que se requiere que proporcione el alma
    // a la junta lomo-web
    k_phi_wg = pi*pi/(L*L) * t*ho*ho*ho/60.0;

    Fd = ( k_phi_fe + k_phi_we + k_phi ) / ( k_phi_fg + k_phi_wg );

    return Fd;

}

