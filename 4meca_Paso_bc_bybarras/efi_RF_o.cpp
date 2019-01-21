//*****************************
// Rutina para Eficiencias
//*****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "efi_RF_o.h"
#include "eficiencias.h"
#include "raros.h"

//************************************************
//************************************************

// Rutinas de rolado en frio
//------------------------------------------------------------------------------
double f2eso(double H, double W, double D, double R, double t, double Fy,
             double E, double LongX, double LongY, double LongT, int tb,
             int ifcz, int ifcy, int ifdz, int ifdy, double mfcz, double mfcy,
             double dfcz, double dfcy, double mfdz, double mfdy, double fxi,
             double fyi, double fzi, double fxf,double fyf, double fzf,
             double mxi, double myi, double mzi, double mxf, double myf,
             double mzf, int ielem)
{

    //printf ( "\nENTRA EFICIENCIAS OMEGAS" );
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
	double Sf;         // Modulo de Seccion.
	double Cs;         // Coeficiente de esfuerzo critico de pandeo.
	double Mnx;        // Momento Nominal en X.
	double Mny;        // Momento Nominal en Y.
    double Mnxt;	   // Momento Nominal de la sección completa = Sftx * Fy
    double Mnyt;	   // Momento Nominal de la sección completa = Sfty * Fy
	double Mny1;       // Momento Nominal en Y1.
	double Mny2;       // Momento Nominal en Y2.
	double Vn;         // Fuerza Cortante Nominal.
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
	double A_aux;
	double xcg_aux;
	double Ix_aux;
	double Iy_aux;
	//double SF, SM;
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

    //Ejemplo de la pag 158 del Hancock
    /*
     H     =  8;
     W = 2.75;
     D = 0.625;
     R = .1875;
     t = 0.060;
     Fy = 55;
     //*/

    //Ejemplo de la pag 399 de Yu
    /*
    H =	8;
    W = 	3;
    D = 	0.115;
    R = 	3.0/16.0;
    t = 	0.105;
    Fy = 	50;
    */


	//-------------------------------------------------------------------------//
	//             PROPIEDADES GEOMETRICAS DE LA SECCION COMPLETA              //
	//                  CALCULADAS POR EL METODO DE LINEAS                     //
	//-------------------------------------------------------------------------//
	// Propiedades de Cuarto de Rondana.
	r = R + 0.5 * t;   // Radio medio.
	u = 1.570 * r;     // Longitud de arco.
	c = 0.637 * r;     // Centroide.



	// Longitudes planas de la seccion.
	h = H - 2.0 * ( R + t );   // Peralte.
	w = W - 2.0 * ( R + t );   // Flange.
	d = D - 1.0 * ( R + t );   // Labio.

	// Area de la seccion completa.
	A = t * ( 2.0 * w + h + 4.0 * u + 2.0 * d);
    //A= ( aa + 2* bb + 2 * cc ) * t;

	// Momento de inercia en X, para sección Omega
	Ix =  	t * (	2.0 * u * ( 0.5 * h + c ) * ( 0.5 * h + c) +
                 2.0 * u * ( 0.5 * H + R - c ) * ( 0.5 * H + R - c )) +
    2 * t * w * ( 0.5 * H - 0.5 * t ) * ( 0.5 * H - 0.5 * t ) +
    2 * t * d * ( 0.5 * H + R + 0.5 * d ) * ( 0.5 * H + R + 0.5 * d ) +
    DOVO * h * h * h * t;

	// Abscisa del centroide de la seccion medido
	// a partir de la parte EXTERNA del peralte.
    //   __
    //  |  |   La parte externa es del lado donde esta el punto.
    // .|
    //  |
    //  |  |
    //   --
	xcg = ( 0.5 * t * h + 2.0 * u * ( t + R - c ) +
		   2.0 * u * ( W - R - t + c ) +
		   w * W + 2.0 * d * ( W - 0.5 * t ) ) * t / A;

	// Momento de inercia en Y.
	Iy = ( 0.25 * t * t * h + 2.0 * u * ( t + R - c ) * ( t + R - c ) +
		  2.0 * u * ( W - R - t + c ) * ( W - R - t + c ) + 0.5 * W * W * w +
		  2.0 * d * ( W - 0.5 * t ) * ( W - 0.5 * t ) ) * t +
	SEXT * w * w * w * t - A * xcg * xcg;

	// Momento polar de inercia.
	// constante torsional de St. Venant, p 196 y 631 Yu
    // jacob.revisar
    J = TERC * A * t * t;

	// Centro de cortante y constantes para evaluar el alabeo.
	// p160 Hancock
    aa = H - t;
	bb = W - t;
	cc = D - 0.5 * t;
	Ms = ( bb * ( 3.0 * aa * aa * bb + cc * ( 6.0 * aa * aa - 8.0 * cc * cc ) ) ) /
    ( aa * aa * aa + 6.0 * aa * aa * bb +
     cc * ( 8.0 * cc * cc - 12.0 * aa * cc + 6.0 * aa * aa ) );

	// Abscisa del centro de cortante de la seccion completa.
	xco = - ( Ms - 0.5 * t ) - xcg;

	// Constante de alabeo.
	Cw = DOVO * aa * aa * bb * bb * t * ( ( 2.0 * aa * aa * aa * bb +
										   003.0 * aa * aa * bb * bb + 48.0 * cc * cc * cc * cc +
                                           112.0 * bb * cc * cc * cc + 8.0 * aa * cc * cc * cc +
                                           48.0 * aa * bb * cc * cc + 12.0 * aa * aa * cc * cc +
                                           12.0 * aa * aa * bb * cc + 6.0 * aa * aa * aa * cc ) /
										 ( 6.0 * aa * aa * bb + ( aa + 2.0 * cc ) * ( aa + 2.0 * cc ) *
										  ( aa + 2.0 * cc ) - 24.0 * aa * cc * cc ) );
    //jess.provisional impresion
    //printf ( "\n Con Metodología y Fórmulas del Hancock");
    //printf("\nA= %f, \t xcg= %f, \t Ix= %f, \t M= %f, \t xc0 = %f, \t Cw= %f",A, xcg, Ix, Ms, xco, Cw );


    // Fórmulas para Cw de una sección C utilizando tabla B.1 de Yu.
    // Sin embargo para secciones C daban resultados de aprox. el
    // Doble de lo que da con formula del Hancock. Comparado con ejercicio
    // del mismo Hancock. Pag. 159.
    /*
     Ix = t / 12.0 * ( aa* aa * aa + 6 * bb * aa * aa + 6 * cc * aa * aa
     - 12.0 * aa * cc * cc
     + 8 * cc * cc * cc);

     Ms = bb * t / (12.0 * Ix) * (  6 * cc * aa * aa
     + 3 * bb * aa * aa
     - 8 * cc * cc * cc);

     Cw = t * t / A * ( xcg * A * aa * aa / t * ( bb * bb / 3.0  + Ms * Ms - Ms * bb)
     + A / (3.0 * t ) * ( Ms * Ms * aa * aa * aa + bb* bb * cc * cc * ( 2 * cc + 3.0 * aa))
     - Ix * Ms * Ms / t * (2.0 * aa * + 4.0 * cc)
     + Ms * cc * cc / 3.0 * ( 8 * bb * bb *cc + 2 * Ms * ( 2 * cc *( cc - aa ) + bb * ( 2 * cc - 3 * aa )))
     + bb * bb * aa * aa / 2.0 * ( ( 3 * cc + bb ) * ( 4 * cc + aa ) - 6 * cc * cc )
     - Ms * Ms * aa * aa * aa * aa / 4.0);

     //*/

    // Para sección Omega
    // El recálculo de A, xcg, Ix y Iy se realiza unicamente para el cálculo
    // de Cw. de acuerdo a las fórmulas de la tabla B.1 pag 629 del Yu.

    A_aux= ( aa + 2* bb + 2 * cc ) * t;
    xcg_aux = bb * t * ( bb + 2*cc) /A_aux;


    Ix_aux = t / 12.0 * ( aa* aa * aa + 6 * bb * aa * aa + 6 * cc * aa * aa
                         + 12.0 * aa * cc * cc
                         + 8 * cc * cc * cc);

    //xco_aux = bb*t*(bb+2*cc)/A_aux + bb * t / (12 * Ix_aux) * (6 * cc * aa * aa + 3* bb * aa * aa
    //                                                           - 8 * cc * cc * cc);

    Iy_aux = t * bb * bb / ( 3.0 * ( aa + 2.0 * bb + 2.0 * cc)) *
    ( 2 * aa * bb + bb * bb + 4.0 * bb * cc + 6 * cc * aa );

    Cw = aa * aa / 4.0 * ( Iy_aux + xcg_aux * xcg_aux * A_aux * ( 1 - aa * aa * A_aux / (4.0 * Ix_aux)))
    + 2.0 * bb * bb * t * cc * cc * cc / 3.0
    - aa * bb * bb * cc * cc * t
    + ( aa * aa * bb * t * cc * cc * cc * xcg * A_aux ) / (3.0 * Ix_aux)
    - ( 4 * bb * bb * t * t * cc*cc*cc * cc*cc*cc ) / ( 9.0 * Ix_aux);



	//jess.provisional impresion
    //printf ( "\n\n Con Metodología y Fórmulas del Yu");
    //printf("\nA= %f, \t xcg= %f, \t Ix= %f, \t M= %f, \t xc0 = %f, \t Cw= %f",A_aux, xcg_aux, Ix_aux, Ms, xco, Cw );



	// Radios de giro.
	// jacob.comentario p161 Hancock
    rx = sqrt( Ix / A );
	ry = sqrt( Iy / A );
	ro = sqrt( ( Ix + Iy ) / A + xco * xco );


	// Distancia del lomo a la coordenada centroidal.
    // jacob.comentario en realidad al parecer es del eje del alma al centro de
    // gravedad, considerando dimensiones como en
    // pag 159 delHancock fig 5.20 (b) ; formula en pag. 645 Appendix C lipped channel.
	xbar = bb * ( bb + 2.0 * cc ) / ( aa + 2.0 * bb + 2.0 * cc );

	// Distancia de los centros de cortante al de gravedad.
	xo = - Ms - xbar;

	// Constantes de evaluacion de propiedad j.
    // jacob.comentario Xbar debe ser considerada negativa por
    // eso signos negativos.
    // Appendix C del Yu
	betaw = - DOVO * t * xbar * aa * aa * aa - t * xbar * xbar * xbar * aa;
	tempo = bb - xbar;
	betaf = 0.5 * t * ( tempo * tempo * tempo * tempo -
					   xbar * xbar * xbar * xbar ) + 0.25 * aa * aa * t *
	( tempo * tempo - xbar * xbar );
	beta1 = 2.0 * cc * t * tempo * tempo * tempo + (2.0/3.0) * t *
	tempo * ( ( aa / 2.0 + cc ) * ( aa / 2.0 + cc ) * ( aa / 2.0 + cc )
             - ( aa * aa * aa / 8.0 ));
	betay = ( ( betaw + betaf + beta1 ) / Iy ) - 2.0 * xo;

	// Propiedad de Pandeo de Torsion-Flexion.
	// Ec. 6.47 pag 375 del Yu
    jtofl = 0.5 * betay;

    //-------------------------------------------------------------------------//
	//      OBTIENE ACCIONES MECANICAS ACTUANTES Y COEFICIENTES DE APOYO       //
	//-------------------------------------------------------------------------//
	oama( &K, &P, &Vx, &Vy, &Mx, &My, &Cb, &Ctf, tb, ifcy, ifcz, ifdy, ifdz, LongT, mfcy, mfcz, dfcy, dfcz, mfdy, mfdz, fxi, fyi, fzi,
		 fxf, fyf, fzf, mxi, myi, mzi, mxf, myf, mzf, ielem, A, Fy );

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
    // El valor de u ( coeficiente de Poisson )  se tomó igual a 0.3 .
    SigmaT = ( E / ( A * ro * ro ) ) *
             ( 0.38461538461538 * J + PICU * Cw / ( K * K * LongT * LongT ) );


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



	//-------------------------------------------------------------------------//

    //-------------------------------------------------------------------------//
	//          DETERMINACION DEL ESFUERZO CRITICO DE TORSION LATERAL          //
	//               CUANDO LA FLEXION ACTUA EN EL EJE X (Fcx)                 //
	//-------------------------------------------------------------------------//
    /*          	Y
                    ^
                    |				  +		^
                    |				  +		|
                    |				  +		d  labio
                    |				  +		|
                    |	ala	   		  +		v
                <- - - - w - - - -> +
                + + + + + + + + + +
              +
     ^		+		|
     |		+		|
     |		+		|
     |		+		|
     |		+		|
     |		+		|
     h pe-  + - - - - - - - - - - - - - - - - - -> X
     | ralte+		|
     |		+		|
     |		+		|
     |		+		|
     v		+		|
              +
                + + + + + + + + + +
                    |			    +
                    |		  		  +
                    |				  +
                    |				  +
                    |				  +
                    |				  +
     */
	// Inciso a) de articula C3.1.2.1 del AISI. Pagina 207 de Yu.              //
	//-------------------------------------------------------------------------//
	// Modulo de seccion elastico de la seccion completa.
	Sf = Ix / ( ( H +  2 * ( R + d ) ) / 2.0);

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
	/*

	                                  X
	                                  ^
	         + + + +                  |                  + + + +
	                 +                |                +
	                   +              |              +
	                   +              |              +
	                   +              |              +
            Y <- - - - + - - - - - - -|- - - - - - - + - - - -
	                   +              |              +
	                   +              |              +
	                   +              |              +
	                     +            |            +
	                       + + + + + +|+ + + + + +
                                      |

	*/
	// Inciso b) de articula C3.1.2.1 del AISI. Pagina 207 de Yu.              //
	//-------------------------------------------------------------------------//
	// Modulo de seccion elatico de la seccion completa.
	Sf = Iy / ( W - xcg );
	Cs = -1.0;

	// Esfuerzo critico de pandeo lateral de torsion.
	FE = ( ( Cs * A * SigmaX) / ( Ctf * Sf)) * ( jtofl + Cs *
                                                sqrt( jtofl * jtofl + ro * ro * ( SigmaT / SigmaX)));
    FE = ( ( Cs * A * SigmaX) / Sf) * ( jtofl + Cs *
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
	/*

	                                  |
	                       + + + + + +|+ + + + + +
	                     +            |            +
	                   +              |              +
	                   +              |              +
	                   +              |              +
               - - - - + - - - - - - -|- - - - - - - + - - - -> Y
	                   +              |              +
	                   +              |              +
	                   +              |              +
                     +                |                +
             + + + +                  |                  + + + +
                                      |
                                      v
                                      X

	*/                                                        //
	// Inciso b) de articula C3.1.2.1 del AISI. Pagina 207 de Yu.
	//-------------------------------------------------------------------------//
	// Modulo de seccion elatico de la seccion completa.
	Sf = Iy / xcg;
	Cs = 1.0;

	// Esfuerzo critico de pandeo lateral de torsion.
	// Antes
	/*
	FE = ( ( Cs * A * SigmaX) / Sf) * ( jtofl + Cs *
                                       sqrt( jtofl * jtofl + ro * ro * ( SigmaT / SigmaX)));
    */
    // Despues
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
	//      CRITICOS DE PANDEO LATERAL DE TORSION ACTUANDO EN AMBOS EJES       //
	//-------------------------------------------------------------------------//
	// Funcion Momentos Nominales a Flexion de Seccion C.

    // Funcion Empleada para sección C:
   	// mnfsc( H, W, D, R, t, E, Fcx, Fcy1, Fcy2, &Mnx, &Mny1, &Mny2 , ielem);

    mnfso ( H, W, D, R, t, E, Fcx, Fcy1, Fcy2, &Mnx, &Mny1, &Mny2 );

	// Obtiene el menor de los momentos permisibles en Y.
	Mny = Mny1;
	if( Mny > Mny2){
		Mny = Mny2;
	}
	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//                   OBTIENE FUERZA CORTANTE PERMISIBLE                    //
	//-------------------------------------------------------------------------//
	// Resistencia Nominal a Cortante cuando hay combinacion con flexion
	// segun capitulo 9 de STRUCTURAL STEEL DESIGNER'S HANDBOOK 4ed.
	// pag 170 del Hancock
    Vn = 0.64 * t * t * sqrt( 5.34 * Fy * E);

	if( (h/t) < sqrt( 5.34 * E / Fy)){
		Vn = 0.6 * Fy * h * t;
	}
	if( (h/t) > ( 1.51 * sqrt( 5.34 * E / Fy))){
		Vn = 0.904 * 5.34 * E * t * t * t / h;
	}
	//-------------------------------------------------------------------------//

	//-------------------------------------------------------------------------//
	//                        EFICIENCIA DE LA SECCION                         //
	//-------------------------------------------------------------------------//
	// Factores de seguridad ASD.
	omega_c = 1.80;   // Compresion.
	omega_b = 1.67;   // Flexion.
	omega_v = 1.67;   // Corte.
    omega_t = 1.67;		//Tension




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
    	EFIFT =   omega_b * ( Mx / Mnxt + My / Mnyt) + omega_t * P / (A * Fy);

        EFIFT2 =  omega_b * ( Mx / Mnx + My / Mny ) - omega_t * P / (A * Fy);

    	if ( EFIFT2 < EFIFT){
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
            ax  = 1.0 - omega_c * P / Pex;
            ay  = 1.0 - omega_c * P / Pey;
        }

        // Eficiencia de la seccion.
        EFIFC = omega_c * P / Pn + omega_b * ( ( Mx / ( ax * Mnx ) ) + ( My / ( ay * Mny ) ) );

        //jacob.modif cambio de eficiencia a 10 cuando ay o ax sea menor a 1
        if(ay<0||ax<0){
            //printf("Elemento %d \n\n\n\n",ielem);
            EFIFC=10;
        }


    }

	// Eficiencia de la seccion a Flexion y Cortante.
	//jacob.modif
    //se le saca raiz a la suma de los cuadrados de las relaciones Actuante/nominal
    // En realidad no importa mucho la reaiz por 0 < efi < 1
    EFIFV = ( omega_b * Mx / Mnx) * ( omega_b * Mx / Mnx) +
	( omega_b * My / Mny) * ( omega_b * My / Mny) +
	( omega_v * Vx / Vn ) * ( omega_v * Vx / Vn );

    EFIFV=sqrt(EFIFV);
	// EFIciencia final de la seccion.
	if( EFIFC < EFIFV){
		EFI = EFIFV;
	}
	else{
		EFI = EFIFC;
	}

    printf ( "\nEFI_OMEGA = %lf", EFI );
	return EFI;
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void mnfso( double H, double W, double D, double R, double t, double E, double Fcx,
            double Fcy1, double Fcy2,  double *Mnx, double *Mny1, double *Mny2 )
{
	/** NOTAS **/
	// Esta funcion determina el MOMENTO NOMINAL que resiste una
	// seccion dada en las cuatro posiciones posibles de esta a la
	// flexiÛn.
	// Este momento se determina en base en lo descrito en el libro
	// Wei-Wen Yu p·g. Inicia en p·g.147.
	// El momento se calcula a partir de INITIATION OF YIELDING.
	// El signo de los  momentos est· determinado de acuerdo a la
	// regla de la mano derecha.
	/** Variables de entrada **/
	// Fy
	// E
	// COMPLETAR

	/** Variables de salida **/
	// COMPLETAR

	int conta, i;
	double Lam_H, Rho_H, h_efe, w_efe_1, w_efe_2, w_inefe;
	double Lam_L, Rho_L, Ku, d_efe, Kw;
	double K, n, Ia, Is, C2, Ka, esf_w, S, h_efe_1, h_efe_2, h_inefe;
	double Lam_W, Rho_W;
	double sum_efe_a, sum_efe_s, AREAS[11], DISTA[11], MILOC[11];
	double Mome, Area, CC, f1, f2, Psi, tempo, Kf, Lam_F, Rho_F;
	double w_efe, esf_h, Cy;
	double Sy2, Iy2, Sy1, Iy1, Ixx, Sx;
	double w, h, d, a, b;
	double ACR, CCR, ICR;
	double DOVO = 1.0 / 12.0;

	/** Ejes propios de la secciÛn **/

	/*				Y
                    ^
                    |				  +		^
                    |				  +		|
                    |				  +		d  labio
                    |				  +		|
                    |	ala	   		  +		v
                <- - - - w - - - -> +
                + + + + + + + + + +
              +
     ^		+		|
     |		+		|
     |		+		|
     |		+		|
     |		+		|
     |		+		|
     h pe-  + - - - - - - - - - - - - - - - - - -> X
     | ralte+		|
     |		+		|
     |		+		|
     |		+		|
     v		+		|
              +
                + + + + + + + + + +
                    |			    +
                    |		  		  +
                    |				  +
                    |				  +
                    |				  +
                    |				  +
     */


	/**********************************
     FLEXI”N SOBRE EJE X
     **********************************/

	/** NOTA **/
	// En este caso se analiza de la misma manera los dos posibles casos,
	// dado que la seccion es simetrica respecto al eje X y por tanto la
	// direccion del momento no afecta.
	/* Caso 1						Caso 2

         |	 +					 	|	+
      + +|+ +					 + +|+ +
     +	 |						+	|
     +	 |						+ 	|
   --+---|----->> Mx	 Mx <<--+---|-----
     +	 |						+ 	|
     +	 |						+	|
      + +|+ + 					 + +|+ +
         |	 +						|	+

     Momento pos. causa com-		Momento neg. causa ten-
     presion en el lado infe-	    sion en el lado inferior
     rior de la seccion.			de la seccion.
     }*/

	/**		1 Propiedades de la seccion		**/

	/** 1.1 Cuarto de rondana, peralte, ala y labio **/

    /** Variables locales **/
    // Para un cuarto de rondana:
    //   r		Radio medio
    //   u		Longitud de arco
    //   c		Centroide ( mÈtodo de lÌneas )
    //   ACR	¡rea
    //   CCR	Centroide
    //	 ICR 	Momento de inercia
    // h		Peralte interno
    // w		Ala interna
    // D		Labio interno

	//r = R + 0.5 * t;
	//u = 1.570 * r;
	//c = 0.637 * r;

	// Los ejes para el c·lculo del centroide y momento de inercia no est·n sobre el
	// centroide. Estos contienen a todo el cuarto de rondana. El archivo de anota-
	// ciones contiene la orientaciÛn de los ejes para estos c·lculos.
	ACR  = 0.7853981633974483 * ( 2.0*R*t + t*t );
	CCR  = 0.4244131815783876 * ( (3.0*R*R + 3.0*R*t + t*t ) / ( 2.0*R + t ) );
	ICR = 0.1963495408493621 * ( 4.0*R*R*R*t + 6.0*R*R*t*t + 4.0*R*t*t*t + t*t*t*t);

	// Longitudes internas
	h = H - 2.0 * ( R + t ); // Peralte
	w = W - 2.0 * ( R + t ); // Ala
	d = D - 1.0 * ( R + t ); // Labio

	/** 1.2 Ancho efectivo del ala ( w ) que esta compresion **/

    /** Variables **/
    // n			Constante
    // Ia			Momento de inercia adecuado del labio
    // Is			Momento de inercia completo del labio
    // C1, C2		Coeficientes
    // Lam_F		Factor de esbeltez del ala a compresiÛn (w)
    // Rho_F		Factor de reduccion del ala a compresiÛn <= 1
    // w_efe		Longitud efectiva del ala en compresion uniforme
    //				En este caso sobre el ala act˙a un esfuerzo
    //				de compresiÛn constante
    // Lam_L     	Factor de esbeltez del labio (d)
    // Rho_L      	Factor de reduccion del labio <= 1
    // d_efe		Longitud efectiva del labio ( posiblemente reducido por C2 )
    //				Sobre el labio act˙an esfuerzos de gradiente
    // K			Coef. de rigidez del labio a compresiÛn o del ala.
    //				Se emplea para ambos.
    // esf_w	Esfuerzo de compresiÛn al que est· sometido el ala
    //				Para el c·lculo se considera que el eje neutro
    //				est· a la mitad de H.
    // S			S = sqrt ( E / f ), de donde f es el esfuerzo al que est·
    //				sometido el ala.

	// Para el labio
	K = 0.43;

	Ku = 0.43;

	// C·lculo del esfuerzo al que est· sometido el ala a compresiÛn
	esf_w = Fcx * H / ( 2 * ( D - t ) + H );

	S = 1.28 * sqrt( E / esf_w );

	// Ancho efectivo del labio
	Lam_L = (1.052/sqrt(K)) * (d/t) * sqrt(Fcx/E);
	if( Lam_L <= 0.673 ){
		Rho_L = 1.0;
	}
	else{
		Rho_L = (1.0 - 0.22/Lam_L) / Lam_L;
	}
	if( Rho_L > 1.0 ){
		Rho_L = 1.0;
	}

	d_efe = Rho_L * d;


	// Determinacion del caso a evaluar segun el articulo B4.2 del AISI.
	if( (w/t) <= (S/3.0) ){ // Caso I
		// En este caso el ala y el labio son
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

		// Longitud efectiva del ala
		w_efe = w;

		// Longitud efectiva del labio
		d_efe = d_efe;
	}
	else{
		if( (w/t) < S){ // Caso II
			// Constante de rigidez del ala.
			n  = 1.0 / 2.0;
			// Momento de inercia requerido del labio
			tempo = ((w/t)/S) - (sqrt(Ku)/2.0);
			tempo = tempo * tempo * tempo;
			Ia = 399.00 * tempo * t * t * t * t;
		}
		else{ // Caso III
			// Constante de rigidez del ala.
			n  = 1.0 / 3.0;
			// Momento de inercia requerido del labio.
			Ia = ((115.0*(w/t)/S) + 5.0) * t * t * t * t;
		}

		// Momento de inercia del labio superior.
		Is = 1.0/12.0 * d * d * d * t;

		// Razon de momento de inercia del atiesador
		// lateral del ala presente respecto al
		// requerido que no puede ser mayor que 1.0.
		C2 = Is / Ia;
		if( C2 > 1.0 ){
			C2 = 1.0;
		}

		// Coeficiente de rigidez del labio respecto
		// al ala que no puede ser mayor que 4.0.
		Ka = 5.25 - 5.0 * (D/w);
		if( Ka > 4.0 ){
			Ka = 4.0;
		}

		// Coeficiente de rigidez del ala
		// que no puede ser mayor que 4.0.
		K  = pow(C2,n) * (Ka - Ku) + Ku;
		//if( K > 4.0 ){
		//	K = 4.0;
		//}

		// Factor de esbeltez del ala.
		Lam_F = (1.052/sqrt(K)) * (w/t) * sqrt(esf_w/E);

		// Factor de Reduccion de longitud efectiva
		// del Flange que no puede ser mayor que 1.0.
		if( Lam_F <= 0.673 ){
			Rho_F = 1.0;
		}
		else{
			Rho_F = (1.0 - 0.22/Lam_F) / Lam_F;
		}
		if( Rho_F > 1.0 ){
			Rho_F = 1.0;
		}

		// Longitud efectiva del ala a comrpesiÛn
		w_efe = Rho_F * w;

		// Longitud efectiva del labio reducida a compresiÛn
		d_efe = Rho_L * C2 * d_efe;
	}



	/** 1.3 LocalizaciÛn del eje neutro **/

    /** Variables **/
    // Lam_W     	Factor de esbeltez del peralte (h)
    // Rho_W      	Factor de reduccion del peralte <= 1
    // h_efe		Longitud efectiva del peralte
    //				El peralte est· sometido a esfuerzos que varÌan
    //				linealmente (de gradiente).
    // h_efe_1		Longitud efectiva del peralte en la parte superior
    // h_efe_2		Longitud efectiva del peralte en la parte inferior
    // h_inefe     	Longitud no efectiva del peralte que no trabaja
    // sum_efe_a	Suma de longitudes efectivas del peralte
    // sum_efe_s  	Suma de longitudes efectivas del peralte
    // conta 		Contador
    // AREAS[]		Vector de ·reas de cada parte de la secciÛn
    // DISTA[]		Distancia del centroide de cada parte de la
    //				secciÛn, medida a partir de la fibra superior
    //				de la secciÛn paralela al eje Y
    // CC			Centroide medido desde la fibra superior del labio
    // Area			¡rea total de la secciÛn
    // f1			Esfuerzo de compresiÛn en la parte superior de w.
    // f2 			Esfuerzo de compresiÛn o tensiÛn en la parte infe-
    //				rior de w.
    // psi			RelaciÛn f1/f2
    // Kw			Coeficiente de pandeo, en este caso del peralte


	h_efe_1 = 0.5 * h;
	h_efe_2 = 0.0;
	h_inefe = 0.0;

	// Inicializacion de controladores del ciclo while.
	sum_efe_a = h;
	sum_efe_s = h_efe_1 + h_efe_2;
	conta = 0;

	while ((fabs((sum_efe_s-sum_efe_a)) > 1.0e-6) && (conta < 20) && (h_inefe >= 0.0)){

		// LocalizaciÛn del eje neutro y calculo del momento de inercia
		// Ixx y mÛdulo de secciÛn Sx respecto al eje X, tomando como
		// referencia el eje que pasa sobre la fibra superior.

		// ¡reas de cada subsecciÛn
		AREAS[0]  = d * t;
		AREAS[1]  = ACR;
		AREAS[2]  = w * t;
		AREAS[3]  = AREAS[1];
		AREAS[4]  = ( h - h_efe_1 - h_inefe ) * t;
		AREAS[5]  = h_efe_1 * t;
		AREAS[6]  = AREAS[1];
		AREAS[7]  = w_efe * t;
		AREAS[8]  = AREAS[1];
		AREAS[9]  = d_efe * t;

		// Distancias a los centroides
		DISTA[0]  = d + h + 4.0 * R + 2.0 * t + 0.5 * d;
		DISTA[1]  = d + h + 4.0 * R + 2.0 * t - CCR;
		DISTA[2]  = d + h + 3.0 * R + 1.5 * t;
		DISTA[3]  = d + h + 2.0 * R + t + CCR;
		DISTA[4]  = d + h + 2.0 * R + t - 0.5 * ( h - h_efe_1- h_inefe );
		DISTA[5]  = d + 2.0 * R + t + 0.5 * h_efe_1;
		DISTA[6]  = d + 2.0 * R + t - CCR;
		DISTA[7]  = d + R + 0.5 * t;
		DISTA[8]  = d + CCR;
		DISTA[9]  = d - 0.5 * d_efe;

		Mome = 0.0;
		for ( i = 0 ; i < 10 ; i++ ){
			Mome = Mome + AREAS[i] * DISTA[i];
		}

		Area = 0.0;
		for ( i = 0 ; i < 10 ; i++ ){
			Area = Area + AREAS[i];
		}

		// Centroide
		CC = Mome / Area;

		// Esfuerzos en la parte plana del peralte (Web).
		f1 = + Fcx * ( ( CC - 2 * R - t - d ) / CC );
		f2 = - Fcx * ( ( h - ( CC - 2 * R - t - d ) ) / CC );

		Psi = f2 / f1;

		tempo = 1.0 - Psi;

		// Coef. de pandeo de la web
		Kw = 4.0 + 2.0 * tempo * tempo * tempo + 2.0 * tempo;

		// Factor de Esbeltez del Web.
		if ( ( h / t ) >= 200 ){
			printf ( "\nError, cociente h/t del elemento pasa el lÌmite de 200" );
		}

		Lam_W = ( 1.052 / sqrt(Kw) ) * ( h/t ) * sqrt( f1/E );

		if( Lam_W <= 0.673 ){
			Rho_W = 1.0;
		}
		else{
			Rho_W = ( 1.0 - 0.22/Lam_W ) / Lam_W;
		}
		if( Rho_W > 1.0 ){
			Rho_W = 1.0;
		}

		// Longitud efectiva del peralte.
		h_efe = Rho_W * h;

		// Longitud efectiva del peralte en la parte superior.
		h_efe_1 = h_efe / (3.0-Psi);

		// Longitud efectiva del peralte en la parte inferior.
		if( Psi <= -0.236 ){
			h_efe_2 = h_efe / 2.0;
		}
		else{
			h_efe_2 = h_efe - h_efe_1;
		}

		// Longitud no efectiva del peralte.
		h_inefe = ( CC - 2 * R - t - d ) - ( h_efe_1 + h_efe_2 );

		// Actualizacion de la suma de longitudes efectivas en el peralte.
		sum_efe_a = sum_efe_s;
		sum_efe_s = h_efe_1 + h_efe_2;

		conta += 1;
	}

	// CorrecciÛn de la longitud no efectiva del peralte
	if( h_inefe < 0.0 ){
		h_inefe = 0.0;
		h_efe_1 = 0.5 * h;
		h_efe_2 = 0.0;
	}

	/** 1.4 Momento de inercia y mÛdulo de seccion **/

    /** Variables **/
    // MILOC[]		Momentos de inercia en X respecto a ejes locales de
    //				cada sub-secciÛn
    // Ixx			Momento de inercia de la secciÛn respecto al eje X
    // Cy 			Distancia m·xima de las fibra m·s lejada al eje
    //				centroidal
    // Sx			MÛdulo de secciÛn con con respecto al eje X
    // Mnx			Momento nominal a la flexiÛn sobre el eje X

	// Areas de las subsecciones
	AREAS[0]  = d * t;
	AREAS[1]  = ACR;
	AREAS[2]  = w * t;
	AREAS[3]  = AREAS[1];
	AREAS[4]  = ( h - h_efe_1 - h_inefe ) * t;
	AREAS[5]  = h_efe_1 * t;
	AREAS[6]  = AREAS[1];
	AREAS[7]  = w_efe * t;
	AREAS[8]  = AREAS[1];
	AREAS[9]  = d_efe * t;

	// Distancias a los centroides
	DISTA[0]  = d + h + 4.0 * R + 2.0 * t + 0.5 * d_efe;
	DISTA[1]  = d + h + 4.0 * R + 2.0 * t - CCR;
	DISTA[2]  = d + h + 3.0 * R + 1.5 * t;
	DISTA[3]  = d + h + 2.0 * R + t + CCR;
	DISTA[4]  = d + h + 2.0 * R + t - 0.5 * ( h - h_efe_1- h_inefe );
	DISTA[5]  = d + 2.0 * R + t + 0.5 * h_efe_1;
	DISTA[6]  = d + 2.0 * R + t - CCR;
	DISTA[7]  = d + R + 0.5 * t;
	DISTA[8]  = d + CCR;
	DISTA[9]  = d - 0.5 * d_efe;

	for ( i = 0 ; i < 10 ; i++ ){
		Mome = Mome + AREAS[i] * DISTA[i];
	}

	for ( i = 0 ; i < 10 ; i++ ){
		Area = Area + AREAS[i];
	}

	// Centroide
	CC = Mome / Area;

	// Momentos de Inercia en X
	MILOC[0] = DOVO * d * d * d * t;
	MILOC[1] = ICR;
	MILOC[2] = DOVO * t * t * t * w;
	MILOC[3] = ICR;
	MILOC[4] = DOVO * ( h - h_efe_1 - h_inefe ) * ( h - h_efe_1 - h_inefe )
    * ( h - h_efe_1 - h_inefe ) * t;
	MILOC[5] = DOVO * h_efe_1 * h_efe_1 * h_efe_1 * t;
	MILOC[6] = ICR;
	MILOC[7] = DOVO * t * t * t * w_efe;
	MILOC[8] = ICR;
	MILOC[9] = DOVO * d_efe * d_efe * d_efe * t;

	// Distancia paralela al eje Y del centroide de la seccion
    // a los ejesinerciales locales de cada sub-seccion.
	DISTA[0]  = d + h + 4.0 * R + 2.0 * t + 0.5 * d_efe - CC;
	DISTA[1]  = d + h + 4.0 * R + 2.0 * t - CC;
	DISTA[2]  = d + h + 3.0 * R + 1.5 * t - CC;
	DISTA[3]  = d + h + 2.0 * R + t - CC;
	DISTA[4]  = d + h + 2.0 * R + t - 0.5 * ( h - h_efe_1- h_inefe ) - CC;
	DISTA[5]  = d + 2.0 * R + t + 0.5 * h_efe_1 - CC;
	DISTA[6]  = d + 2.0 * R + t - CC;
	DISTA[7]  = d + R + 0.5 * t - CC;
	DISTA[8]  = d - CC;
	DISTA[9]  = d - 0.5 * d_efe - CC;


	// Momento de inercia

	Ixx = 0.0;
	for ( i = 0 ; i < 10 ; i++ ){
		Ixx = Ixx + MILOC[i] + DISTA[i] * DISTA[i] * AREAS[i];
	}

	// Distancias maximas de las fibras
	// extremas a los ejes centroidales.
	Cy = CC;
	if( Cy < H + 2 * D - 2 * t - Cy ){
		Cy = H + 2 * D - 2 * t - Cy;
	}

	// Modulo de Seccion calculado con la maxima
	// distancia del centroide a fibras extremas.
	Sx = Ixx / Cy;

	/**		2 Momento nominal y permitido		**/

	*Mnx = Sx * Fcx;

    //jacob.provisional impresion
	//printf ( "\n Mnx = %lf", *Mnx );

    //[][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][]

    //[][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][]

	/**********************************
     FLEXON SOBRE EJE Y
     **********************************/

	/** NOTA **/
	// En este caso hay dos posibles casos a analizar y dado que la seccion
	// no es simetrica, sobre el eje Y, la direccion del momento si importa.
	/* Caso 1					Caso 2

        My
        ^
        ^
        |	+					 	|	+
     + +|+ +					 + +|+ +
     +	|						+	|
     +	|						+ 	|
  ---+--|-----				  --+---|-----
     +	|						+ 	|
     +	|						+	|
     + +|+ + 					 + +|+ +
        |	+						|	+
        v
        v
        My

     Momento pos. causa ten-		Momento neg. causa compre-
     siÛn en el lado izq.		siÛn en el lado izq. de la
     de la secciÛn				secciÛn.

     }*/

    //[][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][]

    //[][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][]

	/** Caso 1 **/

	/**		1 Propiedades de la secciÛn		**/

	/** 1.1 Cuarto de rondana **/
	// Ya fueron calculadas anteriormente

	/** 1.2 Longitud efeciva del labio que est· a compresiÛn **/
    /** Variables **/
    // Lam_ L		Factor de esbeltez del labio
    // Rho_L		Factor de reducciÛn efectiva del labio
    // d_efe		Longitd efectiva de cada labio
    // w_efe_1		Longitud efectiva del ala en la parte superior
    // w_efe_2		Longitud efectiva del ala en la parte inferior
    //				En este caso el ala se halla sobre esfuerzos de
    //				gradiente.
    // w_inefe		Longitud no efectiva del ala que no trabaja


	// Factor de Esbeltez del Labio.
	Ku = 0.43;
	Lam_L = (1.052/sqrt(Ku)) * (d/t) * sqrt(Fcy1/E);

	// Factor de reduccion de longitud efectiva
	// del labio que no puede ser mayor que 1.0
	if( Lam_L <= 0.673 ){
		Rho_L = 1.0;
	}
	else{
		Rho_L = (1.0 - 0.22/Lam_L) / Lam_L;
	}
	if( Rho_L > 1.0 ){
		Rho_L = 1.0;
	}

	// Longitud efectiva de cada labio a compresiÛn
	d_efe = Rho_L * d;

	// Evaluacion de las longitudes efectivas de las alas consider·ndolos
	// como elementos atiesados sujeto a gradiente de esfuerzo segun la
	// metodologÌa del AISI B2.3 (pagina 113 de Yu).

	/** 1.3 LocalizaciÛn del eje neutro **/

    /** Variables **/



	// InicializaciÛn de variables para el proceso iterativo.
	w_efe_1 = 0.5 * w;
	w_efe_2 = 0.0;
	w_inefe = 0.0;

	// InicializaciÛn de controladores del ciclo while.
	sum_efe_a = w;
	sum_efe_s = w_efe_1 + w_efe_2;
	conta = 0;

	// Inicio de iteraciones para encontrar la longitud efectiva
	// del ala, esta restringido a 20 iteraciones.
	while((fabs((sum_efe_s-sum_efe_a)) > 1.e-6) && (conta < 20) && (w_inefe >= 0.0)){

		// LocalizaciÛn del eje neutro y c·lculo del momento de inercia
		// Iy1 y mÛdulo de secciÛn Sy1 respecto al eje Y, tomando como
		// referencia el eje que pasa sobre la fibra superior.

		// Areas de cada sub-secciÛn
		AREAS[0]  = d_efe * t;
		AREAS[1]  = ACR;
		AREAS[2]  = w_efe_1 * t;
		AREAS[3]  = ( w - w_efe_1 - w_inefe ) * t;
		AREAS[4]  = AREAS[1];
		AREAS[5]  = h * t;
		AREAS[6]  = AREAS[4];
		AREAS[7]  = AREAS[3];
		AREAS[8]  = AREAS[2];
		AREAS[9]  = AREAS[1];
		AREAS[10] = AREAS[0];

		// Distancias de centroides de cada
		// sub-secciÛn a la fibra superior
		DISTA[0]  = 0.5 * t;
		DISTA[1]  = t + R - CCR;
		DISTA[2]  = t + R + 0.5 * w_efe_1;
		DISTA[3]  = t + R + w - 0.5 * ( w - w_efe_1 - w_inefe );
		DISTA[4]  = t + R + w + CCR;
		DISTA[5]  = 1.5 * t + 2.0 * R;
		DISTA[6]  = DISTA[4];
		DISTA[7]  = DISTA[3];
		DISTA[8]  = DISTA[2];
		DISTA[9]  = DISTA[1];
		DISTA[10] = DISTA[0];

		// Distancia de fibra superior al eje neutro
		// Momento de area total de la secciÛn
		Mome = 0.0;
		for ( i = 0 ; i < 11 ; i++ ){
			Mome = Mome + AREAS[i] * DISTA[i];
		}

		// ¡rea total de la seccion
		Area = 0.0;
		for ( i = 0 ; i < 11 ; i++ ){
			Area = Area + AREAS[i];
		}

		// Distancia de la fibra superior al eje neutro.
		CC = Mome / Area;

		// Esfuerzos en la parte plana de las alas.
		f1 = + Fcy1 * ( ( CC - R - t ) / CC );     // CompresiÛn en fibra superior
		f2 = - Fcy1 * ( ( w - ( CC - R - t ) ) / CC ); // TensiÛn en fibra inferior

		// Razon de esfuerzos presentes en el ala
		Psi = f2 / f1;

		// Ceoficiente de rigidez del ala
		tempo = 1.0 - Psi;
		Kf = 4.0 + 2.0 * tempo * tempo * tempo + 2.0 * tempo;

		// Factor de esbeltez del ala
		Lam_F = ( 1.052 / sqrt( Kf ) ) * ( w / t ) * sqrt ( f1 / E );

		// Factor de reduccion de longitud efectiva
		// del ala que no puede ser mayor que 1.0.
		if( Lam_F <= 0.673 ){
			Rho_F = 1.0;
		}
		else{
			Rho_F = ( 1.0 - 0.22 / Lam_F ) / Lam_F;
		}
		if( Rho_F > 1.0 ){
			Rho_F = 1.0;
		}

		// Longitud efectiva del ala
		w_efe = Rho_F * w;

		// Longitud efectiva del ala en la parte superior
		w_efe_1 = w_efe / ( 3.0 - Psi );

		// Longitud efectiva del ala en la parte inferior.
		if( Psi <= -0.236 ){
			w_efe_2 = w_efe / 2.0;
		}
		else{
			w_efe_2 = w_efe - w_efe_1;
		}

		// Longitud no efectiva del ala
		w_inefe = ( CC - R - t ) - ( w_efe_1 + w_efe_2 );

		// Actualizacion de la suma de
		// longitudes efectivas del ala.
		sum_efe_a = sum_efe_s;
		sum_efe_s = w_efe_1 + w_efe_2;

		// Actualizacion del contador.
		conta = conta + 1;
	}


	/** 1.4 Momento de inercia y modulo de seccion **/

    /** Variables **/



	// CorrecciÛn de la longitud no efectiva del
	// ala o Ancho que no puede ser negativa.
	if( w_inefe < 0.0 ){
		w_inefe = 0.0;
		w_efe_1 = 0.5 * w;
		w_efe_2 = 0.0;
	}

	// Areas de cada sub-secciÛn
	AREAS[0]  = d_efe * t;
	AREAS[1]  = ACR;
	AREAS[2]  = w_efe_1 * t;
	AREAS[3]  = ( w - w_efe_1 - w_inefe ) * t;
	AREAS[4]  = AREAS[1];
	AREAS[5]  = h * t;
	AREAS[6]  = AREAS[4];
	AREAS[7]  = AREAS[3];
	AREAS[8]  = AREAS[2];
	AREAS[9]  = AREAS[1];
	AREAS[10] = AREAS[0];

	// Distancias de centroides de cada
	// sub-secciÛn a la fibra superior
	DISTA[0]  = 0.5 * t;
	DISTA[1]  = t + R - CCR;
	DISTA[2]  = t + R + 0.5 * w_efe_1;
	DISTA[3]  = t + R + w - 0.5 * ( w - w_efe_1 - w_inefe );
	DISTA[4]  = t + R + w + CCR;
	DISTA[5]  = 1.5 * t + 2.0 * R;
	DISTA[6]  = DISTA[4];
	DISTA[7]  = DISTA[3];
	DISTA[8]  = DISTA[2];
	DISTA[9]  = DISTA[1];
	DISTA[10] = DISTA[0];

	// Distancia de fibra superior al eje neutro.
	// Momento de ·rea total de la secciÛn.
	Mome = 0.0;
	for ( i = 0 ; i < 11 ; i++ ){
		Mome = Mome + AREAS[i] * DISTA[i];
	}

	// Area total de la seccion.
	Area = 0.0;
	for ( i = 0 ; i < 11 ; i++ ){
		Area = Area + AREAS[i];
	}

	// Distancia de la fibra superior al eje neutro.
	CC = Mome / Area;

	// Momentos de inercia en Y respecto a
	// ejes locales de cada sub-secciÛn.
	MILOC[0]  = DOVO * t * t * t * d_efe;
	MILOC[1]  = ICR;
	MILOC[2]  = DOVO * w_efe_1 * w_efe_1 * w_efe_1 * t;
	MILOC[3]  = DOVO * (w - w_efe_1 - w_inefe) * (w - w_efe_1 - w_inefe)
    * (w - w_efe_1 - w_inefe) * t;
	MILOC[4]  = ICR;
	MILOC[5]  = DOVO * t * t * t * h;
	MILOC[6]  = MILOC[4];
	MILOC[7]  = MILOC[3];
	MILOC[8]  = MILOC[2];
	MILOC[9]  = MILOC[1];
	MILOC[10] = MILOC[0];

	// Distancia paralela al eje Y del
	// centroide de la secciÛn a los ejes
	// inerciales locales de cada sub-secciÛn.
	DISTA[0]  = 0.5 * t - CC;
	DISTA[1]  = t + R - CC;
	DISTA[2]  = t + R + 0.5 * w_efe_1 - CC;
	DISTA[3]  = t + R + w - 0.5 * ( w - w_efe_1 - w_inefe ) - CC;
	DISTA[4]  = t + R + w - CC;
	DISTA[5]  = 1.5 * t + 2.0 * R - CC;
	DISTA[6]  = DISTA[4];
	DISTA[7]  = DISTA[3];
	DISTA[8]  = DISTA[2];
	DISTA[9]  = DISTA[1];
	DISTA[10] = DISTA[0];

	// Momento de inercia de la secciÛn respecto al eje Y
	Iy1 = 0.0;
	for ( i = 0 ; i < 11 ; i++ ){
		Iy1 = Iy1 + MILOC[i] + DISTA[i] * DISTA[i] * AREAS[i];
	}

	// Distancias m·ximas de las fibras
	// extremas a los ejes centroidales.
    Cy = CC;
	if( Cy < fabs ( fabs ( W ) - fabs ( Cy ) ) ){
		Cy = fabs ( fabs ( W ) - fabs ( Cy ) );
	}

	// Modulo de seccion calculado con la maxima
	// distancia del centroide a fibras extremas.
	Sy1 = Iy1 / Cy;

	/**		2 Momento nominal y permitido		**/

	*Mny1 = Sy1 * Fcy1;

    //jess.privional impresion
	//printf ( "\n Mny1 = %lf", *Mny1 );

    //[][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][]

    //[][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][]

	/** Caso 2 **/

	/**		1 Propiedades de la seccion		**/

	/** 1.1 Cuarto de rondana **/
	// Ya fueron calculadas anteriormente

	/** 1.2 Longitud efectiva del peralte que est· a compresion y **/
	/** 1.3 Localizacion del eje neutro **/

    /** Variables **/
    // Lam_H		Factor de esbeltez del peralte a compresiÛn (h)
    // Rho_H		Factor de reduccion de longitud efectiva del
    //				peralte a compresiÛn <= 1 (h)
    // h_efe		Longitud efectiva del peralte
    //				En este caso la compresiÛn sobre el peralte
    //				es constante.
    // Lam_F		Factor de esbeltez del ala a compresiÛn (w)
    // Rho_F		Factor de reduccion del ala a compresiÛn <= 1
    // w_efe		Longitud efectiva del ala
    //				En este el ala est· sometida a esfuerzos de
    //				gradiente.
    // K			Coeficiente de rigidez del lado a peralte (h)
    // w_efe_1		Longitud efectiva del ala en la parte superior.
    // w_efe_2		Longitud efectiva del ala en la parte inferior.
    // w_inefe		Longitud no efectiva del ala que no trabaja.


	// a.- Se calcula normalmente la longitud efectiva del peralte a compresiÛn
	//	   considerando un esfuerzo igual a Fy.
	Lam_H = 0.526 * ( h / t ) * sqrt ( Fcy2 / E );

	if( Lam_H <= 0.673 ){
		Rho_H = 1.0;
	}
	else{
		Rho_H = ( 1.0 - 0.22 / Lam_H ) / Lam_H;
	}
	if( Rho_H > 1.0 ){
		Rho_H = 1.0;
	}

	h_efe = Rho_H * h;

	// b.- Se itera para encontrar el eje neutro para que el ala sea total-
	// mente efectiva

	w_efe_1 = 0.5 * h;
	w_efe_2 = 0.0;
	w_inefe = 0.0;

	sum_efe_a = w;
	sum_efe_s = w_efe_1 + w_efe_2;
	conta = 0;

	while((fabs((sum_efe_s-sum_efe_a)) > 1.0e-6) && (conta < 20) && (w_inefe >= 0.0)){

		// LocalizaciÛn del eje neutro y c·lculo del momento de inercia
		// Iy2 y mÛdulo de seccion Sy2 respecto al eje Y, tomando como
		// referencia el eje que pasa sobre la fibra superior.

		// Areas de cada sub-seccion./

		AREAS[0]  = d * t;
		AREAS[1]  = ACR;
		AREAS[2]  = ( w - w_efe_1 - w_inefe ) * t;
		AREAS[3]  = w_efe_1 * t;
		AREAS[4]  = AREAS[1];
		AREAS[5]  = h_efe * t;
		AREAS[6]  = AREAS[1];
		AREAS[7]  = AREAS[3];
		AREAS[8]  = AREAS[2];
		AREAS[9]  = AREAS[1];
		AREAS[10] = AREAS[0];

		// Distancias de centroides de cada
		// sub-seccion a la fibra superior.

		DISTA[0]  = 1.5 * t + 2.0 * R + w;
		DISTA[1]  = t + R + w + CCR;
		DISTA[2]  = t + R + w - 0.5 * ( w - w_efe_1 - w_inefe );
		DISTA[3]  = t + R + 0.5 * w_efe_1;
		DISTA[4]  = t + R - CCR;
		DISTA[5]  = 0.5 * t;
		DISTA[6]  = DISTA[4];
		DISTA[7]  = DISTA[3];
		DISTA[8]  = DISTA[2];
		DISTA[9]  = DISTA[1];
		DISTA[10] = DISTA[0];

		// Distancia de fibra superior al eje neutro.
		// Momento de area total de la seccion.
		Mome = 0.0;
		for ( i = 0 ; i < 11 ; i++ ){
			Mome = Mome + AREAS[i] * DISTA[i];
		}

		// Area total de la seccion.
		Area = 0.0;
		for ( i = 0 ; i < 11 ; i++ ){
			Area = Area + AREAS[i];
		}

		// Distancia de la fibra superior al eje neutro.
		CC = Mome / Area;

		// Esfuerzos en la parte plana de las alas.
		f2 = - Fcy2 * ( w + R + t - CC ) / ( w + 2 * ( R + t ) - CC );
		f1 = + Fcy2 * ( CC - R - t ) / ( w + 2 * ( R + t ) - CC );

		// Razon de esfuerzos presentes en el ala.
		Psi = f2 / f1;

		// Coeficiente de rigidez del ala.
		tempo = 1.0 - Psi;
		Kf = 4.0 + 2.0*tempo*tempo*tempo + 2.0*tempo;

		// Factor de esbeltez del ala.
		Lam_F = (1.052/sqrt(Kf)) * (w/t) * sqrt(f1/E);

		// Factor de reduccion de longitud efectiva
		// del ala que no puede ser mayor que 1.0.
		if( Lam_F <= 0.673 ){
			Rho_F = 1.0;
		}
		else{
			Rho_F = (1.0 - 0.22/Lam_F) / Lam_F;
		}
		if( Rho_F > 1.0 ){
			Rho_F = 1.0;
		}

		// Longitud efectiva del ala.
		w_efe = Rho_F * w;

		// Longitud efectiva del ala en la parte superior.
		w_efe_1 = w_efe / (3.0-Psi);

		// Longitud efectiva del Flange en la parte inferior.
		if( Psi <= -0.236 ){
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

	// c. Determina si el eje neutro est· m·s cerca del lado de compresiÛn
	//	  y en caso de ser necesario manda a recalcular la longitud
	//	  efectiva del peralte.

	if ( CC < 0.5 * w + R + t ){
		esf_h = Fcy2 * ( CC / ( 2.0 * R + 2.0 * t + w ) );
		if ( fabs ( Fcy2 - esf_h ) > 1.0e-2){
			h_efe = recalc_h ( 1, w, h, d, R, t, Fcy2, ACR, CCR, E, &a, &b );
		}
	}

	// Valores calculados en recalc_h
	w_efe_1 = a;
	w_inefe = b;

	// Correccion de la longitud no efectiva del
	// ala o ancho que no puede se negativa.
	if( w_inefe < 0.0 ){
		w_inefe = 0.0;
		w_efe_1 = 0.5 * w;
		w_efe_2 = 0.0;
	}

	// Areas de cada sub-seccion./

	AREAS[0]  = d * t;
	AREAS[1]  = ACR;
	AREAS[2]  = ( w - w_efe_1 - w_inefe ) * t;
	AREAS[3]  = w_efe_1 * t;
	AREAS[4]  = AREAS[1];
	AREAS[5]  = h_efe * t;
	AREAS[6]  = AREAS[1];
	AREAS[7]  = AREAS[3];
	AREAS[8]  = AREAS[2];
	AREAS[9]  = AREAS[1];
	AREAS[10] = AREAS[0];

	// Distancias de centroides de cada
	// sub-seccion a la fibra superior.

	DISTA[0]  = 1.5 * t + 2.0 * R + w;
	DISTA[1]  = t + R + w + CCR;
	DISTA[2]  = t + R + w - 0.5 * ( w - w_efe_1 - w_inefe );
	DISTA[3]  = t + R + 0.5 * w_efe_1;
	DISTA[4]  = t + R - CCR;
	DISTA[5]  = 0.5 * t;
	DISTA[6]  = DISTA[4];
	DISTA[7]  = DISTA[3];
	DISTA[8]  = DISTA[2];
	DISTA[9]  = DISTA[1];
	DISTA[10] = DISTA[0];

	// Distancia de fibra superior al eje neutro.
	// Momento de area total de la seccion.
	Mome = 0.0;
	for ( i = 0 ; i < 11 ; i++ ){
		Mome = Mome + AREAS[i] * DISTA[i];
	}

	// Area total de la seccion.
	Area = 0.0;
	for ( i = 0 ; i < 11 ; i++ ){
		Area = Area + AREAS[i];
	}

	// Distancia de la fibra superior al eje neutro.
	CC = Mome / Area;


	/** 1.4 Momento de inercia y mÛdulo de seccion **/

    /** Variables **/
    // MILOC[]		Momentos de inercia en Y respecto a ejes locales de
    //				cada sub-secciÛn
    // Iy2			Momento de inercia de la secciÛn respecto al eje Y
    // Cy 			Distancia m·xima de las fibra m·s lejada al eje
    //				centroidal
    // Sy2			MÛdulo de secciÛn con con respecto al eje Y
    // Mny2			Momento nominal a la flexiÛn sobre el eje Y

	// Momentos de Inercia en Y respecto a
	// ejes locales de cada sub-seccion.

	MILOC[0]  = DOVO * t * t * t * d;
	MILOC[1]  = ICR;
	MILOC[2]  = DOVO * (w - w_efe_1 - w_inefe) * (w - w_efe_1 - w_inefe)
    * (w - w_efe_1 - w_inefe) * t;
	MILOC[3]  = DOVO * w_efe_1 * w_efe_1 * w_efe_1 * t;
	MILOC[4]  = ICR;
	MILOC[5]  = DOVO * t * t * t * h_efe;
	MILOC[6]  = MILOC[4];
	MILOC[7]  = MILOC[3];
	MILOC[8]  = MILOC[2];
	MILOC[9]  = MILOC[1];
	MILOC[10] = MILOC[0];

	// Distancia paralela al eje Y del
	// centroide de la seccion a los ejes
	// inerciales locales de cada sub-seccion.

	DISTA[0]  = 1.5 * t + 2.0 * R + w - CC;
	DISTA[1]  = t + R + w - CC;
	DISTA[2]  = t + R + w - 0.5 * ( w - w_efe_1 - w_inefe ) - CC;
	DISTA[3]  = t + R + 0.5 * w_efe_1 - CC;
	DISTA[4]  = t + R - CC;
	DISTA[5]  = 0.5 * t - CC;
	DISTA[6]  = DISTA[4];
	DISTA[7]  = DISTA[3];
	DISTA[8]  = DISTA[2];
	DISTA[9]  = DISTA[1];
	DISTA[10] = DISTA[0];

	// Momento de inercia de la seccion respecto al eje Y.
	Iy2 = 0.0;
	for ( i = 0 ; i < 11 ; i++ ){
		Iy2 = Iy2 + MILOC[i] + DISTA[i] * DISTA[i] * AREAS[i];
	}

	// Distancias maximas de las fibras
	// extremas a los ejes centroidales.
	Cy = CC;
	if( Cy < fabs(fabs(W)-fabs(Cy)) ){
		Cy = fabs(fabs(W)-fabs(Cy));
	}

	// Modulo de Seccion calculado con la maxima
	// distancia del centroide a fibras extremas.
	Sy2 = Iy2 / Cy;

	/**		2 Momento nominal y permitido		**/

	*Mny2 = Sy2 * Fcy2;

    //jess.provisional impresionx
	//printf ( "\n Cy = %lf", Cy );
	//printf ( "\n Iy2 = %lf", Iy2 );
	//printf ( "\n Sy2 = %lf", Sy2 );
	//printf ( "\n Mny2 = %lf", *Mny2 );
}

//[][][][][][][][][][][][][][][][][][][][[][][][][][][][][][][][][][][][][][][][]
// PREGUNTAR ACERCA DEL IDENTIFICADOR DE CADA TIPO DE SECCION
// CHECAR LA CONVERGENCIA, EL PASO,
double recalc_h ( int tipo, double w, double h, double d, double R,
                 double t, double Fy, double ACR, double CCR,
                 double E, double *a, double *b )
{
	/** DescrpciÛn **/
	// Recalcula la longitud efectiva del elemento a compresiÛn cuando
	// el eje neutro se encuentra sobre esta zona y no sobre la de tensiÛn.
	// Ahora sola est· implementado lo correspondiente a Omegas.

	/** Variables de entrada **/
	// tipo		Tipo de elemento que se emplea:
	//			1.- Omegas
	//			2.- C con labios
	//			3.- C sin labios
	// w		Longitud interna del ala de la secciÛn
	// h		Longitud interna del peralte
	// d		Longitud internna del labio
	// R		Radio de doblez
	// t		Espesor de la l·mina empleada para el R.F.
	// a		Variable auxiliar para guardar w_efe_1
	// b		VAriable auxiliar para guardar w_inefe



	/** Variables de salida **/
	// h_efe		Longitud efectiva del peralte

	/** Variables locales **/
	// paso		Paso con el que se aumenta o reduce el esfuerzo supesto
	//			sobre el peralte

	int conta, i;
	double paso, Lam_H, Rho_H, h_efe, w_efe_1, w_efe_2, w_inefe;
	double sum_efe_a, sum_efe_s, Mome, Area, CC, f1, f2, Psi;
	double tempo, Kf, Lam_F, Rho_F, AREAS [11], DISTA[11];
	double w_efe, f_real, f_sup;

	f_sup = Fy;

	paso = 1.0;

	do{

		Lam_H = 0.526 * ( h / t ) * sqrt ( f_sup / E );
		if( Lam_H <= 0.673 ){
			Rho_H = 1.0;
        }
		else{
			Rho_H = ( 1.0 - 0.22 / Lam_H ) / Lam_H;
		}
		if( Rho_H > 1.0 ){
			Rho_H = 1.0;
		}

		h_efe = Rho_H * h;

		w_efe_1 = 0.5 * h;
		w_efe_2 = 0.0;
		w_inefe = 0.0;

		sum_efe_a = w;
		sum_efe_s = w_efe_1 + w_efe_2;
		conta = 0;

		while((fabs((sum_efe_s-sum_efe_a)) > 1.0e-6) && (conta < 20) && (w_inefe >= 0.0)){

			// Localizacion del eje neutro y calculo del momento de inecia
			// Iyy y Modulo de seccion Sy respecto al eje Y tomando como
			// referencia el eje que pasa sobre la fibra superior.

			// Areas de cada sub-seccion./
			AREAS[0]  = d * t;
			AREAS[1]  = ACR;
			AREAS[2]  = ( w - w_efe_1 - w_inefe ) * t;
			AREAS[3]  = w_efe_1 * t;
			AREAS[4]  = AREAS[1];
			AREAS[5]  = h_efe * t;
			AREAS[6]  = AREAS[1];
			AREAS[7]  = AREAS[3];
			AREAS[8]  = AREAS[2];
			AREAS[9]  = AREAS[1];
			AREAS[10] = AREAS[0];

			// Distancias de centroides de cada
			// sub-seccion a la fibra superior.

			DISTA[0]  = 1.5 * t + 2.0 * R + w;
			DISTA[1]  = t + R + w + CCR;
			DISTA[2]  = t + R + w - 0.5 * ( w - w_efe_1 - w_inefe );
			DISTA[3]  = t + R + 0.5 * w_efe_1;
			DISTA[4]  = t + R - CCR;
			DISTA[5]  = 0.5 * t;
			DISTA[6]  = DISTA[4];
			DISTA[7]  = DISTA[3];
			DISTA[8]  = DISTA[2];
			DISTA[9]  = DISTA[1];
			DISTA[10] = DISTA[0];

			// Distancia de fibra superior al eje neutro.
			// Momento de area total de la seccion.
			Mome = 0.0;
			for ( i = 0 ; i < 11 ; i++ ){
				Mome = Mome + AREAS[i] * DISTA[i];
			}

			// Area total de la seccion.
			Area = 0.0;
			for ( i = 0 ; i < 11 ; i++ ){
				Area = Area + AREAS[i];
			}

			// Distancia de la fibra superior al eje neutro.
			CC = Mome / Area;

			// Esfuerzos en la parte plana de las alas.
			// CHECAR NUEVAMENTE F1 Y F2 PARA ESTA FUNCION recalc_h
			f2 = - Fy * ( w + R + t - CC ) / ( w + 2 * ( R + t ) - CC );
			f1 = + Fy * ( CC - R - t ) / ( w + 2 * ( R + t ) - CC );

			// Razon de esfuerzos presentes en el ala.
			Psi = f2 / f1;

			// Coeficiente de rigidez del ala.
			tempo = 1.0 - Psi;
			Kf = 4.0 + 2.0*tempo*tempo*tempo + 2.0*tempo;

			// Factor de esbeltez del ala.
			Lam_F = (1.052/sqrt(Kf)) * (w/t) * sqrt(f1/E);

			// Factor de reduccion de longitud efectiva
			// del ala que no puede ser mayor que 1.0.
			if( Lam_F <= 0.673 ){
				Rho_F = 1.0;
			}
			else{
				Rho_F = (1.0 - 0.22/Lam_F) / Lam_F;
			}
			if( Rho_F > 1.0 ){
				Rho_F = 1.0;
			}

			// Longitud efectiva del ala.
			w_efe = Rho_F * w;

			// Longitud efectiva del ala en la parte superior.
			w_efe_1 = w_efe / (3.0-Psi);

			// Longitud efectiva del Flange en la parte inferior.
			if( Psi <= -0.236 ){
				w_efe_2 = w_efe / 2.0;
			}
			else{
				w_efe_2 = w_efe - w_efe_1;
			}
			// Longitud no efectiva del Flange.
			w_inefe = (CC-R-t) - (w_efe_1+w_efe_2);//

			// Actualizacion de la suma de
			// longitudes efectivas en Flange.
			sum_efe_a = sum_efe_s;
			sum_efe_s = w_efe_1+w_efe_2;//

			// Actualizacion del contador.
			conta = conta + 1;
		}

		// Calcula el esfuerzo en el peralte en base a la
		// posiciÛn del eje neutro

		f_real = Fy * CC / ( w + 2 * ( R + t ) - CC );
		//printf ( "\nf_sup = %lf\tf_real = %lf\tCC = %lf\tw_inefe = %lf", f_sup, f_real, CC, w_inefe );
		// Modifica el paso
		if ( f_sup - f_real < paso ){
			paso /=10;
		}
		f_sup -= paso;

	} while ( f_sup - f_real > .001 );

	*a = w_efe_1;
	*b = w_inefe;
	return h_efe;
}



