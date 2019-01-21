#ifndef _efi_RF_o_h_
#define _efi_RF_o_h_

double f2eso(double H, double W, double D, double R, double t, double Fy,
             double E, double LongX, double LongY, double LongT, int tb,
             int ifcz, int ifcy, int ifdz, int ifdy, double mfcz, double mfcy,
             double dfcz, double dfcy, double mfdz, double mfdy, double fxi,
             double fyi, double fzi, double fxf,double fyf, double fzf,
             double mxi, double myi, double mzi, double mxf, double myf,
             double mzf, int ielem);

void mnfso( double H, double W, double D, double R, double t, double E, double Fcx,
            double Fcy1, double Fcy2,  double *Mnx, double *Mny1, double *Mny2 );

double recalc_h ( int tipo, double w, double h, double d, double R,
                 double t, double Fy, double ACR, double CCR,
                 double E, double *a, double *b );



//----------------------------------------------------------------------------//
// Este programa calcula el area efectiva de una seccion Omega sometida a
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

#endif
