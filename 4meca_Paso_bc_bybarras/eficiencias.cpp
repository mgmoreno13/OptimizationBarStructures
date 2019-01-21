//*****************************
// Rutina para Eficiencias
//*****************************
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "eficiencias.h"
#include "raros.h"

//************************************************
void Eficiencias(int nelem, int ndime, int ntipo, int* matnu, int* ntips,
				 double** props, double* fuerc, double* valor, double* xlong,
				 double ***arear, double ***arerf, double ***arerc,
				 double ***arecr, double ****dtcon, int** indfu,
				 double** fuerb, double** fuefl, double* wcarg, int iwrit,
				 int *ubiCat, double *efi_max, int nevab )
{
	int ielem,imats,iefi=1,rolado;
	double valef;
	int iCat;


///////////////////////////////////////////////////////////////////////////////////////////
 //   Cambia_Extremos_Barra(ndime,ntipo,nelem,fuefl);
//Escritura en texto para verificación de entrada de función
    int nprop,iprop, nmats=1, lmats,jelem, nbarr;
    int n_vigas, n_floors, f;
	nprop = (ndime == 2) ? 6 : 9;

   // FILE *efi_prueba, *fuerzas_excel;

	//efi_prueba=fopen("eficiencia.txt","w");
	//fuerzas_excel=fopen("fuerzas.xls","w");

    /*if (efi_prueba == NULL){
    printf("Error opening file!\n");
    }*/

    //Numero de elementos, dimension y tipo de barra
	//fprintf(efi_prueba,"nelem=%d, ndime=%d, ntipo=%d\n", nelem, ndime, ntipo);
    //fprintf(efi_prueba,"\n");

    //Numero de seccion utilizada en cada barra
    //fprintf(efi_prueba,"\nMateriales - Secciones: \n");
    for(ielem=1; ielem<=nelem; ielem++){
        //fprintf(efi_prueba,"matnu[%d]=%d \n", ielem,matnu[ielem]);

        if (matnu[ielem]>=nmats){
            nmats=matnu[ielem];
        }
    }
    // Tipo de barra
    //fprintf(efi_prueba,"\n");
    //fprintf(efi_prueba,"\nTipo de Barra: \n");
  //  for(ielem=1; ielem<=nelem; ielem++)
        //fprintf(efi_prueba,"ntips[%d]=%d \n",ielem, ntips[ielem]);

    //Propiedades de los materiales
    //fprintf(efi_prueba,"\n");
    //fprintf(efi_prueba,"\nPropiedades de Materiales: \n");
    // Lee y escribe propiedades de los materiales:
   // for(imats=1; imats<=nmats; imats++) {
        //fprintf(efi_prueba,"Mat[%d] -",imats);
        // Lee las propiedades antes del tipo de material
     //   for(iprop=1; iprop<=nprop; iprop++)
       //     fprintf(efi_prueba, " %lf \t", props[imats][iprop]);

        //fprintf(efi_prueba,"\n");
    

    //Longitud de cada barra
    //fprintf(efi_prueba,"\nLongitudes de Barra: \n");
    //for(ielem=1; ielem<=nelem; ielem++)
        //fprintf(efi_prueba,"xlong[%d]=%lf \n",ielem, xlong[ielem]);


    //Vector en ceros para guardar cada eficiencia
    //fprintf(efi_prueba,"\n");
    //fprintf(efi_prueba,"\nEspacio para eficiencias: \n");
    for(ielem=1; ielem<=nelem; ielem++){
        if(ntipo ==3) {
            if(ntips[ielem]==1)fuerc[ielem]=-fuefl[ielem][1];
        }
        //fprintf(efi_prueba,"fuerc[%d]=%lf \n ",ielem, fuerc[ielem]);
    }

/*
    //Propiedades de seccion en cada barra
    fprintf(efi_prueba,"\n");
    fprintf(efi_prueba,"\n Seccion de Concreto: \n");
    for(ielem=1; ielem<=nelem; ielem++){
        imats=matnu[ielem];
        iCat = ubiCat[imats];
        if(ndime== 2) lmats=(int)props[matnu[ielem]][6]+0.5;
        if(ndime ==3) lmats=(int)props[matnu[ielem]][9]+0.5;
        fprintf(efi_prueba,"b-arecr[%d][%d][1]=%lf ",ielem,matnu[ielem], arecr[iCat][lmats][1]);
        fprintf(efi_prueba,"h-arecr[%d][%d][1]=%lf ",ielem,matnu[ielem], arecr[iCat][lmats][2]);
        fprintf(efi_prueba,"fc-arecr[%d][%d][1]=%lf ",ielem,matnu[ielem], arecr[iCat][lmats][3]);
        fprintf(efi_prueba,"fy-arecr[%d][%d][1]=%lf \n",ielem,matnu[ielem], arecr[iCat][lmats][4]);
    }

*/
    //Indicador de las fuerzas que existen 1-0
/*
    fprintf(efi_prueba,"\n");
    fprintf(efi_prueba,"\n Indicador de Fuerzas (NF_nodos NF_punt [NFpunt_y] NF_dist): \n");
    for(ielem=1; ielem<=nelem; ielem++){
     for(jelem=1; jelem<=4; jelem++)
        fprintf(efi_prueba,"indfu[%d][%d]=%d ",ielem,jelem, indfu[jelem][ielem]);

    fprintf(efi_prueba,"\n");
    }
*/
    //Magnitud de las fuerzas que existen
    //fprintf(efi_prueba,"\n");
    //fprintf(efi_prueba,"\n Tipo de fuerza asignada a barra: ");
    //fprintf(efi_prueba,"\n 1 - Magnitud de la Fuerza Concentrada (Rev esb Mom Max en Z)");
    //fprintf(efi_prueba,"\n 2 - Magnitud de la Fuerza Concentrada (Rev esb Mom Max en Y)");
    //fprintf(efi_prueba,"\n 3 - Distancia de Aplicacion de Fuerza Concentrada (Rev esb Mom Max en Z)");
    //fprintf(efi_prueba,"\n 4 - Distancia de Aplicacion de Fuerza Concentrada (Rev esb Mom Max en Y)");
    //fprintf(efi_prueba,"\n 5 - Magnitud de la Fuerza Distribuida mfdy (Rev esb Mom Max en Z)");
    //fprintf(efi_prueba,"\n 6 - Magnitud de la Fuerza Distribuida mfdz (Rev esb Mom Max en Y)\n\n");
   /* for(ielem=1; ielem<=nelem; ielem++){
     for(jelem=1; jelem<=6; jelem++)
        //fprintf(efi_prueba,"fuerb[%d][%d]=%lf ",ielem,jelem, fuerb[jelem][ielem]);

    fprintf(efi_prueba,"\n");
    }*/

    //Fuerzas en extremo de barra
    //fprintf(efi_prueba,"\n");
    //fprintf(efi_prueba,"\n Fuerzas locales en extremo de barra: \n");
    /*for(ielem=1; ielem<=nelem; ielem++){
     for(jelem=1; jelem<=nevab; jelem++)
        fprintf(efi_prueba,"fuefl[%d][%d]=%lf ",ielem,jelem, fuefl[ielem][jelem]);

    fprintf(efi_prueba,"\n");
    }*/
    Cambia_Extremos_Barra(ndime,ntipo,nelem,fuefl);
   //fprintf(efi_prueba,"\n Fuerzas locales en extremo de barra: \n");
    /*for(ielem=1; ielem<=nelem; ielem++){
     for(jelem=1; jelem<=12; jelem++)
        fprintf(efi_prueba,"fuefl[%d][%d]=%lf ",ielem,jelem, fuefl[ielem][jelem]);
        fprintf(efi_prueba,"\n");
    }*/
    //Fuerzas en extremo de barra para excel
    //Solo funciona para numero de vigas fijo
//_______________________________________________________
/*
    n_vigas=24;
    n_floors=1;

    fprintf(fuerzas_excel,"\n Fuerzas locales en extremo de barra: \n");
    for(f=1; f<=n_floors; f++){

            for(ielem=1; ielem<=n_vigas; ielem++){

                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][1]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][3]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][2]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][4]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][6]);
                fprintf(fuerzas_excel,"%lf \n", fuefl[ielem][5]);

                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][7]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][9]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][8]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][10]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][12]);
                fprintf(fuerzas_excel,"%lf ", fuefl[ielem][11]);

                fprintf(fuerzas_excel,"\n");

            }

            for(ielem=n_vigas+1; ielem<=nelem; ielem++){

                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][1]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][3]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][2]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][4]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][6]);
                fprintf(fuerzas_excel,"%lf \n", fuefl[ielem][5]);

                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][7]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][9]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][8]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][10]);
                fprintf(fuerzas_excel,"%lf \t", fuefl[ielem][12]);
                fprintf(fuerzas_excel,"%lf ", fuefl[ielem][11]);

                fprintf(fuerzas_excel,"\n");

            }



    }
fclose(fuerzas_excel);
//_______________________________________________________
*/
    //Vector en ceros desconocido
    //fprintf(efi_prueba,"\n");
 //   fprintf(efi_prueba,"\n wcarg? - no se usa en eficiencias: \n");
 //   for(ielem=1; ielem<=2*nelem; ielem++)
 //       fprintf(efi_prueba,"wcarg[%d]=%lf \n",ielem, wcarg[ielem]);
/*
    fprintf(efi_prueba,"\n");
    fprintf(efi_prueba,"\n Información de acero de refuerzo: \n");
    for(ielem=1; ielem<=nelem; ielem++){
        imats=matnu[ielem];
        iCat = ubiCat[imats];
        if(ndime== 2) lmats=(int)props[matnu[ielem]][6]+0.5;
        if(ndime ==3) lmats=(int)props[matnu[ielem]][9]+0.5;
        nbarr=(int)arecr[iCat][lmats][5]+0.5;
        for(iprop=1; iprop<=nbarr; iprop++)
            fprintf(efi_prueba,"dtcon[%d][%d][1]=%lf ",ielem,iprop, dtcon[iCat][lmats][1][iprop]);
        fprintf(efi_prueba,"\n");
        for(iprop=1; iprop<=nbarr; iprop++)
            fprintf(efi_prueba,"dtcon[%d][%d][2]=%lf ",ielem,iprop, dtcon[iCat][lmats][2][iprop]);
        fprintf(efi_prueba,"\n");
        for(iprop=1; iprop<=nbarr; iprop++)
            fprintf(efi_prueba,"dtcon[%d][%d][3]=%lf ",ielem,iprop, dtcon[iCat][lmats][3][iprop]);

    fprintf(efi_prueba,"\n");
    fprintf(efi_prueba,"\n");
    }
 */
    //fprintf(efi_prueba,"\n");
    //for(ielem=1; ielem<=10; ielem++)
     //   fprintf(efi_prueba,"valor[%d]=%lf ",ielem, valor[ielem]);

///////////////////////////////////////////////////////////////////////////////




	for(ielem=1; ielem<=nelem; ielem++){
		imats=matnu[ielem];
		// Identifica el material del elementos
		//if(ndime == 2) rolado=(int)props[imats][5]+0.5;
		//if(ndime == 3) rolado=(int)props[imats][8]+0.5;
        rolado=(int)props[imats][8]+0.5;
		// Identifica el catalogo que le pertence
	 	iCat = ubiCat[imats];
        if(iefi ==1 && (ntipo ==1 || ntipo ==3) && rolado==0){
 //           printf("entro armadura elemento =%d iefi=%d  rolado=%d imats =%d props[imats][8]= %lf \n",ielem,iefi,rolado,imats,props[imats][8]);
			Eficiencia_Armaduras(ielem,ndime,matnu,props,xlong,
								 arear[iCat],fuerc);
        }
        if(iefi ==1 && (rolado ==1 || rolado ==4 )){
 //           printf("entro rolado en frio elemento =%d iefi=%d  rolado=%d  imats =%d props[imats][8]= %lf \n",ielem,iefi,rolado,imats,props[imats][8]);
			Eficiencia_MarcosRF(ielem,nelem,ndime,ntipo,matnu,ntips,props,xlong,
								fuerc,arerf[iCat],indfu,fuerb,fuefl,
								wcarg,rolado);
        }
        if(iefi ==1 && rolado ==2){
//            printf("entro rolado en caliente elemento =%d  iefi=%d  rolado=%d  imats =%d props[imats][8]= %lf\n",ielem,iefi,rolado,imats,props[imats][8]);
			Eficiencia_MarcosRC(ielem,nelem,ndime,ntipo,matnu,ntips,props,xlong,
								fuerc,arerc[iCat],indfu,fuerb,fuefl,wcarg);
        }
        if(iefi ==1 && rolado ==3){
//            printf("entro concreto elemento =%d  iefi=%d  rolado=%d  imats =%d props[imats][8]= %lf\n",ielem,iefi,rolado,imats,props[imats][8]);
			Eficiencia_MarcosCR(ielem,nelem,ndime,ntipo,matnu,ntips,props,xlong,
								fuerc,arecr[iCat],indfu,fuerb,fuefl,wcarg,
								dtcon[iCat]);
        }
	}
	//if((iefi ==1 && ntipo >1 && rolado==0))
	//	fprintf(fp16,"imposible calcular eficiencias\n");
    //else {
		valef=0.0;
		// Obtiene las eficiencias maximas de los elementos
		for(ielem=1; ielem<=nelem; ielem++){
        	if(valef < fabs(fuerc[ielem]))  valef = fabs(fuerc[ielem]);
        	if(efi_max[ielem] < fabs(fuerc[ielem])) efi_max[ielem] = fabs(fuerc[ielem]);
		}
        if(iwrit != 0){
			for(ielem=1; ielem<=nelem; ielem++){
				fprintf(fp16,"elemento = %d eficiencia =%f\n",ielem,fuerc[ielem]);
			}
        }
		//fprintf(fp16,"Eficiecia maxima en la estructura = %lf\n\n",valef);
		//printf("Eficiecia maxima en la estructura = %lf\n\n",valef);
        // Verifica si la eficiencia de este caso de carga supera a la anterior
		if( valef > valor[3] ) valor[3] = valef;
	//}



///////////////////////////////////////////////////////////////////////////////

    //fprintf(efi_prueba,"\n\n");
    //fprintf(efi_prueba,"\nEficiencias ya calculadas: \n");
    /*for(ielem=1; ielem<=nelem; ielem++)
        fprintf(efi_prueba,"fuerc[%d]=%g \n",ielem, fuerc[ielem]);*/

    //fprintf(efi_prueba,"\n");
    //fprintf(efi_prueba,"\n");
	//fclose(efi_prueba);

///////////////////////////////////////////////////////////////////////////////

}
//***********************************************************************************************************
void Cambia_Extremos_Barra(int ndime,int ntipo,int nelem, double **fuefl)
{
    int ielem,ii;
    double temp[13];
    if(ndime ==2)
       for(ielem=1; ielem<=nelem; ielem++){
           if(ntipo==1) {
               temp[1]=fuefl[ielem][1];
               temp[2]=0;
               temp[3]=fuefl[ielem][2];
               temp[4]=0;
               temp[5]=0;
               temp[6]=0;
               temp[7]=fuefl[ielem][3];
               temp[8]=0;
               temp[9]=fuefl[ielem][4];
               temp[10]=0;
               temp[11]=0;
               temp[12]=0;
             }
           if(ntipo==2) {
               temp[1]=fuefl[ielem][1];
               temp[2]=0;
               temp[3]=fuefl[ielem][2];
               temp[4]=0;
               temp[5]=fuefl[ielem][3];
               temp[6]=0;
               temp[7]=fuefl[ielem][4];
               temp[8]=0;
               temp[9]=fuefl[ielem][5];
               temp[10]=0;
               temp[11]=fuefl[ielem][6];
               temp[12]=0;
            }
          if(ntipo==3) {
               temp[1]=fuefl[ielem][1];
               temp[2]=0;
               temp[3]=fuefl[ielem][2];
               temp[4]=0;
               temp[5]=fuefl[ielem][3];
               temp[6]=0;
               temp[7]=fuefl[ielem][4];
               temp[8]=0;
               temp[9]=fuefl[ielem][5];
               temp[10]=0;
               temp[11]=0;
               temp[12]=0;
            }
          for(ii=1; ii<=12; ii++)
            fuefl[ielem][ii]=temp[ii];
    }
}



