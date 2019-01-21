#ifndef _principal_h_
#define _principal_h_

//*****************************
// Rutina Principal
//*****************************


//******************************************************************************
//							Principal
//******************************************************************************
	// Finaliza 0 si no hubo ningun error.
	// Finaliza -1 si se hubo algun error de calculo.
	// Finaliza -2 si falto algun dato es erroneo.
	// Finaliza -3 si no se cuenta con la version completa (o la copia de evaluacion).
	// Finaliza -4 si el archivo no es de datos.
	// Finaliza -5 si no se pudo crear el archivo de resultados.
int Principal(int indso, int ishot, char* input1, char* output1, bool Optimiza,int prueba);

void Actualiza(int nmats, int* matva, int* matvm, double* valor, double* valom);

void Recupera(int nmats, int* matva, int* matvm, double* valor, double* valom);

void Inicializa_Materiales( int nmats, int ndime, double** props, int* matva,
							int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
							int *Con_tamCat, int *ubiCat );


#endif
