#include <stdio.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include "main.h"
#include "principal.h"
#include "raros.h"


int main(int argc, char **argv)
{
	chdir("/home/maria.moreno/meca_Paso_reemp_bc_1000");
	int prueba = strtol(argv[1], NULL, 10);


		//printf("%d\n",prueba );

    int  indsal,indso,ishot;
    char input1[256], output1[256];

    /* Solver 2 funciona bien por el uso del vector stiff
        si se cambia por otro solver se deberá comentar la
        liberación del espacio
    */
    indso = 2;//mg indice de solucion *(manera de resolverlo)
    ishot = 0;

    // Obtiene los nombres del archivo de datos y el de resultados
    fp = fopen("coman.dat","rt") ;
    fscanf(fp,"%s",input1) ;
    fscanf(fp,"%s",output1) ;
    fclose(fp) ;

    // Hace el analisis estatico y/o la optimizacion
	Probini();//inicia conteo de tiempo
	indsal = Principal(indso,ishot,input1, output1, false,prueba);

    // Cierra los flujos a archivos
	if( indsal == 0 ){
		CierraLimpia();
	}
	return(0);

}
//·············································································
void Probini(void)
{
	a= time(NULL);
}
//·············································································
int Prob(void)
{
	b = time(NULL);
	return (int)(difftime(b, a)+0.5);
}
//·············································································
void Comult(int l, int m, int n, double** a, double** b, double** c)
{
	register int i, j, k;

	for(i=1; i<=l; i++) {
		for(j=1; j<=n; j++) {
			a[i][j] = 0.0;

			for(k=1; k<=m; k++) a[i][j] = a[i][j]+b[i][k]*c[k][j];
        }
	}
}

//·············································································
void CierraLimpia(void)
{
	fclose(fp5);
	fclose(fp16);
    //fclose(fp79);

    //	remove(nom[0]);
    //	remove(nom[1]);
    //	remove(nom[2]);
    //	remove(nom[3]);
}
//·············································································
void ConvSegHrs(long int segs, char* hora)
{
	int hrs, mins;

	// Obtiene las horas:
	//hrs = int(segs/3600);
	hrs = (int) segs/3600;

	// Obtiene lo minutos:
	segs-= hrs*3600;
	//mins = int(segs/60);
	mins = (int)segs/60;

	// Obtiene los segundos restantes:
	segs-= mins*60;

	// Completa la cadena de tiempo:
	sprintf(hora, " Tiempo total de ejecucion :   %d:%d:%ld", hrs, mins, segs);
	//W1.printf(10,300,"Termino programa exitosamente en :   %d:%d:%d", hrs, mins, segs);
}
//·············································································
