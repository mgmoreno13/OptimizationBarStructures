#ifndef _Optimizador_h_
#define _Optimizador_h_
#include <vector>
#include <stdio.h>

using namespace std; 
int Optimizador( int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                double** coord, double** carga, double** carpp, double** presc, double** props,
                double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                int nwktl,int prueba);
/*void Genera_Nueva_Poblacion( int nmats, int ndime, int* matva, double** props,
                            int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
                            int *Con_tamCat, int *ubiCat );*/
void Genera_Poblacion_Inicial(int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,//@mg 
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                 int nwktl);
void inicializa_array1D( int n, double *v, double cte );
void escribe_efi_max( int nelem, double *efi_max );
//void Recosido_Simulado(int Num_iter);
void Genera_Nuevo_Individuo( int nmats, int ndime, int* matva, double** props,//@mg 
                             int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
                             int *Con_tamCat, int *ubiCat );
double Funcion_fitness( int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,//@mg 
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl);

void parent_selection();//@mg 
void Crossover();
void mutacion(int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl);//@mg 
int new_pob();    //@mg     
void actualizacion(int* matvm,double* valom);//@mg 
void local_search(int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl);//@mg


struct Individuo //@mg 
{   
    std::vector<int> secciones;//matva_modif
   // std::vector<int> barras;//matnu_modif

    std::vector<double> valores;
    double fit;
    double desp;
    double efi;
            

    bool operator<(Individuo a) {

        return fit < a.fit;
    }

};

double distancia_inicial(std::vector<Individuo> Pob);
double distancia_individuos(Individuo a, Individuo b);


void actualizar_vecino_poblacion(Individuo vecino,int i);//@mg
void DFS_DT(int v, int padre);


#endif
