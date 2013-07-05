#ifdef _SS_OIL_WATER_H
#define SS_OIL_WATER_H
#include <string.h>
#include "phg.h"
#include <math.h>
extern  FLOAT T;
extern  FLOAT stime;    			//time step size
extern  FLOAT ctime;
extern  FLOAT ptime;     		//previous time 
/*some parameters for oil*/
extern  FLOAT MU_O;          //the viscosity
extern  FLOAT B0_O;           		//the init of reservoir volume ratio (bo)
extern  FLOAT C_O;          
/*some parameters for water*/
extern  FLOAT MU_W;
extern  FLOAT B0_W;
extern  FLOAT C_W;
/*some parameters for poros medium*/
extern  FLOAT K;          		//absolute permeability
extern  FLOAT C_R;        
extern  FLOAT PHI0;          		//porosity 
/*some parameters for permeability*/
extern  FLOAT SORW;
extern  FLOAT KRW_SORW;
extern  FLOAT KRO_SWC;
extern  FLOAT SWC;
extern  INT NO;
extern  INT NW;
/*some parameters for well*/
extern  FLOAT PRODUCTION;
/*Initnal Data for unknown variable*/
extern  FLOAT PRESSURE0;           //init pressure
extern  FLOAT SW0;                 //init water saturation

#define USE_PC 1
#define DEBUG 0
/**********************************
	       well.c                    
***********************************/
void Well_init(DOF *source);
void Well_Pressure(DOF *p_h, FLOAT *pressure);
/**********************************
     	  uzawa.c                    
***********************************/
DOF *DivRT(DOF *u, DOF *p);
MAT *BTAB(MAT * A, DOF *B, DOF *u, DOF *p);
void pc_proc(SOLVER *pc_solver, VEC *b0, VEC **x0);
FLOAT phgUzawa(SOLVER *pc_solver, MAT *H0, DOF *u, DOF *p, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g, INT *nits);
/**********************************
	       parameter.c                
**********************************/
void create_kro(DOF *s_w, DOF *kro);
void create_krw(DOF *s_w, DOF *krw);
void update_bo(DOF *p_h, DOF *p0_h, DOF *b_o);
void update_bw(DOF *p_h, DOF *p0_h, DOF *b_w);
void update_phi(DOF *p_h, DOF *p0_h, DOF *phi);
void create_dot_phibo(DOF *p_h, DOF *p0_h, DOF *dot_phibo);
void create_dot_phibw(DOF *p_h, DOF *p0_h, DOF *dot_phibw);
void create_phi_bo(DOF *phi, DOF *b_o, DOF *phi_bo);
void create_phi_bw(DOF *phi, DOF *b_w, DOF *phi_bw);
/**********************************
		 quadfunc.c
**********************************/
FLOAT phgQuadFacefunc(SIMPLEX *e, int face, DOF *u, int order);
FLOAT phgQuadDofTimesDivBas(SIMPLEX *e, DOF *u, DOF *v, int m, int order);
FLOAT L2_norm(DOF *u);
FLOAT L2(DOF *a, DOF *b, int order);
#endif /*end of ss_oil_water.h*//
