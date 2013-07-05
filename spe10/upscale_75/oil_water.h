#ifdef _OIL_WATER_H
#define _OIL_WATER_H
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#endif  /*endof oil_water.h*/
#include "phg.h"
/*Define Some Constants*/
#define USE_UZAWA 0
#define USE_BLOCK 1
#define EXPORT_VTK 0
#define REFINE_WELL 0
#define USE_BLOCK_MATRIX 0
/********Scale Parameters***********/
#define PRESS_FACTOR 0.006894757  //psi-->MPs
#define LEGTH_FACTOR 0.3048       // ft-->m
#define VOLUME_FACTOR 0.159		  // bbl-->m^3
#define VIS_FACTOR     1e-9/(3600*24)  //cp->MPa \cdot day
typedef struct phase PHASE;
struct phase
{
	/*parameters for phase*/
	FLOAT Mu;               /*the viscosity*/
	DOF *B;					/*the reservoir volume ratio*/
	DOF *DB;					/*the reservoir volume ratio*/
	DOF *Kr; 				/*relavant permeability*/
	DOF *DKr;            /*\frac{\paritial kr}{\partial s}*/
	DOF *S;                /*phase saturation*/
	FLOAT S0;
	FLOAT PRESSURE0;
	DOF *P;                /*phase pressue*/
	DOF *P0;
	DOF *Source;              /*phase production*/
	DOF *U;                /*phase velocity*/
};

typedef struct medium MEDIUM;
struct medium
{
	/*parameters for poros medium*/
	DOF *phi;            /*the rock porosity*/
	DOF *phi0;
	DOF *Dphi;
	FLOAT C_R;
	//FLOAT perm;           /*the absolute permeability*/
	DOF *perm;
	INT NX;
	INT NY;
	INT NZ;
	FLOAT DX;
	FLOAT DY;
	FLOAT DZ;
};
typedef struct well WELL;
struct well
{
	FLOAT PROD;
	FLOAT INJE;
	FLOAT MAX_INJE_BHP;
	FLOAT PROD_BHP;
	FLOAT *radius;
	FLOAT *eff_radius;
	FLOAT skin;
	FLOAT *position;
	FLOAT *region;
	INT Number;
};
typedef struct computate_control COM_TRO;
struct computate_control
{
	FLOAT T;
	FLOAT dt;
	FLOAT ct;
	FLOAT pt;
	FLOAT max_dt;
	FLOAT init_dt;
	FLOAT TOL_sys;
	FLOAT TOL_non;
	FLOAT TOL_con;
	FLOAT MAX_DS;
	FLOAT MAX_DP;
};
/**********************************
	       well.c                    
***********************************/
void 
Well_init(PHASE *oil, WELL *well);
void 
Well_Pressure(DOF *p_h, FLOAT *pressure);
void Mark_All_Well(GRID *g, MEDIUM *rock, WELL *well);
/**********************************
     	  uzawa.c                    
***********************************/
DOF *DivRT(DOF *u, DOF *p, COM_TRO *control);
MAT *BTAB(MAT * A, DOF *B, DOF *u, DOF *p);
void pc_proc(SOLVER *pc_solver, VEC *b0, VEC **x0);
FLOAT phgUzawa(SOLVER *pc_solver, MAT *H0, DOF *u, DOF *p, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g, INT *nits);
/**********************************
	       parameter.c                
**********************************/
void create_kr(PHASE *oil, PHASE *water);
void create_dot_kr(PHASE *oil, PHASE *water);
void update_phi(MEDIUM *rock, PHASE *oil);
void update_dot_phi(MEDIUM *rock);
void update_bo(PHASE *oil);
void update_dot_bo(PHASE *oil);
void update_bw(PHASE *water);
void update_dot_bw(PHASE *water);
void
Init_Medium(const char *file, MEDIUM *rock);
/**********************************
		 quad.c
**********************************/
FLOAT phgQuadBasPermBas(SIMPLEX *e, DOF *u, int n, DOF *perm, DOF *v, int m, int order);
FLOAT phgQuadBasKDDDBas(SIMPLEX *e, DOF *u, int n, DOF *D1, DOF *D2, DOF *D3, DOF *v, int m, int order);
FLOAT phgQuadBasAAABas(SIMPLEX *e, DOF *u, int n, DOF *A1, DOF *A2, DOF *A3, DOF *v, int m, int order);
FLOAT *phgQuadDofAAABas(SIMPLEX *e, DOF *u, DOF *A1, DOF *A2, DOF *A3, DOF *v, int n, int order, FLOAT *res);
FLOAT phgQuadDivBasDotBas(SIMPLEX *e, DOF *u, int m, DOF *v, int n, int order);
FLOAT phgQuadFacefunc(SIMPLEX *e, int face, DOF *u, int order);
FLOAT phgQuadDofTimesDivBas(SIMPLEX *e, DOF *u, DOF *v, int m, int order);
FLOAT L2_norm(DOF *u);
FLOAT L2(DOF *a, DOF *b, int order);
/******************************************
    RT1 Basis Function
******************************************/
DOF_TYPE DOF_RT1_;
#define DOF_RT1 (&DOF_RT1_)
# define BasisOrder(u, e, i) (!DofIsHP(u) ? (u)->type->order : \
	(u)->hp->info->types[(u)->hp->elem_order[e->index]]->order)
/****************************************
*	Read Data File
****************************************/
void
read_data(const char *file, PHASE *oil, PHASE *water, WELL *well, COM_TRO *control, MEDIUM *rock);
void
units_convert(PHASE *oil, PHASE *water, WELL *well, MEDIUM *rock);
/****************************************
			export.c
****************************************/
void export_vtkfile(const char *vtkfile, PHASE *oil, PHASE *water, MEDIUM *rock, int flag);
void export_test_data(PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, COM_TRO *control, WELL *well);
void export_test_parameters(PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control, WELL *well);
void
result_export(PHASE *oil, PHASE *water, WELL *well, COM_TRO *control, FLOAT *total_oil);
/*****************************************
			struct.c
*****************************************/
PHASE *phgGetPhase(GRID *g);
MEDIUM *phgGetMedium(GRID *g);
COM_TRO *phgGetComTro(GRID *g);
WELL *phgGetWellParameters(GRID *g);
void phgFreePhase(PHASE *_phase);
void phgFreeMedium(MEDIUM *_medium);
void phgFreeComTro(MEDIUM *_control);
void phgFreeWell(WELL *_well);
