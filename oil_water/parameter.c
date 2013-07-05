#include "phg.h"
#include <string.h>
#include <math.h>
//#include "ss_oil_water.h"
/*paremater function for oil phase*/
extern  FLOAT T= 720;
extern  FLOAT stime = 1;    			//time step size
extern  FLOAT ctime = 0;
extern  FLOAT ptime = 0.0;     		//previous time 
/*some parameters for oil*/
extern  FLOAT MU_O =1e-9 / 3600.;          //the viscosity
extern  FLOAT B0_O = 1.05;           		//the init of reservoir volume ratio (bo)
extern  FLOAT C_O = 2e-3;          
/*some parameters for water*/
extern  FLOAT MU_W = 1e-10 / 3600;
extern  FLOAT B0_W = 1.02;
extern  FLOAT C_W = 1e-3;
/*some parameters for poros medium*/
extern  FLOAT K = 1e-13;          		//absolute permeability
extern  FLOAT C_R = 0.15e-3;        
extern  FLOAT PHI0 = 0.2;          		//porosity 
/*some parameters for permeability*/
extern  FLOAT SORW = 0.;
extern  FLOAT KRW_SORW = 0.6;
extern  FLOAT KRO_SWC = 1;
extern  FLOAT SWC = 0.;
extern  INT NO = 4;
extern  INT NW = 3;
/*some parameters for well*/
extern  FLOAT PRODUCTION;
/*Initnal Data for unknown variable*/
extern  FLOAT PRESSURE0 = 20.;           //init pressure
extern  FLOAT SW0 = 0.4;                 //init water saturation
extern  FLOAT SO0 = 0.6;                 //init water saturation
void
create_kro(DOF *s_w, DOF *kro)
{
	GRID *g = s_w->g;
	SIMPLEX *e;
	FLOAT *p_sw, *p_kro;
	ForAllElements(g, e){
		p_sw = DofElementData(s_w, e->index);
		p_kro = DofElementData(kro, e->index);
		p_kro[0] = KRO_SWC * pow((1 - p_sw[0] - SORW) / (1 - SWC - SORW), NO);
	}
}
void
create_dot_kro(DOF *s_w, DOF *dot_kro)
{
	GRID *g = s_w->g;
	SIMPLEX *e;
	FLOAT *p_sw, *p_dot;
	ForAllElements(g, e){
		p_sw = DofElementData(s_w, e->index);
		p_dot = DofElementData(dot_kro, e->index);
		p_dot[0] = NO * (-1. / (1. - SWC -SORW)) * KRO_SWC * pow((1 - p_sw[0] - SORW) / (1 - SWC - SORW), NO-1);
	}
}
void
update_bo(DOF *p_h, DOF *p0_h, DOF *b_o)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_bo, *p_p0;
	ForAllElements(g, e){
		p_p0 = DofElementData(p0_h, e->index);
		p_p = DofElementData(p_h, e->index);
		p_bo = DofElementData(b_o, e->index);
		p_bo[0] = (1.0 + C_O * (p_p[0] - p_p0[0])) / B0_O;
	}
}
void
update_phi(DOF *p_h, DOF *p0_h, DOF *phi)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_phi, *p_p0;
	ForAllElements(g, e){
		p_p0 = DofElementData(p0_h, e->index);
		p_p = DofElementData(p_h, e->index);
		p_phi = DofElementData(phi, e->index);
		p_phi[0] = PHI0 * (1.0 + C_R * (p_p[0] - p_p0[0]));
	}
}
/*paremater function for water phase*/
void
create_krw(DOF *s_w, DOF *krw)
{
	GRID *g = s_w->g;
	SIMPLEX *e;
	FLOAT *p_sw, *p_krw;
	ForAllElements(g, e){
		p_sw = DofElementData(s_w, e->index);
		p_krw = DofElementData(krw, e->index);
		p_krw[0] = KRW_SORW * pow((p_sw[0] - SWC) / (1 - SWC - SORW), NW);
	}
}
void
create_dot_krw(DOF *s_w, DOF *dot_krw)
{
	GRID *g = s_w->g;
	SIMPLEX *e;
	FLOAT *p_sw, *p_dot;
	ForAllElements(g, e){
		p_sw = DofElementData(s_w, e->index);
		p_dot = DofElementData(dot_krw, e->index);
		p_dot[0] = NW * (1. / (1. - SWC - SORW)) * KRW_SORW * pow((p_sw[0] - SWC) / (1 - SWC - SORW), NW-1);
	}
}
void
update_bw(DOF *p_h, DOF *p0_h, DOF *b_w)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_bw, *p_p0;
	ForAllElements(g, e){
		p_p0 = DofElementData(p0_h, e->index);
		p_p = DofElementData(p_h, e->index);
		p_bw = DofElementData(b_w, e->index);
		p_bw[0] = (1.0 + C_W * (p_p[0] - p_p0[0])) / B0_W;
	}
}
