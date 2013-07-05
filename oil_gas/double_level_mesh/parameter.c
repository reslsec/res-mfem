#include "phg.h"
#include <string.h>
#include <math.h>
#include "oil_gas.h"

static int RANGE = 4;

int 
Interval_NO(FLOAT pressure, BOOLEAN flags)
{
	int i;
	static int interval_no = 0;

	
	if(!flags)
	{
		if(pressure >= PVT[interval_no * SUM_PHY] 
				&& pressure <= PVT[interval_no * SUM_PHY + SUM_PHY])
			return interval_no;
		for(i = 0; i < RANGE; i++)
		{
			if((interval_no + i + 2 > SUM_INTERVAL + 1)
					|| (interval_no - i - 1 < 0))
				break;

			if(pressure <= PVT[interval_no * SUM_PHY + (i + 2) * SUM_PHY]
				&& pressure >= PVT[interval_no * SUM_PHY + (i + 1) * SUM_PHY])
			{
				interval_no = interval_no + i + 1;
				return interval_no;
			}

			if(pressure >= PVT[interval_no * SUM_PHY - (i + 1) * SUM_PHY]
				&& pressure <= PVT[interval_no * SUM_PHY - i * SUM_PHY])
			{
				interval_no = interval_no - i - 1;
				return interval_no;
			}
		}
		flags = TRUE;
	}

	if(flags)
	{
		int max = SUM_INTERVAL, min = 0;
		interval_no = (int) (max + min) / 2;
		while((pressure <= PVT[interval_no * SUM_PHY])
				|| (pressure > PVT[interval_no * SUM_PHY + SUM_PHY]))
		{
			if(pressure > PVT[interval_no * SUM_PHY + SUM_PHY])
				min = interval_no;
			if(pressure <= PVT[interval_no * SUM_PHY])
				max = interval_no;
			interval_no = (int) (max + min) / 2;
			if(interval_no == min)
				break;
		}

		if(pressure > PVT[interval_no * SUM_PHY] 
				&& pressure <= PVT[interval_no * SUM_PHY + SUM_PHY])
			return interval_no;

	}

	phgError(1, "Can not find %lf\n", pressure);
}
/*Unfiorm unit to MPa, m, h*/
void
Unit_Conversion(void)
{
	int i = 0;
	/*For Pressure*/
	for (i = 0; i < 1000; i++){
		PVT[P_INDEX(i)] = PVT[P_INDEX(i)] / 1e6;
	}
	/*For mu_o and mu_g*/
	for (i = 0; i < 1000; i++){
		PVT[MUO_INDEX(i)] = PVT[MUO_INDEX(i)] / (1e6 * 3600.);
		PVT[MUG_INDEX(i)] = PVT[MUG_INDEX(i)] / (1e6 * 3600.);
	}
}
void
Test_PVT_Unit(void)
{
 	int i = 0;
	for (i = 0; i < 20; i++){
		printf("pressure = %lf\n", PVT[P_INDEX(i)]);
	}
	for (i = 0; i < 20; i++){
		printf("mu_o = %le\n", PVT[MUO_INDEX(i)]);
	}
}
void
update_muo(DOF *p_h, DOF *mu_o)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_mu;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_mu = DofElementData(mu_o, e->index);
		NO = Interval_NO(p_p[0], 0);
		FLOAT temp; 
		temp = (p_p[0] - PVT[P_INDEX(NO)]) / (PVT[P_INDEX(NO + 1)] - PVT[P_INDEX(NO)]);
		p_mu[0] = 1. / (PVT[MUO_INDEX(NO)] + temp * (PVT[MUO_INDEX(NO + 1)] - PVT[MUO_INDEX(NO)]));
	}
}
void
update_bo(DOF *p_h, DOF *b_o)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_bo;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_bo = DofElementData(b_o, e->index);
		NO = Interval_NO(p_p[0], 0);
		FLOAT temp; 
		temp = (p_p[0] - PVT[P_INDEX(NO)]) / (PVT[P_INDEX(NO + 1)] - PVT[P_INDEX(NO)]);
//		p_bo[0] = 1.0 / PVT[BO_INDEX(NO)] + temp * (1.0 / PVT[BO_INDEX(NO + 1)] - 1.0 / PVT[BO_INDEX(NO)]);
		p_bo[0] = 1.0 / (PVT[BO_INDEX(NO)] + temp * ( PVT[BO_INDEX(NO + 1)] -  PVT[BO_INDEX(NO)]));
	}
}
void
update_mug(DOF *p_h, DOF *mu_g)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_mu;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_mu = DofElementData(mu_g, e->index);
		NO = Interval_NO(p_p[0], 0);
		FLOAT temp; 
		temp = (p_p[0] - PVT[P_INDEX(NO)]) / (PVT[P_INDEX(NO + 1)] - PVT[P_INDEX(NO)]);
		p_mu[0] = 1. / (PVT[MUG_INDEX(NO)] + temp * (PVT[MUG_INDEX(NO + 1)] - PVT[MUG_INDEX(NO)]));
	}
}
void
update_bg(DOF *p_h, DOF *b_g)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_bg;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_bg = DofElementData(b_g, e->index);
		NO = Interval_NO(p_p[0], 0);
		FLOAT temp; 
		temp = (p_p[0] - PVT[P_INDEX(NO)]) / (PVT[P_INDEX(NO + 1)] - PVT[P_INDEX(NO)]);
	//	p_bg[0] = 1.0 / PVT[BG_INDEX(NO)] + temp * (1.0 / PVT[BG_INDEX(NO + 1)] - 1.0 / PVT[BG_INDEX(NO)]);
		p_bg[0] = 1.0 / (PVT[BG_INDEX(NO)] + temp * (PVT[BG_INDEX(NO + 1)] - PVT[BG_INDEX(NO)]));
	}
}
void
update_Rs(DOF *p_h, DOF *Rs)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_rs;
	int NO = 0;;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_rs = DofElementData(Rs, e->index);
		NO = Interval_NO(p_p[0], 0);
		FLOAT temp; 
		temp = (p_p[0] - PVT[P_INDEX(NO)]) / (PVT[P_INDEX(NO + 1)] - PVT[P_INDEX(NO)]);
		p_rs[0] = PVT[RS_INDEX(NO)] + temp * (PVT[RS_INDEX(NO + 1)] - PVT[RS_INDEX(NO)]);
	}
}
void
update_dot_bo(DOF *p_h, DOF *dot_bo)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dotbo;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dotbo = DofElementData(dot_bo, e->index);
		NO = Interval_NO(p_p[0], 0);
		p_dotbo[0] = (1.0 / PVT[(BO_INDEX(NO + 1))] - 1.0 / PVT[(BO_INDEX(NO))])
			   		/ (PVT[(P_INDEX(NO + 1))] - PVT[(P_INDEX(NO))]);
	}
}
void
update_dot_bg(DOF *p_h, DOF *dot_bg)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dotbg;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dotbg = DofElementData(dot_bg, e->index);
		NO = Interval_NO(p_p[0], 0);
		p_dotbg[0] = (1.0 / PVT[(BG_INDEX(NO + 1))] - 1.0 / PVT[(BG_INDEX(NO))])
			   		/ (PVT[(P_INDEX(NO + 1))] - PVT[(P_INDEX(NO))]);
	}
}
void
update_dot_muo(DOF *p_h, DOF *dot_muo)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dot;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dot = DofElementData(dot_muo, e->index);
		NO = Interval_NO(p_p[0], 0);
		p_dot[0] = (1.0 / PVT[(MUO_INDEX(NO + 1))] - 1.0 / PVT[(MUO_INDEX(NO))])
			   		/ (PVT[(P_INDEX(NO + 1))] - PVT[(P_INDEX(NO))]);
	}
}
void
update_dot_mug(DOF *p_h, DOF *dot_mug)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dot;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dot = DofElementData(dot_mug, e->index);
		NO = Interval_NO(p_p[0], 0);
		p_dot[0] = (1.0 / PVT[(MUG_INDEX(NO + 1))] - 1.0 / PVT[(MUG_INDEX(NO))])
			   		/ (PVT[(P_INDEX(NO + 1))] - PVT[(P_INDEX(NO))]);
	}
}
void
update_dot_Rs(DOF *p_h, DOF *dot_Rs)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dotRs;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dotRs = DofElementData(dot_Rs, e->index);
		NO = Interval_NO(p_p[0], 0);
		p_dotRs[0] = (PVT[RS_INDEX(NO + 1)] - PVT[(RS_INDEX(NO))])
			        / (PVT[(P_INDEX(NO + 1))] - PVT[(P_INDEX(NO))]);
	}
}

void
create_krg(DOF *s_o, DOF *krg)
{
	GRID *g =s_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_krg;
	ForAllElements(g, e){
		p_so = DofElementData(s_o, e->index);
		p_krg = DofElementData(krg, e->index);
		if (e->region_mark == 0 | e->region_mark == 2){
			p_krg[0] = KRG_SWC_1 * pow((1. - p_so[0] - SGC_1) / (1. - SLC_1 - SGC_1), NG_1); 
		}
		else if (e->region_mark == 1 | e->region_mark == 3){
			p_krg[0] = KRG_SWC_2 * pow((1. - p_so[0] - SGC_2) / (1. - SLC_2 - SGC_2), NG_2); 
		}
	}
}
void
create_dot_krg(DOF *s_o, DOF *dot)
{
	GRID *g =s_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_dot;
	ForAllElements(g, e){
		p_so = DofElementData(s_o, e->index);
		p_dot = DofElementData(dot, e->index);
		if (e->region_mark == 0 | e->region_mark == 2){
			p_dot[0] = KRG_SWC_1 * (-1. / (1. - SLC_1 - SGC_1)) * NG_1 * pow((1. - p_so[0] - SGC_1) / (1. - SLC_1 - SGC_1), NG_1-1); 
		}
		else if (e->region_mark == 1 | e->region_mark == 3){
			p_dot[0] = KRG_SWC_2 * (-1. / (1. - SLC_2 - SGC_2)) * NG_2 * pow((1. - p_so[0] - SGC_2) / (1. - SLC_2 - SGC_2), NG_2-1); 
		}
	}
}
void
create_kro(DOF *s_o, DOF *kro)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_kro;
	ForAllElements(g, e){
		p_so = DofElementData(s_o, e->index);
		p_kro = DofElementData(kro, e->index);
		if (e->region_mark == 0 | e->region_mark == 2){
			p_kro[0] = KRO_SGC_1 * pow((p_so[0] - SLC_1) / (1. - SGC_1 - SLC_1), NO_1); 
		}
		else if (e->region_mark == 1 | e->region_mark == 3){
			p_kro[0] = KRO_SGC_2 * pow((p_so[0] - SLC_2) / (1. - SGC_2 - SLC_2), NO_2); 
		}
	}
}
void
create_dot_kro(DOF *s_o, DOF *dot)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_dot;
	ForAllElements(g, e){
		p_so = DofElementData(s_o, e->index);
		p_dot = DofElementData(dot, e->index);
		if (e->region_mark == 0 | e->region_mark == 2){
			p_dot[0] = KRO_SGC_1 * (1. / (1. - SGC_1 - SLC_1)) * NO_1 * pow((p_so[0] - SLC_1) / (1. - SGC_1 - SLC_1), NO_1-1); 
		}
		else if (e->region_mark == 1 | e->region_mark == 3){
			p_dot[0] = KRO_SGC_2 * (1. / (1. - SGC_2 - SLC_2)) * NO_2 * pow((p_so[0] - SLC_2) / (1. - SGC_2 - SLC_2), NO_2-1); 
		}
	}
}
void
update_phi(DOF *p_h, DOF *p0_h, DOF *phi)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_p0, *p_phi;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_p0 = DofElementData(p0_h, e->index);
		p_phi = DofElementData(phi, e->index);
		if (e->region_mark == 0 | e->region_mark == 2){
			p_phi[0] = PHI0_1 * (1. + CR_1 * (p_p[0] - p_p0[0]));
		}
		else if (e->region_mark == 1 | e->region_mark == 3){
			p_phi[0] = PHI0_2 * (1. + CR_2 * (p_p[0] - p_p0[0]));
		}
	}
}
void
update_dot_phi(DOF *p_h, DOF *dot_phi)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dotphi;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dotphi = DofElementData(dot_phi, e->index);
		NO = Interval_NO(p_p[0], 0);
		if (e->region_mark == 0 | e->region_mark == 2){
			p_dotphi[0] = PHI0_1 * CR_1;
		}
		else if (e->region_mark == 1 | e->region_mark == 3){
			p_dotphi[0] = PHI0_2 * CR_2;
		}
	}
}
