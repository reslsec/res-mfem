#include "phg.h"
#include <string.h>
#include <math.h>
#include "black_oil.h"

static int RANGE = 4;

/*some parameters for permeability*/
static FLOAT SORW = 0.;
static FLOAT KRW_SORW = 0.6;
static FLOAT KRO_SWC = 1;
static FLOAT SWC = 0.;
static INT NO = 4;
static INT NW = 3;
static FLOAT KRO_SGC = 1.;
static  FLOAT KRG_SWC = 1;
static  FLOAT SLC = 0.;
static  FLOAT SGC = 0.;
static INT NG = 3;

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
		PVT[MUW_INDEX(i)] = PVT[MUW_INDEX(i)] / (1e6 * 3600.);
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
update_muw(DOF *p_h, DOF *mu_w)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_mu;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_mu = DofElementData(mu_w, e->index);
		NO = Interval_NO(p_p[0], 0);
		FLOAT temp; 
		temp = (p_p[0] - PVT[P_INDEX(NO)]) / (PVT[P_INDEX(NO + 1)] - PVT[P_INDEX(NO)]);
		p_mu[0] = 1. / (PVT[MUW_INDEX(NO)] + temp * (PVT[MUW_INDEX(NO + 1)] - PVT[MUW_INDEX(NO)]));
	}
}
void
update_bw(DOF *p_h, DOF *b_w)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_bw;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_bw = DofElementData(b_w, e->index);
		NO = Interval_NO(p_p[0], 0);
		FLOAT temp; 
		temp = (p_p[0] - PVT[P_INDEX(NO)]) / (PVT[P_INDEX(NO + 1)] - PVT[P_INDEX(NO)]);
//		p_bo[0] = 1.0 / PVT[BO_INDEX(NO)] + temp * (1.0 / PVT[BO_INDEX(NO + 1)] - 1.0 / PVT[BO_INDEX(NO)]);
		p_bw[0] = 1.0 / (PVT[BW_INDEX(NO)] + temp * ( PVT[BW_INDEX(NO + 1)] -  PVT[BW_INDEX(NO)]));
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
		p_dotphi[0] = PHI0 * CR;
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
update_dot_bw(DOF *p_h, DOF *dot_bw)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dotbw;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dotbw = DofElementData(dot_bw, e->index);
		NO = Interval_NO(p_p[0], 0);
		p_dotbw[0] = (1.0 / PVT[(BW_INDEX(NO + 1))] - 1.0 / PVT[(BW_INDEX(NO))])
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
update_dot_muo(DOF *p_h, DOF *dot)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dot;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dot = DofElementData(dot, e->index);
		NO = Interval_NO(p_p[0], 0);
		p_dot[0] = (1.0 / PVT[(MUO_INDEX(NO + 1))] - 1.0 / PVT[(MUO_INDEX(NO))])
			   		/ (PVT[(P_INDEX(NO + 1))] - PVT[(P_INDEX(NO))]);
	}
}	
void
update_dot_muw(DOF *p_h, DOF *dot)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dot;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dot = DofElementData(dot, e->index);
		NO = Interval_NO(p_p[0], 0);
		p_dot[0] = (1.0 / PVT[(MUW_INDEX(NO + 1))] - 1.0 / PVT[(MUW_INDEX(NO))])
			   		/ (PVT[(P_INDEX(NO + 1))] - PVT[(P_INDEX(NO))]);
	}
}	
void
update_dot_mug(DOF *p_h, DOF *dot)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_dot;
	int NO = 0;
	ForAllElements(g, e){
		p_p = DofElementData(p_h, e->index);
		p_dot = DofElementData(dot, e->index);
		NO = Interval_NO(p_p[0], 0);
		p_dot[0] = (1.0 / PVT[(MUG_INDEX(NO + 1))] - 1.0 / PVT[(MUG_INDEX(NO))])
			   		/ (PVT[(P_INDEX(NO + 1))] - PVT[(P_INDEX(NO))]);
	}
}	
/*need to be consisdered*/
#if 0
void
create_kro(DOF *s_o, DOF *kro)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_kro;
	ForAllElements(g, e){
		p_so = DofElementData(s_o, e->index);
		p_kro = DofElementData(kro, e->index);
		int no = 0, i;
		for (i = 0; i < 102; i++)
		{
			if (KRO[no] <= p_so[0]){
				no++;
				if (p_so[0] < KRO[no]){
					no--;
					break;
				}
			}
		}
		p_kro[0] = (KRO[no] + KRO[no + 1 + 1]) / 2;
	//	printf("no = %d,  kro[no] = %lf  kro[no+2] = %lf\n", no, KRO[no], KRO[(no+2)]);
//		printf("s_o = %lf, kro = %lf\n", p_so[0], p_kro[0]);
//		exit(1);
	}
}
void
create_krw(DOF *s_w, DOF *krw)
{
	GRID *g = s_w->g;
	SIMPLEX *e;
	FLOAT *p_sw, *p_krw;
	ForAllElements(g, e){
		p_sw = DofElementData(s_w, e->index);
		p_krw = DofElementData(krw, e->index);
		int nw = 0, i;
		for (i = 0; i < 102; i++)
		{
			if(KRW[nw] <= p_sw[0]){
				nw++;
				if(p_sw[0] < KRW[nw]){
					nw--;
					break;
				}
			}
		}
		p_krw[0] = (KRW[nw] + KRW[nw +2]) / 2;
	//	printf("no = %d,  krw[no] = %lf  krw[no+2] = %lf\n", nw, KRW[nw], KRW[(nw+2)]);
	//	printf("s_w = %lf, krw = %lf\n", p_sw[0],p_krw[0]);
	//	exit(1);
	}
}
void
create_krg(DOF *s_o, DOF *s_w, DOF *krg)
{
	GRID *g =s_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_sw, *p_krg;
	ForAllElements(g, e){
		p_so = DofElementData(s_o, e->index);
		p_sw = DofElementData(s_w, e->index);
		p_krg = DofElementData(krg, e->index);
		static int no = 0;
		FLOAT sg = 1. - p_so[0] - p_sw[0];
		while((sg >= KRG[0]) && (sg <= KRG[51 * 2]))
		{
			if (sg >= KRG[no]){
				no++;
				if (sg < KRG[no]){
					no--;
					break;
				}
			}
		}
		p_krg[0] = (KRG[no +1] + KRG[no + 1 + 1]) / 2;
	//	printf("no = %d,  krg[no] = %lf  krg[no+2] = %lf\n", no, KRG[no], KRG[(no+2)]);
		//printf("krg = %lf\n", p_krg[0]);
	//	exit(1);
	}
}
#endif
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
		p_phi[0] = PHI0 * (1. + CR * (p_p[0] - p_p0[0]));
	}
}
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
update_dot_krw(DOF *s_w, DOF *dsw_krw)
{
	GRID *g = s_w->g;
	SIMPLEX *e;
	FLOAT *p_sw, *p_dswkrw;
	ForAllElements(g, e){
		p_sw = DofElementData(s_w, e->index);
		p_dswkrw = DofElementData(dsw_krw, e->index);
		p_dswkrw[0] = NW * (1. / (1. - SWC - SORW)) * KRW_SORW * pow((p_sw[0] - SWC) / (1 - SWC - SORW), NW-1);
	}
}
void
create_krg(DOF *s_g, DOF *krg)
{
	GRID *g =s_g->g;
	SIMPLEX *e;
	FLOAT *p_sg, *p_krg;
	ForAllElements(g, e){
		p_sg = DofElementData(s_g, e->index);
		p_krg = DofElementData(krg, e->index);
		p_krg[0] = KRG_SWC * pow((p_sg[0] - SGC) / (1. - SLC - SGC), NG); 
	}
}
void
update_dot_krg(DOF *s_o, DOF *s_w, DOF *dso_krg, DOF *dsw_krg)
{
	GRID *g =s_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_sw, *p_dsokrg, *p_dswkrg;
	ForAllElements(g, e){
		p_so = DofElementData(s_o, e->index);
		p_sw = DofElementData(s_w, e->index);
		p_dsokrg = DofElementData(dso_krg, e->index);
		p_dswkrg = DofElementData(dsw_krg, e->index);
		p_dswkrg[0] = NG * (-1./ (1. - SLC- SGC)) * KRG_SWC * pow((1. - p_sw[0] - p_so[0] - SGC) / (1. - SLC - SGC), NG-1); 
		p_dsokrg[0] = NG * (-1./ (1. - SLC- SGC)) * KRG_SWC * pow((1. - p_sw[0] - p_so[0] - SGC) / (1. - SLC - SGC), NG-1); 
	}
}
void
create_kro(DOF *s_w, DOF *s_g, DOF *kro)
{
	GRID *g = kro->g;
	SIMPLEX *e;
	FLOAT Kro_ow, Kro_og, *p_kro, *p_so, *p_sw, *p_sg;
	ForAllElements(g, e){
		p_sw = DofElementData(s_w, e->index);
		p_sg = DofElementData(s_g, e->index);
		p_kro = DofElementData(kro, e->index);
		Kro_ow = KRO_SWC * pow((1.- p_sw[0] - SORW) / (1 - SWC - SORW), NO);
		Kro_og = KRO_SGC * pow((1.- p_sg[0] - SLC) / (1. - SLC - SGC), NO);
		p_kro[0] = Kro_ow * Kro_og;
	}
}
void
update_dot_kro(DOF *s_o, DOF *s_w, DOF *dso_kro, DOF *dsw_kro)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	FLOAT Kro_ow, Kro_og, *p_dsokro, *p_dswkro, *p_so, *p_sw;
	ForAllElements(g, e){
		p_sw = DofElementData(s_w, e->index);
		p_so = DofElementData(s_o, e->index);
		p_dsokro = DofElementData(dso_kro, e->index);
		p_dswkro = DofElementData(dsw_kro, e->index);
		Kro_ow = KRO_SWC * pow((1.- p_sw[0] - SORW) / (1 - SWC - SORW), NO);
		Kro_og = KRO_SGC * pow((p_so[0] + p_sw[0] - SLC) / (1. - SLC - SGC), NO);
		p_dsokro[0] = Kro_ow * NO * (1./ (1. -SLC -SGC)) * KRO_SGC * pow((p_so[0] + p_sw[0] - SLC) / (1. - SLC - SGC), NO-1);
		p_dswkro[0] = Kro_ow * NO * (1./ (1. -SLC -SGC)) * KRO_SGC * pow((p_so[0] + p_sw[0] - SLC) / (1. - SLC - SGC), NO-1) + Kro_og * NO * (-1. / (1 - SWC- SORW)) * KRO_SWC * pow((1.- p_sw[0] - SORW) / (1 - SWC - SORW), NO-1);
	}
}
