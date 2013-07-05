#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "phg.h"
//#include "ss_oil_water.h"
FLOAT PRODUCTION = 20. / 24.;
void
Well_init(DOF *source)
{
	GRID *g = source->g;
	SIMPLEX *e;
	FLOAT *p_source;
	FLOAT vol_1 = 0, vol_2 =0;
	ForAllElements(g, e){
		if(e->region_mark == 2){
			vol_1 += phgGeomGetVolume(g, e);
		}
		else if(e->region_mark == 3)
		{
			vol_2 += phgGeomGetVolume(g, e);
		}
	}
#if USE_MPI
	FLOAT vol_10 = vol_1;
	FLOAT vol_20 = vol_2;
//	memcpy(vol0, vol, sizeof(FLOAT));
	MPI_Allreduce(&vol_10, &vol_1, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&vol_20, &vol_2, 1, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	ForAllElements(g, e){
		if(e->region_mark == 2){
			p_source = DofElementData(source, e->index);
			p_source[0] = - PRODUCTION / (vol_1+vol_2);
		}
		else if(e->region_mark == 3){
			p_source = DofElementData(source, e->index);
			p_source[0] = - PRODUCTION / (vol_2+vol_1);
		}
	}
}
void
Well_Pressure(DOF *p_h, FLOAT *pressure)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, sum = 0;
	INT count = 0, i, well_no, nbas = p_h->type->nbas;

	ForAllElements(g, e)
	{
		*pressure = 0.0;

		if(e->region_mark == 2){
			p_p = DofElementData(p_h, e->index);
			sum += p_p[0];
			count++;
		}
		else if (e->region_mark == 3){
			p_p = DofElementData(p_h, e->index);
			sum += p_p[0];
			count++;
		}
	}

#if USE_MPI
	FLOAT sum0 = sum;
	INT count0 = count;
//	memcpy(sum0, sum, sizeof(FLOAT));
//	memcpy(count0, count, sizeof(INT));
	MPI_Allreduce(&sum0, &sum, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&count0, &count, 1, MPI_INT, MPI_SUM, g->comm);
#endif
	*pressure = sum / count;
}
void 
Peaceman_Model(DOF *p_h, DOF *b_o, DOF *kro, DOF *mu_o, FLOAT *pwf)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_bo, *p_kro, *p_muo;
	int count = 0;
	FLOAT sum_lambda = 0, sum_vol = 0, sum_pwf = 0;
	FLOAT sum_lambda1 = 0, sum_vol1 = 0, sum_pwf1 = 0;
	FLOAT sum_lambda2 = 0, sum_vol2 = 0, sum_pwf2 = 0;
	INT i;
	ForAllElements(g, e){
			p_p = DofElementData(p_h, e->index);
			p_kro = DofElementData(kro, e->index);
			p_bo = DofElementData(b_o, e->index);
			p_muo = DofElementData(mu_o, e->index);
		if (e->region_mark == 2){
			FLOAT K = 1e-13;
			sum_pwf += p_p[0];
			sum_lambda+=  K * p_kro[0] * p_muo[0] * p_bo[0];
			count++;
		}
		else if (e->region_mark == 3){
			FLOAT K = 2e-13;
			sum_pwf += p_p[0];
			sum_lambda+=  K * p_kro[0] * p_muo[0] * p_bo[0];
			count++;
		}
	}
#if USE_MPI
	FLOAT sum_pwf0 = sum_pwf, sum_lambda0 = sum_lambda;
	INT count0 = count;
	MPI_Allreduce(&sum_pwf0, &sum_pwf, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&sum_lambda0, &sum_lambda, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&count0, &count, 1, MPI_INT, MPI_SUM, g->comm);
#endif
	*pwf = sum_pwf / count;
	sum_lambda = sum_lambda / count;
	*pwf = *pwf - PRODUCTION * (log(10. / 0.1) + 3.0 - 0.5) / (2. * M_PI * 10. * sum_lambda);
}
void 
PeacemanModel(DOF *p_h, DOF *b_o, DOF *kro, DOF *mu_o, FLOAT *pwf)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT sum_lambda = 0, sum_vol = 0, sum_pwf = 0;
	FLOAT tmp = 0;
	FLOAT *p_p, *p_bo, *p_kro, *p_muo;
	ForAllElements(g, e){
		if (e->region_mark == 2){
			p_p = DofElementData(p_h, e->index);
			p_kro = DofElementData(kro, e->index);
			p_bo = DofElementData(b_o, e->index);
			p_muo = DofElementData(mu_o, e->index);
			phgQuadDofTimesBas(e, p_h, b_o, 0, QUAD_DEFAULT, &tmp);
			sum_pwf += tmp;
			sum_vol += phgGeomGetVolume(g, e);
			FLOAT K = 1e-13;
			sum_lambda +=  K * p_bo[0] * p_kro[0] * p_muo[0] * phgGeomGetVolume(g, e);
		}
		else if (e->region_mark == 3){
			p_p = DofElementData(p_h, e->index);
			p_kro = DofElementData(kro, e->index);
			p_bo = DofElementData(b_o, e->index);
			p_muo = DofElementData(mu_o, e->index);
			phgQuadDofTimesBas(e, p_h, b_o, 0, QUAD_DEFAULT, &tmp);
			sum_pwf += tmp;
			sum_vol += phgGeomGetVolume(g, e);
			FLOAT K = 2e-13;
			sum_lambda +=  K * p_bo[0] * p_kro[0] * p_muo[0] * phgGeomGetVolume(g, e);
		}
	}
#if USE_MPI
	FLOAT sum_vol0 = 0, pwf0 = 0, sum_lambda0 = 0;
	sum_vol0 = sum_vol, sum_lambda0 = sum_lambda, pwf0 = sum_pwf;
	MPI_Allreduce(&sum_vol0, &sum_vol, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&sum_lambda0, &sum_lambda, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&pwf0, &sum_pwf, 1, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	*pwf = (sum_pwf) / (sum_vol);
	sum_lambda = (sum_lambda) / (sum_vol);
	*pwf = *pwf - PRODUCTION * (log(10. / 0.1) + 3.0 - 0.5) / (2. * M_PI * 10. * sum_lambda);
}
