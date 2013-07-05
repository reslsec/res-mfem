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
	FLOAT vol = 0;
	ForAllElements(g, e){
		if(e->region_mark == 1){
			vol += phgGeomGetVolume(g, e);
		}
	}
#if USE_MPI
	FLOAT vol0 = vol;
//	memcpy(vol0, vol, sizeof(FLOAT));
	MPI_Allreduce(&vol0, &vol, 1, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	ForAllElements(g, e){
		if(e->region_mark == 1){
			p_source = DofElementData(source, e->index);
			p_source[0] = - PRODUCTION / vol;
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

		if(e->region_mark == 1){
			p_p = DofElementData(p_h, e->index);
			for(i = 0; i < nbas; i++){
				*pressure += p_p[i];
			}
			*pressure = *pressure / nbas;
			sum += *pressure;
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
	FLOAT sum_lambda = 0, sum_vol = 0, sum_pwf = 0;
	FLOAT K = 1e-13;
	INT i;
	ForAllElements(g, e){
		if (e->region_mark == 1){
			FLOAT tmp = 0;
			phgQuadDofTimesBas(e, p_h, b_o, 0, QUAD_DEFAULT, &tmp);
			sum_pwf += tmp;
			sum_vol += phgGeomGetVolume(g, e);
			sum_lambda += phgQuadDofDofDof(e, kro, b_o, mu_o, QUAD_DEFAULT);
		}
	}
#if USE_MPI
	FLOAT sum_vol0 = 0, pwf0 = 0, sum_lambda0 = 0;
	sum_vol0 = sum_vol, sum_lambda0 = sum_lambda, pwf0 = sum_pwf;
	MPI_Allreduce(&sum_vol0, &sum_vol, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&sum_lambda0, &sum_lambda, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&pwf0, &sum_pwf, 1, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	*pwf = sum_pwf / sum_vol;
	sum_lambda = sum_lambda / sum_vol;
	*pwf = *pwf - PRODUCTION * (log(10. / 0.1) + 3.0 - 0.5) / (2. * M_PI * K * 10. * sum_lambda);
}
