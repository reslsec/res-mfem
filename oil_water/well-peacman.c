#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "phg.h"
//#include "ss_oil_water.h"
FLOAT PRODUCTION = 10. / 24.;
static FLOAT MESH_HEIGHT = 10.;
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
Well_Pressure(DOF *p_h, DOF *b_o, FLOAT *pressure, FLOAT *pwf)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, sum_vol = 0, sum_lambda = 0, sum_pressure = 0, *p_bo;
	INT count = 0;

	FLOAT	mu = 1e-9 / 3600, K = 1e-13, rw = 0.1, re = 10.; 
	ForAllElements(g, e){
		if(e->region_mark == 1){
			p_p = DofElementData(p_h, e->index);
			p_bo = DofElementData(b_o, e->index);
			sum_pressure += p_p[0] * phgGeomGetVolume(g, e);
			sum_lambda += p_bo[0] * phgGeomGetVolume(g, e); 
			sum_vol += phgGeomGetVolume(g, e);
		}
	}

#if USE_MPI
	FLOAT sum_vol0 = sum_vol, sum_lambda0 = sum_lambda, sum_pressure0 = sum_pressure;
	MPI_Allreduce(&sum_vol0, &sum_vol, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&sum_lambda0, &sum_lambda, 1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(&sum_pressure0, &sum_pressure, 1, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	*pressure = sum_pressure / sum_vol;
	sum_lambda = sum_lambda / sum_vol;
    FLOAT tmp  = (-PRODUCTION * mu * sum_lambda * (log(re/rw) + 3 - 0.5)) / (2 * M_PI * K * MESH_HEIGHT); 
    *pwf = *pressure + tmp;
}
