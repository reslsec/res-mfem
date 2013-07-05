/* Well */
/* use average Kr mu B p, etc. */

#define HAVE_WELL(e) (e->region_mark >= MESH_LEVEL)
#define WELL(e) ((e->region_mark - MESH_LEVEL) / MESH_LEVEL)
#define LEVEL(e) (e->region_mark % MESH_LEVEL)

void
Well_init(DOF *source)
{
	GRID *g = source->g;
	SIMPLEX *e;
	FLOAT *p_source, vol[WELL_NO];
	int well_no;

	bzero(vol, WELL_NO * sizeof(FLOAT));

	ForAllElements(g, e)
	{
#if 0
		for(well_no = 0; well_no < WELL_NO; well_no++){
			if(e->region_mark == well_no + 1){
				vol[well_no] += phgGeomGetVolume(g, e);
			}
		}
#else
		if(HAVE_WELL(e)){
			well_no = WELL(e);
			vol[well_no] += phgGeomGetVolume(g, e);
		}
#endif
	}

#if USE_MPI
	FLOAT vol0[WELL_NO];
	memcpy(vol0, vol, WELL_NO * sizeof(FLOAT));
	MPI_Allreduce(vol0, vol, WELL_NO, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	phgPrintf("well volume: %lf\n", vol[0]);

	FLOAT pro[WELL_NO];
	for(well_no = 0; well_no < WELL_NO; well_no++){
		pro[well_no] = PRODUCTION[well_no] / vol[well_no];
	}

	ForAllElements(g, e)
	{
#if 0
		for(well_no = 0; well_no < WELL_NO; well_no++){
			if(e->region_mark == well_no + 1)
			{
				p_source = DofElementData(source, e->index);
				p_source[0] = - pro[well_no]; // Assume 
			}
		}
#else
		if(HAVE_WELL(e)){
			well_no = WELL(e);
			p_source = DofElementData(source, e->index);
			p_source[0] = - pro[well_no];
		}
#endif
	}
}

void
Well_Pressure(PHASE *oil, FLOAT *pressure)
{
	GRID *g = oil->g;
	SIMPLEX *e;
	FLOAT *p_p, sum_vol[WELL_NO], sum_lambda[WELL_NO];
	INT i, well_no, nbas = oil->pressure->type->nbas;

	bzero(sum_vol, WELL_NO * sizeof(FLOAT));
	bzero(sum_lambda, WELL_NO * sizeof(FLOAT));
	bzero(pressure, WELL_NO * sizeof(FLOAT));

	ForAllElements(g, e)
	{
#if 0
		for(well_no = 0; well_no < WELL_NO; well_no++){
			if(e->region_mark == well_no + 1){
				FLOAT temp;
				phgQuadDofTimesBas(e, oil->pressure, oil->source, 0, QUAD_DEFAULT, &temp);
				pressure[well_no] += temp;
				sum_vol[well_no] += phgGeomGetVolume(g, e);
				sum_lambda[well_no] += phgQuadDofDotDof(e, oil->Kr, oil->B, QUAD_DEFAULT);
			}
		}
#else
		if(HAVE_WELL(e)){
			well_no = WELL(e);
			FLOAT temp;
			phgQuadDofTimesBas(e, oil->pressure, oil->source, 0, QUAD_DEFAULT, &temp);
			pressure[well_no] += temp;
			sum_vol[well_no] += phgGeomGetVolume(g, e);
			sum_lambda[well_no] += phgQuadDofDotDof(e, oil->Kr, oil->B, QUAD_DEFAULT);
		}
#endif
	}

#if USE_MPI
	FLOAT sum_vol0[WELL_NO], pressure0[WELL_NO], sum_lambda0[WELL_NO];
	memcpy(sum_vol0, sum_vol, WELL_NO * sizeof(FLOAT));
	memcpy(sum_lambda0, sum_lambda, WELL_NO * sizeof(FLOAT));
	memcpy(pressure0, pressure, WELL_NO * sizeof(FLOAT));
	MPI_Allreduce(sum_vol0, sum_vol, WELL_NO, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(sum_lambda0, sum_lambda, WELL_NO, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(pressure0, pressure, WELL_NO, MPI_DOUBLE, MPI_SUM, g->comm);
#endif

	for(well_no = 0; well_no < WELL_NO; well_no++){
		pressure[well_no] = pressure[well_no] / sum_vol[well_no];
		sum_lambda[well_no] /= sum_vol[well_no];

		pressure[well_no] = pressure[well_no] - oil->MU * PRODUCTION[well_no] * 
			(log(8.0 / WELL_RADIUS[well_no]) + 3.0 - 0.5) 
			/ (2.0 * M_PI * AB_K[0] * MESH_HEIGHT * sum_lambda[well_no]); 
		FLOAT temp = 
		 - oil->MU * PRODUCTION[well_no] * 
			(log(8.0 / WELL_RADIUS[well_no]) + 3.0 - 0.5) 
			/ (2.0 * M_PI * AB_K[0] * MESH_HEIGHT * sum_lambda[well_no]); 
		phgPrintf("correct: %lf\n", temp);
	}
}

void
Oil_Source2Water_Source(PHASE *oil, PHASE *water)
{
	GRID *g = oil->g;
	SIMPLEX *e;
	FLOAT *p_source_oil, *p_source_water;
	FLOAT lambda_oil, lambda_water;
	int well_no;

	ForAllElements(g, e)
	{
		if(HAVE_WELL(e)){
			p_source_oil = DofElementData(oil->source, e->index);
			p_source_water = DofElementData(water->source, e->index);
			lambda_oil = phgQuadDofDotDof(e, oil->Kr, oil->B, QUAD_DEFAULT) / oil->MU;
			lambda_water = phgQuadDofDotDof(e, water->Kr, water->B, QUAD_DEFAULT) / water->MU;
		//printf("oil = %lf, water = %lf\n", lambda_oil, lambda_water);
			p_source_water[0] = lambda_water * p_source_oil[0] / lambda_oil;
		}
	}
}

void
WOR(PHASE *oil, PHASE *water, FLOAT *OP,FLOAT *WP)
{
	GRID *g = oil->g;
	SIMPLEX *e;
	FLOAT *p_source_oil, *p_source_water;
	FLOAT OOP[WELL_NO], WWP[WELL_NO];
	int well_no;

	bzero(OOP, WELL_NO * sizeof(FLOAT));
	bzero(WWP, WELL_NO * sizeof(FLOAT));

	ForAllElements(g, e)
	{
		if(HAVE_WELL(e)){
			well_no = WELL(e);
			p_source_oil = DofElementData(oil->source, e->index);
			p_source_water = DofElementData(water->source, e->index);
			OOP[well_no] += p_source_oil[0] * phgGeomGetVolume(g, e);
			WWP[well_no] += p_source_water[0] * phgGeomGetVolume(g, e);
		}
	}
#if USE_MPI
	FLOAT OOP0[WELL_NO], WWP0[WELL_NO];
	memcpy(OOP0, OOP, WELL_NO * sizeof(FLOAT));
	memcpy(WWP0, WWP, WELL_NO * sizeof(FLOAT));
	MPI_Allreduce(OOP0, OOP, WELL_NO, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(WWP0, WWP, WELL_NO, MPI_DOUBLE, MPI_SUM, g->comm);
#endif

	for(well_no = 0; well_no < WELL_NO; well_no++){
		OP[well_no] += OOP[well_no];
		WP[well_no] += WWP[well_no];
	}
}
