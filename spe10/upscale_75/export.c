#include "oil_water.h"
#include "phg.h"
void
export_vtkfile(const char *vtkfile, PHASE *oil, PHASE *water, MEDIUM *rock, int flag)
{
	sprintf(vtkfile, "oil_water%00003d.vtk", flag);
	phgExportVTK(oil->P->g, vtkfile, oil->P, oil->S, oil->U, oil->Source, oil->B, oil->DB, oil->Kr, oil->DKr, water->Source, water->B, water->DB, water->Kr, water->DKr, rock->phi, rock->perm, NULL);
}
void
export_test_data(PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, COM_TRO *control, WELL *well)
{
	int i = 0;
	phgPrintf("**************Test Read Data********************\n");
	phgPrintf("------------PHASE DATA:\n");
	phgPrintf("initial pressure  :  %lf, %lf\n", oil->PRESSURE0, water->PRESSURE0);
	phgPrintf("initial saturation:  %lf, %lf\n", oil->S0, water->S0);
	phgPrintf("viscosity:           %le, %le\n", oil->Mu, water->Mu);
	phgPrintf("rock compressibility:     %le\n", rock->C_R);
	phgPrintf("---------------WELL DATA:\n");
	for (i = 0; i < well->Number; i++){
		phgPrintf("well radius:     %lf, eff:   %lf\n", well->radius[i], well->eff_radius[i]);
	}
	for (i = 0; i < well->Number; i++){
		phgPrintf("well position:     %lf  %lf, region_mark:  %lf\n", well->position[2*i], well->position[2*i+1], well->region[i]);
	}
	phgPrintf("well   skin:              %lf\n", well->skin);
	phgPrintf("injection rate:           %lf\n", well->INJE);
	phgPrintf("max inje pwf:             %lf\n", well->MAX_INJE_BHP);
	phgPrintf("production pwf:           %lf\n", well->PROD_BHP);
	phgPrintf("--------MESH INFORMATION:\n");
	phgPrintf("NX, NY, NZ :   %d,  %d,  %d\n", rock->NX, rock->NY, rock->NZ);
	phgPrintf("DX, DY, DZ :   %lf,  %lf,  %lf\n", rock->DX, rock->DY, rock->DZ);
	phgPrintf("------------COMTROL DATA:\n");
	phgPrintf("max_ds     :              %lf\n", control->MAX_DS);
	phgPrintf("init_dt    :              %lf\n", control->init_dt);
	phgPrintf("max_dt     :              %lf\n", control->max_dt);
	phgPrintf("totoal time:              %lf\n", control->T);
	phgPrintf("system  err:              %le\n", control->TOL_sys);
	phgPrintf("conser  err:              %le\n", control->TOL_con);
	phgPrintf("iters   err:              %le\n", control->TOL_non);
	phgPrintf("***********End Test Read Data********************\n");

}
void
export_test_parameters(PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control, WELL *well)
{
	GRID *g = oil->P->g;
	SIMPLEX *e;
	FLOAT *so, *po, *muo, *bo, *dbo, *kro, *dkro, *sw, *pw, *muw, *bw, *dbw, *krw, *dkrw;
	phgPrintf("*********************Test Parameters******************\n");
	ForAllElements(g, e){
		so = DofElementData(oil->S, e->index);
		po = DofElementData(oil->P, e->index);
		bo = DofElementData(oil->B, e->index);
		dbo= DofElementData(oil->DB,e->index);
		kro = DofElementData(oil->Kr,e->index);
		dkro=DofElementData(oil->DKr,e->index);
		sw = DofElementData(water->S, e->index);
		pw = DofElementData(water->P, e->index);
		bw = DofElementData(water->B, e->index);
		dbw= DofElementData(water->DB,e->index);
		krw =DofElementData(water->Kr,e->index);
		dkrw=DofElementData(water->DKr,e->index);
		phgPrintf("oil saturation:    %lf\n", so[0]); 
		phgPrintf("oil 	 pressure:    %lf\n", po[0]);
		phgPrintf("oil  viscosity:    %le\n", oil->Mu);
		phgPrintf("oil         B :    %lf\n", bo[0]);
		phgPrintf("oil        DB :    %lf\n", dbo[0]);
		phgPrintf("oil        Kr :    %lf\n", kro[0]);
		phgPrintf("oil       DKr :    %lf\n", dkro[0]);
		phgPrintf("water saturation:    %lf\n", sw[0]);
		phgPrintf("water   pressure:    %lf\n", pw[0]);
		phgPrintf("water  viscosity:    %le\n", water->Mu);
		phgPrintf("water         B :    %lf\n", bw[0]);
		phgPrintf("water        DB :    %lf\n", dbw[0]);
		phgPrintf("water        Kr :    %lf\n", krw[0]);
		phgPrintf("rock        phi :    %lf\n", *(DofElementData(rock->phi, e->index)));
		phgPrintf("rock       Dphi :    %lf\n", *(DofElementData(rock->Dphi, e->index)));
		phgPrintf("rock       perm :    %e\n", rock->perm);
		phgPrintf("control 	T:         %lf\n", control->T);
		phgPrintf("control dt:         %lf\n", control->dt);
		phgPrintf("control TOL_non:         %e\n", control->TOL_non);
		phgPrintf("control TOL_sys:         %e\n", control->TOL_sys);
		phgPrintf("control TOL_con:         %e\n", control->TOL_con);
		phgPrintf("well production:         %e\n", well->PROD);
		phgPrintf("**************n-1 time step***************\n");
		phgPrintf("oil    B:      %lf\n", *DofElementData(oil_l->B, e->index));
		phgPrintf("oil    S:      %lf\n", *DofElementData(oil_l->S, e->index));
		phgPrintf("water  B:      %lf\n", *DofElementData(water_l->B, e->index));
		phgPrintf("water  S:      %lf\n", *DofElementData(water_l->S, e->index));
		phgPrintf("rock phi:      %lf\n", *DofElementData(rock_l->phi, e->index));
		exit(1);
	}
	phgPrintf("*****************End Of Test Parameters******************\n");
}
void
average_pressure(PHASE *oil, FLOAT *avp)
{
	GRID *g = oil->P->g;
	SIMPLEX *e;
	FLOAT *p_p, sum =0;
	ForAllElements(g, e){
		p_p = DofElementData(oil->P, e->index);
		sum += p_p[0];
	}
#if USE_MPI
	FLOAT sum0 = sum;
	MPI_Allreduce(&sum0, &sum, 1, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	*avp = sum / DofGetDataCountGlobal(oil->P) / PRESS_FACTOR;
}
void
Prod_result(PHASE *oil, PHASE *water, WELL *well, FLOAT *oilrate, FLOAT *waterrate)
{
	GRID *g = oil->Source->g;
	SIMPLEX *e;
	FLOAT oil_rate[well->Number-1] ;
	FLOAT water_rate[well->Number-1];
	bzero(oil_rate, (well->Number-1) * sizeof(FLOAT));
	bzero(water_rate, (well->Number-1) * sizeof(FLOAT));
	FLOAT *p_qo, *p_qw, vol = 0;
	int i = 0;
	ForAllElements(g, e){
		p_qo = DofElementData(oil->Source, e->index);
		p_qw = DofElementData(water->Source, e->index);
		vol = phgGeomGetVolume(g, e);
		for (i = 0; i < well->Number-1; i++){
			if (e->region_mark == well->region[i+1]){
				oil_rate[i] += p_qo[0] * vol;
				water_rate[i] += p_qw[0] * vol;
				break;
			}
		}
	}
#if USE_MPI
	FLOAT or[well->Number-1], wr[well->Number-1];
	memcpy(or, oil_rate, (well->Number-1) * sizeof(FLOAT));
	memcpy(wr, water_rate, (well->Number-1) * sizeof(FLOAT));
	MPI_Allreduce(or, oil_rate, well->Number-1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(wr, water_rate, well->Number-1, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	for (i = 0; i < well->Number-1; i++){
		oilrate[i] = -oil_rate[i] / VOLUME_FACTOR;
		waterrate[i] = -water_rate[i] / VOLUME_FACTOR;
	}
}
void
Bottom_Hole_Pressure(PHASE *oil, WELL *well, FLOAT *press)
{
	GRID *g = oil->P->g;
	SIMPLEX *e;
	FLOAT *p_p, vol =0;
	INT count[well->Number-1], i =0;
	FLOAT sum[well->Number-1];
	bzero(count, (well->Number-1) * sizeof(INT));
	bzero(sum, (well->Number-1) * sizeof(FLOAT));
	ForAllElements(g, e){
		if(e->region_mark == well->region[0]){
			vol += phgGeomGetVolume(g, e);
		}
	}
#if USE_MPI
	FLOAT vol0 = vol;
	MPI_Allreduce(&vol0, &vol, 1, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	ForAllElements(g, e){
		p_p = DofElementData(oil->P, e->index);
		for (i = 0; i < well->Number-1; i++){
			if (e->region_mark == well->region[i+1]){
				sum[i] += p_p[0];
				count[i] = count[i] + 1;
				break;
			}
		}
	}
#if USE_MPI
	FLOAT s[well->Number-1];
	INT c[well->Number-1];
	memcpy(s, sum, (well->Number-1) * sizeof(FLOAT));
	memcpy(c, count, (well->Number-1) * sizeof(INT));
	MPI_Allreduce(s, sum, well->Number-1, MPI_DOUBLE, MPI_SUM, g->comm);
	MPI_Allreduce(c, count, well->Number-1, MPI_INT, MPI_SUM, g->comm);
#endif
	for (i = 0; i < well->Number-1; i++){
		if (count[i] != 0)
			press[i] = sum[i] / count[i] / PRESS_FACTOR;
		else
			press[i] = sum[i] / PRESS_FACTOR;
	}
}
void
result_export(PHASE *oil, PHASE *water, WELL *well, COM_TRO *control, FLOAT *total_oil)
{
	FLOAT avp = 0;
	int i = 0;
	FLOAT *oil_rate = phgAlloc(4* sizeof(FLOAT));
	bzero(oil_rate, 4 * sizeof(FLOAT));
	FLOAT *cut = phgAlloc(4* sizeof(FLOAT));
	bzero(cut, 4 * sizeof(FLOAT));
	FLOAT *bhp = phgAlloc(4* sizeof(FLOAT));
	bzero(bhp, 4 * sizeof(FLOAT));
	average_pressure(oil, &avp);
	Prod_result(oil, water, well, oil_rate, cut);
	Bottom_Hole_Pressure(oil, well, bhp);
	*total_oil = (oil_rate[0] + oil_rate[1] + oil_rate[2] + oil_rate[3]);
	phgPrintf("---------Final Result For Time Step\n");
	phgPrintf("Producer Name-----Time(day)-----oil_rate(bbl/day)----water_rate(bbl/day)----BHP(Psi)\n");
	for (i = 0; i < well->Number-1; i++){
		phgPrintf("Producer%d:   %lf,    %lf,    %lf,    %lf\n", i+1, control->ct, oil_rate[i], cut[i], bhp[i]);
	}
	phgPrintf("Time:   %lf day,  average_press:    %lf Psi\n", control->ct, avp);
//	FLOAT inje_press = 0;
//	Inje_Well_BHP(well, oil->P, &inje_press);
//	phgPrintf("Time:   %lf day,  inje_press:    %lf Psi\n", control->ct, inje_press);
	phgFree(oil_rate);
	phgFree(cut);
	phgFree(bhp);
}
