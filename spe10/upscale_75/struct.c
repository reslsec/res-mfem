#include "oil_water.h"
#include "phg.h"
PHASE *phgGetPhase(GRID *g)
{
	static PHASE *_phase = NULL;
	_phase = phgAlloc(sizeof(PHASE));
	if(_phase != NULL){
		_phase->Kr = phgDofNew(g, DOF_P0, 1, "Kr", DofNoAction);
		_phase->DKr = phgDofNew(g, DOF_P0, 1, "Dot_Kr", DofNoAction);
		_phase->B = phgDofNew(g, DOF_P0, 1, "B", DofNoAction);
		_phase->DB = phgDofNew(g, DOF_P0, 1, "DB", DofNoAction);
		_phase->S = phgDofNew(g, DOF_P0, 1, "saturantion", DofNoAction);
		_phase->P = phgDofNew(g, DOF_P0, 1, "phase pressure", DofNoAction);
		_phase->Source = phgDofNew(g, DOF_P0, 1, "phase source", DofNoAction);
		_phase->U = phgDofNew(g, DOF_RT1, 1, "phase velocity", DofNoAction);
		_phase->P0 = phgDofNew(g, DOF_P0, 1, "Initial pressure", DofNoAction);
		_phase->Mu = 0.;//initial value
		_phase->S0 = 0.;
		_phase->PRESSURE0 = 20.;
		return _phase;
	}
	else
		phgError(1, "Can't create %s\n", "phase");
}

MEDIUM *phgGetMedium(GRID *g)
{
	static MEDIUM *_medium = NULL;
	_medium = phgAlloc(sizeof(MEDIUM));
	if (_medium != NULL){
	//	_medium->perm = 1e-13; 
		_medium->perm = phgDofNew(g, DOF_P0, 3, "perm", DofNoAction);
		_medium->phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
		_medium->phi0 = phgDofNew(g, DOF_P0, 1, "phi0", DofNoAction);
		//_medium->phi0 = 0.2; 
		_medium->Dphi = phgDofNew(g, DOF_P0, 1, "Dphi", DofNoAction);
		_medium->C_R = 0.15e-3;//initial value	
		_medium->NX = 3;
		_medium->NY = 5;
		_medium->NZ = 5;
		_medium->DX = 121.9200;
		_medium->DY = 134.1120;
		_medium->DZ = 10.3632;
		return _medium;
	}
	else
		phgError(1, "Can't create %s\n", "medium");
}

COM_TRO *phgGetComTro(GRID *g)
{
	static COM_TRO *_control = NULL;
	_control = phgAlloc(sizeof(COM_TRO));
	if (_control != NULL){
		_control->T = 30.;//initial value
		_control->init_dt = 1e-4;
		_control->dt = 1e-4;
		_control->ct = 0;
		_control->pt = 0;
		_control->max_dt = 1e+1;
		_control->TOL_sys = 1e-6;
		_control->TOL_non = 1e-4;
		_control->TOL_con = 1e-1;
		_control->MAX_DS = 0.05;
		_control->MAX_DP = 300 * 0.00689;
		return _control;
	}
	else
		phgError(1, "Can't create %s\n", "control");
}

WELL *phgGetWellParameters(GRID *g)
{
	static WELL *_well = NULL;
	_well = phgAlloc(sizeof(WELL));
	if (_well != NULL){
		_well->PROD = 15.9;//initial value
		_well->INJE = 0.;
		_well->MAX_INJE_BHP = 0.;
		_well->PROD_BHP = 0.;
//		_well->radius = 0.;
//		_well->eff_radius = 0.;
//		_well->skin = 3.;
		return _well;
	}
	else
		phgError(1, "Can't create %s\n", "well");

}
void phgFreePhase(PHASE *_phase)
{
	if(!_phase)
		return;
	phgDofFree(&(_phase->B));
	phgDofFree(&(_phase->DB));
	phgDofFree(&(_phase->Kr));
	phgDofFree(&(_phase->DKr));
	phgDofFree(&(_phase->Source));
	phgDofFree(&(_phase->P0));
	phgDofFree(&(_phase->P));
	phgDofFree(&(_phase->S));
	phgDofFree(&(_phase->U));
	phgFree(_phase);
}
void phgFreeMedium(MEDIUM *_medium)
{
	if(!_medium)
		return;
	phgDofFree(&(_medium->phi));
	phgDofFree(&(_medium->phi0));
	phgDofFree(&(_medium->Dphi));
	phgFree(_medium);
}
void phgFreeComTro(MEDIUM *_control)
{
	if(!_control)
		return;
	phgFree(_control);
}
void phgFreeWell(WELL *_well)
{
	if(!_well)
		return;
	phgFree(_well);
}
