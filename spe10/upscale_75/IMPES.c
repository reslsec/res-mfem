 /* This sample code solves the oil_water_phase equation:         *
 * ****************************************************************
 * problem: p_{t} -\Delta{p} = func_f(x, t) \in \Omega X [0, T]   *
 *          p = func_g (x, t) \in \partial\Omega X [0, T]         * 
 *          p = func_p0       \in \Omega t = 0.0                  *
 *WARNING! The unit is h, m, mpa!!                                *
 revision 1.1
 date: 2011/04/20 07:41:10;  author: liuming;  state: Exp;
 ******************************************************************/
#include "phg.h"
#include "oil_water.h"
#include <math.h>
#define rw 0.1
#define re 10
#define prod 18
/* build linear system */
static double
use_time(GRID *g, BOOLEAN flag, double mflops)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;
 
    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    mem = phgMemoryUsage(g, NULL);

    if (flag) {
	if (mflops <= 0)
	    phgPrintf("[%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
	else
	    phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n", mem / (1024.0 * 1024.0),
			 et, mflops*1e-3);
    }
    return et;
}
static void
build_mat_vec(DOF *u_o, DOF *p_h, DOF *s_o, PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g)
{
	GRID *g = u_o->g;
	SIMPLEX *e;
	FLOAT *p_bo, *p_bw, *p_bol, *p_bwl, *p_kro, *p_krw, *p_sol, *p_so, *p_phi, *p_ph, *p_dbw, *p_dbo, *p_dphi;
	INT N = u_o->type->nbas * u_o->dim;
	INT M = p_h->type->nbas * p_h->dim;
	INT I[N], J[M];
	FLOAT mat_A[N][N], mat_TB[N][M], mat_B[M][N], mat_C[M][M], rhs_f[N], rhs_g[M];
	phgVecDisassemble(vec_f);
	phgVecDisassemble(vec_g);
	int i, j, k;
	ForAllElements(g, e){
		p_phi = DofElementData(rock->phi, e->index);
		p_dphi = DofElementData(rock->Dphi, e->index);
		p_bo = DofElementData(oil->B, e->index);
		p_bol = DofElementData(oil_l->B, e->index);
		p_dbo = DofElementData(oil->DB, e->index);
		p_bw = DofElementData(water->B, e->index);
		p_bwl = DofElementData(water_l->B, e->index);
		p_dbw = DofElementData(water->DB, e->index);
		p_so = DofElementData(oil->S, e->index);
		p_sol = DofElementData(oil_l->S, e->index);
		p_kro  = DofElementData(oil->Kr, e->index);
		p_krw  = DofElementData(water->Kr, e->index);
		FLOAT T_o = 0, T_w = 0;
		T_w = p_krw[0] * p_bw[0] / water->Mu;
		T_o = p_kro[0] * p_bo[0] / oil->Mu;
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				mat_A[i][j] = control->dt * oil->Mu * phgQuadBasDotBas(e, u_o, i, u_o, j, QUAD_DEFAULT) / (p_kro[0] * p_bo[0] * rock->perm);
			}
		}
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				mat_TB[i][j] = -control->dt * phgQuadDivBasDotBas(e, u_o, i, p_h, j, QUAD_DEFAULT);
			}
		}
		FLOAT wat_sw = 0, oil_sw = 0, wat_cp = 0, oil_cp = 0, beta = 0;
		for (i = 0; i < M; i++){
			for (j = 0; j < M; j++){
				wat_sw = p_phi[0] * p_bw[0] * phgQuadBasDotBas(e, s_o, i, s_o, j, QUAD_DEFAULT);
				oil_sw = p_phi[0] * p_bo[0] * phgQuadBasDotBas(e, s_o, i, s_o, j, QUAD_DEFAULT);
			}	
		}
		wat_cp = (1. - p_so[0]) * (p_bw[0]*p_dphi[0] + p_phi[0] * p_dbw[0]);
		oil_cp = (p_so[0]) * (p_bo[0] * p_dphi[0] + p_phi[0] * p_dbo[0]);
		beta = 1. / oil_sw + T_w / (T_o * wat_sw);
		FLOAT quad = 0.;
		for (i = 0; i < M; i++){
			for (j = 0; j < M; j++){
			 	quad =  phgQuadBasDotBas(e, p_h, i, p_h, j, QUAD_DEFAULT);
				mat_C[i][j] = (wat_cp / wat_sw + oil_cp / oil_sw) * quad / beta;
			}
		}
		/*create oil rhs*/
		FLOAT quad_phi = 0., quad_phil = 0., quad_qo = 0, quad_qw = 0, quad_pnk=0;
		FLOAT rhs_oil = 0., rhs_wat = 0.;
		for (i = 0; i < M; i++){
			phgQuadDofTimesBas(e, oil->Source, p_h, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, water->Source, p_h, i, QUAD_DEFAULT, &quad_qw);
			phgQuadDofTimesBas(e, rock->phi, p_h, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, p_h, p_h, i, QUAD_DEFAULT, &quad_pnk);
			phgQuadDofTimesBas(e, rock_l->phi, p_h, i, QUAD_DEFAULT, &quad_phil);
			rhs_wat = -control->dt * quad_qw - wat_cp * quad_pnk -  p_bwl[0] *  (1. - p_sol[0]) * quad_phil + p_bw[0] * quad_phi;
			rhs_oil = -control->dt * quad_qo -oil_cp * quad_pnk - p_bol[0] * p_sol[0] * quad_phil;
			rhs_g[i] = (rhs_wat / wat_sw + rhs_oil / oil_sw) / beta;
		}
		for (i = 0; i < N; i++){
			rhs_f[i] = 0.;
		}
		/* Handle Bdry Conditions*/
		for (i = 0; i < N; i++){
			if (phgDofGetElementBoundaryType(u_o, e, i) & (NEUMANN | DIRICHLET)){
				bzero(mat_A[i], N * sizeof(mat_A[i][0]));
				bzero(mat_TB[i], M * sizeof(mat_TB[i][0]));
				for (j = 0; j < N; j++){
					mat_A[j][i] = 0.;
				}
				mat_A[i][i] = 1.;
			}
		}

		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				mat_B[j][i] = mat_TB[i][j];
			}
		}
		for (i = 0; i < N; i++){
			I[i] = phgMapE2L(map_u, 0, e, i);
		}
		for (i = 0; i < M; i++){
			J[i] = phgMapE2L(map_p, 0, e, i);
		}
		phgMatAddEntries(A, N, I, N, I, mat_A[0]);
		phgMatAddEntries(TB, N, I, M, J, mat_TB[0]);
		phgMatAddEntries(B, M, J, N, I, mat_B[0]);
		phgMatAddEntries(C, M, J, M, J, mat_C[0]);
		phgVecAddEntries(vec_f, 0, N, I, rhs_f);
		phgVecAddEntries(vec_g, 0, M, J, rhs_g);
	}
 	phgMatAssemble(A);
  	phgMatAssemble(B);
	phgMatAssemble(TB);			   
	phgMatAssemble(C);
   	phgVecAssemble(vec_f);
	phgVecAssemble(vec_g);
}
static void
Solve_Pressure(DOF *u_o, DOF *p_h, DOF *s_o, PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control)
{
		GRID *g = u_o->g;	
	    MAT *A, *B, *TB, *C;
	    VEC *vec_f, *vec_g;
	    MAP *map_u, *map_p;
		/*Create MAP for Mat and Vec*/
	    map_p = phgMapCreate(p_h, NULL);
		map_u = phgMapCreate(u_o, NULL);
		A = phgMapCreateMat(map_u, map_u);
		B = phgMapCreateMat(map_p, map_u);
		TB = phgMapCreateMat(map_u, map_p);
		C = phgMapCreateMat(map_p, map_p);
		vec_f = phgMapCreateVec(map_u, 1);
		vec_g = phgMapCreateVec(map_p, 1);
			
		phgPrintf("build  p_h:           ");
		use_time(g, FALSE, 0.);
		build_mat_vec(u_o, p_h, s_o, oil, water, oil_l, water_l, rock, rock_l, control, map_u, map_p, A, B, TB,C, vec_f, vec_g);
		use_time(g, TRUE, 0.);
#if USE_UZAWA
		int nits_amg = 0, nits_uzawa = 0;
		use_time(g, FALSE, 0.);
		phgPrintf("Assemble H:             ");
		DOF *B_data = DivRT(p_h, u_o, control);
        MAT *H = BTAB(A, B_data, p_h, u_o);
      	phgDofFree(&B_data);
		use_time(g, TRUE, 0.);
    	/* new implementation */
		phgPrintf("solve p_h:           ");
		use_time(g, FALSE, 0.);
		nits_uzawa = uzawapcg(H, u_o, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
		use_time(g, TRUE, 0.);
		phgMatDestroy(&H);
		phgPrintf("Max iter of AMG---------%d\n", nits_amg);
		phgPrintf("Nits: uzawa:----------------%d\n", nits_uzawa);
#endif
#if USE_BLOCK
		SOLVER *solver;
		MAT *mat[4], *solver_mat;
		FLOAT coef[4];

		mat[0] = A;  mat[1] = TB;
		coef[0] = 1.; coef[1] = 1.;
		mat[2] = B;  mat[3] = C;
		coef[2] = 1.; coef[3] = -1.;
		solver_mat= phgMatCreateBlockMatrix(g, 2, 2, mat, coef, NULL);
		solver = phgSolverCreate(SOLVER_MUMPS, u_o, p_h, NULL);
		phgMatDestroy(&(solver->mat));
		solver->mat=solver_mat;
		solver_mat->refcount++;
		solver->mat->handle_bdry_eqns = FALSE;
		solver->rhs->mat = solver->mat;
   		INT N = vec_f->map->nlocal;
   		INT M = vec_g->map->nlocal;
		memcpy(solver->rhs->data, vec_f->data, sizeof(*vec_f->data) * N);
		memcpy(solver->rhs->data+N, vec_g->data, sizeof(*vec_g->data) * M);
		use_time(g, FALSE, 0.);
		phgPrintf("solve linear system:           ");
		phgSolverSolve(solver, FALSE, u_o, p_h, NULL);
		use_time(g, TRUE, 0.);
		phgSolverDestroy(&solver);
		phgMatDestroy(&(solver_mat));
#endif
		phgMatDestroy(&A);
		phgMatDestroy(&B);
		phgMatDestroy(&TB);
		phgMatDestroy(&C);
		phgMapDestroy(&map_p);
		phgMapDestroy(&map_u);
		phgVecDestroy(&vec_g);
		phgVecDestroy(&vec_f);
}
static void
build_oil_so(SOLVER *solver, DOF *s_o, DOF *div_uo, DOF *p_h, DOF *p_nk, PHASE *oil, PHASE *oil_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	int i,j;
	int N = s_o->type->nbas;
	FLOAT A[N][N], rhs[N];
	INT I[N];
	FLOAT *p_phi, *p_bo, *p_bol, *p_so, *p_sol, *p_phil, *p_ph, *p_dphi, *p_db;
	ForAllElements(g, e){
		p_phi = DofElementData(rock->phi, e->index);
		p_dphi = DofElementData(rock->Dphi, e->index);
		p_bo  = DofElementData(oil->B, e->index);
		p_bol = DofElementData(oil_l->B, e->index);
		p_db = DofElementData(oil->DB, e->index);
		p_so  = DofElementData(oil->S, e->index);
		p_sol  = DofElementData(oil_l->S, e->index);
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for (j = 0; j < N; j++){
				A[i][j] = p_phi[0] * p_bo[0] * phgQuadBasDotBas(e, s_o, i, s_o, j, QUAD_DEFAULT); 
			}
		}
		FLOAT quad_qo = 0, quad_ph = 0, quad_pnk = 0, quad_phil = 0, quad_phi = 0.,quad_divuo = 0;
		FLOAT oil_cp = 0;
		oil_cp = p_so[0] * (p_bo[0] * p_dphi[0] + p_phi[0] * p_db[0]);
		for (i = 0; i < N; i++){
			phgQuadDofTimesBas(e, oil->Source, s_o, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, p_h, s_o, i, QUAD_DEFAULT, &quad_ph);
			phgQuadDofTimesBas(e, p_nk, s_o, i, QUAD_DEFAULT, &quad_pnk);
			phgQuadDofTimesBas(e, rock_l->phi, s_o, i, QUAD_DEFAULT, &quad_phil);
			phgQuadDofTimesBas(e, div_uo, s_o, i, QUAD_DEFAULT, &quad_divuo);
			rhs[i] = control->dt * quad_qo - control->dt * quad_divuo - oil_cp * (quad_ph-quad_pnk) + p_sol[0] * p_bol[0] * quad_phil;
		}
		for (i = 0; i < N; i++){
			phgSolverAddMatrixEntries(solver, 1, I+i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
}
static void
Solver_Oileqn_So(DOF *s_o, DOF *p_h, DOF *u_o, DOF *p_nk, PHASE *oil, PHASE *oil_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control)
{
	GRID *g = s_o->g;
	SOLVER *solver;
	DOF *div_uo;
	div_uo = phgDofDivergence(u_o, NULL, NULL, NULL);
	solver = phgSolverCreate(SOLVER_PCG, s_o, NULL);
	phgPrintf("build s_o:           ");
	use_time(g, FALSE, 0.);
	build_oil_so(solver, s_o, div_uo, p_h, p_nk, oil, oil_l, rock, rock_l, control);
	use_time(g, TRUE, 0.);
	phgPrintf("solve s_o:           ");
	use_time(g, FALSE, 0.);
	phgSolverSolve(solver, TRUE, s_o, NULL);
	use_time(g, TRUE, 0.);
	phgDofFree(&div_uo);
	phgSolverDestroy(&solver);
}
static void
ParametersDump(PHASE *oil, PHASE *water, MEDIUM *rock)
{
	GRID *g = oil->P->g;
	SIMPLEX *e;
	FLOAT *p_bo, *p_bw, *p_phi, *p_kro, *p_krw, *p_dkro, *p_dkrw;
	phgPrintf("parameters value:\n");
	ForAllElements(g, e){
		p_bo = DofElementData(oil->B, e->index);
		p_bw = DofElementData(water->B, e->index);
		p_phi = DofElementData(rock->phi, e->index);
		p_kro = DofElementData(oil->Kr, e->index);
		p_dkro = DofElementData(oil->DKr, e->index);
		p_krw = DofElementData(water->Kr, e->index);
		p_dkrw = DofElementData(water->DKr, e->index);
		phgPrintf("bo:   %lf, bw:   %lf, phi:  %lf, kro:  %lf, dkro:   %lf, krw:   %lf, dkrw   %lf\n", p_bo[0], p_bw[0], p_phi[0], p_kro[0], p_dkro[0], p_krw[0], p_dkrw[0]);
	}
}
static void
PSDump(DOF *p_h, DOF *s_o)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_ph, *p_so;
	ForAllElements(g, e){
		p_ph = DofElementData(p_h, e->index);
		p_so = DofElementData(s_o, e->index);
//		if (e->region_mark == PRODWELL1 | e->region_mark == PRODWELL2 | e->region_mark ==PRODWELL3 | e->region_mark ==PRODWELL4)
		phgPrintf("p_h:   %lf,   s_o:     %lf\n", p_ph[0], p_so[0]);
	}
}
static void 
LastTime(DOF *bo, DOF *bw, DOF *phi, DOF *sw)
{
	GRID *g = bo->g;
	SIMPLEX *e;
	FLOAT *p_bo, *p_bw, *p_phi, *p_sw;
	phgPrintf("Last Time Parameters\n");
	ForAllElements(g, e){
		p_bo = DofElementData(bo, e->index);
		p_bw = DofElementData(bw, e->index);
		p_phi = DofElementData(phi, e->index);
		p_sw = DofElementData(sw, e->index);
		phgPrintf("bol;  %lf,  bwl:   %lf,  phil: %lf,  swl:  %lf\n", p_bo[0], p_bw[0], p_phi[0], p_sw[0]);
	}
}
static void
Find_Min_Max(DOF *s, FLOAT *min, FLOAT *max)
{
	GRID *g = s->g;
	SIMPLEX *e;
	FLOAT *p_s;
	FLOAT Min = 2.0, Max = 1e-16;
	ForAllElements(g, e){
		p_s = DofElementData(s, e->index);
		if (Max <= p_s[0]){
			Max = p_s[0];
		}
		if (Min >= p_s[0]){
			Min = p_s[0];
		}
	}
#if USE_MPI
	FLOAT tmp_min = Min;
	FLOAT tmp_max = Max;
	MPI_Allreduce(&tmp_min, &Min, 1, PHG_MPI_FLOAT, MPI_MIN, g->comm);
	MPI_Allreduce(&tmp_max, &Max, 1, PHG_MPI_FLOAT, MPI_MAX, g->comm);
#endif
	*min = Min;
	*max = Max;
}
static void
Disp_Saturation(DOF *s_o, DOF *p_h, PHASE *oil, PHASE *water, MEDIUM *rock)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	FLOAT so_min = 0, so_max = 0, sw_min = 0, sw_max = 0;
	FLOAT *p_sw, *p_so;
/*	ForAllElements(g, e){
		p_sw = DofElementData(s_w, e->index);
		p_so = DofElementData(s_o, e->index);
		p_sw[0] = 1. - p_so[0];
	}
	Find_Min_Max(s_w, &sw_min, &sw_max);
*/	Find_Min_Max(s_o, &so_min, &so_max);
	phgPrintf("-------------Saturation\n");
	phgPrintf("oil     min;  %le, max:  %le\n", so_min, so_max);
	phgPrintf("water   min:  %le, max:  %le\n", 1-so_max, 1-so_min);
}

static void 
water_fluidity(PHASE *oil, PHASE *water)
{
	GRID *g = oil->Source->g;
	SIMPLEX *e;
	FLOAT lambda_o = 0, lambda_w = 0;
	FLOAT *p_qo, *p_qw, *p_kro, *p_krw, *p_bo, *p_bw;
	ForAllElements(g, e){
			if (e->region_mark == 1){
				p_qw = DofElementData(water->Source, e->index);
				p_qo = DofElementData(oil->Source, e->index);
				p_kro = DofElementData(oil->Kr, e->index);
				p_krw = DofElementData(water->Kr, e->index);
				p_bo = DofElementData(oil->B, e->index);
				p_bw = DofElementData(water->B, e->index);
	//			phgPrintf("lambda_w:   %le,  lambda_o:  %le\n", lambda_w, lambda_o);
				lambda_o = p_kro[0] * p_bo[0] / oil->Mu;
				lambda_w = p_krw[0] * p_bw[0] / water->Mu;
				p_qw[0] = lambda_w / lambda_o * p_qo[0];
//				phgPrintf("warer source:  %lf\n", p_qw[0]);
			}
	}
}
void
update_production(PHASE *oil, PHASE *water, MEDIUM *rock, WELL *well)
{
    GRID *g = oil->P->g;
	SIMPLEX *e;
	FLOAT bar_k = 0, h =0;
	FLOAT *p_qw, *p_p, *p_k, *p_krw, *p_qo, *p_kro, *p_bo, *p_bw;
	ForAllElements(g, e){
		if (e->region_mark == 1){
				p_qw = DofElementData(water->Source, e->index);
				p_qo = DofElementData(oil->Source, e->index);
				p_p = DofElementData(oil->P, e->index);
				p_krw = DofElementData(water->Kr, e->index);
				p_kro = DofElementData(oil->Kr, e->index);
				p_bo = DofElementData(oil->B, e->index);
				p_bw = DofElementData(water->B, e->index);
				h = phgGeomGetDiameter(g, e);
				p_qo[0] = -0.002;
				p_qw[0] = -0.002;
				if (prod >= p_p[0]){
		//			p_qw[0] = -2. * M_PI * rock->perm * p_krw[0] * p_bw[0] * h * (prod - p_p[0]) / ((log(re / rw) + well->skin )* water->Mu);
		//			p_qo[0] = -2. * M_PI * rock->perm * p_kro[0] * p_bo[0] * h * (prod - p_p[0]) / ((log(re / rw) + well->skin )* oil->Mu);
			//		p_qw[0] = 1e-2 * p_qw[0];
			//		p_qo[0] = 1e-2 * p_qo[0];
				}
				if(prod < p_p[0]){
		//			p_qw[0] = 2. * M_PI * rock->perm * p_krw[0] * p_bw[0] * h * (prod - p_p[0]) / ((log(re / rw) + well->skin )* water->Mu);
		//			p_qo[0] = 2. * M_PI * rock->perm * p_kro[0] * p_bo[0] * h * (prod - p_p[0]) / ((log(re / rw) + well->skin )* oil->Mu);
			//		p_qw[0] = 1e-2 * p_qw[0];
			//		p_qo[0] = 1e-2 * p_qo[0];
				}
//				phgPrintf("qo:   %lf,   qw:     %lf\n", p_qo[0], p_qw[0]);
		}
	}
}
static void
Source_Dump(PHASE *oil, PHASE *water)
{
	GRID *g = oil->Source->g;
	SIMPLEX *e;
	FLOAT *p_qo, *p_qw;
	ForAllElements(g, e){
		if (e->region_mark == PRODWELL){
		p_qo = DofElementData(oil->Source, e->index);
		p_qw = DofElementData(water->Source, e->index);
		phgPrintf("qo:  %lf,  qw:  %lf\n", p_qo[0], p_qw[0]);
		}
	}
}
static int
bc_map(int bctype)
{
    switch (bctype) {
	case 2:
	    return NEUMANN;
	default:
	    return UNDEFINED;
    }
}
int
main(int argc, char *argv[])
{
    char *fn = "flow.mesh";
    GRID *g;
    char vtkfile[1000];
    static int pre_refines = 0;
	FLOAT stime = 1, ctime =0, ptime = 0;
	FLOAT T = 720;
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
	phgOptionsRegisterFilename("mesh_file", "Mesh File", (char **)&fn);
	phgOptionsRegisterFloat("step", "Time step", &stime);
	phgOptionsRegisterFloat("T", "computational damain: [0, T]", &T);
    phgInit(&argc, &argv);
	phgOptionsShowUsed(); 
    g = phgNewGrid(-1);
    phgImportSetBdryMapFunc(bc_map);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
//	phgPrintf("Import Mesh Ready\n");	
    phgRefineAllElements(g, pre_refines);
		if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
				(double)g->lif);
	/*******************************************/
		static FLOAT PRESSURE0 = 20.;
		static FLOAT SW0 = 0.4, MUW = 1e-10/3600;
		static FLOAT SO0 = 0.6, MUO = 1e-9/3600;
		PHASE *oil, *oil_l, *water, *water_l;
		MEDIUM *rock, *rock_l;
		COM_TRO *control;
		WELL *well;


		oil = phgGetPhase(g);
		oil_l = phgGetPhase(g);
		water = phgGetPhase(g);
		water_l = phgGetPhase(g);
		rock = phgGetMedium(g);
		rock_l = phgGetMedium(g);
		control = phgGetComTro(g);
		well = phgGetWellParameters(g);

		
		phgPrintf("=======================Computation Information=========================\n");
		phgPrintf("Flow     Model: SPE10 Upscale Small Model\n");
		phgPrintf("Flow     Types: oil-water system\n");
		phgPrintf("Total Elements:	  %d\n",DofGetDataCountGlobal(oil->P));
		phgPrintf("Total     DOFS:    %d\n",DofGetDataCountGlobal(oil->P) + DofGetDataCountGlobal(oil->U));
		phgPrintf("Total     Time:    %lf   Time      step:    %lf\n", control->T, control->dt);
		phgPrintf("***************************DEBUG INFORMATION**********************\n");
		/*struct data copy*/
		oil->S0 = SO0;
		oil_l->S0 = SO0;
		water->S0 = SW0;
		water_l->S0 = SW0;
		oil->Mu = MUO;
		oil_l->Mu = MUO;
		water->Mu = MUW;
		water_l->Mu = MUW;
		phgDofSetDataByValue(oil->S, oil->S0);
		phgDofSetDataByValue(oil_l->S, oil_l->S0);
		phgDofSetDataByValue(water->S, water->S0);
		phgDofSetDataByValue(water_l->S, water_l->S0);
		phgDofSetDataByValue(oil->P0, oil->PRESSURE0);
		phgDofSetDataByValue(oil_l->P0, oil_l->PRESSURE0);
		phgDofSetDataByValue(oil->P, oil->PRESSURE0);
		phgDofSetDataByValue(oil_l->P, oil->PRESSURE0);
		phgDofSetDataByValue(water->P0, water->PRESSURE0);
		phgDofSetDataByValue(water_l->P0, water_l->PRESSURE0);
		phgDofSetDataByValue(water->P, water->PRESSURE0);
		phgDofSetDataByValue(water_l->P, water_l->PRESSURE0);
		control->T = T;
		control->dt = stime;

	/*******************************************/
    DOF *p_h, *p0_h;
    p_h = phgDofNew(g, DOF_P0, 1, "p_h", DofNoAction);
	phgDofSetDataByValue(p_h, oil->PRESSURE0);
    p0_h = phgDofNew(g, DOF_P0, 1, "p0_h", DofNoAction);
	phgDofSetDataByValue(p0_h, oil->PRESSURE0);
    DOF *s_o;
    s_o = phgDofNew(g, DOF_P0, 1, "s_o", DofNoAction);
	phgDofSetDataByValue(s_o, oil->S0);
    DOF *u_o;
    u_o = phgDofNew(g, DOF_RT1, 1, "u_o", DofNoAction);
//	Well_init(oil, well);
//	Injection_Well(water, rock, well, control);
//	sprintf(vtkfile, "%s", "initial.vtk");
//	phgExportVTK(g, vtkfile, oil->Source, oil->S, water->Source, water->S, rock->phi, rock->perm, NULL);
	int flag = 0;
	while (ctime < T -1e-8 ){
		flag++;
		if (stime > 10){
			stime = 10;
		}
		ptime = ctime;
		ctime += stime;
		phgPrintf("==================================================================\n");
		phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)ctime);
		if (ctime > T){ 
			ctime = T;
			stime = T- ptime;
			phgPrintf("current time layer : [%lf, %lf]\n",(double)ptime, (double)ctime);
		}
		phgDofCopy(oil->P, &(oil_l->P), NULL, "pol");
		phgDofCopy(water->P, &(water_l->P), NULL, "pol");
		update_bo(oil_l);
		update_bw(water_l);
		update_phi(rock_l, oil_l);
		phgDofCopy(s_o, &(oil_l->S), NULL, "sol");
		int count = 0;
		while (TRUE){
			count++;
			phgPrintf("------------------------------%d'th\n", count);
			DOF *p_nk, *s_nk;
			p_nk = phgDofCopy(p_h, NULL, NULL, "p_nk");
			s_nk = phgDofCopy(s_o, NULL, NULL, "s_nk");
			/*update the parameters*/
			update_bo(oil);	
			update_dot_bo(oil);
			update_bw(water);
			update_dot_bw(water);
			create_kr(oil, water);
			create_dot_kr(oil, water);
			update_phi(rock, oil);
			update_dot_phi(rock);
			/*end of the parameters*/
#if DEBUG_PARS
			export_test_parameters(oil, water, oil_l, water_l, rock, rock_l, control, well);
			ParametersDump(oil, water, rock);
#endif
			update_production(oil, water, rock, well);
			//Source_Dump(oil, water);
	//		water_fluidity(oil, water);
//			Prod_Well_Rate(oil, water, control);
			phgPrintf("-----------Solve Equations\n");
			Solve_Pressure(u_o, p_h, s_o, oil, water, oil_l, water_l, rock, rock_l, control);
			phgDofCopy(p_h, &(oil->P), NULL, "oil pressure");
			phgDofCopy(p_h, &(water->P), NULL, "water pressure");
			phgDofCopy(u_o, &(oil->U), NULL, "oil velocity");
//			export_vtkfile(vtkfile, oil, water, rock, flag);
			Solver_Oileqn_So(s_o, p_h, u_o, p_nk, oil, oil_l, rock, rock_l, control);
			phgDofCopy(s_o, &(oil->S), NULL, "oil saturation");
//			PSDump(p_h, s_o);
//			exit(1);
//			flag++;
//			Disp_Saturation(oil->S, water->S);
//			Disp_Saturation(s_o, p_h, oil, water, rock);
//			update_bw(water, p_h);
//			create_kr(oil, water);
//			Solve_water_velocity(water, rock, p_h);
//			Test_PROD_BHP(oil, rock, well);
			FLOAT err_p = 0.0, norm_p = 0.0, TOL_p = 0.0;
			FLOAT err_s = 0.0, norm_s = 0.0, TOL_s = 0.0;
			err_p = L2(p_nk, p_h, 5);
			norm_p = (L2_norm(p_nk) > L2_norm(p_h))?L2_norm(p_nk):L2_norm(p_h);
			TOL_p = err_p / norm_p;
			err_s = L2(s_nk, s_o, 5);
			norm_s = (L2_norm(s_nk) > L2_norm(s_o))?L2_norm(s_nk):L2_norm(s_o);
			TOL_s = err_s / norm_s;
			phgPrintf("************************************\n");
			phgDofFree(&p_nk);
			phgDofFree(&s_nk);
			if (count > 20){
				stime = stime / 3;
			}
			if ((TOL_p < control->TOL_non)  & (TOL_s < control->TOL_non)){
	//		if ((oil->con_err < control->TOL_con) && (water->con_err < control->TOL_con)){
				phgPrintf("----------Final Results\n");
				phgPrintf("Pressure     err:            %le\n", TOL_p);
				phgPrintf("Saturation   err:            %le\n", TOL_s);
				phgPrintf("Nonlinear  iters:             %d\n", count);
				break;
			}
		}
		FLOAT pwf = 0;
		Well_Pressure(oil->P, &pwf);
		phgPrintf("time:   %lf,   press:   %lf\n", ctime, pwf);
#if EXPORT_VTK
		export_vtkfile(vtkfile, oil, water, rock, flag);
#endif
    }
    phgDofFree(&p_h);
    phgDofFree(&p0_h);
    phgDofFree(&u_o);
    phgDofFree(&s_o);
 	phgFreePhase(oil);
 	phgFreePhase(oil_l);
 	phgFreePhase(water);
 	phgFreePhase(water_l);
	phgFreeMedium(rock);
	phgFreeMedium(rock_l);
 	phgFreeComTro(control);
 	phgFreeWell(well);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
