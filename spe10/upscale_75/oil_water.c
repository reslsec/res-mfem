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
/* build linear system */
static double
use_time(GRID *g, BOOLEAN flag, double mflops)
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
#if USE_BLOCK_MATRIX
build_mat_vec(DOF *u_o, DOF *dp, DOF *ds, PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g)
#else
build_mat_vec(SOLVER *solver, DOF *u_o, DOF *dp, DOF *ds, PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control)
#endif
{
	GRID *g = u_o->g;
	SIMPLEX *e;
	FLOAT *p_bo, *p_bw, *p_bol, *p_bwl, *p_kro, *p_krw, *p_so, *p_sol, *p_phi, *p_dbw, *p_dbo, *p_dphi, *p_qo, *p_qw, *p_dkro, *p_dkrw, *p_p;
	INT N = u_o->type->nbas * u_o->dim;
	INT M = dp->type->nbas * dp->dim;
	INT I[N+M];
#if USE_BLOCK_MATRIX
	FLOAT mat_A[N][N], mat_TB[N][M], mat_B[M][N], mat_C[M][M], rhs_f[N], rhs_g[M];
	INT I[N], J[M];
	phgVecDisassemble(vec_f);
	phgVecDisassemble(vec_g);
#else
	FLOAT A0[N+M][N+M], rhs[N+M];
#endif
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
		p_dkro = DofElementData(oil->DKr, e->index);
		p_dkrw = DofElementData(water->DKr, e->index);
		p_qo = DofElementData(oil->Source, e->index);
		p_qw = DofElementData(water->Source, e->index);
		p_p = DofElementData(oil->P, e->index);
		FLOAT lambda_o = 0, lambda_w = 0;
		lambda_w = p_krw[0] * p_bw[0] / water->Mu;
		lambda_o = p_kro[0] * p_bo[0] / oil->Mu;
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
#if USE_BLOCK_MATRIX
				mat_A[i][j] = control->dt * oil->Mu * phgQuadBasPermBas(e, u_o, i, rock->perm, u_o, j, QUAD_DEFAULT) / (p_kro[0] * p_bo[0]);
#else
				A0[i][j] = control->dt * oil->Mu * phgQuadBasPermBas(e, u_o, i, rock->perm, u_o, j, QUAD_DEFAULT) / (p_kro[0] * p_bo[0]);
#endif

			}
		}
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
#if USE_BLOCK_MATRIX
				mat_TB[i][j] = -control->dt * phgQuadDivBasDotBas(e, u_o, i, dp, j, QUAD_DEFAULT);
#else
				A0[i][j + N] = -control->dt * phgQuadDivBasDotBas(e, u_o, i, dp, j, QUAD_DEFAULT);
				A0[j + N][i] = -control->dt * phgQuadDivBasDotBas(e, u_o, i, dp, j, QUAD_DEFAULT);
#endif
			}
		}
		for (i = 0; i < N; i++){
#if USE_BLOCK_MATRIX
			rhs_f[i] = control->dt * phgQuadDofTimesDivBas(e, oil->P, u_o, i, 5);
#else

			rhs[i] = control->dt * phgQuadDofTimesDivBas(e, oil->P, u_o, i, 5);
#endif
		}
		FLOAT wat_sw = 0, oil_sw = 0, wat_cp = 0, oil_cp = 0, beta = 0, quad = 0;
			for (i = 0; i < M; i++){
				for (j = 0; j < M; j++){
					wat_sw = p_phi[0] * p_bw[0] * phgQuadBasDotBas(e, ds, i, ds, j, QUAD_DEFAULT);
					oil_sw = p_phi[0] * p_bo[0] * phgQuadBasDotBas(e, ds, i, ds, j, QUAD_DEFAULT);
				}	
			}
			wat_cp = (1. - p_so[0]) * (p_bw[0]*p_dphi[0] + p_phi[0] * p_dbw[0]);
			oil_cp = (p_so[0]) * (p_bo[0] * p_dphi[0] + p_phi[0] * p_dbo[0]);
			beta = 1. / oil_sw + lambda_w / (lambda_o * wat_sw);
			for (i = 0; i < M; i++){
				for (j = 0; j < M; j++){
				 	quad =  phgQuadBasDotBas(e, dp, i, dp, j, QUAD_DEFAULT);
#if USE_BLOCK_MATRIX
					mat_C[i][j] = (wat_cp / wat_sw + oil_cp / oil_sw) * quad / beta;
#else
					A0[i + N][j + N] = -(wat_cp / wat_sw + oil_cp / oil_sw) * quad / beta;

#endif
				}
			}
		/*create oil rhs*/
		FLOAT quad_phi = 0., quad_phil = 0., quad_qo = 0, quad_qw = 0, quad_pnk=0;
		FLOAT rhs_oil = 0., rhs_wat = 0.;
		for (i = 0; i < M; i++){
			phgQuadDofTimesBas(e, oil->Source, dp, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, water->Source, dp, i, QUAD_DEFAULT, &quad_qw);
			phgQuadDofTimesBas(e, rock->phi, dp, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, rock_l->phi, dp, i, QUAD_DEFAULT, &quad_phil);
			rhs_wat = -control->dt * quad_qw + p_bw[0] * (1. - p_so[0]) * quad_phi - p_bwl[0] *  (1. - p_sol[0]) * quad_phil;
			rhs_oil = -control->dt * quad_qo + p_bo[0] * p_so[0] * quad_phi - p_bol[0] * p_sol[0] * quad_phil;
#if USE_BLOCK_MATRIX
			rhs_g[i] = (rhs_wat / wat_sw + rhs_oil / oil_sw) / beta;
#else

			rhs[i + N] = (rhs_wat / wat_sw + rhs_oil / oil_sw) / beta;
#endif
		}
		/* Handle Bdry Conditions*/
		for (i = 0; i < N; i++){
			if (phgDofGetElementBoundaryType(u_o, e, i) & (NEUMANN | DIRICHLET)){
#if USE_BLOCK_MATRIX
				bzero(mat_A[i], N * sizeof(mat_A[i][0]));
				bzero(mat_TB[i], M * sizeof(mat_TB[i][0]));
				mat_A[i][i] = 1.;
				rhs_f[i] = 0.;
#else
				bzero(A0[i], (N + M)* sizeof(A0[i][0]));
				A0[i][i] = 1.;
				rhs[i] = 0.;
#endif
			}
		}
#if USE_BLOCK_MATRIX
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
#else
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
		}
		for (i = 0; i < M; i++){
			I[i + N] = phgSolverMapE2L(solver, 1, e, i);
		}
		phgSolverAddMatrixEntries(solver, N + M, I, N + M, I, A0[0]);
		phgSolverAddRHSEntries(solver, N + M, I, rhs);
#endif
	}
#if USE_BLOCK_MATRIX
 	phgMatAssemble(A);
  	phgMatAssemble(B);
	phgMatAssemble(TB);	
	phgMatAssemble(C);
   	phgVecAssemble(vec_f);
	phgVecAssemble(vec_g);
#endif
}
static void
Solve_Pressure(SOLVER *solver, DOF *u_o, DOF *dp, DOF *ds, PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control)
{
		GRID *g = u_o->g;	
#if USE_BLOCK_MATRIX
	    MAT *A, *B, *TB, *C;
	    VEC *vec_f, *vec_g;
	    MAP *map_u, *map_p;
		/*Create MAP for Mat and Vec*/
	    map_p = phgMapCreate(dp, NULL);
		map_u = phgMapCreate(u_o, NULL);
		A = phgMapCreateMat(map_u, map_u);
		B = phgMapCreateMat(map_p, map_u);
		TB = phgMapCreateMat(map_u, map_p);
		C = phgMapCreateMat(map_p, map_p);
		vec_f = phgMapCreateVec(map_u, 1);
		vec_g = phgMapCreateVec(map_p, 1);
		phgPrintf("build  p_h:           ");
		use_time(g, FALSE, 0.);
		build_mat_vec(u_o, dp, ds, oil, water, oil_l, water_l, rock, rock_l, control, map_u, map_p, A, B, TB,C, vec_f, vec_g);
		use_time(g, TRUE, 0.);
		int nits_amg = 0, nits_uzawa = 0;
		use_time(g, FALSE, 0.);
		phgPrintf("Assemble H:             ");
		DOF *B_data = DivRT(dp, u_o, control);
        MAT *H = BTAB(A, B_data, dp, u_o);
      	phgDofFree(&B_data);
		use_time(g, TRUE, 0.);
    	/* new implementation */
		phgPrintf("solve p_h:           ");
		use_time(g, FALSE, 0.);
		nits_uzawa = uzawapcg(H, u_o, dp, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
		use_time(g, TRUE, 0.);
		phgMatDestroy(&H);
		phgPrintf("Max iter of AMG---------%d\n", nits_amg);
		phgPrintf("Nits: uzawa:----------------%d\n", nits_uzawa);
		/*
#if USE_MUMPS
		MAT *mat[4], *solver_mat;
		FLOAT coef[4];
		mat[0] = A;  mat[1] = TB;
		coef[0] = 1.; coef[1] = 1.;
		mat[2] = B;  mat[3] = C;
		coef[2] = 1.; coef[3] = -1.;
		solver_mat= phgMatCreateBlockMatrix(g, 2, 2, mat, coef, NULL);
		solver = phgSolverCreate(SOLVER_LASPACK, u_o, dp, NULL);
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
		phgPrintf("solve p_h:           ");
		phgSolverSolve(solver, FALSE, u_o, dp, NULL);
		use_time(g, TRUE, 0.);
		phgSolverDestroy(&solver);
		phgMatDestroy(&(solver_mat));
#endif
		*/
		phgMatDestroy(&A);
		phgMatDestroy(&B);
		phgMatDestroy(&TB);
		phgMatDestroy(&C);
		phgMapDestroy(&map_p);
		phgMapDestroy(&map_u);
		phgVecDestroy(&vec_g);
		phgVecDestroy(&vec_f);
#else
		phgPrintf("build  p_h:           ");
		use_time(g, FALSE, 0.);
		build_mat_vec(solver, u_o, dp, ds, oil, water, oil_l, water_l, rock, rock_l, control);
		use_time(g, TRUE, 0.);
		phgPrintf("solve p_h:           ");
		use_time(g, FALSE, 0.);
		phgSolverSolve(solver, FALSE, u_o, dp, NULL);
		use_time(g, TRUE, 0.);
#endif
}
#if 0
static void
build_mat_vec(DOF *u_o, DOF *dp, DOF *ds, PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g)
{
	GRID *g = u_o->g;
	SIMPLEX *e;
	FLOAT *p_bo, *p_bw, *p_bol, *p_bwl, *p_kro, *p_krw, *p_so, *p_sol, *p_phi, *p_dbw, *p_dbo, *p_dphi, *p_qo, *p_qw, *p_dkro, *p_dkrw, *p_p;
	INT N = u_o->type->nbas * u_o->dim;
	INT M = dp->type->nbas * dp->dim;
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
		p_dkro = DofElementData(oil->DKr, e->index);
		p_dkrw = DofElementData(water->DKr, e->index);
		p_qo = DofElementData(oil->Source, e->index);
		p_qw = DofElementData(water->Source, e->index);
		p_p = DofElementData(oil->P, e->index);
		FLOAT lambda_o = 0, lambda_w = 0;
		lambda_w = p_krw[0] * p_bw[0] / water->Mu;
		lambda_o = p_kro[0] * p_bo[0] / oil->Mu;
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				mat_A[i][j] = control->dt * oil->Mu * phgQuadBasPermBas(e, u_o, i, rock->perm, u_o, j, QUAD_DEFAULT) / (p_kro[0] * p_bo[0]);
			}
		}
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				mat_TB[i][j] = -control->dt * phgQuadDivBasDotBas(e, u_o, i, dp, j, QUAD_DEFAULT);
			}
		}
		for (i = 0; i < N; i++){
			rhs_f[i] = control->dt * phgQuadDofTimesDivBas(e, oil->P, u_o, i, 5);
		}
		FLOAT wat_sw = 0, oil_sw = 0, wat_cp = 0, oil_cp = 0, beta = 0, quad = 0;
			for (i = 0; i < M; i++){
				for (j = 0; j < M; j++){
					wat_sw = p_phi[0] * p_bw[0] * phgQuadBasDotBas(e, ds, i, ds, j, QUAD_DEFAULT);
					oil_sw = p_phi[0] * p_bo[0] * phgQuadBasDotBas(e, ds, i, ds, j, QUAD_DEFAULT);
				}	
			}
			wat_cp = (1. - p_so[0]) * (p_bw[0]*p_dphi[0] + p_phi[0] * p_dbw[0]);
			oil_cp = (p_so[0]) * (p_bo[0] * p_dphi[0] + p_phi[0] * p_dbo[0]);
			beta = 1. / oil_sw + lambda_w / (lambda_o * wat_sw);
			for (i = 0; i < M; i++){
				for (j = 0; j < M; j++){
				 	quad =  phgQuadBasDotBas(e, dp, i, dp, j, QUAD_DEFAULT);
					mat_C[i][j] = (wat_cp / wat_sw + oil_cp / oil_sw) * quad / beta;
				}
			}
		/*create oil rhs*/
		FLOAT quad_phi = 0., quad_phil = 0., quad_qo = 0, quad_qw = 0, quad_pnk=0;
		FLOAT rhs_oil = 0., rhs_wat = 0.;
		for (i = 0; i < M; i++){
			phgQuadDofTimesBas(e, oil->Source, dp, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, water->Source, dp, i, QUAD_DEFAULT, &quad_qw);
			phgQuadDofTimesBas(e, rock->phi, dp, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, rock_l->phi, dp, i, QUAD_DEFAULT, &quad_phil);
			rhs_wat = -control->dt * quad_qw + p_bw[0] * (1. - p_so[0]) * quad_phi - p_bwl[0] *  (1. - p_sol[0]) * quad_phil;
			rhs_oil = -control->dt * quad_qo + p_bo[0] * p_so[0] * quad_phi - p_bol[0] * p_sol[0] * quad_phil;
			rhs_g[i] = (rhs_wat / wat_sw + rhs_oil / oil_sw) / beta;
		}
		/* Handle Bdry Conditions*/
		for (i = 0; i < N; i++){
			if (phgDofGetElementBoundaryType(u_o, e, i) & (NEUMANN | DIRICHLET)){
				bzero(mat_A[i], N * sizeof(mat_A[i][0]));
				bzero(mat_TB[i], M * sizeof(mat_TB[i][0]));
				mat_A[i][i] = 1.;
				rhs_f[i] = 0.;
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
Solve_Pressure(DOF *u_o, DOF *dp, DOF *ds, PHASE *oil, PHASE *water, PHASE *oil_l, PHASE *water_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control)
{
		GRID *g = u_o->g;	
	    MAT *A, *B, *TB, *C;
	    VEC *vec_f, *vec_g;
	    MAP *map_u, *map_p;
		/*Create MAP for Mat and Vec*/
	    map_p = phgMapCreate(dp, NULL);
		map_u = phgMapCreate(u_o, NULL);
		A = phgMapCreateMat(map_u, map_u);
		B = phgMapCreateMat(map_p, map_u);
		TB = phgMapCreateMat(map_u, map_p);
		C = phgMapCreateMat(map_p, map_p);
		vec_f = phgMapCreateVec(map_u, 1);
		vec_g = phgMapCreateVec(map_p, 1);
			
		phgPrintf("build  p_h:           ");
		use_time(g, FALSE, 0.);
		build_mat_vec(u_o, dp, ds, oil, water, oil_l, water_l, rock, rock_l, control, map_u, map_p, A, B, TB,C, vec_f, vec_g);
		use_time(g, TRUE, 0.);
#if USE_UZAWA
		int nits_amg = 0, nits_uzawa = 0;
		use_time(g, FALSE, 0.);
		phgPrintf("Assemble H:             ");
		DOF *B_data = DivRT(dp, u_o, control);
        MAT *H = BTAB(A, B_data, dp, u_o);
      	phgDofFree(&B_data);
		use_time(g, TRUE, 0.);
    	/* new implementation */
		phgPrintf("solve p_h:           ");
		use_time(g, FALSE, 0.);
		nits_uzawa = uzawapcg(H, u_o, dp, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
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
		solver = phgSolverCreate(SOLVER_SUPERLU, u_o, dp, NULL);
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
		phgPrintf("solve p_h:           ");
		phgSolverSolve(solver, FALSE, u_o, dp, NULL);
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
#endif
static void
build_oil_so(SOLVER *solver, DOF *ds, DOF *div_uo, DOF *dp, PHASE *oil, PHASE *oil_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control)
{
	GRID *g = ds->g;
	SIMPLEX *e;
	int i,j;
	int N = ds->type->nbas;
	FLOAT A[N][N], rhs[N];
	INT I[N];
	FLOAT *p_phi, *p_bo, *p_bol, *p_so, *p_sol, *p_phil, *p_ph, *p_dphi, *p_db, *p_kro, *p_dkro, *p_p;
	ForAllElements(g, e){
		p_phi = DofElementData(rock->phi, e->index);
		p_dphi = DofElementData(rock->Dphi, e->index);
		p_bo  = DofElementData(oil->B, e->index);
		p_bol = DofElementData(oil_l->B, e->index);
		p_db = DofElementData(oil->DB, e->index);
		p_kro = DofElementData(oil->Kr, e->index);
		p_dkro = DofElementData(oil->DKr, e->index);
		p_so  = DofElementData(oil->S, e->index);
		p_sol  = DofElementData(oil_l->S, e->index);
		p_p  = DofElementData(oil->P, e->index);
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for (j = 0; j < N; j++){
					A[i][j] = p_phi[0] * p_bo[0] * phgQuadBasDotBas(e, ds, i, ds, j, QUAD_DEFAULT); 
			}
		}
		FLOAT oil_cp = 0;
		oil_cp = p_so[0] * (p_bo[0] * p_dphi[0] + p_phi[0] * p_db[0]);
		FLOAT quad_qo = 0, quad_phil = 0, quad_phi = 0., quad_divuo = 0, quad_dp = 0, quad_dpwf=0;
		for (i = 0; i < N; i++){
			phgQuadDofTimesBas(e, oil->Source, ds, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, dp, ds, i, QUAD_DEFAULT, &quad_dp);
			phgQuadDofTimesBas(e, rock->phi, ds, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, rock_l->phi, ds, i, QUAD_DEFAULT, &quad_phil);
			phgQuadDofTimesBas(e, div_uo, ds, i, QUAD_DEFAULT, &quad_divuo);
			rhs[i] = control->dt * quad_qo - control->dt * quad_divuo - oil_cp * quad_dp + p_sol[0] * p_bol[0] * quad_phil - p_bo[0] * p_so[0] * quad_phi;
		}
		for (i = 0; i < N; i++){
			phgSolverAddMatrixEntries(solver, 1, I+i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
}
static void
Solver_Oileqn_So(DOF *ds, DOF *u_o, DOF *dp, PHASE *oil, PHASE *oil_l, MEDIUM *rock, MEDIUM *rock_l, COM_TRO *control)
{
	GRID *g = ds->g;
	SOLVER *solver;
	DOF *div_uo;
	div_uo = phgDofDivergence(u_o, NULL, NULL, NULL);
	solver = phgSolverCreate(SOLVER_PCG, ds, NULL);
	phgPrintf("build s_o:           ");
	use_time(g, FALSE, 0.);
	build_oil_so(solver, ds, div_uo, dp, oil, oil_l, rock, rock_l, control);
	use_time(g, TRUE, 0.);
	phgPrintf("solve s_o:           ");
	use_time(g, FALSE, 0.);
	phgSolverSolve(solver, TRUE, ds, NULL);
	use_time(g, TRUE, 0.);
	phgDofFree(&div_uo);
	phgSolverDestroy(&solver);
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
Disp_Saturation(PHASE *oil, PHASE *water)
{
	GRID *g = oil->S->g;
	SIMPLEX *e;
	FLOAT so_min = 0, so_max = 0, sw_min = 0, sw_max = 0;
	FLOAT *p_sw, *p_so;
	Find_Min_Max(oil->S, &so_min, &so_max);
	phgPrintf("-------------Saturation\n");
	phgPrintf("oil     min:  %le, max:  %le\n", so_min, so_max);
	phgPrintf("water   min:  %le, max:  %le\n", 1-so_max, 1-so_min);
	Find_Min_Max(oil->Kr, &so_min, &so_max);
	phgPrintf("kro: min: %le, max:  %le\n", so_min, so_max);
	Find_Min_Max(water->Kr, &so_min, &so_max);
	phgPrintf("krw: min: %le, max:  %le\n", so_min, so_max);
}
void
Well_EffRadius(GRID *g, WELL *well, MEDIUM *rock)
{
	SIMPLEX *e;
	well->eff_radius = phgAlloc(well->Number * sizeof(FLOAT));
	if (well->eff_radius != NULL)
		printf("radius alloc success\n");
	bzero(well->eff_radius, well->Number * sizeof(FLOAT));
	FLOAT *p_perm, eff[4] = {111110, 111110, 111110, 111110}, tmp = 0;
	int i = 0;
	ForAllElements(g, e){
		p_perm = DofElementData(rock->perm, e->index);
		tmp = 0.28 * Sqrt(Sqrt(p_perm[1]/p_perm[0]) * pow(rock->DX, 2) + Sqrt(p_perm[0]/p_perm[1])*pow(rock->DY, 2))/(pow(p_perm[1]/p_perm[0], 1./4.)+pow(p_perm[0]/p_perm[1], 1./4.));
		for (i = 0; i < well->Number-1; i++){
			if (e->region_mark == well->region[i+1]){
				if(tmp < eff[i]){
					eff[i] = tmp;
					break;
				}
			}
		}
	}
#if USE_MPI
	FLOAT radius[well->Number-1];
	memcpy(radius, eff, (well->Number-1) * sizeof(FLOAT));
	for (i = 0; i < well->Number-1; i++){
		MPI_Allreduce(radius, eff, well->Number-1, PHG_MPI_FLOAT, MPI_MIN, g->comm);
	}
#endif
	for (i = 0; i < well->Number-1; i++){
		well->eff_radius[i+1] = eff[i];
	}
}
void
update_production(PHASE *oil, PHASE *water, MEDIUM *rock, WELL *well)
{
    GRID *g = oil->P->g;
	SIMPLEX *e;
	FLOAT *p_qw, *p_p, *p_perm, *p_krw, *p_qo, *p_kro, *p_bo, *p_bw;
	INT k =0, i = 0;
	FLOAT vol = 0, WI = 0, bar_k = 0, lambda_o = 0, lambda_w = 0;
	ForAllElements(g, e){
		if(e->region_mark == well->region[0]){
			vol += phgGeomGetVolume(g, e);
		}
	}
#if USE_MPI
	FLOAT vol0 = vol;
	MPI_Allreduce(&vol0, &vol, 1, MPI_DOUBLE, MPI_SUM, g->comm);
#endif
	phgPrintf("Producter Name------WI-----lambda_o-------lambda_w------BHP------Block-P-----Q_o-----Q_w\n");
	ForAllElements(g, e){
		for (i = 0; i < well->Number-1; i++){
			if (e->region_mark == well->region[i+1]){
				p_bw = DofElementData(water->B, e->index);
				p_bo = DofElementData(oil->B, e->index);
				p_krw = DofElementData(water->Kr, e->index);
				p_kro = DofElementData(oil->Kr, e->index);
				p_qo = DofElementData(oil->Source, e->index);
				p_qw = DofElementData(water->Source, e->index);
				p_p = DofElementData(oil->P, e->index);
				p_perm = DofElementData(rock->perm, e->index);
				FLOAT eff_radius = 0;
				eff_radius = 0.28 * Sqrt(Sqrt(p_perm[1]/p_perm[0]) * pow(rock->DX, 2) + Sqrt(p_perm[0]/p_perm[1])*pow(rock->DY, 2))/(pow(p_perm[1]/p_perm[0], 1./4.)+pow(p_perm[0]/p_perm[1], 1./4.));
				WI = 2. * M_PI * rock->DZ / (log(eff_radius/well->radius[i+1]));
				bar_k = Sqrt(p_perm[0] * p_perm[1]);
				lambda_o = bar_k * p_kro[0] * p_bo[0] / oil->Mu;
				lambda_w = bar_k * p_krw[0] * p_bw[0] / water->Mu;
		//		phgPrintf("prod %d, eff_radius:   %lf\n", i+1, eff_radius);
				p_qo[0] = WI * lambda_o * (well->PROD_BHP - p_p[0]) / vol;
				p_qw[0] = WI * lambda_w * (well->PROD_BHP - p_p[0]) / vol;
#if 0
				phgPrintf("Producter-%d,  %lf,  %lf,  %lf, %lf, %lf,  %lf,  %le\n", i+1, WI, lambda_o, lambda_w, well->PROD_BHP, p_p[0], p_qo[0], p_qw[0]);
#endif
				break;
			}
		}
	}
}
void
Test_PROD_BHP(PHASE *oil, MEDIUM *rock, WELL *well)
{
	GRID *g = oil->P->g;
	SIMPLEX *e;
	FLOAT *p_p;
	FLOAT bar_k = 0;
	FLOAT *p_qw, *p_perm, *p_krw, *p_qo, *p_kro, *p_bo;
	FLOAT press=0, h =0, vol =0;
	INT count[well->Number-1], i = 0;
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
	FLOAT lambda_o = 0, WI = 0;
	ForAllElements(g, e){
		p_p = DofElementData(oil->P, e->index);
		p_qo = DofElementData(oil->Source, e->index);
		p_perm = DofElementData(rock->perm, e->index);
		p_kro = DofElementData(oil->Kr, e->index);
		p_bo = DofElementData(oil->B, e->index);
		for (i = 0; i < well->Number-1; i++){
			if(e->region_mark == well->region[i+1]){
				bar_k = Sqrt(p_perm[0] * p_perm[1]);
				FLOAT eff_radius = 0;
				eff_radius = 0.28 * Sqrt(Sqrt(p_perm[1]/p_perm[0]) * pow(rock->DX, 2) + Sqrt(p_perm[0]/p_perm[1])*pow(rock->DY, 2))/(pow(p_perm[1]/p_perm[0], 1./ 4.)+pow(p_perm[0]/p_perm[1], 1./ 4.));
				WI = 2. * M_PI * rock->DZ / (log(eff_radius/well->radius[i+1]));
				lambda_o = bar_k * p_kro[0] * p_bo[0] / oil->Mu;
				if (p_p[0] >= well->PROD_BHP)
					press = p_qo[0] * vol / (WI * lambda_o) + p_p[0];
				sum[i] += press;
				count[i] += 1;
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
	phgPrintf("--------Test Compute PROD_BHP\n");
	for (i = 0; i < well->Number-1; i++){
		if (count[i] != 0)
			phgPrintf("Pord Well %d:    %lf Psi     %lf MPa\n", i+1, sum[i] / count[i] / PRESS_FACTOR, sum[i]/count[i]);
		else
			phgPrintf("Pord Well %d:    %lf Psi     %lf MPa\n", i+1, sum[i]  / PRESS_FACTOR, sum[i]);
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
void
Injection_Well(PHASE *water, WELL *well)
{
	GRID *g =water->Source->g;
	SIMPLEX *e;
	FLOAT *p_source;
	FLOAT vol = 0;
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
		if(e->region_mark == well->region[0]){
			p_source = DofElementData(water->Source, e->index);
			p_source[0] = well->INJE / vol;
		}
	}
}
void
pressure_recovery(PHASE *oil, PHASE *water, WELL *well)
{
	GRID *g = oil->P->g;
	SIMPLEX *e;
	FLOAT *p_po, *p_pw;
	INT i = 0;
	ForAllElements(g, e){
		p_po = DofElementData(oil->P, e->index);
		p_pw = DofElementData(water->P, e->index);
		for (i = 0; i < well->Number-1; i++){
			if(e->region_mark == well->region[i+1]){
				if (p_po[0] <= well->PROD_BHP){
					p_po[0] = well->PROD_BHP;
				}
				if (p_pw[0] <= well->PROD_BHP){
					p_pw[0] = well->PROD_BHP;
				}
				break;
			}
		}
	}
}
static void
refine_mesh(WELL *well, DOF *error)
{
    GRID *g = error->g;
    SIMPLEX *e;
	INT i;
    ForAllElements(g, e) {
		for (i = 0; i < well->Number-1; i++){
			if(e->region_mark == well->region[i+1]){
				*DofElementData(error, e->index) = 100;
				break;
			}
		}
    }
}
int
main(int argc, char *argv[])
{
    GRID *g;
    char *fn = "75.mesh";
    char *ff = "spe10_case2_75.inc";
    char *fp = "spe10_case2_75.dat";
    char vtkfile[1000];
    static int pre_refines = 0;
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
	phgOptionsRegisterFilename("input_file", "Input File", (char **)&fp);
	phgOptionsRegisterFilename("medium_file", "Miedum File", (char **)&ff);
    phgInit(&argc, &argv);
	phgOptionsShowUsed(); 
    g = phgNewGrid(-1);
    phgImportSetBdryMapFunc(bc_map);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, pre_refines);
		if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
				(double)g->lif);
	/*******************************************/
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

		/*Import Data*/
		read_data(fp, oil, water, well, control, rock);
		units_convert(oil, water, well, rock);
		/*Mark the well position*/
		Mark_All_Well(g, rock, well);

#if REFINE_WELL
		FLOAT tol= 5e-1;
    	DOF *error;
	    error = phgDofNew(g, DOF_P0, 1, "error", DofNoAction);
		refine_mesh(well, error);
		phgMarkRefine(MARK_DEFAULT, error, Pow(0.8,2), NULL, 0., 6,
				Pow(tol, 2) / g->nleaf_global);
		phgRefineMarkedElements(g);
    	phgDofFree(&error);
#endif
		
		/*Set Init Data*/
		phgDofSetDataByValue(oil->S, oil->S0);
		phgDofSetDataByValue(oil_l->S, oil->S0);
		phgDofSetDataByValue(water->S, water->S0);
		phgDofSetDataByValue(water_l->S, water->S0);
		phgDofSetDataByValue(oil->P0, oil->PRESSURE0);
		phgDofSetDataByValue(oil_l->P0, oil->PRESSURE0);
		phgDofSetDataByValue(oil->P, oil->PRESSURE0);
		phgDofSetDataByValue(oil_l->P, oil->PRESSURE0);
		phgDofSetDataByValue(water->P0, water->PRESSURE0);
		phgDofSetDataByValue(water_l->P0, water->PRESSURE0);
		phgDofSetDataByValue(water->P, water->PRESSURE0);
		phgDofSetDataByValue(water_l->P, water->PRESSURE0);
	Init_Medium(ff, rock);
	Well_EffRadius(g, well, rock);
	export_test_data(oil, water, oil_l, water_l, rock, control, well);
	phgDofCopy(rock->phi0, &(rock->phi), NULL, "phi");
	phgDofCopy(rock->phi0, &(rock_l->phi0), NULL, "phi_l0");
	phgDofCopy(rock->phi0, &(rock_l->phi), NULL, "phi_l");
	export_vtkfile(vtkfile, oil, water, rock, 0);
	/*******************************************/
    DOF *dp;
    dp = phgDofNew(g, DOF_P0, 1, "dp", DofNoAction);
    DOF *ds;
    ds = phgDofNew(g, DOF_P0, 1, "ds", DofNoAction);
    DOF *u_o;
    u_o = phgDofNew(g, DOF_RT1, 1, "u_o", DofNoAction);
	phgPrintf("====================Upcasle SEP10 Samll Model===================\n");
	phgPrintf("Flow     Types:    oil-water system\n");
	phgPrintf("Total Elements:	  %d\n",DofGetDataCountGlobal(oil->P));
	phgPrintf("Total     DOFS:    %d\n",DofGetDataCountGlobal(oil->P) + DofGetDataCountGlobal(oil->U));
	phgPrintf("Total     Time:    %lf   Time      step:    %lf\n", control->T, control->dt);
	phgPrintf("oil Mu:    %le, water Mu:   %le\n", oil->Mu, water->Mu);
	int flag = 0;
	while (control->ct < control->T -1e-8){
		if ((control->dt < control->max_dt)){
			control->dt = 2 * control->dt;
		}
		if (control->dt >= control->max_dt){
			control->dt = control->max_dt;
		}
		control->pt = control->ct;
		control->ct += control->dt;
		phgPrintf("==================================================================\n");
		phgPrintf("current time layer: [%lf, %lf]\n", (double)control->pt, (double)control->ct);
		if (control->ct > control->T){ 
			control->ct = control->T;
			control->dt = control->T- control->pt;
			phgPrintf("current time layer : [%lf, %lf]\n",(double)control->pt, (double)control->ct);
		}
		//Injection_Well(water, well);
		phgDofCopy(oil->P, &(oil_l->P), NULL, "pol");
		phgDofCopy(water->P, &(water_l->P), NULL, "pol");
		phgDofCopy(oil->S, &(oil_l->S), NULL, "sol");
		update_bo(oil_l);
		update_bw(water_l);
		update_phi(rock_l, oil_l);
		int count = 0;
		while (TRUE){
			/*
			if (count >= 20){
				control->dt = control->dt / 3;
				control->ct = control->pt + control->dt;
				phgPrintf("Adjust Time Step\n");
				phgPrintf("current time layer: [%lf, %lf]\n", (double)control->pt, (double)control->ct);
				phgDofCopy(oil_l->P, &(oil->P), NULL, NULL);
				phgDofCopy(oil_l->P, &(water->P), NULL, NULL);
				phgDofCopy(oil_l->S, &(oil->S), NULL, NULL);
				phgDofCopy(tmp, &dp, NULL, NULL);
				count = 0;
			}*/
			count++;
			DOF *p_nk, *s_nk;
			p_nk = phgDofCopy(oil->P, NULL, NULL, NULL);
			s_nk = phgDofCopy(oil->S, NULL, NULL, NULL);
			/*update the parameters*/
			update_bo(oil);	
			update_dot_bo(oil);
			update_bw(water);
			update_dot_bw(water);
			create_kr(oil, water);
			create_dot_kr(oil, water);
			update_phi(rock, oil);
			update_dot_phi(rock);
			update_production(oil, water, rock, well);
			SOLVER *solver;
			solver = phgSolverCreate(SOLVER_DEFAULT, u_o, dp, NULL);
			phgPrintf("-----------Solve Equations\n");
			Solve_Pressure(solver, u_o, dp, ds, oil, water, oil_l, water_l, rock, rock_l, control);
			phgSolverDestroy(&solver);
			phgDofAXPBY(1.0, dp, 1, &(oil->P));
			phgDofAXPBY(1.0, dp, 1, &(water->P));
			phgDofCopy(u_o, &(oil->U), NULL, "oil velocity");
			Test_PROD_BHP(oil, rock, well);
			Solver_Oileqn_So(ds, u_o, dp, oil, oil_l, rock, rock_l, control);
			phgDofAXPBY(1.0, ds, 1, &(oil->S));
			Disp_Saturation(oil, water);
			FLOAT err_p = 0.0, norm_p = 0.0, TOL_p = 0.0;
			FLOAT err_s = 0.0, norm_s = 0.0, TOL_s = 0.0;
			err_p = L2(p_nk, oil->P, 5);
			norm_p = (L2_norm(p_nk) > L2_norm(oil->P))?L2_norm(p_nk):L2_norm(oil->P);
			TOL_p = err_p / norm_p;
			err_s = L2(s_nk, oil->S, 5);
			norm_s = (L2_norm(s_nk) > L2_norm(oil->S))?L2_norm(s_nk):L2_norm(oil->S);
			TOL_s = err_s / norm_s;
			phgPrintf("************************************\n");
			phgDofFree(&p_nk);
			phgDofFree(&s_nk);
			phgPrintf("Relative err:   %le,  %le\n", TOL_p, TOL_s);
			if ((TOL_p < control->TOL_non) && (TOL_s < control->TOL_non)){
				phgPrintf("----------Final Results\n");
				phgPrintf("Pressure     err:            %le\n", TOL_p);
				phgPrintf("Saturation   err:            %le\n", TOL_s);
				phgPrintf("Nonlinear  iters:             %d\n", count);
				break;
			}
		}
		flag++;
		FLOAT total_oil = 0;
		result_export(oil, water, well, control, &total_oil);
#if EXPORT_VTK
		export_vtkfile(vtkfile, oil, water, rock, flag);
#endif
    }
    phgDofFree(&dp);
    phgDofFree(&ds);
    phgDofFree(&u_o);
 	phgFreePhase(water);
 	phgFreePhase(oil);
 	phgFreePhase(oil_l);
 	phgFreePhase(water_l);
	phgFreeMedium(rock);
	phgFreeMedium(rock_l);
 	phgFreeComTro(control);
 	phgFreeWell(well);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
