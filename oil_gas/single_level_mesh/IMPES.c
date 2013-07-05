 /* This sample code solves the oil_water_phase equation:         *
 * ****************************************************************
 * problem: p_{t} -\Delta{p} = func_f(x, t) \in \Omega X [0, T]   *
 *          p = func_g (x, t) \in \partial\Omega X [0, T]         * 
 *          p = func_p0       \in \Omega t = 0.0                  *
 *WARNING! The unit is h, m, mpa!!                                *
 ******************************************************************/
#include "phg.h"
#include <string.h>
#include <math.h>
#include "parameter.c"
#include "well.c"
#include "quadfunc.c"
#include "uzawa.c"
#define USE_UZAWA 1
#define USE_BLOCK 0
#define DEBUG 0
#define USE_OIL_EQN 1
#define USE_GAS_EQN 0
static DOF_TYPE DOF_RT1_;
#define DOF_RT1 (&DOF_RT1_)
# define BasisOrder(u, e, i) (!DofIsHP(u) ? (u)->type->order : \
	(u)->hp->info->types[(u)->hp->elem_order[e->index]]->order)
static double
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
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
/*build_linear_system*/
build_mat_vec(DOF *u_o, DOF *p_h, DOF *p_h_newton, DOF *s_o, DOF *s_o_l, DOF *mu_o, DOF *b_o, DOF *b_o_l, DOF *kro, DOF *phi, DOF *phi_l, DOF *dot_phi, DOF *dot_bo, DOF *Rs, DOF *Rs_l, DOF *dot_Rs, DOF *q_o, DOF *mu_g, DOF *b_g, DOF *b_g_l, DOF *krg, DOF *dot_bg, DOF *q_g, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g)
{
	GRID *g = u_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_sol, *p_bo, *p_bol, *p_kro, *p_phi, *p_phil, *p_dotbo, *p_Rs, *p_dotRs, *p_bg, *p_bgl, *p_krg, *p_dotbg, *p_muo, *p_mug, *p_dotphi, *p_Rsl;
	INT N = u_o->type->nbas * u_o->dim;
	INT M = p_h->type->nbas * p_h->dim;
	INT I[N], J[M];
	FLOAT mat_A[N][N], mat_TB[N][M], mat_B[M][N], mat_C[M][M], rhs_f[N], rhs_g[M];
	phgVecDisassemble(vec_f);
	phgVecDisassemble(vec_g);
	int i, j, k;
	ForAllElements(g, e){
		p_so = DofElementData(s_o, e->index);
		p_sol = DofElementData(s_o_l, e->index);
		p_bo = DofElementData(b_o, e->index);
		p_bol = DofElementData(b_o_l, e->index);
		p_kro = DofElementData(kro, e->index);
		p_phi = DofElementData(phi, e->index);
		p_phil = DofElementData(phi_l, e->index);
		p_dotphi = DofElementData(dot_phi, e->index);
		p_dotbo = DofElementData(dot_bo, e->index);
		p_Rs = DofElementData(Rs, e->index);
		p_Rsl = DofElementData(Rs_l, e->index);
		p_dotRs = DofElementData(dot_Rs, e->index);
		p_bg = DofElementData(b_g, e->index);
		p_bgl = DofElementData(b_g_l, e->index);
		p_krg= DofElementData(krg, e->index);
		p_dotbg = DofElementData(dot_bg, e->index);
		p_muo = DofElementData(mu_o, e->index);
		p_mug = DofElementData(mu_g, e->index);
		/*Create two inverse*/
		FLOAT oil_so = 0, gas_so = 0;
		for (i = 0; i < M; i++){
			for (j = 0; j < M; j++){
				oil_so = p_phi[0] * p_bo[0] * phgQuadBasDotBas(e, s_o, i, s_o, j, QUAD_DEFAULT);
				gas_so = p_phi[0] * (p_bg[0] - p_Rs[0] * p_bo[0]) * phgQuadBasDotBas(e, s_o, i, s_o, j, QUAD_DEFAULT);
			}
		}
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				mat_A[i][j] = stime * phgQuadBasDotBas(e, u_o, i, u_o, j, QUAD_DEFAULT) / (p_kro[0] * p_bo[0] * p_muo[0] * K);
			}
		}
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				mat_TB[i][j] = -stime * phgQuadDivBasDotBas(e, u_o, i, p_h, j, QUAD_DEFAULT);
			}
		}
		FLOAT beta = 0;
		beta = 1. / oil_so + (p_Rs[0] + p_krg[0] * p_bg[0] * p_mug[0] / (p_kro[0] * p_bo[0] * p_muo[0])) / gas_so;
		FLOAT quad = 0., oil_cp = 0., gas_cp = 0;
		oil_cp = p_so[0] * (p_bo[0] * p_dotphi[0] + p_phi[0] * p_dotbo[0]);
		gas_cp = p_so[0] * p_phi[0] * p_bo[0] * p_dotRs[0] + p_Rs[0] * oil_cp + (1. - p_so[0]) * (p_phi[0] * p_dotbg[0] + p_bg[0] * p_dotphi[0]);
		for (i = 0; i < M; i++){
			for (j = 0; j < M; j++){
				quad = phgQuadBasDotBas(e, p_h, i, p_h, j, QUAD_DEFAULT);
				mat_C[i][j] = (oil_cp /oil_so + gas_cp / gas_so) * quad / beta;
			}
		}
		/*Create rhs*/
		FLOAT quad_qo = 0;
		FLOAT quad_qg = 0, quad_phi = 0; 
		FLOAT quad_phil = 0, quad_pnew = 0;
		FLOAT rhs_oil = 0, rhs_gas = 0;
		for (i = 0; i < M; i++){
			phgQuadDofTimesBas(e, q_o, p_h, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, p_h_newton, p_h, i, QUAD_DEFAULT, &quad_pnew);
			phgQuadDofTimesBas(e, phi_l, p_h, i, QUAD_DEFAULT, &quad_phil);
			
			phgQuadDofTimesBas(e, q_g, p_h, i, QUAD_DEFAULT, &quad_qg);
			phgQuadDofTimesBas(e, phi, p_h, i, QUAD_DEFAULT, &quad_phi);
			rhs_oil = -stime * quad_qo - oil_cp * quad_pnew - p_sol[0] * p_bol[0] * quad_phil;
			rhs_gas = -stime * quad_qg + p_bg[0] * quad_phi - gas_cp * quad_pnew - (p_Rsl[0] * p_bol[0] * p_sol[0] + p_bgl[0] * (1. - p_sol[0])) * quad_phil;
			rhs_g[i] = (rhs_oil / oil_so + rhs_gas / gas_so) / beta;
		}
		for (i = 0; i < N; i++){
			rhs_f[i] = 0.0;
		}
		/* Handle Bdry Conditions */
		for (i = 0; i < N; i++){
			if (phgDofGetElementBoundaryType(u_o, e, i) & (NEUMANN | DIRICHLET)){
				bzero(mat_A[i], N *sizeof(mat_A[i][0]));
				bzero(mat_TB[i], M *sizeof(mat_TB[i][0]));
				for (j = 0; j < N; j++){
					mat_A[j][i] = 0.0;
				}
				mat_A[i][i] = 1.0;
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
	phgMatAssemble(TB);
	phgMatAssemble(B);
	phgMatAssemble(C);
	phgVecAssemble(vec_f);
	phgVecAssemble(vec_g);
}
static void
oil_conservation(DOF *p_h, DOF *p0_h, DOF *phi_l, DOF *b_o_l, DOF *s_o, DOF *s_o_l, DOF *q_o, DOF *u_o)
{
	GRID *g = p_h->g;
   	DOF *phi, *b_o;		
	phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
	update_phi(p_h, p0_h, phi);
 	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
	update_bo(p_h, b_o);	
	
	DOF *tmp;
	tmp = phgDofNew(g, DOF_P0, 1, "tmp", DofNoAction);
	phgDofSetDataByValue(tmp, 1.0);
	DOF *div_uo;
	div_uo = phgDofDivergence(u_o, NULL, NULL, NULL);
	SIMPLEX *e;
	FLOAT *p_bo, *p_phi, *p_bol, *p_phil;
	FLOAT v1 =0, v2=0, div_u =0;
	FLOAT input = 0, output = 0, div = 0;
	ForAllElements(g, e){
		p_bo = DofElementData(b_o, e->index);
		p_phi = DofElementData(phi, e->index);
		p_bol = DofElementData(b_o_l, e->index);
		p_phil = DofElementData(phi_l, e->index);
		div_u += phgQuadDofDotDof(e, div_uo, tmp, 5);
		v1 += phgQuadDofDotDof(e, q_o, tmp, 5);
		v2 += p_bo[0] * p_phi[0] * phgQuadDofDotDof(e, s_o, tmp, 5) / stime
			- p_bol[0] * p_phil[0] * phgQuadDofDotDof(e, s_o_l, tmp, 5) / stime;
	}
#if USE_MPI
    output = v2, input = v1, div = div_u;
    MPI_Allreduce(&v1, &input, 1, PHG_MPI_FLOAT, MPI_SUM, p_h->g->comm);
    MPI_Allreduce(&v2, &output, 1, PHG_MPI_FLOAT, MPI_SUM,p_h->g->comm);
    MPI_Allreduce(&div_u, &div, 1, PHG_MPI_FLOAT, MPI_SUM,p_h->g->comm);
#else
    input = v1;
    output = v2;
    div = div_u;
#endif
	phgPrintf("Oil_Conserve: LHS = %le,RHS = %le\n\n", output, input);
	phgPrintf("(Div uo, 1) = %lf\n", div);
	phgPrintf("v2 + (div_uo, 1) = %lf , v1 = %lf\n", output+div, input);
	phgDofFree(&phi);
	phgDofFree(&b_o);
	phgDofFree(&tmp);
	phgDofFree(&div_uo);
}
static void 
gas_fluidity(DOF *q_o, DOF *q_g, DOF *kro, DOF *krg, DOF *b_o, DOF *b_g, DOF *mu_o, DOF *mu_g, DOF *Rs)
{
	GRID *g = q_o->g;
	SIMPLEX *e;
	FLOAT lambda_o = 0, lambda_g = 0;
	FLOAT *p_qo, *p_qg, *p_kro, *p_krg, *p_bo, *p_bg, *p_muo, *p_mug, *p_Rs;
	ForAllElements(g, e){
			p_qg = DofElementData(q_g, e->index);
			p_qo = DofElementData(q_o, e->index);
			p_kro = DofElementData(kro, e->index);
			p_krg = DofElementData(krg, e->index);
			p_bo = DofElementData(b_o, e->index);
			p_bg = DofElementData(b_g, e->index);
			p_muo = DofElementData(mu_o, e->index);
			p_mug = DofElementData(mu_g, e->index);
			p_Rs = DofElementData(Rs, e->index);
			lambda_o = p_kro[0] * p_bo[0] * p_muo[0];
			lambda_g = p_krg[0] * p_bg[0] * p_mug[0];
			p_qg[0] = (lambda_g + p_Rs[0] * lambda_o) / lambda_o * p_qo[0];
	}
}
static void
Solve_Pressure(DOF *u_o, DOF *p_h, DOF *p0_h, DOF *q_o, DOF *s_o, DOF *s_o_l, DOF *q_g, DOF *p_h_newton, DOF *phi_l, DOF *b_o_l, DOF *b_g_l, DOF *Rs_l)
{
		GRID *g = u_o->g;	
	     MAT *A, *B, *TB, *C;
	     VEC *vec_f, *vec_g;
	     MAP *map_u, *map_p;
#if USE_BLOCK
		SOLVER *solver;
		MAT *pmat[4];
		FLOAT coef[4];
#endif
     	/* The parameter DOF */
 	  	DOF *phi, *b_o, *b_g, *kro, *krg, *Rs, *mu_o, *mu_g, *dot_bo, *dot_bg, *dot_Rs, *dot_phi;
		/*Create MAP for Mat and Vec*/
	     map_p = phgMapCreate(p_h, NULL);
		map_u = phgMapCreate(u_o, NULL);
		A = phgMapCreateMat(map_u, map_u);
		A->handle_bdry_eqns = FALSE;
		B = phgMapCreateMat(map_p, map_u);
		B->handle_bdry_eqns = FALSE;
		TB = phgMapCreateMat(map_u, map_p);
		TB->handle_bdry_eqns = FALSE;
		C = phgMapCreateMat(map_p, map_p);
		C->handle_bdry_eqns = FALSE;
		vec_f = phgMapCreateVec(map_u, 1);
		vec_g = phgMapCreateVec(map_p, 1);
			
		phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
		update_phi(p_h, p0_h, phi);
 		b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
		update_bo(p_h, b_o);
 		mu_o = phgDofNew(g, DOF_P0, 1, "mu_o", DofNoAction);
		update_muo(p_h, mu_o);
 		b_g = phgDofNew(g, DOF_P0, 1, "b_g", DofNoAction);
		update_bg(p_h, b_g);
 		mu_g = phgDofNew(g, DOF_P0, 1, "mu_g", DofNoAction);
		update_mug(p_h, mu_g);
 		Rs = phgDofNew(g, DOF_P0, 1, "Rs", DofNoAction);
		update_Rs(p_h, Rs);
		kro = phgDofNew(g, DOF_P0, 1, "kro", DofNoAction);
		create_kro(s_o, kro);
		krg = phgDofNew(g, DOF_P0, 1, "krg", DofNoAction);
		create_krg(s_o, krg);
		dot_bo = phgDofNew(g, DOF_P0, 1, "dot_bo", DofNoAction);
		dot_bg = phgDofNew(g, DOF_P0, 1, "dot_bg", DofNoAction);
		dot_Rs = phgDofNew(g, DOF_P0, 1, "dot_Rs", DofNoAction);
		dot_phi = phgDofNew(g, DOF_P0, 1, "dot_phi", DofNoAction);
		update_dot_bo(p_h, dot_bo);
		update_dot_bg(p_h, dot_bg);
		update_dot_Rs(p_h, dot_Rs);
		update_dot_phi(p_h, dot_phi);
		gas_fluidity(q_o, q_g, kro, krg, b_o, b_g, mu_o, mu_g, Rs);
		build_mat_vec(u_o, p_h, p_h_newton, s_o, s_o_l, mu_o, b_o, b_o_l, kro, phi, phi_l, dot_phi, dot_bo, Rs, Rs_l, dot_Rs, q_o, mu_g, b_g, b_g_l, krg, dot_bg, q_g, map_u, map_p, A, B, TB, C, vec_f, vec_g);
		phgDofFree(&phi);
		phgDofFree(&mu_o);
		phgDofFree(&b_o);
		phgDofFree(&kro);
		phgDofFree(&mu_g);
		phgDofFree(&b_g);
		phgDofFree(&krg);
		phgDofFree(&Rs);
		phgDofFree(&dot_bo);
		phgDofFree(&dot_bg);
		phgDofFree(&dot_Rs);
		phgDofFree(&dot_phi);
#if USE_UZAWA
		int nits_amg = 0;
        int nits_uzawa = 0;
		double time_amg = 0;
		elapsed_time(g, FALSE, 0.);
		DOF *B_data = DivRT(p_h, u_o);
        MAT *H = BTAB(A, B_data, p_h, u_o);
      	phgDofFree(&B_data);
		phgPrintf("Assemble H:             ");
		elapsed_time(g, TRUE, 0.);
		phgPrintf("solve p_h:              \n");
	//	nits_uzawa = phgUzawa(H, u_w, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg, &time_amg);
		elapsed_time(g, FALSE, 0.);
		nits_uzawa = uzawapcg(H, u_o, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
		elapsed_time(g, TRUE, 0.);
		phgPrintf("MAx iter of AMG---------%d\n", nits_amg );
		phgPrintf("Nits: uzawa:----------------%d\n", nits_uzawa);
		phgMatDestroy(&H);
#endif
#if USE_BLOCK
		solver = phgSolverCreate(SOLVER_MUMPS, u_o, p_h, NULL);
		solver->mat->handle_bdry_eqns = FALSE;
		pmat[0] = A;  pmat[1] = TB;
		coef[0] = 1.; coef[1] = 1.;
		pmat[2] = B;  pmat[3] = C;
		coef[2] = 1.; coef[3] = -1.;
		phgMatDestroy(&solver->mat);
		solver->mat = phgMatCreateBlockMatrix(g, 2, 2, pmat, coef, NULL);
		solver->rhs->mat = solver->mat;
    		INT N = vec_f->map->nlocal;
    		INT M = vec_g->map->nlocal;
		memcpy(solver->rhs->data, vec_f->data, sizeof(*vec_f->data) * N);
		memcpy(solver->rhs->data+N, vec_g->data, sizeof(*vec_g->data) * M);
		phgSolverSolve(solver, TRUE, u_o, p_h, NULL);
		phgSolverDestroy(&solver);
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
gas_conservation(DOF *p_h, DOF *p0_h, DOF *phi_l, DOF *b_g_l, DOF *b_o_l, DOF *s_o, DOF *s_o_l, DOF *Rs_l, DOF *q_g)
{
	GRID *g = p_h->g;
   	DOF *phi, *b_g,  *b_o, *Rs;
	phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
	update_phi(p_h, p0_h, phi);
 	b_g = phgDofNew(g, DOF_P0, 1, "b_g", DofNoAction);
	update_bg(p_h, b_g);	
 	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
	update_bo(p_h, b_o);	
 	Rs = phgDofNew(g, DOF_P0, 1, "Rs", DofNoAction);
	update_Rs(p_h, Rs);
	
	DOF *tmp;
	tmp = phgDofNew(g, DOF_P0, 1, "tmp", DofNoAction);
	phgDofSetDataByValue(tmp, 1.0);
	SIMPLEX *e;
	FLOAT *p_bg, *p_bgl, *p_bo, *p_bol, *p_Rs, *p_Rsl, *p_so, *p_sol;
	FLOAT v1 =0, v2=0;
	FLOAT output = 0, input = 0;
	ForAllElements(g, e){
		p_bo = DofElementData(b_o, e->index);
		p_bol = DofElementData(b_o_l, e->index);
		p_so = DofElementData(s_o, e->index);
		p_sol = DofElementData(s_o_l, e->index);
		p_Rs = DofElementData(Rs, e->index);
		p_Rsl = DofElementData(Rs_l, e->index);
		p_bgl = DofElementData(b_g_l, e->index);
		p_bg = DofElementData(b_g, e->index);
		v1 += phgQuadDofDotDof(e, q_g, tmp, 5);
		v2 += (p_Rs[0] * p_bo[0] * p_so[0] + p_bg[0] * (1. - p_so[0])) * phgQuadDofDotDof(e, phi, tmp, 5) / stime
			- (p_Rsl[0] * p_bol[0] * p_sol[0] + p_bgl[0] * (1. -p_sol[0]))* phgQuadDofDotDof(e, phi_l, tmp, 5) / stime;
	}
#if USE_MPI
    input = v1, output = v2;
    MPI_Allreduce(&v1, &input, 1, PHG_MPI_FLOAT, MPI_SUM, p_h->g->comm);
    MPI_Allreduce(&v2, &output, 1, PHG_MPI_FLOAT, MPI_SUM,p_h->g->comm);
#else
    input = v1;
    output = v2;
#endif
	phgPrintf("Gas_Conserve: LHS = %le,RHS = %le\n\n", output, input);
	phgDofFree(&phi);
	phgDofFree(&b_o);
	phgDofFree(&b_g);
	phgDofFree(&Rs);
	phgDofFree(&tmp);
}
#if USE_OIL_EQN
static void
build_oileqn_so(SOLVER *solver, DOF *s_o, DOF *s_o_l, DOF *div_uo, DOF *q_o, DOF *p_h, DOF *p_h_newton, DOF *b_o, DOF *phi, DOF *dot_phi, DOF *dot_bo, DOF *phi_l, DOF *b_o_l)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	int i,j;
	int N = s_o->type->nbas;
	FLOAT A[N][N], rhs[N];
	INT I[N];
	FLOAT *p_phi, *p_bo, *p_so, *p_dotbo, *p_dotphi, *p_phil, *p_bol;
	ForAllElements(g, e){
		p_phi = DofElementData(phi, e->index);
		p_phil = DofElementData(phi_l, e->index);
		p_bo  = DofElementData(b_o, e->index);
		p_bol = DofElementData(b_o_l, e->index);
		p_so  = DofElementData(s_o, e->index);
		p_dotbo = DofElementData(dot_bo, e->index);
		p_dotphi = DofElementData(dot_phi, e->index);
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for (j = 0; j < N; j++){
				A[i][j] =  p_bo[0] * p_phi[0] * phgQuadBasDotBas(e, s_o, i, s_o, j, QUAD_DEFAULT); 
			}
		}
		FLOAT quad_qo = 0, quad_p = 0, quad_pnew = 0, quad_sol = 0, quad_divuo = 0;
		FLOAT oil_cp = 0;
		oil_cp = p_so[0] * (p_bo[0] * p_dotphi[0] + p_phi[0] * p_dotbo[0]);
		for (i = 0; i < N; i++){
			phgQuadDofTimesBas(e, q_o, s_o, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, p_h, s_o, i, QUAD_DEFAULT, &quad_p);
			phgQuadDofTimesBas(e, p_h_newton, s_o, i, QUAD_DEFAULT, &quad_pnew);
			phgQuadDofTimesBas(e, s_o_l, s_o, i, QUAD_DEFAULT, &quad_sol);
			phgQuadDofTimesBas(e, div_uo, s_o, i, QUAD_DEFAULT, &quad_divuo);
			rhs[i] = stime * quad_qo + p_phil[0] * p_bol[0] * quad_sol - oil_cp * (quad_p - quad_pnew) - stime * quad_divuo;
		}
		for (i = 0; i < N; i++){
			phgSolverAddMatrixEntries(solver, 1, I+i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
}
static void
Solver_oileqn_so(DOF *s_o, DOF *u_o, DOF *q_o, DOF *s_o_l, DOF *p_h, DOF *p_h0, DOF *p_h_newton, DOF *phi_l, DOF *b_o_l)
{
	GRID *g = s_o->g;
	SOLVER *solver;
	DOF *div_uo, *dot_phi, *dot_bo, *b_o, *phi;
	div_uo = phgDofDivergence(u_o, NULL, NULL, NULL);
	solver = phgSolverCreate(SOLVER_PCG, s_o, NULL);
	dot_bo = phgDofNew(g, DOF_P0, 1, "dot_bo", DofNoAction);
	dot_phi = phgDofNew(g, DOF_P0, 1, "dot_phi", DofNoAction);
	update_dot_bo(p_h, dot_bo);
	update_dot_phi(p_h, dot_phi);
	phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
	update_phi(p_h, p_h0, phi);
	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
	update_bo(p_h, b_o);
	build_oileqn_so(solver, s_o, s_o_l, div_uo, q_o, p_h, p_h_newton, b_o, phi, dot_phi, dot_bo, phi_l, b_o_l);
	phgSolverSolve(solver, TRUE, s_o, NULL);
	phgDofFree(&div_uo);
	phgDofFree(&dot_bo); 
	phgDofFree(&dot_phi);
	phgDofFree(&phi);
	phgDofFree(&b_o);
	phgSolverDestroy(&solver);
}
#endif
#if USE_GAS_EQN
static void
build_gaseqn_ug(SOLVER *solver, DOF *u_g, DOF *p_h, DOF *krg, DOF *b_g, DOF *mu_g, DOF *kro, DOF *b_o, DOF *mu_o)
{
	int i, j, k;
	GRID *g = p_h->g;
	SIMPLEX *e;
	ForAllElements(g, e){
		int N = DofGetNBas(u_g, e);
		FLOAT A[N][N], rhs[N], *p_krg, *p_bg, *p_mug;
		INT I[N];
		p_krg = DofElementData(krg, e->index);
		p_bg = DofElementData(b_g, e->index);
		p_mug = DofElementData(mu_g, e->index);
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for (j = 0; j < N; j++){
				   A[i][j] = phgQuadBasDotBas(e, u_g, i, u_g, j, QUAD_DEFAULT) * p_mug[0] / (p_krg[0] * p_bg[0] * K);
			}
			rhs[i] = phgQuadDofTimesDivBas(e, p_h, u_g, i, 5);
		     /* handle bdry */
	   		if (phgDofGetElementBoundaryType(u_g, e, i) & (DIRICHLET | NEUMANN)) {
				bzero(A[i], N * sizeof(A[i][0]));
				for (j = 0; j < N; j++){
					A[j][i] = 0.;
				}
				A[i][i] = 1.;
				rhs[i] = 0.;
	   		 }
		}
		for(i = 0; i < N; i++){
			phgSolverAddMatrixEntries(solver, 1, I+i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}

}
static void
build_gaseqn_so(SOLVER *solver, DOF *s_o, DOF *s_o_l, DOF *div_uo, DOF *q_o, DOF *p_h, DOF *p_h_newton, DOF *b_o, DOF *phi, DOF *dot_phi, DOF *dot_bo, DOF *phi_l, DOF *b_o_l)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	int i,j;
	int N = s_o->type->nbas;
	FLOAT A[N][N], rhs[N];
	INT I[N];
	FLOAT *p_phi, *p_bo, *p_so, *p_dotbo, *p_dotphi, *p_phil, *p_bol;
	ForAllElements(g, e){
		p_phi = DofElementData(phi, e->index);
		p_phil = DofElementData(phi_l, e->index);
		p_bo  = DofElementData(b_o, e->index);
		p_bol = DofElementData(b_o_l, e->index);
		p_so  = DofElementData(s_o, e->index);
		p_dotbo = DofElementData(dot_bo, e->index);
		p_dotphi = DofElementData(dot_phi, e->index);
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for (j = 0; j < N; j++){
				A[i][j] =  p_bo[0] * p_phi[0] * phgQuadBasDotBas(e, s_o, i, s_o, j, QUAD_DEFAULT); 
			}
		}
		FLOAT quad_qo = 0, quad_p = 0, quad_pnew = 0, quad_sol = 0, quad_divuo = 0;
		FLOAT oil_cp = 0;
		oil_cp = p_so[0] * (p_bo[0] * p_dotphi[0] + p_phi[0] * p_dotbo[0]);
		for (i = 0; i < N; i++){
			phgQuadDofTimesBas(e, q_o, s_o, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, p_h, s_o, i, QUAD_DEFAULT, &quad_p);
			phgQuadDofTimesBas(e, p_h_newton, s_o, i, QUAD_DEFAULT, &quad_pnew);
			phgQuadDofTimesBas(e, s_o_l, s_o, i, QUAD_DEFAULT, &quad_sol);
			phgQuadDofTimesBas(e, div_uo, s_o, i, QUAD_DEFAULT, &quad_divuo);
			rhs[i] = stime * quad_qo + p_phil[0] * p_bol[0] * quad_sol - oil_cp * (quad_p - quad_pnew) - stime * quad_divuo;
		}
		for (i = 0; i < N; i++){
			phgSolverAddMatrixEntries(solver, 1, I+i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
}
static void
Solver_gaseqn_so(DOF *s_o, DOF *u_o, DOF *q_o, DOF *q_g, DOF *s_o_l, DOF *p_h, DOF *p_h0, DOF *p_h_newton, DOF *phi_l, DOF *b_o_l)
{
	GRID *g = s_o->g;
	SOLVER *solver;
	DOF *div_uo, *dot_phi, *dot_bo, *b_o, *phi;
	div_uo = phgDofDivergence(u_o, NULL, NULL, NULL);
	solver = phgSolverCreate(SOLVER_PCG, s_o, NULL);
	dot_bo = phgDofNew(g, DOF_P0, 1, "dot_bo", DofNoAction);
	dot_phi = phgDofNew(g, DOF_P0, 1, "dot_phi", DofNoAction);
	update_dot_bo(p_h, dot_bo);
	update_dot_phi(p_h, dot_phi);
	phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
	update_phi(p_h, p_h0, phi);
	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
	update_bo(p_h, b_o);
	build_oileqn_so(solver, s_o, s_o_l, div_uo, q_o, p_h, p_h_newton, b_o, phi, dot_phi, dot_bo, phi_l, b_o_l);
	phgSolverSolve(solver, TRUE, s_o, NULL);
	DOF *b_g = phgDofNew(g, DOF_P0, 1, "b_g", DofNoAction);
	update_bg(p_h, b_g);
	DOF *krg = phgDofNew(g, DOF_P0, 1, "krg", DofNoAction);
	create_krg(p_h, krg);
	DOF *kro = phgDofNew(g, DOF_P0, 1, "kro", DofNoAction);
	create_kro(p_h, kro);
	DOF *mu_o = phgDofNew(g, DOF_P0, 1, "mu_o", DofNoAction);
	DOF *mu_g = phgDofNew(g, DOF_P0, 1, "mu_g", DofNoAction);
	update_muo(p_h, mu_o);
	update_mug(p_h, mu_g);
	gas_fluidity(q_o, q_g, kro, krg, b_o, b_g, mu_o, mu_g);
	phgDofFree(&b_g);
	phgDofFree(&krg);
	phgDofFree(&kro);
	phgDofFree(&mu_o);
	phgDofFree(&mu_g);
	phgDofFree(&div_uo);
	phgDofFree(&dot_bo); 
	phgDofFree(&dot_phi);
	phgDofFree(&phi);
	phgDofFree(&b_o);
	phgSolverDestroy(&solver);
}
#endif
static void
test_grid(DOF *p_h)
{
	GRID *g = p_h->g;
	SIMPLEX *e;
	ForAllElements(g, e){
		printf("e->region_mark = %d\n", e->region_mark);
	}
	exit(1);
}
int
main(int argc, char *argv[])
{
    INT mem_max = 10240;
    char *fn = "single_level.mesh";
    GRID *g;
    char vtkfile[1000];
    static int pre_refines = 0;
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgInit(&argc, &argv);
    g = phgNewGrid(-1);

    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
	
    phgRefineAllElements(g, pre_refines); 
		if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
				(double)g->lif);

    /*The pressure variable*/
    DOF *p_h;
    p_h = phgDofNew(g, DOF_P0, 1, "p_h", DofNoAction);
    phgDofSetDataByValue(p_h, PRESSURE0);
    DOF *p0_h;
    p0_h = phgDofNew(g, DOF_P0, 1, "p0_h", DofNoAction);
    phgDofSetDataByValue(p0_h, PRESSURE0);

    DOF *s_o;
    s_o = phgDofNew(g, DOF_P0, 1, "s_o", DofNoAction);
    phgDofSetDataByValue(s_o, SO0);
   /*The velocity variable*/
    DOF *u_o;
    u_o = phgDofNew(g, DOF_RT1, 1, "u_o", DofNoAction);
    DOF *u_g;
    u_g = phgDofNew(g, DOF_RT1, 1, "u_g", DofNoAction);
    /* RHS function */
    DOF *q_g;
    q_g = phgDofNew(g, DOF_P0, 1, "q_g", DofNoAction);

    DOF *q_o;			     
    q_o = phgDofNew(g, DOF_P0, 1, "q_o", DofNoAction);
    Well_init(q_o);
	Unit_Conversion();
//	Test_PVT_uniform();
//	exit(1);
	phgPrintf("the elements is :%d\n",DofGetDataCountGlobal(p_h));
	phgPrintf("the Dof Number is :%d\n",DofGetDataCountGlobal(p_h)+DofGetDataCountGlobal(u_o));
	int flag = 0;
//	test_grid(p_h);
	while (ctime < T -1e-8 ){
		ptime = ctime;
		ctime += stime;
		phgPrintf("==================================================================\n");
		phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)ctime);
		if (ctime > T){ 
			ctime = T;
			stime = T- ptime;
			phgPrintf("current time layer : [%lf, %lf]\n",(double)ptime, (double)ctime);
		}
	 	DOF *phi_l, *b_o_l, *Rs_l, *b_g_l, *s_o_l;
		phi_l = phgDofNew(g, DOF_P0, 1, "phi_l", DofNoAction);
		update_phi(p_h, p0_h, phi_l);
	 	b_o_l = phgDofNew(g, DOF_P0, 1, "b_o_l", DofNoAction);
		update_bo(p_h, b_o_l);
	 	b_g_l = phgDofNew(g, DOF_P0, 1, "b_o_l", DofNoAction);
		update_bg(p_h, b_g_l);
		s_o_l = phgDofCopy(s_o, NULL, NULL, NULL);
	 	Rs_l = phgDofNew(g, DOF_P0, 1, "Rs_l", DofNoAction);
		update_Rs(p_h, Rs_l);
		int count = 0;
		while (TRUE){
			count++;
			DOF *p_h_newton, *s_o_newton;;
			p_h_newton = phgDofCopy(p_h, NULL, NULL, NULL);
			s_o_newton = phgDofCopy(s_o, NULL, NULL, NULL);
			/*  For other formulation 
			DOF *p_h_newton, *b_w, *phi, *dot_phibw, *kro, *krw;
		 	b_w = phgDofNew(g, DOF_P0, 1, "b_w_l", DofNoAction);
			update_bw(p_h, p0_h, b_w);
			phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
			update_phi(p_h, p0_h, phi);
			krw = phgDofNew(g, DOF_P0, 1, "krw", DofNoAction);
			create_krw(s_w, krw);
			dot_phibw = phgDofNew(g, DOF_P0, 1, "dot_phibw", DofNoAction);
			create_dot_phibw(p_h, p0_h, dot_phibw);
			*/
			Solve_Pressure(u_o, p_h, p0_h, q_o, s_o, s_o_l, q_g, p_h_newton, phi_l, b_o_l, b_g_l, Rs_l);
			Solver_oileqn_so(s_o, u_o, q_o, s_o_l, p_h, p0_h, p_h_newton, phi_l, b_o_l);
			FLOAT err = 0.0, norm = 0.0, TOL = 0.0;
			err = L2(p_h_newton, p_h, 4);
			norm = L2_norm(p_h_newton);
			TOL = err / norm;
			phgDofFree(&p_h_newton);
			FLOAT err_s = 0.0, norm_s = 0.0, TOL_s = 0.0;
			err_s = L2(s_o_newton, s_o, 4);
			norm_s = L2_norm(s_o_newton);
			TOL_s = err_s / norm_s;
			phgDofFree(&s_o_newton);
		     /*
			phgDofFree(&b_w);
			phgDofFree(&phi);
			phgDofFree(&dot_phibw);
			phgDofFree(&kro);
			phgDofFree(&krw);
			*/
			phgPrintf("err = %lf, norm = %lf,  TOL = %le\n", err, norm, TOL);
			phgPrintf("err_s = %lf, norm_s = %lf,  TOL_s = %le\n", err_s, norm_s, TOL_s);
			if ((TOL < 1e-6 )& (TOL_s < 1e-6)){
				phgPrintf("Non_ints = %d\n", count);
				break;
			}
		}
		/*Conseravation Test*/
		phgPrintf("====================Conservation Test======================\n");
		oil_conservation(p_h, p0_h, phi_l, b_o_l, s_o, s_o_l, q_o, u_o);
		gas_conservation(p_h, p0_h, phi_l, b_g_l, b_o_l, s_o, s_o_l, Rs_l, q_g);
		phgPrintf("====================End Of ConTest==========================\n");
		phgDofFree(&b_o_l);
		phgDofFree(&b_g_l);
		phgDofFree(&phi_l);
		phgDofFree(&s_o_l);
		phgDofFree(&Rs_l);
		FLOAT pwf = 0;
		Well_Pressure(p_h, &pwf);
		phgPrintf("t = %lf, Pwf = %lf\n", ctime, pwf);
		DOF *b_o, *kro, *mu_o;
	 	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
		update_bo(p_h, b_o);
	 	kro = phgDofNew(g, DOF_P0, 1, "kr_o", DofNoAction);
		create_kro(s_o, kro);
	 	mu_o = phgDofNew(g, DOF_P0, 1, "mu_o", DofNoAction);
		update_bo(p_h, mu_o);
		FLOAT pwf0 = 0;
		Peaceman_Model(p_h, b_o, kro, mu_o, &pwf0);
		phgPrintf("USE Peaceman t = %lf, Pwf = %lf\n", ctime, pwf0);
		phgDofFree(&b_o);
		phgDofFree(&kro);
		phgDofFree(&mu_o);

/*		MAP *map_p = phgMapCreate(p_h, NULL);
		VEC *y = phgMapCreateVec(map_p, 1);
 	  	phgMapDofToLocalData(map_p, 1, &p_h, y->data);
		phgVecDumpMATLAB(y, "p", "p.m");
		phgMapDestroy(&map_p);
		phgVecDestroy(&y);
		MAP *map_s = phgMapCreate(s_o, NULL);
		VEC *x = phgMapCreateVec(map_s, 1);
	 	phgMapDofToLocalData(map_s, 1, &s_o, x->data);
	  	phgVecDumpMATLAB(x, "s", "s.m");
		phgMapDestroy(&map_s);
		phgVecDestroy(&x);
*/		/*Create VTK files*/
/*		DOF *kro = phgDofNew(g, DOF_P0, 1, "kro", DofNoAction);
			create_kro(s_o, kro);
		DOF *krg = phgDofNew(g, DOF_P0, 1, "krg", DofNoAction);
			create_krg(s_o, krg);
		DOF *b_g = phgDofNew(g, DOF_P0, 1, "b_g", DofNoAction);
			update_bg(p_h, b_g);
		DOF *b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
			update_bo(p_h, b_o);
		DOF *phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
			update_phi(p_h, p0_h, phi);
		DOF *mu_o = phgDofNew(g, DOF_P0, 1, "mu_o", DofNoAction);
		DOF *mu_g = phgDofNew(g, DOF_P0, 1, "mu_g", DofNoAction);
			update_muo(p_h, mu_o);
			update_mug(p_h, mu_g);
		flag++;
		sprintf(vtkfile, "oil_gas_%03d.vtk", flag);
		phgExportVTK(g, vtkfile, p_h, s_o, q_o, q_g, u_o, kro, krg, b_g, b_o, phi, mu_o, mu_g, NULL);
		phgDofFree(&kro);
		phgDofFree(&krg);
		phgDofFree(&b_g);
		phgDofFree(&b_o);
		phgDofFree(&mu_o);
		phgDofFree(&mu_g);
		phgDofFree(&phi);
		*/
    }
    phgDofFree(&p_h);
    phgDofFree(&p0_h);
    phgDofFree(&u_o);
    phgDofFree(&u_g);
    phgDofFree(&s_o);
    phgDofFree(&q_o);
    phgDofFree(&q_g);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
/*------------------------------------------------------
 * implementation of RT by Wang Kun
 * ---------------------------------------------------*/
static const FLOAT *
RT1_bas(DOF *dof, SIMPLEX *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of basis functions */
{
    static FLOAT values[NFace][Dim];
    FLOAT (*p)[Dim] = values;
    FLOAT L0 = lambda[0], L1 = lambda[1], L2 = lambda[2], L3 = lambda[3];
    INT i, a, b, c;
    const FLOAT *J, *g1, *g2, *g3, *g4;

    J = phgGeomGetJacobian(dof->g, e);
#if 0
    FLOAT *f1, *f2, *f3, *f4;
    f1 = phgAlloc(3 * sizeof(FLOAT));
    f2 = phgAlloc(3 * sizeof(FLOAT));
    f3 = phgAlloc(3 * sizeof(FLOAT));
    f4 = phgAlloc(3 * sizeof(FLOAT));
#endif
    FLOAT f1[3], f2[3], f3[3], f4[3];

    if (no1 <= 0)
	no1 = NFace;
    assert(no0 < no1);

    switch (no0) {
	case 0:
	    a = e->verts[1];
	    b = e->verts[2];
	    c = e->verts[3];
	    if ((a < b) && (b < c)) {
		g2 = J + 4;
		g3 = J + 8;
		g4 = J + 12;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L1 * *(f2 + i) + L2 * *(f3 + i) + L3 * *(f4 + i);
		}
	    }
	    if ((a < c) && (c < b)) {
		g2 = J + 4;
		g3 = J + 12;
		g4 = J + 8;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L1 * *(f2 + i) + L3 * *(f3 + i) + L2 * *(f4 + i);
		}
	    }
	    if ((b < c) && (c < a)) {
		g2 = J + 8;
		g3 = J + 12;
		g4 = J + 4;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L2 * *(f2 + i) + L3 * *(f3 + i) + L1 * *(f4 + i);
		}
	    }
	    if ((b < a) && (a < c)) {
		g2 = J + 8;
		g3 = J + 4;
		g4 = J + 12;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L2 * *(f2 + i) + L1 * *(f3 + i) + L3 * *(f4 + i);
		}
	    }
	    if ((c < a) && (a < b)) {
		g2 = J + 12;
		g3 = J + 4;
		g4 = J + 8;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L3 * *(f2 + i) + L1 * *(f3 + i) + L2 * *(f4 + i);
		}
	    }
	    if ((c < b) && (b < a)) {
		g2 = J + 12;
		g3 = J + 8;
		g4 = J + 4;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L3 * *(f2 + i) + L2 * *(f3 + i) + L1 * *(f4 + i);
		}
	    }
	    for (i = 0; i < Dim; i++)
		(*p)[i] = (*p)[i] * 2;
	    if (no1 == 1)
		return (FLOAT *)values;
	    p++;

	case 1:
	    a = e->verts[0];
	    b = e->verts[2];
	    c = e->verts[3];
	    if ((a < b) && (b < c)) {
		g1 = J;
		g3 = J + 8;
		g4 = J + 12;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L0 * *(f1 + i) + L2 * *(f3 + i) + L3 * *(f4 + i);
		}
	    }
	    if ((a < c) && (c < b)) {
		g1 = J;
		g3 = J + 12;
		g4 = J + 8;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L0 * *(f1 + i) + L3 * *(f3 + i) + L2 * *(f4 + i);
		}
	    }
	    if ((b < c) && (c < a)) {
		g1 = J + 8;
		g3 = J + 12;
		g4 = J;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L2 * *(f1 + i) + L3 * *(f3 + i) + L0 * *(f4 + i);
		}
	    }
	    if ((b < a) && (a < c)) {
		g1 = J + 8;
		g3 = J;
		g4 = J + 12;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L2 * *(f1 + i) + L0 * *(f3 + i) + L3 * *(f4 + i);
		}
	    }
	    if ((c < b) && (b < a)) {
		g1 = J + 12;
		g3 = J + 8;
		g4 = J;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L3 * *(f1 + i) + L2 * *(f3 + i) + L0 * *(f4 + i);
		}
	    }
	    if ((c < a) && (a < b)) {
		g1 = J + 12;
		g3 = J;
		g4 = J + 8;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L3 * *(f1 + i) + L0 * *(f3 + i) + L2 * *(f4 + i);
		}
	    }
	    for (i = 0; i < Dim; i++)
		(*p)[i] = (*p)[i] * 2;
	    if (no1 == 2)
		return (FLOAT *)values;
	    p++;
	case 2:
	    a = e->verts[0];
	    b = e->verts[1];
	    c = e->verts[3];

	    if ((a < b) && (b < c)) {
		g1 = J;
		g2 = J + 4;
		g4 = J + 12;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L0 * *(f1 + i) + L1 * *(f2 + i) + L3 * *(f4 + i);
		}
	    }
	    if ((a < c) && (c < b)) {
		g1 = J;
		g2 = J + 12;
		g4 = J + 4;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L0 * *(f1 + i) + L3 * *(f2 + i) + L1 * *(f4 + i);
		}
	    }
	    if ((b < c) && (c < a)) {
		g1 = J + 4;
		g2 = J + 12;
		g4 = J;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L1 * *(f1 + i) + L3 * *(f2 + i) + L0 * *(f4 + i);
		}
	    }
	    if ((b < a) && (a < c)) {
		g1 = J + 4;
		g2 = J;
		g4 = J + 12;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L1 * *(f1 + i) + L0 * *(f2 + i) + L3 * *(f4 + i);
		}
	    }
	    if ((c < b) && (b < a)) {
		g1 = J + 12;
		g2 = J + 4;
		g4 = J;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L3 * *(f1 + i) + L1 * *(f2 + i) + L0 * *(f4 + i);
		}
	    }
	    if ((c < a) && (a < b)) {
		g1 = J + 12;
		g2 = J;
		g4 = J + 4;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L3 * *(f1 + i) + L0 * *(f2 + i) + L1 * *(f4 + i);
		}
	    }
	    for (i = 0; i < Dim; i++)
		(*p)[i] = (*p)[i] * 2;
	    if (no1 == 3)
		return (FLOAT *)values;
	    p++;
	case 3:
	    a = e->verts[0];
	    b = e->verts[1];
	    c = e->verts[2];

	    if ((a < b) && (b < c)) {
		g1 = J;
		g2 = J + 4;
		g3 = J + 8;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L0 * *(f1 + i) + L1 * *(f2 + i) + L2 * *(f3 + i);
		}
	    }
	    if ((a < c) && (c < b)) {
		g1 = J;
		g2 = J + 8;
		g3 = J + 4;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L0 * *(f1 + i) + L2 * *(f2 + i) + L1 * *(f3 + i);
		}
	    }
	    if ((b < a) && (a < c)) {
		g1 = J + 4;
		g2 = J;
		g3 = J + 8;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L1 * *(f1 + i) + L0 * *(f2 + i) + L2 * *(f3 + i);
		}
	    }
	    if ((b < c) && (c < a)) {
		g1 = J + 4;
		g2 = J + 8;
		g3 = J;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L1 * *(f1 + i) + L2 * *(f2 + i) + L0 * *(f3 + i);
		}
	    }
	    if ((c < a) && (a < b)) {
		g1 = J + 8;
		g2 = J;
		g3 = J + 4;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L2 * *(f1 + i) + L0 * *(f2 + i) + L1 * *(f3 + i);
		}
	    }
	    if ((c < b) && (b < a)) {
		g1 = J + 8;
		g2 = J + 4;
		g3 = J;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i] =
			L2 * *(f1 + i) + L1 * *(f2 + i) + L0 * *(f3 + i);
		}
	    }
	    for (i = 0; i < Dim; i++)
		(*p)[i] = (*p)[i] * 2;
	    if (no1 == 4)
		return (FLOAT *)values;
    }

//      phgFree(f1);
//      phgFree(f2);
//      phgFree(f3);
//      phgFree(f4);
    return (FLOAT *)values;
}

static const FLOAT *
RT1_grad(DOF *dof, SIMPLEX *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT values[NFace][Dim][Dim + 1];
    FLOAT (*p)[Dim][Dim + 1] = values;
    const FLOAT *J, *g1, *g2, *g3, *g4;
//      FLOAT *f1,*f2,*f3,*f4;
    INT i, j, a, b, c;
    FLOAT f1[3], f2[3], f3[3], f4[3];
    J = phgGeomGetJacobian(dof->g, e);
//      f1 = phgAlloc(3*sizeof(FLOAT));
//      f2 = phgAlloc(3*sizeof(FLOAT));
//      f3 = phgAlloc(3*sizeof(FLOAT));
//      f4 = phgAlloc(3*sizeof(FLOAT));

    if (no1 <= 0)
	no1 = NFace;
    assert(no0 < no1);

    switch (no0) {
	case 0:

	    a = e->verts[1];
	    b = e->verts[2];
	    c = e->verts[3];
	    if ((a < b) && (b < c)) {
		g2 = J + 4;
		g3 = J + 8;
		g4 = J + 12;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = 0;
		    (*p)[i][1] = *(f2 + i);
		    (*p)[i][2] = *(f3 + i);
		    (*p)[i][3] = *(f4 + i);
		}
	    }
	    if ((a < c) && (c < b)) {
		g2 = J + 4;
		g3 = J + 12;
		g4 = J + 8;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = 0;
		    (*p)[i][1] = *(f2 + i);
		    (*p)[i][2] = *(f4 + i);
		    (*p)[i][3] = *(f3 + i);
		}
	    }
	    if ((b < c) && (c < a)) {
		g2 = J + 8;
		g3 = J + 12;
		g4 = J + 4;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = 0;
		    (*p)[i][1] = *(f4 + i);
		    (*p)[i][2] = *(f2 + i);
		    (*p)[i][3] = *(f3 + i);
		}
	    }
	    if ((b < a) && (a < c)) {
		g2 = J + 8;
		g3 = J + 4;
		g4 = J + 12;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = 0;
		    (*p)[i][1] = *(f3 + i);
		    (*p)[i][2] = *(f2 + i);
		    (*p)[i][3] = *(f4 + i);
		}
	    }
	    if ((c < a) && (a < b)) {
		g2 = J + 12;
		g3 = J + 4;
		g4 = J + 8;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = 0;
		    (*p)[i][1] = *(f3 + i);
		    (*p)[i][2] = *(f4 + i);
		    (*p)[i][3] = *(f2 + i);
		}
	    }
	    if ((c < b) && (b < a)) {
		g2 = J + 12;
		g3 = J + 8;
		g4 = J + 4;
		*f2 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f2 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g2 + 2) - *(g4 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g2) - *(g4) * *(g2 + 2);
		*(f3 + 2) = *(g4) * *(g2 + 1) - *(g4 + 1) * *(g2);
		*f4 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f4 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);

		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = 0;
		    (*p)[i][1] = *(f4 + i);
		    (*p)[i][2] = *(f3 + i);
		    (*p)[i][3] = *(f2 + i);
		}
	    }
	    for (i = 0; i < Dim; i++)
		for (j = 0; j < Dim + 1; j++)
		    (*p)[i][j] = (*p)[i][j] * 2;
	    if (no1 == 1)
		return (FLOAT *)values;
	    p++;
	case 1:
	    a = e->verts[0];
	    b = e->verts[2];
	    c = e->verts[3];
	    if ((a < b) && (b < c)) {
		g1 = J;
		g3 = J + 8;
		g4 = J + 12;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f1 + i);
		    (*p)[i][1] = 0;
		    (*p)[i][2] = *(f3 + i);
		    (*p)[i][3] = *(f4 + i);
		}
	    }
	    if ((a < c) && (c < b)) {
		g1 = J;
		g3 = J + 12;
		g4 = J + 8;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f1 + i);
		    (*p)[i][1] = 0;
		    (*p)[i][2] = *(f4 + i);
		    (*p)[i][3] = *(f3 + i);
		}
	    }
	    if ((b < c) && (c < a)) {
		g1 = J + 8;
		g3 = J + 12;
		g4 = J;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f4 + i);
		    (*p)[i][1] = 0;
		    (*p)[i][2] = *(f1 + i);
		    (*p)[i][3] = *(f3 + i);
		}
	    }
	    if ((b < a) && (a < c)) {
		g1 = J + 8;
		g3 = J;
		g4 = J + 12;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f3 + i);
		    (*p)[i][1] = 0;
		    (*p)[i][2] = *(f1 + i);
		    (*p)[i][3] = *(f4 + i);
		}
	    }
	    if ((c < b) && (b < a)) {
		g1 = J + 12;
		g3 = J + 8;
		g4 = J;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f4 + i);
		    (*p)[i][1] = 0;
		    (*p)[i][2] = *(f3 + i);
		    (*p)[i][3] = *(f1 + i);
		}
	    }
	    if ((c < a) && (a < b)) {
		g1 = J + 12;
		g3 = J;
		g4 = J + 8;
		*f1 = *(g3 + 1) * *(g4 + 2) - *(g3 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g3 + 2) * *(g4) - *(g3) * *(g4 + 2);
		*(f1 + 2) = *(g3) * *(g4 + 1) - *(g3 + 1) * *(g4);
		*f3 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f3 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f3 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g3 + 2) - *(g1 + 2) * *(g3 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g3) - *(g1) * *(g3 + 2);
		*(f4 + 2) = *(g1) * *(g3 + 1) - *(g1 + 1) * *(g3);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f3 + i);
		    (*p)[i][1] = 0;
		    (*p)[i][2] = *(f4 + i);
		    (*p)[i][3] = *(f1 + i);
		}
	    }
	    for (i = 0; i < Dim; i++)
		for (j = 0; j < Dim + 1; j++)
		    (*p)[i][j] = (*p)[i][j] * 2;
	    if (no1 == 2)
		return (FLOAT *)values;
	    p++;
	case 2:

	    a = e->verts[0];
	    b = e->verts[1];
	    c = e->verts[3];

	    if ((a < b) && (b < c)) {
		g1 = J;
		g2 = J + 4;
		g4 = J + 12;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f1 + i);
		    (*p)[i][1] = *(f2 + i);
		    (*p)[i][2] = 0;
		    (*p)[i][3] = *(f4 + i);
		}
	    }
	    if ((a < c) && (c < b)) {
		g1 = J;
		g2 = J + 12;
		g4 = J + 4;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f1 + i);
		    (*p)[i][1] = *(f4 + i);
		    (*p)[i][2] = 0;
		    (*p)[i][3] = *(f2 + i);
		}
	    }
	    if ((b < c) && (c < a)) {
		g1 = J + 4;
		g2 = J + 12;
		g4 = J;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f4 + i);
		    (*p)[i][1] = *(f1 + i);
		    (*p)[i][2] = 0;
		    (*p)[i][3] = *(f2 + i);
		}
	    }
	    if ((b < a) && (a < c)) {
		g1 = J + 4;
		g2 = J;
		g4 = J + 12;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f2 + i);
		    (*p)[i][1] = *(f1 + i);
		    (*p)[i][2] = 0;
		    (*p)[i][3] = *(f4 + i);
		}
	    }
	    if ((c < b) && (b < a)) {
		g1 = J + 12;
		g2 = J + 4;
		g4 = J;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f4 + i);
		    (*p)[i][1] = *(f2 + i);
		    (*p)[i][2] = 0;
		    (*p)[i][3] = *(f1 + i);
		}
	    }
	    if ((c < a) && (a < b)) {
		g1 = J + 12;
		g2 = J;
		g4 = J + 4;
		*f1 = *(g2 + 1) * *(g4 + 2) - *(g2 + 2) * *(g4 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g4) - *(g2) * *(g4 + 2);
		*(f1 + 2) = *(g2) * *(g4 + 1) - *(g2 + 1) * *(g4);
		*f2 = *(g4 + 1) * *(g1 + 2) - *(g4 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g4 + 2) * *(g1) - *(g4) * *(g1 + 2);
		*(f2 + 2) = *(g4) * *(g1 + 1) - *(g4 + 1) * *(g1);
		*f4 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f4 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f4 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f2 + i);
		    (*p)[i][1] = *(f4 + i);
		    (*p)[i][2] = 0;
		    (*p)[i][3] = *(f1 + i);
		}
	    }
	    for (i = 0; i < Dim; i++)
		for (j = 0; j < Dim + 1; j++)
		    (*p)[i][j] = (*p)[i][j] * 2;
	    if (no1 == 3)
		return (FLOAT *)values;
	    p++;
	case 3:
	    a = e->verts[0];
	    b = e->verts[1];
	    c = e->verts[2];

	    if ((a < b) && (b < c)) {
		g1 = J;
		g2 = J + 4;
		g3 = J + 8;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f1 + i);
		    (*p)[i][1] = *(f2 + i);
		    (*p)[i][2] = *(f3 + i);
		    (*p)[i][3] = 0;
		}
	    }
	    if ((a < c) && (c < b)) {
		g1 = J;
		g2 = J + 8;
		g3 = J + 4;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);

		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f1 + i);
		    (*p)[i][1] = *(f3 + i);
		    (*p)[i][2] = *(f2 + i);
		    (*p)[i][3] = 0;
		}
	    }
	    if ((b < a) && (a < c)) {
		g1 = J + 4;
		g2 = J;
		g3 = J + 8;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f2 + i);
		    (*p)[i][1] = *(f1 + i);
		    (*p)[i][2] = *(f3 + i);
		    (*p)[i][3] = 0;
		}
	    }
	    if ((b < c) && (c < a)) {
		g1 = J + 4;
		g2 = J + 8;
		g3 = J;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f3 + i);
		    (*p)[i][1] = *(f1 + i);
		    (*p)[i][2] = *(f2 + i);
		    (*p)[i][3] = 0;
		}
	    }
	    if ((c < a) && (a < b)) {
		g1 = J + 8;
		g2 = J;
		g3 = J + 4;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f2 + i);
		    (*p)[i][1] = *(f3 + i);
		    (*p)[i][2] = *(f1 + i);
		    (*p)[i][3] = 0;
		}
	    }
	    if ((c < b) && (b < a)) {
		g1 = J + 8;
		g2 = J + 4;
		g3 = J;
		*f1 = *(g2 + 1) * *(g3 + 2) - *(g2 + 2) * *(g3 + 1);
		*(f1 + 1) = *(g2 + 2) * *(g3) - *(g2) * *(g3 + 2);
		*(f1 + 2) = *(g2) * *(g3 + 1) - *(g2 + 1) * *(g3);
		*f2 = *(g3 + 1) * *(g1 + 2) - *(g3 + 2) * *(g1 + 1);
		*(f2 + 1) = *(g3 + 2) * *(g1) - *(g3) * *(g1 + 2);
		*(f2 + 2) = *(g3) * *(g1 + 1) - *(g3 + 1) * *(g1);
		*f3 = *(g1 + 1) * *(g2 + 2) - *(g1 + 2) * *(g2 + 1);
		*(f3 + 1) = *(g1 + 2) * *(g2) - *(g1) * *(g2 + 2);
		*(f3 + 2) = *(g1) * *(g2 + 1) - *(g1 + 1) * *(g2);
		for (i = 0; i < Dim; i++) {
		    (*p)[i][0] = *(f3 + i);
		    (*p)[i][1] = *(f2 + i);
		    (*p)[i][2] = *(f1 + i);
		    (*p)[i][3] = 0;
		}
	    }
	    for (i = 0; i < Dim; i++)
		for (j = 0; j < Dim + 1; j++)
		    (*p)[i][j] = (*p)[i][j] * 2;
	    if (no1 == 4)
		return (FLOAT *)values;
    }
//      phgFree(f1);
//      phgFree(f2);
//      phgFree(f3);
//      phgFree(f4);

    return (FLOAT *)values;
}

static const FLOAT *
face_normal(GRID *g, SIMPLEX *e, int face)
{
    static FLOAT normal[3 * NFace];
    FLOAT *n = normal + 3 * face;
    FLOAT x0, y0, z0, x1, y1, z1, x2, y2, z2;
    COORD *c;
    int u0, u1, u2;
    if (IsLeaf(e)) {
	memcpy(n, phgGeomGetFaceNormal(g, e, face), 3 * sizeof(*n));
    }
    else {
	GetFaceVertices(e, face, u0, u1, u2);
	x0 = (*(c = g->verts + e->verts[u0]))[0];
	y0 = (*c)[1];
	z0 = (*c)[2];
	x1 = (*(c = g->verts + e->verts[u1]))[0] - x0;
	y1 = (*c)[1] - y0;
	z1 = (*c)[2] - z0;
	x2 = (*(c = g->verts + e->verts[u2]))[0] - x0;
	y2 = (*c)[1] - y0;
	z2 = (*c)[2] - z0;
	n[0] = y1 * z2 - y2 * z1;
	n[1] = z1 * x2 - z2 * x1;
	n[2] = x1 * y2 - x2 * y1;
	/* FIXME: normalize n? */
    }
    /* FIXME: adjust direction of n? */
    return n;
}

void
RT1_init(DOF *dof, SIMPLEX *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
{
    GRID *g = dof->g;
    const FLOAT *n = NULL;
    const FLOAT *m = NULL;
    const FLOAT *p, *w, *bas;
    FLOAT d, d0;
    COORD *c;
    FLOAT x, y, z, *pd;
    FLOAT lambda[Dim + 1];
    INT i, j, k = 0, v0, v1, v2, face, order, nn;
    DOF_TYPE *t;
    QUAD *quad;
    FLOAT values[Dim], A[NFace][NFace], b[NFace];

    i = DofTypeOrder(dof, e);
    order = i + i;
    quad = phgQuadGetQuad2D(order);
    Unused(pdofvalues);

    switch (type) {
	case FACE:
	    t = dof->type;
	    nn = t->np_face;
	    break;
	default:
	    t = NULL;
	    nn = 0;
    }

    if (funcvalues != NULL) {
	memcpy(dofvalues, funcvalues, nn * dof->dim * sizeof(*dofvalues));
	return;
    }

    if (userfunc == NULL && userfunc_lambda == NULL)
	phgError(1, "%s:%d, invalid user function!\n", __FILE__, __LINE__);

    switch (type) {
	case FACE:
	    i = DofTypeOrder(dof, e);
	    order = i + i;
	    quad = phgQuadGetQuad2D(order);

	    n = face_normal(g, e, index);
	    v0 = GetFaceVertex(index, 0);
	    v1 = GetFaceVertex(index, 1);
	    v2 = GetFaceVertex(index, 2);
	    lambda[index] = 0.;
	    p = quad->points;
	    w = quad->weights;
	    d = 0.;
	    for (i = 0; i < quad->npoints; i++) {
		m = n;
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);
		if (userfunc_lambda != NULL) {
		    userfunc_lambda(dof, e, index, lambda, values);
		}
		else {
		    x = y = z = 0.;
		    for (k = 0; k < NVert; k++) {
			c = g->verts + e->verts[k];
			x += (*c)[0] * lambda[k];
			y += (*c)[1] * lambda[k];
			z += (*c)[2] * lambda[k];
		    }
		    userfunc(x, y, z, values);
		}
		pd = values;
		d0 = 0.;
		for (k = 0; k < Dim; k++) {
		    d0 += *(pd++) * *(m++);
		}
		d += d0 * *(w++);
	    }
	    b[index] = d;
	    break;
	default:
	    break;
    }
    switch (type) {
	case FACE:
	    *dofvalues = b[index] * phgGeomGetFaceArea(g, e, index);
	    break;
	default:
	    break;
    }

}

static BYTE RT1_orders[] = { 1, 1, 1, 1 };	/* Right??? */

static DOF_TYPE DOF_RT1_ =
    { DofCache, "RT1", NULL, RT1_orders, DOF_DG0, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, RT1_init, RT1_bas,
	RT1_grad,
    NULL, FE_Hdiv,
    FALSE, FALSE, -1,		/* invariant, free_after_use, id */
    NFace, 1, -1, 3,		/* nbas, order, cont, dim */
    0, 0, 1, 0
};				/* np_vert, np_edge, np_face, np_elem */

/*------------------------------------------------------
 * quad
 * ---------------------------------------------------*/
static FLOAT
phgQuadDivBasDotBas(SIMPLEX *e, DOF *u, int m, DOF *v, int n, int order)
{
    /*
     * computes integration of
     *        (Div of 'm'-th basis function of 'u')
     *        \cdot
     *        ('m'-th basis function of 'v')
     *      on element 'e' using quadrature rule 'quad'.
     */
    int i, j, nvalues = DofTypeDim(v);
    const FLOAT *g1, *g2, *w, *lambda;
    FLOAT d, d0;
    QUAD *quad;

    assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
    // phgPrintf("dim of p_h %d\n",DofTypeDim(u));
    // phgPrintf("dim of u_h %d",DofTypeDim(v));
    if (3 * nvalues != DofTypeDim(u))
	phgError(1, "%s:%d, dimensions mismatch: grad(%s) <==> (%s)\n",
		 __FILE__, __LINE__, u->name, v->name);

    if (order < 0)
	order = BasisOrder(u, e, m) - 1 + BasisOrder(v, e, n);
    if (order < 0)
	order = 0;
    quad = phgQuadGetQuad3D(order);

//    phgPrintf("dof dim %d %d\n", DofTypeDim(u),DofTypeDim(v));
    g1 = phgQuadGetBasisGradient(e, u, m, quad);
    g2 = phgQuadGetBasisValues(e, v, n, quad);
    d = 0.;
    lambda = quad->points;
    w = quad->weights;
    for (i = 0; i < quad->npoints; i++) {
	d0 = 0.;
	for (j = 0; j < nvalues; j++) {
	    d0 += (*(g1) + (*(g1 + 4)) + (*(g1 + 8))) * (*(g2++));
	    g1 = g1 + 9;
	}
	d += d0 * (*(w++));
	lambda += Dim + 1;
    }
    return d * phgGeomGetVolume(u->g, e);
}

static FLOAT
phgQuadFaceDofDotBas_(SIMPLEX *e, int face, DOF *u, DOF_PROJ proj,
		      DOF *v, DOF_PROJ vproj, int N, int order)
{
    GRID *g = u->g;
    int i, j, k, nvalues, v0, v1, v2, dim;
    FLOAT d, d0, lambda[Dim + 1];
    FLOAT *dof, *buffer, *p0, *p1;
    const FLOAT *bas, *p, *w;
    QUAD *quad;
    const FLOAT *n = NULL;	/* out normal vector */
    DOF_TYPE *type;

    assert(!SpecialDofType(v->type));
    assert(face >= 0 && face <= 3);

    type = (DofIsHP(v) ? v->hp->info->types[v->hp->max_order] : v->type);

    if (order < 0) {
	i = DofTypeOrder(u, e);
	j = DofTypeOrder(v, e);
	if (i < 0)
	    i = j;
	order = i + j;
    }
    quad = phgQuadGetQuad2D(order);

    v0 = GetFaceVertex(face, 0);
    v1 = GetFaceVertex(face, 1);
    v2 = GetFaceVertex(face, 2);
    lambda[face] = 0.;

    assert(proj == DOF_PROJ_NONE);

    dim = DofDim(u) / (((proj == DOF_PROJ_DOT)) ? Dim : 1);

    nvalues = DofTypeDim(v) / (((vproj == DOF_PROJ_DOT)) ? Dim : 1);

    if (nvalues != dim)
	phgError(1, "%s: dimensions mismatch\n", __func__);

    buffer = phgAlloc(DofDim(u) * sizeof(*buffer));
    p = quad->points;
    w = quad->weights;

    d = 0.;
    switch (proj) {
	case DOF_PROJ_NONE:
	    for (i = 0; i < quad->npoints; i++) {
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);

		dof = phgDofEval(u, e, lambda, buffer);
		bas = (FLOAT *)type->BasFuncs(v, e, N, N + 1, lambda);
		n = phgGeomGetFaceOutNormal(g, e, face);
		switch (vproj) {
		    case DOF_PROJ_NONE:
			d0 = 0.;
			for (j = 0; j < nvalues; j++) {
			    d0 += *(bas++) * *(dof++);
			}
			d += d0 * *(w++);
			break;

		    case DOF_PROJ_DOT:
			d0 = 0.;
			for (j = 0; j < nvalues; j++) {
			    d0 +=
				(bas[0] * n[0] + bas[1] * n[1] +
				 bas[2] * n[2]) * *(dof++);
			    bas += Dim;
			}
			d += d0 * *(w++);
			break;
		    case DOF_PROJ_CROSS:
			phgError(1, " To DO ......\n");
			d0 = 0.;
			for (j = 0; j < nvalues; j++) {
			    d0 += *(bas++) * *(dof++);
			}
			d += d0 * *(w++);
		}

	    }

	    phgFree(buffer);
	    break;

	case DOF_PROJ_DOT:
	    p0 = phgAlloc(nvalues * Dim * sizeof(*p0));
	    n = phgGeomGetFaceOutNormal(g, e, face);
	    for (i = 0; i < quad->npoints; i++) {
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);

		p1 = phgDofEval(u, e, lambda, p0);
		bas = type->BasFuncs(v, e, N, N + 1, lambda);
		dof = buffer;
		for (k = 0; k < nvalues; k++) {
		    *(dof++) = n[0] * p1[0] + n[1] * p1[1] + n[2] * p1[2];
		    p1 += Dim;
		}
		d0 = 0.;
		dof = buffer;
		for (j = 0; j < nvalues; j++) {
		    d0 += *(bas++) * *(dof++);
		}
		d += d0 * *(w++);
	    }

	    phgFree(buffer);
	    phgFree(p0);
	    break;

	case DOF_PROJ_CROSS:
	    p0 = phgAlloc(nvalues * sizeof(*p0));
	    n = phgGeomGetFaceOutNormal(g, e, face);
	    for (i = 0; i < quad->npoints; i++) {
		lambda[v0] = *(p++);
		lambda[v1] = *(p++);
		lambda[v2] = *(p++);

		p1 = phgDofEval(u, e, lambda, p0);
		bas = type->BasFuncs(v, e, N, N + 1, lambda);
		dof = buffer;
		for (k = 0; k < nvalues / Dim; k++) {
		    *(dof++) = p1[1] * n[2] - p1[2] * n[1];
		    *(dof++) = p1[2] * n[0] - p1[0] * n[2];
		    *(dof++) = p1[0] * n[1] - p1[1] * n[0];
		    p1 += Dim;
		}
		d0 = 0.;
		dof = buffer;
		for (j = 0; j < nvalues; j++) {
		    d0 += *(bas++) * *(dof++);
		}
		d += d0 * *(w++);
	    }

	    phgFree(buffer);
	    phgFree(p0);
	    break;

	default:
	    phgError(1, "%s: unknown projection %d\n", __func__, proj);
    }

    return d * phgGeomGetFaceArea(u->g, e, face);
}

