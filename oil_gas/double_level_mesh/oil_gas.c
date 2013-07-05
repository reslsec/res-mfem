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
#include "rt.c"
#include "write_matlab.c"
#define USE_UZAWA 1
#define USE_BLOCK 0
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
static void
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
				mat_A[i][j] = stime * phgQuadKBasDotBas(e, u_o, i, u_o, j, QUAD_DEFAULT) / (p_kro[0] * p_bo[0] * p_muo[0]);
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
Solve_Pressure(DOF *u_o, DOF *p_h, DOF *p0_h, DOF *q_o, DOF *s_o, DOF *s_o_l, DOF *q_g, DOF *p_h_newton, DOF *phi_l, DOF *b_o_l, DOF *b_g_l, DOF *Rs_l, double *time, INT *uzawa, INT *amg)
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
		elapsed_time(g, FALSE, 0.);
		build_mat_vec(u_o, p_h, p_h_newton, s_o, s_o_l, mu_o, b_o, b_o_l, kro, phi, phi_l, dot_phi, dot_bo, Rs, Rs_l, dot_Rs, q_o, mu_g, b_g, b_g_l, krg, dot_bg, q_g, map_u, map_p, A, B, TB, C, vec_f, vec_g);
		phgPrintf("Build system:  ");
		elapsed_time(g, TRUE, 0.);
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
		phgPrintf("solve p_h:              ");
		elapsed_time(g, FALSE, 0.);
		nits_uzawa = phgUzawa(H, u_o, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
	//	nits_uzawa = uzawapcg(H, u_o, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
		*time = elapsed_time(g, TRUE, 0.);
		*uzawa = nits_uzawa;
		*amg = nits_amg;
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
		elapsed_time(g, FALSE, 0.);
		phgSolverSolve(solver, TRUE, u_o, p_h, NULL);
		phgPrintf("Solve P:   ");
		elapsed_time(g, TRUE, 0.);
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
static void
Get_Sg(DOF *s_o, DOF *s_g)
{
    GRID *g = s_o->g;
    SIMPLEX *e;
    FLOAT *p_sg, *p_so;
    ForAllElements(g, e){
	        p_sg = DofElementData(s_g, e->index);
	        p_so = DofElementData(s_o, e->index);
	        p_sg[0] = 1. - p_so[0];
	    }
}
int
main(int argc, char *argv[])
{
    INT mem_max = 10240;
    char *fn = "two_level.mesh";
	char *ff = "two_level_1T.mesh";
    GRID *g;
    char vtkfile[1000];
    static int pre_refines = 0;
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterFilename("mesh_file", "Mesh File", (char **)&fn);
    phgOptionsRegisterFloat("step", "Time step", &stime);
    phgOptionsRegisterFloat("T", "computational damain: [0, T]", &T);
    phgInit(&argc, &argv);
    g = phgNewGrid(-1);

    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
	
    phgRefineAllElements(g, pre_refines); 
#if 0
    if(!phgExportMedit(g, ff))
			            phgError(1, "can't write file \"%s\".\n",ff);
#endif
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
    DOF *s_g;
    s_g = phgDofNew(g, DOF_P0, 1, "s_g", DofNoAction);
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
	phgPrintf("the Total DOF is :%d\n",DofGetDataCountGlobal(p_h)+DofGetDataCountGlobal(u_o));
	int flag = 0;
	double total_time = 0;
	INT total_uzawa = 0, total_amg = 0, newton =0;
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
			double time = 0;
			INT uzawa = 0, amg = 0;
			DOF *p_h_newton, *s_o_newton;;
			p_h_newton = phgDofCopy(p_h, NULL, NULL, NULL);
			s_o_newton = phgDofCopy(s_o, NULL, NULL, NULL);
			Solve_Pressure(u_o, p_h, p0_h, q_o, s_o, s_o_l, q_g, p_h_newton, phi_l, b_o_l, b_g_l, Rs_l, &time, &uzawa, &amg);
			total_time += time;
			total_uzawa += uzawa;
			total_amg += amg;
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
			if ((TOL < 1e-4 )& (TOL_s < 1e-4)){
                phgPrintf("TOL_p:                    %le\n", TOL);
                phgPrintf("TOL_s:                    %le\n", TOL_s);
                phgPrintf("Non_ints:                 %d\n", count);
				newton += count;
				break;
			}
		}
#if 0
		/*Conseravation Test*/
		phgPrintf("====================Conservation Test======================\n");
		oil_conservation(p_h, p0_h, phi_l, b_o_l, s_o, s_o_l, q_o, u_o);
		gas_conservation(p_h, p0_h, phi_l, b_g_l, b_o_l, s_o, s_o_l, Rs_l, q_g);
		phgPrintf("====================End Of ConTest==========================\n");
#endif
		FLOAT pwf = 0;
		Well_Pressure(p_h, &pwf);
		phgPrintf("NO Peaceman t = %lf, Pwf = %lf\n", ctime, pwf);
		DOF *b_o, *kro, *mu_o;
	 	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
		update_bo(p_h, b_o);
	 	kro = phgDofNew(g, DOF_P0, 1, "kr_o", DofNoAction);
		create_kro(s_o, kro);
	 	mu_o = phgDofNew(g, DOF_P0, 1, "mu_o", DofNoAction);
		update_muo(p_h, mu_o);
		FLOAT pwf0 = 0;
		PeacemanModel(p_h, b_o, kro, mu_o, &pwf0);
//		Peaceman_Model(p_h, b_o, kro, mu_o, &pwf0);
		phgPrintf("USE Peaceman t = %lf, Pwf = %lf\n", ctime, pwf0);
		phgPrintf("Total time  :--------------------%lf\n", total_time);
		phgPrintf("Total uzawa :--------------------%d\n", total_uzawa);
		phgPrintf("Total  amg  :--------------------%d\n", total_amg);
		phgPrintf("Total newton:--------------------%d\n", newton);
		phgDofFree(&b_o);
		phgDofFree(&kro);
		phgDofFree(&mu_o);
#if 0
		MAP *map_p = phgMapCreate(p_h, NULL);
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
#endif
		phgDofFree(&b_o_l);
		phgDofFree(&b_g_l);
		phgDofFree(&phi_l);
		phgDofFree(&s_o_l);
		phgDofFree(&Rs_l);
#if 1
		/*Create VTK files*/
		flag++;
		sprintf(vtkfile, "oil_gas_%03d.vtk", flag);
		phgExportVTK(g, vtkfile, p_h, s_o, q_o, q_g, u_o, NULL);
#endif
    }
	Get_Sg(s_o, s_g);
#	if 0
	char *fso = "fso", *fsg="fsg", *fp="fp";
	MATLAB_Draw3D(fso, s_o, 100, 100, 1);
	MATLAB_Draw3D(fsg, s_g, 100, 100, 1);
	MATLAB_Draw3D(fp,  p_h, 100, 100, 1);
	sprintf(vtkfile, "oil_gas_%s.vtk", "finaltime");
	phgExportVTK(g, vtkfile, p_h, s_o, s_g, q_o, q_g, u_o, u_g, NULL);
#endif
    phgDofFree(&p_h);
    phgDofFree(&p0_h);
    phgDofFree(&u_o);
    phgDofFree(&u_g);
    phgDofFree(&s_o);
    phgDofFree(&s_g);
    phgDofFree(&q_o);
    phgDofFree(&q_g);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
