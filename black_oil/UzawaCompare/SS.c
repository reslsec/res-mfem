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
#define SS 1
#define IMPES 0
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
build_mat_vec(DOF *dp, DOF *d_so, DOF *d_sw, DOF *dot_muo, DOF *dot_muw, DOF *dot_mug, DOF *dso_kro, DOF *dsw_kro, DOF *dsw_krw, DOF *dso_krg, DOF *dsw_krg, DOF *u_o, DOF *p_h, DOF *p_h_newton, DOF *s_o, DOF *s_o_l, DOF *mu_o, DOF *b_o, DOF *b_o_l, DOF *kro, DOF *dot_bo, DOF *q_o, DOF *s_w, DOF *s_w_l, DOF *mu_w, DOF *b_w, DOF *b_w_l, DOF *krw, DOF *dot_bw, DOF *q_w, DOF *phi, DOF *phi_l, DOF *dot_phi, DOF *Rs, DOF *Rs_l, DOF *dot_Rs, DOF *mu_g, DOF *b_g, DOF *b_g_l, DOF *krg, DOF *dot_bg, DOF *q_g, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g)
{
	GRID *g = u_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_sol, *p_bo, *p_bol, *p_kro, *p_phi, *p_phil, *p_dotbo, *p_Rs, *p_dotRs, *p_bg, *p_bgl, *p_krg, *p_dotbg, *p_muo, *p_mug, *p_dotphi, *p_Rsl, *p_sw, *p_swl, *p_bw, *p_bwl, *p_krw, *p_dotbw, *p_muw;
	FLOAT *p_dotmuo, *p_dotmuw, *p_dotmug, *p_dsokro, *p_dswkro, *p_dswkrw, *p_dsokrg, *p_dswkrg, *p_dp, *p_dso, *p_dsw;
	INT N = u_o->type->nbas * u_o->dim;
	INT M = dp->type->nbas * dp->dim;
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
		p_muw = DofElementData(mu_w, e->index);
		p_mug = DofElementData(mu_g, e->index);
		p_sw = DofElementData(s_w, e->index);
		p_swl = DofElementData(s_w_l, e->index);
		p_bw = DofElementData(b_w, e->index);
		p_bwl = DofElementData(b_w_l, e->index);
		p_krw = DofElementData(krw, e->index);
		p_dotbw = DofElementData(dot_bw, e->index);
		/* Add Dofs For SS */
		p_dp = DofElementData(dp, e->index);
		p_dso = DofElementData(d_so, e->index);
		p_dsw = DofElementData(d_sw, e->index);
		p_dotmuo = DofElementData(dot_muo, e->index);
		p_dotmuw = DofElementData(dot_muw, e->index);
		p_dotmug = DofElementData(dot_mug, e->index);
		p_dsokro = DofElementData(dso_kro, e->index);
		p_dswkro = DofElementData(dsw_kro, e->index);
		p_dswkrw = DofElementData(dsw_krw, e->index);
		p_dsokrg = DofElementData(dso_krg, e->index);
		p_dswkrg = DofElementData(dsw_krg, e->index);
		/*Create inverse matrix*/
		FLOAT oil_so = 0., wat_sw = 0., gas_so = 0., gas_sw = 0.;
		for (i = 0; i < M; i++){
			for (j = 0; j < M; j++){
				oil_so = p_phi[0] * p_bo[0] * phgQuadBasDotBas(e, d_so, i, d_so, j, QUAD_DEFAULT);
				wat_sw = p_phi[0] * p_bw[0] * phgQuadBasDotBas(e, d_sw, i, d_sw, j, QUAD_DEFAULT);
				gas_so = p_phi[0] * (p_bg[0] - p_Rs[0] * p_bo[0]) * phgQuadBasDotBas(e, d_so, i, d_so, j, QUAD_DEFAULT);
				gas_sw = p_phi[0] * p_bg[0] * phgQuadBasDotBas(e, d_sw, i, d_sw, j, QUAD_DEFAULT);
			}
		}
#if SS
		FLOAT alpha_o = 0., alpha_w = 0., alpha_g = 0.;
		alpha_o = K * p_kro[0] * p_bo[0] * p_muo[0] + K * p_kro[0] * (p_bo[0] * p_dotmuo[0] + p_dotbo[0] * p_muo[0]) * p_dp[0] \ 
				+ K * p_bo[0] * p_muo[0] * (p_dsokro[0] * p_dso[0] + p_dswkro[0] * p_dsw[0]);
		alpha_w = K * p_krw[0] * p_bw[0] * p_muw[0] + K * p_krw[0] * (p_bw[0] * p_dotmuw[0] + p_dotbw[0] * p_muw[0]) * p_dp[0] \ 
				+ K * p_bw[0] * p_muw[0] * p_dswkrw[0] * p_dsw[0];
		alpha_g = K * p_krg[0] * p_bg[0] * p_mug[0] + K * p_krg[0] * (p_bg[0] * p_dotmug[0] + p_dotbg[0] * p_mug[0]) * p_dp[0] \ 
				+ K * p_bg[0] * p_mug[0] * (p_dsokrg[0] * p_dso[0] + p_dswkrg[0] * p_dsw[0]);
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				mat_A[i][j] = stime * phgQuadBasDotBas(e, u_o, i, u_o, j, QUAD_DEFAULT) / alpha_o;
			}
		}
		FLOAT beta = 0;
		beta = gas_so / oil_so + gas_sw * alpha_w / (wat_sw * alpha_o) + (p_Rs[0] + alpha_g / alpha_o);
#endif
#if IMPES
		FLOAT T_o = 0, T_w = 0, T_g = 0;
		T_o = K * p_kro[0] * p_bo[0] * p_muo[0];
		T_w = K * p_krw[0] * p_bw[0] * p_muw[0];
		T_g = K * p_krg[0] * p_bg[0] * p_mug[0];
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				mat_A[i][j] = stime * phgQuadBasDotBas(e, u_o, i, u_o, j, QUAD_DEFAULT) / T_o;
			}
		}
		FLOAT beta = 0;
		beta = gas_so / oil_so + gas_sw * T_w / (wat_sw * T_o) + (p_Rs[0] + T_g / T_o);
#endif
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				mat_TB[i][j] = -stime * phgQuadDivBasDotBas(e, u_o, i, dp, j, QUAD_DEFAULT);
			}
		}
		FLOAT quad = 0., oil_cp = 0., wat_cp = 0., gas_cp = 0.;
		oil_cp = p_so[0] * (p_bo[0] * p_dotphi[0] + p_phi[0] * p_dotbo[0]);
		wat_cp = p_sw[0] * (p_bw[0] * p_dotphi[0] + p_phi[0] * p_dotbw[0]);
		gas_cp = p_so[0] * p_phi[0] * p_bo[0] * p_dotRs[0] + p_Rs[0] * oil_cp + (1. - p_so[0] - p_sw[0]) * (p_phi[0] * p_dotbg[0] + p_bg[0] * p_dotphi[0]);
		for (i = 0; i < M; i++){
			for (j = 0; j < M; j++){
				quad = phgQuadBasDotBas(e, dp, i, dp, j, QUAD_DEFAULT);
				mat_C[i][j] = (gas_so * oil_cp / oil_so + gas_sw * wat_cp / wat_sw + gas_cp) * quad / beta;
			}
		}
		/*Create rhs*/
		FLOAT quad_qo = 0;
		FLOAT quad_qw = 0.;
		FLOAT quad_qg = 0.;
		FLOAT quad_phi = 0, quad_phil = 0;
		FLOAT rhs_oil = 0, rhs_wat = 0., rhs_gas = 0;
		for (i = 0; i < M; i++){
			phgQuadDofTimesBas(e, q_o, dp, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, q_w, dp, i, QUAD_DEFAULT, &quad_qw);
			phgQuadDofTimesBas(e, q_g, dp, i, QUAD_DEFAULT, &quad_qg);
			phgQuadDofTimesBas(e, phi, dp, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, phi_l, dp, i, QUAD_DEFAULT, &quad_phil);
			rhs_oil = -stime * quad_qo  + p_so[0] * p_bo[0] * quad_phi - p_sol[0] * p_bol[0] * quad_phil;
			rhs_wat = -stime * quad_qw  + p_sw[0] * p_bw[0] * quad_phi - p_swl[0] * p_bwl[0] * quad_phil;
			rhs_gas = -stime * quad_qg  + (p_bg[0] * (1. - p_so[0] - p_sw[0]) + p_Rs[0] * p_so[0] * p_bo[0]) * quad_phi \
					  - (p_Rsl[0] * p_bol[0] * p_sol[0] + p_bgl[0] * (1. - p_sol[0] - p_swl[0])) * quad_phil;
			rhs_g[i] = (gas_so * rhs_oil / oil_so + gas_sw * rhs_wat / wat_sw + rhs_gas) / beta;
		}
		for (i = 0; i < N; i++){
			rhs_f[i] = stime * phgQuadDofTimesDivBas(e, p_h, u_o, i, QUAD_DEFAULT);
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
Update_fluidity(DOF *q_o, DOF *q_w, DOF *q_g, DOF *dp, DOF *d_so, DOF *d_sw, DOF *b_o, DOF *dot_bo, DOF *mu_o, DOF *dot_muo, DOF *kro, DOF *dso_kro, DOF *dsw_kro, DOF *b_w, DOF *dot_bw, DOF *mu_w, DOF *dot_muw, DOF *krw, DOF *dsw_krw, DOF *b_g, DOF *dot_bg, DOF *mu_g, DOF *dot_mug, DOF *krg, DOF *dso_krg, DOF *dsw_krg, DOF *Rs)
{
	GRID *g = q_o->g;
	SIMPLEX *e;
	FLOAT *p_qo, *p_qw, *p_qg, *p_kro, *p_krw, *p_krg, *p_bo, *p_bw, *p_bg, *p_muo, *p_muw, *p_mug, *p_Rs, *p_dotbo, *p_dotbw, *p_dotbg;
	FLOAT *p_dotmuo, *p_dotmuw, *p_dotmug, *p_dsokro, *p_dswkro, *p_dswkrw, *p_dsokrg, *p_dswkrg, *p_dp, *p_dso, *p_dsw;
	ForAllElements(g, e){
		p_qw = DofElementData(q_w, e->index);
		p_qo = DofElementData(q_o, e->index);
		p_qg = DofElementData(q_g, e->index);
		p_kro = DofElementData(kro, e->index);
		p_krw = DofElementData(krw, e->index);
		p_krg = DofElementData(krg, e->index);
		p_bo = DofElementData(b_o, e->index);
		p_bw = DofElementData(b_w, e->index);
		p_bg = DofElementData(b_g, e->index);
		p_dotbo = DofElementData(dot_bo, e->index);
		p_dotbw = DofElementData(dot_bw, e->index);
		p_dotbg = DofElementData(dot_bg, e->index);
		p_muo = DofElementData(mu_o, e->index);
		p_muw = DofElementData(mu_w, e->index);
		p_mug = DofElementData(mu_g, e->index);
		p_Rs = DofElementData(Rs, e->index);
		/* Add Dofs For SS */
		p_dp = DofElementData(dp, e->index);
		p_dso = DofElementData(d_so, e->index);
		p_dsw = DofElementData(d_sw, e->index);
		p_dotmuo = DofElementData(dot_muo, e->index);
		p_dotmuw = DofElementData(dot_muw, e->index);
		p_dotmug = DofElementData(dot_mug, e->index);
		p_dsokro = DofElementData(dso_kro, e->index);
		p_dswkro = DofElementData(dsw_kro, e->index);
		p_dswkrw = DofElementData(dsw_krw, e->index);
		p_dsokrg = DofElementData(dso_krg, e->index);
		p_dswkrg = DofElementData(dsw_krg, e->index);
#if SS
		FLOAT alpha_o = 0., alpha_w = 0., alpha_g = 0.;
		alpha_o = K * p_kro[0] * p_bo[0] * p_muo[0] + K * p_kro[0] * (p_bo[0] * p_dotmuo[0] + p_dotbo[0] * p_muo[0]) * p_dp[0] \ 
				+ K * p_bo[0] * p_muo[0] * (p_dsokro[0] * p_dso[0] + p_dswkro[0] * p_dsw[0]);
		alpha_w = K * p_krw[0] * p_bw[0] * p_muw[0] + K * p_krw[0] * (p_bw[0] * p_dotmuw[0] + p_dotbw[0] * p_muw[0]) * p_dp[0] \ 
				+ K * p_bw[0] * p_muw[0] * p_dswkrw[0] * p_dsw[0];
		alpha_g = K * p_krg[0] * p_bg[0] * p_mug[0] + K * p_krg[0] * (p_bg[0] * p_dotmug[0] + p_dotbg[0] * p_mug[0]) * p_dp[0] \ 
				+ K * p_bg[0] * p_mug[0] * (p_dsokrg[0] * p_dso[0] + p_dswkrg[0] * p_dsw[0]);
		p_qg[0] = (alpha_g + p_Rs[0] * alpha_o) / alpha_o * p_qo[0];
		p_qw[0] = alpha_w / alpha_o * p_qo[0];
#endif
#if IMPES
		FLOAT T_o = 0, T_w = 0, T_g = 0;
		T_o = K * p_kro[0] * p_bo[0] * p_muo[0];
		T_w = K * p_krw[0] * p_bw[0] * p_muw[0];
		T_g = K * p_krg[0] * p_bg[0] * p_mug[0];
		p_qw[0] = T_w / T_o * p_qo[0];
		p_qg[0] = (p_Rs[0] + T_g  / T_o) * p_qo[0];
#endif
	}
}
static void
Solve_Pressure(DOF *dp, DOF *d_so, DOF *d_sw, DOF *u_o, DOF *p_h, DOF *p0_h, DOF *p_h_newton, DOF *q_o, DOF *s_o, DOF *s_o_l, DOF *b_o_l, DOF *q_w, DOF *s_w, DOF *s_w_l, DOF *b_w_l, DOF *q_g, DOF *phi_l, DOF *b_g_l, DOF *Rs_l, DOF *s_g)
{
		GRID *g = u_o->g;	
	     MAT *A, *B, *TB, *C;
	     VEC *vec_f, *vec_g;
	     MAP *map_u, *map_p;
     	/* The parameter DOF */
 	  	DOF *phi, *b_o, *b_w, *b_g, *kro, *krw, *krg, *Rs, *mu_o, *mu_w, *mu_g, *dot_bo, *dot_bw, *dot_bg, *dot_Rs, *dot_phi;
		DOF *dot_muo, *dot_muw, *dot_mug, *dso_kro, *dsw_kro, *dsw_krw, *dso_krg, *dsw_krg;
		/*Create MAP for Mat and Vec*/
	    map_p = phgMapCreate(dp, NULL);
		map_u = phgMapCreate(u_o, NULL);
		A = phgMapCreateMat(map_u, map_u);
		B = phgMapCreateMat(map_p, map_u);
		TB = phgMapCreateMat(map_u, map_p);
		C = phgMapCreateMat(map_p, map_p);
		vec_f = phgMapCreateVec(map_u, 1);
		vec_g = phgMapCreateVec(map_p, 1);
			
		phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
		update_phi(p_h, p0_h, phi);
 		b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
		update_bo(p_h, b_o);
 		mu_o = phgDofNew(g, DOF_P0, 1, "mu_o", DofNoAction);
		update_muo(p_h, mu_o);
 		b_w = phgDofNew(g, DOF_P0, 1, "b_w", DofNoAction);
		update_bw(p_h, b_w);
 		mu_w = phgDofNew(g, DOF_P0, 1, "mu_w", DofNoAction);
		update_muw(p_h, mu_w);
 		b_g = phgDofNew(g, DOF_P0, 1, "b_g", DofNoAction);
		update_bg(p_h, b_g);
 		mu_g = phgDofNew(g, DOF_P0, 1, "mu_g", DofNoAction);
		update_mug(p_h, mu_g);
 		Rs = phgDofNew(g, DOF_P0, 1, "Rs", DofNoAction);
		update_Rs(p_h, Rs);
		/*kro, krw, krg not sure*/
		krw = phgDofNew(g, DOF_P0, 1, "krw", DofNoAction);
		create_krw(s_w, krw);
		krg = phgDofNew(g, DOF_P0, 1, "krg", DofNoAction);
		create_krg(s_g, krg);
		kro = phgDofNew(g, DOF_P0, 1, "kro", DofNoAction);
		create_kro(s_w, s_g, kro);
		dot_bo = phgDofNew(g, DOF_P0, 1, "dot_bo", DofNoAction);
		dot_bw = phgDofNew(g, DOF_P0, 1, "dot_bw", DofNoAction);
		dot_bg = phgDofNew(g, DOF_P0, 1, "dot_bg", DofNoAction);
		dot_Rs = phgDofNew(g, DOF_P0, 1, "dot_Rs", DofNoAction);
		dot_phi = phgDofNew(g, DOF_P0, 1, "dot_phi", DofNoAction);
		update_dot_bo(p_h, dot_bo);
		update_dot_bw(p_h, dot_bw);
		update_dot_bg(p_h, dot_bg);
		update_dot_Rs(p_h, dot_Rs);
		update_dot_phi(p_h, dot_phi);
		/*Add DOFs For SS Methods*/
		dot_muo = phgDofNew(g, DOF_P0, 1, "dot_muo", DofNoAction);
		dot_muw = phgDofNew(g, DOF_P0, 1, "dot_muw", DofNoAction);
		dot_mug = phgDofNew(g, DOF_P0, 1, "dot_mug", DofNoAction);
		dso_kro = phgDofNew(g, DOF_P0, 1, "dso_kro", DofNoAction);
		dsw_kro = phgDofNew(g, DOF_P0, 1, "dsw_kro", DofNoAction);
		dsw_krw = phgDofNew(g, DOF_P0, 1, "dsw_krw", DofNoAction);
		dso_krg = phgDofNew(g, DOF_P0, 1, "dso_krg", DofNoAction);
		dsw_krg = phgDofNew(g, DOF_P0, 1, "dsw_krg", DofNoAction);
		update_dot_muo(p_h, dot_muo);
		update_dot_muw(p_h, dot_muw);
		update_dot_mug(p_h, dot_mug);
		update_dot_kro(s_o, s_w, dso_kro, dsw_kro);
		update_dot_krw(s_w, dsw_krw);
		update_dot_krg(s_o, s_w, dso_krg, dsw_krg);
		/*--------------------------------------------------------*/
		Update_fluidity(q_o, q_w, q_g, dp, d_so, d_sw, b_o, dot_bo, mu_o, dot_muo, kro, dso_kro, dsw_kro, \
						b_w, dot_bw, mu_w, dot_muw, krw, dsw_krw, b_g, dot_bg, mu_g, dot_mug, krg, dso_krg, dsw_krg, Rs);
		build_mat_vec(dp, d_so, d_sw, dot_muo, dot_muw, dot_mug, dso_kro, dsw_kro, dsw_krw, dso_krg, dsw_krg, \
						u_o, p_h, p_h_newton, s_o, s_o_l, mu_o, b_o, b_o_l, kro, dot_bo, q_o, s_w, s_w_l, mu_w, \
						b_w, b_w_l, krw, dot_bw, q_w, phi, phi_l, dot_phi, Rs, Rs_l, dot_Rs, mu_g, b_g, b_g_l, krg, \
						dot_bg, q_g, map_u, map_p, A, B, TB, C, vec_f, vec_g);
		phgDofFree(&phi);
		phgDofFree(&mu_o);
		phgDofFree(&b_o);
		phgDofFree(&kro);
		phgDofFree(&mu_w);
		phgDofFree(&b_w);
		phgDofFree(&krw);
		phgDofFree(&mu_g);
		phgDofFree(&b_g);
		phgDofFree(&krg);
		phgDofFree(&Rs);
		phgDofFree(&dot_bo);
		phgDofFree(&dot_bw);
		phgDofFree(&dot_bg);
		phgDofFree(&dot_Rs);
		phgDofFree(&dot_phi);
		/* Add DOFs For SS Methods*/
		phgDofFree(&dot_muo);
		phgDofFree(&dot_muw);
		phgDofFree(&dot_mug);
		phgDofFree(&dso_kro);
		phgDofFree(&dsw_kro);
		phgDofFree(&dsw_krw);
		phgDofFree(&dso_krg);
		phgDofFree(&dsw_krg);
		/*-----------------------------------------*/
#if USE_UZAWA
		int nits_amg = 0;
        int nits_uzawa = 0;
		double time_amg = 0;
		elapsed_time(g, FALSE, 0.);
		DOF *B_data = DivRT(dp, u_o);
        MAT *H = BTAB(A, B_data, dp, u_o);
      	phgDofFree(&B_data);
		phgPrintf("Assemble H:             ");
		elapsed_time(g, TRUE, 0.);
		phgPrintf("solve p_h:              ");
	//	nits_uzawa = phgUzawa(H, u_w, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg, &time_amg);
		elapsed_time(g, FALSE, 0.);
		nits_uzawa = uzawapcg(H, u_o, dp, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
		elapsed_time(g, TRUE, 0.);
		phgPrintf("MAx iter of AMG---------%d\n", nits_amg );
		phgPrintf("Nits: uzawa:----------------%d\n", nits_uzawa);
		phgMatDestroy(&H);
#endif
#if USE_BLOCK
		SOLVER *solver;
		MAT *pmat[4];
		FLOAT coef[4];
		solver = phgSolverCreate(SOLVER_MUMPS, u_o, dp, NULL);
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
		phgSolverSolve(solver, TRUE, u_o, dp, NULL);
		phgSolverDestroy(&solver);
		elapsed_time(g, FALSE, 0.);
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
build_wateqn_uw(SOLVER *solver, DOF *u_w, DOF *p_h, DOF *krw, DOF *b_w, DOF *mu_w, DOF *dp, DOF *d_sw, DOF *dot_bw, DOF *dot_muw, DOF *dsw_krw)
{
	int i, j, k;
	GRID *g = p_h->g;
	SIMPLEX *e;
	ForAllElements(g, e){
		int N = DofGetNBas(u_w, e);
		FLOAT A[N][N], rhs[N], *p_krw, *p_bw, *p_muw, *p_dp, *p_dsw, *p_dotbw, *p_dotmuw, *p_dswkrw;
		INT I[N];
		p_dp = DofElementData(dp, e->index);
		p_dsw = DofElementData(d_sw, e->index);
		p_krw = DofElementData(krw, e->index);
		p_dswkrw = DofElementData(dsw_krw, e->index);
		p_bw = DofElementData(b_w, e->index);
		p_dotbw = DofElementData(dot_bw, e->index);
		p_muw = DofElementData(mu_w, e->index);
		p_dotmuw = DofElementData(dot_muw, e->index);
		FLOAT alpha_w = K * p_krw[0] * p_bw[0] * p_muw[0] + K * p_krw[0] * (p_bw[0] * p_dotmuw[0] + p_dotbw[0] * p_muw[0]) * p_dp[0] + K * p_bw[0] * p_muw[0] * p_dswkrw[0] * p_dsw[0];
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for (j = 0; j < N; j++){
				   A[i][j] = phgQuadBasDotBas(e, u_w, i, u_w, j, QUAD_DEFAULT) / alpha_w;
			}
			rhs[i] = phgQuadDofTimesDivBas(e, p_h, u_w, i, 5);
		     /* handle bdry */
	   		if (phgDofGetElementBoundaryType(u_w, e, i) & (DIRICHLET | NEUMANN)) {
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
Get_Uw(DOF *dp, DOF *d_so, DOF *d_sw, DOF *u_w, DOF *u_o, DOF *b_w, DOF *dot_bw, DOF *mu_w, DOF *dot_muw, DOF *krw, DOF *dsw_krw, DOF *b_o, DOF *dot_bo, DOF *mu_o, DOF *dot_muo, DOF *kro, DOF *dso_kro, DOF *dsw_kro)
{
	GRID *g = u_w->g;
	SIMPLEX *e;
	INT nbas = u_w->type->nbas;
	int i; 
	FLOAT *p_bw, *p_dotbw, *p_bo, *p_dotbo, *p_muw, *p_muo, *p_kro, *p_krw, *p_dotmuo, *p_dotmuw, *p_dsokro, *p_dswkro, *p_dswkrw, *p_dp, *p_dsw, *p_dso, *p_uo, *p_uw;
	ForAllElements(g, e){
		p_uo = DofElementData(u_o, e->index);
		p_uw = DofElementData(u_w, e->index);
		p_kro = DofElementData(kro, e->index);
		p_krw = DofElementData(krw, e->index);
		p_bo = DofElementData(b_o, e->index);
		p_bw = DofElementData(b_w, e->index);
		p_dotbo = DofElementData(dot_bo, e->index);
		p_dotbw = DofElementData(dot_bw, e->index);
		p_muo = DofElementData(mu_o, e->index);
		p_muw = DofElementData(mu_w, e->index);
		/* Add Dofs For SS */
		p_dp = DofElementData(dp, e->index);
		p_dso = DofElementData(d_so, e->index);
		p_dsw = DofElementData(d_sw, e->index);
		p_dotmuo = DofElementData(dot_muo, e->index);
		p_dotmuw = DofElementData(dot_muw, e->index);
		p_dsokro = DofElementData(dso_kro, e->index);
		p_dswkro = DofElementData(dsw_kro, e->index);
		p_dswkrw = DofElementData(dsw_krw, e->index);
#if SS
		FLOAT alpha_o = 0., alpha_w = 0.;
		alpha_o = K * p_kro[0] * p_bo[0] * p_muo[0] + K * p_kro[0] * (p_bo[0] * p_dotmuo[0] + p_dotbo[0] * p_muo[0]) * p_dp[0] + K * p_bo[0] * p_muo[0] * (p_dsokro[0] * p_dso[0] + p_dswkro[0] * p_dsw[0]);
		alpha_w = K * p_krw[0] * p_bw[0] * p_muw[0] + K * p_krw[0] * (p_bw[0] * p_dotmuw[0] + p_dotbw[0] * p_muw[0]) * p_dp[0] + K * p_bw[0] * p_muw[0] * p_dswkrw[0] * p_dsw[0];
		for (i = 0; i < nbas; i++){
			p_uw[i] = alpha_w / alpha_o * p_uo[i];
		}
#endif
#if IMPES
		FLOAT T_o = 0, T_w = 0, T_g = 0;
		T_o = K * p_kro[0] * p_bo[0] * p_muo[0];
		T_w = K * p_krw[0] * p_bw[0] * p_muw[0];
		for (i = 0; i < nbas; i++){
			p_uw[i] = T_w / T_o * p_uo[i];
		}
#endif
	}
}
static void
build_wateqn_sw(SOLVER *solver, DOF *dp, DOF *d_sw, DOF *s_w, DOF *s_w_l, DOF *div_uw, DOF *q_w, DOF *p_h, DOF *p_h_newton, DOF *b_w, DOF *phi, DOF *dot_phi, DOF *dot_bw, DOF *phi_l, DOF *b_w_l)
{
	GRID *g = d_sw->g;
	SIMPLEX *e;
	int i,j;
	int N = d_sw->type->nbas;
	FLOAT A[N][N], rhs[N];
	INT I[N];
	FLOAT *p_phi, *p_bw, *p_sw, *p_dotbw, *p_dotphi, *p_swl, *p_bwl;
	ForAllElements(g, e){
		p_phi = DofElementData(phi, e->index);
		p_bw  = DofElementData(b_w, e->index);
		p_bwl = DofElementData(b_w_l, e->index);
		p_sw  = DofElementData(s_w, e->index);
		p_swl  = DofElementData(s_w_l, e->index);
		p_dotbw = DofElementData(dot_bw, e->index);
		p_dotphi = DofElementData(dot_phi, e->index);
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for (j = 0; j < N; j++){
				A[i][j] =  p_bw[0] * p_phi[0] * phgQuadBasDotBas(e, d_sw, i, d_sw, j, QUAD_DEFAULT); 
			}
		}
		FLOAT quad_qw = 0, quad_p = 0, quad_phi = 0., quad_phil = 0, quad_divuw = 0;
		FLOAT wat_cp = 0;
		wat_cp = p_sw[0] * (p_bw[0] * p_dotphi[0] + p_phi[0] * p_dotbw[0]);
		for (i = 0; i < N; i++){
			phgQuadDofTimesBas(e, q_w, d_sw, i, QUAD_DEFAULT, &quad_qw);
			phgQuadDofTimesBas(e, dp, d_sw, i, QUAD_DEFAULT, &quad_p);
			phgQuadDofTimesBas(e, phi_l, d_sw, i, QUAD_DEFAULT, &quad_phil);
			phgQuadDofTimesBas(e, phi, d_sw, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, div_uw, d_sw, i, QUAD_DEFAULT, &quad_divuw);
			rhs[i] = stime * quad_qw + p_swl[0] * p_bwl[0] * quad_phil - wat_cp * quad_p  -  p_sw[0] * p_bw[0] * quad_phi - stime * quad_divuw;
		}
		for (i = 0; i < N; i++){
			phgSolverAddMatrixEntries(solver, 1, I+i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
}
static void
Solver_wateqn_sw(DOF *dp, DOF *d_so, DOF *d_sw, DOF *s_w, DOF *s_o, DOF *s_g, DOF *u_w, DOF *u_o, DOF *q_w, DOF *s_w_l, DOF *p_h, DOF *p_h0, DOF *p_h_newton, DOF *phi_l, DOF *b_w_l)
{
	GRID *g = d_sw->g;
	SOLVER *solver;// *solver_uw;
	DOF *div_uw, *dot_phi, *dot_bw, *b_w, *phi, *mu_w, *krw, *kro, *b_o, *mu_o, *dot_bo;
	dot_bw = phgDofNew(g, DOF_P0, 1, "dot_bw", DofNoAction);
	dot_phi = phgDofNew(g, DOF_P0, 1, "dot_phi", DofNoAction);
	update_dot_bw(p_h, dot_bw);
	update_dot_phi(p_h, dot_phi);
	phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
	update_phi(p_h, p_h0, phi);
	b_w = phgDofNew(g, DOF_P0, 1, "b_w", DofNoAction);
	update_bw(p_h, b_w);
	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
	update_bo(p_h, b_o);
	dot_bo = phgDofNew(g, DOF_P0, 1, "dot_bo", DofNoAction);
	update_dot_bo(p_h, dot_bo);
	mu_w = phgDofNew(g, DOF_P0, 1, "mu_w", DofNoAction);
	update_muw(p_h, mu_w);
	mu_o = phgDofNew(g, DOF_P0, 1, "mu_o", DofNoAction);
	update_muo(p_h, mu_o);
	krw = phgDofNew(g, DOF_P0, 1, "krw", DofNoAction);
	create_krw(s_w, krw);
	kro = phgDofNew(g, DOF_P0, 1, "kro", DofNoAction);
	create_kro(s_w, s_g, kro);
		/*Add DOFs For SS Methods*/
	DOF *dot_muo, *dot_muw, *dso_kro, *dsw_kro, *dsw_krw;
		dot_muo = phgDofNew(g, DOF_P0, 1, "dot_muo", DofNoAction);
		dot_muw = phgDofNew(g, DOF_P0, 1, "dot_muw", DofNoAction);
		dso_kro = phgDofNew(g, DOF_P0, 1, "dso_kro", DofNoAction);
		dsw_kro = phgDofNew(g, DOF_P0, 1, "dsw_kro", DofNoAction);
		dsw_krw = phgDofNew(g, DOF_P0, 1, "dsw_krw", DofNoAction);
		update_dot_muo(p_h, dot_muo);
		update_dot_muw(p_h, dot_muw);
		update_dot_kro(s_o, s_w, dso_kro, dsw_kro);
		update_dot_krw(s_w, dsw_krw);

	SOLVER *solver_uw;
	phgOptionsPush();
	phgOptionsSetKeyword("-hypre_solver", "boomeramg");
	solver_uw = phgSolverCreate(SOLVER_DEFAULT, u_w, NULL);
	phgOptionsPop();

	build_wateqn_uw(solver_uw, u_w, p_h, krw, b_w, mu_w, dp, d_sw, dot_bw, dot_muw, dsw_krw);
	phgSolverSolve(solver_uw, TRUE, u_w, NULL);
	phgSolverDestroy(&solver_uw);
	
//	Get_Uw(dp, d_so, d_sw, u_w, u_o, b_w, dot_bw, mu_w, dot_muw, krw, dsw_krw, b_o, dot_bo, mu_o, dot_muo, kro, dso_kro, dsw_kro);
	div_uw = phgDofDivergence(u_w, NULL, NULL, NULL);

	solver = phgSolverCreate(SOLVER_PCG, d_sw, NULL);
	build_wateqn_sw(solver, dp, d_sw, s_w, s_w_l, div_uw, q_w, p_h, p_h_newton, b_w, phi, dot_phi, dot_bw, phi_l, b_w_l);
	phgSolverSolve(solver, TRUE, d_sw, NULL);
	phgDofFree(&div_uw);
	phgDofFree(&dot_bw); 
	phgDofFree(&dot_phi);
	phgDofFree(&phi);
	phgDofFree(&b_w);
	phgDofFree(&mu_w);
	phgDofFree(&krw);
	phgDofFree(&b_o);
	phgDofFree(&dot_bo);
	phgDofFree(&dot_bw);
	phgDofFree(&mu_o);
	phgDofFree(&kro);
		/* Add DOFs For SS Methods*/
		phgDofFree(&dot_muo);
		phgDofFree(&dot_muw);
		phgDofFree(&dso_kro);
		phgDofFree(&dsw_kro);
		phgDofFree(&dsw_krw);
		/*-----------------------------------------*/
	phgSolverDestroy(&solver);
}
static void
build_oileqn_so(SOLVER *solver, DOF *dp, DOF *d_so, DOF *s_o, DOF *s_o_l, DOF *div_uo, DOF *q_o, DOF *p_h, DOF *p_h_newton, DOF *b_o, DOF *phi, DOF *dot_phi, DOF *dot_bo, DOF *phi_l, DOF *b_o_l)
{
	GRID *g = d_so->g;
	SIMPLEX *e;
	int i,j;
	int N = d_so->type->nbas;
	FLOAT A[N][N], rhs[N];
	INT I[N];
	FLOAT *p_phi, *p_bo, *p_so, *p_sol, *p_dotbo, *p_dotphi, *p_bol;
	ForAllElements(g, e){
		p_phi = DofElementData(phi, e->index);
		p_bo  = DofElementData(b_o, e->index);
		p_bol = DofElementData(b_o_l, e->index);
		p_so  = DofElementData(s_o, e->index);
		p_sol  = DofElementData(s_o_l, e->index);
		p_dotbo = DofElementData(dot_bo, e->index);
		p_dotphi = DofElementData(dot_phi, e->index);
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for (j = 0; j < N; j++){
				A[i][j] =  p_bo[0] * p_phi[0] * phgQuadBasDotBas(e, d_so, i, d_so, j, QUAD_DEFAULT); 
			}
		}
		FLOAT quad_qo = 0, quad_p = 0, quad_phi = 0, quad_phil = 0, quad_divuo = 0;
		FLOAT oil_cp = 0;
		oil_cp = p_so[0] * (p_bo[0] * p_dotphi[0] + p_phi[0] * p_dotbo[0]);
		for (i = 0; i < N; i++){
			phgQuadDofTimesBas(e, q_o, d_so, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, dp, d_so, i, QUAD_DEFAULT, &quad_p);
			phgQuadDofTimesBas(e, phi, d_so, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, phi_l, d_so, i, QUAD_DEFAULT, &quad_phil);
			phgQuadDofTimesBas(e, div_uo, s_o, i, QUAD_DEFAULT, &quad_divuo);
			rhs[i] = stime * quad_qo + p_sol[0] * p_bol[0] * quad_phil - oil_cp * quad_p  - p_so[0] * p_bo[0] * quad_phi - stime * quad_divuo;
		}
		for (i = 0; i < N; i++){
			phgSolverAddMatrixEntries(solver, 1, I+i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
}
static void
Solver_oileqn_so(DOF *dp, DOF *d_so, DOF *d_sw, DOF *s_o, DOF *u_o, DOF *q_o, DOF *s_o_l, DOF *p_h, DOF *p_h0, DOF *p_h_newton, DOF *phi_l, DOF *b_o_l)
{
	GRID *g = d_so->g;
	SOLVER *solver;
	DOF *div_uo, *dot_phi, *dot_bo, *b_o, *phi;
	div_uo = phgDofDivergence(u_o, NULL, NULL, NULL);
	solver = phgSolverCreate(SOLVER_PCG, d_so, NULL);
	dot_bo = phgDofNew(g, DOF_P0, 1, "dot_bo", DofNoAction);
	dot_phi = phgDofNew(g, DOF_P0, 1, "dot_phi", DofNoAction);
	update_dot_bo(p_h, dot_bo);
	update_dot_phi(p_h, dot_phi);
	phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
	update_phi(p_h, p_h0, phi);
	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
	update_bo(p_h, b_o);
	build_oileqn_so(solver, dp, d_so, s_o, s_o_l, div_uo, q_o, p_h, p_h_newton, b_o, phi, dot_phi, dot_bo, phi_l, b_o_l);
	phgSolverSolve(solver, TRUE, d_so, NULL);
	phgDofFree(&div_uo);
	phgDofFree(&dot_bo); 
	phgDofFree(&dot_phi);
	phgDofFree(&phi);
	phgDofFree(&b_o);
	phgSolverDestroy(&solver);
}

static void
Get_Sg(DOF *s_o, DOF *s_w, DOF *s_g)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_sw, *p_sg;
	ForAllElements(g, e){
		p_so = DofElementData(s_o, e->index);
		p_sw = DofElementData(s_w, e->index);
		p_sg = DofElementData(s_g, e->index);
		p_sg[0]= 1. - p_sw[0] - p_so[0];
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
Disp_Saturation(DOF *s_o, DOF *s_w, DOF *s_g)
{
	GRID *g = s_o->g;
	SIMPLEX *e;
	FLOAT so_min = 0, so_max = 0, sw_min = 0, sw_max = 0, sg_min = 0, sg_max = 0;
	Get_Sg(s_o, s_w, s_g);
	Find_Min_Max(s_o, &so_min, &so_max);
	Find_Min_Max(s_w, &sw_min, &sw_max);
	Find_Min_Max(s_g, &sg_min, &sg_max);
	phgPrintf("S_O      : Min =  %le,   Max =  %le\n", so_min, so_max);
	phgPrintf("S_W      : Min =  %le,   Max =  %le\n", sw_min, sw_max);
	phgPrintf("S_G      : Min =  %le,   Max =  %le\n", sg_min, sg_max);
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
    phgOptionsRegisterFilename("mesh_file", "Mesh_file", (char **)&fn);
    phgOptionsRegisterFloat("step", "Time step", &stime);
    phgInit(&argc, &argv);
    g = phgNewGrid(-1);

    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
	
    phgRefineAllElements(g, pre_refines); 
		if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
				(double)g->lif);
	phgOptionsShowUsed();
    /*The pressure variable*/
    DOF *p_h;
    p_h = phgDofNew(g, DOF_P0, 1, "p_h", DofNoAction);
    phgDofSetDataByValue(p_h, PRESSURE0);
    DOF *p0_h;
    p0_h = phgDofNew(g, DOF_P0, 1, "p0_h", DofNoAction);
    phgDofSetDataByValue(p0_h, PRESSURE0);
	DOF *dp;
	dp = phgDofNew(g, DOF_P0, 1, "dp", DofNoAction);
	phgDofSetDataByValue(dp, 0.);

    DOF *s_o;
    s_o = phgDofNew(g, DOF_P0, 1, "s_o", DofNoAction);
    phgDofSetDataByValue(s_o, SO0);
	DOF *d_so;
	d_so = phgDofNew(g, DOF_P0, 1, "d_so", DofNoAction);
	phgDofSetDataByValue(d_so, 0.);
    DOF *s_w;
    s_w = phgDofNew(g, DOF_P0, 1, "s_w", DofNoAction);
    phgDofSetDataByValue(s_w, SW0);
	DOF *d_sw;
	d_sw = phgDofNew(g, DOF_P0, 1, "d_sw", DofNoAction);
	phgDofSetDataByValue(d_sw, 0.);
    DOF *s_g;
    s_g = phgDofNew(g, DOF_P0, 1, "s_g", DofNoAction);
    phgDofSetDataByValue(s_g, SG0);
   /*The velocity variable*/
    DOF *u_o;
    u_o = phgDofNew(g, DOF_RT1, 1, "u_o", DofNoAction);
    DOF *u_w;
    u_w = phgDofNew(g, DOF_RT1, 1, "u_w", DofNoAction);
    DOF *u_g;
    u_g = phgDofNew(g, DOF_RT1, 1, "u_g", DofNoAction);
    /* RHS function */
    DOF *q_g;
    q_g = phgDofNew(g, DOF_P0, 1, "q_g", DofNoAction);
    DOF *q_w;			     
    q_w = phgDofNew(g, DOF_P0, 1, "q_w", DofNoAction);
    DOF *q_o;			     
    q_o = phgDofNew(g, DOF_P0, 1, "q_o", DofNoAction);
    Well_init(q_o);
    Unit_Conversion();
	phgPrintf("=======================================================\n");
	phgPrintf("the elements is :%d\n",DofGetDataCountGlobal(p_h));
	phgPrintf("\nThis Program Solve Black Oil Moedel For bubble Point\n \
    	Problem Use SS Method \n");
	phgPrintf("=======================================================\n");
	int flag = 0;
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
	 	DOF *phi_l, *b_o_l, *Rs_l, *b_g_l, *s_o_l, *s_w_l, *b_w_l;
		phi_l = phgDofNew(g, DOF_P0, 1, "phi_l", DofNoAction);
		update_phi(p_h, p0_h, phi_l);
	 	b_o_l = phgDofNew(g, DOF_P0, 1, "b_o_l", DofNoAction);
		update_bo(p_h, b_o_l);
	 	b_g_l = phgDofNew(g, DOF_P0, 1, "b_o_l", DofNoAction);
		update_bg(p_h, b_g_l);
		s_o_l = phgDofCopy(s_o, NULL, NULL, NULL);
		s_w_l = phgDofCopy(s_w, NULL, NULL, NULL);
	 	b_w_l = phgDofNew(g, DOF_P0, 1, "b_w_l", DofNoAction);
		update_bw(p_h, b_w_l);
	 	Rs_l = phgDofNew(g, DOF_P0, 1, "Rs_l", DofNoAction);
		update_Rs(p_h, Rs_l);
		int count = 0;
		while (TRUE){
			count++;
			DOF *p_h_newton, *s_o_newton, *s_w_newton;
			p_h_newton = phgDofCopy(p_h, NULL, NULL, NULL);
			s_o_newton = phgDofCopy(s_o, NULL, NULL, NULL);
			s_w_newton = phgDofCopy(s_w, NULL, NULL, NULL);
			Solve_Pressure(dp, d_so, d_sw, u_o, p_h, p0_h, p_h_newton, q_o, s_o, s_o_l, b_o_l, q_w, s_w, s_w_l, b_w_l, q_g, phi_l, b_g_l, Rs_l, s_g);
			phgDofAXPBY(1.0, dp, 1.0, &p_h);
			Solver_wateqn_sw(dp, d_so, d_sw, s_w, s_o, s_g, u_w, u_o, q_w, s_w_l, p_h, p0_h, p_h_newton, phi_l, b_w_l);
			phgDofAXPBY(1.0, d_sw, 1.0, &s_w);
			Solver_oileqn_so(dp, d_so, d_sw, s_o, u_o, q_o, s_o_l, p_h, p0_h, p_h_newton, phi_l, b_o_l);
			phgDofAXPBY(1.0, d_so, 1.0, &s_o);

			Get_Sg(s_o, s_w, s_g);

			FLOAT err_p = 0.0, norm_p = 0.0, TOL_p = 0.0, norm = 0;
			err_p = L2(p_h_newton, p_h, 4);
			norm_p = L2_norm(p_h_newton);
			norm = L2_norm(p_h);
			norm_p = (norm > norm_p)?norm:norm_p;
			TOL_p = err_p / norm_p;
			phgDofFree(&p_h_newton);
			FLOAT err_so = 0.0, norm_so = 0.0, TOL_so = 0.0;
			err_so = L2(s_o_newton, s_o, 4);
			norm_so = L2_norm(s_o_newton);
			norm = L2_norm(s_o);
			norm_so = (norm > norm_so)?norm:norm_so;
			TOL_so = err_so / norm_so;
			phgDofFree(&s_o_newton);
			FLOAT err_sw = 0.0, norm_sw = 0.0, TOL_sw = 0.0;
			err_sw = L2(s_w_newton, s_w, 4);
			norm_sw = L2_norm(s_w_newton);
			norm = L2_norm(s_w);
			norm_sw = (norm > norm_sw)?norm:norm_sw;
			TOL_sw = err_sw / norm_sw;
			phgDofFree(&s_w_newton);
			if ((TOL_p < 1e-6) && (TOL_so < 1e-6) && (TOL_sw < 1e-6)){
				phgPrintf("err_p =       %le\n", TOL_p);
				phgPrintf("err_so =      %le\n", TOL_so);
				phgPrintf("err_sw =      %le\n", TOL_sw);
				phgPrintf("Non_ints = %d\n", count);
				break;
			}
		}
		Disp_Saturation(s_o, s_w, s_g);
		FLOAT pwf = 0;
		Well_Pressure(p_h, &pwf);
		phgPrintf("t = %lf, Pwf = %lf\n", ctime, pwf);
		phgDofFree(&b_o_l);
		phgDofFree(&b_w_l);
		phgDofFree(&b_g_l);
		phgDofFree(&phi_l);
		phgDofFree(&s_o_l);
		phgDofFree(&s_w_l);
		phgDofFree(&Rs_l);
    }
	phgDofFree(&dp);
	phgDofFree(&d_so);
	phgDofFree(&d_sw);
    phgDofFree(&p_h);
    phgDofFree(&p0_h);
    phgDofFree(&u_o);
    phgDofFree(&u_w);
    phgDofFree(&u_g);
    phgDofFree(&s_o);
    phgDofFree(&s_w);
    phgDofFree(&s_g);
    phgDofFree(&q_o);
    phgDofFree(&q_w);
    phgDofFree(&q_g);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
