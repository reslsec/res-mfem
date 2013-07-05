 /* This sample code solves the oil_water_phase equation:         *
 * ****************************************************************
 * problem: p_{t} -\Delta{p} = func_f(x, t) \in \Omega X [0, T]   *
 *          p = func_g (x, t) \in \partial\Omega X [0, T]         * 
 *          p = func_p0       \in \Omega t = 0.0                  *
 *WARNING! The unit is h, m, mpa!!                                *
 revision 1.2
date: 2011/05/05 02:20:46;  author: liuming;  state: Exp;
 ******************************************************************/
#include "phg.h"
#include <string.h>
#include <math.h>
#include "well.c"
#include "utils.c"
#include "parameter.c"
#include "uzawa.c"
/* build linear system */
static void
build_mat_vec(DOF *u_w, DOF *dp, DOF *p_nk, DOF *ds, DOF *s_w, DOF *s_w_l, DOF *b_o, DOF *b_o_l, DOF *kro, DOF *dot_kro, DOF *q_o, DOF *b_w, DOF *b_w_l, DOF *krw, DOF *dot_krw, DOF *q_w, DOF *phi, DOF *phi_l, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g)
{
	GRID *g = u_w->g;
	SIMPLEX *e;
	FLOAT *p_bo, *p_bw, *p_bol, *p_bwl, *p_kro, *p_dotkro, *p_krw, *p_dotkrw, *p_swl, *p_dp, *p_ds, *p_sw, *p_phi;
	INT N = u_w->type->nbas * u_w->dim;
	INT M = dp->type->nbas * dp->dim;
	INT I[N], J[M];
	FLOAT mat_A[N][N], mat_TB[N][M], mat_B[M][N], mat_C[M][M], rhs_f[N], rhs_g[M];
	phgVecDisassemble(vec_f);
	phgVecDisassemble(vec_g);
	int i, j, k;
	ForAllElements(g, e){
		p_phi = DofElementData(phi, e->index);
		p_bo = DofElementData(b_o, e->index);
		p_bol = DofElementData(b_o_l, e->index);
		p_bw = DofElementData(b_w, e->index);
		p_bwl = DofElementData(b_w_l, e->index);
		p_sw = DofElementData(s_w, e->index);
		p_swl = DofElementData(s_w_l, e->index);
		p_kro  = DofElementData(kro, e->index);
		p_dotkro  = DofElementData(dot_kro, e->index);
		p_krw  = DofElementData(krw, e->index);
		p_dotkrw  = DofElementData(dot_krw, e->index);
		p_dp = DofElementData(dp, e->index);
		p_ds = DofElementData(ds, e->index);
		FLOAT T_o = 0, T_w = 0;
		T_w = K * p_krw[0] * p_bw[0] / MU_W + K * p_krw[0] * C_W / (MU_W * B0_W)  * p_dp[0] + K * p_bw[0] * p_dotkrw[0] / MU_W * p_ds[0]; 
		T_o = K * p_kro[0] * p_bo[0] / MU_O + K * p_kro[0] * C_O / (MU_O * B0_O)  * p_dp[0] + K * p_bo[0] * p_dotkro[0] / MU_O * p_ds[0]; 
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				mat_A[i][j] = stime * phgQuadBasDotBas(e, u_w, i, u_w, j, QUAD_DEFAULT) / T_w;
			}
		}
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				mat_TB[i][j] = -stime * phgQuadDivBasDotBas(e, u_w, i, dp, j, QUAD_DEFAULT);					
			}
		}
		FLOAT wat_sw = 0, oil_sw = 0, wat_cp = 0, oil_cp = 0, beta = 0;
		for (i = 0; i < M; i++){
			for (j = 0; j < M; j++){
				wat_sw = p_phi[0] * p_bw[0] * phgQuadBasDotBas(e, ds, i, ds, j, QUAD_DEFAULT);
				oil_sw = p_phi[0] * p_bo[0] * phgQuadBasDotBas(e, ds, i, ds, j, QUAD_DEFAULT);
			}
		}
		wat_cp = p_sw[0] * (p_bw[0] * PHI0 * C_R + p_phi[0] * C_W / B0_W);
		oil_cp = (1. - p_sw[0]) * (p_bo[0] * PHI0 * C_R + p_phi[0] * C_O / B0_O);
		beta = 1. / wat_sw + T_o / (T_w * oil_sw);
		FLOAT quad = 0.;
		for (i = 0; i < M; i++){
			for (j = 0; j < M; j++){
			 	quad =  phgQuadBasDotBas(e, dp, i, dp, j, QUAD_DEFAULT);
				mat_C[i][j] = (wat_cp / wat_sw + oil_cp / oil_sw) * quad / beta;
			}
		}
		/*create oil rhs*/
		FLOAT quad_phi = 0., quad_phil = 0., quad_qo = 0, quad_qw = 0;
		FLOAT rhs_oil = 0., rhs_wat = 0.;
		for (i = 0; i < M; i++){
			phgQuadDofTimesBas(e, q_o, dp, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, q_w, dp, i, QUAD_DEFAULT, &quad_qw);
			phgQuadDofTimesBas(e, phi, dp, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, phi_l, dp, i, QUAD_DEFAULT, &quad_phil);
			rhs_oil = -stime * quad_qo + p_bo[0] * (1. - p_sw[0]) * quad_phi - p_bol[0] * (1. - p_swl[0]) * quad_phil;
			rhs_wat = -stime * quad_qw + p_bw[0] * p_sw[0] * quad_phi - p_bwl[0] * p_swl[0] * quad_phil;
			rhs_g[i] = (rhs_wat / wat_sw + rhs_oil / oil_sw) / beta;
		}
		for (i = 0; i < N; i++){
			rhs_f[i] = stime * phgQuadDofTimesDivBas(e, p_nk, u_w, i, QUAD_DEFAULT);
		}
		/* Handle Bdry Conditions*/
		for (i = 0; i < N; i++){
			if (phgDofGetElementBoundaryType(u_w, e, i) & (NEUMANN | DIRICHLET)){
				bzero(mat_A[i], N * sizeof(mat_A[i][0]));
				bzero(mat_TB[i], M * sizeof(mat_TB[i][0]));
				for (j = 0; j < N; j++){
					mat_A[j][i] = 0.;
				}
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
Solve_Pressure(DOF *u_w, DOF *dp, DOF *p_nk, DOF *p_h, DOF *p0_h, DOF *q_o, DOF *ds, DOF *s_w, DOF *s_w_l, DOF *q_w, DOF *phi_l, DOF *b_o_l, DOF *b_w_l)
{
		GRID *g = u_w->g;	
	     MAT *A, *B, *TB, *C;
	     VEC *vec_f, *vec_g;
	     MAP *map_u, *map_p;
		SOLVER *pc;
//		SOLVER *solver;
     	/* The parameter DOF */
 	  	DOF *phi, *b_o, *kro, *dot_kro, *b_w, *krw, *dot_krw;
		/*Create MAP for Mat and Vec*/
	     map_p = phgMapCreate(dp, NULL);
		map_u = phgMapCreate(u_w, NULL);
		A = phgMapCreateMat(map_u, map_u);
		B = phgMapCreateMat(map_p, map_u);
		TB = phgMapCreateMat(map_u, map_p);
		C = phgMapCreateMat(map_p, map_p);
		vec_f = phgMapCreateVec(map_u, 1);
		vec_g = phgMapCreateVec(map_p, 1);
			
		phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
		update_phi(p_h, p0_h, phi);
 		b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
		update_bo(p_h, p0_h, b_o);	
 		b_w = phgDofNew(g, DOF_P0, 1, "b_w", DofNoAction);
		update_bw(p_h, p0_h, b_w);	
		kro = phgDofNew(g, DOF_P0, 1, "kro", DofNoAction);
		create_kro(s_w, kro);
		krw = phgDofNew(g, DOF_P0, 1, "krw", DofNoAction);
		create_krw(s_w, krw);
		dot_kro = phgDofNew(g, DOF_P0, 1, "dot_kro", DofNoAction);
		create_dot_kro(s_w, dot_kro);
		dot_krw = phgDofNew(g, DOF_P0, 1, "dot_krw", DofNoAction);
		create_dot_krw(s_w, dot_krw);
		water_fluidity(q_o, q_w, kro, krw, b_o, b_w);
		phgPrintf("build linear system:           ");
		use_time(g, FALSE, 0.);
		build_mat_vec(u_w, dp, p_nk, ds, s_w, s_w_l, b_o, b_o_l, kro, dot_kro, q_o, b_w, b_w_l, krw, dot_krw, q_w, phi, phi_l, map_u, map_p, A, B, TB, C, vec_f, vec_g);
		use_time(g, TRUE, 0.);
		phgDofFree(&phi);
		phgDofFree(&b_o);
		phgDofFree(&kro);
		phgDofFree(&b_w);
		phgDofFree(&krw);
		phgDofFree(&dot_kro);
		phgDofFree(&dot_krw);
    	/* new implementation */
		int nits_amg = 0;
        int nits_uzawa = 0;
		double time_amg = 0;
		use_time(g, FALSE, 0.);
		DOF *B_data = DivRT(dp, u_w);
        MAT *H = BTAB(A, B_data, dp, u_w);
      	phgDofFree(&B_data);
		phgPrintf("Assemble H:             ");
		use_time(g, TRUE, 0.);
		phgPrintf("solve p_h:              \n");
	//	nits_uzawa = phgUzawa(H, u_w, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg, &time_amg);
		nits_uzawa = uzawapcg(H, u_w, dp, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
		phgPrintf("Max iter of AMG---------%d\n", nits_amg);
		phgPrintf("Nits: uzawa:----------------%d\n", nits_uzawa);
		phgMatDestroy(&H);
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
water_fluidity(DOF *q_o, DOF *q_w, DOF *kro, DOF *krw, DOF *b_o, DOF *b_w)
{
	GRID *g = q_o->g;
	SIMPLEX *e;
	FLOAT lambda_o = 0, lambda_w = 0;
	FLOAT *p_qo, *p_qw, *p_kro, *p_krw, *p_bo, *p_bw;
	ForAllElements(g, e){
			p_qw = DofElementData(q_w, e->index);
			p_qo = DofElementData(q_o, e->index);
			p_kro = DofElementData(kro, e->index);
			p_krw = DofElementData(krw, e->index);
			p_bo = DofElementData(b_o, e->index);
			p_bw = DofElementData(b_w, e->index);
			lambda_o = K * p_kro[0] * p_bo[0] / MU_O;
			lambda_w = K * p_krw[0] * p_bw[0] / MU_W;
			p_qw[0] = lambda_w / lambda_o * p_qo[0];
	}
}
static void
build_water_sw(SOLVER *solver, DOF *ds, DOF *s_w, DOF *s_w_l, DOF *div_uw, DOF *q_w, DOF *dp, DOF *b_w, DOF *b_w_l, DOF *phi, DOF *phi_l)
{
	GRID *g = s_w->g;
	SIMPLEX *e;
	int i,j;
	int N = s_w->type->nbas;
	FLOAT A[N][N], rhs[N];
	INT I[N];
	FLOAT *p_phi, *p_bw, *p_bwl, *p_sw, *p_swl, *p_phil;
	ForAllElements(g, e){
		p_phi = DofElementData(phi, e->index);
		p_phil = DofElementData(phi_l, e->index);
		p_bw  = DofElementData(b_w, e->index);
		p_bwl = DofElementData(b_w_l, e->index);
		p_sw  = DofElementData(s_w, e->index);
		p_swl  = DofElementData(s_w_l, e->index);
		for (i = 0; i < N; i++){
			I[i] = phgSolverMapE2L(solver, 0, e, i);
			for (j = 0; j < N; j++){
				A[i][j] = p_phi[0] * p_bw[0] * phgQuadBasDotBas(e, ds, i, ds, j, QUAD_DEFAULT); 
			}
		}
		FLOAT quad_qw = 0, quad_ph = 0, quad_phil = 0, quad_phi = 0.,quad_divuw = 0;
		for (i = 0; i < N; i++){
			phgQuadDofTimesBas(e, q_w, ds, i, QUAD_DEFAULT, &quad_qw);
			phgQuadDofTimesBas(e, dp, ds, i, QUAD_DEFAULT, &quad_ph);
			phgQuadDofTimesBas(e, phi, ds, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, phi_l, ds, i, QUAD_DEFAULT, &quad_phil);
			phgQuadDofTimesBas(e, div_uw, ds, i, QUAD_DEFAULT, &quad_divuw);
			rhs[i] = stime * quad_qw  - stime * quad_divuw - (p_sw[0] * p_bw[0] * quad_phi - p_swl[0] * p_bwl[0] * quad_phil) - p_sw[0] * (p_bw[0] * PHI0 * C_R + p_phi[0] * C_W / B0_W) * quad_ph;
		}
		for (i = 0; i < N; i++){
			phgSolverAddMatrixEntries(solver, 1, I+i, N, I, A[i]);
		}
		phgSolverAddRHSEntries(solver, N, I, rhs);
	}
}
static void
Solver_Wateqn_Sw(DOF *ds, DOF *s_w, DOF *u_w, DOF *q_w, DOF *s_w_l, DOF *dp, DOF *p_h, DOF *p0_h, DOF *phi_l, DOF *b_w_l)
{
	GRID *g = s_w->g;
	SOLVER *solver;
	DOF *div_uw, *phi, *b_w;
	phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
	update_phi(p_h, p0_h, phi);
 	b_w = phgDofNew(g, DOF_P0, 1, "b_w", DofNoAction);
	update_bw(p_h, p0_h, b_w);	
	div_uw = phgDofDivergence(u_w, NULL, NULL, NULL);
	solver = phgSolverCreate(SOLVER_PCG, ds, NULL);
	build_water_sw(solver, ds, s_w, s_w_l, div_uw, q_w, dp, b_w, b_w_l, phi, phi_l);
	phgSolverSolve(solver, TRUE, ds, NULL);
	phgDofFree(&div_uw);
	phgDofFree(&b_w);
	phgDofFree(&phi);
	phgSolverDestroy(&solver);
}
int
main(int argc, char *argv[])
{
    INT mem_max = 10240;
    char *fn = "ow.mesh";
 //   char *fn = "one2.mesh";
    GRID *g;
    char vtkfile[1000];
    static int pre_refines = 0;
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgInit(&argc, &argv);
    g = phgNewGrid(-1);

    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
	
    phgRefineAllElements(g, pre_refines); 
		if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
				(double)g->lif);

    /*The pressure variable*/
    DOF *p_h; /*This is delta p*/
    p_h = phgDofNew(g, DOF_P0, 1, "p_h", DofNoAction);
    phgDofSetDataByValue(p_h, PRESSURE0);
    DOF *p0_h;
    p0_h = phgDofNew(g, DOF_P0, 1, "p0_h", DofNoAction);
    phgDofSetDataByValue(p0_h, PRESSURE0);
    DOF *dp;
    dp = phgDofNew(g, DOF_P0, 1, "dp", DofNoAction);
    phgDofSetDataByValue(dp, 0.);

    DOF *s_w;
    s_w = phgDofNew(g, DOF_P0, 1, "s_w", DofNoAction);
    phgDofSetDataByValue(s_w, SW0);
    DOF *ds;
    ds = phgDofNew(g, DOF_P0, 1, "ds", DofNoAction);
    phgDofSetDataByValue(ds, 0.);
   /*The velocity variable*/
    DOF *u_o;
    u_o = phgDofNew(g, DOF_RT1, 1, "u_o", DofNoAction);
    DOF *u_w;
    u_w = phgDofNew(g, DOF_RT1, 1, "u_w", DofNoAction);
    /* RHS function */
    DOF *q_w;
    q_w = phgDofNew(g, DOF_P0, 1, "q_w", DofNoAction);
    DOF *q_o;			     
    q_o = phgDofNew(g, DOF_P0, 1, "q_o", DofNoAction);
    Well_init(q_o);
 
	phgPrintf("The elements is  :    %d\n",DofGetDataCountGlobal(p_h));
	phgPrintf("The total DOFs is:    %d\n",DofGetDataCountGlobal(p_h) + DofGetDataCountGlobal(u_o));

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
	 	DOF *b_w_l, *b_o_l, *s_w_l, *phi_l;
	 	b_w_l = phgDofNew(g, DOF_P0, 1, "b_w_l", DofNoAction);
		update_bw(p_h, p0_h, b_w_l);
	 	b_o_l = phgDofNew(g, DOF_P0, 1, "b_o_l", DofNoAction);
		update_bo(p_h, p0_h, b_o_l);
		s_w_l = phgDofCopy(s_w, NULL, NULL, NULL);
		phi_l = phgDofNew(g, DOF_P0, 1, "phi_l", DofNoAction);
		update_phi(p_h, p0_h, phi_l);
		int count = 0;
		while (TRUE){
			count++;
			DOF *p_nk, *s_nk;
			p_nk = phgDofCopy(p_h, NULL, NULL, NULL);
			s_nk = phgDofCopy(s_w, NULL, NULL, NULL);
			
			Solve_Pressure(u_w, dp, p_nk, p_h, p0_h, q_o, ds, s_w, s_w_l, q_w, phi_l, b_o_l, b_w_l);
			Solver_Wateqn_Sw(ds, s_w, u_w, q_w, s_w_l, dp, p_h, p0_h, phi_l, b_w_l);
			
			phgDofAXPBY(1.0, dp, 1.0, &p_h);
			phgDofAXPBY(1.0, ds, 1.0, &s_w);

			FLOAT err_p = 0.0, norm_p = 0.0, TOL_p = 0.0;
			FLOAT err_s = 0.0, norm_s = 0.0, TOL_s = 0.0;
			err_p = L2(p_nk, p_h, 4);
			norm_p = (L2_norm(p_nk) > L2_norm(p_h))?L2_norm(p_nk):L2_norm(p_h);
			TOL_p = err_p / norm_p;
			err_s = L2(s_nk, s_w, 4);
			norm_s = (L2_norm(s_nk) > L2_norm(s_w))?L2_norm(s_nk):L2_norm(s_w);
			TOL_s = err_s / norm_s;
			phgDofFree(&p_nk);
			phgDofFree(&s_nk);
			if ((TOL_p < 1e-6)  & (TOL_s < 1e-6)){
				phgPrintf("TOL_p:                    %le\n", TOL_p);
				phgPrintf("TOL_s:                    %le\n", TOL_s);
				phgPrintf("Non_ints:                 %d\n", count);
				break;
			}
		}
		FLOAT pwf = 0;
		Well_Pressure(p_h, &pwf);
		phgPrintf("t = %lf, Pwf = %lf\n", ctime, pwf);
		DOF *b_o, *kro;
	 	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
		update_bo(p_h, p0_h, b_o);	
	 	kro = phgDofNew(g, DOF_P0, 1, "kro", DofNoAction);
		create_kro(s_w, kro);
		FLOAT pwf0 = 0;
		Peaceman_Model(p_h, b_o, kro, &pwf0);
		phgPrintf("USE_Peaceman   t = %lf, Pwf = %lf\n", ctime, pwf0);
#if 0
		MAP *map_p = phgMapCreate(p_h, NULL);
		VEC *y = phgMapCreateVec(map_p, 1);
 	  	phgMapDofToLocalData(map_p, 1, &p_h, y->data);
		phgVecDumpMATLAB(y, "p", "p.m");
		phgMapDestroy(&map_p);
		phgVecDestroy(&y);
		MAP *map_s = phgMapCreate(s_w, NULL);
		VEC *x = phgMapCreateVec(map_s, 1);
	 	phgMapDofToLocalData(map_s, 1, &s_w, x->data);
	  	phgVecDumpMATLAB(x, "s", "s.m");
		phgMapDestroy(&map_s);
		phgVecDestroy(&x);
		flag++;
		sprintf(vtkfile, "new_%03d.vtk", flag);
		/*Create VTK files*/
		phgExportVTK(g, vtkfile, p_h, s_w, q_o, q_w, u_o, u_w, NULL);
#endif
		phgDofFree(&b_o_l);
		phgDofFree(&b_w_l);
		phgDofFree(&phi_l);
		phgDofFree(&s_w_l);
    }
    phgDofFree(&p_h);
    phgDofFree(&dp);
    phgDofFree(&p0_h);
    phgDofFree(&u_o);
    phgDofFree(&u_w);
    phgDofFree(&s_w);
    phgDofFree(&ds);
    phgDofFree(&q_o);
    phgDofFree(&q_w);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
