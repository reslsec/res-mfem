 /* This sample code solves the oil_water_phase equation:         *
 * ****************************************************************
 * problem: p_{t} -\Delta{p} = func_f(x, t) \in \Omega X [0, T]   *
 *          p = func_g (x, t) \in \partial\Omega X [0, T]         * 
 *          p = func_p0       \in \Omega t = 0.0                  *
 *WARNING! The unit is h, m, mpa!!                                *
 revision 1.1
date: 2011/04/29 01:17:21;  author: liuming;  state: Exp;
 ******************************************************************/
#include "phg.h"
#include <string.h>
#include <math.h>
#include "parameter.c"
#include "well.c"
#include "quadfunc.c"
#define USE_PC 1
#define USE_UZAWA 0
#define USE_BLOCK 1
#define DEBUG 0
#define USE_OIL_EQN 1
#define USE_GAS_EQN 0
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
				mat_A[i][j] = stime * p_muo[0] * phgQuadBasDotBas(e, u_o, i, u_o, j, QUAD_DEFAULT) / (p_kro[0] * p_bo[0] * K);
			}
		}
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				mat_TB[i][j] = -stime * phgQuadDivBasDotBas(e, u_o, i, p_h, j, QUAD_DEFAULT);
			}
		}
		FLOAT beta = 0;
		beta = 1. / oil_so + (p_Rs[0] + p_krg[0] * p_bg[0] * p_muo[0] / (p_kro[0] * p_bo[0] * p_mug[0])) / gas_so;
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
/*--------------------------------------------------*/
/*  parellel matmultmat implemented by cheng jie    */
/*    2010.5.25                                     */ 
/*--------------------------------------------------*/
static DOF *
DivRT(DOF *u, DOF *p)
{
    GRID *g = p->g;
    SIMPLEX *e;
    DOF *B;
    FLOAT v;
    int N = p->type->nbas;
    int M = u->type->nbas;
    int i, j;

    B = phgDofNew(g, DOF_P0, N * M, "div bas", DofNoAction);

    ForAllElements(g, e) {
	for (i = 0; i < N; i++) {
	    for (j = 0; j < M; j++) {
		v = -stime * phgQuadDivBasDotBas(e, p, i, u, j, QUAD_DEFAULT);
		*(DofElementData(B, e->index) + i * M + j) = v;
	    }
	}
    }
    return B;
}
static MAT *
BTAB(MAT * A, DOF *B, DOF *u, DOF *p)
{
    GRID *g = B->g;
    SIMPLEX *e, *ne;
    MAP *map = phgMapCreate(u, NULL);
    MAT *S = phgMapCreateMat(map, map);	/* the Schur complement */
    VEC *diagA;
    DOF *diagA2Dof;
    NEIGHBOUR_DATA *nd, *ndu;
    int i, nface;
    INT lindex, gindex;
    FLOAT lvalue, gvalue, a;

    assert(u->type == DOF_P0);
    assert(u->dim == 1);

    if (A->diag != NULL)
	phgFree(A->diag);
    phgMatSetupDiagonal(A);

    diagA = phgMapCreateVec(A->rmap, 1);
    memcpy(diagA->data, A->diag, A->rmap->nlocal * sizeof(*(diagA->data)));
    phgFree(A->diag);
    A->diag = NULL;

    diagA2Dof = phgDofNew(g, p->type, 1, "diagonal data", DofNoAction);
    phgMapVecToDofArrays(diagA->map, diagA, FALSE, &diagA2Dof, NULL);
    phgVecDestroy(&diagA);

    /* to get neighbouring B */
    nd = phgDofInitNeighbourData(B, NULL);
    /* to get neighbouring element global vector index */
    ndu = phgDofInitNeighbourData(u, map);

    ForAllElements(g, e) {
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & (NEUMANN | DIRICHLET))
		continue;

	    a = *(DofFaceData(diagA2Dof, e->faces[i]));

	    /* e <--> e */
	    lindex = phgMapE2L(map, 0, e, 0);
	    lvalue = *(DofElementData(B, e->index) + i) *
		*(DofElementData(B, e->index) + i) / a;
	    phgMatAddEntry(S, lindex, lindex, lvalue);

	    if (!(e->bound_type[i] & INTERIOR))
		continue;

	    /* e <--> e->neighbour[i] */
	    if(e->bound_type[i] & REMOTE){
		/* remote neighbour */
		int ind = (size_t)e->neighbours[i];
		nface = (g->neighbours.list + ind)->op_vertex;
		gindex = *(ndu->index + ind * u->dim);
		gvalue = *(nd->data + ind * B->dim + nface);
		gvalue *= *(DofElementData(B,e->index)+i) / a;
	    }
	    else{ /* local neighbour */
		ne = (SIMPLEX *)e->neighbours[i];
		nface = phgOppositeVertex(g,e,i,ne);
		gindex = phgMapE2G(map, 0, ne, 0);
		gvalue = *(DofElementData(B, e->index) + i) *
		    *(DofElementData(B, ne->index) + nface) / a;
	    }
	    phgMatAddLGEntry(S, lindex, gindex, gvalue);
	}
    }
    phgMapDestroy(&map);
    phgDofFree(&diagA2Dof);
    phgDofReleaseNeighbourData(&nd);
    phgDofReleaseNeighbourData(&ndu);

    return S;
}
static void
pc_proc(SOLVER *pc_solver, VEC *b0, VEC **x0)
{
	FLOAT *old_rhs;
	VEC *b, *x;
	INT n = pc_solver->mat->rmap->nlocal;
	assert(n == b0->map->nlocal);
	phgSolverAssemble(pc_solver);
	b = phgMapCreateVec(pc_solver->mat->rmap, 1);
	x = phgMapCreateVec(pc_solver->mat->rmap, 1);
	/* setup RHS */
	old_rhs = pc_solver->rhs->data;
	pc_solver->rhs->data = b->data;
	pc_solver->rhs->assembled = TRUE;
	pc_solver->rhs_updated = TRUE;
	/* Updated pc_solver */
	memcpy(b->data, b0->data, sizeof(*b->data) * n);
	bzero(x->data, n * sizeof(*x->data));
	phgSolverVecSolve(pc_solver, FALSE, x);
	memcpy((*x0)->data, x->data, sizeof(*x->data) * n);

	pc_solver->rhs->data = old_rhs;

	phgVecDestroy(&b);
	phgVecDestroy(&x);
}
static FLOAT
phgUzawa(SOLVER *pc_solver, MAT *H0, DOF *u, DOF *p, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g, INT *nits)
{
	 VEC *x, *y, *TBy, *d1, *d2, *fr, *gr;
	 MAT *H;
	 SOLVER *solver1, *solver2, *solver3;
	 FLOAT value = 0., b1=0.0, b2=0.0,f=0.0, g=0.0, res1 =0., error = 1e-10, res = 0, res_f = 0, res_g = 0;
	 INT inter_amg = 0, Nits = 0, inter_pcg = 0, inter_max = 0;
	 int i;

	 x = phgMapCreateVec(map_u, 1);
	 TBy = phgMapCreateVec(map_u, 1);
	 d1 = phgMapCreateVec(map_u, 1);
	 fr = phgMapCreateVec(map_u,1);

	 y = phgMapCreateVec(map_p, 1);
	 d2 = phgMapCreateVec(map_p, 1);
	 gr = phgMapCreateVec(map_p, 1);
	 H = phgMapCreateMat(map_p, map_p);

	 phgMapDofToLocalData(map_u, 1, &u, x->data);
      phgMapDofToLocalData(map_p, 1, &p, y->data);

	 H = H0;
	 phgMatAXPBY(1.0, C, 1.0, &H);

	 while (TRUE){
		  // f = f - (AX + TBy)	
		  Nits++;
      	  phgMatVec(0, 1.0, TB, y, 0.0, &TBy);
		  phgMatVec(0, 1.0, A, x, 1.0, &TBy); //  TBy = Ax + TBy
		  phgVecCopy(vec_f, &fr);
		  phgVecAXPBY(-1.0, TBy, 1.0, &fr);
		  res_f = phgVecNorm2(fr, 0, NULL);
		  f = phgVecNorm2(vec_f, 0, NULL);
		  //solve Ad1 = f-(Ax+TBy)
		  solver1 = phgMat2Solver(SOLVER_PCG, A);
//		  printf("Right\n");
		  solver1->mat->handle_bdry_eqns = FALSE;
		  phgVecCopy(fr, &(solver1->rhs));
		  phgSolverVecSolve(solver1, FALSE, d1);
		  phgVecAXPBY(1.0, d1, 1.0, &x);
		  inter_pcg = solver1->nits;
		  phgSolverDestroy(&solver1);
		  // g = Bx-Cy-g
		  phgVecCopy(vec_g, &gr);
		  phgMatVec(0, 1.0, B, x, -1.0, &gr);
		  phgMatVec(0, -1.0, C, y, 1.0, &gr);
		  res_g = phgVecNorm2(gr, 0, NULL);

    		  phgOptionsPush();
		  phgOptionsSetKeyword("-hypre_solver", "boomeramg");
		  //phgOptionsSetKeyword("-hypre_pc", "none");	
		  solver2 = phgMat2Solver(SOLVER_DEFAULT, H);
		  phgOptionsPop();

		  phgVecCopy(gr, &(solver2->rhs));
		  if(USE_PC) {
		 	 phgSolverSetPC(solver2, pc_solver, pc_proc);
		  }
		  phgSolverVecSolve(solver2, FALSE, d2);
		  inter_amg = solver2->nits;
		  if (inter_amg > inter_max){
		  	 inter_max = inter_amg;
		  }

		  phgSolverDestroy(&solver2);
		  res1 = phgVecNorm2(d2, 0, NULL);
		  if (res1 != 0.0 ){
		 		b1 = phgVecDot(d2, 0, gr, 0, NULL);
				for (i = 0; i < TBy->map->nlocal; i++){
					TBy->data[i] = 0.0;
				}
				phgMatVec(0, 1.0, TB, d2, 0.0, &TBy);
				solver3 = phgMat2Solver(SOLVER_PCG, A);
		          phgVecCopy(TBy, &(solver3->rhs));
				for (i = 0; i < d1->map->nlocal; i++){
					d1->data[i] = 0.0;
				}
				phgSolverVecSolve(solver3, FALSE, d1);
				b2 = phgVecDot(d1, 0, TBy, 0, NULL);
				phgSolverDestroy(&solver3);
				value = b1 / b2;
		  }else{
		  		value = 1.0;
		  }
	 	 
		  phgVecAXPBY(7.0/10.0*value, d2, 1.0, &y);
		  //compute the residual ||b-Ax||
		  phgMatVec(0, 1.0, TB, y, 1.0, &TBy);
		  phgVecCopy(vec_f, &fr);
		  phgMatVec(0, 1.0, A, x, 1.0, &TBy);
		  phgVecAXPBY(-1.0, TBy, 1.0, &fr);
		  res_f = phgVecNorm2(fr, 0, NULL);
//		  phgPrintf("res_f = %lf", res_f);
		  f = phgVecNorm2(vec_f, 0, NULL);
		  if (f == 0){
		  	res_f = res_f;
		  }else{
		  	res_f = res_f / f;
		  }
		  g = phgVecNorm2(vec_g, 0, NULL);

	//	  phgPrintf("res_g = %lf\n",res_g);
		  if (g == 0){
		  	res_g = res_g;
		  }else{
		  	res_g = res_g / g;
		  }
		  res = Sqrt(res_f * res_f + res_g * res_g);
		  if (res < error){
		  		phgMapLocalDataToDof(map_u, 1, &u, x->data);
				phgMapLocalDataToDof(map_p, 1, &p, y->data);
				break;
		  }
	 }
	 phgVecDestroy(&d2);
	 phgVecDestroy(&TBy);
	 phgVecDestroy(&d1);
	 phgVecDestroy(&x);
	 phgVecDestroy(&y);
	 phgVecDestroy(&fr);
	 phgVecDestroy(&gr);

	 *nits = inter_max;
	 return Nits;
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
			lambda_o = p_kro[0] * p_bo[0] / p_muo[0];
			lambda_g = p_krg[0] * p_bg[0] / p_mug[0];
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
  	     int nits_amg = 0;
          int nits_uzawa = 0;
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
    		/* new implementation */
		DOF *B_data = DivRT(p_h, u_o);
          MAT *H = BTAB(A, B_data, p_h, u_o);
         	phgDofFree(&B_data);
		/* Set PC solver for schur complement*/
		SOLVER *pc;
		if(USE_PC){
			phgOptionsPush();
			phgOptionsSetKeyword("-hypre_solver", "boomeramg");
			phgOptionsSetKeyword("-hypre_pc", "none");	
			pc = phgMat2Solver(SOLVER_DEFAULT, H);
			pc->mat->handle_bdry_eqns = FALSE;
			pc->maxit = 3;
			pc->rtol = 1e-2;
			pc->warn_maxit = FALSE;
			phgOptionsPop();
		}
		elapsed_time(g, FALSE, 0.);
		nits_uzawa = phgUzawa(pc, H, u_o, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
		phgMatDestroy(&H);
		phgSolverDestroy(&pc);
		phgPrintf("Nits_Uzawa = %d\n",nits_uzawa);
		phgPrintf("Nits_PAMG = %d\n", nits_amg);
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
		elapsed_time(g, TRUE, 0.);
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
		FLOAT pwf = 0;
		Well_Pressure(p_h, &pwf);
		phgPrintf("t = %lf, Pwf = %lf\n", ctime, pwf);
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
		phgDofFree(&b_o_l);
		phgDofFree(&b_g_l);
		phgDofFree(&phi_l);
		phgDofFree(&s_o_l);
		phgDofFree(&Rs_l);
		/*Create VTK files*/
		DOF *kro = phgDofNew(g, DOF_P0, 1, "kro", DofNoAction);
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
