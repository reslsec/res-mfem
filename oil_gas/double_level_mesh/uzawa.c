#include <math.h>
#include "phg.h"
#include <string.h>
//#include "ss_oil_water.h"

/*--------------------------------------------------*/
/*  parellel matmultmat implemented by cheng jie    */
/*    2010.5.25                                     */ 
/*--------------------------------------------------*/
FLOAT stime = 1.;
DOF *
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
MAT *
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
void
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
FLOAT
phgUzawa(MAT *H0, DOF *u, DOF *p, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g, INT *nits_amg)//, double *time_amg)
{
	 VEC *x, *y, *TBy, *d1, *d2, *fr, *gr;
	 MAT *H;
	 SOLVER *solver1, *solver2, *solver3;
	 FLOAT value = 0., b1=0.0, b2=0.0,f=0.0, g=0.0, res1 =0., error = 1e-8, res = 0, res_f = 0, res_g = 0;
	 INT inter_amg = 0, Nits = 0, inter_pcg = 0, inter_max = 0;
	 int i;
	 SOLVER *pc;

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
		  solver1 = phgMat2Solver(SOLVER_DEFAULT, A);
		  solver1->mat->handle_bdry_eqns = FALSE;

		  phgVecCopy(fr, &(solver1->rhs));
//		  phgPrintf("vec_f->nec = %d, sol->rhs->map->ndof = %d\n", vec_f->nvec, solver1->rhs->map->ndof);
		  phgSolverVecSolve(solver1, FALSE, d1);
		  phgVecAXPBY(1.0, d1, 1.0, &x);
		  inter_pcg = solver1->nits;
		  phgSolverDestroy(&solver1);
		  // g = Bx-Cy-g
		  phgVecCopy(vec_g, &gr);
		  phgMatVec(0, 1.0, B, x, -1.0, &gr);
		  phgMatVec(0, -1.0, C, y, 1.0, &gr);
		  res_g = phgVecNorm2(gr, 0, NULL);
		  solver2 = phgMat2Solver(SOLVER_DEFAULT, H);
		  phgVecCopy(gr, &(solver2->rhs));
				phgOptionsPush();
				phgOptionsSetOptions("-solver hypre -hypre_solver boomeramg");
				phgOptionsSetOptions("-hypre_amg_coarsen_type pmis");
				phgOptionsSetOptions("-hypre_amg_coarsest_relax_type gs");
				pc = phgMat2Solver(SOLVER_DEFAULT, H);
				pc->mat->handle_bdry_eqns = FALSE;
				pc->maxit = 2;
				pc->warn_maxit = FALSE;
				phgOptionsPop();
		 		phgSolverSetPC(solver2, pc, pc_proc);
	//	  elapsed_time(u->g, FALSE, 0.);
		  phgSolverVecSolve(solver2, FALSE, d2);
	//	  *time_amg += elapsed_time(u->g, FALSE, 0.);
		  *nits_amg += solver2->nits;

		  phgSolverDestroy(&solver2);
		  phgSolverDestroy(&pc);

		  res1 = phgVecNorm2(d2, 0, NULL);
		  if (res1 != 0.0 ){
		 		b1 = phgVecDot(d2, 0, gr, 0, NULL);
				for (i = 0; i < TBy->map->nlocal; i++){
					TBy->data[i] = 0.0;
				}
				phgMatVec(0, 1.0, TB, d2, 0.0, &TBy);
				solver3 = phgMat2Solver(SOLVER_DEFAULT, A);
			    solver3->mat->handle_bdry_eqns = FALSE;
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
		  f = phgVecNorm2(vec_f, 0, NULL);
		  if (f == 0){
		  	res_f = res_f;
		  }else{
		  	res_f = res_f / f;
		  }
		  g = phgVecNorm2(vec_g, 0, NULL);
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

	 return Nits;
}
FLOAT 
uzawapcg(MAT *H, DOF *u_h, DOF *p_h, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g, INT *nits_amg)
{
		VEC *x, *y, *r, *tempr, *tempp, *Br, *p, *Ap;
		SOLVER *sol_1, *sol_2, *sol_3, *sol_4, *pc;
		FLOAT tol = 1e-8, tolA = 0.5 * 1e-8, ng = 0, err = 0, rho = 0, rho_old = 0, beta = 0, alpha = 0;
		INT k = 1;
		x = phgMapCreateVec(map_u, 1);
		tempr = phgMapCreateVec(map_u, 1);
		tempp = phgMapCreateVec(map_u, 1);

		y = phgMapCreateVec(map_p, 1);
		r = phgMapCreateVec(map_p, 1);
		p = phgMapCreateVec(map_p, 1);

		Br = phgMapCreateVec(map_p, 1);
		Ap = phgMapCreateVec(map_p, 1);
		
		phgMapDofToLocalData(map_u, 1, &u_h, x->data);
		phgMapDofToLocalData(map_p, 1, &p_h, y->data);
		phgMatAXPBY(1.0, C, 1.0, &H);
		/*solve tempr*/
		sol_1 = phgMat2Solver(SOLVER_PCG, A);
		sol_1->mat->handle_bdry_eqns = FALSE;
		phgVecCopy(vec_f, &(sol_1->rhs));
		phgMatVec(0, -1.0, TB, y, 1.0, &(sol_1->rhs));
		phgSolverVecSolve(sol_1, FALSE, tempr);
		phgSolverDestroy(&sol_1);
		/*end of solve temp_r*/
		phgMatVec(0, -1.0, NULL, vec_g, 0.0, &r);
		phgMatVec(0.0, 1.0, B, tempr, 1.0, &r);
		phgMatVec(0.0, -1.0, C, y, 1.0, &r);
		ng = phgVecNorm2(vec_g, 0, NULL);
		err = phgVecNorm2(r, 0, NULL) / ng;
		while (err > tol)
		{
			sol_2 = phgMat2Solver(SOLVER_PCG, H);
			phgVecCopy(r, &(sol_2->rhs));
			phgOptionsPush();
			phgOptionsSetOptions("-solver hypre -hypre_solver boomeramg");
			phgOptionsSetOptions("-hypre_amg_coarsen_type pmis");
			phgOptionsSetOptions("-hypre_amg_coarsest_relax_type gs");
			pc = phgMat2Solver(SOLVER_DEFAULT, H);
			pc->mat->handle_bdry_eqns = FALSE;
			pc->maxit = 2;
			pc->warn_maxit = FALSE;
			phgOptionsPop();
			phgSolverSetPC(sol_2, pc, pc_proc);
			phgSolverVecSolve(sol_2, FALSE, Br);

		//	if (*nits_amg < sol_2->nits){
				*nits_amg += sol_2->nits;
		//	}
			phgSolverDestroy(&sol_2);
			phgSolverDestroy(&pc);
			rho = phgVecDot(Br, 0, r, 0, NULL);
			if (k == 1)
				phgVecCopy(Br, &p);
			else{
				beta = rho / rho_old;
				phgVecAXPBY(1.0, Br, beta, &p);
			}
			sol_3 = phgMat2Solver(SOLVER_PCG, A);
			phgMatVec(0.0, 1.0, TB, p, 0.0, &(sol_3->rhs));
			phgSolverVecSolve(sol_3, FALSE, tempp);
			phgSolverDestroy(&sol_3);
			phgMatVec(0.0, 1.0, B, tempp, 0.0, &Ap);
			phgMatVec(0.0, 1.0, C, p, 1.0, &Ap);
			alpha = rho / phgVecDot(Ap, 0, p, 0, NULL);
			phgVecAXPBY(-alpha, Ap, 1.0, &r);
			phgVecAXPBY(alpha, p, 1.0, &y);
			rho_old = rho;
			k = k + 1;
			err = phgVecNorm2(r, 0, NULL) / ng;
		}
		sol_4 = phgMat2Solver(SOLVER_PCG, A);
		sol_4->mat->handle_bdry_eqns = FALSE;
		phgVecCopy(vec_f, &(sol_4->rhs));
		phgMatVec(0, -1.0, TB, y, 1.0, &(sol_4->rhs));
		phgSolverVecSolve(sol_4, FALSE, x);
		phgSolverDestroy(&sol_4);

		phgMapLocalDataToDof(map_u, 1, &u_h, x->data);
		phgMapLocalDataToDof(map_p, 1, &p_h, y->data);
		phgVecDestroy(&x);
		phgVecDestroy(&y);
		phgVecDestroy(&r);
		phgVecDestroy(&tempr);
		phgVecDestroy(&tempp);
		phgVecDestroy(&Br);
		phgVecDestroy(&p);
		phgVecDestroy(&Ap);
		
		return k;
}
