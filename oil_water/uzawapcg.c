#include "phg.h"
#include "math.h"

FLOAT uzawapcg(MAT *H0, DOF *u_h, DOF *p_h, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g, INT *nits_amg, double *time_amg)
{
		VEC *x, *y, *f, *r, *tempr, *tempp, *Br, *p, *Ap;
		MAT *H;
		SOLVER *sol_1, *sol_2, *sol_3, *sol_4, *pc;
		FLOAT tol = 1e-8, tolA = 0.5 * 1e-8, ng = 0, err = 0, rho = 0, rho_old = 0, beta = 0, alpha = 0;
		INT k = 1;
		x = phgMapCreateVec(map_u, 1);
		f = phgMapCreateVec(map_u, 1);
		tempr = phgMapCreateVec(map_u, 1);
		tempp = phgMapCreateVec(map_u, 1);

		y = phgMapCreateVec(map_p, 1);
		r = phgMapCreateVec(map_p, 1);
		Br = phgMapCreateVec(map_p, 1);
		p = phgMapCreateVec(map_p, 1);
		H = phgMapCreateMat(map_p, map_p);
		
		phgMapDofToLocalData(map_u, 1, &u_h, x->data);
		phgMapDofToLocalData(map_p, 1, &p_h, y->data);
		H = H0;
		phgMatAXPBY(1.0, C, 1.0, &H);
		/*solve tempr*/
		phgVecCopy(vec_f, &f);
		phgMatVec(0, -1.0, TB, y, 1.0, &f);
		sol_1 = phgMat2Solver(SOLVER_PCG, A);
		sol_1->mat->handle_bdry_eqns = FALSE;
		phgVeccopy(f, &(sol_1->rhs));
		phgSolverSolve(sol_1, FALSE, tempr);
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
			use_time(u_h->g, FALSE, 0.);
			phgSolverSolve(sol_2, FALSE, Br);
			*time_amg += use_time(u_h->g, FALSE, 0.);
			*nits_amg += sol_2->nits;
			phgSolverDestroy(&sol_2);
			phgSolverDestroy(&pc);
			rho = phgVecDot(Br, 0, r, 0, NULL);
			if (k ==1)
				phgVecCopy(Br, &p);
			else{
				beta = rho / rho_old;
				phgVecAXPBY(1.0, Br, beta, &p);
			}
			sol_3 = phgMat2Solver(SOLVER_PCG, A);
			phgMatVec(0.0, 1.0, TB, p, 0.0, &(sol_3->rhs));
			phgSolverSolve(sol_3, FALSE, tempp);
			phgSolverDestroy(&sol_3);
			phgMatVec(0.0, 1.0, B, tempp, 0.0, &Ap);
			phgMatVec(0.0, 1.0, C, p, 1.0, &Ap);
			alpha = rho / phgVecDot(Ap, 0, p, 0, NULL);
			phgVecAXPBY(-alpha, Ap, 1.0, r);
			phgVecAXPBY(alpha, p, 1.0, y);
			rho_old = rho;
			k = k+1;
			err = phgVecNorm2(r, 0, NULL) / ng;
		}
		phgVecCopy(vec_f, &f);
		phgMatVec(0, -1.0, TB, y, 1.0, &f);
		sol_4 = phgMat2Solver(SOLVER_PCG, A);
		sol_4->mat->handle_bdry_eqns = FALSE;
		phgVeccopy(f, &(sol_4->rhs));
		phgSolverSolve(sol_4, FALSE, x);
		phgSolverDestroy(&sol_4);

		phgMapLocalDataToDof(map_u, 1, &u_h, x->data);
		phgMapLocalDataToDof(map_p, 1, &p_h, y->data);
		phgVecDestroy(&x);
		phgVecDestroy(&y);
		phgVecDestroy(&f);
		phgVecDestroy(&r);
		phgVecDestroy(&tempr);
		phgVecDestroy(&tempp);
		phgVecDestroy(&Br);
		phgVecDestroy(&p);
		phgVecDestroy(&Ap);
		
		return k;
}
