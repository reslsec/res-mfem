 /*
 * This sample code solves the Singlephase equation:
 * ************************************************
 * problem: p_{t} -\Delta{p} = func_f(x, t) \in \Omega X [0, T]
 *          p = func_g (x, t) \in \partial\Omega X [0, T]
 *          p = func_p0       \in \Omega t = 0.0
 *WARNING! The unit is h, m, mpa!!
 * 
 * ***********************************************
 */
#include "phg.h"
#include <string.h>
#include <math.h>
#include "well.c"
#include "uzawa.c"
#define USE_UZAWA 0
#define USE_BLOCK 1
static DOF_TYPE DOF_RT1_;
#define DOF_RT1 (&DOF_RT1_)
# define BasisOrder(u, e, i) (!DofIsHP(u) ? (u)->type->order : \
	(u)->hp->info->types[(u)->hp->elem_order[e->index]]->order)
/* some parameters */
static FLOAT T =720.;
static FLOAT stime;    			//time step size
static FLOAT ctime = 0.0;      		//current time 
static FLOAT ptime = 0.0;     		//previous time 
/*some parameters for phase flow*/
static FLOAT mu =1e-9 / 3600.;          //the viscosity
static FLOAT b0 = 1.05;           		//the init of reservoir volume ratio (bo)
static FLOAT co = 2e-3;          
static FLOAT K = 1e-13;          		//absolute permeability
static FLOAT cr = 0.15e-3;        
static FLOAT phi0 = 0.2;          		//porosity 
static FLOAT PRESSURE0 = 20.;           //init pressure
static FLOAT rw = 0.1;
static FLOAT flow_rate = 10. / 24.;

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
/* build linear system */
static void
build_mat_vec(DOF *u_h, DOF *p_h, DOF *p_nk, DOF *q_o, DOF *b_o, DOF *phi, DOF *phi_l, DOF *b_o_l, MAP *map_u, MAP *map_p, MAT *A, MAT *B, MAT *TB, MAT *C, VEC *vec_f, VEC *vec_g)

{
	 int i, j, n,k;
	 int N = u_h->type->nbas * u_h->dim;
	 int M = p_h->type->nbas * p_h->dim;
	 GRID *g = p_h->g;
	 SIMPLEX *e;
	 assert(p_h->dim == 1);
	 phgVecDisassemble(vec_f);
	 phgVecDisassemble(vec_g);
	 INT I[N],J[N],L[M];
	 FLOAT A0[N][N], TB0[N][M], B0[M][N],C0[M][M],rhs_f[N],rhs_g[M];
	 FLOAT *p_phi, *p_bo, *p_bol;
	 ForAllElements(g,e) {
		p_phi = DofElementData(phi, e->index);
		p_bo = DofElementData(b_o, e->index);
		p_bol = DofElementData(b_o_l, e->index);
		for (i = 0; i < N; i++){
			 for (j = 0; j < N; j++){
			 	 A0[i][j] = stime * mu * phgQuadBasDotBas(e, u_h, i, u_h, j, QUAD_DEFAULT) / (K * p_bo[0]);
			 }
		}
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				TB0[i][j]= -stime * phgQuadDivBasDotBas(e, u_h, i, p_h, j, QUAD_DEFAULT);
			}
		}       
		FLOAT c_p = p_phi[0] * co / b0 + p_bo[0] * cr * phi0;
		 for (i = 0; i < M; i++){
		 	for (j = 0; j< M; j++){
		  		C0[i][j] = c_p * phgQuadBasDotBas(e, p_h, i, p_h, j, QUAD_DEFAULT);
		  	}
		}
		for(i = 0; i < N; i++){
			rhs_f[i] = 0;
		}
		FLOAT quad_phi = 0, quad_pnk = 0, quad_phil = 0, quad_qo = 0;
		for (i = 0; i < M; i++){
			phgQuadDofTimesBas(e, q_o, p_h, i, QUAD_DEFAULT, &quad_qo);
			phgQuadDofTimesBas(e, phi_l, p_h, i, QUAD_DEFAULT, &quad_phil);
			phgQuadDofTimesBas(e, phi, p_h, i, QUAD_DEFAULT, &quad_phi);
			phgQuadDofTimesBas(e, p_nk, p_h, i, QUAD_DEFAULT, &quad_pnk);
			rhs_g[i] = -quad_qo * stime + p_bo[0] * quad_phi - c_p * quad_pnk - p_bol[0] * quad_phil;
		}
		/* Handle Bdry Conditions */
		for (i = 0; i < N; i++){
			if (phgDofGetElementBoundaryType(u_h, e, i) & (NEUMANN | DIRICHLET)){
				bzero(A0[i], N *sizeof(A0[i][0]));
				bzero(TB0[i], M *sizeof(TB0[i][0]));
				for (j = 0; j < N; j++){
					A0[j][i] = 0.0;
				}
				A0[i][i] = 1.0;
				rhs_f[i] = 0.;
			}
		}
		for (i = 0; i < N; i++){
			for (j = 0; j < M; j++){
				B0[j][i] = TB0[i][j];
			}
		}
		for (i = 0; i < N; i++){
			I[i] = phgMapE2L(map_u, 0, e, i);
		}
		for (i = 0; i < M; i++){
			J[i] = phgMapE2L(map_p, 0, e, i);
		}
		phgMatAddEntries(A, N, I, N, I, A0[0]);
		phgMatAddEntries(TB, N, I, M, J, TB0[0]);
		phgMatAddEntries(B, M, J, N, I, B0[0]);
		phgMatAddEntries(C, M, J, M, J, C0[0]);
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
/* cal L2 norm for DOF_ANALYTIC */
static FLOAT
L2_norm(DOF *u)
{
    SIMPLEX *e;
    FLOAT tmp, result;

    tmp = 0.0;
    ForAllElements(u->g, e){
            tmp += phgQuadDofDotDof(e, u, u, 5);
        }
#if USE_MPI
    MPI_Allreduce(&tmp, &result, 1, PHG_MPI_FLOAT, MPI_SUM, u->g->comm);
#else
    result = tmp;
#endif
    return Sqrt(result);
}

/* L2 error for true and numerical solution  */
static FLOAT
L2(DOF *a, DOF *b, int order)
{
    GRID *g = a->g;
    SIMPLEX *e;
    QUAD *quad = phgQuadGetQuad3D(order);
    FLOAT v;
    int i, n;
    FLOAT *pa, *pb;

    assert(DofDim(a) == DofDim(b));

    FLOAT *vala = phgAlloc(quad->npoints * DofDim(a) * sizeof(*vala));
    FLOAT *valb = phgAlloc(quad->npoints * DofDim(b) * sizeof(*valb));

    FLOAT sum = 0.;
    ForAllElements(g, e) {
	memcpy(vala, phgQuadGetDofValues(e, a, quad),
	       quad->npoints * DofDim(a) * sizeof(*vala));
	memcpy(valb, phgQuadGetDofValues(e, b, quad),
	       quad->npoints * DofDim(b) * sizeof(*valb));

	const FLOAT *w = quad->weights;
	v = 0.;
	pa = vala;
	pb = valb;
	for (n = 0; n < quad->npoints; n++) {
	    FLOAT v0 = 0.;
	    for (i = 0; i < DofDim(a); i++) {
		v0 += (*pa - *pb) * (*pa - *pb);
		pa++;
		pb++;
	    }
	    v += v0 * (*w);
	    w++;
	}

	v *= phgGeomGetVolume(g, e);
	sum += v;
    }

#if USE_MPI
    FLOAT s = sum;
    MPI_Allreduce(&s, &sum, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
#endif

    phgFree(vala);
    phgFree(valb);

    return Sqrt(sum);
}
static void
update_bo(DOF *p_h, DOF *p0_h, DOF *bo)
{

	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_bo, *p_p0;

	ForAllElements(g, e){
		p_p0 = DofElementData(p0_h, e->index);
		p_p = DofElementData(p_h, e->index);
		p_bo = DofElementData(bo, e->index);
		p_bo[0] = (1.0 + co * (p_p[0] - p_p0[0])) / b0;
	}
}

static void
update_phi(DOF *p_h, DOF *p0_h, DOF *phi)
{

	GRID *g = p_h->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_phi, *p_p0;

	ForAllElements(g, e){
		p_p0 = DofElementData(p0_h, e->index);
		p_p = DofElementData(p_h, e->index);
		p_phi = DofElementData(phi, e->index);
		p_phi[0] = phi0 * (1.0 + cr * (p_p[0] - p_p0[0]));
	}
}

static void
conservation(DOF *a, DOF *b, DOF *well)
{
	GRID *g = a->g;
	DOF *tmp, *tmp2;
	tmp = phgDofNew(g, DOF_P0, 1, "tmp", DofNoAction);
	phgDofSetDataByValue(tmp, 1.0);
	tmp2 = phgDofCopy(a, NULL, NULL, NULL);
	phgDofAXPBY(-1.0, b, 1.0, &tmp2);

	SIMPLEX *e;
	FLOAT output = 0, input = 0;
	FLOAT v1 =0, v2=0;
	ForAllElements(g, e){
		v1 += phgQuadDofDotDof(e, well, tmp, 5);
		v2 += phgQuadDofDotDof(e, tmp2, tmp, 5)/stime;
	}
#if USE_MPI
    MPI_Allreduce(&v1, &input, 1, PHG_MPI_FLOAT, MPI_SUM, a->g->comm);
    MPI_Allreduce(&v2, &output, 1, PHG_MPI_FLOAT, MPI_SUM, a->g->comm);
#else
    input = v1;
    output = v2;
#endif
	phgPrintf("----------------Total Conservation---------------------\n");
	phgPrintf("The Output flow is %le, The Input Flow is %le\n", output, input);
	phgDofFree(&tmp);
	phgDofFree(&tmp2);
}
static FLOAT
Get_Diameter(DOF *u)
{
	GRID *g = u->g;
	SIMPLEX *e;
	FLOAT h = 1e-9, h_max = 1e-8;
	ForAllElements(g, e){
		h = phgGeomGetDiameter(g, e);
//		phgPrintf("h = %lf\n", h);
		if (h > h_max){
			h_max = h;
//			phgPrintf("h_max = %lf\n", h_max);
		}
	}
	return h_max;
}
int
main(int argc, char *argv[])
{
    INT mem_max = 10240, N = 0;
    FLOAT tol = 1e-06;
  //  char *fn = "flow.mesh";
   char *fn = "single_level.mesh";
//	char *fn = "one-million.mesh";
    char vtkfile[1000];
    GRID *g;
    MAT *A, *B, *TB, *C;
    VEC *vec_f, *vec_g;
    MAP *map_u, *map_p;
    SOLVER *pc;
    static int pre_refines = 0;
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterFloat("tol", "Tolerance", &tol);
    phgOptionsRegisterFloat("step", "Time step", &stime);
    phgOptionsRegisterFloat("T", "Time domain", &T);
	phgOptionsRegisterFilename("mesh_file", "Mesh File", (char **)&fn);
    phgOptionsRegisterInt("mem_max", "Maximum memory (MB)", &mem_max);
//	phgOptionsRegisterKeyword("hypre_amg_coarsest_relax_type", "Coarsest grid solver",
//	gs-h-forward, &global_params.amg_coarsest_relax_type);

    phgInit(&argc, &argv);
    g = phgNewGrid(-1);
	phgOptionsShowUsed();

    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
		if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
				(double)g->lif);
	phgRefineAllElements(g, pre_refines);
#if 0
    if(!phgExportMedit(g, ff))
			            phgError(1, "can't write file \"%s\".\n",ff);
	exit(1);
#endif 
    /*The pressure variable*/
    DOF *p_h;
    p_h = phgDofNew(g, DOF_P0, 1, "p_h", DofNoAction);
    phgDofSetDataByValue(p_h, PRESSURE0);
    DOF *p0_h;
    p0_h = phgDofNew(g, DOF_P0, 1, "p0_h", DofNoAction);
    phgDofSetDataByValue(p0_h, PRESSURE0);
   /*The velocity variable*/
    DOF *u_h;
    u_h = phgDofNew(g, DOF_RT1, 1, "u_h", DofNoAction);
    /* RHS function */
    DOF *well;

    well = phgDofNew(g, DOF_P0, 1, "well", DofNoAction);
    Well_init(well);
    FLOAT h = 0;
	h = Get_Diameter(p_h);
	int flag = 0;
	double total_time = 0;
	INT total_uzawa =0, total_amg=0, newton=0;
	phgPrintf("the elements is :%d\n",DofGetDataCountGlobal(p_h));
	phgPrintf("the DOF number is :%d\n",DofGetDataCountGlobal(p_h)+ DofGetDataCountGlobal(u_h));
	phgPrintf("The max diameter is: %lf\n", h);
	while (ctime < T -1e-8 ){
		ptime = ctime;
		ctime += stime;
		phgPrintf("\n/*****start new time layer *****/\n");
		phgPrintf("current time layer: [%lf, %lf]\n", (double)ptime, (double)ctime);
		if (ctime > T){ 
			ctime = T;
			stime = T- ptime;
			phgPrintf("current time layer : [%lf, %lf]\n",(double)ptime, (double)ctime);
		}
		/*For Newton Interation*/
		DOF *p_nk, *u_nk;
		int count = 0;
 
//		elapsed_time(g, FALSE, 0.);
		/* The parameter DOF */
 	   	DOF *b_o_l, *phi_l;
		phi_l = phgDofNew(g, DOF_P0, 1, "phi_l", DofNoAction);
		update_phi(p_h, p0_h, phi_l);
 	    b_o_l = phgDofNew(g, DOF_P0, 1, "b_o_l", DofNoAction);
		update_bo(p_h, p0_h, b_o_l);
    	double time = 0; // total time AMG use in T;
		INT NITS_amg = 0; //total NITS use in T;
		while (TRUE)
		{	
			count++;
			
			/*Create MAP for Mat and Vec*/
	     	map_p = phgMapCreate(p_h, 0);
			map_u = phgMapCreate(u_h, 0);
			A = phgMapCreateMat(map_u, map_u);
			B = phgMapCreateMat(map_p, map_u);
			TB = phgMapCreateMat(map_u, map_p);
			C = phgMapCreateMat(map_p, map_p);
			vec_f = phgMapCreateVec(map_u, 1);
			vec_g = phgMapCreateVec(map_p, 1);
			
     		/* The parameter DOF */
 	  	 	DOF *phi, *b_o;
			phi = phgDofNew(g, DOF_P0, 1, "phi", DofNoAction);
			update_phi(p_h, p0_h, phi);
 		   	b_o = phgDofNew(g, DOF_P0, 1, "b_o", DofNoAction);
			update_bo(p_h, p0_h, b_o);	
		
			p_nk = phgDofCopy(p_h, NULL, NULL, NULL);
			u_nk = phgDofCopy(u_h, NULL, NULL, NULL);
		elapsed_time(g, FALSE, 0.);
			build_mat_vec(u_h, p_h, p_nk, well, b_o, phi, phi_l, b_o_l, map_u, map_p, A, B, TB, C, vec_f, vec_g);			
			phgPrintf("Build linear:        ");
			elapsed_time(g, TRUE, 0.);
			phgDofFree(&phi);
			phgDofFree(&b_o);
#if USE_UZAWA
		int nits_amg = 0;
        int nits_uzawa = 0;
		double time_amg = 0;
		elapsed_time(g, FALSE, 0.);
		DOF *B_data = DivRT(p_h, u_h);
        MAT *H = BTAB(A, B_data, p_h, u_h);
      	phgDofFree(&B_data);
		phgPrintf("Assemble H:             ");
		elapsed_time(g, TRUE, 0.);
		phgPrintf("solve p_h:              \n");
//		nits_uzawa = phgUzawa(H, u_h, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg, &time_amg);
		elapsed_time(g, FALSE, 0.);
		nits_uzawa = uzawapcg(H, u_h, p_h, map_u, map_p, A, B, TB, C, vec_f, vec_g, &nits_amg);
		total_time += elapsed_time(g, TRUE, 0.);
		total_uzawa += nits_uzawa;
		total_amg += nits_amg;
		phgPrintf("Total AMG for once newton---------%d\n", nits_amg);
		phgPrintf("Nits: uzawa:----------------%d\n", nits_uzawa);
		phgMatDestroy(&H);
#endif
#if USE_BLOCK
		SOLVER *solver;
		MAT *pmat[4];
		FLOAT coef[4];
		solver = phgSolverCreate(SOLVER_GMRES, u_h, p_h, NULL);
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
				phgOptionsPush();
		   		phgOptionsSetOptions("-solver petsc "
							     "-oem_options -ksp_type gmres "
								 "-oem_options -pc_type asm "													           "-oem_options -pc_asm_overlap 2 "
							     "-oem_options -sub_ksp_type preonly "													   "-oem_options -sub_pc_type lu "															 "-oem_options -sub_pc_factor_mat_solver_package mumps "								   "-oem_options -mat_mumps_sym 2 "
								 "-oem_options -mat_mumps_icntl_4 0");
				pc = phgMat2Solver(SOLVER_DEFAULT, solver->mat);
				pc->mat->handle_bdry_eqns = FALSE;
				phgOptionsPop();
		elapsed_time(g, FALSE, 0.);
		phgSolverSolve(solver, TRUE, u_h, p_h, NULL);
		phgPrintf("Solve P:  iter = %d", solver->nits);
		elapsed_time(g, TRUE, 0.);
		phgSolverDestroy(&solver);
#endif
			/*compute the relative error of ||p^{n,k}-p^{n,k-1}||*/
			FLOAT err2 = 0.0, norm2 = 0.0, err1 =0.0 , norm1 = 0.0, TOL = 0.0;
			err1 = L2(u_nk, u_h, 4);
			err2 = L2(p_nk, p_h, 4);
			norm1 = L2_norm(u_nk);
			norm2 = L2_norm(p_nk);
			TOL = err2 / norm2;// + err1 / norm1;
			phgDofFree(&u_nk);
			phgDofFree(&p_nk);
			phgMatDestroy(&A);
			phgMatDestroy(&B);
			phgMatDestroy(&TB);
			phgMatDestroy(&C);	
			phgMapDestroy(&map_p);
			phgMapDestroy(&map_u);
			phgVecDestroy(&vec_g);		
			phgVecDestroy(&vec_f);
			if(TOL < 1e-4){
				phgPrintf("----------------------------------------------\n");	
				phgPrintf("the nonlinear interation number is %d\n", count);
				phgPrintf("----------------------------------------------\n\n");	
				newton += count;
				break;
			}
		}
//		elapsed_time(g, TRUE, 0.);
		DOF *b_o;
 		b_o = phgDofNew(g, DOF_P0, 1, "bo", DofNoAction);
		update_bo(p_h, p0_h, b_o);
		FLOAT pwf = 0, pressure =0;
		Well_Pressure(p_h, b_o, &pressure, &pwf);
		phgPrintf("t =  %lf,  Pwf =  %lf,  Peaceman = %lf\n", ctime, pressure, pwf);
		phgPrintf("Total time :------------------------%lf\n", total_time);
		phgPrintf("Total uzawa:------------------------%d\n", total_uzawa);
		phgPrintf("total amg  :------------------------%d\n", total_amg);
		phgPrintf("     newton:------------------------%d\n", newton);
		phgDofFree(&phi_l);
		phgDofFree(&b_o_l);
		phgDofFree(&b_o);
	}
	
		sprintf(vtkfile, "single_phase_%s.vtk", "finaltime");
		phgExportVTK(g, vtkfile, p_h, well, u_h, NULL);
    phgDofFree(&p_h);
    phgDofFree(&p0_h);
    phgDofFree(&u_h);
    phgSolverDestroy(&pc);
	phgDofFree(&well);	
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

