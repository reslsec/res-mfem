#include "phg.h"
#include <string.h>
#include <math.h>
//#include "ss_oil_water.h"
FLOAT
phgQuadFacefunc(SIMPLEX *e, int face, DOF *u, int order)
{
    GRID *g = u->g;
    int i,v0, v1, v2;
    FLOAT d,lambda[Dim + 1];
    FLOAT *dof, *buffer;
    const FLOAT *p, *w;
    QUAD *quad;

    assert(face >= 0 && face <= 3);

    if (order < 0) {
	i = DofTypeOrder(u, e);
	order = i;
    }
    quad = phgQuadGetQuad2D(order);

    v0 = GetFaceVertex(face, 0);
    v1 = GetFaceVertex(face, 1);
    v2 = GetFaceVertex(face, 2);
    lambda[face] = 0.;

    buffer = phgAlloc(sizeof(*buffer));
    p = quad->points;
    w = quad->weights;

    d = 0.;
    for (i = 0; i < quad->npoints; i++) {
         lambda[v0] = *(p++);
         lambda[v1] = *(p++);
         lambda[v2] = *(p++);
           
         dof = phgDofEval(u, e, lambda, buffer);
	    d += *(dof)* *(w++);
	}

    phgFree(buffer);

    return d * phgGeomGetFaceArea(g, e, face);
}
FLOAT
phgQuadDofTimesDivBas(SIMPLEX *e, DOF *u, DOF *v, int m, int order)
{
    int i, j, k, nvalues = DofDim(u);
    const FLOAT *g1, *g2, *w;
    FLOAT d, d0;
    QUAD *quad;

    assert(!SpecialDofType(v->type));
    assert(nvalues * Dim == DofTypeDim(v));

    if (order < 0) {
	order = DofTypeOrder(u, e);
    }
    quad = phgQuadGetQuad3D(order);

    d = 0.;
    g1 = phgQuadGetDofValues(e, u, quad);
    g2 = phgQuadGetBasisGradient(e, v, m, quad);
    w = quad->weights;

	for (i = 0; i < quad->npoints; i++) {
	    d0 = 0.;
	    for (j = 0; j < nvalues; j++) {
		d0 += *(g1++) * (*(g2) + *(g2+4) + *(g2+8));
		g2 = g2 +9;
	    }
	    d += d0 * (*(w++));
	}
	return d * phgGeomGetVolume(u->g, e);
}
/* cal L2 norm for DOF_ANALYTIC */
FLOAT
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
FLOAT
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
double
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
