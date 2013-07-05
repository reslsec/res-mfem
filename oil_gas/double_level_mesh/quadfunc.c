#include "phg.h"
#include <string.h>
#include <math.h>
//#include "oil_gas_phase.h"
FLOAT K_1[3] = {1e-13, 1e-13, 1e-13};          		//absolute permeability
FLOAT K_2[3] = {2e-13, 2e-13, 1e-13};          		//absolute permeability
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
FLOAT
phgQuadKBasDotBas(SIMPLEX *e, DOF *u, int n, DOF *v, int m, int order)
{
    int i, j, nvalues = DofTypeDim(u);
    FLOAT *g1, *g2, *w;
    FLOAT d, d0;
    QUAD *quad;
    static FLOAT K[3] = {0.0, 0.0, 0.0};
    if (nvalues != DofTypeDim(v))
	phgError(1, "%s:%d: dimensions mismatch: %s (%d) <==> %s ($d))\n",
		 __FILE__, __LINE__, DofTypeName(u), DofTypeDim(u),
		 DofTypeName(v), DofTypeDim(v));
//    phgPrintf("Dim(u) = %d\n", DofTypeDim(u));
    if (order < 0)
		order = 5;
    if(e->region_mark == 0 | e->region_mark == 2){
    		for (i = 0; i < nvalues; i++){
			K[i] = K_1[i];
		}
    }
    else if(e->region_mark == 1 | e->region_mark == 3){
		for (i = 0; i < nvalues; i++){
			K[i] = K_2[i];
		}
    }
    quad = phgQuadGetQuad3D(order);
    g1 = phgQuadGetBasisValues(e, u, n, quad);
    g2 = phgQuadGetBasisValues(e, v, m, quad);
    d = 0.;
    w = quad->weights;
    for (i = 0; i < quad->npoints; i++) {
	d0 = 0.;
	for (j = 0; j < nvalues; j++) {
	    d0 += *(g1++) * (*(g2++)) / K[j];
	}
	d += d0 * (*(w++));
    }
    return d * phgGeomGetVolume(u->g, e);
}

FLOAT
phgQuadDofDofDof(SIMPLEX *e, DOF *u1, DOF *u2, DOF *u3, int order)
{
    int i, j, nvalues;
    FLOAT d, d0, K = 0;
    const FLOAT *v1, *v2, *v3, *w;
    QUAD *quad;

    nvalues = DofDim(u1);

    if (order < 0) {
	i = DofTypeOrder(u1, e);
	j = DofTypeOrder(u2, e);
	if (i < 0 && j < 0) {
	    phgInfo(-1, "phgQuadDofDotDof: don't use QUAD_DEFAULT when both "
		    "DOF types are analytic.\n");
	    phgError(1, "phgQuadDofDotDof: can't determine quadrature order, "
			"abort.\n");
	}
	if (i < 0)
	    i = j;
	else if (j < 0)
	    j = i;
	order = i + j;
    }
    quad = phgQuadGetQuad3D(order);

    w = quad->weights;
    v1 = phgQuadGetDofValues(e, u1, quad);
    v2 = phgQuadGetDofValues(e, u2, quad);
    v3 = phgQuadGetDofValues(e, u3, quad);

    if (nvalues == 1) {		/* faster code for special case */
	d = 0.;
	for (i = 0; i < quad->npoints; i++)
	    d += *(v1++) * *(v2++) * *(v3++) * *(w++);
    }
    return d * phgGeomGetVolume(u1->g, e);
}
/*this codes for (div u, p)*/
FLOAT
phgQuadDivBasDotBas(SIMPLEX *e, DOF *u, int m, DOF *v, int n, int order)
{
    int i, j, nvalues = DofTypeDim(v);
    const FLOAT *g1, *g2, *w, *lambda;
    FLOAT d, d0;
    QUAD *quad;

    assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
    if (3*nvalues !=  DofTypeDim(u))
	phgError(1, "%s:%d, dimensions mismatch: grad(%s) <==> (%s)\n",
		 __FILE__, __LINE__, u->name, v->name);

	order = 4;
    quad = phgQuadGetQuad3D(order);

    g1 = phgQuadGetBasisGradient(e, u, m, quad);
    g2 = phgQuadGetBasisValues(e, v, n, quad);
    d = 0.;
    lambda = quad->points;
    w = quad->weights;
    for (i = 0; i < quad->npoints; i++) {
	d0 = 0.;
	for (j = 0; j < nvalues; j++) {
	    d0 += (*(g1)+(*(g1+4))+(*(g1+8))) * (*(g2++));
	    g1=g1+9;
	}
	d += d0 * (*(w++));
	lambda += Dim + 1;
    }
    return d * phgGeomGetVolume(u->g, e);
}
