#include "oil_water.h"
#include "phg.h"
/*This quad func for (u,K, A,v) / (D1 * D2)*/
FLOAT K[3] = {1e-13, 1e-13, 1e-13};
FLOAT
phgQuadBasKDDDBas(SIMPLEX *e, DOF *u, int n, DOF *D1, DOF *D2, DOF *D3, DOF *v, int m, int order)
{
    int i, j, k, nvalues = DofTypeDim(u);
    FLOAT *g1, *g2, *coef1, *coef2, *coef3, *w;
    FLOAT d, d0;
    QUAD *quad;

    assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));

    if (nvalues != DofTypeDim(v))
	phgError(1, "%s:%d: dimensions mismatch: %s (%d) <==> %s ($d))\n",
		 __FILE__, __LINE__, DofTypeName(u), DofTypeDim(u),
		 DofTypeName(v), DofTypeDim(v));
	if (order < 0)
		order = 4;
    quad = phgQuadGetQuad3D(order);

    g1 = phgQuadGetBasisValues(e, u, n, quad);
    g2 = phgQuadGetBasisValues(e, v, m, quad);
    w = quad->weights;
	    coef1 = phgQuadGetDofValues(e, D1, quad);
	    coef2 = phgQuadGetDofValues(e, D2, quad);
		coef3 = phgQuadGetDofValues(e, D3, quad);
		for (i = 0; i < quad->npoints; i++) {
		    d0 = 0.;
		    for (j = 0; j < nvalues; j++) {
			d0 += *(g1++) * *(g2++) / K[j];
		    }
		    d += d0 * *(w++)  / (*(coef1++) * *(coef2++) * *(coef3++));
		}
    return d * phgGeomGetVolume(u->g, e);
}
/*This quad func for (u, A1, A2, A3, v) suppose that A1!= NULL*/
FLOAT
phgQuadBasAAABas(SIMPLEX *e, DOF *u, int n, DOF *A1, DOF *A2, DOF *A3, DOF *v, int m, int order)
{
    int i, j, k, nvalues = DofTypeDim(u);
    const FLOAT *g1, *g2, *coef1, *coef2, *coef3,*w;
    FLOAT d, d0;
    QUAD *quad;

    assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));

    if (nvalues != DofTypeDim(v))
	phgError(1, "%s:%d: dimensions mismatch: %s (%d) <==> %s ($d))\n",
		 __FILE__, __LINE__, DofTypeName(u), DofTypeDim(u),
		 DofTypeName(v), DofTypeDim(v));
	order = 4;
    quad = phgQuadGetQuad3D(order);

    g1 = phgQuadGetBasisValues(e, u, n, quad);
    g2 = phgQuadGetBasisValues(e, v, m, quad);
    d = 0.;
    w = quad->weights;

    if (A3 != NULL) {
    coef1 = phgQuadGetDofValues(e, A1, quad);
    coef2 = phgQuadGetDofValues(e, A2, quad);
    coef3 = phgQuadGetDofValues(e, A3, quad);
		for (i = 0; i < quad->npoints; i++) {
		    d0 = 0.;
		    for (j = 0; j < nvalues; j++) {
			d0 += *(g1++) * (*(g2++));
		    }
		    d += d0 * (*(w++)) * *(coef1++) * *(coef2++) * *(coef3++);
		}
    return d * phgGeomGetVolume(u->g, e);
	}
	else if (A2 != NULL){
    coef1 = phgQuadGetDofValues(e, A1, quad);
    coef2 = phgQuadGetDofValues(e, A2, quad);
		for (i = 0; i < quad->npoints; i++) {
		    d0 = 0.;
		    for (j = 0; j < nvalues; j++) {
			d0 += *(g1++) * *(g2++);
		    }
		    d += d0 * *(w++) * *(coef1++) * *(coef2++);
		}
    return d * phgGeomGetVolume(u->g, e);
	}
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
/*this codes for (div (\alpha / \beta u), p)*/
FLOAT
phgQuadDivBasADBas(SIMPLEX *e, DOF *u, int m, DOF *A, DOF *D, DOF *v, int n, int order)
{
    int i, j, nvalues = DofTypeDim(v);
    const FLOAT *g1, *g2, *w, *coef1, *coef2, *lambda;
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
    coef1 = phgQuadGetDofValues(e, A, quad);
    coef2 = phgQuadGetDofValues(e, D, quad);
    d = 0.;
    lambda = quad->points;
    w = quad->weights;
    for (i = 0; i < quad->npoints; i++) {
	d0 = 0.;
	for (j = 0; j < nvalues; j++) {
	    d0 += (*(g1)+(*(g1+4))+(*(g1+8))) * (*(g2++)) * *(coef1++) * *(coef2++);
	    g1=g1+9;
	}
	d += d0 * (*(w++));
	lambda += Dim + 1;
    }
    return d * phgGeomGetVolume(u->g, e);
}
FLOAT *
phgQuadDofAAABas(SIMPLEX *e, DOF *u, DOF *A1, DOF *A2, DOF *A3, DOF *v, int n, int order, FLOAT *res)
/* computes \int 'u' * ('n'-th basis function of 'v') on element 'e'
 * using quadrature rule 'quad'. The results are returned in the user provided
 * array 'res' whose size should be greater than or equal to u->dim.
 */
{
    int i, j, nvalues;
    const FLOAT *bas, *w, *f, *coef1, *coef2, *coef3;
    FLOAT d;
    QUAD *quad;

    assert(v != NULL && !SpecialDofType(v->type) && DofTypeDim(v) == 1);

	order = 4;
    quad = phgQuadGetQuad3D(order);

    nvalues = DofDim(u);
    f = phgQuadGetDofValues(e, u, quad);
    bas = phgQuadGetBasisValues(e, v, n, quad);
    w = quad->weights;

	if (A3 != NULL){
   		coef1 = phgQuadGetDofValues(e, A1, quad);
   		coef2 = phgQuadGetDofValues(e, A2, quad);
   		coef3 = phgQuadGetDofValues(e, A3, quad);
		d = 0.0;
		for (i = 0; i < quad->npoints; i++) {
		    d += *(bas++) * *(w++) * *(f++) * *(coef1++) * *(coef2++) * *(coef3++);
		}
		res[0] = d * phgGeomGetVolume(u->g, e);
	    return res;
	}
	else if (A2 != NULL){
   		coef1 = phgQuadGetDofValues(e, A1, quad);
   		coef2 = phgQuadGetDofValues(e, A2, quad);
		d = 0.0;
		for (i = 0; i < quad->npoints; i++) {
		    d += *(bas++) * *(w++) * *(f++) * *(coef1++) * *(coef2++);
		}
		res[0] = d * phgGeomGetVolume(u->g, e);
	    return res;
	}
	else if (A1 != NULL){
   		coef1 = phgQuadGetDofValues(e, A1, quad);
		d = 0.0;
		for (i = 0; i < quad->npoints; i++) {
		    d += *(bas++) * *(w++) * *(f++) * *(coef1++);
		}
		res[0] = d * phgGeomGetVolume(u->g, e);
 	   return res;
	}
}
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
phgQuadBasPermBas(SIMPLEX *e, DOF *u, int n, DOF *perm, DOF *v, int m, int order)
{
    int i, j, k, nvalues = DofTypeDim(u);
    const FLOAT *g1, *g2, *coef, *w;
    FLOAT d, d0;
    QUAD *quad;

    assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));

    if (nvalues != DofTypeDim(v))
	phgError(1, "%s:%d: dimensions mismatch: %s (%d) <==> %s ($d))\n",
		 __FILE__, __LINE__, DofTypeName(u), DofTypeDim(u),
		 DofTypeName(v), DofTypeDim(v));
	if (order < 0){
		order = DofTypeOrder(u, e);
	}
    quad = phgQuadGetQuad3D(order);

    g1 = phgQuadGetBasisValues(e, u, n, quad);
    g2 = phgQuadGetBasisValues(e, v, m, quad);
    d = 0.;
    w = quad->weights;
    coef = phgQuadGetDofValues(e, perm, quad);
	for (i = 0; i < quad->npoints; i++) {
	    d0 = 0.;
	    for (j = 0; j < nvalues; j++) {
		d0 += *(g1++) * *(g2++) / *(coef++);
	    }
	    d += d0 * *(w++);
	}
    return d * phgGeomGetVolume(u->g, e);
}
