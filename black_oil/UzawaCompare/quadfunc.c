#include "phg.h"
#include <string.h>
#include <math.h>

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
FLOAT
phgQuadFaceDofDotBas_(SIMPLEX *e, int face, DOF *u, DOF_PROJ proj,
		     DOF *v,  DOF_PROJ vproj, int N, int order)
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
    
    nvalues = DofTypeDim(v)/ (((vproj == DOF_PROJ_DOT)) ? Dim : 1);

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
		switch(vproj){
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
                     d0 +=  (bas[0]*n[0]+bas[1]*n[1]+bas[2]*n[2])* *(dof++);
                     bas+=Dim;
		 }
		 d += d0 * *(w++);
             break;
		case DOF_PROJ_CROSS:
		   phgError(1," To DO ......\n");
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
