#include "oil_water.h"
#include "phg.h"
static FLOAT SWC =0.2, SOR = 0.2;
static INT NO =2, NW = 2;
void
create_kr(PHASE *oil, PHASE *water)
{
	GRID *g = water->Kr->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_kro, *p_krw, S =0;
	ForAllElements(g, e){
		p_so = DofElementData(oil->S, e->index);
		p_kro = DofElementData(oil->Kr, e->index);
		p_krw = DofElementData(water->Kr, e->index);
		S = (1. - p_so[0] - SWC) / (1 - SWC - SOR);
		p_kro[0] =pow((1 - S), NO);
		p_krw[0] = pow(S, NW);
	}
}
void
create_dot_kr(PHASE *oil, PHASE *water)
{
	GRID *g = water->DKr->g;
	SIMPLEX *e;
	FLOAT *p_so, *p_dotw, *p_doto;
	assert(water->DKr->type == oil->DKr->type);
	ForAllElements(g, e){
		p_so = DofElementData(oil->S, e->index);
		p_dotw = DofElementData(water->DKr, e->index);
		p_doto = DofElementData(oil->DKr, e->index);
		p_doto[0] = NO * (1. / (1. - SWC - SOR)) * pow((p_so[0] - SOR) / (1 - SWC - SOR), NO-1);
		p_dotw[0] = NW * (-1. / (1. - SWC - SOR)) * pow((1. - p_so[0] - SWC) / (1 - SWC - SOR), NW-1);
	}
}
void
update_bo(PHASE *oil)
{
	GRID *g = oil->B->g;
	SIMPLEX *e;
	FLOAT *p_bo, *p_p;
	FLOAT press_factor = 0.00689;
	FLOAT p1 = 300 * press_factor, p2 = 800 * press_factor, p3 = 8000 * press_factor;
	FLOAT B1 =1.05, B2 =1.02, B3 = 1.01;
		ForAllElements(g, e){
			p_p = DofElementData(oil->P, e->index);
			p_bo = DofElementData(oil->B, e->index);
			if(p_p[0] <= p1)
					p_bo[0] = 1. / B1;
			if(p_p[0]>p1 && p_p[0] <=p2)
					p_bo[0] = 1./ (B1 + (p_p[0] - p1)/ (p2-p1) * (B2 -B1)); 
			if(p_p[0]>p2 && p_p[0] <=p3)
					p_bo[0] = 1./ (B2 + (p_p[0] - p2)/ (p3-p2) * (B3 -B2));
		}
}
void
update_dot_bo(PHASE *oil)
{
	GRID *g = oil->B->g;
	SIMPLEX *e;
	FLOAT *p_db, *p_p;
	FLOAT press_factor = 0.00689;
	FLOAT p1 = 300 * press_factor, p2 = 800 * press_factor, p3 = 8000 * press_factor;
	FLOAT B1 =1.05, B2 =1.02, B3 = 1.01;
		ForAllElements(g, e){
			p_p = DofElementData(oil->P, e->index);
			p_db = DofElementData(oil->DB, e->index);
			if(p_p[0]>=p1 && p_p[0] <=p2)
				p_db[0] = (1./ B2 - 1./B1) / (p2 -p1); 
			if(p_p[0]>p2 && p_p[0] <=p3)
				p_db[0] = (1./ B3 - 1./B2) / (p3 -p2); 
		}
}
void
update_bw(PHASE *water)
{
	GRID *g = water->B->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_bw, *p_p0;
	FLOAT B0 = 1.01, C_W = 0.435e-3;
		ForAllElements(g, e){
			p_p = DofElementData(water->P, e->index);
			p_p0 = DofElementData(water->P0, e->index);
			p_bw = DofElementData(water->B, e->index);
			p_bw[0] = (1. + C_W * (p_p[0] - p_p0[0])) / B0;
		}
}
void
update_dot_bw(PHASE *water)
{
	GRID *g = water->B->g;
	SIMPLEX *e;
	FLOAT *p_db;
	FLOAT B0 = 1.01, C_W = 0.435e-3;
		ForAllElements(g, e){
			p_db = DofElementData(water->DB, e->index);
			p_db[0] = C_W / B0;
		}
}
void
update_phi(MEDIUM *rock, PHASE *oil)
{
	GRID *g = rock->phi->g;
	SIMPLEX *e;
	FLOAT *p_p, *p_phi, *p_p0, *p_phi0;
	INT nbas = rock->phi->type->nbas;
	int i;
	ForAllElements(g, e){
		p_p0 = DofElementData(oil->P0, e->index);
		p_p = DofElementData(oil->P, e->index);
		p_phi = DofElementData(rock->phi, e->index);
		p_phi0 = DofElementData(rock->phi0, e->index);
		for (i = 0; i < nbas; i++){
			p_phi[i] = p_phi0[0] * (1.0 + rock->C_R * (p_p[i] - p_p0[i]));
		}
	}
}
void
update_dot_phi(MEDIUM *rock)
{
	GRID *g = rock->Dphi->g;
	SIMPLEX *e;
	FLOAT *p_phi0, *p_dot;
	INT nbas = rock->phi->type->nbas;
	int i;
	ForAllElements(g, e){
		p_dot = DofElementData(rock->Dphi, e->index);
		p_phi0 = DofElementData(rock->phi0, e->index);
		p_dot[0] = p_phi0[0] * rock->C_R;
	}
}
void
Init_Medium(const char *file, MEDIUM *rock)
{		
		FILE *fp = fopen(file, "r");
		char *ff = "medium.tmp";
		FILE *fn = fopen(ff, "w");
		INT N = rock->NX * rock->NY * rock->NZ; 
		FLOAT *p_phi, *p_perm;
		GRID *g = rock->phi->g;
		SIMPLEX *e;
		char line[500];
		int i;
		while((fgets(line, 500, fp)) != NULL){
			for(i = 0; line[i] != '\0'; i++){
				if(line[i] != '#')
					fprintf(fn, "%c", line[i]);
				else
					break;
			}
		}
		fclose(fp);
		fclose(fn);
		fp = fopen(ff, "r");
		FLOAT *array = phgAlloc(4*N * sizeof(FLOAT));
		if (array != NULL){
			for (i = 0; i < 4*N; i++){
				fscanf(fp, "%lf", array+i);
			}
		}
		COORD *p0 = NULL, *p1 = NULL, *p2 = NULL, *p3 = NULL;
		int  j =0, k = 0, n=0;
		FLOAT x =0, y =0, z=0;
		FLOAT factor = 9.86e-16;
		if (!feof(fp)){
			ForAllElements(g, e){
				p0 = g->verts + e->verts[0];
				p1 = g->verts + e->verts[1];
				p2 = g->verts + e->verts[2];
				p3 = g->verts + e->verts[3];
				x = ((*p0)[0] + (*p1)[0] + (*p2)[0] + (*p3)[0]) / 4;
				y = ((*p0)[1] + (*p1)[1] + (*p2)[1] + (*p3)[1]) / 4;
				z = ((*p0)[2] + (*p1)[2] + (*p2)[2] + (*p3)[2]) / 4;
				i = (x / rock->DX);
				j = (y / rock->DY);
				k = (z / rock->DZ);
				n = i + j * rock->NX + (rock->NZ -1-k) * rock->NX * rock->NY;
				p_phi = DofElementData(rock->phi0, e->index);
				p_phi[0] = array[n];
				p_perm = DofElementData(rock->perm, e->index);
				p_perm[0] = array[n+N] * factor;
				p_perm[1] = array[n+2*N] * factor;
				p_perm[2] = array[n+3*N] * factor;
			}
		}
		/*Compute The Average Preproty*/
		k  = 0;
		while(k < rock->NZ){
			FLOAT aphi = 0, kx = 0, ky =0, kz =0;
			for (i = 0; i < rock->NX * rock->NY; i++){
				aphi = aphi + array[i + k * rock->NX * rock->NY];
				kx  = kx + array[i + k * rock->NX * rock->NY + N];
				ky  = ky + array[i + k * rock->NX * rock->NY + 2*N];
				kz  = kz + array[i + k * rock->NX * rock->NY + 3*N];
			}
			aphi =  aphi / rock->NX / rock->NY;
			kx =  kx / rock->NX / rock->NY;
			ky =  ky / rock->NX / rock->NY;
			kz =  kz / rock->NX / rock->NY;
			phgPrintf("\nlevel %d:  aphi:  %lf, kx:   %le, ky:  %le,  kz:  %le\n", k+1, aphi, kx, ky, kz);
			k++;
		}
		phgFree(array);
		fclose(fp);
	    remove(ff);
}
