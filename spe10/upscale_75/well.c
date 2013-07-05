#include "phg.h"
#include "oil_water.h"
#include "math.h"
void 
Mark_All_Well(GRID *g, MEDIUM *rock, WELL *well)
{
	SIMPLEX *e;
	COORD *p0 = NULL, *p1 = NULL, *p2 = NULL, *p3 = NULL;
	int i=0, j=0, k=0, l =0;
	int N =0, flag=0;
	FLOAT x =0, y =0, z=0;
	ForAllElements(g, e){
		p0 = g->verts + e->verts[0];
		p1 = g->verts + e->verts[1];
		p2 = g->verts + e->verts[2];
		p3 = g->verts + e->verts[3];
		x = ((*p0)[0] + (*p1)[0] + (*p2)[0] + (*p3)[0]) / 4;
		y = ((*p0)[1] + (*p1)[1] + (*p2)[1] + (*p3)[1]) / 4;
		i = (x / rock->DX) + 1;
		j = (y / rock->DY) + 1;
		for (flag=0; flag < well->Number; flag++){
			if ((i == well->position[2*flag])& (j == well->position[2*flag+1]))
			{
				e->region_mark = well->region[flag];
				break;
			}
		}
	}
}
