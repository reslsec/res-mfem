#include <stdio.h>
#include <string.h>
#include <math.h>
#include "phg.h"

void readfile(void)
{
	FILE *fp = fopen("pvt.dat", "r");
	FLOAT p[1000], mu_o[1000], bo[1000], mu_g[1000], bg[1000], Rs[1000];
	int i = 0;
	while(!feof(fp)){
		FLOAT co = 0, mu_w = 0, bw = 0, cw = 0, cg = 0,
		fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &p[i],&mu_o[i], &bo[i], &co, &mu_w, &bw, &cw, &mu_g[i], &bg[i], &cg, &Rs[i]);
		i++;
	}
}
void
Create_pvt_bo()
void
Create_pvt_muo()
void
Create_pvt_rs()
void
Create_pvt_bg()
void
Create_pvt_mug()
