#include "oil_water.h"
#include "phg.h"
#include "math.h"
void
read_data(const char *file, PHASE *oil, PHASE *water, WELL *well, COM_TRO *control, MEDIUM *rock)
{
	char data_name[50];
	FILE *fd, *fw;
	int i, j;

	fd = fopen(file, "r");
	char *f2 = "data.dat";
	fw = fopen(f2, "w");

	char line[500];
	char line2[500];

	/* read a line of the file */
	while((fgets(line, 500, fd)) != NULL){
		for(i = 0; line[i] != '\0'; i++){
			/* line[i] not # ==> print in file */
			if(line[i] != '#')
				fprintf(fw, "%c", line[i]);
			/* line[i] is # ==> delete */
			else
				break;
		}
	}

	fclose(fd);
	fclose(fw);

	fd = fopen("data.dat", "r");
	while(fscanf(fd, "%s", line) == 1){
		if((strcmp(line, "init_pressure")) == 0){
			fscanf(fd, "%lf", &(oil->PRESSURE0));
			fscanf(fd, "%lf", &(water->PRESSURE0));
		}
		/* saturation */
		if((strcmp(line, "saturation")) == 0){
			fscanf(fd, "%lf", &(oil->S0));
			fscanf(fd, "%lf", &(water->S0));
		}
		/*viscosity*/
		if((strcmp(line, "viscosity")) == 0){
			fscanf(fd, "%lf", &(oil->Mu));
			fscanf(fd, "%lf", &(water->Mu));
		}
		/*rock_compressibility*/
		if((strcmp(line, "rock_com")) == 0){
			fscanf(fd, "%lf", &(rock->C_R));
		}
		/* mas_ds */
		if((strcmp(line, "max_ds")) == 0){
			fscanf(fd, "%lf", &(control->MAX_DS));
		}
		if((strcmp(line, "init_dt")) == 0){
			fscanf(fd, "%lf", &(control->init_dt));
		}
		if((strcmp(line, "max_dt")) == 0){
			fscanf(fd, "%lf", &(control->max_dt));
		}
		if((strcmp(line, "time_end")) == 0){
			fscanf(fd, "%lf", &(control->T));
		}
		if((strcmp(line, "err_con")) == 0){
			fscanf(fd, "%le", &(control->TOL_con));
		}
		if((strcmp(line, "err_non")) == 0){
			fscanf(fd, "%le", &(control->TOL_non));
		}
		if((strcmp(line, "err_sys")) == 0){
			fscanf(fd, "%le", &(control->TOL_sys));
		}
		if((strcmp(line, "well_number")) == 0){
			fscanf(fd, "%d", &(well->Number));
		}
		if((strcmp(line, "skin")) == 0){
			fscanf(fd, "%lf", &(well->skin));
		}
		if((strcmp(line, "injection")) == 0){
			fscanf(fd, "%lf", &(well->INJE));
		}
		if((strcmp(line, "max_inje_pwf")) == 0){
			fscanf(fd, "%lf", &(well->MAX_INJE_BHP));
		}
		if((strcmp(line, "prod_pwf")) == 0){
			fscanf(fd, "%lf", &(well->PROD_BHP));
		}
		if((strcmp(line, "mesh_information")) == 0){
			fscanf(fd, "%d", &(rock->NX));
			fscanf(fd, "%d", &(rock->NY));
			fscanf(fd, "%d", &(rock->NZ));
			fscanf(fd, "%lf", &(rock->DX));
			fscanf(fd, "%lf", &(rock->DY));
			fscanf(fd, "%lf", &(rock->DZ));
		}
	}
	fclose(fd);
	well->radius = phgAlloc(well->Number * sizeof(FLOAT));
	bzero(well->radius, well->Number * sizeof(FLOAT));
	well->position = phgAlloc(well->Number * 2 * sizeof(FLOAT));
	bzero(well->position, well->Number * 2 * sizeof(FLOAT));
	well->region = phgAlloc(well->Number * sizeof(FLOAT));
	bzero(well->region, well->Number * sizeof(FLOAT));
	fd = fopen("data.dat", "r");
	while(fscanf(fd, "%s", line) == 1){
		if((strcmp(line, "well_radius")) == 0){
			for (i = 0; i < well->Number; i++){
				fscanf(fd, "%lf", well->radius+i);
			}
		}
		if((strcmp(line, "position")) == 0){
			for (i = 0; i < well->Number * 2; i++){
				fscanf(fd, "%lf", well->position+i);
			}
		}
		if((strcmp(line, "region_mark")) == 0){
			for (i = 0; i < well->Number; i++){
				fscanf(fd, "%lf", well->region+i);
			}
		}
	}
	fclose(fd);
	remove(f2);
}
void
units_convert(PHASE *oil, PHASE *water, WELL *well, MEDIUM *rock)
{
	/* units convert */
		/*
	FLOAT PRESS_FACTOR = 0.006894757; // psi => MPa
	FLOAT LEGTH_FACTOR = 0.3048;//ft->m
	FLOAT VOLUME_FACTOR = 0.159; //bbl->m^3
	FLOAT VIS_FACTOR = 1e-9 /(3600*24);//cp->MPa.day
*/
	oil->PRESSURE0 *= PRESS_FACTOR;
	water->PRESSURE0 *= PRESS_FACTOR;
	well->INJE *= VOLUME_FACTOR;
	well->MAX_INJE_BHP *= PRESS_FACTOR;
	well->PROD_BHP *= PRESS_FACTOR;
	
	oil->Mu *= VIS_FACTOR;
	water->Mu *= VIS_FACTOR;
	
	rock->C_R = rock->C_R / PRESS_FACTOR;

	rock->DX *= LEGTH_FACTOR;
	rock->DY *= LEGTH_FACTOR;
	rock->DZ *= LEGTH_FACTOR;
}

