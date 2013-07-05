#include "phg.h"
#include "math.h"
int
MATLAB_Draw3D(char *file_name, DOF *dof, INT NX, INT NY, INT NZ)
{
	if(NX <= 0 || NY <= 0 || file_name == NULL)
		return;
	FLOAT dx, dy, dz;
	FLOAT MESH_POINT1[3] = {0,0,0};
	FLOAT MESH_POINT2[3] = {400,400,10};

	dx = (MESH_POINT2[0] - MESH_POINT1[0]) / NX;
	dy = (MESH_POINT2[1] - MESH_POINT1[1]) / NY;
	dz = (MESH_POINT2[2] - MESH_POINT1[2]) / NZ;
	COORD points[(NX + 1) * (NY + 1)];
	int i, j, level;
	FLOAT *value;
	value = phgAlloc((NX + 1) * (NY + 1) * sizeof(*value));

	FLOAT zz = MESH_POINT1[2];
	for(level = 0; level < NZ; level++){
		for(i = 0; i < NX + 1; i++){
			for(j = 0; j < NY + 1; j++){
				points[i * (NY + 1) + j][0] = MESH_POINT1[0] + i * dx;
				points[i * (NY + 1) + j][1] = MESH_POINT1[1] + j * dy;
				points[i * (NY + 1) + j][2] = MESH_POINT1[2] + level * dz;;//zz + 0.5 * MESH_DZ[MESH_NZ - 1 - level];
			}
		}
//		zz += MESH_DZ[MESH_NZ - 1 - level];

		phgInterGridDofEval(dof, (NX + 1) * (NY + 1), points, value, 0);

		char file[60];
		if(dof->g->rank == 0){
			FILE *fp;
			sprintf(file, "%s_level%d.m", file_name, level);
			fp = fopen(file, "w");
			fprintf(fp, "function %s\n", file_name);
			fprintf(fp, "x = [\n");
			for(i = 0; i < NX + 1; i++){
				for(j = 0; j < NY + 1; j++){
					fprintf(fp, "%lf ", points[i * (NY + 1) + j][0]);
				}
				fprintf(fp, "\n");
			}
			fprintf(fp, "];\n");

			fprintf(fp, "y = [\n");
			for(i = 0; i < NX + 1; i++){
				for(j = 0; j < NY + 1; j++){
					fprintf(fp, "%lf ", points[i * (NY + 1) + j][1]);
				}
				fprintf(fp, "\n");
			}
			fprintf(fp, "];\n");

			fprintf(fp, "z = [\n");
			for(i = 0; i < NX + 1; i++){
				for(j = 0; j < NY + 1; j++){
					fprintf(fp, "%lf ", value[i * (NY + 1) + j]);
				}
				fprintf(fp, "\n");
			}
			fprintf(fp, "];\n");
	
			fprintf(fp, "surf(x,y,z);\n");
			fclose(fp);
		}
	}

	phgFree(value);

}
