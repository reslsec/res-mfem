#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "phg.h"
int main ()
{
	      FILE *fr, *fw;
	      fr = fopen("krg.dat","r");
	      fw = fopen("krg.txt","w");
	      double kr[2 * 209];
	      int i = 0, j = 0;
	      for (i = 0; i < 2 * 209; i++){
		            kr[i]= 0.; 
		 } 
	      for(i = 0; i < 2*209; i++){
		           fscanf(fr, "%lf", &kr[i]);
		          // printf("%le,  ", kr[i]);
		       }   
	      fclose(fr);
	      for (i = 0;i < 2*209; i++){
		           fprintf(fw, "%lf,", kr[i]);
				 if (i %2 ==1)
			    fprintf(fw, "\n");
		  }
	      fclose(fw);
}

