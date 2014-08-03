#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GeoCross.h"
/*----------------------------------------------------------------------*/
/* #This function is to save the results.				*/
/*									*/
/* #The input varialbes are:						*/
/* FILE *fp: points to the file to save data.				*/
/* particle *parlist: transport in all the particle information.	*/
/* int parnum: the total number of particles.				*/
/* int points: the total points create in memory.			*/
/* int nvar: the dimensions of y, often 4.				*/
/*									*/
/* #This function output a file named "result" with 7 columns:		*/
/* No, time, x, y, vx, vy, energy					*/
/*									*/
/* #The function iMIN() is to return the small integor in the input two.*/
/*----------------------------------------------------------------------*/

int iMIN(int, int);

void geoSave(FILE *fp,
particle *parlist, 
int parnum, 
int points, 
int nvar){
	int No;
	int i, j, k;
	int cstep, length;
	double *xp, **yp, *Ene;
	
	for(i = 0; i < parnum; i++){
		No = parlist[i].No;
		cstep = parlist[i].orbdata.cstep;
		xp = parlist[i].orbdata.xp;
		yp = parlist[i].orbdata.yp;
		Ene = parlist[i].orbdata.Ene;
		length = iMIN(cstep, points);
//		printf("cstep=%d, points=%d, length=%d\n", cstep, points, length);
		for(j = 0; j < length; j++){
			fprintf(fp, "%5d \t", No);
			fprintf(fp, "%10.6f \t", xp[j]);
			for(k = 0; k < nvar; k++){
				fprintf(fp, "%20.15f \t", yp[j][k]);
			}
			fprintf(fp, "%20.15f \n", Ene[j]);
		}
	}
}

int iMIN(int a, int b){
	return (a < b ? a:b);
}
