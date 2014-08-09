#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GeoCross.h"
/*----------------------------------------------------------------------*/
/* #This function is to save the results.				*/
/*									*/
/* #The input varialbes are:						*/
/* particle *parlist: transport in all the particle information.	*/
/* shockhead *sh: the record of particles collisions.			*/
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

int geoSave(particle *parlist, 
shockhead *sh,
int parnum, 
int points, 
int nvar){
	FILE *fp1, *fp2;
	int No, No1, No2;
	int i, j, k;
	int cstep, clength, length;
	int *active;
	double *xp, **yp, *Ene;
	shocklist *sl;
	double t, x, y, *y1, *y2, Erel;
	
	fp1 = fopen("resultPar", "w+");
	if(fp1 == NULL) printf("Cannot open file: fp1\n");
	fp2 = fopen("resultSh", "w+");
	if(fp2 == NULL) printf("Cannot open file: fp2\n");
	
	for(i = 0; i < parnum; i++){
		No = parlist[i].No;
		cstep = parlist[i].orbdata.cstep;
		active = parlist[i].orbdata.active;
		xp = parlist[i].orbdata.xp;
		yp = parlist[i].orbdata.yp;
		Ene = parlist[i].orbdata.Ene;
		length = iMIN(cstep, points);
//		printf("geoSave: cstep = %d\n", cstep);
//		printf("cstep=%d, points=%d, length=%d\n", cstep, points, length);
		for(j = 0; j < length; j++){
			fprintf(fp1, "%5d \t", No);
			fprintf(fp1, "%10d \t", active[j]);
			fprintf(fp1, "%10.6f \t", xp[j]);
			for(k = 0; k < nvar; k++){
				fprintf(fp1, "%20.15f \t", yp[j][k]);
			}
			fprintf(fp1, "%20.15f \t %d\n", Ene[j], j);
		}
	}

	clength = sh->clength;
	if(clength == -1){
		printf("geoSave: No shock has been recorded!\n");
		fprintf(fp2, "%3d\t %3d\t %10.6f\t %20.15f\t %20.15f\t %10.6f\n", -1, -1, -1, -1, -1, -1);
		fclose(fp1);
		fclose(fp2);
		return 0;
	}
	sl = sh->sl;
	length = iMIN((clength + 1), parnum);
//	printf("geoSave: length=%d\n", length);
	for(i = 0; i < length; i++){
		No1 = sl[i].pN1;
		No2 = sl[i].pN2;
		t = sl[i].t;
		x = sl[i].x;
		y = sl[i].y;
		Erel = sl[i].Erel;
		y1 = sl[i].y1;
		y2 = sl[i].y2;
		fprintf(fp2, "%3d \t%3d \t", No1, No2);
		fprintf(fp2, "%10.6f \t", t);
		fprintf(fp2, "%20.15f \t %20.15f \t", x, y);
		fprintf(fp2, "%10.6f \t", Erel);
		for(j = 0; j < nvar; j++){
			fprintf(fp2, "%20.15f \t", y1[j]);
		}
		for(j = 0; j < nvar; j++){
			fprintf(fp2, "%20.15f \t", y2[j]);
		}
		fprintf(fp2, "\n");
	}
	fclose(fp1);
	fclose(fp2);
	return 1;
}

int iMIN(int a, int b){
	return (a < b ? a:b);
}
