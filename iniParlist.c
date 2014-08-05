#include <stdio.h>
#include <stdlib.h>
#include "GeoCross.h"
/*----------------------------------------------------------------------*/
/* #This function is to initialize the particles.			*/
/*									*/
/* #The input varialbes are:						*/
/* int parnum: The total number of particles.				*/
/* int points: The total length required from the memory.		*/
/* int nvar: The dimension of orbital parameters.			*/
/* particle *parlist: Points to the list of particles.			*/
/* double *tbinList: the beginning time list.				*/
/* double **ybinList: the beginning y list.				*/
/* double *EbinList: the beginning energy list.				*/
/*									*/
/* #This function output the parlist with initial values on the 0th	*/
/* element.								*/
/*									*/
/*----------------------------------------------------------------------*/

void iniParlist(
int parnum,
int points,
int nvar,
particle *parlist,
double *tbinList,
double **ybinList,
double *EbinList){
	int i, j;
	double *xp, *Ene;
	double **yp;

	for(i = 0; i < parnum; i++){
		parlist[i].No = i;
		parlist[i].active = 1;
		parlist[i].tini = tbinList[i];
		parlist[i].Eini = EbinList[i];

		xp = (double *)malloc(sizeof(double) * points);
		xp[0] = 0;
		Ene = (double *)malloc(sizeof(double) * points);
		Ene[0] = EbinList[i];
		yp = (double **)malloc(sizeof(double *) * points);
		for(j = 0; j < points; j++)
			yp[j] = (double *)malloc(sizeof(double) * nvar);
		for(j = 0; j < nvar; j++){
			yp[0][j] = ybinList[i][j];
		}

		parlist[i].orbdata.cstep = 0;
		parlist[i].orbdata.xp = xp;
		parlist[i].orbdata.yp = yp;
		parlist[i].orbdata.Ene = Ene;
	}
}
