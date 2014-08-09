#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "GeoCross.h"

void iniTidal_toy(
norpar nordata, 
double zetaMin,
double *mbinList, 
double **ybinList, 
double *EbinList, 
int parnum, 
int nvar){
	int i;
	double vx, vy, E;
	if(parnum < 2){
		printf("The are not enough particles!!!\n");
		exit(0);
	}
	printf("iniTidal_toy: funtion start!\n");
	
	srand(time(NULL));
	//srand(i);
	for(i = 0; i < parnum; i++){
		vx = 2.0 * (rand() - 0.5 * RAND_MAX) / RAND_MAX;
		vy = 2.0 * (rand() - 0.5 * RAND_MAX) / RAND_MAX;
		E = 0.5 * (vx *vx + vy * vy);
		mbinList[i] = 1;
		EbinList[i] = E;	
		ybinList[i][0] = 200.0 * (rand() - 0.5 * RAND_MAX) / RAND_MAX;
		ybinList[i][1] = 200.0 * (rand() - 0.5 * RAND_MAX) / RAND_MAX;
		ybinList[i][2] = vx;
		ybinList[i][3] = vy;
	}
/*	
	mbinList[0] = 1;
	EbinList[0] = 0;	
	ybinList[0][0] = 0;
	ybinList[0][1] = -100.0;
	ybinList[0][2] = 0;
	ybinList[0][3] = 1.0;

	mbinList[1] = 1;
	EbinList[1] = 0;	
	ybinList[1][0] = -100.0;
	ybinList[1][1] = 0;
	ybinList[1][2] = 1.0;
	ybinList[1][3] = 0;

	for(i = 2; i < parnum; i++){
		printf("Extra particles: %d\n", i+1);
		mbinList[i] = 0;
		EbinList[i] = 0;	
		ybinList[i][0] = 0;
		ybinList[i][1] = 0;
		ybinList[i][2] = 0;
		ybinList[i][3] = 0;
	}
*/
}

void F_toy(double t, double *y, double *y1, void *p)
{
	y1[0] = y[2];
	y1[1] = y[3];
	y1[2] = 0;
	y1[3] = 0;
//	printf("-----%f  %f  %f  %f  \n", y[0], y[1], y[2], y[3]);
//	printf("%f  %f  %f  %f  \n", y1[0], y1[1], y1[2], y1[3]);
}

double E_toy(int ndim, double *y, void *p)
{
	double vx, vy;
	vx = y[2];
	vy = y[3];
	return (0.5 * (vx * vx + vy * vy));
}
