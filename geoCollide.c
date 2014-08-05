#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GeoCross.h"

void colliCheck(particle *, particle *, int);
int lineCross(double *, double *, double *, double *);

void geoCollide(particle *parlist, int parnum, int steps){	
	int i, j;
	particle *particle1, *particle2;
	for(i = 0; i < (parnum - 1); i++){
	 for(j = i + 1; j < parnum; j ++){
		particle1 = parlist + i;
		particle2 = parlist + j;
		colliCheck(particle1, particle2, steps);
	 }
	}
}

void colliCheck(particle *particle1, particle *particle2, int steps){
	double *y1p, *y1, *y2p, *y2;
	int cstep, i, j, k;
	cstep = particle1->orbdata.cstep;
	for(i = 0; i < steps; i++){
		j = cstep - steps + i;
		y1p = particle1->orbdata.yp[j - 1];
		y2p = particle2->orbdata.yp[j - 1];
		y1 = particle1->orbdata.yp[j];
		y2 = particle2->orbdata.yp[j];
		if(lineCross(y1p, y1, y2p, y2)){
			printf("colliCheck: Find collide!\n");
			printf("colliCheck: step = %d\n", j);
		}
	}
}
int lineCross(double *y1p, double *y1, double *y2p, double *y2){
	return 0;
}
