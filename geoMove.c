#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GeoCross.h"
/*----------------------------------------------------------------------*/
/* #This function is used to move all particles into the next position.	*/
/*									*/
/* #The input varialbes are:						*/
/* norpar *nordata: transport in al the normalized physical data.	*/
/* particle *parlist: transport in all the particles' information.	*/
/* double tbin: the begin time.						*/
/* double tend: the end time.						*/
/* int parnum: the total number of particles.				*/
/* int steps: the total steps for integration.				*/
/* int haccu: the least substeps for each step.				*/
/* int nvar: the dimensions of y, often 4.				*/
/* fluxfn fn: pointer to the function of force.				*/
/* fluxEn En: pointer to the function of energy.			*/
/*									*/
/* #This function mainly use function geoIntergrator().			*/
/*----------------------------------------------------------------------*/

void geoMove(
norpar *nordata,
particle *parlist, 
double tbin, 
double tend, 
int parnum,
int steps, 
int haccu, 
int nvar, 
fluxfn fn, 
fluxEn En){
	int i;
	int cstep;
	double *ybin;
	orbpar *orbdata;
//	printf("parlist: %p\n", parlist);
//	printf("%d\n", omp_get_num_procs());
//	printf("geoMove: tbin = %e, tend = %e\n", tbin, tend);
//	#pragma omp parallel for
	for(i = 0; i < parnum; i++){
		orbdata = &(parlist[i].orbdata);
		cstep = orbdata->cstep;
		ybin = orbdata->yp[cstep];
//		printf("geoMove: x=%f, y=%f, vx=%f, vy=%f\n", ybin[0], ybin[1], ybin[2], ybin[3]);
		if(orbdata->active[cstep] != -1)
			//printf("geoMove: cstep=%d\n", cstep);
			geoIntegrator(nordata, orbdata, tbin, tend, ybin, steps, haccu, nvar, fn, En);
//		printf("%d\n", cstep);
	}
}
