#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>
#include "PseN.h"
#include "dopri8.h"
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
//		printf("orbdata: %p\n", orbdata);
//		printf("ybin: %p\n", ybin);
		geoIntegrator(nordata, orbdata, tbin, tend, ybin, steps, haccu, nvar, fn, En);
//		printf("%d\n", cstep);
	}
}
