#include <stdio.h>
#include <math.h>
#include "PseN.h"
#include "dopri8.h"
#include "GeoCross.h"
/*----------------------------------------------------------------------*/
/* #This function is used as port to the integration function dopri8().	*/
/*									*/
/* #The input varialbes are:						*/
/* norpar nordata: transport in al the normalized physical data.	*/
/* orbpar orbdata: transport in all the orbital information.		*/
/* double tbin: the begin time.						*/
/* double tend: the end time.						*/
/* double* ybin: the beginning point of general coordinates y.		*/
/* int steps: the total steps for integration.				*/
/* int haccu: the least substeps for each step.				*/
/* int nvar: the dimensions of y, often 4.				*/
/* fluxfn fn: pointer to the function of force.				*/
/* fluxEn En: pointer to the function of energy.			*/
/*									*/
/* #This function output the time series ,orbital information and energy*/
/* to xp, yp and Ene. The initial value, namely, tbin and ybin are NOT	*/
/* saved.								*/
/*----------------------------------------------------------------------*/
void geoIntegrator(
norpar *nordata,
orbpar *orbdata,
double tbin,
double tend,
double* ybin,
int steps,
int haccu,
int nvar,
fluxfn fn,
fluxEn En){
	int cstep;
	int *active;
	double *xp, *Ene;
	double **yp;
	int i, j, k;
	double *y;
	double Ene1, x1, x2, xbin, xend, delX;
	double eps, h0, hmax;
	inData data;
	void *custom_data = &data;
	int flag;

// Initialize-----------------------------------------------------------
	cstep = orbdata->cstep + 1;
	active = &(orbdata->active[cstep]);
	xp = &(orbdata->xp[cstep]);
	yp = &(orbdata->yp[cstep]);
	Ene = &(orbdata->Ene[cstep]);
	x1 = tbin;
	x2 = tend;
	delX = (x2 - x1) / steps;
	eps = 1e-13;
	hmax = delX / haccu;
	h0 = 1e-10;
	data.mass = nordata->bhmass;
	data.lightspeed = nordata->lightspeed;
	data.Newton_G = nordata->newton_g;

	y = (double *)malloc(sizeof(double) * nvar);
	for(i = 0; i < nvar; i++){
		y[i] = ybin[i];
	}
	Ene1 = En(nvar, y, custom_data);
//	printf("geoIntegrator: Initial quantities:\n");

// Run--------------------------------------------------------------
	for(i = 0; i < steps; i++){
		xbin = x1 + delX * i;
		xend = x1 + delX * (i + 1);
		flag = dopri8(fn, nvar, xbin, y, xend, eps, hmax, &h0, NULL, NULL, custom_data);
		if(flag == 0) printf("!!%d\n", i);
		Ene[i] = En(nvar, y, custom_data);
		for(j = 0; j < nvar; j++){
			yp[i][j] = y[j];
		}
		active[i] = active[i - 1];
		xp[i] = xend;
	}
	orbdata->cstep += steps;
//	printf("geoIntegrator: Integration finished!\n");
	free(y);
}
