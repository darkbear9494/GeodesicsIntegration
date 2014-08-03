#include <stdio.h>
#include <math.h>
#include "PseN.h"
#include "dopri8.h"
#include "GeoCross.h"

void dyOvdx(double x, double *y, double *dydx, void *p);

int main(){
	double mstar, rstar, bhmass, rtidal, period, rstar_norm;
	double Mstar, Rstar, BHmass, Rtidal, Period, Tratio, Vcsun;
	double vnorm, Vorig;
	double Ezeta, AngMom, apo, per, Va, Vp, mjAxis, ecc;
	double rIni, vIni, T;
	double Newton_G, lightspeed, R_sch;

	mstar = 1.0;
	rstar = 1.0;
	bhmass = 1.0e6;

// Physics----------------------------------------------------------
	Mstar = mstar * MSUN;
	Rstar = rstar * RSUN;
	BHmass = bhmass * MSUN;
	Rtidal = Rstar * pow((BHmass / Mstar), 1.0 / 3);
	Vorig = pow((NEWTON_G * BHmass / Rtidal), 0.5);
	Vcsun = pow((NEWTON_G * Mstar / Rtidal), 0.5);
	Period = 2.0 * PI * pow(Rtidal, 1.5) * pow((NEWTON_G * BHmass), -0.5);
	Tratio = pow(Rtidal, 1.5) * pow((NEWTON_G * Mstar), -0.5);
	printf("Physical quantities:\n");
	printf("Mstar = %e, Rstar = %e, BHmass = %e\n", Mstar, Rstar, BHmass);
	printf("Rtidal = %e, Vorig = %e, Period = %e\n", Rtidal, Vorig, Period);
	printf("Tratio = %e, Vcsun = %e\n", Tratio, Vcsun);
// Set M_sun, r_tidal and Newton_G as unity-----------------------------
	Newton_G = 1.0;
	lightspeed = LIGHTSPEED / Vcsun;
	R_sch = (Newton_G * bhmass / dSqr(lightspeed));
	rtidal = 1.0;
	rstar_norm = Rstar / Rtidal;
	vnorm = pow((Newton_G * bhmass / rtidal), 0.5);
	period = 2.0 * PI * pow(rtidal, 1.5) * pow((Newton_G * bhmass), -0.5);
	printf("Normalized quantities:\n");
	printf("mstar = %e, rstar_norm = %e, bhmass = %e\n", mstar, rstar_norm, bhmass);
	printf("rtidal = %e, vnorm = %e, period = %e\n", rtidal, vnorm, period);
	printf("Newton_G = %e, lightspeed = %e\n", Newton_G, lightspeed);
	printf("R_sch = %e\n", R_sch);

// Tidal disruption parameters:
	Ezeta = -1.0 * Newton_G * bhmass * rstar_norm * pow(rtidal, -2);
	AngMom = pow((2.0 * Newton_G * bhmass), 0.5);
	mjAxis = -0.5 * Newton_G * bhmass / Ezeta;
	ecc = pow((1.0 - dSqr(AngMom) / (Newton_G * bhmass * mjAxis)), 0.5);
	apo = mjAxis * (1.0 + ecc);
	per = mjAxis * (1.0 - ecc);
	Va = pow(2 * (Ezeta + (Newton_G * bhmass / apo)), 0.5);
	Vp = pow(2 * (Ezeta + (Newton_G * bhmass / per)), 0.5);
	T = 2.0 * PI * pow(mjAxis, 1.5) * pow((Newton_G * bhmass), -0.5);
	printf("TD parameters:\n");
	printf("E = %e, L = %e, a = %e, T = %e\n", Ezeta, AngMom, mjAxis, T);
	printf("e = %e, apo = %e, per = %e, Va = %e, Vp = %e\n", ecc, apo, per, Va, Vp);

// Initialize-----------------------------------------------------------
	rIni = apo;
	vIni = Va;
	printf("Initial quantities:\n");
	printf("rIni = %e, vIni = %e, T = %e\n", rIni, vIni, T);
	printf("Etotal = %e\n",(0.5 * vIni*vIni - Newton_G * bhmass / rIni));

// Integration----------------------------------------------------------
	double xp[STEPS], yp[STEPS][NVAR], Ene[STEPS];
	int i, j, k;
	double y1[NVAR], y[NVAR], Ene1, x1, x2, xbin, xend, delX;
	double eps, h0, hmax;
	inData data;
	void *custom_data = &data;
	int flag, nvar;
	fluxfn fn = F_SG;
	fluxEn En = E_SG;
	FILE *fp;

	x1 = 0;
	x2 = 5 * T;
	delX = (x2 - x1) / STEPS;
	eps = 1e-13;
	hmax = delX / HACCU;
	h0 = 1e-10;
	nvar = NVAR;
	data.mass = bhmass;
	data.lightspeed = lightspeed;
	data.Newton_G = Newton_G;

	for(i = 0; i < NVAR; i++){
		y1[i] = 0;
	}
	y1[1] = rIni;
	y1[2] = -1.0 * vIni;
	for(i = 0; i < NVAR; i++){
		y[i] = y1[i];
	}
	Ene1 = En(nvar, y, custom_data);
// Run--------------------------------------------------------------
	for(i = 0; i < STEPS; i++){
		xbin = x1 + delX * i;
		xend = x1 + delX * (i + 1);
		flag = dopri8(fn, nvar, xbin, y, xend, eps, hmax, &h0, NULL, NULL, custom_data);
		if(flag == 0) printf("!!%d\n", i);
		Ene[i] = En(nvar, y, custom_data);
		for(j = 0; j < NVAR; j++) yp[i][j] = y[j];
		xp[i] = xend;
	}	

	fp = fopen("result","w+");
	if(fp == NULL)
		printf("Cannot open the file!\n");
	
	fprintf(fp, "%10.6f \t", x1);
	for(j = 0; j < NVAR; j++){
		fprintf(fp, "%20.15f \t", y1[j]);
	}
	fprintf(fp, "%20.15f \n", Ene1);
	for(i = 0; i < STEPS; i++){
		fprintf(fp, "%10.6f \t", xp[i]);
		for(j = 0; j < NVAR; j++){
			fprintf(fp, "%20.15f \t", yp[i][j]);
		}
		fprintf(fp, "%20.15f \n", Ene[i]);
	}
}

void dyOvdx(double x, double *y, double *dydx, void *p){
	int i;
	for(i = 0; i < NVAR; i++){
		dydx[i] = 2 * x;
	}
}
