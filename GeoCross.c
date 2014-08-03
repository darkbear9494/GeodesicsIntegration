#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PseN.h"
#include "dopri8.h"
#include "GeoCross.h"

#define MPIRUN
// Parameters for integration.
#define MOVES 100	// Total cycles for the whole simulation.
#define STEPS 2	// Total integration steps
#define HACCU 10000	// The at-least accuracy of h
#define NVAR 4		// The dimension of y
#define FORCE F_GN	// The force function used in geoIntegrator.c
#define ENERGY E_GN	// The energy function used in geoIntegrator.c

// Parameters for simulation. (The units are in SI)
#define PARNUM 2	// Particle number
#define PI 3.141592653589793
#define PARSEC (3.08567758e16)	// parsec
#define ASRUNI 149597870700.0	// AU 
#define LIGHTSPEED 299792458.0	// speed of light
#define NEWTON_G (6.67384e-11)	// gravity constant
#define MSUN (1.9891e30)	// solar mass
#define RSUN (6.955e8)		// solar radius

int main(){
// P1: PREPARE--------------------------------------------------------------------
	double mstar, rstar, bhmass, rtidal, period, rstar_norm;
	double Mstar, Rstar, BHmass, Rtidal, Period, Tratio, Vcsun;
	double vnorm, Vorig;
	double Ezeta, AngMom, apo, per, Va, Vp, mjAxis, ecc, beta;
	double rIni, vIni, T;
	double Newton_G, lightspeed, R_sch;
	int i, j, k;
	phypar phydata;
	norpar nordata;
//	orbpar orbdata;
	tidpar tidata;

	mstar = 1.0;
	rstar = 1.0;
	bhmass = 1.0e6;

// Physics
	Mstar = mstar * MSUN;
	Rstar = rstar * RSUN;
	BHmass = bhmass * MSUN;
	Rtidal = Rstar * pow((BHmass / Mstar), 1.0 / 3);
	Vorig = pow((NEWTON_G * BHmass / Rtidal), 0.5);
	Vcsun = pow((NEWTON_G * Mstar / Rtidal), 0.5);
	Period = 2.0 * PI * pow(Rtidal, 1.5) * pow((NEWTON_G * BHmass), -0.5);
	Tratio = pow(Rtidal, 1.5) * pow((NEWTON_G * Mstar), -0.5);
/*	printf("Physical quantities:\n");
	printf("Mstar = %e, Rstar = %e, BHmass = %e\n", Mstar, Rstar, BHmass);
	printf("Rtidal = %e, Vorig = %e, Period = %e\n", Rtidal, Vorig, Period);
	printf("Tratio = %e, Vcsun = %e\n", Tratio, Vcsun);
*/	
	phydata.Mstar = Mstar;
	phydata.Rstar = Rstar;
	phydata.BHmass = BHmass;
	phydata.Rtidal = Rtidal;
	phydata.Tratio = Tratio;
	phydata.Vcsun = Vcsun;
// Set M_sun, r_tidal and Newton_G as unity
	Newton_G = 1.0;
	lightspeed = LIGHTSPEED / Vcsun;
	R_sch = (Newton_G * bhmass / dSqr(lightspeed));
	rtidal = 1.0;
	rstar_norm = Rstar / Rtidal;
	vnorm = pow((Newton_G * bhmass / rtidal), 0.5);
	period = 2.0 * PI * pow(rtidal, 1.5) * pow((Newton_G * bhmass), -0.5);
/*	printf("Normalized quantities:\n");
	printf("mstar = %e, rstar_norm = %e, bhmass = %e\n", mstar, rstar_norm, bhmass);
	printf("rtidal = %e, vnorm = %e, period = %e\n", rtidal, vnorm, period);
	printf("Newton_G = %e, lightspeed = %e\n", Newton_G, lightspeed);
	printf("R_sch = %e\n", R_sch);
*/
	nordata.mstar = mstar;
	nordata.rstar = rstar;
	nordata.bhmass = bhmass;
	nordata.newton_g = Newton_G;
	nordata.lightspeed = lightspeed;
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
/*	printf("TD parameters:\n");
	printf("E = %e, L = %e, a = %e, T = %e\n", Ezeta, AngMom, mjAxis, T);
	printf("e = %e, apo = %e, per = %e, Va = %e, Vp = %e\n", ecc, apo, per, Va, Vp);
*/	
	tidata.Ezeta = Ezeta;
	tidata.AngMom = AngMom;
	tidata.Rtidal = Rtidal;
	tidata.beta = beta;
// P2: INITIALIZE-----------------------------------------------------------------
	printf("GeoCross: INITIALIZE!\n");
	particle *parlist;
	orbpar *orbdata;
	parlist = (particle *)malloc(sizeof(particle) * PARNUM);
	double *tbinList, **ybinList, *EbinList;
	int parnum, points, moves, steps, haccu, nvar;

	parnum = PARNUM;
	moves = MOVES;
	steps = STEPS;
	points = moves * steps + 1; // The 0th point reserves the initial value.
	haccu = HACCU;
	nvar = NVAR;
	
	tbinList = (double *)malloc(sizeof(double) * parnum);
	ybinList = (double **)malloc(sizeof(double *) * parnum);
	EbinList = (double *)malloc(sizeof(double) * parnum);
	for(i = 0; i < parnum; i++){
		ybinList[i] = (double *)malloc(sizeof(double) * nvar);
		tbinList[i] = 0; // This should be some integrated number.
		EbinList[i] = Ezeta;
		ybinList[i][0] = 0;
		ybinList[i][1] = apo;
		ybinList[i][2] = Va;
		ybinList[i][3] = 0;
	}
	
	iniPar(parnum, points, nvar, parlist, tbinList, ybinList, EbinList);
//	printf("Initialize over!\n");
// P3: INTEGRATION----------------------------------------------------------------
	printf("GeoCross: INTEGRATION!\n");
	double tbin, tend;
//	double *ybin;
	fluxfn fn = FORCE;
	fluxEn En = ENERGY;
//	printf("GeoCross: xp = %p, yp = %p, Ene = %p\n", orbdata.xp, orbdata.yp, orbdata.Ene);
/*	for(i = 0; i < parnum; i++){
		orbdata = &(parlist[i].orbdata);
		ybin = ybinList[i];
		geoIntegrator(&nordata, orbdata, tbin, tend, ybin, steps, haccu, nvar, fn, En);
	}
*/
	for(i = 0; i < moves; i++){
		tbin = i * (5 * T) / moves;
		tend = (i + 1) * (5 * T) / moves;
		geoMove(&nordata, parlist, tbin, tend, parnum, steps, haccu, nvar, fn, En);
//		printf("LOOP%d ok!!!!!\n", i);
	}
// P4: SAVING RESULTS------------------------------------------------------------
	printf("GeoCross: SAVING RESULTS!\n");

	FILE *fp;
	fp = fopen("result","w+");
	if(fp == NULL)
		printf("Cannot open the file!\n");
	geoSave(fp, parlist, parnum, points, nvar);
	fclose(fp);
}
