#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PseN.h"
#include "dopri8.h"
#include "GeoCross.h"

int main(){
// P1: PREPARE--------------------------------------------------------------------
	double mstar, rstar, bhmass, rtidal, rp, period, rstar_norm;
	double Mstar, Rstar, BHmass, Rtidal, Rp, Period, Tratio, Vcsun;
	double vnorm, Vorig;
	double zetaMin, Ezeta, AngMom, apo, per, Va, Vp, mjAxis, ecc, beta;
	double rIni, vIni, T;
	double Newton_G, lightspeed, R_sch;
	int i, j, k;
	phypar phydata;
	norpar nordata;
	tidpar tidata;

	mstar = 1.0;
	rstar = 1.0;
	bhmass = 1.0e6;
	beta = 1.0;
	zetaMin = 0.3;

// Physics
	Mstar = mstar * MSUN;
	Rstar = rstar * RSUN;
	BHmass = bhmass * MSUN;
	Rtidal = Rstar * pow((BHmass / Mstar), 1.0 / 3);
	Rp = Rtidal / beta;
	Vorig = pow((NEWTON_G * BHmass / Rtidal), 0.5);
	Vcsun = pow((NEWTON_G * Mstar / Rtidal), 0.5);
	Period = 2.0 * PI * pow(2 * Rp, 1.5) * pow((NEWTON_G * BHmass), -0.5);
	Tratio = pow(Rtidal, 1.5) * pow((NEWTON_G * Mstar), -0.5);
//	printf("Physical quantities:\n");
//	printf("Mstar = %e, Rstar = %e, BHmass = %e\n", Mstar, Rstar, BHmass);
//	printf("Rtidal = %e, Vorig = %e, Period = %e\n", Rtidal, Vorig, Period);
//	printf("Tratio = %e, Vcsun = %e\n", Tratio, Vcsun);
	
	phydata.Msun = MSUN;
	phydata.Mstar = Mstar;
	phydata.Rstar = Rstar;
	phydata.BHmass = BHmass;
	phydata.Rtidal = Rtidal;
	phydata.Rp = Rp;
	phydata.Tratio = Tratio;
	phydata.Vcsun = Vcsun;
// Set M_sun, r_tidal and Newton_G as unity
	Newton_G = 1.0;
	lightspeed = LIGHTSPEED / Vcsun;
	R_sch = (Newton_G * bhmass / (lightspeed * lightspeed));
	rtidal = 1.0;
	rp = Rp / Rtidal;
	rstar_norm = Rstar / Rtidal;
	vnorm = pow((Newton_G * bhmass / rtidal), 0.5);
	period = 2.0 * PI * pow(2 * rp, 1.5) * pow((Newton_G * bhmass), -0.5);
//	T = 100 * period;
//	printf("Normalized quantities:\n");
//	printf("mstar = %e, rstar_norm = %e, bhmass = %e\n", mstar, rstar_norm, bhmass);
//	printf("rtidal = %e, vnorm = %e, period = %e\n", rtidal, vnorm, period);
//	printf("Newton_G = %e, lightspeed = %e\n", Newton_G, lightspeed);
//	printf("R_sch = %e\n", R_sch);
	
	nordata.mstar = mstar;
	nordata.rstar = rstar_norm;
	nordata.bhmass = bhmass;
	nordata.rtidal = rtidal;
	nordata.rp = rp;
	nordata.newton_g = Newton_G;
	nordata.lightspeed = lightspeed;

//	iniQuant(mstar, rstar, bhmass, beta, zetaMin, &phydata, &nordata);
// Tidal disruption parameters:
	Ezeta = -1.0 * Newton_G * bhmass * rstar_norm * pow(rtidal, -2);
//	Ezeta = -1.0 * Newton_G * bhmass * rstar_norm * pow(rtidal, -2) * zetaMin;
	AngMom = pow((2.0 * Newton_G * bhmass), 0.5);
	mjAxis = -0.5 * Newton_G * bhmass / Ezeta;
	ecc = pow((1.0 - (AngMom * AngMom) / (Newton_G * bhmass * mjAxis)), 0.5);
	apo = mjAxis * (1.0 + ecc);
	per = mjAxis * (1.0 - ecc);
	Va = pow(2 * (Ezeta + (Newton_G * bhmass / apo)), 0.5);
	Vp = pow(2 * (Ezeta + (Newton_G * bhmass / per)), 0.5);
	T = 2.0 * PI * pow(mjAxis, 1.5) * pow((Newton_G * bhmass), -0.5);
	printf("TD parameters:\n");
	printf("E = %e, L = %e, a = %e, T = %e\n", Ezeta, AngMom, mjAxis, T);
//	printf("e = %e, apo = %e, per = %e, Va = %e, Vp = %e\n", ecc, apo, per, Va, Vp);
//	printf("GeoCross: %e, %e\n", (AngMom * AngMom), (Newton_G * bhmass * mjAxis));

	tidata.Ezeta = Ezeta;
	tidata.AngMom = AngMom;
	tidata.Rtidal = Rtidal;
	tidata.beta = beta;

// P2: INITIALIZE-----------------------------------------------------------------
	printf("GeoCross: INITIALIZE!\n");
	particle *parlist;
	orbpar *orbdata;
	shockhead sh;
	double *mbinList, **ybinList, *EbinList;
	int parnum, points, cycle, moves, steps, haccu, nvar;

	parnum = PARNUM;
	cycle = CYCLE;
	moves = MOVES;
	steps = STEPS;
	points = moves * steps + 1; // The 0th point reserves the initial value.
	haccu = HACCU;
	nvar = NVAR;
	
	parlist = (particle *)malloc(sizeof(particle) * parnum);
	sh.clength = -1;
	sh.sl = (shocklist *)malloc(sizeof(shocklist) * parnum);
	mbinList = (double *)malloc(sizeof(double) * parnum);
	ybinList = (double **)malloc(sizeof(double *) * parnum);
	EbinList = (double *)malloc(sizeof(double) * parnum);
	for(i = 0; i < parnum; i++)
		ybinList[i] = (double *)malloc(sizeof(double) * nvar);

	iniTidal(nordata, zetaMin, mbinList, ybinList, EbinList, parnum, nvar);
//	iniTidal_toy(nordata, zetaMin, mbinList, ybinList, EbinList, parnum, nvar);
	iniParlist(parnum, points, nvar, parlist, sh.sl, mbinList, ybinList, EbinList);
//	printf("Initialize over!\n");
// P3: INTEGRATION----------------------------------------------------------------
	printf("GeoCross: INTEGRATION!\n");
	double tbin, tend;
	double time = 10.0 * T;
//	double time = 10.0 * nordata.tmin;
//	double time = 50;
	fluxfn fn = FORCE;
	fluxEn En = ENERGY;
//	printf("GeoCross: xp = %p, yp = %p, Ene = %p\n", orbdata.xp, orbdata.yp, orbdata.Ene);
	for(i = 0; i < moves; i++){
		tbin = (cycle * time) * i / moves;
		tend = (cycle * time) * (i + 1) / moves;
		geoMove(&nordata, parlist, tbin, tend, parnum, steps, haccu, nvar, fn, En);
//		printf("LOOP%d ok!!!!!\n", i);
		geoCollide(&nordata, parlist, &sh, parnum, steps);
	}
// P4: SAVING RESULTS------------------------------------------------------------
	printf("GeoCross: SAVING RESULTS!\n");

	geoSave(parlist, &sh, parnum, points, nvar);

	free(parlist);
	free(mbinList);
	free(ybinList);
	free(EbinList);
}
