#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GeoCross.h"

void iniQuant(
double mstar,
double rstar,
double bhmass,
double beta,
double zetaMin,
phypar *phydata,
norpar * nordata){
	double Mstar, Rstar, BHmass, Rtidal, Rp, Tratio, Vcsun;
	double newton_g, lightspeed, rtidal, rp, tmin, rstar_norm;

	Mstar = mstar * MSUN;
	Rstar = rstar * RSUN;
	BHmass = bhmass * MSUN;
	Rtidal = Rstar * pow((BHmass / Mstar), 1.0 / 3);
	Rp = Rtidal / beta;
	Vcsun = pow((NEWTON_G * Mstar / Rtidal), 0.5);
	Tratio = pow(Rtidal, 1.5) * pow((NEWTON_G * Mstar), -0.5);
	
	phydata->Msun = MSUN;
	phydata->Mstar = Mstar;
	phydata->Rstar = Rstar;
	phydata->BHmass = BHmass;
	phydata->Rtidal = Rtidal;
	phydata->Rp = Rp;
	phydata->Tratio = Tratio;
	phydata->Vcsun = Vcsun;

	newton_g = 1.0;
	lightspeed = LIGHTSPEED / Vcsun;
	rtidal = 1.0;
	rp = Rp / Rtidal;
	rstar_norm = Rstar / Rtidal;
	tmin= 2.0 * PI * pow(2 * rp, 1.5) * pow((newton_g * bhmass), -0.5);
	
	nordata->mstar = mstar;
	nordata->rstar = rstar_norm;
	nordata->bhmass = bhmass;
	nordata->rtidal = rtidal;
	nordata->rp = rp;
	nordata->zetaMin = zetaMin;
	nordata->tmin = tmin;
	nordata->newton_g = newton_g;
	nordata->lightspeed = lightspeed;
}
