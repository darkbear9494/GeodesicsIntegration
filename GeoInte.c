#include <stdio.h>
//#include <math.h>
#include "GeoInte.h"

#define MAXSTEPS 10000
#define NVAR 1
int main(){
	float xp[MAXSTEPS], yp[MAXSTEPS][NVAR];
	int i, j, k;
	float ystart[NVAR], x1, x2;
	float eps, h1, hmin;
	int nvar, nok[MAXSTEPS], nbad[MAXSTEPS], kmax, kount;

	x1 = 0;
	x2 = 10;
	eps = 1e-13;
	h1 = 1e-3;
	hmin = 1e-10;
	nvar = NVAR;
	kmax = 0;
	kount = 0;
	for(i = 0; i < NVAR; i++){
		ystart[i] = 0;
	}
	odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, dyOvdx, bsstep);
}

void dyOvdx(float x, float *y, float *dydx){
	int i;
	for(i = 0; i < NVAR; i++){
		dydx[i] = 2 * x;
	}
}
#undef MAXSTEPS
#undef NVAR
