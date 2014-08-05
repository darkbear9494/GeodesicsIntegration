#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "GeoCross.h"

void orbiPar(norpar , double , double , double *, double *, double *);
void cubiSolver(double *coes, double *roots);
int cmpfunc (const void * a, const void * b);

void iniTidal(
norpar nordata, 
double zetaMin,
double *tbinList, 
double **ybinList, 
double *EbinList, 
int parnum, 
int nvar){
	int i, j;
	double mstar, rstar, bhmass, rp, newton_g, lightspeed, r_g;
	double zeta;
	double Ezeta, AngMom, mjAxis, ecc, apo, per, Va, Vp, tini;

	mstar = nordata.mstar;
	rstar = nordata.rstar;
	bhmass = nordata.bhmass;
	rp= nordata.rp;
	newton_g = nordata.newton_g;
	lightspeed = nordata.lightspeed;
	r_g = newton_g * bhmass / (lightspeed * lightspeed);

	for(i = 0; i < parnum; i++){
// Tidal disruption parameters:
		zeta = 1.0 - (1.0 - zetaMin) * i / parnum;
		Ezeta = -1.0 * newton_g * bhmass * rstar * pow(rp, -2) * zeta;
		AngMom = pow((2.0 * newton_g * bhmass), 0.5);
		orbiPar(nordata, Ezeta, AngMom, &apo, &per, &tini);
		Va = -1.0 * AngMom * (apo - 2.0 * r_g) / (apo * apo * pow((1.0 - 2.0 * Ezeta / (lightspeed * lightspeed)), 0.5));
		Vp = 1.0 * AngMom * (per - 2.0 * r_g) / (per * per * pow((1.0 - 2.0 * Ezeta / (lightspeed * lightspeed)), 0.5));
//		printf("iniTidal: %e, %e, %e\n", apo, per, tini);
// Set initial data:
		tbinList[i] = tini;	// This should be some integrated number.
		EbinList[i] = Ezeta;	// This initial energy value is not the exact energy in general relativistic sense.
		ybinList[i][0] = 0;
		ybinList[i][1] = -1.0 * per;
		ybinList[i][2] = Vp;
		ybinList[i][3] = 0;
	}
}

void orbiPar(
norpar nordata, 
double E, 
double L, 
double *apo, 
double *per, 
double *tini){
	double mstar, rstar, bhmass, newton_g, lightspeed, r_s;
	double ecc, mjAxis;
	double coes[4], roots[3];

	mstar = nordata.mstar;
	rstar = nordata.rstar;
	bhmass = nordata.bhmass;
	newton_g = nordata.newton_g;
	lightspeed = nordata.lightspeed;
	r_s = 2.0 * newton_g * bhmass / (lightspeed * lightspeed);

	coes[0] = 1.0;
	coes[1] = -1.0 * pow((-0.5 * L * L / (E * r_s * r_s)), 1.0/3);
	coes[2] = pow((-0.5 * r_s / (L * E)), 2.0/3) * (lightspeed * lightspeed); 
	coes[3] = -1.0;
//	printf("coes: %f, %f, %f, %f\n", coes[0], coes[1], coes[2], coes[3]);

	cubiSolver(coes, roots);
//	printf("roots: %f, %f, %f\n", roots[0], roots[1], roots[2]);
	
	*apo = pow((-0.5 * L * L / (E * r_s * r_s)), 1.0/3) * r_s / roots[0];
	*per = pow((-0.5 * L * L / (E * r_s * r_s)), 1.0/3) * r_s / roots[1];
	*tini = 0;
}
/*----------------------------------------------------------------------*/
/* #This function is to solve a 3rd order polynominal equation with	*/
/* four coefficients a, b, c, d. The input coefficients comes from 	*/
/* double *coes and the roots provided in double *roots.		*/
/* 									*/
/* #This function is based on Wikipedia: cubic_function.		*/
/*----------------------------------------------------------------------*/
void cubiSolver(double *coes, double *roots){
	double a, b, c , d;
	double delt, delt0, delt1;
	double complex C0, C;
	double complex u[3], croots[3];
	int i;
	a = coes[0];
	b = coes[1];
	c = coes[2];
	d = coes[3];
	if(abs(a) < 1e-12){
		printf("This is not a 3rd order polynominal equation!\n");
		exit(0);
	}
	delt = 18.0 * a * b * c * d - 4.0 * b * b * b * d + b * b * c * c - 4.0 * a * c * c * c - 27.0 * a * a * d * d;
	if(delt < 1e-12){
		printf("The parameters may not be correct!\n");
		exit(0);
	}
	delt0 = b * b - 3.0 * a * c;
	delt1 = 2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d;
	C0 = cpow((delt1 * delt1 - 4.0 * delt0 * delt0 * delt0), 0.5);
	C = cpow((delt1 + C0) / 2.0, 1.0 / 3);
//	printf("%f, %f, %f + %fi, %f + %fi\n", delt0, delt1, creal(C0), cimag(C0), creal(C), cimag(C));
	u[0] = 1.0;
	u[1] = (-1.0 + pow(3, 0.5) * I) * 0.5;
	u[2] = (-1.0 - pow(3, 0.5) * I) * 0.5;
//	printf("%f+i%f\n", creal(u[2]), cimag(u[2]));
	for(i = 0; i < 3; i++){
		croots[i] = -1.0 / (3 * a) * (b + u[i] * C + delt0 / (u[i] * C));
		roots[i] = creal(croots[i]);
//		printf("root %d: %f\n", i, roots[i]);
	}
	qsort(roots, 3, sizeof(double), cmpfunc);
}

int cmpfunc (const void * a, const void * b)
{
   return ( *(double*)a - *(double*)b );
}
