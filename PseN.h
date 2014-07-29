#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct indata{
	double mass;
	double lightspeed;
	double Newton_G;
} inData;
typedef double (*fluxEn)(int ndim, double *y, void *p);
//----------------------------------------------------------------------//
// There are four kind of potential with their corresponding force.
void F_NT(double x, double *y, double *y1, void *p);
double E_NT(int ndim, double *y, void *p);
// Newtonian
void F_GN(double t, double *y, double *y1, void *p);
double E_GN(int ndim, double *y, void *p);
// Generalized Newtonian (Tejeda&Rosswog 2013)
void F_WG(double t, double *y, double *y1, void *p);
double E_WG(int ndim, double *y, void *p);
// Wegg (Wegg 2012)
void F_PW(double t, double *y, double *y1, void *p);
double E_PW(int ndim, double *y, void *p);
//Paczinsky-Witta
void F_WG(double t, double *y, double *y1, void *p);
double E_WG(int ndim, double *y, void *p);
// Wegg (Wegg 2012)

//----------------------------------------------------------------------//
double calR(double x, double y, double z);
double dSqr(double x);
double dCb(double x);
int EPSNijk(int i, int j, int k);
