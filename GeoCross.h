#include "dopri8.h"

// #define MPIRUN
// Parameters for integration.
#define CYCLE 1		// The total cycle of the circular period at 2Rp.
#define MOVES 1000	// Total moves for the whole simulation.
#define STEPS 1		// Total integration steps
#define HACCU 1000	// The at-least accuracy of h
#define NVAR 4		// The dimension of y
#define FORCE F_SG// The force function used in geoIntegrator.c
#define ENERGY E_SG	// The energy function used in geoIntegrator.c

// Parameters for simulation. (The units are in SI)
#define PARNUM 1	// Particle number
#define PI 3.141592653589793
#define PARSEC (3.08567758e16)	// parsec
#define ASRUNI 149597870700.0	// AU 
#define LIGHTSPEED 299792458.0	// speed of light
#define NEWTON_G (6.67384e-11)	// gravity constant
#define MSUN (1.9891e30)	// solar mass
#define RSUN (6.955e8)		// solar radius

// struct for physical parameters.
typedef struct phyPar{
	double Msun;
	double Mstar;
	double Rstar;
	double BHmass;
	double Rtidal;
	double Rp;
	double Tratio;
	double Vcsun;
} phypar;

// normalized parameters.
typedef struct norPar{
	double mstar;
	double rstar;
	double bhmass;
	double rtidal;
	double rp;
	double zetaMin;
	double tmin;
	double newton_g;
	double lightspeed;
} norpar;

// tidal disruption parameters.
typedef struct tidPar{
	double Ezeta; // This ought to be a function.
	double AngMom;
	double Rtidal;
	double beta;
} tidpar;

// orbital parameters.
typedef struct orbPar{
	int cstep;	// The current steps. Initial value of csteps = -1.
	int* active;	// save the colliding times; inactive if active=-1;
	double* xp;	// Point to the vector saving time information.
	double** yp;	// Point to the vectors saving orbital parameters.
	double* Ene;	// Point to the vector saving total energy for conservation check.
} orbpar;

typedef struct pArticle{
	int No;
	double Md_org;
	double Md;
	double Eini;
	double tini;
	orbpar orbdata;
} particle;

typedef struct sHockList{
	int pN1, pN2; // Particle number 1 and 2.
	double t; // The time of particle collision.
	double y1[NVAR], y2[NVAR]; // The position and velocity of the two particles.
	double x, y;
	double Erel;
} shocklist;

typedef struct sHockHead{
	int clength; // The length of the recorded list.
	shocklist *sl; // Pointing to the shock list.
} shockhead;

// The functional of energies.
typedef double (*fluxEn)(int ndim, double *y, void *p);

void iniQuant(double, double, double, double, double, phypar *, norpar *);
void iniTipar(norpar , double , double *, double **, double *, int , int );
void iniTipar_toy(norpar , double , double *, double **, double *, int , int );
void iniParlist(int ,int ,int ,particle *, shocklist *, double *, double **,double *);
void geoIntegrator(norpar *, orbpar *, double , double , double* , int , int , int , fluxfn , fluxEn);
void geoMove(norpar *,particle *, double , double , int ,int , int , int , fluxfn , fluxEn );
void geoMove_toy(norpar *,particle *, double , double , int ,int , int , int , fluxfn , fluxEn );
int geoSave(particle *, shockhead *, int, int, int);
void geoCollide(norpar *, particle *, shockhead *, int , int );

void F_toy(double t, double *y, double *y1, void *p);
double E_toy(int ndim, double *y, void *p);
