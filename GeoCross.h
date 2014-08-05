#include "dopri8.h"

// #define MPIRUN
// Parameters for integration.
#define CYCLE 500	// The total cycle of the circular period at 2Rp.
#define MOVES 100	// Total moves for the whole simulation.
#define STEPS 10 // Total integration steps
#define HACCU 1000	// The at-least accuracy of h
#define NVAR 4		// The dimension of y
#define FORCE F_SG	// The force function used in geoIntegrator.c
#define ENERGY E_SG	// The energy function used in geoIntegrator.c

// Parameters for simulation. (The units are in SI)
#define PARNUM 3	// Particle number
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
	double newton_g;
	double lightspeed;
} norpar;

// orbital parameters.
typedef struct orbPar{
	int cstep;	// The current steps. Initial value of csteps = -1.
	double* xp;	// Point to the vector saving time information.
	double** yp;	// Point to the vectors saving orbital parameters.
	double* Ene;	// Point to the vector saving total energy for conservation check.
} orbpar;

// tidal disruption parameters.
typedef struct tidPar{
	double Ezeta; // This ought to be a function.
	double AngMom;
	double Rtidal;
	double beta;
} tidpar;

typedef struct pArticle{
	int No;
	int active;
	double Md_org;
	double Md;
	double Eini;
	double tini;
	orbpar orbdata;
} particle;

// The functional of energies.
typedef double (*fluxEn)(int ndim, double *y, void *p);

void geoIntegrator(norpar *, orbpar *, double , double , double* , int , int , int , fluxfn , fluxEn);
void iniParlist(int ,int ,int ,particle *, double *, double **,double *);
void geoSave(FILE *, particle *, int, int, int);
void geoMove(norpar *,particle *, double , double , int ,int , int , int , fluxfn , fluxEn );
void iniTipar(norpar , double , double *, double **, double *, int , int );
void geoCollide(particle *, int , int );

