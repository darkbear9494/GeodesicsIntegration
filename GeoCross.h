#include "dopri8.h"

// struct for physical parameters.
typedef struct phyPar{
	double Mstar;
	double Rstar;
	double BHmass;
	double Rtidal;
	double Tratio;
	double Vcsun;
} phypar;

// normalized parameters.
typedef struct norPar{
	double mstar;
	double rstar;
	double bhmass;
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
	double Eini;
	double tini;
	orbpar orbdata;
} particle;

// The functional of energies.
typedef double (*fluxEn)(int ndim, double *y, void *p);

void geoIntegrator(norpar *, orbpar *, double , double , double* , int , int , int , fluxfn , fluxEn);
void iniPar(int ,int ,int ,particle *, double *, double **,double *);
void geoSave(FILE *, particle *, int, int, int);
void geoMove(norpar *,particle *, double , double , int ,int , int , int , fluxfn , fluxEn );

