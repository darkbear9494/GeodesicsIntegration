#include "dopri8.h"

int dopri8(fluxfn calc, int n, double x, double* y, double xend, double eps, double hmax, double* h0, intout out, intfinal final, void* custom_data);

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

