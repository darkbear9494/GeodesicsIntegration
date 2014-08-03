#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "PseN.h"

// double lightspeed = 1;
// double Newton_G = 1;
//----------------------------------------------------------------------------//
// This file provides four kinds of potential, Newtonian, Generalized
// Newtonian, Paczinsky-Witta potential and Wegg potential (see the
// announcement of the functions below).
// The quantity G and M should be provided value.
//----------------------------------------------------------------------------//

// Schwarzschild Geodesics------------------------------------------------------//
void F_SG(double t, double *y, double *y1, void *p)
{// This function is to calculate the spacial acceleration of Schwarzschild geodesics.
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;
	double x[3], dx[3];
	double A, B, C; // Temped variables for calculation.
	double r, r_g, M;
	double part1[3], part2[3], part3[3], part4[3]; // 3 parts of the force equation (Eq.A5 in T&S(2013)).
	int i, j, k; // Temp variables for loops.

	x[0] = y[0];
	x[1] = y[1];
	x[2] = 0;
	dx[0] = y[2];
	dx[1] = y[3];
	dx[2] = 0;
// We only consider two dimensional situatioin at present.
	M = data->mass;
	r = calR(x[0], x[1], x[2]);
	r_g = Newton_G * M / dSqr(lightspeed);

	A = 0; // Calculate A.
	for(i = 0; i < 3; i++)
	{
		A += x[i] * dx[i];
	}
	B = 0; 
	C = 0; // Calculate B and C.
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		for(k = 0; k < 3; k++)
		{
			C += EPSNijk(i, j, k) * x[j] * dx[k];
		}
		B += C * C;
	}
	for(i = 0; i < 3; i++)
	{
		part1[i] = -1.0 * Newton_G * M * x[i] / dCb(r) * (1 - 2.0 * r_g / r);
		part2[i] = 2.0 * r_g * dx[i] / (dSqr(r) * (r - 2 * r_g)) * A;
		part3[i] = r_g * x[i] / (dSqr(dSqr(r)) * (r - 2 * r_g)) * dSqr(A);
		part4[i] = -2.0 * r_g * x[i] / (dSqr(r) * dCb(r)) * B;
	}
	
	y1[0] = dx[0];
	y1[1] = dx[1];
	y1[2] = part1[0] + part2[0] + part3[0] + part4[0];
	y1[3] = part1[1] + part2[1] + part3[1] + part4[1];
}

double E_SG(int ndim, double *y, void *p)
{
	double x[3], dx[3]; // The dimension is less than 6;
	double r, r_g, E, M;
// r is radius of the origin. r_g is the gravitational radius. E is the energy should be returned.
	double A, B, C; // Temped variables for calculate E.
	double t, nt, ad;
// t: kinetic energy; nt: Newtonian potential; ad: additional part. 
	int i, j, k; // Temped variables for loop.
	int n = ndim / 2;
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;

	for(i = 0; i < 3; i++)
	{
		x[i] = 0;
		dx[i] = 0;
	}
	for(i = 0; i < n; i++)
	{
		x[i] = y[i];
		dx[i] = y[n + i];
	}

	M = data->mass;
	r = calR(x[0], x[1], x[2]);
	r_g = Newton_G * M / dSqr(lightspeed);

	A = 0; // Calculate A.
	for(i = 0; i < 3; i++)
	{
		A += x[i] * dx[i];
	}
	B = 0; 
	C = 0; // Calculate B and C.
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		for(k = 0; k < 3; k++)
		{
			C += EPSNijk(i, j, k) * x[j] * dx[k];
		}
		B += C * C;
	}

	t = 0.5 * (dSqr(dx[0]) + dSqr(dx[1]) + dSqr(dx[2]));
	nt = -1 * Newton_G * M / r;
	ad = (2 * r_g / (r - 2 * r_g)) * ((r - r_g) / (r - 2 * r_g) * dSqr(A / r) + 0.5 * B / dSqr(r));

	E = t + nt + ad;
	return E;
	
}
// Generalized Newtonian------------------------------------------------------//
void F_GN(double t, double *y, double *y1, void *p)
{// This function is to calculate the direvative of general coordinates using Generalized Newtonian potential according to appendix A in Tejeda & Rosswog (2013).
//	double *mass = (double *)p;
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;
	double x[3], dx[3];
	double A, B, C; // Temped variables for calculation.
	double r, r_g, M;
	double part1[3], part2[3], part3[3]; // 3 parts of the force equation (Eq.A5 in T&S(2013)).
	int i, j, k; // Temp variables for loops.

	x[0] = y[0];
	x[1] = y[1];
	x[2] = 0;
	dx[0] = y[2];
	dx[1] = y[3];
	dx[2] = 0;
// We only consider two dimensional situatioin at present.
	M = data->mass;
	r = calR(x[0], x[1], x[2]);
	r_g = Newton_G * M / dSqr(lightspeed);

	A = 0; // Calculate A.
	for(i = 0; i < 3; i++)
	{
		A += x[i] * dx[i];
	}
	B = 0; 
	C = 0; // Calculate B and C.
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		for(k = 0; k < 3; k++)
		{
			C += EPSNijk(i, j, k) * x[j] * dx[k];
		}
		B += C * C;
	}
	for(i = 0; i < 3; i++)
	{
		part1[i] = -1.0 * Newton_G * M * x[i] / dCb(r) * dSqr(1 - 2.0 * r_g / r);
		part2[i] = 2.0 * r_g * dx[i] / (dSqr(r) * (r - 2 * r_g)) * A;
		part3[i] = -3.0 * r_g * x[i] / (dSqr(r) * dCb(r)) * B;
	}
	
	y1[0] = dx[0];
	y1[1] = dx[1];
	y1[2] = part1[0] + part2[0] + part3[0];
	y1[3] = part1[1] + part2[1] + part3[1];
}

double E_GN(int ndim, double *y, void *p)
{
	double x[3], dx[3]; // The dimension is less than 6;
	double r, r_g, E, M;
// r is radius of the origin. r_g is the gravitational radius. E is the energy should be returned.
	double A, B, C; // Temped variables for calculate E.
	double t, nt, ad;
// t: kinetic energy; nt: Newtonian potential; ad: additional part. 
	int i, j, k; // Temped variables for loop.
	int n = ndim / 2;
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;

	for(i = 0; i < 3; i++)
	{
		x[i] = 0;
		dx[i] = 0;
	}
	for(i = 0; i < n; i++)
	{
		x[i] = y[i];
		dx[i] = y[n + i];
	}

	M = data->mass;
	r = calR(x[0], x[1], x[2]);
	r_g = Newton_G * M / dSqr(lightspeed);

	A = 0; // Calculate A.
	for(i = 0; i < 3; i++)
	{
		A += x[i] * dx[i];
	}
	B = 0; 
	C = 0; // Calculate B and C.
	for(i = 0; i < 3; i++)
	{
		for(j = 0; j < 3; j++)
		for(k = 0; k < 3; k++)
		{
			C += EPSNijk(i, j, k) * x[j] * dx[k];
		}
		B += C * C;
	}

	t = 0.5 * (dSqr(dx[0]) + dSqr(dx[1]) + dSqr(dx[2]));
	nt = -1 * Newton_G * M / r;
	ad = (2 * r_g / (r - 2 * r_g)) * ((r - r_g) / (r - 2 * r_g) * dSqr(A / r) + 0.5 * B / dSqr(r));

	E = t + nt + ad;
	return E;
	
}

// Wegg-----------------------------------------------------------------------//
void F_WG(double t, double *y, double *y1, void *p)
{
//	double *mass = (double *)p;
	double r, r_g, M; // radius from the origin and the gravitational radius.
	double Alpha, R_x, R_y; // coefficients in Wegg potential.
	double coef; // The coefficient before (xi/r), to calculate force.
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;

	M = data->mass;
	r = calR(y[0], y[1], 0);
	r_g = Newton_G * M / dSqr(lightspeed);
	Alpha = -4.0 / 3.0 * (2.0 + sqrt(6.0));
	R_x = (4.0 * sqrt(6) - 9) * Newton_G * M / dSqr(lightspeed);
	R_y = -4.0 / 3.0 * (2.0 * sqrt(6.0) - 3.0) * Newton_G * M / dSqr(lightspeed);
	coef = Alpha * Newton_G * M / dSqr(r) + ( 1.0 - Alpha) * Newton_G * M / dSqr(r - R_x) + 2 * Newton_G * M * R_y / dCb(r);

	y1[0] = y[2];
	y1[1] = y[3];
	y1[2] = -1.0 * coef * (y[0] / r);
	y1[3] = -1.0 * coef * (y[1] / r);
// We only consider two dimensional situation at present.
}

double E_WG(int ndim, double *y, void *p)
{
	double r, r_g, E, M;
	double Alpha, R_x, R_y; // coefficients in Wegg potential.
	double x[3], dx[3]; // The dimension is less than 6;
	int n = ndim / 2;
	int i, j, k; // Temped variables for loop.
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;

	for(i = 0; i < 3; i++)
	{
		x[i] = 0;
		dx[i] = 0;
	}
	for(i = 0; i < n; i++)
	{
		x[i] = y[i];
		dx[i] = y[n + i];
	}
	M = data->mass;
	r = calR(x[0], x[1], x[2]);
	r_g = Newton_G * M / dSqr(lightspeed);
	Alpha = -4.0 / 3.0 * (2.0 + sqrt(6.0));
	R_x = (4.0 * sqrt(6) - 9) * Newton_G * M / dSqr(lightspeed);
	R_y = -4.0 / 3.0 * (2.0 * sqrt(6.0) - 3.0) * Newton_G * M / dSqr(lightspeed);

	E = 0.5 * (dSqr(dx[0]) + dSqr(dx[1]) + dSqr(dx[2])) - (Alpha * Newton_G * M / r + ( 1.0 - Alpha) * Newton_G * M / (r - R_x) + Newton_G * M * R_y / dSqr(r));
	return E;
}

// Newtonian------------------------------------------------------------------//
void F_NT(double t, double *y, double *y1, void *p)
{
// The dimension of the function is 4 at present.
	double r , M;
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;
	
	r = calR(y[0], y[1], 0);
	M = data->mass;
	y1[0] = y[2];
	y1[1] = y[3];
	y1[2] = -1.0 * Newton_G * M / pow(r, 3) * y[0];
	y1[3] = -1.0 * Newton_G * M / pow(r, 3) * y[1];
//	printf("-----%f  %f  %f  %f  \n", y[0], y[1], y[2], y[3]);
//	printf("%f  %f  %f  %f  \n", y1[0], y1[1], y1[2], y1[3]);
}

double E_NT(int ndim, double *y, void *p)
{
	double r, E, M;
	double x[3], dx[3]; // The dimension is less than 6;
	int n = ndim / 2;
	int i, j, k; // Temped variables for loop.
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;

	for(i = 0; i < 3; i++)
	{
		x[i] = 0;
		dx[i] = 0;
	}
	for(i = 0; i < n; i++)
	{
		x[i] = y[i];
		dx[i] = y[n + i];
	}
	M = data->mass;
	r = calR(x[0], x[1], x[2]);
	E = 0.5 * (dSqr(dx[0]) + dSqr(dx[1]) + dSqr(dx[2])) - Newton_G * M / r; 
	return E;
}
// Paczynski&Witta------------------------------------------------------------
void F_PW(double t, double *y, double *y1, void *p)
{
	double r, r_g, M; // radius from the origin and the gravitational radius.
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;
	
	M = data->mass;
	r = calR(y[0], y[1], 0);
	r_g = Newton_G * M / dSqr(lightspeed);

	y1[0] = y[2];
	y1[1] = y[3];
	y1[2] = -1.0 * Newton_G * M / dSqr(r - 2 * r_g) * (y[0] / r);
	y1[3] = -1.0 * Newton_G * M / dSqr(r - 2 * r_g) * (y[1] / r);
// We only consider two dimensional situation at present.
}

double E_PW(int ndim, double *y, void *p)
{
	double r, r_g, E, M;
	double x[3], dx[3]; // The dimension is less than 6;
	int n = ndim / 2;
	int i, j, k; // Temped variables for loop.
	inData *data = (inData *)p;
	double lightspeed = data->lightspeed;
	double Newton_G = data->Newton_G;

	M = data->mass;
	for(i = 0; i < 3; i++)
	{
		x[i] = 0;
		dx[i] = 0;
	}
	for(i = 0; i < n; i++)
	{
		x[i] = y[i];
		dx[i] = y[n + i];
	}
	r = calR(x[0], x[1], x[2]);
	r_g = Newton_G * M / dSqr(lightspeed);
	E = 0.5 * (dSqr(dx[0]) + dSqr(dx[1]) + dSqr(dx[2])) - Newton_G * M / (r - 2 * r_g); 
	return E;
}

// Accessary functions---------------------------------------------------------//
double calR(double x, double y, double z)
{
	return sqrt(x * x + y * y + z * z);
}

double dSqr(double x)
{// x^2
	return x * x;
}

double dCb(double x)
{// x^3
	return x * x * x;
}

int EPSNijk(int i, int j, int k)
{// This function is to get the Levi-civita epsilon symbol.
	if(!(i - j) || !(i - k) || !(k - j))
		return 0;
	int A, B;
	A = ((i % 2) * 2 - 1); // if i = 2, A = -1. if i = 1 or 3, A = 1;
	B = ((j < k) ? 1 : -1); // if j < k, B = 1. Otherwise, B = -1;
	return (A * B);
}
