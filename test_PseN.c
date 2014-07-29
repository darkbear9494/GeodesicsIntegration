#include <stdio.h>
#include <math.h>
#include "PseN.h"
#include "dopri8.h"
#include <string.h>

#define STEPS 500	// Total integration steps
#define HACCU 1000	// The at-least accuracy of h
#define NVAR 4		// The dimension of y

#define PARNUM 1	// Particle number
#define PI 3.141592653589793
#define PARSEC (3.08567758e16)	// parsec
#define ASRUNI 149597870700.0	// AU 
#define LIGHTSPEED 299792458.0	// speed of light
#define NEWTON_G (6.67384e-11)
#define MSUN (1.9891e30)
#define RSUN (6.955e8)
	
int main(){
	double lightspeed, Newton_G, mass;
	double R_sch; 
	double xp[STEPS], yp[STEPS][NVAR], Ene[STEPS];
	double apo, per, Va, Vp, E, AngMom, Vc, ecc, mjAxis;
	double rIni, vIni, T;
	int i, j, k, loop;
	double y1[NVAR], y[NVAR], Ene1, x1, x2, xbin, xend, delX;
	double eps, h0, hmax;
	inData data;
	void *custom_data = &data;
	int flag, nvar;
	fluxfn fn;
	fluxEn En;
	FILE *fp;
	char *str;
	char str1[] = "result_NT";
	char str2[] = "result_WG";
	char str3[] = "result_PW";
	char str4[] = "result_GN";

// Integration-----------------------------------------------------------
	lightspeed = 10;
	Newton_G = 1;
	mass = 1e6;
	R_sch = 2 * Newton_G * mass / dSqr(lightspeed);

	per = 5.0 * R_sch;
	ecc = 0.998;
	apo = (1 + ecc) / (1 - ecc) * per;
	mjAxis = (per + apo) / 2.0;
	E = -0.5 * Newton_G * mass / mjAxis;
	Va = pow(2 * (E + (Newton_G * mass / apo)), 0.5);
	Vp = pow(2 * (E + (Newton_G * mass / per)), 0.5);

	rIni = apo;
	vIni = Va;
//	rIni = per;
//	vIni = Vp;
	AngMom = rIni * vIni;
	T = 2.0 * PI * pow(mjAxis, 1.5) * pow((Newton_G * mass), -0.5);
	Vc = pow((Newton_G * mass / rIni), 0.5);
	
/*	rIni = 10000.0;
	vIni = 4;
	E = 0.5 * vIni*vIni - (Newton_G * mass / rIni);
	AngMom = rIni * vIni;
	mjAxis = -1.0 * Newton_G * mass / (2 * E);
	ecc = pow((1.0 - dSqr(AngMom) / (Newton_G * mass * mjAxis)), 0.5);
	apo = mjAxis * (1.0 + ecc);
	per = mjAxis * (1.0 - ecc);
	T = 2.0 * PI * pow(mjAxis, 1.5) * pow((Newton_G * mass), -0.5);
	Vc = pow((Newton_G * mass / rIni), 0.5);
*/

for(loop = 1; loop <= 4; loop++){	
	switch (loop){
		case 1:
			fn = F_NT;
			En = E_NT;
			str = &str1[0];
			break;
		case 2:
			fn = F_WG;
			En = E_WG;
			str = &str2[0];
			break;
		case 3:
			fn = F_PW;
			En = E_PW;
			str = &str3[0];
			break;
		case 4:
			fn = F_GN;
			En = E_GN;
			str = &str4[0];
			break;
	}
	printf("%s\n", str);


	x1 = 0;
	x2 = 6.0 * T;
	delX = (x2 - x1) / STEPS;
	eps = 1e-13;
	hmax = delX / HACCU;
	h0 = 1e-10;
	nvar = NVAR;
	data.mass = mass;
	data.lightspeed = lightspeed;
	data.Newton_G = Newton_G;

	for(i = 0; i < NVAR; i++){
		y1[i] = 0;
	}
	y1[1] = -1.0 * rIni;
	y1[2] = vIni;
	for(i = 0; i < NVAR; i++){
		y[i] = y1[i];
	}
	Ene1 = En(nvar, y, custom_data);
// Run--------------------------------------------------------------
	for(i = 0; i < STEPS; i++){
		xbin = x1 + delX * i;
		xend = x1 + delX * (i + 1);
		flag = dopri8(fn, nvar, xbin, y, xend, eps, hmax, &h0, NULL, NULL, custom_data);
		if(flag == 0) printf("!!%d\n", i);
		Ene[i] = En(nvar, y, custom_data);
		for(j = 0; j < NVAR; j++) yp[i][j] = y[j];
		xp[i] = xend;
	}	

	fp = fopen(str,"w+");
	if(fp == NULL)
		printf("Cannot open the file!\n");
	
	fprintf(fp, "%10.6f \t", x1);
	for(j = 0; j < NVAR; j++){
		fprintf(fp, "%20.15f \t", y1[j]);
	}
	fprintf(fp, "%20.15f \n", Ene1);
	for(i = 0; i < STEPS; i++){
		fprintf(fp, "%10.6f \t", xp[i]);
		for(j = 0; j < NVAR; j++){
			fprintf(fp, "%20.15f \t", yp[i][j]);
		}
		fprintf(fp, "%20.15f \n", Ene[i]);
	}
	fclose(fp);
}
	
	printf("END-----------------------------------\n");
	printf("R_sch = %e, Vc = %e\n", R_sch, Vc);
	printf("apo = %e, per = %e, ecc = %e, mjAxis = %e\n", apo, per, ecc, mjAxis);
	printf("Va = %e, Vp = %e\n", Va, Vp);
	printf("Initial quantities:\n");
	printf("rIni = %e, vIni = %e, T = %e\n", rIni, vIni, T);
	printf("E = %e, apo = %e\n", E, apo);
}
