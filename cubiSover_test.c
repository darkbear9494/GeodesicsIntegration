#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

void cubiSolver(double *coes, double *roots);
int cmpfunc (const void * a, const void * b);

int main(){
	double roots[3];
	double coes[4];

	coes[0] = -1;
	coes[1] = 118.11;
	coes[2] = 70301.0/50;
	coes[3] = -15609.0/100;

	cubiSolver(coes, roots);
	printf("%f, %f, %f\n", roots[0], roots[1], roots[2]);
}
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
	printf("%f, %f, %f + %fi, %f + %fi\n", delt0, delt1, creal(C0), cimag(C0), creal(C), cimag(C));
	u[0] = 1.0;
	u[1] = (-1.0 + pow(3, 0.5) * I) * 0.5;
	u[2] = (-1.0 - pow(3, 0.5) * I) * 0.5;
//	printf("%f+i%f\n", creal(u[2]), cimag(u[2]));
	for(i = 0; i < 3; i++){
		croots[i] = -1.0 / (3 * a) * (b + u[i] * C + delt0 / (u[i] * C));
		roots[i] = creal(croots[i]);
		printf("root %d: %f\n", i, roots[i]);
	}
	qsort(roots, 3, sizeof(double), cmpfunc);
	for(i = 0; i < 3; i++){
		printf("root %d: %f\n", i, roots[i]);
	}
}

int cmpfunc (const void * a, const void * b)
{
   return ( *(double*)a - *(double*)b );
}
