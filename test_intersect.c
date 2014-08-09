#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Point{
	double x, y;
} point;

int main(){
	point p1, p2, p3, p4, p5;
	p1.x = 0;
	p1.y = -2;

	p2.x = 0;
	p2.y = 3.09e-13;

	p3.x = -2;
	p3.y = 0;

	p4.x = 3.09e-13;
	p4.y = 0;

	if(intersect(p1, p2, p3, p4, &p5)){
		printf("Find!\n");
		printf("%f, %f\n", p5.x, p5.y);
	}
	else printf("Not intersection found!\n");
}
