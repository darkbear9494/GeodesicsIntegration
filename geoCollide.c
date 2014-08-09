#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GeoCross.h"
#include "PseN.h"

typedef struct Point{
	double x, y;
} point;

int colliCheck(particle *, particle *, point *, int);
int intersect(point , point , point , point , point *);
void strongShock(norpar *, particle *, particle *, point , shockhead *);


void geoCollide(norpar *nordata, particle *parlist, shockhead *sh, int parnum, int steps){
	int i, j;
	int cstep, act1, act2;
	point xx;
	particle *p1, *p2;
	for(i = 0; i < (parnum - 1); i++){
		p1= parlist + i;
		cstep = p1->orbdata.cstep;
		act1 = p1->orbdata.active[cstep - 1];
		if(act1 != -1){
			for(j = i + 1; j < parnum; j ++){
				p2= parlist + j;
//				cstep = p2->orbdata.cstep;
				act2 = p2->orbdata.active[cstep - 1];
				if((act2 != -1) && (colliCheck(p1, p2, &xx, steps))){
				//printf("geoCollide: loop (%d, %d)\n", i, j);
					strongShock(nordata, p1, p2, xx, sh);
				}
			}
		}
	}
}
//------------------------------------------------------------------------------//
int colliCheck(particle *particle1, particle *particle2, point *xx, int steps){
	double *y1p, *y1, *y2p, *y2;
	point aa, bb, cc, dd;
	int cstep, i, j, k;

	cstep = particle1->orbdata.cstep;
	for(i = 0; i < steps; i++){
		j = cstep - steps + i;
		y1p = particle1->orbdata.yp[j];
		y2p = particle2->orbdata.yp[j];
		y1 = particle1->orbdata.yp[j + 1];
		y2 = particle2->orbdata.yp[j + 1];
		
		aa.x = y1p[0];
		aa.y = y1p[1];
		bb.x = y1[0];
		bb.y = y1[1];
		cc.x = y2p[0];
		cc.y = y2p[1];
		dd.x = y2[0];
		dd.y = y2[1];
		if(intersect(aa, bb, cc, dd, xx)){
			printf("colliCheck: Find collide!\n");
			printf("colliCheck: step=%d\n", j+1);
			printf("colliCheck: x=%e, y=%e\n", xx->x, xx->y);
			return 1;
		}
	}

	return 0 ;
}


double determinant(double v1, double v2, double v3, double v4)  // determinant of matrix {v1, v2; v3, v4}.
{
	return (v1*v4-v2*v3);
}

int intersect(point aa, point bb, point cc, point dd, point *xx)
{
	double eps = 1e-10;
	double delta, namda, miu;
	delta = determinant(bb.x-aa.x, cc.x-dd.x, bb.y-aa.y, cc.y-dd.y);
//	printf("delta=%f\n", delta);
	if (fabs(delta) <= eps)  // delta=0，parallel or coincidant.
	{
		return 0;
	}
	namda = determinant(cc.x-aa.x, cc.x-dd.x, cc.y-aa.y, cc.y-dd.y) / delta;
//	printf("namda=%f\n", namda);
	if ( namda>1 || namda<0 )
	{
		return 0;
	}
	miu = determinant(bb.x-aa.x, cc.x-aa.x, bb.y-aa.y, cc.y-aa.y) / delta;
//	printf("miu=%f\n", miu);
	if ( miu>1 || miu<0 )
	{
		return 0;
	}
	xx->x = aa.x + namda * (bb.x - aa.x);
	xx->y = aa.y + namda * (bb.y - aa.y);
	return 1;
}
//------------------------------------------------------------------------------//
void strongShock(norpar *nordata,
particle *particle1,
particle *particle2,
point xx,
shockhead *sh){
	int cstep, clength, No1, No2;
	int i, nvar;
	
	particle *p1, *p2;
	shocklist *sl = sh->sl;
	double Md1, Md2, Md_fin;
	double x, y1[NVAR], y2[NVAR], y_fin[NVAR];
	double E1, E2, E_fin, Erel; // Energy released.
	
	nvar = NVAR;
	p1 = particle1;
	p2 = particle2;
	No1 = p1->No;
	No2 = p2->No;
	Md1 = p1->Md;
	Md2 = p2->Md;
	cstep = p1->orbdata.cstep;
	x = p1->orbdata.xp[cstep];
	E1 = p1->orbdata.Ene[cstep];
	E2 = p2->orbdata.Ene[cstep];
	for(i = 0; i < nvar; i++){
		y1[i] = p1->orbdata.yp[cstep][i];
		y2[i] = p2->orbdata.yp[cstep][i];
	}

	Md_fin = Md1 + Md2;
	y_fin[0] = xx.x;
	y_fin[1] = xx.y;
	y_fin[2] = (Md1 * y1[2] + Md2 * y2[2]) / Md_fin;
	y_fin[3] = (Md1 * y1[3] + Md2 * y2[3]) / Md_fin;
	
	fluxEn En = ENERGY;
	inData data;
	void *custom_data = &data;
	data.mass = nordata->bhmass;
	data.lightspeed = nordata->lightspeed;
	data.Newton_G = nordata->newton_g;
	E_fin = En(nvar, y_fin, custom_data);
	Erel = E1 + E2 - E_fin;
	
	sh->clength += 1;
	clength = sh->clength;
	sl[clength].pN1 = No1;
	sl[clength].pN2 = No2;
	sl[clength].t = x;
	sl[clength].x = xx.x;
	sl[clength].y = xx.y;
	sl[clength].Erel = Erel;
	for(i = 0; i < nvar; i++){
		sl[clength].y1[i] = y1[i];
		sl[clength].y2[i] = y2[i];
	}
	
	p1->Md = Md_fin;
	p1->orbdata.active[cstep] += 1;
	p2->orbdata.active[cstep] = -1;
	p1->orbdata.xp[cstep] = x;
	p1->orbdata.Ene[cstep] = E_fin;
	for(i = 0; i < nvar; i++)
		p1->orbdata.yp[cstep][i] = y_fin[i];
	
}

