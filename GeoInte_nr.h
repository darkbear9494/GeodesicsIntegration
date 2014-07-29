void odeint(float ystart[], int nvar, float x1, float x2,
	float eps, float h1, float hmin, int *nok, int *nbad,
	void (*derivs)(float, float [], float []),
	void (*rkqs)(float [], float [], int, float *, float, float,
	float [], float *, float *, void (*)(float, float [], float [])));
void bsstep(float y[], float dydx[], int nv, float *xx, float htry, float eps,
	float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void mmid(float y[], float dydx[], int nvar, float xs, float htot,
	int nstep, float yout[], void (*derivs)(float, float[], float[]));
void pzextr(int iest, float xest, float yest[], float yz[], float dy[],
	int nv);
void dyOvdx(float x, float *y, float *dydx);
