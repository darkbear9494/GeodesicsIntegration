//#define NUME

// Parameters for integration.
#define STEPS 500	// Total integration steps
#define HACCU 1000	// The at-least accuracy of h
#define NVAR 4		// The dimension of y

// Parameters for simulation. (The units are in SI)
#define PARNUM 1	// Particle number
#define PI 3.141592653589793
#define PARSEC (3.08567758e16)	// parsec
#define ASRUNI 149597870700.0	// AU 
#define LIGHTSPEED 299792458.0	// speed of light
#define NEWTON_G (6.67384e-11)
#define MSUN (1.9891e30)
#define RSUN (6.955e8)

/*
#ifdef NUME
	#define NEWTON_G 1
	#define MSUN 1
	#define RSUN 1
//	#define BHMASS (1e6)	// BH mass
#else 
	#define NEWTON_G (6.67384e-11)
	#define MSUN (1.9891e30)
	#define RSUN (6.955e8)
//	#define BHMASS (1e6 * MSUN)	// BH mass
#endif
*/
