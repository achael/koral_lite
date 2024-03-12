#define MASS 1.//(1./MSUNCM) //so that x-coordinate in centimeters
#define LINEARALFVEN
#define PULSEWIDTH 1.//0.5;

/************************************/
//restart
/************************************/
//#define RESTART 
#define RESTARTNUM -1
#define RESTARTGENERALINDICES

#define DOFIXUPS 1
#define DOU2PRADFIXUPS 0
#define DOU2PMHDFIXUPS 1
#define DORADIMPFIXUPS 1

/************************************/
//magnetic fields
/************************************/
#define MAGNFIELD

#define VELPRIM VELR
//#define ENFORCEENTROPY
//#define NOLOGINS

#define FORCEFREE

#ifdef FORCEFREE

#define HYBRID_FORCEFREE

//#define HYBRID_FORCEFREE_SIGMACUT 50
#define HYBRID_FORCEFREE_XCUT 0.5
#define HYBRID_FORCEFREE_WIDTH 0.1

#define FORCEFREE_SOLVE_PARALLEL
#define FORCEFREE_PARALLEL_COLD
//#define FORCEFREE_PARALLEL_ENTROPY
//#define FORCEFREE_PARALLEL_MHD
//#define SKIPALLFLOORS
#endif

/************************************/
//coordinates / resolution
/************************************/
#define MYCOORDS MINKCOORDS
#define LLL 1
#define MINX -2*LLL
#define MAXX 2*LLL
#define MINY -2*LLL
#define MAXY 2*LLL
#define MINZ 0.
#define MAXZ 1.
#define TNX 256//512
#define TNY 1//256
#define TNZ 1
#define NTX 1
#define NTY 1
#define NTZ 1

#define COPY_XBC
#define PERIODIC_YBC
#define PERIODIC_ZBC

/************************************/
//Output
/************************************/
#define OUTCOORDS MYCOORDS
#define SIMOUTPUT 2
#define OUTPUTINGU
#define DTOUT1 0.01 //res
#define DTOUT2 1.e20 //avg

/************************************/
//reconstruction / stepping
/************************************/
#define INT_ORDER 2 //1
#define TIMESTEPPING RK2HEUN
#define TSTEPLIM .1
#define FLUXLIMITER 1
#define MINMOD_THETA 1.5
#define NOUTSTOP 150
//#define NSTEPSTOP 1
//#define ALLSTEPSOUTPUT 1

/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-15
#define B2RHORATIOMAX 1000.
#define UURHORATIOMAX 1000.
#define B2UURATIOMAX 10000.

#define GAMMAMAXHD 100.//500.
#define B2RHOFLOORFRAME FFFRAME 

#define GAMMAMAXFF 100.

/************************************/
//physics
/************************************/
//#define GAMMA (5./3.)
#define GAMMA (4./3.)
#define HFRAC 1.
#define HEFRAC 0.

/************************************/
//problem parameters
/************************************/

#define SIGMAINIT 250

#ifdef SIGMAINIT
#define THETAINIT 0.1//1.//.25
#else
#define RHOINIT 1.
#define PINIT 1.
#define UUINIT PINIT/(GAMMA-1.)
#endif


