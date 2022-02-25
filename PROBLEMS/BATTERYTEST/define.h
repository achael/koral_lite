#define SKIPRADEVOLUTION
#define SKIPRADSOURCE

#define MASS (1./MSUNCM) //so that x-coordinate in centimeters


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
//radiation
/************************************/
#define RADIATION

/************************************/
//magnetic fields
/************************************/
#define MAGNFIELD
#define BATTERY
#define BATTERYMULTIPLIER 1.e0

/************************************/
//coordinates / resolution
/************************************/
#define MYCOORDS MINKCOORDS
#define LLL 1.
#define MINX -LLL
#define MAXX LLL
#define MINY -LLL
#define MAXY LLL
#define MINZ 0.
#define MAXZ 1.
#define TNX 512
#define TNY 1
#define TNZ 1
#define NTX 1
#define NTY 1
#define NTZ 1

#define PERIODIC_XBC
#define PERIODIC_YBC
#define PERIODIC_ZBC

/************************************/
//Output
/************************************/
#define OUTCOORDS MYCOORDS
#define SIMOUTPUT 2
#define DTOUT1 1.e1 //res
#define DTOUT2 1.e20 //avg

/************************************/
//reconstruction / stepping
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2HEUN
#define TSTEPLIM .5
#define FLUXLIMITER 0
#define MINMOD_THETA  1.5
#define NOUTSTOP 30


/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-10
#define B2RHORATIOMAX 100.
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-15
#define EERHORATIOMAX 1.e15
#define ERADLIMIT 1.e-50
#define GAMMAMAXRAD 20.
#define GAMMAMAXHD 20.

/************************************/
//physics
/************************************/
#define GAMMA (5./3.)
#define HFRAC 1.
#define HEFRAC 0.

/************************************/
//problem params
/************************************/
#define RHOINIT rhoCGS2GU(1.)
#define UUINIT endenCGS2GU(1.*CCC*CCC*1.e-3)
#define ERADINIT UUINIT
#define VGASINIT 0.
#define VRADINIT (2.5*10./9./0.99792*1e-6)
