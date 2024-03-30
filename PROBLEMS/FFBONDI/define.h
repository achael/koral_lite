/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM -1
#define BHDISK_PROBLEMTYPE

/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1    //must be 1 for MAGNFIELD

//#define VECPOTGIVEN // do not use, instead, field defined manually
//#define MONOPOLE_FIELD_CORNERS
//#define INIT_MAGN_CORNERS

/************************************/
//blackhole
/************************************/
#define MASS 1.
#define BHSPIN 0.
#define RHOR 2

/************************************/
//init
/***********************************/
#define MDOTINIT 1.
#define RSONIC 8.
#define GAMMA (4./3.)
#define SIGMACRIT 25.//10.
#define BONDICOORDS BLCOORDS

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 2
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .5
#define FLUXLIMITER 1
#define FLUXMETHOD LAXF_FLUX
#define MINMOD_THETA 1.5

/************************************/
//mhd floors
/************************************/

//#define NOLOGINS // AC TODO -- needed? 
//#define ENFORCEENTROPY
//#define CORRECT_POLARAXIS
//#define NCCORRECTPOLAR 2
//#define SKIPALLFLOORS 

#define FORCEFREE

#if defined(FORCEFREE)
#define HYBRID_FORCEFREE

#define FORCEFREE_SOLVE_PARALLEL
#define FORCEFREE_PARALLEL_COLD
//#define FORCEFREE_PARALLEL_ENTROPY

#if defined(HYBRID_FORCEFREE)
#define HYBRID_FORCEFREE_SIGMACUT 50
#define HYBRID_FORCEFREE_WIDTH 0.05
#endif


//#define B2RHOFLOORFRAME DRIFTFRAME // not used
#define UURHORATIOMIN 0.
#define UURHORATIOMAX 100. 
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 1.e100
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 1.e100

#define GAMMAMAXFF 100//50//100.  //lower than GAMMAMAXHD? 
#define GAMMAMAXHD 100//100. //why can't this be pushed higher on the monopole? 

#else

#define B2RHOFLOORFRAME DRIFTFRAME
#define UURHORATIOMIN 0.
#define UURHORATIOMAX 100.
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 1.e100
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 1.e100

#define GAMMAMAXHD 100.//5. // TODO this is important for keeping things sane at horizon
#endif


/************************************/
//coordinates / resolution
/************************************/
#define myMKS1COORDS
#define METRICAXISYMMETRIC
#define RMIN 0.7*RHOR //1.8<->6 ANDREW
#define RMAX 100.//1000.

#define TNX 64
#define TNY 1
#define TNZ 1

//number of tiles
#define NTX 4//48
#define NTY 2//24
#define NTZ 1

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS MKS1COORDS
#define MKSR0 0.//-4*RHOR
#define MINX (log(RMIN - MKSR0))
#define MAXX (log(RMAX - MKSR0))
#define MINY 0.005*M_PI
#define MAXY M_PI - 0.005*M_PI

#elif defined(myMKS2COORDS) //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MKSR0 -4.*RHOR
#define MKSH0 0.8

#define MINX (log(RMIN - MKSR0))
#define MAXX (log(RMAX - MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)

#else //Schwarzschild
#define MYCOORDS SCHWCOORDS
#define MINX RMIN //(1.1*r_horizon_BL(BHSPIN))
#define MAXX RMAX
#define MINY (0.01*Pi/2.)
#define MAXY (Pi-0.01*Pi/2.)
#endif

#define PHIWEDGE (M_PI/2.)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

#define SPECIFIC_BC
#define PERIODIC_ZBC
  

/************************************/
//output
/************************************/
#define OUTCOORDS BLCOORDS
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10

#define DTOUT1 1
#define DTOUT2 1.

#define NOUTSTOP 100
#define SILOOUTPUT 0
#define PRIMOUTPUT 0
#define SIMOUTPUT 2
#define OUTPUTINGU
//#define PRIMOUTPUTINMYCOORDS
#define SCAOUTPUT 0
#define SILO2D_XZPLANE

/************************************/
//common physics / atmosphere
/************************************/
