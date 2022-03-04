/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM -1

//#define RADIATION


/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1    //must be 1 for MAGNFIELD
#define VECPOTGIVEN
//#define INIT_MAGN_CORNERS

#define FORCEFREE
//#define NOLOGINS 

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1 //TODO why is int_order 1 more stable???
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .9
#define FLUXLIMITER 1
#define FLUXMETHOD LAXF_FLUX
#define MINMOD_THETA 1.5

/************************************/
//rmhd floors
/************************************/

#define B2RHORATIOMAXINIT 500 
#define B2UURATIOMAXINIT 500
//#define SIGMAWCONSTINIT 1.e4

#if defined(FORCEFREE)

//DIFTFRAME not compatible with FORCEFREE yet
#define B2RHOFLOORFRAME FFFRAME // ZAMOFRAME 

//#define CORRECT_POLARAXIS
//#define NCCORRECTPOLAR 2

#define UURHORATIOMIN 0.
#define UURHORATIOMAX 50.
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 1.e10
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 1.e10

#define GAMMAMAXFF 1000. //lower than GAMMAMAXHD? 
#define GAMMAMAXHD 10000. //why can't this be pushed higher on the monopole? 

#define ALLOWENTROPYU2P 0

#else
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2

#define B2RHOFLOORFRAME DRIFTFRAME //ZAMOFRAME 
#define UURHORATIOMIN 0.
#define UURHORATIOMAX 50.
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 2500.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 500.
#define GAMMAMAXHD 50.
#endif

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.9375
#define RHOR (1.+sqrt(1. - BHSPIN*BHSPIN))

/************************************/
//coordinates / resolution
/************************************/
#define myMKS1COORDS
#define METRICAXISYMMETRIC
#define RMIN 0.7*RHOR //1.8<->6 ANDREW
#define RMAX 100.
#define MKS1R0 MKSR0

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS MKS1COORDS
#define MINX (log(RMIN - MKSR0))
#define MAXX (log(RMAX - MKSR0))
#define MINY 0.01*M_PI/2.
#define MAXY M_PI - 0.01*M_PI/2.
#define TNX 128
#define TNY 64//256
#define TNZ 1

#elif defined(myMKS2COORDS) //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MKSR0 -4.*RHOR
#define MKSH0 0.8

#define MINX (log(RMIN - MKSR0))
#define MAXX (log(RMAX - MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)
#define TNX 128
#define TNY 128
#define TNZ 1

#else //Schwarzschild
#define MYCOORDS SCHWCOORDS
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX (25.3)
#define MINY (0.01*Pi/2.)
#define MAXY (Pi-0.01*Pi/2.)
#define TNX 64
#define TNY 32
#define TNZ 1
#endif

#define PHIWEDGE (M_PI/2.)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

#define SPECIFIC_BC
#define PERIODIC_ZBC
  
//number of tiles
#define NTX 4//48
#define NTY 2//24
#define NTZ 1

/************************************/
//output
/************************************/
#define OUTCOORDS KSCOORDS //MKS1COORDS
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 100
#define SILOOUTPUT 1
#define PRIMOUTPUT 1
//#define PRIMOUTPUTINMYCOORDS
#define SCAOUTPUT 0
#define SILO2D_XZPLANE

/************************************/
//common physics / atmosphere
/************************************/
#define GAMMA (4./3.)
#define DTOUT1 1
#define DTOUT2 1.e0
#define RHOATMMIN  1.e-20
#define UINTATMMIN 1.e-20
