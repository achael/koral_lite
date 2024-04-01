/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM -1
#define BHDISK_PROBLEMTYPE

//#define RADIATION


/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1    //must be 1 for MAGNFIELD
#define VECPOTGIVEN
#define MONOPOLE_FIELD_CORNERS
//#define INIT_MAGN_CORNERS

//#define U2PCONV 1.e-14


/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1 //TODO why is int_order 1 more stable???
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .75
#define FLUXLIMITER 1
#define FLUXMETHOD LAXF_FLUX
#define MINMOD_THETA 1.5

/************************************/
//mhd floors
/************************************/
//#define SPLIT_MONOPOLE

#define RHOATMMIN  1.e-20
//#define UINTATMMIN 1.e-8//1.e-20

#define B2RHORATIOMAXINIT 1000
//#define B2UURATIOMAXINIT 100
#define RHO2UURATIOMAXINIT 1000 

//#define NOLOGINS
#define ENFORCEENTROPY

//#define FORCEFREE
//#define HYBRID_FORCEFREE

#if defined(FORCEFREE)
#define FORCEFREE_SOLVE_PARALLEL
#define FORCEFREE_PARALLEL_COLD
//#define FORCEFREE_PARALLEL_ENTROPY

#if defined(HYBRID_FORCEFREE)
#define HYBRID_FORCEFREE_SIGMACUT 100//50
#define HYBRID_FORCEFREE_WIDTH 0.2
#endif //HYBRID_FORCEFREE

//#define SKIPALLFLOORS // TODO seems critical for SOLVE_PARALLEL? 

//TODO
//#define CORRECT_POLARAXIS //messes up omegaH at the axis
//#define NCCORRECTPOLAR 1//2

#define B2RHOFLOORFRAME DRIFTFRAME // not used
#define UURHORATIOMIN 0.
#define UURHORATIOMAX 100. 
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 1.e100
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 1.e100

#define GAMMAMAXFF 2000//50//100.  //lower than GAMMAMAXHD? 
#define GAMMAMAXHD 2000//100. //why can't this be pushed higher on the monopole? 

#else

#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2

#define B2RHOFLOORFRAME DRIFTFRAME
#define UURHORATIOMIN 0.
#define UURHORATIOMAX 100.
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 1000//B2UURATIOMAXINIT
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 1000//B2RHORATIOMAXINIT

#define GAMMAMAXHD 50.//100.
#endif

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.5//9375
#define RHOR (1.+sqrt(1. - BHSPIN*BHSPIN))

/************************************/
//coordinates / resolution
/************************************/
#define myMKS2COORDS
#define METRICAXISYMMETRIC
#define RMIN 0.7*RHOR //1.8<->6 ANDREW
#define RMAX 250. //200.

#define TNX 192//128//256
#define TNY 128//256
#define TNZ 1

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS MKS1COORDS
#define MKSR0 0.
#define MINX (log(RMIN - MKSR0))
#define MAXX (log(RMAX - MKSR0))
#define MINY 0.005*M_PI
#define MAXY M_PI - 0.005*M_PI

#elif defined(myMKS2COORDS) //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MKSR0 0.//-4.*RHOR
#define MKSH0 0.8

#define MINX (log(RMIN - MKSR0))
#define MAXX (log(RMAX - MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)

#else //Schwarzschild
#define MYCOORDS SCHWCOORDS
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX (25.3)
#define MINY (0.01*Pi/2.)
#define MAXY (Pi-0.01*Pi/2.)
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
#define OUTCOORDS2 KSCOORDS
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 100
#define SILOOUTPUT 1
//#define PRIMOUTPUT 1
#define ANAOUT_HDF5
//#define PRIMOUTPUTINMYCOORDS
#define SCAOUTPUT 0
#define SILO2D_XZPLANE

/************************************/
//common physics / atmosphere
/************************************/
#define GAMMA (4./3.)
#define DTOUT1 1
#define DTOUT2 1.e0
