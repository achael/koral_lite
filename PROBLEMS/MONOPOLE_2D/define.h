/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM 44

/************************************/
//radiation choices
/************************************/
//#define RADIATION
//#define SKIPRADSOURCE

/************************************/
//magnetic choices
/************************************/
//#define MAGNFIELD
//#define VECPOTGIVEN
//#define MONOPOLE_FIELD_CORNERS // generates monopole field potential directly on corners
                               // no cell averaging

//#define MAXBETA 0.01 //close to the target pmag/pgas

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .8
#define FLUXLIMITER 0
#define FLUXMETHOD LAXF_FLUX
#define MINMOD_THETA 1.5

/************************************/
//viscosity choices
/************************************/
#define HDVISCOSITY NOVISCOSITY
#define RADVISCOSITY NOVISCOSITY

/************************************/
//rmhd floors
/************************************/
#define UURHORATIOMIN 1.e-15
#define UURHORATIOMAX 50.
//#define EERHORATIOMIN 1.e-20
//#define EERHORATIOMAX 1.e20
//#define EEUURATIOMIN 1.e-20
//#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 2500.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 50.
#define GAMMAMAXRAD 50.

/************************************/
//blackhole
/************************************/
#define MASS 10.
//#define BHSPIN 0.
#define BHSPIN 0.9375
#define RHOR (1.+sqrt(1. - BHSPIN*BHSPIN))
/************************************/
//coordinates / resolution
/************************************/
#define myMKS2COORDS
#define MKSR0 -2
#define MKSH0 0.8
#define METRICAXISYMMETRIC
#define RMIN 0.7*RHOR //1.8<->6 ANDREW
#define RMAX 1000.
#define MKS1R0 MKSR0

#ifdef myMKS2COORDS //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN - MKSR0))
#define MAXX (log(RMAX - MKSR0))
//#define MINY (0.1*Pi/2.)
//#define MAXY (Pi-0.1*Pi/2.)
#define MINY (0.001)
#define MAXY (1.-0.001)
#define TNX 512
#define TNY 40
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
//#define MINZ -1.
//#define MAXZ 1.

#define SPECIFIC_BC
#define PEROIDIC_ZBC
  
//number of tiles
#define NTX 4//48
#define NTY 1//24
#define NTZ 1

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define SILOOUTPUT 1
#define SCAOUTPUT 0
#define SILO2D_XZPLANE

/************************************/
//common physics / atmosphere
/************************************/
#define GAMMA (4./3.)
#define NODONUT 0
#define INFLOWING 0
#define ELL 3.85
#define URIN (0.)
#define KKK 7127.
#define UTPOT .965
#define DTOUT1 1.e0
#define DTOUT2 1.e0
#define RHOATMMIN  1.e-20
#define UINTATMMIN 1.e-20
//#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
//#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)
