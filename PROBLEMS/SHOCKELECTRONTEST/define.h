//#define NONRELMHD

/************************************/
//general
/************************************/
//#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
#define RESTART
//#define RESTARTGENERALINDICES
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/

//#define RADIATION
#define SKIPRADSOURCE    //advective only
#define SKIPRADEVOLUTION

#define RADIMPLICITDAMPINGFACTOR 3.
#define RADIMPLICITMAXENCHANGEDOWN 100.
#define RADIMPLICIMAXENCHANGEUP 100.
#define MAXRADIMPDAMPING 1.e-3
#define RADIMPMAXITER 70
#define RADIMPCONVREL 1.e-6
#define RADIMPCONVRELERR .99
#define RADIMPCONVRELENTR 1.e-4
#define RADIMPCONVRELENTRERR 1.e-2
#define OPDAMPINIMPLICIT 0

#define U2P_SOLVER U2P_SOLVER_WP

#define MAXDIFFTRADS 1.e6
#define U2PCONV 1.e-12
#define RADIMPCONV 1.e-10
#define RADIMPEPS 1.e-6

//#define NCOMPTONIZATION
//#define BREMSSTRAHLUNG
//#define SYNCHROTRON
//#define SYNCHROTRONALLOWEDTRATIO 1.e5

//#define BOUNDFREE
//#define KLEINNISHINA
//#define BOUNDELECTRON
//#define COLDOPACITIES
//#define SKIPFANCYOPACITIES

/************************************/
//electron choices
/************************************/

#define EVOLVEELECTRONS
//#define HEATELECTRONS_USEFIT
#define HEATELECTRONS
#define HEATELECTRONS_DELTA .5

#define NOLOGINS
#define NOLOGINS2

#define MIXENTROPIESPROPERLY
//#define DONOTMIXGASENTROPY
//#define DONOTMIXSPECIESENTROPY
#define UPWINDENTROPYMIXING
//#define PARTLYUPWINDENTROPYMIXING
//#define PARTLYRECONSTRUCTEDUPWINDENTROPYMIXING
//#define LIMITFACTORINMIXING 1
#define DONOTLIMITENTRINMIXING
//#define OLDENTROPYMIXING
//#define MIXENTROPIES_CONST_PRESSURE

//updates entri at end of each timestep to satisfy ui = uu - ue
//doesn't seem helpful
//#define UPDATE_ENTRI_CONSTRAINT 

#define UEUINTMINRATIO 1.e-16
#define UIUINTMINRATIO 1.e-16
#define TEMPEMINIMAL 1.e-4
#define TEMPIMINIMAL 1.e-4
#define TEMPEMINIMALFRACTION 1.e-16
#define TEMPIMINIMALFRACTION 1.e-16
#define TEMPEMAXIMALFRACTION 1.e3
#define TEMPIMAXIMALFRACTION 1.e3

/************************************/
//magnetic choices
/************************************/
//#define MAGNFIELD
//#define GDETIN 1

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .8
#define FLUXLIMITER 0
#define FLUXMETHOD LAXF_FLUX
#define MINMOD_THETA 1.8
#define SHUFFLELOOPS 0
#define DOFIXUPS 1
#define DORADFIXUPS 1

/************************************/
//viscosity choices
/************************************/
#define HDVISCOSITY NOVISCOSITY
#define RADVISCOSITY NOVISCOSITY

/************************************/
//rmhd floors
/************************************/
//#define CORRECT_POLARAXIS
//#define POLARAXISAVGIN3D
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-14
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 50.
#define GAMMAMAXRAD 50.
#define GAMMAMAXHD 50.
#define RHOFLOOR 1.e-50

/************************************/
//blackhole
/************************************/
#define MASS (1./MSUNCM*CCC_CGS) //so that time in seconds!
#define BHSPIN 0.

/************************************/
//initial setup
/************************************/
#define TESTNO 100

#if (TESTNO==0) //Ressler+15
#define GAMMA (5./3.)
#define GAMMAE (4./3.) //gamma of electrons
#define GAMMAI (5./3.) //gamma of ions
#define HEATELECTRONSATENDRK2
#define DISSIPATIONFROMGASONLY 
#define NSTEPSTOP 50000.
#define DTOUT1 (1.e1)  //dt for output in seconds
#define RHO_INIT 1.0
#define TE_INIT 30.
#define TI_INIT 30.
#define TEINITFACTOR 1.
#define VEL 1.e-3
#define FORCEGAMMAGASFIXED  
#endif

#if (TESTNO==1) //modified Ressler+15
#define GAMMA (5./3.)
#define GAMMAE (4./3.) //gamma of electrons
#define GAMMAI (5./3.) //gamma of ions
#define HEATELECTRONSATENDRK2
#define DISSIPATIONFROMGASONLY
#define NSTEPSTOP 50000.
#define NOUTSTOP 5
#define DTOUT1 (100)  //dt for output in seconds
#define RHO_INIT 1.0
#define TE_INIT .1
#define TI_INIT 10.
#define TEINITFACTOR 1.
#define FORCEGAMMAGASFIXED
#define VEL 1.e-3
#endif

#if (TESTNO==10) //to test electron heating 
//#define CONSISTENTGAMMA
//#define DISSIPATIONFROMGASONLY 
#define GAMMA (5./3.)
#define GAMMAE (4./3.) //gamma of electrons
#define GAMMAI (5./3.) //gamma of ions
#define NSTEPSTOP 50000.
#define NOUTSTOP 5
#define DTOUT1 (100)  //dt for output in seconds
#define RHO_INIT 1.0
#define TE_INIT .1
#define TI_INIT 10.
#define TEINITFACTOR 1.
#define VEL 1.e-3
#define FORCEGAMMAGASFIXED
#endif


#if (TESTNO==20) //to test electron heating 
#define CONSISTENTGAMMA
#define FIXEDGAMMASPECIES
//#define DISSIPATIONFROMGASONLY 
#define GAMMA (5./3.)
#define GAMMAE (4./3.) //gamma of electrons
#define GAMMAI (5./3.) //gamma of ions
#define NSTEPSTOP 50000.
#define NOUTSTOP 5
#define DTOUT1 (100)  //dt for output in seconds
#define RHO_INIT 1.
#define TE_INIT .1
#define TI_INIT 10.
#define TEINITFACTOR 1.
#define VEL 1.e-3
#endif


#if (TESTNO==100) //to test electron heating 
#define FORCEGAMMAGASFIXED
//#define CONSISTENTGAMMA
//#define FIXEDGAMMASPECIES

#define DISSIPATIONFROMGASONLY
#define HEATELECTRONSATENDRK2

#define GAMMA (5./3.)
//#define GAMMA 1.5
#define GAMMAE (4./3.) //gamma of electrons
#define GAMMAI (5./3.) //gamma of ions
#define NSTEPSTOP 50000.
#define NOUTSTOP 5
#define DTOUT1 (100.)  //dt for output in seconds
#define RHO_INIT 1.
#define TE_INIT .1
#define TI_INIT 10.
#define TEINITFACTOR 1.
#define VEL 1.e-3
#endif

/************************************/
//coordinates / resolution
/************************************/
#define MYCOORDS MINKCOORDS
#define PHIWEDGE (M_PI/2.)

#define MINX 0.
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

//total resolution
#define TNX 2000//128
#define TNY 1
#define TNZ 1

//number of tiles
#define NTX 1
#define NTY 1
#define NTZ 1

//#define SPECIFIC_BC
#define PERIODIC_ZBC
#define PERIODIC_XBC
#define PERIODIC_YBC

/************************************/
//output
/************************************/
#define OUTCOORDS MINKCOORDS
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define SILOOUTPUT 0
#define OUTOUTPUT 0
#define SIMOUTPUT 0
//#define SIMOUTPUTINTERNAL
#define RADOUTPUT 0
#define SCAOUTPUT 0
#define AVGOUTPUT 0
#define DTOUT2 DTOUT1

/************************************/
//common physics
/************************************/
//#define NOTALLOWFORNEGATIVEHEATING
#define HFRAC 1. //mass fraction of the hydrogen X
#define HEFRAC 0. //mass fraction of helium Y
