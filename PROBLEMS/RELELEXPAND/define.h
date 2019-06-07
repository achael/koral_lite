/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
//#define RADIATION
#define RADIMPSTOPWHENFAIL
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

#define MAXDIFFTRADS 1.e6
#define U2PCONV 1.e-12
#define RADIMPCONV 1.e-10
#define RADIMPEPS 1.e-6

//#define NCOMPTONIZATION
#define SCALE_JACOBIAN
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
#define RELELECTRONS

//#define RELEL_SYN
//#define RELEL_SYN_ART
//#define RELEL_FF
//#define RELEL_FF_ART
//#define RELEL_IC
//#define RELEL_IC_ART
//#define RELEL_COUL
//#define RELEL_COUL_ART

#define HEATELECTRONS
//#define HEATELECTRONS_USEFIT
//#define HEATELECTRONS_DELTA 0.01
//#define TEMPEDAMPFRACTION 0.99
#define UEUINTMINRATIO 1e-2
#define UIUINTMINRATIO 1e-2
#define TEMPEMINIMAL 1e2
#define TEMPIMINIMAL 1e2
#define TEMPEMINIMALFRACTION 1e-6
#define TEMPIMINIMALFRACTION 1e-6
#define TEMPEMAXIMALFRACTION 1e3
#define TEMPIMAXIMALFRACTION 1e3

/************************************/
//magnetic choices
/************************************/
//#define MAGNFIELD
#define GDETIN 1

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2HEUN//IMEX
#define TSTEPLIM .5
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 1
#define DORADFIXUPS 1

/************************************/
//viscosity choices
/************************************/

/************************************/
//rmhd floors
/************************************/
//#define CORRECT_POLARAXIS
//#define POLARAXISAVGIN3D
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-10
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
//#define MASS 4.e6 //Sgr A*
#define BHSPIN 0.


/************************************/
//initial setup
/************************************/
#define TEST 100

#if (TEST==100) // pure Coulomb
#define CONSISTENTGAMMA
#define RELELENTROPY
#define CONSISTENTGAMMAINT
#define NSTEPSTOP 1.e20
#define SKIPRADSOURCE
#define SIZE 1.e-3 //scales time step

#define DTOUT_LOG  1
#define DTOUT1_LOG_INIT -2 

#define TMAX (10000)    //stop in seconds, if takes too long, increase  SIZE
#define DTOUT1 (10000)  //dt for output in seconds

#define NTH_SPEC_IX 0
#define NTH_SPEC_IY 0
#define NTH_SPEC_IZ 0

#define RHO_INIT rhoCGS2GU(10.0)
#define NPH_INIT numdensCGS2GU(2.0287e25)
#define TE_INIT tempCGS2GU(1.e11)
#define TI_INIT tempCGS2GU(1.e11)
#define TR_INIT tempCGS2GU(1.e10)


//For relativistic electrons
#define NRELBIN 128
#define NSCALARS (NRELBIN + 12)
#define RELGAMMAMIN 1.e2
#define RELGAMMAMAX 1.e12

#define RELEL_ADIAB_ART
//#define SKIPRELELEXPANSION
#define RELEL_DIV_ART -0.005 

#define RELEL_INTEGRATE_LOGSPACE

#define RELEL_ADIAB_LOGSPACE_LF
//#define RELEL_ADIAB_LOGSPACE
//#define RELEL_ADIAB_RK2
#define RELEL_MINMOD_THETA 1.25

#define RELEL_INIT_NORM 1000.
//#define RELEL_INIT_MAXWELL
//#define RELEL_INIT_THETA 100. 
//

#define RELEL_NRELEL 0.0
//#define RELEL_INIT_PLAW
#define RELEL_INIT_INDEX 3.5
//#define RELEL_INIT_MIN 49
//#define RELEL_INIT_MAX 500001

#define RELEL_HEAT_ART
//#define SKIPRELELHEATING
#define RELEL_HEAT_NORM 1000. 
#define RELEL_HEAT_FIX_INDEX
#define RELEL_HEAT_INDEX 3.5
#define RELEL_HEAT_FIX_LIMITS
#define RELEL_INJ_MIN 5.e2
#define RELEL_INJ_MAX 5.e5

#endif

/************************************/
//coordinates / resolution
/************************************/
#define MYCOORDS MINKCOORDS

#define MINX 0.
#define MAXX SIZE
#define MINY 0.
#define MAXY SIZE
#define MINZ 0.
#define MAXZ SIZE

//total resolution
#define TNX 4//256 //28*9
#define TNY 1//256 //26*9
#define TNZ 1 //2*8
//number of tiles
#define NTX 16
#define NTY 1
#define NTZ 1

#define PERIODIC_ZBC
#define PERIODIC_XBC
#define PERIODIC_YBC

/************************************/
//output
/************************************/
#define OUTCOORDS MINKCOORDS
#define TIMEINSCALARSINSEC
#define ALLSTEPSOUTPUT 0
#define NOUTSTOP 5000
#define SILOOUTPUT 0
#define OUTOUTPUT 0
#define SIMOUTPUT 0
#define RADOUTPUT 0
#define SCAOUTPUT 1
#define RELELSPECTRUMOUTPUT 1
#define AVGOUTPUT 0
//#define DTOUT2 1.e100

/************************************/
//common physics
/************************************/
#define GAMMA (5./3.)
#define GAMMAE (5./3.) //gamma of electrons
#define GAMMAI (5./3.) //gamma of ions
#define HFRAC 1. //mass fraction of the hydrogen X
#define HEFRAC 0. //mass fraction of helium Y

