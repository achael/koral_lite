/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
#define RADIATION
#define RADIMPSTOPWHENFAIL
#define RADIMPLICITDAMPINGFACTOR 3.
#define RADIMPLICITMAXENCHANGEDOWN 1.e5
#define RADIMPLICITMAXENCHANGEUP 1.e5
#define MAXRADIMPDAMPING 1.e-3
#define RADIMPMAXITER 70
#define RADIMPCOVREL 1.e-6
#define RADIMPCONVRELERR .99
#define RADIMPCONVRELENTR 1.e-4
#define RADIMPCONVRELENTRERR 1.e-2
#define OPDAMPINIMPLICIT 0

#define MAXDIFFTRADS 1.e6
#define U2PCONV 1.e-12
#define RADIMPCONV 1.e-10
#define RADIMPEPS 1.e-6

//#define EVOLVEPHOTONNUMBER
#define SCALE_JACOBIAN
//#define BREMSSTRAHLUNG
//#define SYNCHROTRON
//#define SYNCHROTRON_OLD_FIT
//#define SYNCHROTRONALLOWEDTRATIO 1.e5

//#define BOUNDFREE
//#define KLEINNISHINA
//#define BOUNDELECTRON
//#define COLDOPACITIES
//#define SKIPFANCYOPACITIES
#define SKIPCOULOMBCOUPLING
#define NO_COMPTONIZATION

/************************************/
//electron choices
/************************************/

#define EVOLVEELECTRONS
#define RELELECTRONS
#define HEATELECTRONS
//#define HEATELECTRONS_USEFIT
#define HEATELECTRONS_DELTA 1.
//#define TEMPEDAMPFRACTION 0.99
#define UEUINTMINRATIO 1e-2
#define UIUINTMINRATIO 1e-2
#define TEMPEMINIMAL 1e2
#define TEMPIMINIMAL 1e2
#define TEMPEMINIMALFRACTION 1e-6
#define TEMPIMINIMALFRACTION 1e-6
#define TEMPEMAXIMALFRACTION 1e2
#define TEMPIMAXIMALFRACTION 1e2
//#define NOTALLOWFORNEGATIVEHEATING

#define MAX_RELEL_FRAC_N 0.999
#define MAX_RELEL_FRAC_U 0.999
#define MAX_RELEL_FRAC_P 0.999
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
#define TSTEPLIM .95
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 1
#define DORADFIXUPS 1

/************************************/
//thermal electrons
/************************************/
#define CONSISTENTGAMMA
#define GAMMAINTCONSISTENTWITHCV //Ramesh's routine for inverting gammaint
//#define ELECTRONIONHEATTYPE ELECTRONIONHEATTYPE_THROUGHENTROPY
#define RELELENTROPY
//#define SKIPRADSOURCE
//#define ENFORCE_HEATING_ENERGY_SUM

/************************************/
//rmhd floors
/************************************/
//#define CORRECT_POLARAXIS
//#define POLARAXISAVGIN3D
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 0.//1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 0.//1.e-20
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

#define BHSPIN 0.

/************************************/
// test problem setup
/************************************/
// common setup
#define RHO_INIT rhoCGS2GU(1.e-8)
#define NPH_INIT numdensCGS2GU(2.0287e12)
#define TE_INIT tempCGS2GU(1.e12) 
#define TI_INIT tempCGS2GU(1.e12) 
#define TR_INIT tempCGS2GU(1.e4)

//#define TEST 100 // synchrotron only
//#define TEST 200 // sychrotron + IC
//#define TEST 300 // Wong+Zhdankin turbulent diffusion-advection
#define TEST 400

// Synchrotron Test
#if (TEST == 100)
#define MASS (1./MSUNCM*CCC_CGS) //so that time in seconds!

#define SIZE 1.e-2 //scales time step
#define NSTEPSTOP 100000000

#define DTOUT_LOG  1        //logarithmic output in time
#define DTOUT1_LOG_INIT -3  //begin at 1.e-3 s, then 1.e-2, etc

#define TMAX  10000 //3 //stop in seconds, if takes too long, increase  SIZE
#define DTOUT1 .01 //.001  //dt for output in seconds

#define NRELBIN 100//64
#define NSCALARS (NRELBIN + 12)
#define RELGAMMAMIN 1.001
#define RELGAMMAMAX 1.e10

#define RELEL_IMPLICIT_LOGSPACE_LF
#define RELEL_MINMOD_THETA 1.5
#define RELEL_NRELEL 0.0000001 //initial condition
#define SKIPRELELEXPANSION
//#define SKIPRELELHEATING

#define RELEL_SYN
#define RELEL_SYN_ART
#define RELEL_B_ART 200.
//#define RELEL_NION_ART 3.7e15 //cm^-3
//#define RELEL_NETH_ART 3.7e15 //cm^-3
//#define RELEL_TRAD_ART 3e4 //K
//#define RELEL_ERAD_ART 10e8 //erg cm^-3
//#define RELEL_ADIAB_ART
//#define RELEL_DIV_ART 0. 

#define RELEL_HEAT_FIX_FRAC
#define RELEL_HEAT_FRAC 0.2
#define RELEL_HEAT_ART
#define RELEL_HEAT_NORM 1 //cm^-3  
#define RELEL_HEAT_FIX_INDEX
#define RELEL_HEAT_INDEX 3.5
#define RELEL_HEAT_FIX_LIMITS

#define RELEL_INJ_MIN 1000.
#define RELEL_INJ_MAX 1.e6
#endif

//Synchrotron + IC -- linear time output
#if (TEST == 200)
#define MASS (1./MSUNCM*CCC_CGS)*3155692597470 //so that time in 10^5 yr!
#define SIZE 1.e-3 //scales time step
#define NSTEPSTOP 100000000
#define TMAX  100. //stop in 10^5 yr
#define DTOUT1 1. //dt for output in 10^5 yr

#define NRELBIN 100
#define NSCALARS (NRELBIN + 12)
#define RELGAMMAMIN 1.001
#define RELGAMMAMAX 1.e10

#define RELEL_IMPLICIT_LOGSPACE_LF
#define RELEL_MINMOD_THETA 1.5
#define RELEL_NRELEL 0.0000001 //initial condition
#define SKIPRELELEXPANSION
//#define SKIPRELELHEATING

#define RELEL_SYN
#define RELEL_SYN_ART
#define RELEL_IC
#define RELEL_IC_ART
//#define RELEL_ADIAB_ART
//#define RELEL_DIV_ART 0. 

#define RELEL_B_ART 10.e-6
#define RELEL_TRAD_ART 30000. //K 
#define RELEL_ERAD_ART 7.95775e-10 //erg cm^-3 

#define RELEL_HEAT_ART
#define RELEL_HEAT_NORM 0.001 //0.001 //cm^-3
#define RELEL_HEAT_INDEX 2.0000001
#define RELEL_INJ_MIN 1.e2
#define RELEL_INJ_MAX 1.e9
#define RELEL_HEAT_FRAC 0.2
#define RELEL_HEAT_FIX_LIMITS
#define RELEL_HEAT_FIX_FRAC
#define RELEL_HEAT_FIX_INDEX
#endif

// Diffusion-Advection Test
#if (TEST == 300)
#define MASS (1.) // TODO -- figure out time scaling later

#define SIZE 1.e-4 //scales time step
#define NSTEPSTOP 100000000
#define TMAX  10 
#define DTOUT1 1 

#define NRELBIN 100//64
#define NSCALARS (NRELBIN + 12)
#define RELGAMMAMIN 1.001
#define RELGAMMAMAX 1.e10

#define RELEL_IMPLICIT_LOGSPACE_LF
#define RELEL_MINMOD_THETA 1.5
#define RELEL_NRELEL 0.0000001 //initial condition
#define RELEL_INIT_MAXWELL
#define RELEL_INIT_THETA 100.
#define RELEL_INIT_NORM 1.

#define SKIPRELELEXPANSION
#define SKIPRELELHEATING

#define RELEL_TURB_ADVECT // TODO -- figure out parameter dependence
#define RELEL_DIFFUSE
#define RELEL_TURB_DIFFUSE

#define RELEL_HEAT_INDEX 3.5
#endif

// Diffusion-Advection Test, steady-state
#if (TEST == 400)
#define MASS (1.) // TODO -- figure out time scaling later

#define SIZE 1.e-6 //scales time step
#define NSTEPSTOP 100000000
#define TMAX  .3
#define DTOUT1 .1 

#define NRELBIN 64
#define NSCALARS (NRELBIN + 12)
#define RELGAMMAMIN 1.001
#define RELGAMMAMAX 1.e6

#define RELEL_IMPLICIT_LOGSPACE_LF
#define RELEL_MINMOD_THETA 1.5
#define RELEL_NRELEL 0.00000000001 //initial condition
#define RELEL_INIT_MAXWELL
#define RELEL_INIT_THETA 100.
#define RELEL_INIT_NORM 1.

#define SKIPRELELEXPANSION
#define SKIPRELELHEATING

#define RELEL_ADVECT_Z19
#define RELEL_GAMMA0 300. 
#define RELEL_RATEA -3.5
#define RELEL_RATEH -5.0

#define RELEL_DIFFUSE
#define RELEL_DIFFUSE_Z19
#define RELEL_RATE2 1.4
#define RELEL_TAUC 1.0

//#define RELEL_TURB_ADVECT // TODO -- figure out parameter dependence
//#define RELEL_TURB_DIFFUSE
#define RELEL_HEAT_INDEX 3.5

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
#define NTX 4
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
#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define SIMOUTPUT 0 //2
#define OUTPUTINGU
#define RADOUTPUT 1
#define SCAOUTPUT 1
#define RELELSPECTRUMOUTPUT 1
#define AVGOUTPUT 0
#define DTOUT2 1.e100

#define NTH_SPEC_IX 0 // which bin to print output from
#define NTH_SPEC_IY 0
#define NTH_SPEC_IZ 0

/************************************/
//common physics
/************************************/
#define GAMMA (5./3.)
#define GAMMAE (5./3.) //gamma of electrons
#define GAMMAI (5./3.) //gamma of ions
#define HFRAC 1. //mass fraction of the hydrogen X
#define HEFRAC 0. //mass fraction of helium Y

