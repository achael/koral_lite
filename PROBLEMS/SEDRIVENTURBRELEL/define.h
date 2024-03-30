#pragma once

#define FFT
#define SRANDSEED 2142352

/************************************/
//general
/************************************/
//#define BHDISK_PROBLEMTYPE
//#define ENFORCEENTROPY

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM 80
#define NOUTSTOP 500

/************************************/
//radiation choices
/************************************/
//#define RADIATION
//#define SKIPRADSOURCE
#define RADIMPSTOPWHENFAIL
#define MAXDIFFTRADS 1.e6
#define U2PCONV 1.e-12
#define RADIMPCONV 1.e-10
#define RADIMPEPS 1.e-6
#define RADIMPMAXITER 50

//no opacities!
//#define NCOMPTONIZATION 
//#define BREMSSTRAHLUNG
//#define SYNCHROTRON
//#define SYNCHROTRONALLOWEDTRATIO 1.e5

//#define BOUNDFREE
//#define KLEINNISHINA
//#define BOUNDELECTRON
//#define COLDOPACITIES
//#define SKIPFANCYOPACITIES

#define EVOLVEELECTRONS
#define CONSISTENTGAMMA
//#define ELECTRONHEATTYPE ELECTRONHEATTYPE_THROUGHUINT
#define GAMMAINTCONSISTENTWITHCV //Ramesh's routine for inverting gammaint
#define RELELENTROPY

#define MIXENTROPIESPROPERLY //ANDREW Entropy mixing
#define UPWINDENTROPYMIXING
#define DONOTLIMITENTRINMIXING

//#define HEATELECTRONSATENDRK2
#define ACCUMULATEVISCHEATING
//#define SKIPCOULOMBCOUPLING
#define HEATELECTRONS
//#define DISSIPATIONFROMGASONLY
//#define DISSIPATIONFROMGASONLYU2P
//#define HEATELECTRONS_USEFIT
#define HEATELECTRONS_DELTA .1

/************************************/
//Relativistic electron choices
/************************************/
#define RELELSPECTRUMOUTPUT 1
#define NTH_SPEC_IX NTX/2
#define NTH_SPEC_IY NTY/2
#define NTH_SPEC_IZ 0

#define SCALE_JACOBIAN

#define RELELECTRONS
//#define RELEL_INTEGRATE_SIMPSON
#define NRELBIN 65
#define NSCALARS (NRELBIN + 12)
#define RELGAMMAMIN 500.
#define RELGAMMAMAX 5.e6
//#define SKIPRELELEXPANSION
#define RELEL_ADIAB_DT
#define RELEL_ADIAB_LOGSPACE

#define NO_RELEL_ADIAB_SCALE
//#define RELEL_ADIAB_LOGSPACE_2NDORDER
//#define RELEL_ADIAB_PLAW_INTERP
//#define RELEL_ADIAB_RK2

//#define RELEL_ADIAB_ART
//#define RELEL_SYN
//#define RELEL_SYN_ART
//#define RELEL_FF
//#define RELEL_FF_ART
//#define RELEL_IC
//#define RELEL_IC_ART
//#define RELEL_COUL
//#define RELEL_COUL_ART

//#define RELEL_B_ART 5.0e2
//#define RELEL_NION_ART 3.7e15 //cm^-3
//#define RELEL_NETH_ART 3.7e15 //cm^-3
//#define RELEL_TRAD_ART 3e4 //K
//#define RELEL_ERAD_ART 10e8 //erg cm^-3

//#define SKIPRELELHEATING
//#define NORELELNEGHEAT
//#define RELEL_HEAT_ART
//#define RELEL_HEAT_NORM .001
#define RELEL_HEAT_FIX_FRAC //ANDREW
#define RELEL_HEAT_FRAC 0.05
#define RELEL_HEAT_INDEX 3.5
#define RELEL_INJ_MIN 499.9
#define RELEL_INJ_MAX 5.001e6
//#define RELEL_CONST_INJPARAMS
#define RELEL_HEAT_FIX_LIMITS

#define MAX_RELEL_FRAC_N 0.9
#define MAX_RELEL_FRAC_U 0.9
#define MAX_RELEL_FRAC_P 0.95
/************************************/
//magnetic choices
/************************************/
//#define MAGNFIELD
#define GDETIN 1

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2HEUN
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
#define MASS (1./MSUNCM)
#define BHSPIN 0.


/************************************/
//coordinates / resolution
/************************************/
#define MYCOORDS MINKCOORDS
#define MINX 0.
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

//total resolution
#define TNX 128//512 //28*9
#define TNY 128//512 //26*9
#define TNZ 1//2*8
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
#define TESTCOMPTINSCALARS4FLAT
//#define TIMEINSCALARSINSEC
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e12

#define SILOOUTPUT 1
#define PRINTVISCHEATINGTOSILO

#define OUTOUTPUT 0
#define SIMOUTPUT 0
#define RADOUTPUT 0
#define SCAOUTPUT 1
#define AVGOUTPUT 0
#define DTOUT1 10//100.//100.

#define DTOUT2 1000.

/************************************/
//common physics
/************************************/
#define GAMMA (5./3.)
#define GAMMAE (5./3.) //gamma of electrons
#define GAMMAI (5./3.) //gamma of ions

//#define TEINIT 5.e9 //if not defined - start from CS0_INIT 
//#define TIINIT 5.e9
//default 5.e9

#ifdef TEINIT
#undef DTOUT1
#define DTOUT1 (.1*sqrt(TEINIT/5.e9))
#endif
//#define FORCEGAMMAGASFIXED
//#define DISSIPATIONFROMGASONLY
#define HFRAC 1. //mass fraction of the hydrogen X
#define HEFRAC 0. //mass fraction of helium 

/************************************/
//initial setup
/************************************/
#define RHO_INIT 1.
#define CS0_INIT 8.6e-4
#define VX_INIT 0.0
#define VY_INIT 0.0
#define VZ_INIT 0.0
#define TR_INIT tempCGS2GU(1.e7)
//#define EE_INIT RHO_INIT
#define NPH_INIT numdensCGS2GU(2.0287e3)

#define RELEL_NRELEL_INIT 0.0 //ANDREW

//these don't work with implicit!
/*
#define RHO_INIT 1.
#define CS0_INIT 8.6e-4
#define VX_INIT 0.0
#define VY_INIT 0.0
#define VZ_INIT 0.
//#define EE_INIT RHO_INIT
#define TR_INIT tempCGS2GU(1.e3)
#define NPH_INIT numdensCGS2GU(2.0287e3)
*/

#define BX_INIT (RHO_INIT/1.e2)
#define BY_INIT (RHO_INIT/1.e2)
#define BZ_INIT 0.
#define BETA_INIT 6.

//#define HUBBLEBACKGROUND
//#define HUBBLEMAGN 1.e-2

#define PERT_SPECTRUM
#define PERT_SPECTRUM_MODE 2.
#define DELTAV_NORM 3.e-7
//#define DELTAV_NORM 0.0e-8
//3.e-8 - fiducial from the paper test
//1.e-5 - fiducial for these test
//1.e-4 - strong
//1.e-6 - weak
//#define PERT_RANDOM
//#define PERT_INIT
//#define VPERT 0.01
//#define VERTCOLLISION
//#define DBLVERTCOLLISION
