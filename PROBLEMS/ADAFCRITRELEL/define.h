/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
//#define NORESTART
//#define RESTARTFROMNORELEL
#define RESTARTGENERALINDICES
#define RESTARTNUM -1

//#define RESCALEDENSITY 3

/************************************/
//radiation choices
/************************************/
//#define RADIATION
//#define SKIPRADSOURCE    //advective only
//#define SKIPRADEVOLUTION //keeps initial values

#ifdef RADIATION
//#define BASICRADIMPLICIT
//#define RADIMPSTOPWHENFAIL
//#define RADIMP_START_WITH_BISECT
//#define BALANCEENTROPYWITHRADIATION
#define EVOLVEPHOTONNUMBER
#define OPDAMPINIMPLICIT 0
#define SCALE_JACOBIAN
#endif

#define U2PCONV 1.e-12
#define CHECKENTROPYAFTEREXPLICIT 
#define CHECKENTROPYAFTEREXPLICITFACTOR 1.// - accept everything
#define RADIMPCONVREL 1.e-6
#define RADIMPCONVRELERR (1.e-4)
#define RADIMPCONVRELENTR 1.e-6
#define RADIMPCONVRELENTRERR 1.e-4
#define RADIMPCONV 1.e-8
#define RADIMPENTRCONV 1.e-5
#define RADIMPEPS 1.e-8
#define RADIMPMAXITER 100
#define RADIMPLICITDAMPINGFACTOR 5.
#define RADIMPLICITMAXNPHCHANGE 100.
#define RADIMPLICITMAXENCHANGEDOWN 10.
#define RADIMPLICITMAXENCHANGEUP 10.
#define RADIMPLICITMAXTECHANGE 5.
#define MAXRADIMPDAMPING 1.e-6

//opacities
#define SCATTERING
#define BREMSSTRAHLUNG
#define SYNCHROTRON
#define MAXDIFFTRADS 1.e4
#define MAXDIFFTRADSNEARBH 1.e2

//#define BOUNDFREE
//#define KLEINNISHINA
//#define BOUNDELECTRON
//#define COLDOPACITIES
//#define SKIPFANCYOPACITIES
//#define SKIPCOULOMBCOUPLING

/************************************/
//electron choices
/************************************/
/*
//#define EVOLVEELECTRONS

#define CONSISTENTGAMMA
#define GAMMAINTCONSISTENTWITHCV //Ramesh's routine for inverting gammaint
#define RELELENTROPY

//#define HEATELECTRONSATENDRK2
//#define DISSIPATIONFROMGASONLY
//#define FORCEGAMMAGASFIXED
#define MIXENTROPIESPROPERLY
#define UPWINDENTROPYMIXING
#define DONOTLIMITENTRINMIXING
//#define MIXENTROPIES_CONST_PRESSURE

#define PRINTVISCHEATINGTOSILO
#define PRINTCOULOMBTOSILO
#define CALCVISCHEATING
#define HEATELECTRONS
#define EXPLICITELECTRONHEATING

#define HEATELECTRONS_HOWES
//#define HEATELECTRONS_ROWAN2

#define UEUINTMINRATIO 1.e-3
#define UIUINTMINRATIO 1.e-3
#define TEMPEMINIMAL 1.e2
#define TEMPIMINIMAL 1.e2
#define TEMPEMINIMALFRACTION 1.e-6
#define TEMPIMINIMALFRACTION 1.e-6
#define TEMPEMAXIMALFRACTION 1.e3
#define TEMPIMAXIMALFRACTION 1.e3

*/
/************************************/
//Relativistic electron choices
/************************************/
//#define RELELSPECTRUMOUTPUT 2
//#define RELELECTRONS

//#define NRELBIN 48
//#define NSCALARS (NRELBIN + 12)
//#define RELGAMMAMIN 100.
//#define RELGAMMAMAX 1.e8

//#define SKIPRELELEXPANSION
//#define SKIPRELELINJECTION

//#define RELEL_ADIAB_LOGSPACE_LF
//#define RELEL_IMPLICIT_LOGSPACE_LF
//#define RELEL_MINMOD_THETA 1.8
//#define ZERO_NONTHERMAL_LOWGAMMA

//#define RELEL_SYN
//#define RELEL_SYN_ART
//#define RELEL_B_ART 2.0e3

//#define RELEL_FF
//#define RELEL_FF_ART
//#define RELEL_NION_ART 3.7e15 //cm^-3

//#define RELEL_IC
//#define RELEL_IC_ART
//#define RELEL_NETH_ART 3.7e15 //cm^-3

//#define RELEL_COUL
//#define RELEL_COUL_ART
//#define RELEL_TRAD_ART 3e4 //K
//#define RELEL_ERAD_ART 10e8 //erg cm^-3

//#define NORELELHEATATBH

//#define RELEL_HEAT_FIX_FRAC
//#define RELEL_HEAT_FIX_INDEX
//#define RELEL_HEAT_FRAC 0.015 
//#define RELEL_HEAT_INDEX 3.5

//#define RELEL_HEAT_FIX_LIMITS
//#define RELEL_INJ_MIN 100.
//#define RELEL_INJ_MAX 1.e8

//#define MAX_RELEL_FRAC_N 0.9
//#define MAX_RELEL_FRAC_U 0.9
//#define MAX_RELEL_FRAC_P 0.95

/************************************/
//magnetic choices
/************************************/

//#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN
#define INIT_MAGN_CORNERS //for consistent vec potential    
#define MAXBETA .01 //target pmag/pgas int the midplane
//#define MAXBETA_SEPARATE

/************************************/
//dynamo choices
/************************************/

#define MIMICDYNAMO
#define CALCHRONTHEGO
#define THETAANGLE 0.25
#define ALPHAFLIPSSIGN                                                        

#define ALPHADYNAMO 0.314
#define DAMPBETA
#define BETASATURATED 0.1
#define ALPHABETA 6.28

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .9
#define FLUXMETHOD LAXF_FLUX

#define FLUXLIMITER 0
#define MINMOD_THETA 1.9
#define SHUFFLELOOPS 0

#define DOFIXUPS 1
#define DOU2PRADFIXUPS 0 //ANDREW -- investigate why this leads to failures!
#define DOU2PMHDFIXUPS 1
#define DORADIMPFIXUPS 1

/************************************/
//viscosity choices
/************************************/
#ifdef RADIATION
#define RADVISCOSITY SHEARVISCOSITY
#define ACCELRADVISCOSITY
#define RADVISCMFPSPHMAX 10.
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL 0.05
#endif

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-8
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-30
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 1.e5
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 100.
#define GAMMAMAXRAD 20.
#define GAMMAMAXHD 20.
#define RHOFLOOR 1.e-50

/************************************/
//blackhole
/************************************/
#define MASS 4.e6
#define BHSPIN 0.0

/************************************/
//coordinates / resolution
/************************************/
#define myMKS3COORDS
#define METRICAXISYMMETRIC
#define RMIN 1.5 
#define RMAX 1000.
#define MKSR0 -2.
#define MKSH0 0.7
#define MKSMY1 0.002
#define MKSMY2 0.02
#define MKSMP0 1.3

//#ifdef myMKS2COORDS //modified Kerr-Shild
//#define MYCOORDS MKS2COORDS
//#define MINX (log(RMIN-MKSR0))
//#define MAXX (log(RMAX-MKSR0))
//#define MINY (0.001)
//#define MAXY (1.-0.001)
//#endif

#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY 0.
#define MAXY 1.
#endif

#define PHIWEDGE (2*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

//total resolution
#define TNX 512//256 //16*16 //128 <->256
#define TNY 64//192 //320 //14*9 //92 <->256
#define TNZ 1//96

//number of tiles
#define NTX 32
#define NTY 16
#define NTZ 1//8

#define SPECIFIC_BC
#define PERIODIC_ZBC
//#define PERIODIC_XBC
//#define PERIODIC_YBC

/************************************/
//output
/************************************/
#define DTOUT1 10 //res
#define DTOUT2 1000 //avg
#define DTOUT3 100 //box,var

#define OUTCOORDS BLCOORDS//KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
//#define RADOUTPUTINZAMO
#define NSTEPSTOP 1.e20
#define NOUTSTOP 2000

//#define OUTPUTPERCORE
#define COORDOUTPUT 0
#define SILOOUTPUT 0
#define CGSOUTPUT
#define OUTOUTPUT 0
#define SIMOUTPUT 0
#define SIMOUTPUT_GILAB2FF
//#define SIMOUTPUTWITHINDTHETA .05
#define GRTRANSSIMOUTPUT

#define RADOUTPUT 0
//#define RADOUTPUTWITHINDTHETA .05
#define SCAOUTPUT 0
#define AVGOUTPUT 0
#define NORELELAVGS
#define THOUTPUT 0
#define THPROFRADIUS 30
#define SILO2D_XZPLANE
#define CBAUTOSCALE

#define BOXOUTPUT 0
#define BOXR1 10.
#define BOXR2 15.
#define BOXITH 30 //distance from eq.plane in cells  

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)
#define GAMMAI (5./3.)
#define GAMMAE (4./3.)
#define NTORUS 4
#define HFRAC 1. //mass fraction of the hydrogen X
#define HEFRAC 0. //mass fraction of helium Y

#if(NTORUS==4) //a=0 SANE, no rad, denser loops
#define EXPECTEDHR 0.4
#define LT_KAPPA (6.e11)
#define LT_XI 0.708
#define LT_R1 42.
#define LT_R2 1000.
#define LT_GAMMA 5./3.
#define LT_RIN 10.
#undef MAXBETA
#define MAXBETA (1./100.) //target pmag/pgas inside torus
#define BETANORMFULL
#endif

//#if(NTORUS==14) //a=0 SANE, no rad, fewer loops
//#define EXPECTEDHR 0.4
//#define LT_KAPPA 1.e16
//#define LT_XI 0.708
//#define LT_R1 42.
//#define LT_R2 1000.
//#define LT_GAMMA 5./3.
//#define LT_RIN 10.
//#undef MAXBETA
//#define MAXBETA (1./30.) //target pmag/pgas inside torus
//#define BETANORMFULL
//#endif

#define RHOATMMIN  1.e-27
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN,0,0,0))
#define ATMTRADINIT 3 //1.e-2
#define ERADATMMIN  calc_LTE_EfromT(ATMTRADINIT) //((calc_LTE_EfromT(3.e6)/10)/1.e19)

