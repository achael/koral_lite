/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE
//#define PRECOMPUTE_MY2OUT

/************************************/
//restart
/************************************/
//#define NORESTART
//#define RESTARTFROMNORELEL
#define RESTARTGENERALINDICES
#define RESTARTNUM -1

/************************************/
//blackhole
/************************************/
#define MASS 6.2e9
#define BHSPIN 0.9375

/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1 //must be 1 for MAGNFIELD
#define VECPOTGIVEN
#define INIT_MAGN_CORNERS    
#define MAXBETA .01 //target pmag/pgas int the midplane

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
//radiation choices
/************************************/
#define RADIATION
//#define SKIPRADSOURCE    //advective only
//#define SKIPRADEVOLUTION //keeps initial values

#ifdef RADIATION
//#define BASICRADIMPLICIT
//#define RADIMPSTOPWHENFAIL
//#define RADIMP_START_WITH_BISECT
//#define BALANCEENTROPYWITHRADIATION
#define ALLOWRADCEILINGINIMPLICIT
#define EVOLVEPHOTONNUMBER
#define OPDAMPINIMPLICIT 1
#define SCALE_JACOBIAN
#endif

#define U2PCONV 1.e-10
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
#define RADIMPLICITMAXENCHANGEDOWN 100.
#define RADIMPLICITMAXENCHANGEUP 10.
#define RADIMPLICITMAXTECHANGE 5.
#define MAXRADIMPDAMPING 1.e-6

//opacities
#define SCATTERING
#define BREMSSTRAHLUNG
#define SYNCHROTRON
#define NO_SYNCHROTRON_BRIDGE_FUNCTIONS
#define MAXDIFFTRADS 1.e3
#define MAXDIFFTRADSNEARBH 1.e2
//#define KLEINNISHINA


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
#define MAXRADVISCVEL 0.1
#endif

/************************************/
//electron choices
/************************************/
//#define EVOLVEELECTRONS

#ifdef EVOLVEELECTRONS
#define CONSISTENTGAMMA
#define GAMMAINTCONSISTENTWITHCV //Ramesh's routine for inverting gammaint
#define RELELENTROPY

//#define HEATELECTRONSATENDRK2
//#define DISSIPATIONFROMGASONLY
//#define FORCEGAMMAGASFIXED
#define MIXENTROPIESPROPERLY
#define UPWINDENTROPYMIXING
#define DONOTLIMITENTRINMIXING

#define PRINTVISCHEATINGTOSILO
#define PRINTCOULOMBTOSILO
#define CALCVISCHEATING
#define HEATELECTRONS
#define EXPLICITELECTRONHEATING

#define HEATELECTRONS_HOWES
//#define HEATELECTRONS_ROWAN
//#define HEATELECTRONS_ROWAN2

#define UEUINTMINRATIO 1.e-3
#define UIUINTMINRATIO 1.e-3
#define TEMPEMINIMAL 1.e2
#define TEMPIMINIMAL 1.e2
#define TEMPEMINIMALFRACTION 1.e-6
#define TEMPIMINIMALFRACTION 1.e-6
#define TEMPEMAXIMALFRACTION 1.e3
#define TEMPIMAXIMALFRACTION 1.e3
#endif

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 2
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .5
#define FLUXMETHOD LAXF_FLUX

#define FLUXLIMITER 0
#define MINMOD_THETA 1.9
#define SHUFFLELOOPS 0

#define DOFIXUPS 1
#define DOU2PRADFIXUPS 0 //ANDREW -- investigate why this leads to failures!
#define DOU2PMHDFIXUPS 1
#define DORADIMPFIXUPS 1

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-8
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
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
//resolution
/************************************/
//total resolution
#define TNX 64//320//288 //256 //16*16 //128 <->256
#define TNY 64//256//224 //320 //14*9 //92 <->256
#define TNZ 1//128//96

//number of tiles
#define NTX 32
#define NTY 16//32
#define NTZ 8

/************************************/
//coordinates
/************************************/
//#define myJETCOORDS
//#define myMKS2COORDS
#define myMKS3COORDS
#define METRICAXISYMMETRIC
#define RMIN 1.
#define RMAX 1.e5

#ifdef myMKS2COORDS //modified Kerr-Shild
#define METRICNUMERIC
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)

#define MKSR0 -1.35
#define MKSH0 0.7
#define MKSMY1 0.002
#define MKSMY2 0.02
#define MKSMP0 1.3
#endif

#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY 0.
#define MAXY 1.

#define MKSR0 -1.35
#define MKSH0 0.7
#define MKSMY1 0.002
#define MKSMY2 0.02
#define MKSMP0 1.3
#endif

#ifdef myJETCOORDS //concentrate resolution in jet and disk zones
#define MYCOORDS JETCOORDS

#define METRICNUMERIC
#define DERIVS_NOGSL // use a faster numeric derivative for coordinate transformations
#define MINX 0
#define MAXX 1.
#define MINY -(1.-1.e-8) //10^-8 away from the poles seems to be the last safe point
#define MAXY 1.-1.e-8    //TODO -- can we fix this!??!?!?

#define MKSR0 -1.35
#define HYPRBRK 5000
#define FJET 0.4
#define FDISK 0.5
#define RUNI RMIN
#define RCOLL_JET 1000
#define RDECOLL_JET 2*RMIN
#define RCOLL_DISK 5*RMIN
#define RDECOLL_DISK 2*RMIN
#define ALPHA_1 1
#define ALPHA_2 0.375

#define CYLINDRIFY
#define RCYL 20
#define NCYL 1
#endif

#define PHIWEDGE (2*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

#define COPY_XBC     //simpler outflowing boundary conditions
//#define COPY_YBC
#define SPECIFIC_BC
#define PERIODIC_ZBC
//#define PERIODIC_XBC
//#define PERIODIC_YBC

/************************************/
//output
/************************************/
#define DTAVG .1 //how frequently to compute quantities in avg
#define DTOUT1 10 //res
#define DTOUT2 100 //avg
#define DTOUT3 100 //box,var

//#define DUMPS_READ_HDF5
//#define DUMPS_WRITE_HDF5

#define OUTCOORDS BLCOORDS
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e20
#define NOUTSTOP 1000

#define COORDOUTPUT 1
#define GRIDOUTPUT 1
#define SILOOUTPUT 1
#define CGSOUTPUT
#define OUTOUTPUT 0
#define SIMOUTPUT 0//2
//#define SIMOUTPUT_GILAB2FF
#define GRTRANSSIMOUTPUT

#define RADOUTPUT 0
#define SCAOUTPUT 0
#define AVGOUTPUT 1
#define NORELELAVGS
#define THOUTPUT 0
#define THPROFRADIUS 30
#define SILO2D_XZPLANE

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)
#define GAMMAI (5./3.)
#define GAMMAE (4./3.)
#define NTORUS 0
#define HFRAC 1. //mass fraction of the hydrogen X
#define HEFRAC 0. //mass fraction of helium Y


#define RHOATMMIN  1.e-30
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN,0,0,0))
#define ATMTRADINIT 3 //1.e-2
#define ERADATMMIN  calc_LTE_EfromT(ATMTRADINIT) //((calc_LTE_EfromT(3.e6)/10)/1.e19)

#if(NTORUS==0) //a=0 SANE, no rad, denser loops
#define RHOMAX_VPOT_CGS 1.88e-18//4.725e-17
#define RMAX_VPOT 23
#define RHO_VPOT_CUT .1
#define EXPECTEDHR 0.4
#define LT_KAPPA 3.e8//3.e10//6.e11
#define LT_XI 0.7135
#define LT_R1 42.
#define LT_R2 800
#define LT_GAMMA 5./3.
#define LT_RIN 10.5
#undef MAXBETA
#define MAXBETA (1./100.) //target pmag/pgas inside torus
#define BETANORMFULL
#endif

#if(NTORUS==14 || NTORUS==10) //a=0 SANE, no rad, fewer loops
#define EXPECTEDHR 0.4
#define LT_KAPPA 1.e16
#define LT_XI 0.708
#define LT_R1 42.
#define LT_R2 1000.
#define LT_GAMMA 5./3.
#define LT_RIN 10.
#undef MAXBETA
#define MAXBETA (1./100.) //target pmag/pgas inside torus
#define BETANORMFULL
#endif

