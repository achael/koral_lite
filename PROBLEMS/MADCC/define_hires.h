/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
//#define NORESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1
//#define PERTURBUINT 0.02

/************************************/
//blackhole
/************************************/
#define MASS 6.5e9
#define BHSPIN 0.9375

/************************************/
//U2P-related choices
/************************************/
//#define U2P_EQS U2P_EQS_JONS
//#define U2P_SOLVER U2P_SOLVER_WP
#define U2PCONV 1.e-10

/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1    //must be 1 for MAGNFIELD
#define VECPOTGIVEN
#define INIT_MAGN_CORNERS //initialize magnetic field on corners/edges (which?)

/************************************/
//dynamo choices --  2D ONLY
/************************************/
/*
#define MIMICDYNAMO
#define CALCHRONTHEGO
#define THETAANGLE 0.25
#define ALPHAFLIPSSIGN                                                        
#define ALPHADYNAMO 0.314
#define DAMPBETA
#define BETASATURATED 0.1
#define ALPHABETA 6.28
*/
/************************************/
//radiation choices
/************************************/
//#define RADIATION

//#define SKIPRADSOURCE    //advective only
//#define SKIPRADEVOLUTION //keeps initial values in place

#ifdef RADIATION

//#define BASICRADIMPLICIT
//#define RADIMPSTOPWHENFAIL
//#define RADIMP_START_WITH_BISECT
//#define BALANCEENTROPYWITHRADIATION

//#define ALLOWRADCEILINGINIMPLICIT
//#define EVOLVEPHOTONNUMBER
#define OPDAMPINIMPLICIT 0
#define SCALE_JACOBIAN
#endif

//opacities
#define SCATTERING
#define BREMSSTRAHLUNG
//#define KLEINNISHINA
#define SYNCHROTRON
#define NO_SYNCHROTRON_BRIDGE_FUNCTIONS
#define MAXDIFFTRADS 1.e3
#define MAXDIFFTRADSNEARBH 1.e2

//implicit convergence
#define RADIMPCONV 1.e-14
#define RADIMPCONVREL 1.e-6
#define RADIMPCONVRELERR (1.e-4)  
#define RADIMPCONVRELENTR 1.e-6
#define RADIMPCONVRELENTRERR 1.e-4
#define RADIMPENTRCONV 1.e-5
#define RADIMPEPS 1.e-8
#define RADIMPMAXITER 100
#define RADIMPLICITDAMPINGFACTOR 5.
#define RADIMPLICITMAXNPHCHANGE 100.
#define RADIMPLICITMAXENCHANGEDOWN 100.
#define RADIMPLICITMAXENCHANGEUP 10.
#define RADIMPLICITMAXTECHANGE 2.
#define IMPLICITMAXTGASCHANGE 2.
#define MAXRADIMPDAMPING 1.e-7

//Radiation Viscosity choices
#ifdef RADIATION

#define RADVISCOSITY SHEARVISCOSITY
#define ACCELRADVISCOSITY
#define RADVISCMFPSPHMAX 10.
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL 0.1

//#define SKIPHDEVOLUTION
//#define SKIPRADEVOLUTION
//#define SKIPEVOLUTION
//#define SKIPRADSOURCE
//#define SKIPCOULOMBCOUPLING
//#define RADIMPSTOPWHENFAIL
#define DAMPCOMPTONIZATIONATBH
//#define NO_COMPTONIZATION
//#define SKIPFANCYOPACITIES
//#define ENFORCEENTROPY
//#define GASRADCOUPLEDWAVESPEEDS

#endif

/************************************/
//electron choices
/************************************/
//#define EVOLVEELECTRONS

#ifdef EVOLVEELECTRONS
#define CONSISTENTGAMMA
#define GAMMAINTCONSISTENTWITHCV //Ramesh's routine for inverting gamma_int
#define RELELENTROPY

//heating
#define HEATELECTRONS
//#define HEATELECTRONS_HOWES
//#define HEATELECTRONS_ROWAN
#define HEATELECTRONS_ROWAN2
//#define HEATELECTRONS_ROWAN3

#define NOHEATATBH

//#define HEATELECTRONSATENDRK2
//#define DISSIPATIONFROMGASONLY
#define FORCEGAMMAGASFIXED

//entropy mixing
//#define MIXENTROPIESPROPERLY
//#define UPWINDENTROPYMIXING
//#define DONOTLIMITENTRINMIXING

//silo output
//#define PRINTVISCHEATINGTOSILO
//#define PRINTCOULOMBTOSILO

//floors
#define UEUINTMINRATIO 1.e-3
#define UIUINTMINRATIO 1.e-3
#define TEMPEMINIMAL 1.e2
#define TEMPIMINIMAL 1.e2
#define TEMPEMINIMALFRACTION 1.e-6
#define TEMPIMINIMALFRACTION 1.e-6
#define TEMPEMAXIMALFRACTION 1.e2
#define TEMPIMAXIMALFRACTION 1.e2
#endif

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 2
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .9
#define FLUXMETHOD LAXF_FLUX

#define FLUXLIMITER 0
#define MINMOD_THETA 1.9
#define SHUFFLELOOPS 0

#define DOFIXUPS 1
#define DOU2PRADFIXUPS 0 //ANDREW -- these fixups lead to failures!
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
#define B2UURATIOMAX 1.e3
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 100.
#define GAMMAMAXRAD 20.
#define GAMMAMAXHD 20.
#define RHOFLOOR 1.e-40

/************************************/
//resolution
/************************************/
//total resolution
#define TNX 336//256//32//128//312 
#define TNY 192//192//32//128//200 
#define TNZ 192//192

//number of tiles
#define NTX 16 //28
#define NTY 8 //16
#define NTZ 8 //16

/************************************/
//coordinates
/************************************/
//#define myJETCOORDS
#define myMKS2COORDS
//#define myMKS3COORDS

#define METRICAXISYMMETRIC

#define RH 1.3479
#define RMIN 0.88*RH
#define RMAX 5000.//1.e5

#ifdef myMKS2COORDS //modified Kerr-Shild
#define METRICNUMERIC
#define MYCOORDS MKS2COORDS

#define MKSR0 0.//35
#define MKSH0 0.7
#define MKSMY1 0.002
#define MKSMY2 0.02
#define MKSMP0 1.3

#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY (0.01)
#define MAXY (1.-0.01)

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
#define PRECOMPUTE_MY2OUT // precompute transformation matrices from BL <--> JET

#define MINX 0
#define MAXX 1.
#define Y_OFFSET 0.009
#define MINY -(1.-Y_OFFSET) //10^-8 away from the poles seems to be the last safe point
#define MAXY 1. - Y_OFFSET  

#define MKSR0 0 //-1.35 // Should probably be 0 for jetcoords! (issue with ix=-2 metric)
#define HYPRBRK 1000
#define FJET 0.25
#define FDISK 0.4

//#define RUNI RMIN
//#define RCOLL_JET 1000
//#define RDECOLL_JET 2*RMIN//2*RMIN
//#define RCOLL_DISK 20*RMIN//5*RMIN
//#define RDECOLL_DISK 2*RMIN//2*RMIN

#define RUNI RMIN
#define RCOLL_JET 1000
#define RDECOLL_JET 2*RMIN
#define RCOLL_DISK 5*RMIN
#define RDECOLL_DISK 2*RMIN

#define ALPHA_1 1
#define ALPHA_2 0.375

#define CYLINDRIFY
#define RCYL 20.//4.//10.//10
#define NCYL 1.//0.5//1.
#endif

#define PHIWEDGE (2*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

/************************************/
//boundary conditions
/************************************/

#define SPECIFIC_BC
#define PERIODIC_ZBC
//#define COPY_XBC     //simpler outflowing boundary conditions in  radius

/************************************/
//common physics 
/************************************/
#define GAMMA (13./9.) //(5./3)
#define GAMMAI (5./3.)
#define GAMMAE (4./3.)

#define HFRAC 1. //mass fraction of the hydrogen X
#define HEFRAC 0. //mass fraction of helium Y

/************************************/
//Initial Torus/Atmosphere
/************************************/

#define NTORUS 1

// Fishbone-Moncrief
#define FM_rin 20.
#define FM_rmax 41.
#define FM_rho0 1.//rhoCGS2GU(1.e-18) //density normalization

// atmosphere
#define ATMTYPE 0
#define RHOATMMIN  1.e-4*FM_rho0*pow(2.,-1.5) //rhoCGS2GU(1.e-30)
#define UINTATMMIN  1.e-6*pow(2.,-2.5)/(GAMMA-1.)/3. //(calc_PEQ_ufromTrho(1.e9,RHOATMMIN,0,0,0)) //1.e-6
#define ATMTRADINIT 1.
#define ERADATMMIN  calc_LTE_EfromT(ATMTRADINIT) 

// B-field
#if(NTORUS==0) // single loop v1
#define FM_Aphi_cut 0.2
#undef MAXBETA
#define MAXBETA (100.) //target pgas/pmag inside torus
#define BETANORMFULL
#endif

#if(NTORUS==1) // single loop v2
#define FM_rcut 400.
#define FM_Aphi_cut 0.2
#undef MAXBETA
#define MAXBETA (100.) //target pgas/pmag inside torus
#define BETANORMFULL
#endif

#if(NTORUS==2) // dipolar loops
#undef MAXBETA
#define MAXBETA (100.) //target pgas/pmag inside torus
#define BETANORMFULL
#endif

#if(NTORUS==3) // quadrupolar loops
#undef MAXBETA
#define MAXBETA (100.) //target pgas/pmag inside torus
#define BETANORMFULL
#endif

/************************************/
//output
/************************************/
//#define DTAVG .1   //how frequently to compute quantities included in avg
#define DTOUT1 10  //res
#define DTOUT2 100 //avg
#define DTOUT3 1000 //box,var

//stopping condition

#define NOUTSTOP 1000//2000

//#define DUMPS_READ_HDF5
//#define DUMPS_WRITE_HDF5

#define OUTCOORDS KSCOORDS//BLCOORDS
#define OUTVEL VEL4

#define COORDOUTPUT 0
#define GRIDOUTPUT 0
#define SILOOUTPUT 0
#define OUTOUTPUT 0
#define SIMOUTPUT 5  //needs to be nonzero for phiavg or phicorr
#define SIMOUTPUT_PHIAVG
//#define SIMOUTPUT_PHICORR
//#define CGSOUTPUT
//#define GRTRANSSIMOUTPUT_2

#define RADOUTPUT 0
//#define RADOUTPUTWITHINDTHETA (M_PI/6.)

#define SCAOUTPUT 0
#define AVGOUTPUT 0
#define NORELELAVGS
#define THOUTPUT 0
//#define THPROFRADIUS 30
//#define SILO2D_XZPLANE

