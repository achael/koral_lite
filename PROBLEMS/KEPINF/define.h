//************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM -1
//#define PERTURBRHO 0.01

/************************************/
//radiation choices
/************************************/
#define RADIATION
#ifdef RADIATION
//#define EVOLVEPHOTONNUMBER
//#define NO_COMPTONIZATION
#ifndef EVOLVEPHOTONNUMBER
#ifndef NO_COMPTONIZATION
#define COMPTONIZATION
#endif
#endif
//#define RADIMPSTOPWHENFAIL
#define RADIMP_START_WITH_BISECT
//#define BALANCEENTROPYWITHRADIATION
#define ALLOWRADCEILINGINIMPLICIT
#define ALLOWFORENTRINF4DPRIM
#define OPDAMPINIMPLICIT 1
#endif
#define U2PCONV 1.e-12
#define RADIMPCONVREL 1.e-6
#define RADIMPCONVRELERR 1.e-2
#define RADIMPCONVRELENTR 1.e-6
#define RADIMPCONVRELENTRERR 1.e-2
#define RADIMPCONV 1.e-8
#define RADIMPENTRCONV 1.e-8
#define RADIMPEPS 1.e-8
#define RADIMPMAXITER 50
#define RADIMPLICITDAMPINGFACTOR 5.
#define RADIMPLICITMAXNPHCHANGE 100.
#define RADIMPLICITMAXENCHANGEDOWN 100.
#define RADIMPLICITMAXENCHANGEUP 100.
#define MAXRADIMPDAMPING 1.e-4


/************************************/
//opacities
/************************************/

//opacities
#define SCATTERING
//#define NO_COMPTONIZATION
#define BREMSSTRAHLUNG
#define SYNCHROTRON
#define NO_SYNCHROTRON_BRIDGE_FUNCTIONS
#define SUTHERLAND_DOPITA_LAMBDA
//#define BOUNDFREE
//#define KLEINNISHINA
//#define BOUNDELECTRON
//#define COLDOPACITIES
//#define SKIPFANCYOPACITIES


/************************************/
//magnetic choices
/************************************/
//if we want a magnetic field, uncomment MAGNFIELD
#define MAGNFIELD
#define GDETIN 1
//#define MIMICDYNAMO
//#define EXPECTEDHR 1.
//#define CALCHRONTHEGO
#define THETAANGLE 0.25
#define ALPHAFLIPSSIGN                                                        
#define ALPHADYNAMO 0.314
#define DAMPBETA
#define BETASATURATED 0.1
#define ALPHABETA 6.28
//#define STREAMDAMPFAC 0.1 //more damping in stream? (damp propto R^-(3/2)) - Brandon

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 2
#ifdef RADIATION
#define TIMESTEPPING RK2IMEX
#else
#define TIMESTEPPING RK2HEUN
#endif
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.0
#define DOFIXUPS 0
#define DOU2PRADFIXUPS 0
#define DOU2PMHDFIXUPS 0
#define DORADIMPFIXUPS 0

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY SHEARVISCOSITY
//#define ACCELRADVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define ALPHARADVISC .1
#define RADVISCMAXVELDAMP
#define MAXRADVISCVEL 1.
//#define RADVISCSTARTTIME 10000.

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2
#define RHOFLOOR 1.e-50
#define UURHORATIOMIN 1.e-8 //1.e-10
#define UURHORATIOMAX 1.e0 //1.e2
#define EERHORATIOMIN 1.e-20 //1.e-20
#define EERHORATIOMAX 1.e20 //1.e20
#define EEUURATIOMIN 1.e-20 //1.e-20
#define EEUURATIOMAX 1.e20 //1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000. //100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 50.//50.//100.
#define GAMMAMAXRAD 20.//30.
#define GAMMAMAXHD 20.//20.

/************************************/
//blackhole
/************************************/
#define MASS 1.e6
#define BHSPIN 0.0

/************************************/
//coordinates / resolution
/************************************/
#define myMKS2COORDS
#define ROUT 1.e4
#define RMIN 1.8

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MKSR0 0.
#define MYCOORDS MKS1COORDS
#define MINX (log(3.6-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.0025*Pi/2.)
#define MAXY (Pi-0.0025*Pi/2.)
#endif

#ifdef myMKS2COORDS //modified Kerr-Shild with more cells towards the eq.plane
#define MKSR0 0.
#define MKSH0 0.7 //makes cells smaller towards equatorial plane
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(ROUT-MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)
#endif

//total resolution
#define TNX 160
#define TNY 160//96
#define TNZ 1
//number of tiles
#define NTX 16
#define NTY 16//4//16
#define NTZ 1//16//4


#define PHIWEDGE (2.*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

#define PERTMAGN 1.e-2
#define SPECIFIC_BC
#define PERIODIC_ZBC

//special definitions for new grid
#define RMIN_DISK 40.//100. //ROUT //special definition used in BC
#define SPECIAL_BC_CHECK
#define STREAM_IX 58
#define STREAM_IYT 73
#define STREAM_IYB 86
#define STREAM_IZT 63
#define STREAM_IZB 64
//#define SEARCH_STREAM_BOUNDARY //need to run in serial first to use this

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
//#define RADOUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 150//1000//120

#define SIMOUTPUT 0
//#define RAD_INTEGRATION //strictly for getting tau surface in simoutput

#define COORDOUTPUT 0
#define SCAOUTPUT 1
#define SILOOUTPUT 1
#define RADOUTPUT 0
#define AVGOUTPUT 0
#define SILO2D_XZPLANE
#define PRINTXGC_RIGHT
//#define PRINTEACHT

#define DTOUT1 100.//100.
#define DTOUT2 1000.

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)
#define HFRAC HFRAC_SUN
#define HEFRAC HEFRAC_SUN
#define MFRAC (1. - HFRAC - HEFRAC)

//parameters for inflow
//#define EPSILON_TIME_DEPENDENT 
//#define MDOT_TIME_DEPENDENT //Make inflow rate decay following T^(-5/3)
#define USE_TPLUSTFB

//Inputs for stream
#define EPSILON0 -7.5e-3 // not used yet...
#define TVISC 10.*pow(RMIN_DISK,1.5) // viscous delay from injection point
#define TFB0 5000. //Decay rate during injection (if MDOT_TIME_DEPENDENT is on)
#define TSTARTDECAY 25000. //TFB0
//#define TSHUTOFF 20. //forced shutoff time (debugging only)
//#define SHUTOFF_AFTER_TSHUTOFF
#define DISKHR 0.05
#define DISKHCGS (DISKHR*RMIN_DISK*MASS*147700.) 
#define MDOTEDD0 1000.//81.5 //initial inflow rate in units of Eddington accretion rate (~Mdot_target/10)
//#define MDOTIN (MDOTEDDIN*2.48e18*MASS)
//#define DISKSIGMA surfdensCGS2GU(-MDOTIN/(2.*M_PI*RMIN_DISK*147700.*MASS*3.e10*DISKVR))
#define DISKRCIR RMIN_DISK //20. //circularization radius
#define DISKTEMP 1.e9//1.e4//DISKHR*DISKHR*(5./7.)/RMIN_DISK/(1.5e-13)
//#define COPY_DOMAIN_UINT //copy disk uint from domain (after TVISC)
//#define COPY_DOMAIN_EHAT //copy ehat from inner cell
//#define VERTBTIME 1000.
#define VERTICALB

///#define MAXRADIUS4DYNAMO DISKRCIR

#define DISKH (DISKHR*RMIN_DISK)
//#define DISKRHO (DISKSIGMA/2./DISKH)

#define MAGNOMEGA 0.//(2.*M_PI/500.)//0.
//#define MAGNOMEGA 0.*(1.5e-2*pow(200./50.,-1.5)) //omega should follow the free fall time?
#define MAGBETA 3.e-2
#define SANEFIELD

//atmosphere
#define RHOATMMIN  1.e-21
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e8,RHOATMMIN,0,0,0))
#define ATMTRADINIT 1.e5
#define ERADATMMIN  (calc_LTE_EfromT(ATMTRADINIT))

