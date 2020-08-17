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
//#define SYNCHROTRON
#define NO_SYNCHROTRON_BRIDGE_FUNCTIONS
//#define SUTHERLAND_DOPITA_LAMBDA
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
#define MINMOD_THETA 1.5
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
#define UURHORATIOMIN 1.e-10 //1.e-10
#define UURHORATIOMAX 1.e2 //1.e2
#define EERHORATIOMIN 1.e-20 //1.e-20
#define EERHORATIOMAX 1.e20 //1.e20
#define EEUURATIOMIN 1.e-20 //1.e-20
#define EEUURATIOMAX 1.e20 //1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 1e4. //100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 10.//100.
#define GAMMAMAXRAD 50.
#define GAMMAMAXHD 10.

/************************************/
//blackhole
/************************************/
#define MASS 1.e6
#define BHSPIN 0.0

/************************************/
//coordinates / resolution
/************************************/
#define myMKS2COORDS
#define ROUT 1.e3
#define RMIN 1.3

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
#define TNX 64
#define TNY 64
#define TNZ 1
//number of tiles
#define NTX 8
#define NTY 8
#define NTZ 4


#define PHIWEDGE (2.*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

#define PERTMAGN 1.e-2
#define SPECIFIC_BC
#define PERIODIC_ZBC

//special definitions for new grid
#define RMIN_DISK 80.//100. //ROUT //special definition used in BC
#define SPECIAL_BC_CHECK
#define STREAM_IX 40//212//149//55
#define STREAM_IYT 29//131//73//73//22 //YBCHI
#define STREAM_IYB 34//86//86//41 //YBCLO
#define STREAM_IZT 62
#define STREAM_IZB 65
//#define SEARCH_STREAM_BOUNDARY //need to run in serial first to use this
//#define FIX_INFLOW_RATE
//#define CELL_BY_CELL_RESCALE //better for arbitrary grid resolutions with MPI
#define STR_USETRADATM //copy trad from inner cell. Don't split recalc Trad from Tgas_init

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 10//1000//120

//#define SIMOUTPUT 2
//#define RAD_INTEGRATION //strictly for getting tau surface in simoutput

#define COORDOUTPUT 0
#define SCAOUTPUT 0
#define SILOOUTPUT 1
#define RADOUTPUT 0
#define AVGOUTPUT 0
#define SILO2D_XZPLANE
#define PRINTXGC_RIGHT
//#define PRINTEACHT

#define DTOUT1 1.//100.
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
#define MDOT_TIME_DEPENDENT
#define USE_TPLUSTFB

#define EPSILON0 -2.15e-4*(1.)*pow((MASS/1e6),1./3.) //Base epsilon for model, defines t_fb
#define TFB0 2.*M_PI*pow((-1./(2.*EPSILON0)),1.5) //Fallback time of most bound material
//#define DISKVR (-1.*sqrt(2.*0.1/RMIN_DISK))
#define DISKHR 0.05
#define DISKHCGS (DISKHR*RMIN_DISK*MASS*147700.) 
#define MDOTEDD0 100.// initial inflow rate in units of Eddington accretion rate
//#define MDOTIN (MDOTEDDIN*2.48e18*MASS)
//#define DISKSIGMA surfdensCGS2GU(-MDOTIN/(2.*M_PI*RMIN_DISK*147700.*MASS*3.e10*DISKVR))
#define DISKRCIR 20. //circularization radius
#define DISKTEMP 1.e4//1.e4//DISKHR*DISKHR*(5./7.)/RMIN_DISK/(1.5e-13)
//#define VERTBTIME 1000.
#define VERTICALB

///#define MAXRADIUS4DYNAMO DISKRCIR

#define DISKH (DISKHR*RMIN_DISK)
//#define DISKRHO (DISKSIGMA/2./DISKH)

#define MAGNOMEGA 0.//(2.*M_PI/2000.)//0.
//#define MAGNOMEGA 0.*(1.5e-2*pow(200./50.,-1.5)) //omega should follow the free fall time?
#define MAGBETA 3.e-2
//#define SANEFIELD

//atmosphere
#define RHOATMMIN  1.e-21
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e8,RHOATMMIN,0,0,0))
#define ATMTRADINIT 1.e5
#define ERADATMMIN  (calc_LTE_EfromT(ATMTRADINIT))

