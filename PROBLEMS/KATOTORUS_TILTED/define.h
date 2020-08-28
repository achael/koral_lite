
/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/ 
#define RESTART
//#define RESTARTGENERALINDICES
#define RESTARTNUM -1

/************************************/
//radiation-related choices
/************************************/
//#define RADIATION
//#define RADIMPSTOPWHENFAIL
//#define SKIPRADSOURCE
#define BALANCEENTROPYWITHRADIATION
#define EVOLVEPHOTONNUMBER
#define RADIMP_START_WITH_BISECT
//#define NPH_START_WITH_BISECT
#define ALLOWRADCEILINGINIMPLICIT
#define ALLOWFORENTRINF4DPRIM

#define RADIMPCONV 1.e-8
#define RADIMPEPS 1.e-8
#define RADIMPMAXITER 50
#define RADIMPCONVREL 1.e-6
#define RADIMPCONVRELERR 1.e-2
#define RADIMPCONVRELENTR 1.e-6
#define RADIMPCONVRELENTRERR 1.e-2
#define RADIMPLICITTHRESHOLD 1.e-2 // Is this even in use anymore?
#define RADIMPLICITDAMPINGFACTOR 5.
#define RADIMPLICITMAXNPHCHANGE 100.
#define RADIMPLICITMAXENCHANGEDOWN 100.
#define RADIMPLICITMAXENCHANGEUP 100.
#define RADIMPLICITMAXTECHANGE 5.
#define MAXRADIMPDAMPING 1.e-4
#define OPDAMPINIMPLICIT 0
#define OPDAMPMAXLEVELS 3
#define OPDAMPFACTOR 10.

//opacities
#define SCATTERING
#define BREMSSTRAHLUNG
#define SYNCHROTRON
#define NO_SYNCHROTRON_BRIDGE_FUNCTIONS
//#define SYNCHROTRON_OLD_FIT
#define SUTHERLAND_DOPITA_LAMBDA
//#define METALLICITY 0.

//metallicity
#define HFRAC HFRAC_SUN
#define HEFRAC HEFRAC_SUN
#define MFRAC (1.-HFRAC-HEFRAC)

/************************************/
//rad. viscosity-related choices
/************************************/
#ifdef RADIATION
#define RADVISCOSITY SHEARVISCOSITY
#endif
#define ACCELRADVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL 0.1

/************************************/
//MHD-related choices
/************************************/
//#define U2P_EQS U2P_EQS_JONS
//#define U2P_SOLVER U2P_SOLVER_WP
#define U2PCONV 1.e-12

/************************************/
//non-rel MHD?
/************************************/
//#define NONRELMHD
//#define NONRELMHDENTROPYCUT 1.e-10 // Tcut = 3e12*this number

/************************************/
//magnetic-related choices
/************************************/
#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN
//#define CURL_IN_SPHERICAL
#define MAXBETA .01 //target initial pmag/pgas 
//#define INIT_MAGN_CORNERS

/************************************/
//dynamo-related choices
/************************************/
//#define MIMICDYNAMO
//#define CALCHRONTHEGO
#define THETAANGLE 0.25
#define ALPHAFLIPSSIGN                                                        
#define ALPHADYNAMO 0.314
#define DAMPBETA
#define BETASATURATED 0.1
#define ALPHABETA 6.28

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 2
#define TIMESTEPPING RK2IMEX 
#define TSTEPLIM 0.9
#define FLUXLIMITER 0
#define MINMOD_THETA 1.9
#define SHUFFLELOOPS 0
#define DOFIXUPS 0
#define DORADIMPFIXUPS 0

/************************************/
//blackhole
/************************************/
#define MASS 1.e6
#define BHSPIN 0.0

/************************************/
//coordinates 
/************************************/
#define myMKS2COORDS
//#define mySPHCOORDS
//#define myKSCOORDS
//#define myCYLCOORDS
#define RMIN 1.8
#define RMAX 1.e2
#define MKS_THCUT 0.005
#define MKSR0 0.
#define MKSH0 0.01
#define MKSMY1 0.001
#define MKSMY2 0.2
#define MKSMP0 1.5
#define METRICAXISYMMETRIC

#ifdef myMSPH1COORDS //modified Kerr-Shild
#define PWPOTENTIAL
#define MYCOORDS MSPH1COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(100.-MKSR0))
#define MINY (0.45)
#define MAXY (M_PI-MINY)
#endif

#ifdef mySPHCOORDS //modified Kerr-Shild
#define PWPOTENTIAL
#define MYCOORDS SPHCOORDS
#define MINX RMIN
#define MAXX 100.
#define MINY (0.05)
#define MAXY (M_PI-MINY)
#endif

#ifdef myKSCOORDS //modified Kerr-Shild
//#define PWPOTENTIAL
#define MYCOORDS KSCOORDS
#define MINX RMIN
#define MAXX 100.
#define MINY (0.05)
#define MAXY (M_PI-MINY)
#endif


#ifdef myCYLCOORDS //modified Kerr-Shild
#define PWPOTENTIAL
#define MYCOORDS CYLCOORDS
#define MINX RMIN
#define MAXX 100.
#define MINY (-60.)
#define MAXY 60.
#endif

#ifdef myMKS2COORDS //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY MKS_THCUT//(0.001)
#define MAXY (1.-MKS_THCUT)//(1.-0.001)
#endif

#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS
#define MINX (log(1.85-MKSR0))
#define MAXX (log(100.-MKSR0))
#define MINY 0.
#define MAXY 1.
#endif

#define PHIWEDGE (2.*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

/************************************/
//resolution 
/************************************/
//total resolution
#define TNX 32//28*9
#define TNY 64//26*9
#define TNZ 32 //2*8
//number of tiles
#define NTX 4
#define NTY 4
#define NTZ 4

/************************************/
//boundary conditions 
/************************************/
#define SPECIFIC_BC
#define PERIODIC_ZBC
//#define PERIODIC_XBC
//#define PERIODIC_YBC
#define TRANSMITTING_YBC

/************************************/
//output
/************************************/
#define DTOUT1 1. //res - files
#define DTOUT2 1000.  //avg - files
#define DTOUT3 1. //box,var - files
#define OUTCOORDS BLCOORDS//KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define OUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 1//5000
#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define RADOUTPUT 0
#define SCAOUTPUT 0
#define AVGOUTPUT 0
#define SIMOUTPUT 0
#if(TNZ==1)
#define SILO2D_XZPLANE
//#endif
#else
#define FULLPHI
#endif


#define BOXVERTOUTPUT 0
#if (BOXVERTOUTPUT==1)
#define BOXR1 35.//(pow(pow(19.,3.5)-10.5*(global_time-TSTART),1./3.5))//12. //v=-3./R^2.5
#define BOXR2 40.//(pow(pow(20.,3.5)-10.5*(global_time-TSTART),1./3.5))//14.
#define BOXITH (320/2-31) //distance from eq.plane in cells  
#endif

/************************************/
//physics 
/************************************/
#define GAMMA (4./3.)

/************************************/
//initial torus 
/************************************/
#define NTORUS 1

#if(NTORUS==1) //TDE disk with single poloidal loop (scales with rho)
#define KT_RMIN 5.
#define KT_A 0.0
#define KT_R0 10.
#define KT_RHO0 (3.*rhoCGS2GU(1.e-3)*6.62/MASS)
#define KT_T0 1.0e10
#define KT_DROT 90.*(M_PI/180.)//M_PI/12.
#define RADIUS_POWER 3
#define SIN_POWER 0.5
#define SIN_POWER_THETA 8. //Sin term to damp magnetic field at torus boundary.
#define RHO_CUT_FACTOR 0.
#define MAXBETA_SEPARATE
#undef MAXBETA
#define MAXBETA 0.015
#endif

#if(NTORUS==2) //TDE disk with single poloidal loop (scales with rho)
#define KT_RMIN 5.
#define KT_A 0.0
#define KT_R0 10.
#define KT_RHO0 (3.*rhoCGS2GU(1.e-3)*6.62/MASS)
#define KT_T0 1.0e10
#define KT_DROT 0.*(M_PI/180.)//M_PI/12.
#define RADIUS_POWER 3
#define SIN_POWER 0.5
#define SIN_POWER_THETA 8. //Sin term to damp magnetic field at torus boundary.
#define RHO_CUT_FACTOR 0.
#define MAXBETA_SEPARATE
#undef MAXBETA
#define MAXBETA 0.015
#endif

#define BETANORMFULL

/************************************/
//initial atmosphere 
/************************************/
#define RHOATMMIN  (KT_RHO0*1.e-4)
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN,0,0,0))
#define ERADATMMIN  (calc_LTE_EfromT(3.e4)/10*6.62/MASS) //If this is set too high, the torus is not in hydrostatic equilibrium since torus will take on atm value

/************************************/
//rmhd floors
/************************************/
#define RHOFLOOR 1.e-50
#define UURHORATIOMIN 1.e-8//1K: uu/rho = 7.259162e+12
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 50.
#define GAMMAMAXRAD 20.
#define GAMMAMAXHD 20.

/************************************/
//polar axis
/************************************/
//#ifndef myCYLCOORDS
//#define CORRECT_POLARAXIS
//#endif
//#define POLARAXISAVGIN3D
//#define NCCORRECTPOLAR 2
