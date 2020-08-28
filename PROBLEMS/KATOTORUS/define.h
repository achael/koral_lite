
/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/ 
#define RESTART
//#define RESTARTGENERALINDICES
//#define RESCALEDENSITY 2.
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
//#define ACCELRADVISCOSITY
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
//#define MAGNFIELD
#define GDETIN 1
//#define VECPOTGIVEN
#define MAXBETA .01 //target initial pmag/pgas 

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
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX 
#define TSTEPLIM 0.9
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
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
#define myJETCOORDS
//#define myMKS2COORDS
//#define myMKS3COORDS

#define METRICAXISYMMETRIC

#define RH 2.//1.348
#define RMIN 1.3
#define RMAX 1.e4//1.e5

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
#define MKSR0 0.
#define MKSH0 0.7
#define MKSMY1 0.001
#define MKSMY2 0.2
#define MKSMP0 1.5
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)
#endif

#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS
#define MINX (log(1.85-MKSR0))
#define MAXX (log(100.-MKSR0))
#define MINY 0.
#define MAXY 1.
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
#define HYPRBRK 10000
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
#define ALPHA_2 0.25 //0.375

#define CYLINDRIFY
#define RCYL 20.//4.//10.//10
#define NCYL 1.//0.5//1.
#endif //End of myJETCOORDS

#define PHIWEDGE (0.5*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

/************************************/
//resolution 
/************************************/
//total resolution
#define TNX 64//28*9
#define TNY 128//128//26*9
#define TNZ 1//32//32 //2*8
//number of tiles
#define NTX 16
#define NTY 16
#define NTZ 1//4

/************************************/
//boundary conditions 
/************************************/
#define SPECIFIC_BC
#define PERIODIC_ZBC
//#define PERIODIC_XBC
//#define PERIODIC_YBC

/************************************/
//output
/************************************/
#define DTOUT1 100. //res - files
#define DTOUT2 1000.  //avg - files
#define DTOUT3 1. //box,var - files
#define OUTCOORDS BLCOORDS//KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define OUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 0//5000
#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define RADOUTPUT 0
#define SCAOUTPUT 0
#define AVGOUTPUT 1
#define SIMOUTPUT 0
#define COORDOUTPUT 0
#if(TNZ==1)
#define SILO2D_XZPLANE
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
#define NTORUS 2

#if(NTORUS==1) //Jiang+14
#define KT_A 0.0
#define KT_R0 20.
#define KT_RMIN 5.
#define KT_RHO0 10.*(3.*rhoCGS2GU(1.e-2)*6.62/MASS)
#define KT_T0 2.e10
#undef MAXBETA
#define MAXBETA 0.01 //eq.plane
#endif

#if(NTORUS==2) //TDE disk with single poloidal loop (scales with rho)
#define KT_A 0.0
#define KT_RMIN 5.
#define KT_R0 20.
#define KT_RIN 18. //used for B field setup
#define KT_ROUT 100. //used for B field setup
#define KT_RHO0 10.*(3.*rhoCGS2GU(1.e-3)*6.62/MASS)
#define KT_T0 2.0e10
#define RADIUS_POWER 3
#define SIN_POWER 0.5
#define RHO_CUT_FACTOR 0.

#define BATMZ
#define BATMZ_B0 1.e-22//1.e-18
#define BATMZ_MINR RMIN
#define BATMZ_MAXR 200.//RMAX
#define BATMZ_TORUSRESCALE 1e3 //rescale outside of torus

#define MAXBETA_SEPARATE
#undef MAXBETA
#define MAXBETA 0.003
#endif

#if(NTORUS==3) //TDE disk with single poloidal loop. LRTORUS setup
#define KT_A 0.0
#define KT_R0 40.
#define KT_RHO0 (3.*rhoCGS2GU(1.e-3)*6.62/MASS)
#define KT_T0 2.0e10
#define RADIUS_POWER 3
#define RHO_CUT_FACTOR 0.
//#define MAXBETA_SEPARATE
#undef MAXBETA
#define MAXBETA 0.03
#endif

#if(NTORUS==4) //TDE disk with Quadrupolar loops. LRTORUS setup
#define KT_A 0.
#define KT_RMIN 5.
#define KT_R0 20.
#define KT_RHO0 (3.*rhoCGS2GU(1.e-3)*6.62/MASS)
#define KT_T0 3.8e10
//#define MAXBETA_SEPARATE
#undef MAXBETA
#define MAXBETA 0.03
#endif

#if(NTORUS==5) //TDE disk with Quadrupolar loops. Wavelength of loops is broken power law. 
#define KT_A 0.
#define KT_RMIN 5.
#define KT_R0 20.
#define KT_RHO0 (3.*rhoCGS2GU(1.e-3)*6.62/MASS)
#define KT_T0 2.0e10
//#define MAXBETA_SEPARATE
#undef MAXBETA
#define MAXBETA 0.03
#endif


#define BETANORMFULL

/************************************/
//initial atmosphere 
/************************************/
#define RHOATMMIN  (KT_RHO0*1.e-4)
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN,0,0,0))
#define ERADATMMIN  (calc_LTE_EfromT(3.e4)/10*6.62/MASS)

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
#ifndef myCYLCOORDS
#define CORRECT_POLARAXIS
#endif
//#define POLARAXISAVGIN3D
#define NCCORRECTPOLAR 2
