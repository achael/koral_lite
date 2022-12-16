//unfortunately this is what it was
#define MU_GAS 1.
#define MU_I 2.
#define MU_E 2.


/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM 0

/************************************/
//radiation choices
/************************************/
#define RADIATION
// does this work? 
#define BALANCEENTROPYWITHRADIATION

#define COMPTONIZATION
#define U2PCONV 1.e-10
#define RADIMPLICITDAMPINGFACTOR 3.
#define RADIMPLICITMAXENCHANGEDOWN 10.
#define RADIMPLICIMAXENCHANGEUP 10.
#define MAXRADIMPDAMPING 1.e-3
#define RADIMPCONV 1.e-10
#define RADIMPEPS 1.e-6
#define RADIMPMAXITER 40
#define RADIMPCONVREL 1.e-8
#define RADIMPCONVRELERR 1.e-4
#define RADIMPCONVRELENTR 1.e-4
#define RADIMPCONVRELENTRERR 1.e-2
// how does the opacity damping work?
#define OPDAMPINIMPLICIT 0 

//no bremsstrahlung in AGN
#define BREMSSTRAHLUNG
#define SYNCHROTRON
//Should be used, the bridge function is adding a NR component to synchotron opacity
#define NO_SYNCHROTRON_BRIDGE_FUNCTIONS
#define KLEINNISHINA

// try - from Brandon
#define RADIMP_START_WITH_BISECT
#define NPH_START_WITH_BISECT
#define ALLOWRADCEILINGINIMPLICIT
#define ALLOWFORENTRINF4DPRIM

#define RADIMPLICITMAXNPHCHANGE 100 // default is 1e10 - goes to nph floors??? strange, check

//#define SKIPFANCYOPACITIES

//#define EVOLVEPHOTONNUMBER
//we need this for evolvephotonnumber
#define SCALE_JACOBIAN

/************************************/
//magnetic choices
/************************************/
#define MIMICDYNAMO
#define CALCHRONTHEGO
#define THETAANGLE 0.25
#define ALPHAFLIPSSIGN                                                        

#define ALPHADYNAMO (2.*0.314)
#define DAMPBETA
#define BETASATURATED 0.1
#define ALPHABETA (2.*6.28)
#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN
#define MAXBETA .01 //target pmag/pgas int the midplane of the initial torus, overwritten later

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 2 // 1 - lets try better (bfrom Brandon's problem)
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .5
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define DOFIXUPS 1
#define DOU2PRADFIXUPS 0
#define DOU2PMHDFIXUPS 1
#define DORADIMPFIXUPS 0 // 1 let's try without ... 

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY SHEARVISCOSITY
//#define ACCELRADVISCOSITY // not great ? better calc. each step - if the rad. velocities are not changed much since the last step, it will use the old values
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL 0.1

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
//#define POLARAXISAVGIN3D
#define NCCORRECTPOLAR 2

// we are hitting 1e-8
#define UURHORATIOMIN 1.e-8 // 1.e-10 new from Brandon
#define UURHORATIOMAX 1.e0 // 1.e2 new from Brandon

#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1e4 //1.e20 - can go even lower after regridding to force it to relax

#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20


#define B2UURATIOMAX 50. //around the same as uu to rho
// maybe try this (after it is implemented correctly) to avoid breaking of the simulations - can be also usefull after regridding 
// when we want to force it to run until it relaxes a bit (eg after regridding)
//#define B2RHOFLOOR_BACKUP_FFFRAME
#define B2RHOFLOORFRAME DRIFTFRAME
//#define DUINTFROMDRHO - in ZAMOFRAME, maybe better than duint=0


#define B2RHORATIOMAX 50.

#define GAMMAMAXRAD 10.
#define GAMMAMAXHD 10.
//maybe? 
#define RHOFLOOR 1.e-30
// #define UUEERATIOMAX 10. // David - trying to remove the extra uu from ehat

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.

/************************************/
//coordinates / resolution
/************************************/
#define myMKS2COORDS
#define MKSR0 0.1
#define MKSH0 0.9
#define MKSMY1 0.001
#define MKSMY2 0.2
#define MKSMP0 1.5
#define METRICAXISYMMETRIC
#define RMIN 1.85
#define RMAX 500.

#ifdef myMKS2COORDS //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)
#endif

#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY 0.
#define MAXY 1.
#endif

#define PHIWEDGE (M_PI/2.)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

//total resolution
#define TNX 384 //32*10
#define TNY 360 //64*5
#define TNZ 1 //32 //2*8
//number of tiles
#define NTX 32
#define NTY 30
#define NTZ 1

#define SPECIFIC_BC
#define PERIODIC_ZBC

/************************************/
//output
/************************************/
#define DTOUT1 100.
#define DTOUT2 500.


#define BOXOUTPUT 0


#define OUTCOORDS BLCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NOUTSTOP 5000
#define SILOOUTPUT 0
#define OUTOUTPUT 0
#define COORDOUTPUT 0
#define SIMOUTPUT 0
#define RADOUTPUT 0 
#define SCAOUTPUT 1
#define AVGOUTPUT 0 
#define THOUTPUT 0
#define THPROFRAD1US 30
#define SILO2D_XZPLANE

#define PRINT_FIXUPS_TO_SILO

/************************************/
//common physics / torus / atmosphere
/************************************/
#define RHOATMMIN  1.e-24
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN,0,0,0))
#define ATMTRADINIT 3.e6
//define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10) //(calc_LTE_EfromT(ATMTRADINIT)/10)
// From Brandon's collapsing sims
#define ERADATMMIN  (calc_LTE_EfromT(3.e5)/10*6.62/MASS) //If this is set too high, the torus is not in hydrostatic equilibrium since torus will take on atm value

#define GAMMA (5./3.)

#define QUADLOOPS

#define RESTORETORUS
//#define RESCALEDENSITY 0.9
//#define RESCALETORUS 0.8

#define LT_KAPPA 6.e1
#define EXPECTEDHR 0.3
#define LT_XI 0.995
#define LT_R1 20.
#define LT_R2 350.
#define LT_GAMMA 4./3.
#define LT_RIN 35.
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (1./20.) 
