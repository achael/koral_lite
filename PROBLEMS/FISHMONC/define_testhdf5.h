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
#define MAXBETA .01 //target initial pmag/pgas 
#define MAXBETA_SEPARATE // find maxima of ptot and pmag independently and take ratio
#define INIT_MAGN_CORNERS // for consistent definition of vector potential

/************************************/
//dynamo-related choices
/************************************/
//#define MIMICDYNAMO  // should be used only for 2D
//#define CALCHRONTHEGO
//#define THETAANGLE 0.25
//#define ALPHAFLIPSSIGN
//#define ALPHADYNAMO 0.314
//#define DAMPBETA
//#define BETASATURATED 0.1
//#define ALPHABETA 6.28

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 2
#define TIMESTEPPING RK2IMEX 
#define TSTEPLIM 0.8
//#define TSTEPLIM 0.45
#define FLUXLIMITER 0
#define MINMOD_THETA 1.9
#define SHUFFLELOOPS 0
#define DOFIXUPS 0
#define DORADIMPFIXUPS 0

/************************************/
//blackhole
/************************************/
#define MASS 1.e1
#define BHSPIN 0.9375

/************************************/
//coordinates 
/************************************/
//#define myMKS2COORDS
#define myMKS3COORDS
//#define mySPHCOORDS
//#define myKSCOORDS
//#define myCYLCOORDS
#define RMIN 1.15
#define RMAX 50.
#define MKSR0 0.
#define MKSH0 0.6
#define MKSMY1 0.0025
#define MKSMY2 0.025
#define MKSMP0 1.2
#define METRICAXISYMMETRIC

#ifdef myMSPH1COORDS //modified Kerr-Shild
#define PWPOTENTIAL
#define MYCOORDS MSPH1COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY (0.45)
#define MAXY (M_PI-MINY)
#endif

#ifdef mySPHCOORDS //modified Kerr-Shild
#define PWPOTENTIAL
#define MYCOORDS SPHCOORDS
#define MINX RMIN
#define MAXX RMAX
#define MINY (0.05)
#define MAXY (M_PI-MINY)
#endif

#ifdef myKSCOORDS //modified Kerr-Shild
//#define PWPOTENTIAL
#define MYCOORDS KSCOORDS
#define MINX RMIN
#define MAXX RMAX
#define MINY (0.05)
#define MAXY (M_PI-MINY)
#endif


#ifdef myCYLCOORDS //modified Kerr-Shild
#define PWPOTENTIAL
#define MYCOORDS CYLCOORDS
#define MINX RMIN
#define MAXX RMAX
#define MINY (-60.)
#define MAXY 60.
#endif

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

#define PHIWEDGE (2. * M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

/************************************/
//resolution 
/************************************/
//total resolution
#define TNX 8 // 96
#define TNY 8 // 96
#define TNZ 8 // 96
//number of tiles
#define NTX 1
#define NTY 1
#define NTZ 1

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
#define DTOUT1 10. //res - files
#define DTOUT2 1000.  //avg - files
#define DTOUT3 1. //box,var - files

#define OUTCOORDS BLCOORDS//KERRCOORDS
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define OUTPUTINZAMO
#define NSTEPSTOP 1.e12
#define NOUTSTOP 10 // 1000

#define FOLDER_HDF5 "./hdf5"
//#define COORDOUTPUT_HDF5
#define DUMPS_READ_HDF5
#define DUMPS_WRITE_HDF5

//#define SILOOUTPUT 1
//#define OUTOUTPUT 0
//#define RADOUTPUT 1
//#define SCAOUTPUT 1
#define AVGOUTPUT 1
//#define SIMOUTPUT 2
//#define COORDOUTPUT 1
#define INIT_GU_OUTPUT
#if(TNZ==1)
#define SILO2D_XZPLANE
#else
#define FULLPHI
#endif


//#define BOXVERTOUTPUT 1
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
//initial torus: Fishbone-Moncrief
/************************************/
#define NTORUS 1

#if(NTORUS==1) 
#define FM_rin 6.
#define FM_rmax 12.
#define FM_rho0 1.
#define FM_Aphi_cut 0.2
#endif

//#define PERTURBVZ 0.05
//#define PERTURBUINT 0.02
#define BETANORMFULL  // compute BETA over the whole torus, not just the equatorial plane

#define RADOUTPUTWITHINDTHETA (M_PI / 6.)

/************************************/
//initial atmosphere (constants correspond to rout=2)
/************************************/
#define RHOATMMIN (3.5355e-6)
#define UINTATMMIN (1.7678e-8)
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10*6.62/MASS)

/************************************/
//rmhd floors
/************************************/

// For density floor, use only one of the following: RHOFLOOR, RHOFLOOR_BH, RHOFLOOR_INIT. If using either of the latter, set RHOFLOOR to a very low value
#define RHOFLOOR 1.e-50
//#define RHOFLOOR_INIT

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
