/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

//***********************************/
// Main problem specification
// SPIN, SANE vs MAD, Coordinates, Torus
// Modified from Ramesh's survey
//***********************************/
#define SPINp5     //or: SPINp9 SPINp7 SPINp5 SPIN0 SPINm5 SPINm7 SPINm9
#define NTORUS 1   //3 for SANE, 1 for MAD

#define myMKS2COORDS //or: myMKS2COORDS, myMKS3COORDS
#define FISHMONCTORUS   // LIMOTORUS or FISHMONCTORUS

/************************************/
//restart
/************************************/
//#define NORESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1
#define PERTURBUINT 0.02
//#define PERTURBAFTERINIT //will not perturb after init if this is not defined

/************************************/
//blackhole
/************************************/
#define MASS 6.5e9   

// this defines the spin and sets the grid inner radius and torus parameters

#if defined(SPINm9)
#define BHSPIN -0.9
#define RH 1.43589   //a=-0.9 1.43589

#elif defined(SPINm75)
#define BHSPIN -0.75
#define RH 1.66144   //a=-.75 1.66144

#elif defined(SPINm7)
#define BHSPIN -0.7
#define RH 1.71414   //a=-0.7 1.71414

#elif defined(SPINm5)
#define BHSPIN -0.5
#define RH 1.86603   //a=-.5 1.86603

#elif defined(SPINm3)
#define BHSPIN -0.3
#define RH 1.95394   //a=-0.3 1.95393

#elif defined(SPINm25)
#define BHSPIN -0.25
#define RH 1.96825   //a=-.25 1.96825

#elif defined(SPIN0)
#define BHSPIN 0.0
#define RH 2.00000   //a=0.0 2.00000

#elif defined(SPINp25)
#define BHSPIN 0.25
#define RH 1.96825   //a=.25 1.96825

#elif defined(SPINp3)
#define BHSPIN 0.3
#define RH 1.95394   //a=0.3 1.95393

#elif defined(SPINp5)
#define BHSPIN 0.5
#define RH 1.86603   //a=.5 1.86603

#elif defined(SPINp7)
#define BHSPIN 0.7
#define RH 1.71414   //a=0.7 1.71414

#elif defined(SPINp75)
#define BHSPIN 0.75
#define RH 1.66144   //a=.75 1.66144

#elif defined(SPINp9)
#define BHSPIN 0.9
#define RH 1.43589   //a=0.9 1.43589

#endif

/************************************/
//U2P-related choices
/************************************/
//#define U2P_EQS U2P_EQS_JONS
//#define U2P_SOLVER U2P_SOLVER_WP
#define U2PCONV 1.e-12

/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1    //must be 1 for MAGNFIELD
#define VECPOTGIVEN
#define INIT_MAGN_CORNERS //initialize magnetic field on corners/edges (which?)
#define MAXBETA_SEPARATE // find maxima of ptot and pmag independently and take ratio
#define MAXBETA .01 //target initial pmag/pgas int the midplane


//#define FORCEFREE
//#define SKIPALLFLOORS

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
//reconstruction / Courant
/************************************/
#define INT_ORDER 2
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .8
#define FLUXMETHOD LAXF_FLUX

#define FLUXLIMITER 1
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0

#define DOFIXUPS 1
#define DOU2PRADFIXUPS 0 //ANDREW -- these fixups lead to failures!
#define DOU2PMHDFIXUPS 1
#define DORADIMPFIXUPS 1

/************************************/
//rmhd floors
/************************************/
#if defined(FORCEFREE)

//#define REDUCEORDERWHENNEEDED
//#define REDUCEORDERATBH
//#define REDUCEORDERFF

//#define ENFORCEENTROPY

#define HYBRID_FORCEFREE
#define HYBRID_FORCEFREE_SIGMACUT 50
#define HYBRID_FORCEFREE_WIDTH 0.25
#define FORCEFREE_SOLVE_PARALLEL

#define FORCEFREE_PARALLEL_COLD
//#define FORCEFREE_PARALLEL_ENTROPY
//#define NO_FORCEFREE_PARALLEL_SOURCETERM
#define FORCEFREE_NO_PARALLEL_ATBH

//#define NOLOGINS // TODO do we need this in the new version? 

#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-6
#define UURHORATIOMAX 10. //1.e3
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 1.e10
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 1.e10
#define GAMMAMAXFF 1000.
#define GAMMAMAXHD 1000.
#define RHOFLOOR 1.e-14
#define B2RHOFLOORFRAME FFFRAME //DRIFTFRAME
//#define RHOFLOOR_INIT

#else
//#define NOLOGINS
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-6
#define UURHORATIOMAX 10.
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 50.
#define GAMMAMAXHD 25.
#define RHOFLOOR 1.e-14
#define B2RHOFLOORFRAME ZAMOFRAME
#endif

/************************************/
//resolution
/************************************/
//total resolution
#define TNX 128 //320
#define TNY 128 //192
#define TNZ 1//64 //192

#define SILO2D_XZPLANE
//#define CURRENTTIMEDERIV

//number of tiles
#define NTX 8 
#define NTY 8 
#define NTZ 4 

/************************************/
//coordinates
/************************************/

#define METRICAXISYMMETRIC

#define RMIN 0.7*RH
#define RMAX 1.e4

#ifdef myMKS1COORDS //modified Kerr-Shild
#define METRICNUMERIC
#define MYCOORDS MKS1COORDS

#define MKSR0 0. //-1.35
#define MKSH0 0.7

#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY (0.01)*M_PI //(0.005)*M_PI
#define MAXY (1. - 0.01)*M_PI //(1. - 0.005)*M_PI
#endif

#ifdef myMKS2COORDS //modified Kerr-Shild
#define METRICNUMERIC
#define MYCOORDS MKS2COORDS

#define MKSR0 0. //-1.35
#define MKSH0 0.7

#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY (0.005)
#define MAXY (1. - 0.005)
#endif

#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS

#define MKSR0 0. //-1.35
#define MKSH0 0.7
#define MKSMY1 0.002
#define MKSMY2 0.02
#define MKSMP0 1.3

#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
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
#define Y_OFFSET 1.e-5
#define MINY -(1.-Y_OFFSET) //10^-8 away from the poles seems to be the last safe point
#define MAXY 1.-Y_OFFSET  

#define MKSR0 0 //-1.35 // Should be 0 for jetcoords! (issue with ix=-2 metric)
#define HYPRBRK 1.e4 //RMAX
#define FJET 0.3 // 0.4
#define FDISK 0.4

#define RUNI RMIN
#define RCOLL_JET 1000     //10000
#define RDECOLL_JET 2.*RH  //2*RMIN
#define RCOLL_DISK 20.*RH  //5*RMIN
#define RDECOLL_DISK 2.*RH //2*RMIN

#define ALPHA_1 1
#define ALPHA_2 0.25 //0.2

#define CYLINDRIFY
#define RCYL 30. //20.
#define NCYL 1.  //0.5
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
#define GAMMA  (13./9.)
#define GAMMAI (5./3.)
#define GAMMAE (4./3.)

#define HFRAC 1.  //mass fraction of the hydrogen X
#define HEFRAC 0. //mass fraction of helium Y

/************************************/
//Initial Torus/Atmosphere
/************************************/

#ifdef LIMOTORUS
// Penna et all 2013 limotorus
#define LT_KAPPA 0.008   //This density normalization gives peak density ~1
                         // in code units
#define LT_MAXDENS 1000  //This maximum density value will remove issues
                         // with high density points outside the disk
                         // we see in retrograde esp.
                         // be careful with scale,esp with LT_KAPPA!

//LT_XI sets the torus size -- these values chosen to give outer radius ~500
//TODO fill in for intermediate values
//a=0.75 0.752, a=0.5 0.757, a=0.25 0.761
// a=-0.25 0.771, a=-.5 0.776, a=-.75 0.782
#if defined(SPINm75)
#define LT_XI 0.782
#elif defined(SPINm5)
#define LT_XI 0.776
#elif defined(SPINm25)
#define LT_XI 0.771
#elif defined(SPINp25)
#define LT_XI 0.761
#elif defined(SPINp5)
#define LT_XI 0.757
#elif defined(SPINp75)
#define LT_XI 0.752
#endif

#define LT_R1 42.        //42.
#define LT_R2 800.       //500.
#define LT_GAMMA 13./9.
#define LT_RIN 12.       //10.5

#else

// Fishbone-Moncrief torus
#define FM_rin 20.
#define FM_rho0 1. //density normalization

#if defined(SPINm9)
#define FM_rmax 43.06
#elif defined(SPINm7)
#define FM_rmax 42.90
#elif defined(SPINm5)
#define FM_rmax 42.75
#elif defined(SPINm3)
#define FM_rmax 42.62
#elif defined(SPIN0)
#define FM_rmax 42.43
#elif defined(SPINp3)
#define FM_rmax 42.25
#elif defined(SPINp5)
#define FM_rmax 42.15
#elif defined(SPINp7)
#define FM_rmax 42.05
#elif defined(SPINp9)
#define FM_rmax 41.96
#endif


#endif

// atmosphere //TODO better scalings!
#define ATMTYPE 0
#define RHOATMMIN   1.e-5*pow(2.,-1.5)*FM_rho0 
#define UINTATMMIN  1.e-7*pow(2.,-2.5)*FM_rho0/(GAMMA-1.)/3. 
#define ATMTRADINIT 2.7
#define ERADATMMIN  calc_LTE_EfromT(ATMTRADINIT) 

// B-field
#if(NTORUS==0) // single loop v1
#define Aphi_rho_cut 0.2
#undef MAXBETA
#define MAXBETA (100.) //target pgas/pmag inside torus
#define BETANORMFULL
#endif

#if(NTORUS==1) // single loop v2
#define Aphi_r_cut 400.
#define Aphi_rho_cut 0.2
#undef MAXBETA
#define MAXBETA (100.) //target pgas/pmag inside torus
#define BETANORMFULL
#endif

#if(NTORUS==2) // dipolar loops
#define Aphi_lambda 0.5
#define Aphi_rchop 400.
#undef MAXBETA
#define MAXBETA (100.) //target pgas/pmag inside torus
#define BETANORMFULL
#endif

#if(NTORUS==3) // quadrupolar loops
#define Aphi_lambda 1.
#define Aphi_rchop 400.
#undef MAXBETA
#define MAXBETA (100.) //target pgas/pmag inside torus
#define BETANORMFULL
#endif

/************************************/
//output
/************************************/
//#define DTAVG .1  //how frequently to compute quantities included in avg
#define DTOUT1 10   //res
#define DTOUT2 1000 //avg
#define DTOUT3 1000 //box,var

//stopping condition

#define NOUTSTOP 0//1000 //2000

//#define DUMPS_READ_HDF5
//#define DUMPS_WRITE_HDF5
//#define COORDOUTPUT_HDF5

//#define ANAOUT_HDF5
#define OUTCOORDS2 KSCOORDS

#define OUTCOORDS KSCOORDS
#define OUTVEL VEL4

//#define CGSOUTPUT
#define COORDOUTPUT 0
#define GRIDOUTPUT 0
#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define SIMOUTPUT 0  //needs to be nonzero for phiavg or phicorr
//#define SIMOUTPUT_PHIAVG
//#define SIMOUTPUT_PHICORR
//#define GRTRANSSIMOUTPUT_2

#define RADOUTPUT 0
//#define RADOUTPUTWITHINDTHETA (M_PI/6.)

#define SCAOUTPUT 0
#define AVGOUTPUT 0
#define NORELELAVGS
#define THOUTPUT 0
//#define THPROFRADIUS 30


