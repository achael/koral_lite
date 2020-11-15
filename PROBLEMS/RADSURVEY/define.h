/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE


//***********************************/
// Main problem specification
// SPIN, SANE vs MAD, Coordinates, Torus
//***********************************/
#define SPINp25     //or: SPINp5 SPINp75 SPINm25 SPINm5 SPINm75
#define NTORUS 3    //3 for SANE, 1 for MAD

#define myJETCOORDS //or: myMKS2COORDS, myMKS3COORDS
#define LIMOTORUS   //or: FISHMONCTORUS

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

#if defined(SPINm75)
#define BHSPIN -0.75
#define RH 1.66144   //a=.75 1.66144, a=.5 1.86603, a=.25 1.96825

#elif defined(SPINm5)
#define BHSPIN -0.5
#define RH 1.86603   //a=.75 1.66144, a=.5 1.86603, a=.25 1.96825

#elif defined(SPINm25)
#define BHSPIN -0.25
#define RH 1.96825   //a=.75 1.66144, a=.5 1.86603, a=.25 1.96825

#elif defined(SPINp25)
#define BHSPIN 0.25
#define RH 1.96825   //a=.75 1.66144, a=.5 1.86603, a=.25 1.96825

#elif defined(SPINp5)
#define BHSPIN 0.5
#define RH 1.86603   //a=.75 1.66144, a=.5 1.86603, a=.25 1.96825

#elif defined(SPINp75)
#define BHSPIN 0.75
#define RH 1.66144   //a=.75 1.66144, a=.5 1.86603, a=.25 1.96825
#endif

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
#define TSTEPLIM .8
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
#define GAMMAMAXRAD 25.
#define GAMMAMAXHD 25.
#define RHOFLOOR 1.e-40

/************************************/
//resolution
/************************************/
//total resolution
#define TNX 320 //336//256//32//128//312 
#define TNY 192 //256 //192//192//32//128//200 
#define TNZ 1//192 //192

#define SILO2D_XZPLANE

//number of tiles
#define NTX 8 //28
#define NTY 6 //16
#define NTZ 6 //8 //16

/************************************/
//coordinates
/************************************/

#define METRICAXISYMMETRIC

#define RMIN 0.825*RH
#define RMAX 1.e5

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
#define Y_OFFSET 0.01
#define MINY -(1.-Y_OFFSET) //10^-8 away from the poles seems to be the last safe point
#define MAXY 1.-Y_OFFSET  

#define MKSR0 0 //-1.35 // Should probably be 0 for jetcoords! (issue with ix=-2 metric)
#define HYPRBRK RMAX
#define FJET 0.4
#define FDISK 0.4

#define RUNI RMIN
#define RCOLL_JET 10000
#define RDECOLL_JET 2.//2*RMIN
#define RCOLL_DISK 10.//5*RMIN
#define RDECOLL_DISK 2.//2*RMIN

#define ALPHA_1 1
#define ALPHA_2 0.2//0.375

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
#define GAMMA  (5./3.) //(13./9.)
#define GAMMAI (5./3.)
#define GAMMAE (4./3.)

#define HFRAC 1. //mass fraction of the hydrogen X
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
#define LT_GAMMA 5./3.
#define LT_RIN 12.       //10.5

#else

// Fishbone-Moncrief torus
#define FM_rin 20.
#define FM_rmax 41.
#define FM_rho0 1.//rhoCGS2GU(1.e-18) //density normalization

#endif

// atmosphere //TODO better scalings!
#define ATMTYPE 0
#define RHOATMMIN  1.e-4*pow(2.,-1.5) //*FM_rho0 
#define UINTATMMIN  1.e-6*pow(2.,-2.5)/(GAMMA-1.)/3. 
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
#define DTOUT2 100  //avg
#define DTOUT3 1000 //box,var

//stopping condition

#define NOUTSTOP 0 //1000//2000

//#define DUMPS_READ_HDF5
//#define DUMPS_WRITE_HDF5

#define OUTCOORDS BLCOORDS //KSCOORDS
#define OUTVEL VEL4

#define COORDOUTPUT 0
#define GRIDOUTPUT 0
#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define SIMOUTPUT 0  //needs to be nonzero for phiavg or phicorr
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


