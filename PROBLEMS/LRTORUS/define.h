/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM 30

/************************************/
//radiation choices
/************************************/
#define RADIATION
#ifdef RADIATION
#define NO_COMPTONIZATION 
//#define EVOLVEPHOTONNUMBER
//#define RADIMPSTOPWHENFAIL
#define RADIMP_START_WITH_BISECT
//#define BALANCEENTROPYWITHRADIATION
#define ALLOWRADCEILINGINIMPLICIT
#define ALLOWFORENTRINF4DPRIM
#endif

//#define DAMPRADWAVESPEEDNEARAXIS
//#define DAMPRADWAVESPEEDNEARAXISNCELLS 6
//#define BASICRADIMPLICIT
#define OPDAMPINIMPLICIT 1
#define OPDAMPMAXLEVELS 3
#define OPDAMPFACTOR 10.
#define RADIMP_START_WITH_BISECT
//#define ALLOWRADCEILINGINIMPLICIT

//#define RADIMPSTOPWHENFAIL
#define MAXDIFFTRADS 1.e4
#define MAXDIFFTRADSNEARBH 1.e2
//#define DAMPCOMPTONIZATIONATBH

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
#define RADIMPLICITMAXTECHANGE 5.
#define MAXRADIMPDAMPING 1.e-4

//opacities
#define SCATTERING
#define BREMSSTRAHLUNG
//#define SYNCHROTRON
//#define NO_SYNCHROTRON_BRIDGE_FUNCTIONS
//#define KLEINNISHINA

/************************************/
//magnetic choices
/************************************/
//#define MIMICDYNAMO
#ifdef MIMICDYNAMO
#define CALCHRONTHEGO
#define THETAANGLE 0.25
#define ALPHAFLIPSSIGN                                                        
#define ALPHADYNAMO 0.314
#define DAMPBETA
#define BETASATURATED 0.1
#define ALPHABETA 6.28
#endif

#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#ifdef RADIATION
#define TIMESTEPPING RK2IMEX
#else
#define TIMESTEPPING RK2HEUN
#endif
#define TSTEPLIM .8
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 1
#define DOU2PRADFIXUPS 1
#define DOU2PMHDFIXUPS 1
#define DORADIMPFIXUPS 1

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY SHEARVISCOSITY
//#define ACCELRADVISCOSITY
#define RADVISCMFPSPH
//#define RADVISCMFPSPHRMIN 10.
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL 0.3

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 100.
#define GAMMAMAXRAD 50.
#define GAMMAMAXHD 50.

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.9

/************************************/
//coordinates / resolution
/************************************/
//#define myMKS2COORDS
#define myMKS3COORDS
#define ROUT 1000.
#define MKSR0 0.
#define MKSH0 0.6
#define MKSMY1 0.0025
#define MKSMY2 0.025
#define MKSMP0 1.2
#define METRICAXISYMMETRIC
#define RMIN 1.3

#ifdef myMKS2COORDS //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(ROUT-MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)
#endif

#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(ROUT-MKSR0))
#define MINY 0.
#define MAXY 1.
#endif

#define PHIWEDGE (M_PI/2.)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

//total resolution

#define TNX 160 //320 //16*17
#define TNY 160 //320 //16*12
#define TNZ 1//32//128 //16*8

//number of tiles

#define NTX 8//32//16//32//16//17
#define NTY 8//32//16//32//8//12
#define NTZ 1//2//8

#define SPECIFIC_BC
#define PERIODIC_ZBC
//#define PERIODIC_XBC
//#define PERIODIC_YBC

/************************************/
//output
/************************************/
//#define OUTPUTPERCORE
#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
//#define RADOUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 50 //5000

#define BOXOUTPUT 0
#define BOXR1 15.
#define BOXR2 20.
#define BOXITH 20 
#define VAROUTPUT 0
#define VARRADIUS 100.
#define NVARCUTS 20

#define DTOUT3 1.
#define DTOUT4 1.

#define COORDOUTPUT 0
#define SILOOUTPUT 1
//#define SIMOUTPUT 2
#define OUTOUTPUT 0
#define RADOUTPUT 0
#define SCAOUTPUT 0
#define AVGOUTPUT 0
#define SILO2D_XZPLANE
#define CBAUTOSCALE
#define DTOUT1 100.
#define DTOUT2 1000.

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)

#define NTORUS 17

#if(NTORUS==81) //
#define LT_KAPPA 2.e2
#define LT_XI 0.705
#define LT_R1 40.
#define LT_R2 1000.
#define UINT_FACTOR 0.25
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 10.
#define BETANORMEQPLANE
//#define BETANORMFACTOR 1.e-04
#undef MAXBETA
#define MAXBETA (.1) 
#endif


#if(NTORUS==79 || NTORUS==80) //
#define LT_KAPPA 5.e2
#define LT_XI 0.705
#define LT_R1 40.
#define LT_R2 1000.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 10.
#define BETANORMEQPLANE
//#define BETANORMFACTOR 1.e-04
#undef MAXBETA
#define MAXBETA (.005) 
#endif

#if(NTORUS==78) //
#define LT_KAPPA 5.e2
#define LT_XI 0.96
#define LT_R1 14.
#define LT_R2 400.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 10.
#define BETANORMEQPLANE
#undef MAXBETA
#define MAXBETA (.1) 
#endif

#if(NTORUS==77) //flat sigma, single poloidal loop
#define LT_KAPPA 5.e2
#define LT_XI 0.975
#define LT_R1 30.
#define LT_R2 200.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 22.
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (.05) 
#endif


#if(NTORUS==17) //flat sigma, ~10Mdot_edd, single loop -> MAD
#define LT_KAPPA 5.e1 //5.e2
#define EXPECTEDHR 0.3
#define LT_XI 0.975
#define LT_R1 30.
#define LT_R2 200.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 22.
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (1./20.)
#endif

#if(NTORUS==7) //flat sigma, ~10Mdot_edd
#define LT_KAPPA 5.e2
#define EXPECTEDHR 0.3
#define LT_XI 0.975
#define LT_R1 30.
#define LT_R2 200.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 22.
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (1./20.)
#endif

#if(NTORUS==6) //for not-so-hyper
#define LT_KAPPA 1.5e3
#define EXPECTEDHR 0.4
#define LT_XI 0.95
#define LT_R1 16.
#define LT_R2 200.
#define LT_GAMMA 4./3.
//#define LT_RIN 10.25
#define LT_RIN 10.6
#undef MAXBETA
#define MAXBETA (1./25.) //eq.plane
#endif

#if(NTORUS==5) //single toroidal loop

#define EXPECTEDHR 0.4
#define LT_KAPPA 1.e-2
#define LT_XI 0.708
#define LT_R1 42.
#define LT_R2 1000.
#define LT_GAMMA 5./3.
#define LT_RIN 10.
#undef MAXBETA
#define MAXBETA (1./30.) //target pmag/pgas inside torus
#define BETANORMFULL
//#define BETANORMFACTOR 2.e-10
#endif

#if(NTORUS==4) //a=0 SANE, no rad, denser loops
#define EXPECTEDHR 0.4
#define LT_KAPPA 1.e-2
#define LT_XI 0.708
#define LT_R1 42.
#define LT_R2 1000.
#define LT_GAMMA 5./3.
#define LT_RIN 10.
#undef MAXBETA
#define MAXBETA (1./30.) //target pmag/pgas inside torus
#define BETANORMFULL
#endif

#if(NTORUS==3) //a=0 SANE, no rad!
#define EXPECTEDHR 0.4
#define LT_KAPPA 1.e-2
#define LT_XI 0.708
#define LT_R1 42.
#define LT_R2 1000.
#define LT_GAMMA 5./3.
#define LT_RIN 10.
#undef MAXBETA
#define MAXBETA (1./30.) //target pmag/pgas inside torus
#define BETANORMFULL
#endif

#if(NTORUS==1) //original (2nd koral paper)
#define LT_KAPPA 1.5e3
#define LT_XI 0.9
#define LT_R1 31.75
#define LT_R2 200.
#define LT_GAMMA 4./3.
#define LT_RIN 15.
#endif

#if(NTORUS==2) //for Yucong?
#define LT_KAPPA 2.e3
#define LT_XI 0.95
#define LT_R1 16.
#define LT_R2 200.
#define LT_GAMMA 4./3.
#define LT_RIN 10.
#endif

#define RHOATMMIN  1.e-22
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN,0,0,0))
#define ERADATMMIN  (calc_LTE_EfromT(1.e6))
