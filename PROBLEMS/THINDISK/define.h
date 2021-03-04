//unfortunately this is what it was:
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
//#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
#define RADIATION
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
#define OPDAMPINIMPLICIT 0
#define MAXRADIMPDAMPING 1.e-3

/************************************/
//magnetic choices
/************************************/
#define MIMICDYNAMO
//#define CALCHRONTHEGO
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
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX//HEUN//IMEX
#define TSTEPLIM .5
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define DOFIXUPS 1
#define DOU2PRADFIXUPS 0
#define DOU2PMHDFIXUPS 1
#define DORADIMPFIXUPS 1

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY SHEARVISCOSITY
#define ACCELRADVISCOSITY
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
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 50.
#define GAMMAMAXRAD 10.
#define GAMMAMAXHD 10.

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.

/************************************/
//coordinates / resolution
/************************************/
#define myMKS2COORDS
#define MKSR0 0.
#define MKSH0 0.875
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
#define TNX 320 //32*10
#define TNY 320 //64*5
#define TNZ 1 //32 //2*8
//number of tiles
#define NTX 20
#define NTY 40
#define NTZ 2

#define SPECIFIC_BC
#define PERIODIC_ZBC

/************************************/
//output
/************************************/
#define DTOUT1 50.
#define DTOUT2 100.


#define BOXOUTPUT 0
#define BOXVERTOUTPUT 0
//#undef CALCHRONTHEGO
#define TSTART 64000.
#define BOXR1 6.//(pow(pow(19.9,3.5)-10.5*(global_time-TSTART),1./3.5))//12. //v=-3./R^2.5
#define BOXR2 7.//(pow(pow(20.,3.5)-10.5*(global_time-TSTART),1./3.5))//14.
//y for testing the runonthego feature of ./ana
//#define NSTEPSTOP 2.
#define PRINTEACHT
#define BOXITH (320/2-31) //distance from eq.plane in cells  
//#define BOXITH (320/2) //distance from eq.plane in cells  


#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
//#define RADOUTPUTINZAMO
//#define RADOUTPUTWITHINDISK
#define NOUTSTOP 5000
#define SILOOUTPUT 0
#define OUTOUTPUT 0
#define COORDOUTPUT 0
#define SIMOUTPUT 4 //1
//#define GRTRANSSIMOUTPUT
#define RADOUTPUT 0 
#define SCAOUTPUT 1
#define AVGOUTPUT 0 
#define THOUTPUT 0
#define THPROFRAD1US 30
#define SILO2D_XZPLANE

/************************************/
//common physics / torus / atmosphere
/************************************/
#define RHOATMMIN  1.e-24
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN,0,0,0))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)
#define GAMMA (5./3.)

#define NTORUS 8

#define RESTORETORUS
//#define RESCALEDENSITY 0.9
#define RESCALETORUS 0.8

#if(NTORUS==8) //thinnish, colder disk
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
#endif



#if(NTORUS==7) //flat sigma, r001 in 3d paper
#define LT_KAPPA 5.e2
#define EXPECTEDHR 0.3
#define LT_XI 0.975
#define LT_R1 30.
#define LT_R2 200.
#define LT_GAMMA 4./3.
#define LT_RIN 22.
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (1./10.) //eq.plane
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

