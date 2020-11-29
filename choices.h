/*! \file ko.h
 \brief default choices
*/

#define RESTARTGENERALINDICES

//small and big numbers
#define SMALL 1.e-80 
#define BIG (1./SMALL) 
#define NUMEPSILON DBL_EPSILON

// hydrogen fraction
#ifndef HFRAC 
#define HFRAC 1.
#endif

// helium fraction
#ifndef HEFRAC 
#define HEFRAC 0.
#endif

// metallicity
#ifndef MFRAC 
#define MFRAC (1.-HFRAC-HEFRAC) 
#endif

// weighted average of inverse of atomic weight, default solar
#ifndef A_INV_MEAN 
#define A_INV_MEAN 0.0570490404519 
#endif

// weighted average of nuclear charge squared divided by atomic weight, default solar
#ifndef Z2divA_MEAN 
#define Z2divA_MEAN 5.14445150955 
#endif

// Gas mean molecular weight
#ifndef MU_GAS
#define MU_GAS (1./(0.5 + 1.5*HFRAC + 0.25*HEFRAC + A_INV_MEAN*MFRAC))
#endif

// Ion mean molecular weight
#ifndef MU_I
#define MU_I (1./(HFRAC + 0.25*HEFRAC + A_INV_MEAN*MFRAC))
#endif

// Electron mean molecular weight 
#ifndef MU_E
#define MU_E (2./(1. + HFRAC))
#endif

// Adiabatic index of mixed gas
#ifndef GAMMA
#define GAMMA (5./3.) 
#endif
#define GAMMAM1 (GAMMA-1.) 

// Electron Adiabatic index
#ifndef GAMMAE
#define GAMMAE (5./3.) 
#endif

// Ion Adiabatic index
#ifndef GAMMAI
#define GAMMAI (5./3.) 
#endif

// Black hole mass
#ifndef MASS
#define MASS 1./MSUNCM 
#endif

#ifndef MASSCM
#define MASSCM (MASS*MSUNCM) 
#endif

#ifndef BHSPIN
#define BHSPIN 0.
#endif

/****** Frames/coordinates  ******/

#ifndef MYCOORDS2
#define MYCOORDS2 MYCOORDS
#endif

#ifndef OUTCOORDS
#define OUTCOORDS MYCOORDS
#endif

#ifndef RADCLOSURECOORDS
#define RADCLOSURECOORDS MYCOORDS
#endif

#ifndef SILOCOORDS
#define SILOCOORDS MINKCOORDS
#endif

// Velocity choices
#ifndef VELPRIM
#define VELPRIM VELR
#endif

#ifndef VELPRIMRAD
#define VELPRIMRAD VELR
#endif

#ifndef OUTVEL
#define OUTVEL VEL3
#endif

/********Simulation run ***********/
//maximum number of outputs
#ifndef NOUTSTOP
#define NOUTSTOP 1e50
#endif

//maximum  number of steps
#ifndef NSTEPSTOP
#define NSTEPSTOP 1e50 
#endif

//maximum time
#ifndef TMAX
#define TMAX 1.e50 
#endif

//verbose level
#ifndef VERBOSE0
#define VERBOSE0 0 
#endif

//time stepping
#ifndef TIMESTEPPING
#define TIMESTEPPING RK2IMEX 
#endif

//flux method
#ifndef FLUXMETHOD
#define FLUXMETHOD LAXF_FLUX
#endif

//integration  order
#ifndef INT_ORDER
#define INT_ORDER 1 
#endif

//flux-limiter: Use basic Minmod as the default
#ifndef FLUXLIMITER
#define FLUXLIMITER 0
#endif

//midway between MinMod and MC slope limiter
#ifndef MINMOD_THETA
#define MINMOD_THETA 1.5
#endif

//whether to include metric determinant into the fluxes
#ifndef GDETIN
#define GDETIN 1 
#endif

//radiation source term
#ifndef EXPLICIT_LAB_RAD_SOURCE
#ifndef EXPLICIT_SUBSTEP_RAD_SOURCE
#ifndef IMPLICIT_LAB_RAD_SOURCE
#define IMPLICIT_LAB_RAD_SOURCE
#endif
#endif
#endif

//Radiation closure (no choices beyond M1 anymore)
#ifndef RADCLOSURE
#define RADCLOSURE M1CLOSURE
#endif

//whether to allow reducing implicit_lab to explicit
#ifndef ALLOW_EXPLICIT_RAD_SOURCE
#define ALLOW_EXPLICIT_RAD_SOURCE 0 
#endif

//electron-ion heating
#ifndef ELECTRONIONHEATTYPE
#define ELECTRONIONHEATTYPE ELECTRONIONHEATTYPE_THROUGHUINT
#endif

//number of input arguments in the command line
#ifndef NUM_INPUTARG
#define NUM_INPUTARG 0 
#endif

//shuffle loops
#ifndef SHUFFLELOOPS
#define SHUFFLELOOPS 0
#endif


/*******number of evolved variables********/
//number of nth electron bins
#ifndef NRELBIN
#define NRELBIN 0
#endif

//number of hydro variables
#ifndef EVOLVEELECTRONS
#define NVHD (6)

#else

#ifdef RELELECTRONS
#define NVHD (6+2+NRELBIN)
#else
#define NVHD (6+2)
#endif
#endif

//number of magneto-hydro variables
#ifdef MAGNFIELD
#define NVMHD (NVHD+3)
#else
#define NVMHD (NVHD)
#endif


#ifdef RADIATION

//number of radiative quantities per fluid
#ifndef EVOLVEPHOTONNUMBER
#define NRADVAR 4
#else
#define NRADVAR 5
#endif

//number of total variables
#define NV (NVMHD+NRADVAR)

#else //no RADIATION

//number of total variables
#define NV (NVMHD)
#define NRADVAR 4 //not used
#endif

//number of averaged  quantities
#ifndef NAVGVARS
#define NAVGVARS (180+3*NV) 
#endif

/****** Initialization *********/

//atmosphere initialization parameters
#ifndef RHOATMMIN
#define RHOATMMIN 1.
#endif

#ifndef ERADATMMIN
#define ERADATMMIN 1.
#endif

#ifndef UINTATMMIN
#define UINTATMMIN 1.e-2
#endif

/****** Fixups *****************/

#ifndef NCCORRECTPOLAR
#define NCCORRECTPOLAR 2
#endif

#ifndef ALLOWENTROPYU2P
#define ALLOWENTROPYU2P 1
#endif

#ifndef DOFIXUPS
#define DOFIXUPS 1
#endif

#ifndef DOU2PRADFIXUPS
#define DOU2PRADFIXUPS 0
#endif

#ifndef DOU2PMHDFIXUPS
#define DOU2PMHDFIXUPS 0
#endif

#ifndef DORADIMPFIXUPS
#define DORADIMPFIXUPS 0
#endif


/****** Floors *******/

//min uint over rho 
#ifndef UURHORATIOMIN
#define UURHORATIOMIN 1.e-10
#endif

//uint over rho for u2p_cold
#ifndef UURHORATIOU2PCOLD
#define UURHORATIOU2PCOLD 1.e-10
#endif

//max uint over rho
#ifndef UURHORATIOMAX 
#define UURHORATIOMAX 1.e2
#endif

//min Erad over rho
#ifndef EERHORATIOMIN
#define EERHORATIOMIN 1.e-30
#endif

//max Erad over rho
#ifndef EERHORATIOMAX 
#define EERHORATIOMAX 1.e30
#endif

//min Erad over uint
#ifndef EEUURATIOMIN
#define EEUURATIOMIN 1.e-30
#endif

//max Erad over uint
#ifndef EEUURATIOMAX 
#define EEUURATIOMAX 1.e30
#endif

//min absolute Erad
#ifndef ERADFLOOR
#define ERADFLOOR (10.*SMALL)
#endif

//min B^2 over uint
#ifndef B2UURATIOMIN 
#define B2UURATIOMIN 0.
#endif

//max B^2 over uint
#ifndef B2UURATIOMAX 
#define B2UURATIOMAX 100.
#endif

//max B^2 over Ehat
#ifndef B2EERATIOMAX
#define B2EERATIOMAX 100.
#endif

//min B^2 over rho 
#ifndef B2RHOFLOORFRAME
#define B2RHOFLOORFRAME ZAMOFRAME
#endif

#ifndef B2RHORATIOMIN 
#define B2RHORATIOMIN 0.
#endif

//max B^2 over rho 
#ifndef B2RHORATIOMAX 
#define B2RHORATIOMAX 100.
#endif

//absolute density floor
#ifndef RHOFLOOR
#define RHOFLOOR 1.e-50
#endif

//absolute energy floor
#ifndef UUFLOOR
#define UUFLOOR 1.e-50
#endif

//absolute rad energy floor
#ifndef EEFLOOR
#define EEFLOOR 1.e-50
#endif

//maximum lorentz factor of fluid frame
#ifndef GAMMAMAXHD
#define GAMMAMAXHD 100.
#endif

//maximum lorentz factor or radiation  frame
#ifndef GAMMAMAXRAD
#define GAMMAMAXRAD 100.
#endif

//at what point above which assume gamma^2=1.0
#ifndef GAMMASMALLLIMIT
#define GAMMASMALLLIMIT (1.0-1E-10) 
#endif

//max Te/Tgas
#ifndef TEMPEMAXIMALFRACTION
#define TEMPEMAXIMALFRACTION 1000.
#endif

//max Ti/Tgas
#ifndef TEMPIMAXIMALFRACTION
#define TEMPIMAXIMALFRACTION 1000.
#endif

//minimal ratio of electron to gas energy
#ifndef UEUINTMINRATIO
#define UEUINTMINRATIO 0.01
#endif

//minimal ratio of ion to gas energy
#ifndef UIUINTMINRATIO
#define UIUINTMINRATIO 0.01
#endif

//minimal temperature for electrons
#ifndef TEMPEMINIMAL
#define TEMPEMINIMAL 1.e2
#endif

//minimal ratio of electrons to gas temperature
#ifndef TEMPEMINIMALFRACTION
#define TEMPEMINIMALFRACTION 1.e-5
#endif

//minimal temperature for ions
#ifndef TEMPIMINIMAL
#define TEMPIMINIMAL 1.e2
#endif

//minimal ratio of ions to gas temperature
#ifndef TEMPIMINIMALFRACTION
#define TEMPIMINIMALFRACTION 1.e-5
#endif

//entropy mixing limiter
#ifndef LIMITFACTORINMIXING
#define LIMITFACTORINMIXING 1.
#endif

#ifndef NONRELMHDENTROPYCUT
#define NONRELMHDENTROPYCUT 1.e-7
#endif

//whether to check if the advection operator keeps entropy increasing,
//if not invert with the independently evolved entropy
#ifndef VERIFYENTROPYAFTERADVECTION
#define VERIFYENTROPYAFTERADVECTION 0
#endif

#ifndef CHECKENTROPYAFTEREXPLICITFACTOR
#define CHECKENTROPYAFTEREXPLICITFACTOR 1.e-5
#endif

//maximum fraction of relativistic electron energy to total
#ifndef MAX_RELEL_FRAC_U
#define MAX_RELEL_FRAC_U 0.95
#endif

//maximum fraction of relativistic electron pressure to total
#ifndef MAX_RELEL_FRAC_P
#define MAX_RELEL_FRAC_P 0.95
#endif

//maximum fraction of relativistic electrons to total
#ifndef MAX_RELEL_FRAC_N
#define MAX_RELEL_FRAC_N 0.95
#endif

/******* output ********/

#ifndef OUTOUTPUT
#define OUTOUTPUT 0
#endif

#ifndef SCAOUTPUT
#define SCAOUTPUT 0
#endif

#ifndef RADOUTPUT
#define RADOUTPUT 0
#endif

#ifndef AVGOUTPUT
#define AVGOUTPUT 0
#endif

#ifndef SILOOUTPUT
#define SILOOUTPUT 0
#endif

#ifndef GRIDOUTPUT
#define GRIDOUTPUT 0
#endif

#ifndef SIMOUTPUT
#define SIMOUTPUT 0
#endif

#ifndef ALLSTEPSOUTPUT
#define ALLSTEPSOUTPUT 0
#endif

#ifndef DTOUT2
#define DTOUT2 DTOUT1
#endif

#ifndef DTOUT3
#define DTOUT3 DTOUT1
#endif

#ifndef DTOUT4
#define DTOUT4 DTOUT1
#endif

//number of output scalars
#ifndef NSCALARS
#define NSCALARS 13
#endif

//number of profiles
#ifndef NRADPROFILES
#define NRADPROFILES (68-1)
#endif

#ifndef NTHPROFILES
#define NTHPROFILES 8
#endif

//number of box scalars
#ifndef NBOXSCALARS
#define NBOXSCALARS 44
#endif

#ifndef NBOXCORRSCALARS
#define NBOXCORRSCALARS 10
#endif

#ifndef NBOXVERTSCALARS
#define NBOXVERTSCALARS 26
#endif

#ifndef VARRADIUS
#define VARRADIUS 100.
#endif

#ifndef NVARCUTS
#define NVARCUTS 20
#endif

#ifndef NVARVARSPERCUT
#define NVARVARSPERCUT 2
#endif

#ifndef NVARSCALARS
#define NVARSCALARS (NVARCUTS*NVARVARSPERCUT)
#endif

#ifndef NANARELRADPROFILES
#define NANARELRADPROFILES 4
#endif

#ifndef IMAGETYPE
#define IMAGETYPE "gif"
#endif

/*********U2P implicit solver *****/
#ifndef U2PCONV
#define U2PCONV 1.e-10
#endif

#ifndef U2P_EQS
#define U2P_EQS U2P_EQS_NOBLE
#endif

#ifndef U2P_SOLVER
#define U2P_SOLVER U2P_SOLVER_W
#endif

/********Rad implicit solver *******/
#ifndef FORCEUEQPINIMPLICIT
#define FORCEUEQPINIMPLICIT 1
#endif

#ifndef OPDAMPINIMPLICIT
#define OPDAMPINIMPLICIT 0
#endif

#ifndef OPDAMPMAXLEVELS
#define OPDAMPMAXLEVELS 5.
#endif

#ifndef OPDAMPFACTOR
#define OPDAMPFACTOR 3.
#endif

#ifndef RADIMPLICITTHRESHOLD
#define RADIMPLICITTHRESHOLD 1.e-2
#endif

#ifndef RADIMPLICITDAMPINGFACTOR
#define RADIMPLICITDAMPINGFACTOR 3.
#endif

#ifndef RADIMPLICITMAXNPHCHANGE
#define RADIMPLICITMAXNPHCHANGE 1.e10
#endif

#ifndef RADIMPLICITMAXNECHANGE
#define RADIMPLICITMAXNECHANGE 1.e5
#endif

#ifndef RADIMPLICITMAXTRADCHANGE
#define RADIMPLICITMAXTRADCHANGE 10.
#endif

#ifndef IMPLICITMAXTGASCHANGE
#define IMPLICITMAXTGASCHANGE 10.
#endif

#ifndef RADIMPLICITMAXENCHANGEDOWN
#define RADIMPLICITMAXENCHANGEDOWN 10.
#endif

#ifndef RADIMPLICITMAXENCHANGEUP
#define RADIMPLICITMAXENCHANGEUP 10.
#endif

#ifndef RADIMPLICITMAXTECHANGE
#define RADIMPLICITMAXTECHANGE 10.
#endif

#ifndef MAXRADIMPDAMPING
#define MAXRADIMPDAMPING 1.e-5
#endif

#ifndef RADIMPCONVREL
#define RADIMPCONVREL 1.e-6
#endif

#ifndef RADIMPCONVRELERR
#define RADIMPCONVRELERR (1.e-5)
#endif

#ifndef RADIMPCONVRELENTR
#define RADIMPCONVRELENTR 1.e-6
#endif

#ifndef RADIMPCONVRELENTRERR
#define RADIMPCONVRELENTRERR 1.e-4
#endif

#ifndef RADIMPCONV
#define RADIMPCONV 1.e-10
#endif

#ifndef RADIMPENTRCONV
#define RADIMPENTRCONV 1.e-5
#endif

#ifndef RADIMPEPS
#define RADIMPEPS 1.e-6
#endif

#ifndef RADIMPMAXITER
#define RADIMPMAXITER 50
#endif


/*********metric********/
#ifndef MKSR0
#define MKSR0 0.
#endif

#ifndef MODYFIKUJKRZYSIE
#if (GDETIN==1)
#define MODYFIKUJKRZYSIE 1
#else
#define MODYFIKUJKRZYSIE 0
#endif
#endif

/*********viscosity********/

#ifndef HDVISCOSITY
#define HDVISCOSITY NOVISCOSITY
#endif

#ifndef RADVISCOSITY
#define RADVISCOSITY NOVISCOSITY
#endif

#ifndef MAXRADVISCVEL
#define MAXRADVISCVEL 1.
#endif

/*********synchrotron*********/
#ifdef SYNCHROTRON
// TODO -- do we want synchrotron bridge functions on by default?
//#ifndef NO_SYNCHROTRON_BRIDGE_FUNCTIONS
//#define USE_SYNCHROTRON_BRIDGE_FUNCTIONS
//#endif

#ifndef MAXDIFFTRADS
#define MAXDIFFTRADS 1.e2
#endif
#endif 

#ifndef SYNCHROTRONALLOWEDTRATIO_RELEL
#define SYNCHROTRONALLOWEDTRATIO_RELEL 1.e30
#endif

#ifndef SYNCHROTRONALLOWEDTRATIO
#define SYNCHROTRONALLOWEDTRATIO 1.e30
#endif

/*********comptonization*********/

// comptonization is on by default
#ifdef RADIATION
#ifndef NO_COMPTONIZATION
#define COMPTONIZATION
#define COMPTONIZATIONFLAG
#endif
#endif


/*******battery and dynamo********/
#ifndef DYNAMORADIUS
#define DYNAMORADIUS rISCOBL
#endif

#ifndef BATTERYMULTIPLIER
#define BATTERYMULTIPLIER 1.
#endif

#ifndef COUPLING_MEASURE_THRESHOLD
#define COUPLING_MEASURE_THRESHOLD 0.01
#endif

#ifndef COUPLING_MEASURE_THRESHOLD_SHARPNESS 
#define COUPLING_MEASURE_THRESHOLD_SHARPNESS 0.5
#endif

/*******nonthermal electrons********/
#ifdef RELELECTRONS
#ifndef RELEL_MINMOD_THETA
#define RELEL_MINMOD_THETA 1.8 //midway between MinMod and MC slope limiter
#endif

#ifndef RELEL_INJ_MIN
#define RELEL_INJ_MIN RELGAMMAMIN
#endif

#ifndef RELEL_INJ_MAX
#define RELEL_INJ_MAX RELGAMMAMAX
#endif
#endif


/**********MPI***********/
#ifndef MPIMSGBUFSIZE
#define MPIMSGBUFSIZE 12
#endif

#ifndef MPI4CORNERS //required by MAGNFIELD
#ifdef MAGNFIELD
#define MPI4CORNERS
#endif
#endif

#ifdef MPI4CORNERS
#undef MPIMSGBUFSIZE
#if (TNX>1 && TNY>1 && TNY>1) //3d
#define MPIMSGBUFSIZE 52
#else //2d
#define MPIMSGBUFSIZE 20
#endif
#endif

#ifndef MPI
#ifndef OUTPUTPERCORE
#define OUTPUTPERCORE
#endif
#endif
