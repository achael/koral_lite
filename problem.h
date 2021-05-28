//KORAL - problem.h
//choice of the problem plus some definitions

//available problems:

//1 RADBEAM2D - beam of light 
//2 RADINFALL - radial inflow
//3 DONUT - 2d Polish donut
//4 GEODESICINFALL - geodesic infall with blobs or not
//5 EDDINFALL - infall with flux from inside
//6 RADTUBE- radiative shock tubes as in Farris et al 09
//7 BONDI - like in Fragile's paper
//8 HDTUBE - relativistic shock tube
//9 HDTUBE2D - in 2d
//10 RADPULSE - radiative blob spreading around
//11 RADSHADOW - radiative shadow
//12 RADATM - atmosphere enlighted
//13 DONUTOSC - 2d Polish donut oscillating
//14 RADWAVEBC - 1d linear rad wave imposed on boundary
//15 RADWAVE - 1d linear rad wave with periodic BC
//16 RADPULSE3D - radiative blob spreading around
//17 RADDBLSHADOW - radiative shadow with two beams inclined
//18 ATMSTATIC - hydro atmosphere 
//19 RADBEAM2DKS - beam of light in KS coordinates
//20 ATMKS - radial atmosphere infalling in KS
//21 DONUTKS - 2d Polish donut in KS
//22 DONUTMKS1 - 2d Polish donut in MKS1
//23 ATMMKS1 - radial atmosphere infalling in MKS1
//24 RADBEAMFLAT - beam of light in Cartesian 
//25 RDONUT - 2d radiative Polish donut in KS
//26 RADBEAM2DKSVERT - 2d radiative beam in r,theta plane
//27 RADFLATNESS - flat but with non-zero four-force
//28 BOWSHOCK - bow shock hydro test
//29 RADWALL - flat with wall
//30 RADNT - emission from midplane
//31 FLATDISK - emission from flat disk
//32 CYLBEAM - beam towards the axis in cylindrical
//33 RADDOT - radiating dots
//34 MFPULSE - multi fluid pulse
//35 MFBEAMS - multi fluid colliding beams
//36 MFDOTS - multi fluid radiating dots
//37 MFCYLBEAM - beam towards the axis in cylindrical with multifluids
//40 CYLBEAMCART - similar to discrete CYLBEAM but in cartesian 
//41 FLATDOT - dot which may be bigger than a dot
//42 RVDONUT - radiative and viscous dougnut
//43 RVDONUTIN - radiative and viscous dougnut inflowing
//44 RADNTCYL - emission from midplane in cylindrical
//45 MFRADNTCYL - multi fluid emission from midplane in cylindrical
//46 MFRADNTSPH - multi fluid emission from midplane in spherical
//47 MFDONUT - inflowing donut with radiation/viscosity and multi-fluids
//48 RVDTEST - testing the entropy near the axis
//49 RVRING - 1d donut ring
//50 HDWAVE - hydro wave for testing entropy
//51 HUBBLE - hubble-type test as in the WHAM paper
//52 SPHFLAT - spherical flat to test Christoffels
//53 1DDONUT - 1d donut equatorial plane structure
//54 RVDISK - viscous rad. disk damped from larger radius
//55 G2STATIC - cloud through static atmosphere
//56 EDDSPH - 1d Edd Sphere
//57 LIMOTORUS - torus from Bob's/Akshay's paper
//58 LUKE - beam hitting a net gas flow / the simplest 3d
//59 BBBLOB - rad blobs
//60 JONTESTS - testing jon's fails.dat
//61 MAGTESTS - simple mag tests
//62 MAGTUBES - magnetic tubes
//63 ORSZAG - Orszag-Tang vortex
//64 MAGDONUT - MHD donut with poloidal magnetic fields
//65 MRTORUS - RMHD torus with radiation
//66 KTORUS - RMHD Newtonian torus with radiation
//67 LRTORUS - RMHD limo torus
//68 RMHDWAVE - radiation modified linear magnetosonic waves
//69 INJDISK - magnetized gas damped from outer boundary
//70 M1SHOCK - light beam penetrating absorbing medium, Monreal & Frank 08
//71 MBBBLOB - rad blobs with magn field
//72 BOWSHOCK - bow shock in wind tunnel
//73 VISCZERO - tests of radiative diffusion to compare with ZERO solution
//74 ERICTEST - tests for Eric
//75 LBULLET - bullet through limo torus
//76 ZEBRA - Eric's tidal disks
//77 EMPTY - empty problem
//78 FUNNY - for Ronaldo & Grzesiek & not for Bhupendra
//79 SOFTBALL - Grzegorz's oscillating torus
//80 GENCLOSURE - tests of closures
//81 TWOBEAM - test of VET with two beams
//82 VETTHIN - thin disk emission for VET tests
//83 RADFLATTEST - tests of radiation solvers in flat space
//84 AXISBLOB - blob near 3d polar axis
//85 NSTEST - tests of axisymmetric NS
//86 TFLAT - tests of TFLATCOORDS
//87 INFDISK - gas dumped to circularize
//88 VETPHIBEAMS - beams in phi
//89 RADTORUS - radiating fixed torus
//90 VETSHADOW - vet shadow test
//91 TDEMILIO - tidal disruption with Emilio's input
//92 ADVDIFF - diffusion in scattering moving blob
//93 KATOTORUS - radiative torus initiated like in Kato+2004
//94 TDEMILIO3D - tidal disruption with Emilio's input
//95 XINYI - two half-planes colliding in cartesian
//96 BATTERYTORUS - tests of Contopoulos battery
//97 HUBBLE - Hubble like in WENO paper
//98 HUBBLE2D - Hubble like in WENO paper
//99 THINDISK - starting from a torus, adjusting torus
//100 NSTEST1 - first tests for nuetron stars
//101 NSDISK - NS with accretion disk
//102 ADAFCRIT - critical ADAFs
//103 KATOLANL - KATOTORUS that uses opacity table
//104 NSHOVER - Hovering atmosphere
//105 COULOMBETC - uniform gas interacting with radiation
//106 SEDRIVENTURB - driven turbulence to test electron heating
//107 RELELTEST - test relativistic electron evolution
//108 SYNCTEST - tests of synchrotron emission
//109 SHEARINGBOX - simple shearing box
//110 SHEARMRI - development of MRI test
//111 EDDOSC - Oscillation of Levitating Atmospheres
//112 TDEMILIOSPIN
//113 ULXBUBBLE - Bubble expansion
//114 RELELEXPAND - cooling of relativistic electrons in Hubble expansion
//115 SHOCKELECTRONTEST - 1D (Noh) Shock test from Ressler+15
//116 MHDSHOCKS - relativistic MHD shocks
//117 BATTERYTEST - Antonio's test of radiation battery
//118 SEDDRIVENTURBRELEL - Driven turbulence with rel. electrons
//119 ADAFCRITRELEL - ADAFCRIT + non thermal electrons
//120 RADSHOCKS - to measure convergence of radiative shocks
//121 AGNBUBBLE - Bubble expansion
//122 TURBKNL - 118 modified for KNL testing
//123 HIGHZBHS - super-Eddington accretion on high z BHs
//124 NSTORUS -  super-Eddington accretion on a hard surface
//125 M87 - importance of radiation for M87 simulations
//126 SHOCKENTROPY - 1D hydro shock to study entropy generation
//127 TDESURFACE - shell falling on a BH-surface
//128 TEST_CONV_VEL - testing the conv_vel routines
//129 TEST_WAVESPEEDS - testing the wavespeeds routines
//130 TEST_EXPLICIT_IMPLICIT - testing one step of the explicit-implicit cycle
//131 TEST_EXPLICIT_IMPLICIT_FULL - testing multiple steps of explicit-implicit
//132 MONOPOLE_2D - copy of HARMPI 2D monopole test problem
//133 TEST_EXPLICIT_IMPLICIT_2TEMP_RELEL - testing one step of the explicit-implicit cycle with adafcritrelel
//134 FISHMONC - Fishbone & Moncrief torus model: ApJ, 207, 962, 1976
//135 WINDACCRETION
//136 CONTACTELECTRONTEST - Contact discontinuity electron heating test
//137 KATOTORUS_TILTED - radiative torus initiated like in Kato+2004 with a tilt from BH spin axis
//138 JETCOORDS -- test jet coordinate system
//139 MADCC -- MAD code comparison
//140 RADSURVEY -- radiative parameter survey
//141 KEPINF -- INFDISK modified for injecting keplerian 
//142 PARTIALTDE -- INFDISK modified for partial TDE binding energy distribution non-constant
//143 GPUTEST -- disk for gpu hackathon test problem

#define PROBLEM 143

#if(PROBLEM==143)

#define PR_DEFINE "PROBLEMS/GPUTEST/define.h"
#define PR_BC "PROBLEMS/GPUTEST/bc.c"
#define PR_INIT "PROBLEMS/GPUTEST/init.c"
#define PR_POSTINIT "PROBLEMS/GPUTEST/postinit.c"
#define PR_TOOLS "PROBLEMS/RADSURVEY/tools.c"

#endif

#if(PROBLEM==142)

#define PR_DEFINE "PROBLEMS/PARTIALTDE/define.h"
#define PR_BC "PROBLEMS/PARTIALTDE/bc.c"
#define PR_BC_SPECIAL "PROBLEMS/PARTIALTDE/bc_special.c"
#define PR_BC_SPECIAL_LOOP "PROBLEMS/PARTIALTDE/loop_alloc_special.c"
#define PR_INIT "PROBLEMS/PARTIALTDE/init.c"
//#define PR_KAPPA "PROBLEMS/KEPINF/kappa.c"
#define PR_KAPPAES "PROBLEMS/PARTIALTDE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/PARTIALTDE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/PARTIALTDE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/PARTIALTDE/dump.c"
#define PR_TOOLS "PROBLEMS/PARTIALTDE/tools.c"

#endif

#if(PROBLEM==141)

#define PR_DEFINE "PROBLEMS/KEPINF/define.h"
#define PR_BC "PROBLEMS/KEPINF/bc.c"
#define PR_BC_SPECIAL "PROBLEMS/KEPINF/bc_special.c"
#define PR_BC_SPECIAL_LOOP "PROBLEMS/KEPINF/loop_alloc_special.c"
#define PR_INIT "PROBLEMS/KEPINF/init.c"
//#define PR_KAPPA "PROBLEMS/KEPINF/kappa.c"
#define PR_KAPPAES "PROBLEMS/KEPINF/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/KEPINF/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/KEPINF/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/KEPINF/dump.c"
#define PR_TOOLS "PROBLEMS/KEPINF/tools.c"

#endif

#if(PROBLEM==140)

#define PR_DEFINE "PROBLEMS/RADSURVEY/define.h"
#define PR_BC "PROBLEMS/RADSURVEY/bc.c"
#define PR_INIT "PROBLEMS/RADSURVEY/init.c"
#define PR_KAPPAES "PROBLEMS/RADSURVEY/kappaes.c"
#define PR_TOOLS "PROBLEMS/RADSURVEY/tools.c"
#define PR_POSTINIT "PROBLEMS/RADSURVEY/postinit.c"

#endif

#if(PROBLEM==139)

#define PR_DEFINE "PROBLEMS/MADCC/define.h"
#define PR_BC "PROBLEMS/MADCC/bc.c"
#define PR_INIT "PROBLEMS/MADCC/init.c"
#define PR_KAPPAES "PROBLEMS/MADCC/kappaes.c"
#define PR_TOOLS "PROBLEMS/MADCC/tools.c"
#define PR_POSTINIT "PROBLEMS/MADCC/postinit.c"

#endif

#if(PROBLEM==138)

#define PR_DEFINE "PROBLEMS/JETCOORDS/define.h"
#define PR_BC "PROBLEMS/JETCOORDS/bc.c"
#define PR_INIT "PROBLEMS/JETCOORDS/init.c"
#define PR_KAPPAES "PROBLEMS/JETCOORDS/kappaes.c"
#define PR_TOOLS "PROBLEMS/JETCOORDS/tools.c"
#define PR_POSTINIT "PROBLEMS/JETCOORDS/postinit.c"

#endif

#if(PROBLEM==137)

#define PR_DEFINE "PROBLEMS/KATOTORUS_TILTED/define.h"
#define PR_BC "PROBLEMS/KATOTORUS_TILTED/bc.c"
#define PR_INIT "PROBLEMS/KATOTORUS_TILTED/init.c"
#define PR_POSTINIT "PROBLEMS/KATOTORUS_TILTED/postinit.c"
//#define PR_KAPPA "PROBLEMS/KATOTORUS_TILTED/kappa.c"
#define PR_KAPPAES "PROBLEMS/KATOTORUS_TILTED/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/KATOTORUS_TILTED/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/KATOTORUS_TILTED/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/KATOTORUS_TILTED/dump.c"
#define PR_TOOLS "PROBLEMS/KATOTORUS_TILTED/tools.c"

#endif

#if(PROBLEM==136)
#define PR_DEFINE "PROBLEMS/CONTACTELECTRONTEST/define.h"
#define PR_BC "PROBLEMS/CONTACTELECTRONTEST/bc.c"
#define PR_INIT "PROBLEMS/CONTACTELECTRONTEST/init.c"
#define PR_KAPPAES "PROBLEMS/CONTACTELECTRONTEST/kappaes.c"
#define PR_DUMP "PROBLEMS/CONTACTELECTRONTEST/dump.c"
#define PR_TOOLS "PROBLEMS/CONTACTELECTRONTEST/tools.c"
#define PR_POSTINIT "PROBLEMS/CONTACTELECTRONTEST/postinit.c"
#endif

#if(PROBLEM==135)

#define PR_DEFINE "PROBLEMS/WINDACCRETION/define.h"
#define PR_BC "PROBLEMS/WINDACCRETION/bc.c"
#define PR_INIT "PROBLEMS/WINDACCRETION/init.c"
//#define PR_KAPPA "PROBLEMS/WINDACCRETION/kappa.c"
#define PR_KAPPAES "PROBLEMS/WINDACCRETION/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/WINDACCRETION/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/WINDACCRETION/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/WINDACCRETION/dump.c"
#define PR_TOOLS "PROBLEMS/WINDACCRETION/tools.c"
#endif

#if(PROBLEM==134)
#define PR_DEFINE "PROBLEMS/FISHMONC/define.h"
#define PR_BC "PROBLEMS/FISHMONC/bc.c"
#define PR_INIT "PROBLEMS/FISHMONC/init.c"
//#define PR_KAPPA "PROBLEMS/FISHMONC/kappa.c"
#define PR_KAPPAES "PROBLEMS/FISHMONC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/FISHMONC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/FISHMONC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/FISHMONC/dump.c"
#define PR_TOOLS "PROBLEMS/FISHMONC/tools.c"
#define PR_POSTINIT "PROBLEMS/FISHMONC/postinit.c"
#define PR_WRITE_OUTPUT "PROBLEMS/FISHMONC/write_output.c"
#endif

#if(PROBLEM==133)

#define PR_DEFINE "PROBLEMS/ADAFCRITRELEL/define_explicit_implicit.h"
#define PR_BC "PROBLEMS/ADAFCRITRELEL/bc.c"
#define PR_INIT "PROBLEMS/ADAFCRITRELEL/init.c"
//#define PR_KAPPA "PROBLEMS/ADAFCRITRELEL/kappa.c"
#define PR_KAPPAES "PROBLEMS/ADAFCRITRELEL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ADAFCRITRELEL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ADAFCRITRELEL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ADAFCRITRELEL/dump.c"
#define PR_TOOLS "PROBLEMS/ADAFCRITRELEL/tools.c"
#define PR_POSTINIT "PROBLEMS/ADAFCRITRELEL/postinit.c"
#endif

#if(PROBLEM==132)

#define PR_DEFINE "PROBLEMS/MONOPOLE_2D/define.h"
#define PR_BC "PROBLEMS/MONOPOLE_2D/bc.c"
#define PR_INIT "PROBLEMS/MONOPOLE_2D/init.c"
//#define PR_KAPPA "PROBLEMS/MONOPOLE_2D/kappa.c"
//#define PR_KAPPAES "PROBLEMS/MONOPOLE_2D/kappaes.c"
//#define PR_OUT2GIF_2D "PROBLEMS/MONOPOLE_2D/out2gif_2d.c"
//#define PR_OUT2GIF_1D "PROBLEMS/MONOPOLE_2D/out2gif_1d.c"
//#define PR_DUMP "PROBLEMS/MONOPOLE_2D/dump.c"
//#define PR_TOOLS "PROBLEMS/MONOPOLE_2D/tools.c"
//#define PR_POSTINIT "PROBLEMS/MAGDONUT/postinit.c"
#endif

#if(PROBLEM==131)

#define PR_DEFINE "PROBLEMS/HIGHZBHS/define_explicit_implicit.h"
#define PR_BC "PROBLEMS/HIGHZBHS/bc.c"
#define PR_INIT "PROBLEMS/HIGHZBHS/init.c"
//#define PR_KAPPA "PROBLEMS/HIGHZBHS/kappa.c"
#define PR_KAPPAES "PROBLEMS/HIGHZBHS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HIGHZBHS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HIGHZBHS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HIGHZBHS/dump.c"
#define PR_TOOLS "PROBLEMS/HIGHZBHS/tools.c"
#define PR_POSTINIT "PROBLEMS/HIGHZBHS/postinit.c"

#endif

#if(PROBLEM==130)

#define PR_DEFINE "PROBLEMS/HIGHZBHS/define_explicit_implicit.h"
#define PR_BC "PROBLEMS/HIGHZBHS/bc.c"
#define PR_INIT "PROBLEMS/HIGHZBHS/init.c"
//#define PR_KAPPA "PROBLEMS/HIGHZBHS/kappa.c"
#define PR_KAPPAES "PROBLEMS/HIGHZBHS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HIGHZBHS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HIGHZBHS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HIGHZBHS/dump.c"
#define PR_TOOLS "PROBLEMS/HIGHZBHS/tools.c"
#define PR_POSTINIT "PROBLEMS/HIGHZBHS/postinit.c"

#endif

#if(PROBLEM==129)

#define PR_DEFINE "PROBLEMS/HIGHZBHS/define_wavespeeds.h"
#define PR_BC "PROBLEMS/HIGHZBHS/bc.c"
#define PR_INIT "PROBLEMS/HIGHZBHS/init.c"
//#define PR_KAPPA "PROBLEMS/HIGHZBHS/kappa.c"
#define PR_KAPPAES "PROBLEMS/HIGHZBHS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HIGHZBHS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HIGHZBHS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HIGHZBHS/dump.c"
#define PR_TOOLS "PROBLEMS/HIGHZBHS/tools.c"
#define PR_POSTINIT "PROBLEMS/HIGHZBHS/postinit.c"

#endif

#if(PROBLEM==128)

#define PR_DEFINE "PROBLEMS/HIGHZBHS/define_conv_vel.h"
#define PR_BC "PROBLEMS/HIGHZBHS/bc.c"
#define PR_INIT "PROBLEMS/HIGHZBHS/init.c"
//#define PR_KAPPA "PROBLEMS/HIGHZBHS/kappa.c"
#define PR_KAPPAES "PROBLEMS/HIGHZBHS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HIGHZBHS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HIGHZBHS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HIGHZBHS/dump.c"
#define PR_TOOLS "PROBLEMS/HIGHZBHS/tools.c"
#define PR_POSTINIT "PROBLEMS/HIGHZBHS/postinit.c"

#endif

#if(PROBLEM==127)

#define PR_DEFINE "PROBLEMS/TDESURFACE/define.h"
#define PR_BC "PROBLEMS/TDESURFACE/bc.c"
#define PR_INIT "PROBLEMS/TDESURFACE/init.c"
//#define PR_KAPPA "PROBLEMS/TDESURFACE/kappa.c"
#define PR_KAPPAES "PROBLEMS/TDESURFACE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/TDESURFACE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/TDESURFACE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/TDESURFACE/dump.c"
#define PR_TOOLS "PROBLEMS/TDESURFACE/tools.c"
#define PR_POSTINIT "PROBLEMS/TDESURFACE/postinit.c"

#endif

#if(PROBLEM==126)

#define PR_DEFINE "PROBLEMS/SHOCKENTROPY/define.h"
#define PR_BC "PROBLEMS/SHOCKENTROPY/bc.c"
#define PR_INIT "PROBLEMS/SHOCKENTROPY/init.c"
#define PR_KAPPA "PROBLEMS/SHOCKENTROPY/kappa.c"
#define PR_KAPPAES "PROBLEMS/SHOCKENTROPY/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/SHOCKENTROPY/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/SHOCKENTROPY/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/SHOCKENTROPY/dump.c"
#define PR_TOOLS "PROBLEMS/SHOCKENTROPY/tools.c"
#define PR_POSTINIT "PROBLEMS/SHOCKENTROPY/postinit.c"

#endif

#if(PROBLEM==125)

#define PR_DEFINE "PROBLEMS/M87/define.h"
#define PR_BC "PROBLEMS/M87/bc.c"
#define PR_INIT "PROBLEMS/M87/init.c"
//#define PR_KAPPA "PROBLEMS/M87/kappa.c"
#define PR_KAPPAES "PROBLEMS/M87/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/M87/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/M87/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/M87/dump.c"
#define PR_TOOLS "PROBLEMS/M87/tools.c"
#define PR_POSTINIT "PROBLEMS/M87/postinit.c"

#endif

#if(PROBLEM==124)

#define PR_DEFINE "PROBLEMS/NSTORUS/define.h"
#define PR_BC "PROBLEMS/NSTORUS/bc.c"
#define PR_INIT "PROBLEMS/NSTORUS/init.c"
//#define PR_KAPPA "PROBLEMS/NSTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/NSTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/NSTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/NSTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/NSTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/NSTORUS/tools.c"
#define PR_POSTINIT "PROBLEMS/NSTORUS/postinit.c"

#endif

#if(PROBLEM==123)

#define PR_DEFINE "PROBLEMS/HIGHZBHS/define.h"
#define PR_BC "PROBLEMS/HIGHZBHS/bc.c"
#define PR_INIT "PROBLEMS/HIGHZBHS/init.c"
//#define PR_KAPPA "PROBLEMS/HIGHZBHS/kappa.c"
#define PR_KAPPAES "PROBLEMS/HIGHZBHS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HIGHZBHS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HIGHZBHS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HIGHZBHS/dump.c"
#define PR_TOOLS "PROBLEMS/HIGHZBHS/tools.c"
#define PR_POSTINIT "PROBLEMS/HIGHZBHS/postinit.c"

#endif

#if(PROBLEM==122)
#define PR_DEFINE "PROBLEMS/TURBKNL/define.h"
#define PR_BC "PROBLEMS/TURBKNL/bc.c"
#define PR_INIT "PROBLEMS/TURBKNL/init.c"
//#define PR_KAPPA "PROBLEMS/TURBKNL/kappa.c"
#define PR_KAPPAES "PROBLEMS/TURBKNL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/TURBKNL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/TURBKNL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/TURBKNL/dump.c"
#define PR_TOOLS "PROBLEMS/TURBKNL/tools.c"
#define PR_POSTINIT "PROBLEMS/TURBKNL/postinit.c"
#define PR_FINGER "PROBLEMS/TURBKNL/finger.c"
#endif

#if(PROBLEM==121)
#define PR_DEFINE "PROBLEMS/AGNBUBBLE/define.h"
#define PR_BC "PROBLEMS/AGNBUBBLE/bc.c"
#define PR_INIT "PROBLEMS/AGNBUBBLE/init.c"
#define PR_POSTINIT "PROBLEMS/AGNBUBBLE/postinit.c"
#define PR_KAPPA "PROBLEMS/AGNBUBBLE/kappa.c"
#define PR_KAPPAES "PROBLEMS/AGNBUBBLE/kappaes.c"
#define PR_TOOLS "PROBLEMS/AGNBUBBLE/tools.c"
#endif


#if(PROBLEM==120)
#define PR_DEFINE "PROBLEMS/RADSHOCKS/define.h"
#define PR_BC "PROBLEMS/RADSHOCKS/bc.c"
#define PR_INIT "PROBLEMS/RADSHOCKS/init.c"
#define PR_POSTINIT "PROBLEMS/RADSHOCKS/postinit.c"
#define PR_KAPPA "PROBLEMS/RADSHOCKS/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADSHOCKS/kappaes.c"
#define PR_TOOLS "PROBLEMS/RADSHOCKS/tools.c"
#endif

#if(PROBLEM==119)

#define PR_DEFINE "PROBLEMS/ADAFCRITRELEL/define.h"
#define PR_BC "PROBLEMS/ADAFCRITRELEL/bc.c"
#define PR_INIT "PROBLEMS/ADAFCRITRELEL/init.c"
//#define PR_KAPPA "PROBLEMS/ADAFCRITRELEL/kappa.c"
#define PR_KAPPAES "PROBLEMS/ADAFCRITRELEL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ADAFCRITRELEL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ADAFCRITRELEL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ADAFCRITRELEL/dump.c"
#define PR_TOOLS "PROBLEMS/ADAFCRITRELEL/tools.c"
#define PR_POSTINIT "PROBLEMS/ADAFCRITRELEL/postinit.c"

#endif

#if(PROBLEM==118)
#define PR_DEFINE "PROBLEMS/SEDRIVENTURBRELEL/define.h"
#define PR_BC "PROBLEMS/SEDRIVENTURBRELEL/bc.c"
#define PR_INIT "PROBLEMS/SEDRIVENTURBRELEL/init.c"
//#define PR_KAPPA "PROBLEMS/SEDRIVENTURBRELEL/kappa.c"
#define PR_KAPPAES "PROBLEMS/SEDRIVENTURBRELEL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/SEDRIVENTURBRELEL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/SEDRIVENTURBRELEL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/SEDRIVENTURBRELEL/dump.c"
#define PR_TOOLS "PROBLEMS/SEDRIVENTURBRELEL/tools.c"
#define PR_POSTINIT "PROBLEMS/SEDRIVENTURBRELEL/postinit.c"
#define PR_FINGER "PROBLEMS/SEDRIVENTURBRELEL/finger.c"
#endif

#if(PROBLEM==117)
#define PR_DEFINE "PROBLEMS/BATTERYTEST/define.h"
#define PR_BC "PROBLEMS/BATTERYTEST/bc.c"
#define PR_INIT "PROBLEMS/BATTERYTEST/init.c"
//#define PR_KAPPA "PROBLEMS/BATTERYTEST/kappa.c"
#define PR_KAPPAES "PROBLEMS/BATTERYTEST/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BATTERYTEST/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BATTERYTEST/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BATTERYTEST/dump.c"
#define PR_FINGER "PROBLEMS/BATTERYTEST/finger.c"
#define PR_TOOLS "PROBLEMS/BATTERYTEST/tools.c"
#define PR_PREPINIT "PROBLEMS/BATTERYTEST/prepinit.c"
#define PR_POSTINIT "PROBLEMS/BATTERYTEST/postinit.c"
#endif

#if(PROBLEM==116)
#define PR_DEFINE "PROBLEMS/MHDSHOCKS/define.h"
#define PR_BC "PROBLEMS/MHDSHOCKS/bc.c"
#define PR_INIT "PROBLEMS/MHDSHOCKS/init.c"
#define PR_POSTINIT "PROBLEMS/MHDSHOCKS/postinit.c"
#define PR_KAPPA "PROBLEMS/MHDSHOCKS/kappa.c"
#define PR_KAPPAES "PROBLEMS/MHDSHOCKS/kappaes.c"
#define PR_TOOLS "PROBLEMS/MHDSHOCKS/tools.c"
#endif

#if(PROBLEM==115)
#define PR_DEFINE "PROBLEMS/SHOCKELECTRONTEST/define.h"
#define PR_BC "PROBLEMS/SHOCKELECTRONTEST/bc.c"
#define PR_INIT "PROBLEMS/SHOCKELECTRONTEST/init.c"
//#define PR_KAPPA "PROBLEMS/SHOCKELECTRONTEST/kappa.c"
//#define PR_KAPPAES "PROBLEMS/SHOCKELECTRONTEST/kappaes.c"
//#define PR_DUMP "PROBLEMS/SHOCKELECTRONTEST/dump.c"
//#define PR_TOOLS "PROBLEMS/SHOCKELECTRONTEST/tools.c"
//#define PR_POSTINIT "PROBLEMS/SHOCKELECTRONTEST/postinit.c"
#endif

#if(PROBLEM==114)
#define PR_DEFINE "PROBLEMS/RELELEXPAND/define.h"
#define PR_BC "PROBLEMS/RELELEXPAND/bc.c"
#define PR_INIT "PROBLEMS/RELELEXPAND/init.c"
//#define PR_KAPPA "PROBLEMS/RELELEXPAND/kappa.c"
#define PR_KAPPAES "PROBLEMS/RELELEXPAND/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RELELEXPAND/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RELELEXPAND/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RELELEXPAND/dump.c"
#define PR_TOOLS "PROBLEMS/RELELEXPAND/tools.c"
#define PR_POSTINIT "PROBLEMS/RELELEXPAND/postinit.c"
#endif

#if(PROBLEM==113)
#define PR_DEFINE "PROBLEMS/ULXBUBBLE/define.h"
#define PR_BC "PROBLEMS/ULXBUBBLE/bc.c"
#define PR_INIT "PROBLEMS/ULXBUBBLE/init.c"
#define PR_POSTINIT "PROBLEMS/ULXBUBBLE/postinit.c"
#define PR_KAPPA "PROBLEMS/ULXBUBBLE/kappa.c"
#define PR_KAPPAES "PROBLEMS/ULXBUBBLE/kappaes.c"
#define PR_TOOLS "PROBLEMS/ULXBUBBLE/tools.c"
#endif

#if(PROBLEM==112)
#define PR_PREPINIT "PROBLEMS/TDEMILIOSPIN/prepinit.c"
#define PR_DEFINE "PROBLEMS/TDEMILIOSPIN/define.h"
#define PR_BC "PROBLEMS/TDEMILIOSPIN/bc.c"
#define PR_INIT "PROBLEMS/TDEMILIOSPIN/init.c"
#define PR_POSTINIT "PROBLEMS/TDEMILIOSPIN/postinit.c"
#define PR_KAPPA "PROBLEMS/TDEMILIOSPIN/kappa.c"
#define PR_KAPPAES "PROBLEMS/TDEMILIOSPIN/kappaes.c"
#define PR_TOOLS "PROBLEMS/TDEMILIOSPIN/tools.c"
#endif


#if(PROBLEM==111)
#define PR_DEFINE "PROBLEMS/EDDOSC/define.h"
#define PR_BC "PROBLEMS/EDDOSC/bc.c"
#define PR_INIT "PROBLEMS/EDDOSC/init.c"
#define PR_KAPPA "PROBLEMS/EDDOSC/kappa.c"
#define PR_KAPPAES "PROBLEMS/EDDOSC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/EDDOSC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/EDDOSC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/EDDOSC/dump.c"
#define PR_TOOLS "PROBLEMS/EDDOSC/tools.c"
#define PR_POSTINIT "PROBLEMS/EDDOSC/postinit.c"

#endif


#if(PROBLEM==110)
#define PR_DEFINE "PROBLEMS/SHEARMRI/define.h"
#define PR_BC "PROBLEMS/SHEARMRI/bc.c"
#define PR_INIT "PROBLEMS/SHEARMRI/init.c"
#define PR_POSTINIT "PROBLEMS/SHEARMRI/postinit.c"
//#define PR_KAPPA "PROBLEMS/SHEARMRI/kappa.c"
#define PR_KAPPAES "PROBLEMS/SHEARMRI/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/SHEARMRI/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/SHEARMRI/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/SHEARMRI/dump.c"
#define PR_TOOLS "PROBLEMS/SHEARMRI/tools.c"
#define PR_POSTINIT "PROBLEMS/SHEARMRI/postinit.c"

#endif



#if(PROBLEM==109)
#define PR_DEFINE "PROBLEMS/SHEARINGBOX/define.h"
#define PR_BC "PROBLEMS/SHEARINGBOX/bc.c"
#define PR_INIT "PROBLEMS/SHEARINGBOX/init.c"
//#define PR_KAPPA "PROBLEMS/SHEARINGBOX/kappa.c"
#define PR_KAPPAES "PROBLEMS/SHEARINGBOX/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/SHEARINGBOX/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/SHEARINGBOX/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/SHEARINGBOX/dump.c"
#define PR_TOOLS "PROBLEMS/SHEARINGBOX/tools.c"
#define PR_POSTINIT "PROBLEMS/SHEARINGBOX/postinit.c"

#endif



#if(PROBLEM==108)
#define PR_DEFINE "PROBLEMS/SYNCTEST/define.h"
#define PR_BC "PROBLEMS/SYNCTEST/bc.c"
#define PR_INIT "PROBLEMS/SYNCTEST/init.c"
//#define PR_KAPPA "PROBLEMS/SYNCTEST/kappa.c"
#define PR_KAPPAES "PROBLEMS/SYNCTEST/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/SYNCTEST/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/SYNCTEST/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/SYNCTEST/dump.c"
#define PR_TOOLS "PROBLEMS/SYNCTEST/tools.c"
#define PR_FINGER "PROBLEMS/SYNCTEST/finger.c"

#endif

#if(PROBLEM==107)
#define PR_DEFINE "PROBLEMS/RELELTEST/define.h"
#define PR_BC "PROBLEMS/RELELTEST/bc.c"
#define PR_INIT "PROBLEMS/RELELTEST/init.c"
//#define PR_KAPPA "PROBLEMS/RELELTEST/kappa.c"
#define PR_KAPPAES "PROBLEMS/RELELTEST/kappaes.c"
//#define PR_OUT2GIF_2D "PROBLEMS/RELELTEST/out2gif_2d.c"
//#define PR_OUT2GIF_1D "PROBLEMS/RELELTEST/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RELELTEST/dump.c"
#define PR_TOOLS "PROBLEMS/RELELTEST/tools.c"
#define PR_POSTINIT "PROBLEMS/RELELTEST/postinit.c"
#endif


#if(PROBLEM==106)

#define PR_DEFINE "PROBLEMS/SEDRIVENTURB/define.h"
#define PR_BC "PROBLEMS/SEDRIVENTURB/bc.c"
#define PR_INIT "PROBLEMS/SEDRIVENTURB/init.c"
//#define PR_KAPPA "PROBLEMS/SEDRIVENTURB/kappa.c"
#define PR_KAPPAES "PROBLEMS/SEDRIVENTURB/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/SEDRIVENTURB/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/SEDRIVENTURB/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/SEDRIVENTURB/dump.c"
#define PR_TOOLS "PROBLEMS/SEDRIVENTURB/tools.c"
#define PR_POSTINIT "PROBLEMS/SEDRIVENTURB/postinit.c"
#define PR_FINGER "PROBLEMS/SEDRIVENTURB/finger.c"


#endif

#if(PROBLEM==105)

#define PR_DEFINE "PROBLEMS/COULOMBETC/define.h"
#define PR_BC "PROBLEMS/COULOMBETC/bc.c"
#define PR_INIT "PROBLEMS/COULOMBETC/init.c"
//#define PR_KAPPA "PROBLEMS/COULOMBETC/kappa.c"
#define PR_KAPPAES "PROBLEMS/COULOMBETC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/COULOMBETC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/COULOMBETC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/COULOMBETC/dump.c"
#define PR_TOOLS "PROBLEMS/COULOMBETC/tools.c"
#define PR_POSTINIT "PROBLEMS/COULOMBETC/postinit.c"
#endif

#if(PROBLEM==104)

#define PR_DEFINE "PROBLEMS/NSHOVER/define.h"
#define PR_BC "PROBLEMS/NSHOVER/bc.c"
#define PR_INIT "PROBLEMS/NSHOVER/init.c"
#define PR_KAPPA "PROBLEMS/NSHOVER/kappa.c"
#define PR_KAPPAES "PROBLEMS/NSHOVER/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/NSHOVER/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/NSHOVER/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/NSHOVER/dump.c"
#define PR_TOOLS "PROBLEMS/NSHOVER/tools.c"
#define PR_POSTINIT "PROBLEMS/NSHOVER/postinit.c"
#endif

#if(PROBLEM==103)
#define PR_DEFINE "PROBLEMS/KATOLANL/define.h"
#define PR_BC "PROBLEMS/KATOLANL/bc.c"
#define PR_INIT "PROBLEMS/KATOLANL/init.c"
#define PR_KAPPA "PROBLEMS/KATOLANL/kappa.c"
#define PR_KAPPAES "PROBLEMS/KATOLANL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/KATOLANL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/KATOLANL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/KATOLANL/dump.c"
#define PR_TOOLS "PROBLEMS/KATOLANL/tools.c"
#define PR_POSTINIT "PROBLEMS/KATOLANL/postinit.c"
#endif

#if(PROBLEM==102)

#define PR_DEFINE "PROBLEMS/ADAFCRIT/define.h"
#define PR_BC "PROBLEMS/ADAFCRIT/bc.c"
#define PR_INIT "PROBLEMS/ADAFCRIT/init.c"
//#define PR_KAPPA "PROBLEMS/ADAFCRIT/kappa.c"
#define PR_KAPPAES "PROBLEMS/ADAFCRIT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ADAFCRIT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ADAFCRIT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ADAFCRIT/dump.c"
#define PR_TOOLS "PROBLEMS/ADAFCRIT/tools.c"
#define PR_POSTINIT "PROBLEMS/ADAFCRIT/postinit.c"

#endif

#if(PROBLEM==101)

#define PR_DEFINE "PROBLEMS/NSDISK/define.h"
#define PR_BC "PROBLEMS/NSDISK/bc.c"
#define PR_INIT "PROBLEMS/NSDISK/init.c"
#define PR_KAPPA "PROBLEMS/NSDISK/kappa.c"
#define PR_KAPPAES "PROBLEMS/NSDISK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/NSDISK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/NSDISK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/NSDISK/dump.c"
#define PR_TOOLS "PROBLEMS/NSDISK/tools.c"
#define PR_POSTINIT "PROBLEMS/NSDISK/postinit.c"
#endif

#if(PROBLEM==100)

#define PR_DEFINE "PROBLEMS/NSTEST1/define.h"
#define PR_BC "PROBLEMS/NSTEST1/bc.c"
#define PR_INIT "PROBLEMS/NSTEST1/init.c"
#define PR_KAPPA "PROBLEMS/NSTEST1/kappa.c"
#define PR_KAPPAES "PROBLEMS/NSTEST1/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/NSTEST1/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/NSTEST1/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/NSTEST1/dump.c"
#define PR_TOOLS "PROBLEMS/NSTEST1/tools.c"
#define PR_POSTINIT "PROBLEMS/NSTEST1/postinit.c"

#endif

#if(PROBLEM==99)

#define PR_DEFINE "PROBLEMS/THINDISK/define.h"
#define PR_BC "PROBLEMS/THINDISK/bc.c"
#define PR_INIT "PROBLEMS/THINDISK/init.c"
#define PR_KAPPA "PROBLEMS/THINDISK/kappa.c"
#define PR_KAPPAES "PROBLEMS/THINDISK/kappaes.c"
#define PR_DUMP "PROBLEMS/THINDISK/dump.c"
#define PR_TOOLS "PROBLEMS/THINDISK/tools.c"
#define PR_POSTINIT "PROBLEMS/THINDISK/postinit.c"
#define PR_PREPINIT "PROBLEMS/THINDISK/prepinit.c"
#define PR_FINGER "PROBLEMS/THINDISK/finger.c"

#endif


#if(PROBLEM==98)
#define PR_DEFINE "PROBLEMS/HUBBLE2D/define.h" //most-important, definitions, choices
#define PR_PREPINIT "PROBLEMS/HUBBLE2D/prepinit.c" //auxiliary, if you need to precalculate something
#define PR_INIT "PROBLEMS/HUBBLE2D/init.c" //initial state
#define PR_POSTINIT "PROBLEMS/HUBBLE2D/postinit.c" //auxilary, if you need to e.g. normalize something (magn.field?)
#define PR_BC "PROBLEMS/HUBBLE2D/bc.c" //boundary conditions (if not periodic)
#define PR_KAPPA "PROBLEMS/HUBBLE2D/kappa.c" //absorption opacity
#define PR_KAPPAES "PROBLEMS/HUBBLE2D/kappaes.c" //scattering opacity
#define PR_TOOLS "PROBLEMS/HUBBLE2D/tools.c" //auxiliary
#endif

#if(PROBLEM==97)
#define PR_DEFINE "PROBLEMS/HUBBLE/define.h" //most-important, definitions, choices
#define PR_PREPINIT "PROBLEMS/HUBBLE/prepinit.c" //auxiliary, if you need to precalculate something
#define PR_INIT "PROBLEMS/HUBBLE/init.c" //initial state
#define PR_POSTINIT "PROBLEMS/HUBBLE/postinit.c" //auxilary, if you need to e.g. normalize something (magn.field?)
#define PR_BC "PROBLEMS/HUBBLE/bc.c" //boundary conditions (if not periodic)
#define PR_KAPPA "PROBLEMS/HUBBLE/kappa.c" //absorption opacity
#define PR_KAPPAES "PROBLEMS/HUBBLE/kappaes.c" //scattering opacity
#define PR_TOOLS "PROBLEMS/HUBBLE/tools.c" //auxiliary
#endif


#if(PROBLEM==96)
#define PR_DEFINE "PROBLEMS/BATTERYTORUS/define.h"
#define PR_BC "PROBLEMS/BATTERYTORUS/bc.c"
#define PR_INIT "PROBLEMS/BATTERYTORUS/init.c"
#define PR_KAPPA "PROBLEMS/BATTERYTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/BATTERYTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BATTERYTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BATTERYTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BATTERYTORUS/dump.c"
#define PR_FINGER "PROBLEMS/BATTERYTORUS/finger.c"
#define PR_TOOLS "PROBLEMS/BATTERYTORUS/tools.c"
#define PR_PREPINIT "PROBLEMS/BATTERYTORUS/prepinit.c"
#define PR_POSTINIT "PROBLEMS/BATTERYTORUS/postinit.c"
#endif

#if(PROBLEM==95)
#define PR_DEFINE "PROBLEMS/XINYI/define.h" //most-important, definitions, choices
#define PR_PREPINIT "PROBLEMS/XINYI/prepinit.c" //auxiliary, if you need to precalculate something
#define PR_INIT "PROBLEMS/XINYI/init.c" //initial state
#define PR_POSTINIT "PROBLEMS/XINYI/postinit.c" //auxilary, if you need to e.g. normalize something (magn.field?)
#define PR_BC "PROBLEMS/XINYI/bc.c" //boundary conditions (if not periodic)
#define PR_KAPPA "PROBLEMS/XINYI/kappa.c" //absorption opacity
#define PR_KAPPAES "PROBLEMS/XINYI/kappaes.c" //scattering opacity
#define PR_TOOLS "PROBLEMS/XINYI/tools.c" //auxiliary
#endif


#if(PROBLEM==94)
#define PR_PREPINIT "PROBLEMS/TDEMILIO3D/prepinit.c"
#define PR_DEFINE "PROBLEMS/TDEMILIO3D/define.h"
#define PR_BC "PROBLEMS/TDEMILIO3D/bc.c"
#define PR_INIT "PROBLEMS/TDEMILIO3D/init.c"
#define PR_POSTINIT "PROBLEMS/TDEMILIO3D/postinit.c"
#define PR_KAPPA "PROBLEMS/TDEMILIO3D/kappa.c"
#define PR_KAPPAES "PROBLEMS/TDEMILIO3D/kappaes.c"
#define PR_TOOLS "PROBLEMS/TDEMILIO3D/tools.c"
#endif


#if(PROBLEM==93)

#define PR_DEFINE "PROBLEMS/KATOTORUS/define.h"
#define PR_BC "PROBLEMS/KATOTORUS/bc.c"
#define PR_INIT "PROBLEMS/KATOTORUS/init.c"
#define PR_POSTINIT "PROBLEMS/KATOTORUS/postinit.c"
//#define PR_KAPPA "PROBLEMS/KATOTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/KATOTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/KATOTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/KATOTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/KATOTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/KATOTORUS/tools.c"

#endif


#if(PROBLEM==92)

#define PR_DEFINE "PROBLEMS/ADVDIFF/define.h"
#define PR_BC "PROBLEMS/ADVDIFF/bc.c"
#define PR_INIT "PROBLEMS/ADVDIFF/init.c"
#define PR_KAPPA "PROBLEMS/ADVDIFF/kappa.c"
#define PR_KAPPAES "PROBLEMS/ADVDIFF/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ADVDIFF/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ADVDIFF/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ADVDIFF/dump.c"
#define PR_TOOLS "PROBLEMS/ADVDIFF/tools.c"

#endif

#if(PROBLEM==91)

#define PR_PREPINIT "PROBLEMS/TDEMILIO/prepinit.c"
#define PR_FINGER "PROBLEMS/TDEMILIO/finger.c"
#define PR_DEFINE "PROBLEMS/TDEMILIO/define.h"
#define PR_BC "PROBLEMS/TDEMILIO/bc.c"
#define PR_INIT "PROBLEMS/TDEMILIO/init.c"
#define PR_KAPPA "PROBLEMS/TDEMILIO/kappa.c"
#define PR_KAPPAES "PROBLEMS/TDEMILIO/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/TDEMILIO/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/TDEMILIO/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/TDEMILIO/dump.c"
#define PR_TOOLS "PROBLEMS/TDEMILIO/tools.c"

#endif


#if(PROBLEM==90)

#define PR_DEFINE "PROBLEMS/VETSHADOW/define.h"
#define PR_BC "PROBLEMS/VETSHADOW/bc.c"
#define PR_INIT "PROBLEMS/VETSHADOW/init.c"
#define PR_KAPPA "PROBLEMS/VETSHADOW/kappa.c"
#define PR_KAPPAES "PROBLEMS/VETSHADOW/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/VETSHADOW/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/VETSHADOW/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/VETSHADOW/dump.c"
#define PR_TOOLS "PROBLEMS/VETSHADOW/tools.c"

#endif


#if(PROBLEM==89)

#define PR_DEFINE "PROBLEMS/RADTORUS/define.h"
#define PR_BC "PROBLEMS/RADTORUS/bc.c"
#define PR_INIT "PROBLEMS/RADTORUS/init.c"
#define PR_KAPPA "PROBLEMS/RADTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/RADTORUS/tools.c"
#define PR_FINGER "PROBLEMS/RADTORUS/finger.c"

#endif


#if(PROBLEM==88)

#define PR_DEFINE "PROBLEMS/VETPHIBEAMS/define.h"
#define PR_BC "PROBLEMS/VETPHIBEAMS/bc.c"
#define PR_INIT "PROBLEMS/VETPHIBEAMS/init.c"
#define PR_KAPPA "PROBLEMS/VETPHIBEAMS/kappa.c"
#define PR_KAPPAES "PROBLEMS/VETPHIBEAMS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/VETPHIBEAMS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/VETPHIBEAMS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/VETPHIBEAMS/dump.c"
#define PR_TOOLS "PROBLEMS/VETPHIBEAMS/tools.c"

#endif

#if(PROBLEM==87)

#define PR_DEFINE "PROBLEMS/INFDISK/define.h"
#define PR_BC "PROBLEMS/INFDISK/bc.c"
#define PR_BC_SPECIAL "PROBLEMS/INFDISK/bc_special.c"
#define PR_BC_SPECIAL_LOOP "PROBLEMS/INFDISK/loop_alloc_special.c"
#define PR_INIT "PROBLEMS/INFDISK/init.c"
//#define PR_KAPPA "PROBLEMS/INFDISK/kappa.c"
#define PR_KAPPAES "PROBLEMS/INFDISK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/INFDISK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/INFDISK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/INFDISK/dump.c"
#define PR_TOOLS "PROBLEMS/INFDISK/tools.c"

#endif

#if(PROBLEM==86)

#define PR_DEFINE "PROBLEMS/TFLAT/define.h"
#define PR_BC "PROBLEMS/TFLAT/bc.c"
#define PR_INIT "PROBLEMS/TFLAT/init.c"
#define PR_KAPPA "PROBLEMS/TFLAT/kappa.c"
#define PR_KAPPAES "PROBLEMS/TFLAT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/TFLAT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/TFLAT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/TFLAT/dump.c"
#define PR_TOOLS "PROBLEMS/TFLAT/tools.c"

#endif


#if(PROBLEM==85)


#define PR_DEFINE "PROBLEMS/NSTEST/define.h"
//#define PR_POSTINIT "PROBLEMS/NSTEST/postinit.c"
#define PR_PREPINIT "PROBLEMS/NSTEST/prepinit.c"
#define PR_BC "PROBLEMS/NSTEST/bc.c"
#define PR_INIT "PROBLEMS/NSTEST/init.c"
#define PR_KAPPA "PROBLEMS/NSTEST/kappa.c"
#define PR_KAPPAES "PROBLEMS/NSTEST/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/NSTEST/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/NSTEST/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/NSTEST/dump.c"
#define PR_TOOLS "PROBLEMS/NSTEST/tools.c"

#endif




#if(PROBLEM==84)

#define PR_DEFINE "PROBLEMS/AXISBLOB/define.h"
#define PR_BC "PROBLEMS/AXISBLOB/bc.c"
#define PR_INIT "PROBLEMS/AXISBLOB/init.c"
#define PR_KAPPA "PROBLEMS/AXISBLOB/kappa.c"
#define PR_KAPPAES "PROBLEMS/AXISBLOB/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/AXISBLOB/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/AXISBLOB/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/AXISBLOB/dump.c"
#define PR_TOOLS "PROBLEMS/AXISBLOB/tools.c"

#endif


#if(PROBLEM==83)

#define PR_DEFINE "PROBLEMS/RADFLATTEST/define.h"
#define PR_BC "PROBLEMS/RADFLATTEST/bc.c"
#define PR_INIT "PROBLEMS/RADFLATTEST/init.c"
#define PR_KAPPA "PROBLEMS/RADFLATTEST/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADFLATTEST/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADFLATTEST/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADFLATTEST/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADFLATTEST/dump.c"
#define PR_TOOLS "PROBLEMS/RADFLATTEST/tools.c"

#endif

#if(PROBLEM==82)

#define PR_DEFINE "PROBLEMS/VETTHIN/define.h"
#define PR_BC "PROBLEMS/VETTHIN/bc.c"
#define PR_INIT "PROBLEMS/VETTHIN/init.c"
#define PR_KAPPA "PROBLEMS/VETTHIN/kappa.c"
#define PR_KAPPAES "PROBLEMS/VETTHIN/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/VETTHIN/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/VETTHIN/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/VETTHIN/dump.c"
#define PR_TOOLS "PROBLEMS/VETTHIN/tools.c"

#endif

#if(PROBLEM==81)

#define PR_DEFINE "PROBLEMS/TWOBEAM/define.h"
#define PR_BC "PROBLEMS/TWOBEAM/bc.c"
#define PR_INIT "PROBLEMS/TWOBEAM/init.c"
#define PR_KAPPA "PROBLEMS/TWOBEAM/kappa.c"
#define PR_KAPPAES "PROBLEMS/TWOBEAM/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/TWOBEAM/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/TWOBEAM/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/TWOBEAM/dump.c"
#define PR_TOOLS "PROBLEMS/TWOBEAM/tools.c"
#define PR_FINGER "PROBLEMS/TWOBEAM/finger.c"

#endif


#if(PROBLEM==80)

#define PR_DEFINE "PROBLEMS/GENCLOSURE/define.h"
#define PR_BC "PROBLEMS/GENCLOSURE/bc.c"
#define PR_INIT "PROBLEMS/GENCLOSURE/init.c"
#define PR_KAPPA "PROBLEMS/GENCLOSURE/kappa.c"
#define PR_KAPPAES "PROBLEMS/GENCLOSURE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/GENCLOSURE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/GENCLOSURE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/GENCLOSURE/dump.c"
#define PR_TOOLS "PROBLEMS/GENCLOSURE/tools.c"
#define PR_FINGER "PROBLEMS/GENCLOSURE/finger.c"

#endif

#if(PROBLEM==79)

#define PR_DEFINE "PROBLEMS/SOFTBALL/define.h"
#define PR_BC "PROBLEMS/SOFTBALL/bc.c"
#define PR_INIT "PROBLEMS/SOFTBALL/init.c"
#define PR_KAPPA "PROBLEMS/SOFTBALL/kappa.c"
#define PR_KAPPAES "PROBLEMS/SOFTBALL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/SOFTBALL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/SOFTBALL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/SOFTBALL/dump.c"
#define PR_TOOLS "PROBLEMS/SOFTBALL/tools.c"
#define PR_PREPINIT "PROBLEMS/SOFTBALL/prepinit.c"
#define PR_POSTINIT "PROBLEMS/SOFTBALL/postinit.c"

#endif

#if(PROBLEM==78)

#define PR_DEFINE "PROBLEMS/FUNNY/define.h"
#define PR_BC "PROBLEMS/FUNNY/bc.c"
#define PR_INIT "PROBLEMS/FUNNY/init.c"
#define PR_KAPPA "PROBLEMS/FUNNY/kappa.c"
#define PR_KAPPAES "PROBLEMS/FUNNY/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/FUNNY/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/FUNNY/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/FUNNY/dump.c"
#define PR_TOOLS "PROBLEMS/FUNNY/tools.c"

#endif

#if(PROBLEM==77)

#define PR_DEFINE "PROBLEMS/EMPTY/define.h"
#define PR_BC "PROBLEMS/EMPTY/bc.c"
#define PR_INIT "PROBLEMS/EMPTY/init.c"
#define PR_KAPPA "PROBLEMS/EMPTY/kappa.c"
#define PR_KAPPAES "PROBLEMS/EMPTY/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/EMPTY/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/EMPTY/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/EMPTY/dump.c"
#define PR_TOOLS "PROBLEMS/EMPTY/tools.c"

#endif

#if(PROBLEM==76)

#define PR_DEFINE "PROBLEMS/ZEBRA/define.h"
#define PR_BC "PROBLEMS/ZEBRA/bc.c"
#define PR_INIT "PROBLEMS/ZEBRA/init.c"
#define PR_KAPPA "PROBLEMS/ZEBRA/kappa.c"
#define PR_KAPPAES "PROBLEMS/ZEBRA/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ZEBRA/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ZEBRA/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ZEBRA/dump.c"
#define PR_TOOLS "PROBLEMS/ZEBRA/tools.c"

#endif

#if(PROBLEM==75)

#define PR_DEFINE "PROBLEMS/LBULLET/define.h"
#define PR_BC "PROBLEMS/LBULLET/bc.c"
#define PR_INIT "PROBLEMS/LBULLET/init.c"
#define PR_FINGER "PROBLEMS/LBULLET/finger.c"
#define PR_KAPPA "PROBLEMS/LBULLET/kappa.c"
#define PR_KAPPAES "PROBLEMS/LBULLET/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/LBULLET/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/LBULLET/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/LBULLET/dump.c"
#define PR_TOOLS "PROBLEMS/LBULLET/tools.c"
#define PR_PREPINIT "PROBLEMS/LBULLET/prepinit.c"

#endif

#if(PROBLEM==74)

#define PR_DEFINE "PROBLEMS/ERICTEST/define.h"
#define PR_BC "PROBLEMS/ERICTEST/bc.c"
#define PR_INIT "PROBLEMS/ERICTEST/init.c"
#define PR_KAPPA "PROBLEMS/ERICTEST/kappa.c"
#define PR_KAPPAES "PROBLEMS/ERICTEST/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ERICTEST/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ERICTEST/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ERICTEST/dump.c"
#define PR_TOOLS "PROBLEMS/ERICTEST/tools.c"

#endif

#if(PROBLEM==73)

#define PR_DEFINE "PROBLEMS/VISCZERO/define.h"
#define PR_BC "PROBLEMS/VISCZERO/bc.c"
#define PR_INIT "PROBLEMS/VISCZERO/init.c"
#define PR_KAPPA "PROBLEMS/VISCZERO/kappa.c"
#define PR_KAPPAES "PROBLEMS/VISCZERO/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/VISCZERO/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/VISCZERO/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/VISCZERO/dump.c"
#define PR_TOOLS "PROBLEMS/VISCZERO/tools.c"

#endif

#if(PROBLEM==72)

#define PR_DEFINE "PROBLEMS/BOWSHOCK/define.h"
#define PR_BC "PROBLEMS/BOWSHOCK/bc.c"
#define PR_INIT "PROBLEMS/BOWSHOCK/init.c"
#define PR_KAPPA "PROBLEMS/BOWSHOCK/kappa.c"
#define PR_KAPPAES "PROBLEMS/BOWSHOCK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BOWSHOCK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BOWSHOCK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BOWSHOCK/dump.c"
#define PR_TOOLS "PROBLEMS/BOWSHOCK/tools.c"

#endif

#if(PROBLEM==71)

#define PR_DEFINE "PROBLEMS/MBBBLOB/define.h"
#define PR_BC "PROBLEMS/MBBBLOB/bc.c"
#define PR_INIT "PROBLEMS/MBBBLOB/init.c"
#define PR_KAPPA "PROBLEMS/MBBBLOB/kappa.c"
#define PR_KAPPAES "PROBLEMS/MBBBLOB/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MBBBLOB/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MBBBLOB/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MBBBLOB/dump.c"
#define PR_TOOLS "PROBLEMS/MBBBLOB/tools.c"

#endif

#if(PROBLEM==70)

#define PR_DEFINE "PROBLEMS/M1SHOCK/define.h"
#define PR_BC "PROBLEMS/M1SHOCK/bc.c"
#define PR_INIT "PROBLEMS/M1SHOCK/init.c"
#define PR_KAPPA "PROBLEMS/M1SHOCK/kappa.c"
#define PR_KAPPAES "PROBLEMS/M1SHOCK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/M1SHOCK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/M1SHOCK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/M1SHOCK/dump.c"
#define PR_TOOLS "PROBLEMS/M1SHOCK/tools.c"

#endif

#if(PROBLEM==69)

#define PR_DEFINE "PROBLEMS/INJDISK/define.h"
#define PR_BC "PROBLEMS/INJDISK/bc.c"
#define PR_INIT "PROBLEMS/INJDISK/init.c"
#define PR_KAPPA "PROBLEMS/INJDISK/kappa.c"
#define PR_KAPPAES "PROBLEMS/INJDISK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/INJDISK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/INJDISK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/INJDISK/dump.c"
#define PR_TOOLS "PROBLEMS/INJDISK/tools.c"
#define PR_POSTINIT "PROBLEMS/INJDISK/postinit.c"

#endif

#if(PROBLEM==68)

#define PR_DEFINE "PROBLEMS/RMHDWAVE/define.h"
#define PR_BC "PROBLEMS/RMHDWAVE/bc.c"
#define PR_INIT "PROBLEMS/RMHDWAVE/init.c"
#define PR_KAPPA "PROBLEMS/RMHDWAVE/kappa.c"
#define PR_KAPPAES "PROBLEMS/RMHDWAVE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RMHDWAVE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RMHDWAVE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RMHDWAVE/dump.c"
#define PR_TOOLS "PROBLEMS/RMHDWAVE/tools.c"

#endif

#if(PROBLEM==67)

#define PR_DEFINE "PROBLEMS/LRTORUS/define.h"
#define PR_BC "PROBLEMS/LRTORUS/bc.c"
#define PR_INIT "PROBLEMS/LRTORUS/init.c"
#define PR_KAPPA "PROBLEMS/LRTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/LRTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/LRTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/LRTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/LRTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/LRTORUS/tools.c"
#define PR_POSTINIT "PROBLEMS/LRTORUS/postinit.c"

#endif

#if(PROBLEM==66)

#define PR_DEFINE "PROBLEMS/KTORUS/define.h"
#define PR_BC "PROBLEMS/KTORUS/bc.c"
#define PR_INIT "PROBLEMS/KTORUS/init.c"
#define PR_KAPPA "PROBLEMS/KTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/KTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/KTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/KTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/KTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/KTORUS/tools.c"
#define PR_POSTINIT "PROBLEMS/KTORUS/postinit.c"

#endif

#if(PROBLEM==65)

#define PR_DEFINE "PROBLEMS/MRTORUS/define.h"
#define PR_BC "PROBLEMS/MRTORUS/bc.c"
#define PR_INIT "PROBLEMS/MRTORUS/init.c"
#define PR_KAPPA "PROBLEMS/MRTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/MRTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MRTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MRTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MRTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/MRTORUS/tools.c"
#define PR_POSTINIT "PROBLEMS/MRTORUS/postinit.c"

#endif

#if(PROBLEM==64)

#define PR_DEFINE "PROBLEMS/MAGDONUT/define.h"
#define PR_BC "PROBLEMS/MAGDONUT/bc.c"
#define PR_INIT "PROBLEMS/MAGDONUT/init.c"
#define PR_KAPPA "PROBLEMS/MAGDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/MAGDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MAGDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MAGDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MAGDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/MAGDONUT/tools.c"
//#define PR_POSTINIT "PROBLEMS/MAGDONUT/postinit.c"

#endif

#if(PROBLEM==63)

#define PR_DEFINE "PROBLEMS/ORSZAG/define.h"
#define PR_BC "PROBLEMS/ORSZAG/bc.c"
#define PR_INIT "PROBLEMS/ORSZAG/init.c"
#define PR_KAPPA "PROBLEMS/ORSZAG/kappa.c"
#define PR_KAPPAES "PROBLEMS/ORSZAG/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ORSZAG/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ORSZAG/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ORSZAG/dump.c"
#define PR_TOOLS "PROBLEMS/ORSZAG/tools.c"

#endif

#if(PROBLEM==62)

#define PR_DEFINE "PROBLEMS/MAGTUBES/define.h"
#define PR_BC "PROBLEMS/MAGTUBES/bc.c"
#define PR_INIT "PROBLEMS/MAGTUBES/init.c"
#define PR_KAPPA "PROBLEMS/MAGTUBES/kappa.c"
#define PR_KAPPAES "PROBLEMS/MAGTUBES/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MAGTUBES/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MAGTUBES/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MAGTUBES/dump.c"
#define PR_TOOLS "PROBLEMS/MAGTUBES/tools.c"

#endif

#if(PROBLEM==61)

#define PR_DEFINE "PROBLEMS/MAGTESTS/define.h"
#define PR_BC "PROBLEMS/MAGTESTS/bc.c"
#define PR_INIT "PROBLEMS/MAGTESTS/init.c"
#define PR_KAPPA "PROBLEMS/MAGTESTS/kappa.c"
#define PR_KAPPAES "PROBLEMS/MAGTESTS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MAGTESTS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MAGTESTS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MAGTESTS/dump.c"
#define PR_TOOLS "PROBLEMS/MAGTESTS/tools.c"

#endif

#if(PROBLEM==60)

#define PR_DEFINE "PROBLEMS/JONTESTS/define.h"
#define PR_BC "PROBLEMS/JONTESTS/bc.c"
#define PR_INIT "PROBLEMS/JONTESTS/init.c"
#define PR_KAPPA "PROBLEMS/JONTESTS/kappa.c"
#define PR_KAPPAES "PROBLEMS/JONTESTS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/JONTESTS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/JONTESTS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/JONTESTS/dump.c"
#define PR_TOOLS "PROBLEMS/JONTESTS/tools.c"

#endif

#if(PROBLEM==59)

#define PR_DEFINE "PROBLEMS/BBBLOB/define.h"
#define PR_BC "PROBLEMS/BBBLOB/bc.c"
#define PR_INIT "PROBLEMS/BBBLOB/init.c"
#define PR_KAPPA "PROBLEMS/BBBLOB/kappa.c"
#define PR_KAPPAES "PROBLEMS/BBBLOB/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BBBLOB/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BBBLOB/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BBBLOB/dump.c"
#define PR_TOOLS "PROBLEMS/BBBLOB/tools.c"

#endif

#if(PROBLEM==58)

#define PR_DEFINE "PROBLEMS/LUKE/define.h"
#define PR_BC "PROBLEMS/LUKE/bc.c"
#define PR_INIT "PROBLEMS/LUKE/init.c"
#define PR_KAPPA "PROBLEMS/LUKE/kappa.c"
#define PR_KAPPAES "PROBLEMS/LUKE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/LUKE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/LUKE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/LUKE/dump.c"
#define PR_TOOLS "PROBLEMS/LUKE/tools.c"

#endif

#if(PROBLEM==57)

#define PR_DEFINE "PROBLEMS/LIMOTORUS/define.h"
#define PR_BC "PROBLEMS/LIMOTORUS/bc.c"
#define PR_INIT "PROBLEMS/LIMOTORUS/init.c"
#define PR_KAPPA "PROBLEMS/LIMOTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/LIMOTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/LIMOTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/LIMOTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/LIMOTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/LIMOTORUS/tools.c"

#endif

#if(PROBLEM==56)

#define PR_DEFINE "PROBLEMS/EDDSPH/define.h"
#define PR_BC "PROBLEMS/EDDSPH/bc.c"
#define PR_INIT "PROBLEMS/EDDSPH/init.c"
#define PR_KAPPA "PROBLEMS/EDDSPH/kappa.c"
#define PR_KAPPAES "PROBLEMS/EDDSPH/kappaes.c"
//#define PR_OUT2GIF_2D "PROBLEMS/EDDSPH/out2gif_2d.c"
//#define PR_OUT2GIF_1D "PROBLEMS/EDDSPH/out2gif_1d.c"
//#define PR_DUMP "PROBLEMS/EDDSPH/dump.c"
#define PR_TOOLS "PROBLEMS/EDDSPH/tools.c"

#endif

#if(PROBLEM==55)

#define PR_PREPINIT "PROBLEMS/G2STATIC/prepinit.c"
#define PR_DEFINE "PROBLEMS/G2STATIC/define.h"
#define PR_BC "PROBLEMS/G2STATIC/bc.c"
#define PR_INIT "PROBLEMS/G2STATIC/init.c"
#define PR_KAPPA "PROBLEMS/G2STATIC/kappa.c"
#define PR_KAPPAES "PROBLEMS/G2STATIC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/G2STATIC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/G2STATIC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/G2STATIC/dump.c"
#define PR_TOOLS "PROBLEMS/G2STATIC/tools.c"
#define PR_DEFS "PROBLEMS/G2STATIC/defs.h"
#define PR_FINGER "PROBLEMS/G2STATIC/finger.c"
#endif

#if(PROBLEM==54)

#define PR_DEFINE "PROBLEMS/RVDISK/define.h"
#define PR_BC "PROBLEMS/RVDISK/bc.c"
#define PR_INIT "PROBLEMS/RVDISK/init.c"
#define PR_KAPPA "PROBLEMS/RVDISK/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVDISK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVDISK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVDISK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVDISK/dump.c"
#define PR_TOOLS "PROBLEMS/RVDISK/tools.c"

#endif

#if(PROBLEM==53)

#define PR_DEFINE "PROBLEMS/1DDONUT/define.h"
#define PR_BC "PROBLEMS/1DDONUT/bc.c"
#define PR_INIT "PROBLEMS/1DDONUT/init.c"
#define PR_KAPPA "PROBLEMS/1DDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/1DDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/1DDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/1DDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/1DDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/1DDONUT/tools.c"

#endif

#if(PROBLEM==52)

#define PR_DEFINE "PROBLEMS/SPHFLAT/define.h"
#define PR_BC "PROBLEMS/SPHFLAT/bc.c"
#define PR_INIT "PROBLEMS/SPHFLAT/init.c"
#define PR_KAPPA "PROBLEMS/SPHFLAT/kappa.c"
#define PR_KAPPAES "PROBLEMS/SPHFLAT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/SPHFLAT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/SPHFLAT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/SPHFLAT/dump.c"
#define PR_TOOLS "PROBLEMS/SPHFLAT/tools.c"

#endif

#if(PROBLEM==51)

#define PR_DEFINE "PROBLEMS/HUBBLE/define.h"
#define PR_BC "PROBLEMS/HUBBLE/bc.c"
#define PR_INIT "PROBLEMS/HUBBLE/init.c"
#define PR_KAPPA "PROBLEMS/HUBBLE/kappa.c"
#define PR_KAPPAES "PROBLEMS/HUBBLE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HUBBLE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HUBBLE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HUBBLE/dump.c"
#define PR_TOOLS "PROBLEMS/HUBBLE/tools.c"

#endif

#if(PROBLEM==50)

#define PR_DEFINE "PROBLEMS/HDWAVE/define.h"
#define PR_BC "PROBLEMS/HDWAVE/bc.c"
#define PR_INIT "PROBLEMS/HDWAVE/init.c"
#define PR_KAPPA "PROBLEMS/HDWAVE/kappa.c"
#define PR_KAPPAES "PROBLEMS/HDWAVE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HDWAVE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HDWAVE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HDWAVE/dump.c"
#define PR_TOOLS "PROBLEMS/HDWAVE/tools.c"

#endif


#if(PROBLEM==49)

#define PR_DEFINE "PROBLEMS/RVRING/define.h"
#define PR_BC "PROBLEMS/RVRING/bc.c"
#define PR_INIT "PROBLEMS/RVRING/init.c"
#define PR_KAPPA "PROBLEMS/RVRING/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVRING/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVRING/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVRING/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVRING/dump.c"
#define PR_TOOLS "PROBLEMS/RVRING/tools.c"

#endif

#if(PROBLEM==48)

#define PR_DEFINE "PROBLEMS/RVDTEST/define.h"
#define PR_BC "PROBLEMS/RVDTEST/bc.c"
#define PR_INIT "PROBLEMS/RVDTEST/init.c"
#define PR_KAPPA "PROBLEMS/RVDTEST/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVDTEST/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVDTEST/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVDTEST/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVDTEST/dump.c"
#define PR_TOOLS "PROBLEMS/RVDTEST/tools.c"

#endif

#if(PROBLEM==47)

#define PR_DEFINE "PROBLEMS/MFDONUT/define.h"
#define PR_BC "PROBLEMS/MFDONUT/bc.c"
#define PR_INIT "PROBLEMS/MFDONUT/init.c"
#define PR_KAPPA "PROBLEMS/MFDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/MFDONUT/tools.c"

#endif

#if(PROBLEM==46)

#define PR_DEFINE "PROBLEMS/MFRADNTSPH/define.h"
#define PR_BC "PROBLEMS/MFRADNTSPH/bc.c"
#define PR_INIT "PROBLEMS/MFRADNTSPH/init.c"
#define PR_KAPPA "PROBLEMS/MFRADNTSPH/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFRADNTSPH/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFRADNTSPH/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFRADNTSPH/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFRADNTSPH/dump.c"
#define PR_TOOLS "PROBLEMS/MFRADNTSPH/tools.c"

#endif

#if(PROBLEM==45)

#define PR_DEFINE "PROBLEMS/MFRADNTCYL/define.h"
#define PR_BC "PROBLEMS/MFRADNTCYL/bc.c"
#define PR_INIT "PROBLEMS/MFRADNTCYL/init.c"
#define PR_KAPPA "PROBLEMS/MFRADNTCYL/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFRADNTCYL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFRADNTCYL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFRADNTCYL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFRADNTCYL/dump.c"
#define PR_TOOLS "PROBLEMS/MFRADNTCYL/tools.c"

#endif

#if(PROBLEM==44)

#define PR_DEFINE "PROBLEMS/RADNTCYL/define.h"
#define PR_BC "PROBLEMS/RADNTCYL/bc.c"
#define PR_INIT "PROBLEMS/RADNTCYL/init.c"
#define PR_KAPPA "PROBLEMS/RADNTCYL/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADNTCYL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADNTCYL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADNTCYL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADNTCYL/dump.c"
#define PR_TOOLS "PROBLEMS/RADNTCYL/tools.c"

#endif

#if(PROBLEM==43)

#define PR_DEFINE "PROBLEMS/RVDONUTIN/define.h"
#define PR_BC "PROBLEMS/RVDONUTIN/bc.c"
#define PR_INIT "PROBLEMS/RVDONUTIN/init.c"
#define PR_KAPPA "PROBLEMS/RVDONUTIN/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVDONUTIN/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVDONUTIN/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVDONUTIN/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVDONUTIN/dump.c"
#define PR_TOOLS "PROBLEMS/RVDONUTIN/tools.c"

#endif

#if(PROBLEM==42)

#define PR_DEFINE "PROBLEMS/RVDONUT/define.h"
#define PR_BC "PROBLEMS/RVDONUT/bc.c"
#define PR_INIT "PROBLEMS/RVDONUT/init.c"
#define PR_KAPPA "PROBLEMS/RVDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/RVDONUT/tools.c"

#endif

#if(PROBLEM==41)

#define PR_DEFINE "PROBLEMS/FLATDOT/define.h"
#define PR_BC "PROBLEMS/FLATDOT/bc.c"
#define PR_INIT "PROBLEMS/FLATDOT/init.c"
#define PR_KAPPA "PROBLEMS/FLATDOT/kappa.c"
#define PR_KAPPAES "PROBLEMS/FLATDOT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/FLATDOT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/FLATDOT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/FLATDOT/dump.c"
#define PR_TOOLS "PROBLEMS/FLATDOT/tools.c"
#define PR_FINGER "PROBLEMS/FLATDOT/finger.c"

#endif

#if(PROBLEM==40)

#define PR_DEFINE "PROBLEMS/CYLBEAMCART/define.h"
#define PR_BC "PROBLEMS/CYLBEAMCART/bc.c"
#define PR_INIT "PROBLEMS/CYLBEAMCART/init.c"
#define PR_KAPPA "PROBLEMS/CYLBEAMCART/kappa.c"
#define PR_KAPPAES "PROBLEMS/CYLBEAMCART/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/CYLBEAMCART/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/CYLBEAMCART/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/CYLBEAMCART/dump.c"
#define PR_TOOLS "PROBLEMS/CYLBEAMCART/tools.c"

#endif

#if(PROBLEM==37)

#define PR_DEFINE "PROBLEMS/MFCYLBEAM/define.h"
#define PR_BC "PROBLEMS/MFCYLBEAM/bc.c"
#define PR_INIT "PROBLEMS/MFCYLBEAM/init.c"
#define PR_KAPPA "PROBLEMS/MFCYLBEAM/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFCYLBEAM/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFCYLBEAM/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFCYLBEAM/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFCYLBEAM/dump.c"
#define PR_TOOLS "PROBLEMS/MFCYLBEAM/tools.c"

#endif

#if(PROBLEM==36)

#define PR_DEFINE "PROBLEMS/MFDOTS/define.h"
#define PR_BC "PROBLEMS/MFDOTS/bc.c"
#define PR_INIT "PROBLEMS/MFDOTS/init.c"
#define PR_KAPPA "PROBLEMS/MFDOTS/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFDOTS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFDOTS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFDOTS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFDOTS/dump.c"
#define PR_TOOLS "PROBLEMS/MFDOTS/tools.c"
#define PR_FINGER "PROBLEMS/MFDOTS/finger.c"

#endif

#if(PROBLEM==35)

#define PR_DEFINE "PROBLEMS/MFBEAMS/define.h"
#define PR_BC "PROBLEMS/MFBEAMS/bc.c"
#define PR_INIT "PROBLEMS/MFBEAMS/init.c"
#define PR_KAPPA "PROBLEMS/MFBEAMS/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFBEAMS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFBEAMS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFBEAMS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFBEAMS/dump.c"
#define PR_TOOLS "PROBLEMS/MFBEAMS/tools.c"

#endif

#if(PROBLEM==34)

#define PR_DEFINE "PROBLEMS/MFPULSE/define.h"
#define PR_BC "PROBLEMS/MFPULSE/bc.c"
#define PR_INIT "PROBLEMS/MFPULSE/init.c"
#define PR_KAPPA "PROBLEMS/MFPULSE/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFPULSE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFPULSE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFPULSE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFPULSE/dump.c"
#define PR_TOOLS "PROBLEMS/MFPULSE/tools.c"

#endif

#if(PROBLEM==33)

#define PR_DEFINE "PROBLEMS/RADDOT/define.h"
#define PR_BC "PROBLEMS/RADDOT/bc.c"
#define PR_INIT "PROBLEMS/RADDOT/init.c"
#define PR_KAPPA "PROBLEMS/RADDOT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADDOT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADDOT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADDOT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADDOT/dump.c"
#define PR_TOOLS "PROBLEMS/RADDOT/tools.c"
#define PR_FINGER "PROBLEMS/RADDOT/finger.c"

#endif

#if(PROBLEM==32)

#define PR_DEFINE "PROBLEMS/CYLBEAM/define.h"
#define PR_BC "PROBLEMS/CYLBEAM/bc.c"
#define PR_INIT "PROBLEMS/CYLBEAM/init.c"
#define PR_KAPPA "PROBLEMS/CYLBEAM/kappa.c"
#define PR_KAPPAES "PROBLEMS/CYLBEAM/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/CYLBEAM/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/CYLBEAM/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/CYLBEAM/dump.c"
#define PR_TOOLS "PROBLEMS/CYLBEAM/tools.c"

#endif

#if(PROBLEM==31)

#define PR_DEFINE "PROBLEMS/FLATDISK/define.h"
#define PR_BC "PROBLEMS/FLATDISK/bc.c"
#define PR_INIT "PROBLEMS/FLATDISK/init.c"
#define PR_KAPPA "PROBLEMS/FLATDISK/kappa.c"
#define PR_KAPPAES "PROBLEMS/FLATDISK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/FLATDISK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/FLATDISK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/FLATDISK/dump.c"
#define PR_TOOLS "PROBLEMS/FLATDISK/tools.c"

#endif

#if(PROBLEM==30)

#define PR_DEFINE "PROBLEMS/RADNT/define.h"
#define PR_BC "PROBLEMS/RADNT/bc.c"
#define PR_INIT "PROBLEMS/RADNT/init.c"
#define PR_KAPPA "PROBLEMS/RADNT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADNT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADNT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADNT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADNT/dump.c"
#define PR_TOOLS "PROBLEMS/RADNT/tools.c"

#endif

#if(PROBLEM==1)

#define PR_DEFINE "PROBLEMS/RADBEAM2D/define.h"
#define PR_BC "PROBLEMS/RADBEAM2D/bc.c"
#define PR_INIT "PROBLEMS/RADBEAM2D/init.c"
#define PR_KAPPA "PROBLEMS/RADBEAM2D/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADBEAM2D/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADBEAM2D/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADBEAM2D/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADBEAM2D/dump.c"
#define PR_TOOLS "PROBLEMS/RADBEAM2D/tools.c"

#endif

#if(PROBLEM==2)

#define PR_DEFINE "PROBLEMS/RADINFALL/define.h"
#define PR_BC "PROBLEMS/RADINFALL/bc.c"
#define PR_INIT "PROBLEMS/RADINFALL/init.c"
#define PR_KAPPA "PROBLEMS/RADINFALL/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADINFALL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADINFALL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADINFALL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADINFALL/dump.c"
#define PR_TOOLS "PROBLEMS/RADINFALL/tools.c"

#endif

#if(PROBLEM==3)

#define PR_DEFINE "PROBLEMS/DONUT/define.h"
#define PR_BC "PROBLEMS/DONUT/bc.c"
#define PR_INIT "PROBLEMS/DONUT/init.c"
#define PR_KAPPA "PROBLEMS/DONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/DONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/DONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/DONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/DONUT/dump.c"
#define PR_TOOLS "PROBLEMS/DONUT/tools.c"

#endif

#if(PROBLEM==4)

#define PR_DEFINE "PROBLEMS/GEODESICINFALL/define.h"
#define PR_BC "PROBLEMS/GEODESICINFALL/bc.c"
#define PR_INIT "PROBLEMS/GEODESICINFALL/init.c"
#define PR_KAPPA "PROBLEMS/GEODESICINFALL/kappa.c"
#define PR_KAPPAES "PROBLEMS/GEODESICINFALL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/GEODESICINFALL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/GEODESICINFALL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/GEODESICINFALL/dump.c"
#define PR_TOOLS "PROBLEMS/GEODESICINFALL/tools.c"

#endif

#if(PROBLEM==5)

#define PR_DEFINE "PROBLEMS/EDDINFALL/define.h"
#define PR_BC "PROBLEMS/EDDINFALL/bc.c"
#define PR_INIT "PROBLEMS/EDDINFALL/init.c"
#define PR_KAPPA "PROBLEMS/EDDINFALL/kappa.c"
#define PR_KAPPAES "PROBLEMS/EDDINFALL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/EDDINFALL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/EDDINFALL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/EDDINFALL/dump.c"
#define PR_TOOLS "PROBLEMS/EDDINFALL/tools.c"

#endif

#if(PROBLEM==6)

#define PR_DEFINE "PROBLEMS/RADTUBE/define.h"
#define PR_BC "PROBLEMS/RADTUBE/bc.c"
#define PR_INIT "PROBLEMS/RADTUBE/init.c"
#define PR_KAPPA "PROBLEMS/RADTUBE/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADTUBE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADTUBE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADTUBE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADTUBE/dump.c"
#define PR_TOOLS "PROBLEMS/RADTUBE/tools.c"

#endif

#if(PROBLEM==7)

#define PR_DEFINE "PROBLEMS/BONDI/define.h"
#define PR_PREPINIT "PROBLEMS/BONDI/prepinit.c"
//#define PR_POSTINIT "PROBLEMS/BONDI/postinit.c"
#define PR_BC "PROBLEMS/BONDI/bc.c"
#define PR_INIT "PROBLEMS/BONDI/init.c"
#define PR_KAPPA "PROBLEMS/BONDI/kappa.c"
#define PR_KAPPAES "PROBLEMS/BONDI/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BONDI/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BONDI/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BONDI/dump.c"
#define PR_TOOLS "PROBLEMS/BONDI/tools.c"

#endif

#if(PROBLEM==8)

#define PR_DEFINE "PROBLEMS/HDTUBE/define.h"
#define PR_BC "PROBLEMS/HDTUBE/bc.c"
#define PR_INIT "PROBLEMS/HDTUBE/init.c"
#define PR_KAPPA "PROBLEMS/HDTUBE/kappa.c"
#define PR_KAPPAES "PROBLEMS/HDTUBE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HDTUBE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HDTUBE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HDTUBE/dump.c"
#define PR_TOOLS "PROBLEMS/HDTUBE/tools.c"

#endif

#if(PROBLEM==9)

#define PR_DEFINE "PROBLEMS/HDTUBE2D/define.h"
#define PR_BC "PROBLEMS/HDTUBE2D/bc.c"
#define PR_INIT "PROBLEMS/HDTUBE2D/init.c"
#define PR_KAPPA "PROBLEMS/HDTUBE2D/kappa.c"
#define PR_KAPPAES "PROBLEMS/HDTUBE2D/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HDTUBE2D/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HDTUBE2D/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HDTUBE2D/dump.c"
#define PR_TOOLS "PROBLEMS/HDTUBE2D/tools.c"

#endif

#if(PROBLEM==10)

#define PR_DEFINE "PROBLEMS/RADPULSE/define.h"
#define PR_BC "PROBLEMS/RADPULSE/bc.c"
#define PR_INIT "PROBLEMS/RADPULSE/init.c"
#define PR_KAPPA "PROBLEMS/RADPULSE/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADPULSE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADPULSE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADPULSE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADPULSE/dump.c"
#define PR_TOOLS "PROBLEMS/RADPULSE/tools.c"

#endif

#if(PROBLEM==11)

#define PR_DEFINE "PROBLEMS/RADSHADOW/define.h"
#define PR_BC "PROBLEMS/RADSHADOW/bc.c"
#define PR_INIT "PROBLEMS/RADSHADOW/init.c"
#define PR_KAPPA "PROBLEMS/RADSHADOW/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADSHADOW/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADSHADOW/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADSHADOW/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADSHADOW/dump.c"
#define PR_TOOLS "PROBLEMS/RADSHADOW/tools.c"

#endif

#if(PROBLEM==12)

#define PR_DEFINE "PROBLEMS/RADATM/define.h"
#define PR_BC "PROBLEMS/RADATM/bc.c"
#define PR_INIT "PROBLEMS/RADATM/init.c"
#define PR_KAPPA "PROBLEMS/RADATM/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADATM/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADATM/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADATM/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADATM/dump.c"
#define PR_TOOLS "PROBLEMS/RADATM/tools.c"

#endif

#if(PROBLEM==13)

#define PR_DEFINE "PROBLEMS/DONUTOSC/define.h"
#define PR_BC "PROBLEMS/DONUTOSC/bc.c"
#define PR_INIT "PROBLEMS/DONUTOSC/init.c"
#define PR_KAPPA "PROBLEMS/DONUTOSC/kappa.c"
#define PR_KAPPAES "PROBLEMS/DONUTOSC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/DONUTOSC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/DONUTOSC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/DONUTOSC/dump.c"
#define PR_TOOLS "PROBLEMS/DONUTOSC/tools.c"

#endif

#if(PROBLEM==14)

#define PR_DEFINE "PROBLEMS/RADWAVEBC/define.h"
#define PR_BC "PROBLEMS/RADWAVEBC/bc.c"
#define PR_INIT "PROBLEMS/RADWAVEBC/init.c"
#define PR_KAPPA "PROBLEMS/RADWAVEBC/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADWAVEBC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADWAVEBC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADWAVEBC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADWAVEBC/dump.c"
#define PR_TOOLS "PROBLEMS/RADWAVEBC/tools.c"

#endif

#if(PROBLEM==15)

#define PR_DEFINE "PROBLEMS/RADWAVE/define.h"
#define PR_BC "PROBLEMS/RADWAVE/bc.c"
#define PR_INIT "PROBLEMS/RADWAVE/init.c"
#define PR_KAPPA "PROBLEMS/RADWAVE/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADWAVE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADWAVE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADWAVE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADWAVE/dump.c"
#define PR_TOOLS "PROBLEMS/RADWAVE/tools.c"

#endif

#if(PROBLEM==16)

#define PR_DEFINE "PROBLEMS/RADPULSE3D/define.h"
#define PR_BC "PROBLEMS/RADPULSE3D/bc.c"
#define PR_INIT "PROBLEMS/RADPULSE3D/init.c"
#define PR_KAPPA "PROBLEMS/RADPULSE3D/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADPULSE3D/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADPULSE3D/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADPULSE3D/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADPULSE3D/dump.c"
#define PR_TOOLS "PROBLEMS/RADPULSE3D/tools.c"

#endif

#if(PROBLEM==17)

#define PR_DEFINE "PROBLEMS/RADDBLSHADOW/define.h"
#define PR_BC "PROBLEMS/RADDBLSHADOW/bc.c"
#define PR_INIT "PROBLEMS/RADDBLSHADOW/init.c"
#define PR_KAPPA "PROBLEMS/RADDBLSHADOW/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADDBLSHADOW/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADDBLSHADOW/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADDBLSHADOW/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADDBLSHADOW/dump.c"
#define PR_TOOLS "PROBLEMS/RADDBLSHADOW/tools.c"

#endif

#if(PROBLEM==18)

#define PR_DEFINE "PROBLEMS/ATMSTATIC/define.h"
#define PR_BC "PROBLEMS/ATMSTATIC/bc.c"
#define PR_INIT "PROBLEMS/ATMSTATIC/init.c"
#define PR_KAPPA "PROBLEMS/ATMSTATIC/kappa.c"
#define PR_KAPPAES "PROBLEMS/ATMSTATIC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ATMSTATIC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ATMSTATIC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ATMSTATIC/dump.c"
#define PR_TOOLS "PROBLEMS/ATMSTATIC/tools.c"

#endif

#if(PROBLEM==19)

#define PR_DEFINE "PROBLEMS/RADBEAM2DKS/define.h"
#define PR_BC "PROBLEMS/RADBEAM2DKS/bc.c"
#define PR_INIT "PROBLEMS/RADBEAM2DKS/init.c"
#define PR_KAPPA "PROBLEMS/RADBEAM2DKS/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADBEAM2DKS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADBEAM2DKS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADBEAM2DKS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADBEAM2DKS/dump.c"
#define PR_TOOLS "PROBLEMS/RADBEAM2DKS/tools.c"

#endif

#if(PROBLEM==20)

#define PR_DEFINE "PROBLEMS/ATMKS/define.h"
#define PR_BC "PROBLEMS/ATMKS/bc.c"
#define PR_INIT "PROBLEMS/ATMKS/init.c"
#define PR_KAPPA "PROBLEMS/ATMKS/kappa.c"
#define PR_KAPPAES "PROBLEMS/ATMKS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ATMKS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ATMKS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ATMKS/dump.c"
#define PR_TOOLS "PROBLEMS/ATMKS/tools.c"

#endif

#if(PROBLEM==21)

#define PR_DEFINE "PROBLEMS/DONUTKS/define.h"
#define PR_BC "PROBLEMS/DONUTKS/bc.c"
#define PR_INIT "PROBLEMS/DONUTKS/init.c"
#define PR_KAPPA "PROBLEMS/DONUTKS/kappa.c"
#define PR_KAPPAES "PROBLEMS/DONUTKS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/DONUTKS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/DONUTKS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/DONUTKS/dump.c"
#define PR_TOOLS "PROBLEMS/DONUTKS/tools.c"

#endif

#if(PROBLEM==22)

#define PR_DEFINE "PROBLEMS/DONUTMKS1/define.h"
#define PR_BC "PROBLEMS/DONUTMKS1/bc.c"
#define PR_INIT "PROBLEMS/DONUTMKS1/init.c"
#define PR_KAPPA "PROBLEMS/DONUTMKS1/kappa.c"
#define PR_KAPPAES "PROBLEMS/DONUTMKS1/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/DONUTMKS1/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/DONUTMKS1/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/DONUTMKS1/dump.c"
#define PR_TOOLS "PROBLEMS/DONUTMKS1/tools.c"

#endif



#if(PROBLEM==23)

#define PR_DEFINE "PROBLEMS/ATMMKS1/define.h"
#define PR_BC "PROBLEMS/ATMMKS1/bc.c"
#define PR_INIT "PROBLEMS/ATMMKS1/init.c"
#define PR_KAPPA "PROBLEMS/ATMMKS1/kappa.c"
#define PR_KAPPAES "PROBLEMS/ATMMKS1/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ATMMKS1/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ATMMKS1/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ATMMKS1/dump.c"
#define PR_TOOLS "PROBLEMS/ATMMKS1/tools.c"

#endif

#if(PROBLEM==24)

#define PR_DEFINE "PROBLEMS/RADBEAMFLAT/define.h"
#define PR_BC "PROBLEMS/RADBEAMFLAT/bc.c"
#define PR_INIT "PROBLEMS/RADBEAMFLAT/init.c"
#define PR_KAPPA "PROBLEMS/RADBEAMFLAT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADBEAMFLAT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADBEAMFLAT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADBEAMFLAT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADBEAMFLAT/dump.c"
#define PR_TOOLS "PROBLEMS/RADBEAMFLAT/tools.c"

#endif

#if(PROBLEM==25)

#define PR_DEFINE "PROBLEMS/RDONUT/define.h"
#define PR_BC "PROBLEMS/RDONUT/bc.c"
#define PR_INIT "PROBLEMS/RDONUT/init.c"
#define PR_KAPPA "PROBLEMS/RDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/RDONUT/tools.c"

#endif

#if(PROBLEM==26)

#define PR_DEFINE "PROBLEMS/RADBEAM2DKSVERT/define.h"
#define PR_BC "PROBLEMS/RADBEAM2DKSVERT/bc.c"
#define PR_INIT "PROBLEMS/RADBEAM2DKSVERT/init.c"
#define PR_KAPPA "PROBLEMS/RADBEAM2DKSVERT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADBEAM2DKSVERT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADBEAM2DKSVERT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADBEAM2DKSVERT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADBEAM2DKSVERT/dump.c"
#define PR_TOOLS "PROBLEMS/RADBEAM2DKSVERT/tools.c"

#endif

#if(PROBLEM==27)

#define PR_DEFINE "PROBLEMS/RADFLATNESS/define.h"
#define PR_BC "PROBLEMS/RADFLATNESS/bc.c"
#define PR_INIT "PROBLEMS/RADFLATNESS/init.c"
#define PR_KAPPA "PROBLEMS/RADFLATNESS/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADFLATNESS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADFLATNESS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADFLATNESS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADFLATNESS/dump.c"
#define PR_TOOLS "PROBLEMS/RADFLATNESS/tools.c"

#endif

#if(PROBLEM==28)

#define PR_DEFINE "PROBLEMS/BOWSHOCK/define.h"
#define PR_BC "PROBLEMS/BOWSHOCK/bc.c"
#define PR_INIT "PROBLEMS/BOWSHOCK/init.c"
#define PR_KAPPA "PROBLEMS/BOWSHOCK/kappa.c"
#define PR_KAPPAES "PROBLEMS/BOWSHOCK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BOWSHOCK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BOWSHOCK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BOWSHOCK/dump.c"
#define PR_TOOLS "PROBLEMS/BOWSHOCK/tools.c"

#endif


#if(PROBLEM==29)

#define PR_DEFINE "PROBLEMS/RADWALL/define.h"
#define PR_BC "PROBLEMS/RADWALL/bc.c"
#define PR_INIT "PROBLEMS/RADWALL/init.c"
#define PR_KAPPA "PROBLEMS/RADWALL/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADWALL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADWALL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADWALL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADWALL/dump.c"
#define PR_TOOLS "PROBLEMS/RADWALL/tools.c"

#endif



/*********************/
//including problem specific definitions from PROBLEMS/XXX/define.h
/*********************/

#ifdef _OPENMP
#define OMP
#endif

#include PR_DEFINE
