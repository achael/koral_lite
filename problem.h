#pragma once
//KORAL - problem.h
//choice of the problem plus some definitions

// WARNING: old problem definitions are listed but have been removed as they have not been tested on the updated koral_lite code
// Please contact achael@princeton.edu if you would like to request that they be restored to the code

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
//143 FFTESTS -- Force Free tests
//144 FFALFVEN -- hybrid mhd/ff alfven wave test from Komissarov 99
//145 FFTORUSTEST -- test of force-free torus
//146 FFBONDI -- bondi flow for ff test
//147 PUFFY


#define PROBLEM 145

#if(PROBLEM==147)
#define PR_DEFINE "PROBLEMS/PUFFY/define.h"
#define PR_BC "PROBLEMS/PUFFY/bc.c"
#define PR_INIT "PROBLEMS/PUFFY/init.c"
#define PR_PREPINIT "PROBLEMS/PUFFY/prepinit.c"
#define PR_POSTINIT "PROBLEMS/PUFFY/postinit.c"
#define PR_KAPPAES "PROBLEMS/PUFFY/kappaes.c"
#define PR_TOOLS "PROBLEMS/PUFFY/tools.c"
#endif


#if(PROBLEM==146)
#define PR_DEFINE "PROBLEMS/FFBONDI/define.h"
#define PR_BC "PROBLEMS/FFBONDI/bc.c"
#define PR_INIT "PROBLEMS/FFBONDI/init.c"
#define PR_KAPPAES "PROBLEMS/FFBONDI/kappaes.c"
#define PR_TOOLS "PROBLEMS/FFBONDI/tools.c"
#define PR_POSTINIT "PROBLEMS/FFBONDI/postinit.c"
#endif


#if(PROBLEM==145)
#define PR_DEFINE "PROBLEMS/FFTORUSTEST/define.h"
#define PR_BC "PROBLEMS/FFTORUSTEST/bc.c"
#define PR_INIT "PROBLEMS/FFTORUSTEST/init.c"
#define PR_KAPPAES "PROBLEMS/FFTORUSTEST/kappaes.c"
#define PR_TOOLS "PROBLEMS/FFTORUSTEST/tools.c"
#define PR_POSTINIT "PROBLEMS/FFTORUSTEST/postinit.c"
#endif

#if(PROBLEM==144)
#define PR_DEFINE "PROBLEMS/FFALFVEN/define.h"
#define PR_BC "PROBLEMS/FFALFVEN/bc.c"
#define PR_INIT "PROBLEMS/FFALFVEN/init.c"
#define PR_KAPPAES "PROBLEMS/FFALFVEN/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/FFALFVEN/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/FFALFVEN/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/FFALFVEN/dump.c"
#define PR_FINGER "PROBLEMS/FFALFVEN/finger.c"
#define PR_TOOLS "PROBLEMS/FFALFVEN/tools.c"
#define PR_PREPINIT "PROBLEMS/FFALFVEN/prepinit.c"
#define PR_POSTINIT "PROBLEMS/FFALFVEN/postinit.c"
#endif

#if(PROBLEM==143)
#define PR_DEFINE "PROBLEMS/FFTESTS/define.h"
#define PR_BC "PROBLEMS/FFTESTS/bc.c"
#define PR_INIT "PROBLEMS/FFTESTS/init.c"
#define PR_KAPPAES "PROBLEMS/FFTESTS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/FFTESTS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/FFTESTS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/FFTESTS/dump.c"
#define PR_FINGER "PROBLEMS/FFTESTS/finger.c"
#define PR_TOOLS "PROBLEMS/FFTESTS/tools.c"
#define PR_PREPINIT "PROBLEMS/FFTESTS/prepinit.c"
#define PR_POSTINIT "PROBLEMS/FFTESTS/postinit.c"
#endif

#if(PROBLEM==142)
#define PR_DEFINE "PROBLEMS/PARTIALTDE/define.h"
#define PR_BC "PROBLEMS/PARTIALTDE/bc.c"
#define PR_BC_SPECIAL "PROBLEMS/PARTIALTDE/bc_special.c"
#define PR_BC_SPECIAL_LOOP "PROBLEMS/PARTIALTDE/loop_alloc_special.c"
#define PR_INIT "PROBLEMS/PARTIALTDE/init.c"
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

#if(PROBLEM==134)
#define PR_DEFINE "PROBLEMS/FISHMONC/define.h"
#define PR_BC "PROBLEMS/FISHMONC/bc.c"
#define PR_INIT "PROBLEMS/FISHMONC/init.c"
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
#endif

#if(PROBLEM==119)
#define PR_DEFINE "PROBLEMS/ADAFCRITRELEL/define.h"
#define PR_BC "PROBLEMS/ADAFCRITRELEL/bc.c"
#define PR_INIT "PROBLEMS/ADAFCRITRELEL/init.c"
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
#define PR_KAPPAES "PROBLEMS/BATTERYTEST/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BATTERYTEST/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BATTERYTEST/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BATTERYTEST/dump.c"
#define PR_FINGER "PROBLEMS/BATTERYTEST/finger.c"
#define PR_TOOLS "PROBLEMS/BATTERYTEST/tools.c"
#define PR_PREPINIT "PROBLEMS/BATTERYTEST/prepinit.c"
#define PR_POSTINIT "PROBLEMS/BATTERYTEST/postinit.c"
#endif

#if(PROBLEM==115)
#define PR_DEFINE "PROBLEMS/SHOCKELECTRONTEST/define.h"
#define PR_BC "PROBLEMS/SHOCKELECTRONTEST/bc.c"
#define PR_INIT "PROBLEMS/SHOCKELECTRONTEST/init.c"
#endif

#if(PROBLEM==114)
#define PR_DEFINE "PROBLEMS/RELELEXPAND/define.h"
#define PR_BC "PROBLEMS/RELELEXPAND/bc.c"
#define PR_INIT "PROBLEMS/RELELEXPAND/init.c"
#define PR_KAPPAES "PROBLEMS/RELELEXPAND/kappaes.c"
#define PR_DUMP "PROBLEMS/RELELEXPAND/dump.c"
#define PR_TOOLS "PROBLEMS/RELELEXPAND/tools.c"
#define PR_POSTINIT "PROBLEMS/RELELEXPAND/postinit.c"
#endif

#if(PROBLEM==107)
#define PR_DEFINE "PROBLEMS/RELELTEST/define.h"
#define PR_BC "PROBLEMS/RELELTEST/bc.c"
#define PR_INIT "PROBLEMS/RELELTEST/init.c"
#define PR_KAPPAES "PROBLEMS/RELELTEST/kappaes.c"
#define PR_DUMP "PROBLEMS/RELELTEST/dump.c"
#define PR_TOOLS "PROBLEMS/RELELTEST/tools.c"
#define PR_POSTINIT "PROBLEMS/RELELTEST/postinit.c"
#endif


#if(PROBLEM==93)
#define PR_DEFINE "PROBLEMS/KATOTORUS/define.h"
#define PR_BC "PROBLEMS/KATOTORUS/bc.c"
#define PR_INIT "PROBLEMS/KATOTORUS/init.c"
#define PR_POSTINIT "PROBLEMS/KATOTORUS/postinit.c"
#define PR_KAPPAES "PROBLEMS/KATOTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/KATOTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/KATOTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/KATOTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/KATOTORUS/tools.c"
#endif


#if(PROBLEM==87)
#define PR_DEFINE "PROBLEMS/INFDISK/define.h"
#define PR_BC "PROBLEMS/INFDISK/bc.c"
#define PR_BC_SPECIAL "PROBLEMS/INFDISK/bc_special.c"
#define PR_BC_SPECIAL_LOOP "PROBLEMS/INFDISK/loop_alloc_special.c"
#define PR_INIT "PROBLEMS/INFDISK/init.c"
#define PR_KAPPAES "PROBLEMS/INFDISK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/INFDISK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/INFDISK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/INFDISK/dump.c"
#define PR_TOOLS "PROBLEMS/INFDISK/tools.c"
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



/*********************/
//including problem specific definitions from PROBLEMS/XXX/define.h
/*********************/

#ifdef _OPENMP
#define OMP
#endif

#include PR_DEFINE
