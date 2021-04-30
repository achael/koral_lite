/*! \file ko.h
 \brief constants, variables, and function definitions
*/

//small and big numbers
#define SMALL 1.e-80 
#define BIG (1./SMALL) 
#define NUMEPSILON DBL_EPSILON

// physical constants
#define GGG0 (6.674e-8)
#define CCC0 (2.9979246e10)
#define CCC_CGS (CCC0)
#define gTILDA (1.)
#define cTILDA (1.)
#define GGG (GGG0/gTILDA)
#define CCC (CCC0/cTILDA)
#define MSUNCM 147700. // GM_sun/C^2
#define PARSECCGS (3.086e18)

#define K_BOLTZ_CGS (1.3806488e-16)
#define M_ELECTR_CGS (9.1094e-28)
#define M_PROTON_CGS (1.67262158e-24)
#define SIGMA_RAD_CGS (5.670367e-5)
#define H_CGS (6.6260755e-27)
#define SIGMATH_CGS (6.6524e-25)
#define Pi (3.141592654)
#define HFRAC_SUN 0.7381 //Solar H, He, and Metal abundances
#define HEFRAC_SUN 0.2485
#define MFRAC_SUN (1. - HFRAC_SUN - HEFRAC_SUN) //0.0134

//conversions
#define GMC2 (MSUNCM*MASS)
#define GMC3 (GMC2/CCC0)

#define tempCGS2GU(x)     (x)
#define tempGU2CGS(x)     (x)
#define lenCGS2GU(x)      (x/MASSCM)
#define lenGU2CGS(x)      (x*MASSCM)
#define numdensCGS2GU(x)  (x*(MASSCM*MASSCM*MASSCM))
#define numdensGU2CGS(x)  (x/(MASSCM*MASSCM*MASSCM))
#define timeCGS2GU(x)     (x*(CCC/MASSCM))
#define timeGU2CGS(x)     (x*(MASSCM/CCC))
#define velCGS2GU(x)      (x/CCC)
#define velGU2CGS(x)      (x*CCC)
#define rhoCGS2GU(x)      (x*((GGG*MASSCM*MASSCM)/(CCC*CCC)))
#define rhoGU2CGS(x)      (x*((CCC*CCC)/(GGG*MASSCM*MASSCM)))
#define surfdensCGS2GU(x) (x*((GGG*MASSCM)/(CCC*CCC)))
#define surfdensGU2CGS(x) (x*((CCC*CCC)/(GGG*MASSCM)))
#define massCGS2GU(x)     (x*(GGG/(MASSCM*CCC*CCC)))
#define massGU2CGS(x)     (x*((MASSCM*CCC*CCC)/GGG))
#define kappaCGS2GU(x)    (x*((CCC*CCC)/(GGG*MASSCM)))
#define kappaGU2CGS(x)    (x*((GGG*MASSCM)/(CCC*CCC)))
#define endenCGS2GU(x)    (x*((GGG*MASSCM*MASSCM)/(CCC*CCC*CCC*CCC)))
#define endenGU2CGS(x)    (x*((CCC*CCC*CCC*CCC)/(GGG*MASSCM*MASSCM)))
#define heatcoolCGS2GU(x) (endenCGS2GU(timeGU2CGS(x)))
#define heatcoolGU2CGS(x) (endenGU2CGS(timeCGS2GU(x)))
#define fluxCGS2GU(x)     (x*((GGG*MASSCM*MASSCM)/(CCC*CCC*CCC*CCC*CCC)))
#define fluxGU2CGS(x)     (x*((CCC*CCC*CCC*CCC*CCC)/(GGG*MASSCM*MASSCM)))
#define ergCGS2GU(x)      (x*(GGG/(MASSCM*CCC*CCC*CCC*CCC)))
#define ergGU2CGS(x)      (x*((MASSCM*CCC*CCC*CCC*CCC)/GGG))
#define chargeCGS2GU(x)   (x*sqrt(GGG/MASSCM/CCC/CCC/CCC/CCC)*sqrt(1./MASSCM))
#define chargeGU2CGS(x)   (x/sqrt(GGG/MASSCM/CCC/CCC/CCC/CCC)/sqrt(1./MASSCM))
#define crossCGS2GU(x)    (x/(MASSCM*MASSCM))
#define crossGU2CGS(x)    (x*(MASSCM*MASSCM))

//constants in GU
#define K_BOLTZ (K_BOLTZ_CGS* GGG / CCC / CCC / CCC / CCC / MASSCM)
#define M_PROTON massCGS2GU(M_PROTON_CGS)
#define M_ELECTR massCGS2GU(M_ELECTR_CGS)
#define MPE_RATIO (M_PROTON/M_ELECTR)
#define SIGMA_RAD (SIGMA_RAD_CGS * GGG / CCC / CCC / CCC / CCC / CCC * MASSCM * MASSCM)
#define A_RAD (4.*SIGMA_RAD)
#define Z_RATIO (1.0)
#define KAPPA_ES_COEFF (kappaCGS2GU(0.2*(1. + HFRAC)))
#define E_CHARGE chargeCGS2GU(4.80320425e-10)
#define SIGMA_THOMPSON crossCGS2GU(0.665e-24)

/***** external libraries ******/
#include "problem.h"
#include "mnemonics.h"
#include "choices.h"
#include "mdefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>

#ifdef FFT
#include <fftw3.h>
#endif

#ifdef MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

#if defined COORDOUTPUT_HDF5 || defined DUMPS_CONVERT_HDF5 || defined DUMPS_READ_HDF5 || defined DUMPS_WRITE_HDF5 || defined ANAOUT_HDF5
#include "hdf5.h"
#endif

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>  
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_lambert.h>

#ifdef PR_DEFS
#include PR_DEFS
#endif

#define ldouble double
//#define int long long int
#define FTYPE ldouble
#define gSIZE 20 //size of metric arrays = 16 + 1 (gdet) + 3 (dlgdet)

/***** global variables ******/

//physical constants
ldouble tempcgs2gu, tempgu2cgs, lencgs2gu, lengu2cgs, numdenscgs2gu, numdensgu2cgs,
  timecgs2gu, timegu2cgs, velcgs2gu, velgu2cgs, rhocgs2gu, rhogu2cgs, surfdenscgs2gu,
  surfdensgu2cgs, masscgs2gu, massgu2cgs, kappacgs2gu, kappagu2cgs, endencgs2gu, endengu2cgs,
  heatcoolcgs2gu, heatcoolgu2cgs, fluxcgs2gu, fluxgu2cgs, ergcgs2gu, erggu2cgs, chargecgs2gu,
  chargegu2cgs, crosscgs2gu, crossgu2cgs,
  k_boltz_cgs, k_boltz_cgs_inv, k_boltz, k_boltz_inv, m_proton_cgs, m_proton,
  m_electr_cgs, m_electr, mpe_ratio, sigma_rad_cgs, sigma_rad, h_cgs, sigmath_cgs,
  a_rad, z_ratio, e_charge, sigma_thomson,
  fourpi, fourmpi, per_volume_per_time_cgs2gu, sigma_rad_over_pi, four_sigmarad,
  one_over_four_sigmarad, k_over_mecsq, mugas_mp_over_kB, kB_over_mugas_mp, kB_over_mp,
  kB_over_me, one_over_kB_2p70118, one_over_mugas_mp, one_over_mui_mp, one_over_mue_mp,
  kB_over_mui_mp, kB_over_mue_mp,mui_over_mue,
  four_third, five_third, one_third, two_third, log_2p6, one_over_log_2p6;

ldouble global_double1,global_double2; //for general use

//time stepping
ldouble global_time;
ldouble global_dt;
ldouble global_expdt, global_impdt,global_expdt2;
ldouble global_tstepdenmax,global_tstepdenmin;
ldouble start_time, end_time, mid1_time, mid2_time, maxmp_time;
ldouble max_u2ptime,  min_u2ptime, start_u2ptime, end_u2ptime;
int max_u2ptime_loc,min_u2ptime_loc;
ldouble avgtime,dt;
ldouble *avgselftime;
int nstep;
int global_int_slot[NGLOBALINTSLOT];
ldouble max_ws[3],max_dt,ttm1,ttm2,max_ws_ph;
ldouble tstepdenmax,tstepdenmin;
ldouble min_dx,min_dy,min_dz;

//precalculated metric parameters
ldouble rhorizonBL,rISCOBL,rmboundBL,rphotonBL,etaNT;

// jet coordinate specific
#if (MYCOORDS==JETCOORDS)
ldouble hypx1in,hypx1out,hypx1brk;
#endif
#ifdef CYLINDRIFY
ldouble thetaCYL;
ldouble thetaAX;
ldouble x2cyl;
ldouble sinthetaCYL;
ldouble sinthetaAX;
ldouble rmidcyl;
#endif

//tile specific
int TI,TJ,TK; //tile order
int TOI,TOJ,TOK; //indices of the tile origin
int PROCID;
int NPROCS;

//relativistic electron bins
#ifdef RELELECTRONS
ldouble log10binspace; 
ldouble logbinspace, logbinspace_inv;
ldouble binspace, binspace_inv;
ldouble relel_gammas[NRELBIN];
ldouble relel_gammas_e[NRELBIN+1];
ldouble relel_gammas_inv[NRELBIN];
ldouble relel_gammas_e_inv[NRELBIN+1];

#if defined(RELEL_HEAT_FIX_FRAC) && defined(RELEL_HEAT_FIX_INDEX)
ldouble gamma_inj_fixed_ratio;
#endif

#if defined(RELEL_HEAT_FIX_LIMITS) && defined(RELEL_HEAT_FIX_INDEX)
ldouble relel_injpow[NRELBIN];
ldouble xx_relel_inj;
#endif
#endif

/***** mpi-spec ******/

#if !defined(MPI) && !defined(OMP)
#undef NTX
#undef NTY
#undef NTZ
#define NTX 1 //number of tiles in X
#define NTY 1
#define NTZ 1
#endif

#ifndef NTX
#define NTX 1 //number of tiles in X
#endif

#ifndef NTY
#define NTY 1
#endif

#ifndef NTZ
#define NTZ 1
#endif

#ifdef MPI
#ifndef NX
#define NX (TNX/NTX)
#endif 

#ifndef NY
#define NY (TNY/NTY)
#endif
 
#ifndef NZ
#define NZ (TNZ/NTZ)
#endif 

#else //OMP or single core

#ifndef NX
#define NX (TNX)
#endif 

#ifndef NY
#define NY (TNY)
#endif
 
#ifndef NZ
#define NZ (TNZ)
#endif 
#endif

//number of ghost cells
#ifndef NG
#if (INT_ORDER==0)
#define NG 2 
#endif
#if (INT_ORDER==1)
#define NG 2 
#endif
#if (INT_ORDER==2)
#define NG 3
#endif
#endif

//size of 3d arrays
#define SX ((NX)+2*(NG))
#define NGCX (NG)
#define iX(ix) (ix)

#if(NY>1)
#define NGCY (NG)
#define SY ((NY)+2*(NG))
#define iY(iy) (iy)
#else
#define NGCY 0
#define SY 1
#define iY(iy) (0)
#endif

#if(NZ>1)
#define NGCZ (NG)
#define SZ ((NZ)+2*(NG))
#define iZ(iz) (iz)
#else
#define NGCZ 0
#define SZ 1
#define iZ(iz) (0)
#endif

#define SZMET (SZ)
#define NGCZMET (NGCZ)
#define iZMET(iz) iZ(iz)

//for METRICAXISYMMETRIC, store only single phi slice
#ifdef METRICAXISYMMETRIC
#undef SZMET
#undef NGCZMET
#undef iZMET
#define SZMET 1
#define NGCZMET 0
#define iZMET(iz) (0)
#endif

/*****  arrays   *****/

ldouble **msgbufs;
ldouble *scent, *savg;
ldouble *u,*x,*ucent,*xb,*xout, *xbout_xface, *xbout_yface, *xbout_zface,
  *du,*ut1,*ut2,*ut3,*ut4,*ut0,*u_bak_fixup,*p_bak_fixup,
  *u_step1,*u_step2,*gammagas,
  *vischeating,*vischeatingnegebalance,*vischeatingnegibalance,*vischeatingtimesdeltae,
  *u_step3,*u_step4,*ahdx,*ahdy,*ahdz,
  *aradx,*arady,*aradz,
  *cell_tstepden,*cell_dt,
  *dut0,*dut1,*dut2,*uforget,*drt0,*drt1,*drt2,
  *ahdxl,*ahdyl,*ahdzl,*aradxl,
  *aradyl,*aradzl,  *ahdxr,*ahdyr,*ahdzr,*aradxr,*aradyr,*aradzr,
  *p,*pinit,*pproblem1,*pproblem2,*emf,*ptemp1,*pvecpot,*ppostimplicit,
  *ppreexplicit,*upreexplicit,
  *ptm1,*pt0,*px,*py,*pz,*s,*pavg,
  *dxdx_my2out,*dxdx_out2my,
  *g,*gbx,*gby,*gbz,
  *G,*Gbx,*Gby,*Gbz,
  *a_array, *ax_array, *ay_array, *az_array,
  *pbLx,*pbRx,*pbLy,*pbRy,*pbLz,*pbRz,*sbLx,*sbRx,*sbLy,*sbRy,*sbLz,*sbRz,
  *flbx,*flby,*flbz,*flLx,*flRx,*flLy,*flRy,*flLz,*flRz,
  *flbx2,*flby2,*flbz2,*flLx2,*flRx2,*flLy2,*flRy2,*flLz2,*flRz2,
  *gKr;
  //*emuup,*emulo,*emuupbx,*emulobx,*emuupby,*emuloby,*emuupbz,*emulobz,
  //*emuup2,*emulo2,*emuupbx2,*emulobx2,*emuupby2,*emuloby2,*emuupbz2,*emulobz2,
  //*tmuup,*tmulo,*tmuupbx,*tmulobx,*tmuupby,*tmuloby,*tmuupbz,*tmulobz,
  //*tmuup2,*tmulo2,*tmuupbx2,*tmulobx2,*tmuupby2,*tmuloby2,*tmuupbz2,*tmulobz2;

int *cellflag;

//rad-viscosity specific
#ifdef RADIATION
#if (RADVISCOSITY==SHEARVISCOSITY)
ldouble *Rijviscprev,*radvisclasttime,*Rijviscglobal;
#endif
#endif

//arrays for on the go averaging throughout the domain
ldouble sigma_otg[TNX];
ldouble sigma_otg_temp[TNX];
ldouble scaleth_otg[TNX];
ldouble scaleth_otg_temp[TNX];
ldouble axis1_primplus[NV+2][TNX];
ldouble axis2_primplus[NV+2][TNX];
ldouble axis1_primplus_temp[NV+2][TNX];
ldouble axis2_primplus_temp[NV+2][TNX];

ldouble scalars[NSCALARS];
int doingavg;
int doingpostproc;
int doingpostproc_avg;
ldouble Kr_tmp[4][4][4],g_tmp[4][4];
ldouble inputarg[10];
int **gcidx;
FILE *fout1,*fout_scalars,*fout_radprofiles,*fout_fail,*fhandle_problem1;
int nfout1,nfout2;

/*****  loops   ******/

int Nloop_0,Nloop_1,Nloop_2,Nloop_02,Nloop_3,Nloop_4,Nloop_5,Nloop_6;
int **loop_0,**loop_02,**loop_1,**loop_2,**loop_3,**loop_4,**loop_5,**loop_6;

/* loop over all primitives */
#define PLOOP(j) for(j=0;j<NV;j++)
/* loop over all Dimensions; second rank loop */
#define DLOOP(j,k) for(j=0;j<4;j++)for(k=0;k<4;k++)
/* loop over all Dimensions; second rank loop */
#define DLOOPB(i,j,k) for(i=0;i<4;i++)for(j=0;j<4;j++)for(k=0;k<4;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA(j) for(j=0;j<4;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP(j,k) for(j=1;j<4;j++)for(k=1;k<4;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA(j) for(j=1;j<4;j++)

/*****  structures   ******/

struct opacities
{
  ldouble kappaGasAbs;
  ldouble kappaRadAbs;

  ldouble kappaGasNum;
  ldouble kappaRadNum;

  ldouble totEmissivity;

  ldouble kappaGasRoss;
  ldouble kappaRadRoss;
};

struct geometry
{
  int coords;
  int ix,iy,iz;
  int ifacedim; //-1 - cell center, 0 - x face, 1 - y face, 2 - z face
  ldouble xxvec[4];
  ldouble xx,yy,zz;
  ldouble gg[4][5];
  ldouble GG[4][5];
  ldouble gdet;
  ldouble alpha;
  //ldouble tlo[4][4];
  //ldouble tup[4][4];
  //ldouble elo[4][4];
  //ldouble eup[4][4];
  int par; //some parameter to be used by user	
  ldouble gttpert; //perturbed part of g_tt
};

struct struct_of_state
{
  ldouble rho,uint,gamma,pgas,dV,Tgas,Te,Ti,Sgas,STe,STi,ucon[4],ucov[4],entr,betamag,betarad,cs;
  ldouble bcon[4],bcov[4],Bcon[4],bsq,pmag;
  ldouble urfcon[4],urfcov[4],relgamma,Ehat,prad,Nph;
  ldouble kappa,kappaes;
  struct opacities opac;
  ldouble Rij[4][4],Gi[4],Gic[4],Giff[4],Gicff[4],radentr,Trad,TradBB;
  ldouble gammae,gammai;
  ldouble ne,ue,pe;
  ldouble ni,ui,pi;
  ldouble nenth, uenth, penth;
};

// CHIANTI and OPAL opacities definitions
#define ONEOVERLOGTEN 0.43429448190325182765
struct OpTable {
  int NLOGT; 
  int NLOGRHO; 
  ldouble *logTgrid; 
  ldouble *logRhogrid;
  ldouble **table;
};

/*****  function wrappers  ******/

// coordinates
#define get_x(ic,idim) (idim==0 ? x[ic+NG] : \
		       (idim==1 ? x[ic+NG + NX+2*NG] : \
		       (idim==2 ? x[ic+NG + NX+2*NG + NY+2*NG ] : 0.)))

#define get_xb(ic,idim) (idim==0 ? xb[ic+NG] : \
			(idim==1 ? xb[ic+NG + NX+2*NG + 1] : \
			(idim==2 ? xb[ic+NG + NX+2*NG +1 + NY+2*NG +1 ] : 0.)))

// emf
#define get_emf(iv,ix,iy,iz) (emf[iv-1 + (ix)*3 + (iy)*((NX)+1)*3 + (iz)*((NY)+1)*((NX)+1)*3])
#define set_emf(iv,ix,iy,iz,val) emf[iv-1 + (ix)*3 + (iy)*((NX)+1)*3 + (iz)*((NY)+1)*((NX)+1)*3] = val

// flags
#define get_cflag(iflag,ix,iy,iz)  cellflag[iflag + (iX(ix)+(NGCX))*NFLAGS + \
					            (iY(iy)+(NGCY))*(SX)*NFLAGS + \
						    (iZ(iz)+(NGCZ))*(SY)*(SX)*NFLAGS]
				  
#define set_cflag(iflag,ix,iy,iz,val) cellflag[iflag + (iX(ix)+(NGCX))*NFLAGS + \
					               (iY(iy)+(NGCY))*(SX)*NFLAGS + \
						       (iZ(iz)+(NGCZ))*(SY)*(SX)*NFLAGS] = val

// outcoords at centers 
//note as in get_x idim is spatial (0=x, 1=y, 2=z)
#define get_xout(idim,ix,iy,iz) xout[idim + (iX(ix)+(NGCX))*3 + \
				            (iY(iy)+(NGCY))*(SX)*3 + \
					    (iZ(iz)+(NGCZ))*(SY)*(SX)*3]
#define set_xout(idim,ix,iy,iz,val) xout[idim + (iX(ix)+(NGCX))*3 +	\
				                (iY(iy)+(NGCY))*(SX)*3 + \
					        (iZ(iz)+(NGCZ))*(SY)*(SX)*3]=val

//outcoords at faces
#define get_xbout_xface(idim,ix,iy,iz) xbout_xface[idim + (iX(ix)+(NGCX))*3 +	     \
				                          (iY(iy)+(NGCY))*(SX+1)*3 + \
					                  (iZ(iz)+(NGCZ))*(SY)*(SX+1)*3]
#define set_xbout_xface(idim,ix,iy,iz,val) xbout_xface[idim + (iX(ix)+(NGCX))*3 +	\
				                              (iY(iy)+(NGCY))*(SX+1)*3 + \
					                      (iZ(iz)+(NGCZ))*(SY)*(SX+1)*3]=val

#define get_xbout_yface(idim,ix,iy,iz) xbout_yface[idim + (iX(ix)+(NGCX))*3 + \
				                          (iY(iy)+(NGCY))*(SX)*3 + \
					                  (iZ(iz)+(NGCZ))*(SY+1)*(SX)*3]
#define set_xbout_yface(idim,ix,iy,iz,val) xbout_yface[idim + (iX(ix)+(NGCX))*3 +	\
				                              (iY(iy)+(NGCY))*(SX)*3 + \
					                      (iZ(iz)+(NGCZ))*(SY+1)*(SX)*3]=val

#define get_xbout_zface(idim,ix,iy,iz) xbout_zface[idim + (iX(ix)+(NGCX))*3 + \
				                          (iY(iy)+(NGCY))*(SX)*3 + \
					                  (iZ(iz)+(NGCZ))*(SY)*(SX)*3]
#define set_xbout_zface(idim,ix,iy,iz,val) xbout_zface[idim + (iX(ix)+(NGCX))*3 +	\
				                              (iY(iy)+(NGCY))*(SX)*3 + \
					                      (iZ(iz)+(NGCZ))*(SY)*(SX)*3]=val

//primitive and flux arrays
#define get_u(uarr,iv,ix,iy,iz) uarr[iv + (iX(ix)+(NGCX))*NV + \
				          (iY(iy)+(NGCY))*(SX)*NV + \
					  (iZ(iz)+(NGCZ))*(SY)*(SX)*NV]

#define set_u(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+(NGCX))*NV + \
					      (iY(iy)+(NGCY))*(SX)*NV + \
					      (iZ(iz)+(NGCZ))*(SY)*(SX)*NV] = val

#define NVAVG ((NV+NAVGVARS))
#define SXNVAVG ((long long)(SX*NVAVG))
#define SXSYNVAVG ((long long)(SY*SXNVAVG))
#define get_uavg(uarr,iv,ix,iy,iz) (uarr[iv + (iX(ix)+(NGCX))*(NVAVG) + \
					      (iY(iy)+(NGCY))*(SXNVAVG) + \
					      (iZ(iz)+(NGCZ))*(SXSYNVAVG)])
#define set_uavg(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+(NGCX))*(NVAVG) + \
					         (iY(iy)+(NGCY))*(SXNVAVG) + \
						 (iZ(iz)+(NGCZ))*(SXSYNVAVG)]=val
  
#define get_u_scalar(uarr,ix,iy,iz) uarr[(iX(ix)+(NGCX)) + \
					  (iY(iy)+(NGCY))*(SX) + \
					  (iZ(iz)+(NGCZ))*(SY)*(SX)]
  
#define set_u_scalar(uarr,ix,iy,iz,val) uarr[iX(ix)+(NGCX) + \
					    (iY(iy)+(NGCY))*(SX) + \
					    (iZ(iz)+(NGCZ))*(SY)*(SX)] = val
  
#define set_ubx(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+(NGCX))*NV + \
					        (iY(iy)+(NGCY))*(SX+1)*NV + \
					        (iZ(iz)+(NGCZ))*(SY)*(SX+1)*NV] = val
  
#define set_uby(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+(NGCX))*NV + \
					        (iY(iy)+(NGCY))*(SX)*NV + \
					        (iZ(iz)+(NGCZ))*(SY+1)*(SX)*NV] = val
  
#define set_ubz(uarr,iv,ix,iy,iz,val) uarr[iv + (iX(ix)+(NGCX))*NV + \
					        (iY(iy)+(NGCY))*(SX)*NV + \
					        (iZ(iz)+(NGCZ))*(SY)*(SX)*NV] = val

#define get_ub(uarr,iv,ix,iy,iz,idim) (idim==0 ? uarr[iv + (iX(ix)+(NGCX))*NV + \
						           (iY(iy)+(NGCY))*(SX+1)*NV + \
					                   (iZ(iz)+(NGCZ))*(SY)*(SX+1)*NV] : \
				      (idim==1 ? uarr[iv + (iX(ix)+(NGCX))*NV + \
						           (iY(iy)+(NGCY))*(SX)*NV + \
						           (iZ(iz)+(NGCZ))*(SY+1)*(SX)*NV] : \
				      (idim==2 ? uarr[iv + (iX(ix)+(NGCX))*NV + \
						           (iY(iy)+(NGCY))*(SX)*NV + \
						           (iZ(iz)+(NGCZ))*(SY)*(SX)*NV] : \
				      0.)))

#define get_ub_ptr(uarr,iv,ix,iy,iz,idim) (idim==0 ? &uarr[iv + (iX(ix)+(NGCX))*NV + \
							        (iY(iy)+(NGCY))*(SX+1)*NV + \
					                        (iZ(iz)+(NGCZ))*(SY)*(SX+1)*NV] : \
					  (idim==1 ? &uarr[iv + (iX(ix)+(NGCX))*NV + \
							        (iY(iy)+(NGCY))*(SX)*NV + \
							        (iZ(iz)+(NGCZ))*(SY+1)*(SX)*NV] : \
					  (idim==2 ? &uarr[iv + (iX(ix)+(NGCX))*NV + \
							        (iY(iy)+(NGCY))*(SX)*NV + \
							        (iZ(iz)+(NGCZ))*(SY)*(SX)*NV] : \
					  NULL)))
//saved jacobians

#define get_dxdx(dxdxarr,i,j,ix,iy,iz) dxdxarr[i*4+j + (iX(ix)+(NGCX))*16 + \
				                   (iY(iy)+(NGCY))*(SX)*16 + \
					           (iZMET(iz)+(NGCZMET))*(SY)*(SX)*16]

#define set_dxdx(dxdxarr,i,j,ix,iy,iz,val) dxdxarr[i*4+j + (iX(ix)+(NGCX))*16 + \
				                       (iY(iy)+(NGCY))*(SX)*16 + \
					               (iZMET(iz)+(NGCZMET))*(SY)*(SX)*16] = val

//metric specific

#define get_g(uarr,i,j,ix,iy,iz) uarr[i*5+j + (iX(ix)+(NGCX))*gSIZE + \
				              (iY(iy)+(NGCY))*(SX)*gSIZE + \
					      (iZMET(iz)+(NGCZMET))*(SY)*(SX)*gSIZE]

#define set_g(uarr,i,j,ix,iy,iz,val) uarr[i*5+j + (iX(ix)+(NGCX))*gSIZE +	\
				                  (iY(iy)+(NGCY))*(SX)*gSIZE + \
					          (iZMET(iz)+(NGCZMET))*(SY)*(SX)*gSIZE] = val

#define get_T(uarr,i,j,ix,iy,iz) uarr[i*4+j + (iX(ix)+(NGCX))*16 + \
				              (iY(iy)+(NGCY))*(SX)*16 + \
					      (iZMET(iz)+(NGCZMET))*(SY)*(SX)*16]

#define set_T(uarr,i,j,ix,iy,iz,val) uarr[i*4+j + (iX(ix)+(NGCX))*16 +	\
				                  (iY(iy)+(NGCY))*(SX)*16 + \
					          (iZMET(iz)+(NGCZMET))*(SY)*(SX)*16] = val

#define get_Tfull(uarr,i,j,ix,iy,iz) uarr[i*4+j + (iX(ix)+(NGCX))*16 + \
					          (iY(iy)+(NGCY))*(SX)*16 + \
					          (iZ(iz)+(NGCZ))*(SY)*(SX)*16]

#define set_Tfull(uarr,i,j,ix,iy,iz,val) uarr[i*4+j + (iX(ix)+(NGCX))*16 + \
				     	              (iY(iy)+(NGCY))*(SX)*16 + \
					              (iZ(iz)+(NGCZ))*(SY)*(SX)*16] = val

#define get_Tb(uarr,i,j,ix,iy,iz,idim) (idim==0 ? uarr[i*4+j + (iX(ix)+(NGCX))*16 + \
						               (iY(iy)+(NGCY))*(SX+1)*16 + \
					                       (iZMET(iz)+(NGCZMET))*(SY)*(SX+1)*16] : \
				       (idim==1 ? uarr[i*4+j + (iX(ix)+(NGCX))*16 + \
						               (iY(iy)+(NGCY))*(SX)*16 + \
						               (iZMET(iz)+(NGCZMET))*(SY+1)*(SX)*16] : \
				       (idim==2 ? uarr[i*4+j + (iX(ix)+(NGCX))*16 + \
						               (iY(iy)+(NGCY))*(SX)*16 + \
						               (iZMET(iz)+(NGCZMET))*(SY)*(SX)*16] : \
				       0.)))
  
#define get_gb(uarr,i,j,ix,iy,iz,idim) (idim==0 ? uarr[i*5+j + (iX(ix)+(NGCX))*gSIZE + \
						               (iY(iy)+(NGCY))*(SX+1)*gSIZE + \
					                       (iZMET(iz)+(NGCZMET))*(SY)*(SX+1)*gSIZE] : \
				       (idim==1 ? uarr[i*5+j + (iX(ix)+(NGCX))*gSIZE + \
						               (iY(iy)+(NGCY))*(SX)*gSIZE + \
						               (iZMET(iz)+(NGCZMET))*(SY+1)*(SX)*gSIZE] : \
				       (idim==2 ? uarr[i*5+j + (iX(ix)+(NGCX))*gSIZE + \
						               (iY(iy)+(NGCY))*(SX)*gSIZE + \
						               (iZMET(iz)+(NGCZMET))*(SY)*(SX)*gSIZE] : \
				       0.)))
  
#define get_gKr(i,j,k,ix,iy,iz) gKr[i*4*4+j*4+k + (iX(ix)+(NGCX))*64 + \
				                  (iY(iy)+(NGCY))*(SX)*64 + \
						  (iZMET(iz)+(NGCZMET))*(SY)*(SX)*64]
				       
#define set_gKr(i,j,k,ix,iy,iz,val) gKr[i*4*4+j*4+k + (iX(ix)+(NGCX))*64 + \
					              (iY(iy)+(NGCY))*(SX)*64 + \
						      (iZMET(iz)+(NGCZMET))*(SY)*(SX)*64] = val

//other wrappers
#define delta(i,j) (i==j ? 1 : 0)
#define kron(i,j) (i == j ? 1. : 0.)
#define dot(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3])

#ifdef NONRELMHD
#define dotB(A,B) (A[1]*B[1]+A[2]*B[2]+A[3]*B[3]) 
#else
#define dotB(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2]+A[3]*B[3])
#endif

#define dot3(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define dot3nr(A,B) (A[3]*B[3]+A[1]*B[1]+A[2]*B[2])

//some more macros
#define my_max(x,y) (x>y?x:y)
#define my_max3(x,y,z) (x>my_max(y,z)?x:my_max(y,z))
#define mybilinear(x, y, i, j, xx, yy, zz)  (( zz[i][j]*(xx[i+1]-x)*(yy[j+1]-y) + \
					       zz[i+1][j]*(x-xx[i])*(yy[j+1]-y) + \
					       zz[i][j+1]*(xx[i+1]-x)*(y-yy[j]) + \
					       zz[i+1][j+1]*(x-xx[i])*(y-yy[j])) / \
					     ((xx[i+1]-xx[i])*(yy[j+1]-yy[j])) )
				       
#define mylinear(x,i,xx,yy) (yy[i] + (yy[i+1] - yy[i])/(xx[i+1] - xx[i])*(x- xx[i])   )

///////////////////////////////////////////////////////////////
// ko.c ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int main(int argc, char **argv);

///////////////////////////////////////////////////////////////
// problem.c //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int solve_the_problem(ldouble tstart, char* folder);
int my_finger(ldouble t);

///////////////////////////////////////////////////////////////
// misc.c /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void initialize_constants();
int print_scalings();
int set_initial_profile();
void am_i_sane();
int initialize_arrays();
int free_arrays();
void init_pointers();
int fill_arrays_at_init();
int inverse_matrix(ldouble *a, ldouble *ia, int N);
int multiply_44matrices(ldouble T1[][4],ldouble T2[][4],ldouble Tout[][4]);
int inverse_44matrix(ldouble a[][4], ldouble ia[][4]);
ldouble determinant_44matrix(ldouble a[][4]);
ldouble step_function(ldouble x,ldouble x9);
int my_err(char *message);
int my_warning(char *message);
int getch();
ldouble my_min(ldouble a, ldouble b);
ldouble my_min_N(ldouble *v,int N);
ldouble my_max_N(ldouble *v,int N);
ldouble my_sign(ldouble x);
ldouble my_atan2(ldouble y, ldouble x);
void shuffle_loop(int **array, size_t n);
int print_tensor(ldouble T[][4]);
int print_NNtensor(ldouble **T,int N);
int print_metric(ldouble T[][5]);
int print_4vector(ldouble v[4]);
int print_NVvector(ldouble *v);
int print_primitives(ldouble *p);
int print_conserved(ldouble *u);
int print_Nvector(ldouble v[4],int N);
int my_clock_gettime(void* tsptr);
int decompose_vels(ldouble *pp,int velidx, ldouble v[4],void *ggg,  void *gggBL);
int get_cellsize_out(int ix,int iy,int iz,ldouble dx[3]);
int get_cellsize_arb(int ix,int iy,int iz,ldouble dx[3],int COORDS);

int pertsign_3d(int a,int b,int c);
int pertsign_4d(int a,int b,int c, int d);
int epsilon_4d(int a,int b,int c, int d);
ldouble epsilon_4d_tensor(int a,int b,int c, int d, void *ggg);
ldouble epsilon_3d_tensor(int a,int b,int c, void *ggg);
int convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t);
int convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t);

// pointers in misc.c which point to functions in physics.c
ldouble (*calc_SefromrhoT)(ldouble rho, ldouble temp,int type);
ldouble (*calc_Sefromrhou)(ldouble rho, ldouble temp,int type);
ldouble (*calc_Sefromrhop)(ldouble rho, ldouble p,int type);
ldouble (*calc_TfromSerho)(ldouble S2,ldouble rho,int type,int,int,int);
ldouble (*calc_ufromSerho)(ldouble S2,ldouble rho,int type,int,int,int);

///////////////////////////////////////////////////////////////
// mpi.c //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int mpi_exchangedata();

#ifdef MPI
MPI_Group mpi_all_group;
#define MPI_LDOUBLE MPI_DOUBLE
int mpi_recvdata(MPI_Request *reqs, int *nreqs);
int mpi_senddata(MPI_Request *reqs, int *nreqs);
int mpi_savedata();
#endif

int mpi_hasBC();
int mpi_isitBC(int BCtype);
int mpi_isitBC_forcorners(int BCtype);
void mpi_synchtiming(ldouble *time);
void mpi_myinit(int argc, char *argv[]);
void mpi_myfinalize();
//void mpi_tileorigin(long long ti, long long tj, long long tk, long long* toi, long long* toj, long long* tok);
//void mpi_global2localidx(long long gix, long long giy, long long giz, long long *lix, long long *liy, long long *liz);
//void mpi_local2globalidx(long long lix, long long liy, long long liz, long long *gix, long long *giy, long long *giz);

void mpi_tileorigin(int ti, int tj, int tk, int* toi, int* toj, int* tok);
void mpi_global2localidx(int gix,int giy, int giz, int *lix, int *liy, int *liz);
void mpi_local2globalidx(int lix,int liy, int liz, int *gix, int *giy, int *giz);

void mpi_procid2tile(int procid, int* tilei, int* tilej, int* tilek);
int mpi_tile2procid(int tilei, int tilej, int tilek);
int omp_myinit();
int calc_avgs_throughout();

///////////////////////////////////////////////////////////////
// metric.c ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
int calc_metric();
ldouble r_horizon_BL(ldouble a);
ldouble r_ISCO_BL(ldouble ac);
ldouble r_mbound_BL(ldouble a);
ldouble r_photon_BL(ldouble a);

ldouble calc_gdet(ldouble *xx);
ldouble calc_gdet_arb(ldouble *xx,int COORDS);
ldouble calc_gdet_arb_ana(ldouble *xx,int coords);
ldouble calc_dlgdet(ldouble *xx, int idim);
ldouble calc_dlgdet_arb(ldouble *xx, int idim,int coords);
int calc_g(ldouble *xx, ldouble g[][5]);
int calc_g_arb(ldouble *xx, ldouble g[][5],int COORDS);
int calc_g_arb_ana(ldouble *xx, ldouble g[][5],int coords);
int calc_G(ldouble *xx, ldouble G[][5]);
int calc_G_arb(ldouble *xx, ldouble G[][5],int COORDS);
int calc_G_arb_ana(ldouble *xx, ldouble G[][5],int coords);
int calc_Krzysie(ldouble *xx, ldouble Krzys[][4][4]);
int calc_Krzysie_arb(ldouble *xx, ldouble Krzys[][4][4],int COORDS);
int calc_Krzysie_at_center(int ix,int iy,int iz, ldouble Krzys[][4][4]);
int calc_Krzysie_arb_ana(ldouble *xx, ldouble Krzys[][4][4],int coords);

int fill_geometry(int ix,int iy,int iz,void *geom);
int fill_geometry_arb(int ix,int iy,int iz,void *geom,int COORDS);
int fill_geometry_face(int ix,int iy,int iz,int idim, void *geom);
int fill_geometry_face_arb(int ix,int iy,int iz,int idim, void *geom,int COORDS);

//deprecated
//int calc_tetrades(ldouble g[][5], ldouble tmuup[][4], ldouble tmulo[][4],int coords);
//int calc_ZAMOes(ldouble g[][5], ldouble emuup[][4], ldouble emulo[][4], int coords);

int coco_N(ldouble *x1, ldouble *x2,int CO1, int CO2);
int coco_BL2KS(ldouble *xBL, ldouble *xKS);
int coco_KS2BL(ldouble *xKS, ldouble *xBL);
int coco_KS2MKS1(ldouble *xKS, ldouble *xMKS1);
int coco_KS2MKS2(ldouble *xKS, ldouble *xMKS1);
int coco_KS2MKS3(ldouble *xKS, ldouble *xMKS);
int coco_KS2TKS3(ldouble *xKS, ldouble *xMKS);
int coco_MINK2TFLAT(ldouble *xKS, ldouble *xMKS);
int coco_MKS12KS(ldouble *xMKS1, ldouble *xKS);
int coco_MKS22KS(ldouble *xMKS1, ldouble *xKS);
int coco_MKS32KS(ldouble *xMKS, ldouble *xKS);
int coco_TKS32KS(ldouble *xMKS, ldouble *xKS);
int coco_TFLAT2MINK(ldouble *xMKS, ldouble *xKS);
int coco_MCYL12CYL(ldouble *xMCYL1, ldouble *xCYL);
int coco_CYL2MCYL1(ldouble *xCYL, ldouble *xMCYL1);
int coco_CYL2SPH(ldouble *xCYL, ldouble *xSPH);
int coco_SPH2CYL(ldouble *xCYL, ldouble *xSPH);
int coco_MSPH12SPH(ldouble *xMSPH1, ldouble *xSPH);
int coco_SPH2MSPH1(ldouble *xSPH, ldouble *xMSPH1);
int coco_MKER12KER(ldouble *xMKER1, ldouble *xKER);
int coco_KER2MKER1(ldouble *xKER, ldouble *xMKER1);
int coco_SPH2MINK(ldouble *xSPH, ldouble *xMINK);
int coco_MINK2SPH(ldouble *xMINK, ldouble *xSPH);
int coco_JET2KS(ldouble *xJET, ldouble *xKS);
int coco_KS2JET(ldouble *xKS, ldouble *xJET);
  
int dxdx_BL2KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KS2BL(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KS2MKS1(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KS2MKS3(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KS2TKS3(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MINK2TFLAT(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KS2MKS2(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MKS12KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MKS22KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MKS32KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_TKS32KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_TFLAT2MINK(ldouble *xx, ldouble dxdx[][4]);
int dxdx_CYL2MCYL1(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MCYL12CYL(ldouble *xx, ldouble dxdx[][4]);
int dxdx_CYL2SPH(ldouble *xx, ldouble dxdx[][4]);
int dxdx_SPH2CYL(ldouble *xx, ldouble dxdx[][4]);
int dxdx_SPH2MSPH1(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MSPH12SPH(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KER2MKER1(ldouble *xx, ldouble dxdx[][4]);
int dxdx_MKER12KER(ldouble *xx, ldouble dxdx[][4]);
int dxdx_JET2KS(ldouble *xx, ldouble dxdx[][4]);
int dxdx_KS2JET(ldouble *xx, ldouble dxdx[][4]);
ldouble dxdxF(ldouble x, void* params);
int calc_dxdx_arb_num(ldouble *xx, ldouble dxdx[][4], int CO1, int CO2);
int calc_dxdx_arb(ldouble *xx, ldouble dxdx[][4], int CO1, int CO2);

int calc_g_arb_num(ldouble *xx, ldouble gout[][5],int COORDS);
int calc_G_arb_num(ldouble *xx, ldouble Gout[][5],int COORDS);
ldouble fg (double x, void * params);
int calc_Krzysie_arb_num(ldouble *xx, ldouble Krzys[][4][4],int COORDS);
ldouble calc_gdet_arb_num(ldouble *xx,int COORDS);
ldouble calc_gttpert(ldouble *xx);
ldouble calc_gttpert_arb(double *xx, int COORDS);

int if_coords_sphericallike(int coords);
int if_coords_cylindricallike(int coords);
int print_Krzysie(ldouble g[][4][4]);
int print_g(ldouble g[][5]);
int test_metric();

ldouble psi_smooth(ldouble x);
ldouble theta_smooth(ldouble x);
ldouble minn(ldouble a, ldouble b, ldouble df);
ldouble maxx(ldouble a, ldouble b, ldouble df);
ldouble wjet(ldouble x2, ldouble fdisk, ldouble fjet);
ldouble theta_disk_or_jet(ldouble r, ldouble x2, ldouble rdecoll, ldouble rcoll, ldouble runi, ldouble a1, ldouble a2);
ldouble theta_diskjet(ldouble r, ldouble x2, void *params);
ldouble jetcoords_theta(ldouble r, ldouble x2, void *params);
ldouble jetcoords_theta_root(ldouble x2, void *params);
ldouble jetcoords_theta_inv(ldouble r, ldouble theta, void *params);
ldouble hyperexp_func(ldouble x1, void *params);
ldouble hyperexp_func_root(ldouble x1, void *params);
ldouble hyperexp_func_inv(ldouble r, void *params);
ldouble hyperexp_x1max(ldouble rmax, ldouble rbrk, ldouble r0);

ldouble to1stquad(ldouble  x2);
ldouble calc_theta0(ldouble rcyl, ldouble x2cyl);

ldouble sinth0(ldouble r, ldouble x2, void* params);
ldouble sinth1(ldouble r, ldouble x2, void* params);
ldouble sinth2(ldouble r, ldouble x2, void* params);

int set_cyl_params();
//ldouble sinth2(ldouble r, ldouble theta, ldouble theta2); 
ldouble sinth2(ldouble r, ldouble x2, void* params); 

//ldouble f2func(ldouble r, ldouble x2, ldouble theta, void* params);
ldouble f2func(ldouble r, ldouble x2, void* params);
ldouble cylindrify(ldouble r, ldouble x2, void* params);
///////////////////////////////////////////////////////////////
// frames.c ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int trans_pall_coco(ldouble *pp1, ldouble *pp2, int CO1,int CO2, ldouble *xxvec, void* ggg1, void* ggg2);
int trans_pmhd_coco(ldouble *ppin, ldouble *ppout, int CO1,int CO2, ldouble *xxvec, void* ggg1,void* ggg2);
int trans_prad_coco(ldouble *ppin, ldouble *ppout, int CO1,int CO2, ldouble *xxvec, void* ggg1, void* ggg2);

//deprecated
int prad_ff2lab(ldouble *pp1, ldouble *pp2, void* ggg);
//int prad_lab2ff(ldouble *pp1, ldouble *pp2, void *ggg);
//int prad_on2lab(ldouble *pp1, ldouble *pp2, void* ggg);
//int prad_lab2on(ldouble *pp1, ldouble *pp2, void *ggg);
//int prad_ff2zamo(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble GG[][5], ldouble eup[][4]);
//int prad_zamo2ff(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble GG[][5], ldouble eup[][4]);
//int boost2_zamo2ff(ldouble* A1,ldouble* A2,ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4]);
//int boost2_ff2zamo(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4]);
//int boost22_zamo2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4]);
//int boost22_ff2zamo(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4]);
//int trans22_zamo2lab(ldouble T1[][4],ldouble T2[][4],ldouble elo[][4]);
//int trans22_lab2zamo(ldouble T1[][4],ldouble T2[][4],ldouble eup[][4]);
//int trans2_lab2zamo(ldouble *u1,ldouble *u2,ldouble eup[][4]);
//int trans2_zamo2lab(ldouble *u1,ldouble *u2,ldouble elo[][4]);
//int trans22_on2cc(ldouble T1[][4],ldouble T2[][4],ldouble tlo[][4]);
//int trans22_cc2on(ldouble T1[][4],ldouble T2[][4],ldouble tup[][4]);
//int trans2_cc2on(ldouble *u1,ldouble *u2,ldouble tup[][4]);
//int trans2_on2cc(ldouble *u1,ldouble *u2,ldouble tlo[][4]);

int calc_Lorentz_lab2ff(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble L[][4]);
int calc_Lorentz_lab2ff_4vel(ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble L[][4], ldouble ucon[4], ldouble ucov[4]);
int calc_Lorentz_ff2lab(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble L[][4]);

int boost22_lab2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5]);
int boost22_ff2lab(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5]);
int boost22_ff2lab_with_alpha(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5], ldouble alpha);
int boost22_rf2lab(ldouble T1[][4],ldouble T2[][4],ldouble *pp0,ldouble gg[][5],ldouble GG[][5]);
int boost22_lab2rf(ldouble T1[][4],ldouble T2[][4],ldouble *pp0,ldouble gg[][5],ldouble GG[][5]);
int boost2_lab2ff(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5]);
int boost2_lab2ff_4vel(ldouble A1[4], ldouble A2[4], ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble ucon[4], ldouble ucov[4]);
int boost2_lab2rf(ldouble A1[4],ldouble A2[4],ldouble *pp0,ldouble gg[][5],ldouble GG[][5]);
int boost2_ff2lab(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5]);

int multiply11(ldouble T1[][4],ldouble T2[][4],ldouble A[][4]);
int multiply22(ldouble T1[][4],ldouble T2[][4],ldouble A[][4]);
int multiply2(ldouble *u1,ldouble *u2,ldouble A[][4]);

int coco_3vector(ldouble A1[3],ldouble A2[3],int CO1,int CO2,void* ggg);
int trans2_coco(ldouble *xx,ldouble *u1,ldouble *u2,int CO1, int CO2);
int trans22_coco(ldouble *xx,ldouble T1[][4],ldouble T2[][4],int CO1, int CO2);

int indices_1122(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5]);
int indices_2211(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5]);
int indices_2122(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5]);
int indices_1121(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5]);
int indices_2221(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5]);
int indices_12(ldouble A1[4],ldouble A2[4],ldouble GG[][5]);
int indices_21(ldouble A1[4],ldouble A2[4],ldouble gg[][5]);

int trans_pall_coco_my2out(ldouble *pp1, ldouble *pp2, void* ggg1, void* ggg2) ;
int trans_pall_coco_out2my(ldouble *pp1, ldouble *pp2, void* ggg1, void* ggg2) ;

int trans_pmhd_coco_my2out(ldouble *ppin, ldouble *ppout, void* ggg1,void* ggg2);
int trans_prad_coco_my2out(ldouble *ppin, ldouble *ppout, void* ggg1, void* ggg2);
int trans_pmhd_coco_out2my(ldouble *ppin, ldouble *ppout, void* ggg1,void* ggg2);
int trans_prad_coco_out2my(ldouble *ppin, ldouble *ppout, void* ggg1, void* ggg2);
int trans_pmhd_coco_precompute(ldouble *ppin, ldouble *ppout, void* ggg1,void* ggg2, int which);
int trans_prad_coco_precompute(ldouble *ppin, ldouble *ppout, void* ggg1, void* ggg2, int which);

int trans2_coco_my2out(ldouble *u1, ldouble *u2, int ix, int iy, int iz);
int trans22_coco_my2out(ldouble T1[][4], ldouble T2[][4], int ix, int iy, int iz);
int trans2_coco_out2my(ldouble *u1, ldouble *u2, int ix, int iy, int iz);
int trans22_coco_out2my(ldouble T1[][4], ldouble T2[][4], int ix, int iy, int iz);
int trans2_coco_precompute(ldouble *u1, ldouble *u2, int ix, int iy, int iz, int which);
int trans22_coco_precompute(ldouble T1[][4], ldouble T2[][4], int ix, int iy, int iz, int which);

///////////////////////////////////////////////////////////////
// relele.c ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int calc_ucon_ucov_from_prims(ldouble *pr, void *ggg, ldouble *ucon, ldouble *ucov);
int calc_urcon_urcov_from_prims(ldouble *pr, void *ggg, ldouble *urcon, ldouble *urcov);
int conv_vels(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
int conv_vels_ut(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
int conv_vels_both(ldouble *u1,ldouble *u2con,ldouble *u2cov,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
int conv_vels_core(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5],int);
ldouble calc_alpgam(ldouble *u1, ldouble gg[][5], ldouble GG[][5]);

int fill_utinvel3(ldouble *u1,double gg[][5],ldouble GG[][5]);
int fill_utinucon(ldouble *u1,double gg[][5],ldouble GG[][5]);
int conv_velsinprims(ldouble *pp,int which1, int which2,ldouble gg[][5],ldouble GG[][5]);

//TODO: remove these? 
int calc_normalobs_ncon(ldouble GG[][5], ldouble alpha, ldouble *ncon);
int calc_normalobs_ncov_ncon(ldouble GG[][5], ldouble alpha, ldouble *ncov, ldouble *ncon);
int calc_normalobs_4vel(ldouble GG[][5], ldouble *ncon);
int calc_normalobs_relvel(ldouble GG[][5], ldouble *ncon);
//

int set_hdatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype);
int set_radatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype);

int pick_Tb(ldouble *arr,int ix,int iy,int iz,int,ldouble T[][4]);
int pick_T(ldouble *arr,int ix,int iy,int iz,ldouble T[][4]);
int pick_g(int ix,int iy,int iz,ldouble gg[][5]);
ldouble pick_gdet(int ix,int iy,int iz);
int pick_G(int ix,int iy,int iz,ldouble gg[][5]);
int pick_gb(int ix,int iy,int iz,int,ldouble gg[][5]);
int pick_Gb(int ix,int iy,int iz,int,ldouble gg[][5]);

int print_p(ldouble *p);
int print_u(ldouble *p);


///////////////////////////////////////////////////////////////
// p2u.c //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int p2u(ldouble *p, ldouble *u, void *ggg);
int p2u_mhd(ldouble *p, ldouble *u, void *ggg);
ldouble calc_utp1(ldouble *vcon, ldouble *ucon, void *ggg);
int p2u_mhd_nonrel(ldouble *p, ldouble *u, void *ggg);
int p2u_rad(ldouble *pp,ldouble *uu,void *ggg);
int p2avg(int ix,int iy,int iz,ldouble *avg);
int test_maginv();

///////////////////////////////////////////////////////////////
// u2p.c //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int calc_primitives(int ix,int iy,int iz,int type,int setflags);
int u2p(ldouble *uu0, ldouble *pp,void *ggg,int corrected[3],int fixups[2],int type);
int check_floors_mhd(ldouble *pp, int whichvel,void *ggg);

int u2p_solver(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose);
int u2p_solver_nonrel(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose);
int u2p_solver_Wp(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose);
int u2p_solver_W(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose);
int u2p_solver_Bonly(ldouble *uu, ldouble *pp, void *ggg);

int count_entropy(int *n, int *n2);
int copy_entropycount();
int update_entropy();

int test_inversion();
int test_inversion_nonrel();
int test_inversion_5d();

///////////////////////////////////////////////////////////////
// u2prad.c ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int u2p_rad(ldouble *uu, ldouble *pp, void *ggg, int *corrected);
int u2p_rad_urf(ldouble *uu, ldouble *pp,void* ggg, int *corrected);
int check_floors_rad(ldouble *pp, int whichvel,void *ggg);

///////////////////////////////////////////////////////////////
// finite.c ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int reduce_order_check(ldouble *pm2,ldouble *pm1,ldouble *p0,
		       ldouble *pp1,ldouble *pp2,
		       int ix,int iy,int iz);
ldouble reduce_minmod_theta(ldouble *pm2,ldouble *pm1,ldouble *p0,
			    ldouble *pp1,ldouble *pp2,
			    int ix,int iy,int iz);
int avg2point(ldouble *um2,ldouble *um1,ldouble *u0,ldouble *up1,ldouble *up2,
	      ldouble *ul,ldouble *ur,ldouble dxm2,ldouble dxm1,ldouble dx0,
	      ldouble dxp1,ldouble dxp2,int param,ldouble theta);
int save_wavespeeds(int ix,int iy,int iz, ldouble *aaa);
int save_timesteps();
int calc_u2p(int type,int setflags);
int calc_wavespeeds();
int do_correct();

int op_explicit(ldouble t, ldouble dtin);
int op_intermediate(ldouble t, ldouble dt);
int apply_dynamo(ldouble t, ldouble dt);
int op_implicit(ldouble t, ldouble dtin);

ldouble f_calc_fluxes_at_faces(int ix,int iy,int iz);

int set_grid(ldouble *mindx,ldouble *mindy, ldouble *mindz, ldouble *maxdtfac);
int set_grid_outcoords();

int alloc_loops();
int print_grid(ldouble min_dx, ldouble min_dy, ldouble min_dz);

ldouble get_size_x(int ic, int idim);
int get_xx(int ix,int iy,int iz,ldouble *xx);
int get_xxout(int ix,int iy,int iz,ldouble *xx);
int get_xx_arb(int ix,int iy,int iz,ldouble *xx,int COORDSOUT);
int set_x(int ic, int idim,ldouble val);
int set_xb(int ic, int idim,ldouble val);
ldouble calc_xb(int i,int idim);
int calc_bc(int ix,int iy,int iz,ldouble t, ldouble *uu,ldouble *pp,int ifinit,int BCtype);
//int set_g(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value);
//int set_T(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value);
//int set_Tfull(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value);
int set_ub(ldouble* uarr,int iv,int ix,int iy,int iz,ldouble value,int idim);
int set_gb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim);
int set_Tb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim);

int copy_u_core(ldouble factor,ldouble *uu1,ldouble* uu2, long long N);
int copy_u(ldouble factor,ldouble *uu1,ldouble* uu2 );
int copyi_u(ldouble factor,ldouble *uu1,ldouble* uu2);
int add_u_core(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble *uu3, long long N);
int add_u(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble *uu3);
int addi_u(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble *uu3);
int add_u_core_3(ldouble f1, ldouble* uu1,
		 ldouble f2, ldouble *uu2,
		 ldouble f3, ldouble *uu3,
		 ldouble *uu4,long long N);
int add_u_3(ldouble f1, ldouble* uu1,
	    ldouble f2, ldouble *uu2,
	    ldouble f3, ldouble *uu3,
	    ldouble *uu4);
int addi_u_3(ldouble f1, ldouble* uu1,
	     ldouble f2, ldouble *uu2,
	     ldouble f3, ldouble *uu3,
	     ldouble *uu4);

int if_indomain(int ix,int iy,int iz);
int if_outsidegc(int ix,int iy,int iz);

int set_bc_core(int ix,int iy,int iz,double t,ldouble *uval,ldouble *pval,int ifinit,int BCtype);
int set_bc(ldouble t,int ifinit);

int cell_fixup(int type);

int f_implicit_metric(const gsl_vector * x, void *paramsp, gsl_vector * f);
int print_state_metric (int iter, gsl_multiroot_fsolver * s);
int solve_implicit_metric(int ix,int iy,int iz,ldouble dt,ldouble *ubase);

int smooth_polaraxis();
int correct_nssurface();
int correct_polaraxis();
int correct_polaraxis_3d();
int is_cell_corrected_polaraxis(int ix, int iy, int iz);
int is_cell_active(int ix, int iy, int iz);
int skip_cell_implicit(int ix, int iy, int iz);

int get_factors_entropies_following_gas(int ix,int iy,int iz,ldouble *f0,
					ldouble *fxl,ldouble *fxr,
					ldouble* fyl, ldouble *fyr,
					ldouble *fzl, ldouble *fzr,
					ldouble dt, int iv);
int mix_entropies(ldouble dt);


///////////////////////////////////////////////////////////////
// magn.c /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

void calc_bcon_bcov_bsq_from_4vel(ldouble *pr, ldouble *ucon, ldouble *ucov, void* ggg,
				  ldouble *bcon, ldouble *bcov, double *bsq);
void calc_bcon_4vel(double *pr, double *ucon, double *ucov, double *bcon);
void calc_Bcon_4vel(double *pr, double *ucon, double *bcon, double *Bcon);
void calc_bcon_prim(double *pp, double *bcon, void* ggg);
void calc_Bcon_prim(double *pp, double *bcon,double *Bcon, void* ggg);

int fl_x(int i);
int fl_y(int i);
int fl_z(int i);
int flx_x(int i);
int flx_y(int i);
int flx_z(int i);
int fly_x(int i);
int fly_y(int i);
int fly_z(int i);
int flz_x(int i);
int flz_y(int i);
int flz_z(int i);
int flux_ct();
int adjust_fluxcttoth_emfs();

int calc_BfromA(ldouble* pinput, int ifoverwrite);
int calc_BfromA_core();
ldouble calc_divB(int ix,int iy,int iz);

int calc_Qthetaphi(int ix, int iy, int iz,ldouble *Qtheta,ldouble *Qphi);
int calc_angle_brbphibsq(int ix, int iy, int iz, ldouble *brbphi, ldouble *bsq,
			 ldouble *bcon,ldouble *bcov);
int calc_angle_bpbphibsq(int ix, int iy, int iz, ldouble *bpbphi, ldouble *bsq,
			 ldouble *bcon, ldouble *bcov);
int calc_curl(ldouble *ptu, int ix, int iy, int iz, void* ggg, ldouble *curl);
int mimic_dynamo(ldouble dtin);

///////////////////////////////////////////////////////////////
// rad.c //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int implicit_lab_rad_source_term(int ix,int iy, int iz,ldouble dt);
int solve_implicit_lab(int ix,int iy,int iz,ldouble dt,int verbose);
int solve_implicit_lab_4dprim(ldouble *uu00,ldouble *pp00,void *ggg,
			      ldouble dt,int verbose,int *params,ldouble *pp);
int implicit_apply_constraints(ldouble *pp, ldouble *uu, ldouble *uu0,  void* ggg, int whichprim);
int f_implicit_lab_4dprim_with_state(ldouble *uuin, ldouble *ppin, void *sss,
				     ldouble *uu0, ldouble *pp0, void *sss0,
				     ldouble dt, void* ggg,
				     ldouble *f, int *params, ldouble *err0);
int print_state_implicit_lab_4dprim (int iter, ldouble *x, ldouble *f,ldouble err,int N);
int free_solve_implicit_lab_4dprim(ldouble** J, ldouble** iJ, ldouble *tJ, ldouble *tiJ,
				   ldouble *f1, ldouble *f2, ldouble *f3,
				   ldouble *xxx, ldouble *xxxbest,int N);

int solve_implicit_lab_1dprim(ldouble *uu0, ldouble *pp0, void *sss0, void *ggg,
			      ldouble dt,  int verbose, ldouble *ppout, void *sss);
ldouble f_implicit_1dprim_err(ldouble xxx, ldouble *uu0, ldouble *pp0, void *sss0, void *sss,
			      ldouble dtau, void *ggg, int *params,
			      ldouble totenergy, ldouble ratio, int verbose);
int solve_implicit_nphoton_rad_1dprim(ldouble *uu0, ldouble *pp0, void *sss, void *ggg,
				      ldouble dt, int verbose, ldouble *ppout, void *sssout);
ldouble f_implicit_photon_rad_1dprim_err(ldouble xxx, ldouble *uu0, ldouble *pp0, void *sss,
					 ldouble dtaurad, void *ggg);

int explicit_rad_source_term(int ix,int iy, int iz,ldouble dt);
int solve_explicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int verbose);
int solve_explicit_lab_core(ldouble *uu,ldouble *pp,void* ggg,ldouble dt,ldouble* deltas,int verbose);
int apply_rad_source_del4(int ix,int iy,int iz,ldouble *del4);

ldouble calc_Gi(ldouble *pp, void *ggg, ldouble Gi[4], ldouble relel_dudtau, int type, int reltype);
ldouble calc_Gi_with_state(ldouble *pp, void *sss, void *ggg,
			   ldouble Gi[4], ldouble relel_dudtau, int type, int reltype);
ldouble calc_all_Gi_with_state(ldouble *pp, void *sss, void* ggg,
			       ldouble Gitot_ff[4], ldouble Gitot_lab[4],
			       ldouble Gith_ff[4], ldouble Gith_lab[4],
			       ldouble relel_dudtau, int reltype);
ldouble calc_Gi_nonrel_with_state(ldouble *pp, void *sss, void *ggg, ldouble Gi[4],int labframe);
int calc_Compt_Gi(ldouble *pp, void* ggg, ldouble *Gic,
		  ldouble Ehatrad, ldouble Te, ldouble kappaes, ldouble *ucon);
int calc_Compt_Gi_with_state(ldouble *pp, void *sss, void* ggg, ldouble *Gic, ldouble *ucon_frame);

ldouble calc_CoulombCoupling(ldouble *pp,void *ggg);
ldouble calc_CoulombCoupling_with_state(ldouble *pp,void *sss,void *ggg);

void calc_Ehat_from_Rij_ucov(double Rij[4][4], double uffcov[4], ldouble *Ehat);
int calc_Rij(ldouble *pp, void* ggg, ldouble Rij[][4]);
int calc_Rij_M1_ff(ldouble *pp, ldouble Rij[][4]);
int calc_Rij_M1_from_4vel(ldouble *pp, void* ggg, ldouble *urfcon, ldouble Rij[][4]);
int calc_Rij_M1(ldouble *pp, void* ggg, ldouble Rij[][4]);

ldouble calc_Tnfromn(ldouble n);
ldouble calc_NFfromT(ldouble T);
ldouble calc_NFfromE(ldouble E);
ldouble calc_LTE_EfromT(ldouble T);
ldouble calc_LTE_TfromE(ldouble E );
ldouble calc_LTE_Efromurho(ldouble u,ldouble rho);
int calc_ff_Rtt(ldouble *pp,ldouble *Rttret, ldouble* ucon,void* ggg);
int calc_ff_Ehat(ldouble *pp,ldouble *Ehat, ldouble* ucon,void* ggg);

ldouble calc_nsource(ldouble *pp, void* ggg);
ldouble calc_nsource_with_state(ldouble *pp, void *sss, void* ggg);
ldouble calc_ncompt_Thatrad_4vel(ldouble *pp, void* ggg, ldouble Ehatrad,
				 ldouble* urfcon, ldouble* uffcov);
ldouble calc_ncompt_Ehatrad(ldouble Tradhat, ldouble nphhat);
ldouble calc_ncompt_Thatrad_fromEN(ldouble Ehat, ldouble nphhat);
ldouble calc_ncompt_Thatrad_full(ldouble *pp, void* ggg);
ldouble calc_ncompt_Thatrad(ldouble *pp, void* ggg, ldouble Ehatrad);
ldouble calc_ncompt_Thatrad_4vel(ldouble *pp, void* ggg, ldouble Ehatrad,
				 ldouble* urfcon, ldouble* uffcov);
ldouble calc_ncompt_Thatrad_nphhat(ldouble nphhat, ldouble Ehatrad);
ldouble calc_ncompt_nphlab(ldouble *pp, void* ggg);
ldouble calc_ncompt_nphhat(ldouble *pp, void* ggg);

int calc_rad_wavespeeds(ldouble *pp,void *ggg,ldouble tautot[3],ldouble *aval,int verbose);
int f_flux_prime_rad_total(ldouble *pp, void *ggg,
			   ldouble Rij[][4],ldouble Rij0[][4], ldouble Rijvisc[][4]);
int f_flux_prime_rad( ldouble *pp, int idim, void *ggg,ldouble *ff);

int calc_rad_shearviscosity(ldouble *pp,void* ggg,ldouble shear[][4],ldouble *nuret,int *derdir);
int calc_shear_lab(ldouble *pp0, void* ggg,ldouble S[][4],ldouble *div, int hdorrad,int *derdir);
int calc_fluid_div_lab(ldouble *pp0, void* ggg, ldouble dt, ldouble *div, int hdorrad, int *derdir);
int calc_rad_visccoeff(ldouble *pp,void *ggg,ldouble *nuret,ldouble *mfpret,ldouble *mindxret);
int calc_Rij_visc_total();
int calc_Rij_visc(ldouble *pp, void* ggg, ldouble Rvisc[][4], int *derdir);
void reset_radviscaccel();

int estimate_Bgrowth_battery(int ix,int iy,int iz,ldouble dBdt[4]);
int calc_batteryflux(ldouble *pp, void* ggg, ldouble *eterm,int idim,ldouble *ucov);
int calc_Efield_battery(ldouble *pp,void *ggg,ldouble econ[4]);

ldouble estimate_gas_radiation_coupling(ldouble *pp, void *ggg);

int test_solve_implicit_lab();
int test_Gi();
int test_solve_implicit_lab_file();
int test_jon_solve_implicit_lab();
int test_Ccoupling();
int test_heatfit();
int test_opacities();
int test_calckappa();
int test_Tsynch();


///////////////////////////////////////////////////////////////
//physics.c
///////////////////////////////////////////////////////////////

int fill_struct_of_state(ldouble *pp,void* ggg,void* sss);
int copy_state_to_state(void *sss1, void *sss2);
int update_state_for_nphoton(ldouble *pp, void *ggg, void *sss);
ldouble calc_thermal_ne(ldouble *pp);

int calc_wavespeeds_lr(int ix, int iy, int iz,ldouble *aaa);
int calc_wavespeeds_lr_pure(ldouble *pp,void *ggg,ldouble *aaa);
int calc_wavespeeds_lr_core(ldouble *ucon, ldouble GG[][5], ldouble *aret,
			    ldouble wspeed2x, ldouble wspeed2y, ldouble wspeed2z);

int f_metric_source_term_arb(ldouble *pp,void *ggg,ldouble *ss);
int f_general_source_term_arb(ldouble *pp,void *ggg,ldouble *ss);
int f_metric_source_term(int ix, int iy, int iz,ldouble *ss);
int f_general_source_term(int ix, int iy, int iz,ldouble *ss);

int f_flux_prime(ldouble *pp, int idim, int ix, int iy, int iz,ldouble *ff,int lr);
int calc_Tij(ldouble *pp, void* ggg, ldouble T[][4]);

ldouble calc_ufromS(ldouble S,ldouble rho,int ix,int iy,int iz);
ldouble calc_TfromS(ldouble S,ldouble rho,int ix,int iy,int iz);
ldouble calc_Sfromu(ldouble rho,ldouble u,int ix,int iy,int iz);
ldouble calc_SfromT(ldouble rho,ldouble u,int ix,int iy,int iz);

ldouble calc_S2fromrhou(ldouble rho, ldouble uint,int type);
ldouble calc_S2fromrhoT(ldouble rho, ldouble temp,int type);
ldouble calc_TfromS2rho(ldouble S2,ldouble rho, int type,int ix,int iy,int iz);
ldouble calc_ufromS2rho(ldouble S2,ldouble rho,int type,int ix,int iy,int iz);

ldouble calc_S3fromrhoT(ldouble rho, ldouble temp,int type);
ldouble calc_S3fromrhou(ldouble rho, ldouble uint,int type);
ldouble calc_TfromS3rho(ldouble S3,ldouble rho, int type,int ix,int iy,int iz);
ldouble calc_ufromS3rho(ldouble S3,ldouble rho,int type,int ix,int iy,int iz);

int update_entropy_cell(int ix,int iy,int iz,int u2pflag);
ldouble entri_from_entre_energy_cons(ldouble *pp, int ix, int iy, int iz);
ldouble calc_PEQ_Teifrompp(ldouble* pp,ldouble* Te, ldouble* Ti,int ix, int iy, int iz);
ldouble calc_PEQ_ugasfromrhoTei(double rho,ldouble Te,ldouble Ti,int ix,int iy,int iz);
ldouble calc_PEQ_ufromTrho(ldouble T,ldouble rho,int ix,int iy,int iz);
ldouble calc_PEQ_Tfromurho(ldouble u,ldouble rho,int ix,int iy,int iz);
ldouble calc_PEQ_Tfromprho(ldouble p,ldouble rho,int ix,int iy,int iz);

int set_gammagas(int type);
ldouble pick_gammagas(int ix,int iy,int iz);
ldouble calc_gammagas(ldouble* pp,int ix,int iy,int iz);
ldouble calc_gammaintfromTei(ldouble Te,ldouble Ti);
ldouble calc_gammaintfromtheta(ldouble theta);
ldouble calc_meanlorentz(ldouble theta);
ldouble calc_gammaintfromtemp(ldouble temp,int type);
ldouble solve_Teifromnmu(ldouble n, ldouble m, ldouble u, int species);
ldouble solve_Teifromnmu_inconsistent(ldouble n, ldouble m, ldouble u);
ldouble solve_Teifromnmu_consistent(ldouble n, ldouble m, ldouble u);
ldouble uint_function(ldouble n, ldouble m, ldouble u,ldouble theta);
ldouble duint_dtheta(ldouble n, ldouble m, ldouble theta);

ldouble heat_electronions_with_state(ldouble dtin);
ldouble calc_ViscousElectronHeatingFraction(ldouble *pp,void *ggg);
ldouble calc_ViscousElectronHeatingFraction_from_state(ldouble *pp,void *sss,void* ggg);
ldouble apply_du_dn_2_species(ldouble *pp, ldouble u, ldouble n, ldouble du, ldouble dn,
			      void* ggg, int type, ldouble* Sreturn);
ldouble pick_ViscousHeating(int ix,int iy,int iz);

int test_gammagas();
int test_calcgamma();

///////////////////////////////////////////////////////////////
// nonthermal.c ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int set_relel_gammas();
ldouble integrate_relel_q(ldouble *q);
ldouble calc_relel_ne(ldouble *pp);
ldouble calc_relel_uint(ldouble *pp);
ldouble calc_relel_p(ldouble *pp);
ldouble calc_gammainj_max_syncool(ldouble bsq, ldouble dtau);
ldouble calc_gammainj_min_jointhermal(ldouble theta, ldouble delta_nth, ldouble p_index, ldouble gammamax);
int reconnection_plaw_params_from_state(ldouble *pp, void *ggg, void *sss, ldouble* delta_back, ldouble* pindex_back);
//ldouble calc_S4fromnT(ldouble n, ldouble temp, int type);
//ldouble calc_S4fromnu(ldouble n, ldouble uint,int type);
//ldouble calc_TfromS4n(ldouble S4,ldouble n, int type,int ix,int iy,int iz);
//ldouble calc_ufromS4n(ldouble S4,ldouble n,int type,int ix,int iy,int iz);
ldouble calc_gammaint_relel(ldouble* pp, ldouble Te, ldouble Ti);
ldouble calc_PEQ_ugasfrom_Tei_relel(ldouble *pp, ldouble Te,ldouble Ti);
ldouble chemical_potential_short(ldouble theta, ldouble neth);
ldouble chemical_potential_long(ldouble p, ldouble u, ldouble sn, ldouble neth, ldouble Te);
int apply_relel_visc_heating(ldouble pp[NV], ldouble durelel, ldouble p_index, ldouble injmin, ldouble injmax, ldouble dtau);
ldouble relel_adiab(ldouble dt);
ldouble calc_gammainj(ldouble theta, ldouble delta_nth, ldouble p_index);
int relel_adiab_derivs(ldouble theta, ldouble dtau, ldouble *ngammas, ldouble *q);
int relel_adiab_derivs_logspace(ldouble theta, ldouble dtau, ldouble *x, ldouble *q);
int relel_adiab_derivs_logspace_lf(ldouble theta, ldouble dtau, ldouble *x, ldouble *q);
ldouble remove_lowgamma_electrons(ldouble thetae, ldouble neth, ldouble *pp);
ldouble gdot_syn(ldouble gamma, ldouble bsq_cgs);
ldouble gdot_ff(ldouble gamma, ldouble nion_cgs);
ldouble gdot_ic(ldouble gamma, ldouble trad_cgs, ldouble eradhat_cgs);
ldouble gdot_coul(ldouble gamma, ldouble neth_cgs);
ldouble gdot_turb_advect(ldouble gamma);
ldouble d_turb_diffuse(ldouble gamma);

int calc_relel_f_and_fmag_from_state(ldouble *pp, void *sss, ldouble *pp0, void *ggg, ldouble dtau,ldouble *frel,ldouble *frelmag);
int calc_relel_cooling_from_state(ldouble *pp, void *sss, ldouble *pp0, ldouble dtau, ldouble *qcool);
int calc_relel_cooling_lf_from_state(ldouble *pp, void *sss, ldouble *pp0, ldouble dtau, ldouble *qcool);
ldouble calc_relel_G0_fluidframe_from_state(ldouble *pp, void *sss, void *ggg, ldouble relel_dudtau, int type);
ldouble calc_relel_G0_fluidframe_direct_from_state(ldouble *pp, void *sss, int radtype);
ldouble calc_relel_photon_ndot_from_state(ldouble *pp,  void* sss, int radtype);
ldouble calc_relel_cool_dq_from_state(ldouble *pp, void *sss);
ldouble calc_relel_cool_dn_from_state(ldouble *pp, void *sss);
ldouble calc_relel_CoulombCoupling_from_state(ldouble *pp, void *sss);

//ANDREW Needed still for silo.c ... 
ldouble calc_relel_G0_fluidframe(ldouble *pp, void *ggg, ldouble relel_dudtau, int type);
ldouble calc_relel_G0_fluidframe_direct(ldouble *pp, void *ggg, int radtype);

int fprint_relel_spectrum(ldouble t, int ix, int iy, int iz, int nfile, char* folder, char* prefix, int doingavg);
int fprint_relel_avg_spectrum(ldouble t, int jx, int jy, int jz, int nfile, char* folder, char* prefix,int doingavg);

///////////////////////////////////////////////////////////////
// postproc.c /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

int calc_radialprofiles(ldouble profiles[][NX]);
int calc_thetaprofiles(ldouble profiles[][NY]);
int calc_scalars(ldouble *scalars,ldouble t);

ldouble calc_totalmass();
ldouble calc_mdotEdd();
ldouble calc_lumEdd();
int calc_local_lum(int ix,int iy,int iz,ldouble *radlum, ldouble *totallum);
int calc_lum(ldouble radius,int type,ldouble *radlum, ldouble *totallum);
ldouble calc_exitlum();
ldouble calc_resmri(ldouble radius);
ldouble calc_meantemp(ldouble radius);
ldouble calc_scaleheight(ldouble radius);
ldouble calc_photloc(int ix);
ldouble calc_mdot(ldouble radius,int type);
ldouble calc_Edot(ldouble radius);
ldouble calc_lum_proxy(ldouble radius, ldouble theta_min, ldouble theta_max);
ldouble calc_Ldot(ldouble radius);
int calc_Bflux(ldouble radius,int type,ldouble *Bflux, ldouble* Bfluxquad);

int calc_lum_tausurface(ldouble taumax,ldouble *radlum);

///////////////////////////////////////////////////////////////
// fileop.c ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//deleted OUTOUTPUT BOXOUTPUT BOXCORROUTPUT VAROUTPUT, BOXVERTOUTPUT ANARELRADOUTPUT SLICEOUTPUT

int save_avg(ldouble dtin);
int fprint_openfiles(char* folder);
int fprint_closefiles();
int fprint_gridfile(char* folder);
int fprint_scalars(ldouble t, ldouble *scalars, int nscalars);
int fprint_radprofiles(ldouble t, int nfile, char* folder, char* prefix);
int fprint_thprofiles(ldouble t, int nfile, char* folder, char* prefix);
int fprint_restartfile(ldouble t, char* folder);
int fprint_restartfile_mpi(ldouble t, char* folder);
int fprint_restartfile_bin(ldouble t, char* folder);
int fread_restartfile(int nout1, char* folder,ldouble *t);
int fread_restartfile_bin(int nout1, char *folder, ldouble *t);
int fread_restartfile_mpi(int nout1, char *folder, ldouble *t);
int fprint_avgfile(ldouble t, char* folder,char* prefix);
int fprint_avgfile_mpi(ldouble t, char* folder, char* prefix);
int fprint_avgfile_bin(ldouble t, char* folder,char *prefix);
int fread_avgfile(int nout1, char* base,ldouble *pavg, ldouble *dt,ldouble *t);
int fread_avgfile_bin(int nout1, char *base,ldouble *pavg, ldouble *dt,ldouble *t);
int fprint_coordfile(char* folder,char* prefix);
int fprint_coordBL(char* folder,char* prefix);
int fprint_simplefile(ldouble t, int nfile, char* folder,char* prefix);
int fprint_simplecart(ldouble t, int nfile, char* folder,char* prefix);
int fprint_simplesph(ldouble t, int nfile, char* folder,char* prefix);
int fprint_simple_phiavg(ldouble t, int nfile, char* folder,char* prefix);
int fprint_simple_phicorr(ldouble t, int nfile, char* folder,char* prefix);

int fprint_restartfile_mpi_hdf5(ldouble t, char* folder);
int fprint_restartfile_serial_hdf5(ldouble t, char* folder);
int fread_restartfile_mpi_hdf5(int nout1, char *folder, ldouble *t);
int fread_restartfile_serial_hdf5(int nout1, char *folder, ldouble *t);
int fprint_anaout_hdf5(ldouble t, char* folder,char* prefix);

void get_prim_name(char* prim_name, int iv);


///////////////////////////////////////////////////////////////
// silo.c /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#ifndef NOSILO
int fprint_silofile(ldouble time, int num, char* folder, char* prefix);
#endif

///////////////////////////////////////////////////////////////
// opacities.c ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

ldouble calc_kappa(ldouble *pp,void *ggg,void *op);
ldouble calc_kappa_from_state(ldouble *pp, void *sss, void *ggg, void *op);
ldouble calc_kappaes_with_temperatures(ldouble rho, ldouble Tgas, ldouble Te, ldouble Ti, ldouble Trad);
ldouble calc_kappaes(ldouble *pp,void *ggg);
ldouble calc_chi(ldouble *pp,void *ggg);
int calc_tautot(ldouble *pp, void *ggg, ldouble *dx, ldouble *tautot);
int calc_tauabs(ldouble *pp, void *ggg, ldouble *dx, ldouble *tauabs);
ldouble calc_opacities_from_state(ldouble *pp, void *ggg, void *sss, void *op);
int init_OpTable(void *optab0, char *filename);
int init_all_kappa_table();

#ifdef USE_CHIANTI_ISM_TABLE
int init_ChiantiISMTable(void);
int ChiantiISMTableLength;
ldouble **ChiantiISMTable;
ldouble return_ChiantiISMTableOpacity(ldouble T, ldouble rho);
#endif

#ifdef USE_PLANCK_TABLE
struct OpTable PlanckTable;
ldouble *chiantilogkappa;
ldouble *chiantilogT;
#define NCHIANTI 81
ldouble return_Chianti(ldouble logTin, ldouble logrhoin);
ldouble return_PlanckOpacity_from_table(ldouble logTin, ldouble logrhoin);
#endif

#ifdef USE_ROSS_TABLE
struct OpTable RossTable;
ldouble return_RossOpacity_from_table(ldouble logTin, ldouble logrhoin);
#endif
