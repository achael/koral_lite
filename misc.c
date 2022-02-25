/*! \file misc.c
 \brief misceleanous routines
*/

#include "ko.h"

//**********************************************************************
/*! \fn int initialize_consants()
 \brief Set constants to be used throughout
*/
//**********************************************************************
void initialize_constants()
{
  // Conversion Factors  
  tempcgs2gu = 1.;
  tempgu2cgs = 1.;
  lencgs2gu = 1. / MASSCM;
  lengu2cgs = MASSCM;
  numdenscgs2gu = MASSCM*MASSCM*MASSCM;
  numdensgu2cgs = (1./MASSCM/MASSCM/MASSCM);
  timecgs2gu = 1. /MASSCM*CCC;
  timegu2cgs = 1. / timecgs2gu;
  velcgs2gu = 1. / CCC;
  velgu2cgs = CCC;
  rhocgs2gu = GGG/CCC/CCC*MASSCM*MASSCM;
  rhogu2cgs = 1./GGG*CCC*CCC/MASSCM/MASSCM;
  surfdenscgs2gu = GGG/CCC/CCC*MASSCM;
  surfdensgu2cgs = 1. / surfdenscgs2gu;
  masscgs2gu = GGG/CCC/CCC/MASSCM;
  massgu2cgs = 1. / masscgs2gu;
  kappacgs2gu = 1. /GGG*CCC*CCC/MASSCM;
  kappagu2cgs = 1. / kappacgs2gu;
  endencgs2gu = GGG*MASSCM*MASSCM/CCC/CCC/CCC/CCC;
  endengu2cgs = 1. / endencgs2gu;
  heatcoolcgs2gu =  endencgs2gu * timegu2cgs;
  heatcoolgu2cgs = 1. / heatcoolcgs2gu;
  fluxcgs2gu = GGG*MASSCM*MASSCM/CCC/CCC/CCC/CCC/CCC;
  fluxgu2cgs = 1. / fluxcgs2gu;
  ergcgs2gu = GGG/MASSCM/CCC/CCC/CCC/CCC;
  erggu2cgs = 1. / ergcgs2gu;
  chargecgs2gu = sqrt(GGG/MASSCM/CCC/CCC/CCC/CCC)*sqrt(1./MASSCM);
  chargegu2cgs = 1. / chargecgs2gu;
  crosscgs2gu = 1. / MASSCM/MASSCM;
  crossgu2cgs = 1. / crosscgs2gu;
  
  // Physical constants  
  k_boltz_cgs = K_BOLTZ_CGS * 1.;
  k_boltz_cgs_inv = 1. / k_boltz_cgs;
  k_boltz = K_BOLTZ;
  k_boltz_inv = (1. / K_BOLTZ);
  m_proton_cgs = M_PROTON_CGS;
  m_proton = M_PROTON;
  m_electr_cgs = M_ELECTR_CGS;
  m_electr = M_ELECTR;
  mpe_ratio = MPE_RATIO;
  sigma_rad_cgs = SIGMA_RAD_CGS;
  sigma_rad = SIGMA_RAD;
  h_cgs = H_CGS;
  sigmath_cgs = SIGMATH_CGS;
  a_rad = A_RAD;
  z_ratio = Z_RATIO;
  e_charge = E_CHARGE;
  sigma_thomson = SIGMA_THOMPSON;

  // Miscellaneous coefficients  
  fourpi = 4. * Pi;
  fourmpi = 4. * M_PI;
  per_volume_per_time_cgs2gu = 1. / (lencgs2gu * lencgs2gu * lencgs2gu * timecgs2gu);
  sigma_rad_over_pi = SIGMA_RAD / Pi;
  four_sigmarad = 4. * sigma_rad;
  one_over_four_sigmarad = 1. / four_sigmarad;
  k_over_mecsq = 1.69e-10;
  kB_over_mui_mp = (K_BOLTZ / MU_I / M_PROTON);
  kB_over_mue_mp = (K_BOLTZ / MU_E / M_PROTON);
  kB_over_mugas_mp = (K_BOLTZ / MU_GAS / M_PROTON);
  mugas_mp_over_kB = 1. / kB_over_mugas_mp;
  kB_over_mp = (K_BOLTZ / M_PROTON);
  kB_over_me = (K_BOLTZ / M_ELECTR);
  one_over_kB_2p70118 = k_boltz_inv / 2.70118;
  one_over_mugas_mp = (1. / MU_GAS / M_PROTON);
  one_over_mui_mp = (1. / MU_I / M_PROTON);
  one_over_mue_mp =  (1. / MU_E / M_PROTON);
  mui_over_mue = (MU_I/MU_E);
    
  four_third = 4. / 3.;
  one_third = 1. / 3.;
  two_third = 2. / 3.;
  log_2p6 = log(2.6);
  one_over_log_2p6 = 1. / log_2p6;

  /*  
  printf("Testing gammainterp\n");
  printf("%e %e %e \n",calc_meanlorentz(1.e-3),calc_meanlorentz(1.1e-3),calc_meanlorentz(5.e-3));
  printf("%e %e %e \n",calc_meanlorentz(0.74),calc_meanlorentz(1.23),calc_meanlorentz(7.54));
  printf("%e %e %e \n",calc_meanlorentz(14.3),calc_meanlorentz(77.3),calc_meanlorentz(144.));
  printf("%e %e %e \n",calc_meanlorentz(555.),calc_meanlorentz(999.),calc_meanlorentz(1000.));
  exit(-1);
  */
  // Coordinate specific factors
  #if (MYCOORDS==JETCOORDS)
  //printf("Finding hypx1out\n");
  hypx1in = log(RMIN-MKSR0);
  hypx1brk= log(HYPRBRK-MKSR0);
  hypx1out= hyperexp_x1max(RMAX, HYPRBRK, MKSR0);

  //printf("hyperx1in %e | hyperx2brk %e | hyperx1out %e\n",hypx1in,hypx1brk,hypx1out);
  #ifdef CYLINDRIFY
  set_cyl_params();
  #endif

  
  //ANDREW -- diagnostics for new jet coordinates/metric

  /*
  //printf("%.7f %.7f %.7f\n",hypx1in,hypx1brk,hypx1out);

  ldouble x0[4] = {0, .15, 0.14863, 1};
  ldouble x1[4], x2[4];

  printf("%.7f %.7f %.7f %.7f\n",x0[0],x0[1],x0[2],x0[3]);

  struct timespec temp_clock;
  my_clock_gettime(&temp_clock);    
  start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
  
  coco_JET2KS(x0,x1);
  
  my_clock_gettime(&temp_clock);    
  end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
  printf("Transform Time %.7f\n",end_time-start_time);
  
  printf("%.7f %.7f %.7f %.7f\n",x1[0],x1[1],x1[2],x1[3]);  

  coco_KS2JET(x1,x2);
  printf("%.7f %.7f %.7f %.7f\n",x2[0],x2[1],x2[2],x2[3]);

  
  ldouble dxdx[4][4], dxdxinv[4][4], dxdxinv2[4][4];
  int i,j,tmp;

  my_clock_gettime(&temp_clock);    
  start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

  tmp = dxdx_arb_num(x0, dxdx, MYCOORDS, KSCOORDS);

  my_clock_gettime(&temp_clock);    
  end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
  printf("DXDX Time %.7f\n",end_time-start_time);

  tmp = dxdx_arb_num(x1, dxdxinv, KSCOORDS, MYCOORDS);
  inverse_44matrix(dxdx,dxdxinv2);    

  print_tensor(dxdx);
  print_tensor(dxdxinv);
  //print_tensor(dxdxinv2);
  printf("\n\n====================\n\n");
  ldouble gtmp[4][5], g[4][4], G[4][4], ginv[4][4], Ginv[4][4];

  my_clock_gettime(&temp_clock);    
  start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

  calc_g_arb_num(x0, gtmp, MYCOORDS);
  my_clock_gettime(&temp_clock);    
  end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
  printf("Metric Time %.7f\n",end_time-start_time);
  
  DLOOP(i,j) g[i][j]=gtmp[i][j];
  calc_G_arb_num(x0, gtmp, MYCOORDS);
  DLOOP(i,j) G[i][j]=gtmp[i][j];
  inverse_44matrix(g,ginv);
  inverse_44matrix(G,Ginv);

  //DLOOP(i,j) {if(fabs(g[i][j])<1.e-9) g[i][j]=0;}
  //DLOOP(i,j) {if(fabs(G[i][j])<1.e-9) G[i][j]=0;}
  //DLOOP(i,j) {if(fabs(ginv[i][j])<1.e-9) ginv[i][j]=0;}
  //DLOOP(i,j) {if(fabs(Ginv[i][j])<1.e-9) Ginv[i][j]=0;}
  
  print_tensor(g);
  //print_tensor(Ginv);
  //DLOOP(i,j) Ginv[i][j] = (g[i][j]-Ginv[i][j])/(g[i][j]+SMALL);
  //print_tensor(Ginv);
  
  printf("\n\n====================\n\n");
  print_tensor(G);
  //print_tensor(ginv);
  //DLOOP(i,j) ginv[i][j] = (G[i][j]-ginv[i][j])/(G[i][j]+SMALL);
  //print_tensor(ginv);

 
  exit(-1);
  */
  #endif
  
  return;
}


//**********************************************************************
/*! \fn int print_scalings()
 \brief Print out key parameters and scalings of the simulation
*/
//**********************************************************************
int
print_scalings()
{
  printf("\n ***************************************\n\n");
  printf("BH mass: %.6f\nspin: %.6f\n\nscalings  (GU->CGS):\nrho: %.6e\nmdot: %.6e\nsigma: %.6e\nlen: %.6e\ntime: %.6e\nenden:"
         "%.6e\nflux: %.6e\nT(1,1): %.6e\nkbt: %.6e\nkb/me: %.6e\nsigma_rad: %.6e\nkappa: %.6e\nGt: %.6e\nmass: %.6e\n\n"
         "rhorizonBL: %.6f\nrISCOBL: %.6f\netaNT: %.6f\n\n->mdotEdd: %.6e\n->lumEdd: %.6e\n\nGMc2: %.6e\nGMc3: %.6e\n"
         "Mdot_Edd: %.6e g/s\nLdot_Edd: %.6e erg/s\n",
         MASS,BHSPIN,
         rhoGU2CGS(1.),
         rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.),
         rhoGU2CGS(1.)*lenGU2CGS(1.),
         lenGU2CGS(1.),
         timeGU2CGS(1.),
         endenGU2CGS(1.),
         fluxGU2CGS(1.),
         calc_PEQ_Tfromurho(1.,1.,0,0,0),
         K_BOLTZ/MU_GAS/M_PROTON,
         K_BOLTZ/M_ELECTR,
         SIGMA_RAD,
         kappaGU2CGS(1.),
         kappaGU2CGS(1.)*rhoGU2CGS(1.)*endenGU2CGS(1.)*CCC,
         massGU2CGS(1.),
         rhorizonBL,
         rISCOBL,
         etaNT,
         rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)/calc_mdotEdd(),
         rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)*CCC0*CCC0/calc_lumEdd(),
         GMC2,
         GMC3,
         calc_mdotEdd(),
         calc_lumEdd()
         );
  printf("\n ***************************************\n\n");
  return 0;
}

//**********************************************************************
/*! \fn int set_initial_profile()
 \brief Sets the initial distributions of quantities
 
 Details are to be supplied in the file: PROBLEMS/XXX/init.c
*/
//**********************************************************************
int
set_initial_profile()
{
  if(PROCID==0) {printf("Initializing problem... \n");fflush(stdout);}
  int ix,iy,iz;
  
#pragma omp parallel for private(ix,iy,iz) schedule (static)
  for(ix=0;ix<NX;ix++)
  {
    for(iy=0;iy<NY;iy++)
    {
      for(iz=0;iz<NZ;iz++)
      {
        
#include PR_INIT
        
      }
    }
  }

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
#ifdef OMP
#pragma omp barrier
#endif
  
  if(PROCID==0) printf("done!\n");
  
  return 0;
} 

//**********************************************************************
/*! \fn void am_i_sane()
 \brief Checks that there are no conflicts in the chosen settings
 */
//**********************************************************************
void
am_i_sane()
{

  if(!(MYCOORDS==SPHCOORDS || MYCOORDS==CYLCOORDS || MYCOORDS==MINKCOORDS ||
     MYCOORDS==KERRCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==MKS1COORDS ||
     MYCOORDS==MKS2COORDS || MYCOORDS==MCYL1COORDS || MYCOORDS==MSPH1COORDS ||
     MYCOORDS==MKER1COORDS))
  {
    if(GDETIN==0)
    {
      printf("GDETIN==1 does not work with this coordinate system!\n");
      exit(-1);
    }
     #ifndef METRICNUMERIC
    //printf("this coordinate system requires METRICNUMERIC!\n");
    //exit(-1);
     #endif
  }

  if(HFRAC+HEFRAC+MFRAC != 1.)
  {
    printf("Chosen H, He, and Metal abundance does not sum to 1.\n");
    exit(-1);
  }

  if (INT_ORDER > 2)
  {
    printf("Only INT_ORDER=1 and INT_ORDER=2 supported!\n");
    exit(-1);
  }

#ifdef CORRECT_POLARAXIS
  if (TNY == 1)
  {
    printf("\nERROR!! Using CORRECT_POLARAXIS on a 2D problem in r-phi. No evolution will occur! Switch off polar axis correction and recompile!\n\n");
    exit(-1);
  }
#endif

#ifdef TRANSMITTING_YBC
  if(PROCID==0) 
  {
    printf("Note: Using TRANSMITTING_YBC. Recommended MINY >= 0.005. Difussion around pole dominates transmission otherwise.\n");
  }


if (TNZ % 2 != 0)
{
  printf("\nERROR!! Using transmitting y boundary. TNZ = %d must be divisible by 2!\n\n", TNZ);
  exit(-1);
}
if (NTZ % 2 != 0)
{
  printf("\nERROR!! Using transmitting y boundary. NTZ = %d must be divisible by 2!\n\n", NTZ);
  exit(-1);
}
#endif  

#ifdef PRECOMPUTE_MY2OUT
  if (PROCID==0){
  printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  printf("\nWarning -- using precomputed MYCOORDS --> OUTCOORDS for floors, averages, and boundary conditions\n");
  printf("Check your bcs.c!!\n");
#ifdef BHDISK_PROBLEMTYPE
  if(!((OUTCOORDS == BLCOORDS) || (OUTCOORDS == KSCOORDS))) {
    printf("For BHDISK_PROBLEMTYPE PRECOMPUTE_MY2OUT currently only works with OUTCOORDS=BLCOORDS!\n");
    printf("(It should also work with KSCOORDS, but just being safe for now...)\n");
      exit(-1);
  }
#endif
 printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
}
#endif
  
#ifdef PWPOTENTIAL
  printf("PWPOTENTIAL has been removed!\n");
  exit(-1);
#endif

#ifdef NCOMPTONIZATION
  printf("NCOMPTONIZATION has been replaced by EVOLVEPHOTONNUMBER!\n");
  exit(-1);
#endif

#ifdef RADIATION  
#if defined(COMPTONIZATIONFLAG)
  if (PROCID == 0)
  {
    printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\nAutomatically switching on COMPTONIZATION\n");
    printf("Please define NO_COMPTONIZATION if Comptonization should remain off\n");
#ifdef EVOLVEPHOTONNUMBER
    printf("EVOLVEPHOTONNUMBER is ON\n");
#else
    printf("EVOLVEPHOTONNUMBER is OFF\n");
#endif
    printf("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
  }
#endif
#endif

#ifdef COLDOPACITIES
  printf("COLDOPACITIES no longer defined!\n");
  exit(-1);
#endif
	 
#ifdef EVOLVEELECTRONS
#ifdef FIXEDTETIRATIO
  printf("FIXEDTETIRATIO can't be used together with EVOLVEELECTRONS\n");
  exit(-1);
#endif

#if defined(RADIATION) && !defined(SKIPRADSOURCE)
#ifdef HEATELECTRONSATENDRK2
  printf("HEATELECTRONSATENDRK2 only works without RADIATION!");
  exit(-1);
#endif
#endif  
#endif
  
#ifdef SHEARINGBOX
#ifdef MPI
  printf("SHEARINGBOX BC not implemented into MPI\n");
  exit(-1);
#endif
  
#ifndef NONRELMHD
  printf("SHEARINGBOX requires NONRELMHD.\n");
  exit(-1);
#endif
  
#if (MYCOORDS!=MINKCOORDS)
  printf("SHEARINGBOX requires  MINKCOORDS.\n");
  exit(-1);
#endif
#endif  // ifdef SHEARINGBOX
  
#ifdef NONRELMHD  
  if(MYCOORDS!=SPHCOORDS && MYCOORDS!=CYLCOORDS && MYCOORDS!=MINKCOORDS)
  {
    printf("NONRELMHD implemented only for flat coords so far.\n");
    exit(-1);
  }

#ifdef EVOLVEPHOTONNUMBER
  printf("NONRELMHD not implemented for EVOLVEPHOTONNUMBER.\n");
  exit(-1);
#endif
#endif  // ifdef NONRELMHD
  
#ifdef CONSISTENTGAMMA
#ifndef EVOLVEELECTRONS
  printf("CONSISTENTGAMMA works only with EVOLVEELECTRONS\n");
  exit(-1);
#endif
#endif

#ifdef RELELECTRONS
#ifndef EVOLVEELECTRONS
  printf("RELELECTRONS requires EVOLVEELECTRONS\n");
  exit(-1);
#endif
#ifndef NRELBIN
  printf("RELELECTRONS requires NRELBIN\n");
  exit(-1);
#endif
#ifndef RELGAMMAMIN
  printf("RELELECTRONS requires RELGAMMAMIN\n");
  exit(-1);
#endif
#ifndef RELGAMMAMAX
  printf("RELELECTRONS requires RELGAMMAMAX\n");
  exit(-1);
#endif
#endif
  
#ifdef PWPOTENTIAL
  if(MYCOORDS!=SPHCOORDS && MYCOORDS!=CYLCOORDS)
  {
    printf("PWPOTENTIAL implemented only for SPHCOORDS or CYLCOORDS so far.\n");
    exit(-1);
  }
#endif
  
  if (FLUXMETHOD == HLLC_FLUX)
  {
    printf("\n HLLC_FLUX is not supported anymore!\n");
  }
  
#ifdef SUBZONES 
  printf("\nSUBZONES is not supported anymore!!\n\n");
  exit(0);
#endif

#ifdef EVOLVEINTENSITIES 
  printf("\nEVOLVEINTENSITIES not supported anymore!!\n\n");
  exit(0);
#endif

#if (RADCLOSURE==VETCLOSURE) 
  printf("\nVETCLOSURE is not supported anymore!!\n\n");
  exit(0);
#endif
  
#ifdef METRICTIMEDEPENDENT 
  printf("\nMETRICTIMEDEPENDENT is not currently supported!!\n\n");
  exit(0);
#endif

#ifdef NUMRADWAVESPEEDS
  printf("NUMRADWAVESPEEDS is not currently supported!\n");
  exit(-1);
#endif

#ifdef RADOUTPUTINFF
  printf("RADOUTPUTINFF no longer supported!\n");
  exit(-1);
#endif

#ifdef RADOUTPUTINZAMO
  printf("RADOUTPUTINZAMO no longer supported!\n");
  exit(-1);
#endif

#ifdef RESTARTFROMMHD
  if(PROCID==0) printf("RESTARTFROMMHD\n");
  if(PROCID==0) printf("urad/uu: %e ue/uu: %e\n", INITURADFRAC, INITUEFRAC);
#endif
  return;
}


//**********************************************************************
// allocate arrays
//**********************************************************************

int
initialize_arrays()
{
  long long i,j,k;

  /********************* Basic arrays ***************************/

  long long Ngrid=(SX)*(SY)*(SZ);
  long long GridSize=Ngrid*sizeof(ldouble);
  long long Nprim=Ngrid*NV;
  long long PrimSize=Nprim*sizeof(ldouble);
  long long Navg=Ngrid*(NV+NAVGVARS);
  long long AvgSize=Navg*sizeof(ldouble);


  long long Ngridmet=(SX)*(SY)*(SZMET);
  long long Nmet=Ngridmet*gSIZE;
  long long MetSize=Nmet*sizeof(ldouble);
  long long Nkris=(SX)*(SY)*(SZMET)*64;
  long long KrisSize=Nkris*sizeof(ldouble);
  
  long long Ntensor=Ngrid*16;
  long long TensorSize=Ntensor*sizeof(ldouble);
  
  //grid
  if((x=(ldouble*)malloc((NX+NY+NZ+6*NG)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  if((xb=(ldouble*)malloc((NX+1+NY+1+NZ+1+6*NG)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");

#ifdef PRECOMPUTE_MY2OUT
  //arrays for MYCOORDS->OUTCOORDS transformation
  if((xout=(ldouble*)malloc((Ngrid*3)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  if((dxdx_my2out=(ldouble*)malloc((Ngridmet*16)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  if((dxdx_out2my=(ldouble*)malloc((Ngridmet*16)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");

  //outcoords at x-faces // only if postproc==1 ??
  long long Ngrid_xface=(SX+1)*(SY)*(SZ);
  long long Ngrid_yface=(SX)*(SY+1)*(SZ);
  long long Ngrid_zface=(SX)*(SY)*(SZ+1);
  if((xbout_xface=(ldouble*)malloc((Ngrid_xface*3)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  if((xbout_yface=(ldouble*)malloc((Ngrid_yface*3)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  if((xbout_zface=(ldouble*)malloc((Ngrid_zface*3)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
#endif
  
  //primitives at cell centers
  if((p=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

  //quantities to average in time
  if((pavg=(ldouble*)malloc(AvgSize))==NULL) my_err("malloc err.\n");
  
  //conserved averages
  if((u=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

  //flags at cell centers
  if((cellflag=(int*)malloc(Ngrid*NFLAGS*sizeof(int)))==NULL) my_err("malloc err.\n");
 
  //metric at cell centers
  if((g=(ldouble*)malloc(MetSize))==NULL) my_err("malloc err.\n");
  if((G=(ldouble*)malloc(MetSize))==NULL) my_err("malloc err.\n");

  //Kristofels at cell centers
  if((gKr=(ldouble*)malloc(KrisSize))==NULL) my_err("malloc err.\n");

  //primitives at cell centers at initial state - used for fixed boundary conditions
  if((pinit=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
 
  //primitives at cell centers at the beginning and end of explicit operator
  if((upreexplicit=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
  if((ppreexplicit=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
   
  //primitives at cell centers at the end of implicit operator
  if((ppostimplicit=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

  //arrays for temporary use (e.g., vector potential, mimic_dynamo)
  if((pproblem1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
  if((ptemp1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
  if((pvecpot=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

  //arrays for avg time
  if((avgselftime=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");  

  //arrays for radiation tensor
#ifdef RADIATION
#if (RADVISCOSITY==SHEARVISCOSITY)
  if((Rijviscprev=(ldouble*)malloc(TensorSize))==NULL) my_err("malloc err.\n");
  if((Rijviscglobal=(ldouble*)malloc(TensorSize))==NULL) my_err("malloc err.\n");
  if((radvisclasttime=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
#endif
#endif

  //arrays for viscous heating
  if((vischeating=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
  if((vischeatingnegebalance=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
  if((vischeatingnegibalance=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
  if((vischeatingtimesdeltae=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");

  //these will aggregate over time so should start with zeros
  for(i=0;i<Ngrid;i++)
    vischeatingnegebalance[i]=vischeatingnegibalance[i]=0.;

  //gamma of gas at the beginnning of timestep
  if((gammagas=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
 
  /****************** extra arrays, used only for time evolution **********************/
  //we might need some of these arrays in postproc (if doingpostproc_avg==1)
  #ifdef DIVIDEVISCHEATBYDT
  if(1.)
  #else
  if(doingpostproc==0 || doingpostproc_avg==1)
  #endif
  {     
     //buffer for sending/receiving messages
#ifdef MPI
     if((msgbufs=(ldouble**)malloc(MPIMSGBUFSIZE*sizeof(ldouble*)))==NULL) my_err("malloc err.\n");
     for(i=0;i<MPIMSGBUFSIZE;i++)
       if((msgbufs[i]=(ldouble*)malloc(my_max3(NX*NY*NV*NG,NY*NZ*NV*NG,NZ*NX*NV*NG)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
#endif

     long long NMetX = (SX+1)*(SY)*(SZMET)*gSIZE;
     long long MetXSize=NMetX*sizeof(ldouble);
     long long NMetY = (SX)*(SY+1)*(SZMET)*gSIZE;
     long long MetYSize=NMetY*sizeof(ldouble);
     long long NMetZ = (SX)*(SY)*(SZMET+1)*gSIZE;
     long long MetZSize=NMetZ*sizeof(ldouble);
     
     long long NMetVec = (SX)*(SY)*(SZMET)*16;
     long long MetVecSize=NMetVec*sizeof(ldouble);
            
     //metric at cell x-faces
     if((gbx=(ldouble*)malloc(MetXSize))==NULL) my_err("malloc err.\n");
     if((Gbx=(ldouble*)malloc(MetXSize))==NULL) my_err("malloc err.\n");

     //metric at cell y-faces
     if((gby=(ldouble*)malloc(MetYSize))==NULL) my_err("malloc err.\n");
     if((Gby=(ldouble*)malloc(MetYSize))==NULL) my_err("malloc err.\n");

     //metric at cell z-faces
     if((gbz=(ldouble*)malloc(MetZSize))==NULL) my_err("malloc err.\n");
     if((Gbz=(ldouble*)malloc(MetZSize))==NULL) my_err("malloc err.\n");
      
     //LNRF basis one-forms and vectors
     //if((emuup=(ldouble*)malloc(MetVecSize))==NULL) my_err("malloc err.\n");
     //if((emulo=(ldouble*)malloc(MetVecSize))==NULL) my_err("malloc err.\n");

     //tetrad one-forms and vectors
     //if((tmuup=(ldouble*)malloc(MetVecSize))==NULL) my_err("malloc err.\n");
     //if((tmulo=(ldouble*)malloc(MetVecSize))==NULL) my_err("malloc err.\n");

     //Fluxes and wavespeeds
     long long NfluxX = (SX+1)*(SY)*(SZ)*NV;
     long long fluxXSize = NfluxX*sizeof(double);
     long long NfluxY = (SX)*(SY+1)*(SZ)*NV;
     long long fluxYSize = NfluxY*sizeof(double);
     long long NfluxZ = (SX)*(SY)*(SZ+1)*NV;
     long long fluxZSize = NfluxZ*sizeof(double);

     //wavespeeds hd and rad - max(al,ar)
     if((ahdx=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdy=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdz=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradx=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((arady=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradz=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     
     //wavespeeds hd and rad - leftgoing
     if((ahdxl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdyl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdzl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradxl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradyl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradzl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
  
     //wavespeeds hd and rad - rightgoing
     if((ahdxr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdyr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdzr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradxr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradyr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradzr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
       
     //left-interpolated primitives at cell faces
     if((pbLx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((pbLy=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((pbLz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //right-interpolated primitives at cell faces
     if((pbRx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((pbRy=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((pbRz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //corrected flux at faces
     if((flbx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((flby=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((flbz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //flux based on left-interpolated conserved at cell faces
     if((flLx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((flLy=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((flLz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //flux based on right-interpolated conserved at cell faces
     if((flRx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((flRy=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((flRz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //auxiliary primitive arrays
     if((pproblem2=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ptm1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ut0=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ut1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ut2=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ut3=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((dut0=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((dut1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((dut2=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((drt0=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((drt1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((drt2=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((uforget=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((u_bak_fixup=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((p_bak_fixup=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

     //timesteps required by each cell
     if((cell_tstepden=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((cell_dt=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");

#ifdef MAGNFIELD
     //electromotive force at corners
     long long Nemf = (NX+1)*(NY+1)*(NZ+1)*3;
     long long EmfSize = Nemf*sizeof(ldouble);
     if((emf=(ldouble*)malloc(EmfSize))==NULL) my_err("malloc err.\n");
#endif
  
  }

  init_all_kappa_table();

  return 0;
}

//**********************************************************************
//Free arrays at end
//**********************************************************************
int
free_arrays()
{
  free(cellflag);
  free(x);
  free(xb);
  free(p);
  free(pavg);
  free(pinit);
  free(upreexplicit);
  free(ppreexplicit);
  free(ppostimplicit);
  free(pproblem1);
  free(pproblem2);
  free(avgselftime);
  free(vischeating);
  free(vischeatingnegebalance);
  free(vischeatingnegibalance);
  free(vischeatingtimesdeltae);
  free(ptemp1);
  free(pvecpot);
#ifdef RADIATION
#if (RADVISCOSITY==SHEARVISCOSITY)
  free(Rijviscprev);
  free(Rijviscglobal);
  free(radvisclasttime);
#endif
#endif  
#ifdef MAGNFIELD
  free(emf);
#endif

  free(ptm1);
#ifdef MPI
  int i;
  for(i=0;i<MPIMSGBUFSIZE;i++)
    free(msgbufs[i]);
  free(msgbufs);
#endif
  
  free(px);
  free(py);
  free(pz);
  free(u);
  free(g);
  free(G);
  free(gKr);
  //free(emuup);
  //free(emulo);
  //free(emuupbx);
  //free(emulobx);
  //free(emuupby);
  //free(emuloby);
  //free(emuupbz);
  //free(emulobz);
  //free(tmuup);
  //free(tmulo);
  //free(tmuupbx);
  //free(tmulobx);
  //free(tmuupby);
  //free(tmuloby);
  //free(tmuupbz);
  //free(tmulobz);

  free(pbLx);
  free(pbRx);
  free(pbLy);
  free(pbRy);
  free(pbLz);
  free(pbRz);
 
  free(flbx);
  free(flby);
  free(flbz);

  free(flLx);
  free(flRx);
  free(flLy);
  free(flRy);
  free(flLz);
  free(flRz);

  free(gbx);
  free(gby);
  free(gbz);
  free(Gbx);
  free(Gby);
  free(Gbz);
 
  free(ut0);
  free(ut1);
  free(ut2);
  free(ut3);
  free(dut0);
  free(dut1);
  free(dut2);
  free(drt0);
  free(drt1);
  free(drt2);
  free(uforget);
  free(u_bak_fixup);
  free(p_bak_fixup);

  free(aradx);
  free(arady);
  free(aradz);
  free(ahdx);
  free(ahdy);
  free(ahdz);
  free(aradxl);
  free(aradyl);
  free(aradzl);
  free(ahdxl);
  free(ahdyl);
  free(ahdzl);
  free(aradxr);
  free(aradyr);
  free(aradzr);
  free(ahdxr);
  free(ahdyr);
  free(ahdzr);

  free(cell_dt);
  free(cell_tstepden);
  free(gammagas);

  #ifdef USE_PLANCK_TABLE
  free(chiantilogkappa);
  free(chiantilogT);
  #endif
  #ifdef USE_CHIANTI_ISM_TABLE
  for(i=0,i<=ChiantiISMTableLength,i++)
    free(ChiantiISMTable[idx]);
  free(ChiantiISMTable);
  #endif
  //#ifdef SUTHERLAND_DOPITA_LAMBDA
  //free(temperaturelog);
  //free(Lambdalog);
  //#endif
  
  return 0;
}

//**********************************************************************
// Initialize pointers to entropy functions 
//**********************************************************************
void
init_pointers()
{

#if defined(CONSISTENTGAMMA) && !defined(FIXEDGAMMASPECIES)
  calc_SefromrhoT=&calc_S3fromrhoT;
  calc_Sefromrhou=&calc_S3fromrhou;
  calc_TfromSerho=&calc_TfromS3rho;
  calc_ufromSerho=&calc_ufromS3rho;
#else
  calc_SefromrhoT=&calc_S2fromrhoT;
  calc_Sefromrhou=&calc_S2fromrhou;
  calc_TfromSerho=&calc_TfromS2rho;
  calc_ufromSerho=&calc_ufromS2rho;
#endif

}

//**********************************************************************
// initialize gammagas at init 
//**********************************************************************

int
fill_arrays_at_init()
{
  int ii;
  for(ii=0;ii<Nloop_02;ii++) 
    { 
      int ix,iy,iz;
      ix=loop_02[ii][0];
      iy=loop_02[ii][1];
      iz=loop_02[ii][2];
      set_u_scalar(gammagas, ix, iy, iz, GAMMA);
    }
  return 0;
}


//**********************************************************************
//inverse of a general matrix using gsl
//**********************************************************************

int
inverse_matrix(ldouble *a, ldouble *ia, int N)
{
  gsl_matrix *m
    = gsl_matrix_alloc (N, N);
  gsl_matrix *im
    = gsl_matrix_alloc (N, N);
  int i,j;
  for(i=0;i<N;i++)
      for(j=0;j<N;j++)
	gsl_matrix_set(m,i,j,a[i*N+j]);

  gsl_permutation * p = gsl_permutation_alloc (N);

  int s;

  gsl_linalg_LU_decomp (m, p, &s);

  gsl_linalg_LU_invert (m, p, im);
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      ia[i*N+j]=gsl_matrix_get(im,i,j);

  gsl_matrix_free(m);
  gsl_matrix_free(im);
  gsl_permutation_free(p);
  
  return 0;
}

//**********************************************************************
//multiply 4by4 matrices
//**********************************************************************

int
multiply_44matrices(ldouble T1[][4],ldouble T2[][4],ldouble Tout[][4])
{
  int i,j,k;

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tout[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      Tout[i][j] += T1[i][k] * T2[k][j]; 
	    }
	}
    }
  return 0;
}

//**********************************************************************
//inverse 4by4 matrix
//**********************************************************************

int
inverse_44matrix(ldouble a[][4], ldouble ia[][4])
{

  ldouble mat[16],dst[16];
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      mat[i*4+j]=a[i][j];

  ldouble tmp[12];
  ldouble src[16];
  ldouble det;
  
  // transpose matrix
  for (i = 0; i <4; i++)
    {
      src[i]=mat[i*4];
      src[i+4]=mat[i*4+1];
      src[i+8]=mat[i*4+2];
      src[i+12]=mat[i*4+3];
    }
  
  // calculate pairs for first 8 elements (cofactors)
  tmp[0] = src[10] * src[15];
  tmp[1] = src[11] * src[14];
  tmp[2] = src[9] * src[15];
  tmp[3] = src[11] * src[13]; 
  tmp[4] = src[9] * src[14]; 
  tmp[5] = src[10] * src[13];
  tmp[6] = src[8] * src[15];
  tmp[7] = src[11] * src[12];
  tmp[8] = src[8] * src[14];
  tmp[9] = src[10] * src[12];
  tmp[10] = src[8] * src[13];
  tmp[11] = src[9] * src[12];
  
  // calculate first 8 elements (cofactors)
  dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7]; 
  dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7]; 
  dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7]; 
  dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7]; 
  dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6]; 
  dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6]; 
  dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3]; 
  dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3]; 
  dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3]; 
  dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3]; 
  dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];

  // calculate pairs for second 8 elements (cofactors)
  tmp[0] = src[2]*src[7]; 
  tmp[1] = src[3]*src[6];
  tmp[2] = src[1]*src[7];
  tmp[3] = src[3]*src[5]; 
  tmp[4] = src[1]*src[6];
  tmp[5] = src[2]*src[5];
  tmp[6] = src[0]*src[7];
  tmp[7] = src[3]*src[4];
  tmp[8] = src[0]*src[6];
  tmp[9] = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];

  // calculate second 8 elements (cofactors)
  dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15]; 
  dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15]; 
  dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15]; 
  dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15]; 
  dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14]; 
  dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10]; 
  dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10]; 
  dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8]; 
  dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8]; 
  dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9]; 
  dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9]; 
  dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];

  // calculate determinant
  det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

  // calculate matrix inverse
  det = 1/det; 

  if(isnan(det))
    //my_err("det in inverse 4x4 zero\n");
    return -1;

  for (j = 0; j < 16; j++)
    dst[j] *= det;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	ia[i][j]= dst[i*4+j];
	if(isnan(ia[i][j])) return -1;
      }

  return 0;
}

//**********************************************************************
//determinant of 4by4 matrix
//**********************************************************************

ldouble
determinant_44matrix(ldouble a[][4])
{

  ldouble mat[16],dst[16];
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      mat[i*4+j]=a[i][j];

  ldouble tmp[12];
  ldouble src[16];
  ldouble det;

  // transpose matrix
  for (i = 0; i <4; i++)
    {
      src[i]=mat[i*4];
      src[i+4]=mat[i*4+1];
      src[i+8]=mat[i*4+2];
      src[i+12]=mat[i*4+3];
    }
  
  // calculate pairs for first 8 elements (cofactors)
  tmp[0] = src[10] * src[15];
  tmp[1] = src[11] * src[14];
  tmp[2] = src[9] * src[15];
  tmp[3] = src[11] * src[13]; 
  tmp[4] = src[9] * src[14]; 
  tmp[5] = src[10] * src[13];
  tmp[6] = src[8] * src[15];
  tmp[7] = src[11] * src[12];
  tmp[8] = src[8] * src[14];
  tmp[9] = src[10] * src[12];
  tmp[10] = src[8] * src[13];
  tmp[11] = src[9] * src[12];
  
  // calculate first 8 elements (cofactors)
  dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7]; 
  dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
  dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7]; 
  dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7]; 
  dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
  dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7]; 
  dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6]; 
  dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6]; 
  dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3]; 
  dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3]; 
  dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3]; 
  dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
  dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3]; 
  dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
  dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
  dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
  
  // calculate pairs for second 8 elements (cofactors)
  tmp[0] = src[2]*src[7]; 
  tmp[1] = src[3]*src[6];
  tmp[2] = src[1]*src[7];
  tmp[3] = src[3]*src[5]; 
  tmp[4] = src[1]*src[6];
  tmp[5] = src[2]*src[5];
  tmp[6] = src[0]*src[7];
  tmp[7] = src[3]*src[4];
  tmp[8] = src[0]*src[6];
  tmp[9] = src[2]*src[4];
  tmp[10] = src[0]*src[5];
  tmp[11] = src[1]*src[4];
  
  // calculate second 8 elements (cofactors)
  dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15]; 
  dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
  dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15]; 
  dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15]; 
  dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
  dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15]; 
  dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
  dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14]; 
  dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
  dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10]; 
  dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10]; 
  dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8]; 
  dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8]; 
  dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9]; 
  dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9]; 
  dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
  
  // calculate determinant
  det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

  return det;
}

//**********************************************************************
//prints error message and gets chars
//**********************************************************************

int
my_err(char *message)
{
  if(PROCID==0)
    {
      char bufor[200];
      sprintf(bufor,"|err| : %s\n",message);
      printf("%s",bufor);
      getchar();
    }
  return 0;  
}

//**********************************************************************
//prints warning message and returns
//**********************************************************************

int
my_warning(char *message)
{
  if(PROCID==0)
    {
      char bufor[200];
      sprintf(bufor,"|warning| : %s\n",message);
      printf("%s",bufor);
    }
  return 0;  
}

//**********************************************************************
//getch to pause until character is entered
//**********************************************************************

int
getch()
{
  getchar();
  return 0;
}

//*********************************************************************
// minimum of 2 numbers
//**********************************************************************

ldouble my_min(ldouble a, ldouble b)
{
  if(a<b) return a;
  else return b;
}

//**********************************************************************
// minimum of a length N array
//**********************************************************************
ldouble my_min_N(ldouble *v,int N)
{
  ldouble min=v[0];
  int i;
  for(i=1;i<N;i++)
    if(v[i]<min) min=v[i];
  return min;
}

//**********************************************************************
// Maximum of a length N array
//**********************************************************************

ldouble my_max_N(ldouble *v,int N)
{
  ldouble max=v[0];
  int i;
  for(i=1;i<N;i++)
    if(v[i]>max) max=v[i];
  return max;
}

//**********************************************************************
// sign of a double
//**********************************************************************

ldouble my_sign(ldouble x)
{
  if(x>0.) return 1.;
  if(x<0.) return -1.;
  if(x==0.) return 1.;
  return 0;
}

//**********************************************************************
//* atan2 in [0,2pi] 
//**********************************************************************

ldouble
my_atan2(ldouble y, ldouble x)
{
  ldouble res=atan2(y,x);
  if(res<0.) res+=2.*M_PI;   
  return res;
}

//**********************************************************************
//Heaviside step function around 0
//x9 determines the sharpness and says where step function equals 0.95
//**********************************************************************

ldouble
step_function(ldouble x,ldouble x9)
{
  ldouble k=1.47222/x9;
  return 1./(1.+exp(-2.*k*x));
}

//**********************************************************************
// Arrange the N elements of ARRAY in random order.
//**********************************************************************

void shuffle_loop(int **array, size_t n)
{
  if (n > 1) {
    size_t i;
    for (i = 0; i < n - 1; i++) {
      size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
      int t[3] = {array[j][0],array[j][1],array[j][2]};
      array[j][0] = array[i][0];
      array[j][1] = array[i][1];
      array[j][2] = array[i][2];
      array[i][0] = t[0];
      array[i][1] = t[1];
      array[i][2] = t[2];
    }
  }
}

//****************************************************************
// prints tensor to screen
//****************************************************************

int
print_tensor(ldouble T[][4])
{
  int i;
  printf("============\n");
  for(i=0;i<4;i++)
    printf("%10.16e %10.16e %10.16e %10.16e\n",T[0][i],T[1][i],T[2][i],T[3][i]);
  printf("============\n");
  return 0;  
}

//****************************************************************
// prints tensor to screen 
//****************************************************************

int
print_NNtensor(ldouble **T,int N)
{
  int i,j;
  printf("============\n");
  for(i=0;i<N;i++)
    {
      for(j=0;j<N;j++)
	{
	  printf("%14.6e ",T[i][j]);
	}
      printf("\n");
    }
  printf("============\n");
  return 0;  
}


//****************************************************************
// prints metric to screen 
//****************************************************************

int
print_metric(ldouble T[][5])
{
  int i;
  printf("============\n");
  for(i=0;i<4;i++)
    printf("%10.3e %10.3e %10.3e %10.3e\n",T[i][0],T[i][1],T[i][2],T[i][3]);
  printf("============\n");
  return 0;  
}


//*****************************************************************
// prints 4vector to screen
//*****************************************************************

int
print_4vector(ldouble v[4])
{
  int i;
  printf("============\n");
  printf("%10.8e %10.8e %10.8e %10.8e\n",v[0],v[1],v[2],v[3]);
  printf("============\n");
  return 0;  
}


//****************************************************************
// prints Nvvector to screen 
//****************************************************************

int
print_NVvector(ldouble *v)
{
  print_Nvector(v,NV);
  return 0;
}


//****************************************************************
// prints primitives to screen 
//****************************************************************

int
print_primitives(ldouble *p)
{
  printf("\n");
  printf("rho = %.15e\n",p[RHO]);
  printf("ugas = %.15e\n",p[UU]);
  printf("u^1 = %.15e\n",p[VX]);
  printf("u^2 = %.15e\n",p[VY]);
  printf("u^3 = %.15e\n",p[VZ]);
  printf("S = %.15e\n",p[ENTR]);
#ifdef MAGNFIELD
  printf("B^1 = %.15e\n",p[B1]);
  printf("B^2 = %.15e\n",p[B2]);
  printf("B^3 = %.15e\n",p[B3]);
#endif
#ifdef RADIATION
  printf("Erf = %.15e\n",p[EE0]);
  printf("ur^1 = %.15e\n",p[FX0]);
  printf("ur^2 = %.15e\n",p[FY0]);
  printf("ur^3 = %.15e\n",p[FZ0]);
#ifdef EVOLVEPHOTONNUMBER
  printf("nph = %.15e\n",p[NF0]);
#endif
#endif
#ifdef EVOLVEELECTRONS
  printf("entre = %.15e\n",p[ENTRE]);
  printf("entri = %.15e\n",p[ENTRI]);
#endif
#ifdef RELELECTRONS
  int ib;
  for(ib=0;ib<NRELBIN;ib++)
    printf("rele(%d) = %.15e\n",ib,p[NEREL(ib)]);
#endif

   printf("\n");

  return 0;
}

//****************************************************************
// prints conserveds to screen
//****************************************************************

int
print_conserved(ldouble *u)
{
  printf("\n");
  printf("rho u^t = %.15e\n",u[RHO]);
  printf("T^t_t + rho u^t = %.15e\n",u[UU]);
  printf("T^t_1 = %.15e\n",u[VX]);
  printf("T^t_2 = %.15e\n",u[VY]);
  printf("T^t_3 = %.15e\n",u[VZ]);
  printf("S u^t = %.15e\n",u[ENTR]);
#ifdef MAGNFIELD
  printf("B^1 = %.15e\n",u[B1]);
  printf("B^2 = %.15e\n",u[B2]);
  printf("B^3 = %.15e\n",u[B3]);
#endif
#ifdef RADIATION
  printf("R^t_t = %.15e\n",u[EE0]);
  printf("R^t_1 = %.15e\n",u[FX0]);
  printf("R^t_2 = %.15e\n",u[FY0]);
  printf("R^t_3 = %.15e\n",u[FZ0]);
#ifdef EVOLVEPHOTONNUMBER
  printf("nph u^t = %.15e\n",u[NF0]);
#endif
#endif
#ifdef EVOLVEELECTRONS
  printf("entre u^t = %.15e\n",u[ENTRE]);
 printf("entri u^t = %.15e\n",u[ENTRI]);
#endif

#ifdef RELELECTRONS
  int ib;
  for(ib=0;ib<NRELBIN;ib++)
    printf("rele(%d) = %.15e\n",ib,u[NEREL(ib)]);
#endif
  printf("\n");
 
  return 0;
}


//****************************************************************
// prints Nvector to screen
//****************************************************************

int
print_Nvector(ldouble v[4],int N)
{
  int i;
  printf("============\n");
  for(i=0;i<N;i++)
  printf("%14.6e ",v[i]);
  printf("\n============\n");
  return 0;  
}


//*****************************************************************
// get clock time 
//*****************************************************************

int
my_clock_gettime(void* tsptr)
{
struct timespec *ts
    = (struct timespec *) tsptr;

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;

#else
  clock_gettime(CLOCK_REALTIME, ts);
#endif
  return 0;
}

//**********************************************************************
//decomposes velocities into cartesian from spherical
//**********************************************************************

int
decompose_vels(ldouble *pp,int velidx, ldouble v[4],void *ggg,  void *gggBL)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  struct geometry *geomBL
    = (struct geometry *) gggBL;
  int iv;
  int velprim=VELPRIM; if(velidx==FX) velprim=VELPRIMRAD;
  ldouble r,th,ph,ucon[4];
  DLOOPA(iv) ucon[iv]=pp[velidx-1+iv];
  ucon[0]=0.;
  r=geomBL->xx;	     th=geomBL->yy;	     ph=geomBL->zz;
  //to ortonormal cartesian
  ucon[1]*=sqrt(geom->gg[1][1]);
  ucon[2]*=sqrt(geom->gg[2][2]);
  ucon[3]*=sqrt(geom->gg[3][3]);

  v[1] = sin(th)*cos(ph)*ucon[1] + cos(th)*cos(ph)*ucon[2] - sin(ph)*ucon[3];
  v[2] = sin(th)*sin(ph)*ucon[1] + cos(th)*sin(ph)*ucon[2] + cos(ph)*ucon[3];
  v[3] = cos(th)*ucon[1] - sin(th)*ucon[2];
  return 0;
}

//**********************************************************************
// get cell size in x,y,z
//**********************************************************************


// ANDREW TODO DEBUG
// get the size of a cell dx in 3 dimensions in OUTCOOORDS
int get_cellsize_out(int ix, int iy, int iz, ldouble dx[3])
{
   ldouble xx1[4], xx2[4];
#ifdef PRECOMPUTE_MY2OUT
   int ii;
   xx1[0] = 0.; xx2[0] = 0.;
   for(ii=0;ii<3;ii++)
   {
     xx1[ii+1] = get_xbout_xface(ii,ix,iy,iz);
     xx2[ii+1] = get_xbout_xface(ii,ix+1,iy,iz);
   }

   dx[0] = fabs(xx2[1] - xx1[1]);

   for(ii=0;ii<3;ii++)
   {
     xx1[ii+1] = get_xbout_yface(ii,ix,iy,iz);
     xx2[ii+1] = get_xbout_yface(ii,ix,iy+1,iz);
   }
   
   dx[1] = fabs(xx2[2] - xx1[2]);

   for(ii=0;ii<3;ii++)
   {
     xx1[ii+1] = get_xbout_zface(ii,ix,iy,iz);
     xx2[ii+1] = get_xbout_zface(ii,ix,iy,iz+1);
   }
   
   dx[2] = fabs(xx2[3] - xx1[3]);
          
#else
   xx1[0] = 0.; xx1[1] = get_xb(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_x(iz, 2);
   xx2[0] = 0.; xx2[1] = get_xb(ix+1, 0);xx2[2] = get_x(iy, 1); xx2[3] = get_x(iz, 2);
   coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
   coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);

   dx[0] = fabs(xx2[1] -xx1[1]);

   xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_xb(iy, 1); xx1[3] = get_x(iz, 2);
   xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_xb(iy+1, 1); xx2[3] = get_x(iz, 2);
   coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
   coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);

   dx[1] = fabs(xx2[2] - xx1[2]);

   xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_xb(iz, 2);
   xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_x(iy, 1); xx2[3] = get_xb(iz+1, 2);
   coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
   coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);

   dx[2] = fabs(xx2[3] - xx1[3]);

#endif
   return 0;

}

int
get_cellsize_arb(int ix,int iy,int iz,ldouble dx[3],int COORDS)
{
  
  ldouble xx1[4],xx2[4]; 
  xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_x(iz,2);
  xx2[0]=0.;xx2[1]=get_xb(ix+1,0);xx2[2]=get_x(iy,1);xx2[3]=get_x(iz,2);
  coco_N(xx1,xx1,MYCOORDS,COORDS);
  coco_N(xx2,xx2,MYCOORDS,COORDS);
  dx[0]=fabs(xx2[1]-xx1[1]);
  xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_x(iz,2);
  xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_xb(iy+1,1);xx2[3]=get_x(iz,2);
  coco_N(xx1,xx1,MYCOORDS,COORDS);
  coco_N(xx2,xx2,MYCOORDS,COORDS);
  dx[1]=fabs(xx2[2]-xx1[2]);
  xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_xb(iz,2);
  xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_x(iy,1);xx2[3]=get_xb(iz+1,2);
  coco_N(xx1,xx1,MYCOORDS,COORDS);
  coco_N(xx2,xx2,MYCOORDS,COORDS);
  dx[2]=fabs(xx2[3]-xx1[3]);
  if(NZ==1) dx[2]=2.*M_PI;

  return 0;
}

//**********************************************************************
//sign of a 3 perturbation {abc}
//**********************************************************************

int
pertsign_3d(int a,int b,int c)
{
  a+=1;
  b+=1;
  c+=1;
  if(a==1 && b==2 && c==3) return 1;
  if(a==1 && b==3 && c==2) return -1;
  if(a==2 && b==1 && c==3) return -1;
  if(a==2 && b==3 && c==1) return 1;
  if(a==3 && b==1 && c==2) return 1;
  if(a==3 && b==2 && c==1) return -1;

  return 0;
}

//**********************************************************************
//sign of a 4 perturbation {abcd}
//**********************************************************************

int
pertsign_4d(int a,int b,int c, int d)
{
  a+=1;
  b+=1;
  c+=1;
  d+=1;
  if(a==1 && b==2 && c==3 && d==4)  return 1;
  if(a==2 && b==1 && c==3 && d==4)  return -1;
  if(a==1 && b==3 && c==2 && d==4)  return -1;
  if(a==3 && b==1 && c==2 && d==4)  return 1;
  if(a==2 && b==3 && c==1 && d==4)  return 1;
  if(a==3 && b==2 && c==1 && d==4)  return -1;
  if(a==1 && b==2 && c==4 && d==3)  return -1;
  if(a==2 && b==1 && c==4 && d==3)  return 1;
  if(a==1 && b==4 && c==2 && d==3)  return 1;
  if(a==4 && b==1 && c==2 && d==3)  return -1;
  if(a==2 && b==4 && c==1 && d==3)  return -1;
  if(a==4 && b==2 && c==1 && d==3)  return 1;
  if(a==1 && b==3 && c==4 && d==2)  return 1;
  if(a==3 && b==1 && c==4 && d==2)  return -1;
  if(a==1 && b==4 && c==3 && d==2)  return -1;
  if(a==4 && b==1 && c==3 && d==2)  return 1;
  if(a==3 && b==4 && c==1 && d==2)  return 1;
  if(a==4 && b==3 && c==1 && d==2)  return -1;
  if(a==2 && b==3 && c==4 && d==1)  return -1;
  if(a==3 && b==2 && c==4 && d==1)  return 1;
  if(a==2 && b==4 && c==3 && d==1)  return 1;
  if(a==4 && b==2 && c==3 && d==1)  return -1;
  if(a==3 && b==4 && c==2 && d==1)  return -1;
  if(a==4 && b==3 && c==2 && d==1)  return 1;

  return 0;
}

//**********************************************************************
//Levi-Civita symbol eta_{abcd}
//**********************************************************************

int
epsilon_4d(int a,int b,int c, int d)
{
  if(a==b || a==c || a==d) return 0;
  if(b==c || b==d) return 0;
  if(c==d) return 0;

  return pertsign_4d(a,b,c,d);
}

//**********************************************************************
//Levi-Civita tensor: (-g)^−1/2  [ijk]
//**********************************************************************
ldouble
epsilon_4d_tensor(int a,int b,int c, int d, void *ggg)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  if(a==b || a==c || a==d) return 0;
  if(b==c || b==d) return 0;
  if(c==d) return 0;

  return 1./geom->gdet*(ldouble)pertsign_4d(a,b,c,d);
}

//**********************************************************************
//Levi-Civita tensor: (-g)^−1/2  [ijk]
//**********************************************************************

ldouble
epsilon_3d_tensor(int a,int b,int c, void *ggg)
{
   struct geometry *geom
     = (struct geometry *) ggg;

  if(a==b || a==c || b==c) return 0;

  return 1./geom->gdet*(ldouble)pertsign_3d(a,b,c);
}


//**********************************************************************
//calls gnuplot
//**********************************************************************

int
convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
{
  char buf1[50],buf2[50];  
  sprintf(buf1,"plot.gp.%d",PROCID);
  FILE *fgnu=fopen(buf1,"w");

#ifdef PR_OUT2GIF_2D
  #include PR_OUT2GIF_2D   //PROBLEMS/XXX/out2gid_2d.c
#endif

  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  sprintf(buf2,"gnuplot %s",buf1);
  int i=system(buf2);
  return 0;
}

int
convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t)
{
  FILE *fgnu=fopen("plot.gp","w");
  char bufor[50];

#ifdef PR_OUT2GIF_1D
  #include PR_OUT2GIF_1D   //PROBLEMS/XXX/out2gid_2d.c
#endif
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
}

