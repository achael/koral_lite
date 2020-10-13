/*! \file u2p.c
 \brief MHD Conserved to primitives conversion
 */

#include "ko.h"

static FTYPE dpdWp_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq,FTYPE gamma);
static FTYPE compute_idwmrho0dp(FTYPE wmrho0,FTYPE gamma);
static FTYPE compute_idrho0dp(FTYPE wmrho0);
static int f_u2p_hot(ldouble Wp, ldouble* cons,ldouble *f,ldouble *df,ldouble *err,ldouble pgamma);
static FTYPE pressure_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma);
static FTYPE compute_inside_entropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma);
static FTYPE compute_specificentropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma);
static FTYPE compute_dspecificSdwmrho0_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma);
static FTYPE compute_dspecificSdrho_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0, FTYPE gamma);
static int f_u2p_entropy(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df, ldouble *err,ldouble pgamma);

//**********************************************************************
//calculates primitives in given cell basing on global array u[]
// type: not used
// int setflags -- is always set to 1 in the current code
//**********************************************************************

int
calc_primitives(int ix,int iy,int iz,int type,int setflags)
{
  int verbose=0;
  int iv,u2pret,u2pretav;
  ldouble uu[NV],uuav[NV],pp[NV],ppav[NV];
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  //temporary local arrays
  gg=geom.gg;
  GG=geom.GG;
  gdet=geom.gdet;gdetu=gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif
  
  int corrected[3]={0,0,0}, fixups[2]={0,0};
  for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,ix,iy,iz);
    pp[iv]=get_u(p,iv,ix,iy,iz);
  }
  
  if(setflags)
  {
    set_cflag(ENTROPYFLAG,ix,iy,iz,0);
    set_cflag(ENTROPYFLAG2,ix,iy,iz,0);
  }

  //u to p inversion is done here
  if(is_cell_corrected_polaraxis(ix,iy,iz))
  {
    u2p_solver_Bonly(uu,pp,&geom); // invert only the magnetic field, the rest will be overwritten
  }
  else
  {
    u2p(uu,pp,&geom,corrected,fixups,type); // regular inversion
  }

  //set flags for entropy solver
  if(corrected[0]==1 && setflags) //hd correction - entropy solver
  {
    set_cflag(ENTROPYFLAG,ix,iy,iz,1);
  }
  
  if(corrected[2]==1 && setflags) //borrowing energy from radiation didn't work
  {  
    set_cflag(ENTROPYFLAG2,ix,iy,iz,1);
  }
  
  //check hd floors
  int floorret=0;
  
  if(is_cell_active(ix,iy,iz) && !is_cell_corrected_polaraxis(ix,iy,iz))
  {
    floorret=check_floors_mhd(pp,VELPRIM,&geom);
  }
  
  if(floorret<0.)
  {
    corrected[0]=1;
  }

  //check rad floors
#ifdef RADIATION
  floorret=0;
  if(is_cell_active(ix,iy,iz) &&  !is_cell_corrected_polaraxis(ix,iy,iz))
  {
    floorret=check_floors_rad(pp,VELPRIMRAD,&geom);
  }
  
  if(floorret<0.)
  {
    corrected[1]=1;
  }
#endif
  
  //set new primitives and conserved
  for(iv=0;iv<NV;iv++)
  { 
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }
  
  //set flags for fixups of unsuccessful cells
  if(setflags)
  {
    if(fixups[0]>0)
    {
      set_cflag(HDFIXUPFLAG,ix,iy,iz,1);
      global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]++;
    }
    else
      set_cflag(HDFIXUPFLAG,ix,iy,iz,0);
    
    if(fixups[1]>0)
    {
      set_cflag(RADFIXUPFLAG,ix,iy,iz,-1);
      global_int_slot[GLOBALINTSLOT_NTOTALRADFIXUPS]++;
    }
    else
      set_cflag(RADFIXUPFLAG,ix,iy,iz,0); 
  }
  
  return 0;
} 


//**********************************************************************
//high-level u2p solver
// type: not used
//**********************************************************************

int
u2p(ldouble *uu0, ldouble *pp, void *ggg, int corrected[3], int fixups[2], int type)
{
  struct geometry *geom
  = (struct geometry *) ggg;
  
  ldouble uu[NV];
  int iv;
  PLOOP(iv) uu[iv]=uu0[iv];
  
  ldouble (*gg)[5], (*GG)[5], gdet, gdetu, gdetu_inv;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif
  gdetu_inv = 1. / gdetu;
  
  int verbose=0;
  int hdcorr=0;
  int radcor=0;
  corrected[0]=corrected[1]=0;
  fixups[0]=fixups[1]=0;
  
  int u2pret,u2pentrret,ret;
  ldouble ppbak[NV];
  for(u2pret=0;u2pret<NV;u2pret++)
    ppbak[u2pret]=pp[u2pret];

  //************************************
  //magneto-hydro part
  //************************************
  ldouble u0=pp[1];
  
  //************************************
  //hot hydro - conserving energy
  ret=0;
  u2pret=-1;
  
  //test
  ldouble ppold[NV];
  ldouble ppentr[NV];
  PLOOP(iv)
  {
    ppold[iv]=pp[iv];
    ppentr[iv]=-1.; // negative value indicates that entropy inversion not yet calculated
  }
  
  //negative uu[0] = rho u^t
  if(uu[0] * gdetu_inv < 0.)
  {
    int gix,giy,giz;
    mpi_local2globalidx(geom->ix,geom->iy,geom->iz,&gix,&giy,&giz);
    if(verbose) printf("%4d > %4d %4d %4d > NEGUU  > neg uu[0] - requesting fixup\n",PROCID,gix,giy,giz);
    pp[0]=RHOFLOOR; //used when not fixing up
    uu[0]=RHOFLOOR*gdetu;
    ret=-2;    //to request fixup
               //ANDREW -- but ret=-1 if energy inversion failes but entropy inversion does not!
               //ANDREW -- do we always want a fixup if we have negative uu[0] ? 
    u2pret=-1; // indicates that inversion is needed
    
#ifndef SWAPPAPC
    global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]++;  //but count as fixup
#endif
  }
    
  if(u2pret!=0)  // u2pret=-1 at this stage, so this is always satisfied
  {
#ifdef ENFORCEENTROPY  
    u2pret=-1;  //skip hot energy-conserving inversion and go to entropy inversion
#else
    u2pret = u2p_solver(uu,pp,ggg,U2P_HOT,0);  // invert using the hot energy equation    
#endif //ENFORCEENTROPY
  }
  
  if(ALLOWENTROPYU2P)  // Inversion with entropy equation -- on by default (see choices.h)
  {
    if(u2pret<0)  // true if energy equation failed, or if energy equation was not required (because ENFORCEENTROPY is defined)
    {
      ret=-1;
      
      if(verbose>2 )
      {
        printf("u2p_entr     >>> %d %d <<< %d >>> %e > %e\n",geom->ix + TOI, geom->iy + TOJ,u2pret,u0,pp[1]);
      }
      
      //************************************
      //entropy solver - conserving entropy
      if(ppentr[RHO]<0.) //if not yet calculated
      {
        u2pret=u2p_solver(uu,pp,ggg,U2P_ENTROPY,0);  // invert using entropy equation
      }
      
      if(u2pret<0)
      {
        ret=-2;
        
        if(verbose>1)
        {
          printf("u2p_entr err No. %4d > %e %e %e > %e %e > %4d %4d %4d\n",u2pret,uu[0],uu[1],uu[5],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
        }
	
      } // if(u2pret<0) // second time -- entropy eqn
    } // if(u2pret<0) // first time -- energy eqn
  }  // if(ALLOWENTROPYU2P)
  
  if(u2pret<0)  // entropy equation also failed
  {
 
    //leaving primitives unchanged - should not happen
    if(verbose>1 || 1)
    {
      printf("%4d > %4d %4d %4d > MHDU2PFAIL > u2p prim. unchanged > %d \n",PROCID,geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,u2pret);
    }
    ret=-3;
    for(u2pret=0;u2pret<NV;u2pret++)
      pp[u2pret]=ppbak[u2pret];
  }
  
  if(ret<0) //to update conserved
    hdcorr=1;  
  if(ret<-1) //request fixup when entropy failed
    fixups[0]=1;
  else
    fixups[0]=0;
  
  //************************************
  //radiation part
  //************************************
  
  corrected[2]=0;
  
#ifdef RADIATION  
#ifdef BALANCEENTROPYWITHRADIATION
  
  //trying to balance gain of energy because of entropy inversion
  //by borrowing from the radiation field
  if(ret==-1) //entropy u2p was used in MHD part
  {
    ldouble uunew[NV],ppnew[NV];
    PLOOP(iv) { uunew[iv]=uu[iv]; ppnew[iv]=pp[iv]; }
    p2u_mhd(pp,uunew,geom);
    ldouble dugas = uunew[UU] - uu[UU];  //this much energy was introduced
    if(fabs(dugas)<0.1*fabs(uunew[EE0])) //correction relatively small - is this general enough?
    {
      uunew[EE0]-=dugas; //balancing with radiation
      u2p_rad(uunew,ppnew,geom,&radcor);
    }
    else
      radcor=1;
    
    if(radcor==0) //there was enough energy to borrow from and uunew inverts with hot
    {
      PLOOP(iv)
      uu[iv]=uunew[iv];
      //printf("entropy correction worked at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
    }
    else
    {
      corrected[2]=1; //entropy correction didn't work
      //printf("entropy correction didn't work at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
    }
  }
#endif //BALANCEENTROPYWITHRADIATION

  //Do the radiative inversion from u2p_rad.c
  u2p_rad(uu,pp,geom,&radcor);

#endif // RADIATION

  //************************************  
  //output
  //************************************
  
  //rad fixups only for critical failure in implicit
  if(radcor>0)     
    fixups[1]=1;
  else
    fixups[1]=0;
    
  if(hdcorr>0) corrected[0]=1;
  if(radcor>0) corrected[1]=1;
  
  return ret;
} 


//**********************************************************************
//checks if hydro primitives make sense
//**********************************************************************

int
check_floors_mhd(ldouble *pp, int whichvel,void *ggg)
{

  int verbose=0;
  int ret=0;
  int iv;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble uu[NV],pporg[NV];
  p2u_mhd(pp,uu,ggg);
 
  //**********************************************************************
  //rho too small
  if(pp[0]<RHOFLOOR) 
  {
      for(iv=0;iv<NVMHD;iv++)
      {
	pporg[iv]=pp[iv];
      }

      if(verbose ) printf("hd_floors CASE 1 at %d %d %d | %d %d %d (%e) | tijk: %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,geom->ix,geom->iy,geom->iz,pp[0],TI,TJ,TK);
      pp[0]=RHOFLOOR; 
     
      ret=-1; 
  }

  //**********************************************************************
  //rho too small, use floor scaling as r^-3/2
#ifdef RHOFLOOR_BH
  ldouble xxBL[4];
  #ifdef PRECOMPUTE_MY2OUT
  get_xxout(geom->ix, geom->iy, geom->iz, xxBL);
  #else
  coco_N(geom->xxvec,xxBL,MYCOORDS,BLCOORDS);
  #endif

  ldouble rr = xxBL[1] / rhorizonBL;
  ldouble rhofloor = RHOFLOOR_BH_NORM / sqrt(rr*rr*rr);
  if(pp[0]<rhofloor) 
  {
    for(iv=0;iv<NVMHD;iv++)
    {
      pporg[iv]=pp[iv];
    }
      
      if(verbose ) printf("hd_floors BH CASE 1 at %d %d (%e)\n",geom->ix+TOI,geom->iy+TOJ,pp[0]);
      pp[0]=rhofloor;

      ret=-1; 
  }
#endif

  //***********************************************************************
  //rho too small, BH-disk like: Here we use the initial atmosphere as the floor on both density and pressure
#ifdef RHOFLOOR_INIT
  ldouble xxBL[4], rout = 2.;

  #ifdef PRECOMPUTE_MY2OUT
  get_xxout(geom->ix, geom->iy, geom->iz, xxBL);
  #else
  coco_N(geom->xxvec,xxBL,MYCOORDS,BLCOORDS);
  #endif

  ldouble rr = xxBL[1] / rout;
  ldouble rhofloor = RHOATMMIN / sqrt(rr*rr*rr);
  ldouble uintfloor = UINTATMMIN / sqrt(rr*rr*rr*rr*rr);
  if(pp[0] < rhofloor || pp[1] < uintfloor)
  {
    for(iv = 0; iv < NVMHD; iv++)
    {
      pporg[iv] = pp[iv];
    }
    
    if(verbose ) printf("hd_floors BH CASE INIT at %d %d (%e)\n",geom->ix+TOI, geom->iy+TOJ, pp[0]);
    pp[0] = rhofloor;
    pp[1] = uintfloor;
    
    ret=-1;
  }
#endif
  
  //**********************************************************************
  //too cold
  if(pp[1]<UURHORATIOMIN*pp[0]) 
  {
    for(iv=0;iv<NVMHD;iv++)
    {
      pporg[iv]=pp[iv];
    }

    if(verbose) {printf("hd_floors CASE 2 at (%d,%d,%d | %d,%d,%d): %e %e | tijk: %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,geom->ix,geom->iy,geom->iz,pp[0],pp[1],TI,TJ,TK);}//getchar();}
      pp[1]=UURHORATIOMIN*pp[0]; //increasing uint

              
    ret=-1;
  }

  //**********************************************************************
  //too hot
  if(pp[1]>UURHORATIOMAX*pp[0]) 
  {
    for(iv=0;iv<NVMHD;iv++)
    {
      pporg[iv]=pp[iv];
    }

    pp[1]=UURHORATIOMAX*pp[0]; //decreasing uint

    ret=-1;      
    if(verbose ) printf("hd_floors CASE 3 at (%d,%d,%d): %e %e\n",geom->ix+TOI,geom->iy+TOJ,geom->iz,pp[0],pp[1]);
  }
  
  //**********************************************************************
  //too magnetized
#ifdef MAGNFIELD
  ldouble ucond[4],ucovd[4];
  ldouble bcond[4],bcovd[4],bsq,magpre;
  ldouble etacon[4],etarel[4];
  for(iv=1;iv<4;iv++)
    ucond[iv]=pp[1+iv];
  calc_ucon_ucov_from_prims(pp, geom, ucond, ucovd);
  calc_bcon_bcov_bsq_from_4vel(pp, ucond, ucovd, geom, bcond, bcovd, &bsq);
  magpre = 0.5 * bsq;
  
  calc_normalobs_ncon(GG, geom->alpha, etacon);
  conv_vels_ut(etacon,etarel,VEL4,VELPRIM,gg,GG);

  if(magpre>B2RHORATIOMAX*pp[RHO]) 
  {
    if(verbose) printf("mag_floors CASE 2 at (%d,%d,%d): %e %e\n",geom->ix+TOI,geom->iy+TOJ,geom->iz,pp[RHO],magpre);
    ldouble f=magpre/(B2RHORATIOMAX*pp[RHO]);

    for(iv=0;iv<NVMHD;iv++)
    {
      pporg[iv]=pp[iv];
    }

#if (B2RHOFLOORFRAME==ZAMOFRAME) //new mass in ZAMO

      ldouble dpp[NV],duu[NV];
      ldouble drho=pp[RHO]*(f-1.);
   
      for(iv=0;iv<NVMHD;iv++)
	dpp[iv]=0.0;

      dpp[RHO]=drho;
      //do not inject energy - just density
      dpp[UU]=0.;
      dpp[VX] = etarel[1];
      dpp[VY] = etarel[2];
      dpp[VZ] = etarel[3];
      dpp[ENTR] = 0.;
      dpp[B1] = dpp[B2] = dpp[B3] = 0.;

      p2u_mhd(dpp,duu,geom);
 
      for(iv=0;iv<NVMHD;iv++)
      {
	uu[iv]+=duu[iv];
      }

      int rettemp=0;
      rettemp=u2p_solver(uu,pp,geom,U2P_HOT,0); 
      if(rettemp<0)
	rettemp=u2p_solver(uu,pp,geom,U2P_ENTROPY,0); 
      
      if(rettemp<0) 
      {
#ifdef BHDISK_PROBLEMTYPE
	//if(geom->ix+TOI>5) //report only outside horizon
#endif
	  printf("u2p failed after imposing bsq over rho floors at %d %d %d with gamma=%f\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,get_u_scalar(gammagas,geom->ix,geom->iy,geom->iz));

#ifdef B2RHOFLOOR_BACKUP_FFFRAME
       // Backup bsq/rho floor -- if zamo frame fails, do fluid frame instead of crashing 
       for(iv=0;iv<NVMHD;iv++)
         pp[iv]=pporg[iv];
       pp[RHO]*=f;
       pp[UU]*=f;
#else
       print_primitives(pp);
       exit(-1);
#endif
      }
    
#elif(B2RHOFLOORFRAME==FFFRAME) //new mass in fluid frame
      pp[RHO]*=f;
      pp[UU]*=f;

#endif //B2RHOFLOORFRAME==ZAMOFRAME

#ifdef EVOLVEELECTRONS

      //keep energy density in ions and electrons fixed after modifying B
      ldouble Tg,Te,Ti,ptot,uint,theta;

      //get temperatures after explicit
      Tg=calc_PEQ_Teifrompp(pporg,&Te,&Ti,geom->ix,geom->iy,geom->iz);
    
      ldouble ne=calc_thermal_ne(pporg); //thermal only
      ldouble pe=K_BOLTZ*ne*Te;
      ldouble gammae=GAMMAE;
#ifdef CONSISTENTGAMMA
#ifndef FIXEDGAMMASPECIES
      gammae=calc_gammaintfromtemp(Te,ELECTRONS);
#endif
#endif
      ldouble ue=pe/(gammae-1.);  
      ldouble ni=pporg[RHO]/MU_I/M_PROTON;
      ldouble pi=K_BOLTZ*ni*Ti;
      ldouble gammai=GAMMAI;
#ifdef CONSISTENTGAMMA
#ifndef FIXEDGAMMASPECIES
      gammai=calc_gammaintfromtemp(Ti,IONS);
#endif
#endif
      ldouble ui=pi/(gammai-1.);

      //calculate new entropy using pp[]
      ldouble mass, gamma, Tenew,Tinew, Senew, Sinew;
      ldouble  n, pnew;

      //electrons
      mass = M_ELECTR;
      gamma = GAMMAE;
      n = calc_thermal_ne(pp);
#ifdef CONSISTENTGAMMA
      Tenew=solve_Teifromnmu(n, mass, ue,ELECTRONS); //solves in parallel for gamma and temperature
      theta=K_BOLTZ*Tenew/mass;  
#ifndef FIXEDGAMMASPECIES
      gamma=calc_gammaintfromtheta(theta); //the same gamma as just solved     
#endif
#endif
      pnew=(ue)*(gamma-1.);
      Tenew=pnew/K_BOLTZ/n;

      ldouble rhoe=n*MU_E*M_PROTON;
      Senew=calc_SefromrhoT(rhoe,Tenew,ELECTRONS);

      pp[ENTRE]=Senew;

      //ions
      mass = M_PROTON;
      gamma = GAMMAI;
      n = calc_thermal_ne(pp);
#ifdef CONSISTENTGAMMA
      Tinew=solve_Teifromnmu(n, mass, ui,IONS); //solves in parallel for gamma and temperature
      theta=K_BOLTZ*Tinew/mass;
#ifndef FIXEDGAMMASPECIES
      gamma=calc_gammaintfromtheta(theta); //the same gamma as just solved     
#endif
#endif
      
      pnew=(ui)*(gamma-1.);
      Tinew=pnew/K_BOLTZ/n;
      Sinew=calc_SefromrhoT(pp[RHO],Tinew,IONS);
      
      pp[ENTRI]=Sinew;
#endif  //EVOLVEELECTRONS

      ret=-1;      
  } //if(magpre>B2RHORATIOMAX*pp[RHO]) 
#endif //MAGNFIELD

  //**********************************************************************
  //too fast
  if(VELPRIM==VELR) 
  {
      ldouble qsq=0.;
      int i,j;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=pp[UU+i]*pp[UU+j]*gg[i][j];
      ldouble gamma2=1.+qsq;
      if(gamma2>GAMMAMAXHD*GAMMAMAXHD)
      {
	  ldouble qsqmax=GAMMAMAXHD*GAMMAMAXHD-1.;
	  ldouble A=sqrt(qsqmax/qsq);
	  for(j=1;j<4;j++)
	    pp[UU+j]*=A;
	  ret=-1;
	  if(verbose )
	  {
	      printf("hd_floors CASE 4 at (%d,%d,%d): %e",geom->ix+TOI,geom->iy+TOJ,geom->iz,sqrt(gamma2));
	      qsq=0.;
	      for(i=1;i<4;i++)
		for(j=1;j<4;j++)
		  qsq+=pp[UU+i]*pp[UU+j]*gg[i][j];
	      gamma2=1.+qsq;
	      printf(" -> %e\n",sqrt(gamma2));
	  }
      }
  }
  
  //TODO: implement checks for other VELPRIM

  //**********************************************************************  
  //Species temperature floors/ceilings
#ifdef EVOLVEELECTRONS
  ldouble mue,mui;
  mui=MU_I;
  mue=MU_E;
  ldouble Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],geom->ix,geom->iy,geom->iz);

  //Electrons
  ldouble Teloc,Teloc0;
  ldouble neth=calc_thermal_ne(pp);
  ldouble rhoeth=MU_E*M_PROTON*neth;
  Teloc=calc_TfromSerho(pp[ENTRE],rhoeth,ELECTRONS,geom->ix,geom->iy,geom->iz);
  Teloc0=Teloc;
  
  // absolute floor
  if(Teloc<TEMPEMINIMAL)
  {
    Teloc=TEMPEMINIMAL;
  }

  // relative floor
  ldouble Teminimal=TEMPEMINIMALFRACTION*Tgas;
  if(Teloc<Teminimal)
  {  
    Teloc=Teminimal;
  }
  
  // ceiling
  ldouble Temaximal=TEMPEMAXIMALFRACTION*Tgas;
  if(Teloc>Temaximal)
  {
    Teloc=Temaximal;
  }
  
  //Ion Temperature
  ldouble Tiloc,Tiloc0;
  Tiloc=calc_TfromSerho(pp[ENTRI],pp[RHO],IONS,geom->ix,geom->iy,geom->iz);
  Tiloc0=Tiloc;
  
  // absolute floor
  if(Tiloc<TEMPIMINIMAL)
  {
    Tiloc=TEMPIMINIMAL;
  }
 
  // relative floor
  ldouble Timinimal=TEMPIMINIMALFRACTION*Tgas;
  if(Tiloc<Timinimal)
  {
      Tiloc=Timinimal;
  }

  // celing 
  ldouble Timaximal=TEMPIMAXIMALFRACTION*Tgas;
  if(Tiloc>Timaximal)
  {
    Tiloc=Timaximal;
  }
  
  if(Teloc!=Teloc0) //update temperature of electrons
  {  
    pp[ENTRE]=calc_SefromrhoT(rhoeth,Teloc,ELECTRONS);
    ret=-1;
  }

  if(Tiloc!=Tiloc0) //update temperature of ioms
  { 
    pp[ENTRI]=calc_SefromrhoT(pp[RHO],Tiloc,IONS);
    ret=-1;
  }
      
#ifdef RELELECTRONS
  int ie;
  ldouble ne_relel,ne_tot,uint_relel,uint_tot,p_relel,p_tot;
 
  //No negative rel. electron numbers
  for (ie=0; ie<NRELBIN; ie++)
  {
    if (pp[NEREL(ie)] < 0.0) 
    {
      pp[NEREL(ie)] = 0.0;
    }
  }

  ldouble relfracn, relfracu,relfracp;
  
  //Not too many rel. electrons
  ne_relel = calc_relel_ne(pp);
  ne_tot = pp[RHO]/MU_E/M_PROTON;
  relfracn = ne_relel/ne_tot;
  if (relfracn > MAX_RELEL_FRAC_N) 
  { 
    for (ie=0; ie<NRELBIN; ie++)
    {  
      pp[NEREL(ie)] *= (MAX_RELEL_FRAC_N/relfracn);
    }
  }

  //Not too much rel electron energy
  uint_relel = calc_relel_uint(pp);
  uint_tot = pp[UU];
  relfracu = uint_relel/uint_tot;
  if (relfracu > MAX_RELEL_FRAC_U) 
  {
    for (ie=0; ie<NRELBIN; ie++)
    {  
      pp[NEREL(ie)] *= (MAX_RELEL_FRAC_U/relfracu);
    }
  }
  
#endif //RELELECTRONS
#endif //EVOLVEELECTRONS

  //ANDREW TODO do we want this? Is this inconsistent with keeping entropy as a backup until the very end of time step? 
  //updates entropy after floor corrections
  if(ret<0)
    pp[5]=calc_Sfromu(pp[RHO],pp[UU],geom->ix,geom->iy,geom->iz);

  return ret;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//routines for various types of inversions
//**********************************************************************
//**********************************************************************
//**********************************************************************

//********************************************
//Harm u2p_hot
//********************************************

static FTYPE
dpdWp_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq, FTYPE gamma)
{
  FTYPE W=Wp+D;
  return( (gamma - 1.) * (1. - vsq) /  gamma ) ;
}

// 1 / (d(u+p)/dp)
static FTYPE
compute_idwmrho0dp(FTYPE wmrho0, FTYPE gamma)
{
  return((gamma-1.)/gamma);
}


// 1 / (drho0/dp) holding wmrho0 fixed
static FTYPE
compute_idrho0dp(FTYPE wmrho0)
{
  return(0.0);
}

static int
f_u2p_hot(ldouble Wp, ldouble* cons,ldouble *f,ldouble *df,ldouble *err,ldouble pgamma)
{

  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Qdotnp=cons[6];
  
  ldouble W=Wp+D;

  FTYPE W3,X3,Ssq,Wsq,X,X2,Xsq; 
  FTYPE Qtsq = Qt2;
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X2 = X*X;
  Xsq = X2;
  X3 = X2*X;

  ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
  ldouble gamma2 = 1./(1.-v2);
  ldouble gamma = sqrt(gamma2);
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = Wp/gamma2 - D*v2/(1.+gamma);
  ldouble u = wmrho0 / pgamma;
  ldouble p = (pgamma-1)*u;

  //original:
#if (U2P_EQS==U2P_EQS_NOBLE)
   *f = Qn + W - p + 0.5*Bsq*(1.+v2) - QdotBsq/2./Wsq;
 *err = fabs(*f) / (fabs(Qn) + fabs(W) + fabs(p) + fabs(0.5*Bsq*(1.+v2)) + fabs(QdotBsq/2./Wsq));
#endif

  //JONS:
#if (U2P_EQS==U2P_EQS_JON)
   *f = Qdotnp + Wp - p + 0.5*Bsq + (Bsq*Qtsq - QdotBsq)/X2;
 *err = fabs(*f) / (fabs(Qdotnp) + fabs(Wp) + fabs(p) + fabs(0.5*Bsq) + fabs((Bsq*Qtsq - QdotBsq)/X2));
#endif

  // dp/dWp = dp/dW + dP/dv^2 dv^2/dW  
  ldouble dvsq=(-2.0/X3 * ( Qtsq  +  QdotBsq * (3.0*W*X + Bsq*Bsq)/W3));
  ldouble dp1 = dpdWp_calc_vsq(Wp, D, v2 ,pgamma); // vsq can be unphysical

  ldouble idwmrho0dp=compute_idwmrho0dp(wmrho0,pgamma);
  ldouble dwmrho0dvsq = (D*(gamma*0.5-1.0) - Wp);

  ldouble drho0dvsq = -D*gamma*0.5; // because \rho = D/\gamma
  ldouble idrho0dp = compute_idrho0dp(wmrho0);

  ldouble dp2 =   drho0dvsq *idrho0dp  +   dwmrho0dvsq *idwmrho0dp;

  ldouble dpdW = dp1  + dp2*dvsq; // dp/dW = dp/dWp

  //original:
  #if (U2P_EQS==U2P_EQS_NOBLE)
  *df=1.-dpdW + QdotBsq/(Wsq*W) + 0.5*Bsq*dvsq;
  #endif

  //JONS:
  #if (U2P_EQS==U2P_EQS_JON)
  *df=1. -dpdW + (Bsq*Qtsq - QdotBsq)/X3*(-2.0);
  #endif

  return 0;  
}

//********************************************
//Harm u2p_entropy
//********************************************

// p(rho0, w-rho0 = u+p)
static FTYPE
pressure_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma)
{
  ldouble igammar = (gamma-1.)/gamma;
  return(igammar*wmrho0) ;
}

// local aux function
static FTYPE
compute_inside_entropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma)
{
  FTYPE pressure,indexn,insideentropy;

  pressure=pressure_wmrho0_idealgas(rho0,wmrho0,gamma);
  indexn=1.0/(gamma-1.);  
  insideentropy=pow(pressure,indexn)/pow(rho0,indexn+1.0);

  return(insideentropy);
}


// specific entropy as function of rho0 and internal energy (u)
// Ss(rho0,\chi=u+p)
// specific entropy = \ln( p^n/\rho^{n+1} )
static FTYPE
compute_specificentropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma)
{
  FTYPE insideentropy,specificentropy;

  insideentropy=compute_inside_entropy_wmrho0_idealgas(rho0, wmrho0,gamma);
  specificentropy=log(insideentropy);

  return(specificentropy);

}

// used for utoprim_jon when doing entropy evolution
// dSspecific/d\chi
static FTYPE
compute_dspecificSdwmrho0_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma)
{
  FTYPE dSdchi;

  dSdchi = 1.0/((gamma-1.)*wmrho0);
  // Again, GAMMA->1 means dSdchi->\infty unless \chi->0 or rho0->0

  return(dSdchi);

}

// dSspecific/drho0
static FTYPE
compute_dspecificSdrho_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0, FTYPE gamma)
{
  FTYPE dSdrho;
  
  dSdrho=gamma/((1.0-gamma)*rho0);

  return(dSdrho);
}

static int
f_u2p_entropy(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df, ldouble *err,ldouble pgamma)
{
  ldouble Qn=cons[0];
  ldouble Qt2=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Sc=cons[5];
 
  ldouble W=Wp+D;

  FTYPE W3,X3,Ssq,Wsq,X,X2,Xsq; 
  FTYPE Qtsq = Qt2;
  X = Bsq + W;
  Wsq = W*W;
  W3 = Wsq*W ;
  X2 = X*X;
  Xsq = X2;
  X3 = X2*X;

  ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
  ldouble gamma2 = 1./(1.-v2);
  ldouble gamma = sqrt(gamma2);
  ldouble rho0 = D/gamma;
  ldouble wmrho0 = Wp/gamma2 - D*v2/(1.+gamma);
  ldouble u = wmrho0 / pgamma;
  ldouble p = (pgamma-1)*u;

  ldouble Ssofchi=compute_specificentropy_wmrho0_idealgas(rho0,wmrho0,pgamma);

  *f= -Sc + D*Ssofchi;

  *err = fabs(*f) / (fabs(Sc) + fabs(D*Ssofchi));

  FTYPE dSsdW,dSsdvsq,dSsdWp,dScprimedWp,dSsdrho,dSsdchi;
  FTYPE dvsq,dwmrho0dW,drho0dW;
  FTYPE dwmrho0dvsq,drho0dvsq;

  dSsdrho=compute_dspecificSdrho_wmrho0_idealgas(rho0,wmrho0,pgamma);
  dSsdchi=compute_dspecificSdwmrho0_wmrho0_idealgas(rho0,wmrho0,pgamma);

  dwmrho0dW = 1.0/gamma; // holding utsq fixed
  drho0dW = 0.0; // because \rho=D/\gamma and holding utsq fixed
  dwmrho0dvsq = (D*(gamma*0.5-1.0) - Wp); // holding Wp fixed
  drho0dvsq = -D*gamma*0.5; // because \rho=D/\gamma and holding Wp fixed

  dvsq=(-2.0/X3 * ( Qtsq  +  QdotBsq * (3.0*W*X + Bsq*Bsq)/W3));

  dSsdW =   drho0dW   *dSsdrho +   dwmrho0dW   *dSsdchi; // dSs/dW' holding utsq fixed
  dSsdvsq = drho0dvsq *dSsdrho +   dwmrho0dvsq *dSsdchi;
  dSsdWp = dSsdW  + dSsdvsq*dvsq; // dSs/dW = dSs/dWp [total derivative]

  dScprimedWp = D*dSsdWp;

  *df = dScprimedWp;
  
  return 0;
 
}

//**********************************************************************
// solver wrapper
//**********************************************************************

int
u2p_solver(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{
  #ifdef NONRELMHD
  return u2p_solver_nonrel(uu,pp,ggg,Etype,verbose);
  #endif

  int (*solver)(ldouble*,ldouble*,void*,int,int);
  struct geometry *geom
    = (struct geometry *) ggg;
  
#if (U2P_SOLVER==U2P_SOLVER_WP)
  solver = & u2p_solver_Wp;
#endif
  
#if (U2P_SOLVER==U2P_SOLVER_W)  // this is the default
  solver = & u2p_solver_W;
#endif

  return (*solver)(uu,pp,ggg,Etype,verbose);
} 
 
//**********************************************************************
//non-relativistic, analytical, solver
//**********************************************************************

int
u2p_solver_nonrel(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{
  //prepare geometry
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5], (*GG)[5], gdet, gdetu, gdetu_inv;
  gg=geom->gg;  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  gdetu_inv = 1. / gdetu;

  //density
  ldouble rho=uu[RHO] * gdetu_inv;
  pp[RHO]=rho;

  //velocities u_i
  ldouble ucov[4],ucon[4],vcov[4];
  ucov[0]=-1.;
  ucov[1]=uu[VX] * gdetu_inv / rho;
  ucov[2]=uu[VY] * gdetu_inv / rho;
  ucov[3]=uu[VZ] * gdetu_inv / rho;

  indices_12(ucov,ucon,GG);

  ucon[0]=1.;

  ldouble v2=dot3nr(ucon,ucov);

  pp[VX]=ucon[1];
  pp[VY]=ucon[2];
  pp[VZ]=ucon[3];

  ldouble bsq=0.;

#ifdef MAGNFIELD
 
  ldouble bcon[4],bcov[4],Bcon[4];

  Bcon[0]=0.;
  Bcon[1]=uu[B1] * gdetu_inv;
  Bcon[2]=uu[B2] * gdetu_inv ;
  Bcon[3]=uu[B3] * gdetu_inv ;

  pp[B1]=Bcon[1];
  pp[B2]=Bcon[2];
  pp[B3]=Bcon[3];

  calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, geom, bcon, bcov, &bsq);

#endif

  
  ldouble uint;
  if(Etype==U2P_HOT)
    {
      uint = -uu[UU] * gdetu_inv - bsq/2. - rho*v2/2.;

      if(uint<NONRELMHDENTROPYCUT*rho || !isfinite(uint)) 
	{
	  //if(!isfinite(uint)) printf("%d %d > %e %e %e %e %e\n",geom->ix+TOI,geom->iy+TOI,uint,uu[UU]/gdetu,bsq/2.,rho*v2/2.,rho); 
	  //printf("%e %e\n",uint,rho);
	  //if(geom->ix>50) getch();
	  return -1;
	}

      pp[UU]=uint;
    }
  else if(Etype==U2P_ENTROPY)
    {
      ldouble S=uu[ENTR] * gdetu_inv ;
      uint= calc_ufromS(S,rho,geom->ix,geom->iy,geom->iz);
      pp[UU]=uint;
    }

  //pure entropy evolution - updated only in the end of RK2
  pp[ENTR]= uu[ENTR] * gdetu_inv;


  #ifdef EVOLVEELECTRONS
  ldouble Se=uu[ENTRE] * gdetu_inv ;
  pp[ENTRE]=Se;
 
  ldouble Si=uu[ENTRI] * gdetu_inv ;
  pp[ENTRI]=Si;
  #endif

  #ifdef RELELECTRONS
   int ib;
   for(ib=0;ib<NRELBIN;ib++)
     pp[NEREL(ib)]=uu[NEREL(ib)] * gdetu_inv ;

  #endif


  return 0;
} //u2p_solver_nonrel


//**********************************************************************
//Newton-Raphson solver 
//upgraded - uses Wp instead of W
//Etype == 0 -> hot inversion (uses D,Ttt,Tti)
//Etype == 1 -> entropy inversion (uses D,S,Tti)
//**********************************************************************

int
u2p_solver_Wp(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{
  int i,j,k;
  ldouble rho,uint,w,W,Wp,Wpprev,alpha,D,Sc,alphasq,betasqoalphasq;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],Qconp[4],Qcovp[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn,Qdotnp;
  ldouble QdotB,QdotBsq,Bcon[4],Bcov[4],Bsq;

  /****************************/
  //prepare geometry
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble pgamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  pgamma=pick_gammagas(geom->ix,geom->iy,geom->iz);
  #endif
  ldouble pgammam1=pgamma-1.;

  ldouble (*gg)[5], (*GG)[5], gdet, gdetu, gdetu_inv;
  gg=geom->gg;  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  gdetu_inv = 1. / gdetu;
                
  /****************************/
  //equations choice
  int (*f_u2p)(ldouble,ldouble*,ldouble*,ldouble*,ldouble*,ldouble);
 if(Etype==U2P_HOT) 
   f_u2p=&f_u2p_hot;
 if(Etype==U2P_ENTROPY) 
   f_u2p=&f_u2p_entropy;
  /****************************/
   
  if(verbose>1)
  {
    printf("********************\n");
    print_conserved(uu);
    print_primitives(pp);
  }

  //conserved quantities etc
  
  //alpha
  alpha=geom->alpha;
  alphasq=alpha*alpha;

  //D
  D=uu[0] * gdetu_inv * alpha; //uu[0]=gdetu rho ut

  //conserved entropy "gdet S u^t"
  Sc=uu[5] * gdetu_inv * alpha;

  //Qp_mu=alpha T^t_mu 
  Qcovp[0]=uu[1] * gdetu_inv * alpha;
  Qcovp[1]=uu[2] * gdetu_inv * alpha;
  Qcovp[2]=uu[3] * gdetu_inv * alpha;
  Qcovp[3]=uu[4] * gdetu_inv * alpha;

  //Qp^mu
  indices_12(Qcovp,Qconp,GG);

  //Q_mu=alpha (T^t_mu - rho u^t delta(t,mu)) - avoid this one
  Qcov[0]=(uu[1] * gdetu_inv -uu[0] * gdetu_inv ) * alpha;
  Qcov[1]=uu[2] * gdetu_inv * alpha;
  Qcov[2]=uu[3] * gdetu_inv * alpha;
  Qcov[3]=uu[4] * gdetu_inv * alpha;

  //Q^mu
  indices_12(Qcov,Qcon,GG);

#ifdef MAGNFIELD
  //curly B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1] * gdetu_inv * alpha;
  Bcon[2]=uu[B2] * gdetu_inv * alpha;
  Bcon[3]=uu[B3] * gdetu_inv * alpha;

  //B_mu
  indices_21(Bcon,Bcov,gg);

  Bsq = dot(Bcon,Bcov);

  QdotB = dot(Qcov,Bcon);

  QdotBsq = QdotB*QdotB;
  
#else
  Bsq=QdotB=QdotBsq=0.;
  Bcon[0]=Bcon[1]=Bcon[2]=Bcon[3]=0.;
#endif  

  //n_mu = (-alpha, 0, 0, 0)
  ncov[0]=-alpha;
  ncov[1]=ncov[2]=ncov[3]=0.;
  
  //n^mu
  indices_12(ncov,ncon,GG);

  //Q_mu n^mu = Q^mu n_mu = -alpha*Q^t
  Qn=Qcon[0] * ncov[0];

  //\beta^i \beta_i / \alpha^2 = g^{ti} g_{ti}
  betasqoalphasq=gg[0][1]*GG[0][1] + gg[0][2]*GG[0][2] + gg[0][3]*GG[0][3]; 
  
  //Qdotnp=-E'=-E+D
  ldouble Dfactor = (-geom->gttpert + alphasq*betasqoalphasq)/(alphasq+alpha);
  Qdotnp = Qconp[0]*ncov[0] + D*(Dfactor) ; // -Qdotn-W = -Qdotnp-Wp

  //j^mu_nu=delta^mu_nu +n^mu n_nu
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      jmunu[i][j] = delta(i,j) + ncon[i]*ncov[j];

  //Qtilda^nu = j^nu_mu Q^mu
  for(i=0;i<4;i++)
  {
    Qtcon[i]=0.;
    for(j=0;j<4;j++)
      Qtcon[i]+=jmunu[i][j]*Qcon[j];
  }

  //Qtilda_nu
  indices_21(Qtcon,Qtcov,gg);

  //Qt2=Qtilda^mu Qtilda_mu
  Qt2=dot(Qtcon,Qtcov);
  FTYPE Qtsq = Qt2;

  //initial guess for Wp = w gamma**2 based on primitives
  rho=pp[0];
  uint=pp[1];
  utcon[0]=0.;
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  
  if (VELPRIM != VELR)
  {
    conv_vels(utcon,utcon,VELPRIM,VELR,gg,GG);
  }

  ldouble qsq=0.;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      qsq+=utcon[i]*utcon[j]*gg[i][j];
  ldouble gamma2=1.+qsq;
  ldouble gamma=sqrt(gamma2);

  //Wp
  Wp=(pgamma*uint)*gamma2;

  ldouble Wpinit, Winit;
  if(verbose>1) printf("initial Wp:%e\n",Wp);
  Wpinit=Wp;
  Winit=Wpinit+D;

  //test if does not provide reasonable gamma2
  // Make sure that W is large enough so that v^2 < 1 and w-rho > 0 : 
  int i_increase = 0;
  ldouble f0,f1,dfdW,err;
  ldouble CONV=U2PCONV; 
  ldouble cons[7]={Qn,Qt2,D,QdotBsq,Bsq,Sc,Qdotnp};
  int iter=0,fu2pret;
  
  do
    {
      W=Wp+D;
      f0=dfdW=0.;

      FTYPE Wsq,Xsq,X; 
      X = Bsq + W;
      Xsq = X*X;
      Wsq = W*W;

      ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
      ldouble gamma2 = 1./(1.-v2);
      ldouble gamma = sqrt(gamma2);
      ldouble rho0 = D/gamma;
      ldouble wmrho0 = Wp/gamma2 - D*v2/(1.+gamma);

      (*f_u2p)(Wp,cons,&f0,&dfdW,&err,pgamma);
      
      if((gamma2<0. || Wp<0. || wmrho0<0.|| !isfinite(f0) || !isfinite(dfdW)) && (i_increase < 50))
	{
	  if(verbose>0) printf("init Wp : %e - %e %e %e %e\n",Wp,v2,wmrho0,f0,dfdW);
	  Wp *= 2.;
	  i_increase++;
	  continue;
	}
      else
	break;    
    }
  while(1);

  if(i_increase>=50)
    {
      if(verbose>0) 
	{printf("failed to find initial W for Etype: %d\n",Etype);
	  printf("at %d %d\n",geom->ix+TOI,geom->iy+TOJ);}
      return -150;

      print_NVvector(uu);
      print_NVvector(pp);
      //getchar();
    }

  //1d Newton solver 
  do
    {
      Wpprev=Wp;
      iter++;
     
      fu2pret=(*f_u2p)(Wp,cons,&f0,&dfdW,&err,pgamma);

      if(verbose>1) printf("%d %e %e %e %e\n",iter,Wp,f0,dfdW,err);
 
      //convergence test
      //if(err<CONV)
      //break;
      
      if(dfdW==0.)
      {
	Wp*=1.1;
	continue;
      }

      ldouble Wpnew=Wp-f0/dfdW;
      
      //don't jump over the zero
      Wpnew = my_max(Wpnew, Wp/100.);

      ldouble Wnew=Wpnew+D;
      int idump=0;
      ldouble dumpfac=1.;

      //test if goes out of bounds and damp step if so
      int itmaxdamp=50;
      do
	{
	  ldouble f0tmp,dfdWtmp,errtmp;
	  Wnew=Wpnew+D;
	  f0tmp=dfdWtmp=0.;

	  FTYPE Wsq,Xsq,X; 
	  X = Bsq + Wnew;
	  Xsq = X*X;
	  Wsq = Wnew*Wnew;

	  ldouble v2=( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*Wnew)) / (Wsq*Xsq);
	  ldouble gamma2 = 1./(1.-v2);
	  ldouble gamma = sqrt(gamma2);
	  ldouble rho0 = D/gamma;
	  ldouble wmrho0 = Wpnew/gamma2 - D*v2/(1.+gamma);

	  if(verbose>1) printf("sub (%d) :%d %e %e %e %e %e %e\n",idump,iter,Wpnew,f0tmp,dfdWtmp,v2,gamma2,wmrho0);

	  //if((gamma2<0. || Wpnew<0. || wmrho0<0. || !isfinite(f0tmp) || !isfinite(dfdWtmp)) && (idump<itmaxdamp))
	  if((gamma2<0. || Wpnew<0. || wmrho0<=0.) && (idump<itmaxdamp))
	    {
	      idump++;
	      dumpfac/=2.;
	      Wpnew=Wp-dumpfac*f0/dfdW;
	      continue;
	    }
	  else
	    break;
	}
      while(1);

      if(idump>=itmaxdamp) 
	{
	  if(verbose>0) printf("damped unsuccessfuly\n");
	  return -101;
	}

      if(fabs(W)>BIG) 
	{
	  if(verbose>1) printf("W has gone out of bounds at %d,%d,%d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz); 
	  return -103;
	}
	
      Wp=Wpnew; 

      //convergence test:
      if(err<CONV || (Etype==U2P_HOT && fabs((Wp-Wpprev)/Wpprev)<CONV && err<(sqrt(CONV)))
	 || (Etype==U2P_ENTROPY && fabs((Wp-Wpprev)/Wpprev)<CONV && err<0.99)) break;
      //if(err<CONV || fabs((Wp-Wpprev)/Winit)<CONV) break;
    }
  while(iter<50);
   
  if(iter>=50)
    {
      if(verbose>0) printf("iter exceeded in u2p_solver with Etype: %d\n",Etype); //getchar();
      return -102;
    }

  if(!isfinite(Wp) || !isfinite(Wp)) {if(verbose) printf("nan/inf W in u2p_solver with Etype: %d\n",Etype); return -103;}
 
  if(verbose>1) 
    {
      fu2pret=(*f_u2p)(Wp,cons,&f0,&dfdW,&err,pgamma);
      printf("end: %d %e %e %e %e\n",iter,Wp,f0,dfdW,err);
    }

  //W found, let's calculate v2 and the rest
  W=Wp+D;
  ldouble Wsq,Xsq,v2,wmrho0,entr;
	
  Wsq = W*W ;
  Xsq = (Bsq + W) * (Bsq + W);  
  v2 = ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);

  gamma2=1./(1.-v2);
  gamma=sqrt(gamma2);
  rho=D/gamma;
  entr=Sc/gamma;
  // w-\rho_0 = (u+p) = W'/\gamma^2 - D v^2/(1+\gamma)
  wmrho0 = Wp/gamma2 - D*v2/(1.+gamma);
  uint=1./pgamma*wmrho0;
  utcon[0]=0.;
  utcon[1]=gamma/(W+Bsq)*(Qtcon[1]+QdotB*Bcon[1]/W);
  utcon[2]=gamma/(W+Bsq)*(Qtcon[2]+QdotB*Bcon[2]/W);
  utcon[3]=gamma/(W+Bsq)*(Qtcon[3]+QdotB*Bcon[3]/W);

  if(!isfinite(utcon[1]))
    {
      //print_4vector(utcon);
      return -120;
    }


  if(uint<0. || gamma2<0. ||isnan(Wp) || !isfinite(Wp)) 
    {
      if(verbose>0) printf("neg u in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);//getchar();
      return -104;
    }

  //converting to VELPRIM
  if (VELR != VELPRIM)
  {
    conv_vels(utcon,utcon,VELR,VELPRIM,gg,GG);
  }
  
  //returning new primitives
  pp[RHO]=rho;
  pp[UU]=uint;
  pp[VX]=utcon[1];
  pp[VY]=utcon[2];
  pp[VZ]=utcon[3];

  if(rho<0.) 
    {
      if(verbose>0) printf("neg rho in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);//getchar();
      return -105;
    }

  //pure entropy evolution - updated only in the end of RK2
  pp[ENTR]=entr;

#ifdef MAGNFIELD
  //magnetic conserved=primitives
  pp[B1]=uu[B1] * gdetu_inv;
  pp[B2]=uu[B2] * gdetu_inv;
  pp[B3]=uu[B3] * gdetu_inv;
#endif

#ifdef EVOLVEELECTRONS
  conv_vels(utcon,utcon,VELPRIM,VEL4,gg,GG);
  ldouble Se=uu[ENTRE] * gdetu_inv / utcon[0];
  pp[ENTRE]=Se;
  ldouble Si=uu[ENTRI] * gdetu_inv / utcon[0];
  pp[ENTRI]=Si;

#ifdef RELELECTRONS
  int ib;
  for(ib=0;ib<NRELBIN;ib++)
    pp[NEREL(ib)]=uu[NEREL(ib)] * gdetu_inv / utcon[0];
#endif
#endif
  
  if(verbose) print_primitives(pp);

  //verify uunew against uuoriginal
  //seems unnecessary - if this point reached, it already did its best and its best to stay like that
  if(Etype==U2P_HOT)
    {
      ldouble uunew[NV];
      p2u(pp,uunew,geom);
      
      ldouble errinv,maxerrinv=-1.;
      int iv;
      
      //do we recover rho properly
      iv=RHO;
      errinv = fabs((uunew[iv]-uu[iv])/uu[iv]);
      if(errinv > maxerrinv) maxerrinv=errinv;
      //internal energy
      if(Etype==U2P_HOT)
	{
	  iv=UU;
	  errinv = fabs((uunew[iv]-uu[iv])/uu[iv]);
	  if(errinv > maxerrinv) maxerrinv=errinv;
	}
      if(Etype==U2P_ENTROPY)
	{
	  iv=ENTR;
	  errinv = fabs((uunew[iv]-uu[iv])/uu[iv]);
	  if(errinv > maxerrinv) maxerrinv=errinv;
	}
      
      double inverr=1.e-2;
      
      if(Etype==U2P_ENTROPY) inverr=0.999;
      
      if(Etype==U2P_HOT) inverr=0.1;
      
      if(maxerrinv>inverr)// && verbose>0) 
    {
      
      if(Etype==U2P_ENTROPY && 0) { 
	printf("verify u2p (%d) failed: %e || ",Etype,maxerrinv);
	printf("%e %e | %e %e | %e %e\n",uunew[RHO],uu[RHO],uunew[ENTR],uu[ENTR],uunew[VX],uu[VX]);
	print_conserved(uu);
	print_conserved(uunew);
	print_primitives(pp);
	getch();
      }
      if(Etype==U2P_HOT && 0) { 
	printf("verify u2p (%d) failed: %e || ",Etype,maxerrinv);
	printf("%e %e | %e %e | %e %e\n",uunew[RHO],uu[RHO],uunew[UU],uu[UU],uunew[VX],uu[VX]);
	print_conserved(uu);
	print_conserved(uunew);
	print_primitives(pp);
	getch();
      }
      return -200;
    }
    }

  if(verbose>0) printf("u2p_solver returns 0\n");
  return 0; //ok

}

//**********************************************************************
//old Newton-Raphson solver 
//iterates W, not Wp
//Etype == 0 -> hot inversion (uses D,Ttt,Tti)
//Etype == 1 -> entropy inversion (uses D,S,Tti)
//**********************************************************************

int
u2p_solver_W(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{
  int i,j,k;
  ldouble rho,uint,w,W,alpha,D,Sc;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],Qconp[4],Qcovp[4],jmunu[4][4],Qtcon[4],Qtcov[4],Qt2,Qn;
  ldouble QdotB,QdotBsq,Bcon[4],Bcov[4],Bsq;
  
  /*******************************************************/
  //prepare geometry
  struct geometry *geom
  = (struct geometry *) ggg;
  
  ldouble pgamma=GAMMA;
#ifdef CONSISTENTGAMMA
  pgamma=pick_gammagas(geom->ix,geom->iy,geom->iz);
#endif
  ldouble pgammam1=pgamma-1.;
  
  ldouble (*gg)[5], (*GG)[5], gdet, gdetu, gdetu_inv;
  gg=geom->gg; GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  gdetu_inv = 1. / gdetu;
  
  /****************************/
  //equations choice
  int (*f_u2p)(ldouble,ldouble*,ldouble*,ldouble*,ldouble*,ldouble);
  if(Etype==U2P_HOT)
    f_u2p=&f_u2p_hot;
  if(Etype==U2P_ENTROPY)
    f_u2p=&f_u2p_entropy;
  /****************************/
  
  if(verbose>1)
  {
    printf("********************\n");
    print_conserved(uu);
    print_primitives(pp);
  }
  
  //conserved quantities etc
  
  //alpha
  //alpha=sqrt(-1./GG[0][0]);
  alpha = geom->alpha;
  
  //D
  D = uu[0] * gdetu_inv * alpha; //uu[0]=gdetu rho ut
  
  //conserved entropy "S u^t"
  Sc = uu[5] * gdetu_inv * alpha;
  
  //Q_mu=alpha T^t_mu
  Qcov[0] = (uu[1] * gdetu_inv - uu[0] * gdetu_inv) * alpha;
  Qcov[1] = uu[2] * gdetu_inv * alpha;
  Qcov[2] = uu[3] * gdetu_inv * alpha;
  Qcov[3] = uu[4] * gdetu_inv * alpha;
  
  //Qp_mu=alpha T^t_mu
  Qcovp[0] = uu[1] * gdetu_inv *alpha;
  Qcovp[1] = uu[2] * gdetu_inv *alpha;
  Qcovp[2] = uu[3] * gdetu_inv *alpha;
  Qcovp[3] = uu[4] * gdetu_inv *alpha;
  
  //Qp^mu
  indices_12(Qcovp,Qconp,GG);
  
  //Q^mu
  indices_12(Qcov,Qcon,GG);
  
#ifdef MAGNFIELD
  //curly B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1] * gdetu_inv * alpha;
  Bcon[2]=uu[B2] * gdetu_inv * alpha;
  Bcon[3]=uu[B3] * gdetu_inv * alpha;
  
  //B_mu
  indices_21(Bcon,Bcov,gg);
  
  Bsq = dot(Bcon,Bcov);
  QdotB = dot(Qcov,Bcon);
  QdotBsq = QdotB*QdotB;

#else
  Bsq=QdotB=QdotBsq=0.;
  Bcon[0]=Bcon[1]=Bcon[2]=Bcon[3]=0.;
#endif  // MAGNFIELD

  //normal observer velocity
  //n_mu = (-alpha, 0, 0, 0)
  ncov[0]=-alpha;
  ncov[1]=ncov[2]=ncov[3]=0.;
  
  //n^mu
  indices_12(ncov,ncon,GG);
  
  //Q_mu n^mu = Q^mu n_mu = -alpha*Q^t
  Qn = Qcon[0] * ncov[0];
  
  //j^mu_nu = delta^mu_nu +n^mu n_nu
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      jmunu[i][j] = delta(i,j) + ncon[i]*ncov[j];
  
  //Qtilda^nu = j^nu_mu Q^mu
  for(i=0;i<4;i++)
  {
    Qtcon[i]=0.;
#ifdef APPLY_OMP_SIMD
    //#pragma omp simd
#endif
    for(j=0;j<4;j++)
      Qtcon[i]+=jmunu[i][j]*Qcon[j];
  }
  
  //Qtilda_nu
  indices_21(Qtcon,Qtcov,gg);
  
  //Qt2=Qtilda^mu Qtilda_mu
  Qt2=dot(Qtcon,Qtcov);
  FTYPE Qtsq = Qt2;
  
  //\beta^i \beta_i / \alpha^2 = g^{ti} g_{ti}
  ldouble betasqoalphasq=gg[0][1]*GG[0][1] + gg[0][2]*GG[0][2] + gg[0][3]*GG[0][3];
  ldouble alphasq=alpha*alpha;
  //Qdotnp=-E'=-E+D
  ldouble Dfactor = (-geom->gttpert + alphasq*betasqoalphasq)/(alphasq+alpha);
  ldouble Qdotnp = Qconp[0]*ncov[0] + D*(Dfactor) ; // -Qdotn - W = -Qdotnp-Wp
  
  //initial guess for W = w gamma**2 based on current primitives
  rho=pp[0];
  uint=pp[1];
  utcon[0]=0.;
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  
  if (VELPRIM != VELR)
  {
    conv_vels(utcon,utcon,VELPRIM,VELR,gg,GG);
  }
  
  ldouble qsq=0.;
  for(i=1;i<4;i++)
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
    for(j=1;j<4;j++)
      qsq+=utcon[i]*utcon[j]*gg[i][j];
  ldouble gamma2=1.+qsq;
  ldouble gamma=sqrt(gamma2);
  
  //W
  W=(rho+pgamma*uint)*gamma2;
  
  if(verbose>1) printf("initial W:%e\n",W);
  
  // test if does not provide reasonable gamma2
  // Make sure that W is large enough so that v^2 < 1 :
  int i_increase = 0;
  ldouble f0,f1,dfdW,err;
  ldouble CONV=U2PCONV;
  ldouble EPS=1.e-4;
  ldouble Wprev=W;
  ldouble cons[7]={Qn,Qt2,D,QdotBsq,Bsq,Sc,Qdotnp};
  
  do
  {
    f0=dfdW=0.;
    
    //if(Etype!=U2P_HOT) //entropy-like solvers require this additional check
    //now invoked for all solvers:
    (*f_u2p)(W-D,cons,&f0,&dfdW,&err,pgamma);
    
    if( ((( W*W*W * ( W + 2.*Bsq )
           - QdotBsq*(2.*W + Bsq) ) <= W*W*(Qtsq-Bsq*Bsq))
         || !isfinite(f0) || !isfinite(f0)
         || !isfinite(dfdW) || !isfinite(dfdW))
       && (i_increase < 50))
    {
      if(verbose>0) printf("init W : %e -> %e (%e %e)\n",W,100.*W,f0,dfdW);
      W *= 10.;
      i_increase++;
      continue;
    }
    else
      break;
  }
  while(1);
  
  if(i_increase>=50)
  {
    return -150;
    printf("failed to find initial W for Etype: %d\n",Etype);
    printf("at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
    print_NVvector(uu);
    print_NVvector(pp);
    //getchar();
  }
  
  //1d Newton solver
  int iter=0, fu2pret;
  do
  {
    Wprev=W;
    iter++;

    fu2pret=(*f_u2p)(W-D,cons,&f0,&dfdW,&err,pgamma);
    
    if(verbose>1) printf("%d %e %e %e %e\n",iter,W,f0,dfdW,err);
    
    //convergence test
    if(err<CONV)
      break;
    
    if(dfdW==0.) {W*=1.1; continue;}
    
    ldouble Wnew = W-f0/dfdW;
    int idump=0;
    ldouble dumpfac=1.;
    
    //test if goes out of bounds and damp solution if so
    do
    {
      ldouble f0tmp,dfdWtmp,errtmp;
      f0tmp=dfdWtmp=0.;
      //now for all solvers
      //if(Etype!=U2P_HOT) //entropy-like solvers require this additional check
      (*f_u2p)(Wnew-D,cons,&f0tmp,&dfdWtmp,&errtmp,pgamma);
      if(verbose>1) printf("sub (%d) :%d %e %e %e %e\n",idump,iter,Wnew,f0tmp,dfdWtmp,errtmp);
      if( ((( Wnew*Wnew*Wnew * ( Wnew + 2.*Bsq )
             - QdotBsq*(2.*Wnew + Bsq) ) <= Wnew*Wnew*(Qtsq-Bsq*Bsq))
           || !isfinite(f0tmp) || !isfinite(f0tmp)
           || !isfinite(dfdWtmp) || !isfinite(dfdWtmp))
         && (idump<100))
      {
        idump++;
        dumpfac/=2.;
        Wnew=W-dumpfac*f0/dfdW;
        continue;
      }
      else
        break;
    }
    while(1);
    
    if(idump>=100)
    {
      if(verbose>0) printf("damped unsuccessfuly\n");
      return -101;
    }
    
    W=Wnew;
    
    if(fabs(W)>BIG)
    {
      if(verbose>1) printf("W has gone out of bounds at %d,%d,%d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz);
      return -103;
    }
    
    if(fabs((W-Wprev)/Wprev)<CONV && err<1.e-1) break;
    //if(fabs((W-Wprev)/Wprev)<CONV && err<sqrt(CONV)) break;
  }
  while(iter<50);
  
  if(iter>=50)
  {
    if(verbose>0) printf("iter exceeded in u2p_solver with Etype: %d\n",Etype); //getchar();
    return -102;
  }
  
  if(!isfinite(W))
  {
    if(verbose) printf("nan/inf W in u2p_solver with Etype: %d\n",Etype);
    return -103;
  }
  
  if(verbose>1)
  {
    fu2pret=(*f_u2p)(W-D,cons,&f0,&dfdW,&err,pgamma);
    printf("end: %d %e %e %e %e\n",iter,W,f0,dfdW,err);
  }
  
  //W found, let's calculate v2 and the rest
  ldouble Wsq,Xsq,v2,entr;
  
  Wsq = W*W ;
  Xsq = (Bsq + W) * (Bsq + W);
  v2 = ( Wsq * Qtsq  + QdotBsq * (Bsq + 2.*W)) / (Wsq*Xsq);
  
  gamma2=1./(1.-v2);
  gamma=sqrt(gamma2);
  rho=D/gamma;
  entr=Sc/gamma;
  uint=1./pgamma*(W/gamma2-rho);
  utcon[0]=0.;
  utcon[1]=gamma/(W+Bsq)*(Qtcon[1]+QdotB*Bcon[1]/W);
  utcon[2]=gamma/(W+Bsq)*(Qtcon[2]+QdotB*Bcon[2]/W);
  utcon[3]=gamma/(W+Bsq)*(Qtcon[3]+QdotB*Bcon[3]/W);

  if(verbose>1)
    printf("end2: %e %e %e %e %e %e\n",W,D,pgamma,gamma2,rho,uint);  
  if(!isfinite(utcon[1]))
  {
    return -120;
  }
  
  if(uint<0. || gamma2<0. || isnan(W) || !isfinite(W))
  {
    if(verbose>0) printf("neg u in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);//getchar();
    return -104;
  }
  
  //converting to VELPRIM
  if (VELR != VELPRIM)
  {
    conv_vels(utcon,utcon,VELR,VELPRIM,gg,GG);
  }
  
  //returning new primitives
  pp[RHO]=rho;
  pp[UU]=uint;
  pp[VX]=utcon[1];
  pp[VY]=utcon[2];
  pp[VZ]=utcon[3];
  
  if(rho<0.)
  {
    if(verbose>0) printf("neg rho in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);
    //getchar();
    return -105;
  }
  
  //entropy based on Etype  
  //pure entropy evolution - updated only in the end of RK2
  pp[ENTR]=entr;
  
#ifdef MAGNFIELD
  //magnetic conserved=primitives
  pp[B1]=uu[B1] * gdetu_inv;
  pp[B2]=uu[B2] * gdetu_inv;
  pp[B3]=uu[B3] * gdetu_inv;
#endif
  
#ifdef EVOLVEELECTRONS
  conv_vels(utcon,utcon,VELPRIM,VEL4,gg,GG);
  ldouble Se=uu[ENTRE] * gdetu_inv / utcon[0];
  pp[ENTRE]=Se;
  ldouble Si=uu[ENTRI] * gdetu_inv / utcon[0];
  pp[ENTRI]=Si;
  
#ifdef RELELECTRONS
  int ib;
  for(ib=0;ib<NRELBIN;ib++)
    pp[NEREL(ib)]=uu[NEREL(ib)] * gdetu_inv / utcon[0];
#endif
#endif
  
  if(verbose) print_primitives(pp);
  
  if(verbose>0) printf("u2p_solver returns 0\n");
  return 0; //ok
}


//**********************************************************************
//recovers only magnetic field primitives - used when correcting polar axis
//**********************************************************************

int
u2p_solver_Bonly(ldouble *uu, ldouble *pp, void *ggg)
{
  //prepare geometry
  struct geometry *geom
  = (struct geometry *) ggg;
  
  ldouble gdet, gdetu, gdetu_inv;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  gdetu_inv = 1. / gdetu;
  
#ifdef MAGNFIELD
  //magnetic conserved=primitives
  pp[B1]=uu[B1] * gdetu_inv;
  pp[B2]=uu[B2] * gdetu_inv;
  pp[B3]=uu[B3] * gdetu_inv;
#endif
  
  return 0; //always succeeds
}  // int u2p_solver_Bonly


//**********************************************************************
//count the number of entropy inversions
//**********************************************************************

int count_entropy(int *n, int *n2)
{
  int nentr=0,nentrloc=0,ii,ix,iy,iz;
  int nentr2=0,nentrloc2=0;

  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  nentrloc+=get_cflag(ENTROPYFLAG,ix,iy,iz); 
	  nentrloc2+=get_cflag(ENTROPYFLAG2,ix,iy,iz); 
	}

#ifdef MPI
  MPI_Allreduce(&nentrloc, &nentr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  
  MPI_Allreduce(&nentrloc2, &nentr2, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  
#else
  nentr=nentrloc;
  nentr2=nentrloc2;
#endif
    
  *n = nentr;
  *n2 = nentr2;
  return 0;
}

//**********************************************************************
//backups entropy count to spit it into a silo file
//**********************************************************************

int copy_entropycount()
{
  int ii,ix,iy,iz;
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2];
      ldouble val=get_cflag(ENTROPYFLAG,ix,iy,iz);
      set_cflag(ENTROPYFLAG3,ix,iy,iz,val);
    }

  return 0;
}

//**********************************************************************
//calculates entropy corresponding to given rho and uint
//**********************************************************************

int
update_entropy()
{
  int ii;
  //#pragma omp parallel for private(ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      int ix,iy,iz;
      ldouble rho,uint,entr;
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);
     
      rho=get_u(p,RHO,ix,iy,iz);
      uint=get_u(p,UU,ix,iy,iz);
      entr=calc_Sfromu(rho,uint,ix,iy,iz);
     
      //printf("%d %d > %e %e > %e\n",ix,iy,rho,uint,entr); getch();

      set_u(p,ENTR,ix,iy,iz,entr);

      p2u_mhd(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);

    }
  return 0;

}


//**********************************************************************
//tests
//**********************************************************************

int
test_inversion()
{
  ldouble pp[NV],pp2[NV],uu[NV],ucon[4]={0.,0.,0.,0.};
  struct geometry geom,geomBL;
  int iv;

  fill_geometry(NX-2,0,0,&geom);
  fill_geometry_arb(NX-2,0,0,&geomBL,BLCOORDS);
  ucon[1]=1.e-12;
  conv_vels(ucon,ucon,VEL4,VEL4,geomBL.gg,geomBL.GG);
  trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
  conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
  pp[RHO]=1.;
  pp[UU]=1.;
  pp[VX]=ucon[1];
  pp[VY]=pp[VZ]=0.;

  print_primitives(pp);
  p2u(pp,uu,&geom);

  pp[VX]*=100000.*M_PI;
  pp[RHO]*=0.001245325124;
  pp[UU]*=23.124124214421124;
  
  PLOOP(iv) pp2[iv]=pp[iv];
  

  
  print_conserved(uu);
  printf("gdet = %e\n",geom.gdet);
  u2p_solver(uu,pp,&geom,U2P_HOT,0); 
  print_primitives(pp);

  return 0;

}

int
test_inversion_nonrel()
{
  ldouble pp[NV],pp2[NV],uu[NV],ucon[4]={0.,0.,0.,0.};
  struct geometry geom,geomBL;
  int iv;

  fill_geometry(NX-2,NY/2,0,&geom);

  print_metric(geom.gg);

  ucon[1]=0.1;
  ucon[2]=0.001;
  ucon[3]=0.001;
  conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);

  pp[RHO]=100.;
  pp[UU]=1.e-5;
  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],geom.ix,geom.iy,geom.iz);
  pp[VX]=ucon[1];
  pp[VY]=ucon[2];
  pp[VZ]=ucon[3];

#ifdef MAGNFIELD
  pp[B1]=pp[B2]=pp[B3]=0.;
  pp[B1]=1.e-1;
  pp[B2]=1.e-4;
  pp[B3]=1.e-3;
#endif

#ifdef RADIATION
  pp[EE]=pp[UU];

  pp[FX]=pp[FY]=pp[FZ]=0.;
#endif


  print_primitives(pp);
  p2u(pp,uu,&geom);

  print_conserved(uu);
  printf("gdet = %e\n",geom.gdet);

  u2p_solver(uu,pp,&geom,U2P_HOT,0); 
  print_primitives(pp);

  return 0;

}
