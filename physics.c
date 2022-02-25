/*! \file physics.c
 \brief Problem independent physics
*/

#include "ko.h"

//********************************************************************
//calculates the state corresponding to the given primitives
//**********************************************************************

int
fill_struct_of_state(ldouble *pp, void* ggg, void* sss)
{
  struct geometry *geom
  = (struct geometry *) ggg;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;
  
  int i,j;
  ldouble rho, uint, gamma, pgas, Tgas, Te, Ti, Sgas, STe, STi;
  ldouble ucon[4], ucov[4], bcon[4], bcov[4], Bcon[4], bsq;
  ldouble urfcon[4], urfcov[4], relgamma, Rij[4][4], Ehat, Trad, TradBB, kappaes;
  struct opacities opac;
  
  rho = pp[RHO];
  uint = pp[UU];
  
  state->rho = rho;
  state->uint = uint;
  gamma = GAMMA;
#ifndef FORCEGAMMAGASFIXED
#ifdef CONSISTENTGAMMA
  gamma = pick_gammagas(geom->ix,geom->iy,geom->iz);
  //ANDREW -- can we do this consistently instead?
  //gamma = 1. + (state->pi+state->pe+state->penth)/(state->ue+state->ui+state->uenth);
#endif
#endif
  pgas = (gamma - 1.) * uint;
  state->gamma = gamma;
  state->pgas = pgas;

  //ANDREW this is in different units than in pp[ENTR]!
  Sgas = kB_over_mugas_mp*calc_Sfromu(rho, uint, geom->ix,geom->iy,geom->iz);

  Tgas = calc_PEQ_Teifrompp(pp,&Te,&Ti,geom->ix,geom->iy,geom->iz);
  
  state->Tgas = Tgas;
  state->Sgas = Sgas;

  ldouble netot=one_over_mue_mp*rho;  
  ldouble nnth=0.;
  ldouble pnth=0.;
  ldouble unth=0.;
#ifdef RELELECTRONS
  nnth = calc_relel_ne(pp);
  //if(nnth/netot > MAX_RELEL_FRAC_N) nnth = MAX_RELEL_FRAC_N*netot; 
 
  unth = calc_relel_uint(pp);
  //if (unth/uint > MAX_RELEL_FRAC_U) unth = MAX_RELEL_FRAC_U*uint;

  pnth = calc_relel_p(pp);
  //if (pnth/pgas > MAX_RELEL_FRAC_P) pnth = MAX_RELEL_FRAC_P*pgas;
#endif
  state->nenth=nnth;
  state->uenth=unth;
  state->penth=pnth;

  state->ne = netot - nnth;
  state->ni = one_over_mui_mp*rho;

  state->Te = Te;
  state->Ti = Ti;

  state->pe = (state->ne)*K_BOLTZ*Te;
  state->pi = (state->ni)*K_BOLTZ*Ti;

  ldouble gammae=GAMMA;
  ldouble gammai=GAMMA;
  #ifdef GAMMAE
  gammae=GAMMAE;
  #endif
  #ifdef GAMMAI
  gammai=GAMMAI;
  #endif
  #if defined(CONSISTENTGAMMA) && !defined(FIXEDGAMMASPECIES)
  gammae=calc_gammaintfromtemp(Te,ELECTRONS);
  gammai=calc_gammaintfromtemp(Ti,IONS);
  #endif
  state->gammae=gammae;
  state->gammai=gammai;

  state->ue=(state->pe)/(gammae-1.);
  state->ui=(state->pi)/(gammai-1.);
  
  ldouble rhoeth = (state->ne)*MU_E*M_PROTON;
  STe = calc_SefromrhoT(rhoeth, Te, ELECTRONS);
  STi = calc_SefromrhoT(rho, Ti, IONS);
  
  state->STe = STe;
  state->STi = STi;

  calc_ucon_ucov_from_prims(pp, geom, ucon, ucov);
  
  DLOOPA(i)
  {
    state->ucon[i] = ucon[i];
    state->ucov[i] = ucov[i];
  }

#ifdef MAGNFIELD
  
  calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, geom, bcon, bcov, &bsq);
  calc_Bcon_4vel(pp, ucon, bcon, Bcon);
  
  DLOOPA(i)
  {
    state->bcon[i] = bcon[i];
    state->bcov[i] = bcov[i];
    state->Bcon[i] = Bcon[i];
  }
  
  state->bsq = bsq;
  
#else
  
  state->bsq = 0.;
  
#endif
  
#ifdef RADIATION
  
  calc_urcon_urcov_from_prims(pp, geom, urfcon, urfcov);
  calc_Rij_M1_from_4vel(pp, geom, urfcon, Rij);
  calc_Ehat_from_Rij_ucov(Rij, ucov, &Ehat);
  
  DLOOPA(i)
  {
    DLOOPA(j)
    {
      state->Rij[i][j] = Rij[i][j];
    }
    
    state->urfcon[i] = urfcon[i];
    state->urfcov[i] = urfcov[i];
  }
  
  state->Ehat = Ehat;
  
  TradBB = calc_LTE_TfromE(Ehat);
  state->TradBB = TradBB;

#ifdef EVOLVEPHOTONNUMBER 
  relgamma = urfcon[0]*ucov[0] + urfcon[1]*ucov[1] +urfcon[2]*ucov[2] +urfcon[3]*ucov[3];
  state->relgamma = relgamma;
  ldouble nphhat = -relgamma * pp[NF];
  Trad = calc_ncompt_Thatrad_nphhat(nphhat, Ehat);
  
#ifdef MAXDIFFTRADS
  ldouble maxfac=MAXDIFFTRADS;
  if (Trad > maxfac * TradBB)
  {
    pp[NF] *= (Trad / (maxfac * TradBB)); 
    Trad = maxfac * TradBB;
  }
#ifndef SYNCHROTRON
  if (Trad < TradBB)
  {
    pp[NF] *= (Trad / TradBB);    
    Trad = TradBB;
  }
#else
  if (Trad < TradBB / maxfac)
  {
    pp[NF] *= (Trad * maxfac / TradBB);    
    Trad = TradBB / maxfac;
  }
#endif //SYNCHROTRON
#endif //MAXDIFFTRADS
#else //EVOLVEPHOTONNUMBER
  relgamma = -1.;
  Trad = TradBB;
#endif
  
  state->relgamma = relgamma;
  state->Trad = Trad;
  
  kappaes = calc_kappaes_with_temperatures(rho, Tgas, Te, Ti, Trad);
  state->kappaes = kappaes;
  
  ldouble kappaGasAbs, kappaRadAbs, kappaGasNum, kappaRadNum, kappaGasRoss, kappaRadRoss;
  ldouble kappa = calc_kappa_from_state(pp, state, geom, &opac);
  
  kappaGasAbs = opac.kappaGasAbs;
  kappaRadAbs = opac.kappaRadAbs;
  kappaGasNum = opac.kappaGasNum;
  kappaRadNum = opac.kappaRadNum;
  kappaGasRoss = opac.kappaGasRoss;
  kappaRadRoss = opac.kappaRadRoss;
  
  state->opac.kappaGasAbs = kappaGasAbs;
  state->opac.kappaRadAbs = kappaRadAbs;
  state->opac.kappaGasNum = kappaGasNum;
  state->opac.kappaRadNum = kappaRadNum;
  state->opac.kappaGasRoss = kappaGasRoss;
  state->opac.kappaRadRoss = kappaRadRoss;

#endif  // ifdef RADIATION
  
  return 0;
}


//**********************************************************************
//Copies one struct_of_state to another
//**********************************************************************

int copy_state_to_state(void *sss1, void *sss2)
{
  int i, j;
  
  struct struct_of_state *state1
  = (struct struct_of_state *) sss1;
  struct struct_of_state *state2
  = (struct struct_of_state *) sss2;
  
  state2->rho = state1->rho;
  state2->uint = state1->uint;
  state2->gamma = state1->gamma;
  state2->pgas = state1->pgas;
  
  state2->Tgas = state1->Tgas;
  state2->Te = state1->Te;
  state2->Ti = state1->Ti;
  
  state2->Sgas = state1->Sgas;
  state2->STe = state1->STe;
  state2->STi = state1->STi;

  state2->ne = state1->ne;
  state2->pe = state1->pe;
  state2->ue = state1->ue;
  state2->ni = state1->ni;
  state2->ui = state1->ui;
  state2->pi = state1->pi;

  state2->nenth = state1->nenth;
  state2->uenth = state1->uenth;
  state2->penth = state1->penth;
  
  for (i = 0; i < 4; i++)
  {
    state2->ucon[i] = state1->ucon[i];
    state2->ucov[i] = state1->ucov[i];
  }
  
#ifdef MAGNFIELD
  
  for (i = 0; i < 4; i++)
  {
    state2->bcon[i] = state1->bcon[i];
    state2->bcov[i] = state1->bcov[i];
    state2->Bcon[i] = state1->Bcon[i];
  }
  
  state2->bsq = state1->bsq;
  
#endif
  
#ifdef RADIATION
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      state2->Rij[i][j] = state1->Rij[i][j];
    }
    
    state2->urfcon[i] = state1->urfcon[i];
    state2->urfcov[i] = state1->urfcov[i];
  }
  
  state2->Ehat = state1->Ehat;
  state2->TradBB = state1->TradBB;
  state2->relgamma = state1->relgamma;
  state2->Trad = state1->Trad;
  state2->kappaes = state1->kappaes;
  
  state2->opac.kappaGasAbs = state1->opac.kappaGasAbs;
  state2->opac.kappaRadAbs = state1->opac.kappaRadAbs;
  state2->opac.kappaGasNum = state1->opac.kappaGasNum;
  state2->opac.kappaRadNum = state1->opac.kappaRadNum;
  state2->opac.kappaGasRoss = state1->opac.kappaGasRoss;
  state2->opac.kappaRadRoss = state1->opac.kappaRadRoss;
#endif  // ifdef RADIATION
  
  return 0;
}

//**********************************************************************
// Updates the state for a new value of nphoton
//**********************************************************************

int update_state_for_nphoton(ldouble *pp, void *ggg, void *sss)
{
  struct geometry *geom
  = (struct geometry *) ggg;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;

  // Trad changes when nphoton is changed and this also affects the opacities  
#ifdef EVOLVEPHOTONNUMBER
  
  ldouble nphhat = -state->relgamma * pp[NF];
  state->Trad = calc_ncompt_Thatrad_nphhat(nphhat, state->Ehat);
  
  ldouble kappaes = calc_kappaes_with_temperatures(state->rho, state->Tgas, state->Te, state->Ti, state->Trad);
  state->kappaes = kappaes;

  struct opacities opac;
  ldouble kappaGasAbs, kappaRadAbs, kappaGasNum, kappaRadNum, kappaGasRoss, kappaRadRoss;
  ldouble kappa=calc_kappa_from_state(pp, state, geom, &opac);
  
  kappaGasAbs = opac.kappaGasAbs;
  kappaRadAbs = opac.kappaRadAbs;
  kappaGasNum = opac.kappaGasNum;
  kappaRadNum = opac.kappaRadNum;
  kappaGasRoss = opac.kappaGasRoss;
  kappaRadRoss = opac.kappaRadRoss;
  
  state->opac.kappaGasAbs = kappaGasAbs;
  state->opac.kappaRadAbs = kappaRadAbs;
  state->opac.kappaGasNum = kappaGasNum;
  state->opac.kappaRadNum = kappaRadNum;
  state->opac.kappaGasRoss = kappaGasRoss;
  state->opac.kappaRadRoss = kappaRadRoss;
  
#endif  //  no change otherwise

  return 0;
}

//**********************************************************************
//calculates fluid quantities for thermal electrons
//**********************************************************************


//Total number density of thermal electrons  
ldouble calc_thermal_ne(ldouble *pp)
{
  ldouble ne_relel=0.0;
  ldouble ne_tot;
  ldouble ne_th;
  
  ne_tot = one_over_mue_mp * pp[RHO];

#ifdef RELELECTRONS
  ne_relel=calc_relel_ne(pp);
  
  //if(ne_relel/ne_tot > MAX_RELEL_FRAC_N)
  //  ne_relel = MAX_RELEL_FRAC_N*ne_tot; 
#endif

  return ne_tot - ne_relel;
}


//*****************************************************
//calculates left and right wave speeds at cell center
//******************************************************
// July 8, 17, Ramesh: This version computes wavespeeds in all three directions simultaneously,
// and save a little by not repeating certain direction-independent quantities.
// Also, computes wavespeeds only for the relevant directions.
//*****************************************************
//calculates left and right wave speeds at cell center
//*****************************************************

int
calc_wavespeeds_lr(int ix, int iy, int iz, ldouble *aaa)
{
  ldouble (*gg)[5],(*GG)[5];

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  //temporary using local arrays
  gg=geom.gg;
  GG=geom.GG;

  ldouble pp[NV];
  int iv;

  //picking up primitives 
  for(iv=0;iv<NV;iv++)
    pp[iv]=get_u(p,iv,ix,iy,iz);

  calc_wavespeeds_lr_pure(pp,&geom,aaa);

  return 0;
}

//*************************************************
//calculates left and right wave speeds at cell center
//*************************************************

int
calc_wavespeeds_lr_pure(ldouble *pp,void *ggg,ldouble *aaa)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  gg=geom->gg;
  GG=geom->GG;

  int iv;
  
  ldouble axhdl,axhdr,ayhdl,ayhdr,azhdl,azhdr;
  ldouble axl0,axr0,ayl0,ayr0,azl0,azr0;
  ldouble axl,axr,ayl,ayr,azl,azr;
  axl=axr=ayl=ayr=azl=azr=1.;
  
  ldouble utcon[4],ucon[4],ucov[4],cst1,cst2,cst3,cst4;
  ldouble bcon[4],bcov[4],bsq;
  ldouble cs2,va2,EF,EEloc; 
  ldouble rho,uu,pre,prerad;

  //**********************************************************************
  //***** four velocity **************************************************
  //**********************************************************************

  for(iv=1;iv<4;iv++)
    utcon[iv]=pp[1+iv];
  utcon[0]=0.;
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

  //**********************************************************************
  //***** hydro: speed of sound ******************************************
  //**********************************************************************
	      
  rho=pp[RHO];
  uu=pp[UU];
 
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(geom->ix,geom->iy,geom->iz);
  #endif
  ldouble gammam1=gamma-1.;
  //gas pressure
  pre=(gamma-1.)*uu;

  ldouble preeff=pre;
  ldouble uueff=uu;
#ifdef GASRADCOUPLEDWAVESPEEDS
#ifdef RADIATION  // gas-radiation coupling is needed only with radiation
  //radiation pressure
  ldouble temp[4],Ehat;
  calc_ff_Ehat(pp,&Ehat,temp,geom);
  prerad=one_third*Ehat;
  
  ldouble fcouple;  //coupling coefficient
  
  // Choose which version of gas-radiation coupling to use
#ifdef GASRADCOUPLEDWAVESPEEDS_SIMPLE
  // Simple scheme just uses the total pressure to estimate hydro wave speed
  // Still doubtful if this works
  fcouple = 1.;
  preeff += prerad;
  uueff += Ehat;
#else
  // More sophisticated method, which estimates the degree of coupling between gas and radiation (fcouple) and estimates the effective pressure for hydro wave speed
  fcouple=estimate_gas_radiation_coupling(pp,ggg);
  if (isnan(fcouple)) fcouple = 1.; // Ramesh: needed sometimes near the inner boundary
  preeff += fcouple*prerad;
  uueff += fcouple*Ehat;
#endif  // ifdef GASRADCOUPLEDWAVESPEEDS_SIMPLE
  
#endif  // ifdef RADIATION  
#endif  // ifdef GASRADCOUPLEDWAVESPEEDS

  //Sound speed
  //ANDREW changed the denominator
  //cs2=gamma*preeff/(rho+uu+preeff); //old version
  cs2=gamma*preeff/(rho+uueff+preeff); //new version 
  
  if(cs2>0.95) cs2=0.95;

  //**********************************************************************
  //***** magn: alfvenic speed ****** *************************************
  //**********************************************************************
	      
  va2=0.;
#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, geom, bcon, bcov, &bsq);
  EF = rho + uu + pre;
  EEloc = bsq + EF ;
  va2 = bsq/EEloc ;
  if(va2<0.) va2=0.;
#endif

  //**********************************************************************
  //***** mhd: fast magnetosonic speed ***********************************
  //**********************************************************************

  ldouble vtot2; //total characteristic velocity
  vtot2=cs2 + va2 - cs2*va2;

#ifdef NONRELMHD //non-rel version
  
  ldouble vx,vy,vz,cs,csx,csy,csz;
  vx=pp[VX];
  vy=pp[VY];
  vz=pp[VZ];
  cs=sqrt(vtot2);
  csx=cs/sqrt(gg[1][1]);
  csy=cs/sqrt(gg[2][2]);
  csz=cs/sqrt(gg[3][3]);
  
  axhdr=vx+csx;
  axhdl=vx-csx;

  ayhdr=vy+csy;
  ayhdl=vy-csy;
  
  azhdr=vz+csz;
  azhdl=vz-csz;

#else //fully relativistic

  //**********************************************************************
  //algorithm from HARM to transform the fluid frame wavespeed into lab frame
  //**********************************************************************

  ldouble aret[6];
  int ret;
  
  ret=calc_wavespeeds_lr_core(ucon,GG,aret,vtot2,vtot2,vtot2);
  if(ret<0) {printf("error occurred at %d | %d | %d\n",geom->ix,geom->iy,geom->iz);}
  axhdl=aret[0];
  axhdr=aret[1];
  ayhdl=aret[2];
  ayhdr=aret[3];
  azhdl=aret[4];
  azhdr=aret[5];

#endif
 
#ifdef RADIATION
  
  //**********************************************************************
  //***** radiation: characteristic wave speed ***************************
  //**********************************************************************

  ldouble aval[18];
  int verbose=0;

  //physical size of the cell
  ldouble dx[3];
  ldouble xx[4]={0.,geom->xx,geom->yy,geom->zz};
  
  //ix,iy,iz could be the indices of a face, so the depth taken from left/right
  dx[0]=my_max(get_size_x(geom->ix,0)*sqrt(gg[1][1]),get_size_x(geom->ix+1,0)*sqrt(gg[1][1]));
  dx[1]=my_max(get_size_x(geom->iy,1)*sqrt(gg[2][2]),get_size_x(geom->iy+1,1)*sqrt(gg[2][2]));
  dx[2]=my_max(get_size_x(geom->iz,2)*sqrt(gg[3][3]),get_size_x(geom->iz+1,2)*sqrt(gg[3][3]));
  ldouble tautot[3];
  calc_tautot(pp,geom,dx,tautot);

  //compute radiative wavespeeds
  calc_rad_wavespeeds(pp,geom,tautot,aval,verbose);

  //unlimited by optical depth
  axl0=aval[0];
  axr0=aval[1];
  ayl0=aval[2];
  ayr0=aval[3];
  azl0=aval[4];
  azr0=aval[5];

  //affected by optical depth - by default scaled as 1/tau
  axl=aval[6+0];
  axr=aval[6+1];
  ayl=aval[6+2];
  ayr=aval[6+3];
  azl=aval[6+4];
  azr=aval[6+5];

  //in the other approach - radiative wavespeeds unlimited in optically thin medium (fcouple==0)
  //and equal to gas wavespeeds in optically thick medium (fcouple==1)
#ifdef GASRADCOUPLEDWAVESPEEDS
  axl=fcouple*axhdl+(1.-fcouple)*axl;
  axr=fcouple*axhdr+(1.-fcouple)*axr;
  ayl=fcouple*ayhdl+(1.-fcouple)*ayl;
  ayr=fcouple*ayhdr+(1.-fcouple)*ayr;
  azl=fcouple*azhdl+(1.-fcouple)*azl;
  azr=fcouple*azhdr+(1.-fcouple)*azr;
#endif
  
#endif //RADIATION

#ifdef OVERWRITERADWAVESPEEDSWITHHD
  axl=axl0=axl2=axhdl;
  axr=axr0=axr2=axhdr;
  ayl=ayl0=ayl2=ayhdl;
  ayr=ayr0=ayr2=ayhdr;
  azl=azl0=azl2=azhdl;
  azr=azr0=azr2=azhdr;
#endif

  //zeroing 'co-going' velocities
  if(axhdl>0.) axhdl=0.;
  if(axhdr<0.) axhdr=0.;
  if(ayhdl>0.) ayhdl=0.;
  if(ayhdr<0.) ayhdr=0.;
  if(azhdl>0.) azhdl=0.;
  if(azhdr<0.) azhdr=0.;

  if(axl>0.) axl=0.;
  if(axr<0.) axr=0.;
  if(ayl>0.) ayl=0.;
  if(ayr<0.) ayr=0.;
  if(azl>0.) azl=0.;
  if(azr<0.) azr=0.;

  if(axl0>0.) axl0=0.;
  if(axr0<0.) axr0=0.;
  if(ayl0>0.) ayl0=0.;
  if(ayr0<0.) ayr0=0.;
  if(azl0>0.) azl0=0.;
  if(azr0<0.) azr0=0.;

  //saving and passing up
  //hd:
  aaa[0]=axhdl;
  aaa[1]=axhdr;
  aaa[2]=ayhdl;
  aaa[3]=ayhdr;
  aaa[4]=azhdl;
  aaa[5]=azhdr;

  //rad:
  //unlimited by optical depth - used for calculation of timestep
  aaa[6]=axl0;
  aaa[7]=axr0;
  aaa[8]=ayl0;
  aaa[9]=ayr0;
  aaa[10]=azl0;
  aaa[11]=azr0;

  //affected by optical depth
  aaa[6+6]=axl;
  aaa[6+7]=axr;
  aaa[6+8]=ayl;
  aaa[6+9]=ayr;
  aaa[6+10]=azl;
  aaa[6+11]=azr;

  return 0;
}

int
calc_wavespeeds_lr_core(ldouble *ucon, ldouble GG[][5], ldouble *aret,
			ldouble wspeed2x, ldouble wspeed2y, ldouble wspeed2z)
{
  int ierr = 0;
  ldouble Acov[4], Acon[4], Bcov[4], Bcon[4];
  ldouble Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, discr, cst1, cst2;
  
  // Compute direction-independent quantities first
  
  Bcov[0] = 1.;
  Bcov[1] = 0.;
  Bcov[2] = 0.;
  Bcov[3] = 0.;
  
  indices_12(Bcov, Bcon, GG);
  Bsq = dot(Bcon, Bcov);
  Bu = dot(Bcov, ucon);
  Bu2 = Bu * Bu;

  // Now work on the relevant directions
  if (TNX > 1)  // x-direction is needed
  {
    Acov[0] = 0.;
    Acov[1] = 1.;
    Acov[2] = 0.;
    Acov[3] = 0.;
    
    indices_12(Acov, Acon, GG);
    
    Asq = dot(Acon, Acov);
    Au = dot(Acov, ucon);
    AB = dot(Acon, Bcov);
    Au2 = Au * Au;
    AuBu = Au * Bu;
    
    B = -2. * (AuBu * (1.0 - wspeed2x)  - AB * wspeed2x);
    A = Bu2 * (1.0 - wspeed2x) - Bsq * wspeed2x;
    discr = 4.0 * wspeed2x * ((AB * AB - Asq * Bsq) * wspeed2x + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2x - 1.0));
    
    if(discr < 0.) {printf("discr in x-wavespeeds lt 0\n"); ierr = -1;}
    discr = sqrt(discr);
    cst1 = (-B + discr) / (2. * A);
    cst2 = (-B - discr) / (2. * A);
    if(cst2 > cst1)
    {
      aret[0] = cst1;  aret[1] = cst2;
    }
    else
    {
      aret[0] = cst2;  aret[1] = cst1;
    }
  }
  
  if (TNY > 1)  // y-direction is needed
  {
    Acov[0] = 0.;
    Acov[1] = 0.;
    Acov[2] = 1.;
    Acov[3] = 0.;
    
    indices_12(Acov, Acon, GG);
    
    Asq = dot(Acon, Acov);
    Au = dot(Acov, ucon);
    AB = dot(Acon, Bcov);
    Au2 = Au * Au;
    AuBu = Au * Bu;
    
    B = -2. * (AuBu * (1.0 - wspeed2y)  - AB * wspeed2y);
    A = Bu2 * (1.0 - wspeed2y) - Bsq * wspeed2y;
    discr = 4.0 * wspeed2y * ((AB * AB - Asq * Bsq) * wspeed2y + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2y - 1.0));
    
    if(discr < 0.) {printf("discr in y-wavespeeds lt 0\n"); ierr = -1;}
    discr = sqrt(discr);
    cst1 = (-B + discr) / (2. * A);
    cst2 = (-B - discr) / (2. * A);
    if(cst2 > cst1)
    {
      aret[2] = cst1;  aret[3] = cst2;
    }
    else
    {
      aret[2] = cst2;  aret[3] = cst1;
    }
  }
  
  if (TNZ > 1)  // z-direction is needed
  {
    Acov[0] = 0.;
    Acov[1] = 0.;
    Acov[2] = 0.;
    Acov[3] = 1.;
    
    indices_12(Acov, Acon, GG);
    
    Asq = dot(Acon, Acov);
    Au = dot(Acov, ucon);
    AB = dot(Acon, Bcov);
    Au2 = Au * Au;
    AuBu = Au * Bu;
    
    B = -2. * (AuBu * (1.0 - wspeed2z)  - AB * wspeed2z);
    A = Bu2 * (1.0 - wspeed2z) - Bsq * wspeed2z;
    discr = 4.0 * wspeed2z * ((AB * AB - Asq * Bsq) * wspeed2z + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2z - 1.0));
    
    if(discr < 0.) {printf("discr in z-wavespeeds lt 0\n"); ierr = -1;}
    discr = sqrt(discr);
    cst1 = (-B + discr) / (2. * A);
    cst2 = (-B - discr) / (2. * A);
    if(cst2 > cst1)
    {
      aret[4] = cst1;  aret[5] = cst2;
    }
    else
    {
      aret[4] = cst2;  aret[5] = cst1;
    }
  }
  
  if (ierr == 0)
  {
    return 0;
  }
  else
  {
    return -1;
  }
}




//*************************************************************
//returns geometrical source terms for all conserved quantities
//*************************************************************

int f_metric_source_term_arb(ldouble *pp,void *ggg,ldouble *ss)
{
  int i;

  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  int ix,iy,iz;
  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;
 
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;
  gdetu=gdet;

  #if (GDETIN==0) //no metric determinant inside derivatives
  gdetu=1.;
  #endif

  ldouble dlgdet[3];
  dlgdet[0]=gg[0][4]; //D[gdet,x1]/gdet
  dlgdet[1]=gg[1][4]; //D[gdet,x2]/gdet
  dlgdet[2]=gg[2][4]; //D[gdet,x3]/gdet
  
  ldouble ut;
  ldouble T[4][4];

  //calculating stress energy tensor components
  calc_Tij(pp,geom,T);
  indices_2221(T,T,gg);

  int ii, jj;
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
      {
	if(isnan(T[ii][jj])) 
	  {
	    printf("%d %d %e\n",ii,jj,T[ii][jj]);
	    my_err("nan in metric_source_terms\n");
	  }
      }
 
  ldouble rho=pp[RHO];
  ldouble u=pp[UU];
  ldouble vcon[4],ucon[4];
  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];
  ldouble S=pp[5];

  //converting to 4-velocity
  conv_vels(vcon,ucon,VELPRIM,VEL4,gg,GG);
  
  int k,l,iv;
  for(iv=0;iv<NV;iv++)
    ss[iv]=0.;  // zero out all source terms initially

#ifdef RADIATION

  ldouble Rij[4][4];
  calc_Rij(pp,geom,Rij); //R^ij
  indices_2221(Rij,Rij,gg); //R^i_j
  
  //terms with Christoffels
  for(k=0;k<4;k++)
    for(l=0;l<4;l++)
      {
	ss[1]+=gdetu*T[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[2]+=gdetu*T[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[3]+=gdetu*T[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[4]+=gdetu*T[k][l]*get_gKr(l,3,k,ix,iy,iz);
	ss[EE0]+=gdetu*Rij[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[FX0]+=gdetu*Rij[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[FY0]+=gdetu*Rij[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[FZ0]+=gdetu*Rij[k][l]*get_gKr(l,3,k,ix,iy,iz);
      }

#if (GDETIN==0)   //terms with dloggdet if gdet not inside the derivatives
#ifdef EVOLVEPHOTONNUMBER
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);
#endif

  for(l=1;l<4;l++)
    {
      ss[0]+=-dlgdet[l-1]*rho*ucon[l];
      ss[1]+=-dlgdet[l-1]*(T[l][0]+rho*ucon[l]);
      ss[2]+=-dlgdet[l-1]*(T[l][1]);
      ss[3]+=-dlgdet[l-1]*(T[l][2]);
      ss[4]+=-dlgdet[l-1]*(T[l][3]);
      ss[5]+=-dlgdet[l-1]*S*ucon[l];
      ss[EE0]+=-dlgdet[l-1]*(Rij[l][0]);
      ss[FX0]+=-dlgdet[l-1]*(Rij[l][1]);
      ss[FY0]+=-dlgdet[l-1]*(Rij[l][2]);
      ss[FZ0]+=-dlgdet[l-1]*(Rij[l][3]);
      
      #ifdef EVOLVEPHOTONNUMBER
      ss[NF0]+=-dlgdet[l-1]*pp[NF0]*urfcon[l];
      #endif

      #ifdef EVOLVEELECTRONS
      ss[ENTRE] += -dlgdet[l-1]*pp[ENTRE]*ucon[l];
      ss[ENTRI] += -dlgdet[l-1]*pp[ENTRI]*ucon[l];
      #ifdef RELELECTRONS
      int ie;
      for (ie=0; ie<NRELBIN; ie++)
	{
	  ss[NEREL(ie)]+=-dlgdet[l-1]*pp[NEREL(ie)]*ucon[l];
	}  
      #endif
      #endif
    }
#endif //GDETIN

#else //ndef RADIATION , pure hydro

  //terms with Christoffels
  for(k=0;k<4;k++)
    for(l=0;l<4;l++)
      {
	ss[1]+=gdetu*T[k][l]*get_gKr(l,0,k,ix,iy,iz);
	ss[2]+=gdetu*T[k][l]*get_gKr(l,1,k,ix,iy,iz);
	ss[3]+=gdetu*T[k][l]*get_gKr(l,2,k,ix,iy,iz);
	ss[4]+=gdetu*T[k][l]*get_gKr(l,3,k,ix,iy,iz);
      }

  //terms with dloggdet  
#if (GDETIN==0)
  for(l=1;l<4;l++)
    {
      ss[0]+=-dlgdet[l-1]*rho*ucon[l];
      ss[1]+=-dlgdet[l-1]*(T[l][0]+rho*ucon[l]);
      ss[2]+=-dlgdet[l-1]*(T[l][1]);
      ss[3]+=-dlgdet[l-1]*(T[l][2]);
      ss[4]+=-dlgdet[l-1]*(T[l][3]);
      ss[5]+=-dlgdet[l-1]*S*ucon[l];
    }   
#endif
#endif //RADIATION

#ifdef SHEARINGBOX 
  //signs the same despite rho u^i u_t evolved because in koral source terms on lhs
  
  //-2 rho (Omega \hat z) x \vec v
  ss[VX]+=gdet*(2.*rho*SHEAROM*pp[VY]);
  ss[VY]+=gdet*(-2.*rho*SHEAROM*pp[VX]);
  ss[VZ]+=0.;

  //2 q Omega^2 x \hat x
  ss[VX]+=gdet*(2.*SHEARQ*rho*SHEAROM*SHEAROM*geom->xx);
  ss[VY]+=0.;
  ss[VZ]+=0.;

  //-rho Omega^2 z \hat z
  ss[VX]+=0.;
  ss[VY]+=0.;
  #ifndef SHEARINGBOXUNSTRATIFIED
  ss[VZ]+=-gdet*rho*SHEAROM*SHEAROM*geom->zz;
  #endif

  //Omega^2 rho \vec v 2 q x \hat x
  ss[UU]+=gdet*SHEAROM*SHEAROM*rho*2.*SHEARQ*geom->xx*pp[VX];

  //Omega^2 rho \vec v z \hat z
  ss[UU]+=gdet*SHEAROM*SHEAROM*rho*geom->zz*pp[VZ];

#endif

  return 0;
}


int f_metric_source_term(int ix, int iy, int iz,ldouble *ss)
{
  int i;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  f_metric_source_term_arb(&get_u(p,0,ix,iy,iz), &geom, ss);

  return 0;
}


//**********************************************************
//returns general source terms for all conserved quantities
//**********************************************************

int f_general_source_term_arb(ldouble *pp,void *ggg,ldouble *ss)
{
  int i,j;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble gdet=geom->gdet, gdetu=gdet;
  
  #if (GDETIN==0) //no metric determinant inside derivatives
  gdetu=1.;
  #endif

  int ix,iy,iz,iv;
  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;
 
  PLOOP(iv) ss[iv]=0.;  // zero out all source terms initially

/***************************************************/
//artificial heating of gas at constant rate per unit mass
#if defined(HEATINGRATEPERMASS) || defined(HEATINGRATEPERMASSSQ) || defined(HEATINGCONSTANT)
  #ifdef EVOLVEELECTRONS
  my_err("artifical heating not implemented for Electron Evolution!");
  #endif

  //gas velocity
  ldouble ucon[4];
  ucon[1]=pp[VX];
  ucon[2]=pp[VY];
  ucon[3]=pp[VZ];
  conv_vels(ucon,ucon,VELPRIM,VEL4,geom->gg,geom->GG);

  //gas density
  ldouble rho=pp[RHO];

  //four-vector of heating to be applied to gas, H^\mu
  ldouble Hmu[4]={0.,0.,0.,0.};
  ldouble Hthat =0.;

  #if defined(HEATINGRATEPERMASS)
  Hthat += HEATINGRATEPERMASS*rho;
  #elif defined(HEATINGRATEPERMASSSQ)
  Hthat += HEATINGRATEPERMASSSQ*rho*rho;
  #endif

  #ifdef HEATINGLIMITINGRHO
  if(rho<HEATINGLIMITINGRHO) Hthat = 0.;
  #endif

  //can have ADDITIONAL constant, density independent heating
  #if defined(HEATINGCONSTANT)
  Hthat += HEATINGCONSTANT;
  #endif
  
  for(i=0;i<4;i++)
    Hmu[i]=Hthat*ucon[i];
  indices_21(Hmu,Hmu,geom->gg); //H_\mu

  //source terms
  ss[UU]+=gdetu*Hmu[0];
  ss[VX]+=gdetu*Hmu[1];
  ss[VY]+=gdetu*Hmu[2];
  ss[VZ]+=gdetu*Hmu[3];
#endif

/***************************************************/
//artificial heating of gas towards prescribed entropy to make disk thin
#if defined(COOLINGTOWARDSENTROPY)
  #ifdef EVOLVEELECTRONS
  my_err("COLINGTOWARDSENTROPY not implemented for Electron Evolution!");
  #endif

  //gas velocity
  ldouble ucon[4];
  ucon[1]=pp[VX];
  ucon[2]=pp[VY];
  ucon[3]=pp[VZ];
  conv_vels(ucon,ucon,VELPRIM,VEL4,geom->gg,geom->GG);

  //local Keplerian orbital time
  ldouble xx[4];
  get_xx_arb(geom->ix,geom->iy,geom->iz,xx,BLCOORDS);
  ldouble r=xx[1];
  ldouble Omk=1./(BHSPIN+sqrt(r*r*r));
  ldouble Pk = 2.*M_PI/Omk;
  ldouble taucool = Pk;

  //current entropy
  ldouble Sloc = log(GAMMAM1*pp[UU]/pow(pp[RHO],GAMMA)); //log(p/rho^Gamma)
  //target entropy
  ldouble Star = TARGETLOGENTROPY;

  //vertical limit for artificial cooling
  ldouble theta = xx[2];
  ldouble thnocool = THETANOCOOL;
  ldouble stheta = exp(-(theta - M_PI/2.)*(theta - M_PI/2.)/2./thnocool/thnocool);

  //comoving cooling rate
  ldouble dudtau;
  if(Sloc>Star) dudtau=-pp[UU]*(Sloc-Star)*stheta/taucool;

  //sanity checks
  if(fabs(theta-M_PI/2.)>(0.9*M_PI/2.)) stheta=0.; //nothing at the axis                                                                                        
  if(dudtau>0.) dudtau=0.; //no heating
  //if(global_dt>taucool) dudtau=0.; //not to overshoot - following Penna
  if(fabs(dudtau) > 0.5*pp[UU]) dudtau=-0.5*pp[UU]; //not to overshoot - mine
  dudtau *= step_function(r - 0.5*(rISCOBL + rhorizonBL),.1*rhorizonBL); //to avoid BH
  if(r<rhorizonBL) dudtau=0.;

  //four-vector of cooling to be applied to gas, C^\mu
  ldouble Cmu[4]={0.,0.,0.,0.};
  for(i=0;i<4;i++)
    Cmu[i]=dudtau*ucon[i];
  indices_21(Cmu,Cmu,geom->gg);  //C_\mu

  //source terms
  ss[UU]+=gdetu*Cmu[0];
  ss[VX]+=gdetu*Cmu[1];
  ss[VY]+=gdetu*Cmu[2];
  ss[VZ]+=gdetu*Cmu[3];
#endif

  return 0;
}


int f_general_source_term(int ix, int iy, int iz,ldouble *ss)
{
  int i;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  f_general_source_term_arb(&get_u(p,0,ix,iy,iz),&geom,ss);

  return 0;
}


//***************************************************************
// calculates fluxes at faces
//***************************************************************
int f_flux_prime(ldouble *pp, int idim, int ix, int iy, int iz,ldouble *ff,int lr)
{  

  int iv;
  for(iv=0;iv<NV;iv++) 
    ff[iv]=0.;

  //picking up metric from a cell face  
  struct geometry geom;
  fill_geometry_face(ix,iy,iz,idim,&geom);

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom.gg;
  GG=geom.GG;
  gdet=geom.gdet;
  gdetu=gdet;

  #if (GDETIN==0) //no metric determinant inside derivative
  gdetu=1.;
  #endif

  //calculating Tij
  ldouble T[4][4];
  calc_Tij(pp,&geom,T);
  indices_2221(T,T,gg);//T^ij --> T^i_j

  //primitives
#ifdef EVOLVEELECTRONS
  ldouble Se=pp[ENTRE]; //entropy of electrons
  ldouble Si=pp[ENTRI]; //entropy of ions
#endif
  ldouble rho=pp[RHO];
  ldouble u=pp[UU];

  ldouble vcon[4],ucon[4],ucov[4],bcon[4],bcov[4],bsq=0.;

  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];
  ldouble S=pp[5];

  //converting to 4-velocity
  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

#ifdef NONRELMHD
  ucon[0]=1.;
  ucov[0]=-1.;
#endif

#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, &geom, bcon, bcov, &bsq);
#endif

  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(ix,iy,iz);
  #endif
  ldouble gammam1=gamma-1.;

  ldouble pre=(gamma-1.)*u; 
  ldouble w=rho+u+pre;
  ldouble eta=w+bsq;
  ldouble etap = u+pre+bsq; //eta-rho

  int ii, jj, irf;
  for(ii=0;ii<4;ii++)
    for(jj=0;jj<4;jj++)
    {
	if(isnan(T[ii][jj])) 
	{
	    printf("%d > nan tmunu: %d %d %e at %d %d %d\n",PROCID,ii,jj,T[ii][jj],ix+TOI,iy+TOJ,iz+TOK);
	    printf("%d > nan tmunu: %e %e %e %e\n",PROCID,gamma,pre,w,eta);
	    print_4vector(ucon);
	    print_metric(geom.gg);
	    print_Nvector(pp,NV);
	    my_err("nan in flux_prime\n");
	    exit(1);
	}
    }

  ldouble utp1=calc_utp1(vcon,ucon,&geom);

  //***************************************
  //fluxes
  //***************************************
  //hydro fluxes
  ff[0]= gdetu*rho*ucon[idim+1];
  
  //ff[1]= gdetu*(T[idim+1][0]+rho*ucon[idim+1]);
  //to avoid slow cancellation:
  ff[1]= gdetu*(etap*ucon[idim+1]*ucov[0] + rho*ucon[idim+1]*utp1);
#ifdef MAGNFIELD
  ff[1]+= -gdetu*bcon[idim+1]*bcov[0];
#endif

  ff[2]= gdetu*(T[idim+1][1]);
  ff[3]= gdetu*(T[idim+1][2]); 
  ff[4]= gdetu*(T[idim+1][3]);
  ff[5]= gdetu*S*ucon[idim+1];

#ifdef NONRELMHD
  ff[1]= gdetu*T[idim+1][0];
#endif

#ifdef EVOLVEELECTRONS
  ff[ENTRE]= gdetu*Se*ucon[idim+1]; 
  ff[ENTRI]= gdetu*Si*ucon[idim+1]; 

#ifdef RELELECTRONS
  int ie;
  for (ie=0; ie < NRELBIN ; ie++)
    ff[NEREL(ie)] = gdetu*pp[NEREL(ie)]*ucon[idim+1];
#endif
#endif

  //magnetic fluxes
#ifdef MAGNFIELD
  ff[B1]=gdetu*(bcon[1]*ucon[idim+1] - bcon[idim+1]*ucon[1]);
  ff[B2]=gdetu*(bcon[2]*ucon[idim+1] - bcon[idim+1]*ucon[2]);
  ff[B3]=gdetu*(bcon[3]*ucon[idim+1] - bcon[idim+1]*ucon[3]);

#ifdef BATTERY //radiation battery
#ifdef RADIATION
  ldouble eterm[4]={-1.,0.,0.,0.};
  calc_batteryflux(pp,&geom,eterm,idim,ucov);
  ldouble suppfac=1.;
 
  ff[B1]+=suppfac*gdetu*eterm[1];
  ff[B2]+=suppfac*gdetu*eterm[2];
  ff[B3]+=suppfac*gdetu*eterm[3];
#endif
#endif
#endif

  //radiation fluxes
#ifdef RADIATION
  f_flux_prime_rad(pp,idim,&geom,ff);
#endif

  return 0;
}

//********************************************************************************************
//calculates energy-momentum tensor components basing on vector of primitivies p and metric g
//returns T^munu
//********************************************************************************************

int
calc_Tij(ldouble *pp, void* ggg, ldouble T[][4])
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  int iv,i,j;
  ldouble rho=pp[RHO];
  ldouble uu=pp[UU];
  ldouble utcon[4],ucon[4],ucov[4];  
  ldouble bcon[4],bcov[4],bsq=0.;
  
  //converts to 4-velocity
  for(iv=1;iv<4;iv++)
    utcon[iv]=pp[1+iv];
  utcon[0]=0.;
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

#ifdef NONRELMHD
  ucon[0]=1.;
  ucov[0]=-1.;
#endif


#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, geom, bcon, bcov, &bsq);
#else
  bcon[0]=bcon[1]=bcon[2]=bcon[3]=0.;
  bsq=0.;
#endif
  
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(geom->ix,geom->iy,geom->iz);
  #endif
  ldouble gammam1=gamma-1.;

  ldouble p=(gamma-1.)*uu; 
  ldouble w=rho+uu+p;
  ldouble eta=w+bsq;
  ldouble ptot=p+0.5*bsq;

#ifndef NONRELMHD  
  for(i=0;i<4;i++)
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
    for(j=0;j<4;j++)
      T[i][j]=eta*ucon[i]*ucon[j] + ptot*GG[i][j] - bcon[i]*bcon[j];

#else
  
  ldouble v2=dot3nr(ucon,ucov);

  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      T[i][j]=(rho)*ucon[i]*ucon[j] + ptot*GG[i][j] - bcon[i]*bcon[j];

  T[0][0]=uu + bsq/2. + rho*v2/2.;

  
   for(i=1;i<4;i++)
     T[0][i]=T[i][0]=(T[0][0] + ptot) *ucon[i]*ucon[0] + ptot*GG[i][0] - bcon[i]*bcon[0];

#endif  // ifndef NONRELMHD

  return 0;
}

//**********************************************************************
//entropy-related routines
//**********************************************************************
//S1
//Gas only, uses gamma from pick_gammagas
//Note the units are scaled by K_BOLTZ/MU_GAS/M_PROTON compared to S2,S3
ldouble
calc_Sfromu(ldouble rho,ldouble u,int ix,int iy,int iz)
{
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(ix,iy,iz);
  #endif
  ldouble gammam1=gamma-1.;
  ldouble indexn=1.0/gammam1;
  ldouble ret;
  ldouble pre=gammam1*u;
  #ifdef NOLOGINS
  ret= rho*u / pow(rho,gamma);
  #else
  ret = rho*log(pow(pre,indexn)/pow(rho,indexn+1.));
  #endif

  return ret;
}

ldouble
calc_SfromT(ldouble rho,ldouble temp,int ix,int iy,int iz)
{
   
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(ix,iy,iz);
  #endif
  ldouble gammam1=gamma-1.;


  ldouble n=rho*one_over_mugas_mp;
  ldouble p=K_BOLTZ*n*temp;
  ldouble u=p/gammam1;

  ldouble ret=calc_Sfromu(rho,u,ix,iy,iz);

  return ret;
}

ldouble
calc_ufromS(ldouble S,ldouble rho,int ix,int iy,int iz)
{						
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(ix,iy,iz);
  #endif
  ldouble gammam1=gamma-1.;
  ldouble indexn=1.0/gammam1;
  ldouble ret;
  
  #ifdef NOLOGINS
  ret = pow(rho,gamma)*(S/rho);
  #else
  ret = pow((pow(rho,indexn+1.)*exp(S/rho)), gammam1)/(gammam1);
  #endif
  
return ret;
}


ldouble
calc_TfromS(ldouble S,ldouble rho,int ix,int iy,int iz)
{						
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(ix,iy,iz);
  #endif
  ldouble gammam1=gamma-1.;

  ldouble n = rho*one_over_mugas_mp;  
  ldouble u = calc_ufromS(S,rho,ix,iy,iz);
  ldouble p = u*gammam1;
  ldouble ret=p/(K_BOLTZ*n);
      
return ret;
}

///////////////////////////////////////////////////////////////
//S2
//appropriate for electrons & ions OR gas, but uses fixed gamma

ldouble
calc_S2fromrhou(ldouble rho, ldouble uint, int type)
{
 ldouble gamma,mu;
  if(type==IONS) 
    {
      gamma=GAMMAI;
      mu=MU_I;
    }
  else if(type==ELECTRONS)
    {
      gamma=GAMMAE;
      mu=MU_E;
    }
  else if(type==GAS) //hardly used
    {
      gamma=GAMMA;
      mu=MU_GAS;
    }
  ldouble gammam1=gamma-1.;
  ldouble indexn=1.0/gammam1;

  ldouble n=rho/(mu*M_PROTON); 
  ldouble pre=uint*gammam1;
  ldouble S2;
  
  #ifdef NOLOGINS2
  S2 = n*uint / pow(n,gamma);
  #else
  S2 = n*K_BOLTZ*log(pow(pre,indexn)/pow(n,indexn+1.)); 
  #endif
  return S2;

}

ldouble
calc_S2fromrhoT(ldouble rho, ldouble temp,int type)
{
  ldouble gamma,mu;
  if(type==IONS) 
    {
      gamma=GAMMAI;
      mu=MU_I;
    }
  else if(type==ELECTRONS)
    {
      gamma=GAMMAE;
      mu=MU_E;
    }
  else if(type==GAS) //hardly used
    {
      gamma=GAMMA;
      mu=MU_GAS;
    }

  ldouble n=rho/(mu*M_PROTON);
  ldouble p=K_BOLTZ*n*temp;
  ldouble u=p/(gamma-1.);

  return calc_S2fromrhou(rho,u,type);
}

ldouble
calc_ufromS2rho(ldouble S2,ldouble rho,int type,int ix,int iy,int iz)
{
  ldouble gamma,mu;
  if(type==IONS) 
    {
      gamma=GAMMAI;
      mu=MU_I;
    }
  else if(type==ELECTRONS)
    {
      gamma=GAMMAE;
      mu=MU_E;
    }
  else if(type==GAS) //hardly used
    {
      gamma=GAMMA;
      mu=MU_GAS;
    }

  ldouble gammam1=gamma-1.;
  ldouble indexn=1.0/gammam1;
  
  ldouble n=rho/(mu*M_PROTON); //number density, thermal only
  ldouble ret;
  #ifdef NOLOGINS2
  ret= (S2/n)*pow(n,gamma);
  #else
  ret= pow(( pow(n,indexn+1.)* exp(S2/(n*K_BOLTZ)) ), gammam1) / (gammam1); 
  #endif
 
 return ret;
}

ldouble
calc_TfromS2rho(ldouble S2,ldouble rho, int type,int ix,int iy,int iz)
{
  ldouble gamma,mu,T;
  if(type==IONS) 
    {
      gamma=GAMMAI;
      mu=MU_I;
    }
  else if(type==ELECTRONS)
    {
      gamma=GAMMAE;
      mu=MU_E;
    }
  else if(type==GAS)
    {
      gamma=GAMMA;
      mu=MU_GAS;
    }

  ldouble n=rho/(mu*M_PROTON); //number density, thermal only
  ldouble u=calc_ufromS2rho(S2,rho,type,ix,iy,iz);
  T=u*(gamma-1.)/(K_BOLTZ*n);

  return T;
}

//////////////////////////////////////////////////////////////////
//S3 - entropy which smoothly transitions from theta=0 to theta=1
//used with CONSISTENTGAMMA
ldouble
calc_S3fromrhoT(ldouble rho, ldouble temp,int type)
{  
  ldouble gamma,mu,theta,n, rhospec;

  if(type==IONS) 
    {
      mu=MU_I;
      theta = kB_over_mui_mp * temp;
    }
  else if(type==ELECTRONS)
    {
      mu=MU_E;
      theta = kB_over_me * temp;
    }
  else if(type==GAS) 
    {
      my_err("calc_S3fromrhoT with type==GAS is inconsistent!");
    }
 
  n=rho/(mu*M_PROTON);
  ldouble S3= n*K_BOLTZ*log(sqrt(theta*theta*theta*(theta+0.4)*(theta+0.4)*(theta+0.4))/n); 
  return S3;
}

ldouble
calc_S3fromrhou(ldouble rho, ldouble uint,int type)
{  
  ldouble mass, gamma, T;
  ldouble  n, p, u;
  ldouble tmin, tminfrac, S;
  int  enttype;
  if (type==IONS)
    {
      mass = M_PROTON*MU_I;
      n = one_over_mui_mp * rho;
    }
  else if (type==ELECTRONS)
    {
      mass = M_ELECTR;
      n = one_over_mue_mp * rho;
    }
  else if(type==GAS)
    {
      my_err("calc_S3fromrhou with type==GAS is inconsistent!");
    }
   
  T=solve_Teifromnmu(n, mass, uint,type); //solves in parallel for gamma and temperature
  S=calc_S3fromrhoT(rho, T, type); 

  return S;
}

//Eq A9 from Sadowski, Wielgus, Narayan, Abarca, McKinney 2016
ldouble
calc_TfromS3rho(ldouble S3,ldouble rho, int type,int ix,int iy,int iz)
{
  ldouble gamma,mu,mass,T;
  if(type==IONS) 
    {
      mu=MU_I;
      mass=M_PROTON*MU_I;
    }
  else if(type==ELECTRONS)
    {
      mu=MU_E;
      mass=M_ELECTR;
    }  
  else if(type==GAS)
    {
      my_err("calc_TfromS3rho with type==GAS is inconsistent!");
    }

  ldouble n=rho/(mu*M_PROTON);
  ldouble ex=exp(S3/(n*K_BOLTZ));

  if(isfinite(ex)) //well defined temperature
    {
      ldouble rhs=cbrt(n*ex*n*ex);
      ldouble theta=0.2*(-1.+sqrt(1.+25.*rhs));
      T =  theta*mass/K_BOLTZ;
    }
  else
    {
      T = BIG; //ceiling put in calc_PEQ_Tei and in floors_mhd
    }

  if(!isfinite(T))
    T = BIG;

  return T;
}

ldouble
calc_ufromS3rho(ldouble S3,ldouble rho,int type,int ix,int iy,int iz)
{
  ldouble uint=0.0;
  ldouble T=calc_TfromS3rho(S3,rho,type,ix,iy,iz);
  if(type==ELECTRONS)
    {
      ldouble ne = one_over_mue_mp * rho; //thermal only
      ldouble pe=K_BOLTZ*ne*T;
      ldouble gammae=GAMMAE;
      #ifdef CONSISTENTGAMMA
      gammae=calc_gammaintfromtemp(T,ELECTRONS);
      #endif
      uint=pe/(gammae-1.);
    }
 else if(type==IONS)
    {
      ldouble ni = one_over_mui_mp * rho; //thermal only
      ldouble pi=K_BOLTZ*ni*T;
      ldouble gammai=GAMMAI;
      #ifdef CONSISTENTGAMMA
      gammai=calc_gammaintfromtemp(T,IONS);
      #endif
      uint=pi/(gammai-1.);
   }
 else if(type==GAS)
   {
      my_err("calc_ufromS3rho with type==GAS is inconsistent!");
   }
    
 return uint;
}


//**********************************************************************
//updates entropy in the specified cell (p[5]) basing on new primitives
//or stays with the old one if entropy u2p solver was involved
// Ramesh: This function is called by init.c in HIGHZBHS;
// we should change all the others
//**********************************************************************
int
update_entropy_cell(int ix,int iy,int iz,int u2pflag)
{
  ldouble gg[4][5],GG[4][5];
  pick_G(ix,iy,iz,GG);
  pick_g(ix,iy,iz,gg);
  ldouble gdet=gg[3][4];
  ldouble gdetu=gdet;
  #if (GDETIN==0) //no metric determinant inside derivatives
  gdetu=1.;
  #endif

  ldouble ucon[4],ut,S,Sut,rho,uu;
  int iv;

  //density and energy density
  rho=get_u(p,0,ix,iy,iz);
  uu=get_u(p,1,ix,iy,iz);

  //converts to 4-velocity
  for(iv=1;iv<4;iv++)
    ucon[iv]=get_u(p,iv+1,ix,iy,iz);
  ucon[0]=0.;
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  ut=ucon[0];

  //u2p_hot worked
  if(u2pflag==0 && uu>0. && rho>0.)
    {
      S=calc_Sfromu(rho,uu,ix,iy,iz);      
      set_u(p,5,ix,iy,iz,S);
      set_u(u,5,ix,iy,iz,S*ut*gdetu); 
    }

  return 0;
}

//********************************************************************************
//Reset ion entropy from electron entropy, number density, and energy conservation
//********************************************************************************
ldouble
entri_from_entre_energy_cons(ldouble *pp, int ix, int iy, int iz)
{
 ldouble ue, ui, uur, rhoeth, neth, entri;
 //ldouble nith = one_over_mui_mp * pp[RHO];

 #ifdef RELELECTRONS
 uur = calc_relel_uint(pp);
 neth = calc_thermal_ne(pp);
 rhoeth = (MU_E * M_PROTON) * neth; 
 #else
 uur = 0.0;
 rhoeth=pp[RHO];
 #endif

 ue = calc_ufromSerho(pp[ENTRE],rhoeth,ELECTRONS,ix,iy,iz); 

 ui = pp[UU] - ue - uur;

 if (ui < UIUINTMINRATIO*pp[UU]) 
 {
  ui = UIUINTMINRATIO*pp[UU];
 }

 entri = calc_Sefromrhou(pp[RHO],ui,IONS);

 return entri;
}

//********************************************************************************
//Calculate both species temperatures from entropies
//And total gas temperature from internal energy density
//********************************************************************************
ldouble calc_PEQ_Teifrompp(ldouble* pp, ldouble* Te, ldouble* Ti,int ix, int iy, int iz)
{
  ldouble Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);
  
  ldouble Tiloc,Teloc;
#ifndef EVOLVEELECTRONS
  Tiloc=Teloc=Tgas;

#ifdef FIXEDTETIRATIO
  ldouble factor=FIXEDTETIRATIO;
  Teloc=  (factor*MU_E*MU_I*Tgas)/(MU_GAS * (MU_E + factor * MU_I));
  Tiloc = Teloc/factor;
#endif

#else //EVOLVEELECTRONS

  ldouble neth = calc_thermal_ne(pp);
  ldouble rhoeth = MU_E * M_PROTON * neth;
  Teloc=calc_TfromSerho(pp[ENTRE],rhoeth,ELECTRONS,ix,iy,iz);
  Tiloc=calc_TfromSerho(pp[ENTRI],pp[RHO],IONS,ix,iy,iz);

#endif //EVOLVEELECTRONS
  
  //!AC -- modified for  issues with gammagas in  corners
  //!AC -- this will effectively make all corner gammagas 4/3
  // I don't think that matters though...
  if(!isfinite(Teloc)) Teloc=BIG;
  if(!isfinite(Tiloc)) Tiloc=BIG;
  
  if(Teloc<TEMPEMINIMAL) Teloc=TEMPEMINIMAL;
  if(Tiloc<TEMPIMINIMAL) Tiloc=TEMPIMINIMAL;

  *Te=Teloc;
  *Ti=Tiloc;

  return Tgas;
}

//********************************************************************************
//other calc_PEQ functions for thermodynamic quantities
//********************************************************************************
//Gas energy from temperature
ldouble
calc_PEQ_ufromTrho(ldouble T,ldouble rho,int ix,int iy,int iz)
{
   ldouble gamma=GAMMA;
   #ifdef CONSISTENTGAMMA
   gamma=pick_gammagas(ix,iy,iz); 
   #endif

   ldouble p = kB_over_mugas_mp * rho * T;
   ldouble u=p/(gamma-1.);
   return u;
}

//Total gas temperature from energy density (uses precalculated gamma)
ldouble
calc_PEQ_Tfromurho(ldouble u,ldouble rho,int ix,int iy,int iz)
{
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(ix,iy,iz); 
  #endif
  
  ldouble p = u*(gamma-1.);
  ldouble T = mugas_mp_over_kB * p / rho;

  return T;
}

//Total gas temperature from pressure
ldouble
calc_PEQ_Tfromprho(ldouble p,ldouble rho,int ix,int iy,int iz)
{
  ldouble T = mugas_mp_over_kB * p / rho;

  return T;
}

//Total gas energy density from species temperatures
//ANDREW this does NOT work with rel electrons, but is only used in some init.c files
ldouble
calc_PEQ_ugasfromrhoTei(double rho,ldouble Te,ldouble Ti,int ix,int iy,int iz)
{
  ldouble gamma=GAMMA;
  ldouble gammae=GAMMAE;
  ldouble gammai=GAMMAI;
  #ifdef CONSISTENTGAMMA
  //gamma=pick_gammagas(ix,iy,iz);
  gamma=calc_gammaintfromTei(Te,Ti); //ANDREW in case gammagas not set yet in init.c
  #ifndef FIXEDGAMMASPECIES
  gammae=calc_gammaintfromtemp(Te,ELECTRONS);
  gammai=calc_gammaintfromtemp(Ti,IONS);
  #endif
  #endif

  ldouble mue,mui,mug;
  mui=MU_I;
  mue=MU_E;
  mug=MU_GAS;

  ldouble ue = kB_over_mue_mp * rho * Te / (gammae-1.);
  ldouble ui = kB_over_mui_mp * rho * Ti / (gammai-1.);
  ldouble u = ue + ui;
  return u;
}

//********************************************************************************
//Calculate gamma_gas for all cells and fill gammagas array
//********************************************************************************
//type: 0 - domain + ghost cells + corners
//type 1 - domain only
//type 2 - domain + ghost cells
//type 3 - domain + ghost cells + corners
int
set_gammagas(int type)
{
#ifdef CONSISTENTGAMMA
  int ii;
  if(type==0) //domain + ghost cells + corners
    {
#pragma omp parallel for private(ii) schedule (static)
      for(ii=0;ii<Nloop_5;ii++) 
	{ 
	  int ix,iy,iz;
	  ix=loop_5[ii][0];      iy=loop_5[ii][1];      iz=loop_5[ii][2];
	  ldouble gamma = calc_gammagas(&get_u(p,0,ix,iy,iz),ix,iy,iz);
	  set_u_scalar(gammagas,ix,iy,iz,gamma);  
	}
    }
  if(type==2) //ghost cells
    {
#pragma omp parallel for private(ii) schedule (static)
      for(ii=0;ii<Nloop_02;ii++) 
	{ 
	  int ix,iy,iz;
	  ix=loop_02[ii][0];      iy=loop_02[ii][1];      iz=loop_02[ii][2];
	  if(if_indomain(ix,iy,iz)==1) continue;
	  ldouble gamma=calc_gammagas(&get_u(p,0,ix,iy,iz),ix,iy,iz);
	  set_u_scalar(gammagas,ix,iy,iz,gamma);
	}
    }
  if(type==3) //domain + ghost cells + corners
    {
#pragma omp parallel for private(ii) schedule (static)
      for(ii=0;ii<Nloop_5;ii++) 
	{ 
	  int ix,iy,iz;
	  ix=loop_5[ii][0];      iy=loop_5[ii][1];      iz=loop_5[ii][2];
	  if(if_indomain(ix,iy,iz)==1) continue;
	  ldouble gamma=calc_gammagas(&get_u(p,0,ix,iy,iz),ix,iy,iz);
	  set_u_scalar(gammagas,ix,iy,iz,gamma);
	}
    }
  if(type==1) //domain only
    {
#pragma omp parallel for private(ii) schedule (static)
      for(ii=0;ii<Nloop_0;ii++) 
	{ 
	  int ix,iy,iz;
	  ix=loop_0[ii][0];      iy=loop_0[ii][1];      iz=loop_0[ii][2];
	  ldouble gamma=calc_gammagas(&get_u(p,0,ix,iy,iz),ix,iy,iz);
	  set_u_scalar(gammagas,ix,iy,iz,gamma);
	}
    }

#endif //CONSISTENTGAMMA
  return 0;
}

//********************************************************************************
//Get adiabatic index of gas from memory
//********************************************************************************
ldouble pick_gammagas(int ix,int iy,int iz)
{
  ldouble gamma;

#ifdef CONSISTENTGAMMA
  gamma=get_u_scalar(gammagas,ix,iy,iz);
#else
  gamma = GAMMA;
#endif

  return gamma;
}

//********************************************************************************
// Calculate gamma_int for total gas
//********************************************************************************
ldouble
calc_gammagas(ldouble* pp,int ix,int iy,int iz)
{
  ldouble gamma=GAMMA;

#ifndef FORCEGAMMAGASFIXED
#ifdef CONSISTENTGAMMA
  ldouble Te,Ti;
  calc_PEQ_Teifrompp(pp,&Te,&Ti,ix,iy,iz);

  #ifdef RELELECTRONS
  gamma=calc_gammaint_relel(pp,Te,Ti); 
  #else
  gamma=calc_gammaintfromTei(Te,Ti);
  #endif

  //!AC
  /*
  ldouble Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);
  if(!isfinite(Te) || !isfinite(Ti) || !isfinite(Tgas))
    {
      int gix, giy, giz;
      mpi_local2globalidx(ix, iy, iz, &gix, &giy, &giz);
      printf("Te/Ti not finite >>> %d: > %d %d %d > %e %e %e %e %e %e\n",PROCID,gix,giy,giz,Te,Ti,Tgas,pp[RHO],pp[UU],gamma); 
      fflush(stdout);
      //getch();
      //exit(-1);
    }
  */
  
  if(!isfinite(gamma))
    {
      ldouble Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);
      int gix, giy, giz;
      mpi_local2globalidx(ix, iy, iz, &gix, &giy, &giz);
      printf("gammagas not finite >>> %d: > %d %d %d > %e %e %e %e %e %e\n",PROCID,gix,giy,giz,Te,Ti,Tgas,pp[RHO],pp[UU],gamma); 
      fflush(stdout);
      //getch();
      //exit(-1);
    }
#endif
#endif

  return gamma;
}

//Eq A14 from Sadowski, Wielgus, Narayan, Abarca, McKinney 2016
//ANDREW with rel electrons, use calc_gammaint_relel
ldouble
calc_gammaintfromTei(ldouble Te, ldouble Ti)
{
  ldouble tratio=Ti/Te;
  ldouble the = kB_over_me * Te;
  ldouble thi = kB_over_mui_mp * Ti;
  ldouble gammae=GAMMAE;
  ldouble gammai=GAMMAI;
  #ifndef FIXEDGAMMASPECIES
  gammae=calc_gammaintfromtheta(the);
  gammai=calc_gammaintfromtheta(thi);
  #endif
  ldouble gamma;

  //the original formulae (mui==mue)
  //gamma= 1.+(gammae-1.)*(gammai-1.)*(1.+Ti/Te)/
  //((gammai-1.)+(gammae-1.)*Ti/Te);

  //Andrew's generalized one:
  gamma =  1.+((gammae-1.)*(gammai-1.)*(mui_over_mue + tratio))/
              ((gammae-1.)*tratio + (gammai-1.)*mui_over_mue);
  
  return gamma;
}

//Eq A14 from Sadowski, Wielgus, Narayan, Abarca, McKinney 2016
//used for the equation of state
ldouble
calc_gammaintfromtemp(ldouble temp, int type)
{
  ldouble theta;
  ldouble gint;
  if(type==IONS) 
    {
      theta = kB_over_mui_mp * temp;
    }
  if(type==ELECTRONS)
    {
      theta = kB_over_me * temp;
    }
  if(type==GAS)
    {
      theta = kB_over_mugas_mp * temp;
    }
  gint = calc_gammaintfromtheta(theta);
  return gint;
}

ldouble
calc_gammaintfromtheta(ldouble theta)
{
  ldouble gamma;
  
  #ifdef GAMMAINTCONSISTENTWITHCV
  gamma = 1.+theta/(3.*theta - 0.6 * log1p(2.5*theta));
  #else
  gamma = (10. + 20. * theta) / (6. + 15. * theta);
  #endif

  return gamma;
}

ldouble calc_meanlorentz(ldouble theta)
{
  ldouble gamma;
  if(theta<1.e-3)
  {
    gamma = 1.5*theta + 1.;
  }
  else if(theta>1.e3)
  {
    gamma = 3*theta;
  }
  else //Interpolate in the transition zone
  {
    ldouble xval = log10(theta);
    //printf("theta: %e, xval: %e\n",theta, xval);
    int n=61;
    ldouble xarr[61] = {-3., -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2., \
	    -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1., -0.9, \
	    -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, \
	     0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, \
	     1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.};

    ldouble gammarr[61] = {1.0015, 1.00189, 1.00238, 1.003, 1.00378, 1.00476, 1.006, 1.00756, \
	       1.00954, 1.01203, 1.01519, 1.01918, 1.02424, 1.03066, 1.03883, \
	       1.04925, 1.06257, 1.07966, 1.10165, 1.13008, 1.16699, 1.2151, 1.2781, \
	       1.36087, 1.46995, 1.6139, 1.80386, 2.05412, 2.38271, 2.81217, \
	       3.37044, 4.09199, 5.01943, 6.20561, 7.71639, 9.63422, 12.0626, \
	       15.1319, 19.006, 23.8917, 30.0494, 37.8071, 47.5782, 59.8828, \
	       75.3764, 94.8841, 119.445, 150.366, 189.295, 238.305, 300.005, \
	       377.682, 475.471, 598.581, 753.568, 948.685, 1194.32, 1503.56, \
	       1892.87, 2382.99, 3000.};
    //locate
    int ipos;
    int jl,ju,jm;
    jl = 0;
    ju = n+1;
    //jm = (ju+jl)/2;
    while((ju-jl) > 1)
    {
      jm = (ju+jl)/2;
      //printf("%d %d %d \n",jl,ju,jm);
      if(xval>=xarr[jm-1])
	jl = jm;
      else
	ju = jm;
    }
    if(xval==xarr[0]) ipos = 0;
    else if(xval==xarr[n-1]) ipos=n-2;
    else if(xval<xarr[0]) ipos=-1;
    else if(xval>xarr[n-1]) ipos = n-1;
    else ipos = jl-1;
    
    //printf("\njl: %d , ipos: %d\n",jl, ipos);

    //interpolate
    if(ipos<(n-1) && ipos>-1)
    {
       ldouble frac = (xval - xarr[ipos])/(xarr[ipos+1] - xarr[ipos]);
       gamma = gammarr[ipos] + frac*(gammarr[ipos+1] - gammarr[ipos]);
    }
    else if(ipos==n-1)
    {
      gamma = gammarr[n-1];
    }
    else if (ipos==-1)
    {
      gamma = gammarr[0];
    }
    else
    {
      gamma = theta;
      printf("error in interp gammaavg!\n");
    }
  }
  return gamma;
}

//********************************************************************************
//Solves for species temperature from number density, mass, and energy
//********************************************************************************

//species=1 IONS
//species=2 ELECTRONS
ldouble
solve_Teifromnmu(ldouble n, ldouble m, ldouble u, int species)
{
  ldouble T;
  #if defined(FIXEDGAMMASPECIES)
  ldouble gamma;
  if(species==IONS) gamma=GAMMAI;
  else if(species==ELECTRONS) gamma=GAMMAE;
  T = (gamma-1.)*u/(n*K_BOLTZ);
  #elif defined(GAMMAINTCONSISTENTWITHCV)
  T = solve_Teifromnmu_consistent(n, m, u);
  #else
  T = solve_Teifromnmu_inconsistent(n, m, u);
  #endif
  return T;
}

//solves for gamma_int from the equation of state with the inconsistent gamma_int equation
ldouble
solve_Teifromnmu_inconsistent(ldouble n, ldouble m, ldouble u)
{
  ldouble k=K_BOLTZ;
  ldouble T;
  
  if (u > m * n)
  // regular form of the quadratic solution
  {
    T = (-6.*k*m*n + 5.*k*u + k*sqrt(36.*pow(m,2)*pow(n,2) + 180.*m*n*u + 25.*pow(u,2))) / (30.*pow(k,2)*n); 
  }
  else
  // this form of the quadratic solution is better for very small values of u
  {
    T = 8* m*u / (6.*k*m*n - 5.*k*u + k*sqrt(36.*pow(m,2)*pow(n,2) + 180.*m*n*u + 25.*pow(u,2))); 
  }
  
  return T;
}

//Solves for the temperature corresponding to a given u using the self-consistent formula for u
//Uses Ramesh's Newton-Raphson routine
ldouble
solve_Teifromnmu_consistent(ldouble n, ldouble m, ldouble u)
{
  // Obtain a first approximation for the temperature using Olek's routine
  ldouble T = solve_Teifromnmu_inconsistent(n, m, u);
  
  // If necessary, improve the solution using the self-consistent formula
  ldouble factor = K_BOLTZ / m;
  ldouble theta = factor * T;
  
  if (theta < 6.e-7) // Olek's routine is already accurate enough
  {
    return T;
  }
  else
  {
    // imax is the number of iterations of Newton-Raphson iterations needed for 12 digit accuracy.
    // The ranges below have been calibrated with tests.
    int imax;
    if (theta < 0.0035 || theta > 1.5e4)
    {
      imax = 1;
    }
    else if (theta < 0.1 || theta > 17.)
    {
      imax = 2;
    }
    else
    {
      imax = 3;
    }
    
    int i;
    for (i = 0; i < imax; i++)
    {
      theta -= uint_function(n, m, u, theta) / duint_dtheta(n, m, theta);
      //printf ("i, theta, temperature: %d %e %.12e\n", i, theta, theta / factor);
    }
    
    return theta / factor;
  }
}

// Self-consistent formula for u(theta)
// sam equation as in calc_gammaintfromtheta(ldouble theta)
ldouble
uint_function(ldouble n, ldouble m, ldouble u,ldouble theta)
{
  // Use log1p function instead of log for better accuracy at low temperatures
  ldouble uint = n*m * ( 3. * theta - 0.6 * log1p(2.5 * theta) );

  //printf("theta, uint: %.12e %.12e\n", theta, uint);
  return uint - u;
}


// du/dtheta corresponding to the self-consistent formula for u(theta)
ldouble
duint_dtheta(ldouble n, ldouble m, ldouble theta)
{
  ldouble gamma_CV = (5. + 20. * theta) / (3. + 15. * theta);
  ldouble duintdtheta = n * m / (gamma_CV - 1.);

  //printf("theta, gamma_CV, duintdtheta: %.12e %.12e %.12e\n", theta, gamma_CV, duintdtheta);
  return duintdtheta;
}

//********************************************************************************
//********************************************************************************
//viscous heating of electrons and ions
//applied after the explicit operator
//compares adiabatic and non-adiabatic evolution, 
//identifies viscous heating as the difference of internal energies between the two,
//and applies the surplus to electrons and ions according to delta_e
//the very heating is positive only on average!
//ANDREW Modified heavily for rel. electrons
//********************************************************************************
//********************************************************************************

ldouble
heat_electronions_with_state(ldouble dtin)
{
  int ii;

#ifdef EVOLVEELECTRONS
  
  #pragma omp parallel for schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      ldouble dt, dtau;
      int ix,iy,iz,iv;
      ldouble vcon[4];
      ldouble ucon[4];
      ldouble ucov[4];
      ldouble pp0[NV],pp[NV],uu0[NV];

      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      if(is_cell_active(ix,iy,iz)==0) 
	continue;

      //timestep
      dt=dtin;

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom); 
      ldouble gdet=geom.gg[3][4], gdetu=gdet;

      #if (GDETIN==0) //no metric determinant inside derivatives
      gdetu=1.;
      #endif

      //read primitives and conserved from memory
      PLOOP(iv)
      {
	pp[iv]=pp0[iv]=get_u(p, iv, ix, iy, iz);
	uu0[iv]=get_u(u, iv, ix, iy, iz);
      }

      int itergamma=0;  
      const int maxitergamma=1;
      struct struct_of_state state;
      
      //iterations of (inversion + heating) with updated gammagas
      for(itergamma=0; itergamma<maxitergamma; itergamma++)
      {
	
	  //invert locally, in case gammagas has changed
	  //ANDREW -- u2p with no fixups causes crashes with cylindrified coordinates!!!	
	  int corrected[3]={0,0,0}, fixups[2]={0,0};
	  if(itergamma>0)
	    u2p(uu0,pp,&geom,corrected,fixups,0);

	  //ANDREW -- skip heating operator if u2p failed
	  if(fixups[0]>0 || fixups[1]>0)
	  {
            PLOOP(iv) pp[iv] = pp0[iv];
	    break;
	  }
	  

#ifdef HEATELECTRONS
          fill_struct_of_state(pp, &geom, &state);
	  dtau = dt/(state.ucon[0]);

	  if(!isfinite(pp[ENTRE]))
	  {
	    printf("start: not finite ENTRE at  %d %d %d\n",ix,iy,iz);
	    //getch();
	  }
	  /**************/
	  //Total gas
	  /**************/

	  ldouble rho=pp[RHO];
	  ldouble uint=pp[UU];

	  ldouble Tg=state.Tgas;
	  
	  ldouble gamma=state.gamma;
	  ldouble gammam1=gamma-1.;
	  ldouble ptot=state.pgas;

	  /**************/
	  //electrons
	  /**************/
	  ldouble Te=state.Te;
          ldouble theta_e = kB_over_me * Te;	  
	  ldouble ne=state.ne; 
	  ldouble pe=state.pe;
	  ldouble gammae=state.gammae;
	  ldouble ue=state.ue;
	  ldouble ueorg=ue;

	  /**************/
	  //ions
	  /**************/
	  ldouble Ti=state.Ti;	  
          ldouble ni = state.ni;
	  ldouble pi= state.pi;
	  ldouble gammai=state.gammai;
	  ldouble ui=state.ui;
      
	  /**************/
	  //nonthermal
	  /**************/
	  ldouble utotrelel = state.uenth;
	  ldouble ptotrelel = state.penth;
 
	  /**************/
	  //total species
	  /**************/
	  ldouble ptotentr = pe+pi+ptotrelel;
	  ldouble utotentr = ue+ui+utotrelel;

	  /***********************/
          //calculate dissipation
	  /***********************/

	  
	  ldouble du = uint - utotentr;
	   
#ifdef DISSIPATIONFROMGASONLY
	  ldouble ut=state.ucon[0];
	  ldouble spost=get_u(u,ENTR,ix,iy,iz)/ut/gdetu;
	  ldouble utotentrgas=calc_ufromS(spost,pp[RHO],ix,iy,iz);
	  du = uint - utotentrgas;	  
#endif

          #ifdef NOTALLOWFORNEGATIVEHEATING
	  if(du<0.) du=0.;
          #endif

	  //track viscous heating
	  set_u_scalar(vischeating,ix,iy,iz,du/dtau);
	  
	  /**********************************************/
	  //determine fraction going into electrons 
	  /***********************************************/
	 
          ldouble fe=calc_ViscousElectronHeatingFraction_from_state(pp,&state,&geom);

	  //to calcuate avg delta_e later
	  set_u_scalar(vischeatingtimesdeltae,ix,iy,iz,du*fe/dtau);

	  /************************************/
	  //Heat the nonthermal distribution
	  /************************************/
	  ldouble frel=0.0; //fraction of electron energy injection into relel
	  ldouble p_index=10;
	  ldouble gamma_injmin, gamma_injmax;
	  ldouble utotrelel2=0.0;
	  ldouble durelel=0.0;
	  ldouble due2=0.0;
          ldouble dne=0.0;

#ifdef RELELECTRONS
#ifndef SKIPRELELHEATING
	  //Determine injection power law slope and fraction
	  
	  #ifdef RELEL_HEAT_RECONNECTION //David's reconnection  fits
          reconnection_plaw_params_from_state(pp, &geom, &state, &frel, &p_index);
          #else //fixed frel & p_index

          #ifdef RELEL_HEAT_FIX_FRAC 
          frel=RELEL_HEAT_FRAC; 
          #endif
	  
	  #ifdef RELEL_HEAT_FIX_INDEX
	  p_index=RELEL_HEAT_INDEX;
	  #endif

          #endif //RELEL_HEAT_RECONNECTION

	  //Determine gamma_inj_min & gamma_inj_max
	  #ifndef RELEL_HEAT_FIX_LIMITS
	  
	  ldouble bsq_cgs = 0.;

          #ifdef RELEL_SYN_ART
	  bsq_cgs = RELEL_B_ART*RELEL_B_ART;
          #else
	  bsq_cgs = 4.*M_PI*endenGU2CGS(state.bsq);
          #endif
	  
	  ldouble dtau_cgs = timeGU2CGS(dtau);
	  gamma_injmax= calc_gammainj_max_syncool(bsq_cgs, dtau_cgs); //from synchrotron cooling
	  gamma_injmin= calc_gammainj_min_jointhermal(theta_e, frel, p_index, gamma_injmax); //from continuity w/ thermal

          #else // fix injection minimum and maximum 
	  gamma_injmax=RELEL_INJ_MAX;	  
	  gamma_injmin=RELEL_INJ_MIN;
          #endif
	  
          //energy density injected into relativistic distribution

	  durelel = frel*fe*du;  
          
          //ceilings and floors to limit heating
	  if((durelel+utotrelel)>uint*MAX_RELEL_FRAC_U)
	    durelel=uint*MAX_RELEL_FRAC_U - utotrelel;
          
	  else if((durelel+utotrelel)<0.)
	    durelel=0.0;

          #ifdef NORELELHEATATBH // no rel heating inside bh
          ldouble xxBLh[4];
          #ifdef PRECOMPUTE_MY2OUT
          get_xxout(geom.ix, geom.iy, geom.iz, xxBLh);
          #else
          coco_N(geom.xxvec,xxBLh,MYCOORDS,BLCOORDS);
          #endif
	  
          if(xxBLh[1]<=rhorizonBL)
	    durelel=0.;
          #endif     

          #ifdef NORELELNEGHEAT // no negative rel heating
	  if(durelel<0.)
	    durelel=0.;
          #endif

	  //heating parameters for test problems
	  //RELEL_HEAT_NORM is total injected NUMBER DENSITY (cgs) 
          #ifdef RELEL_HEAT_ART
	  ldouble xxx;
	  
	  xxx = (pow(gamma_injmax, 2.0-p_index) - pow(gamma_injmin, 2.0-p_index))/(2.0-p_index);
	  xxx *= (1.0-p_index)/(pow(gamma_injmax, 1.0-p_index) - pow(gamma_injmin, 1.0-p_index));
          xxx -= 1.;
	  
	  du = endenCGS2GU(xxx * RELEL_HEAT_NORM * M_ELECTR_CGS * CCC_CGS * CCC_CGS)*timeGU2CGS(dtau);
          fe=1.;
	  frel=1.; //ANDREW override??
          durelel = du;

	  //printf("%e %e %e %e\n",durelel, p_index, gamma_injmin, gamma_injmax);
	  //exit(-1);
          #endif    //RELEL_HEAT_ART

	  //Add electrons to nonthermal population
          apply_relel_visc_heating(pp, durelel, p_index, gamma_injmin, gamma_injmax, dtau);
           
          // Ensure energy is conserved by recomputing the *actual* energy fraction put into nonthermal
          utotrelel2 = calc_relel_uint(pp);
	  durelel = utotrelel2 - utotrelel;
          if (du!=0.0 && fe!=0.0)
	    frel = durelel / (du*fe); 

	  // Optional: Remove nonthermal electrons below the thermal peak 
          // TODO: are we doing anything with due2? 
          #ifdef ZERO_NONTHERMAL_LOWGAMMA
	  due2 = remove_lowgamma_electrons(theta_e, ne, pp);
          #endif

          //The (total) change of thermal number density
	  dne = calc_thermal_ne(pp) - ne; 

	  // Increase the total gas energy density if heating is artificial
          #ifdef RELEL_HEAT_ART
	  uint += durelel;
	  pp[UU] = pp[UU] + durelel;
          #endif
	  
#endif //SKIPRELELHEATING
#endif //RELELECTRONS

	  /************************************************************/
	  //apply the rest of the heating to thermal electrons and ions
	  /************************************************************/
	  ldouble due=fe*(1.-frel)*du;
	  ldouble dui=(1.- fe)*du;

          #ifdef NOHEATATBH
          ldouble xxBL[4];

          #ifdef PRECOMPUTE_MY2OUT
          get_xxout(geom.ix, geom.iy, geom.iz, xxBL);
          #else
          coco_N(geom.xxvec,xxBL,MYCOORDS,BLCOORDS);
          #endif
	  
          if(xxBL[1]<=rhorizonBL)
	    {
	      due=0.;
	      dui=0.;
	    }
          #endif
	  
	  //Apply floors
	  //electrons
	  if((ue+due) < UEUINTMINRATIO*uint) 
	  {  
	      due=UEUINTMINRATIO*uint - ue;
	  }  
	  //ions
	  if((ui+dui) < UIUINTMINRATIO*uint) 
	  {
	      dui=UIUINTMINRATIO*uint - ui;
	  }
	  
	  
	  //Update electron and ion entropies to refelect heating
	  ldouble dni = 0.0; // Never any ion number change
	  ldouble Senew, Sinew;
	 	  
          if (!isfinite(due) || !isfinite(dui)) {
            printf("due or dui not finite!!");
	    printf("due: %e dui: %e du: %e frel: %e fe: %e \n",due,dui,du,frel,fe);
            due=0.;
	    dui=0.;
	    //getch();
          }

	  ldouble ue2, ui2;
	  ue2 = apply_du_dn_2_species(pp, ue, ne, due, dne, &geom, ELECTRONS, &Senew);    
	  ui2 = apply_du_dn_2_species(pp, ui, ni, dui, dni, &geom, IONS, &Sinew);
          pp[ENTRE] = Senew;
          pp[ENTRI] = Sinew;

#endif //HEATELECTRONS

	  //updates the gamma of the gas
	  //ANDREW -- currently not important since maxitergamma=1
	  ldouble gammanew=calc_gammagas(pp, ix, iy, iz);
	  set_u_scalar(gammagas, ix, iy, iz, gammanew);
	  
      } //rerun the inversion and heating application with the new gammagas
      
      //Save final pp to memory
      PLOOP(iv)
	set_u(p,iv,ix,iy,iz,pp[iv]);

      //perform p2u
      p2u_mhd(pp,&get_u(u,0,ix,iy,iz),&geom);
    }

#endif //EVOLVEELECTRONS

  return 0;
}


//**************************************************************
//calculates fraction of dissipative heat going into electrons
//**************************************************************

ldouble calc_ViscousElectronHeatingFraction(ldouble *pp,void *ggg)
{
  struct geometry *geom
    = (struct geometry*) ggg;
  struct struct_of_state state;
  fill_struct_of_state(pp,geom,&state);

  ldouble delta = calc_ViscousElectronHeatingFraction_from_state(pp,&state,geom);

  return delta;
}

ldouble calc_ViscousElectronHeatingFraction_from_state(ldouble *pp,void *sss, void *ggg)
{
  ldouble delta=0.;

 struct geometry *geom
  = (struct geometry *) ggg;
 struct struct_of_state *state
  = (struct struct_of_state *) sss;

#ifdef HEATELECTRONS 
  delta=0.;
   
#if defined(HEATELECTRONS_DELTA)
  delta=HEATELECTRONS_DELTA;

#elif defined(HEATELECTRONS_HOWES)
  ldouble gamma=state->gamma;
  ldouble gammam1=gamma-1.;


  // get beta = gas pressure / magn pressure
  ldouble pion = state->pi;  
  ldouble bsq=state->bsq;
  ldouble Ti = state->Ti;
  ldouble Te = state->Te;
  ldouble beta=0.;
  
  #ifdef MAGNFIELD
  
  //ANDREW change beta -> beta_i in Howes
  ldouble betaion = 2.*pion/bsq;
  beta=betaion;
  
  if(!isfinite(beta) || beta>1.e20) beta=1.e20; //no magnetic field 
  #endif

  ldouble me = M_ELECTR_CGS;
  ldouble mp = M_PROTON_CGS;
  ldouble mi = MU_I*mp;
  //ldouble mi = mp*(HFRAC*1. + HEFRAC*4.);

  ldouble tratio=Ti/Te;
  ldouble mratio=mi/me;
  ldouble logtratio = log10(tratio);
  
  //fraction of total viscous heating going to electrons
  ldouble c1 = 0.92;
  ldouble c2 = 1.6*Te/Ti;
  if(Ti<Te) c2 = 1.2*Te/Ti; //non-continuous!
  ldouble c3 = 18. + 5.*logtratio;
  
  if (Ti < Te) c3 = 18.;
  
  ldouble QpQeRatio;
  if(beta>1.e10) 
    {
      QpQeRatio = c1*sqrt(mratio*tratio);
    }
  else
  QpQeRatio = c1*(c2*c2 + pow(beta,2. - 0.2*logtratio ) )/(c3*c3 + pow(beta, 2. - 0.2*logtratio ) )*sqrt(mratio*tratio)*exp(-1./beta);
 
  ldouble factorVEHeating = 1./(1. + QpQeRatio); //fe from Ressler et al. 2015 following Howes 2010
  delta=factorVEHeating;

  if(!isfinite(delta))
    {
      //printf("problem with delta fit: %d %d f : %e %e %e %e %e\n",geom->ix,geom->iy,delta,QpQeRatio,beta,Te,Ti);
      //print_primitives(pp);
    }

#elif defined(HEATELECTRONS_ROWAN1)
  //Michael's fit to the reconnection electron heating fraction
  //fit for data generated at Te=Ti as a function of sigma_w and beta
  //and mion = mproton
  
  // get beta = gas pressure / magn pressure
  // and sigma = magn energy density / enthalpy
  ldouble beta=0.;
  ldouble sigmaw=0.;
  ldouble betamax=0.;
  ldouble betanorm=0.;

  ldouble rho=state->rho;
  ldouble bsq=state->bsq;
  ldouble Ti = state->Ti;
  ldouble Te = state->Te;
  ldouble ue = state->ue;
  ldouble ui = state->ui;
  ldouble pion = state->pi;
  ldouble gammae = state->gammae;
  ldouble gammai = state->gammai;
    
  ldouble enth_tot = rho + gammai*ui + gammae*ue;

  #ifdef MAGNFIELD
  beta = 2.*pion/bsq;
  sigmaw = bsq/enth_tot;
  betamax = 0.25/sigmaw;

  if(!isfinite(sigmaw)) sigmaw=1.e10; //magnetic field dominates, delta->.5
  if(!isfinite(beta)) beta=1.e10; //no magnetic field, delta->.5
  if(!isfinite(betamax)) betamax=1.e10; //no magnetic field, delta->.4

  betanorm = beta/betamax;
  if(betanorm>1.) betanorm=1.;
  if(betanorm<0.) betanorm=0.;
  #endif

  //Michael's fit to simulation numbers

  delta = fabs(0.5*exp(-pow(1.-betanorm,3.3)/(1.+1.2*pow(sigmaw,0.7))));
  
  if(!isfinite(delta))
    {
      //printf("problem with Michael delta fit: %d %d f : %e %e %e %e %e\n",geom->ix,geom->iy,delta,beta,sigma,Te,Ti);
      //print_primitives(pp);
      delta=.5;
    }

#elif defined(HEATELECTRONS_ROWAN2)
  //Michael's fit to the reconnection electron heating fraction
  //Fit to outflow data
  //fit for data generated at Te=Ti as a function of sigma_w and beta
  //and mion = mproton
  
  // get beta = gas pressure / magn pressure
  // and sigma = magn energy density / enthalpy
  ldouble beta=0.;
  ldouble sigmaw=0.;
  ldouble betamax=0.;
  ldouble betanorm=0.;

  ldouble rho=state->rho;
  ldouble bsq=state->bsq;
  ldouble Ti = state->Ti;
  ldouble Te = state->Te;
  ldouble ue = state->ue;
  ldouble ui = state->ui;
  ldouble pion = state->pi;
  ldouble gammae = state->gammae;
  ldouble gammai = state->gammai;
    
  ldouble enth_tot = rho + gammai*ui + gammae*ue;

  #ifdef MAGNFIELD
  beta = 2.*pion/bsq;
  sigmaw = bsq/enth_tot;

  if(beta<1.e-10) beta=1.e-10;
  if(sigmaw<1.e-10) sigmaw=1.e-10;
  if(!isfinite(sigmaw) || sigmaw>1.e10) sigmaw=1.e10; //magnetic field dominates, delta->.5
  if(!isfinite(beta) || sigmaw>1.e10) beta=1.e10; //no magnetic field, delta->.5
  
  betamax = 0.25/sigmaw;

  if(betamax<1.e-10) betamax=1.e-10;
  if(!isfinite(betamax)) betamax=1.e10; //no magnetic field, delta->.4

  betanorm = beta/betamax;
  if(!isfinite(betanorm) || betanorm>1.) betanorm=1.;
  if(betanorm<0.) betanorm=0.;
  #endif

  //Fit to simulation numbers
  delta = fabs(0.5*exp(-(1.-betanorm)/(0.8+pow(sigmaw,0.5))));
  
  if(!isfinite(delta) || delta>.5)
    {
      //printf("problem with Michael delta fit: %d %d f : %e %e %e %e %e\n",geom->ix,geom->iy,delta,beta,sigma,Te,Ti);
      //print_primitives(pp);
      delta=.5;
    }
  else if(delta<0.)
    delta=0.;
        
#elif defined(HEATELECTRONS_ROWAN3)
  //Michael's fit to the reconnection electron heating fraction with guide field
  //Rowan 2019 eqn 19
  
  // get beta = gas pressure / magn pressure
  // and sigma = magn energy density / enthalpy

  #ifndef GUIDE_RATIO
  #define GUIDE_PREF -0.512 //guide_ratio = 0.33
  #else
  #define GUIDE_PREF 1.7*(tanh(0.33*GUIDE_RATIO)-0.4)
  #endif
  
  ldouble beta=0.;
  ldouble sigmaw=0.;
  ldouble betamax=0.;
  ldouble betanorm=0.;

  ldouble rho=state->rho;
  ldouble bsq=state->bsq;
  ldouble Ti = state->Ti;
  ldouble Te = state->Te;
  ldouble ue = state->ue;
  ldouble ui = state->ui;
  ldouble pion = state->pi;
  ldouble gammae = state->gammae;
  ldouble gammai = state->gammai;
    
  ldouble enth_tot = rho + gammai*ui + gammae*ue;
  ldouble tratio = Te/Ti;
  if(!isfinite(tratio) || tratio>1.e10) tratio=1.e10;
  
  #ifdef MAGNFIELD
  beta = 2.*pion/bsq;
  sigmaw = bsq/enth_tot;

  if(beta<1.e-10) beta=1.e-10;
  if(sigmaw<1.e-10) sigmaw=1.e-10;
  if(!isfinite(sigmaw) || sigmaw>1.e10) sigmaw=1.e10; //magnetic field dominates, delta->.5
  if(!isfinite(beta) || sigmaw>1.e10) beta=1.e10; //no magnetic field, delta->.5
  
  //betamax = 0.25/sigmaw;
  betamax = 0.5/(sigmaw*(1.+tratio));

  if(betamax<1.e-10) betamax=1.e-10;
  if(!isfinite(betamax)) betamax=1.e10; 

  betanorm = beta/betamax;
  if(!isfinite(betanorm) || betanorm>1.) betanorm=1.;
  if(betanorm<0.) betanorm=0.;
  #endif

  //OLD Fit to simulation numbers
  ldouble numer = sqrt(1-betanorm) * (1-betanorm);
  ldouble denom = pow(sigmaw,0.3) * (0.42 + tratio);

  delta = 0.5*(GUIDE_PREF*tanh(numer/denom) + 1.); //GUIDE_PREF = 1.7*(tanh(0.33bg)-0.4)
  
  if(!isfinite(delta))
    {
      //delta=1.;
      delta=0.5; //default to equal heating
    }
  else if(delta<0.)
    delta=0.;

#elif defined(HEATELECTRONS_ZHDANKIN)
  //based on ratio of gyroradi Zhdankin 2019 
  
  ldouble Ti = state->Ti;
  ldouble Te = state->Te;
  ldouble the = kB_over_me * Te;
  ldouble thi = kB_over_mui_mp * Ti;
  //ldouble gammae=GAMMAE;
  //ldouble gammai=GAMMAI;
  //#ifndef FIXEDGAMMASPECIES
  //gammae=calc_gammaintfromtheta(the);
  //gammai=calc_gammaintfromtheta(thi);
  //#endif
  ldouble gmeane = calc_meanlorentz(the);
  ldouble gmeani = calc_meanlorentz(thi);
  ldouble gyroratio = 1836.15267 * sqrt((gmeani*gmeani - 1.)/(gmeane*gmeane - 1.));
  ldouble uioue = pow(gyroratio,2./3.);
  delta = 1./(uioue + 1.);

  if(!isfinite(delta) || delta>1.)
    {
      delta=1.;
    }
  else if(delta<0.)
    delta=0.;
#endif
#endif //HEATELECTRONS

  return delta;
}

//****************************************************************************
//applies change in energy and number to the thermal electrons
//****************************************************************************

ldouble
apply_du_dn_2_species(ldouble *pp, ldouble u, ldouble n,ldouble du, ldouble dn,
		      void* ggg, int type, ldouble* Sreturn)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  int ix=geom->ix;
  int iy=geom->iy;
  int iz=geom->iz;

  //ANDREW only need Tgas for floors below
  ldouble Tg=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);
  
  //Compute the starting temperature and adiabatic index
  ldouble mass, gamma;
  ldouble tmin, tminfrac;
  int  enttype;
  ldouble Snew;
  ldouble mu;
  if (type==IONS) {
    mu = MU_I;
    mass = MU_I*M_PROTON;
    gamma = GAMMAI;
    tmin = TEMPIMINIMAL;
    tminfrac = TEMPIMINIMALFRACTION;
    enttype = ENTRI;
  }
  if (type==ELECTRONS) {
    mu = MU_E;
    mass = M_ELECTR;
    gamma = GAMMAE;
    tmin = TEMPEMINIMAL;
    tminfrac = TEMPEMINIMALFRACTION;
    enttype = ENTRE;
  } 

//Update energy and number directly and compute entropy
#if (ELECTRONIONHEATTYPE==ELECTRONIONHEATTYPE_THROUGHUINT)
    
  u += du;
  #ifdef RELELECTRONS
  n += dn;
  #endif

  ldouble Tnew;
  #ifdef CONSISTENTGAMMA
  Tnew=solve_Teifromnmu(n, mass, u,type); //solve for gamma and Tnew

  //if(type==ELECTRONS)
  if(!isfinite(Tnew))
  {
      int iv;
      PLOOP(iv) printf("%d: %e ",iv,pp[iv]);
      printf("\ndn: %e nfinal: %e\n",dn,n);
      printf("du: %e ufinal: %e\n",du,u);
      printf("mid: not finite Tnew in solve_Tei (%d) at  %d %d %d > %e %e %e %e %e\n",type,ix+TOI,iy+TOJ,iz+TOK,Tnew,n,u,du,dn);
      //getch();
  }
  
  if(type==ELECTRONS && Tnew<TEMPEMINIMAL)
    Tnew=TEMPEMINIMAL;
  if(type==IONS && Tnew<TEMPIMINIMAL)
    Tnew=TEMPIMINIMAL;

  #ifndef FIXEDGAMMASPECIES
  ldouble theta=K_BOLTZ*Tnew/mass;  
  gamma=calc_gammaintfromtheta(theta); //the same gamma, consistent with Tnew    
  #endif
  #endif //CONSISTENTGAMMA
  
  ldouble Tnew1, Tnew2,Tnew3;
  Tnew1 = Tnew;

  ldouble pnew=u*(gamma-1.);
  Tnew2=pnew/(K_BOLTZ*n); 

  //floors on Tnew2
  if(Tnew2<tmin) 
    pnew=n*K_BOLTZ*tmin;
  if(Tnew2<tminfrac*Tg) 
    pnew=n*K_BOLTZ*tminfrac*Tg;
      
  //recompute temperature using floor on pressure
  Tnew3=pnew/(K_BOLTZ*n);
  Tnew=Tnew3;
  
  ldouble rho=n*mu*M_PROTON;
  Snew=calc_SefromrhoT(rho,Tnew,type);


  if(type==ELECTRONS && !isfinite(Snew))
  {
     printf("end: not finite ENTRE at  %d %d %d > %e %e %e %e\n",ix,iy,iz,Tnew,gamma,u,du);
     Snew = pp[ENTRE];
     //getch();
  }
  if(type==IONS && !isfinite(Snew))
  {
     printf("end: not finite ENTRI at  %d %d %d > %e %e %e %e\n",ix,iy,iz,Tnew,gamma,u,du);
     Snew = pp[ENTRI];
     //getch();
  }

  //ANDREW -- systematic bias in slightly overestimating  final energy
  //how can we fix this?
  /*
  if(type==ELECTRONS){
    ldouble u2=calc_ufromS4n(Snew,n,type,ix,iy,iz);
    Tnew=Tnew3*(1+.1*(1.-u2/u)*((rand()%2)*2-1));
    Snew=calc_S4fromnT(n,Tnew,type);
    ldouble u3=calc_ufromS4n(Snew,n,type,ix,iy,iz);
    //printf("%e  %e %e\n",u,u3,1.-u3/u);
  }
  */
#endif //(ELECTRONIONHEATTYPE==ELECTRONIONHEATTYPE_THROUGHUINT)

// or, Update species entropy directly
#if (ELECTRONIONHEATTYPE==ELECTRONIONHEATTYPE_THROUGHENTROPY)

  ldouble Sold=pp[enttype];
  ldouble cpot=0.0;

  //Current temperature
  ldouble Told;
  ldouble rho=n*mu*M_PROTON;
  Told=calc_TfromSerho(pp[enttype],rho,type,ix,iy,iz);

  //calculate chemical potential if there is a change of number density
  #ifdef RELELECTRONS
  ldouble theta= K_BOLTZ*Told/mass;
  cpot = chemical_potential_short(theta, n);
  #endif

  //1st law of thermodynamics
  Snew=Sold + (du - cpot*dn)/Told;

  #ifdef RELELECTRONS
  n += dn;
  #endif
  u += du;
#endif

   ldouble u2=0.;
#ifdef ENFORCE_HEATING_ENERGY_SUM

   u2=calc_ufromSerho(Snew,rho,type,ix,iy,iz);  

#endif
   
   // return the new entropy
   *Sreturn = Snew;

   return u2;
}

ldouble
pick_ViscousHeating(int ix,int iy,int iz)
{
  return get_u_scalar(vischeating,ix,iy,iz);
}



//*********************************************************************************
//tests
//*********************************************************************************


/*! int test_gammagas()
 \brief Test calculation of gamma from temperature
 */
int
test_gammagas()
{
  ldouble Te,Ti,gammagas;
  for(Te=1.;Te<15.;Te*=1.01)
  {
    for(Ti=1.;Ti<15.;Ti*=1.01)
    {
      gammagas=calc_gammaintfromTei(pow(10.,Te),pow(10.,Ti));
      printf("%e %e %e\n",Te,Ti,gammagas);
      
    }
    printf("\n");
  }
  return 0;
}


/*! int test_calcgamma()
 \brief Test calculation of gamma from primitives
 */
int
test_calcgamma()
{
  struct geometry geom;
  fill_geometry(NX/4,0,0,&geom);
  int i;

  ldouble pp[NV];
  ldouble rhocgs=1.e-10;
  ldouble Tgas=1.e10;
  ldouble Te=Tgas;

  pp[RHO]=rhoCGS2GU(rhocgs);
  pp[UU]=calc_PEQ_ufromTrho(Tgas,pp[RHO],NX/4,0,0);

  //nothing
  pp[VX]=pp[VY]=pp[VZ]=0.;
  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],NX/4,0,0);
  
  #ifdef MAGNFIELD
  pp[B1]=pp[B2]=pp[B3]=0.;
  #endif

  #ifdef RADIATION
  //some radiation field
  pp[FX]=pp[FY]=pp[FZ]=0.;
  pp[EE]=pp[UU];
  pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/calc_PEQ_Tfromurho(pp[UU],pp[RHO],NX/4,0,0);
  #endif
  
  #ifdef EVOLVEELECTRONS
  //electron entropy (test)
  ldouble Se=calc_SefromrhoT(pp[RHO],Te,ELECTRONS);
  ldouble S2=calc_S2fromrhoT(pp[RHO],Te,ELECTRONS);
  ldouble S3=calc_S3fromrhoT(pp[RHO],Te,ELECTRONS);

  pp[ENTRE]=Se;

  //ldouble theta= K_BOLTZ*Te/M_ELECTR;
  ldouble theta = kB_over_me * Te;

  ldouble Tee= calc_TfromSerho(pp[ENTRE],pp[RHO],ELECTRONS,NX/4,0,0); //test_calcgamma
  ldouble Te2= calc_TfromS2rho(S2,pp[RHO],ELECTRONS,NX/4,0,0);
  ldouble Te3= calc_TfromS3rho(S3,pp[RHO],ELECTRONS,NX/4,0,0);

  printf("out Te2: %e Te3: %e Tee: %e\n",Te2,Te3,Tee);
  #endif
  
  return 0;
}
