/*! \file p2u.c
 \brief Primitives to conserved conversion
 */


#include "ko.h"


//**********************************************************************
//primitive to conserved converter
//**********************************************************************

int
p2u(ldouble *p, ldouble *u, void *ggg)
{
  p2u_mhd(p,u,ggg);

#ifdef RADIATION
  p2u_rad(p,u,ggg);
#endif

  return 0;
}


//**********************************************************************
//primitive to conserved converter -- hydro
//**********************************************************************


int
p2u_mhd(ldouble *p, ldouble *u, void *ggg)
{

#ifdef NONRELMHD
  p2u_mhd_nonrel(p,u,ggg);
  return 0;
#endif

  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  gdet=geom->gdet;
  GG=geom->GG;
  gdetu=gdet;

  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif

  ldouble rho=p[0];
  ldouble uu=p[1];
  ldouble vcon[4],vcov[4],ucon[4],ucov[4];
  ldouble bcon[4]={0.,0.,0.,0.}, bcov[4]={0.,0.,0.,0.}, bsq=0.;
  vcon[0]=0.;
  vcon[1]=p[2];
  vcon[2]=p[3];
  vcon[3]=p[4];
  ldouble S=p[5];

  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel(p, ucon, ucov, geom, bcon, bcov, &bsq);
#endif

  //************************************
  //hydro part
  //************************************
 
  ldouble ut=ucon[0];
  ldouble rhout = rho*ut;
  ldouble Sut = S*ut;

  //ANDREW Do we need to ensure gamma is consistent here?
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(geom->ix,geom->iy,geom->iz);
  #endif
  ldouble gammam1=gamma-1.;

  ldouble pre=(gamma-1.)*uu; 
  ldouble w=rho+uu+pre;
  ldouble eta=w+bsq;
  ldouble etap = uu+pre+bsq; //eta-rho
  ldouble ptot=pre+0.5*bsq;
  
  //this computes utp1=1+u_t
  ldouble utp1;
  utp1=calc_utp1(vcon,ucon,geom);

  ldouble Tttt =etap*ucon[0]*ucov[0] + rho*ucon[0]*utp1 + ptot - bcon[0]*bcov[0];
  ldouble Ttr  =eta*ucon[0]*ucov[1] - bcon[0]*bcov[1];
  ldouble Ttth =eta*ucon[0]*ucov[2] - bcon[0]*bcov[2];
  ldouble Ttph =eta*ucon[0]*ucov[3] - bcon[0]*bcov[3];

  u[0]=gdetu*rhout;
  u[1]=gdetu*Tttt;
  u[2]=gdetu*Ttr;
  u[3]=gdetu*Ttth;
  u[4]=gdetu*Ttph;
  u[5]=gdetu*Sut;

#ifdef EVOLVEELECTRONS
  ldouble Seut=p[ENTRE]*ut;
  u[ENTRE]= gdetu*Seut;
  ldouble Siut=p[ENTRI]*ut;
  u[ENTRI]= gdetu*Siut;
#endif

#ifdef RELELECTRONS
  int ib;
  for(ib=0;ib<NRELBIN;ib++)
    {
      u[NEREL(ib)]=gdetu*p[NEREL(ib)]*ut;
    }    
#endif

  //************************************
  //magnetic part
  //************************************ 
#ifdef MAGNFIELD
  u[B1]=gdetu*p[B1];
  u[B2]=gdetu*p[B2];
  u[B3]=gdetu*p[B3];
#endif

  return 0.;
}


//**********************************************************************
//this computes utp1=1+u_t , which for nonrelativistic cases is ~0.
//If computed as 1+u_t, then if residual is small there will be a large error.
//**********************************************************************

ldouble
calc_utp1(ldouble *vcon, ldouble *ucon, void *ggg)
{
   struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  gdet=geom->gdet;
  GG=geom->GG;
  gdetu=gdet;

  ldouble utp1;
  if(VELPRIM==VELR) //based on VELR
  {
      int i,j;
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=vcon[i]*vcon[j]*gg[i][j];
      ldouble gamma2=(1.+qsq);
      ldouble alpha = geom->alpha;
      ldouble alphasq = alpha*alpha;
      ldouble alpgam=sqrt(alphasq*gamma2);

      //\beta^i \beta_i / \alpha^2 = g^{ti} g_{ti}
      ldouble betasqoalphasq=gg[0][1]*GG[0][1] + gg[0][2]*GG[0][2] + gg[0][3]*GG[0][3];
      
      ldouble ud0tilde = 0.0;
      SLOOPA(j) ud0tilde += vcon[j]*gg[0][j]; // \tilde{u}_t = \tilde{u}^i g_{ti} since \tilde{u}^t=0
      utp1= ud0tilde + (geom->gttpert - alphasq*(betasqoalphasq + qsq))/(1.0+alpgam);
  }
  else //based on ucon[]
  {
      int i,j,k;

      // 3-velocity in coordinate basis
      ldouble vconp[4];
      SLOOPA(j) vconp[j]=ucon[j]/ucon[0];

      ldouble plus1gv00=geom->gttpert;
      ldouble vsq=geom->gttpert;
      SLOOPA(j) vsq+=2.0*geom->gg[0][j]*vconp[j];
      SLOOP(j,k) vsq+=geom->gg[j][k]*vconp[j]*vconp[k];

      ldouble gvtt=geom->gg[0][0];
      
      ldouble alpha=0.0;
      SLOOPA(j) alpha+=geom->gg[j][0]*ucon[j];

      ldouble uu0 = ucon[0];

      utp1 = alpha + ((1.0-gvtt)*plus1gv00 - uu0*uu0*vsq*gvtt*gvtt)/(1.0-gvtt*uu0);
  }
  return utp1;
}


//**********************************************************************
//primitive to conserved converter - non-relativistic!
//**********************************************************************

int
p2u_mhd_nonrel(ldouble *p, ldouble *u, void *ggg)
{

  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  gdet=geom->gdet;
  GG=geom->GG;
  gdetu=gdet;

  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif

  ldouble rho=p[0];
  ldouble uu=p[1];
  ldouble vcon[4],vcov[4],ucon[4],ucov[4];
  ldouble bcon[4]={0.,0.,0.,0.},bcov[4]={0.,0.,0.,0.},bsq=0.;
  vcon[1]=p[2];
  vcon[2]=p[3];
  vcon[3]=p[4];
  vcon[0]=0.;
  ldouble S=p[5];

  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

  #ifdef NONRELMHD
  ucon[0]=1.;
  ucov[0]=-1.;
  #endif

  ldouble v2=dot3nr(ucon,ucov);

#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel(p, ucon, ucov, geom, bcon, bcov, &bsq);
#endif

  //************************************
  //hydro part
  //************************************
 
  ldouble Ttt=-(uu + bsq/2. + rho*v2/2.);
  ldouble Tttt=Ttt;
  ldouble Ttr =rho*ucov[1];
  ldouble Ttth =rho*ucov[2];
  ldouble Ttph =rho*ucov[3];

  u[0]=gdetu*rho;
  u[1]=gdetu*Tttt;
  u[2]=gdetu*Ttr;
  u[3]=gdetu*Ttth;
  u[4]=gdetu*Ttph;
  u[5]=gdetu*S;

  //************************************
  //magnetic part
  //************************************
 
#ifdef MAGNFIELD
  u[B1]=gdetu*p[B1];
  u[B2]=gdetu*p[B2];
  u[B3]=gdetu*p[B3];
#endif

  return 0.;
}

/********************************************************/
/**** converts radiative primitives *********************/
/********************************************************/

int p2u_rad(ldouble *pp, ldouble *uu, void *ggg)
{
  int i, j, irf;
  
  struct geometry *geom
  = (struct geometry *) ggg;
  
  ldouble (*gg)[5], (*GG)[5],  gdet, gdetu;
  gg = geom->gg;
  GG = geom->GG;
  gdet = geom->gdet;
  gdetu = gdet;
#if (GDETIN == 0) //gdet out of derivatives
  gdetu = 1.;
#endif
  
  //Only M1 closure supported!
  {
    ldouble Erf = pp[EE];
    
    //relative four-velocity
    ldouble urf[4];
    urf[0] = 0.;
    urf[1] = pp[FX];
    urf[2] = pp[FY];
    urf[3] = pp[FZ];
    
    //converting to lab four-velocity
    conv_vels(urf, urf, VELPRIMRAD, VEL4, gg, GG);
    
    ldouble Rtopp[4];    
    for (i = 0; i < 4; i++)  // write like this to help vectorization
    {
      Rtopp[i] = four_third * Erf * urf[0] * urf[i] + one_third * Erf * GG[0][i];
    }
    
    indices_21(Rtopp, Rtopp, gg); //R^t_mu
    
    for (i = 0; i < 4; i++)  // to help vectorization
    {
      uu[EE+i] = gdetu * Rtopp[i];
    }
    
#ifdef EVOLVEPHOTONNUMBER
    uu[NF] = gdetu * pp[NF] * urf[0];
#endif
  }
  
  return 0;
}

//**********************************************************************
//takes primitives and calculates quantities that are averaged for the avg file
//**********************************************************************

int
p2avg(int ix,int iy,int iz,ldouble *avg)
{

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  struct geometry geomout;
  fill_geometry_arb(ix,iy,iz,&geomout,OUTCOORDS);

  // ANDREW now that we have set avg[AVGRHOURDIFF] = 0, we don't need these
  /*
  struct geometry geoml;
  fill_geometry_face(ix,iy,iz,0,&geoml);
  struct geometry geomoutl;
  fill_geometry_face_arb(ix,iy,iz,0,&geomoutl,OUTCOORDS);
  */
 
  int iv,iv2;ldouble pp[NV],uu[NV];
  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz); //conserved 
      pp[iv]=get_u(p,iv,ix,iy,iz); //primitives 
      avg[iv]=pp[iv]; //first NV slots in pavg are regular primitives in MYCOORDS!
   }

  for(iv=NV;iv<NAVGVARS;iv++)
    avg[iv]=0.;

  //primitives to OUTCOORDS
#ifdef PRECOMPUTE_MY2OUT
  trans_pall_coco_my2out(pp,pp,&geom,&geomout);
#else      
  trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomout);
#endif
  
  ldouble (*gg)[5],(*GG)[5];
  gg=geomout.gg;
  GG=geomout.GG;

  //four-vectors etc
  ldouble rho=pp[0];
  ldouble uint=pp[1];
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(ix,iy,iz);
  #endif
  ldouble gammam1=gamma-1.;

  ldouble pregas=(gamma-1.)*uint;
 
  ldouble Ti,Te,trad;
  ldouble Tgas=calc_PEQ_Teifrompp(pp,&Te,&Ti,ix,iy,iz);
  ldouble vcon[4],vcov[4],ucon[4],ucov[4];
  ldouble bcon[4]={0.,0.,0.,0.},bcov[4]={0.,0.,0.,0.},bsq=0.;
  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];
  vcon[0]=0.;
  ldouble S=pp[5];
  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG); 

#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, &geomout, bcon, bcov, &bsq);
#endif
  
  //hydro stress-energy
  ldouble Tij[4][4];
  calc_Tij(pp,&geomout,Tij);
  indices_2221(Tij,Tij,gg);

  //avg already in OUTCOORDS
  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGTIJ(iv,iv2)]=Tij[iv][iv2];
  avg[AVGBSQ]=bsq;
  avg[AVGPGAS]=pregas;
  avg[AVGTGAS]=Tgas;
  for(iv=0;iv<4;iv++)
    avg[AVGUCON(iv)]=ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGUCOV(iv)]=ucov[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGBCON(iv)]=bcon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGBCOV(iv)]=bcov[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGRHOUCON(iv)]=rho*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGRHOUCOV(iv)]=rho*ucov[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGUUUCON(iv)]=uint*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGUUCOV(iv)]=uint*ucov[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGBSQUCON(iv)]=bsq*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGBSQUCOV(iv)]=bsq*ucov[iv];
  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGRHOUCONUCOV(iv,iv2)]=rho*ucon[iv]*ucov[iv2];
  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGUUUCONUCOV(iv,iv2)]=uint*ucon[iv]*ucov[iv2];
  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGBSQUCONUCOV(iv,iv2)]=bsq*ucon[iv]*ucov[iv2];
  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGBCONBCOV(iv,iv2)]=bcon[iv]*bcov[iv2];
  for(iv=0;iv<4;iv++)
    avg[AVGWUCON(iv)]=(rho+uint+bsq/2)*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGWUCON(iv)]=(rho+uint+bsq/2)*ucon[iv];

  //fluxes at faces, including the diffusive part
  //conserved fluxes at left face in MYCOORDS
  for(iv=0;iv<NV;iv++)
    {
      avg[AVGFLUXXL(iv)]=get_ub(flbx,iv,ix,iy,iz,0);
      avg[AVGFLUXYL(iv)]=get_ub(flby,iv,ix,iy,iz,1);
      avg[AVGFLUXZL(iv)]=get_ub(flbz,iv,ix,iy,iz,2);
    }

  // ANDREW - unclear why this is necessary
  // and it has a bug! (iv != 0 in vector[1] = get_ub)
  avg[AVGRHOURDIFF] = 0;
/*
  //converting rest-mass flux to BLCOORDS 
  ldouble vector[4];

  //primitives and conserved at left faces - used to fill missing time-component
  ldouble uface[NV],pface[NV],fd_uLl[NV],fd_uRl[NV];
  int i;
  for(i=0;i<NV;i++)
    {
      fd_uLl[i]=get_ub(pbLx,i,ix,iy,iz,0);
      fd_uRl[i]=get_ub(pbRx,i,ix,iy,iz,0);
      pface[i]=.5*(fd_uLl[i]+fd_uRl[i]);
    }
  p2u(pface,uface,&geoml);

  //rest-mass flux in radius
  vector[0]=uface[RHO]; //rho ut gdet
  vector[1]=get_ub(flbx,iv,ix,iy,iz,0); //rho ur gdet
  vector[2]=0.; //unimportant within Kerr-Shield
  vector[3]=0.;
  trans2_coco(geoml.xxvec,vector,vector,MYCOORDS, OUTCOORDS); //this does not work?
  avg[AVGRHOURDIFF]=vector[1];
*/
    
#ifdef EVOLVEELECTRONS
  //electrons
  ldouble ne=calc_thermal_ne(pp);
  ldouble pree=K_BOLTZ*ne*Te;

  //ions
  ldouble ni=rho/MU_I/M_PROTON; 
  ldouble prei=K_BOLTZ*ni*Ti;
  avg[AVGPE]=pree;
  avg[AVGPI]=prei;

#ifdef RELELECTRONS
  ldouble nerelel=calc_relel_ne(pp);
  ldouble prelel=calc_relel_p(pp);
  ldouble uintrelel=calc_relel_uint(pp);
  avg[AVGNRELEL]=nerelel;
  avg[AVGPRELEL]=prelel;
  avg[AVGURELEL]=uintrelel;
#endif
#endif
   
#ifdef RADIATION
  ldouble Rtt,Ehat,ugas[4];

  ldouble Rij[4][4];
  calc_Rij(pp,&geomout,Rij);
  indices_2221(Rij,Rij,geomout.gg);

  //Ehat calculation from Rij, gas velocity still in ucon, ucov
  Rtt=0.;
  int i1,i2;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
    Rtt+=-Rij[i1][i2]*ucon[i2]*ucov[i1];

  Ehat=-Rtt;       

  vcon[1]=pp[FX];
  vcon[2]=pp[FY];
  vcon[3]=pp[FZ];
  vcon[0]=0.;  
  ldouble urcon[4],urcov[4];
  conv_vels_both(vcon,urcon,urcov,VELPRIM,VEL4,gg,GG); 

  //radiation four-fource
  ldouble Gi[4],Giff[4]={0.,0.,0.,0.};
  ldouble Gic[4],Gicff[4]={0.,0.,0.,0.};
  
  calc_Gi(pp,&geomout,Giff,0.,0,0); //fluid frame 

#if defined(COMPTONIZATION)
  //directly in ff
  ldouble uconff[4];
  uconff[1]=uconff[2]=uconff[3]=0.;
  uconff[0]=1.;
  ldouble kappaes=calc_kappaes(pp,&geomout);
  calc_Compt_Gi(pp,&geomout,Gicff,Ehat,Te,kappaes,uconff);
#endif 

  //radiation temperature
  ldouble Thatrad,nphhat,ThatradBB;
  ThatradBB=calc_LTE_TfromE(Ehat);
  Thatrad=ThatradBB;
  #ifdef EVOLVEPHOTONNUMBER //number of photons conserved
  Thatrad = calc_ncompt_Thatrad(pp,&geomout,Ehat);  

  //Trad Limiter
#ifdef MAXDIFFTRADS
  ldouble maxfac=MAXDIFFTRADS;
  
  //Ceiling: Ramesh: we adjust the photon number when changing the temperature
  if (Thatrad > maxfac * ThatradBB)
  {
    pp[NF] *= (Thatrad / (maxfac * ThatradBB));
    Thatrad = maxfac * ThatradBB;
  }
  
  // Floor: Ramesh: different floors, depending on whether or not we have synchrotron. Again, adjust photon number
  #ifndef SYNCHROTRON
  if (Thatrad < ThatradBB)
  {
    pp[NF] *= (Thatrad / ThatradBB);
    Thatrad = ThatradBB;
  }
  #else
  if (Thatrad < ThatradBB / maxfac)
  {
    pp[NF] *= (Thatrad * maxfac / ThatradBB);
    Thatrad = ThatradBB / maxfac;
  }
  #endif
#endif //MAXDIFFTRADS

  avg[AVGNFHAT]=calc_ncompt_nphhat(pp, &geomout);

#endif //EVOLVEPHOTONNUMBER
  for(iv=0;iv<4;iv++)
    avg[AVGURFCON(iv)]=urcon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGURFCOV(iv)]=urcov[iv];
   
  avg[AVGEHAT]=Ehat;


  for(iv=0;iv<4;iv++)
    for(iv2=0;iv2<4;iv2++)
      avg[AVGRIJ(iv,iv2)]=Rij[iv][iv2];

  for(iv=0;iv<4;iv++)
    avg[AVGEHATUCON(iv)]=Ehat*ucon[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGEHATUCOV(iv)]=Ehat*ucov[iv];

  for(iv=0;iv<4;iv++)
    avg[AVGGHAT(iv)]=Giff[iv];
  for(iv=0;iv<4;iv++)
    avg[AVGGHATCOMPT(iv)]=Gicff[iv];

#endif //RADIATION
if(doingpostproc_avg)
{
  struct struct_of_state state;
  fill_struct_of_state(pp,&geomout,&state);
  ldouble fe=calc_ViscousElectronHeatingFraction_from_state(pp,&state,&geomout);
  
  avg[AVGVISCHEATING]=1.;
  avg[AVGVISCHEATINGNEGE]=0.;
  avg[AVGVISCHEATINGNEGI]=0.;
  avg[AVGVISCHEATINGTIMESDELTAE]=fe;
}
else
{
  avg[AVGVISCHEATING]=get_u_scalar(vischeating,ix,iy,iz);
  avg[AVGVISCHEATINGNEGE]=get_u_scalar(vischeatingnegebalance,ix,iy,iz);
  avg[AVGVISCHEATINGNEGI]=get_u_scalar(vischeatingnegibalance,ix,iy,iz);
  avg[AVGVISCHEATINGTIMESDELTAE]=get_u_scalar(vischeatingtimesdeltae,ix,iy,iz);
}
 
  return 0.;
}

//**********************************************************************
/*! int test_maginv()
 \brief Test p2u and u2p for MHD fluid
*/
//**********************************************************************

int
test_maginv()
{
  ldouble uu[NV];
  ldouble pp[NV];
  
  //geometries
  struct geometry geom;
  fill_geometry(0,0,0,&geom);
  
  print_metric(geom.gg);
  
  pp[RHO]=1.;
  pp[UU]=0.001;
  pp[VX]=pp[VY]=pp[VZ]=0.;
  
  pp[VX]=0.0;
  pp[VY]=0.0;
  pp[VZ]=0.0;
  
#ifdef MAGNFIELD
  pp[B1]=pp[B2]=pp[B3]=0.;
  
  pp[B1]=0.e-2;
  pp[B2]=0.e-4;
  pp[B3]=1.e-0;
#endif
  
  //entropy
  pp[5]=calc_Sfromu(pp[RHO],pp[UU],0,0,0);
  
  print_NVvector(pp);
  p2u(pp,uu,&geom);
  print_NVvector(uu);
  
  int aa[3],bb[2],ret;
  pp[UU]*=1.1;
  ret=u2p(uu,pp,&geom,aa,bb,0);
  printf("u2p ret: %d\n",ret);
  print_NVvector(pp);
  
  return 0;
}
