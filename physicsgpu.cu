extern "C" {

#include "ko.h"

}

#include "kogpu.h"

//**********************************************************************//
// calculate total gas entropy from density & energy density
//**********************************************************************//
__device__ __host__ ldouble calc_Sfromu_device(ldouble rho,ldouble u,int ix,int iy,int iz)
{
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  //gamma=pick_gammagas(ix,iy,iz); //TODO
  #endif
  ldouble gammam1=gamma-1.;
  ldouble indexn=1.0/gammam1;
  ldouble pre=gammam1*u;
  #ifdef NOLOGINS
  ldouble ret = rho*u / pow(rho,gamma);
  #else
  ldouble ret = rho*log(pow(pre,indexn)/pow(rho,indexn+1.));
  #endif

  return ret;
}


//**********************************************************************//
// calculate stress energy tensor
//**********************************************************************//
__device__ __host__ int calc_Tij_device(ldouble *pp, void* ggg, ldouble T[][4])
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble utcon[4],ucon[4],ucov[4];  
  ldouble bcon[4],bcov[4],bsq=0.;
  
  //converts to 4-velocity
  utcon[0]=0.;
  for(int iv=1;iv<4;iv++)
    utcon[iv]=pp[1+iv];
  conv_vels_both_device(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

#ifdef NONRELMHD
  ucon[0]=1.;
  ucov[0]=-1.;
#endif

#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel_device(pp, ucon, ucov, geom, bcon, bcov, &bsq); 
#else
  bcon[0]=bcon[1]=bcon[2]=bcon[3]=0.;
  bsq=0.;
#endif
  
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  //gamma=pick_gammagas(geom->ix,geom->iy,geom->iz); //TODO
  #endif

  ldouble rho=pp[RHO];
  ldouble uu=pp[UU];  
  ldouble p=(gamma-1.)*uu; 
  ldouble w=rho+uu+p;
  ldouble eta=w+bsq;
  ldouble ptot=p+0.5*bsq;

#ifndef NONRELMHD  
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      T[i][j]=eta*ucon[i]*ucon[j] + ptot*GG[i][j] - bcon[i]*bcon[j];
#else
  
  ldouble v2=dot3nr(ucon,ucov); //TODO
  for(int i=1;i<4;i++)
    for(int j=1;j<4;j++)
      T[i][j]=(rho)*ucon[i]*ucon[j] + ptot*GG[i][j] - bcon[i]*bcon[j];

  T[0][0]=uu + bsq/2. + rho*v2/2.;
  for(int i=1;i<4;i++)
    T[0][i]=T[i][0]=(T[0][0] + ptot) *ucon[i]*ucon[0] + ptot*GG[i][0] - bcon[i]*bcon[0];

#endif  // ifndef NONRELMHD

  return 0;
}

//***************************************************************//
// calculates fluxes at faces
//***************************************************************//
__device__ int f_flux_prime_device(ldouble *pp, int idim, ldouble *ff, void* ggg)
{  

  struct geometry *geom
    = (struct geometry *) ggg;

  // zero fluxes initially
  for(int iv=0;iv<NV;iv++)
  {
    ff[iv]=0.;
  }
  
  ldouble (*gg)[5],(*GG)[5],gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdetu=geom->gdet;

  #if (GDETIN==0) //no metric determinant inside derivative
  gdetu=1.;
  #endif

  //calculating Tij
  ldouble T[4][4];
  calc_Tij_device(pp,geom,T);
  indices_2221_device(T,T,gg);//T^ij --> T^i_j

  //primitives
#ifdef EVOLVEELECTRONS
  ldouble Se=pp[ENTRE]; //entropy of electrons
  ldouble Si=pp[ENTRI]; //entropy of ions
#endif
  ldouble rho=pp[RHO];
  ldouble u=pp[UU];
  ldouble S=pp[5];
  
  ldouble vcon[4],ucon[4],ucov[4];
  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];

  //converting to 4-velocity
  conv_vels_both_device(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

#ifdef NONRELMHD
  ucon[0]=1.;
  ucov[0]=-1.;
#endif

  ldouble bsq=0.;
#ifdef MAGNFIELD
  ldouble bcon[4],bcov[4];
  calc_bcon_bcov_bsq_from_4vel_device(pp, ucon, ucov, geom, bcon, bcov, &bsq);
#endif

  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  //gamma=pick_gammagas(geom->ix,geom->iy,geom->iz); // TODO 
  #endif

  ldouble pre=(gamma-1.)*u; 
  ldouble w=rho+u+pre;
  ldouble eta=w+bsq;
  ldouble etap = u+pre+bsq; //eta-rho

  for(int ii=0;ii<4;ii++)
  {
    for(int jj=0;jj<4;jj++)
    {
	if(isnan(T[ii][jj])) 
	{
	  printf("%d %d %e\n",ii,jj,T[ii][jj]);
	  printf("nan in flux_prime \n");
	  // TODO print and my_err
	  //printf("%d > nan tmunu: %d %d %e at %d %d %d\n",PROCID,ii,jj,T[ii][jj],geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK);
	  //  printf("%d > nan tmunu: %e %e %e %e\n",PROCID,gamma,pre,w,eta);
	  //  print_4vector(ucon);
	  //  print_metric(geom->gg);
	  //  print_Nvector(pp,NV);
	  //  my_err("nan in flux_prime\n");
	  //  exit(1);
	}
    }
  }
  
  ldouble utp1=calc_utp1_device(vcon,ucon,geom);

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

  //magnetic fluxes
#ifdef MAGNFIELD
  ff[B1]=gdetu * (bcon[1]*ucon[idim+1] - bcon[idim+1]*ucon[1]);
  ff[B2]=gdetu * (bcon[2]*ucon[idim+1] - bcon[idim+1]*ucon[2]);
  ff[B3]=gdetu * (bcon[3]*ucon[idim+1] - bcon[idim+1]*ucon[3]);

  //TODO TODO -- what to do about fused multiply add that prevents exact cancellation above?
  // turned off in compiler with -fmad=false, but in general more precision/perfomance from fmad is supposed to be good
  // what to do? 
  /*
  if(geom->ix==ixTEST && geom->iy==iyTEST && geom->iz==izTEST && idim==0)
  {
    printf("in flux_prime, idim %d: ff[B1] %e \n",idim,ff[B1]);
  }
  */
#endif

  //radiation fluxes // TODO 
#ifdef RADIATION
  f_flux_prime_rad(pp,idim,geom,ff);
#ifdef MAGNFIELD
#ifdef BATTERY
  ldouble eterm[4]={-1.,0.,0.,0.};
  calc_batteryflux(pp,geom,eterm,idim,ucov);
  ldouble suppfac=1.;
 
  ff[B1]+=suppfac*gdetu*eterm[1];
  ff[B2]+=suppfac*gdetu*eterm[2];
  ff[B3]+=suppfac*gdetu*eterm[3];
#endif
#endif
#endif

  // electron species fluxes
#ifdef EVOLVEELECTRONS
  ff[ENTRE]= gdetu*Se*ucon[idim+1]; 
  ff[ENTRI]= gdetu*Si*ucon[idim+1]; 

#ifdef RELELECTRONS
  for(int ie=0; ie < NRELBIN ; ie++)
    ff[NEREL(ie)] = gdetu*pp[NEREL(ie)]*ucon[idim+1];
#endif
#endif

      if(geom->ix==ixTEST && geom->iy==iyTEST && geom->iz==izTEST && idim==0)
      {
	//	printf("in flux_prime, idim %d: %e %e %e %e %e %e %e %e %e\n",idim,ff[0],ff[1],ff[2],ff[3],ff[4],ff[5],ff[6],ff[7],ff[8]);
      }

  return 0;
}


//***************************************************************//
// calculates metric source term at centers
//***************************************************************//

// TODO: deleted RADIATION and SHEARINGBOX parts
__device__ __host__ int f_metric_source_term_device(ldouble *pp, ldouble *ss, void* ggg,ldouble* gKr_arr)
{
  struct geometry *geom 
    = (struct geometry *) ggg;

  ldouble (*gg)[5],gdetu;
  gg=geom->gg;
  gdetu=geom->gdet;

  #if (GDETIN==0) //no metric determinant inside derivative
  gdetu=1.;
  #endif

  int ix,iy,iz;
  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;
 
  ldouble T[4][4];
  //calculating stress energy tensor components
  //calc_Tij_device(pp,&geom,T);
  calc_Tij_device(pp,geom,T); 
  indices_2221_device(T,T,gg);

  for(int i=0;i<4;i++)
  {
    for(int j=0;j<4;j++)
    {
	if(isnan(T[i][j])) 
	{
	    printf("%d %d %e\n",i,j,T[i][j]);
	    printf("nan in metric_source_terms\n");
	    //my_err("nan in metric_source_terms\n");//TODO
	}
    }
  }
  
  // zero out all source terms initially
  for(int iv=0;iv<NV;iv++)
  {
    ss[iv]=0.;  
  }
  //terms with Christoffels
  for(int k=0;k<4;k++)
  {
    for(int l=0;l<4;l++)
    {
      ss[1]+=gdetu*T[k][l]*get_gKr(gKr_arr,l,0,k,ix,iy,iz);
      ss[2]+=gdetu*T[k][l]*get_gKr(gKr_arr,l,1,k,ix,iy,iz);
      ss[3]+=gdetu*T[k][l]*get_gKr(gKr_arr,l,2,k,ix,iy,iz);
      ss[4]+=gdetu*T[k][l]*get_gKr(gKr_arr,l,3,k,ix,iy,iz);       
    }
  }


#if (GDETIN==0)
  ldouble (*GG)[5];
  GG=ggg->GG;
    
  //gdet derivatives
  ldouble dlgdet[3];
  dlgdet[0]=gg[0][4]; //D[gdet,x1]/gdet
  dlgdet[1]=gg[1][4]; //D[gdet,x2]/gdet
  dlgdet[2]=gg[2][4]; //D[gdet,x3]/gdet

  //get 4-velocity
  ldouble vcon[4],ucon[4];
  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];  
  conv_vels_device(vcon,ucon,VELPRIM,VEL4,gg,GG); 

  //terms with dloggdet  
  for(int l=1;l<4;l++)
  {
    ss[0]+=-dlgdet[l-1]*pp[RHO]*ucon[l];
    ss[1]+=-dlgdet[l-1]*(T[l][0]+pp[RHO]*ucon[l]);
    ss[2]+=-dlgdet[l-1]*(T[l][1]);
    ss[3]+=-dlgdet[l-1]*(T[l][2]);
    ss[4]+=-dlgdet[l-1]*(T[l][3]);
    ss[5]+=-dlgdet[l-1]*pp[ENTR]*ucon[l];
  }   
#endif
  
  return 0;
}


//************************************************************************//
//calculates left and right wave speeds at cell center
//************************************************************************//

__device__ __host__ int calc_wavespeeds_lr_pure_device(ldouble *pp,void *ggg,ldouble *aaa)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  
  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
 
  ldouble axhdl,axhdr,ayhdl,ayhdr,azhdl,azhdr;
  ldouble axl0,axr0,ayl0,ayr0,azl0,azr0;
  ldouble axl,axr,ayl,ayr,azl,azr;
  axl0=axr0=ayl0=ayr0=azl0=azr0=1.;
  axl=axr=ayl=ayr=azl=azr=1.;
 
  //**********************************************************************//
  //***** four velocity **************************************************//
  //**********************************************************************//
  ldouble utcon[4],ucon[4],ucov[4];
  for(int iv=1;iv<4;iv++)
    utcon[iv]=pp[1+iv];
  utcon[0]=0.;
  conv_vels_both_device(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

  //**********************************************************************//
  //***** hydro: speed of sound ******************************************//
  //**********************************************************************//
	      
  ldouble rho=pp[RHO];
  ldouble uu=pp[UU];
 
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA 
  //gamma=pick_gammagas(geom->ix,geom->iy,geom->iz); //TODO
  #endif

  //gas pressure
  ldouble pre=(gamma-1.)*uu;
  ldouble preeff=pre;
  ldouble uueff=uu;
  
#ifdef GASRADCOUPLEDWAVESPEEDS
#ifdef RADIATION  // gas-radiation coupling is needed only with radiation

  //radiation pressure
  ldouble temp[4],Ehat;
  calc_ff_Ehat_device(pp,&Ehat,temp,geom); // TODO
  ldouble prerad=one_third*Ehat;

  //coupling coefficient
  ldouble fcouple;  
  // Choose which version of gas-radiation coupling to use
#ifdef GASRADCOUPLEDWAVESPEEDS_SIMPLE
  // Simple scheme just uses the total pressure to estimate hydro wave speed
  fcouple = 1.;
  preeff += prerad;
  uueff += Ehat;
#else
  // More sophisticated method, which estimates the degree of coupling between gas and radiation
  // and estimates the effective pressure for hydro wave speed
  fcouple = estimate_gas_radiation_coupling_device(pp,ggg); // TODO
  if (isnan(fcouple))
    fcouple = 1.; 
  preeff += fcouple*prerad;
  uueff += fcouple*Ehat;
#endif  
#endif  
#endif  // ifdef GASRADCOUPLEDWAVESPEEDS

  //Sound speed
  ldouble cs2=gamma*preeff/(rho+uueff+preeff);

  //maximum is 95% light
  if(cs2>0.95) cs2=0.95; 

  //********************************************************************//
  //***** magn: alfvenic speed *****************************************//
  //********************************************************************//
	      
  ldouble va2=0.;
#ifdef MAGNFIELD
  ldouble bcon[4],bcov[4],bsq;
  calc_bcon_bcov_bsq_from_4vel_device(pp, ucon, ucov, geom, bcon, bcov, &bsq);
  ldouble EF = rho + uu + pre;
  ldouble EEloc = bsq + EF ;
  va2 = bsq/EEloc ;
  if(va2<0.) va2=0.;
#endif

  //********************************************************************//
  //***** mhd: fast magnetosonic speed *********************************//
  //********************************************************************//

#ifdef NONRELMHD // non-rel version
  
  ldouble vx=pp[VX];
  ldouble vy=pp[VY];
  ldouble vz=pp[VZ];
  ldouble cs=sqrt(vtot2);
  ldouble csx=cs/sqrt(gg[1][1]);
  ldouble csy=cs/sqrt(gg[2][2]);
  ldouble csz=cs/sqrt(gg[3][3]);
  
  axhdr=vx+csx;
  axhdl=vx-csx;

  ayhdr=vy+csy;
  ayhdl=vy-csy;
  
  azhdr=vz+csz;
  azhdl=vz-csz;

#else // relativistic version

  ldouble vtot2=cs2 + va2 - cs2*va2; //total characteristic velocity
  
  //********************************************************************//
  //algorithm from HARM to transform the fluid fr wavespeed into lab fr
  //********************************************************************//

  ldouble aret[6];
  int ret;
  ret = calc_wavespeeds_lr_core_device(ucon,GG,aret,vtot2,vtot2,vtot2);
  if(ret<0)
  {
    printf("error in calc_wavespeeds_lr_core_device at %d | %d | %d\n",geom->ix,geom->iy,geom->iz);
  }
  
  axhdl=aret[0];
  axhdr=aret[1];
  ayhdl=aret[2];
  ayhdr=aret[3];
  azhdl=aret[4];
  azhdr=aret[5];

#endif // ifdef NORELMHD
 
#ifdef RADIATION
  
  //********************************************************************//
  //***** radiation: characteristic wave speed *************************//
  //********************************************************************//

  // physical size of the cell
  // ix,iy,iz could be the indices of a face, so the optical depth taken from left/right
  //TODO pass the xb_arr to get_size_x
  ldouble dx[3];
  dx[0]=my_max(get_size_x(geom->ix,0)*sqrt(gg[1][1]),get_size_x(geom->ix+1,0)*sqrt(gg[1][1]));
  dx[1]=my_max(get_size_x(geom->iy,1)*sqrt(gg[2][2]),get_size_x(geom->iy+1,1)*sqrt(gg[2][2]));
  dx[2]=my_max(get_size_x(geom->iz,2)*sqrt(gg[3][3]),get_size_x(geom->iz+1,2)*sqrt(gg[3][3]));

  // calculate optical depth
  ldouble tautot[3];
  calc_tautot_device(pp,geom,dx,tautot); // TODO

  // compute radiative wavespeeds
  ldouble aval[18];
  int verbose=0;
  calc_rad_wavespeeds_device(pp,geom,tautot,aval,verbose); // TODO 

  // wavespeeds unlimited by optical depth
  axl0=aval[0];
  axr0=aval[1];
  ayl0=aval[2];
  ayr0=aval[3];
  azl0=aval[4];
  azr0=aval[5];

  // wavespeeds scaled as 1/tau
  axl=aval[6+0];
  axr=aval[6+1];
  ayl=aval[6+2];
  ayr=aval[6+3];
  azl=aval[6+4];
  azr=aval[6+5];

  // coupled approach: radiative wavespeeds unlimited in optically thin medium (fcouple==0)
  // and equal to gas wavespeeds in optically thick medium (fcouple==1)
#ifdef GASRADCOUPLEDWAVESPEEDS
  axl=fcouple*axhdl+(1.-fcouple)*axl;
  axr=fcouple*axhdr+(1.-fcouple)*axr;
  ayl=fcouple*ayhdl+(1.-fcouple)*ayl;
  ayr=fcouple*ayhdr+(1.-fcouple)*ayr;
  azl=fcouple*azhdl+(1.-fcouple)*azl;
  azr=fcouple*azhdr+(1.-fcouple)*azr;
#endif
#endif //RADIATION

  // zero out 'co-going' velocities
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

  // save and pass up wavespeed
  // hd
  aaa[0]=axhdl;
  aaa[1]=axhdr;
  aaa[2]=ayhdl;
  aaa[3]=ayhdr;
  aaa[4]=azhdl;
  aaa[5]=azhdr;

  // rad unlimited by optical depth
  aaa[6]=axl0;
  aaa[7]=axr0;
  aaa[8]=ayl0;
  aaa[9]=ayr0;
  aaa[10]=azl0;
  aaa[11]=azr0;

  // rad affected by optical depth
  aaa[12]=axl;
  aaa[13]=axr;
  aaa[14]=ayl;
  aaa[15]=ayr;
  aaa[16]=azl;
  aaa[17]=azr;

  return 0;
}


//********************************************************************//
// Transform the fluid frame wavespeed into lab frame
//********************************************************************//

__device__ __host__ int calc_wavespeeds_lr_core_device(ldouble *ucon, ldouble GG[][5], ldouble *aret,
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
  
  indices_12_device(Bcov, Bcon, GG);
  Bsq = dot(Bcon, Bcov);
  Bu = dot(Bcov, ucon);
  Bu2 = Bu * Bu;

  // Now work on the relevant directions
  if(TNX > 1)  // x-direction 
  {
    Acov[0] = 0.;
    Acov[1] = 1.;
    Acov[2] = 0.;
    Acov[3] = 0.;
    
    indices_12_device(Acov, Acon, GG);
    
    Asq = dot(Acon, Acov);
    Au = dot(Acov, ucon);
    AB = dot(Acon, Bcov);
    Au2 = Au * Au;
    AuBu = Au * Bu;
    
    B = -2. * (AuBu * (1.0 - wspeed2x)  - AB * wspeed2x);
    A = Bu2 * (1.0 - wspeed2x) - Bsq * wspeed2x;
    discr = 4.0 * wspeed2x * ((AB * AB - Asq * Bsq) * wspeed2x + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2x - 1.0));
    
    if(discr < 0.)
    {
      printf("discr in x-wavespeeds lt 0\n");
      ierr = -1;
    }
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
  
  if(TNY > 1)  // y-direction is needed
  {
    Acov[0] = 0.;
    Acov[1] = 0.;
    Acov[2] = 1.;
    Acov[3] = 0.;
    
    indices_12_device(Acov, Acon, GG);
    
    Asq = dot(Acon, Acov);
    Au = dot(Acov, ucon);
    AB = dot(Acon, Bcov);
    Au2 = Au * Au;
    AuBu = Au * Bu;
    
    B = -2. * (AuBu * (1.0 - wspeed2y)  - AB * wspeed2y);
    A = Bu2 * (1.0 - wspeed2y) - Bsq * wspeed2y;
    discr = 4.0 * wspeed2y * ((AB * AB - Asq * Bsq) * wspeed2y + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2y - 1.0));
    
    if(discr < 0.)
    {
      printf("discr in y-wavespeeds lt 0\n");
      ierr = -1;
    }
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
  
  if(TNZ > 1)  // z-direction is needed
  {
    Acov[0] = 0.;
    Acov[1] = 0.;
    Acov[2] = 0.;
    Acov[3] = 1.;
    
    indices_12_device(Acov, Acon, GG);
    
    Asq = dot(Acon, Acov);
    Au = dot(Acov, ucon);
    AB = dot(Acon, Bcov);
    Au2 = Au * Au;
    AuBu = Au * Bu;
    
    B = -2. * (AuBu * (1.0 - wspeed2z)  - AB * wspeed2z);
    A = Bu2 * (1.0 - wspeed2z) - Bsq * wspeed2z;
    discr = 4.0 * wspeed2z * ((AB * AB - Asq * Bsq) * wspeed2z + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2z - 1.0));
    
    if(discr < 0.)
    {
      printf("discr in z-wavespeeds lt 0\n");
      ierr = -1;
    }
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
  
  if(ierr == 0)
  {
    return 0;
  }
  else
  {
    return -1;
  }
}

