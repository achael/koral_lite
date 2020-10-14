
int diskatboundary(ldouble *pp, void *ggg, void *gggBL)
{
  struct geometry *geom
     = (struct geometry *) ggg;
   
  int ix=geom->ix;
  int iy=geom->iy;
  int iz=geom->iz;

   struct geometry *geomBL
     = (struct geometry *) gggBL;
 
  ldouble r=geomBL->xxvec[1];
  ldouble th=geomBL->xxvec[2];
  ldouble phi=geomBL->xxvec[3];

  ldouble theq=fabs(M_PI/2.-th);
  ldouble theq2=M_PI/2.-th;
  ldouble thmax=DISKHR*M_PI/2.;

  ldouble phieq= fabs(phi);
  ldouble phieq2= phi;
  ldouble phimax=DISKHR*M_PI/2.;

  //for 3D only
  ldouble thp, thpmax;

  //ambient
  set_hdatmosphere(pp,geom->xxvec,geom->gg,geom->GG,0);
#ifdef RADIATION
  set_radatmosphere(pp,geom->xxvec,geom->gg,geom->GG,0);
#endif

  if(TNY>1 && TNZ==1)
  {
    if(theq>thmax)
      return -1; //outside disk; //This approach makes one cell atmosphere but prevents -rho 
  }
  else if(TNY==1 && TNZ>1)
  {
    if(phieq>phimax)
      return -1; //outside disk; //This approach makes one cell atmosphere but prevents -rho 
  }
  else if(TNY>1 && TNZ>1)
  {
    #ifdef STREAM_RING
    {
      if(theq>thmax)
        return -1; //outside disk; //This approach makes one cell atmosphere but prevents -rho 
    }
    #else
    if(theq>thmax || phieq>phimax)
      return -1; //outside disk; //This approach makes one cell atmosphere but prevents -rho 
    #endif
  }
  #ifdef SPECIAL_BC_CHECK // don't put stream fluid at outer cell
  int gix = ix + TOI;
  int giy = iy + TOJ;
  int giz = iz + TOK;
  //if(gix > (STREAM_IX+1))
  if(gix > (STREAM_IX+3))
    return -1;
  #endif

  int iiy;
  ldouble EPSILON,MDOTEDDIN,DISKVR,MDOTIN,DISKSIGMA,DISKRHO,TFB_STR;
  ldouble rho,uint,uphi,Om,Omk,ucon[4],temp,pre0,uint0,temp0,ell;
  ldouble meanuint,newT,mdotR,mdotL,mdotIN,rhorescale,rhoL,vL,rhoR,vR,rhoB,vB,mdotB,xx[4],dx[3];

  TFB_STR = TFB0;

  #ifndef USE_TPLUSTFB
  #ifdef EPSILON_TIME_DEPENDENT
  if(global_time <= TFB_STR) EPSILON = EPSILON0;
  else if(global_time > TFB_STR) EPSILON = EPSILON0*pow( (global_time/TFB_STR) , -2./3.);
  #else
  EPSILON = EPSILON0;
  #endif

  #ifdef MDOT_TIME_DEPENDENT
  MDOTEDDIN = MDOTEDD0/(1. + pow((global_time/TFB_STR),5./3.) );
  #else
  MDOTEDDIN = MDOTEDD0;
  #endif
  #endif

  #ifdef USE_TPLUSTFB
  ldouble stream_time = global_time + TFB_STR;

  #ifdef EPSILON_TIME_DEPENDENT
  EPSILON = EPSILON0*pow( (stream_time/TFB_STR) , -2./3.);
  #else
  EPSILON = EPSILON0;
  #endif

  #ifdef MDOT_TIME_DEPENDENT
  MDOTEDDIN = MDOTEDD0*pow((stream_time/TFB_STR),-5./3.);
  #else
  MDOTEDDIN = MDOTEDD0;
  #endif
  #endif

  DISKVR = -1.*sqrt(2.*(1./RMIN_DISK + EPSILON));
  MDOTIN = (MDOTEDDIN*2.48e18*MASS);

  //gas density
  if(TNY>1 && TNZ==1) //2D r-theta
  {
    DISKSIGMA = surfdensCGS2GU(-MDOTIN/(2.*M_PI*RMIN_DISK*147700.*MASS*3.e10*DISKVR));
    DISKRHO = (DISKSIGMA/2./DISKH);
  }
  else if(TNY==1 && TNZ>1) //2D r-phi
  {
    DISKSIGMA = surfdensCGS2GU(-MDOTIN/(2.*M_PI*RMIN_DISK*147700.*MASS*3.e10*DISKVR));
    DISKRHO = (DISKSIGMA/2./DISKH);
  }
  else if(TNY>1 && TNZ>1) //3D
  {
    #ifdef STREAM_RING
    DISKSIGMA = surfdensCGS2GU(-MDOTIN/(2.*M_PI*RMIN_DISK*147700.*MASS*3.e10*DISKVR));
    DISKRHO = (DISKSIGMA/2./DISKH);
    #else
    DISKRHO = rhoCGS2GU(-MDOTIN/(DISKHR*DISKHR*RMIN_DISK*RMIN_DISK*147700.*147700.*MASS*MASS*3.e10*DISKVR));
    #endif
  }

  ldouble rhopref = 1.;
  if(TNY>1 && TNZ==1) rhopref=pow(1.-pow(theq/thmax,2.),3.);
  if(TNY==1 && TNZ>1) rhopref=pow(1.-pow(phieq/phimax,2.),3.);
  if(TNY>1 && TNZ>1) 
  {
     #ifdef STREAM_RING
     rhopref=pow(1.-pow(theq/thmax,2.),3.);
     #else
     thp = sqrt( theq*theq + phieq*phieq );
     thpmax = 1.5*DISKHR*sqrt(2.);
     rhopref= pow(1.-pow(thp/thpmax,2.),3.);
     #endif
  }

  rho=DISKRHO*rhopref; //initial guess at what density should be


  if(rhopref < 0) rho = pp[0]; //Set rho to atm if happens to go negative

  //check the mass flux across STREAMIX-1/STREAMIX cell boundary (domain/ghost cell boundary)
  //rescale if need be. Note mdotL is a global variable set in ko.h and computed in ***
  newT = 0.;
  rhorescale = 1.;

  #ifdef FIX_INFLOW_RATE
  if(gix == STREAM_IX && giy >= STREAM_IYT && giy <= STREAM_IYB)
  {
    mdotL = 0.;
    meanuint = 0.;

    #ifndef CELL_BY_CELL_RESCALE
    for(iiy=0; iiy<(NY+1); ++iiy)
    {
      if((iiy+TOJ) >= STREAM_IYT && (iiy+TOJ) < (STREAM_IYB+1))
      {
        get_xx(ix-1, iiy, iz, xx);
        dx[0] = get_size_x(ix-1, 0);
        dx[1] = get_size_x(iiy, 1);
        dx[2] = get_size_x(iz, 2);
        if(NZ==1)
        {
          dx[2] = 2.* M_PI;
        }
        else
        {
          #ifdef PHIWEDGE
          dx[2]*=(2. * M_PI / PHIWEDGE);
          #endif
        }

        rhoR = rho;
        vR = -DISKVR;

        rhoL = get_u(p, 0, ix-1, iiy, iz);
        vL = my_max( get_u(p, 2, ix-1, iiy, iz)*(-1), 0.);
      
        rhoB = 0.5*(rhoL + rhoR);
        vB = 0.5*(vL + vR);

        mdotL += rhoL*vL*dx[1]*dx[2];
        meanuint += get_u(p, 1, ix-1, iiy, iz);
      }
    }
    #endif
  
    #ifdef CELL_BY_CELL_RESCALE
    get_xx(ix-1, iy, iz, xx);
    dx[0] = get_size_x(ix-1, 0);
    dx[1] = get_size_x(iy, 1);
    dx[2] = get_size_x(iz, 2);
    if(NZ==1)
    {
      dx[2] = 2.* M_PI;
    }
    else
    {
      #ifdef PHIWEDGE
      dx[2]*=(2. * M_PI / PHIWEDGE);
      #endif
    }

    rhoL = get_u(p, 0, ix-1, iy, iz);
    vL = get_u(p, 2, ix-1, iy, iz);//vL = my_max( get_u(p, 2, ix-1, iy, iz)*(-1), 0.);
    
    rhoR = rho;
    vR = DISKVR;
      
    rhoB = 0.5*(rhoL + rhoR);
    vB = 0.5*(vL + vR);

    mdotL = rhoL*vL*dx[1]*dx[2];
    meanuint = get_u(p, 1, ix-1, iy, iz);
    #endif

    mdotIN = MDOTIN*GGG/(CCC*CCC*CCC); //in GU
    mdotR = DISKRHO*(32./35.)*DISKHR*M_PI*M_PI*r*r*(-DISKVR); //in GU

    if( mdotR < (mdotIN+mdotL) || mdotR > (mdotIN+mdotL) )
      rhorescale= ((mdotIN+mdotL)/mdotR);

    //Rescale pressure to match that outside of stream? (this will lead to a hot stream, but maybe it should be heated later as disk forms and interacts with it?)
    #ifndef CELL_BY_CELL_RESCALE
    meanuint = meanuint/(STREAM_IYB + 1 - STREAM_IYT);
    #endif

    newT=calc_PEQ_Tfromurho(meanuint,DISKRHO*rhorescale,ix,iy,iz);

    rho = rho*rhorescale;
    rhoR = (-2.*mdotIN/vB/(DISKHR*M_PI*M_PI*r*r)) - rhoL;
    //if(rhoR > 0.) rho = rhoR;
  }
  else if(gix == (STREAM_IX+1))
  {
    mdotIN = MDOTIN*GGG/(CCC*CCC*CCC); //in GU
    mdotR = DISKRHO*(32./35.)*DISKHR*M_PI*M_PI*r*r*(-DISKVR); //in GU

    if( mdotR < mdotIN || mdotR > mdotIN )
      rhorescale= (mdotIN/mdotR);

    rho = rho*rhorescale;
  }
  #endif

  temp=DISKTEMP;

  //switching off temp fixing for now
  //#ifdef FIX_INFLOW_RATE
  //if(newT > temp) temp = newT;
  //#endif

  uint=calc_PEQ_ufromTrho(temp,rho,ix,iy,iz);

  Omk=1./sqrt(r*r*r);

  // pre0 = DISKSIGMA * Omk * Omk * DISKH * DISKH / 2. / DISKH ;  //P=2Hp, P/S = Om^2 H^2 
  //uint0 = pre0 / GAMMAM1;
  //temp0 = calc_PEQ_Tfromurho(uint0, DISKRHO);

  //ang.momentum corresponding to DISKRCIR
  ell=sqrt(DISKRCIR);
  Om=ell/r/r;
  

  ucon[1]=DISKVR;
  ucon[2]=0.;
  ucon[3]=Om;
    
  conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL->gg,geomBL->GG);
 
  pp[0]=my_max(rho,pp[0]);
  pp[1]=my_max(uint,pp[1]);
  pp[2]=ucon[1]; 
  pp[3]=ucon[2];
  pp[4]=ucon[3];

  //print_primitives(pp);getchar();

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
  pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

#ifdef RADIATION

#ifdef STR_USETRADATM

#ifdef SPECIAL_BC_CHECK
  if(gix == (STREAM_IX)) pp[EE0]=get_u(p,EE0,ix-1,iy,iz);
  if(gix == (STREAM_IX+1)) pp[EE0]=get_u(p,EE0,ix-2,iy,iz);

  //set flux to zero in fluid frame : Assumes opt. thick?
  pp[FX0]=0.;
  pp[FY0]=0.;
  pp[FZ0]=0.;

  #ifdef EVOLVEPHOTONNUMBER
    pp[NF0]=calc_NFfromE(pp[EE0]);
  #endif
#endif

#else 
  //split pressure into rad and gas assuming LTE
  ldouble Eatm=pp[EE0];
  ldouble E,Fx, Fy,Fz;
  //distributing pressure
  ldouble P,aaa,bbb;
  P=GAMMAM1*pp[1];
  //solving for T satisfying P=pgas+prad=bbb T + aaa T^4
  aaa=4.*SIGMA_RAD/3.;
  bbb=K_BOLTZ*pp[0]/MU_GAS/M_PROTON;
  ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));
  ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;
    
  E=calc_LTE_EfromT(T4);
  Fx=Fy=Fz=0.;
  uint=calc_PEQ_ufromTrho(T4,rho,ix,iy,iz);

  pp[UU]=uint;//my_max(uint,uintatm);
  pp[EE0]=E;//my_max(E,Eatm);

  pp[FX0]=Fx;
  pp[FY0]=Fy;
  pp[FZ0]=Fz;

  //#ifdef NCOMPTONIZATION
    //pp[NF0]=calc_NFfromE(pp[EE0]);
    //pp[NF]=calc_NFfromT(calc_PEQ_Tfromurho(pp[UU],pp[RHO]));
    //pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
    //pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/T4;
  //#endif
  #ifdef EVOLVEPHOTONNUMBER
    pp[NF0]=calc_NFfromE(pp[EE0]);
  #endif

#endif //ifndef STR_USETRADATM : Set flux to zero in fluid frame, but it should be non-zero in labVVV

  //transforming from BL lab radiative primitives to code non-ortonormal primitives
  prad_ff2lab(pp,pp,geomBL);
#endif //ifdef RADIATION

  //transforming primitives from BL to MYCOORDS
  trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL->xxvec,gggBL,ggg);
    
 

#ifdef MAGNFIELD 
  //artificially impose poloidal magnetic field
  ldouble Pgas,Prad;
  Pgas=GAMMAM1*uint;
#ifdef RADIATION
  Prad=pp[EE0]/3.;
#else
  Prad=0.;
#endif
  /*old B2
  pp[B2]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom->gg[2][2])*cos(MAGNOMEGA*global_time);
  */
  
  //pp[B2]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom->gg[2][2])*sin(theq2/thmax*M_PI)*cos(MAGNOMEGA*global_time);

  //radial / vertical
#ifdef BDIPOLAR
  pp[B1]=sqrt(DISKRHO*MAGBETA)/sqrt(geom->gg[1][1])*cos(theq2/thmax*M_PI)*cos(MAGNOMEGA*global_time);
  if(global_time < VERTBTIME)
    pp[B2]=sqrt(DISKRHO*MAGBETA)/sqrt(geom->gg[2][2])*sin(theq2/thmax*M_PI)*cos(MAGNOMEGA*global_time);
  /* pp[B1]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom->gg[1][1])*cos(theq2/thmax*M_PI)*cos(MAGNOMEGA*global_time);
  if(global_time < VERTBTIME)
  pp[B2]=sqrt(2.*(Pgas+Prad)*MAGBETA/2.)/sqrt(geom->gg[2][2])*cossin(theq2/thmax*M_PI)*cos(MAGNOMEGA*global_time);*/
  #endif

#ifdef VERTICALB
  //pp[B1]=sqrt(2.*(Pgas+Prad)*MAGBETA)/sqrt(geom->gg[1][1])*sin(theq2/thmax*M_PI*2.)*cos(MAGNOMEGA*global_time); //why 2 here?
  //if(global_time < VERTBTIME)
  #ifdef SANEFIELD
    pp[B2]=sqrt(2.*(Pgas+Prad)*MAGBETA/2.)/sqrt(geom->gg[2][2])*sin(theq2/thmax*M_PI)*my_sign(cos(MAGNOMEGA*global_time));

  #else
    pp[B2]=sqrt(2.*(Pgas+Prad)*MAGBETA/2.)/sqrt(geom->gg[2][2])*cos(theq2/thmax*M_PI)*my_sign(cos(MAGNOMEGA*global_time));
  #endif
#endif

  //Debug lines - Brandon
  /*
  //MYCOORDS vector potential to calculate B's
  ldouble Acov[4];
  Acov[0]=Acov[1]=Acov[2]=0.;
  Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/4.e-20,2.)-0.02,0.)*sqrt(1.e-23)*pow(sin(fabs(geomBL.yy)),4.);

  pp[B1]=0.;
  pp[B2]=0.;
  pp[B3]=Acov[3];
  */
#endif


#ifdef PERTMAGN //perturb to break axisymmetry
pp[UU]*=1.+PERTMAGN*sin(10.*2.*M_PI*(MAXZ-geomBL->zz)/(MAXZ-MINZ));
#endif

//Calculate entropy?
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);
  
//Special BC. Done here for simplicity instead of bc.c. No inflow in X into stream region
#ifdef SPECIAL_BC_CHECK
  int iv = 0; 
  int shutoffatinner=0;

  //Reflecting boundary in front after shutoff
  #ifdef SHUTOFF_AFTER_TSHUTOFF
  #ifdef TSHUTOFF
  if(global_time > TSHUTOFF) shutoffatinner = 1; //turn off inflow completely
  #else
  ldouble TLB;
  if(EPSILON_MAX < 0.) //Only use dynamics based model if all mass returns in finite time
  {
    TLB = 2.*M_PI*pow((-1./(2.*EPSILON_MAX)),1.5);
    if(global_time > (TLB-TFB0)) shutoffatinner = 1; //apply reflecting boundary as in XBCHI?
  }
  #endif

  if(shutoffatinner)
  {
    if(gix == (STREAM_IX))
    {
      for(iv=0; iv<NV; ++iv)
      {
        #ifdef MAGNFIELD
        if(iv==VX || iv==B1 || iv==FX0)
        #else
        if(iv==VX || iv==FX0)
        #endif
          pp[iv]=-get_u(p,iv,ix-1,iy,iz);
        else
          pp[iv]=get_u(p,iv,ix-1,iy,iz);
      }
    }

    if(gix == (STREAM_IX+1))
    {
      for(iv=0; iv<NV; ++iv)
      {
        #ifdef MAGNFIELD
        if(iv==VX || iv==B1 || iv==FX0)
        #else
        if(iv==VX || iv==FX0)
        #endif
          pp[iv]=-get_u(p,iv,ix-3,iy,iz);
        else
          pp[iv]=get_u(p,iv,ix-3,iy,iz);
      }
    }
  }
  #endif

  //Reflecting boundary behind stream
  if(gix == (STREAM_IX+3))
  {
    for(iv=0; iv<NV; ++iv)
    {
      #ifdef MAGNFIELD
      if(iv==VX || iv==B1 || iv==FX0)
      #else
      if(iv==VX || iv==FX0)
      #endif
        pp[iv]=-get_u(p,iv,ix+1,iy,iz);
      else
        pp[iv]=get_u(p,iv,ix+1,iy,iz);
    }

  #ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.; 
  #endif
  }
  if(gix == (STREAM_IX+2))
  {
    for(iv=0; iv<NV; ++iv)
    {
      #ifdef MAGNFIELD
      if(iv==VX || iv==B1 || iv==FX0)
      #else
      if(iv==VX || iv==FX0)
      #endif
        pp[iv]=-get_u(p,iv,ix+3,iy,iz);
      else
        pp[iv]=get_u(p,iv,ix+3,iy,iz);
    }

  #ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.; 
  #endif
  }
#endif

  return 0;
}


int if_indisk(int ix, int iy, int iz, int ixs, int iyst, int iysb)
{
  struct geometry geomBL;
  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
 
  ldouble r=geomBL.xxvec[1];
  ldouble th=geomBL.xxvec[2];
  ldouble phi=geomBL.xxvec[3];

  ldouble theq=fabs(M_PI/2.-th);
  ldouble theq2=M_PI/2.-th;
  ldouble theqmax=DISKHR*M_PI/2.;
  ldouble thmin = -theqmax;
  ldouble thmax = theqmax;

  //based around phieq = [-pi,+pi]
  ldouble phieq= fabs(phi);
  ldouble phieq2= phi;
  ldouble phieqmax=DISKHR*M_PI/2.;
  ldouble phimin = -phieqmax;
  ldouble phimax = phieqmax;

  #ifdef SPECIAL_BC_CHECK 
  //find ix_stream, iy_stream_top, iy_stream_bot
  int iix,iiy,iiz;
  int ix_str, iy_str_top, iy_str_bot, iz_str_top, iz_str_bot;
  ix_str = STREAM_IX; 
  iy_str_top = STREAM_IYT;
  iy_str_bot = STREAM_IYB;
  iz_str_top = STREAM_IZT;
  iz_str_bot = STREAM_IZB;

  #ifdef SEARCH_STREAM_BOUNDARY //ONLY USE IN SERIAL!
  for(iix=0;iix<TNX-1;iix++)
  {
    ldouble rrr = 0.;
    ldouble rrr2 = 0.;
    fill_geometry_arb(iix,0,0,&geomBL,KERRCOORDS);
    rrr = geomBL.xxvec[1];
    fill_geometry_arb(iix+1,0,0,&geomBL,KERRCOORDS);
    rrr2 = geomBL.xxvec[1];

    if(rrr < RMIN_DISK && rrr2 >= RMIN_DISK)
    {
      ix_str = iix+1;
      break;
    }
    //MPI case 1: disk is in domain but not the edge
    //MPI case 2: disk is in domain but only one cell (should try to divide cells so this doesn't happen)
    //MPI case 3: disk is not in domain
    /*
    if(iix == NX-2)
    {
      if(rrr2 < RMIN_DISK) 
      {
        ix_str = -1;
        break;
      }
    }
    */
  }

  //Theta boundary
  if(TNY > 1)
  {
    for(iiy=0;iiy<TNY-1;iiy++)
    {
      ldouble thl = 0.;
      ldouble thl2 = 0.;
      fill_geometry_arb(ix_str,iiy,0,&geomBL,KERRCOORDS);
      thl = geomBL.xxvec[2];
      fill_geometry_arb(ix_str,iiy+1,0,&geomBL,KERRCOORDS);
      thl2 = geomBL.xxvec[2];
      thl = M_PI/2. - thl;
      thl2 = M_PI/2. - thl2;

      if(thl >= thmin && thl2 < thmax)
      {
        iy_str_top = iiy;
        break;
      }
    }

    for(iiy=TNY-1;iiy>0;iiy--)
    {
      ldouble thl = 0.;
      ldouble thl2 = 0.;
      fill_geometry_arb(ix_str,iiy,0,&geomBL,KERRCOORDS);
      thl = geomBL.xxvec[2];
      fill_geometry_arb(ix_str,iiy-1,0,&geomBL,KERRCOORDS);
      thl2 = geomBL.xxvec[2];
      thl = M_PI/2. - thl;
      thl2 = M_PI/2. - thl2;

      if(thl <= thmax && thl2 > thmin)
      {
        iy_str_bot = iiy;
        break;
      }
    }
  }

  //Phi boundary
  #ifndef STREAM_RING
  if(TNZ > 1)
  {
    for(iiz=0;iiz<TNZ-1;iiz++)
    {
      ldouble phil = 0.;
      ldouble phil2 = 0.;
      fill_geometry_arb(ix_str,0,iiz,&geomBL,KERRCOORDS);
      phil = geomBL.xxvec[3];
      fill_geometry_arb(ix_str,0,iiz+1,&geomBL,KERRCOORDS);
      phil2 = geomBL.xxvec[3];
      phil = phil;// - M_PI;
      phil2 = phil2;// - M_PI;

      if(phil >= phimin && phil2 < phimax)
      {
        iz_str_top = iiz;
        break;
      }
    }

    for(iiz=TNZ-1;iiz>0;iiz--)
    {
      ldouble phil = 0.;
      ldouble phil2 = 0.;
      fill_geometry_arb(ix_str,0,iiz,&geomBL,KERRCOORDS);
      phil = geomBL.xxvec[3];
      fill_geometry_arb(ix_str,0,iiz-1,&geomBL,KERRCOORDS);
      phil2 = geomBL.xxvec[3];
      phil = phil;// - M_PI;
      phil2 = phil2;// - M_PI;


      if(phil <= phimax && phil2 > phimin)
      {
        iz_str_bot = iiz;
        break;
      }
    }
  }
  #endif

  printf("ixs iyst iysb izst izsb : %i %i %i %i %i\n",ix_str,iy_str_top+1,iy_str_bot-1,iz_str_top+1,iz_str_bot-1);
  //#define STREAM_IX ix_str
  //#define STREAM_IYT iy_str_top
  //#define STREAM_IYB iy_str_bot
  //#undef SPECIAL_BC_CHECK
  #endif  
  if(ix < TNX) fill_geometry_arb(ix-1,iy,iz,&geomBL,KERRCOORDS);
  ldouble r2=geomBL.xxvec[1];

  //if(ix >= NX) return 1; //apply XBCHI?

  if(TNY>1 && TNZ==1) //r-theta 2D
  {
    if(ix >= ix_str && iy >= iy_str_top && iy <= iy_str_bot) //inside stream region?
    {
      //if(ix == ix_str || ix == (ix_str+1))
      if(ix == ix_str || ix == (ix_str+1) || ix == (ix_str+2) || ix == (ix_str+3))
      {
        return 1; //gc with XBCHI
      }
      else if(ix < TNX)
      {
        return 0;
      }
      else if(ix >= TNX)
      {
        return 1;
      }
    }
    else
    {
      return 0;
    }
  }
  else if(TNY==1 && TNZ>1) //r-phi 2D
  {
    if(ix >= ix_str && iz >= iz_str_top && iz <= iz_str_bot) //inside stream region?
    {
      if(ix == ix_str || ix == (ix_str+1) || ix == (ix_str+2) || ix == (ix_str+3))
      {
        return 1;
      }
      else if(ix < TNX)
      {
        return 0;
      }
      else if(ix >= TNX)
      {
        return 1;
      }
    }
    else
    {
      return 0;
    }
  }
  else if(TNY>1 && TNZ>1) //3D
  {
    #ifdef STREAM_RING
    if(ix >= ix_str && iy >= iy_str_top && iy <= iy_str_bot) //inside stream region?
    #else
    if(ix >= ix_str  && iy >= iy_str_top && iy <= iy_str_bot && iz >= iz_str_top && iz <= iz_str_bot) //inside stream region?
    #endif
    {
      if(ix == ix_str || ix == (ix_str+1) || ix == (ix_str+2) || ix == (ix_str+3))
      {
        return 1;
      }
      else if(ix < TNX)
      {
        return 0;
      }
      else if(ix >= TNX)
      {
        return 1;
      }
    }
    else
    {
      return 0;
    }
  }

  /***
  if(theq>thmax)
  {
    return 0; //outside disk, do nothing special;
  }
  else
  {
    if(r >= RMIN_DISK && r2 < RMIN_DISK) //use XBCHI only at innermost cell
    {
      return 1;
    }
    else if(r2 >= RMIN_DISK && r > RMIN_DISK)
    {
      if(theq2 > 0) //Use YBCHI
      {
        return 2;
      }
      else if(theq2 <= 0) //Use YBCLOW
      {
        return 3;
      }
    }
    else if(r < RMIN_DISK)
    {
      return 0;
    }
  }
  ***/
  #endif

  return 0;
}


//Only used for XBC/YBC checks
int if_indisk_bc_check(int ix, int iy, int iz)
{
  struct geometry geomBL;
  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
 
  ldouble r=geomBL.xxvec[1];
  ldouble th=geomBL.xxvec[2];
  ldouble phi=geomBL.xxvec[3];

  ldouble theq=fabs(M_PI/2.-th);
  ldouble theq2=M_PI/2.-th;
  ldouble theqmax=DISKHR*M_PI/2.;
  ldouble thmin = -theqmax;
  ldouble thmax = theqmax;

  //based around phieq = [-pi,+pi]
  ldouble phieq= fabs(phi);
  ldouble phieq2= phi;
  ldouble phieqmax=DISKHR*M_PI/2.;
  ldouble phimin = -phieqmax;
  ldouble phimax = phieqmax;

  #ifdef SPECIAL_BC_CHECK 
  //find ix_stream, iy_stream_top, iy_stream_bot
  int iix,iiy,iiz;
  int ix_str, iy_str_top, iy_str_bot, iz_str_top, iz_str_bot;
  ix_str = STREAM_IX; 
  iy_str_top = STREAM_IYT;
  iy_str_bot = STREAM_IYB;
  iz_str_top = STREAM_IZT;
  iz_str_bot = STREAM_IZB;

  if(ix < TNX) fill_geometry_arb(ix-1,iy,iz,&geomBL,KERRCOORDS);
  ldouble r2=geomBL.xxvec[1];

  //if(ix >= NX) return 1; //apply XBCHI?

  if(TNY>1 && TNZ==1)
  {
    if(ix >= ix_str && iy >= iy_str_top && iy <= iy_str_bot) //inside stream region?
    {
      //if(ix == ix_str || ix == (ix_str+1) )
      if(ix == ix_str || ix == (ix_str+1) || ix == (ix_str+2) || ix == (ix_str+3) )
      {
        if(iy == iy_str_top) return 0;//-1; // gc with either special corner bc or YBCHI
        else if(iy == iy_str_bot) return 0;//1;  // gc with either special corner bc or YBCHI
        else if(iy > iy_str_top && iy < iy_str_top) return 0; //gc with stream inflow
      }
      else if( ix == TNX || ix == (TNX+1) )
      {
        return 0; //XBCHI, no disk material
      }
      else
      {
        if(iy == iy_str_top) return -1; //1st gc with YBCHI
        else if(iy == (iy_str_top+1)) return -3; //2nd gc with YBCHI
        else if(iy == iy_str_bot) return 1; //1st gc with YBCLO
        else if(iy == (iy_str_bot-1) ) return 3; //2nd gc with YBCLO
        else return 0;
      }
    }
    else
    {
      return 0;
    }
  }
  if(TNY==1 && TNZ>1)
  {
    if(ix >= ix_str && iz >= iz_str_top && iz <= iz_str_bot) //inside stream region?
    {
      //if(ix == ix_str || ix == (ix_str+1) )
      if(ix == ix_str || ix == (ix_str+1) || ix == (ix_str+2) || ix == (ix_str+3) )
      {
        if(iz == iz_str_top) return 0;//-1; // gc with either special corner bc or YBCHI
        else if(iz == iz_str_bot) return 0;//1;  // gc with either special corner bc or YBCHI
        else if(iz > iz_str_top && iz < iz_str_top) return 0; //gc with stream inflow
      }
      else if( ix == TNX || ix == (TNX+1) )
      {
        return 0; //XBCHI, no disk material
      }
      else
      {
        if(iz == iz_str_top) return -1; //1st gc with ZBCHI
        else if(iz == (iz_str_top+1)) return -3; //2nd gc with ZBCHI
        else if(iz == iz_str_bot) return 1; //1st gc with ZBCLO
        else if(iz == (iz_str_bot-1) ) return 3; //2nd gc with ZBCLO
        else return 0;
      }
    }
    else
    {
      return 0;
    }
  }
  if(TNY>1 && TNZ>1)
  {
    #ifdef STREAM_RING
    if(ix >= ix_str && iy >= iy_str_top && iy <= iy_str_bot) //inside stream region?
    {
      //if(ix == ix_str || ix == (ix_str+1) )
      if(ix == ix_str || ix == (ix_str+1) || ix == (ix_str+2) || ix == (ix_str+3) )
      {
        if(iy == iy_str_top) return 0;//-1; // gc with either special corner bc or YBCHI
        else if(iy == iy_str_bot) return 0;//1;  // gc with either special corner bc or YBCHI
        else if(iy > iy_str_top && iy < iy_str_top) return 0; //gc with stream inflow
      }
      else if( ix == TNX || ix == (TNX+1) )
      {
        return 0; //XBCHI, no disk material
      }
      else
      {
        if(iy == iy_str_top) return -1; //1st gc with YBCHI
        else if(iy == (iy_str_top+1)) return -3; //2nd gc with YBCHI
        else if(iy == iy_str_bot) return 1; //1st gc with YBCLO
        else if(iy == (iy_str_bot-1) ) return 3; //2nd gc with YBCLO
        else return 0;
      }
    }
    else
    {
      return 0;
    }
    #else
    if(ix >= ix_str && iy >= iy_str_top && iy <= iy_str_bot && iz >= iz_str_top && iz <= iz_str_bot) //inside stream region?
    {
      //if(ix == ix_str || ix == (ix_str+1) )
      if(ix == ix_str || ix == (ix_str+1) || ix == (ix_str+2) || ix == (ix_str+3) )
      {
        if(iy == iy_str_top) return 0;//-1; // gc with either special corner bc or YBCHI
        else if(iy == iy_str_bot) return 0;//1;  // gc with either special corner bc or YBCHI
        else if(iy > iy_str_top && iy < iy_str_top) return 0; //gc with stream inflow

        if(iz == iz_str_top) return 0;//-1; // gc with either special corner bc or YBCHI
        else if(iz == iz_str_bot) return 0;//1;  // gc with either special corner bc or YBCHI
        else if(iz > iz_str_top && iz < iz_str_top) return 0; //gc with stream inflow
      }
      else if( ix == TNX || ix == (TNX+1) )
      {
        return 0; //XBCHI, no disk material
      }
      else //really not used currently
      {
        if(iy == iy_str_top) return -1; //1st gc with YBCHI
        else if(iy == (iy_str_top+1)) return -3; //2nd gc with YBCHI
        else if(iy == iy_str_bot) return 1; //1st gc with YBCLO
        else if(iy == (iy_str_bot-1) ) return 3; //2nd gc with YBCLO
        else if(iz == iz_str_top) return -1; //1st gc with ZBCHI
        else if(iz == (iz_str_top+1)) return -3; //2nd gc with ZBCHI
        else if(iz == iz_str_bot) return 1; //1st gc with ZBCLO
        else if(iz == (iz_str_bot-1) ) return 3; //2nd gc with ZBCLO
        else return 0;
      }
    }
    else
    {
      return 0;
    }
    #endif
  }
  #endif

  return 0;
}
