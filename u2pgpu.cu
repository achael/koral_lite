/*! \file u2pgpu.cu
 \brief MHD Conserved to primitives conversion
 */

extern "C" {

#include "ko.h"

}

#include "kogpu.h"

#define ERRCONV 1.e-1

// reminder of default settings
//#define U2PCONV 1.e-1039

//#define U2P_EQS U2P_EQS_NOBLE
//#define U2P_SOLVER U2P_SOLVER_W

__device__ __host__ static ldouble dpdWp_calc_vsq(ldouble Wp, ldouble D, ldouble vsq,ldouble gamma);
__device__ __host__ static ldouble compute_idwmrho0dp(ldouble wmrho0,ldouble gamma);
__device__ __host__ static ldouble compute_idrho0dp(ldouble wmrho0);
__device__ __host__ static int f_u2p_hot(ldouble Wp, ldouble* cons,ldouble *f,ldouble *df,ldouble *err,ldouble pgamma);
__device__ __host__ static ldouble compute_specificentropy_wmrho0_idealgas(ldouble rho0, ldouble wmrho0,ldouble gamma);
__device__ __host__ static ldouble compute_dspecificSdwmrho0_wmrho0_idealgas(ldouble rho0, ldouble wmrho0,ldouble gamma);
__device__ __host__ static ldouble compute_dspecificSdrho_wmrho0_idealgas(ldouble rho0, ldouble wmrho0, ldouble gamma);
__device__ __host__ static int f_u2p_entropy(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df, ldouble *err,ldouble pgamma);

/*
__device__ __host__ int set_cflag_device(int *cellflag_arr, int iflag,int ix,int iy,int iz, int val)
{
   cellflag_arr[iflag + (iX(ix)+(NGCX))*NFLAGS +  \
			(iY(iy)+(NGCY))*(SX)*NFLAGS + \
		        (iZ(iz)+(NGCZ))*(SY)*(SX)*NFLAGS] = val;
   return 0;
}
*/

//**********************************************************************
//primitive to conserved converter 
//**********************************************************************

__device__ __host__ int p2u_device(ldouble *p, ldouble *u, void *ggg)
{
  p2u_mhd_device(p,u,ggg);

  /* //TODO 
#ifdef RADIATION
  p2u_rad_device(p,u,ggg);
#endif
  */
  
  return 0;
}


__device__ __host__ int p2u_mhd_device(ldouble *pp, ldouble *uu, void *ggg)
{

#ifdef NONRELMHD
  //p2u_mhd_nonrel(p,u,ggg); //TODO nonrel
  return 0;
#endif

  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdetu;
  gg=geom->gg;
  GG=geom->GG;
  
  gdetu=geom->gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif

  ldouble rho = pp[0];
  ldouble uint= pp[1];
  ldouble S   = pp[5];
  
  ldouble vcon[4],ucon[4],ucov[4];
  vcon[0]=0.;
  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];
  conv_vels_both_device(vcon,ucon,ucov,VELPRIM,VEL4,gg,GG);
  
  ldouble bcon[4]={0.,0.,0.,0.}, bcov[4]={0.,0.,0.,0.}, bsq=0.;
#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel_device(pp, ucon, ucov, geom, bcon, bcov, &bsq);
#endif

  //************************************
  //hydro part
  //************************************
 
  ldouble ut=ucon[0];
  ldouble rhout = rho*ut;
  ldouble Sut = S*ut;

  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  //gamma=pick_gammagas(geom->ix,geom->iy,geom->iz); //TODO
  #endif
 
  ldouble pre=(gamma-1.)*uint; 
  ldouble eta  = rho+uint+pre+bsq;
  ldouble etap = uint+pre+bsq; //eta-rho
  ldouble ptot = pre+0.5*bsq;
  
  //this computes utp1=1+u_t
  ldouble utp1 = calc_utp1_device(vcon,ucon,geom); 

  ldouble Tttt =etap*ucon[0]*ucov[0] + rho*ucon[0]*utp1 + ptot - bcon[0]*bcov[0];
  ldouble Ttr  =eta *ucon[0]*ucov[1] - bcon[0]*bcov[1];
  ldouble Ttth =eta *ucon[0]*ucov[2] - bcon[0]*bcov[2];
  ldouble Ttph =eta *ucon[0]*ucov[3] - bcon[0]*bcov[3];

  uu[0]=gdetu*rhout;
  uu[1]=gdetu*Tttt;
  uu[2]=gdetu*Ttr;
  uu[3]=gdetu*Ttth;
  uu[4]=gdetu*Ttph;
  uu[5]=gdetu*Sut;

#ifdef EVOLVEELECTRONS
  uu[ENTRE]= gdetu*pp[ENTRE]*ut;
  uu[ENTRI]= gdetu*pp[ENTRI]*ut;
#endif

#ifdef RELELECTRONS
  int ib;
  for(int ib=0;ib<NRELBIN;ib++)
    uu[NEREL(ib)]=gdetu*pp[NEREL(ib)]*ut;    
#endif

  //************************************
  //magnetic part
  //************************************ 
#ifdef MAGNFIELD
  uu[B1]=gdetu*pp[B1];
  uu[B2]=gdetu*pp[B2];
  uu[B3]=gdetu*pp[B3];
#endif

  return 0.;
}


//**********************************************************************
// Compute utp1=1+u_t , which for nonrelativistic cases is ~0.
// if computed directly as 1+u_t, then if the residual is small there will be a large error.
//**********************************************************************
__device__ __host__ ldouble calc_utp1_device(ldouble *vcon, ldouble *ucon, void *ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  ldouble utp1;
  if(VELPRIM==VELR) //based on VELR
  {
  
    ldouble qsq=0.;
    for(int i=1;i<4;i++)
      for(int j=1;j<4;j++)
	 qsq+=vcon[i]*vcon[j]*gg[i][j];

    ldouble gamma2=(1.+qsq);
    ldouble alpha = geom->alpha;
    ldouble alphasq = alpha*alpha;
    ldouble alpgam=sqrt(alphasq*gamma2);

    //\beta^i \beta_i / \alpha^2 = g^{ti} g_{ti}
    ldouble betasqoalphasq=gg[0][1]*GG[0][1] + gg[0][2]*GG[0][2] + gg[0][3]*GG[0][3];

    // \tilde{u}_t = \tilde{u}^i g_{ti} since \tilde{u}^t=0
    ldouble ud0tilde = 0.0;
    for(int j=1;j<4;j++)
      ud0tilde += vcon[j]*gg[0][j]; 

    utp1= ud0tilde + (geom->gttpert - alphasq*(betasqoalphasq + qsq))/(1.0+alpgam);
  }
  else //based on ucon[]
  {
    
    // 3-velocity in coordinate basis
    ldouble vconp[4];
    for(int j=1;j<4;j++)
      vconp[j]=ucon[j]/ucon[0];

    // 3 velcocity squared
    ldouble vsq=geom->gttpert;
    for(int j=1;j<4;j++)
      vsq+=2.0*geom->gg[0][j]*vconp[j];

    for(int j=1;j<4;j++)
      for(int k=1;k<4;k++)
	vsq+=geom->gg[j][k]*vconp[j]*vconp[k];
  
    ldouble alpha=0.0;
    for(int j=1;j<4;j++)
      alpha+=geom->gg[j][0]*ucon[j];

    ldouble plus1gv00=geom->gttpert;
    ldouble uu0 = ucon[0];
    ldouble gvtt=geom->gg[0][0];

    utp1 = alpha + ((1.0-gvtt)*plus1gv00 - uu0*uu0*vsq*gvtt*gvtt)/(1.0-gvtt*uu0);
  }
  
  return utp1;
}


//**********************************************************************
//high-level u2p solver
//**********************************************************************

__device__ int u2p_device(ldouble *uu0, ldouble *pp, void *ggg,
			  int corrected[3], int fixups[2], int* int_slot_arr)
{
  int verbose=0;
  
  int ret=0;
  int u2pret=-1;
  int hdcorr=0;
  int radcor=0;
  corrected[0]=corrected[1]=corrected[2]=0;
  fixups[0]=fixups[1]=0;

  struct geometry *geom
  = (struct geometry *) ggg;

  ldouble gdetu, gdetu_inv;
  gdetu=geom->gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif
  gdetu_inv = 1. / gdetu;
  
  ldouble uu[NV];

  for(int iv=0;iv<NV;iv++)
    uu[iv]=uu0[iv];
    
  ldouble ppbak[NV];
  for(int iv=0;iv<NV;iv++)
    ppbak[iv]=pp[iv];

  //************************************
  //magneto-hydro part
  //************************************
  ldouble u0=pp[1]; // save initial energy to print if something goes wrong
  
  //************************************
  //hot hydro - conserving energy
  
  //negative uu[0] = rho u^t
  if(uu[0] * gdetu_inv < 0.)
  {
    //TODO -- print statement with mpi data
    /*
    if verbose
    {
	int gix,giy,giz;
	mpi_local2globalidx(geom->ix,geom->iy,geom->iz,&gix,&giy,&giz);
	printf("%4d > %4d %4d %4d > NEGUU  > neg uu[0] - requesting fixup\n",PROCID,gix,giy,giz);
    }
    */
    pp[0]=RHOFLOOR; //used when not fixing up
    uu[0]=RHOFLOOR*gdetu;
    ret=-2;    //to request fixup
               //TODO what do we want here?
               //ANDREW -- but ret goes back to -1 if energy inversion fails but entropy inversion does not!
               //ANDREW -- do we always want a fixup if we have negative uu[0] ? 
    
    atomicAdd(&int_slot_arr[GLOBALINTSLOT_NTOTALMHDFIXUPS],1); //TODO right??
    //global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]++;  //this counts as a fixup
  }

  if(u2pret!=0)  // u2pret=-1 at this stage, so this is always satisfied
  {
#ifdef ENFORCEENTROPY  
    u2pret=-1;  //skip hot energy-conserving inversion and go to entropy inversion
#else
    u2pret = u2p_solver_device(uu,pp,ggg,U2P_HOT,0);  // invert using the hot energy equation    
#endif //ENFORCEENTROPY
  }

  //************************************
  //entropy solver - conserving entropy
 
  if(ALLOWENTROPYU2P)  // allow entropy equation inversion -- on by default (see choices.h)
  {
    if(u2pret<0)  // true if energy equation failed, or if ENFORCEENTROPY is defined
    {
      ret=-1;
      
      if(verbose>2 )
      {
	//TODO -- print with mpi variables
	printf("u2p_entr     >>> %d %d <<< %d >>> %e > %e\n",geom->ix, geom->iy,u2pret,u0,pp[1]);
        //printf("u2p_entr     >>> %d %d <<< %d >>> %e > %e\n",geom->ix + TOI, geom->iy + TOJ,u2pret,u0,pp[1]);
      }
      
      u2pret=u2p_solver_device(uu,pp,ggg,U2P_ENTROPY,0);  // invert using entropy equation
      
      if(u2pret<0)
      {
        ret=-2;
        
        if(verbose>1)
        {
          printf("u2p_entr err No. %4d > %e %e %e > %e %e > %4d %4d %4d\n",
		 u2pret,uu[0],uu[1],uu[5],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
        }
	
      } // if(u2pret<0) // second time -- entropy eqn
    } // if(u2pret<0) // first time -- energy eqn
  }  // if(ALLOWENTROPYU2P)


  //************************************
  // both energy and entropy inversion failed

  if(u2pret<0)
  {
 
    //leaving primitives unchanged - should not happen
    if(verbose>1 || 1)
    {
      //TODO print with mpi variables
      //printf("%4d > %4d %4d %4d > MHDU2PFAIL > u2p prim. unchanged > %d \n",PROCID,geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,u2pret);
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
  //TODO
  /*
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
  */
  
  //************************************  
  //output
  //************************************
  
  //rad fixups only for critical failure in implicit
  if(radcor>0)     
    fixups[1]=1;
  else
    fixups[1]=0;
    
  if(hdcorr>0)
    corrected[0]=1;
  if(radcor>0)
    corrected[1]=1;

  //if(geom->ix==ixTEST && geom->iy==iyTEST && geom->iz==izTEST) printf("End of u2p_device\n");
  
  return ret;
} 


//**********************************************************************
// solve only for magnetic field 
//**********************************************************************

__device__ __host__ int u2p_solver_Bonly_device(ldouble *uu, ldouble *pp, void *ggg)
{
  //prepare geometry
  struct geometry *geom
  = (struct geometry *) ggg;
  
  ldouble gdetu, gdetu_inv;
  gdetu=geom->gdet;
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
  
  return 0; 
}  


//**********************************************************************
// solver wrapper
//**********************************************************************

__device__ __host__ int u2p_solver_device(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{

#ifdef NONRELMHD
  //return u2p_solver_nonrel_device(uu,pp,ggg,Etype,verbose); //TODO -- u2p_solver nonrel
#endif
    
  int (*solver)(ldouble*,ldouble*,void*,int,int);
    
#if (U2P_SOLVER==U2P_SOLVER_W)  // this is the default
  solver = & u2p_solver_W_device;
#endif

  /*
#if (U2P_SOLVER==U2P_SOLVER_WP)
  solver = & u2p_solver_Wp_device; 
#endif
  */
  
  return (*solver)(uu,pp,ggg,Etype,verbose);
} 


//**********************************************************************
//old Newton-Raphson solver 
//iterates W, not Wp
//Etype == 0 -> hot inversion (uses D,Ttt,Tti)
//Etype == 1 -> entropy inversion (uses D,S,Tti)
//**********************************************************************

__device__ __host__ int u2p_solver_W_device(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{
  //prepare geometry
  struct geometry *geom
  = (struct geometry *) ggg;

  //if(geom->ix==ixTEST && geom->iy==iyTEST && geom->iz==izTEST) printf("In u2p_solver_W_device!\n");
  
  ldouble rho,uint,W,alpha,D,Sc;
  ldouble utcon[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],Qconp[4],Qcovp[4],Qtcon[4],Qtcov[4],Bcon[4],Bcov[4];
  ldouble jmunu[4][4];
  ldouble Qtsq,Qn,QdotB,QdotBsq,Bsq;
  
  
  ldouble pgamma=GAMMA;
  //#ifdef CONSISTENTGAMMA
  //pgamma=pick_gammagas(geom->ix,geom->iy,geom->iz); //TODO
  //#endif
  
  ldouble (*gg)[5], (*GG)[5];
  ldouble gdetu, gdetu_inv;
  gg=geom->gg; GG=geom->GG;
  gdetu=geom->gdet; 
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  gdetu_inv = 1. / gdetu;

  
  //equations choice
  int (*f_u2p)(ldouble,ldouble*,ldouble*,ldouble*,ldouble*,ldouble);
  if(Etype==U2P_HOT)
    f_u2p=&f_u2p_hot;
  if(Etype==U2P_ENTROPY)
    f_u2p=&f_u2p_entropy;
  
  
  //TODO -- print statements
  /*
  if(verbose>1)
  {
    printf("********************\n");
    print_conserved(uu);
    print_primitives(pp);
  }
  */
  
  //conserved quantities etc
  
  //alpha
  alpha = geom->alpha; //alpha=sqrt(-1./GG[0][0]);
  
  //D
  D = uu[0] * gdetu_inv * alpha; //uu[0]=gdetu rho ut
  
  //conserved entropy "S u^t"
  Sc = uu[5] * gdetu_inv * alpha; //uu[5]=gdetu S ut
  
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
  indices_12_device(Qcovp,Qconp,GG);
  
  //Q^mu
  indices_12_device(Qcov,Qcon,GG);
  
#ifdef MAGNFIELD
  //curly B^mu
  Bcon[0]=0.;
  Bcon[1]=uu[B1] * gdetu_inv * alpha;
  Bcon[2]=uu[B2] * gdetu_inv * alpha;
  Bcon[3]=uu[B3] * gdetu_inv * alpha;
  
  //B_mu
  indices_21_device(Bcon,Bcov,gg);
  
  Bsq = dot(Bcon,Bcov);  //NOTE: dot() is a macro
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
  indices_12_device(ncov,ncon,GG);
  
  //Q_mu n^mu = Q^mu n_mu = -alpha*Q^t
  Qn = Qcon[0] * ncov[0];
  
  //j^mu_nu = delta^mu_nu +n^mu n_nu
  for(int i=0;i<4;i++)
  {
    for(int j=0;j<4;j++)
    {
      jmunu[i][j] = delta(i,j) + ncon[i]*ncov[j];
    }
  }
  //Qtilda^nu = j^nu_mu Q^mu
  for(int i=0;i<4;i++)
  {
    Qtcon[i]=0.;
    for(int j=0;j<4;j++)
    {
      Qtcon[i]+=jmunu[i][j]*Qcon[j];
    }
  }
  
  //Qtilda_nu
  indices_21_device(Qtcon,Qtcov,gg);
  
  //Qtsq=Qtilda^mu Qtilda_mu
  Qtsq=dot(Qtcon,Qtcov);
  
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
    conv_vels_device(utcon,utcon,VELPRIM,VELR,gg,GG);
  }
  
  ldouble qsq=0.;
  for(int i=1;i<4;i++)
  {
    for(int j=1;j<4;j++)
    {
      qsq+=utcon[i]*utcon[j]*gg[i][j];
    }
  }

  ldouble gamma2=1.+qsq;
  ldouble gamma=sqrt(gamma2);
  
  //W
  W=(rho+pgamma*uint)*gamma2;
  
  if(verbose>1)
    printf("initial W:%e\n",W);
  
  // test if does not provide reasonable gamma2
  // Make sure that W is large enough so that v^2 < 1 :
  
  ldouble f0,dfdW,err;
  ldouble cons[7]={Qn,Qtsq,D,QdotBsq,Bsq,Sc,Qdotnp};


  // check initial guess of W
  int i_increase = 0;
  do
  {
    f0=dfdW=0.;
    
    (*f_u2p)(W-D,cons,&f0,&dfdW,&err,pgamma);
    
    if( ((( W*W*W * ( W + 2.*Bsq )
           - QdotBsq*(2.*W + Bsq) ) <= W*W*(Qtsq-Bsq*Bsq))
         || !isfinite(f0) || !isfinite(f0)
         || !isfinite(dfdW) || !isfinite(dfdW))
       && (i_increase < 50))
    {
      if(verbose>0) printf("init W : %e -> %e (%e %e)\n",W,10.*W,f0,dfdW);
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
    
    printf("failed to find initial W for Etype: %d\n",Etype);
    //TODO -- print statements
    //printf("at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
    //print_NVvector(uu);
    //print_NVvector(pp);
    return -150;
  }

    
  //1d Newton solver
  int iter=0, fu2pret;
  do
  {
    //if(geom->ix==ixTEST && geom->iy==iyTEST && geom->iz==izTEST) printf("Hello from Newton-Raphson loop\n");
    
    ldouble Wprev=W;
    iter++;

    fu2pret=(*f_u2p)(W-D,cons,&f0,&dfdW,&err,pgamma);
    
    if(verbose>1) printf("%d %e %e %e %e\n",iter,W,f0,dfdW,err);
    
    // convergence test
    if(err<U2PCONV)
      break;

    // nudge if at local max/min
    if(dfdW==0.)
    {
      W*=1.1;
      continue;
    }

    // apply update
    ldouble Wnew = W-f0/dfdW;
    
    
    // test if goes out of bounds and damp solution if so
    int idamp=0;
    ldouble dampfac=1.;
    do
    {
      ldouble f0tmp,dfdWtmp,errtmp;
      f0tmp=dfdWtmp=0.;
      (*f_u2p)(Wnew-D,cons,&f0tmp,&dfdWtmp,&errtmp,pgamma);
      if(verbose>1)
	printf("sub (%d) :%d %e %e %e %e\n",idamp,iter,Wnew,f0tmp,dfdWtmp,errtmp);
      if( ((( Wnew*Wnew*Wnew * ( Wnew + 2.*Bsq )
             - QdotBsq*(2.*Wnew + Bsq) ) <= Wnew*Wnew*(Qtsq-Bsq*Bsq))
           || !isfinite(f0tmp) || !isfinite(f0tmp)
           || !isfinite(dfdWtmp) || !isfinite(dfdWtmp))
         && (idamp<100))
      {
        idamp++;
        dampfac/=2.;
        Wnew=W-dampfac*f0/dfdW; // damp to bring in bounds
        continue;
      }
      else
        break;
    }
    while(1);
    
    if(idamp>=100)
    {
      if(verbose>0) printf("damped unsuccessfuly\n");
      return -101;
    }

    // apply update after damping
    W=Wnew;

    // check if W is too large
    if(fabs(W)>BIG)
    {
      //TODO -- print statement
      //if(verbose>1) printf("W has gone out of bounds at %d,%d,%d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz);
      return -103;
    }

    // check for convergence
    if(fabs((W-Wprev)/Wprev)<U2PCONV && err<ERRCONV)
      break;
  }
  while(iter<50);
  
  
  if(iter>=50)
  {
    if(verbose>0) printf("iter exceeded in u2p_solver with Etype: %d\n",Etype); 
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
    printf("end: %d %e %e %e %e %d \n",iter,W,f0,dfdW,err,fu2pret);
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
    if(verbose>0)
      printf("neg u in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);
    return -104;
  }
  
  //converting to VELPRIM
  if (VELR != VELPRIM)
  {
    conv_vels_device(utcon,utcon,VELR,VELPRIM,gg,GG);
  }
  
  //returning new primitives
  pp[RHO]=rho;
  pp[UU]=uint;
  pp[VX]=utcon[1];
  pp[VY]=utcon[2];
  pp[VZ]=utcon[3];
  
  if(rho<0.)
  {
    if(verbose>0)
      printf("neg rho in u2p_solver %e %e %e %e\n",rho,uint,gamma2,W);
    return -105;
  }
  
  //entropy based on Etype  
  //pure entropy evolution - updated only in the end of RK2
  pp[ENTR]=entr; // TODO why is this entr instead of Sc * gdetu_inv / utcon[0]? 
  
#ifdef MAGNFIELD
  //magnetic conserved=primitives
  pp[B1]=uu[B1] * gdetu_inv;
  pp[B2]=uu[B2] * gdetu_inv;
  pp[B3]=uu[B3] * gdetu_inv;
#endif
  
#ifdef EVOLVEELECTRONS
  conv_vels_device(utcon,utcon,VELPRIM,VEL4,gg,GG);
  ldouble Se=uu[ENTRE] * gdetu_inv / utcon[0];
  pp[ENTRE]=Se;
  ldouble Si=uu[ENTRI] * gdetu_inv / utcon[0];
  pp[ENTRI]=Si;
  
#ifdef RELELECTRONS
  for(int ib=0;ib<NRELBIN;ib++)
    pp[NEREL(ib)]=uu[NEREL(ib)] * gdetu_inv / utcon[0];
#endif
#endif

  //TODO -- print statement
  //if(verbose) print_primitives(pp);
  
  if(verbose>0)
    printf("u2p_solver returns 0\n");

  //if(geom->ix==ixTEST && geom->iy==iyTEST && geom->iz==izTEST) printf("End of u2p_solver_W_device\n");
  
  return 0; 
}


//********************************************
//Harm u2p_hot
//********************************************

__device__ __host__ static ldouble dpdWp_calc_vsq(ldouble Wp, ldouble D, ldouble vsq, ldouble gamma)
{
  ldouble W=Wp+D;
  return( (gamma - 1.) * (1. - vsq) /  gamma ) ;
}

// 1 / (d(u+p)/dp)
__device__ __host__ static ldouble compute_idwmrho0dp(ldouble wmrho0, ldouble gamma)
{
  return((gamma-1.)/gamma);
}


// 1 / (drho0/dp) holding wmrho0 fixed
__device__ __host__ static ldouble compute_idrho0dp(ldouble wmrho0)
{
  return(0.0);
}

__device__ __host__ static int f_u2p_hot(ldouble Wp, ldouble* cons,ldouble *f,ldouble *df,
					 ldouble *err,ldouble pgamma)
{

  //printf("hi from f_u2p_hot\n");
  ldouble Qn=cons[0];
  ldouble Qtsq=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  //ldouble Qdotnp=cons[6]; //not used unles U2P_EQS_JON
  
  ldouble W=Wp+D;
  ldouble X = Bsq + W;
  ldouble Wsq = W*W;
  ldouble W3 = Wsq*W ;
  ldouble Xsq = X*X;
  ldouble X3 = Xsq*X;

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
 ldouble Qdotnp=cons[6];
   *f = Qdotnp + Wp - p + 0.5*Bsq + (Bsq*Qtsq - QdotBsq)/Xsq;
 *err = fabs(*f) / (fabs(Qdotnp) + fabs(Wp) + fabs(p) + fabs(0.5*Bsq) + fabs((Bsq*Qtsq - QdotBsq)/Xsq));
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

// specific entropy as function of rho0 and internal energy (u)
// Ss(rho0,\chi=u+p)
// specific entropy = \ln( p^n/\rho^{n+1} )
__device__ __host__ static ldouble compute_specificentropy_wmrho0_idealgas(ldouble rho0, ldouble wmrho0,ldouble gamma)
{
  ldouble insideentropy,specificentropy;
  ldouble pressure,indexn;

  pressure=((gamma-1.)/gamma)*wmrho0;
  indexn=1.0/(gamma-1.);  
  insideentropy=pow(pressure,indexn)/pow(rho0,indexn+1.0);
  specificentropy=log(insideentropy);

  return(specificentropy);

}

// used for utoprim_jon when doing entropy evolution
// dSspecific/d\chi
__device__ __host__ static ldouble compute_dspecificSdwmrho0_wmrho0_idealgas(ldouble rho0, ldouble wmrho0,ldouble gamma)
{
  ldouble dSdchi;

  dSdchi = 1.0/((gamma-1.)*wmrho0);

  return(dSdchi);

}

// dSspecific/drho0
__device__ __host__ static ldouble compute_dspecificSdrho_wmrho0_idealgas(ldouble rho0, ldouble wmrho0, ldouble gamma)
{
  ldouble dSdrho;
  
  dSdrho=gamma/((1.0-gamma)*rho0);

  return(dSdrho);
}


__device__ __host__ static int f_u2p_entropy(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df,
					     ldouble *err,ldouble pgamma)
{
  //ldouble Qn=cons[0]; //not used
  ldouble Qtsq=cons[1];
  ldouble D=cons[2];
  ldouble QdotBsq=cons[3];
  ldouble Bsq=cons[4];
  ldouble Sc=cons[5];
 
  ldouble W=Wp+D;
  ldouble X = Bsq + W;
  ldouble Wsq = W*W;
  ldouble W3 = Wsq*W ;
  ldouble Xsq = X*X;
  ldouble X3 = Xsq*X;

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

  ldouble dSsdW,dSsdvsq,dSsdWp,dScprimedWp,dSsdrho,dSsdchi;
  ldouble dvsq,dwmrho0dW,drho0dW;
  ldouble dwmrho0dvsq,drho0dvsq;

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
//checks if hydro primitives make sense
//**********************************************************************

// TODO deleted some floors
// deleted RHOFLOOR_BH, RHOFLOORINIT
// deleted EVOLVEELECTRONS corrections in B2RHO ceiling
// deleted EVOLVEELECTRONS, RELELECTRONS floors/ceilings
__device__ __host__ int check_floors_mhd_device(ldouble *pp, int whichvel,void *ggg)
{

  int verbose=0;
  int ret=0;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
 
  //**********************************************************************
  //rho too small -- set equal to floor valu
  if(pp[0] < RHOFLOOR) 
  {

      //TODO print with MPI
      //if(verbose ) printf("hd_floors CASE 1 at %d %d %d | %d %d %d (%e) | tijk: %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,geom->ix,geom->iy,geom->iz,pp[0],TI,TJ,TK);

      pp[0]=RHOFLOOR; 
     
      ret=-1; 
  }

  //**********************************************************************
  //too cold -- raise gas energy to minimum density fraction
  if(pp[1] < UURHORATIOMIN*pp[0]) 
  {

    //TODO print with MPI
    //if(verbose) {printf("hd_floors CASE 2 at (%d,%d,%d | %d,%d,%d): %e %e | tijk: %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,geom->ix,geom->iy,geom->iz,pp[0],pp[1],TI,TJ,TK);}
    
    pp[1]=UURHORATIOMIN*pp[0]; 

    ret=-1;
  }

  //**********************************************************************
  //too hot -- lower gas energy to maximum density fraction 
  if(pp[1]>UURHORATIOMAX*pp[0]) 
  {
    //TODO print with MPI
    //if(verbose ) printf("hd_floors CASE 3 at (%d,%d,%d | %d,%d,%d): %e %e | tijk: %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,geom->ix,geom->iy,geom->iz,pp[0],pp[1],TI,TJ,TK);
 
    pp[1]=UURHORATIOMAX*pp[0]; 

    ret=-1;      
    
  }
  
  //**********************************************************************
  //too magnetized
  
#ifdef MAGNFIELD
  ldouble ucond[4],ucovd[4];
  ldouble bcond[4],bcovd[4],bsq,magpre;
  ldouble etacon[4],etarel[4];

  // get magnetic energy
  for(int iv=1;iv<4;iv++)
    ucond[iv]=pp[1+iv];
  calc_ucon_ucov_from_prims_device(pp, geom, ucond, ucovd);
  calc_bcon_bcov_bsq_from_4vel_device(pp, ucond, ucovd, geom, bcond, bcovd, &bsq);
  magpre = 0.5 * bsq;
  
  calc_normalobs_ncon_device(GG, geom->alpha, etacon);
  conv_vels_ut_device(etacon,etarel,VEL4,VELPRIM,gg,GG);

  if(magpre > B2RHORATIOMAX*pp[RHO]) 
  {
    //TODO -- print with MPI
    //if(verbose) printf("mag_floor at (%d,%d,%d): %e %e\n",geom->ix+TOI,geom->iy+TOJ,geom->iz,pp[RHO],magpre);

    ldouble f=magpre/(B2RHORATIOMAX*pp[RHO]); //correction factor

#ifdef B2RHOFLOOR_BACKUP_FFFRAME
    ldouble pporg[NV];
    for(int iv=0;iv<NVMHD;iv++)
    {
      pporg[iv]=pp[iv];
    }
#endif
    
    ldouble uu[NV];
    p2u_mhd_device(pp,uu,ggg); 

#if (B2RHOFLOORFRAME==ZAMOFRAME) // add new mass in ZAMO by default

    ldouble dpp[NV],duu[NV]; 
   
    for(int iv=0;iv<NVMHD;iv++)
       dpp[iv]=0.;

    // compute delta p
    dpp[RHO]=pp[RHO]*(f-1.);;
    dpp[UU]=0.; //do not inject energy - just density
    dpp[VX] = etarel[1]; //new mass moves at ZAMO velocity
    dpp[VY] = etarel[2];
    dpp[VZ] = etarel[3];
    dpp[ENTR] = 0.;
    dpp[B1] = dpp[B2] = dpp[B3] = 0.;

    // translate to conserved and update
    p2u_mhd_device(dpp,duu,geom); 
    for(int iv=0;iv<NVMHD;iv++)
    {
      uu[iv]+=duu[iv];
    }

    // perform a u2p on the updated conserved
    int rettemp=0;
    rettemp=u2p_solver_device(uu,pp,geom,U2P_HOT,0); 
    if(rettemp<0)
       rettemp=u2p_solver_device(uu,pp,geom,U2P_ENTROPY,0); 
      
    if(rettemp<0) 
    {
      //TODO print with MPI
      //printf("u2p failed after imposing bsq over rho floors at %d %d %d with gamma=%f\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,get_u_scalar(gammagas,geom->ix,geom->iy,geom->iz));

#ifdef B2RHOFLOOR_BACKUP_FFFRAME
      // Backup bsq/rho floor -- if zamo frame fails, do fluid frame instead of crashing 
      for(int iv=0;iv<NVMHD;iv++)
         pp[iv]=pporg[iv];
      pp[RHO]*=f; // inject both energy and density (TODO is this correct?)
      pp[UU]*=f;
#else
      //TODO print statement
      //print_primitives(pp);
      // TODO: note you shouldn't call __host__ exit(...) from __host__ __device__
      //exit(-1);
#endif
    }
    
#elif(B2RHOFLOORFRAME==FFFRAME) // add new mass in fluid frame
    pp[RHO]*=f;
    pp[UU]*=f;

#endif //B2RHOFLOORFRAME==ZAMOFRAME
    
    ret=-1;      
  } //if(magpre>B2RHORATIOMAX*pp[RHO]) 
#endif //MAGNFIELD

  //**********************************************************************
  //too fast
  if(VELPRIM==VELR) 
  {
    ldouble qsq=0.;
    for(int i=1;i<4;i++)
    {
      for(int j=1;j<4;j++)
      {
	qsq+=pp[UU+i]*pp[UU+j]*gg[i][j];
      }
    }

    ldouble gamma2=1.+qsq; // squared lorentz factor
    if(gamma2>GAMMAMAXHD*GAMMAMAXHD)
    {
      ldouble qsqmax=GAMMAMAXHD*GAMMAMAXHD-1.;
      ldouble A=sqrt(qsqmax/qsq);
      for(int j=1;j<4;j++)
        pp[UU+j]*=A;

      ret=-1;

      if(verbose )
      {
        //TODO print with MPI
	//printf("hd_floors CASE 4 at (%d,%d,%d): %e",geom->ix+TOI,geom->iy+TOJ,geom->iz,sqrt(gamma2));
	//ldouble qsqnew=0.;
	//for(int i=1;i<4;i++)
	//   for(int j=1;j<4;j++)
        //       qsqnew+=pp[UU+i]*pp[UU+j]*gg[i][j];
	//printf(" -> %e\n",sqrt(1+qsq));
      }
    }
  }
  
  //TODO do we want this? Is this inconsistent with keeping entropy as a backup until the very end of time step? 
  //updates entropy after floor corrections
  if(ret<0)
    pp[5] = calc_Sfromu_device(pp[RHO],pp[UU],geom->ix,geom->iy,geom->iz);

  return ret;
}

//**********************************************************************
// u2p kernel
//**********************************************************************

__global__ void calc_primitives_kernel(int Nloop_0, 
				       int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
                                       ldouble *x_arr, ldouble *g_arr, ldouble *G_arr,
				       ldouble *u_arr, ldouble *p_arr,
				       int setflags, int* cellflag_arr, int int_slot_arr[NGLOBALINTSLOT])
				       
{
  
  // get index for this thread
  // Nloop_0 is number of cells to update;
  // usually Nloop_0=NX*NY*NZ, but sometimes there are weird bcs inside domain 
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_0) return;
    
  // get indices from 1D arrays
  int ix=loop_0_ix[ii];
  int iy=loop_0_iy[ii];
  int iz=loop_0_iz[ii]; 
  
  //skip if cell is passive
  if(!is_cell_active_device(ix,iy,iz))
    return;
  
  struct geometry geom;
  fill_geometry_device(ix,iy,iz,&geom,x_arr,g_arr, G_arr);

  ldouble uu[NV],pp[NV];
  
  int corrected[3]={0,0,0}, fixups[2]={0,0};
  for(int iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u_arr,iv,ix,iy,iz);
    pp[iv]=get_u(p_arr,iv,ix,iy,iz);
  }
 
  if(setflags)
  {
    set_cflag(cellflag_arr,ENTROPYFLAG,ix,iy,iz,0);
    set_cflag(cellflag_arr,ENTROPYFLAG2,ix,iy,iz,0);
  }
  
  //u to p inversion is done here
  if(is_cell_corrected_polaraxis_device(ix,iy,iz)) 
  {
    // invert only the magnetic field, the rest will be overwritten
    u2p_solver_Bonly_device(uu,pp,&geom); 
  }
  else
  {
    // regular inversion
    u2p_device(uu,pp,&geom,corrected,fixups,int_slot_arr); 
  }
  
  //set flags for entropy solver
  if(corrected[0]==1 && setflags) //hd correction - entropy solver
  {
    set_cflag(cellflag_arr,ENTROPYFLAG,ix,iy,iz,1);
  }
  
  if(corrected[2]==1 && setflags) //borrowing energy from radiation didn't work
  {  
    set_cflag(cellflag_arr,ENTROPYFLAG2,ix,iy,iz,1);
  }

#ifndef NOFLOORS
  //check hd floors
  int floorret=0;

  if(is_cell_active_device(ix,iy,iz) && !is_cell_corrected_polaraxis_device(ix,iy,iz))
  {
    floorret=check_floors_mhd_device(pp,VELPRIM,&geom);
  }
  
  if(floorret<0.)
  {
    corrected[0]=1;
  }

  //TODO 
  /*
  //check rad floors
#ifdef RADIATION
  floorret=0;
  if(is_cell_active_device(ix,iy,iz) &&  !is_cell_corrected_polaraxis_device(ix,iy,iz))
  {
    floorret=check_floors_rad(pp,VELPRIMRAD,&geom);
  }
  
  if(floorret<0.)
  {
    corrected[1]=1;
  }
#endif
  */
#endif
  
  //set new primitives
  for(int iv=0;iv<NV;iv++)
  { 
    set_u(p_arr,iv,ix,iy,iz,pp[iv]);
  }
  
  //set flags for fixups of unsuccessful cells
  if(setflags)
  {
    if(fixups[0]>0)
    {
      set_cflag(cellflag_arr,HDFIXUPFLAG,ix,iy,iz,1);
      atomicAdd(&int_slot_arr[GLOBALINTSLOT_NTOTALMHDFIXUPS],1); //TODO right??
      //global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]++; 
    }
    else
      set_cflag(cellflag_arr,HDFIXUPFLAG,ix,iy,iz,0);
    
    if(fixups[1]>0)
    {
      set_cflag(cellflag_arr,RADFIXUPFLAG,ix,iy,iz,-1);
      atomicAdd(&int_slot_arr[GLOBALINTSLOT_NTOTALRADFIXUPS],1); //TODO right??
      //global_int_slot[GLOBALINTSLOT_NTOTALRADFIXUPS]++;
    }
    else
      set_cflag(cellflag_arr,RADFIXUPFLAG,ix,iy,iz,0); 
  }
} 


__global__ void cell_fixup_kernel(int Nloop_0, 
				  int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
				  ldouble* x_arr, ldouble* g_arr, ldouble* G_arr,
				  int* cellflag_arr, int type,
				  ldouble* u_arr, ldouble* p_arr,
				  ldouble* u_bak_fixup_arr, ldouble* p_bak_fixup_arr)

{

  int verbose=1;

  // get index for this thread
  // Nloop_0 is number of cells to update;
  // usually Nloop_0=NX*NY*NZ, but sometimes there are weird bcs inside domain 
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_0) return;
    
  // get indices from 1D arrays
  int ix=loop_0_ix[ii];
  int iy=loop_0_iy[ii];
  int iz=loop_0_iz[ii]; 

  // copy current primitives/conserved to output arrays
  // TODO could use up_copy_kernel
  for(int iv=0;iv<NV;iv++)
  {
    set_u(u_bak_fixup_arr,iv,ix,iy,iz,get_u(u_arr,iv,ix,iy,iz));
    set_u(p_bak_fixup_arr,iv,ix,iy,iz,get_u(p_arr,iv,ix,iy,iz));
  }
  
  //do not correct if overwritten later on
  if(is_cell_corrected_polaraxis_device(ix,iy,iz))
    return;

  // check if we want to fixup
  if(((get_cflag(cellflag_arr,HDFIXUPFLAG,ix,iy,iz)!=0 && type==FIXUP_U2PMHD) ||
      (get_cflag(cellflag_arr,RADFIXUPFLAG,ix,iy,iz)!=0 && type==FIXUP_U2PRAD) ||
      (get_cflag(cellflag_arr,RADIMPFIXUPFLAG,ix,iy,iz)!=0 && type==FIXUP_RADIMP)) &&
       is_cell_active_device(ix,iy,iz))
  {

    // fill geometry
    struct geometry geom;
    fill_geometry_device(ix,iy,iz,&geom, x_arr,g_arr,G_arr);

    //looking around for good neighbors
    ldouble ppn[6][NV],pp[NV],uu[NV];

    //should care about global but at the stage where it is called knowns not about the boundaries
    //so fixups idividually in tiles but later exchanged 

    int in=0; //number of successful neighbors

    //********************//
    // x neighbors
    //********************//

    if(ix-1>=0 &&  
       ((get_cflag(cellflag_arr,HDFIXUPFLAG,ix-1,iy,iz)==0 && type==FIXUP_U2PMHD) ||
	(get_cflag(cellflag_arr,RADFIXUPFLAG,ix-1,iy,iz)==0 && type==FIXUP_U2PRAD) ||
	(get_cflag(cellflag_arr,RADIMPFIXUPFLAG,ix-1,iy,iz)==0 && type==FIXUP_RADIMP)))	      
    {
      // TODO -- MPI variables
      #ifdef SPECIAL_BC_CHECK //make sure that ix-1 is not a stream cell
      if(TNY>1 && TNZ==1)
      {
        if((iy+TOJ) >= STREAM_IYT && (iy+TOJ) <= STREAM_IYB && (ix+TOI) > (STREAM_IX+1))
        {
          if((ix-1+TOI) > (STREAM_IX+3)) //not a stream cell in ix-1
          {
            in++;
            for(int iv=0;iv<NV;iv++)
              ppn[in-1][iv]=get_u(p_arr,iv,ix-1,iy,iz);
          }
        }
        else
        {
          in++;
          for(int iv=0;iv<NV;iv++)
            ppn[in-1][iv]=get_u(p_arr,iv,ix-1,iy,iz);
        }
      }
      #else

      in++;
      for(int iv=0;iv<NV;iv++)
        ppn[in-1][iv]=get_u(p_arr,iv,ix-1,iy,iz);
      #endif
    }

    if(ix+1<NX &&
      ((get_cflag(cellflag_arr,HDFIXUPFLAG,ix+1,iy,iz)==0 && type==FIXUP_U2PMHD) ||
       (get_cflag(cellflag_arr,RADFIXUPFLAG,ix+1,iy,iz)==0 && type==FIXUP_U2PRAD) ||
       (get_cflag(cellflag_arr,RADIMPFIXUPFLAG,ix+1,iy,iz)==0 && type==FIXUP_RADIMP)))
    {
      // TODO -- MPI variables
      #ifdef SPECIAL_BC_CHECK //make sure that ix-1 is not a stream ghost cell
      if(TNY>1 && TNZ==1)
      {
	if((iy+TOJ) >= STREAM_IYT && (iy+TOJ) <= STREAM_IYB && (ix+TOI) < STREAM_IX)
	{
	  if((ix+1+TOI) < STREAM_IX) //not a stream cell in ix+1
	  {
	    in++;
	    for(int iv=0;iv<NV;iv++)
	      ppn[in-1][iv]=get_u(p_arr,iv,ix+1,iy,iz);
	  }
	}
	else //not in danger of using a stream ghost to do a fixup
	{
	  in++;
	  for(int iv=0;iv<NV;iv++)
	    ppn[in-1][iv]=get_u(p_arr,iv,ix+1,iy,iz);
	}
      }
      #else
      in++;
      for(int iv=0;iv<NV;iv++)
	ppn[in-1][iv]=get_u(p_arr,iv,ix+1,iy,iz);
      #endif
    }

    //********************//
    // y neighbors
    //********************//
    if(iy-1>=0 &&
      ((get_cflag(cellflag_arr,HDFIXUPFLAG,ix,iy-1,iz)==0 && type==FIXUP_U2PMHD) ||
       (get_cflag(cellflag_arr,RADFIXUPFLAG,ix,iy-1,iz)==0 && type==FIXUP_U2PRAD) ||
       (get_cflag(cellflag_arr,RADIMPFIXUPFLAG,ix,iy-1,iz)==0 && type==FIXUP_RADIMP)))
    {
      // TODO -- MPI variables
      #ifdef SPECIAL_BC_CHECK //make sure that ix-1 is not a stream ghost cell
      if(TNY>1 && TNZ==1)
      {
	if((ix+TOI) >= STREAM_IX && (ix+TOI) <= (STREAM_IX+3) && (iy+TOJ) > STREAM_IYB)
	{
	  if((iy-1+TOJ) > STREAM_IYB) //not a stream cell in iy-1
	  {
	    in++;
	    for(int iv=0;iv<NV;iv++)
	      ppn[in-1][iv]=get_u(p_arr,iv,ix,iy-1,iz);
	  }
	}
	else //not in danger of using a stream ghost to do a fixup
	{
	  in++;
	  for(int iv=0;iv<NV;iv++)
	    ppn[in-1][iv]=get_u(p_arr,iv,ix,iy-1,iz);
	}
      }
      #else
      in++;
      for(int iv=0;iv<NV;iv++)
	ppn[in-1][iv]=get_u(p_arr,iv,ix,iy-1,iz);
      #endif
    }

    if(iy+1<NY &&
       ((get_cflag(cellflag_arr,HDFIXUPFLAG,ix,iy+1,iz)==0 && type==FIXUP_U2PMHD) ||
	(get_cflag(cellflag_arr,RADFIXUPFLAG,ix,iy+1,iz)==0 && type==FIXUP_U2PRAD) ||
	(get_cflag(cellflag_arr,RADIMPFIXUPFLAG,ix,iy+1,iz)==0 && type==FIXUP_RADIMP)))
    {
      // TODO -- MPI variables
      #ifdef SPECIAL_BC_CHECK //make sure that iy+1 is not a stream ghost cell
      if(TNY>1 && TNZ==1)
      {
	if((ix+TOI) >= STREAM_IX && (ix+TOI) <= (STREAM_IX+3) && (iy+TOJ) < STREAM_IYT)
	{
	  if((iy+1+TOJ) < STREAM_IYT) //not a stream cell in iy+1
	  {
	    in++;
	    for(int iv=0;iv<NV;iv++)
	      ppn[in-1][iv]=get_u(p_arr,iv,ix,iy+1,iz);
	  }
	}
	else //not in danger of using a stream ghost to do a fixup
	{
	  in++;
	  for(int iv=0;iv<NV;iv++)
	    ppn[in-1][iv]=get_u(p_arr,iv,ix,iy+1,iz);
	}
      }
      #else
      in++;
      for(int iv=0;iv<NV;iv++)
	ppn[in-1][iv]=get_u(p_arr,iv,ix,iy+1,iz);
      #endif
    }


    //********************//
    // z neighbors
    //********************//

    if(iz-1>=0 &&
       ((get_cflag(cellflag_arr,HDFIXUPFLAG,ix,iy,iz-1)==0 && type==FIXUP_U2PMHD) ||
	(get_cflag(cellflag_arr,RADFIXUPFLAG,ix,iy,iz-1)==0 && type==FIXUP_U2PRAD) ||
	(get_cflag(cellflag_arr,RADIMPFIXUPFLAG,ix,iy,iz-1)==0 && type==FIXUP_RADIMP)))
    {
      in++;
      for(int iv=0;iv<NV;iv++)
        ppn[in-1][iv]=get_u(p_arr,iv,ix,iy,iz-1);
    }

    if(iz+1<NZ  &&
       ((get_cflag(cellflag_arr,HDFIXUPFLAG,ix,iy,iz+1)==0 && type==FIXUP_U2PMHD) ||
	(get_cflag(cellflag_arr,RADFIXUPFLAG,ix,iy,iz+1)==0 && type==FIXUP_U2PRAD) ||
	(get_cflag(cellflag_arr,RADIMPFIXUPFLAG,ix,iy,iz+1)==0 && type==FIXUP_RADIMP)))
    {
      in++;
      for(int iv=0;iv<NV;iv++)
        ppn[in-1][iv]=get_u(p_arr,iv,ix,iy,iz+1);
    }

    // TODO do we want to set all flags to zero when we are done?
    // would need a separate kernel after synchronization

    //*******************************************//
    // do fixup if there are sufficient neighbors
    //*******************************************//

    if((NZ==1 && NY==1 && in>=1) ||
       (NZ==1 && in>=2) || (NY==1 && in>=2) ||
       (in>=3))
    {
      for(int iv=0;iv<NV;iv++)
      {
        int fixthis=0;

	if(type==FIXUP_U2PMHD &&  iv!=RHO && iv<B1) 
	  fixthis=1; //skip correcting magnetic field so we don't disrupt div B

	if(type==FIXUP_U2PRAD &&  iv!=RHO && iv>=EE && iv<=FZ) 
	  fixthis=1; //fix only radiative quantites

	if(type==FIXUP_RADIMP &&  iv!=RHO && (iv<B1 || (iv>=EE && iv<=FZ))) 
	  fixthis=1; //fix both mhd and rad but not magn field

	#ifdef EVOLVEPHOTONNUMBER
	if(type==FIXUP_RADIMP &&  iv==NF) 
	  fixthis=1; //fix number of photons
	#endif

	#ifdef EVOLVEELECTRONS
	if(type==FIXUP_RADIMP &&  (iv==ENTRE || iv==ENTRI)) 
	  fixthis=1; //fix electron/ion entropy
	#endif

	// do the fixup
	if(fixthis==1)
	{
	  pp[iv]=0.;
	  for(int iii=0;iii<in;iii++)
	    pp[iv]+=ppn[iii][iv];
	  pp[iv]/=(ldouble)in;  
	}
	else //leave as was
	{
	  pp[iv]=get_u(p_arr,iv,ix,iy,iz); 
	}
      }

      // do p2u to get consistent u 
      p2u_device(pp,uu,&geom);

      // TODO MPI prints
      if(verbose>1) 
      {
	//if(type==FIXUP_U2PMHD) printf("%4d > %4d %4d %4d > U2PMHD > fixing up with %d neighbors\n",PROCID,ix+TOI,iy+TOJ,iz+TOK,in);
	//if(type==FIXUP_U2PRAD) printf("%4d > %4d %4d %4d > U2PRAD > fixing up with %d neighbors\n",PROCID,ix+TOI,iy+TOJ,iz+TOK,in);
	//if(type==FIXUP_RADIMP) printf("%4d > %4d %4d %4d > RADIMP > fixing up with %d neighbors\n",PROCID,ix+TOI,iy+TOJ,iz+TOK,in);
      }

      //save to updated arrays memory
      for(int iv=0;iv<NV;iv++)
      {
	set_u(u_bak_fixup_arr,iv,ix,iy,iz,uu[iv]);
	set_u(p_bak_fixup_arr,iv,ix,iy,iz,pp[iv]);
      }
    }
    else
    {
      //TODO MPI prints
      //if(type==FIXUP_U2PMHD) printf("%4d > %4d %4d %4d > U2PMHD > didn't manage to hd fixup\n",PROCID,ix+TOI,iy+TOJ,iz+TOK);
      //if(type==FIXUP_U2PRAD) printf("%4d > %4d %4d %4d > U2PRAD > didn't manage to hd fixup\n",PROCID,ix+TOI,iy+TOJ,iz+TOK);
      //if(type==FIXUP_RADIMP) printf("%4d > %4d %4d %4d > RADIMP > didn't manage to hd fixup\n",PROCID,ix+TOI,iy+TOJ,iz+TOK);
    }
  }
}

//**********************************************************************
//Call the kernels in u2p
//**********************************************************************

ldouble calc_u2p_gpu(int setflags)
{

  cudaError_t err = cudaSuccess;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // NOTE: size of xb is copied from initial malloc in misc.c
  // these need to be long long if the grid is on one tile and large (~256^3)
  long long Ncellflag = (SX)*(SY)*(SZ)*NFLAGS;
  long long Nprim     = (SX)*(SY)*(SZ)*NV;

  err = cudaMemcpy(d_cellflag_arr, cellflag, sizeof(int)*Ncellflag, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_int_slot_arr, global_int_slot, sizeof(int)*NGLOBALINTSLOT, cudaMemcpyHostToDevice);

  // launch calc_primitives_kernel

  int threadblocks = (Nloop_0 / TB_SIZE) + ((Nloop_0 % TB_SIZE)? 1:0);

  cudaEventRecord(start);
  calc_primitives_kernel<<<threadblocks, TB_SIZE>>>(Nloop_0, 
						    d_loop0_ix, d_loop0_iy, d_loop0_iz,
						    d_x, d_gcov, d_gcon,
						    d_u_arr, d_p_arr,
						    setflags, d_cellflag_arr, d_int_slot_arr);
              
  cudaEventRecord(stop);
  err = cudaPeekAtLastError(); // TODO ?? 
  cudaDeviceSynchronize(); 

  cudaEventSynchronize(stop);
  float tms = 0.;
  cudaEventElapsedTime(&tms, start,stop);
  printf("gpu u2p time: %0.2f \n",tms);

#ifdef CPUKO
  // print results
  ldouble* p_tmp;
  if((p_tmp=(ldouble*)malloc(sizeof(ldouble)*Nprim))==NULL) my_err("malloc err.\n");
  err = cudaMemcpy(p_tmp, d_p_arr, sizeof(ldouble)*Nprim, cudaMemcpyDeviceToHost);
  if(err != cudaSuccess) printf("failed cudaMemcpy of d_p_arr to p_tmp\n");
  printf("gpu u2p pp[NV]: ");
  for(int iv=0;iv<NV;iv++)
    printf("%e ", get_u(p_tmp, iv, ixTEST, iyTEST, izTEST));
  printf("\n");
  free(p_tmp);
#endif

  // fixups

  int do_fixups=1;
  int type=FIXUP_U2PMHD;
  if(DOFIXUPS==0)
    do_fixups=0;

  if(DOU2PMHDFIXUPS==0 && type==FIXUP_U2PMHD)
    do_fixups=0;

  if(DOU2PRADFIXUPS==0 && type==FIXUP_U2PRAD)
    do_fixups=0;

  if(DORADIMPFIXUPS==0 && type==FIXUP_RADIMP)
    do_fixups=0;

  // TODO need a test for fixups
  if(do_fixups)
  {
    // launch cell_fixup kernel 
    cell_fixup_kernel<<<threadblocks, TB_SIZE>>>(Nloop_0, 
						 d_loop0_ix, d_loop0_iy, d_loop0_iz,
						 d_x, d_gcov, d_gcon,
						 d_cellflag_arr, type, 
						 d_u_arr, d_p_arr,
						 d_u_fixup_arr, d_p_fixup_arr);
    err = cudaPeekAtLastError(); // TODO ?? 
    cudaDeviceSynchronize(); 

  #ifdef RADIATION
    //TODO
    cell_fixup_kernel(FIXUP_U2PRAD);
  #endif

    // copy fixed up p&u to p,u
    cudaEventRecord(start);
    up_copy_kernel<<<threadblocks, TB_SIZE>>>(Nloop_0,
					      d_loop0_ix, d_loop0_iy, d_loop0_iz,
					      d_u_fixup_arr, d_p_fixup_arr,
                                              d_u_arr, d_p_arr);
    cudaEventRecord(stop);
    err = cudaPeekAtLastError(); // TODO ?? 
    cudaDeviceSynchronize(); 


    cudaEventSynchronize(stop);
    float tms = 0.;
    cudaEventElapsedTime(&tms, start,stop);
    printf("gpu u2p fixup time: %0.2f \n",tms);

  #ifdef CPUKO
    // print results
    if((p_tmp=(ldouble*)malloc(sizeof(ldouble)*Nprim))==NULL) my_err("malloc err.\n");
    err = cudaMemcpy(p_tmp, d_p_arr, sizeof(ldouble)*Nprim, cudaMemcpyDeviceToHost);
    if(err != cudaSuccess) printf("failed cudaMemcpy of d_p_arr to p_tmp\n");
    printf("gpu u2p fixup pp[NV]: ");
    for(int iv=0;iv<NV;iv++)
      printf("%e ", get_u(p_tmp, iv, ixTEST, iyTEST, izTEST));
    printf("\n");
    free(p_tmp);
  #endif
  }

  return (ldouble)tms;
}



