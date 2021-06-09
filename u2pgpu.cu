/*! \file u2pgpu.cu
 \brief MHD Conserved to primitives conversion
 */

extern "C" {

#include "ko.h"

}

#include "kogpu.h"

// reminder of default settings
//#define U2PCONV 1.e-10
//#define U2P_EQS U2P_EQS_NOBLE
//#define U2P_SOLVER U2P_SOLVER_W

__device__ __host__ static FTYPE dpdWp_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq,FTYPE gamma);
__device__ __host__ static FTYPE compute_idwmrho0dp(FTYPE wmrho0,FTYPE gamma);
__device__ __host__ static FTYPE compute_idrho0dp(FTYPE wmrho0);
__device__ __host__ static int f_u2p_hot(ldouble Wp, ldouble* cons,ldouble *f,ldouble *df,ldouble *err,ldouble pgamma);
__device__ __host__ static FTYPE pressure_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma);
__device__ __host__ static FTYPE compute_inside_entropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma);
__device__ __host__ static FTYPE compute_specificentropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma);
__device__ __host__ static FTYPE compute_dspecificSdwmrho0_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma);
__device__ __host__ static FTYPE compute_dspecificSdrho_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0, FTYPE gamma);
__device__ __host__ static int f_u2p_entropy(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df, ldouble *err,ldouble pgamma);

// todo deleted type
__global__ void calc_primitives_kernel(int Nloop_0, int setflags,
				       int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
				       ldouble *u_arr, ldouble *p_arr,
				       ldouble *x_arr, ldouble *g_arr, ldouble *G_arr)
{

  int verbose=0;
  
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
  //fill_geometry(ix,iy,iz,&geom);
  fill_geometry_device(ix,iy,iz, x_arr,&geom,g_arr, G_arr);

   
  int u2pret,u2pretav;
  ldouble uu[NV],pp[NV];
  
  int corrected[3]={0,0,0}, fixups[2]={0,0};
  for(int iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u_arr,iv,ix,iy,iz);
    pp[iv]=get_u(p_arr,iv,ix,iy,iz);
  }

  //TODO -- put in flags
  /*
  if(setflags)
  {
    set_cflag(ENTROPYFLAG,ix,iy,iz,0);
    set_cflag(ENTROPYFLAG2,ix,iy,iz,0);
  }
  */

  //TODO -- put in check for corrected_polaraxis
  //u to p inversion is done here
  //if(is_cell_corrected_polaraxis(ix,iy,iz))
  //{
  //  u2p_solver_Bonly(uu,pp,&geom); // invert only the magnetic field, the rest will be overwritten
  //}
  //else
  //{
    u2p_device(uu,pp,&geom,corrected,fixups); // regular inversion
    //}

    //TODO
  /*
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

  /*
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
  */
    
  //set new primitives and conserved
  for(int iv=0;iv<NV;iv++)
  { 
    set_u(p_arr,iv,ix,iy,iz,pp[iv]);
  }

  //TODO
  /*
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
    */

} 




//**********************************************************************
//high-level u2p solver
//**********************************************************************

__device__ __host__ int u2p_device(ldouble *uu0, ldouble *pp, void *ggg, int corrected[3], int fixups[2])
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
  ldouble u0=pp[1];
  
  //************************************
  //hot hydro - conserving energy
  
  
  //negative uu[0] = rho u^t
  if(uu[0] * gdetu_inv < 0.)
  {
    //TODO
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
               //ANDREW -- but ret=-1 if energy inversion failes but entropy inversion does not!
               //ANDREW -- do we always want a fixup if we have negative uu[0] ? 
    u2pret=-1; // indicates that inversion is needed
    
#ifndef SWAPPAPC
    //TODO -- global array
    //global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]++;  //but count as fixup
#endif
  }

  if(u2pret!=0)  // u2pret=-1 at this stage, so this is always satisfied
  {
#ifdef ENFORCEENTROPY  
    u2pret=-1;  //skip hot energy-conserving inversion and go to entropy inversion
#else
    //TODO
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
	//TODO -- need to pass in mpi variables TOI,TOJ
        //printf("u2p_entr     >>> %d %d <<< %d >>> %e > %e\n",geom->ix + TOI, geom->iy + TOJ,u2pret,u0,pp[1]);
      }
      
      //TODO
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
      //TODO need to pass TOI, TOJ, TOK, PROCID
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
  
  return ret;
} 


//**********************************************************************
// solver wrapper
//**********************************************************************

__device__ __host__ int u2p_solver_device(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
{

  //TODO
  /*
#ifdef NONRELMHD
  return u2p_solver_nonrel_device(uu,pp,ggg,Etype,verbose);
#endif
  */
  
  int (*solver)(ldouble*,ldouble*,void*,int,int);
  struct geometry *geom
    = (struct geometry *) ggg;

    
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

  if(geom->ix==ixTEST && geom->iy==iyTEST && geom->iz==izTEST)
    printf("In u2p_solver_W_device!\n");
  
  ldouble rho,uint,w,W,alpha,D,Sc;
  ldouble ucon[4],ucov[4],utcon[4],utcov[4],ncov[4],ncon[4];
  ldouble Qcon[4],Qcov[4],Qconp[4],Qcovp[4],Qtcon[4],Qtcov[4],Bcon[4],Bcov[4];
  ldouble jmunu[4][4];
  ldouble Qtsq,Qn,QdotB,QdotBsq,Bsq;
  
  
  ldouble pgamma=GAMMA;
  //TODO
  //#ifdef CONSISTENTGAMMA
  //pgamma=pick_gammagas(geom->ix,geom->iy,geom->iz);
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
    f_u2p=&f_u2p_hot_device;
  if(Etype==U2P_ENTROPY)
    f_u2p=&f_u2p_entropy_device;
  
  
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
  int i_increase = 0;
  ldouble f0,f1,dfdW,err;
  ldouble EPS=1.e-4;
  ldouble Wprev=W;
  ldouble cons[7]={Qn,Qtsq,D,QdotBsq,Bsq,Sc,Qdotnp};

  /*
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
    if(err<U2PCONV)
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
    
    if(fabs((W-Wprev)/Wprev)<U2PCONV && err<1.e-1) break;
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
  */
  
  return 0; //ok
}


//********************************************
//Harm u2p_hot
//********************************************

__device__ __host__ static FTYPE dpdWp_calc_vsq(FTYPE Wp, FTYPE D, FTYPE vsq, FTYPE gamma)
{
  FTYPE W=Wp+D;
  return( (gamma - 1.) * (1. - vsq) /  gamma ) ;
}

// 1 / (d(u+p)/dp)
__device__ __host__ static FTYPE compute_idwmrho0dp(FTYPE wmrho0, FTYPE gamma)
{
  return((gamma-1.)/gamma);
}


// 1 / (drho0/dp) holding wmrho0 fixed
__device__ __host__ static FTYPE compute_idrho0dp(FTYPE wmrho0)
{
  return(0.0);
}

__device__ __host__ static int f_u2p_hot(ldouble Wp, ldouble* cons,ldouble *f,ldouble *df,
					 ldouble *err,ldouble pgamma)
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
__device__ __host__ static FTYPE pressure_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma)
{
  ldouble igammar = (gamma-1.)/gamma;
  return(igammar*wmrho0) ;
}

// local aux function
__device__ __host__ static FTYPE compute_inside_entropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma)
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
__device__ __host__ static FTYPE compute_specificentropy_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma)
{
  FTYPE insideentropy,specificentropy;

  insideentropy=compute_inside_entropy_wmrho0_idealgas(rho0, wmrho0,gamma);
  specificentropy=log(insideentropy);

  return(specificentropy);

}

// used for utoprim_jon when doing entropy evolution
// dSspecific/d\chi
__device__ __host__ static FTYPE compute_dspecificSdwmrho0_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0,FTYPE gamma)
{
  FTYPE dSdchi;

  dSdchi = 1.0/((gamma-1.)*wmrho0);
  // Again, GAMMA->1 means dSdchi->\infty unless \chi->0 or rho0->0

  return(dSdchi);

}

// dSspecific/drho0
__device__ __host__ static FTYPE compute_dspecificSdrho_wmrho0_idealgas(FTYPE rho0, FTYPE wmrho0, FTYPE gamma)
{
  FTYPE dSdrho;
  
  dSdrho=gamma/((1.0-gamma)*rho0);

  return(dSdrho);
}


__device__ __host__ static int f_u2p_entropy(ldouble Wp, ldouble* cons, ldouble *f, ldouble *df,
					     ldouble *err,ldouble pgamma)
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
//Call the kernel
//**********************************************************************

int calc_u2p_gpu(int setflags)
{

  ldouble *d_u_arr, *d_p_arr;

  cudaError_t err = cudaSuccess;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // NOTE: size of xb is copied from initial malloc in misc.c
  // these need to be long long if the grid is on one tile and large (~256^3)
  long long Nprim  = (SX)*(SY)*(SZ)*NV;

  // allocate and sync prims and cons to device 
  err = cudaMalloc(&d_p_arr,    sizeof(ldouble)*Nprim);  
  err = cudaMalloc(&d_u_arr,    sizeof(ldouble)*Nprim);

  err = cudaMemcpy(d_p_arr, p, sizeof(ldouble)*Nprim, cudaMemcpyHostToDevice);  // is this used to seed?
  err = cudaMemcpy(d_u_arr, u, sizeof(ldouble)*Nprim, cudaMemcpyHostToDevice);
  

  // launch calc_primitives_kernel

  int threadblocks = (Nloop_0 / TB_SIZE) + ((Nloop_0 % TB_SIZE)? 1:0);

  cudaEventRecord(start);
  calc_primitives_kernel<<<threadblocks, TB_SIZE>>>(Nloop_0, setflags, 
                                                    d_loop0_ix, d_loop0_iy, d_loop0_iz,
                                                    d_u_arr, d_p_arr,
                                                    d_x, d_gcov, d_gcon);

  cudaEventRecord(stop);
  err = cudaPeekAtLastError();
  cudaDeviceSynchronize(); //TODO: do we need this, does cudaMemcpy synchronize?

  cudaEventSynchronize(stop);
  float tms = 0.;
  cudaEventElapsedTime(&tms, start,stop);
  printf("gpu u2p time: %0.2f \n",tms);

  // ======= TODO
  // Free Device Memory
  cudaFree(d_u_arr);
  cudaFree(d_p_arr);

  return 0;
}


