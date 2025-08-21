/*! \file u2p.c
 \brief MHD Conserved to primitives conversion
 */

#include "ko.h"
static ldouble u2p_solver_ff_parallel(ldouble Wguess,ldouble* cons, int Etype,int verbose);


static int f_u2p_parallel_hot(ldouble W, ldouble* cons,ldouble *f,ldouble *df,ldouble *err);
static int f_u2p_parallel_entropy(ldouble W, ldouble* cons,ldouble *f,ldouble *df,ldouble *err);


//**********************************************************************
//force-free solver 
//**********************************************************************
int
u2p_solver_ff(ldouble *uu, ldouble *pp, void *ggg, int verbose)
{

#ifdef FORCEFREE
  int i,j,k;
  ldouble alpha,alphasq;
  ldouble rho,entr,uint,Bsq,Bmag;
  ldouble gamma2_perp,gamma_perp,gamma2_par,gamma_par,gamma2,gamma;
  ldouble betavec[4];
  ldouble Econ[4],Ecov[4],Bcon[4],Bcov[4];
  ldouble Qcov[4],Qcon[4],Qtildecov[4],Qtildecon[4];
  ldouble vcon_perp[4],vcon_par[4],ucon_tot[4];
  int Etype;
  int rhofloored=0; //did we apply absolute floor to density or not? 
  
  // convert VELR to VELPRIM if necessary
  if(VELPRIM != VELR)
  {
    printf("u2p_solver_ff requires VELPRIM==VELR!!\n");
    exit(-1);
  }
  
  /****************************/
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

  ldouble pgamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  pgamma=pick_gammagas(geom->ix,geom->iy,geom->iz);
  #endif

  //alpha
  alpha   = geom->alpha;
  alphasq = alpha*alpha;

  //beta vector components
  betavec[0] = 0.;
  betavec[1] = alphasq * GG[0][1];
  betavec[2] = alphasq * GG[0][2];
  betavec[3] = alphasq * GG[0][3];

  //////////////////////////////////////////////////////////////////////////////////////////
  // Force-Free EM part
  /////////////////////////////////////////////////////////////////////////////////////////
  
  //McKinney defn of B^mu = B^mu HARM * alpha 
  Bcon[0]=0.;
  Bcon[1]=uu[B1] * alpha * gdetu_inv;
  Bcon[2]=uu[B2] * alpha * gdetu_inv;
  Bcon[3]=uu[B3] * alpha * gdetu_inv;

  indices_21(Bcon,Bcov,gg);
  Bsq = dot(Bcon,Bcov);
  Bmag = sqrt(Bsq);

  //Q_mu= alpha*T^t_mu 
  Qcov[0] = 0.; 
  Qcov[1] = uu[VXFF] * alpha * gdetu_inv;
  Qcov[2] = uu[VYFF] * alpha * gdetu_inv;
  Qcov[3] = uu[VZFF] * alpha * gdetu_inv;
  
  //Qtilde_mu = Q_mu + n_mu (n.Q)
  Qtildecov[0] = Qcov[1]*betavec[1] + Qcov[2]*betavec[2] + Qcov[3]*betavec[3];
  Qtildecov[1] = Qcov[1];
  Qtildecov[2] = Qcov[2];
  Qtildecov[3] = Qcov[3];
  indices_12(Qtildecov,Qtildecon,GG);  

  // get three velocity
  vcon_perp[0] = 0.;
  vcon_perp[1] = Qtildecon[1]/Bsq;
  vcon_perp[2] = Qtildecon[2]/Bsq;
  vcon_perp[3] = Qtildecon[3]/Bsq;

  // get Lorentz factor
  int ll,mm; 
  ldouble vsq_perp=0.;
  for(ll=1;ll<4;ll++)
  {
    for(mm=1;mm<4;mm++)
    {
      vsq_perp += gg[ll][mm]*vcon_perp[ll]*vcon_perp[mm];
    }
  }
  
  gamma2_perp = 1./(1-vsq_perp);
  gamma_perp = sqrt(gamma2_perp);

  // ANDREW moved vsq_perp limiter here 
  ldouble vsqmax = 1. - 1./(GAMMAMAXFF*GAMMAMAXFF);  
  int hitgammaceil=0;

  if(!isfinite(gamma_perp) || vsq_perp>vsqmax)
  {
    //verbose=1;
    if(verbose>0) printf("gamma_perp > gamma_max in u2p_solver_ff %d %d %d | %e %e %e %e| %e\n",
			 geom->ix,geom->iy,geom->iz,
			 gg[1][1]*vcon_perp[1]*vcon_perp[1],
			 gg[2][2]*vcon_perp[2]*vcon_perp[2],
			 gg[3][3]*vcon_perp[3]*vcon_perp[3],
			 2*gg[1][3]*vcon_perp[1]*vcon_perp[3],
			 vsq_perp);

    ldouble vfac = sqrt(vsqmax/vsq_perp);
    vcon_perp[1] *= vfac;
    vcon_perp[2] *= vfac;
    vcon_perp[3] *= vfac;
    vsq_perp = vsqmax;
    gamma2_perp = 1./(1-vsq_perp);
    gamma_perp = sqrt(gamma2_perp);
    hitgammaceil=1;
    
    //getchar();
    //return -200;
  }


  //////////////////////////////////////////////////////////////////////////////////////////
  // Hydro part -- solve for parallel velocity
  /////////////////////////////////////////////////////////////////////////////////////////

  ldouble gammam2_perp=1./gamma2_perp;
  ldouble afac=(pgamma-1.)/pgamma;

  // MHD conserveds
  ldouble D  = uu[RHO] * alpha * gdetu_inv; // uu[RHO] = gdet rho ut, so D = gamma * rho
  ldouble Sc = uu[ENTR] * alpha * gdetu_inv; // uu[ENTR] = gdet S ut, so Sc = gamma * S

  int ret=-1;
  int whicheqs_parallel = 1;
  ldouble vpar = 0, vsq_par=0;
  ldouble W = -1.;
  
  // which equations to start with for parallel solver? 
  #if !defined(FORCEFREE_SOLVE_PARALLEL)  // no parallel solver at all
  whicheqs_parallel = 5; 
  #elif defined(FORCEFREE_PARALLEL_COLD)  // cold momentum conservation
  whicheqs_parallel = 4;
  #elif defined(FORCEFREE_PARALLEL_ENTROPY) // adiabatic specific enthalpy
  whicheqs_parallel = 3;
  #elif defined(FORCEFREE_PARALLEL_MHD) // hot inversion from MHD
  whicheqs_parallel = 1;
  #else
  #error Incorrect option chosen for FORCEFREE_SOLVE_PARALLEL!
  #endif

  #ifdef FORCEFREE_NO_PARALLEL_ATBH // no parallel velocity under BH
  ldouble xx[4],xxBL[4];

  #if (defined(PRECOMPUTE_MY2OUT) && (OUTCOORDS==BLCOORDS || OUTCOORDS==KSCOORDS))
  get_xxout(geom->ix, geom->iy, geom->iz, xxBL);
  #else
  get_xx(geom->ix,geom->iy,geom->iz,xx);      
  coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
  #endif
  
  if(xxBL[1] < rhorizonBL)
    whicheqs_parallel=5;
      
  #endif // FORCEFREE_NO_PARALLEL_ATBH

  //***************************************
  // parallel solver
  while(whicheqs_parallel<6)
  {
    
    // no parallel velocity
    if(whicheqs_parallel == 5)
    {
      vpar = 0.; // no parallel velocity
      Etype = U2P_ENTROPY;
    }
    
    // cold parallel equation, conservation of b^0
    else if(whicheqs_parallel == 4) 
    {
    
       Etype = U2P_ENTROPY;

       ldouble w_s = 1.;

       // specific enthalpy is included in the conserved quantity
       // if we don't define  FORCEFREE_PARALLEL_COLD
       // so these equations are inconsistent
       // as a guess, divide out the CURRENT specific enthalpy
       #ifndef FORCEFREE_PARALLEL_COLD
       printf("whicheqs_parallel=4 but FORCEFREE_PARALLEL_COLD not defined!\n");
       exit(-1);
       //w_s = 1 + pgamma*pp[UU]/pp[RHO]; // specific enthalpy
       #endif
       
       ldouble upar = uu[UUFF] * alpha * gdetu_inv / Bmag / w_s;
  
       gamma2 = gamma2_perp*(1 + upar*upar);
       vpar = my_sign(upar) * sqrt(gammam2_perp*upar*upar / (1 + upar*upar));
       vsq_par = vpar*vpar;
    }

    // Numerical inversion using hot or entropy
    else //(whicheqs_parallel == 3,2,1)
    {
  
      // initial guess for W is based on current primitives
      // recall that force free requires VELPRIM = VELR

      ldouble Y=0;
      ldouble Z=0;
      
      ldouble utcon_prev[4];
      utcon_prev[0]=0.;
      utcon_prev[1]=pp[VX];
      utcon_prev[2]=pp[VY];
      utcon_prev[3]=pp[VZ];
  
      ldouble qsq=0.;
      for(i=1;i<4;i++)
        for(j=1;j<4;j++)
          qsq+=utcon_prev[i]*utcon_prev[j]*gg[i][j];

      ldouble gamma2_prev=1.+qsq;
      ldouble Wguess = (pp[RHO] + pgamma*pp[UU])*gamma2_prev;

      // inversion based on adiabatically evolved specific enthalpy
      if(whicheqs_parallel == 3) 
      {
  
        ldouble M = uu[UUFF] * alpha * gdetu_inv / Bmag;
        Y = M * D;
        Z = 0.;
        Etype = U2P_ENTROPY; 
      }
      // inversion constants based on full MHD conserveds
      else if(whicheqs_parallel == 2 || whicheqs_parallel == 1) 
      {
      
        Qcov[0] = uu[UU] * alpha * gdetu_inv - uu[RHO] * alpha * gdetu_inv;
        Qcov[1] = uu[VX] * alpha * gdetu_inv;
        Qcov[2] = uu[VY] * alpha * gdetu_inv;
        Qcov[3] = uu[VZ] * alpha * gdetu_inv;
  
        indices_12(Qcov,Qcon,GG);  

        ldouble QdotB = dot(Bcon,Qcov);
        ldouble QdotEta = -alpha * Qcon[0];

        Y = QdotB / Bmag;
        Z = QdotEta + 0.5*Bsq*(1. + vsq_perp);

        if(whicheqs_parallel == 1)
	  Etype=U2P_HOT;
        else if(whicheqs_parallel == 2)
	  Etype=U2P_ENTROPY;
      }
    
      // solver constants
      ldouble cons[6];
      cons[0] = D;
      cons[1] = Y;
      cons[2] = Z;
      cons[3] = gammam2_perp;
      cons[4] = afac;
      cons[5] = Sc;
  
      // solve for W (whicheqs==3,2,1)
      W = u2p_solver_ff_parallel(Wguess,cons,Etype,verbose);
      vpar = Y/W;
    } //(whicheqs_parallel == 3,2,1)

    
    // get parallel 3-velocity from vpar solution
    vcon_par[0] = 0.;
    vcon_par[1] = vpar*(Bcon[1]/Bmag);
    vcon_par[2] = vpar*(Bcon[2]/Bmag);
    vcon_par[3] = vpar*(Bcon[3]/Bmag);

    vsq_par = vpar*vpar;

    // parallel lorentz factor
    gamma2_par = 1./(1. - vsq_par);
    gamma_par = sqrt(gamma2_par);

    // total lorentz factor
    gamma2 = 1./(1. - vsq_perp - vsq_par);
    gamma = sqrt(gamma2);

    // determine the density
    rho = D/gamma;

    //floor on density to prevent going negative
    //ANDREW - removed RHOFLOOR, now only in check_floors_mhd
    //if(rho<RHOFLOOR || !isfinite(rho))
    if(rho<0 || !isfinite(rho))       
    {
      if(verbose>0)
      {
        //printf("rho (%e) <RHOFLOOR (%e) in u2p_solver_ff \n",rho,RHOFLOOR); //getchar();
        printf("rho (%e) < 0 in u2p_solver_ff \n",rho); //getchar();	
      }
      rho=.99*RHOFLOOR; // slightly below RHOFLOOR so this will be caught in floors
      rhofloored=1;      
    }
    
    // determine the entropy and energy density
    entr = Sc/gamma;
	
    if(Etype==U2P_HOT)
    {
      uint = (W/gamma2 - rho)/pgamma;
    }
    else //(Etype==U2P_ENTROPY)
    {
      uint = calc_ufromS(entr,rho,geom->ix,geom->iy,geom->iz);
    }

    // is the solution for uint acceptable?
    // if not, continue with next set of equations
    if(!isfinite(uint) || uint < 0) 
    {
      //verbose=1;
      int whicheqs_parallel_old = whicheqs_parallel;
      if(whicheqs_parallel<5)
      {
	// skip cold parallel evolution if adiabatic evolution fails
	// we are not conserving the right quantity for cold evolution
	if(whicheqs_parallel==3)
	{
	  whicheqs_parallel=5; // go straight to no parallel velocity evolution
        }
	else
	{
	  whicheqs_parallel += 1;
        }
	if(verbose>0)
	{
	  printf("neg uint in parallel solver: whicheqs %d -> %d \n",
			   whicheqs_parallel_old, whicheqs_parallel);
        }
	continue;
      }
      else //fail!
      {
	if(verbose>0)
	{
	  printf("neg uint in parallel solver: whicheqs 5 -> fail \n");
        }
	break;
      }
    }
    else // acceptable solution for uint
    {
      break;
    }
  } //while(whicheqs_parallel<6)

  //***************************************
  // Apply final floors/ceilings
  
  // apply separate limiters on perp and parallel lorentz factors  
  // ANDREW removed limiter here, only above on vsq_perp
  // and on total lorentz factor in check_floors_mhd() 
  /*
  if(vsq_perp>vsqmax)
  {
    hitgammaceil=1;
    ldouble vfac = sqrt(vsqmax/vsq_perp);
    vcon_perp[1] *= vfac;
    vcon_perp[2] *= vfac;
    vcon_perp[3] *= vfac;
    vsq_perp = vsqmax;
    gamma2_perp = 1./(1-vsq_perp);
    gamma_perp = sqrt(gamma2_perp);
  }
    
  if(vsq_par>vsqmax)
  {
    hitgammaceil=1;
    ldouble vfac = sqrt(vsqmax/vsq_par);
    vcon_par[1] *= vfac;
    vcon_par[2] *= vfac;
    vcon_par[3] *= vfac;
    vsq_par = vsqmax;
    gamma2_par = 1./(1. - vsq_par);
    gamma_par = sqrt(gamma2_par);
  }
  
  // total lorentz factor (recomputed in case we hit ceiling)
  if(hitgammaceil)
  {
    gamma2 = 1./(1. - vsq_perp - vsq_par);
    gamma = sqrt(gamma2);
  }
  */
  
  // apply floor on density
  // ANDREW removed floor, only in check_floors_mhd() now
  /*
  if(rho<RHOFLOOR || !isfinite(rho)) 
  {
    verbose=1;
    if(verbose>0)
    {
      printf("final rho (%e) <RHOFLOOR (%e) in u2p_solver_ff \n",rho,RHOFLOOR); //getchar();
    }
    rho=RHOFLOOR;
    rhofloored=1;
  }  
  */
  
  // apply floor/ceiling on internal energy
  // ANDREW removed floor, only in check_floors_mhd() now
  /*
  if(uint>UURHORATIOMAX*rho || (!isfinite(uint) && isfinite(entr) && isfinite(rho)))
  {
    uint = UURHORATIOMAX*rho;
  }
  if(uint<UURHORATIOMIN*rho)
  {
    uint = UURHORATIOMIN*rho;
  }
  */
  
  //***************************************
  // Final checks 
  // is the final solution for rho acceptable? 
  if(rho<0. || !isfinite(rho)) 
  {
    //verbose=1;
    if(verbose>0)
    {
      printf("neg rho in u2p_solver_ff %e %e %e \n",rho,D,gamma); //getchar();
      printf("gamma_perp %e gamma_par %e gamma %e\n",gamma_perp,gamma_par,gamma);
    }
    ret = -206;
  }
  // is the final solution for uint acceptable? 
  else if(uint<0. || !isfinite(uint)) 
  {
    //verbose=1;
    if(verbose>0)
    {
      printf("neg uint in u2p_solver_ff %d | %e %e %e (%d) %e | ",whicheqs_parallel, uint,entr,rho,rhofloored,D);
      printf("gamma_perp %e gamma_par %e gamma %e\n",gamma_perp,gamma_par,gamma);
      getchar();
    }
    ret = -207;
  }
  else
  {
    ret=0;
  }
  
  // total velocity (VELR)
  ucon_tot[0] = 0.;
  ucon_tot[1]=(vcon_perp[1] + vcon_par[1])*gamma;
  ucon_tot[2]=(vcon_perp[2] + vcon_par[2])*gamma;
  ucon_tot[3]=(vcon_perp[3] + vcon_par[3])*gamma;
  
  //***************************************
  //final new primitives

  pp[B1]=uu[B1] * gdetu_inv; 
  pp[B2]=uu[B2] * gdetu_inv;
  pp[B3]=uu[B3] * gdetu_inv;

  pp[VX] = ucon_tot[1];
  pp[VY] = ucon_tot[2];
  pp[VZ] = ucon_tot[3];
  pp[RHO] = rho;
  pp[UU] = uint;
  pp[ENTR] = entr;

  pp[UUFF] = vpar*gamma_par; 
  pp[VXFF] = vcon_perp[1]*gamma_perp;
  pp[VYFF] = vcon_perp[2]*gamma_perp;
  pp[VZFF] = vcon_perp[3]*gamma_perp;
  
  return ret; 
#endif
}

//*********************************************************************
//parallel momentum part of the force-free solver
//*********************************************************************
// TODO solve for W? Wp? vparallel^2?
// TODO exact solution?
// TODO entropy version? 
// TODO what to do about variable adiabatic index? 
static ldouble
u2p_solver_ff_parallel(ldouble Wguess,ldouble* cons, int Etype, int verbose)
{
  
  ldouble CONV=U2PCONV;
  
  ldouble W=Wguess;
  ldouble f0,f1,dfdW,err,vpar2guess;

  int (*f_u2p_parallel)(ldouble,ldouble*,ldouble*,ldouble*,ldouble*);
  if(Etype==U2P_HOT) 
   f_u2p_parallel=&f_u2p_parallel_hot;
  if(Etype==U2P_ENTROPY) 
   f_u2p_parallel=&f_u2p_parallel_entropy;
  
  // constants
  ldouble D = cons[0];
  ldouble Y = cons[1];
  ldouble Z = cons[2];
  ldouble gammam2_perp = cons[3];
  ldouble afac = cons[4];
  ldouble Sc = cons[5];

  ldouble Ysq=Y*Y;
  ldouble Wsq=W*W;
  if(verbose>1)
    {
      printf("In parallel entropy solver\n");
      printf("%e %e %e %e %e %e\n",D,Y,Z,gammam2_perp,afac,Sc);
    }
  
  int i_increase = 0;
  ldouble increase_fac=2.;
  int dampmax=500;
  
  // Make sure that W is large enough so that v^2 < 1 :
  // Test if initial guess is out of bounds and damp if so
  do
  {
    f0=dfdW=0.;
    
    (*f_u2p_parallel)(W,cons,&f0,&dfdW,&err);

    //vpar2guess = (Y/W)*(Y/W);
    vpar2guess = Ysq/(W*W);
    if(((gammam2_perp - vpar2guess) < 0
       || !isfinite(f0) || !isfinite(dfdW))
       && (i_increase < dampmax))
    {
      if(verbose>1)
	{
	  printf("error in init W : %e -> %e (%e %e %e)\n",W,increase_fac*W,
		 gammam2_perp-vpar2guess,f0,dfdW);
	  printf("D %e Y %e Z %e gammam2 %e a %e Sc %e\n",D,Y,Z,gammam2_perp,afac,Sc);
          getch();
	}
      W *= increase_fac;
      i_increase++;
      continue;
    }
    else
      break;
  }
  while(1);
  
  if(i_increase>=dampmax)
  {
    if(verbose>1)
      {
        printf("failed to find initial W in parallel solver\n");
        printf("W : %e->%e | (%e %e %e)\n",Wguess,W,gammam2_perp-vpar2guess,f0,dfdW);
      }
    return -250;
  }

  (*f_u2p_parallel)(W,cons,&f0,&dfdW,&err);  
  if(verbose>1) printf("\ninitial W:%e\n",W);
  if(verbose>1) printf("initial f:%e\n",f0);
  if(verbose>1) printf("initial err:%e\n\n",err);

  //1d Newton solver
  int iter=0, fu2pret;
  do
  {
    ldouble Wprev=W;
    iter++;

    fu2pret=(*f_u2p_parallel)(W,cons,&f0,&dfdW,&err);
    
    if(verbose>1) printf("%d parallel solver: %e %e %e %e\n",iter,W,f0,dfdW,err);
    
    //convergence test (moved below)
    //if(err<CONV) break;

    // get out of fixed point
    if(dfdW==0.)
    {
      W*=1.1;
      continue;
    }
    
    ldouble Wnew = W-f0/dfdW;

    int idump=0;
    ldouble dumpfac=1.;
    
    //test if goes out of bounds and damp solution if so
    do
    {
      ldouble f0tmp,dfdWtmp,errtmp;
      f0tmp=dfdWtmp=0.;
      (*f_u2p_parallel)(Wnew,cons,&f0tmp,&dfdWtmp,&errtmp);
    
      vpar2guess = Ysq/(Wnew*Wnew);
      if(((gammam2_perp - vpar2guess) < 0
         || !isfinite(f0tmp) || !isfinite(dfdWtmp) || !isfinite(dfdWtmp))
         && (idump < 100))
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
    
    if(idump>=dampmax)
    {
      if(verbose>0) printf("damped unsuccessfuly in parallel solver\n");
      //getchar();
      return -201;
    }
    
    W=Wnew;
    
    if(fabs(W)>BIG)
    {
      if(verbose>0) printf("W has gone out of bounds in parallel solver\n");
      //getchar();
      return -202;
    }
    
    if(err<CONV || (fabs((W-Wprev)/Wprev)<CONV && err<sqrt(CONV))) break;
    
  }
  while(iter<100);
  
  if(iter>=100)
  {
    if(verbose>0)
    {
      printf("iter exceeded in parallel u2p_solver \n");
      printf("W %e f0 %e dfdW %e\n",W,f0,dfdW);
      printf("D %e Y %e Z %e gm2 %e a %e Sc %e\n",D,Y,Z,gammam2_perp,afac,Sc);
    }
    //getchar();
    return -203;
  }
  
  if(!isfinite(W))
  {
    if(verbose>0) printf("nan/inf W in parallel u2p_solver: \n");
    //getchar();
    return -204;
  }

  vpar2guess = Ysq/(W*W);
  if((gammam2_perp - vpar2guess) < 0)
  {
    if(verbose>0) printf("final W out of bounds in parallel solver\n");
    //getchar();
    return -205;
  }

  if(verbose>1)
  {
    printf("final W %e\n",W);
    (*f_u2p_parallel)(W,cons,&f0,&dfdW,&err);
    printf("final f %e\n",f0);
    printf("final err %e\n\n",err);
  }
  return W;
}

static int
f_u2p_parallel_hot(ldouble W, ldouble* cons,ldouble *f,ldouble *df,ldouble *err)
{

  // constants
  ldouble D = cons[0];
  ldouble Y = cons[1];
  ldouble Z = cons[2];
  ldouble gammam2_perp = cons[3];
  ldouble afac = cons[4];
  ldouble Sc = cons[5];
  
  ldouble Ysq = Y*Y;
  ldouble Wsq = W*W;

  ldouble vparsq = Ysq/Wsq;
  ldouble gammam2 = gammam2_perp - vparsq;
  ldouble gammam1 = sqrt(gammam2);

  ldouble rho = D*gammam1;
  ldouble pgas = afac*W*gammam2 - afac*rho;

  ldouble drhodW = D*vparsq/(W*gammam1);
  ldouble dpgasdW = afac*(gammam2_perp + vparsq - drhodW);
  
  // root function

  *f = Z + W - pgas;
  
  // error
  *err = fabs(*f) / (fabs(pgas) + fabs(W) + fabs(Z));

  // derivitive
  *df = 1. - dpgasdW;
  
  return 0;  
}

static int
f_u2p_parallel_entropy(ldouble W, ldouble* cons, ldouble *f, ldouble *df, ldouble *err)
{
  // constants
  ldouble D = cons[0];
  ldouble Y = cons[1];
  //ldouble Z = cons[2];
  ldouble gammam2_perp = cons[3];
  ldouble afac = cons[4];
  ldouble Sc = cons[5];
  
  ldouble Ysq = Y*Y;
  ldouble Wsq = W*W;

  ldouble vparsq = Ysq/Wsq;
  ldouble gammam2 = gammam2_perp - vparsq; //1/gamma^2
  ldouble gammam1 = sqrt(gammam2);
  ldouble rho = D*gammam1;
  ldouble pgas = afac*(W*gammam2 - rho);

  // ANDREW TODO option for variable gamma
  ldouble pgamma = 1./(1.-afac); //=GAMMA; (should be consistent);
  
  ldouble drhodW = D*vparsq/(W*gammam1);
  ldouble dpgasdW = afac*(gammam2_perp + vparsq - drhodW);

  ldouble scalc,dscalcdrho,dscalcdpgas;

  #ifdef NOLOGINS

  // specific entropy (without leading factor of rho)
  scalc = pgas / pow(rho,pgamma) / (pgamma-1.);
  dscalcdrho = -pgamma*pgas/pow(rho,pgamma+1.) / (pgamma-1.);
  dscalcdpgas = 1. / pow(rho,pgamma) / (pgamma-1.);

  #else
  ldouble indexn = 1./(pgamma-1.);
  
  scalc = log(pow(pgas,indexn)/pow(rho,indexn+1.));

  dscalcdrho = -1.*(indexn+1.) / rho;
  dscalcdpgas = indexn/pgas;
  
  #endif

  // root function
  *f = scalc*D - Sc; 
  
  // error
  *err = fabs(*f) / (fabs(scalc*D) + fabs(Sc));

  *df = D*(dscalcdrho*drhodW + dscalcdpgas*dpgasdW); 
  
  return 0;  

}

//**********************************************************************
//ensure force free variables are consistent
//project force free velocity perpindicular to B 
//**********************************************************************
int
fill_ffprims()
{
#ifdef FORCEFREE
  int ii;

  //  #pragma omp parallel for private(ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      int ix,iy,iz,iv;
      
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);

      ldouble pp[NV],uu[NV];
      for(iv=0;iv<NV;iv++)
      {
        pp[iv]=get_u(p,iv,ix,iy,iz);
        uu[iv]=get_u(u,iv,ix,iy,iz);	
      }
      
      // update primitives in a single cell
      fill_ffprims_cell(pp,&geom);

      // determine inversion mixing factor f(sigma)
      
      int fflag = 1;
      ldouble ffval= 1.;
	
#ifdef HYBRID_FORCEFREE
      
#ifdef HYBRID_FORCEFREE_XCUT //cut based on x-domain for test problems
      //if(geom.xx < HYBRID_FORCEFREE_XCUT) ffval = 0.;
      //else ffval = 1.;

      ffval = calc_ffinv_val_x(geom.xx);

#else // cut based on local sigma
      ldouble ucon[4],ucov[4],bcon[4],bcov[4];
      ldouble bsq, sigma;
      calc_ucon_ucov_from_prims(pp, &geom, ucon, ucov);
      calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, &geom, bcon, bcov, &bsq);
      sigma = bsq/pp[RHO];

      ffval = calc_ffinv_val(sigma);
      #endif
      #endif //HYBRID_FORCEFREE
      
      for(iv=0;iv<NV;iv++)
      {
        set_u(p,iv,ix,iy,iz,pp[iv]);
 	//set_u(u,iv,ix,iy,iz,uu[iv]); //ANDREW -- do not need to modify conserverds
      }

      // determine local inversion flags
      int ffflag, mhdflag;
      if(ffval<=0.) ffflag=0; // ff inversion is NOT required
      else ffflag = 1; // ff inversion is required

      if(ffval>=1.) mhdflag=0; // mhd inversion is NOT required
      else mhdflag = 1; // mhd inversion is required

      //printf("%d %d %d %e\n",ix,iy,iz,ffval);
      set_u_scalar(ffinvarr, ix, iy, iz, ffval);
      set_cflag(FFINVFLAG, ix,iy,iz,ffflag);
      set_cflag(MHDINVFLAG, ix,iy,iz,mhdflag);
    }

#endif  
  return 0;
}


int
fill_ffprims_cell(ldouble *pp, void *ggg)
{
#ifdef MAGNFIELD
#ifdef FORCEFREE

  //prepare geometry
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5], (*GG)[5], gdet, gdetu, gdetu_inv,alpha;
  gg=geom->gg;  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
  alpha   = geom->alpha;


  //********************************
  // make force-free velocity prims consistent
  //********************************
  ldouble uperp[4];
  ldouble vpar;
  vpar = get_driftvel_velr(pp,uperp,geom);
  ldouble gammapar = 1./sqrt(1-vpar*vpar);
  
  // save perpindicular relative 4-velocity
  pp[VXFF] = uperp[1];
  pp[VYFF] = uperp[2];
  pp[VZFF] = uperp[3];
  pp[UUFF] = gammapar*vpar; 
  
#endif
#endif
  return 0.;
}

ldouble
get_driftvel_velr(ldouble *pp,ldouble *velff,void* ggg)
{
  //prepare geometry
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5], (*GG)[5], gdet, gdetu, gdetu_inv,alpha;
  gg=geom->gg;  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
  alpha   = geom->alpha;

  // get magn field
  // Note: using McKinney06 defn of B with alpha
  ldouble Bcon[4],Bcov[4],Bsq,Bmag;
  
  Bcon[0]=0.;
  Bcon[1]=pp[B1] * alpha;
  Bcon[2]=pp[B2] * alpha;
  Bcon[3]=pp[B3] * alpha;

  indices_21(Bcon,Bcov,gg);
  Bsq = dot(Bcon,Bcov);
  Bmag = sqrt(Bsq);
  if(Bmag<SMALL) Bmag=SMALL;
  
  // get velocity (recall FORCEFREE requires VELR)
  ldouble utcon[4];
  utcon[0]=0.;
  utcon[1]=pp[VX];
  utcon[2]=pp[VY];
  utcon[3]=pp[VZ];

  conv_vels(utcon,utcon,VELPRIM,VELR,gg,GG);
  
  ldouble qsq=0.;
  int i,j;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      qsq+=utcon[i]*utcon[j]*gg[i][j];
  ldouble gamma=sqrt(1.+qsq);
  ldouble gammam2 = 1./(1+qsq); // 1/gamma^2
  
  // take dot product and find perpindicular gamma
  ldouble udotB = dot(utcon,Bcov); // = gamma*|B|*vpar
  ldouble vpar = udotB/(gamma*Bmag);

  //TODO apply floor to ensure |vpar|<=1 here? where?
  if(vpar*vpar>1)
  {
    //printf("vpar^2 > 1 in make_ffprims_consistent!");
    vpar=vpar/sqrt(vpar*vpar);   
  }
  
  ldouble gammaperpm2 = gammam2  + vpar*vpar; 
  ldouble gammaperp = 1./sqrt(gammaperpm2);

  // find perpindicular 4-velocity
  
  velff[0] = 0;
  velff[1] = gammaperp*((utcon[1]/gamma) - vpar*Bcon[1]/Bmag);
  velff[2] = gammaperp*((utcon[2]/gamma) - vpar*Bcon[2]/Bmag);
  velff[3] = gammaperp*((utcon[3]/gamma) - vpar*Bcon[3]/Bmag);

  return vpar;
}


// get the flag for whether or not to treat a cell FACE as force free
// OR operation on neighbor cells
int get_ffinv_flag_face(int ix, int iy, int iz,int ifacedim)
{
  int flag;
  int flag_p=0;
  int flag_m=0;
  //x-face
  if(ifacedim==0)
  {
    flag_m = get_cflag(FFINVFLAG, ix, iy, iz);
    flag_p = get_cflag(FFINVFLAG, ix+1, iy, iz);
  }
  if(ifacedim==1)
  {
    flag_m = get_cflag(FFINVFLAG, ix, iy, iz);
    flag_p = get_cflag(FFINVFLAG, ix, iy+1, iz);
  }
  if(ifacedim==2)
  {
    flag_m = get_cflag(FFINVFLAG, ix, iy, iz);
    flag_p = get_cflag(FFINVFLAG, ix, iy, iz+1);
  }

  if(flag_p==0 && flag_m==0)
    flag = 0;
  else
    flag = 1;

  return flag;
}


// calculate f(sigma), the value of the mixing function for hybrid force-free MHD
ldouble
calc_ffinv_val(ldouble sigma)
{
  ldouble ffval = 0.;
#ifdef FORCEFREE
  ffval = 1.; 

#if defined(HYBRID_FORCEFREE) && !defined(HYBRID_FORCEFREE_XCUT)
  ldouble tanhwidth = HYBRID_FORCEFREE_WIDTH;
  
  if(tanhwidth<=0) // step function transition
  {
     if(sigma >= HYBRID_FORCEFREE_SIGMACUT) ffval = 1.;
     else ffval = 0.;
  }
  else //finite width transition
  {
     if(sigma < ffinv_lower_cutoff) ffval=0.;
     else if(sigma > ffinv_upper_cutoff) ffval=1.;
     else
     {
       ldouble sigcut = HYBRID_FORCEFREE_SIGMACUT;
       ldouble normsigma = sigma/sigcut;
       ldouble sigpow = pow(normsigma, 2./tanhwidth);
       ffval = sigpow / (1. + sigpow);
     }
  }
  
  if(isnan(ffval) || ffval <0. || ffval>1. || isnan(sigma))
  {
    ffval = 0.; // default MHD inversion 
  }
#endif
#endif
  return ffval;
}

// calculate f(x), the value of the mixing function for hybrid force free MHD
// (1D tests)
ldouble
calc_ffinv_val_x(ldouble x)
{
  ldouble ffval = 0.;
#ifdef FORCEFREE
  ffval = 1.; 

#if defined(HYBRID_FORCEFREE) && defined(HYBRID_FORCEFREE_XCUT)
  ldouble tanhwidth = HYBRID_FORCEFREE_WIDTH;
  
  if(tanhwidth<=0) // step function transition
  {
     if(x >= HYBRID_FORCEFREE_XCUT) ffval = 1.;
     else ffval = 0.;
  }
  else //finite width transition
  {
     ldouble xcut = HYBRID_FORCEFREE_XCUT;
     ffval = 0.5 + 0.5*tanh((x-xcut)/tanhwidth);
     if(ffval<(1./64.)) ffval=0.;
     if(ffval>(63./64.)) ffval=1.;
  }
  
  if(isnan(ffval) || ffval <0. || ffval>1. || isnan(x))
  {
    ffval = 0.; // default MHD inversion 
  }
#endif
#endif
  return ffval;
}


//**********************************************************************
//count the number of ff and mhd inversions
//**********************************************************************

int count_ff(int *n, int *n2, int *n3)
{
  int nff=0,nmhd=0,nrhofloor=0;
  int ii,ix,iy,iz;
  int nffloc=0,nmhdloc=0,nrhofloorloc=0;

  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  nffloc+=get_cflag(FFINVFLAG,ix,iy,iz); 
	  nmhdloc+=get_cflag(MHDINVFLAG,ix,iy,iz);
	  nrhofloorloc+=get_cflag(RHOFLOORFLAG,ix,iy,iz);
	}

#ifdef MPI
  MPI_Allreduce(&nffloc, &nff, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  
  MPI_Allreduce(&nmhdloc, &nmhd, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  
  MPI_Allreduce(&nrhofloorloc, &nrhofloor, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);  									  
#else
  nff=nffloc;
  nmhd=nmhdloc;
  nrhofloor=nrhofloorloc;
#endif
    
  *n = nff;
  *n2 = nmhd;
  *n3 = nrhofloor;
  return 0;
}

// Source term for adiabatic evolution of mu b^0
// Based on calc_shear_lab in rad.c
ldouble
calc_uuff_source(ldouble *pp0, void* ggg,int *derdir)
{

  int i,j,k,iv;

  struct geometry *geomin
    = (struct geometry *) ggg;

  int ix,iy,iz;
  ix=geomin->ix;
  iy=geomin->iy;
  iz=geomin->iz;

  //lets get geometry again, in case we aren't in internal coordinates
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  ldouble (*gg)[5],(*GG)[5];
  gg=geom.gg;
  GG=geom.GG;

  // derivatives of entropy/particle
  ldouble ds[3]; //d_i (s)
  for(i=0;i<3;i++)
  {
    ds[i] = 0.;
  }
  
  ldouble ppm1[NV],ppp1[NV],pp[NV];
  ldouble ggm1[4][5],GGm1[4][5];
  ldouble ggp1[4][5],GGp1[4][5];
  ldouble xxvecm1[4],xxvec[4],xxvecp1[4];
  int idim;
  
  //primitives at cell basing on pp[]
  get_xx(ix,iy,iz,xxvec);
  for(iv=0;iv<NV;iv++)
  {
    pp[iv]=pp0[iv];
  }
  
  ldouble sc=pp[ENTR]/pp[RHO];
  
  //derivatives
  for(idim=1;idim<4;idim++)
  {
      if(idim==1)
	{
	  get_xx(ix-1,iy,iz,xxvecm1);
	  get_xx(ix+1,iy,iz,xxvecp1);	  	   
	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix-1,iy,iz);
	      ppp1[iv]=get_u(p,iv,ix+1,iy,iz);
	    }
	  pick_g(ix-1,iy,iz,ggm1);  pick_G(ix-1,iy,iz,GGm1);
	  pick_g(ix+1,iy,iz,ggp1);  pick_G(ix+1,iy,iz,GGp1);
	}

      if(idim==2)
	{
	  get_xx(ix,iy-1,iz,xxvecm1);
	  get_xx(ix,iy+1,iz,xxvecp1);	  
	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix,iy-1,iz);
	      ppp1[iv]=get_u(p,iv,ix,iy+1,iz);
	    }
	  pick_g(ix,iy-1,iz,ggm1);  pick_G(ix,iy-1,iz,GGm1);
	  pick_g(ix,iy+1,iz,ggp1);  pick_G(ix,iy+1,iz,GGp1);	    
	}

     if(idim==3)
       {
	 get_xx(ix,iy,iz-1,xxvecm1);
	 get_xx(ix,iy,iz+1,xxvecp1);
	 for(iv=0;iv<NV;iv++)
	   {
	     ppm1[iv]=get_u(p,iv,ix,iy,iz-1);
	     ppp1[iv]=get_u(p,iv,ix,iy,iz+1);
	   }
	 pick_g(ix,iy,iz-1,ggm1);  pick_G(ix,iy,iz-1,GGm1);
	 pick_g(ix,iy,iz+1,ggp1);  pick_G(ix,iy,iz+1,GGp1);
       }

     //calculate entropy per particle
     ldouble sm1=ppm1[ENTR]/ppm1[RHO];
     ldouble sp1=ppp1[ENTR]/ppp1[RHO];
     
     ldouble dl,dr,dc;
     ldouble dl2,dr2,dc2;

     dc=(sp1-sm1) / (xxvecp1[idim] - xxvecm1[idim]);
     dr=(sp1-sc ) / (xxvecp1[idim] - xxvec[idim]);
     dl=(sc -sm1) / (xxvec[idim] - xxvecm1[idim]);
	 
     //to avoid corners
     if((ix<0 && iy==0 && iz==0 && idim!=1) ||
	(iy<0 && ix==0 && iz==0 && idim!=2) ||
	(iz<0 && ix==0 && iy==0 && idim!=3))
     {
	     ds[idim-1]=dr;
     }
     else if((ix<0 && iy==NY-1 && iz==NZ-1 && idim!=1) ||
	    (iy<0 && ix==NX-1 && iz==NZ-1 && idim!=2) ||
	    (iz<0 && ix==NX-1 && iy==NY-1 && idim!=3))
     {
	     ds[idim-1]=dl;
     }
     else if((ix>=NX && iy==0 && iz==0 && idim!=1) ||
	    (iy>=NY && ix==0 && iz==0 && idim!=2) ||
	    (iz>=NZ && ix==0 && iy==0 && idim!=3))
     {
	     ds[idim-1]=dr;

     }
     else if((ix>=NX && iy==NY-1 && iz==NZ-1 && idim!=1) ||
	    (iy>=NY && ix==NX-1 && iz==NZ-1 && idim!=2) ||
	    (iz>=NZ && ix==NX-1 && iy==NY-1 && idim!=3))
     {
	     ds[idim-1]=dl;
     }
     else
     {
        //choice of 1st order derivative
	if(derdir[idim-1]==0)
	{
	  ds[idim-1]=dc;
	}
	if(derdir[idim-1]==1)
	{
	  ds[idim-1]=dl;
	}
	if(derdir[idim-1]==2)
	{
	  ds[idim-1]=dr;
	}
      }
  } //for(idim=1;idim<4;idim++)

  // calculate source term
  ldouble utcon[4],ucon[4],ucov[4];
  utcon[1]=pp[VX];  utcon[2]=pp[VY];  utcon[3]=pp[VZ];
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);  

  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(ix,iy,iz);
  #endif
  ldouble pgas=(gamma-1.)*pp[UU];
  ldouble rho=pp[RHO];

  ldouble BdotdS = pp[B1]*ds[0] + pp[B2]*ds[1] + pp[B3]*ds[2];
  ldouble source = BdotdS * (pgas/rho) * (1./ucon[0]);

  if(!isfinite(source) || isnan(source))
  {
    printf("mub^0 source not finite in FF! %d %d %d | %.2e %.2e | %.2e %.2e %.2e | %.2e %.2e %.2e",
	   ix,iy,iz,pgas/rho,ucon[0],pp[B1],pp[B2],pp[B3],ds[0],ds[1],ds[2]);
    source=0.;
    getch();
  }
  return source;
}
