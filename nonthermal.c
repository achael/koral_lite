/*! \file nonthermal.c
 \brief Routines relating to nonthermal relativistic electrons
 */


#include "ko.h"


//*********************************************/
//* set arrays for gamma faces and centers
//*********************************************/

//ANDREW changed so RELGAMMAMIN and RELGAMMAMAX are the lower and upper bin edges

int set_relel_gammas()
{
#ifdef RELELECTRONS
  
  int ie;
  
  ldouble relel_log_gammas[NRELBIN];
  ldouble relel_log_gammas_e[NRELBIN+1];

#ifdef RELEL_INTEGRATE_SIMPSON   // Simpson's rule in log space -- no bin edges
  if(NRELBIN % 2 == 0) {printf("NRELBIN MUST BE ODD TO USE SIMPSON"); getch(); return -1;}

  log10binspace = (log10(RELGAMMAMAX) - log10(RELGAMMAMIN))/((double)(NRELBIN-1.));
  logbinspace = log10binspace * log(10.);
  logbinspace_inv = 1./logbinspace;
  binspace = pow(10., log10binspace);
  binspace_inv = 1./binspace;
  if(PROCID==0) printf("nonthermal bin spacing: %e %e \n",logbinspace,binspace); 

  relel_gammas_e[0] = RELGAMMAMIN;
  relel_log_gammas_e[0] = log(RELGAMMAMIN);
  relel_log_gammas[0] = log(RELGAMMAMIN);
  relel_gammas[0] = RELGAMMAMIN;

#else // rectangular bins with separate edges
  log10binspace = (log10(RELGAMMAMAX) - log10(RELGAMMAMIN))/((double)(NRELBIN));
  logbinspace = log10binspace * log(10.);
  logbinspace_inv = 1./logbinspace;
  binspace = pow(10., log10binspace);
  binspace_inv = 1./binspace;
  if(PROCID==0) printf("nonthermal bin spacing: %e %e \n",logbinspace,binspace); 

  relel_gammas_e[0] = RELGAMMAMIN;
  relel_log_gammas_e[0] = log(RELGAMMAMIN);
  
  relel_log_gammas[0] = relel_log_gammas_e[0] + 0.5*logbinspace;
  relel_gammas[0] = exp(relel_log_gammas[0]);
#endif

  relel_gammas_inv[0] = 1./relel_gammas[0];
  relel_gammas_e_inv[0] = 1./relel_gammas_e[0];
  
  // gammas of bin centers and and edges
  for(ie=1; ie<NRELBIN; ie++)
  {
    relel_log_gammas_e[ie] = relel_log_gammas_e[ie-1] + logbinspace;
    relel_log_gammas[ie] = relel_log_gammas[ie-1] + logbinspace;
    relel_gammas_e[ie] = exp(relel_log_gammas_e[ie]);
    relel_gammas[ie] = exp(relel_log_gammas[ie]);
    
    relel_gammas_inv[ie] = 1./relel_gammas[ie];
    relel_gammas_e_inv[ie] = 1./relel_gammas_e[ie];
  }

  relel_log_gammas_e[NRELBIN] = log(RELGAMMAMAX);
  relel_gammas_e[NRELBIN] = RELGAMMAMAX; //equal to relel_gammas_e[NRELBIN-1] for simpson
  relel_gammas_e_inv[NRELBIN] = 1./RELGAMMAMAX;
  
#ifdef RELEL_CONST_INJPARAMS
  my_err("RELEL_CONST_INJPARAMS is deprecated, use RELEL_HEAT_FIX_LIMITS instead!")
#endif

  //precompute injection arrays    
#if defined(RELEL_HEAT_FIX_LIMITS) && defined(RELEL_HEAT_FIX_INDEX)
  ldouble p_index = RELEL_HEAT_INDEX;
  ldouble injmax = RELEL_INJ_MAX;
  ldouble injmin = RELEL_INJ_MIN;

  //normalization
  if(p_index < 0)
    if(PROCID==0)
      printf("WARNING: RELEL_HEAT_INDEX < 0 corresponds to positive power law!"); 

  if(p_index!=2 && p_index!=1)
  {
    xx_relel_inj = (pow(injmax, 2.-p_index) - pow(injmin, 2.-p_index))/(2.-p_index) - (pow(injmax, 1.-p_index) - pow(injmin, 1.-p_index))/(1.-p_index);
  }
  else if (p_index == 2)
  {
    xx_relel_inj = 1./injmax - 1./injmin + log(injmax/injmin);
  }
  else if (p_index == 1)
  {
    xx_relel_inj = injmax - injmin - log(injmax/injmin);
  }

  //power array
  for(ie=0; ie<NRELBIN; ie++) relel_injpow[ie] = pow(relel_gammas[ie], -p_index);
#endif

#if defined(RELEL_HEAT_FIX_FRAC) && defined(RELEL_HEAT_FIX_INDEX)
  ldouble eta= RELEL_HEAT_FRAC/(1.-RELEL_HEAT_FRAC);
  ldouble lhs = 6.*(RELEL_HEAT_INDEX-2.)*eta;
  ldouble a = 0.25 * pow(lhs, 0.25);
  gamma_inj_fixed_ratio = -4.*gsl_sf_lambert_Wm1(-a);

#endif
#endif
  
  return 0;
}


//*******************************************************
//calculates fluid quantities for relativistic electrons
//*******************************************************
//Numerically integrate an array q of size NRELBIN
ldouble integrate_relel_q(ldouble *q)
{
  
  int ie;
  ldouble dgamma;
  ldouble qint=0.;
#ifdef RELELECTRONS
#if defined(RELEL_INTEGRATE_SIMPSON) // Simpson's rule in log space
  dgamma = one_third*logbinspace;
  for (ie=1; ie<=(NRELBIN-1)/2; ie++)
  {  
    qint += relel_gammas[2*ie-2]*q[2*ie-2];
    qint += 4.*relel_gammas[2*ie-1]*q[2*ie-1];
    qint += relel_gammas[2*ie]*q[2*ie];
  }
  qint *= dgamma;

#elif defined(RELEL_INTEGRATE_LOGSPACE) // rectangular integration in log space
  dgamma = logbinspace;
  for (ie=0; ie<NRELBIN; ie++)
  {
    qint += dgamma*relel_gammas[ie]*q[ie];
  }
  
#else // rectangular integration in regular space
  for (ie=0; ie<NRELBIN; ie++) 
  {
    dgamma = relel_gammas_e[ie+1]-relel_gammas_e[ie];
    qint += (dgamma*q[ie]);
  }
#endif
#endif
  return qint;
}


//Total number density of relativistic electrons  
ldouble calc_relel_ne(ldouble *pp)
{
  ldouble ne_relel=0.;
#ifdef  RELELECTRONS
  ldouble ne_tot=pp[RHO]*one_over_mue_mp;
  ldouble q[NRELBIN];
  ldouble dgamma;
  int ie;
  
  for (ie=0; ie<NRELBIN; ie++) 
  {
    q[ie] = pp[NEREL(ie)];
  }
  
  ne_relel = integrate_relel_q(q);
  if (ne_relel < 0.)
    ne_relel=0.;
#endif
  return ne_relel;
}


//Total energy density of relativistic electrons  
ldouble calc_relel_uint(ldouble *pp)
{
  ldouble u_relel=0.;
#ifdef RELELECTRONS
  ldouble g2, g1;
  ldouble dgamma;
  int ie;
  ldouble q[NRELBIN];
  
  for (ie=0; ie<NRELBIN; ie++) 
  {
    q[ie] = pp[NEREL(ie)]*(relel_gammas[ie]-1.);
  }

  u_relel = integrate_relel_q(q);  
  u_relel *= M_ELECTR;
  if (u_relel < 0.)
    u_relel=0.;
#endif
  return u_relel;
}


//Total pressure of relativistic electrons  
ldouble calc_relel_p(ldouble *pp)
{
  ldouble p_relel=0.;
#ifdef RELELECTRONS
  ldouble p2, p1;
  ldouble dgamma;
  int ie;
  ldouble q[NRELBIN];
  
  for (ie=0; ie<NRELBIN; ie++) 
  {
    q[ie] = pp[NEREL(ie)]*(relel_gammas[ie]-relel_gammas_inv[ie]);  
  }
  p_relel = integrate_relel_q(q);
  p_relel *= (M_ELECTR/3.);
  
  if (p_relel < 0.)
    p_relel=0.;
#endif
  return p_relel;
}

//*********************************************/
//* Determine gamma_inj_min and gamma_inj_max
//*********************************************/
//this is an approximate form for calculating gamma_inj_min
//it is exact when gamma_inj_max ->infinity and theta>1 
ldouble calc_gammainj_min_jointhermal(ldouble theta, ldouble delta_nth, ldouble p_index, ldouble gammamax)
{
  ldouble ratio=1., gamma_min=0;

#ifdef RELELECTRONS
#ifdef RELEL_HEAT_FIX_LIMITS
  gamma_min=RELEL_INJ_MIN;
#else

#if defined CALC_GAMMAMIN_FIT
  // fitting function to exact solution, works will all p>1
  // see mathematica notebook plaw_tail_full_simple.nb
 
  ldouble p0 = 2.2;
  ldouble y0 = 5.1;
  ldouble y1 = 2.01;
  ldouble al = 2.4;
  ldouble be = 0.6;
  ldouble de = 0.04;
  ldouble eta= delta_nth/(1.-delta_nth);
  ldouble gamma_pk = 2*theta;

  ldouble y;
  if(p_index < 1)
    my_err("in calc_gammainj_min, injection power law index must be >= 1!");
  else if(p_index < p0)
    y = y0 * pow(gammamax, de*(p0-p_index));
  else
  {
    ldouble pmax = al*pow(eta,-be);
    if(p_index<pmax)
    {
      ldouble aa = y1/y0;
      ldouble bb = (p0 - p_index) / (p0 - pmax);
      y = y0 * pow(aa, bb);
    }
    else
      y = y1;
  }

  gamma_min = y * gamma_pk;

#elif defined(CALC_GAMMAMIN_FAST) // assumes p > 2
  if(p_index <= 2.) my_err("in calc_gammainj_min, injection power law index must be > 2!");

  #if defined(RELEL_HEAT_FIX_FRAC) && defined(RELEL_HEAT_FIX_INDEX)
  ratio = gamma_inj_fixed_ratio;
  #else
  
  ldouble eta= delta_nth/(1.-delta_nth);
  ldouble lhs = 6.*(p_index-2.)*eta;
  ldouble a = 0.25 * pow(lhs, 0.25);
  ratio = -4.*gsl_sf_lambert_Wm1(-a);
  #endif
  
  gamma_min = ratio*theta;

#else 
  // iterative -- assumes p > 2
  ldouble ep=1.e-8;
  if(fabs(p_index-2.0) < ep)
    p_index = 2.0 + ep;

  ldouble eta= delta_nth/(1.-delta_nth);
  ldouble gmin_init = 10*theta;
  ldouble gamma_pk = 2*theta;
  ldouble gmin_guess = gmin_init;
  ldouble converge_crit = 0.05;
  int itermax = 25;
      
  ldouble pp2=2.+p_index;
  ldouble pm2=2.-p_index;
  ldouble prefac = 1./(6*eta*pm2*theta*theta*theta*theta);
  ldouble gmaxp2 =  pow(gammamax, pm2); 
  ldouble gminp4;
  
  ldouble rhs, gmin, relchange, deriv;
  int i;
  int sol_found=0;
  for(i=0; i<itermax; i++)
  {
    gminp4 = gmin_guess*gmin_guess*gmin_guess*gmin_guess;
    rhs = prefac*(gmaxp2*pow(gmin_guess,pp2) - gminp4);
    gmin = theta*log(rhs);
    relchange =  fabs(1.-gmin/gmin_guess);
    gmin_guess = gmin;

    if(gmin < 2*gamma_pk)
    {
      gmin_guess = gmin_init*2;
      gmin_init = gmin_guess;
    }
    
    if(relchange < converge_crit)
    {
      sol_found = 1;
      break;
    }
  }

  if(sol_found==1)
    gamma_min = gmin_guess;
  else
    gamma_min = gamma_pk;
  
#endif //CALC_GAMMAMIN_FIT
#endif //RELEL_HEAT_FIX_LIMITS
  
  //Floor & ceiling on gammamin
  if(gamma_min < RELEL_INJ_MIN)
    gamma_min = RELEL_INJ_MIN;
  if(gamma_min > RELEL_INJ_MAX)
    gamma_min = RELEL_INJ_MAX;

#endif //RELELECTRONS
  

  return gamma_min;
}

//calculate gamma_inj_max from synchrotron cooling timescale
ldouble calc_gammainj_max_syncool(ldouble bsq_cgs, ldouble dtau_cgs)
{

  ldouble gamma_max=1.;
#ifdef RELELECTRONS
  gamma_max = RELEL_INJ_MAX;
  
#ifndef RELEL_HEAT_FIX_LIMITS
#ifdef RELEL_SYN
  
  ldouble bsq=bsq_cgs;

  //synchrotron timescale set by equating cooling timescale to dtau
  gamma_max = 7.83994e8 / (bsq * dtau_cgs);
  if (gamma_max > RELEL_INJ_MAX)
     gamma_max = RELEL_INJ_MAX;

#endif //RELEL_SYN
#endif //RELEL_HEAT_FIX_LIMITS
#endif
  
  return gamma_max;
}


int reconnection_plaw_params_from_state(ldouble *pp, void *ggg, void *sss, ldouble* delta_back, ldouble* pindex_back)
{
  ldouble delta=0.;
  ldouble pindex=10.;

   struct geometry *geom
    = (struct geometry *) ggg;
   struct struct_of_state *state
    = (struct struct_of_state *) sss;

  
#ifdef RELELECTRONS
#if defined(RELEL_HEAT_RECONNECTION)
  //David's fit to the reconnection electron heating fraction
  //fit for data generated at Te=Ti as a function of sigma_w and beta
  //and mion = mproton
  
  // get beta = gas pressure / magn pressure
  // and sigma = magn energy density / enthalpy
  ldouble beta=0.;
  ldouble sigma=0.;
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
  sigma = bsq/rho;//enth_tot;
  betamax = 0.25/sigma;

  //ANDREW TODO LIMITS??
  if(sigma<1.e-10) sigma=1.e-10;
  if(!isfinite(enth_tot)) sigma=1.e-10;
  if(!isfinite(beta)) beta=1.e10; //no magnetic field
  if(!isfinite(betamax)) betamax=1.e10; //no magnetic field

  betanorm = beta/betamax;
  if(betanorm>1.) betanorm=1.;
  if(betanorm<0.) betanorm=0.;
  #endif

  //Andrew's fit to David's fit to simulation numbers
  //TODO find a better fit!!
  
  ldouble ap = 1.8 + 0.7/sqrt(sigma);
  ldouble bp = 3.7*pow(sigma, -0.19);
  ldouble cp = 23.4*pow(sigma, 0.26);
  pindex = ap + bp*tanh(cp*beta);

  ldouble ae = 1-1./(4.2*pow(sigma,0.55) + 1);
  ldouble be = 0.63*pow(sigma,0.07);
  ldouble ce = -68*pow(sigma,0.13);
  delta = ae + be*tanh(ce*beta);

  //TODO is this fraction of the total or the ratio?
  //delta = 1./((1+115.994*beta)*(1+0.2838673/sqrt(sigma)));
  //ldouble rtsigi = 1./sqrt(sigma);
  //pindex = 1. + 2*(1+3.65*rtsigi)/(1+pow(betanorm,0.422)) + 0.435*rtsigi;


  //ANDREW TODO CHECK FLOORS                                                                                                                                                                  
  if(!isfinite(delta))
    {
      //printf("problem with David delta fit: %d %d %d : %e %e %e %e %e\n",geom->ix,geom->iy,geom->iz,delta,beta,sigma,Te,Ti);                                                                    
      delta = 0;
      //print_primitives(pp);                                                                                                                                                                     
    }
  if(delta>1)
    delta=1.;
  if(delta<0)
    delta=0.;
  if(!isfinite(pindex));
   {
     //printf("problem with David pindex fit: %d %d %d : %e %e %e %e %e\n",geom->ix,geom->iy,geom->iz, delta,beta,sigma,Te,Ti);                                                                   
     pindex = 10.;
     delta=0.;
   }
  if(pindex>10)
    pindex=10;
  if(pindex<1.01)
    pindex=1.01;
#endif
#endif

  *pindex_back=pindex;
  *delta_back=delta;
  return 0;
}

//*********************************************/
//* Thermal Physics
//*********************************************/

//Calculate total gas adiabatic index with nonthermal population 
ldouble calc_gammaint_relel(ldouble* pp, ldouble Te, ldouble Ti)
{  
  ldouble gamma;
#ifndef RELELECTRONS
  gamma = calc_gammaintfromTei(Te,Ti);
#else
  ldouble nur=calc_relel_ne(pp);
  ldouble pur=calc_relel_p(pp);
  ldouble uur=calc_relel_uint(pp);
  if (nur==0. || pp[RHO]==0.) 
  {
    gamma = calc_gammaintfromTei(Te,Ti);
  }
  else 
  {
    ldouble the= kB_over_me*Te;
    ldouble gammae=GAMMAE;
    #ifndef FIXEDGAMMASPECIES
    gammae=calc_gammaintfromtheta(the);
    #endif
    ldouble peth=calc_thermal_ne(pp)*K_BOLTZ*Te;
    ldouble ueth=peth/(gammae-1.);

    ldouble thi= kB_over_mui_mp*Ti;
    ldouble gammai=GAMMAI;
    #ifndef FIXEDGAMMASPECIES
    gammai=calc_gammaintfromtheta(thi);
    #endif
    ldouble pi=kB_over_mui_mp*Ti*pp[RHO];
    ldouble ui=pi/(gammai-1.);

    gamma = 1. + (peth+pi+pur)/(ueth+ui+uur);
  }
#endif
  return gamma;
}


//ANDREW -- do without pp, or precompute uur and nur?
//calc ugas from Te and Ti using info in pp for density and energy density
//pp[RHO] and pp[NEREL(ie)] must be already set to use this function!
ldouble calc_PEQ_ugasfrom_Tei_relel(ldouble *pp, ldouble Te,ldouble Ti)
{
  ldouble uur = calc_relel_uint(pp);
  ldouble gammae=GAMMAE;
  ldouble gammai=GAMMAI;
  #ifndef FIXEDGAMMASPECIES
  gammae=calc_gammaintfromtemp(Te,ELECTRONS);
  gammai=calc_gammaintfromtemp(Ti,IONS);
  #endif
  ldouble ueth = K_BOLTZ*calc_thermal_ne(pp)*Te/(gammae - 1.);
  ldouble uith = kB_over_mui_mp*pp[RHO]* Ti/(gammai - 1.);
  ldouble u = uur + ueth + uith;
  return u;
}


// Calculate and return the chemical potential for a given theta and thermal number density
// Expression  consistent with the entropy S3. 
ldouble
chemical_potential_short(ldouble theta, ldouble neth)
{
  ldouble mu;
  mu = M_ELECTR * (1. - 0.6*log1p(2.5*theta) + theta*(4-1.5*log(theta*(theta+0.4)) + log(neth)));
  return mu;
}


// Calculate chemical potential from entropy density (s*n)
// pressure, energy density, number density, temperature
ldouble
chemical_potential_long(ldouble p, ldouble u, ldouble sn, ldouble neth, ldouble Te)
{
  return ((p + u - Te*sn)/neth);
}


//*********************************************/
//* Injection
//*********************************************/

//Apply viscous heating / particle injection for relativistic distribution in a single cell
int
apply_relel_visc_heating(ldouble pp[NV], ldouble durelel, ldouble p_index, ldouble injmin, ldouble injmax, ldouble dtau)
{
#ifdef RELELECTRONS
   int ie;
   ldouble xx,pref, gamma, neur, scale;
   ldouble qheat_app_tot;   
   ldouble qheat[NRELBIN];
   ldouble nheat[NRELBIN];
   
   // total relel energy before heating
   ldouble uint0 = calc_relel_uint(pp); 
   ldouble nenth = calc_relel_ne(pp);

   //Skip heating if it would give negative net energy in the distribution, or if already at maximum # 
   if((uint0 + durelel) > 0. && nenth/(pp[RHO]*one_over_mue_mp) < MAX_RELEL_FRAC_N)
   {     
      // Prefactor (precomputed if constant using constant heating params)
      // ANDREW TODO case where 1<p_index<2
      #if defined(RELEL_HEAT_FIX_LIMITS) && defined(RELEL_HEAT_FIX_INDEX)
      xx=xx_relel_inj;
      #else
      if (p_index!=2 && p_index!=1)
      {
      xx = (pow(injmax, 2.-p_index) - pow(injmin, 2.-p_index))/(2.-p_index) - (pow(injmax, 1.-p_index) - pow(injmin, 1.-p_index))/(1.-p_index);
      }
      else if (p_index == 2)
      {
      xx = 1./injmax - 1./injmin + log(injmax/injmin);
      }
      else if (p_index == 1)
      {
      xx = injmax - injmin - log(injmax/injmin);
      }

      
      #endif
      pref = durelel / (xx * M_ELECTR);
      //printf("%e %e\n",numdensGU2CGS(pref), numdensGU2CGS(pref)/timeGU2CGS(dtau));
      //exit(-1);
      //Calculate and heating rate in all bins, and the total energy added
      qheat_app_tot=0.;
      for (ie=0; ie<NRELBIN; ie++) 
      {
        gamma = relel_gammas[ie];

        if ((gamma < injmin) || (gamma > injmax)) 
        {
         nheat[ie] = 0.;
        }
        else 
        {
	 #if defined(RELEL_HEAT_FIX_LIMITS) && defined(RELEL_HEAT_FIX_INDEX)
	 nheat[ie] = pref * relel_injpow[ie];
	 #else 
         nheat[ie] = pref * pow(gamma, -p_index);
         #endif
	}
	
	qheat[ie] = (gamma-1.)*nheat[ie];
      }
      qheat_app_tot = M_ELECTR*integrate_relel_q(qheat);

      //Normalization
      //Scale so that durelel is total numerically integrated energy added
      //And apply the change in each bin

#ifndef NORELELINJSCALE
      if (fabs(qheat_app_tot) > 0.) 
      {
        scale = durelel / qheat_app_tot;
        if(scale < 0.) scale = 0.;
        for (ie=0; ie<NRELBIN; ie++) 
        {
	  pp[NEREL(ie)] += scale*nheat[ie];
        }
      }
      //printf("%e\n",scale);
#endif
      
      // Ceiling: Not too many rel. electrons
      nenth = calc_relel_ne(pp);
      if (nenth/(pp[RHO]*one_over_mue_mp) > MAX_RELEL_FRAC_N) 
      { 
        for (ie=0; ie<NRELBIN; ie++)
     	{
	  ldouble scale2 = MAX_RELEL_FRAC_N/(nenth/(pp[RHO]*one_over_mue_mp));
       	  pp[NEREL(ie)] *= scale2;
       	}
      }
      
      // Floor: No negative rel. electron bin numbers      
      for (ie=0; ie<NRELBIN; ie++)
      {
	if (pp[NEREL(ie)] < 0.) 
          pp[NEREL(ie)] = 0.;
      }
   }
#endif  
  return 0;
}


//*********************************************/
//* Adiabatic Cooling
//*********************************************/

//Apply adiabatic compression/expansion to relativistic electron distribution
//in all bins
ldouble
relel_adiab(ldouble dt)
{
#ifdef RELELECTRONS
  int ii;

  #pragma omp parallel for private(ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      int ix,iy,iz,iv;
      int i,ie,ig;
      ldouble pp[NV];
      ldouble pp0[NV];

      ldouble uu[NV];

      ldouble ngammas[NRELBIN];
      ldouble ngammas0[NRELBIN];
      ldouble q[NRELBIN]; //derivative 
      ldouble napp[NRELBIN];
      ldouble qapp[NRELBIN];
      ldouble theta;

#if defined(RELEL_ADIAB_LOGSPACE) || defined(RELEL_ADIAB_LOGSPACE_LF) 
      ldouble x[NRELBIN]; //gamma*n(gamma)
      ldouble x0[NRELBIN];
#endif
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      if(is_cell_active(ix,iy,iz)==0) 
	continue;

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom); 
      
      //Calculate proper time step
      ldouble vcon[4];
      ldouble ucon[4];
      ldouble ucov[4];
      ldouble dtau;
      vcon[1]=get_u(p, VX, ix, iy, iz);
      vcon[2]=get_u(p, VY, ix, iy, iz);
      vcon[3]=get_u(p, VZ, ix, iy, iz);
      vcon[0]=0.;
      conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,geom.gg,geom.GG);      
      dtau = dt/ucon[0];
      
      //Fill vectors of primitives & conserved quantities
      for (i=0; i < NV; i++) 
      {
        pp0[i] = get_u(p, i, ix, iy, iz);
        pp[i] = pp0[i];    
        uu[i] = get_u(u, i, ix, iy, iz);
      }
      
      //Calculate thermal electron number and energy before adiabatic compression/expansion
      ldouble ueth, neth;
      
      neth = calc_thermal_ne(pp);
      //#ifdef RELELENTROPY
      //ueth=calc_ufromS4n(pp[ENTRE],neth,ELECTRONS,ix,iy,iz);
      //#else
      ueth=calc_ufromSerho(pp[ENTRE],neth*MU_E*M_PROTON,ELECTRONS,ix,iy,iz);
      //#endif
      
      //Calculate expansion four-divergence
      ldouble div;
#ifdef RELEL_ADIAB_ART
      div = RELEL_DIV_ART / timeCGS2GU(1.);
#else
      int derdir[3] = {0,0,0};
      calc_fluid_div_lab(pp, &geom, dt, &div, MHD, derdir);
      //ldouble shear[4][4];      
      //calc_shear_lab(pp, &geom, shear, &div, MHD, derdir);
#endif

      //Constant velocity in gamma space, ignore the 1/gamma (gamma >> 1)
      theta = one_third*div;

      // Make temporary arrays
      for (ie=0; ie<NRELBIN; ie++) 
      {  
	ngammas0[ie] = pp0[NEREL(ie)];
	ngammas[ie] = ngammas0[ie];
#if defined(RELEL_ADIAB_LOGSPACE) || defined(RELEL_ADIAB_LOGSPACE_LF)
        x0[ie] = relel_gammas[ie]*ngammas0[ie];
        x[ie] = x0[ie];
#endif 
        //bin number change expected
        qapp[ie] = -1*dtau*M_ELECTR*theta*(relel_gammas[ie]-relel_gammas_inv[ie])*ngammas0[ie];
      }

      //compute expected total energy loss or gain
      ldouble de1 = 0.;
      ldouble de2 = 0.;
      de1 = integrate_relel_q(qapp);
	
      //Compute derivative with finite diff and apply
#if defined(RELEL_ADIAB_LOGSPACE)
      relel_adiab_derivs_logspace(theta, dtau, x, q);  
      for(ie=0; ie<NRELBIN; ie++) 
	 x[ie] = x0[ie] + dtau*q[ie];
      
#elif defined(RELEL_ADIAB_LOGSPACE_LF)
      relel_adiab_derivs_logspace_lf(theta, dtau, ngammas, q);
      for(ie=0; ie<NRELBIN; ie++)
	x[ie] = x0[ie] +  dtau*q[ie];
#else
      
      relel_adiab_derivs(theta, dtau, ngammas, q);
      for(ie=0; ie<NRELBIN; ie++)
	 ngammas[ie] = ngammas0[ie] + dtau*q[ie];
#endif 

      /*
      //2nd order RK Heun
#ifdef RELEL_ADIAB_RK2
      ldouble q2[NRELBIN];
#if defined(RELEL_ADIAB_LOGSPACE)
      relel_adiab_derivs_logspace(theta, dtau, x, q2);
      for (ie=0; ie<(NRELBIN); ie++) 
        x[ie] = x0[ie] + 0.5 * dtau * (q[ie] + q2[ie]);
#elif defined(RELEL_ADIAB_LOGSPACE_LF)
      for (ie=0; ie<(NRELBIN); ie++) 
         ngammas[ie] = x[ie]*relel_gammas_inv[ie];
      relel_adiab_derivs_logspace_lf(theta, dtau, ngammas, q2);
      for (ie=0; ie<(NRELBIN); ie++) 
         x[ie] = x0[ie] + 0.5 * dtau * (q[ie] + q2[ie]);
#else
      relel_adiab_derivs(theta, dtau, ngammas, q2);
      for (ie=0; ie<(NRELBIN); ie++) 
        ngammas[ie] = ngammas0[ie] + 0.5 * dtau * (q[ie] + q2[ie]);
#endif
#endif
      */
      
      //Determine final heating rate and total energy change
      for (ie=0; ie<NRELBIN; ie++)
      {
	
#if defined(RELEL_ADIAB_LOGSPACE) || defined(RELEL_ADIAB_LOGSPACE_LF)
	napp[ie] = x[ie]*relel_gammas_inv[ie] - pp0[NEREL(ie)];
#else
	napp[ie] = ngammas[ie] - pp0[NEREL(ie)];  
#endif
	qapp[ie] = M_ELECTR*(relel_gammas[ie]-1.)*napp[ie];
      }

      de2 = integrate_relel_q(qapp);
      
      //Compute thermal energy gain from bottom or top bin
      ldouble dueth=0.;
      ldouble dneth=0.;
      ldouble Senew;
      ldouble gdot;

      //Expansion
      if (relel_gammas[0] != 1. && div > 0)
      {
        gdot = theta*(relel_gammas_e[0] - relel_gammas_e_inv[0]);
        dueth = fabs(M_ELECTR*gdot*ngammas0[0]*(relel_gammas_e[0] - 1.)*dtau);
      }
      //Compression
      else if(div < 0.)
      {
        gdot = theta*(relel_gammas_e[NRELBIN] - relel_gammas_e_inv[NRELBIN]);
        dueth = fabs(M_ELECTR*gdot*ngammas0[NRELBIN-1]*(relel_gammas_e[NRELBIN] - 1.)*dtau);
      }

      //NORMALIZATION
#ifndef NO_RELEL_ADIAB_SCALE      
      ldouble frac= de1/de2;
      if(calc_relel_uint(pp)>0. && fabs(de2)>0.)
      {
	for(ie=0;ie<NRELBIN;ie++)
	  napp[ie] *= frac;
      }
#endif
		
      //Apply the adiabatic change to the bins
      for(ie=0;ie<NRELBIN;ie++)
      {
	pp[NEREL(ie)] = pp0[NEREL(ie)] + napp[ie];

        //Floor 
        if(pp[NEREL(ie)]<0.) pp[NEREL(ie)]=0.;
      }
      
      //Apply the combination of edge heating and zeroed bins to thermal electrons
      dneth = calc_thermal_ne(pp) - neth; //total change in thermal particle number
      apply_du_dn_2_species(pp, ueth, neth, dueth, dneth, &geom, ELECTRONS, &Senew); 
      pp[ENTRE] = Senew;
      
      //Convert primitive to conserved and save
      p2u_mhd(pp, uu, &geom);
      for (i=0; i<NV; i++) 
      { 
        set_u(u, i, ix, iy, iz, uu[i]);
        set_u(p, i, ix, iy, iz, pp[i]);
      } 
  }
#endif
  return 0;
}

//Zero out nonthermal bins below thermal peak and move energy to thermal.
//pre-calculate thetae and neth (from state)
ldouble remove_lowgamma_electrons(ldouble thetae, ldouble neth, ldouble *pp)
{
  ldouble du=0.;
#ifdef RELELECTRONS
  //ldouble Te,thetae;
  ldouble pref;
  ldouble uintold;
  int i,ie;

  //energy in rel electrons before zeroing
  uintold = calc_relel_uint(pp);
 
  //zero out all bins below the approximate location of the thermal peak
  ldouble pkgamma = 2.*thetae; //high theta approx
  int pkbin = NRELBIN;
  int zeroed_bins = 0;  
  for (ie=0;ie<NRELBIN;ie++)
  {
    if(relel_gammas[ie] >= pkgamma) //approx location of the peak of the thermal distribution function
      {
	pkbin = ie;
	break;
      }
    else
      {
      if (pp[NEREL(ie)] != 0.)
        {
	  pp[NEREL(ie)] = 0.;
	  zeroed_bins++;
        }
      }
  }

  //above the thermal peak, zero out bins where n_th > n_neth
  if (pkbin<NRELBIN)
  {
    ldouble gamma,ngamma_th,ngamma_nth;
    ldouble pref = 1./(2.*thetae*thetae*thetae) + 1./(8.*thetae*thetae*thetae*thetae*thetae);
    for(ie=pkbin;ie<NRELBIN;ie++)
    {
	gamma = relel_gammas[ie];
	ngamma_th = neth*pref*gamma*gamma*exp(-gamma/thetae); //assuming theta >> 1
	ngamma_nth = pp[NEREL(ie)];
	if (ngamma_nth > ngamma_th)
	{
	  break;
	}
	else
	{
	  pp[NEREL(ie)] = 0.;
	  zeroed_bins++;
	}
    }
  }

  //calculate the missing energy to add to the thermal electron distribution
  if (zeroed_bins>0)
  {
    du = fabs(uintold - calc_relel_uint(pp));
  }
#endif //RELELECTRONS
  
  return du;
  
}

///////////////////////////////////////////////////////////////////////////
//Derivative methods for explicit relel adiabatic evolution
///////////////////////////////////////////////////////////////////////////
//Compute dn/dtau using given array of ngamma values and put derivatives in q
int
relel_adiab_derivs(ldouble theta, ldouble dtau, ldouble *ngammas, ldouble *q)
{
#ifdef RELELECTRONS
  
  int ie; 
  ldouble ngdown, ngup, ngup2;
  ldouble gdown,gup;

  for (ie=0; ie<(NRELBIN); ie++) 
  {
    // Current bin is always downwind 
    ngdown = theta*(relel_gammas[ie]-relel_gammas_inv[ie])*ngammas[ie];
    gdown = relel_gammas[ie];
   
    // Expansion
    if (theta > 0.) 
    {
      if (ie==NRELBIN-1) 
      {
        ngup = 0.;
	gup = relel_gammas[NRELBIN-1]*binspace;
      } 
      else 
      {
        ngup = theta*(relel_gammas[ie+1]-relel_gammas_inv[ie+1])*ngammas[ie+1];
        gup = relel_gammas[ie+1];
      }
    }
    //Compression
    else if (theta < 0.)
    {
      if (ie==0) 
      {
        ngup=0.;
	gup=relel_gammas[0]/binspace;
      } 
      else 
      {
        ngup = theta*(relel_gammas[ie-1]-relel_gammas_inv[ie-1])*ngammas[ie-1];
        gup = relel_gammas[ie-1];
      }
    }
   
    //Derivative
    q[ie] = (ngup - ngdown) / (gup - gdown); 
  }
  
#endif
  
  return 0;
}

// Compute dx/dtau with the given gdot = u^mu_;mu/3 and put the answer in q
// Takes advantage of equal log spacing
// Uses Lax-Frederichs instead of upwind finite diff
int
relel_adiab_derivs_logspace_lf(ldouble theta, ldouble dtau, ldouble *n, ldouble *q)
{
#ifdef RELELECTRONS
  int ie;
  ldouble n_l[NRELBIN]; //primitives at cell left boundaries
  ldouble n_r[NRELBIN]; //primitives at cell right boundaries
  ldouble f[NRELBIN+1]; //fluxes at cell boundaries
  
  //linear interpolation using minmod theta flux limiter to cell walls
  ldouble mmtheta=RELEL_MINMOD_THETA;
  ldouble np1, nm1, n0;
  ldouble slope, deltam, deltap;
  ldouble uL,uR,fL, fR, wspeed;
  
  for (ie=0; ie<NRELBIN; ie++) //Loop over cells and interpolate n
  {
    n0 = n[ie];
    
    //Ghost cells -- all values zero  
    if (ie==0) nm1=0.;
    else nm1=n[ie-1];

    if (ie==NRELBIN-1) np1=0.;
    else np1 = n[ie+1];

    //left and right slopes
    deltam = n0-nm1;
    deltap = np1-n0;
     
    if (deltam * deltap <= 0.)
    {
      // We are at a local maximum or minimum. Use zero slope (i.e., donor cell)
      n_l[ie] = n0;
      n_r[ie] = n0;
    }
    else
    {
      if (deltam > 0.)
      {
      // Slopes are positive. 
         slope = my_min(my_min(mmtheta*deltam, 0.5*(deltam+deltap)), mmtheta*deltap); // theta=1 is MinMod, theta=2 is MC
      }
      else
      {
      // Slopes are negative.
         slope = my_max(my_max(mmtheta*deltam, 0.5*(deltam+deltap)), mmtheta*deltap); // theta=1 is MinMod, theta=2 is MC
      }
        
      n_r[ie] = n0 + 0.5*slope;
      n_l[ie] = n0 - 0.5*slope;

      if (n_r[ie] < 0.) n_r[ie]=0.;
      if (n_l[ie] < 0.) n_l[ie]=0.;
			
    }
  }

  //calculate LF flux on the cell boundaries 
  for(ie=0;ie<NRELBIN+1;ie++)
  {
    if(ie==0)
    {
      fL=0.;
      uL=0.;
    }
    else
    {
      uL = n_r[ie-1]*relel_gammas_e[ie]; //logspace conserved is gamma*n
      //fL = -theta*uL; //flux on the left of boundary is from the cell right face
      fL = -theta*uL*(1. - relel_gammas_e_inv[ie]*relel_gammas_e_inv[ie]);
    }
    
    if(ie==NRELBIN)
    {
      fR=0.;
      uR=0.;
    }
    else
    {
      uR = n_l[ie]*relel_gammas_e[ie]; //logspace conserved is gamma*n
      //fR = -theta*uR; //flux on the right of boundary is from cell left face
      fR = -theta*uR*(1. - relel_gammas_e_inv[ie]*relel_gammas_e_inv[ie]);
    }

    //which type  of LF do we do? 
    wspeed= fabs(theta); //reverts to upwind in appropriate direction
    //wspeed = logbinspace/dtau; //bounded by Courant
    f[ie] = 0.5*(fL + fR) - 0.5*wspeed*(uR - uL);
    //printf("%d %e %e %e \n",ie,fL,fR,f[ie]);
  }

  //calculates derivatives at cell
  for(ie=0; ie<NRELBIN; ie++)
  {
    q[ie] = (f[ie] - f[ie+1]) * logbinspace_inv;
  }
#endif
  return 0; 

}

// Compute dx/dtau with the given gdot = u^mu_;mu/3 and put the answer in q
// Takes advantage of equal log spacing
int
relel_adiab_derivs_logspace(ldouble theta, ldouble dtau, ldouble *x, ldouble *q)
{
  #ifdef RELELECTRONS
  int ie;
  ldouble xdown,xup,sign,gamma,xddown;
  
  for (ie=0; ie<(NRELBIN); ie++) 
  {
    // Current bin is always downwind 
    xdown = x[ie];
    gamma = relel_gammas[ie];
    
    //Expansion
    if (theta > 0.) 
    {
     sign = 1.;
     //Boundary conditions/ghost cells
     if (ie==NRELBIN-1) 
     {
       xup=0.; 
     }
     else if (ie==NRELBIN-2)
     {
       xup=x[ie+1];       
     }
     else 
     {
       xup = x[ie+1];
     }
    }

    //Compression
    else if (theta < 0.) 
    {
      sign = -1;
      //Boundary conditions/ghost cells
      if (ie==0) 
      {
        xup=0.;
      } 
      else if (ie==1)
      {
        xup=x[ie-1];
      }
      else 
      {
        xup = x[ie-1];
      }
    }
    //Derivative
    q[ie] = sign *(theta*logbinspace_inv)*(xup - xdown);

   } 
#endif
   return 0; 

}


//*********************************************/
//* Nonthermal Cooling Rates
//*********************************************/

ldouble
gdot_syn(ldouble gamma, ldouble bsq_cgs)
{
  ldouble bs;
  ldouble g2b2;
  g2b2 = gamma*gamma - 1.;
  bs = 1.292e-9 * bsq_cgs;
  return bs * g2b2 * timegu2cgs;
}


ldouble
gdot_ff(ldouble gamma, ldouble nion_cgs)
{
   ldouble bff;
   bff = 1.37e-16 * nion_cgs;
   return bff * gamma * (log(gamma) + 0.36) * timegu2cgs;
}


ldouble
gdot_ic(ldouble gamma, ldouble trad_cgs, ldouble eradhat_cgs)
{
  ldouble bic1;
  ldouble bic2;
  ldouble g2b2;
  ldouble out;
  g2b2 = gamma*gamma - 1.;
  bic1 = 3.248e-8 * eradhat_cgs;
  //printf("%e %e \n",eradhat_cgs, bic1);
  //exit(-1);
  bic2 = 11.2 * K_BOLTZ_CGS * trad_cgs / (M_ELECTR_CGS * CCC_CGS * CCC_CGS);
  out = bic1 * g2b2 * pow((1. + bic2 * gamma), -1.5) * timegu2cgs;
  return out;
} 


ldouble
gdot_coul(ldouble gamma, ldouble neth_cgs)
{
  ldouble bc1;
  ldouble bc2;
  
  bc1 = 1.491e-14 * neth_cgs;
  if (neth_cgs > 0)
    bc2 = -log(neth_cgs) + 74.7; //ANDREW is the 74.7 correct?
  else 
    bc2 = 0.;
  return bc1 * (log(gamma) + bc2) * timegu2cgs;
} 

// Advection term from Wong, Zhdankin, + 2020
// Cannot use in upwind, since sign changes.
// TODO dependence on other parameters!
ldouble gdot_turb_advect(ldouble gamma)
{
  ldouble a = 545. - 172.*log(1.+gamma/33.);
  return -1*a; // sign convention 
}
// advection from zhdankin 2019
ldouble gdot_z19_advect(ldouble gamma)
{
#ifdef RELEL_ADVECT_Z19
  ldouble D  =RELEL_RATE2*(RELEL_GAMMA0*RELEL_GAMMA0 + gamma*gamma);
  ldouble Ap =RELEL_RATEH*RELEL_GAMMA0 + RELEL_RATEA*gamma;
  ldouble Acool= -gamma*gamma / RELEL_GAMMA0;
  return -1*(2.*D/gamma + Ap + Acool)/RELEL_TAUC;
#else
  return 0.;
#endif
}


// Diffusion term from Wong, Zhdankin, + 2020
// TODO dependence on other parameters!
ldouble d_turb_diffuse(ldouble gamma)
{
  ldouble d = 920.*pow(gamma,2./3.) + 0.055*gamma*gamma;
  return d;
}
ldouble d_z19_diffuse(ldouble gamma)
{
#ifdef RELEL_DIFFUSE_Z19
  ldouble D=RELEL_RATE2*(RELEL_GAMMA0*RELEL_GAMMA0 + gamma*gamma);
  return D/RELEL_TAUC;
#else
  return 0.;
#endif
}
//*********************************************/
//* Nonthermal Error Function
//*********************************************/
int
calc_relel_f_and_fmag_from_state(ldouble *pp, void *sss, ldouble *pp0, void *ggg, ldouble dtau,ldouble *frel,ldouble *frelmag)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  struct struct_of_state *state
    = (struct struct_of_state *) sss;
#ifdef RELELECTRONS  
  ldouble qcool[NRELBIN];
  int fac=1;
  int ie;
  int done=-1;
  
#ifdef RELEL_IMPLICIT_LOGSPACE_LF 
  done = calc_relel_cooling_lf_from_state(pp, state, pp0, dtau, qcool); //Lax-Freidrichs
#endif

  if (done==-1)
  {
    printf("fallback to UPWIND!");
    calc_relel_cooling_from_state(pp, state, pp0, dtau, qcool); //simple upwind 
  }
  
#ifdef NORELELCOOLATBH
  ldouble xxBL[4];
  #ifdef PRECOMPUTE_MY2OUT
  get_xxout(geom->ix, geom->iy, geom->iz, xxBL);
  #else
  coco_N(geom->xxvec,xxBL,MYCOORDS,BLCOORDS);
  #endif
  if(xxBL[1]<=rhorizonBL) fac=0;
#endif     

  // compute f and fmag 
  for (ie=0; ie < NRELBIN; ie++)
    {  
      frel[ie] = pp[NEREL(ie)] - pp0[NEREL(ie)] - dtau*fac*qcool[ie];
      frelmag[ie] = fabs(pp[NEREL(ie)]) + fabs(pp0[NEREL(ie)]) + fabs(dtau*(fac*qcool[ie]));
    }
#endif //RELELECTRONS  
  return 0;
}

//calculate relel cooling rate from upwind finite differencing 
//struct_of_state should already be computed
int
calc_relel_cooling_from_state(ldouble *pp, void *sss, ldouble *pp0, ldouble dtau, ldouble *qcool)
{
  struct struct_of_state *state
    = (struct struct_of_state *) sss;
#ifdef RELELECTRONS
  ldouble gdots[NRELBIN];

  //all in cgs
  ldouble bsq_cgs=0.; //gauss
  ldouble nion_cgs=0.; //cm^-3
  ldouble neth_cgs=0.; //cm^-3
  ldouble trad_cgs=0.; //K
  ldouble erad_cgs=0.; //erg cm^-3 // in fluid frame!

  // synchrotron cooling rate
  #ifdef RELEL_SYN
  #ifdef RELEL_SYN_ART
  bsq_cgs = RELEL_B_ART*RELEL_B_ART;
  #else
  bsq_cgs = (state->bsq)*fourmpi*endengu2cgs;
  #endif
  #endif

  // free-free cooling      
  #ifdef RELEL_FF
  #ifdef RELEL_FF_ART
  nion_cgs = RELEL_NION_ART;
  #else
  nion_cgs = (state->ni) * numdensgu2cgs;
  #endif
  #endif

  // inverse compton cooling
  #ifdef RELEL_IC
  #ifdef RELEL_IC_ART
  trad_cgs = RELEL_TRAD_ART;
  erad_cgs = RELEL_ERAD_ART;
  #else 
  trad_cgs = state->Trad;
  erad_cgs = (state->Ehat) * endengu2cgs;
  #endif
  #endif
  
  // coloumb cooling
  #ifdef RELEL_COUL
  #ifdef RELEL_COUL_ART
  neth_cgs = RELEL_NETH_ART;
  #else 
  neth_cgs = (state->ne)*numdensgu2cgs;
  #endif
  #endif
  
  // Calculate gamma dots
  ldouble gdot_bin, gamma,gdottmp; 
  int ie;
  
  for (ie=0; ie<NRELBIN; ie++) 
  {
      gamma = relel_gammas[ie];
      gdot_bin = 0.;

      #ifdef RELEL_SYN
      gdot_bin += gdot_syn(gamma,bsq_cgs);
      #endif
      #ifdef RELEL_FF
      gdot_bin += gdot_ff(gamma,nion_cgs);
      #endif
      #ifdef RELEL_IC
      gdot_bin += gdot_ic(gamma, trad_cgs, erad_cgs);
      #endif
      #ifdef RELEL_COUL
      gdot_bin += gdot_coul(gamma, neth_cgs);
      #endif
      
      gdots[ie] = gdot_bin;
  }
  
  ldouble gdot_plus_1, n_plus_1;
  //Direct Upwind in Log Space
  for (ie=0; ie<(NRELBIN); ie++) 
  {      
    
    //Rightmost Edge
    if (ie == NRELBIN-1)
    {
      gdot_plus_1 = 0.; 
      n_plus_1=0.;
    } 
    else
    {
      gdot_plus_1=gdots[ie+1];
      n_plus_1=pp[NEREL(ie+1)];
    }

    qcool[ie] = (gdot_plus_1*n_plus_1 - gdots[ie]*pp[NEREL(ie)]) / (logbinspace * relel_gammas[ie]); //* logbinspace_inv * relel_gammas_inv[ie]; 

  }
#endif
  return 0;
}


//calculate relel cooling rate from upwind lax-freidrichs 
//struct_of_state should already be computed
int
calc_relel_cooling_lf_from_state(ldouble *pp, void *sss, ldouble *pp0, ldouble dtau, ldouble *qcool)
{
  struct struct_of_state *state
    = (struct struct_of_state *) sss;

  int out = 0;
  
  //all in cgs
  ldouble bsq_cgs=0.; //gauss
  ldouble nion_cgs=0.; //cm^-3
  ldouble neth_cgs=0.; //cm^-3
  ldouble trad_cgs=0.; //K
  ldouble erad_cgs=0.; //erg cm^-3 // in fluid frame!
  
#ifdef RELELECTRONS
  // synchrotron cooling rate
  #ifdef RELEL_SYN
  #ifdef RELEL_SYN_ART
  bsq_cgs = RELEL_B_ART*RELEL_B_ART;
  #else
  bsq_cgs = (state->bsq)*fourmpi*endengu2cgs;
  #endif
  #endif

  // free-free cooling      
  #ifdef RELEL_FF
  #ifdef RELEL_FF_ART
  nion_cgs = RELEL_NION_ART;
  #else
  nion_cgs = (state->ni) * numdensgu2cgs;
  #endif
  #endif

  // inverse compton cooling
  #ifdef RELEL_IC
  #ifdef RELEL_IC_ART
  trad_cgs = RELEL_TRAD_ART;
  erad_cgs = RELEL_ERAD_ART;
  #else 
  trad_cgs = state->Trad;
  erad_cgs = (state->Ehat) * endengu2cgs;
  #endif
  #endif
  
  // coloumb cooling
  #ifdef RELEL_COUL
  #ifdef RELEL_COUL_ART
  neth_cgs = RELEL_NETH_ART;
  #else 
  neth_cgs = (state->ne)*numdensgu2cgs;
  #endif
  #endif

  int ie;
  ldouble n_l[NRELBIN+1]; //primitives at cell left boundaries
  ldouble n_r[NRELBIN+1]; //primitives at cell right boundaries
  ldouble f[NRELBIN+1]; //fluxes at cell boundaries
  
  //linear interpolation using minmod theta flux limiter to cell walls
  ldouble mmtheta = RELEL_MINMOD_THETA;
  ldouble np1, nm1, n0, nr, nl;
  ldouble slope, deltam, deltap;
  ldouble uL,uR,fL, fR, wspeed,wspeed2;
  ldouble gdot_wall, gamma,aedge;

  // Calculates flux-limited bin values at the cell walls
  for (ie=0; ie<NRELBIN; ie++)
  {
    n0 = pp[NEREL(ie)]*relel_gammas[ie]; // the quantity we interpolate is n(gamma)*gamma
    
    //Ghost cells -- all values zero  
    if (ie==0)
      nm1=0.;
    else
      nm1=pp[NEREL(ie-1)]*relel_gammas[ie-1];

    if (ie==NRELBIN-1)
      np1=0.;
    else
      np1=pp[NEREL(ie+1)]*relel_gammas[ie+1];

    //left and right slopes
    deltam = n0 - nm1;
    deltap = np1 - n0;
     
    if (deltam * deltap <= 0.)
    {
      // We are at a local maximum or minimum. Use zero slope (i.e., donor cell)
      nl = n0;
      nr = n0;
    }
    else
    {
      if (deltam > 0.)
      {
         // Slopes are positive. 
         //slope = my_min(my_min(mmtheta*deltam, 0.5*(deltam+deltap)), mmtheta*deltap); // theta=1 is MinMod, theta=2 is MC
	 slope = my_max(my_min(2*deltam, deltap), my_min(deltam, 2*deltap)); // SUPERBEE
      }
      else
      {
         // Slopes are negative.
         //slope = my_max(my_max(mmtheta*deltam, 0.5*(deltam+deltap)), mmtheta*deltap); // theta=1 is MinMod, theta=2 is MC
	 slope = my_min(my_max(2*deltam, deltap), my_max(deltam, 2*deltap)); // SUPERBEE
      }
        
      nr = n0 + 0.5*slope;
      nl = n0 - 0.5*slope;			
    }

    n_r[ie] = nl; //left wall from cell perspective is right wall from edge perspective
    n_l[ie+1] = nr;

    if (n_r[ie] < 0.) n_r[ie]=0.;
    if (n_l[ie] < 0.) n_l[ie]=0.;
  }

  // Wall values at the right and left boundaries
  n_r[NRELBIN] = n_l[NRELBIN]; //TODO -- just copy rightmost wall? 
  n_l[0] = n_r[0]; // TODO -- just copy leftmost wall? 
  
  // Calculate LF fluxes on all the walls 
  for(ie=0; ie<NRELBIN+1; ie++)
  {
    gamma = relel_gammas_e[ie];
    gdot_wall = 0.;

    #ifdef RELEL_SYN
    gdot_wall += gdot_syn(gamma,bsq_cgs);
    #endif
    #ifdef RELEL_FF
    gdot_wall += gdot_ff(gamma,nion_cgs);
    #endif
    #ifdef RELEL_IC
    gdot_wall += gdot_ic(gamma, trad_cgs, erad_cgs);
    #endif
    #ifdef RELEL_COUL
    gdot_wall += gdot_coul(gamma, neth_cgs);
    #endif
    #ifdef RELEL_TURB_ADVECT
    gdot_wall = gdot_turb_advect(gamma); // TODO -- dependence on parameters
    #endif
    #ifdef RELEL_ADVECT_Z19
    gdot_wall = gdot_z19_advect(gamma);
    #endif
    
    aedge = -gdot_wall / gamma;      

    fR = n_r[ie]*aedge;
    fL = n_l[ie]*aedge;

    f[ie] = 0.5 * (fR + fL - fabs(aedge) * (n_r[ie] - n_l[ie]));

    // Diffusion part
    #ifdef RELEL_DIFFUSE
    ldouble D_edge, nplus, nminus;
    D_edge=0.;
    #ifdef RELEL_TURB_DIFFUSE
    D_edge = d_turb_diffuse(gamma) / gamma; // TODO -- dependence on parameters
    #endif
    #ifdef RELEL_DIFFUSE_Z19
    D_edge = d_z19_diffuse(gamma) / gamma;
    #endif
    
    if(ie==0)
    {
      nplus = pp[NEREL(0)]; // no gamma term in diffusion part
      nminus = 0.;
    }
    else if(ie==NRELBIN)
    {
      nplus=0;
      nminus=pp[NEREL(NRELBIN-1)];
    }
    else
    {
      nplus=pp[NEREL(ie)];
      nminus=pp[NEREL(ie-1)];
    }

    ldouble f_diffusion = -D_edge * (nplus - nminus) * logbinspace_inv;
    f[ie] += f_diffusion;
    #endif // RELEL_DIFFUSE 

  }

  
  // Apply Neumann BCs
  // TODO -- dirichlet? 
  f[0] = 0.;
  f[NRELBIN] = 0.;
  
  //printf("f[0] , f[1] %e %e\n",f[0],f[1]);
  // Calculates derivatives at cell
  for(ie=0; ie<NRELBIN; ie++)
  {
    qcool[ie] = - (f[ie+1] - f[ie]) * logbinspace_inv * relel_gammas_inv[ie];
    if isnan(qcool[ie]) out = -1;
  }
  //printf("q[0] %e %e %e \n",qcool[0],relel_gammas[1],relel_gammas_inv[1]);

#endif
  //getch();
  //if (out==-1) printf("\n \n -1 -1 -1 \n\n\n");
  return out;
}


//*********************************************/
//* Nonthermal Radiation Power
//*********************************************/
// fluid frame total energy loss rate from relel. to radiation 
// type = 0: compute with integral over gammadots
// type = 1: compute using relel_dudtau minus the coloumb & particle contributions to conserve energy
ldouble
calc_relel_G0_fluidframe_from_state(ldouble *pp, void *sss, void *ggg, ldouble relel_dudtau, int type)
{
 
  struct geometry *geom
    = (struct geometry *) ggg;
  struct struct_of_state *state
    = (struct struct_of_state *) sss;

  ldouble G0=0.;

#ifdef RELELECTRONS   
  if (type == 1) 
  {
    // Relel Coulomb Coupling rate
    ldouble CC_relel=calc_relel_CoulombCoupling_from_state(pp, state);

    //ANDREW should dn/dtau be computed this way or passed directly?
    ldouble cool_relel_dq=calc_relel_cool_dq_from_state(pp, state); //Cooling dq into the thermal (positive)
  
    //ANDREW - I don't think the chemical potential term should be here
    //ldouble Te=state->Te;
    //ldouble neth=state->ne;
    //ldouble theta= kB_over_me*Te;
    //ldouble mu = chemical_potential_short(theta, neth); //usually negative
    //ldouble cool_relel_dn=calc_relel_cool_dn_from_state(pp, state, geom); //Cooling dn into the thermal (positive)

    ldouble Qe = CC_relel + cool_relel_dq; // - mu*cool_relel_dn; //ANDREW should the chemical potential term be here? 

    G0 = relel_dudtau + Qe; //Qe should be positive, so should decrease the abs value of the negative relel_dudtau
  }

  else if (type==0) 
  {
    G0 = calc_relel_G0_fluidframe_direct_from_state(pp, state, 0); // radtype=0 returns G0 for all types of radiation
  }

  if (G0 > 0.) G0 = 0.; //No absorption, so G0_rel should always be negative.

#ifdef NORELELCOOLATBH
  ldouble xxBL[4];
  #ifdef PRECOMPUTE_MY2OUT
  get_xxout(geom->ix, geom->iy, geom->iz, xxBL);
  #else
  coco_N(geom->xxvec,xxBL,MYCOORDS,BLCOORDS);
  #endif
  if(xxBL[1]<=rhorizonBL) G0=0.;
#endif     
#endif
  return G0;
}

// fluid frame total energy loss rate from relel. to radiation 
// radtype = 0: all
// radtype = 1: synchrotron only
// radtype = 2: free-free only
// radtype = 3: ic only
// radtype = 4: synchrotron + free-free
ldouble
calc_relel_G0_fluidframe_direct_from_state(ldouble *pp, void *sss, int radtype)
{
 
  struct struct_of_state *state
    = (struct struct_of_state *) sss;

  ldouble G0=0.;
#ifdef RELELECTRONS
  ldouble g1,  gdot1;
  ldouble qG0[NRELBIN];
  int ie;
  
  //all in cgs
  ldouble bsq_cgs=0.; //gauss
  ldouble nion_cgs=0.; //cm^-3
  ldouble neth_cgs=0.; //cm^-3
  ldouble trad_cgs=0.; //K
  ldouble erad_cgs=0.; //erg cm^-3 // in fluid frame!

  // synchrotron
  #ifdef RELEL_SYN
  #ifdef RELEL_SYN_ART
  bsq_cgs = RELEL_B_ART*RELEL_B_ART;
  #else
  bsq_cgs = (state->bsq)*fourmpi*endengu2cgs;
  #endif
  #endif

  // free-free       
  #ifdef RELEL_FF
  #ifdef RELEL_FF_ART
  nion_cgs = RELEL_NION_ART;
  #else
  nion_cgs = (state->ni) * numdensgu2cgs;
  #endif
  #endif

  // inverse compton
  #ifdef RELEL_IC
  #ifdef RELEL_IC_ART
  trad_cgs = RELEL_TRAD_ART;
  erad_cgs = RELEL_ERAD_ART;
  #else 
  trad_cgs = state->Trad;
  erad_cgs = (state->Ehat) * endengu2cgs;
  #endif
  #endif
      
  for (ie=0; ie<NRELBIN; ie++) 
  {
     g1 = relel_gammas[ie];
     gdot1 = 0.;

     if (radtype==0 || radtype==1 || radtype==4) 
     {
       gdot1 += gdot_syn(g1, bsq_cgs);
     } 
     if (radtype==0 || radtype==2 || radtype==4) 
     {
       gdot1 += gdot_ff(g1, nion_cgs);
     }
     if (radtype==0 || radtype==3) 
     {
       gdot1 += gdot_ic(g1, trad_cgs, erad_cgs);
     }

     qG0[ie] = pp[NEREL(ie)]*gdot1;
     
  }
  G0 = -M_ELECTR*integrate_relel_q(qG0);
  
#endif
  return G0;
}


// fluid frame rate of relel photon emission 
// radtype=1; synchrotron only
// radtype=2; free-free only
// radtype=4; both free-free and synchrotron
ldouble
calc_relel_photon_ndot_from_state(ldouble *pp, void *sss, int radtype)
{

  struct struct_of_state *state
    = (struct struct_of_state *) sss;
  
  ldouble ndot=0.;
#ifdef RELELECTRONS
  //all in cgs
  ldouble bsq_cgs=0.; //gauss
  ldouble nion_cgs=0.; //cm^-3

  ldouble g1, gdot1;
  ldouble qndot[NRELBIN];
  int ie;
  
  //Synchrotron
  if (radtype==1 || radtype==4)
  {
    #ifdef RELEL_SYN
    #ifdef RELEL_SYN_ART
    bsq_cgs = RELEL_B_ART*RELEL_B_ART;
    #else
    bsq_cgs = (state->bsq)*fourmpi*endengu2cgs;
    #endif
    #endif

    ldouble ne_cgs=(state->nenth)*numdensgu2cgs; 
    ldouble ndotsyncgs = 1.44e5*sqrt(bsq_cgs)*ne_cgs;
    ndot += (ndotsyncgs*numdenscgs2gu*timegu2cgs);
  }

  //Free-Free
  if (radtype==2 || radtype==4)
  { 
    #ifdef RELEL_FF
    #ifdef RELEL_FF_ART
    nion_cgs = RELEL_NION_ART;
    #else
    nion_cgs = (state->ni) * numdensgu2cgs;
    #endif
    #endif

    for (ie=0; ie<NRELBIN; ie++) 
    {
      gdot1 = gdot_ff(relel_gammas[ie], nion_cgs);
      qndot[ie] = pp[NEREL(ie)]*gdot1*relel_gammas_inv[ie];
    }
    ndot += integrate_relel_q(qndot);
  
  }
#endif
  return ndot;
}


// fluid frame rate of electrons cooling back to thermal distribution (all but adiabatic)
ldouble
calc_relel_cool_dq_from_state(ldouble *pp, void *sss)
{

  struct struct_of_state *state
    = (struct struct_of_state *) sss;
  ldouble qCool=0.;
#ifdef RELELECTRONS  
  ldouble nCool=calc_relel_cool_dn_from_state(pp, state);
  qCool=nCool*M_ELECTR*(RELGAMMAMIN-1.);
#endif  
  return qCool;
}


// fluid frame rate of electrons cooling back to thermal distribution (all but adiabatic)
ldouble
calc_relel_cool_dn_from_state(ldouble *pp, void *sss)
{

  struct struct_of_state *state
    = (struct struct_of_state *) sss;

  ldouble nCool=0.;

#ifdef RELELECTRONS
  ldouble gamma,gdot;


  //all in cgs
  ldouble bsq_cgs=0.; //gauss
  ldouble nion_cgs=0.; //cm^-3
  ldouble neth_cgs=0.; //cm^-3
  ldouble trad_cgs=0.; //K
  ldouble erad_cgs=0.; //erg cm^-3 // in fluid frame!

  // synchrotron cooling rate
  #ifdef RELEL_SYN
  #ifdef RELEL_SYN_ART
  bsq_cgs = RELEL_B_ART*RELEL_B_ART;
  #else
  bsq_cgs = (state->bsq)*fourmpi*endengu2cgs;
  #endif
  #endif

  // free-free cooling      
  #ifdef RELEL_FF
  #ifdef RELEL_FF_ART
  nion_cgs = RELEL_NION_ART;
  #else
  nion_cgs = (state->ni) * numdensgu2cgs;
  #endif
  #endif

  // inverse compton cooling
  #ifdef RELEL_IC
  #ifdef RELEL_IC_ART
  trad_cgs = RELEL_TRAD_ART;
  erad_cgs = RELEL_ERAD_ART;
  #else 
  trad_cgs = state->Trad;
  erad_cgs = (state->Ehat) * endengu2cgs;
  #endif
  #endif
  
  // coloumb cooling
  #ifdef RELEL_COUL
  #ifdef RELEL_COUL_ART
  neth_cgs = RELEL_NETH_ART;
  #else 
  neth_cgs = (state->ne)*numdensgu2cgs;
  #endif
  #endif
    
  // Calculate escape energy rate from lowest bin
  gamma = RELGAMMAMIN;
  gdot = 0.;
  if (gamma > 1.) 
  {
    #ifdef RELEL_SYN
    gdot += gdot_syn(gamma, bsq_cgs); 
    #endif
    #ifdef RELEL_FF
    gdot += gdot_ff(gamma, nion_cgs);
    #endif
    #ifdef RELEL_IC
    gdot += gdot_ic(gamma, trad_cgs, erad_cgs);
    #endif
    #ifdef RELEL_COUL
    gdot += gdot_coul(gamma, neth_cgs); 
    #endif
  }
  
  nCool = gdot * pp[NEREL(0)];  
#endif 

  return nCool;
}

// relelectron Couloub energy exchange with thermal electrons
ldouble calc_relel_CoulombCoupling_from_state(ldouble *pp, void *sss)
{

  struct struct_of_state *state
    = (struct struct_of_state *) sss;
  
  ldouble CoulombC=0.;

#ifdef RELELECTRONS
#ifdef RELEL_COUL
  
  ldouble gdot1;
  ldouble qC[NRELBIN];
  ldouble neth_cgs=0.; //cm^-3
  int ie;

  #ifdef RELEL_COUL
  #ifdef RELEL_COUL_ART
  neth_cgs = RELEL_NETH_ART;
  #else 
  neth_cgs = (state->ne)*numdensgu2cgs;
  #endif
  #endif
  
  for (ie=0; ie<NRELBIN; ie++) 
  {
    gdot1 = gdot_coul(relel_gammas[ie], neth_cgs);
    qC[ie] = pp[NEREL(ie)]*gdot1;
  }
  
  CoulombC = M_ELECTR*integrate_relel_q(qC);

#endif
#endif  
  return CoulombC;
}

// fluid frame total energy loss rate from relel. to radiation 
// type = 0: compute with integral over gammadots
// type = 1: compute using relel_dudtau minus the coloumb & particle contributions to conserve energy
ldouble
calc_relel_G0_fluidframe(ldouble *pp, void *ggg, ldouble relel_dudtau, int type)
{
 
  struct geometry *geom
    = (struct geometry*) ggg;
  struct struct_of_state state;
  fill_struct_of_state(pp,geom,&state);

  ldouble out = calc_relel_G0_fluidframe_from_state(pp,&state,geom,relel_dudtau,type);

  return out;
}

// fluid frame total energy loss rate from relel. to radiation 
// radtype = 0: all
// radtype = 1: synchrotron only
// radtype = 2: free-free only
// radtype = 3: ic only
// radtype = 4: synchrotron + free-free
ldouble
calc_relel_G0_fluidframe_direct(ldouble *pp, void *ggg, int radtype)
{
 
  struct geometry *geom
    = (struct geometry*) ggg;
  struct struct_of_state state;
  fill_struct_of_state(pp,geom,&state);

  ldouble out = calc_relel_G0_fluidframe_direct_from_state(pp,&state,radtype);

  return out;
}

//*********************************************/
//* prints spectra of relativistic electrons at given coordinate
//*********************************************/

int
fprint_relel_spectrum(ldouble t, int ix, int iy, int iz, int nfile, char* folder, char* prefix, int doingavg)
{
  char bufor[50],bufor2[50];
  sprintf(bufor,"%s/%s_%d_%d_%d_%04d.dat",folder,prefix,ix,iy,iz,nfile);

#ifdef RELELECTRONS
  FILE *fout_spectrum=fopen(bufor,"w");
  int ie;
  ldouble n=0.;
  for (ie=0 ; ie<NRELBIN; ie++)
    {
      if(doingavg==0) n = numdensGU2CGS(get_u(p,NEREL(ie),ix,iy,iz));
      else if(doingavg==1) n = numdensGU2CGS(get_uavg(pavg,NEREL(ie),ix,iy,iz));
      fprintf(fout_spectrum,"%e %e\n",relel_gammas[ie], n);
    }
  fclose(fout_spectrum);
#endif
  
  return 0;
}


///////////////////////////////////////////////////////////////
int
fprint_relel_avg_spectrum(ldouble t, int jx, int jy, int jz, int nfile, char* folder, char* prefix,int doingavg)
{
  char bufor[50],bufor2[50];
  sprintf(bufor,"%s/%s_%d_%d_%d_%04d.dat",folder,prefix,jx,jy,jz,nfile);
  
  ldouble pp[NV];

  int ix,iy,iz,iv,ie;
  int bix1,bix2,biy2,biy1;

  //Limits of the average
  //ANDREW SET MANUALLY
  bix1 = jx - 5;
  bix2 = jx + 5; 
  biy1 = jy - 5; 
  biy2 = jy + 5; 
  
  if(bix1<3) bix1=3;
  if(biy1<3) biy1=3;
  if(bix2>(NX-3)) bix2=NX-3;
  if(biy2>(NY-3)) biy2=NY-3;
 
  ldouble volume=0.0;
  ldouble spe[NRELBIN];
  for(ie=0 ; ie<NRELBIN; ie++) spe[ie]=0.0;
   
#ifdef RELELECTRONS
  for(ix=bix1;ix<bix2;ix++)
    for(iy=biy1;iy<biy2;iy++)
      for(iz=0;iz<NZ;iz++)
      {
          struct geometry geom;
          fill_geometry(ix,iy,iz,&geom);
	
          struct geometry geomBL;
          fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	
	  //coordinate
	  ldouble dx[3];
	  get_cellsize_out(ix,iy,iz,dx);
	  //get_cellsize_arb(ix,iy,iz,dx,OUTCOORDS); 

          ldouble dvol=dx[0]*dx[1]*dx[2]*geomBL.gdet;
	  volume+=dvol;

	  //primitives at the cell - either averaged or original, in BL or MYCOORDS
	  for(iv=0 ; iv<NV ; iv++)
	    if(doingavg==0) pp[iv]=get_u(p,iv,ix,iy,iz);
	    else if(doingavg==1) pp[iv]=get_uavg(pavg,iv,ix,iy,iz);
	    else pp[iv]=0.;

          //add spectrum contribution from cell
          for(ie=0 ; ie<NRELBIN ; ie++)
            spe[ie] += pp[NEREL(ie)]*dvol;
      } 

  // Save spectrum file
  FILE *fout_spectrum=fopen(bufor,"w");
  ldouble n;
  for (ie=0 ; ie<NRELBIN ; ie++)
    {
      n = numdensGU2CGS((spe[ie]/volume));
      fprintf(fout_spectrum, "%e %e\n", relel_gammas[ie], n);
    }
  fclose(fout_spectrum);
#endif
  
  return 0;
}


