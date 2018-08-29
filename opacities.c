/*! \file opacities.c
 \brief Opacity-related functions
*/

#include "ko.h"
#include <stdio.h>

// Define global quantities for CHIANTI and Sutherland_Dopita opacities
int CHIANTI_LAMBDA_READ = 0;
int SD_LAMBDA_READ = 0;
int N_TEMPERATURE;

ldouble temperature_min_log, temperature_max_log, delta_temperaturelog, delta_temperaturelog_inv;
ldouble *temperaturelog;
ldouble *Lambdalog;

//**********************************************************************
//******* opacities ****************************************************
//**********************************************************************

//Calculate opacities
ldouble
calc_kappa(ldouble *pp,void *ggg,void *op)
{

  struct geometry *geom
    = (struct geometry *) ggg;
  struct struct_of_state state;
  fill_struct_of_state(pp,geom,&state);

  ldouble out = calc_kappa_from_state(pp,&state,geom,op);

  return out;
}

//Calculate opacities from state structure
ldouble
calc_kappa_from_state(ldouble *pp, void *sss, void *ggg, void *op)
{
  
  struct geometry *geom
  = (struct geometry *) ggg;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;
  
  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  struct opacities *opac
  = (struct opacities *) op;
  
  //default, to mark not assigned opacities
  ldouble kappa=-1.;

  //calculate total emissivity direclty
  ldouble Trad = state->Trad;
  ldouble Tgas = state->Tgas;
  ldouble Te = state->Te;
  ldouble B = sigma_rad_over_pi * Te * Te * Te * Te;
  ldouble rho=state->rho;
  
  //Calculate Trad from full ncompt
#ifdef EVOLVEPHOTONNUMBER
  ldouble ThatradBB = state->TradBB;
  ldouble Thatrad = state->Trad;
#endif //EVOLVEPHOTONNUMBER
  
  //Call to the main opacity code goes here  
  //*******************************************************
  #ifdef PR_KAPPA
  #include PR_KAPPA //defined in PROBLEMS/XXX/kappa.c
  #else //default
  kappa = calc_opacities_from_state(pp,  state, geom, opac);  
  #endif
  //*******************************************************
  
  
  if(opac->kappaGasAbs>=0.)
    opac->totEmissivity=(opac->kappaGasAbs * fourpi * B);
  else
    opac->totEmissivity=(kappa * fourpi * B);
  
  if(kappa<0. && opac->kappaGasAbs<0.)
  {
    printf("negative kappa at local indices: %d %d %d \n primitives: \n",geom->ix,geom->iy,geom->iz);
    int iv; PLOOP(iv) printf("%d: %e \n",iv,pp[iv]);
    printf("Trad: %e\n",Trad);
    printf("Te: %e\n",Te);
  }
  if(opac->kappaGasAbs<0.) //no distinction between opacities
  {
    opac->kappaGasNum=opac->kappaRadNum=kappa;
    opac->kappaGasAbs=opac->kappaRadAbs=kappa;
  }
  if(opac->kappaGasRoss<0.) //no Rosseland singled out
  {
    opac->kappaGasRoss=opac->kappaRadRoss=opac->kappaGasAbs;
  }
  
  return kappa;
}


//scattering opacity using precomputed temperatures
//ANDREW should the effective Te for KN cross section change with rel. electrons?
ldouble
calc_kappaes_with_temperatures(ldouble rho, ldouble Tgas, ldouble Te, ldouble Ti, ldouble Trad)
{

  //Call to the scattering opacity code goes here
  //********************************************
  #ifdef PR_KAPPAES
  #include PR_KAPPAES //defined in PROBLEMS/XXX/kappaes.c
  #else //default
  return 0.;
  #endif

  //********************************************
}


//scattering opacity
ldouble
calc_kappaes(ldouble *pp,void *ggg)
{  
  struct geometry *geom
    = (struct geometry *) ggg;
  
  ldouble rho=pp[RHO];
  ldouble Ti,Te;
  ldouble Tgas=calc_PEQ_Teifrompp(pp,&Te,&Ti,geom->ix,geom->iy,geom->iz);
  ldouble Trad=Te;
  
  //********************************************
  //Call to the scattering opacity code goes here
  #ifdef PR_KAPPAES
  #include PR_KAPPAES //defined in PROBLEMS/XXX/kappaes.c
  #else //default
  return 0.;
  #endif
  //********************************************

}

//Calculates total opacity per dx[]
ldouble
calc_chi(ldouble *pp,void *ggg)
{
  struct opacities opac;
  ldouble kappa=calc_kappa(pp,ggg,&opac);
  ldouble chi=kappa+calc_kappaes(pp,ggg);

  return chi;
}

//Calculates total opacity over dx[]
int
calc_tautot(ldouble *pp, void *ggg, ldouble *dx, ldouble *tautot)
{
  ldouble chi=calc_chi(pp,ggg);

  tautot[0]=chi*dx[0];
  tautot[1]=chi*dx[1];
  tautot[2]=chi*dx[2];

  return 0;
}


//Calculates abs opacity over dx[]
int
calc_tauabs(ldouble *pp, void *ggg, ldouble *dx, ldouble *tauabs)
{
  struct opacities opac;
  ldouble kappa=calc_kappa(pp,ggg,&opac);

  tauabs[0]=kappa*dx[0];
  tauabs[1]=kappa*dx[1];
  tauabs[2]=kappa*dx[2];

  return 0;
}

//Opacities from state, put together by Maciek
ldouble calc_opacities_from_state(ldouble *pp, void *sss, void *ggg, void *op)
{
  ldouble kappa;
  int i, j;
  
  struct geometry *geom
  = (struct geometry *) ggg;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;
  struct opacities *opac
  = (struct opacities *) op;
  
  ldouble rho = state->rho;
  ldouble Te = state->Te;
  ldouble Ti = state->Ti;
  ldouble Tgas = state->Tgas;
  ldouble Trad = state->Trad;
  ldouble TradBB = state->TradBB;
  ldouble rtTe = sqrt(Te);
  ldouble rtTgas = sqrt(Tgas);

  opac->kappaGasAbs=0.; //Planck emission at T_r --> T_e 
  opac->kappaRadAbs=0.; //Planck absorption at T_r and T_e

  opac->kappaGasRoss=0.; //Rosseland at T_r --> T_e 
  opac->kappaRadRoss=0.; //Rosseland at T_r and Te

  opac->kappaGasNum=0.; //Number opacity at T_r --> T_e
  opac->kappaRadNum=0.; //Number opacity at T_r and Te
  
  ldouble k_boltCGS = k_boltz_cgs;
  ldouble h_planckCGS = h_cgs;
  ldouble sigmaCGS = sigma_rad_cgs;
  ldouble mpcgs = m_proton_cgs;
  ldouble rhocgs = rhogu2cgs * rho;
  ldouble nethcgs = state->ne * numdensgu2cgs;
  ldouble Bmagcgs = 0.0;

  if(nethcgs<0) printf("negative neth %e\n",nethcgs);

  ldouble zeta = Trad/Te;
  ldouble zeta_inv = 1. / zeta, zeta_inv_3 = zeta_inv * zeta_inv * zeta_inv;
  ldouble zetaRoot5 = pow(zeta,0.2);
  ldouble zetaRoot5_inv_4 = 1. / (zetaRoot5 * zetaRoot5 * zetaRoot5 * zetaRoot5);
  ldouble zetaRoot5_inv_3 = zetaRoot5 * zetaRoot5_inv_4;
  ldouble scalePlaRosFF;
  
#ifdef MAGNFIELD
  ldouble ucon[4], ucov[4], bcon[4],bcov[4], bsq, bsqcgs;
  
  for (i = 0; i < 4; i++)
  {
    ucon[i] = state->ucon[i];
    ucov[i] = state->ucov[i];
    bcon[i] = state->bcon[i];
    bcov[i] = state->bcov[i];
  }

  bsq = state->bsq;
  bsqcgs = fourmpi *endenGU2CGS(bsq);
  Bmagcgs=sqrt(bsqcgs);
#endif
  
#ifndef SKIPFANCYOPACITIES
  ldouble BBenergy = 4.*sigmaCGS*Te*Te*Te*Te;
  ldouble kappagasff = 0.;
  ldouble kapparadff = 0.;
  ldouble kappagasffross = 0.;
  ldouble kapparadffross = 0.;
  ldouble kappagasbe = 0.;
  ldouble kapparadbe = 0.;
  ldouble kappagassyn = 0.;
  ldouble kapparadsyn = 0.;
  ldouble kapparadnumsyn = 0.;
  ldouble kappagasbf = 0.;
  ldouble kapparadbf = 0.;
  ldouble kappamh = 0.;
  ldouble kapparadsynross = 0.;
  ldouble kappagassynross = 0.;
  ldouble kappagasdc = 0.;
  ldouble kapparaddc = 0.;
  ldouble kapparadnumdc = 0.;
  ldouble kappagasnumdc = 0.;
  
#if defined(SUTHERLAND_DOPITA_LAMBDA)
  
  // When we use SUTHERLAND_DOPITA_LAMBDA opacities
  // We want to switch off CHIANTI, BREMSSTRAHLUNG, BOUNDELECTRON, BOUNDFREE
  #undef CHIANTI_LAMBDA
  #undef BREMSSTRAHLUNG
  #undef BOUNDELECTRON
  #undef BOUNDFREE
  
  if (SD_LAMBDA_READ == 0)  // Read in table of Sutherland_Dopita Lambda values
  {
    // The metallicity is log( (Z/X)/(Z/X)_Sun )
#ifdef METALLICITY //BRANDON - [Fe/H] as defined in define file
    ldouble metallicity = METALLICITY;
#else //BRANDON - [Fe/H] calculated from user defined mass abundances (HFRAC, HEFRAC, MFRAC). This is more consistent.
    ldouble metallicity = -4.; // Z=0 metallicity by default

    if (MFRAC > 0.)
    {
      metallicity = log10( (MFRAC/HFRAC)/(MFRAC_SUN/HFRAC_SUN) );
    }
#endif
    
    if (metallicity > 0.5)
    {
      printf("\nERROR!!  metallicity = %e is outside allowed range!!\n\n", metallicity);
      exit(1);
    }
    
    if (metallicity == 0.)  // solar metallicity
    {
      if(PROCID==0) printf("\nUsing Sutherland-Dopita (1993) opacities for metallicity: %e\n\n", metallicity);
      
      FILE *read_SD_data = fopen("./OPACITIES/Sutherland_Dopita_1993/SD_Z0.dat", "r");
      fscanf(read_SD_data,"%d", &N_TEMPERATURE);
      fscanf(read_SD_data,"%lf %lf %lf", &temperature_min_log, &temperature_max_log, &delta_temperaturelog);
      delta_temperaturelog_inv = 1. / delta_temperaturelog;
      
      if ( ( temperaturelog = (ldouble *) malloc(N_TEMPERATURE*sizeof(ldouble)) ) == NULL )
        my_err("malloc err\n");
      if ( ( Lambdalog = (ldouble *) malloc(N_TEMPERATURE*sizeof(ldouble))) == NULL )
        my_err("malloc err\n");
      ldouble ignore;
      
      int itemperature;
      for (itemperature = 0; itemperature < N_TEMPERATURE; itemperature++)
      {
        fscanf(read_SD_data,"%lf %lf %lf", &temperaturelog[itemperature], &Lambdalog[itemperature], &ignore);
      }
      fclose(read_SD_data);
    }
    else if (metallicity < -3.)  // Use Z=0 metallicity
    {
      if(PROCID==0) printf("\nUsing Sutherland-Dopita (1993) opacities for zero metallicity\n\n");
      
      FILE *read_SD_data = fopen("./OPACITIES/Sutherland_Dopita_1993/SD_Zzero.dat", "r");
      fscanf(read_SD_data,"%d", &N_TEMPERATURE);
      fscanf(read_SD_data,"%lf %lf %lf", &temperature_min_log, &temperature_max_log, &delta_temperaturelog);
      delta_temperaturelog_inv = 1. / delta_temperaturelog;
      
      if ( ( temperaturelog = (ldouble *) malloc(N_TEMPERATURE*sizeof(ldouble)) ) == NULL )
        my_err("malloc err\n");
      if ( ( Lambdalog = (ldouble *) malloc(N_TEMPERATURE*sizeof(ldouble))) == NULL )
        my_err("malloc err\n");
      ldouble ignore;
      
      int itemperature;
      for (itemperature = 0; itemperature < N_TEMPERATURE; itemperature++)
      {
        fscanf(read_SD_data,"%lf %lf %lf", &temperaturelog[itemperature], &Lambdalog[itemperature], &ignore);
      }
      fclose(read_SD_data);
    }
    else  // Interpolate between nearest two metallicities
    {
      if(PROCID==0) printf("\nUsing Sutherland-Dopita (1993) opacities for metallicity: %e\n\n", metallicity);
      
      char *SD_data1;
      char *SD_data2;
      ldouble met1, met2;
      
      if (metallicity > 0. && metallicity <= 0.5)
      {
        SD_data1 = "./OPACITIES/Sutherland_Dopita_1993/SD_Z0.dat";
        SD_data2 = "./OPACITIES/Sutherland_Dopita_1993/SD_Z.5.dat";
        met1 = 0.;
        met2 = 0.5;
      }
      else if (metallicity > -0.5 && metallicity <= 0.)
      {
        SD_data1 = "./OPACITIES/Sutherland_Dopita_1993/SD_Zm.5.dat";
        SD_data2 = "./OPACITIES/Sutherland_Dopita_1993/SD_Z0.dat";
        met1 = -0.5;
        met2 = 0.;
      }
      else if (metallicity > -1. && metallicity <= -0.5)
      {
        SD_data1 = "./OPACITIES/Sutherland_Dopita_1993/SD_Zm1.dat";
        SD_data2 = "./OPACITIES/Sutherland_Dopita_1993/SD_Zm.5.dat";
        met1 = -1.;
        met2 = -0.5;
      }
      else if (metallicity > -1.5 && metallicity <= -1.)
      {
        SD_data1 = "./OPACITIES/Sutherland_Dopita_1993/SD_Zm1.5.dat";
        SD_data2 = "./OPACITIES/Sutherland_Dopita_1993/SD_Zm1.dat";
        met1 = -1.5;
        met2 = -1.;
      }
      else if (metallicity > -2. && metallicity <= -1.5)
      {
        SD_data1 = "./OPACITIES/Sutherland_Dopita_1993/SD_Zm2.dat";
        SD_data2 = "./OPACITIES/Sutherland_Dopita_1993/SD_Zm1.5.dat";
        met1 = -2.;
        met2 = -1.5;
      }
      else if (metallicity >= -3. && metallicity <= -2.)
      {
        SD_data1 = "./OPACITIES/Sutherland_Dopita_1993/SD_Zm3.dat";
        SD_data2 = "./OPACITIES/Sutherland_Dopita_1993/SD_Zm2.dat";
        met1 = -3.;
        met2 = -2.;
      }
      else  // Error in metallicity
      {
        printf("\nERROR!!  metallicity = %e is outside allowed range!!\n\n", metallicity);
        exit(1);
      }
      
      FILE *read_SD_data1 = fopen(SD_data1, "r");
      FILE *read_SD_data2 = fopen(SD_data2, "r");
      fscanf(read_SD_data1,"%d", &N_TEMPERATURE);
      fscanf(read_SD_data1,"%lf %lf %lf", &temperature_min_log, &temperature_max_log, &delta_temperaturelog);
      fscanf(read_SD_data2,"%d", &N_TEMPERATURE);
      fscanf(read_SD_data2,"%lf %lf %lf", &temperature_min_log, &temperature_max_log, &delta_temperaturelog);
      delta_temperaturelog_inv = 1. / delta_temperaturelog;
      
      if ( ( temperaturelog = (ldouble *) malloc(N_TEMPERATURE*sizeof(ldouble)) ) == NULL )
        my_err("malloc err\n");
      if ( ( Lambdalog = (ldouble *) malloc(N_TEMPERATURE*sizeof(ldouble))) == NULL )
        my_err("malloc err\n");
      ldouble ignore;
      
      int itemperature;
      ldouble wt1 = (met2 - metallicity) / (met2 - met1);
      ldouble wt2 = (metallicity - met1) / (met2 - met1);
      ldouble Lambdalog1, Lambdalog2;
      for (itemperature = 0; itemperature < N_TEMPERATURE; itemperature++)
      {
        fscanf(read_SD_data1,"%lf %lf %lf", &temperaturelog[itemperature], &Lambdalog1, &ignore);
        fscanf(read_SD_data2,"%lf %lf %lf", &ignore, &Lambdalog2, &ignore);
        Lambdalog[itemperature] = wt1 * Lambdalog1 + wt2 * Lambdalog2;
      }
      fclose(read_SD_data1);
      fclose(read_SD_data2);
      
    }
    
    SD_LAMBDA_READ = 1;  // Read the data only once
    
  }  // end of SD_LAMBDA_READ
  
  // Compute SD lambda for local temperature
  ldouble logTe_local = log10(Te);
  int ilow = (logTe_local - temperature_min_log) * delta_temperaturelog_inv;
  if (ilow < 0) ilow = 0;
  if (ilow > N_TEMPERATURE - 2) ilow = N_TEMPERATURE - 2;
  int ihigh = ilow + 1;
  
  ldouble Lambdalog_local;
  Lambdalog_local =  Lambdalog[ilow] * (temperaturelog[ihigh] - logTe_local);
  Lambdalog_local += Lambdalog[ihigh] * (logTe_local - temperaturelog[ilow]);
  Lambdalog_local *= delta_temperaturelog_inv;
  
  // Convert Lambda to kappagas and then convert to gravitational units
  kappagasff = lengu2cgs * ( pow(10, Lambdalog_local) * (nethcgs * nethcgs / BBenergy) );
  
  // Compute other kappas from kappagasff using the same conversions as in BREMSSTRAHLUNG
  scalePlaRosFF = (14.12/(432.7 - 106.8 * zetaRoot5_inv_3 + 43.17 * zetaRoot5_inv_4 + 57.88 * zeta_inv)	);
  kappagasffross = kappagasff * 0.0330;

#ifdef GRAY_BREMSS
  kapparadff = kappagasff;
  kapparadffross = kappagasffross;
#else
  kapparadff = kappagasff * log1p(1.6*zeta) * one_over_log_2p6 * zeta_inv_3;
  kapparadffross = kappagasff * scalePlaRosFF * zeta_inv_3;
#endif
  
#elif defined(CHIANTI_LAMBDA)  // end of SUTHERLAND_DOPITA_LAMBDA, start of CHIANTI_LAMBDA
  
  // When we use CHIANTI opacities, we want to switch off BREMSSTRAHLUNG, BOUNDELECTRON, BOUNDFREE
  #undef SUTHERLAND_DOPITA_LAMBDA
  #undef BREMSSTRAHLUNG
  #undef BOUNDELECTRON
  #undef BOUNDFREE

  if (CHIANTI_LAMBDA_READ == 0)  // Read in table of CHIANTI Lambda values
  {
    FILE *read_CHIANTI_data = fopen("./OPACITIES/CHIANTI/opacity_CHIANTI.dat", "r");
    fscanf(read_CHIANTI_data,"%d", &N_TEMPERATURE);
    fscanf(read_CHIANTI_data,"%lf %lf %lf", &temperature_min_log, &temperature_max_log, &delta_temperaturelog);
    delta_temperaturelog_inv = 1. / delta_temperaturelog;
    
    if ( ( temperaturelog = (ldouble *) malloc(N_TEMPERATURE*sizeof(ldouble)) ) == NULL )
      my_err("malloc err\n");
    if ( ( Lambdalog = (ldouble *) malloc(N_TEMPERATURE*sizeof(ldouble))) == NULL )
      my_err("malloc err\n");
    ldouble ignore;
    
    int itemperature;
    for (itemperature = 0; itemperature < N_TEMPERATURE; itemperature++)
    {
      fscanf(read_CHIANTI_data,"%lf %lf %lf", &temperaturelog[itemperature], &Lambdalog[itemperature], &ignore);
    }
    
    fclose(read_CHIANTI_data);
    CHIANTI_LAMBDA_READ = 1;  // Read data only once
    
    if(PROCID==0) printf("\nUsing CHIANTI opacities for solar metallicity\n\n");
  }  // end of CHIANTI_LAMBDA_READ
  
  // Compute Lambda for local temperature
  ldouble logTe_local = log10(Te);
  int ilow = (logTe_local - temperature_min_log) * delta_temperaturelog_inv;
  if (ilow < 0) ilow = 0;
  if (ilow > N_TEMPERATURE - 2) ilow = N_TEMPERATURE - 2;
  int ihigh = ilow + 1;

  ldouble Lambdalog_local;
  Lambdalog_local =  Lambdalog[ilow] * (temperaturelog[ihigh] - logTe_local);
  Lambdalog_local += Lambdalog[ihigh] * (logTe_local - temperaturelog[ilow]);
  Lambdalog_local *= delta_temperaturelog_inv;

  // Convert Lambda to kappagas and then convert to gravitational units
  kappagasff = lengu2cgs * ( pow(10, Lambdalog_local) * (nethcgs * nethcgs / BBenergy) );

  // Compute other kappas from kappagasff using the same conversions as in BREMSSTRAHLUNG
  scalePlaRosFF = (14.12/(432.7 - 106.8 * zetaRoot5_inv_3 + 43.17 * zetaRoot5_inv_4 + 57.88 * zeta_inv)	);
  kappagasffross = kappagasff * 0.0330;

#ifdef GRAY_BREMSS
  kapparadff = kappagasff;
  kapparadffross = kappagasffross;
#else
  kapparadff = kappagasff * log1p(1.6*zeta) * one_over_log_2p6 * zeta_inv_3;
  kapparadffross = kappagasff * scalePlaRosFF * zeta_inv_3;
#endif

#elif defined(BREMSSTRAHLUNG) // end of CHIANTI_LAMBDA, start of BREMSSTRAHLUNG 
  ldouble emisFF;

  //BRANDON - added extra term in n_avg for general metallicity
  // n_avg = (X + Y + <Z_j^2/A_j>*Z)*rho/mp  
  emisFF = 1.4e-27*nethcgs*((HFRAC + HEFRAC + Z2divA_MEAN*MFRAC)/mpcgs)*rtTe*1.2*(1.+4.4e-10*Te); 
  
  //this is standard emisivity divided by rhocgs, with Gaunt factor 1.2 and a relativistic correction
  kappagasff = kappacgs2gu*(emisFF/BBenergy)*rho;

  scalePlaRosFF = (14.12/(432.7 - 106.8 * zetaRoot5_inv_3 + 43.17 * zetaRoot5_inv_4 + 57.88 * zeta_inv)	);
  kappagasffross = kappagasff * 0.0330;

#ifdef GRAY_BREMSS
  kapparadff = kappagasff;
  kapparadffross = kappagasffross;
#else
  kapparadff = kappagasff * log1p(1.6*zeta) * one_over_log_2p6 * zeta_inv_3;
  kapparadffross = kappagasff * scalePlaRosFF * zeta_inv_3;
#endif  

#endif //SUTHERLAND_DOPITA_LAMBDA or CHIANTI_LAMBDA or BREMSSTRAHLUNG
  
#ifdef BOUNDELECTRON
  // taken from McKinney et al. 2015
  kappagasbe = kappacgs2gu * 4.0e34*(1-HFRAC-HEFRAC) * 50. * pow(Te,-4.7)*(nethcgs*mpcgs/rhocgs))*rho;
  kapparadbe = kappagasbe * zeta_inv_3;
#endif //BOUNDELECTRON
  
#ifdef BOUNDFREE
  // this is fine for solar composition HFRAC = 0.70, HEFRAC = 0.28, taken from McKinney et al. 2015
  // ANDREW DID NOT change for nonthermal electrons
  kappagasbf = kappacgs2gu * (3.0e25*(1.+HFRAC + 0.75*HEFRAC)*(1.-HFRAC- HEFRAC)*rhocgs/ (Te * Te * Te * rtTe) * log_2p6) * rho;
  kapparadbf = kappagasbf*log(1. + 1.6*zeta) * one_over_log_2p6 * zeta_inv_3;
#endif //BOUNDFREE
  
#ifdef DOUBLECOMPTON
  ldouble theta_e = k_boltCGS*Te/M_ELECTR_CGS/CCC0/CCC0;
  ldouble theta_rad = k_boltCGS*Trad/M_ELECTR_CGS/CCC0/CCC0;
  ldouble pthe = 1./(1.+ theta_e)/(1.+ theta_e)/(1.+ theta_e);
  if (theta_e > 1) pthe = 0.;

  kapparaddc = kappaCGS2GU(7.36e-46*Trad*Trad)*nethcgs*mpcgs*pthe/(1./6.83 + pow(theta_rad, 3.63)/0.0374 + pow(theta_rad, 3.63/3)/0.134 );
  kappagasdc = kapparaddc*zeta*zeta*zeta*zeta;
  
  kapparadnumdc = kappaCGS2GU(7.36e-46*Trad*Trad)*nethcgs*mpcgs*pthe/(1./116. + pow(theta_rad, 3.03)/1.34 + pow(theta_rad, 3.03/3)/4.72 );
  kappagasnumdc = kapparadnumdc*zeta*zeta*zeta;
#endif //DOUBLECOMPTON

#ifdef SYNCHROTRON
#ifndef SYNCHROTRON_NEW_FIT // use old fit by default

  ldouble  nu_MBsyn, zetaBsyngas, zetaBsynrad, zetaAdenomNum, zetaAdenomrad, IaByBrad, IaByBnum, emisSynchro;
  ldouble Tc_n;

  Tc_n=5.07783e9; //Crossover temperature for number opacity
  
  ldouble Trad_lim = pow(TradBB,1.333333333333)/(pow(Te,0.333333333333));
  ldouble nph_lim = pp[NF]*(Trad/Trad_lim);

  #ifdef USE_SYNCHROTRON_BRIDGE_FUNCTIONS // TradBB should not exceed Te BB
  if(Trad < Trad_lim)
  {
    Trad=Trad_lim;
  }
  #endif

  emisSynchro = 3.61e-34*(nethcgs/rhocgs)*Te*Te*Bmagcgs*Bmagcgs;
  nu_MBsyn = 1.19e-13*Te*Te;
  zetaBsynrad = k_boltCGS*Trad/h_planckCGS/nu_MBsyn;
  zetaBsyngas =  k_boltCGS*Te/h_planckCGS/nu_MBsyn;
  
  //ANDREW - this causes many failures when I try turning it back on.
  /*
   ldouble zetas=zetaBsynrad/Bmagcgs;
   ldouble zetasmin=1.;
   ldouble zetasmax=1.e40;
   if(zetas<zetasmin)
   zetaBsynrad*=zetasmin/zetas;
   if(zetas>zetasmax)
   zetaBsynrad/=zetas/zetasmax;
  */
  
  zetaAdenomrad =
  1.79*cbrt(zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*Bmagcgs*Bmagcgs*Bmagcgs*Bmagcgs)
  + 1.35*cbrt(zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*Bmagcgs*Bmagcgs)
  + 0.248*zetaBsynrad*zetaBsynrad*zetaBsynrad;
  IaByBrad = Bmagcgs*Bmagcgs/zetaAdenomrad;
  
  kapparadsyn = kappacgs2gu * (2.13e39*(nethcgs/rhocgs)/ (Te * Te * Te * Te * Te) *IaByBrad)*rho;  
  kappagassyn = kappacgs2gu * (emisSynchro/BBenergy) * rho;

  //converge to Non-relativistic opacities at low Te
  //Adding NR component, becomes dominant for Te < T_crossover=10^9 K roughly
  #ifdef USE_SYNCHROTRON_BRIDGE_FUNCTIONS 
  kapparadsyn += kappacgs2gu * (2.35869e-21*(nethcgs/rhocgs)*(Bmagcgs*Bmagcgs)/(Trad*Trad*Trad))*rho; 
  kappagassyn += kappacgs2gu * (2.35869e-21*(nethcgs/rhocgs)*(Bmagcgs*Bmagcgs)/(Te*Te*Te))*rho; //Adding NR component
  #endif
  
  //number-of-photons averaged opacities
  zetaAdenomNum =
  (0.025*cbrt(zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*Bmagcgs*Bmagcgs)
   + 0.169*cbrt(zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*Bmagcgs)
   + 0.287*zetaBsynrad*zetaBsynrad);
  IaByBnum = (Bmagcgs/zetaAdenomNum);

  kapparadnumsyn = kappacgs2gu * (2.13e39 * (nethcgs/rhocgs)/ (Te * Te * Te * Te * Te) *IaByBnum)*rho;

  #ifdef USE_SYNCHROTRON_BRIDGE_FUNCTIONS 
  kapparadnumsyn *= (Te/Tc_n)/(1.+(Te/Tc_n));
  #endif
  
  //Rosseland
  ldouble IaByBrossRad, IaByBrossGas, Bmb;
  
  //Bmb is pow(Bmagcgs,-0.463), but we need to account for the case of B = 0 here
  if (Bmagcgs > 0. )
    Bmb = pow(Bmagcgs,-0.463);
  else
    Bmb = 0.; // regardless of Bmb, if it is not NaN, we will get kappaRoss
  
  IaByBrossRad = 0.13*pow(Bmagcgs,0.69)/ (pow(zetaBsynrad,1.69) * exp(1.6*pow(zetaBsynrad,0.463)*Bmb)); 
  
  kapparadsynross = kappacgs2gu * (2.13e39*(nethcgs/rhocgs)/ (Te * Te * Te * Te * Te) *IaByBrossRad)*rho; 
  
  IaByBrossGas = 0.13*pow(Bmagcgs,0.69)/ (pow(zetaBsyngas,1.69) * exp(1.6*pow(zetaBsyngas,0.463)*Bmb));
  
  kappagassynross = kappacgs2gu * (2.13e39*(nethcgs/rhocgs)/ (Te * Te * Te * Te * Te) *IaByBrossGas)*rho; 

  //Ramesh: suppress synchrotron opacity at nonrelativistic temperatures -- avoids numerical problems at low temperatures
  #ifndef USE_SYNCHROTRON_BRIDGE_FUNCTIONS
  ldouble Terel = Te * k_over_mecsq;  // this is kT/mec^2
  ldouble Terelfactor = (Terel * Terel) / (1. + Terel * Terel);  // suppression factor
  kapparadsyn *= Terelfactor;
  kappagassyn *= Terelfactor;
  kapparadnumsyn *= Terelfactor;
  kapparadsynross *= Terelfactor;
  kappagassynross *= Terelfactor;
  #endif
  
#else //SYNCHROTRON_NEW_FIT - this is a failure!!

ldouble  Te_inv, B_033, B_inv, zetaSy, zetaSy_m033, zetaB, zetaB_m033;
ldouble  emisSynchro, fit_rad_syn, fit_num_syn, fit_Ross_rad_syn, fit_Ross_gas_syn;

Te_inv = 1./Te;
B_033 = cbrt(Bmagcgs);

if (Bmagcgs > 0. )
  B_inv = 1./Bmagcgs;
else
  B_inv = 0.; 

kappagassyn = 1.5916*1e-30*Bmagcgs*Bmagcgs*Te_inv*Te_inv*nethcgs*kappacgs2gu*rhocgs2gu; //kappa_{P,e}^{sy} in GU
zetaSy = 1.75e23*Trad*Te_inv*Te_inv;
zetaSy_m033 = 1./cbrt(zetaSy);
zetaB = zetaSy*B_inv;

//limit for zetaB needed because number opacity diverges for low B
ldouble zetaBmax = 1e5;
if (zetaB > zetaBmax )
  zetaB = zetaBmax; 

zetaB_m033 = zetaSy_m033*B_033;

//zetaB is: k_boltCGS*T_{rad,e}/h_planckCGS/nu_MBsyn with nu_MBsyn = 1.19e-13*Te*Te
//some messing with B_inv is to avoid division by zero

fit_rad_syn = 1.005*zeta_inv_3/(1. + 5.444*zetaB_m033*zetaB_m033 + 7.218*zetaB_m033*zetaB_m033*zetaB_m033*zetaB_m033);
kapparadsyn = kappagassyn*fit_rad_syn; //kappa_{P,a}^{sy}

fit_num_syn = zeta_inv_3*0.868*zetaB/(1. + 0.589*zetaB_m033 + 0.087*zetaB_m033*zetaB_m033);
kapparadnumsyn = kappagassyn*fit_num_syn; //kappa_{n}^{sy}

fit_Ross_rad_syn = zeta_inv_3*(3.24*1.e-2)*pow(zetaB,1.31)*exp(-1.6*pow(zetaB,0.463));
kapparadsynross = kappagassyn*fit_Ross_rad_syn; // kappa_{R,a}^{sy}

fit_Ross_gas_syn = zeta_inv_3*(3.24*1.e-2)*pow(zetaB*zeta_inv,1.31)*exp(-1.6*pow(zetaB*zeta_inv,0.463));
kappagassynross = kappagassyn*fit_Ross_gas_syn; // kappa_{R,e}^{sy}

emisSynchro = 3.61e-34*(nethcgs/rhocgs)*Te*Te*Bmagcgs*Bmagcgs; // standard emisivity divided by density

  // Ramesh: suppress synchrotron opacity at nonrelativistic temperatures -- avoids numerical problems at low temperatures 
  ldouble Terel = Te * k_over_mecsq;  // this is kT/mec^2
  ldouble Terelfactor = (Terel * Terel) / (1. + Terel * Terel);  // suppression factor
  kapparadsyn *= Terelfactor;
  kappagassyn *= Terelfactor;
  kapparadnumsyn *= Terelfactor;
  kapparadsynross *= Terelfactor;
  kappagassynross *= Terelfactor;

#endif  //SYNCHROTRON_NEW_FIT
#endif  //SYNCHROTRON

  // sum up all the absorption opacities
  opac->kappaGasAbs=kappagasff+kappagassyn+kappagasbe+kappagasbf;
  opac->kappaRadAbs=kapparadff+kapparadsyn+kapparadbe+kapparadbf;
  opac->kappaGasNum=kappagasff+kappagasbe+kappagasbf; //synchrotron emission number opacity applied separately
  opac->kappaRadNum=kapparadff+kapparadnumsyn+kapparadbe+kapparadbf;
  opac->kappaGasRoss=kappagasffross+kappagassynross+kappagasbe+kappagasbf;
  opac->kappaRadRoss=kapparadffross+kapparadsynross+kapparadbe+kapparadbf;

  // by default return kappa
  kappa=opac->kappaGasAbs;
  if(kappa<0){
    printf("negative kappa %e %e  %e %e\n",kappagasff,kappagassyn,kappagasbe,kappagasbf);
    //getch();
  }

#else //SKIPFANCYOPACITIES
  //the very simplest opacities - free-free only
  
  kappa=kappacgs2gu * ((6.6e-24/(mpcgs*mpcgs))*rhocgs/ (Tgas * Tgas * Tgas * rtTgas))*rho*(1.+4.4e-10*Tgas);
  opac->kappaGasAbs=kappa;
  opac->kappaRadAbs=kappa;
  opac->kappaGasNum=kappa;
  opac->kappaRadNum=kappa;
  opac->kappaGasRoss=kappa;
  opac->kappaRadRoss=kappa;

  if(kappa<0)   
   printf("SKIPFANCYOPACITIES KAPPA: Tgas: %e Trad: %e Kappa: %e \n",Tgas,Trad,kappa);

#endif
  
  return kappa;
}

//**********************************************************************
// CHIANTI opacity table functions *************************************
//**********************************************************************

int init_all_kappa_table()
{

#ifdef USE_CHIANTI_ISM_TABLE
  init_ChiantiISMTable();
#endif

#ifdef USE_PLANCK_TABLE
  init_OpTable(&PlanckTable, PLANCK_FILE_NAME);
#endif

#ifdef USE_ROSS_TABLE
  init_OpTable(&RossTable, ROSS_FILE_NAME);
#endif

#ifdef USE_PLANCK_TABLE
  ldouble  chiantilogT0[NCHIANTI] = {
    4.0 , 4.05 , 4.1 , 4.15 , 4.2 , 4.25 , 4.3 , 4.35 , 4.4 , 
    4.45 , 4.5 , 4.55 , 4.6 , 4.65 , 4.7 , 4.75 , 4.8 , 4.85 ,
    4.9 , 4.95 , 5.0 , 5.05 , 5.1 , 5.15 , 5.2 , 5.25 , 5.3 , 
    5.35 , 5.4 , 5.45 , 5.5 , 5.55 , 5.6 , 5.65 , 5.7 , 5.75 , 
    5.8 , 5.85 , 5.9 , 5.95 , 6.0 , 6.05 , 6.1 , 6.15 , 6.2 , 
    6.25 , 6.3 , 6.35 , 6.4 , 6.45 , 6.5 , 6.55 , 6.6 , 6.65 , 
    6.7 , 6.75 , 6.8 , 6.85 , 6.9 , 6.95 , 7.0 , 7.05 , 7.1 , 
    7.15 , 7.2 , 7.25 , 7.3 , 7.35 , 7.4 , 7.45 , 7.5 , 7.55 , 
    7.6 , 7.65 , 7.7 , 7.75 , 7.8 , 7.85 , 7.9 , 7.95 , 8.0
  };


  // this is log10(kappa/rho) for Chianti
  // to get kappa, take value below + logRho, then do pow10
  ldouble chiantilogkappa0[NCHIANTI] = {
    11.75717 ,   12.08387 ,   12.34371 ,   12.50191 ,   12.49125 ,   
    12.30155 ,   12.01424 ,   11.70515 ,   11.41899 ,   11.17827 ,   
    10.97482 ,   10.79819 ,   10.64821 ,   10.52418 ,   10.42249 ,   
    10.33559 ,   10.25433 ,   10.16331 ,    10.0429 ,   9.892752 ,   
    9.718066 ,   9.510048 ,   9.291102 ,   9.094951 ,   8.917125 ,   
    8.738897 ,   8.549819 ,   8.354955 ,   8.147079 ,   7.889193 ,   
    7.561365 ,   7.220804 ,   6.930424 ,    6.69094 ,   6.473525 ,   
    6.249058 ,   6.007913 ,   5.772385 ,    5.55806 ,   5.355904 ,   
    5.152081 ,   4.941354 ,   4.728707 ,   4.515663 ,   4.294598 ,   
    4.054255 ,   3.784555 ,   3.487911 ,   3.185084 ,   2.890325 ,   
    2.610972 ,   2.353685 ,   2.119692 ,   1.905097 ,   1.704218 ,   
    1.512006 ,   1.324008 ,    1.13596 ,   0.943657 ,   0.743129 ,  
    0.5308462 ,  0.3044171 , 0.06315171 , -0.1874624 , -0.4359885 , 
    -0.6732784 , -0.8971446 ,   -1.10934 ,  -1.312718 ,  -1.509582 ,  
    -1.701588 ,  -1.889981 ,  -2.075686 ,  -2.259378 ,  -2.441602 ,   
    -2.62273 ,    -2.8031 ,  -2.982901 ,   -3.16231 ,  -3.341428 ,  
    -3.520309
  };

  if ( ( chiantilogkappa = (ldouble *) malloc(NCHIANTI*sizeof(ldouble)) ) == NULL )
    my_err("malloc err\n");
  if ( ( chiantilogT = (ldouble *) malloc(NCHIANTI*sizeof(ldouble))) == NULL )
    my_err("malloc err\n");


  int i;

  for (i = 0; i<NCHIANTI; i++) {
    chiantilogkappa[i] = chiantilogkappa0[i];
    chiantilogT[i] = chiantilogT0[i];
  }
  
  //replace the first row of the Planck table
  //to allow smooth transition between Chianti and Planck
  for (i = 0; i<PlanckTable.NLOGT; i++) {
    PlanckTable.table[i][0] = log( return_Chianti( PlanckTable.logTgrid[i], PlanckTable.logRhogrid[0] ) )*ONEOVERLOGTEN;
  }
  
#endif  

  return 0;
}

int init_OpTable(void *optab0, char *filename)
{
  FILE *OpFile;
  int i,j;
  struct OpTable *optab = (struct OpTable *) optab0;


  OpFile = fopen(filename, "r");
  if (OpFile == NULL) {
    my_err("Error Reading opacity table!");
  }
  
  
  fscanf(OpFile, "%d %d", &(optab->NLOGT), &(optab->NLOGRHO) );

  if ( ( optab->logTgrid = (ldouble *) malloc(optab->NLOGT*sizeof(ldouble)) ) == NULL ) {
    my_err("malloc err\n");
  }
  
  if ( ( optab->logRhogrid = (ldouble *) malloc(optab->NLOGRHO*sizeof(ldouble)) ) == NULL ) {
    my_err("malloc err\n");
  }

  if ( ( optab->table = (ldouble **) malloc(optab->NLOGT*sizeof(ldouble*)) )== NULL) 
    my_err("malloc err\n");
  
  for (i = 0; i<optab->NLOGT; i++) {
    if ( (optab->table[i] = (ldouble *)malloc(optab->NLOGRHO*sizeof(ldouble)) ) == NULL ) 
      my_err("malloc err\n");
  }
    
  for (i = 0; i < optab->NLOGT; i++) {
    fscanf(OpFile, "%lf", &(optab->logTgrid[i]) );
  }  

  for (i = 0; i < optab->NLOGRHO; i++) {
    fscanf(OpFile, "%lf", &(optab->logRhogrid[i]) );
  }  

  for (i=0; i<optab->NLOGT; i++) {
    for (j=0; j<optab->NLOGRHO; j++) {
      fscanf(OpFile, "%lf", &(optab->table[i][j]) );
    }
  }
  
  fclose(OpFile);

  return 0;
}

#ifdef USE_CHIANTI_ISM_TABLE
//reads  PROBLEMS/ULXBUBBLE/cooling_function.txt to memory
int
init_ChiantiISMTable(void)
{

  ChiantiISMTable = (ldouble **)malloc(sizeof(ldouble*));
  ChiantiISMTable[0]= (ldouble *)malloc(2*sizeof(ldouble));

  FILE *fin = fopen("PROBLEMS/ULXBUBBLE/cooling_function.txt","r");
  if(fin==NULL)
    my_err("PROBLEMS/ULXBUBBLE/cooling_function.txt file missing.\n");

  ldouble temp,lambda;
  int idx=0;
  while(fscanf(fin,"%lf %lf\n",&temp,&lambda)!=EOF)
    {
      ChiantiISMTable[idx][0]=temp;
      ChiantiISMTable[idx][1]=lambda;
      idx++;
      ChiantiISMTable = (ldouble **)realloc(ChiantiISMTable,(idx+1)*sizeof(ldouble*));
      ChiantiISMTable[idx]= (ldouble *)malloc(2*sizeof(ldouble));
    }
  ChiantiISMTableLength = idx;
  return 0;
  
}

ldouble
return_ChiantiISMTableOpacity(ldouble Tgas, ldouble rhocgs)
{

  int i;
  ldouble stef_bol_cgs = 5.6704e-5;
  ldouble lambda_interp;
  ldouble kappa=-1;
  
  //temperatures in ChiantiISMTable[i][0]
  //cooling fucntions in ChiantiISMTable[i][1]
  ldouble ne=rhocgs/MU_E/M_PROTON_CGS;
  ldouble ni=rhocgs/MU_I/M_PROTON_CGS;
  ldouble nh=ni;

  //do linear interpolation of cooling function
  const ldouble MAXCHIANTITEMP=1.e8;

  if (Tgas <= ChiantiISMTable[0][0])
    {
      lambda_interp = ChiantiISMTable[0][1];
      kappa = (lambda_interp*ne*nh/(4.*rhocgs*stef_bol_cgs*Tgas*Tgas*Tgas*Tgas));
    }      

  else if(Tgas >= MAXCHIANTITEMP)
    {
      lambda_interp  = (2.4093687697571752e-27)*sqrt(Tgas)*(1.+4.4e-10*Tgas);
      kappa = (lambda_interp*ne*nh/(4.*rhocgs*stef_bol_cgs*Tgas*Tgas*Tgas*Tgas));
    }

  else

    {

      for(i=1;i<ChiantiISMTableLength;i++)
	{

	if(Tgas >  ChiantiISMTable[i][0])
	  {
	    continue;
	  }
    
	if(Tgas <= ChiantiISMTable[i][0])
	  {
	    lambda_interp = ChiantiISMTable[i-1][1]+((ChiantiISMTable[i][1]-ChiantiISMTable[i-1][1])/(ChiantiISMTable[i][0]-ChiantiISMTable[i-1][0]))*(Tgas-ChiantiISMTable[i-1][0]);
	    kappa = (lambda_interp*ne*nh/(4.*rhocgs*stef_bol_cgs*pow(Tgas, 4.)));
	    break;
	  }
  
      }
    }

  if(kappa<0.)
    my_err("something went wrong in return_ChiantiISMTableOpacity()\n");
  
  return kappa;
}

#endif  

