/*! \file rad.c
 \brief Radiation-related functions
*/

#include "ko.h"
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>


/************************************************************************/
/******* implicit radiative source term in lab frame ********************/
/******* solves numerically in 4D ***************************************/
/************************************************************************/
int
implicit_lab_rad_source_term(int ix,int iy, int iz,ldouble dt)
{
  int iv;
  int ret;
  int verbose=1; //set to 2 to print out the whole failed iterations

#ifdef RADIMPSTOPWHENFAIL
  verbose=2;
#endif
  
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEIMPLICITLAB);

  int vb2=0;//set to 2 to print out the whole failed iterations
  //if(ix==0 && iy==0) vb2=1;
  //attempt to solve the implicit
  ret = solve_implicit_lab(ix,iy,iz,dt,vb2);
    
  if(ret<0) //numerical implicit in 4D did not work
  {
      set_cflag(RADIMPFIXUPFLAG,ix,iy,iz,-1);
      if(verbose) 
      {
          #ifdef REPORTEVERYRADIMPFAILURE
	  printf("%4d > %4d %4d %4d > RADIMPFAIL > unchanged / fixup \n",PROCID,ix+TOI,iy+TOJ,iz+TOK);
          #endif

	  global_int_slot[GLOBALINTSLOT_NTOTALRADIMPFAILURES]++;

	  if(verbose>1)
	  {
              #ifdef REPORTRADIMPFAILUREMINIX
	      if(ix>REPORTRADIMPFAILUREMINIX)
              #endif
	      {
		//try again with verbose set to 2
		solve_implicit_lab(ix,iy,iz,dt,2);
		exit(0);
	      }
	  }
      }
  }
  else //success
  {
      set_cflag(RADIMPFIXUPFLAG,ix,iy,iz,0); 
  }

  return 0;
}


//**********************************************************************
//* wrapper for implicit ***********************************************
//**********************************************************************
int
solve_implicit_lab(int ix,int iy,int iz,ldouble dt,int verbose)
{

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  set_cflag(RADFIXUPFLAG,ix,iy,iz,0);

  int iv;ldouble pp[NV],uu[NV];
  for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,ix,iy,iz); //conserved after advection and geometry and source terms
    pp[iv]=get_u(p,iv,ix,iy,iz); //some reasonable estimate of primitives 
  }

  //otherwise assumes *u and *p consistent! 
  #if (FORCEUEQPINIMPLICIT==1)  //this is the default in choices.h
  p2u(pp,uu,&geom);
  #endif
  
  int corr[3],fixup[2],params[8],ret;
  ldouble pp0[NV],pp00[NV],uu0[NV],uu00[NV];

  PLOOP(iv) 
  {
    pp0[iv]=pp[iv];
    uu0[iv]=uu[iv];
    pp00[iv]=pp[iv];
    uu00[iv]=uu[iv];
  }
  
  //calculate Ehat to determine if RAD or MHD dominates
  ldouble ugas0[4],Rtt0,Tgas0,Trad0,Ehat;
  calc_ff_Rtt(pp0,&Rtt0,ugas0,&geom);
  Ehat=-Rtt0;
  
  //**** set solver parameters *********
  //threshold for working on RAD vs MHD
  ldouble enratiotreshold = RADIMPLICITTHRESHOLD;  
  
  params[3]=0; //no overshooting check
  params[5]=0; //no photons
  #ifdef EVOLVEPHOTONNUMBER
  params[5]=1; //photons

  #ifdef RADIMPFIXNPH 
  params[5]=0; //no photons
  #endif
  #endif

  params[6]=0; //no electrons
  #ifdef EVOLVEELECTRONS
  params[6]=1; //electrons
  #endif

  params[7]=0; //no relativistic electrons
  #ifdef RELELECTRONS
  params[7]=1; //relativistic electrons
  #endif

  int opdamplevel, opdampmaxlevel;
  #if (OPDAMPINIMPLICIT==1)
  opdampmaxlevel=OPDAMPMAXLEVELS;
  #else
  opdampmaxlevel=0;  // this is the default
  #endif

 //run the implicit solver
 //post explicit state in pp0[], initial guess is in pp[]
 for(opdamplevel=0;opdamplevel<=opdampmaxlevel;opdamplevel++)
 {
    ret=-1; 
    params[4]=opdamplevel; //current damping level of opacities

    int startwith;
    if(Ehat<enratiotreshold*pp0[UU]) startwith=RAD;
    else startwith=MHD;

    //*********** 1st ************
    if(ret!=0)
    {
      PLOOP(iv) { pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
      params[1]=RADIMPLICIT_ENERGYEQ;
      params[2]=RADIMPLICIT_LAB;
      #ifdef RADIMPLICIT_STARTWITHFF
      params[2]=RADIMPLICIT_FF;
      #endif

      params[0]=startwith;
      ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,verbose,params,pp); 

      if(verbose>1) getch();
    }
    
    //********** 1.5st  ***********
    //same quantity, different frame
    #ifndef BASICRADIMPLICIT
    if(ret!=0)
    { 
      PLOOP(iv) { pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
      params[2]=RADIMPLICIT_FF;
      #ifdef RADIMPLICIT_STARTWITHFF
      params[2]=RADIMPLICIT_LAB;
      #endif

      ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,verbose,params,pp);
    }      
    #endif
    
    //*********** 2nd *************
    //switch params MHD <--> RAD
  
    if(ret!=0)
    {
      PLOOP(iv) { pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
      params[1]=RADIMPLICIT_ENERGYEQ;
      params[2]=RADIMPLICIT_LAB;
      #ifdef RADIMPLICIT_STARTWITHFF
      params[2]=RADIMPLICIT_FF;
      #endif

      if(params[0]==RAD) params[0]=MHD;
      else params[0]=RAD;
      ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,verbose,params,pp);
    }
    //************2.5nd ************
    //same quantity, different frame
    #ifndef BASICRADIMPLICIT
    if(ret!=0)
    {
      PLOOP(iv) { pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
      params[2]=RADIMPLICIT_FF;
      #ifdef RADIMPLICIT_STARTWITHFF
      params[2]=RADIMPLICIT_LAB;
      #endif
      ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,verbose,params,pp);
    }  
    #endif
   
    //*********** 3rd *************
    //Use entropy fluid frame equations
    #ifndef BASICRADIMPLICIT
    if(ret!=0)
    {
      PLOOP(iv) { pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
      params[1]=RADIMPLICIT_ENTROPYEQ;
      params[2]=RADIMPLICIT_FF;

      params[0]=startwith;
      ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,verbose,params,pp);
    }

   //*********** 4th ************
   //use entropy fluid frame eqns
   //and switch MHD <--> RAD 
   if(ret!=0)
   {
     PLOOP(iv) { pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
     params[1]=RADIMPLICIT_ENTROPYEQ;
     params[2]=RADIMPLICIT_FF;
     
     if(params[0]==RAD) params[0]=MHD;
     else params[0]=RAD;
     ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,verbose,params,pp);
   }
   #endif //BASICRADIMPLICIT

   if(ret==0) //success!
   {
     set_cflag(RADFIXUPFLAG,ix,iy,iz,0);
     break;
   }
   else if(ret!=0)  
   {
     set_cflag(RADFIXUPFLAG,ix,iy,iz,-1); // flag cell for fixups
     global_int_slot[GLOBALINTSLOT_NTOTALRADIMPFIXUPS]++;      
     continue; //with lower, damped opacities
   }
 } //the end of the opdamp loop

 //report failure
 if(ret<0) 
 {
   //if(ix==0 && iy==0) printf("\n implicit solver failed at %i %i %i\n",ix,iy,iz);
   return -1;
 }

 //else, succeeded!
 if(opdamplevel>0 && opdamplevel<=opdampmaxlevel)
    printf("%4d > %4d %4d %4d > opdamp =  %d helped in implicit\n",PROCID,ix+TOI,iy+TOJ,iz+TOK,opdamplevel);

 //solution given in pp[], get consistent uu
 p2u(pp,uu,&geom);

 //save to memory consistent pp & uu
 PLOOP(iv)
 {
#ifdef SKIPHDEVOLUTION
    if(iv>=NVMHD)
#endif

#ifdef RADIATION
#ifdef SKIPRADEVOLUTION
#ifdef EVOLVEPHOTONNUMBER
    if(iv!=EE && iv!=FX && iv!=FY && iv!=FZ && iv!=NF)
#else
    if(iv!=EE && iv!=FX && iv!=FY && iv!=FZ)
#endif
#endif
#endif

#ifdef SKIPHDBUTENERGY
    if(iv>=NVMHD || iv==UU)
#endif
    {
      set_u(p,iv,ix,iy,iz,pp[iv]);
      set_u(u,iv,ix,iy,iz,uu[iv]);
    }
 }

  //calculate stats on average number of succesful iterations
  int whichprim=params[0];
  int whicheq=params[1];
  int whereeq=params[2];
  if(whichprim==RAD && whicheq==RADIMPLICIT_ENERGYEQ && whereeq==RADIMPLICIT_LAB)
  {
    #pragma omp critical
    global_int_slot[GLOBALINTSLOT_NIMPENERRAD]+=1;
  }
  if(whichprim==RAD && whicheq==RADIMPLICIT_ENERGYEQ && whereeq==RADIMPLICIT_FF)
  {
    #pragma omp critical
    global_int_slot[GLOBALINTSLOT_NIMPENERRADFF]+=1;
  }
  if(whichprim==MHD && whicheq==RADIMPLICIT_ENERGYEQ && whereeq==RADIMPLICIT_LAB)
  {
    #pragma omp critical
    global_int_slot[GLOBALINTSLOT_NIMPENERMHD]+=1; 
  }
  if(whichprim==MHD && whicheq==RADIMPLICIT_ENERGYEQ && whereeq==RADIMPLICIT_FF)
  {
    #pragma omp critical
    global_int_slot[GLOBALINTSLOT_NIMPENERMHDFF]+=1; 
  }
  if(whichprim==RAD && whicheq==RADIMPLICIT_ENTROPYEQ)
  {
    #pragma omp critical
    global_int_slot[GLOBALINTSLOT_NIMPENTRRAD]+=1;
  }
  if(whichprim==MHD && whicheq==RADIMPLICIT_ENTROPYEQ)
  {
    #pragma omp critical
    global_int_slot[GLOBALINTSLOT_NIMPENTRMHD]+=1;
  }
  if(whicheq==RADIMPLICIT_LTEEQ)
  {
    #pragma omp critical
    global_int_slot[GLOBALINTSLOT_NIMPLTE]+=1;
  }
  
  return 0;  
}


//**********************************************************************
//**********************************************************************
//******* 4+D implicit solver working on primitives*********************
//**********************************************************************
//**********************************************************************
int
solve_implicit_lab_4dprim(ldouble *uu00,ldouble *pp00,void *ggg,ldouble dt,int verbose,int *params,ldouble *pp)
{
#ifdef RADIATION
  int i1,i2,i3,iv,i,j,ie, ret_fill, ret_copy;
  int u2pret=0;
  ldouble pp0[NV],ppp[NV],uu[NV],uu0[NV],uup[NV];
  ldouble ucon00[4], ucov00[4], urfcon00[4], urfcov00[4];
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  struct struct_of_state state00;  // corresponds to pp00: post-explicit state
  struct struct_of_state state0;   // corresponds to pp0: initial guess
  struct struct_of_state state;    // corresponds to pp: current guess

  struct geometry *geom
    = (struct geometry *) ggg;

  int ix = geom->ix;
  int iy = geom->iy;
  int iz = geom->iz;

  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;
  gdetu=gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif
  
  // initial conditions
  for(iv=0;iv<NV;iv++)
  {
    uu0[iv]=uu00[iv]; 
    pp0[iv]=pp00[iv]; 
    uu[iv]=uu0[iv]; 
    pp[iv]=pp0[iv];     
  }

  // initial states
  ret_fill = fill_struct_of_state(pp00, geom, &state00);
  ret_copy = copy_state_to_state(&state00, &state0);

  // initial state parameters
  ldouble ugas00[4],Rtt00,Tgas00,Trad00,Te00,Ti00,Trad00BB;
  ldouble gamma=state00.gamma;
  ldouble gammam1=gamma-1.;
  
  Rtt00 = -state00.Ehat;
  
  for(i = 0; i < 4; i++)
  {
    ugas00[i] = state00.ucon[i];
    ucon00[i] = state00.ucon[i];
    ucov00[i] = state00.ucov[i];
    urfcon00[i] = state00.urfcon[i];
    urfcov00[i] = state00.urfcov[i];
  }

  Tgas00 = state00.Tgas;
  Te00 = state00.Te;
  Ti00 = state00.Ti;
  Trad00BB = state00.TradBB;
  Trad00=Trad00BB;
  
  // dump a problematic case to file 
  if(verbose) 
  {
      printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
      printf("\n@@@@@@@@ IMPLICIT @@@@@@@@@@@@");
      printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n");

     
      FILE *out = fopen("imp.problem.dat","w");
      for (i1=0;i1<NV;i1++)
	fprintf(out,"%.20e ",uu0[i1]);
      for (i1=0;i1<NV;i1++)
	fprintf(out,"%.20e ",pp0[i1]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  fprintf(out,"%.20e ",gg[i1][i2]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  fprintf(out,"%.20e ",GG[i1][i2]);
      fprintf(out,"%.20e \n",dt);
      fprintf(out,"%.20e \n",geom->alpha);
      fprintf(out,"%.20e \n",geom->gdet);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<4;i2++)
	  fprintf(out,"%.20e ",geom->tup[i1][i2]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<4;i2++)
	  fprintf(out,"%.20e ",geom->tlo[i1][i2]);
      fprintf(out,"%.20e ",geom->gttpert);
      for (i1=0;i1<4;i1++)
	fprintf(out,"%.20e ",geom->xxvec[i1]);
      fprintf(out,"\n");
      fclose(out);
      printf("dumped problematic case to imp.problem.dat: NV = %d\n", NV);
  }

  // changes to initial conditions & initial states
  #ifdef RADIMPSTARTLOWTEMP
  pp0[UU]=calc_PEQ_ufromTrho(1.e4,pp0[RHO]);
  #endif

  // fix the temperature at the start of implicit
#ifdef RADIMPSTARTFIXTEMP
  pp0[UU]=calc_PEQ_ufromTrho(RADIMPSTARTFIXTEMP,pp0[RHO],ix,iy,iz);
#endif
  
  // solve for initial rad guess with bisection
#ifdef RADIMP_START_WITH_BISECT
  solve_implicit_lab_1dprim(uu00, pp00, &state00, geom, dt,  verbose, pp0, &state0);  
#endif
  
  // solve for initial photon number with bisection
#ifdef EVOLVEPHOTONNUMBER
#ifdef NPH_START_WITH_BISECT
  solve_implicit_nphoton_rad_1dprim(uu00, pp00, &state00, geom, dt, verbose, pp0, &state0);
#endif
#endif

  // fix the temperature to the radiation temperature at the start of implicit
#ifdef RADIMPSTARTTRAD
  pp0[UU] = calc_PEQ_ufromTrho(Trad00,pp0[RHO]);
#endif

  //Now that we have initial guess pp0, get corresponding uu0 and state0
  p2u(pp0,  uu0, geom);
  ret_fill = fill_struct_of_state(pp0, geom, &state0);

  // Calculate opacities and print problems to file
  ldouble TLTE=0.;
  struct opacities opac;
  ldouble kappa,kappaes,chi,xi1,xi2;
  ldouble ucon[4];
  if(verbose) 
  {
      kappa=calc_kappa(pp0,geom,&opac);
      kappaes=calc_kappaes(pp0,geom);
      chi=kappa+kappaes;
      xi1=kappa*dt*(1.+16.*SIGMA_RAD*pow(Tgas00,4.)/pp0[UU]);
      xi2=chi*dt*(1.+(-Rtt00)/(pp0[RHO]+gamma*pp0[UU]));
      
      printf("\nparams: \n  whichprim (RAD 1, MHD 2): %d\n  whicheq (ENERGYEQ 0, ENTROPYEQ 1, LTEEQ 2): %d\n  whichframe (LAB 0, FF 1): %d\n  opdamplevel (0 to OPDAMPMAXLEVELS): %d\n  ifncompt (no EVOLVEPHOTONNUMBER 0, with EVOLVEPHOTONNUMBER 1): %d\n  ifentre (no EVOLVEELECTRONS 0, with EVOLVEELECTRONS 1): %d\n  ifrelel (no RELELECTRONS 0, with RELELECTRONS 1): %d\n\n",params[0],params[1],params[2],params[4],params[5],params[6],params[7]);

      if(params[1]==RADIMPLICIT_ENERGYEQ)
	printf("energy eq.\n");
      if(params[1]==RADIMPLICIT_ENTROPYEQ)
	printf("entropy eq.\n");
      if(params[1]==RADIMPLICIT_LTEEQ)
	printf("LTE eq. - no longer supported\n");
      if(params[2]==RADIMPLICIT_FF)
	printf("fluid\n");
      if(params[2]==RADIMPLICIT_LAB)
	printf("lab\n");

      printf("\n===========\nxi1: %e xi2: %e\n===========\n\n",xi1,xi2);
      ldouble Rtt;
      calc_ff_Rtt(pp0,&Rtt,ucon,geom);
      printf("Ehat: %e\n\n",-Rtt);

      ldouble qsq=0.;
      int i,j;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=pp0[UU+i]*pp0[UU+j]*geom->gg[i][j];
      ldouble gamma2=1.+qsq;
      printf("Lorentz gamma gas: %e\n\n",sqrt(gamma2));

      ldouble ucon[4], ucov[4], bcon[4], bcov[4], bsq;
      calc_ucon_ucov_from_prims(pp0, geom, ucon, ucov);
      calc_bcon_bcov_bsq_from_4vel(pp0, ucon, ucov, geom, bcon, bcov, &bsq);
      printf("bsq: %e\n\n",bsq);

      qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=pp0[EE0+i]*pp0[EE0+j]*geom->gg[i][j];
      gamma2=1.+qsq;
      printf("Lorentz gamma rad: %e\n\n",sqrt(gamma2));

      ldouble Gi00[4],Gihat00[4];
      calc_Gi(pp0, geom,Gi00, 0.0, 1, 0); 
      indices_21(Gi00,Gi00,geom->gg);

      printf("viscous heating: %e\n",0.);
	
      boost2_lab2ff_with_alpha(Gi00, Gihat00, pp0, geom->gg, geom->GG, geom->alpha);
      for(iv=0;iv<4;iv++)
      {
	Gi00[iv]*=dt*gdetu;
	Gihat00[iv]*=dt*gdetu;
      }

      ldouble Trad=calc_LTE_TfromE(-Rtt);
      TLTE=Trad;
      ldouble Tgas=calc_PEQ_Tfromurho(pp0[UU],pp0[RHO],geom->ix,geom->iy,geom->iz);
      
      printf("\n===========\nkappa: %e chi: %e\n===========\n",kappa,kappa+kappaes);
      printf("\nxi1: %e xi2: %e\n",xi1,xi2);
      printf("gas temp: %e\n",Tgas00);      
      printf("ele temp: %e\n",Te00);      
      printf("ion temp: %e\n",Ti00);      
      printf("rad temp: %e\n",Trad00); 
      printf("LTE temp: %e\n\n",TLTE) ;  
      printf("gammagas: %f\n",get_u_scalar(gammagas,geom->ix,geom->iy,geom->iz));
      
      #ifdef EVOLVEELECTRONS
      ldouble ue;
      ldouble ui;
      ldouble rhoe=calc_thermal_ne(pp0)*MU_E*M_PROTON;
      ue=calc_ufromSerho(pp0[ENTRE],rhoe,ELECTRONS,geom->ix,geom->iy,geom->iz);
      ui=calc_ufromSerho(pp0[ENTRI],pp0[RHO],IONS,geom->ix,geom->iy,geom->iz);
      printf("ue: %e ui: %e (%e)\n",ue,ui,(ue+ui)/pp0[UU]);
      #endif
  }
    
  //make the choice of which primitives to evolve
  int whichprim,whicheq,do_mom_over,ifncompt,ifentre,ifrelel;
  if(params[5]==0)
  {
    ifncompt=0;
  }
  else
  {
    ifncompt=1; //evolve photon number
  }

  if(params[6]==0)
  {
    ifentre=0;
  }
  else
  {
    ifentre=1; //evolve thermal electron entropy
  }

  if(params[7]==0)
  {
    ifrelel=0;
  }
  else
  {
    ifrelel=1; //evolve relativistic electrons
  }

  if(verbose) printf("\nifncompt/ifentre/ifrelel : %d/%d/%d\n",ifncompt,ifentre,ifrelel);

  //choose the sub-dominant quantity to evolve
  if(params[0]==-1) 
  {
    if(-Rtt00<1.e-3*pp0[UU]) //rad sub-dominant
      whichprim=RAD;
    else                     //gas sub-dominant
      whichprim=MHD; 
    params[0]=whichprim;  
  }
  else
  {
    if(params[0]==RAD)
      whichprim=RAD;
    if(params[0]==MHD)
      whichprim=MHD;
  }
 
  whicheq=params[1];
  do_mom_over = params[3];

  int sh;
  if(whichprim==MHD) 
    sh=UU;  //solving in hydro primitives
  else if(whichprim==RAD) 
    sh=EE0; //solving in rad primitives
  
  if(verbose && whichprim==MHD) printf("\n\nWorking on MHD\n\n");
  if(verbose && whichprim==RAD) printf("\n\nWorking on RAD\n\n");

  //check if one can compare gas & rad velocities
  if(VELPRIM!=VELPRIMRAD) 
    my_err("implicit solver assumes VELPRIM == VELPRIMRAD\n");
 
  //other parameters
  ldouble err,errbase,errbest;
  ldouble EPS = RADIMPEPS;
  ldouble CONV = RADIMPCONV;  
  if(params[1]==RADIMPLICIT_ENTROPYEQ)
     CONV = RADIMPENTRCONV;  
  ldouble MAXITER = RADIMPMAXITER;
  int corr[2],fixup[2];
  ldouble frdt = 1.0;
  int iter=0;
  int failed=0;
  errbest=BIG;

  if(verbose) 
  {
      printf("=== ix, iy, iz: %d %d %d\n\n",ix,iy,iz);
      printf("Conserveds:\n");
      print_conserved(uu0); //ANDREW this uu0 may not be consistent with pp0 above, haven't applied p2u yet
      printf("Primitives:\n");
      print_primitives(pp0);
      printf("Metric gg:\n");
      print_metric(gg);
      printf("\n===\n Trying imp lab 4d prim with dt : %e \n",dt);
  }

  //number of equations solved
  int np=4,inph=-9,ientre=-9,irelel=-9;
  if(ifncompt)
  {
    inph=4;
    np++;
  }
  if(ifentre)
  {
    np++;
    if(ifncompt)
      ientre=5;
    else
      ientre=4;
  }
  if(ifrelel)
  {
    np+=NRELBIN;
    irelel=ientre+1; //nonthermal electron slots are after thermal electron entropy
  }

  //allocating arrays
  int N=np;
  int ib;
  ldouble *tJ, *tiJ,**J,**iJ;
  ldouble *f1,*f2,*f3,*xxx,*xxxbest;
  f1=(ldouble*)malloc(N*sizeof(ldouble));
  f2=(ldouble*)malloc(N*sizeof(ldouble));
  f3=(ldouble*)malloc(N*sizeof(ldouble));
  xxx=(ldouble*)malloc(N*sizeof(ldouble));
  xxxbest=(ldouble*)malloc(N*sizeof(ldouble));
  J=(ldouble**)malloc(N*sizeof(ldouble*));
  iJ=(ldouble**)malloc(N*sizeof(ldouble*));
  tJ=(ldouble*)malloc(N*N*sizeof(ldouble));
  tiJ=(ldouble*)malloc(N*N*sizeof(ldouble));
  for(ib=0;ib<N;ib++)
  {
    J[ib]=(ldouble*)malloc(N*sizeof(ldouble));
    iJ[ib]=(ldouble*)malloc(N*sizeof(ldouble));
  }
  
  /*************************************************************/
  //Solve the implicit problem
  /*************************************************************/
  // Copy initial guess from pp0 to pp
  for (iv = 0; iv < NV; iv++) 
  {
    pp[iv] = pp0[iv];
  }

  //apply constraints and fill state  
  int fret = implicit_apply_constraints(pp, uu, uu00, geom, whichprim);
  fill_struct_of_state(pp,geom,&state);
  
  //main solver loop
  do 
  {	 
      failed=0;
      iter++;
      
      //************************************************
      //Step 0: copy base state and determine its error
      for(i=0;i<NV;i++)
      {
	ppp[i]=pp[i];
      }	
      for(i=0;i<4;i++)
      {
	xxx[i]=ppp[i+sh];
      }  
      if(ifncompt)
	xxx[inph]=ppp[NF0];
      if(ifentre)
	xxx[ientre]=ppp[ENTRE];
      if(ifrelel)
	for(i=0;i<NRELBIN;i++) xxx[irelel+i]=ppp[NEREL(i)];

      //print initial state
      if(verbose>0)
      {
          int ret=f_implicit_lab_4dprim_with_state(uu0, pp, &state, uu0, pp0, &state0, dt, geom, f1, params, &err);
	  print_state_implicit_lab_4dprim(iter-1,xxx,f1,err,N); 
	  printf("\n kappaCGS: %e ", kappaGU2CGS( calc_kappa(pp,geom,&opac)/pp[RHO]) );
	  printf("rhoCGS: %e ", rhoGU2CGS(pp[RHO]));
	  ldouble Ti0,Te0;
          ldouble Tgas0=state0.Tgas;
	  ldouble Trad0=state0.Trad;
	  printf("T: %e Te: %e ", Tgas0,Te0);
	  printf("Trad: %e ", Trad0);
	  printf("Gi^t %e ", state0.Gi[0]);
	  printf("Gic^t %e \n\n", state0.Gic[0]);
	  if(ret<0) printf("f_lab_4dprim ret: %d\n",ret);
      }

      //error vector at base state
      //uses fret from previous loop iteration of f_implicit_lab_4dprim_with_state in check
      //pp should already be consistent -- i.e. implicit_apply_constraints applied at bottom of loop
      if(fret<-1)
      {
	if(verbose>0) printf("base state\n");
	free_solve_implicit_lab_4dprim(J, iJ, tJ, tiJ, f1, f2, f3, xxx, xxxbest,N);
	return -1;	  
      }
      
      f_implicit_lab_4dprim_with_state(uu, pp, &state, uu00, pp00, &state00, dt, geom, f1, params, &err);
      errbase=err;	  

      //has the error converged? 
      if(err<CONV)
      {
	if(verbose) printf("\n === success (error) ===\n");
	if(verbose==2) getchar();
	break;
      }
      if(err<errbest)
      {
	errbest=err;
	for(j=0;j<np;j++)
	  xxxbest[j]=xxx[j];          
      }

      //*******************************************
      //Step 1: calculate approximate Jacobian    
      ldouble del;
      if (verbose == 2)
      {
        printf("Details of Jacobian calculation:\n============\n");
      }
      
      for(j=0;j<np;j++) //loop over Jacobian columns
      {
	  //one-way derivatives
	  //try both signs
	  ldouble sign=-1.;  //minus avoids u2p_mhd errors when working on radiative
	  for(;;)
	  {
	    if(j==0) //rad or mhd energy
	    {
	      del=sign*EPS*ppp[sh];
	      pp[j+sh]=ppp[j+sh]+del;
	    }
	    else if(ifncompt && j==inph) //photon number
	    {
	      del=sign*EPS*ppp[NF0];
	      pp[NF0]=ppp[NF0]+del;
	    }	   
	    else if(ifentre && j==ientre) //electron entropy
	    {
	      del=sign*EPS*ppp[ENTRE];
	      pp[ENTRE]=ppp[ENTRE]+del;
	    }	    
	    else if(ifrelel && j>=irelel && j<irelel+NRELBIN) //nonthermal electrons
	    {
	      //use opposite initial sign for this case since we often start with nrelel=0
	      ie = j-irelel; 
	      del=-sign*EPS*my_max(1.e-20*ppp[RHO]/MU_GAS/M_PROTON, ppp[NEREL(ie)]);
	      pp[NEREL(ie)]=ppp[NEREL(ie)]+del;
	    }
	    else //rad or mhd velocities
	    {
	      ldouble veleps = EPS*my_sign(ppp[j+sh])*my_max(1.e-8/sqrt(geom->gg[j][j]), fabs(ppp[j+sh]));
	      if(ppp[j+sh]>=0.)
		del=sign*veleps; 
	      else
		del=-sign*veleps; 
	      pp[j+sh]=ppp[j+sh]+del;
	    }	    
	    
            //calculate error functions for the perturbed state
	   
	    fret = implicit_apply_constraints(pp, uu, uu00, geom, whichprim);
	    fill_struct_of_state(pp,geom,&state);
	    
	    if(fret<-1)
	    {
	      if(sign>0.) //already switched signs
	      {	      
		if(verbose) printf("Jac mid-state (%d) both signs lead to trouble\n",j);
		failed=1;
		break;
	      }
	      else //try switching signs
	      {
		if(verbose) printf("Jac mid-state (%d) trying the other one\n",j);
		pp[j+sh]=ppp[j+sh];
		sign*=-1.;
		continue;
	      }
	    }

  	    f_implicit_lab_4dprim_with_state(uu, pp, &state, uu00, pp00, &state00, dt, geom, f2, params, &err); 

	    //Jacobian matrix rows for this f2 vector
	    for(i=0;i<np;i++)
	    {

	      J[i][j] = (f2[i] - f1[i]) / del;
	      
              if (verbose==2)
              {
                printf("parameter: %d  equation: %d  Delta p: %e  Delta f: %e\n", j, i, del, f2[i] - f1[i]);
              }
	    }

	    //recover original values of pp
	    if(j<4) 
	      pp[j+sh]=ppp[j+sh];
	    else if(ifncompt && j==inph)
	      pp[NF0]=ppp[NF0];
	    else if(ifentre && j==ientre) 
	    {
	      pp[ENTRE]=ppp[ENTRE];
              pp[ENTRI]=ppp[ENTRI];
	    }
            else if(ifrelel && j>=irelel && j<irelel+NRELBIN)
	    {
	      ie = j-irelel; 
	      pp[NEREL(ie)]=ppp[NEREL(ie)];
	    }
	    break; 
	  } //end loop to try both signs of perturbation for jacobian

	  if(failed==1)
	    break;
      } //end loop over jacobian columns
   
      if(failed==1) break;

     
#ifdef SCALE_JACOBIAN  //scale jacobian for better inversion 
      ldouble scalevec[N];
      ldouble invscalevec[N];
      ldouble nrelel=1.0;
      ldouble neth=calc_thermal_ne(pp);
      if(ifrelel) nrelel=calc_relel_ne(pp);
      
      //Find the scale factors
      for(j=0;j<np;j++)
      {
	if(j==0) 
          scalevec[j]=pp[j+sh];
        else if(j>0 && j<4)
          scalevec[j]=1.0;
	else if(ifncompt && j==inph)
	  scalevec[j]=pp[NF0];
	else if(ifentre && j==ientre) 
	  scalevec[j]=neth*K_BOLTZ;
        else if(ifrelel && j>=irelel && j<irelel+NRELBIN)
          scalevec[j]=nrelel;
	  
        if(scalevec[j] == 0.0) scalevec[j]=1.0; // ANDREW -- use different small cutoff? 
        invscalevec[j] = 1.0 / scalevec[j];
      }

      //if(verbose) {printf("SCALEVEC\n"); print_Nvector(scalevec,np); print_Nvector(invscalevec,np);}
      
      //Scale the Jacobian columns and rows
      for(i=0;i<N;i++)
      {
	for(j=0;j<N;j++)
	{  
	  if(i!=j) J[i][j] = J[i][j] * (scalevec[j] * invscalevec[i]);
        }
      }
#endif  // SCALE_JACOBIAN

      //*******************************************
      //Step 2: Invert the Jacobian
      int ret,ib;
      #pragma omp critical //for some reason gsl-based inverse does not work on my mac with openmp
      {
	for(i=0;i<N;i++)
	  for(j=0;j<N;j++)
	    tJ[i*N+j]=J[i][j]; 

	ret=inverse_matrix(tJ,tiJ,N);
 

	for(i=0;i<N;i++)
	  for(j=0;j<N;j++)
	    iJ[i][j]=tiJ[i*N+j];
      }

      if(ret<0) //jacobian inversion  failed
      {
	failed=1;
	if(verbose || 1) 
	  printf("Jacobian NxN inversion failed\n");//getchar();
	break;
      }
      if(verbose==2)
      {
        printf("\nJacobian:\n");
        print_NNtensor(J,N);
        printf("\nInverse Jacobian:\n");
        print_NNtensor(iJ,N);
      }

      //*******************************************
      //Step 3: Apply inverse Jacobian corrections
      ldouble xiapp=1.;
      do //loop to damp applied inverse jacobian 
      {
	//inital state vector -- need to refill every stage of opdamp loop
	for(i=0;i<np;i++)
	{
	  if(i<4)
	    xxx[i]=ppp[i+sh];
	  else if(ifncompt && i==inph)
	    xxx[inph]=ppp[NF0];		
	  else if(ifentre && i==ientre)
	    xxx[ientre]=ppp[ENTRE];		
	  else if(ifrelel && i>=irelel && i<irelel+NRELBIN)
	  {
	    ie = i-irelel; 
	    xxx[i]=ppp[NEREL(ie)];
	  }
	}
        
	//updating xxx with computed inverse Jacobian
	for(i=0;i<np;i++)
	{
	  for(j=0;j<np;j++)
	  {
	    //apply the inverse Jacobian scaling 
            #ifdef SCALE_JACOBIAN
            xxx[i] -= xiapp*(iJ[i][j]*scalevec[i]*invscalevec[j])*f1[j]; 
            #else
            xxx[i] -= xiapp*iJ[i][j]*f1[j];
            #endif
	  }
	}
         
	if(verbose)
	{
	  printf("\nsub> trying with xi=%e\n",xiapp);
	  print_Nvector(xxx,np);
	}

	//update primitives in ppp from xxx
	for(i=0;i<np;i++)
	{
	  if(i<4)
	    pp[i+sh]=xxx[i];
	  else if(ifncompt && i==inph)
	    pp[NF0]=xxx[inph];
	  else if(ifentre && i==ientre)
	    pp[ENTRE]=xxx[ientre];
	  else if(ifrelel && i>=irelel && i<irelel+NRELBIN)
          {
	    ie = i-irelel;
	    pp[NEREL(ie)]=xxx[irelel+ie];
	    if (pp[NEREL(ie)] < 0.) 
            {
	      pp[NEREL(ie)] = 0.0; // apply a hard floor on nonthermal here
	    }
	  }
	}

	//apply constraints to pp based on primitives not evolved  
	fret = implicit_apply_constraints(pp, uu, uu00, geom, whichprim);
	fill_struct_of_state(pp,geom,&state);  

	  //check limits on the change per step
	  //ANDREW -- could use state here (will need another stateppp!)
	  int okcheck=0;

	  //minimal energy density that you can go down to in one step
	  ldouble edenmin=ppp[sh]/RADIMPLICITMAXENCHANGEDOWN;
	  
	  //maximal energy density that you can go up to in one step
	  ldouble edenmax=ppp[sh]*RADIMPLICITMAXENCHANGEUP;

	  //check if energy density positive and the inversion worked using U2P_HOT
	  if(xxx[0]>edenmin && xxx[0]<edenmax && fret>=-1) okcheck=1;
	  if(okcheck==0) 
	  {
	    if(verbose==2) 
	    {
	      if(fret<-1)  
		printf("inverison failed\n");
	      else
		printf("change in energy density too large\n");
              }
	  }
	  
	  //ANDREW added Tgas limiter 
	  ldouble Tgasold=calc_PEQ_Tfromurho(ppp[UU],ppp[RHO],ix,iy,iz);
	  ldouble Tgasnew=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);
	  if(okcheck==1 && (Tgasnew>Tgasold*IMPLICITMAXTGASCHANGE || Tgasnew<Tgasold/IMPLICITMAXTGASCHANGE)) 
	  {
	    okcheck=0;
	    if(verbose==2) printf("change in gas temperature too large\n");
	  }		
          
	  if(ifncompt)
	  {  
	    ldouble nphmin=ppp[NF0]/RADIMPLICITMAXNPHCHANGE;
	    ldouble nphmax=ppp[NF0]*RADIMPLICITMAXNPHCHANGE;
	    if(okcheck==1 && !(xxx[inph]>nphmin && xxx[inph]<nphmax)) 
	    {
	      okcheck=0;
	      if(verbose==2) printf("change in photon number too large\n");
	     }		

            ldouble Tradnew=calc_ncompt_Thatrad_full(pp,geom);
            ldouble Tradold=calc_ncompt_Thatrad_full(ppp,geom);	    
	    if(okcheck==1 && (Tradnew>Tradold*RADIMPLICITMAXTRADCHANGE || Tradnew<Tradold/RADIMPLICITMAXTRADCHANGE)) 
	    {
	      okcheck=0;
	      if(verbose==2) printf("change in rad temperature too large\n");
	    }		
	  }
	  
	  //ANDREW -- added new limiter for Trad in both cases. 
          else
          {
	    ldouble Tradnew=calc_LTE_TfromE(pp[EE0]);
	    ldouble Tradold=calc_LTE_TfromE(ppp[EE0]);
	    if(okcheck==1 && (Tradnew>Tradold*RADIMPLICITMAXTRADCHANGE || Tradnew<Tradold/RADIMPLICITMAXTRADCHANGE)) 
	    {
	      okcheck=0;
	      if(verbose==2) printf("change in rad temperature too large\n");
	    }		
	  }
	  
	  if(ifrelel)
	  {
            ldouble neurnew=calc_relel_ne(pp); 
            ldouble neurold=calc_relel_ne(ppp);

            if(okcheck==1 && (neurnew>neurold*RADIMPLICITMAXNECHANGE || neurnew<neurold/RADIMPLICITMAXNECHANGE))
	    {
	       okcheck=0;
               if(verbose==2) printf("change in nonthermal number too large\n");
	    }
	  }

	  if(ifentre)
	  {  
	    ldouble Sechange=fabs(xxx[ientre]-ppp[ENTRE]);
            ldouble nethold=calc_thermal_ne(ppp);
            ldouble nethnew=calc_thermal_ne(pp); 
            ldouble rhoethnew=nethnew*MU_I*M_PROTON;

	    //ANDREW - why does using nethold here crash???
            ldouble Tenew = calc_TfromSerho(xxx[ientre],rhoethnew,ELECTRONS,geom->ix,geom->iy,geom->iz);
	    ldouble Teold = calc_TfromSerho(ppp[ENTRE],rhoethnew,ELECTRONS,geom->ix,geom->iy,geom->iz);

	    if(okcheck==1 && (Tenew>Teold*RADIMPLICITMAXTECHANGE || Tenew<Teold/RADIMPLICITMAXTECHANGE || isnan(Teold/Tenew) || isnan(Tenew/Teold)))
	    {
	      okcheck=0;
	      if(verbose==2) printf("change in electron temperature too large\n");
	      if(verbose==2) printf("Sechange: %e %e > %e > %d\n",xxx[ientre],ppp[ENTRE],Sechange/fabs(ppp[ENTRE]),okcheck);
	      if(verbose==2) printf("Techange: %e %e > %e > %d\n",Tenew,Teold,Tenew/Teold,okcheck);
	    }

	  }

	  //all checks are ok
	  if(okcheck==1)
	    break;
	  
	  //if not okcheck  decrease the applied Jacobian fraction
	  if(xxx[0]<=0. && 1)
	  {
	    xiapp*=ppp[sh]/(ppp[sh]+fabs(xxx[0]));
	    xiapp*=1.e-1; //not to land too close to zero, but sometimes prevents from finding proper solution
	  }	  
	  else //u2p error only or too large change
	  {
	    xiapp/=RADIMPLICITDAMPINGFACTOR; 
	  }

	  if(xiapp<MAXRADIMPDAMPING) 
	  {
	    if(verbose) printf("damped unsuccesfully in implicit_4dprim\n");
	    failed=1;
	    break;
	  }
      }
      while(1); //end loop to apply inverse Jacobian corrections 

      //*******************************************      
      //Step 5: Check new solution pp for convergence
      if(failed==0)
      {
	//compute  relative change of iterated quantities
	f3[0]=fabs((pp[sh]-ppp[sh])/ppp[sh]);
	for(i=1;i<np;i++)
	{
	  if(i<4)
	  {
	    f3[i]=pp[i+sh]-ppp[i+sh];
	    f3[i]=fabs(f3[i]/my_max(EPS,fabs(ppp[i+sh])));	
	  }
	  if(ifncompt && i==inph)
	  {
	    f3[inph]=pp[NF0]-ppp[NF0];
	    f3[inph]=fabs(f3[inph]/fabs(ppp[NF0]));	
	  }
	  if(ifentre && i==ientre)
	  {
	    f3[ientre]=pp[ENTRE]-ppp[ENTRE];
	    f3[ientre]=fabs(f3[ientre]/fabs(ppp[ENTRE]));	
	  }
	  if(ifrelel && i>=irelel && i<irelel+NRELBIN)
	  {
	    ie = i-irelel; //bin number
	    f3[i]=pp[NEREL(ie)]-ppp[NEREL(ie)];
	    f3[i]=fabs(f3[i])/fabs(my_max(1.e-25*ppp[RHO]/MU_GAS/M_PROTON, ppp[NEREL(ie)]));
          }
	}

	if(verbose)
	{
	      printf("rel. change:  \n");
	      print_Nvector(f3,np);
	      getch();
	}

	//convergence criterion
	int convrelcheck=1;
	ldouble CONVREL=RADIMPCONVREL;
	ldouble CONVRELERR=RADIMPCONVRELERR;
	if(params[1]==RADIMPLICIT_ENTROPYEQ)
	{
	  CONVREL=RADIMPCONVRELENTR;
	  CONVRELERR=RADIMPCONVRELENTRERR;
	}
  	for (ib=0;ib<np;ib++)
	{
	  if(f3[ib]>CONVREL) 
	  {
	     convrelcheck=0;
	     break;
	  }
        }
	if(convrelcheck || errbase<=CONVRELERR) 
	{
	  if(verbose) printf("\n === success (rel.change) ===\n");
	  break;
	}     
      }
   
      if(iter>MAXITER || failed==1)
      break;
  }
  while(1); //end of main solver loop
  
/*************************************************************/
/*************************************************************/

  //failure
  if(iter>MAXITER || failed==1)
  {
    if(verbose)
    {
      printf("iter (%d) or failed in solve_implicit_lab_4dprim() for frdt=%f (%e)\n",iter,dt,errbest);	  
    }
      
    free_solve_implicit_lab_4dprim(J, iJ, tJ, tiJ, f1, f2, f3, xxx, xxxbest,N);
    return -1;       
  }

//**************************************
  //success
  if(u2pret<-1) 
  {
    if(verbose>1)
      printf("final solution rejected\n");
    free_solve_implicit_lab_4dprim(J, iJ, tJ, tiJ, f1, f2, f3, xxx, xxxbest,N);
    return -1;
  }
  
  if(verbose)
  {
      int corr[2],fixup[2];
  
      print_conserved(uu);
      print_primitives(pp);

      calc_ff_Rtt(pp,&Rtt00,ugas00,geom);
      Tgas00=calc_PEQ_Tfromurho(pp[UU],pp[RHO],geom->ix,geom->iy,geom->iz);
      Trad00=calc_LTE_TfromE(-Rtt00);
      printf("Tgas: %e\n",Tgas00);
      printf("Trad: %e\n",Trad00);
      int ret=f_implicit_lab_4dprim_with_state(uu, pp, &state, uu00, pp00, &state00, dt, geom, f2, params, &err);
      printf("err %e",err);
  }

  //report on average number of iterations
  if(whichprim==RAD && whicheq==RADIMPLICIT_ENERGYEQ)
  {
      #pragma omp critical //to make the counting precise
      global_int_slot[GLOBALINTSLOT_ITERIMPENERRAD]+=iter;
  }
  if(whichprim==MHD && whicheq==RADIMPLICIT_ENERGYEQ)
  {
      #pragma omp critical
      global_int_slot[GLOBALINTSLOT_ITERIMPENERMHD]+=iter;
  }
  if(whichprim==RAD && whicheq==RADIMPLICIT_ENTROPYEQ)
  {
      #pragma omp critical
      global_int_slot[GLOBALINTSLOT_ITERIMPENTRRAD]+=iter;
  }
  if(whichprim==MHD && whicheq==RADIMPLICIT_ENTROPYEQ)
  {
      #pragma omp critical
      global_int_slot[GLOBALINTSLOT_ITERIMPENTRMHD]+=iter;
  }
  if(whicheq==RADIMPLICIT_LTEEQ)
  {
      #pragma omp critical
      global_int_slot[GLOBALINTSLOT_ITERIMPLTE]+=iter;
  }

  //succeeded
  //primitives vector returned to pp[]  
  free_solve_implicit_lab_4dprim(J, iJ, tJ, tiJ, f1, f2, f3, xxx, xxxbest,N);

#endif //RADIATION
  return 0;
}


//**********************************************************************
//******* Apply constraint equations to make non-iterated **************
//******* quantities in pp, uu consistent, and compute state ***********
//**********************************************************************
int
implicit_apply_constraints(ldouble *pp, ldouble *uu, ldouble *uu0, void* ggg, int whichprim)
{

  struct geometry *geom
  = (struct geometry *) ggg;
  
  ldouble (*gg)[5], (*GG)[5], gdet, gdetu,gdetu_inv;
  ldouble ucon[4],ucov[4];
  ldouble urcon[4],urcov[4];

  int ie;
  int u2pret=0.;
  ldouble entrebak,nfbak;
  int corr[3]={0,0,0}, fixup[2]={0,0};

  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;
  gdetu=gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif
  gdetu_inv=1./gdetu;
  
  if(whichprim==MHD) // use MHD primitives in implicit
  {
    //Get 4 velocity
    #ifdef NONRELMHD
    ucon[0]=1.; ucov[0]=-1.;
    #else    
    calc_ucon_ucov_from_prims(pp, geom, ucon, ucov);
    #endif
    
    //rho may be inconsistent on input if iterating MHD primitives
    //correct it from the conserved quantity in uu0
    pp[RHO] = uu0[RHO]*gdetu_inv/ucon[0];

    //ANDREW Should we even update ENTR? What if it's needed for fixups or
    //what if it is the quantity being iterated?    
    //recompute total gas entropy to be consistent with the energy and density
    pp[ENTR] = calc_Sfromu(pp[RHO],pp[UU],geom->ix,geom->iy,geom->iz);

    //apply energy conservation to get ENTRI
    #ifdef EVOLVEELECTRONS
    pp[ENTRI] = entri_from_entre_energy_cons(pp, geom->ix, geom->iy, geom->iz); 
    #endif

    //save photon number primitive to nfbak
    #ifdef EVOLVEPHOTONNUMBER
    nfbak = pp[NF0];
    #endif
    
    //calculate conserved MHD quantities
    p2u_mhd(pp,uu,geom);

    //make opposite changes in the RAD conserved quantities
    uu[EE0] = uu0[EE0] - (uu[1]-uu0[1]);
    uu[FX0] = uu0[FX0] - (uu[2]-uu0[2]);
    uu[FY0] = uu0[FY0] - (uu[3]-uu0[3]);
    uu[FZ0] = uu0[FZ0] - (uu[4]-uu0[4]);

    //and invert back to primitives
    u2pret=u2p_rad(uu,pp,geom,corr);
    
    if(corr[0]>0)
    {
      #ifdef ALLOWRADCEILINGINIMPLICIT
      u2pret=-1;
      #else
      u2pret=-2; //don't allow hitting radiation ceiling in rad when doing iterations
      #endif
    }
    
    //the u2p_rad won't have the right photon number,
    //because the conserved quantity for photon number could be garbage
    //recover photon number from nfbak
    #ifdef EVOLVEPHOTONNUMBER
    #ifdef NONRELMHD
    ufcon[0]=1.; urcov[0]=-1.;
    #else
    calc_urcon_urcov_from_prims(pp, geom, urcon, urcov);
    #endif

    pp[NF0]=nfbak;
    uu[NF0]=pp[NF0]*urcon[0]*gdetu;
    #endif
  }
  else if(whichprim==RAD) // use radiation primitives in implicit
  {
    //save electron  entropy  primitive to entrebak
    //and rel electron primitives to nrelelbak
    #ifdef EVOLVEELECTRONS
    entrebak = pp[ENTRE];
    #ifdef RELELECTRONS
    ldouble nrelelbak[NRELBIN];
    for(ie=0;ie<NRELBIN;ie++) nrelelbak[ie] = pp[NEREL(ie)];
    #endif
    #endif
    
    //calculate conserved RAD quantities
    p2u_rad(pp,uu,geom);
    
    //make opposite changes in the MHD conserved quantities
    uu[RHO]=uu0[RHO];
    uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
    uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
    uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
    uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);

    //copy B field to uu
    #ifdef MAGNFIELD
    uu[B1]=gdetu*pp[B1];
    uu[B2]=gdetu*pp[B2];
    uu[B3]=gdetu*pp[B3];
    #endif
   
    //and invert back to primitives
    int rettemp=0;
    rettemp=u2p_solver(uu,pp,geom,U2P_HOT,0); 
    #ifdef ALLOWFORENTRINF4DPRIM
    if(rettemp<0)
      rettemp=u2p_solver(uu,pp,geom,U2P_ENTROPY,0);
    #endif
    //u2pret=rettemp;
    if(rettemp<0) u2pret=-2; 
    else u2pret=0;
    
    //get gas 4-velocity
    #ifdef NONRELMHD
    ucon[0]=1.; ucov[0]=-1.;
    #else    
    calc_ucon_ucov_from_prims(pp, geom, ucon, ucov);
    #endif

    //ANDREW Should we even update ENTR? What if it's needed for fixups/
    //what if it is the quantity being iterated?    
    //recompute total gas entropy to be consistent with the energy and density
    pp[ENTR] = calc_Sfromu(pp[RHO],pp[UU],geom->ix,geom->iy,geom->iz);
    uu[ENTR] = pp[ENTR]*gdetu*ucon[0];
      
    //the u2p_mhd won't give the right electron entropy or relel counts
    //because the corresponding conserved quantities could be garbage
    //recover from entrebak & nrelelbak
    #ifdef EVOLVEELECTRONS 
    pp[ENTRE]=entrebak;
    uu[ENTRE]=pp[ENTRE]*gdetu*ucon[0];

    #ifdef RELELECTRONS
    for(ie=0;ie<NRELBIN;ie++)
    {
      pp[NEREL(ie)]=nrelelbak[ie];
      uu[NEREL(ie)]=nrelelbak[ie]*gdetu*ucon[0];
    }
    #endif //RELELECTRONS

    //apply energy conservation to get ENTRI
    pp[ENTRI] = entri_from_entre_energy_cons(pp, geom->ix, geom->iy, geom->iz); 
    uu[ENTRI] = pp[ENTRI]*gdetu*ucon[0];
    #endif //EVOLVEELECTRONS    
  }
  else
    my_err("whichprim can only be RAD or MHD! \n");

  if(u2pret>=0 && (corr[0]!=0 || corr[1]!=0)) 
    u2pret=1;
    
  return u2pret;
}

//**********************************************************************
//******* return error vector for implicit solver **********************
//******* solves implicitly four-force source terms ********************
//******* in the lab frame  working on primitives    *******************
//******* rad or hydro (whichprim) *************************************
//**********************************************************************

//ANDREW: - before calling f_implicit_lab_4dprim_with_state, must make sure that the quantities in pp and sss are CONSISTENT
//with the conserved quantities uu0 by calling implicit_apply_constraints and fill_struct_of_state

int
f_implicit_lab_4dprim_with_state(ldouble *uu, ldouble *pp, void *sss, ldouble *uu0, ldouble *pp0, void *sss0, ldouble dt, void* ggg,
				 ldouble *f, int *params, ldouble *err0)
{
  int ret=0,i;
  
  struct geometry *geom
  = (struct geometry *) ggg;
  struct struct_of_state *state0
  = (struct struct_of_state *) sss0;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;

  ldouble (*gg)[5], (*GG)[5], gdet, gdetu, gdetu_inv;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;
  gdetu=gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif
  gdetu_inv = 1. / gdetu;
  
  int whichprim=params[0];
  int whicheq=params[1];
  int whichframe=params[2];
  int opdamplevel=params[4];
  int ifncompt=params[5];
  int ifentre=params[6];
  int ifrelel=params[7];
  
  ldouble opdamp = 1.;
  for (i=0;i<opdamplevel;i++)
    opdamp/=OPDAMPFACTOR;
   
  // equation numbers for non GRRMHD quantities
  int incompt, ientre, irelel, neq=4;
  if(ifncompt)
  {
    incompt=4;
    neq=5;
  }
  if(ifentre)
  {
    if(ifncompt)
    {
      ientre=5;
      neq=6;
    }
    else
    {
      ientre=4;
      neq=5;
    }
  }
  if(ifrelel)
  {
    irelel=ientre+1;
    neq+=NRELBIN;
  }
  
  //error vector initialized to zero
  ldouble *err=(ldouble*) malloc(neq*sizeof(ldouble));
  int ib;
  for (ib=0;ib<neq;ib++)
    err[ib]=0.;
  
  //species temperatures
  ldouble Tgas0, Tgas, Ti0, Ti, Te0, Te;
  Tgas0 = state0->Tgas;
  Ti0 = state0->Ti;
  Te0 = state0->Te;  
  Tgas = state->Tgas;
  Ti = state->Ti;
  Te = state->Te;
  
  //Ehat etc
  ldouble uconf[4], Rtt0, Rtt, B, Ehat0, Ehat, dtau;
  B = sigma_rad_over_pi * Te * Te * Te * Te;
  Ehat0 = state0->Ehat;
  Rtt0 = -Ehat0;
  Ehat = state->Ehat;
  Rtt = -Ehat;
  
  for (i = 0; i < 4; i++)
  {
    uconf[i] = state->ucon[i];
  }
  dtau=dt/uconf[0];
  
  //radiative four-force
  ldouble Gi[4],Giff[4];
  ldouble Gith[4], Githff[4];
  ldouble Hi[4]={0.,0.,0.,0.}, Hiff[4]={0.,0.,0.,0.};
  
  //ANDREW calculate change in nonthermal energy density from end of explicit
  ldouble relel_dudtau=0.;
  #ifdef RELELECTRONS
  //relel_dudtau = (state->uenth - state0->uenth)/dtau;
  relel_dudtau = (calc_relel_uint(pp) - calc_relel_uint(pp0)) / dtau;
  #endif
  
  //calculate radiation - gas coupling
  int relel_g0_type=1;
  ldouble one_over_scale = calc_all_Gi_with_state(pp, state, ggg, Giff, Gi, Githff, Gith, relel_dudtau, relel_g0_type);

  for(i=0;i<4;i++)
  {
    Gi[i]*=opdamp;
    Giff[i]*=opdamp;
    Gith[i]*=opdamp;
    Githff[i]*=opdamp;
  }
  
  // lower indices for lab frame only
  indices_21(Gi,Gi,gg);
  indices_21(Gith, Gith, gg);
  
  //artificial heating
#ifdef HEATINIMP_CONSTANT
  Hiff[0] = HEATINIMP_CONSTANT;
  for(i=0;i<4;i++)
    Hi[i]=Hiff[0]*uconf[i];
  indices_21(Hi,Hi,gg);
#endif

  //********************************************************************  
  //********************************************************************
  //implicit equations
  //********************************************************************
  //********************************************************************
  
  //photon number error function -- always in lab frame
#ifdef EVOLVEPHOTONNUMBER
  if(ifncompt)
  {
    ldouble nsource=calc_nsource_with_state(pp, state, ggg);
    nsource*=opdamp;
    f[incompt] = uu[NF0] - uu0[NF0] - dt * gdetu * nsource;
    if(fabs(f[incompt])>SMALL) err[incompt]=fabs(f[incompt])/(fabs(uu[NF0])+fabs(uu0[NF0])+fabs(dt*gdetu*nsource));
    else err[incompt]=0.;
  }
#endif
  
  //thermal electron error function -- always in fluid frame
#ifdef EVOLVEELECTRONS
  if(ifentre)
  {
    //emission - absorption, thermal only
    ldouble EA=Githff[0];
    
    //Thermal coloumb coupling
    ldouble CC=0.;
#ifndef  SKIPCOULOMBCOUPLING
    CC=calc_CoulombCoupling(pp,ggg);
#endif
    
    //Nonthermal coloumb coupling + nonthermal cooling
    ldouble CC_relel=0.0;
    ldouble cool_relel_dn=0.0;
    ldouble cool_relel_dq=0.0;
    ldouble mu=0.0;
    
#ifdef RELELECTRONS
    
    CC_relel= calc_relel_CoulombCoupling_from_state(pp, state); //Coulomb coupling between thermal and nonthermal
    cool_relel_dn = calc_relel_cool_dn_from_state(pp, state); //Cooling dn into the thermal
    cool_relel_dq = calc_relel_cool_dq_from_state(pp, state); //Cooling dq into the thermal (positive)
    
    //Chemical potential
    ldouble theta= kB_over_me * Te;
    ldouble neth=state->ne;
    mu = chemical_potential_short(theta, neth);
#endif

    //ANDREW - in the entropy equation, should there be a chemical potential term due to change in volume?
    //is this equation without relel consistent with equation in terms of entropy per particle? 
    ldouble Qe = CC + EA + CC_relel + cool_relel_dq - mu*cool_relel_dn;
    f[ientre]=Te*(pp[ENTRE] - pp0[ENTRE]) - dtau * Qe;
    if(fabs(f[ientre])>SMALL) err[ientre]=fabs(f[ientre])/(fabs(pp[ENTRE]*Te)+fabs(pp0[ENTRE]*Te)+fabs(dtau*Qe)); else err[ientre]=0.;
  }
  
  // relativistic electron error function -- in fluid frame
#ifdef RELELECTRONS
  if(ifrelel)
  {
    int ibin;
    ldouble frel[NRELBIN], frelmag[NRELBIN];
    
    calc_relel_f_and_fmag_from_state(pp,state,pp0,geom,dtau,frel,frelmag);
    
    for(ibin=0;ibin<NRELBIN;ibin++)
    {
      f[irelel+ibin]=frel[ibin];
      err[irelel+ibin]=fabs(frel[ibin]/my_max(frelmag[ibin], 1.e-20*pp[RHO]*one_over_mue_mp));
    }
  }
#endif //RELELECTRONS
#endif //EVOLVEELECTRONS
  
  //update four fource to include artificial heating
  for(i=0;i<4;i++)
    Gi[i]+=Hi[i];
  Giff[0]+=Hiff[0];
  
  //errors in momenta -- always in lab frame
  if(whichprim==MHD) //mhd-primitives
  {
    f[1] = uu[2] - uu0[2] - dt * gdetu * Gi[1];
    f[2] = uu[3] - uu0[3] - dt * gdetu * Gi[2];
    f[3] = uu[4] - uu0[4] - dt * gdetu * Gi[3];
    
    if(fabs(f[1])>SMALL) err[1]=fabs(f[1])/(1.e-20*uu[UU]+fabs(uu[2])+fabs(uu0[2])+fabs(dt*gdetu*Gi[1])); else err[1]=0.;
    if(fabs(f[2])>SMALL) err[2]=fabs(f[2])/(1.e-20*uu[UU]+fabs(uu[3])+fabs(uu0[3])+fabs(dt*gdetu*Gi[2])); else err[2]=0.;
    if(fabs(f[3])>SMALL) err[3]=fabs(f[3])/(1.e-20*uu[UU]+fabs(uu[4])+fabs(uu0[4])+fabs(dt*gdetu*Gi[3])); else err[3]=0.;
  }
  else if(whichprim==RAD) //rad-primitives
  {
    f[1] = uu[FX0] - uu0[FX0] + dt * gdetu * Gi[1];
    f[2] = uu[FY0] - uu0[FY0] + dt * gdetu * Gi[2];
    f[3] = uu[FZ0] - uu0[FZ0] + dt * gdetu * Gi[3];
    
    if(fabs(f[1])>SMALL) err[1]=fabs(f[1])/(1.e-20*uu[EE]+fabs(uu[FX0])+fabs(uu0[FX0])+fabs(dt*gdetu*Gi[1])); else err[1]=0.;
    if(fabs(f[2])>SMALL) err[2]=fabs(f[2])/(1.e-20*uu[EE]+fabs(uu[FY0])+fabs(uu0[FY0])+fabs(dt*gdetu*Gi[2])); else err[2]=0.;
    if(fabs(f[3])>SMALL) err[3]=fabs(f[3])/(1.e-20*uu[EE]+fabs(uu[FZ0])+fabs(uu0[FZ0])+fabs(dt*gdetu*Gi[3])); else err[3]=0.;
  }
  
  //errors in energy/entropy - in either frame
  /***** LAB FRAME ENERGY/ENTROPY EQS *****/
  if(whichframe==RADIMPLICIT_LAB)
  {
    if(whichprim==RAD) //rad-primitives
    {
      if(whicheq==RADIMPLICIT_ENERGYEQ)
      {
        f[0] = uu[EE0] - uu0[EE0] + dt * gdetu * Gi[0];
        if(fabs(f[0])>SMALL) err[0] = fabs(f[0])/(fabs(uu[EE0]) + fabs(uu0[EE0]) + fabs(dt*gdetu*Gi[0]));
        else err[0]=0.;
      }
      else if(whicheq==RADIMPLICIT_ENTROPYEQ)
      {
        my_err("not implemented RADIMPLICIT_ENTROPYEQ\n");
      }
      else
        my_err("not implemented 3\n");
    }
    else if(whichprim==MHD) //hydro-primitives
    {
      if(whicheq==RADIMPLICIT_ENERGYEQ)
      {
        f[0] = uu[UU] - uu0[UU] - dt * gdetu * Gi[0];
        if(fabs(f[0])>SMALL) err[0] = fabs(f[0])/(fabs(uu[UU]) + fabs(uu0[UU]) + fabs(dt * gdetu * Gi[0]));
        else err[0]=0.;
      }
      else if(whicheq==RADIMPLICIT_ENTROPYEQ)
      {
        my_err("not implemented RADIMPLICIT_ENTROPYEQ\n");
      }
      else
        my_err("not implemented 4\n");
    }
  }
  
  /***** FLUID FRAME ENERGY/ENTROPY EQS ******/
  else if(whichframe==RADIMPLICIT_FF)
  {
    
    //fluid frame energy equation:
    if(whichprim==RAD) //rad-primitives
    {
      if(whicheq==RADIMPLICIT_ENERGYEQ)
      {
        f[0]=Ehat - Ehat0 + dtau * Giff[0];
        err[0]=fabs(f[0])/(fabs(Ehat) + fabs(Ehat0) + fabs(dtau * Giff[0]));
      }
      else if(whicheq==RADIMPLICIT_ENTROPYEQ)
      {
        ldouble S0 = state0->Sgas;
        ldouble S = state->Tgas;
        
        f[0]=Tgas*(S-S0) - dtau * Giff[0];
        err[0]=fabs(f[0])/(fabs(Tgas*S) + fabs(Tgas*S0) + fabs(dtau * Giff[0]));
      }
      else
        my_err("not implemented 2\n");
    }
    else if(whichprim==MHD) //mhd-primitives
    {
      if(whicheq==RADIMPLICIT_ENERGYEQ)
      {
        f[0]=pp[UU] - pp0[UU] - dtau * Giff[0];
        err[0]=fabs(f[0])/(fabs(pp[UU]) + fabs(pp0[UU]) + fabs(dtau * Giff[0]));
      }
      else if(whicheq==RADIMPLICIT_ENTROPYEQ)
      {
        ldouble S0 = state0->Sgas;
        ldouble S = state->Tgas;
        
        f[0]=Tgas*(S-S0) - dtau * Giff[0];
        err[0]=fabs(f[0])/(fabs(Tgas*S) + fabs(Tgas*S0) + fabs(dtau * Giff[0]));
      }
      else
        my_err("not implemented 2\n");
    }
  }
  
  //Return error vector
  int isinf=0;
  *err0=0.;
  for (ib=0;ib<neq;ib++)
  {
    if(!isfinite(f[0])) isinf=1;
    if(err[ib]>*err0)
      *err0=err[ib];
  }
  
  free(err);
  
  return ret;
}


//**********************************************************************
//Print solver state
//**********************************************************************
int
print_state_implicit_lab_4dprim (int iter, ldouble *x, ldouble *f,ldouble err,int N)
{
  printf ("\n============\niter = %3d\n============\n",iter);
  int i;
  printf("x   : ");
  for(i=0;i<N;i++)
    {
      printf("%.8e ",x[i]);
    }
  printf("\n============\nf(x): ");
   for(i=0;i<N;i++)
    {
      printf("%.8e ",f[i]);
    }
   printf("\n============\n");
   printf("err : %e\n============\n",err);
   return 0;
}


//**********************************************************************
//free memory used in the solver
//**********************************************************************
int
free_solve_implicit_lab_4dprim(ldouble** J, ldouble** iJ, ldouble *tJ, ldouble *tiJ,
			       ldouble *f1, ldouble *f2, ldouble *f3, ldouble *xxx, ldouble *xxxbest,int N)
{
  int ib;
  for(ib=0;ib<N;ib++)
    {
      free(J[ib]);
      free(iJ[ib]);
    }
	    
  free(f1);
  free(f2);
  free(f3);
  free(xxx);
  free(xxxbest);
  free(tJ);
  free(tiJ);
  free(J);
  free(iJ);
  return 0;
}

///*************************************************************************/
///******* 1D solver in energy densities ***********************************/
///******* totenergy is the total energy density in ff *********************/
///******* ratio is the radio of energy denisities between ff and rad.rest frame */
///*************************************************************************/
///*********** always in FF frame! ************/
/*! \fn ldouble f_implicit_1dprim_err(ldouble xxx,ldouble *uu0,ldouble *pp0,ldouble dtau,void *ggg,int *params,ldouble totenergy,ldouble ratio, int verbose)
 \brief Calculates error for 1D bisection routine
 
 \param[in] xxx either pp[UU] or pp[EE0], depending on which energy equation is being considered
 \param[in] *uu0 pre-implicit conserveds
 \param[in] *pp0 pre-implicit primitives
 \param[in] dtau time step in fluid frame
 \param[in] *ggg geometry structure for current cell
 \param[in] totenergy total energy: Ehat + ugas
 \param[in] ratio radiation frame to fluid frame: Ehat/pp0[EE0]
 \param[in] verbose >0 if details are to be printed, 0 otherwise
 
 \returns Error in the appropriate equation:\n
   MHD energy equation: err = ugas - pp0[UU] - dtau * Gi[0]\n
   MHD entropy equation: err = Tgas * (S - S0) - dtau * Gi[0]\n
   Radiation energy equation: err = Ehat - pp0[EE]/ratio + dtau * Gi[0]\n
   Radiation entropy equation: err = Tgas * (S - S0) - dtau * Gi[0]\n
 
 \todo We only need Gi[0], not the whole vector. Check calc_Gi to see if we could simplify the computation.
 */

// Ramesh-modified version of f_implicit_1dprim_err
ldouble
f_implicit_1dprim_err(ldouble xxx, ldouble *uu0, ldouble *pp0, void *sss0, void *sss, ldouble dtau, void *ggg, int *params, ldouble totenergy, ldouble ratio, int verbose)
{  
  int ret_update, ret_fill;
  
  struct geometry *geom
  = (struct geometry *) ggg;
  struct struct_of_state *state0
  = (struct struct_of_state *) sss0;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;
  
  int whichprim = params[0];
  int whicheq = params[1];
  
  // printf("input: %e %e %e %e\n",xxx,pp0[ENTR],dt,totenergy);
  
  ldouble ugas = 0., Ehat = 0., err=0.;
#ifdef RADIATION  
  if(whichprim == RAD)
  {
    Ehat = xxx * ratio;
    ugas = totenergy - Ehat;
  }
  
  if(whichprim == MHD)
  {
    ugas = xxx;
    Ehat = totenergy - ugas;
  }
  
  //update the new state
  ldouble pp[NV];
  int i;
  for (i = 0; i< NV; i++) {
    pp[i] = pp0[i];
  }
  pp[UU] = ugas;
  pp[EE] = Ehat/ratio;
  ret_fill = fill_struct_of_state(pp, geom, state);
  
  ldouble Gi[4];
  
  //fluid frame
  //calculate local rate of change in relel energy density
  ldouble relel_dudtau=0.0;
  #ifdef RELELECTRONS
  relel_dudtau = ((state->uenth) - (state0->uenth))/dtau;
  #endif
  
  calc_Gi_with_state(pp, state, ggg, Gi, relel_dudtau, 0, 1);
  
  if (whichprim == RAD)
  {
    if (whicheq == RADIMPLICIT_ENERGYEQ)
    {
      err = Ehat - pp0[EE] / ratio + dtau * Gi[0];
      
      if (verbose)
      {
        // This was included to print out details during tests
        ldouble Tgas0, Tgas, Ti, Te, TradBB;
        TradBB=calc_LTE_TfromE(Ehat);
        Tgas0 = state0->Tgas;
        Tgas = state->Tgas;
        Ti = state->Ti;
        Te = state->Te;

        printf("\nEhat, pp0/ratio, rho, pp[UU], Tgas, pp0[UU], Tgas0, TradBB, Gi0, err: %e %e %e %e %e %e %e %e %e %e\n", Ehat, pp0[EE]/ratio, pp0[RHO], pp[UU], Tgas, pp0[UU], Tgas0, TradBB, Gi[0], err);
      }
    }  // if (whicheq == RADIMPLICIT_ENERGYEQ)
    else if (whicheq == RADIMPLICIT_ENTROPYEQ)
    {
      ldouble Tgas0, Tgas, Ti, Te, S0, S;
      Tgas = state->Tgas;
      S0 = state0->Sgas;
      S = state->Sgas;
	     
      err = Tgas * (S - S0) - dtau * Gi[0];
    }  // else if (whicheq == RADIMPLICIT_ENTROPYEQ)
    else
    {
      my_err("not implemented 2\n");
    }
  }  // if (whichprim == RAD)
  
  if(whichprim == MHD)
  {
    if (whicheq == RADIMPLICIT_ENERGYEQ)
    {
      err = pp[UU] - pp0[UU] - dtau * Gi[0];
    }
    else if (whicheq == RADIMPLICIT_ENTROPYEQ)
    {
      ldouble Tgas0, Tgas, Ti, Te, S0, S;
      Tgas = state->Tgas;
      S0 = state0->Sgas;
      S = state->Sgas;
      
      err = Tgas * (S - S0) - dtau * Gi[0];
    }
    else
      my_err("not implemented 2\n");
  }  // if(whichprim == MHD)

#endif  
  return err;

}  // f_implicit_1dprim_err


///////////////////////////////////////////////////////////////
/*! \fn int solve_implicit_lab_1dprim(ldouble *uu0,ldouble *pp0,void *ggg,ldouble dt,int verbose,ldouble *ppout)
 \brief 1D solver of energy equation via the bisection method
 
 \param[in] *uu0 pre-implicit conserveds
 \param[in] *pp0 pre-implicit primitives
 \param[in] *ggg geometry structure for current cell
 \param[in] dt time step in lab frame
 \param[in] verbose >0 if details are to be printed, 0 otherwise
 \param[out] *ppout primitives after solving the energy equation via bisection
 
 The routine decides whether to solve the MHD or radiation equation. In each case, it chooses between the energy equation and the entropy equation. All computations are done in the fluid frame.
*/
// Ramesh-modified version

int
solve_implicit_lab_1dprim(ldouble *uu0, ldouble *pp0, void *sss0, void *ggg, ldouble dt,  int verbose, ldouble *ppout, void *sss)
{
  int i1,i2,i3,iv,i,j,iter;
  ldouble pp[NV],uu[NV];
  ldouble errlo, errhi, errmid;
  ldouble xlo,xhi,xmid;
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  
  struct geometry *geom
  = (struct geometry *) ggg;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;
  struct struct_of_state *state0
  = (struct struct_of_state *) sss0;
  
  int ix=geom->ix;
  int iy=geom->iy;
  int iz=geom->iz;
  
  //temporary using local arrays
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
  
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  
  //ff rad energy density
  ldouble Rtt0, Ehat, ugas0[4];
  Ehat = state->Ehat;
  ugas0[0] = state->ucon[0];
  
  //time step in fluid frame
  ldouble dtau;
  dtau=dt/ugas0[0];
    
  //total ff energy
  ldouble totenergy;
  totenergy = Ehat + pp0[UU];
  
  if(verbose) printf("starting 1D \n");
  
  //residual function parameters
  int params[3],whichprim;
  if(Ehat<1.e-3*pp0[UU]) //hydro preffered
    whichprim=RAD;
  else
    whichprim=MHD;
  
#ifdef BISECT_USING_MHD_ONLY
  whichprim=MHD;
#endif
  
  //solve in rad or mhd temperature?
  params[0]=whichprim;
  //energy or entropy equation to solve
  params[1]=RADIMPLICIT_ENERGYEQ;
  //frame for energy/entropy equation to solve
  //for now, only FF, no lab frame implemented
  params[2]=RADIMPLICIT_FF;
  
  //local vectors
  for(iv=0;iv<NV;iv++)
  {
    uu[iv]=uu0[iv];
    pp[iv]=pp0[iv];
  }

  // Calculate temperature to identify problem cases
  ldouble Ti,Te, ppUU = pp0[UU], ppEE = pp0[EE0];
  ldouble Tgas = state->Tgas;
  ldouble ratio = Ehat/pp0[EE0]; //ratio of energy densities between the ff and rad rest frame
  ldouble FAC=3., FAC1=FAC+1.;

  // Flag unusually large gas temperature. Reduce gas internal energy to a more reasonable value  
  if (Tgas > 1.e14)
  {
    ppUU = pp0[UU] * 1.e14 / Tgas;
    ppEE = (totenergy - ppUU) / ratio;
  }
  
  if(verbose) print_NVvector(pp0);
  if(verbose) printf("enden: %e %e %e\n",Ehat,pp0[UU],totenergy);

  // Choose xlo and xhi
  if(whichprim==MHD)
  {
    xlo = ppUU / 1.05;
    xhi = my_min(ppUU*1.05, (ppUU+FAC*totenergy)/FAC1);  // avoid exceeding totenergy
    if(verbose) printf("working on MHD\n");
  }
  
  if(whichprim==RAD)
  {
    xlo = ppEE / 1.05;
    xhi = my_min(ppEE*1.05, (ppEE+FAC*totenergy/ratio)/FAC1);  // avoid exceeding totenergy
    if(verbose) printf("working on RAD\n");
  }
  
  errlo=f_implicit_1dprim_err(xlo,uu0,pp0,state0,state,dtau,geom,params,totenergy,ratio,verbose);
  errhi=f_implicit_1dprim_err(xhi,uu0,pp0,state0,state,dtau,geom,params,totenergy,ratio,verbose);
  
  if(verbose) printf("%d >>> [%e : %e] > ",0,xlo,xhi);
  if(verbose) printf("[%e : %e]\n",errlo,errhi);
  
  if(verbose) printf("searching for valid brackets\n");
  
  //1dprim
  //first try low the lower bracket
  int MAXITER=50;
  ldouble CONV=1.e-3; // no need for too much precision in bisection
  
  iter=0.;
    
  while(errlo*errhi>0.)
  {
    //by virtue of being in the loop, errlo, errhi have the same sign
    iter++;
    
    if (errlo > 0.)
    { //if f>0 , then extend in low end,
      //remember the old lower bracket and error
      xhi = xlo;
      errhi = errlo;
      
      xlo/=FAC;
      errlo=f_implicit_1dprim_err(xlo,uu0,pp0,state0,state,dtau,geom,params,totenergy,ratio,verbose);
    }
    else
    { // if errlo < 0, expand on the high end
      //remember the old upper bracket and error
      xlo = xhi;
      errlo = errhi;
      
      // Make sure xhi is below the total energy
      if(whichprim==RAD)
        xhi = my_min(xhi*FAC, (xhi+FAC*totenergy/ratio)/FAC1);
      if(whichprim==MHD)
        xhi=my_min(xhi*FAC, (xhi+FAC*totenergy)/FAC1);
      
      errhi=f_implicit_1dprim_err(xhi,uu0,pp0,state0,state,dtau,geom,params,totenergy,ratio,verbose);
    }
    
    if(verbose) printf("%d (%d) >>> [%e : %e] > ",0,iter,xlo,xhi);
    if(verbose) printf("[%e : %e]\n",errlo,errhi);
    if(isnan(errlo) || isnan(errhi)) iter=MAXITER;
    if(iter>=MAXITER) break;
  }
  
  if(iter>=MAXITER)
  {
    if(verbose) {printf("brackets not found at %d %d %d!\n",ix+TOI,iy+TOJ,iz+TOK);}
    PLOOP(i) ppout[i]=pp0[i];
    return -1;
  }
  
  if(verbose) printf("brackets found!\n");
  
  // Bracket has been found. Now bisect to improve the solution
  do
  {
    iter++;
    xmid=0.5*(xlo+xhi); //new estimate
    errmid=f_implicit_1dprim_err(xmid,uu0,pp0,state0,state,dtau,geom,params,totenergy,ratio,verbose); //new error
    
    if(errmid*errlo>0.) //same sign as the lower bracket
    {
      xlo=xmid;
      errlo=errmid;
    }
    else
    {
      xhi=xmid;
      errhi=errmid;
    }
    if(verbose) printf("%d >>> [%e : %e : %e] > ",iter,xlo,xmid,xhi);
    if(verbose) printf("[%e : %e : %e]\n",errlo,errmid,errhi);
    
    if(isnan(xmid) || isnan(errmid))
    {
      my_err("nan in 1dprim\n");
    }
  }
  while(fabs((xhi-xlo)/xmid)>CONV);
  
  if(verbose) printf("solution found: %e\n",xmid);
  
  // Set ppout equal to the new set of primitives, and update energy
  PLOOP(i) ppout[i]=pp0[i];
  if(whichprim==RAD)
  {
    ppout[EE0] = xmid;
    ppout[UU] = totenergy - xmid * ratio;
  }
  else
  {
    ppout[UU] = xmid;
    ppout[EE0] = (totenergy - xmid) / ratio;
  }
  
  // Flag unusually large gas temperature. 
  if (Tgas*ppout[UU]/pp0[UU] > 1.e14)
  {
    printf("\nUNPHYSICALLY HIGH TEMPERATURE!: Initial: ix, iy, iz, pp0[EE0], Ehat, pp0[UU], Tgas, Trad, TradBB: %d %d %d %e %e %e %e %e %e\n", ix, iy, iz, pp0[EE0], pp0[EE0] * ratio, pp0[UU], Tgas, state->Trad, state->TradBB);
    printf("Final: ppout[EE0], Ehat, ppout[UU], Tgas: %e %e %e %e\n", ppout[EE0], ppout[EE0] * ratio, ppout[UU], Tgas * ppout[UU] / pp0[UU]);
  }
  
  if(verbose) print_NVvector(ppout);
  //if(verbose) getchar();
  
  return 0;
}


/*************************************************************************/
/************ 1D solver in photon number *********************************/
/*********** always in RADIATION frame! **********************************/
/*************************************************************************/

ldouble
f_implicit_photon_rad_1dprim_err(ldouble xxx, ldouble *uu0, ldouble *pp0, void *sss, ldouble dtaurad, void *ggg)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;
  
  int ret_state;
  ldouble err, nsource;
  
  //update the new state
  ldouble pp[NV];
  int i;
  for (i = 0; i< NV; i++) {
	  pp[i] = pp0[i];
  }
  
  // Set current nphoton and update state
  pp[NF0] = xxx;
  ret_state = update_state_for_nphoton(pp, geom, state);

  nsource = calc_nsource_with_state(pp, state, geom);
  err = pp[NF0] - pp0[NF0] - dtaurad * nsource;
    
  return err;
}

int
solve_implicit_nphoton_rad_1dprim(ldouble *uu0,ldouble *pp0, void *sss, void *ggg, ldouble dt, int verbose, ldouble *ppout, void *sssout)
{
  int i1,i2,i3,iv,i,j,iter;
  ldouble pp[NV],uu[NV]; 
  ldouble errlo, errhi, errmid;
  ldouble errloold, errhiold;
  ldouble xlo_old, xhi_old;
  ldouble xlo,xhi,xmid;
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  //verbose=1; //ANDREW
  
  struct geometry *geom
    = (struct geometry *) ggg;
  struct struct_of_state *state0
  = (struct struct_of_state *) sss; // this is the input state
  struct struct_of_state *state
  = (struct struct_of_state *) sssout;  // this is updated during the iterations

  int ix=geom->ix;
  int iy=geom->iy;
  int iz=geom->iz;

  //temporary using local arrays
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
  
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  
  //time step in rad frame
  ldouble dtau = dt / state->urfcon[0];
  
  if(verbose) printf("starting 1D NPHOTON \n");
  if(verbose) print_NVvector(pp0);
  
  //local vectors
  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=uu0[iv]; 
      pp[iv]=pp0[iv];     
    }
  
  xmid=pp0[NF0];      
  xlo=pp0[NF0];
  xhi=pp0[NF0];
  
  //test signs
  errmid=f_implicit_photon_rad_1dprim_err(xmid, uu0, pp0, state, dtau, geom);
  errlo = errmid;
  errhi = errmid;

  if(verbose) printf("%d >>> [%e : %e : %e] > ",0,xlo,xmid,xhi);
  if(verbose) printf("[%e : %e : %e]\n",errlo,errmid,errhi);  
  if(verbose) printf("searching for valid brackets\n");

  //1dprim
  //first try low the lower bracket
  int MAXITER=50;
  ldouble FAC=2.;
  ldouble CONV=1.e-4; //cannot be too low!
  
  iter=0.;
 
  // First bracket the solution  
  xlo_old = xlo;
  xhi_old = xhi;
  errloold = errlo;
  errhiold = errhi;
  while(errlo*errhi>0.){ 
  //by virtue of being in the loop, errlo, errhi have the same sign
      iter++;
      
      if (errlo > 0.) { //if f>0 , then extend in low end,
	
      //rename the old lower bracket as the upper bracket
      xhi = xlo;
      errhi = errlo;
      
      xlo/=FAC;
        errlo = f_implicit_photon_rad_1dprim_err(xlo, uu0, pp0, state, dtau, geom);
      }
      else { // if errlo < 0, expand on the high end
      
      //rename the old upper bracket as the lower bracket
      xlo = xhi;
      errlo = errhi;
       
      xhi*=FAC;

        errhi = f_implicit_photon_rad_1dprim_err(xhi, uu0, pp0, state, dtau, geom);
      }
       
      if(verbose) printf("%d (%d) >>> [%e : %e : %e] > ",0,iter,xlo,xmid,xhi);
      if(verbose) printf("[%e : %e : %e]\n",errlo,errmid,errhi);  
      if(isnan(errlo) || isnan(errhi)) iter=MAXITER;
      if(iter>=MAXITER) break;
    }

  if(iter>=MAXITER)
    {
      if(verbose) {printf("brackets not found at %d %d %d!\n",ix+TOI,iy+TOJ,iz+TOK); getchar();}
      PLOOP(i) ppout[i]=pp0[i];
      return -1;
    }      
    
  if(verbose) printf("brackets found!\n");
    
  //Now solve by bisection
  do
    {
      iter++;
      xmid=0.5*(xlo+xhi); //new estimate
      errmid = f_implicit_photon_rad_1dprim_err(xmid, uu0, pp0, state, dtau, geom); //new error
      
      if(errmid*errlo>0.) //same sign as the lower bracket
	{
	  xlo=xmid;
	  errlo=errmid;
	}
      else
	{
	  xhi=xmid;
	  errhi=errmid;
	}
      if(verbose) printf("%d >>> [%e : %e : %e] > ",iter,xlo,xmid,xhi);
      if(verbose) printf("[%e : %e : %e]\n",errlo,errmid,errhi);  

      if(isnan(xmid) || isnan(errmid))
      {
        my_err("nan in 1dprim\n");
      }
    }
  while(fabs((xhi-xlo)/xmid)>CONV);

  if(verbose) printf("solution found: %e\n",xmid);

  //converting to new set of primitives
  PLOOP(i) pp[i]=pp0[i];
  pp[NF0] = xmid;
  p2u(pp,uu,geom);

  //updating entropy
  //pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],geom->ix,geom->iy,geom->iz);

  //returning the new set of primitives
  PLOOP(i) ppout[i]=pp[i];

  if(verbose) print_NVvector(ppout);
  if(verbose) getchar();

  return 0; 
}


/************************************************************************/
/******* explicit radiative source term  ********************************/
/************************************************************************/

int explicit_rad_source_term(int ix,int iy, int iz,ldouble dt)
{
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEEXPLICIT);   
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  ldouble del4[NRADVAR],delapl[NV];
  int iv;

  //applied explicitly directly in lab frame
  solve_explicit_lab(ix,iy,iz,dt,del4,0);

  apply_rad_source_del4(ix,iy,iz,del4);

  set_cflag(RADIMPFIXUPFLAG,ix,iy,iz,0); 

  return 0;
}

//**********************************************************************
//******* solves explicitly gas - radiation interaction  ***************
//******* in the lab frame, returns vector of deltas *******************
//**********************************************************************

int
solve_explicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int verbose)
{
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  int iv,ret;
  ldouble pp[NV],uu[NV];
  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz); //primitives corresponding to zero-state  
      uu[iv]=get_u(u,iv,ix,iy,iz);  
    }

  ret= solve_explicit_lab_core(uu,pp,&geom,dt,deltas,verbose);

  return ret;
}

int
solve_explicit_lab_core(ldouble *uu,ldouble *pp,void* ggg,ldouble dt,ldouble* deltas,int verbose)
{
#ifdef RADIATION
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;
  gdetu=gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif

  ldouble Gi[4];
  calc_Gi(pp,geom,Gi,0.0,1,0);
  indices_21(Gi,Gi,geom->gg);

  deltas[0]=-Gi[0]*dt*gdetu;
  deltas[1]=-Gi[1]*dt*gdetu;
  deltas[2]=-Gi[2]*dt*gdetu;
  deltas[3]=-Gi[3]*dt*gdetu;

  #ifdef EVOLVEPHOTONNUMBER
  ldouble nsource=calc_nsource(pp,geom);
  deltas[5]=dt*gdetu*nsource;
  #endif

  if(verbose)
    {
      ldouble delapl[NV];
      ldouble uu0[NV],pp0[NV];
      
      int iv;
      for(iv=0;iv<NV;iv++)
	{
	  delapl[iv]=0.;
	  uu0[iv]=uu[iv];
	  pp0[iv]=pp[iv];
	}
      
      delapl[1]=-deltas[0];
      delapl[ENTR]=-deltas[0];
      delapl[2]=-deltas[1];
      delapl[3]=-deltas[2];
      delapl[4]=-deltas[3];
      delapl[EE0]=deltas[0];
      delapl[FX0]=deltas[1];
      delapl[FY0]=deltas[2];
      delapl[FZ0]=deltas[3];

      for(iv=0;iv<NV;iv++)
	{
	  uu0[iv]+=delapl[iv];
	}

      int corr[3],fixup[2];

      printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
      printf("\n@@@@@@@@ EXPLICIT @@@@@@@@@@@@");
      printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n");
  
      print_Nvector(uu,NV);
      print_Nvector(pp,NV);

      u2p(uu0,pp0,geom,corr,fixup,0);
      printf("%d %d\n",corr[0],corr[1]);

      print_Nvector(uu0,NV);
      print_Nvector(pp0,NV);
      print_4vector(deltas);
      getchar();
    }
#endif   
  return 0;

}


//************************************************************************/
//******* aplying radiative four velocity-related changes in conserved ***/
//************************************************************************/
int
apply_rad_source_del4(int ix,int iy,int iz,ldouble *del4)
{
  ldouble delapl[NV];
  int iv;  
  ldouble pp0[NV];
#ifdef RADIATION
  for(iv=0;iv<NV;iv++)
    {
      delapl[iv]=0.;
      pp0[iv]=get_u(p,iv,ix,iy,iz);
    }

  delapl[RHO]=0.;
  delapl[UU]=-del4[0];
  delapl[VX]=-del4[1];
  delapl[VY]=-del4[2];
  delapl[VZ]=-del4[3];
  delapl[ENTR]=-del4[0];
#ifdef MAGNFIELD
  delapl[B1]=0.;
  delapl[B2]=0.;
  delapl[B3]=0.;
#endif
  delapl[EE0]=del4[0];
  delapl[FX0]=del4[1];
  delapl[FY0]=del4[2];
  delapl[FZ0]=del4[3];

#ifdef EVOLVEPHOTONNUMBER
  delapl[NF]=del4[5];
#endif

  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz, get_u(u,iv,ix,iy,iz)+delapl[iv] );
    }
#endif
  return 0;
}

//****************************************************************************
//****** wrapper to calculate contravariant four-force ***********************
//****************************************************************************

ldouble
calc_Gi(ldouble *pp, void *ggg, ldouble Gi[4], ldouble relel_dudtau, int type, int reltype)
{
  // type = 0: thermal + nonthermal in fluid frame
  // type = 1: thermal + nonthermal in lab frame
  // type = 2: thermal ONLY in fluid frame
  // type = 3: thermal ONLY in lab frame
  
  // reltype = 1: compute G0_rel with relel_dudtau value directly to conserve energy
  // reltype = 0: compute G0_rel with integral over gamma dots

  struct geometry *geom
    = (struct geometry*) ggg;
  struct struct_of_state state;
  fill_struct_of_state(pp,geom,&state);

  ldouble out = calc_Gi_with_state(pp,&state,geom,Gi,relel_dudtau,type,reltype);

  return out;
}

ldouble
calc_Gi_with_state(ldouble *pp, void *sss, void *ggg, ldouble Gi[4], ldouble relel_dudtau, int type, int reltype)
{
  // type = 0: thermal + nonthermal in fluid frame
  // type = 1: thermal + nonthermal in lab frame
  // type = 2: thermal ONLY in fluid frame
  // type = 3: thermal ONLY in lab frame
  
  // reltype = 1: compute G0_rel with relel_dudtau value directly to conserve energy
  // reltype = 0: compute G0_rel with integral over gamma dots

#ifdef NONRELMHD
  return calc_Gi_nonrel_with_state(pp, sss, ggg, Gi, type);
#endif
  
  ldouble Gi0[4], Gi1[4], Gi2[4], Gi3[4];
  ldouble out;
  
  out = calc_all_Gi_with_state(pp, sss, ggg, Gi0, Gi1, Gi2, Gi3, relel_dudtau, reltype);
  if (type==0) {
    Gi[0] = Gi0[0];    
    Gi[1] = Gi0[1];    
    Gi[2] = Gi0[2];    
    Gi[3] = Gi0[3];    
  } else if (type==1) {
    Gi[0] = Gi1[0];    
    Gi[1] = Gi1[1];    
    Gi[2] = Gi1[2];    
    Gi[3] = Gi1[3]; 
  } else if (type==2) {
    Gi[0] = Gi2[0];    
    Gi[1] = Gi2[1];    
    Gi[2] = Gi2[2];    
    Gi[3] = Gi2[3]; 
  } else {
    Gi[0] = Gi3[0];    
    Gi[1] = Gi3[1];    
    Gi[2] = Gi3[2];    
    Gi[3] = Gi3[3]; 
  }
 
  return out;
}

ldouble
calc_all_Gi_with_state(ldouble *pp, void *sss, void* ggg, ldouble Gitot_ff[4], ldouble Gitot_lab[4], ldouble Gith_ff[4], ldouble Gith_lab[4], ldouble relel_dudtau, int reltype)
{
  // reltype = 1: compute G0_rel with relel_dudtau value directly to conserve energy
  // reltype = 0: compute G0_rel with integral over gamma dots
  
  int i,j,k;
  struct geometry *geom
  = (struct geometry *) ggg;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;
  
  ldouble (*gg)[5], (*GG)[5];
  gg = geom->gg;
  GG = geom->GG;
  
  //radiative stress tensor in the lab frame
  ldouble Rij[4][4];
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      Rij[i][j] = state->Rij[i][j];
    }
  }
  
  //the contrav. four-velocity of fluid in fluid frame
  ldouble ucon_ff[4];
  ucon_ff[1] = ucon_ff[2] = ucon_ff[3] = 0.;
  ucon_ff[0] = 1.;
  
  //the four-velocity of fluid in lab frame
  ldouble ucon[4],  ucov[4];
  for (i = 0; i < 4; i++)
  {
    ucon[i] = state->ucon[i];
    ucov[i] = state->ucov[i];
  }
  
  //gas properties
  ldouble Tgas, Te, Ti;
  Tgas = state->Tgas;
  Te = state->Te;
  ldouble B = sigma_rad_over_pi * Te * Te * Te * Te;
  
  //struct opacities opac;
  ldouble kappa, kappaGasAbs, kappaRadAbs, kappaGasRoss, kappaRadRoss;
  kappaGasAbs = state->opac.kappaGasAbs;
  kappaRadAbs = state->opac.kappaRadAbs;
  kappaGasRoss = state->opac.kappaGasRoss;
  kappaRadRoss = state->opac.kappaRadRoss;
  ldouble kappaes = state->kappaes;
  
  //contravariant four-force in the lab frame
  ldouble Ehatrad = state->Ehat;
  ldouble Ruu = Ehatrad;
   
  //Calculate G thermal in lab frame
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
  for(i = 0; i < 4; i++)
  {
    ldouble Ru = 0.;
    for(j = 0; j < 4; j++)
      Ru += Rij[i][j] * ucov[j];
    Gith_lab[i] = -(kappaRadRoss + kappaes)*Ru - ((kappaRadRoss + kappaes - kappaRadAbs) * Ruu + kappaGasAbs * fourpi * B) * ucon[i];
  }
  
  // Boost to fluid frame
  boost2_lab2ff_4vel_with_alpha(Gith_lab, Gith_ff, pp, gg, GG, geom->alpha, ucon, ucov);
  
  //rewrite the time component directly
  Gith_ff[0] = -kappaGasAbs * fourpi * B + kappaRadAbs * Ehatrad;
  
  //Comptonization
#if defined(COMPTONIZATION) || defined(EVOLVEPHOTONNUMBER)
  ldouble Gic_ff[4];
  ldouble Gic_lab[4];
  calc_Compt_Gi_with_state(pp, state, ggg, Gic_lab, ucon);
  calc_Compt_Gi_with_state(pp, state, ggg,Gic_ff, ucon_ff);
  
  ldouble fac=1.;
#ifdef DAMPCOMPTONIZATIONATBH
  ldouble xxBL[4];
  coco_N(geom->xxvec,xxBL,MYCOORDS,BLCOORDS);
  fac=step_function(xxBL[1]-1.2*rhorizonBL,.1*rhorizonBL);
  if(xxBL[1]<rhorizonBL) fac=0.;
#endif
  
  for(i = 0; i < 4; i++)
  {
    Gith_lab[i] += fac * Gic_lab[i];
    Gith_ff[i] += fac * Gic_ff[i];
  }
#endif
  
  // Set total Gi equal to thermal
  for(i = 0; i < 4; i++)
  {
    Gitot_lab[i] = Gith_lab[i];
    Gitot_ff[i] = Gith_ff[i];
  }
  
  // Add the component from nonthermal electrons
#ifdef RELELECTRONS
  ldouble G0_rel_ff=0.;
  G0_rel_ff = calc_relel_G0_fluidframe_from_state(pp, state, geom, relel_dudtau, reltype);
  
  // lab frame
  Gitot_lab[0] += G0_rel_ff * ucon[0];
  Gitot_lab[1] += G0_rel_ff * ucon[1];
  Gitot_lab[2] += G0_rel_ff * ucon[2];
  Gitot_lab[3] += G0_rel_ff * ucon[3];
  
  // fluid frame
  Gitot_ff[0] += G0_rel_ff;
#endif
  
  return 1. / (Tgas * Tgas * kappa);
}


//*****************************************************************************
//****** takes radiative stress tensor and gas primitives *********************
//****** and calculates nonrel contravariant four-force as in Jiang+14 ********
//*****************************************************************************
ldouble
calc_Gi_nonrel_with_state(ldouble *pp, void *sss, void *ggg, ldouble Gi[4], int labframe)  // needs to be worked on
{
  int i,j,k;
  struct geometry *geom
  = (struct geometry *) ggg;
  
  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  ldouble gamma=GAMMA;
#ifdef CONSISTENTGAMMA
  gamma=pick_gammagas(geom->ix,geom->iy,geom->iz);
#endif
  ldouble gammam1=gamma-1.;
  
  //radiative stress tensor in the lab frame
  ldouble Rij[4][4];
  calc_Rij_M1(pp,ggg,Rij);
  
  //the four-velocity of fluid in lab frame
  ldouble ucon[4],utcon[4],ucov[4],vpr[3],vgas[3];
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);
  vgas[0]=utcon[1];
  vgas[1]=utcon[2];
  vgas[2]=utcon[3];
  
  //gas properties
  ldouble rho=pp[RHO];
  ldouble u=pp[1];
  ldouble p= (gamma-1.)*(ldouble)u;
  ldouble Ti,Te;
  ldouble Tgas=calc_PEQ_Teifrompp(pp,&Te,&Ti,geom->ix,geom->iy,geom->iz);
  ldouble B = SIGMA_RAD*pow(Te,4.)/Pi;
  struct opacities opac;
  ldouble kappaGasAbs,kappaRadAbs,kappaGasNum,kappaRadNum,kappaGasRoss,kappaRadRoss;
  ldouble kappa=calc_kappa(pp,geom,&opac);
  ldouble kappaes=calc_kappaes(pp,geom);
  kappaGasAbs=opac.kappaGasAbs;
  kappaRadAbs=opac.kappaRadAbs;
  kappaGasNum=opac.kappaGasNum;
  kappaRadNum=opac.kappaRadNum;
  kappaGasRoss=opac.kappaGasRoss;
  kappaRadRoss=opac.kappaRadRoss;
  
  //radiation energy and fluxes
  ldouble Er,Fr[3],Pr[3][3];
  Er=Rij[0][0];
  Fr[0]=Rij[1][0];
  Fr[1]=Rij[2][0];
  Fr[2]=Rij[3][0];
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      Pr[i-1][j-1]=Rij[i][j];
  
  //printf("nr %e %e %e\n",kappa,Er, B);
  
  //fluid frame: Gi[0]=-kappagasAbs*4.*Pi*B + kappaRadAbs*Ehatrad;
  Gi[0]=-kappa*(4.*Pi*B - Er) - (kappa - kappaes) *
  (
   vgas[0]*(Fr[0]-(vgas[0]*Er + vgas[0]*Pr[0][0] + vgas[1]*Pr[0][1] + vgas[2]*Pr[0][2])) +
   vgas[1]*(Fr[1]-(vgas[1]*Er + vgas[0]*Pr[1][0] + vgas[1]*Pr[1][1] + vgas[2]*Pr[1][2])) +
   vgas[2]*(Fr[2]-(vgas[2]*Er + vgas[0]*Pr[2][0] + vgas[1]*Pr[2][1] + vgas[2]*Pr[2][2]))
   );
  
  for(i=0;i<3;i++)
  {
    Gi[i+1]=
    (kappa + kappaes)*(Fr[i]-(vgas[i]*Er + vgas[0]*Pr[i][0] + vgas[1]*Pr[i][1] + vgas[2]*Pr[i][2])) + vgas[i]*kappa*(4.*Pi*B - Er);
  }
  
  if(labframe==0)
  {
    //this is approximate and inconsistent with Jiang - to be used as a backup only!
    //R^ab u_a u_b = Erad in fluid frame
    ldouble Ruu=0.;
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        Ruu+=Rij[i][j]*ucov[i]*ucov[j];
    ldouble Ehatrad = Ruu;
    boost2_lab2ff_with_alpha(Gi, Gi, pp, gg, GG, geom->alpha);
    
    //rewrite the time component directly
    Gi[0]=-kappaGasAbs*4.*Pi*B + kappaRadAbs*Ehatrad;
  }
  
  return 1./(Tgas*Tgas*kappa);
}


/************************************************************/
/***** Thermal Comptonization component of the **************/
/***** energy transfer through G^t **************************/
/************************************************************/
int
calc_Compt_Gi(ldouble *pp, void* ggg, ldouble *Gic, ldouble Ehatrad, ldouble Te, ldouble kappaes, ldouble *ucon)
{
  
  struct geometry *geom
   = (struct geometry *) ggg;
  struct struct_of_state state;

  fill_struct_of_state(pp,geom,&state);

  int out = calc_Compt_Gi_with_state(pp,&state,geom,Gic,ucon);
  return out;
}

/************************************************************/
/***** Thermal Comptonization component of the **************/
/***** energy transfer through G^t **************************/
/************************************************************/
int
calc_Compt_Gi_with_state(ldouble *pp, void *sss, void* ggg, ldouble *Gic, ldouble *ucon_frame)
{
#ifdef RADIATION
  int i;
  
  struct geometry *geom
  = (struct geometry *) ggg;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;
  
  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  ldouble Ehatrad = state->Ehat;
  ldouble Te = state->Te;
  ldouble kappaes = state->kappaes;
  ldouble ThatradBB = state->TradBB;
  ldouble Thatrad = state->Trad;
  
  ldouble urfcon[4], uffcov[4];
  for (i = 0; i < 4; i++)
  {
    urfcon[i] = state->urfcon[i];
    uffcov[i] = state->ucov[i];
  }

  //ANDREW correction factor for the nonthermal electrons present in kappa_es
  ldouble relel_corr = 1.0;
#ifdef RELELECTRONS
  relel_corr = (state->ne)/((state->rho)*one_over_mue_mp);
#endif
  
  ldouble thetae = kB_over_me * Te;
  ldouble coeff = relel_corr * kappaes * Ehatrad * (4. * kB_over_me * (Thatrad - Te)) * (1. + 3.683 * thetae + 4. * thetae * thetae) / (1. + 4. * thetae);
  for(i=0;i<4;i++)
  {
    Gic[i]= coeff * ucon_frame[i];
  }
  
#endif
  return 0;
  
}

//**********************************************************************
//Calculate Thermal Couloumb coupling rate
//**********************************************************************
ldouble
calc_CoulombCoupling(ldouble *pp,void *ggg)
{

  struct geometry *geom
    = (struct geometry *) ggg;
  struct struct_of_state state;

  fill_struct_of_state(pp,geom,&state);

  ldouble out = calc_CoulombCoupling_with_state(pp,&state,geom);

  return out;
}

ldouble
calc_CoulombCoupling_with_state(ldouble *pp,void *sss,void *ggg)
{

  struct geometry *geom
    = (struct geometry *) ggg;
  struct struct_of_state *state
  = (struct struct_of_state *) sss;

  ldouble CoulombC=0.;
  ldouble rho=pp[RHO];
  ldouble Ti=state->Ti;
  ldouble Te=state->Te;
  ldouble rhocgs = rhoGU2CGS(rho);

  ldouble k_boltCGS = K_BOLTZ_CGS;
  ldouble me = M_ELECTR_CGS;
  ldouble mp = M_PROTON_CGS;

  //ANDREW -- is this the right mass ratio
  ldouble mratio=M_ELECTR/(M_PROTON*MU_I);

  ldouble cc =  CCC_CGS;
  ldouble sigmaT = SIGMATH_CGS;
  //BRANDON - n_e*n_avg for general metallicity: n_avg = (X + Y + <Z_j^2/A_j>*Z)*rho/mp
  ldouble Z2nine = 0.5*(1+HFRAC)*(HFRAC + HEFRAC + Z2divA_MEAN*MFRAC)*rhocgs*rhocgs/(mp*mp); 
  ldouble theta_e = kB_over_me*Te;
  ldouble theta_i = kB_over_mui_mp*Ti;
  ldouble theta_mean = theta_e*theta_i/(theta_e + theta_i);
  ldouble K0, K1, K2e, K2i;

  //bessel functions
  K2e = gsl_sf_bessel_Kn(2, 1./theta_e);
  K2i = gsl_sf_bessel_Kn(2, 1./theta_i);
  K0 = gsl_sf_bessel_Kn(0, 1./theta_mean);
  K1 = gsl_sf_bessel_Kn(1, 1./theta_mean);

  double theta_min = 1.e-2;
  if (theta_i < theta_min && theta_e < theta_min) //asymptotic expansions
   {
    CoulombC = 1.5*sigmaT*cc*k_boltCGS*mratio*Z2nine*(Ti -
						       Te)/sqrt(0.5*M_PI*(theta_e+theta_i))*20.*( (2.*(theta_e + theta_i)*(theta_e +
												  theta_i) + 1. + 2.*(theta_e + theta_i) )/(theta_e + theta_i));
   }
  else if (theta_i < theta_min)
   {
    CoulombC =  1.5*sigmaT*cc*k_boltCGS*mratio*Z2nine*(Ti -
							Te)/K2e*20./exp(1./theta_e)*sqrt(theta_mean/theta_i)*( (2.*(theta_e + theta_i)*(theta_e +
   													       theta_i) + 1. + 2.*(theta_e + theta_i) )/(theta_e + theta_i));
   }
  else if (theta_e < theta_min)
   {
    CoulombC =  1.5*sigmaT*cc*k_boltCGS*mratio*Z2nine*(Ti -
							Te)/K2i*20./exp(1./theta_i)*sqrt(theta_mean/theta_e)*( (2.*(theta_e + theta_i)*(theta_e +
   													       theta_i) + 1. + 2.*(theta_e + theta_i) )/(theta_e + theta_i));
   }
  else //general
    CoulombC =  1.5*sigmaT*cc*k_boltCGS*mratio*Z2nine*(Ti -
							Te)/K2e/K2i*20.*( (2.*(theta_e + theta_i)*(theta_e + theta_i) +
									   1.)/(theta_e + theta_i)*K1 + 2.*K0 ); //standard formula
 
 if(!isfinite(CoulombC)) CoulombC=0.;
    if(!isfinite(CoulombC) && 0)
    {
      printf("Coulomb coupling: ccc %e Se: %e\n",CoulombC,pp[ENTRE]);
      printf("### %e %e %e %e\n",theta_e,theta_i,theta_mean,-1.);
      printf(">>> %e %e %e %e\n",K2e,K2i,K0,K1);
      printf("ooo %e %e %e\n",Te,Ti,Te-Ti);
      
      ldouble rho=pp[RHO];
      ldouble S3=pp[ENTRE];      
      ldouble mu=MU_E;
      ldouble mass=M_ELECTR;
      ldouble n=rho/mu/M_PROTON;
      ldouble ex=exp(S3/n/K_BOLTZ);
      ldouble rhs=cbrt(rho*ex*rho*ex);
      ldouble theta=1./5.*(-2.+sqrt(4.+25.*rhs));
      ldouble T =  theta*mass/K_BOLTZ;

      printf("exp: %e %e\n",S3/n/K_BOLTZ,ex);
      printf("TfromS3: %e %e %e %e %e %e\n",S3,rho,ex,rhs,theta,T);	
    }

    return kappacgs2gu*rhocgs2gu*endencgs2gu*velcgs2gu*CoulombC;
}


/****************************************************/
/****** Calculate Ehatrad in fluid frame ************/
/****************************************************/
void
calc_Ehat_from_Rij_ucov(double Rij[4][4], double uffcov[4], ldouble *Ehat)
{
  int i, j;
  
  //R^ab u_a u_b = Erad in the fluid frame
  
  *Ehat = 0.;
  for(i=0;i<4;i++)
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
    for(j=0;j<4;j++)
      *Ehat+=Rij[i][j]*uffcov[i]*uffcov[j];
  
  return;
}

/****************************************************/
/***** radiative stress energy tensor ***************/
/***** M1 (or other closure) only here **************/
/****************************************************/
int
calc_Rij(ldouble *pp, void* ggg, ldouble Rij[][4])
{
#ifdef RADIATION
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  int verbose=0;
  int i,j;
  
  calc_Rij_M1(pp,ggg,Rij); // M1 closure
  
#endif //RADIATION
  return 0;
}

//**********************************************************************
//******* takes E and F^i from primitives (artificial) *****************
//******* and calculates radiation stress ******************************
//******* tensor R^ij in fluid frame using M1 closure scheme ***********
//**********************************************************************
int
calc_Rij_M1_ff(ldouble *pp, ldouble Rij[][4])
{
  int irf=0;
  ldouble E=pp[EE];
  ldouble F[3]={pp[FX],pp[FY],pp[FZ]};

  ldouble nx,ny,nz,nlen,f;

  nx=F[0]/E;
  ny=F[1]/E;
  nz=F[2]/E;

  nlen=sqrt(nx*nx+ny*ny+nz*nz);
  
  if(nlen>=1.)
    f=1.;	
  else
    f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  
  
  if(nlen>0) 
    {
      nx/=nlen;
      ny/=nlen;
      nz/=nlen;
    }
 
  Rij[0][0]=E;
  Rij[0][1]=Rij[1][0]=F[0];
  Rij[0][2]=Rij[2][0]=F[1];
  Rij[0][3]=Rij[3][0]=F[2];

  Rij[1][1]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nx*nx);
  Rij[1][2]=E*(.5*(3.*f - 1.)*nx*ny);
  Rij[1][3]=E*(.5*(3.*f - 1.)*nx*nz);

  Rij[2][1]=E*(.5*(3.*f - 1.)*ny*nx);
  Rij[2][2]=E*(.5*(1.-f) + .5*(3.*f - 1.)*ny*ny);
  Rij[2][3]=E*(.5*(3.*f - 1.)*ny*nz);

  Rij[3][1]=E*(.5*(3.*f - 1.)*nz*nx);
  Rij[3][2]=E*(.5*(3.*f - 1.)*nz*ny);
  Rij[3][3]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nz*nz);

  return 0;
}

/****************************************************/
/***** radiative stress energy tensor ***************/
/***** from radiation four-velocity *****************/
/***** pure M1 only *********************************/
/****************************************************/

int
calc_Rij_M1_from_4vel(ldouble *pp, void* ggg, ldouble *urfcon, ldouble Rij[][4])
{
#ifdef RADIATION
  struct geometry *geom
  = (struct geometry *) ggg;
  
  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  int verbose=0;
  int i,j;
  
  //radiative energy density in the radiation rest frame
  ldouble Erf = pp[EE0]; 

  //lab frame stress energy tensor:
  for(i = 0; i < 4; i++)
  {
    for(j = 0; j < 4; j++)
    {
      Rij[i][j] = four_third * Erf * urfcon[i] * urfcon[j] + one_third * Erf * GG[i][j];
    }
  }
#endif
  
  return 0;
}


/****************************************************/
/***** radiative stress energy tensor ***************/
/***** pure M1 only here ****************************/
/****************************************************/

int
calc_Rij_M1(ldouble *pp, void* ggg, ldouble Rij[][4])
{
#ifdef RADIATION
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  int verbose=0;
  int i,j;

  //covariant formulation of M1
  ldouble Erf; //radiative energy density in the radiation rest frame
  ldouble urfcon[4];
  Erf=pp[EE0];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];

  //converting to lab four-velocity
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);

  //lab frame stress energy tensor:
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Rij[i][j] = four_third * Erf * urfcon[i] * urfcon[j] + one_third * Erf * GG[i][j];

#endif
  return 0;
}

//*************************************************
//suplementary routines for radiation conversions
//*************************************************

// calculates LTE temperature from photon number 
ldouble
calc_Tnfromn(ldouble n)
{
  return cbrt(2.70118*K_BOLTZ*n/4./SIGMA_RAD);
}

//calculates LTE number of photons from temp
ldouble
calc_NFfromT(ldouble T)
{
  return A_RAD*T*T*T/2.70118/K_BOLTZ;
}

//calculates LTE number of photons from energy
ldouble
calc_NFfromE(ldouble E)
{
  ldouble temp=calc_LTE_TfromE(E);
  return calc_NFfromT(temp);
}

//calculates LTE energy from temperature
ldouble
calc_LTE_EfromT(ldouble T)
{
  return four_sigmarad * T * T * T * T;
}

//calculates LTE temeprature from energy
ldouble
calc_LTE_TfromE(ldouble E )
{
  return sqrt(sqrt((one_over_four_sigmarad * E)));
}

//radiation energy corresponding to gas temperature corresponding with u and rho
//inconsistent with CONSISTENTGAMMA !!
ldouble
calc_LTE_Efromurho(ldouble u,ldouble rho)
{
  ldouble p=(GAMMA-1.)*(u); 
  ldouble T=p*MU_GAS*M_PROTON/K_BOLTZ/rho;

  return calc_LTE_EfromT(T);
}

//calculates fluid frame radiative energy density (-R^t_t = Ehat)
//and lab-frame four-velocity of gas
int
calc_ff_Rtt(ldouble *pp,ldouble *Rttret, ldouble* ucon,void* ggg)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble ucov[4],utcon[4];
  utcon[0]=0.;
  utcon[1]=pp[VX];
  utcon[2]=pp[VY];
  utcon[3]=pp[VZ];

  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,geom->gg,geom->GG);
  
  ldouble Rij[4][4],Rtt;
  calc_Rij_M1(pp,ggg,Rij);
 
  indices_2221(Rij,Rij,geom->gg);
  Rtt=0.;
  int i1,i2;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
      Rtt+=-Rij[i1][i2]*ucon[i2]*ucov[i1];

  *Rttret = Rtt;

  return 0;
}

//calculates fluid frame radiative energy density
int
calc_ff_Ehat(ldouble *pp,ldouble *Ehat, ldouble* ucon,void* ggg)
{
  calc_ff_Rtt(pp,Ehat,ucon,ggg);
  *Ehat *= -1.;
  return 0;
}

//**********************************************************************
//**********************************************************************
// Photon number routines
//**********************************************************************
//**********************************************************************

//**********************************************************************
//* calculates the source term for number of photons
//**********************************************************************
ldouble
calc_nsource(ldouble *pp, void* ggg)
{

  struct geometry *geom
    = (struct geometry *) ggg;
  struct struct_of_state state;
  fill_struct_of_state(pp,geom,&state);

  ldouble out = calc_nsource_with_state(pp,&state,geom);
  return out;
}

ldouble
calc_nsource_with_state(ldouble *pp, void *sss, void* ggg)
{
  int i,j, i1;
  ldouble uffcon[4], uffcov[4], bsq, urfcon[4], urfcov[4],Ehatrad;
  
#ifdef EVOLVEPHOTONNUMBER
  
  struct struct_of_state *state
  = (struct struct_of_state *) sss;
  
  struct geometry *geom
  = (struct geometry *) ggg;
  
  int gix = geom->ix + TOI;
  int giy = geom->iy + TOJ;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;
  gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  
  for (i1 = 0; i1 < 4; i1++)
  {
    uffcon[i1] = state->ucon[i1];
    uffcov[i1] = state->ucov[i1];
    urfcon[i1] = state->urfcon[i1];
    urfcov[i1] = state->urfcov[i1];
  }

  ldouble Thatrad = state->Trad;
  
  //gas properties & opacities
  ldouble Tgas = state->Tgas;
  ldouble Te = state->Te;
  ldouble Ti = state->Ti;
  ldouble TradBB = state->TradBB;
  ldouble Trad=Thatrad; 
  ldouble B = sigma_rad_over_pi * Te * Te * Te * Te;
  ldouble Tc_n=5.07783e9; //cross over temperature for non-relativistic to ultra-relativistic number opacity

  ldouble Trad_lim = pow(TradBB,1.333333333333)/(pow(Te,0.333333333333));
  ldouble nph_lim = pp[NF]*(Trad/Trad_lim);

  #ifdef USE_SYNCHROTRON_BRIDGE_FUNCTIONS // TradBB should not exceed Te BB
  if(Trad < Trad_lim)
  {
    Trad=Trad_lim;
  }
  #endif

  struct opacities opac;
  ldouble kappaGasAbs,kappaRadAbs,kappaGasNum,kappaRadNum,kappaGasRoss,kappaRadRoss;
  kappaGasAbs = state->opac.kappaGasAbs;
  kappaRadAbs = state->opac.kappaRadAbs;
  kappaGasNum = state->opac.kappaGasNum;
  kappaRadNum = state->opac.kappaRadNum;
  kappaGasRoss = state->opac.kappaGasRoss;
  kappaRadRoss = state->opac.kappaRadRoss;
  
  //Thermal Synchrotron
  ldouble ndotffSynch=0.;
  
#ifdef SYNCHROTRON
  bsq = state->bsq;
  ldouble Bmagcgs=0.;
  
  //calculate B field strength
#ifdef MAGNFIELD
  ldouble bsqcgs = fourpi * endengu2cgs * bsq;
  Bmagcgs = sqrt(bsqcgs);
#endif
  
  ldouble ne = numdensgu2cgs * state->ne; //ANDREW thermal number density in CGS
  ldouble ndotffSynchCgs = 1.44e5 * Bmagcgs * ne;
  ndotffSynch = ndotffSynchCgs * per_volume_per_time_cgs2gu;
  
  #ifdef USE_SYNCHROTRON_BRIDGE_FUNCTIONS
  ndotffSynch *= (Te/Tc_n)/(1.+(Te/Tc_n)); //Bridge Function: makes ndot -> non-relativistic solution for low temperatures
  #else
  // Ramesh: suppress photon generation rate at nonrelativistic temperatures
  // avoids numerical problems at low temperatures
  // Replaced with above bridge function
  ldouble Terel = Te * k_over_mecsq;  // this is kT/mec^2
  ldouble Terelfactor = (Terel * Terel) / (1. + Terel * Terel);  // suppression factor
  ndotffSynch *= Terelfactor;
  #endif
  
  //estimate the synchrotron radiation temperature produced
  ldouble TradSynch = 2.05e-18 * Bmagcgs * Te * Te * 1.e-20 *k_boltz_cgs_inv / 2.70188; //Ramesh's fit
  ldouble Tratio=TradSynch/Te;
  
  //if too low/too high against electron temp - modify the photon generation rate
  ldouble allowance=SYNCHROTRONALLOWEDTRATIO;  
  if(Tratio<1./allowance)
    ndotffSynch*=Tratio*allowance;
  else if(Tratio>allowance)
    ndotffSynch*=Tratio/allowance;
  
#endif //SYNCHROTRON
  
  //Contribution to photon production from relativistic free-free & synchrotron
  ldouble ndot_relel_syn=0.0;
  ldouble ndot_relel_ff=0.0;
  
#ifdef RELELECTRONS
  
  ndot_relel_syn = calc_relel_photon_ndot_from_state(pp, state, 1); //Synchrotron
  ndot_relel_ff = calc_relel_photon_ndot_from_state(pp, state, 2); //Free-Free
  
#ifdef RELEL_SYN
  if (ndot_relel_syn != 0.0 && ndotffSynch!=0.0)
  {
    //Nonthermal synchrotron temperature
    ldouble TradSynch_relel = fabs(calc_relel_G0_fluidframe_direct_from_state(pp, state, 1)) / (ndot_relel_syn * K_BOLTZ);
    ldouble Tratio_relel=TradSynch_relel/Trad; 
    
    //Synchrotron adjustment for large temperature ratios.
    ldouble allowance_relel = SYNCHROTRONALLOWEDTRATIO_RELEL;
    if(Tratio_relel<1./allowance_relel)
      ndot_relel_syn*=Tratio_relel*allowance_relel;
    else if(Tratio_relel>allowance_relel)
      ndot_relel_syn*=Tratio_relel/allowance_relel;
  }
#endif //RELEL_SYN
#endif //RELELECTRONS
  
  //remaining contributions to the thermal photon production rate are from opacity
  Ehatrad=state->Ehat;
  ldouble ndotff_thermal = (0.370209) * k_boltz_inv *(kappaGasNum* fourmpi *B/Te - kappaRadNum * Ehatrad / Thatrad);
  
  //the rate of change of number of photons is invariant
  ldouble ndotff = ndotff_thermal + ndotffSynch + ndot_relel_syn + ndot_relel_ff;
  ldouble ndotrf = ndotff;
  
  return ndotrf;

#else //EVOLVEPHOTONNUMBER not defined
  return 0;
#endif  // EVOLVEPHOTONNUMBER
}


//**********************************************************************
//* calculates rad. en. den. in rad. rest frame
//*  based on temperature and comoving(!) number of photons
//**********************************************************************
ldouble
calc_ncompt_Ehatrad(ldouble Tradhat, ldouble nphhat)
{ 
  //diluted BB
  ldouble Ehatrad=Tradhat*nphhat*2.70118*K_BOLTZ; 
  return Ehatrad;
}

ldouble
calc_ncompt_Thatrad_fromEN(ldouble Ehat, ldouble nphhat)
{
  //diluted BB
  ldouble Tradhat = Ehat/(nphhat*2.70118*K_BOLTZ);
  return Tradhat;
}


//**********************************************************************
//* calculates radiation color temperature in the fluid frame
//* and calls more
//**********************************************************************
ldouble
calc_ncompt_Thatrad_full(ldouble *pp, void* ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  //radiative stress tensor in the lab frame
  ldouble Rij[4][4];
  calc_Rij(pp,ggg,Rij);

  //the four-velocity of fluid 
  ldouble uffcon[4],utcon[4],uffcov[4],vpr[3];
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  conv_vels_both(utcon,uffcon,uffcov,VELPRIM,VEL4,gg,GG);
  
  //R^ab u_a u_b = Erad in fluid frame
  ldouble Ruu=0.;
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Ruu+=Rij[i][j]*uffcov[i]*uffcov[j];
  ldouble Ehatrad = Ruu;
  
  //the four-velocity of radiation in lab frame
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom->gg,geom->GG);

  return calc_ncompt_Thatrad_4vel(pp,ggg,Ehatrad,urfcon,uffcov);
}

//**********************************************************************
//* calculates the radiation temperature using
//* radiative energy density and the number of photons
//**********************************************************************
ldouble
calc_ncompt_Thatrad(ldouble *pp, void* ggg, ldouble Ehatrad)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  //the four-velocity of fluid in lab frame
  ldouble uffcon[4],utcon[4],uffcov[4],vpr[3];
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  conv_vels_both(utcon,uffcon,uffcov,VELPRIM,VEL4,gg,GG);

  //the four-velocity of radiation in lab frame
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom->gg,geom->GG);

  return calc_ncompt_Thatrad_4vel(pp,ggg,Ehatrad,urfcon,uffcov);
}

//**********************************************************************
//* calculates the radiation color temperature using
//* radiative energy density and the number of photons
//* assumes four-velocities known
//**********************************************************************
ldouble
calc_ncompt_Thatrad_4vel(ldouble *pp, void* ggg, ldouble Ehatrad, ldouble* urfcon, ldouble* uffcov)
{
  //number of photons in rad rest frame
  ldouble nphrf = pp[NF];
  
  //relative gamma rad-fluid rest frames
  ldouble relgamma = urfcon[0]*uffcov[0] + urfcon[1]*uffcov[1] +urfcon[2]*uffcov[2] +urfcon[3]*uffcov[3];
  
  //ff quantities
  ldouble nphhat = - relgamma*nphrf;
    
  //dilluted BB
  ldouble Thatrad = one_over_kB_2p70118 * Ehatrad / nphhat;
    
  return Thatrad;
}


//**********************************************************************
//* calculates the radiation color temperature using
//* radiative energy density and the number of photons
//* assumes relgamma is known
//**********************************************************************
ldouble
calc_ncompt_Thatrad_nphhat(ldouble nphhat, ldouble Ehatrad)
{
  //diluted BB
  ldouble Thatrad = one_over_kB_2p70118 * Ehatrad / nphhat;

  return Thatrad;
}


//**********************************************************************
//* calculates the number of photons in the lab frame
//**********************************************************************
ldouble
calc_ncompt_nphlab(ldouble *pp, void* ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  //four velocity of the lab frame
  ldouble ucov[4], ucon[4];
  calc_normalobs_ncov_ncon(GG, geom->alpha, ucov, ucon);

  //the four-velocity of radiation in lab frame
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom->gg,geom->GG);

  //number of photons in rad rest frame
  ldouble nphrf = pp[NF];

  //relative gamma rad:-lab frames
  ldouble relgamma = urfcon[0]*ucov[0] + urfcon[1]*ucov[1] +urfcon[2]*ucov[2] +urfcon[3]*ucov[3]; 

  //lab
  ldouble nphlab = - relgamma*nphrf;

  return nphlab;
}


//**********************************************************************
//* calculates the number of photons in the fluid frame
//**********************************************************************
ldouble
calc_ncompt_nphhat(ldouble *pp, void* ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  //the four-velocity of fluid in lab frame
  ldouble uffcon[4],utcon[4],uffcov[4],vpr[3];
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  conv_vels_both(utcon,uffcon,uffcov,VELPRIM,VEL4,gg,GG);

  //the four-velocity of radiation in lab frame
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom->gg,geom->GG);

  //number of photons in rad rest frame
  ldouble nphrf = pp[NF];

  //relative gamma rad-fluid rest frames
  ldouble relgamma = urfcon[0]*uffcov[0] + urfcon[1]*uffcov[1] +urfcon[2]*uffcov[2] +urfcon[3]*uffcov[3]; 

  //ff quantities
  ldouble nphhat = - relgamma*nphrf;

  return nphhat;
}

//******************************************************************************/
//******* calculates wavespeeds in the lab frame taking 1/sqrt(3) in ***********/
//******* radiative rest frame and boosting it to lab frame ********************/
//******* using the HARM algorithm - with taul limiter *************************/
//******* or calculating numerically at cell centers ***************************/
//******************************************************************************/
// July 8, 17, Ramesh: Modified to make use of calc_wavespeeds_lr_core_new

int
calc_rad_wavespeeds(ldouble *pp,void *ggg,ldouble tautot[3],ldouble *aval,int verbose)
{
  verbose=0;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;
  gdetu=gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif

  int i,j;
  
  //transforming 1/sqrt(3) from radiation rest frame + limiter based on optical depth
  //the regular, non viscous, approach
  ldouble urfcon[4]; 
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);

  //square of radiative wavespeed in radiative rest frame
  ldouble rv2rad = one_third;
  ldouble rv2, rv2x, rv2y, rv2z, rv2tau;
  ldouble rv2dim[3];

  int dim, ret;
  ldouble cst1, cst2;
  ldouble aret[6];

  //first unlimited rad wavespeeds used for time step determination
  rv2 = one_third;
  ret = calc_wavespeeds_lr_core_new(urfcon,GG,aret,rv2,rv2,rv2);
  
  for (i = 0; i < 6; i++)
  {
    aval[i]=aret[i];
  }

  //now damped by the optical depth
  //characterisitic limiter based on the optical depth (Sadowski+13a)
#if defined(SKIPRADWAVESPEEDLIMITER) //radiative wavespeeds as in optically thin, full diffusion
  rv2x=rv2rad;
  rv2y=rv2rad;
  rv2z=rv2rad;
#elif defined(FULLRADWAVESPEEDS) //super-diffusive
  rv2x=1.;
  rv2y=1.;
  rv2z=1.;
#else
  for (dim = 0; dim < 3; dim++)
  {
    if(tautot[dim]<=0.)
      rv2dim[dim]=rv2rad;
    else //damp by optical depth
    {
#if defined(DAMPRADWAVESPEEDSQRTTAU)
      rv2tau=four_third/tautot[dim];
#elif defined(MULTIPLYRADWAVESPEEDCONSTANT)
      rv2tau=((four_third*four_third)/(tautot[dim]*tautot[dim]))*MULTIPLYRADWAVESPEEDCONSTANT*MULTIPLYRADWAVESPEEDCONSTANT;
#elif defined(DAMPRADWAVESPEEDLIMITERINSIDE)
      ldouble xxBL[4];
      coco_N(geom->xxvec,xxBL,MYCOORDS,BLCOORDS);
      ldouble radius=xxBL[1];
      ldouble damp=step_function(radius-DAMPRADWAVESPEEDLIMITERRADIUS,DAMPRADWAVESPEEDLIMITERFRAC*DAMPRADWAVESPEEDLIMITERRADIUS);
      rv2tau=(four_third*four_third)/(damp*tautot[dim]*tautot[dim]);
      
      //to limit the maximal wavespeed to cs
      ldouble temp[4],Ehat,rho,cs2tot;
      rho=pp[RHO];
      calc_ff_Ehat(pp,&Ehat,temp,ggg);
      cs2tot=one_third*Ehat/(rho*my_min(1.,tautot[dim]));
      if (rv2tau>cs2tot) rv2tau=cs2tot;
      
#else
      rv2tau=(four_third*four_third)/(tautot[dim]*tautot[dim]);
#endif // DAMPRADWAVESPEED
      
      rv2dim[dim]=my_min(rv2rad,rv2tau);
      
#ifdef DAMPRADWAVESPEEDNEARAXIS //increases numerical diffusion near axis to promote stability
      int giy;
      giy=geom->iy+TOJ;
      int limit = DAMPRADWAVESPEEDNEARAXISNCELLS;
      if(giy<limit || giy> (TNY-limit-1))
        rv2dim[dim]=rv2rad;
#endif
    }
  }
#endif  // SKIPRADWAVESPEEDLIMITER etc.
  
  rv2x = rv2dim[0];
  rv2y = rv2dim[1];
  rv2z = rv2dim[1];
  
  ret = calc_wavespeeds_lr_core_new(urfcon,GG,aret,rv2x,rv2y,rv2z);
  
  for (i = 0; i < 6; i++)
  {
    aval[i+6]=aret[i];
  }

  return 0;
}

//****************************************************
// calculates radiative tensor at faces or centers
// returns total, pure and visc ~ R^i_j
//****************************************************
int
f_flux_prime_rad_total(ldouble *pp, void *ggg,ldouble Rij[][4],ldouble Rij0[][4], ldouble Rijvisc[][4])
{    

#ifdef RADIATION
  int i,j;

  ldouble uu[NV];
  p2u(pp,uu,ggg);

  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble Rvisc1[4][4],Rvisc2[4][4];
  int ix,iy,iz; int iix,iiy,iiz; 
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  calc_Rij(pp,ggg,Rij0); //regular 0 R^ij
  indices_2221(Rij0,Rij0,gg); //R^i_j
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Rij[i][j]=Rij0[i][j];
	Rijvisc[i][j]=0.;
      }

  
#if (RADVISCOSITY==SHEARVISCOSITY)
  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;
  struct geometry geomcent;
      
  iix=ix;iiy=iy;iiz=iz;

  //when face and shear viscosity put the face primitives at both cell centers and average the viscous stress tensor
  if(geom->ifacedim>-1)
    //face fluxes
    {      
      //left (because ix as a face index corresponds to face between ix-1 and ix cells)
      if(geom->ifacedim==0)
	iix=ix-1;
      if(geom->ifacedim==1)
	iiy=iy-1;
      if(geom->ifacedim==2)
	iiz=iz-1;
      
      //using the face interpolated primitives

      int derdir[3]={0,0,0}; //by default centered derivatives in calc_shear
      
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  {
	    Rijvisc[i][j]=.5*(get_Tfull(Rijviscglobal,i,j,ix,iy,iz)+get_Tfull(Rijviscglobal,i,j,iix,iiy,iiz));
	  }

    }
  else
    //cell centered fluxes for char. wavespeed evaluation
    {
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  {
	    Rijvisc[i][j]=get_Tfull(Rijviscglobal,i,j,ix,iy,iz);
	  }     
    }  

  /**********************************/
  //damping the viscous term based on char. viscous velocity
  /**********************************/
#ifdef RADVISCMAXVELDAMP
  //basing on radiative Reynolds number Re = diffusiv flux of conserved quantity / conserved quantity
  ldouble dampfac=1.;
  int idim;
  ldouble vel[3]={0.,0.,0.},maxvel=-1.;
 
  //face centers, using given uu
  for(idim=0;idim<3;idim++)
    for(i=1;i<4;i++)
      {
	if(i==2 && NY==1) continue;
	if(i==3 && NZ==1) continue;
	if(fabs(uu[EE0+i])<1.e-10 * fabs(uu[EE0])) continue;
	vel[idim]=Rijvisc[idim+1][i]/(uu[EE0+i]/gdetu)*sqrt(gg[idim+1][idim+1]);
      }
  if(geom->ifacedim>-1)
    maxvel=fabs(vel[geom->ifacedim]);
  else
    maxvel=sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);

  //adjust:
  if(maxvel>MAXRADVISCVEL)
    {
      //todo: choose best prescription
      //dampfac= 1. / (1.+pow(maxvel/MAXRADVISCVEL,2.));
      //dampfac = 1.-step_function(-1.+maxvel/MAXRADVISCVEL,.5);    
      dampfac=MAXRADVISCVEL/maxvel;
      //printf("damping at %d (%e)\n",geom->ix,maxvel);
    }
  else
    dampfac=1.;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Rijvisc[i][j]*=dampfac;
      }
  
#endif
  //adding up to Rij
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Rij[i][j]+=Rijvisc[i][j];
      }

#endif //SHEARVISCOSITY
#endif //RADIATION

  return 0;
}


//***************************************
// calculates radiative fluxes at faces or centers
//***************************************

int f_flux_prime_rad( ldouble *pp, int idim, void *ggg,ldouble *ff)
{  
#ifdef RADIATION
  int i,j;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  ldouble Rij[4][4],Rij0[4][4],Rijvisc[4][4];
  f_flux_prime_rad_total(pp,ggg,Rij,Rij0,Rijvisc);

  //fluxes to ff[EE0+]

  ff[EE0]= gdetu*(Rij[idim+1][0]);
      
  ff[FX0]= gdetu*(Rij[idim+1][1]);
      
  ff[FY0]= gdetu*(Rij[idim+1][2]);
      
  ff[FZ0]= gdetu*(Rij[idim+1][3]);

#ifdef EVOLVEPHOTONNUMBER
  ldouble nphrad=pp[NF0];
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);
  ff[NF0]= gdetu*nphrad*urfcon[idim+1];
#endif

#endif
  return 0;
}


/*******************************************/
/* rad viscosity and shear at cell centers */
/*******************************************/
int
calc_rad_shearviscosity(ldouble *pp,void* ggg,ldouble shear[][4],ldouble *nuret,int *derdir)
{  
  
#if(RADVISCOSITY==SHEARVISCOSITY) //full shear tensor
  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble mfp,nu,mindx;

  if(geom->ix<-1) 
    {
      printf("rad_shear too far at %d %d\n",geom->ix,geom->iy);
      //getchar();
    }

  ldouble div=0.0;
  calc_shear_lab(pp,ggg,shear,&div,RAD,derdir);  
  indices_1122(shear,shear,geom->GG);
 
  //calculating the mean free path
  calc_rad_visccoeff(pp,ggg,&nu,&mfp,&mindx);
  
  *nuret=nu;

#else //no rad.viscosity
  int i1,i2;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
      shear[i1][i2]=0.;
  *nuret=0.;
#endif
  return 0;  
}

//***********************************************************************************
//calculates shear tensor sigma_ij in the lab frame at cell centers only!
//hdorrad == MHD -> using gas velocity
//hdorrad == RAD -> using radiative velocity
//derdir[] determines the type of derivative in each dimension (left,right,centered)
//************************************************************************************
int
calc_shear_lab(ldouble *pp0, void* ggg,ldouble S[][4],ldouble *div, int hdorrad,int *derdir)
{
  int i,j,k,iv;

  struct geometry *geomin
    = (struct geometry *) ggg;

  int ix,iy,iz;
  ix=geomin->ix;
  iy=geomin->iy;
  iz=geomin->iz;

  //lets get geometry again, just in case geomin is not in internal coordinates, like in postprocessing
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  ldouble (*gg)[5],(*GG)[5];
  gg=geom.gg;
  GG=geom.GG;

  //let's start with derivatives
  ldouble du[4][4]; //du_i,j
  ldouble du2[4][4]; //du^i,j

  int istart,whichvel;
  if(hdorrad==MHD)
    {
      whichvel=VELPRIM;
      istart=VX;
    }
  else if(hdorrad==RAD)
    {
      whichvel=VELPRIMRAD;
      istart=FX;
    }

  //instead:
  for(i=0;i<4;i++)
    {
      //force d/dt = 0 in shear
      du[i][0] = 0.;
      du2[i][0] = 0.;
    }

  ldouble ppm1[NV],ppp1[NV],pp[NV];
  ldouble ggm1[4][5],GGm1[4][5];
  ldouble ggp1[4][5],GGp1[4][5];
  ldouble xxvecm1[4],xxvec[4],xxvecp1[4];
  ldouble uconm1[4],uconp1[4],utconm1[4],utconp1[4],utcon[4],ucon[4];
  ldouble ucovm1[4],ucovp1[4],ucov[4];
  ldouble enl,enr;
  int idim;

  //four-velocity at cell basing on pp[]
  //xxvec=geom->xxvec;
  get_xx(ix,iy,iz,xxvec);
  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=pp0[iv];
    }
  utcon[1]=pp[istart];  utcon[2]=pp[istart+1];  utcon[3]=pp[istart+2];
  conv_vels_both(utcon,ucon,ucov,whichvel,VEL4,gg,GG);  

  //derivatives
  for(idim=1;idim<4;idim++)
    {
      if(idim==1)
	{
	  get_xx(ix-1,iy,iz,xxvecm1);
	  get_xx(ix+1,iy,iz,xxvecp1);
	  if(hdorrad==RAD) 
	    {
	      enl=get_u(u,EE0,ix-1,iy,iz);
	      enr=get_u(u,EE0,ix+1,iy,iz);
	    }
	  else
	    {
	      enl=get_u(u,UU,ix-1,iy,iz);
	      enr=get_u(u,UU,ix+1,iy,iz);
	    }
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
	  if(hdorrad==RAD) 
	    {
	      enl=get_u(u,EE0,ix,iy-1,iz);
	      enr=get_u(u,EE0,ix,iy+1,iz);
	    }
	  else
	    {
	      enl=get_u(u,UU,ix,iy-1,iz);
	      enr=get_u(u,UU,ix,iy+1,iz);
	    }
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
	 if(hdorrad==RAD) 
	    {
	      enl=get_u(u,EE0,ix,iy,iz-1);
	      enr=get_u(u,EE0,ix,iy,iz+1);
	    }
	  else
	    {
	      enl=get_u(u,UU,ix,iy,iz-1);
	      enr=get_u(u,UU,ix,iy,iz+1);
	    }
	 for(iv=0;iv<NV;iv++)
	   {
	     ppm1[iv]=get_u(p,iv,ix,iy,iz-1);
	     ppp1[iv]=get_u(p,iv,ix,iy,iz+1);
	   }
	 pick_g(ix,iy,iz-1,ggm1);  pick_G(ix,iy,iz-1,GGm1);
	 pick_g(ix,iy,iz+1,ggp1);  pick_G(ix,iy,iz+1,GGp1);
       }

     //calculating four velocity

     utconm1[1]=ppm1[istart];  utconm1[2]=ppm1[istart+1];  utconm1[3]=ppm1[istart+2];
     utconp1[1]=ppp1[istart];  utconp1[2]=ppp1[istart+1];  utconp1[3]=ppp1[istart+2];

     conv_vels_both(utconm1,uconm1,ucovm1,whichvel,VEL4,ggm1,GGm1);
     conv_vels_both(utconp1,uconp1,ucovp1,whichvel,VEL4,ggp1,GGp1);
 
     ldouble dl,dr,dc;
     ldouble dl2,dr2,dc2;
     for(i=0;i<4;i++)
       {
	 dc=(ucovp1[i]-ucovm1[i]) / (xxvecp1[idim] - xxvecm1[idim]);
	 dr=(ucovp1[i]-ucov[i]) / (xxvecp1[idim] - xxvec[idim]);
	 dl=(ucov[i]-ucovm1[i]) / (xxvec[idim] - xxvecm1[idim]);

	 dc2=(uconp1[i]-uconm1[i]) / (xxvecp1[idim] - xxvecm1[idim]);
	 dr2=(uconp1[i]-ucon[i]) / (xxvecp1[idim] - xxvec[idim]);
	 dl2=(ucon[i]-uconm1[i]) / (xxvec[idim] - xxvecm1[idim]);
	 
	 //TODO : whether MPI handles this properly
	 //to avoid corners
	 if((ix<0 && iy==0 && iz==0 && idim!=1) ||
	    (iy<0 && ix==0 && iz==0 && idim!=2) ||
	    (iz<0 && ix==0 && iy==0 && idim!=3))
	   {
	     du[i][idim]=dr;
	     du2[i][idim]=dr2;
	   }
	 else if((ix<0 && iy==NY-1 && iz==NZ-1 && idim!=1) ||
	    (iy<0 && ix==NX-1 && iz==NZ-1 && idim!=2) ||
	    (iz<0 && ix==NX-1 && iy==NY-1 && idim!=3))
	   {
	     du[i][idim]=dl;
	     du2[i][idim]=dl2;
	   }
	 else if((ix>=NX && iy==0 && iz==0 && idim!=1) ||
	    (iy>=NY && ix==0 && iz==0 && idim!=2) ||
	    (iz>=NZ && ix==0 && iy==0 && idim!=3))
	   {
	     du[i][idim]=dr;
	     du2[i][idim]=dr2;
	   }
	 else if((ix>=NX && iy==NY-1 && iz==NZ-1 && idim!=1) ||
	    (iy>=NY && ix==NX-1 && iz==NZ-1 && idim!=2) ||
	    (iz>=NZ && ix==NX-1 && iy==NY-1 && idim!=3))
	   {
	     du[i][idim]=dl;
	     du2[i][idim]=dl2;
	   }
	 else
	   {
	     //choice of 1st order derivative
	     if(derdir[idim-1]==0)
	       {
		 du[i][idim]=dc;
		 du2[i][idim]=dc2;
	       }
	     if(derdir[idim-1]==1)
	       {
		 du[i][idim]=dl;
		 du2[i][idim]=dl2;
	       }
	     if(derdir[idim-1]==2)
	       {
		 du[i][idim]=dr;
		 du2[i][idim]=dr2;
	       }
	   }

	 if(isnan(du[i][idim]) && !doingavg) {
	   printf("nan in shear_lab : %d %d %d %d\n",ix,iy,iz,idim);
	   print_4vector(ucovm1);
	   print_4vector(ucov);
	   print_4vector(ucovp1);
	   print_4vector(xxvecm1);
	   print_4vector(xxvec);	   
	   print_4vector(xxvecp1);
	   //getchar();
	 }	 	 
       }       
    }

  //covariant derivative tensor du_i;j
  ldouble dcu[4][4];
  
  //covariant derivative tensor du^i;j - only for expansion
  ldouble dcu2[4][4];
  ldouble Krsum;

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Krsum=0.;
	  for(k=0;k<4;k++)
	    Krsum+=get_gKr(k,i,j,ix,iy,iz)*ucov[k];

	  dcu[i][j] = du[i][j] - Krsum;
	}
    
      //only diagonal terms for expansion
      Krsum=0.;
      for(k=0;k<4;k++)
	Krsum+=get_gKr(i,i,k,ix,iy,iz)*ucon[k];
    
      dcu2[i][i] = du2[i][i] + Krsum; 
    }

  //expansion u^{\mu}_{;\mu}
  ldouble theta=0.;
  for(i=0;i<4;i++)
    theta += dcu2[i][i];
  *div = theta;

  //projection tensors P11=P_ij, P21=P^i_j
  ldouble P11[4][4],P21[4][4];
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	P11[i][j] = gg[i][j] + ucov[i]*ucov[j];
	P21[i][j] = delta(i,j) + ucon[i]*ucov[j];
      }

  //the shear tensor sigma_ij - only spatial components
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      {
	ldouble sum1,sum2;
	sum1=sum2=0.;
	for(k=0;k<4;k++)
	  {
	    sum1+=dcu[i][k]*P21[k][j];
	    sum2+=dcu[j][k]*P21[k][i];
	  }
	S[i][j] = 0.5*(sum1+sum2) - 1./3.*theta*P11[i][j];

     }

  //filling the time component from u^mu sigma_munu = 0 
  //(zero time derivatives in the comoving frame - no need for separate comoving routine)
  for(i=1;i<4;i++)
    S[i][0]=S[0][i]=-1./ucon[0]*(ucon[1]*S[1][i]+ucon[2]*S[2][i]+ucon[3]*S[3][i]);

  S[0][0]=-1./ucon[0]*(ucon[1]*S[1][0]+ucon[2]*S[2][0]+ucon[3]*S[3][0]);

  
  return 0;
}

//calculate the divergence of  4-velocity using u^\mu_;\mu = (sqrt(g)*u^\mu)_,\mu / sqrt(g)
int
calc_fluid_div_lab(ldouble *pp0, void* ggg, ldouble dt, ldouble *div, int hdorrad, int *derdir)
{
  int i,j,k,iv;

  struct geometry *geomin
    = (struct geometry *) ggg;

  int ix,iy,iz;
  ix=geomin->ix;
  iy=geomin->iy;
  iz=geomin->iz;

  //get geometry again, just in case geomin is not in internal coordinates, like in postprocessing
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  ldouble (*gg)[5],(*GG)[5];
  ldouble gdet;
  gg=geom.gg;
  GG=geom.GG;
  gdet=geom.gdet;
  
  //start with derivatives
  ldouble du2[4]; //(g*du^i),i

  int istart,whichvel;
  if(hdorrad==MHD)
    {
      whichvel=VELPRIM;
      istart=VX;
    }
  else if(hdorrad==RAD)
    {
      whichvel=VELPRIMRAD;
      istart=FX;
    }
  
  ldouble ppm1[NV],ppp1[NV],pp[NV];
  ldouble gdetm1,gdetp1;
  ldouble ggm1[4][5],GGm1[4][5];
  ldouble ggp1[4][5],GGp1[4][5];

  ldouble xxvecm1[4],xxvec[4],xxvecp1[4];
  ldouble uconm1[4],uconp1[4],utconm1[4],utconp1[4],utcon[4],ucon[4];
  ldouble ucovm1[4],ucovp1[4],ucov[4];
  ldouble enl,enr;
  int idim;

  //four-velocity at cell basing on pp[]
  get_xx(ix,iy,iz,xxvec);
  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=pp0[iv];
    }
  
  utcon[1]=pp[istart];  utcon[2]=pp[istart+1];  utcon[3]=pp[istart+2];
  conv_vels_both(utcon,ucon,ucov,whichvel,VEL4,gg,GG);  

  //time derivative
  du2[0] = 0.;
  
#ifdef TIMEINEXPANSION
  uconm1[0]=0.; //time component will be calculated
  uconm1[1]=get_u(ptm1,istart,  ix,iy,iz);
  uconm1[2]=get_u(ptm1,istart+1,ix,iy,iz);
  uconm1[3]=get_u(ptm1,istart+2,ix,iy,iz); 
  conv_vels(uconm1,uconm1,whichvel,VEL4,geom.gg,geom.GG);
  du2[0] = gdet * (ucon[0] - uconm1[0])/dt; //ANDREW assume metric time-independent. 
#endif
  
  //spatial derivatives
  for(idim=1;idim<4;idim++)
  {
      if(idim==1)
	{
	  get_xx(ix-1,iy,iz,xxvecm1);
	  get_xx(ix+1,iy,iz,xxvecp1);

	  if(hdorrad==RAD) 
	    {
	      enl=get_u(u,EE0,ix-1,iy,iz);
	      enr=get_u(u,EE0,ix+1,iy,iz);
	    }
	  else
	    {
	      enl=get_u(u,UU,ix-1,iy,iz);
	      enr=get_u(u,UU,ix+1,iy,iz);
	    }
       
	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix-1,iy,iz);
	      ppp1[iv]=get_u(p,iv,ix+1,iy,iz);
	    }

	  pick_g(ix-1,iy,iz,ggm1);  pick_G(ix-1,iy,iz,GGm1);
	  pick_g(ix+1,iy,iz,ggp1);  pick_G(ix+1,iy,iz,GGp1);

	  gdetm1=pick_gdet(ix-1,iy,iz);  
	  gdetp1=pick_gdet(ix+1,iy,iz); 
	}

      if(idim==2)
	{
	  get_xx(ix,iy-1,iz,xxvecm1);
	  get_xx(ix,iy+1,iz,xxvecp1);	  
	  if(hdorrad==RAD) 
	    {
	      enl=get_u(u,EE0,ix,iy-1,iz);
	      enr=get_u(u,EE0,ix,iy+1,iz);
	    }
	  else
	    {
	      enl=get_u(u,UU,ix,iy-1,iz);
	      enr=get_u(u,UU,ix,iy+1,iz);
	    }

	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix,iy-1,iz);
	      ppp1[iv]=get_u(p,iv,ix,iy+1,iz);
	    }

	  pick_g(ix,iy-1,iz,ggm1);  pick_G(ix,iy-1,iz,GGm1);
	  pick_g(ix,iy+1,iz,ggp1);  pick_G(ix,iy+1,iz,GGp1);
	  
	  gdetm1=pick_gdet(ix,iy-1,iz); 
	  gdetp1=pick_gdet(ix,iy+1,iz); 
	    
	}

     if(idim==3)
       {
	 get_xx(ix,iy,iz-1,xxvecm1);
	 get_xx(ix,iy,iz+1,xxvecp1);
	 enl=get_u(u,UU,ix,iy,iz-1);
	 enr=get_u(u,UU,ix,iy,iz+1);
	 
	 for(iv=0;iv<NV;iv++)
	   {
	     ppm1[iv]=get_u(p,iv,ix,iy,iz-1);
	     ppp1[iv]=get_u(p,iv,ix,iy,iz+1);
	   }
	 if(hdorrad==RAD) 
	    {
	      enl=get_u(u,EE0,ix,iy,iz-1);
	      enr=get_u(u,EE0,ix,iy,iz+1);
	    }
	 else
	    {
	      enl=get_u(u,UU,ix,iy,iz-1);
	      enr=get_u(u,UU,ix,iy,iz+1);
	    }

	 pick_g(ix,iy,iz-1,ggm1);  pick_G(ix,iy,iz-1,GGm1);
	 pick_g(ix,iy,iz+1,ggp1);  pick_G(ix,iy,iz+1,GGp1);
	 
	 gdetm1=pick_gdet(ix,iy,iz-1); 
	 gdetp1=pick_gdet(ix,iy,iz+1); 
       }

     //calculating four velocity

     utconm1[1]=ppm1[istart];  utconm1[2]=ppm1[istart+1];  utconm1[3]=ppm1[istart+2];
     utconp1[1]=ppp1[istart];  utconp1[2]=ppp1[istart+1];  utconp1[3]=ppp1[istart+2];

     conv_vels_both(utconm1,uconm1,ucovm1,whichvel,VEL4,ggm1,GGm1);
     conv_vels_both(utconp1,uconp1,ucovp1,whichvel,VEL4,ggp1,GGp1);
 
     ldouble dl2,dr2,dc2;

	 dc2=(gdetp1*uconp1[idim]-gdetm1*uconm1[idim]) / (xxvecp1[idim] - xxvecm1[idim]);
	 dr2=(gdetp1*uconp1[idim]-gdet*ucon[idim]) / (xxvecp1[idim] - xxvec[idim]);
	 dl2=(gdet*ucon[idim]-gdetm1*uconm1[idim]) / (xxvec[idim] - xxvecm1[idim]);
	 
	 //TODO : whether MPI handles this properly
	 //to avoid corners
	 if((ix<0 && iy==0 && iz==0 && idim!=1) ||
	    (iy<0 && ix==0 && iz==0 && idim!=2) ||
	    (iz<0 && ix==0 && iy==0 && idim!=3))
	   {
	     du2[idim]=dr2;
	   }
	 else if((ix<0 && iy==NY-1 && iz==NZ-1 && idim!=1) ||
	    (iy<0 && ix==NX-1 && iz==NZ-1 && idim!=2) ||
	    (iz<0 && ix==NX-1 && iy==NY-1 && idim!=3))
	   {
	     du2[idim]=dl2;
	   }
	 else if((ix>=NX && iy==0 && iz==0 && idim!=1) ||
	    (iy>=NY && ix==0 && iz==0 && idim!=2) ||
	    (iz>=NZ && ix==0 && iy==0 && idim!=3))
	   {
	     du2[idim]=dr2;
	   }
	 else if((ix>=NX && iy==NY-1 && iz==NZ-1 && idim!=1) ||
	    (iy>=NY && ix==NX-1 && iz==NZ-1 && idim!=2) ||
	    (iz>=NZ && ix==NX-1 && iy==NY-1 && idim!=3))
	   {
	     du2[idim]=dl2;
	   }
	 else
	   {
	     //choice of 1st order derivative
	     #ifdef NUMRADWAVESPEEDS
	     if(fabs(enl)>fabs(enr))
	       {
		 du2[idim]=dl2;
	       }
	     else
	       {
		 du2[idim]=dr2;
	       }
	     #else

	     if(derdir[idim-1]==0)
	       {
		 du2[idim]=dc2;
	       }
	     if(derdir[idim-1]==1)
	       {
		 du2[idim]=dl2;
	       }
	     if(derdir[idim-1]==2)
	       {
		 du2[idim]=dr2;
	       }
	     #endif
	   }

	 if(isnan(du2[idim]) && !doingavg) {
	   printf("nan in calc_div : %d %d %d %d\n",ix,iy,iz,idim);
	   print_4vector(ucovm1);
	   print_4vector(ucov);
	   print_4vector(ucovp1);
	   print_4vector(xxvecm1);
	   print_4vector(xxvec);	   
	   print_4vector(xxvecp1);
	   //getchar();
	  }	 	 
       }       

   //Calculate div
   
   ldouble theta=0.;
   for (i=0; i<4 ; i++) theta+= du2[i];

   theta = theta/gdet;
   *div = theta;

   return 0;
}

//**********************************************************************
//calculates mean free path used in the viscosity coefficient
//**********************************************************************
int
calc_rad_visccoeff(ldouble *pp,void *ggg,ldouble *nuret,ldouble *mfpret,ldouble *mindxret)
{
#if (RADVISCOSITY==SHEARVISCOSITY)
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble chi=calc_chi(pp,geom);
  ldouble mfp = 1./chi; // dr / dtau
  #ifdef SKIPRADSOURCE
  mfp=1.e50;
  #endif
  ldouble mindx,nu;

  //limiting in opt.thin region

  //minimal cell size
  ldouble dx[3]={get_size_x(geom->ix,0)*sqrt(geom->gg[1][1]),   
		 get_size_x(geom->iy,1)*sqrt(geom->gg[2][2]),
		 get_size_x(geom->iz,2)*sqrt(geom->gg[3][3])};

  if(NY==1 && NZ==1) mindx = dx[0];
  else if(NZ==1) mindx = my_min(dx[0],dx[1]);
  else if(NY==1) mindx = my_min(dx[0],dx[2]);
  else mindx = my_min(dx[0],my_min(dx[1],dx[2]));

//choice of mean free path limiter for tau<<1
#ifdef RADVISCMFPCONST

  if(mfp>RADVISCMFPCONST || chi<SMALL) mfp=RADVISCMFPCONST;

#elif defined(RADVISCMFPCYL)

  ldouble xxBL[4]; 
  coco_N(geom->xxvec,xxBL,MYCOORDS, BLCOORDS);
  ldouble mfplim=xxBL[1]*sin(xxBL[2]); //Rcyl = Rsph * sin(th)
  if(mfp>mfplim || chi<SMALL) mfp=mfplim; //Rcyl = Rsph
  if(mfp<0. || !isfinite(mfp)) mfp=0.;
 
#elif defined(RADVISCMFPSPH)

  ldouble xxBL[4];
  #ifdef RADVISCMFPSPHRMIN
  ldouble rmin=RADVISCMFPSPHRMIN;
  #else
  ldouble rmin=1.2*rhorizonBL;
  #endif
  
  coco_N(geom->xxvec,xxBL,MYCOORDS, BLCOORDS);
  ldouble mfplim=xxBL[1];
  if(mfp>mfplim || chi<SMALL) mfp=mfplim; //Rcyl = Rsph
  if(mfp<0. || !isfinite(mfp)) mfp=0.;
  mfp*=step_function(xxBL[1]-rmin,0.2*rmin);
  if(xxBL[1]<=1.*rhorizonBL) mfp=0.;

#ifdef RADVISCMFPSPHMAX
  if(mfp>RADVISCMFPSPHMAX) 
    mfp=RADVISCMFPSPHMAX;
#endif

#elif defined(RADVISCMFPNOLIMIT)

  ldouble xxBL[4];
  ldouble rhor=rhorizonBL;
  coco_N(geom->xxvec,xxBL,MYCOORDS, BLCOORDS);
  ldouble mfplim=100.*xxBL[1];
  if(mfp>mfplim || chi<SMALL) mfp=mfplim; //Rcyl = Rsph
  if(mfp<0. || !isfinite(mfp)) mfp=0.;
  mfp*=step_function(xxBL[1]-1.1*rhor,0.01*rhor);
  if(xxBL[1]<=1.*rhor) mfp=0.;

#else

  if(mfp>mindx || chi<SMALL) mfp=mindx;

#endif

  nu = ALPHARADVISC * mfp;

  /**********************************/
  //damping the viscous term based on maximal diffusion coefficient for given dt ans cell size
  /**********************************/

#ifdef RADVISCNUDAMP
  ldouble nulimit = mindx*mindx / 2. / global_dt / 2.;
  ldouble fac=nu/nulimit;
  if(nu>nulimit)
    {
      //printf("famping at %d %d | %e | %e %e | %e\n",geom->ix,geom->iy,nu,mfp,mindx,dt); getchar();
      nu=nulimit;
    }
#endif


  *nuret=nu;
  *mfpret=mfp;
  *mindxret=mindx;
#endif
  return 0;
}


//**********************************************************************
//*** calculates Rij visc for all cells ********************************
//**********************************************************************
int
calc_Rij_visc_total()
{
#ifdef RADIATION
#if (RADVISCOSITY==SHEARVISCOSITY)

  int ii;
//#pragma omp parallel for private(ii)
  for(ii=0;ii<Nloop_1;ii++) //domain plus some ghost cells
  {
    int ix, iy, iz, i, j;
    struct geometry geomcent;
    int derdir[3]={0,0,0};
    ldouble Rvisc[4][4];

    ix=loop_1[ii][0];
    iy=loop_1[ii][1];
    iz=loop_1[ii][2];
    
    if(if_outsidegc(ix,iy,iz)==1) continue; //avoid corners
    
    fill_geometry(ix,iy,iz,&geomcent);
    calc_Rij_visc(&get_u(p,0,ix,iy,iz),&geomcent,Rvisc,derdir);
    indices_2221(Rvisc,Rvisc,geomcent.gg); //R^i_j
    
    //save to memory
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
      {
        set_Tfull(Rijviscglobal,i,j,ix,iy,iz,Rvisc[i][j]);
      }
  }
  
#endif //SHEARVISCOSITY
#endif //RADVISCOSITY
   return 0;
}

//**********************************************************************
//*** calculates Rij visc for single cell ******************************
//**********************************************************************
int
calc_Rij_visc(ldouble *pp, void* ggg, ldouble Rvisc[][4], int *derdir)
{
  int i,j,ix,iy,iz;
  
  struct geometry *geom
   = (struct geometry *) ggg;

  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;
  
  for(i=0;i<4;i++)
  {
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
    for(j=0;j<4;j++)
    {
      Rvisc[i][j]=0.;
    }
  }

#ifdef RADIATION  
#if (RADVISCOSITY==SHEARVISCOSITY)
  //verify if recalculating shear velocity needed
  int recalcvisc=1;

  #ifdef ACCELRADVISCOSITY
  ldouble lasttime=get_u_scalar(radvisclasttime,ix,iy,iz);
  if(lasttime>0.)
    {
      ldouble dtlast = global_time - lasttime;
      ldouble dtmin = global_dt;
      ldouble sizevec[3] = {get_size_x(ix,0)*sqrt(get_g(g,1,1,ix,iy,iz)),
			    get_size_x(iy,1)*sqrt(get_g(g,2,2,ix,iy,iz)),
			    get_size_x(iz,2)*sqrt(get_g(g,3,3,ix,iy,iz))};
      ldouble minsize;
      if(TNZ==1 && TNY==1)
	minsize=sizevec[0];
      else if(TNZ==1)
	minsize=my_min(sizevec[0],sizevec[1]);
      else
	minsize=my_min_N(sizevec,3);
            
      //if(dtlast/dtmin < radsizefrac)
      if(dtlast < 0.5*minsize)
	{
	  recalcvisc=0.;
	  /*
	  if(iy==0 || 1)
	    {
	      printf("skipping recalculating Rijvisc at %d %d with radsizefrac = %e | %e %e %e %e\n",
		     ix,iy,minsize,global_time,radvisclasttime[ix+NGCX][iy+NGCY][iz+NGCZ],dtmin,minsize);
	    }
	  */
	}
    }
  #endif //ACCELRADVISCOSITY

  if(recalcvisc==1)
    {
      ldouble shear[4][4],nu;
      calc_rad_shearviscosity(pp,ggg,shear,&nu,derdir);
      
      //radiation rest frame energy density :
      ldouble Erad;
      Erad=pp[EE0]; 
      
      //multiply by viscosity to get viscous tensor
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  {
	    Rvisc[i][j]= -2. * nu * Erad * shear[i][j];

	    set_Tfull(Rijviscprev,i,j,ix,iy,iz,Rvisc[i][j]);
	  }

      set_u_scalar(radvisclasttime,ix,iy,iz,global_time);
    }
  else
    {
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  {
	    Rvisc[i][j]=get_Tfull(Rijviscprev,i,j,ix,iy,iz);
	  }
    }

#endif //SHEARVISCOSITY
#endif //RADIATION
  return 0;

}

///////////////////////////////////////////////////////////////
//resetting accelerating arrays
void
reset_radviscaccel()
{
#ifdef RADIATION
#if (RADVISCOSITY==SHEARVISCOSITY)

  int ix,iy,iz,ii;

  for(ii=0;ii<Nloop_02;ii++) //domain and ghost cells 
    {
      ix=loop_02[ii][0];
      iy=loop_02[ii][1];
      iz=loop_02[ii][2]; 
      //radvisclasttime[ix][iy][iz]=-1.;
      set_u_scalar(radvisclasttime,ix,iy,iz,-1);
    }

  #endif
#endif
}


///////////////////////////////////////////////////////////////
//estimates dBi/dt from the battery term
int
estimate_Bgrowth_battery(int ix,int iy,int iz,ldouble dBdt[4])
{
  int coords=BLCOORDS; 
  //this routine used only when exporting to silo

  ldouble fxr[4],fxl[4],fyr[4],fyl[4],fzr[4],fzl[4];
  int iv;
  for(iv=0;iv<4;iv++)
    fxl[iv]=fxr[iv]=fyl[iv]=fyr[iv]=fzl[iv]=fzr[iv]=0.;

  struct geometry geom,geomxr,geomxl,geomyr,geomyl,geomzr,geomzl;
  struct geometry geomBL,geomBLxr,geomBLxl,geomBLyr,geomBLyl,geomBLzr,geomBLzl;
 
  fill_geometry(ix,iy,iz,&geom);
  fill_geometry(ix-1,iy,iz,&geomxl);
  fill_geometry(ix+1,iy,iz,&geomxr);
  fill_geometry(ix,iy-1,iz,&geomyl);
  fill_geometry(ix,iy+1,iz,&geomyr);
  fill_geometry(ix,iy,iz-1,&geomzl);
  fill_geometry(ix,iy,iz+1,&geomzr);
   
  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
  fill_geometry_arb(ix-1,iy,iz,&geomBLxl,BLCOORDS);
  fill_geometry_arb(ix+1,iy,iz,&geomBLxr,BLCOORDS);
  fill_geometry_arb(ix,iy-1,iz,&geomBLyl,BLCOORDS);
  fill_geometry_arb(ix,iy+1,iz,&geomBLyr,BLCOORDS);
  fill_geometry_arb(ix,iy,iz-1,&geomBLzl,BLCOORDS);
  fill_geometry_arb(ix,iy,iz+1,&geomBLzr,BLCOORDS);
 
  struct struct_of_state state;
  ldouble pp[NV];
  ldouble dx,dy,dz,xxvec1[4],xxvec2[4];

  if(coords!=BLCOORDS)
    {
      fill_struct_of_state(&get_u(p,0,ix-1,iy,iz),&geomxl,&state);
      calc_batteryflux(&get_u(p,0,ix-1,iy,iz),&geomxl,fxl,0,state.ucov);

      fill_struct_of_state(&get_u(p,0,ix+1,iy,iz),&geomxr,&state);
      calc_batteryflux(&get_u(p,0,ix+1,iy,iz),&geomxr,fxr,0,state.ucov);

      fill_struct_of_state(&get_u(p,0,ix,iy-1,iz),&geomyl,&state);
      calc_batteryflux(&get_u(p,0,ix,iy-1,iz),&geomyl,fyl,1,state.ucov);

      fill_struct_of_state(&get_u(p,0,ix,iy+1,iz),&geomyr,&state);
      calc_batteryflux(&get_u(p,0,ix,iy+1,iz),&geomyr,fyr,1,state.ucov);

      fill_struct_of_state(&get_u(p,0,ix,iy,iz-1),&geomzl,&state);
      calc_batteryflux(&get_u(p,0,ix,iy,iz-1),&geomzl,fzl,2,state.ucov);

      fill_struct_of_state(&get_u(p,0,ix,iy,iz+1),&geomzr,&state);
      calc_batteryflux(&get_u(p,0,ix,iy,iz+1),&geomzr,fzr,2,state.ucov);

      dx = get_x(ix+1,0) - get_x(ix-1,0);
      dy = get_x(iy+1,1) - get_x(iy-1,1);
      dz = get_x(iz+1,2) - get_x(iz-1,2);
  
      dBdt[1] = ((fxr[1]*geomxr.gdet - fxl[1]*geomxl.gdet)/dx + (fyr[1]*geomyr.gdet - fyl[1]*geomyl.gdet)/dy + (fzr[1]*geomzr.gdet - fzl[1]*geomzl.gdet)/dz)/geom.gdet;
      dBdt[2] = ((fxr[2]*geomxr.gdet - fxl[2]*geomxl.gdet)/dx + (fyr[2]*geomyr.gdet - fyl[2]*geomyl.gdet)/dy + (fzr[2]*geomzr.gdet - fzl[2]*geomzl.gdet)/dz)/geom.gdet;
      dBdt[3] = ((fxr[3]*geomxr.gdet - fxl[3]*geomxl.gdet)/dx + (fyr[3]*geomyr.gdet - fyl[3]*geomyl.gdet)/dy + (fzr[3]*geomzr.gdet - fzl[3]*geomzl.gdet)/dz)/geom.gdet;
  
    }
  else if(coords==BLCOORDS)
    {

      PLOOP(iv) pp[iv]=get_u(p,iv,ix-1,iy,iz);
      trans_pall_coco(pp,pp,MYCOORDS,BLCOORDS,geomxl.xxvec,&geomxl,&geomBLxl);
      fill_struct_of_state(pp,&geomBLxl,&state);
      calc_batteryflux(pp,&geomBLxl,fxl,0,state.ucov);
  
      PLOOP(iv) pp[iv]=get_u(p,iv,ix+1,iy,iz);
      trans_pall_coco(pp,pp,MYCOORDS,BLCOORDS,geomxr.xxvec,&geomxr,&geomBLxr);
      fill_struct_of_state(pp,&geomBLxr,&state);
      calc_batteryflux(pp,&geomBLxr,fxr,0,state.ucov);
  
      PLOOP(iv) pp[iv]=get_u(p,iv,ix,iy-1,iz);
      trans_pall_coco(pp,pp,MYCOORDS,BLCOORDS,geomyl.xxvec,&geomyl,&geomBLyl);
      fill_struct_of_state(pp,&geomBLyl,&state);
      calc_batteryflux(pp,&geomBLyl,fyl,1,state.ucov);
  
      PLOOP(iv) pp[iv]=get_u(p,iv,ix,iy+1,iz);
      trans_pall_coco(pp,pp,MYCOORDS,BLCOORDS,geomyr.xxvec,&geomyr,&geomBLyr);
      fill_struct_of_state(pp,&geomBLyr,&state);
      calc_batteryflux(pp,&geomBLyr,fyr,1,state.ucov);
  
      PLOOP(iv) pp[iv]=get_u(p,iv,ix,iy,iz-1);
      trans_pall_coco(pp,pp,MYCOORDS,BLCOORDS,geomzl.xxvec,&geomzl,&geomBLzl);
      fill_struct_of_state(pp,&geomBLzl,&state);
      calc_batteryflux(pp,&geomBLzl,fzl,2,state.ucov);
  
      PLOOP(iv) pp[iv]=get_u(p,iv,ix,iy,iz+1);
      trans_pall_coco(pp,pp,MYCOORDS,BLCOORDS,geomzr.xxvec,&geomzr,&geomBLzr);
      fill_struct_of_state(pp,&geomBLzr,&state);
      calc_batteryflux(pp,&geomBLzr,fzr,2,state.ucov);
    
  
      get_xx_arb(ix+1,iy,iz,xxvec2,BLCOORDS);
      get_xx_arb(ix-1,iy,iz,xxvec1,BLCOORDS);
      dx=xxvec2[1]-xxvec1[1];
      get_xx_arb(ix,iy+1,iz,xxvec2,BLCOORDS);
      get_xx_arb(ix,iy-1,iz,xxvec1,BLCOORDS);
      dy=xxvec2[2]-xxvec1[2];
      get_xx_arb(ix,iy,iz+1,xxvec2,BLCOORDS);
      get_xx_arb(ix,iy,iz-1,xxvec1,BLCOORDS);
      dz=xxvec2[3]-xxvec1[3];
    
      dBdt[1] = (-(fxr[1]*geomBLxr.gdet - fxl[1]*geomBLxl.gdet)/dx - (fyr[1]*geomBLyr.gdet - fyl[1]*geomBLyl.gdet)/dy - (fzr[1]*geomBLzr.gdet - fzl[1]*geomBLzl.gdet)/dz)/geomBL.gdet;
      dBdt[2] = (-(fxr[2]*geomBLxr.gdet - fxl[2]*geomBLxl.gdet)/dx - (fyr[2]*geomBLyr.gdet - fyl[2]*geomBLyl.gdet)/dy - (fzr[2]*geomBLzr.gdet - fzl[2]*geomBLzl.gdet)/dz)/geomBL.gdet;
      dBdt[3] = (-(fxr[3]*geomBLxr.gdet - fxl[3]*geomBLxl.gdet)/dx - (fyr[3]*geomBLyr.gdet - fyl[3]*geomBLyl.gdet)/dy - (fzr[3]*geomBLzr.gdet - fzl[3]*geomBLzl.gdet)/dz)/geomBL.gdet;
    
  
    }

  return 0;
}

int
calc_batteryflux(ldouble *pp, void* ggg, ldouble *eterm,int idim,ldouble *ucov)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble econ[4],ecov[4];
  calc_Efield_battery(pp,geom,econ);
  indices_21(econ, ecov, geom->gg);


  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	eterm[1]+=-epsilon_4d_tensor(1,idim+1,i,j,ggg)*ecov[i]*ucov[j];
	eterm[2]+=-epsilon_4d_tensor(2,idim+1,i,j,ggg)*ecov[i]*ucov[j];
	eterm[3]+=-epsilon_4d_tensor(3,idim+1,i,j,ggg)*ecov[i]*ucov[j];

     }
  return 0;
}

//calculates electric field in the fluid frame due to the Poynting-Robertson drag
int
calc_Efield_battery(ldouble *pp,void *ggg,ldouble econ[4])
{
  #ifdef BATTERY
  int i,j,k;

  econ[0]=econ[1]=econ[2]=econ[3]=0.;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  //four-velocity
  ldouble ucon[4],ucov[4];
  ldouble bcon[4],bcov[4];
  calc_ucon_ucov_from_prims(pp, geom, ucon, ucov);
  calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, geom, bcon, bcov, &bsq);
  ldouble b2rho=bsq/pp[RHO];
  ldouble suppfac = 1.;

  //projection tensors h^mu_nu
  ldouble hmunu[4][4];
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {	
	hmunu[i][j] = delta(i,j) + ucon[i]*ucov[j];
      }

  //radiative four-force
  ldouble Gi[4],Giff[4];
  ldouble Gith[4], Githff[4];
  calc_Gi(pp, ggg, Gi, 0.0, 1, 0); //ANDREW no nonthermal in battery problem

  //electric field
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	//E^i = sigma_Th F^i / e / c = kappa_es mp F^i / e / c = kappa_es rho F^i mp / rho / e / c = G^i M / rho / e / c;
	econ[i]+=suppfac*BATTERYMULTIPLIER*hmunu[i][j]*Gi[j]/pp[RHO]*M_PROTON/E_CHARGE;
      }

  //limit econ < bcon
#ifdef BATTERYMAXESQOVERRHO
  ldouble ecov[4];
  indices_21(econ,ecov,gg); 
  ldouble esq = dotB(econ,ecov); 
  ldouble e2rho=esq/pp[RHO];
  if(e2rho>BATTERYMAXESQOVERRHO)
    {
      for(i=0;i<4;i++)
	econ[i]*=sqrt(BATTERYMAXESQOVERRHO/e2rho);
    }
#endif
	       
#ifdef LIMITBATTERYTOPHI //take only azimuthal fluid frame flux into consideration
  econ[0]=econ[1]=econ[2]=0.;
#endif

#endif
  return 0;  
}


///////////////////////////////////////////////////////////////
//estimates how much gas is coupled to radiation (or how large is the optical depth)
//by looking at how close is the gas comoving radiation stress energy tensor to the
//Eddington tensor - which is the limit for optically thick gas
//this is done by calculating the difference between gas and radiation rest frame velocity
//if the two velocities agree - this is the optically thick limit, if they are far apart - this
//corresponds to optically thin limit
//returns a number in range 0...1 with 0 - optically thin, 1 - completely optically thick

ldouble
estimate_gas_radiation_coupling(ldouble *pp, void *ggg)
{
  //coordinate system to compute the four-velocities in
  int coords=BLCOORDS;
  if(MYCOORDS==MINKCOORDS)
    coords=MINKCOORDS;

  struct geometry *geom
    = (struct geometry *) ggg;
  struct geometry geomBL;
  fill_geometry_arb(geom->ix,geom->iy,geom->iz,&geomBL,coords);

  //**********************************************************************
  //***** compute gas four velocity ************************************
  //**********************************************************************

  ldouble utcon[4]={0,pp[VX],pp[VY],pp[VZ]};
  ldouble ucon[4], uconbl[4];
  conv_vels(utcon,ucon,VELPRIM,VEL4,geom->gg,geom->GG);
  trans2_coco(geom->xxvec,ucon,uconbl,geom->coords,coords);

  //**********************************************************************
  //***** compute radiation rest frame four velocity *******************
  //**********************************************************************

  ldouble urftcon[4]={0,pp[FX],pp[FY],pp[FZ]};
  ldouble urfcon[4], urfconbl[4];
  conv_vels(urftcon,urfcon,VELPRIMRAD,VEL4,geom->gg,geom->GG);
  trans2_coco(geom->xxvec,urfcon,urfconbl,geom->coords,coords);

  //**********************************************************************
  //***** making the velocities proper and ortonormal **************
  //**********************************************************************
  
  uconbl[1]*=sqrt(geomBL.gg[1][1])/uconbl[0];
  uconbl[2]*=sqrt(geomBL.gg[2][2])/uconbl[0];
  uconbl[3]*=sqrt(geomBL.gg[3][3])/uconbl[0];
  urfconbl[1]*=sqrt(geomBL.gg[1][1])/urfconbl[0];
  urfconbl[2]*=sqrt(geomBL.gg[2][2])/urfconbl[0];
  urfconbl[3]*=sqrt(geomBL.gg[3][3])/urfconbl[0];

  //**********************************************************************
  //***** computing the difference of these proper velocities **********
  //**********************************************************************
  
  ldouble dv[3],dvmag;
  dv[0]=uconbl[1]-urfconbl[1];
  dv[1]=uconbl[2]-urfconbl[2];
  dv[2]=uconbl[3]-urfconbl[3];
  dvmag=sqrt(dot3(dv,dv));

  //**********************************************************************
  //***** computing the coupling measure *****************************
  //**********************************************************************
  
  ldouble factor;
  factor = 1.0 - step_function(dvmag - COUPLING_MEASURE_THRESHOLD, COUPLING_MEASURE_THRESHOLD_SHARPNESS*COUPLING_MEASURE_THRESHOLD);

  //printf("%d %d > %e %e %e > %e %e %e > %e > %f\n",geom->ix,geom->iy,uconbl[1],uconbl[2],uconbl[3],urfconbl[1],urfconbl[2],urfconbl[3],dvmag,factor); 

  return factor;
}


//**********************************************************************
//**********************************************************************
//returns rad primitives for an atmosphere
//**********************************************************************
//**********************************************************************

int
set_radatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype)
{
#ifdef RADIATION   
  if(atmtype==0) //fixed Erf, urf of normal observer
    {
      pp[EE0]=ERADATMMIN; 
      ldouble ucon[4];
      calc_normalobs_relvel(GG,ucon);
      if (VELR != VELPRIMRAD)
      {
        conv_vels(ucon,ucon,VELR,VELPRIMRAD,gg,GG);
      }
 
      pp[FX0]=ucon[1]; 
      pp[FY0]=ucon[2];
      pp[FZ0]=ucon[3];
    }
  if(atmtype==1) //fixed Erf, urf 0 in lab frame
    {
      ldouble ucon[4];
      ldouble xx2[4];
      ldouble GGBL[4][5];

      // BL coords
      coco_N(xx,xx2,MYCOORDS,BLCOORDS);
      calc_G_arb(xx2,GGBL,BLCOORDS);

      // normal observer in BL = stationary observer
      calc_normalobs_4vel(GGBL,ucon);
     
      // to MYCOORDS
      trans2_coco(xx2,ucon,ucon,BLCOORDS,MYCOORDS);
     
      // to VELPRIMRAD
      conv_vels_ut(ucon,ucon,VEL4,VELPRIMRAD,gg,GG);
     
      pp[FX0]=ucon[1];
      pp[FY0]=ucon[2];
      pp[FZ0]=ucon[3];

      // print_4vector(ucon); getchar();
      pp[EE0]=ERADATMMIN; 
     }
  if(atmtype==2) //optically thin atmosphere, scalings from numerical solution of radiall influx
    {
      ldouble gammamax=10.;
      ldouble rout=2.; //normalize at r_BL=2

      ldouble xxBL[4];
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      ldouble r=xxBL[1];
     
      pp[EE0]=ERADATMMIN*(rout/r)*(rout/r)*(rout/r)*(rout/r);

      ldouble ut[4]={0.,-gammamax*pow(r/rout,1.),0.,0.};

      ldouble ggBL[4][5],GGBL[4][5];
      calc_g_arb(xxBL,ggBL,KERRCOORDS);
      calc_G_arb(xxBL,GGBL,KERRCOORDS);

      conv_vels(ut,ut,VELR,VEL4,ggBL,GGBL);

      trans2_coco(xxBL,ut,ut,KERRCOORDS,MYCOORDS);

      conv_vels_ut(ut,ut,VEL4,VELPRIM,gg,GG);
      
      pp[FX0]=ut[1];      
      pp[FY0]=ut[2];      
      pp[FZ0]=ut[3];

    }
  if(atmtype==3) //LTE, normal observer
    {
      ldouble Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],0,0,0); //indices inconsistent with CONSISTENTGAMMA
      pp[EE0]=calc_LTE_EfromT(Tgas);

      ldouble ucon[4];
      calc_normalobs_relvel(GG,ucon);
      if (VELR != VELPRIMRAD)
      {
        conv_vels(ucon,ucon,VELR,VELPRIMRAD,gg,GG);
      }
 
      pp[FX0]=ucon[1]; 
      pp[FY0]=ucon[2];
      pp[FZ0]=ucon[3];
    }

#ifdef EVOLVEPHOTONNUMBER
  pp[NF0]=calc_NFfromE(pp[EE0]);
#endif
#endif //RADIATION
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//************** tests *************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************

int
test_solve_implicit_lab()
{
#ifdef RADIATION
  int i1,i2,iv;
  ldouble uu0[NV],pp0[NV],pp[NV],uu[NV],dt;
  struct geometry geom;
  fill_geometry(TNX-1,TNY/2,0,&geom);

  pp0[0]=1.e-21;
  pp0[1]=1.e-23;//calc_PEQ_ufromTrho(1.e9,pp0[0]);
  pp0[2]=0.00001;
  pp0[3]=0.0;
  pp0[4]=0.0;
  pp0[5]=calc_Sfromu(pp0[0],pp0[1],TNX-1,TNY/2,0);

  pp0[B1]=pp0[B2]=pp0[B3]=0.;
  pp0[EE]=1.e-23;//calc_LTE_Efromurho(pp0[1],pp0[0]);//1.e-21;
  pp0[FX]=0.00001;
  pp0[FY]=0.0;
  pp0[FZ]=0.0;

#ifdef NCOMPTONNIZATION
  pp0[NF]=1./2.70118/K_BOLTZ * pp0[EE]/calc_PEQ_Tfromurho(pp0[UU],pp0[RHO],geom.ix,geom.iy,geom.iz);
#endif

#ifdef EVOLVEELECTRONS
  ldouble rhogas=pp0[RHO];
  ldouble Tgas=calc_PEQ_Tfromurho(pp0[UU],pp0[RHO],geom.ix,geom.iy,geom.iz);
  ldouble Te=Tgas;
  pp0[ENTRE]=calc_SefromrhoT(rhogas,Te,ELECTRONS); //test_solve_implicit_lab
#endif
  
  p2u(pp0,uu0,&geom);

  dt=1.e0;

  PLOOP(iv) pp[iv]=pp0[iv];

  print_primitives(pp0);
  printf("gas temp: %e\nrad temp: %e\n",calc_PEQ_Tfromurho(pp0[1],pp0[0],geom.ix,geom.iy,geom.iz),calc_LTE_TfromE(pp0[EE]));

  ldouble del4[NRADVAR];
  int verbose=0;
  int params[7];

  int solver=0;
  
  //explicit
  if(0)
    {
      solve_explicit_lab_core(uu0,pp,&geom,dt,del4,verbose);
      ldouble delapl[NV];
      delapl[RHO]=0.;
      delapl[UU]=-del4[0];
      delapl[VX]=-del4[1];
      delapl[VY]=-del4[2];
      delapl[VZ]=-del4[3];
      delapl[ENTR]=-del4[0];
#ifdef MAGNFIELD
      delapl[B1]=0.;
      delapl[B2]=0.;
      delapl[B3]=0.;
#endif
      delapl[EE0]=del4[0];
      delapl[FX0]=del4[1];
      delapl[FY0]=del4[2];
      delapl[FZ0]=del4[3];

      for(iv=0;iv<NV;iv++)
	uu[iv]=uu0[iv]+delapl[iv];

      int corr[3],fixup[2];
      u2p(uu,pp,&geom,corr,fixup,0);

      if(corr[0]!=0 || corr[1]!=0) printf("corr: %d %d\n",corr[0],corr[1]);
      print_primitives(pp);
    }

  verbose=1;
  //full implicit
  PLOOP(iv) pp[iv]=pp0[iv];
  if(1)
    { 
  
      params[0]=MHD;
      params[1]=RADIMPLICIT_ENERGYEQ;
      params[2]=RADIMPLICIT_LAB;
      params[4]=0; //opacity damp
      params[5]=0; //ifncompt
      params[6]=0; //ifelectron
      solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,verbose,params,pp);
      print_primitives(pp);
      //printf("gas temp: %e\nrad temp: %e\n",calc_PEQ_Tfromurho(pp[1],pp[0]),calc_LTE_TfromE(pp[EE]));
    }
#endif
  exit(1);
  
}

int
test_Gi()
{
#ifdef RADIATION
  int i1,i2,iv,i;
  ldouble uu0[NV],pp0[NV],pp[NV],uu[NV],dt;
  struct geometry geom;
  fill_geometry_arb(57,20,0,&geom,BLCOORDS);
  printf("radius: %f\n",geom.xx);

  ldouble rho=rhoCGS2GU(9.253e-9);
  ldouble kappaes=kappaCGS2GU(0.34)*rho;
  
  ldouble Thatrad=7.965e+7;
  ldouble Tgas=2.478e+8;
  ldouble Ehatrad=endenCGS2GU(7.963e+13);
  
  ldouble ucon[4];
  ucon[1]=ucon[2]=ucon[3]=0.;
  ucon[0]=1.;

  ldouble Gic[4];
  for(i=0;i<4;i++)
    Gic[i]=kappaes * Ehatrad * (4.*K_BOLTZ*(Thatrad - Tgas)/M_ELECTR) 
      * (1. + 3.683 * K_BOLTZ * Tgas / M_ELECTR + 4. * K_BOLTZ * Tgas / M_ELECTR * K_BOLTZ * Tgas / M_ELECTR)
      / (1. + 4.*K_BOLTZ*Tgas/M_ELECTR) * ucon[i]; 

  ldouble conv=kappaGU2CGS(1.)*rhoGU2CGS(1.)*endenGU2CGS(1.)*CCC; //because (cE-4piB) in non-geom
	       
  printf("Giff[0] = %.3e\n",Gic[0]*conv);
  exit(1);

  pp0[0]=9.253e-9;
  pp0[1]=calc_PEQ_ufromTrho(2.478e9,pp0[RHO],geom.ix,geom.iy,geom.iz);
  pp0[VX]=-1.04792e-01; 
  pp0[VY]=2.69291e-03;
  pp0[VZ]=1.62996e-02;
  pp0[5]=calc_Sfromu(pp0[0],pp0[1],geom.ix,geom.iy,geom.iz);
  pp0[EE]=1.0*calc_LTE_Efromurho(pp0[UU],pp0[RHO]);
  pp0[FX]=0.01;
  pp0[FY]=0.;
  pp0[FZ]=0.0;
  #ifdef EVOLVEPHOTONNUMBER
  pp0[NF0]=calc_NFfromE(pp0[EE0]);
  #endif

  print_primitives(pp0);

  ldouble Rij[4][4];
  calc_Rij_M1(pp0,&geom,Rij);
  print_4vector(&Rij[0][0]);

  //relativistic radiative four-force (turn off NONRELMHD)
  ldouble Gi[4],Giff[4],Giff2[4];
  calc_Gi(pp0,&geom,Gi, 0.0, 0, 0); //ANDREW no rel electrons in test problem
  print_4vector(Gi);
#endif
  return 0;

}

int
test_solve_implicit_lab_file()
{
  FILE *in = fopen("imp.problem.0","r");

  int i1,i2,iv;
  ldouble uu0[NV],pp0[NV],pp[NV],dt;
  struct geometry geom;

  for (i1=0;i1<NV;i1++)
    iv=fscanf(in,"%lf",&uu0[i1]);
  for (i1=0;i1<NV;i1++)
    iv=fscanf(in,"%lf",&pp0[i1]);
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<5;i2++)
      iv=fscanf(in,"%lf ",&geom.gg[i1][i2]);
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<5;i2++)
      iv=fscanf(in,"%lf ",&geom.GG[i1][i2]);
  iv=fscanf(in,"%lf ",&dt);
  iv=fscanf(in,"%lf ",&geom.alpha);
  iv=fscanf(in,"%lf ",&geom.gdet);
  
  //for imp.problems > 21
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<4;i2++)
      iv=fscanf(in,"%lf ",&geom.tup[i1][i2]);
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<4;i2++)
      iv=fscanf(in,"%lf ",&geom.tlo[i1][i2]);
 fscanf(in,"%lf ",&geom.gttpert);

  for (i1=0;i1<4;i1++)
      iv=fscanf(in,"%lf ",&geom.xxvec[i1]);

  fclose(in);

  print_4vector(geom.xxvec);

  geom.ix=geom.iy=geom.iz=0;  

  ldouble Rttcov[4]={uu0[EE0],uu0[FX0],uu0[FY0],uu0[FZ0]};
  ldouble Rttcon[4];
  indices_12(Rttcov,Rttcon,geom.GG);
  //print_4vector(Rttcov);
  //print_4vector(Rttcon);
  ldouble vcon[4],ucon[4],ucov[4];

  //converting to 4-velocity
  vcon[1]=pp0[2];
  vcon[2]=pp0[3];
  vcon[3]=pp0[4];
  vcon[0]=0.;  
  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,geom.gg,geom.GG);
  print_4vector(ucon);
  print_4vector(ucov);

  ldouble deltas[NRADVAR];
  int verbose=2;
  int params[8];
  
  params[0]=RAD;
  params[1]=RADIMPLICIT_ENERGYEQ;
  params[2]=RADIMPLICIT_LAB;
  params[3]=0; //overshooting check - not used
  params[4]=0; //damping level of opacities
  #ifdef EVOLVEPHOTONNUMBER
  params[5]=1; //ifncompt
#else
  params[5]=0; 
#endif  
#ifdef EVOLVEELECTRONS
  params[6]=1; //ifncompt
#else
  params[6]=0; 
#endif
 params[7]=0; //no relativistic electrons
#ifdef RELELECTRONS
  params[7]=1; //relativistic electrons
#endif

  solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,verbose,params,pp);

  exit(1);
  
}

//Jon test wrapper
int
test_jon_solve_implicit_lab()
{
  //NOGDET please

  FILE *in = fopen("jon.problem.pre","r");
  int i1,i2,iv,ifile;
  ldouble uu[NV],pp[NV],pp0[NV],dt;
  struct geometry geom;

  for(ifile=1;ifile<=165;ifile++)
    {
      printf("\n          &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
      printf("            &&&&&&&&&&&& case %4d &&&&&&&&&&&&&\n",ifile);
      printf("            &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n");

      for (i1=0;i1<NV;i1++)
	iv=fscanf(in,"%lf",&uu[i1]);
      for (i1=0;i1<NV;i1++)
	iv=fscanf(in,"%lf",&pp[i1]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  iv=fscanf(in,"%lf ",&geom.gg[i1][i2]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  iv=fscanf(in,"%lf ",&geom.GG[i1][i2]);
      iv=fscanf(in,"%lf ",&dt);
      iv=fscanf(in,"%lf ",&geom.alpha);
      iv=fscanf(in,"%lf ",&geom.gdet);
      
      //fill missing parts
      ldouble ucon[4];
      ucon[1]=pp[VX];
      ucon[2]=pp[VY];
      ucon[3]=pp[VZ];
      conv_vels(ucon,ucon,VEL4,VEL4,geom.gg,geom.GG);
      geom.alpha=sqrt(-1./geom.GG[0][0]);
      pp[5]=calc_Sfromu(pp[0],pp[1],geom.ix,geom.iy,geom.iz);

      //destroy magn field
      //uu[B1]=uu[B2]=uu[B3]=pp[B1]=pp[B2]=pp[B3]=0.;

      printf("\n...........................\nJon's input:\n\n");
      print_Nvector(uu,NV);
      print_Nvector(pp,NV);
      //print_metric(geom.gg);
      //print_metric(geom.GG);
      //printf("%e %e %e\n",dt,geom.alpha,geom.gdet);
      //printf("ut: %e\n",ucon[0]);
      int corr[2],fixup[2],u2pret,radcor;
     
      //printf("inverting...\n");
      u2pret=u2p_solver(uu,pp,&geom,U2P_HOT,0); //hd
      if(u2pret<0) printf("u2pret mhd: (%d)\n",u2pret);
      u2p_rad(uu,pp,&geom,&radcor); //rad
      if(radcor!=0) printf("u2pcor rad: (%d)\n",radcor);
      printf("\n..........................\nafter u2p_HOT:\n\n");
      print_Nvector(pp,NV);

      //compare change in entropy
      ucon[1]=pp[VX];
      ucon[2]=pp[VY];
      ucon[3]=pp[VZ];
      conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
      ldouble s1=exp(uu[ENTR]/ucon[0]/pp[RHO]);
      ldouble s2=exp(pp[ENTR]/pp[RHO]);

      printf("\n..........................\nchange in entropy:\n\n");
      printf("s(adv) | s(inv): %e | %e\n",s1,s2); 
     
      if(s2/s1 < 0.9 | u2pret<0.)
	{ 
	  printf("\n PROBLEM DETECTED IN ENTROPY OR U2P_HOT DID NOT SUCCEED!\n");
	  u2pret=u2p_solver(uu,pp,&geom,U2P_ENTROPY,0); //hd
	  if(u2pret<0) printf("u2pret mhd: (%d)\n",u2pret);
	  printf("\n..........................\nafter u2p_ENTROPY:\n\n");
	  print_Nvector(pp,NV);
	}      
      
      printf("\n..........................\nafter p2u:\n\n");
      p2u(pp,uu,&geom);
      for (i1=0;i1<NV;i1++)
	pp0[i1]=pp[i1];
      print_Nvector(uu,NV);
      print_Nvector(pp,NV);

      getchar();
   
      ldouble deltas[NRADVAR];
      int verbose=1;
      int params[4];
      
      solve_explicit_lab_core(uu,pp,&geom,dt,deltas,verbose);
      params[1]=RADIMPLICIT_ENERGYEQ;
      params[2]=RADIMPLICIT_LAB;
      params[3]=1;
      ldouble ppret[NV];
      solve_implicit_lab_4dprim(uu,pp,&geom,dt,verbose,params,ppret);

      getchar();
    }

  fclose(in);

  return 0;
}

///////////////////////////////////////////////////////////////
int test_calckappa()
{
  char bufor[50];
  sprintf(bufor,"kappatest.dat");
  fout1 = fopen(bufor, "w");

  int ix, iy,iz,iv;
  ldouble pp[NV];
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
	      for (iv=0;iv<NV;iv++){
		pp[iv] = get_u(p, iv, ix, iy, iz);
	      }
	      struct geometry geom, geomBL;
	      fill_geometry(ix,iy,iz,&geom);
	      fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	      struct opacities opac;
	      ldouble kappa=calc_kappa(pp,&geom,&opac);

	      ldouble Ti,Te;
	      ldouble Tgas=calc_PEQ_Teifrompp(pp,&Te,&Ti,ix,iy,iz);

	      ldouble r=geomBL.xx;
	      ldouble th=geomBL.yy;
	      ldouble ph=geomBL.zz;
	     
	     
	      //gas prims
	      fprintf(fout1,"%d %d %d ",ix,iy,iz); //(1-3)

	      fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph); //(4-6)
	      
	      fprintf(fout1,"%.5e %.5e ", rhoGU2CGS(pp[RHO]), Tgas);
	      
	      fprintf(fout1,"%.5e %.5e ", kappaGU2CGS(opac.kappaGasAbs/pp[RHO]), kappaGU2CGS(opac.kappaGasRoss/pp[RHO]));
	      
	      fprintf(fout1,"\n");
	      
		
	    }
	}
    }
  fflush(fout1);
  fclose(fout1);
  return 0;
}


//prints opacities
int
test_opacities()
{
  ldouble temp,rho,rhocgs,pp[NV],k1,k2,k3,k4,totEmissivity,kappa;
  struct geometry geom;
  fill_geometry(NX/2,0,0,&geom);

  rhocgs=MU_GAS*M_PROTON_CGS;
  rhocgs=1.029306e-24;
  rho=rhoCGS2GU(rhocgs);

  
  
  pp[RHO]=rho;

  temp = 9.900092e+03;

  pp[UU]= calc_PEQ_ufromTrho(temp,rho,geom.ix,geom.iy,geom.iz);
  
  for(temp=1.e3;temp<1.e13;temp*=1.01)
    {
      pp[UU]= calc_PEQ_ufromTrho(temp,rho,geom.ix,geom.iy,geom.iz);
      struct opacities opac;
      kappa=calc_kappa(pp,&geom,&opac);
      printf("%e %e\n",temp,kappa);
    }

  
  return 0;
}

//test heating prescription fit
int
test_heatfit()
{
  ldouble me = M_ELECTR_CGS;
  ldouble mp = M_PROTON_CGS;
  ldouble mi = mp*(HFRAC*1. + HEFRAC*4.);
  ldouble c1 = 0.92;
  ldouble Te,Ti,beta;

  Ti=1.e10;
  
  for(Te=1.e2;Te<1.e20;Te*=2.)
    {
    for(beta=1.e-15;beta<1.e15;beta*=2.)
      {
 

	ldouble c2 = 1.6*Te/Ti;
	if(Ti<Te) c2 = 1.2*Te/Ti; //non-continuous!
	ldouble c3 = 18. + 5.*log10(Ti/Te);
	if (Ti < Te) c3 = 18.;
  
	ldouble QpQeRatio;
	if(beta>BIG) 
	  {
	    QpQeRatio = c1*sqrt(mi*Ti/me/Te);
	  }
	else
	  QpQeRatio = c1*(c2*c2 + pow(beta,2. - 0.2*log10(Ti/Te) ) )/(c3*c3 + pow(beta, 2. - 0.2*log10(Ti/Te) ) )*sqrt(mi*Ti/me/Te)*exp(-1./beta);
 
	ldouble factorVEHeating = 1./(1. + QpQeRatio); //fe from Ressler et al. 2015
	ldouble delta=factorVEHeating;
	printf("%e %e %e\n",Te/Ti,beta,delta);
      }
  printf("\n");
}
  

  return 0;
}

//test Coulomb coupling solver
int
test_Ccoupling()
{
#ifdef RADIATION
#ifdef EVOLVEPHOTONNUMBER
#ifdef EVOLVEELECTRONS
  struct geometry geom;
  fill_geometry(NX/4,0,0,&geom);
  int i;

  ldouble pp[NV];
  ldouble rhocgs=1.e5;
  ldouble Tgas=1.e10;
  pp[RHO]=rhoCGS2GU(rhocgs);
  pp[UU]=calc_PEQ_ufromTrho(Tgas,pp[RHO],geom.ix,geom.iy,geom.iz);

  //nothing
  pp[VX]=pp[VY]=pp[VZ]=0.;
  pp[FX]=pp[FY]=pp[FZ]=0.;
  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],geom.ix,geom.iy,geom.iz);
  pp[B1]=pp[B2]=pp[B3]=0.*sqrt(1.*pp[RHO]);

  //negligible radiation field
  pp[EE]=calc_LTE_EfromT(Tgas)/1.e30;
  pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/calc_PEQ_Tfromurho(pp[UU],pp[RHO],geom.ix,geom.iy,geom.iz);

  ldouble Te=Tgas;
  ldouble Te1,Te2=-1.;
  Te2=Tgas;
  Te1=1.;
  Te=0.5*(Te1+Te2);

  while(fabs(Te2-Te1)>1.e-8*Te)
    {
      Te=0.5*(Te1+Te2);
      pp[ENTRE]=calc_SefromrhoT(pp[RHO],Te,ELECTRONS); //test_ccoupling

      //radiative four-force
      ldouble Gi[4],Giff[4],Hi[4]={0.,0.,0.,0.},Hiff[4]={0.,0.,0.,0.};

      //recalculate the ff time component of four-force
      calc_Gi(pp,&geom,Giff,0.0, 0,0);  // Ghat^0 = - kappa 4 Pi B + ... //
      ldouble CC=0.;
      CC=calc_CoulombCoupling(pp,&geom); //Coulomb coupling

      ldouble EA=Giff[0]; //emission - absorption, noe heating here but below
      ldouble Qe=CC+EA;

      //print_primitives(pp);
      ldouble Ti,Te;
      Tgas=calc_PEQ_Teifrompp(pp,&Te,&Ti,geom.ix,geom.iy,geom.iz);
      printf("Temperatures: %e %e %e     ",Tgas,Ti,Te);
      printf("CC: %e EA: %e\n",CC,EA);

      if(CC+EA>0.) 
	{
	  Te1=Te;
	}
      else
	Te2=Te;
    }
#endif
#endif
#endif
  
  return 0;
}

// test synchrotron temperature
int
test_Tsynch()
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
  
  pp[VX]=pp[VY]=pp[VZ]=0.;
  
  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],NX/4,0,0);
#ifdef MAGNFIELD
  pp[B1]=pp[B2]=pp[B3]=0.;
  pp[B1]=sqrt(2.*0.1*GAMMAM1*pp[UU]);
#endif
  
  //some radiation field
#ifdef RADIATION
  pp[FX]=pp[FY]=pp[FZ]=0.;
  pp[EE]=pp[UU];
  pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/calc_PEQ_Tfromurho(pp[UU],pp[RHO],NX/4,0,0);
  
  //electron entropy
#ifdef EVOLVEELECTRONS
  pp[ENTRE]=calc_S3fromrhoT(pp[RHO],Te,ELECTRONS);
  pp[ENTRI]=calc_S3fromrhoT(pp[RHO],Te,IONS);
#endif
  
#endif
  
  //magnetic pressure
  ldouble ucon[4], ucov[4], bcon[4], bcov[4], bsq, bsqcgs, Bmagcgs, pmagcgs;
  calc_ucon_ucov_from_prims(pp, &geom, ucon, ucov);
#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, &geom, bcon, bcov, &bsq);
  pmagcgs = endenGU2CGS(bsq/2.);
  Bmagcgs=sqrt(pmagcgs*8.*M_PI);
#endif
  
  printf("density: %e (GU) = %e (cgs, g/cm3)\n",pp[RHO],rhoGU2CGS(pp[RHO]));
  printf("rest mass en. den: %e (GU) = %e (cgs, erg/cm3)\n",pp[RHO],endenGU2CGS(pp[RHO]));
  printf("ptot: %e (GU) = %e (cgs, erg/cm3)\n",GAMMAM1*pp[UU],endenGU2CGS(GAMMAM1*pp[UU]));
#ifdef MAGNFIELD
  printf("pmag: %e (GU) = %e (cgs, erg/cm3)\n",bsq/2,pmagcgs);
  printf("B: %e (GU) = %e (cgs, G)\n", pp[B1], Bmagcgs);
#endif
  
  
  //Synchrotron emission
  
  
  ldouble k_boltCGS = K_BOLTZ_CGS;
  ldouble h_planckCGS = H_CGS;
  ldouble sigmaCGS = SIGMA_RAD_CGS;
  ldouble mpcgs=M_PROTON_CGS;
  
  ldouble BBenergy = 4.*sigmaCGS*Te*Te*Te*Te;
  
  ldouble B = SIGMA_RAD*pow(Te,4.)/Pi;
  
  
  ldouble rho=pp[RHO];
  ldouble  nu_MBsyn, zetaBsyngas, zetaBsynrad, zetaAdenomNum, zetaAdenomrad, IaByBrad, IaByBnum, emisSynchro;
  emisSynchro = 3.61e-34/mpcgs*0.5*(1. + HFRAC)*Te*Te*Bmagcgs*Bmagcgs; // standard emisivity divided by density
  
  ldouble Trad=Te;
  nu_MBsyn = 1.19e-13*Te*Te;
  zetaBsynrad = k_boltCGS*Trad/h_planckCGS/nu_MBsyn;
  ldouble zeta=zetaBsynrad/Bmagcgs;
  zetaAdenomrad =
  1.79*cbrt(zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*Bmagcgs*Bmagcgs*Bmagcgs*Bmagcgs)
  + 1.35*cbrt(zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*Bmagcgs*Bmagcgs)
  + 0.248*zetaBsynrad*zetaBsynrad*zetaBsynrad;
  
  IaByBrad = Bmagcgs*Bmagcgs/zetaAdenomrad;
  ldouble kapparadsyn = kappaCGS2GU( (2.13e39/mpcgs)*0.5*(1. + HFRAC)/Te/Te/Te/Te/Te*IaByBrad)*rho;
  ldouble kappagassyn = kappaCGS2GU(emisSynchro/BBenergy)*rho;
  
  //number-of-photons averaged opacities
  
  zetaAdenomNum =
  (0.326*cbrt(zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*Bmagcgs*Bmagcgs)
   + 0.146*zetaBsynrad*zetaBsynrad*cbrt(Bmagcgs) +
   0.015*cbrt(zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad*zetaBsynrad));
  
  IaByBnum = cbrt(Bmagcgs*Bmagcgs*Bmagcgs*Bmagcgs)/zetaAdenomNum;
  
  ldouble kapparadnumsyn = kappaCGS2GU( (2.13e39/mpcgs)*0.5*(1. + HFRAC)/Te/Te/Te/Te/Te*IaByBnum)*rho;
  
  // sum up the absorption opacities
  ldouble kappaGasAbs=kappagassyn;
  ldouble kappaRadAbs=kapparadsyn;
  ldouble kappaGasNum=0.;
  ldouble kappaRadNum=kapparadnumsyn;
  
  
  ldouble emission=fabs(kappaGasAbs*4.*Pi*B);
  
  //emission of photons
  ldouble ndotffSynch=0.;
  ldouble mue=2./(1.+HFRAC);
  ldouble ne=rhocgs/mue/M_PROTON_CGS; //number density of photons and electrons
  ldouble ndotffSynchCgs = 1.44e5*Bmagcgs*ne;
  ndotffSynch = ndotffSynchCgs/lenCGS2GU(1.)/lenCGS2GU(1.)/lenCGS2GU(1.)/timeCGS2GU(1.);
  //this is genereal, with or without synchro
  ldouble ndotff = ndotffSynch;
  
  
  //estimate the radiation temperature produced
  ldouble TradSynch=2.05e-18*Bmagcgs*Te*Te/1.e10/1.e10/2.70188/K_BOLTZ_CGS;
  
  
  ldouble TradSynchCode = emission/ndotff/2.70188/K_BOLTZ;
  
  
  //emission
  printf("============\n");
  
  ldouble emcgs=heatcoolGU2CGS(emission);
  ldouble ndotcgs=numdensGU2CGS(timeCGS2GU(ndotff));
  
  printf("emission: %e (GU) = %e (cgs, erg/cm3/s)\n",
         emission,emcgs);
  printf("photon gen. rate: %e (GU) = %e (cgs, 1/cm3/s)\n",
         ndotff,ndotcgs);
  
  
  printf("============\n");
  printf("Tsynch: %e (Ramesh estimate)\n",2.05e-18*Bmagcgs*(Te/1.e10)*(Te/1.e10)/(2.70188*K_BOLTZ_CGS));
  printf("Tsynch: %e (Maciek estimate)\n",1.82*1e-23*Bmagcgs*(Te)*(Te)/2.);
  printf("Tsynch: %e (code estimate)\n",TradSynch);
  printf("Tsynch: %e (code exact)\n",TradSynchCode);
  printf("Tsynch: %e (cgs ratio)\n",emcgs/ndotcgs/(2.70188*K_BOLTZ_CGS));
  
  printf("zeta: %e\n",zeta);
  
  return 0;
}
