ldouble uu[NV], pp[NV];

//last time when recalculated perturbation
global_double1=-1.;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

ldouble gamma=GAMMA;
#ifdef CONSISTENTGAMMA
#ifdef TEINIT
gamma=calc_gammaintfromTei(TEINIT,TIINIT);
#endif
//printf("%f\n",gamma);
#ifdef FORCEGAMMAGASFIXED
set_u_scalar(gammagas,ix,iy,iz,GAMMA);
#else
set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif
#endif

//Calculate initial temperatures from sound speed
pp[RHO]=RHO_INIT;
ldouble cs0=CS0_INIT; 
//cs=sqrt(gamma*p/rho)
//cs2 rho = gamma*p
//cs2 rho / gamma / gamma+1 = u
pp[UU]=cs0*cs0*pp[RHO]/gamma/(gamma-1);
ldouble Teinit,Tiinit;

#ifdef TEINIT
Teinit=TEINIT;
Tiinit=TIINIT;
pp[UU]=calc_PEQ_ugasfromrhoTei(pp[RHO],TEINIT,TIINIT,ix,iy,iz); //placeholder - need to correct for relel later
#endif

pp[VX]=VX_INIT;
pp[VY]=VY_INIT;
pp[VZ]=VZ_INIT;

#ifdef VERTCOLLISION
if(iy>TNY/2.) pp[VY]=-fabs(pp[VY]);
 else pp[VY]=fabs(pp[VY]);
#endif

#ifdef DBLVERTCOLLISION
if(ix>TNX/2.) pp[VX]=-fabs(pp[VX]);
 else pp[VX]=fabs(pp[VX]);
#endif

#ifdef HUBBLEBACKGROUND
pp[VX]=-HUBBLEMAGN*(geom.xx-0.5);
#endif

#ifdef PERT_INIT
 pp[VX]+=VPERT*((double)rand()/(double)RAND_MAX-0.5);
 pp[VY]+=VPERT*((double)rand()/(double)RAND_MAX-0.5);
#endif

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
pp[B1]=BX_INIT;
pp[B2]=BY_INIT;
pp[B3]=BZ_INIT;

#ifdef PERT_SPECTRUM
pp[B1]=sqrt(pp[UU]*(gamma-1.)/BETA_INIT);
pp[B2]=pp[B3]=0.;
#endif
#endif

#ifdef RADIATION
pp[EE]=calc_ncompt_Ehatrad(TR_INIT,NPH_INIT);
pp[FX]=pp[FY]=pp[FZ]=0.;
   
#ifdef EVOLVEPHOTONNUMBER
pp[NF]=NPH_INIT;
#endif
#endif 

//#ifdef RADIATION
//pp[EE]=EE_INIT;
//pp[FX]=pp[FY]=pp[FZ]=0.;
//#endif
   
#ifdef EVOLVEELECTRONS
ldouble ue=1./100.*pp[UU];
ldouble ui=(1.-1./100.)*pp[UU];

pp[ENTRE]=calc_Sefromrhou(pp[RHO],ue,ELECTRONS);
pp[ENTRI]=calc_Sefromrhou(pp[RHO],ui,IONS);
//pp[ENTRE]=calc_SefromrhoT(RHO_INIT,Teinit,ELECTRONS);
//pp[ENTRI]=calc_SefromrhoT(RHO_INIT,Tiinit,IONS);

#ifdef CONSISTENTGAMMA
//ldouble Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz); //ANDREW gas temperature if NO relel
//Teinit=Tiinit=Tgas;
//pp[UU]=calc_PEQ_ugasfromrhoTei(pp[RHO],Teinit,Tiinit,ix,iy,iz); //placeholder - need to correct for relel later
//gamma=calc_gammaintfromTei(Teinit,Tiinit); //good faith estimate
//set_u_scalar(gammagas,ix,iy,iz,gamma);

Tiinit=solve_Teifromnmu(pp[RHO]/MU_I/M_PROTON, M_PROTON, ui,IONS); //solves in parallel for gamma and temperature
Teinit=solve_Teifromnmu(pp[RHO]/MU_E/M_PROTON, M_ELECTR, ue,ELECTRONS); //solves in parallel for gamma and temperature
gamma=calc_gammaintfromTei(Teinit,Tiinit); //good faith estimate
set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif

#ifdef RELELECTRONS
int ie;
for (ie=0; ie < NRELBIN; ie++){
   pp[NEREL(ie)] = RELEL_NRELEL_INIT;
 }

ldouble neur=calc_relel_ne(pp);
ldouble neth=calc_thermal_ne(pp);
if(neur/(pp[RHO]/M_PROTON/MU_E) > MAX_RELEL_FRAC_N){
  printf("neur/netot > MAX_RELEL_FRAC_N !!");
  getch();
 }

//recompute gammagas and total u with relel contribution
gamma = calc_gammaint_relel(pp, Teinit, Tiinit);
pp[UU]=calc_PEQ_ugasfrom_Tei_relel(pp, Teinit,Tiinit); 

//Recompute entropies and total gas energy with correct entropy


pp[ENTRE]=calc_SefromrhoT(neth*MU_E*M_PROTON, Teinit, ELECTRONS);
pp[ENTRI]=calc_SefromrhoT(pp[RHO], Tiinit, IONS);

//ldouble ueth=calc_ufromSerho(pp[ENTRE],neth*MU_E*M_PROTON,ELECTRONS,ix,iy,iz);
//ldouble uith=calc_ufromSerho(pp[ENTRI],pp[RHO],IONS,ix,iy,iz);
//ldouble uur=calc_relel_uint(pp);
//pp[UU] = ueth+uith+uur;
#endif
#endif

//printf("ENTRE: %e ENTRI %e \n", pp[ENTRE], pp[ENTRI]);

//printf("T: %e %e \n", Teinit, Tiinit);
//printf("U: %e %e %e %e \n", ueth, uith, uur, pp[UU]);
//getch();

//Total Gas Entropy
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);

//set gammagas!
set_u_scalar(gammagas,ix,iy,iz,gamma);

p2u(pp,uu,&geom);



/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//print_primitives(pp); getch();
set_cflag(0,ix,iy,iz,0);
