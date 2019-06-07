ldouble uu[NV], pp[NV];

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);


pp[RHO]=RHO_INIT;
pp[UU]=calc_PEQ_ugasfromrhoTei(RHO_INIT,TE_INIT,TI_INIT,ix,iy,iz);
pp[VX]=pp[VY]=pp[VZ]=0.;

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

#ifdef RADIATION
pp[EE]=calc_ncompt_Ehatrad(TR_INIT,NPH_INIT);
pp[FX]=pp[FY]=pp[FZ]=0.;
   
#ifdef NCOMPTONIZATION
pp[NF]=NPH_INIT;
#endif
#endif 

#ifdef EVOLVEELECTRONS
pp[ENTRE]=calc_SefromrhoT(RHO_INIT,TE_INIT,ELECTRONS);
pp[ENTRI]=calc_SefromrhoT(RHO_INIT,TI_INIT,IONS);

#ifdef RELELECTRONS
int ie;
ldouble g;
ldouble b;
//!AC this is an UNNORMALIZED maxwell-juttner distribution
//printf("NBIN: %d \n", NRELBIN);

for (ie=0; ie < NRELBIN; ie++){
   pp[NEREL(ie)] = RELEL_NRELEL;
 }

#ifdef RELEL_INIT_MAXWELL
for (ie=0; ie < NRELBIN; ie++)
{
   g = relel_gammas[ie];
   b = pow((1. - 1./(g*g)), 0.5);
   pp[NEREL(ie)] = numdensCGS2GU(RELEL_INIT_NORM) * g * g * b * exp(-g/RELEL_INIT_THETA);
}
#endif
#ifdef RELEL_INIT_PLAW
   ldouble pref = numdensCGS2GU(RELEL_INIT_NORM)*(RELEL_INIT_INDEX-1.0)/(pow(RELEL_INIT_MIN,1.0-RELEL_INIT_INDEX)-pow(RELEL_INIT_MAX,1.0-RELEL_INIT_INDEX));
  
ldouble alpha=1.e-5;
ldouble beta=1.e-22;
ldouble g0=RELEL_INIT_MIN - RELEL_INIT_INDEX/(2*alpha*RELEL_INIT_MIN);
ldouble g1=RELEL_INIT_MAX - RELEL_INIT_INDEX/(2*beta*RELEL_INIT_MAX);
ldouble aa=pref*pow(RELEL_INIT_MIN,-1*RELEL_INIT_INDEX);
ldouble bb=pref*pow(RELEL_INIT_MAX,-1*RELEL_INIT_INDEX);
printf("%e %e \n",aa,bb);

for (ie=0; ie < NRELBIN; ie++)
{
   g = relel_gammas[ie];
   if (g <= RELEL_INIT_MAX && g >= RELEL_INIT_MIN) pp[NEREL(ie)] = pref* pow(g, -1.*RELEL_INIT_INDEX);
   else if(g<RELEL_INIT_MIN) pp[NEREL(ie)] = aa*exp(-alpha*(g-RELEL_INIT_MIN)*(g-2*g0+RELEL_INIT_MIN));
   else if(g>RELEL_INIT_MAX) pp[NEREL(ie)] = bb*exp(-beta*(g-RELEL_INIT_MAX)*(g-2*g1+RELEL_INIT_MAX));
   else pp[NEREL(ie)] = 1.e-15;
}

//bump at end
//pp[NEREL(NRELBIN-1)] = 10000*pref*pow(g,-1*RELEL_INIT_INDEX);
//flat at end
#endif

ldouble neth = calc_thermal_ne(pp);
if(neth<0)
{
  printf("neth < 0 - change RELEL_NRELEL or RHO_INIT!");
  getch();
 }

//Recompute entropies
pp[ENTRE]=calc_SefromrhoT(neth*MU_E*M_PROTON,TE_INIT,ELECTRONS);
pp[ENTRI]=calc_SefromrhoT(RHO_INIT,TI_INIT,IONS);
#endif
#endif

//entropy
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);

//to conserved
p2u(pp,uu,&geom);



/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);
