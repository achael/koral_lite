ldouble uu[NV], pp[NV];

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

pp[UU]=RHO_INIT*K_BOLTZ*TGAS_INIT/(M_PROTON*MU_GAS*(GAMMA-1.));
//pp[UU]=calc_PEQ_ugasfromrhoTei(RHO_INIT,TE_INIT,TI_INIT,ix,iy,iz);
pp[VX]=pp[VY]=pp[VZ]=0.;
pp[VX] = VEL;
if(geom.xx<0.5)
  {
pp[RHO]=RHO_INIT;
  }
else
  {
pp[RHO]=RHO_FAC*RHO_INIT;
  }

ldouble mach = VEL/sqrt(GAMMA*(GAMMA-1)*pp[UU]/pp[RHO]);
if (ix==1)
  printf("%d MACH %e \n", ix, mach);
if (ix==(NX-1))
  printf("%d MACH %e \n", ix, mach);

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

#ifdef EVOLVEELECTRONS
pp[ENTRE]=calc_Sefromrhou(pp[RHO],pp[UU]*ELEC_U_FRAC,ELECTRONS);
pp[ENTRI]=calc_Sefromrhou(pp[RHO],pp[UU]*(1-ELEC_U_FRAC),IONS);
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
