ldouble uu[NV], pp[NV];

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

pp[RHO]=RHO_INIT;
pp[UU]=calc_PEQ_ugasfromrhoTei(RHO_INIT,TE_INIT,TI_INIT,ix,iy,iz);
pp[VX]=pp[VY]=pp[VZ]=0.;

if(geom.xx<0.5) pp[VX]=VEL;
 else pp[VX]=-VEL;

ldouble cs = sqrt(GAMMA*(GAMMA-1)*pp[UU]/pp[RHO]);
ldouble mach = VEL/cs;
if (ix==NX/2){
  //printf("%d %d %d %d %d %d\n",RHO, UU, VX, VY, VZ, ENTR);
  printf("%e %e\n",M_ELECTR_CGS,pp[RHO]);
  printf("%d CS: %.5e MACH: %e \n", ix, cs, mach);
 }

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

#ifdef EVOLVEELECTRONS
pp[ENTRE]=calc_SefromrhoT(RHO_INIT,TE_INIT/TEINITFACTOR,ELECTRONS);
pp[ENTRI]=calc_SefromrhoT(RHO_INIT,TI_INIT,IONS);
//pp[UU] = (calc_ufromSerho(pp[ENTRE],RHO_INIT,ELECTRONS,ix,iy,iz) + calc_ufromSerho(pp[ENTRI],RHO_INIT,IONS,ix,iy,iz));
#endif
//entropy
pp[5]=calc_Sfromu(pp[0],pp[1],geom.ix,geom.iy,geom.iz);

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
update_entropy_cell(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);
