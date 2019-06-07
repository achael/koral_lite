ldouble uu[NV], pp[NV];

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

pp[RHO]=RHO_INIT;

// with RELEL
pp[UU]=calc_PEQ_ugasfromrhoTei(RHO_INIT,TE_INIT,TI_INIT,ix,iy,iz); //placeholder - need to correct for relel
pp[VX]=pp[VY]=pp[VZ]=0.;
#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

#ifdef RADIATION
pp[EE]=calc_ncompt_Ehatrad(TR_INIT, NPH_INIT);
pp[FX]=pp[FY]=pp[FZ]=0.;
   
#ifdef EVOLVEPHOTONNUMBER
pp[NF]=NPH_INIT;
#endif
#endif 

#ifdef EVOLVEELECTRONS

#ifdef RELELECTRONS
int ie;
for (ie=0; ie < NRELBIN; ie++){
   pp[NEREL(ie)] = RELEL_NRELEL;
 }
#endif

ldouble neth=calc_thermal_ne(pp);
if(neth<0){
  printf("neth < 0 - change RELEL_NRELEL or RHO_INIT!");
  getch();
 }

//#ifdef RELELENTROPY
//pp[ENTRE]=calc_S4fromnT(neth, TE_INIT, ELECTRONS);
//pp[ENTRI]=calc_S4fromnT(pp[RHO]/MU_I/M_PROTON, TI_INIT, IONS);
//#else
pp[ENTRE]=calc_SefromrhoT(neth*MU_E*M_PROTON,TE_INIT,ELECTRONS);
pp[ENTRI]=calc_SefromrhoT(pp[RHO],TI_INIT,IONS);
//#endif

ldouble ueth=calc_ufromSerho(pp[ENTRE], neth*MU_E*M_PROTON, ELECTRONS,ix,iy,iz);
ldouble uith=calc_ufromSerho(pp[ENTRI], pp[RHO], IONS,ix,iy,iz);
ldouble uur=calc_relel_uint(pp);

pp[UU] = ueth+uith+uur;
printf("%i %e %e %e %e \n",ix, endenGU2CGS(ueth), endenGU2CGS(uith), endenGU2CGS(uur),  (ueth+uith+uur-pp[UU])/pp[UU]);


//set gammagas
ldouble gamma=calc_gammagas(pp,ix,iy,iz); //good faith estimate
set_u_scalar(gammagas,ix,iy,iz,gamma);
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

