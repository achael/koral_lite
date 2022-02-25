//sets the initial conditions
//called from a loop going over ix,iy,iz

/***********************************************/
//structure of geometry
//holds coordinates and metric, see ko.h
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);


/***********************************************/
//vectors of primitives and conserved
ldouble pp[NV],uu[NV];

/***********/
pp[RHO]=RHOINIT;
pp[UU]=UUINIT;
pp[VX]=pp[VZ]=0.;
pp[VY]=VGASINIT;

pp[B1]=pp[B2]=pp[B3]=0.;

pp[EE]=ERADINIT;
pp[FX]=pp[FZ]=0.;
if(geom.xx<-3./4.*LLL || geom.xx>3./4.*LLL)
  pp[FY]=0.;
if(geom.xx<1./4.*LLL && geom.xx>-1./4.*LLL)
  pp[FY]=-VRADINIT;
if(geom.xx>-3./4.*LLL && geom.xx<-1./4*LLL)
  pp[FY]=-VRADINIT*(geom.xx+3.*LLL/4.)/(LLL/2.);
if(geom.xx>1./4.*LLL && geom.xx<3./4.*LLL)
  pp[FY]=VRADINIT*(geom.xx-3.*LLL/4.)/(LLL/2.);

//calculate entropy from rho & uint
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);

/***********************************************/
//convert primitives to conserved
p2u(pp,uu,&geom);
/***********************************************/

//save to memory
int iv;
for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

