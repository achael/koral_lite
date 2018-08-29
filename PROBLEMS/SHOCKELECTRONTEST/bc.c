//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,
//	ldouble *uu,ldouble *pp,int ifinit,int BCtype)

/**********************/
//geometries
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv;  	  

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

/**********************/

//radius
if(BCtype==XBCHI || BCtype==XBCLO) //outflow in magn, atm in rad., atm. in HD
  {
    
pp[RHO]=RHO_INIT;
pp[UU]=calc_PEQ_ugasfromrhoTei(RHO_INIT,TE_INIT,TI_INIT,ix,iy,iz);
pp[VX]=pp[VY]=pp[VZ]=0.;

if(geom.xx<0.5) pp[VX]=VEL;
 else pp[VX]=-VEL;

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

  
#ifdef EVOLVEELECTRONS
pp[ENTRE]=calc_S2fromrhoT(RHO_INIT,TE_INIT/TEINITFACTOR,ELECTRONS);
pp[ENTRI]=calc_S2fromrhoT(RHO_INIT,TI_INIT,IONS);
#endif

//entropy
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);

//to conserved
p2u(pp,uu,&geom);

    return 0;  
  }

  
//periodic in phi:
iiz=iz;
iiy=iy;
iix=ix;
if(BCtype==ZBCLO) iiz=iz+NZ;
if(BCtype==ZBCHI) iiz=iz-NZ;
if(BCtype==YBCLO) iiy=iy+NY;
if(BCtype==YBCHI) iiy=iy-NY;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

