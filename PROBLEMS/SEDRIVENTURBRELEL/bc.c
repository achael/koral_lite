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


//periodic in x:
iiz=iz;
iiy=iy;
iix=ix;
if(BCtype==XBCLO) iix=ix+NX;
if(BCtype==XBCHI) iix=ix-NX;

   

//periodic in theta:
if(BCtype==YBCLO) iiy=iy+NY;
if(BCtype==YBCHI) iiy=iy-NY;

   
//periodic in phi:
if(BCtype==ZBCLO) iiz=iz+NZ;
if(BCtype==ZBCHI) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

#ifdef HUBBLEBACKGROUND
pp[VX]=-HUBBLEMAGN*(geom.xx-0.5);
#endif

p2u(pp,uu,&geom);

//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

