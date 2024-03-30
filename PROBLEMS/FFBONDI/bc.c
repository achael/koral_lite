// TODO analytics

//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,
//	ldouble *uu,ldouble *pp,int ifinit,int BCtype)

/**********************/
//geometries
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv;  	  

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

ldouble gg[4][5],GG[4][5];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

/**********************/

//radius

// boundary is analytic solution

if(BCtype==XBCHI || BCtype==XBCLO)
{
   
  init_bondi_full(pp, ix, iy, iz);
  p2u(pp,uu,&geom);
  return 0;
}


/*
  if(BCtype==XBCHI) //outflow at outer radius
  {
    
    iix=NX-1;
    iiy=iy;
    iiz=iz;
    
   //copying everything
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
	ppback[iv]=pp[iv];
      }
    
    p2u(pp,uu,&geom);
    return 0;  
  }


 else if(BCtype==XBCLO) //outflow near BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

      for(iv=0;iv<NV;iv++)
       {
	 pp[iv]=get_u(p,iv,iix,iiy,iiz);
	 ppback[iv] = pp[iv];
       }

     p2u(pp,uu,&geom);
     return 0;
   }
*/

//reflections/outflow in theta 
//in 3D polar cells overwritten with #define CORRECT_POLARAXIS_3D
if(BCtype==YBCLO) //upper spin axis 
  {      
    iiy=-1*iy-1; 
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	//v_theta
#ifdef FORCEFREE
        if(iv==VY || iv==B2 || iv==VYFF || iv==FY0)
#else
	if(iv==VY || iv==B2 || iv==FY0)
#endif
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }
     

    p2u(pp,uu,&geom);
    return 0;
  }

if(BCtype==YBCHI) //lower spin axis
  {
    iiy=NY-(iy-NY)-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
  	  
    for(iv=0;iv<NV;iv++)
      {
#ifdef FORCEFREE
        if(iv==VY || iv==B2 || iv==VYFF || iv==FY0)
#else
	if(iv==VY || iv==B2 || iv==FY0)
#endif
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	
      }
 
    p2u(pp,uu,&geom); 
    return 0; 
  }
   
//periodic in phi:
iiz=iz;
iiy=iy;
iix=ix;
if(BCtype==ZBCLO) iiz=iz+NZ;
if(BCtype==ZBCHI) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

