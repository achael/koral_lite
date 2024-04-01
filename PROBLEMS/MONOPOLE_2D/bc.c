//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,
//	ldouble *uu,ldouble *pp,int ifinit,int BCtype)

#define BCCOORDS KSCOORDS

/**********************/
//geometries
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv;  	  

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BCCOORDS);

ldouble gg[4][5],GG[4][5];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

/**********************/

//radius
  if(BCtype==XBCHI) //outflow in magn, atm in rad., atm. in HD
  {
    iix=NX-1;
    iiy=iy;
    iiz=iz;
    
    //copying everything
    //   for(iv=0;iv<NV;iv++)
    // {
    //	pp[iv]=get_u(p,iv,iix,iiy,iiz);
    //  }
   //copying everything
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
	ppback[iv]=pp[iv];
      }

 
    //ANDREW TODO what to do about outer BC
    /*
    //!! begin rescale
    //first transform to BL
    trans_pmhd_coco(pp, pp, MYCOORDS,BCCOORDS, geom.xxvec,&geom,&geomBL);

    // boundary cell and scale factors
    struct geometry geombdBL;
    fill_geometry_arb(iix,iiy,iiz,&geombdBL,BCCOORDS);
    ldouble rghost = geomBL.xx;
    ldouble rbound = geombdBL.xx;
    ldouble deltar = rbound - rghost;
    ldouble scale1 =  rbound/rghost; //r^-1
    ldouble scale2 =  rbound*rbound/rghost/rghost; //r^-2
    ldouble scale3 =  rbound*rbound*rbound/rghost/rghost/rghost; //r^-3 -- inital atm

    //  bsq = br*br + r^2bth^2 + r^2*sin^2(th)*bph^2
      
    pp[RHO]*=scale3;//scale2;
    pp[UU] *=scale3;//scale2;
    #ifdef MAGNFIELD
    pp[B1] *=scale2;
    pp[B2] *=scale1;
    pp[B3] *=scale1;
    #endif
    pp[VX] *=1.;
    pp[VY] *=scale1;
    pp[VZ] *=scale1;
    #ifdef RADIATION
    pp[EE0] *=scale2;
    pp[FX0] *=1.;
    pp[FY0] *=scale2;
    pp[FZ0] *=scale2;
    #endif
    //transform back after rescaling
    trans_pmhd_coco(pp, pp,BCCOORDS, MYCOORDS, geomBL.xxvec,&geomBL, &geom);
    */
    /*
    //!! end rescale
    /*
    //check for gas inflow
    ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};    
    conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
    trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,BCCOORDS);
    if(ucon[1]<0.) //inflow, reset to zero velocity
      {
	//atmosphere in rho,uint and velocities and zero magn. field
	//set_hdatmosphere(pp,xxvec,gg,GG,4);
	ucon[1]=0.;
	trans2_coco(geomBL.xxvec,ucon,ucon,BCCOORDS,MYCOORDS);
	conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	pp[VX]=ucon[1];
	pp[VY]=ucon[2];
	pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }
    */
    
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
       }

     p2u(pp,uu,&geom);
     return 0;
   }

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

