//returns problem specific BC
//int calc_bc(int ix,int iy,int iz,ldouble t,
//	ldouble *uu,ldouble *pp,int ifinit,int BCtype)

/**********************/
//geometries
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv;  	  

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

/**********************/

#ifdef COPY_XBC //simple inflow/outflow
if(BCtype==XBCHI || BCtype==XBCLO)
{
  iiy=iy;
  iiz=iz;
  if(BCtype==XBCHI) iix=NX-1;
  else if(BCtype==XBCLO) iix=0;
  
  for(iv=0;iv<NV;iv++)
  {
    pp[iv]=get_u(p,iv,iix,iiy,iiz);
  }
  p2u(pp,uu,&geom);
  return 0;
}
#else

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

ldouble gg[4][5],GG[4][5],ggsrc[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

//radius
if(BCtype==XBCHI) //outflow in magn, atm in rad., atm. in HD
                  //this can take a long time bc of coord transforms!
  {
    iix=NX-1;
    iiy=iy;
    iiz=iz;
    
    //copying everything
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

    //!! begin rescale

    // AA -- example of how to use precomputed my2out and out2my for faster coco
    //first transform to BL
    #ifdef PRECOMPUTE_MY2OUT
    trans_pmhd_coco_my2out(pp,pp,&geom,&geomBL);
    #else      
    trans_pmhd_coco(pp, pp, MYCOORDS,BLCOORDS, geom.xxvec,&geom,&geomBL);
    #endif

    
    struct geometry geombdBL;
    fill_geometry_arb(iix,iiy,iiz,&geombdBL,BLCOORDS);
    ldouble rghost = geomBL.xx;
    ldouble rbound = geombdBL.xx;
    ldouble deltar = rbound - rghost;
    ldouble scale1 =  rbound*rbound/rghost/rghost;
    ldouble scale2 = rbound/rghost;

    pp[RHO]*=scale1;
    pp[UU] *=scale1;
    #ifdef MAGNFIELD
    pp[B1] *=scale1;
    pp[B2] *=scale2;
    pp[B3] *=scale2;
    #endif
    //pp[VX] *=1.;
    pp[VY] *=scale1;
    pp[VZ] *=scale1;
    #ifdef RADIATION
    pp[EE0] *=scale1;
    //pp[FX0] *=1.;
    pp[FY0] *=scale1;
    pp[FZ0] *=scale1;
    #endif
    //transform back after rescaling

    #ifdef PRECOMPUTE_MY2OUT
    trans_pmhd_coco_out2my(pp,pp,&geomBL,&geom);
    #else      
    trans_pmhd_coco(pp, pp, BLCOORDS,MYCOORDS, geomBL.xxvec,&geomBL,&geom);
    #endif
    
    //!! end rescale

//checking for the gas inflow
    ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]}; 
    conv_vels(ucon, ucon, VELPRIM, VEL4, geom.gg, geom.GG);
    if(MYCOORDS!=CYLCOORDS)
      {
      #ifdef PRECOMPUTE_MY2OUT
      trans2_coco_my2out(ucon, ucon, geom.ix, geom.iy, geom.iz);
      #else      
      trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
      #endif
      }
    if(ucon[1]<0.) //inflow, resseting to atmosphere
      {
	//atmosphere in rho,uint and velocities and zero magn. field
	//set_hdatmosphere(pp,xxvec,gg,GG,4);
	//set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);

	//!! begin reset to floor        
	//	pp[UU] = pp[RHO]*UURHORATIOMIN*3.;	
	//!! end reset to floor

	ucon[1]=0.;
	/*
	#ifdef MAGNFIELD
	pp[B2]=pp[B3]=0.;
	#endif
	*/
	if(MYCOORDS!=CYLCOORDS)
	{
         #ifdef PRECOMPUTE_MY2OUT
	 trans2_coco_out2my(ucon, ucon, geomBL.ix, geomBL.iy, geomBL.iz);
         #else      
	 trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
         #endif
        }
	conv_vels_ut(ucon, ucon, VEL4, VELPRIM, geom.gg, geom.GG); // AA -- we know ut after coco
	//conv_vels(ucon, ucon, VEL4, VELPRIM, geom.gg, geom.GG);
	
	pp[VX]=ucon[1];
	pp[VY]=ucon[2];
	pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }
    
#ifdef RADIATION
    ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};    
    conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);

    #ifdef PRECOMPUTE_MY2OUT
    trans2_coco_my2out(urfcon, urfcon, geom.ix, geom.iy, geom.iz);
    #else      
    trans2_coco(geom.xxvec,urfcon,urfcon,MYCOORDS,BLCOORDS);
    #endif
    
    if(urfcon[1]<0.) //inflow, resseting to atmosphere
      {
	//atmosphere in radiation
	//set_radatmosphere(pp,xxvec,gg,GG,0);
	urfcon[1]=0.;

        #ifdef PRECOMPUTE_MY2OUT
        trans2_coco_out2my(urfcon, urfcon, geomBL.ix, geomBL.iy, geomBL.iz);
        #else      
	trans2_coco(geomBL.xxvec,urfcon,urfcon,BLCOORDS,MYCOORDS);
        #endif
	conv_vels_ut(urfcon,urfcon,VEL4,VELPRIMRAD,geom.gg,geom.GG); // AA -- we know ut after coco
	//conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,geom.gg,geom.GG);
	
	pp[FX0]=urfcon[1];
	pp[FY0]=urfcon[2];
	pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }
#endif

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
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       }

     p2u(pp,uu,&geom);
     return 0;
   }
#endif

//reflections/outflow in theta 
//in 3D polar cells overwritten with #define CORRECT_POLARAXIS_3D
if(BCtype==YBCLO) //upper spin axis 
  {      
    
    iiy=-iy-1;
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	//v_theta
	if(iv==VY || iv==B2 || iv==FY0)
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
	//v_theta	
	if(iv==VY || iv==B2 || iv==FY0)
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
	
      }
 
    p2u(pp,uu,&geom); 
    return 0; 
  }

if(BCtype==ZBCHI || BCtype==ZBCLO) //periodic in phi:
{
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
  return 0;
}
//testing if interpolated primitives make sense
//check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section
 
return 0;

