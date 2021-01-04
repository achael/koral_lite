int diskatboundary(ldouble *pp, void *ggg, void *gggBL);
int if_indisk_bc_check(int ix, int iy, int iz);

/**********************/

//definitions
ldouble gdet_src,gdet_bc,Fx,Fy,Fz,ppback[NV];
int iix,iiy,iiz,iv;  	  
ldouble xxvec[4],xxvecBL[4],xx,yy,zz;
#ifdef SPECIAL_BC_CHECK
ldouble pp_x,pp_y;
#endif

//coordinates
get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

/**********************/

//outer edge, outflows with velocity check
if(BCtype==XBCHI)
  {

#ifdef SPECIAL_BC_CHECK
    int bc_check = if_indisk_bc_check(ix+TOI,iy+TOJ,iz+TOK);
    int disk_bc_check = diskatboundary(pp, &geom, &geomBL);
    if(ix > NX) bc_check = 0; //ignore bc check value if cell is a gc within tile as it is
    if(disk_bc_check < 0 && bc_check == 0)
#else
    if(diskatboundary(pp, &geom, &geomBL)<0)
#endif
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
    //first transform to BL
    trans_pmhd_coco(pp, pp, MYCOORDS,BLCOORDS, geom.xxvec,&geom,&geomBL);
    
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
    trans_pmhd_coco(pp, pp,BLCOORDS, MYCOORDS, geomBL.xxvec,&geomBL, &geom);
    //!! end rescale

    //checking for the gas inflow
    ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]};    
    conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
    if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,ucon,ucon,MYCOORDS,BLCOORDS);
    if(ucon[1]<0.) //inflow, resseting to atmosphere
      {
	//atmosphere in rho,uint and velocities and zero magn. field
	//set_hdatmosphere(pp,xxvec,gg,GG,4);
	//set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);

	//!! begin reset to floor        
	//	pp[UU] = pp[RHO]*UURHORATIOMIN*3.;
	
	//!! end reset to floor

	ucon[1]=0.;
	/**
	#ifdef MAGNFIELD
	pp[B2]=pp[B3]=0.;
	#endif
	**/
	if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
	conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	pp[VX]=ucon[1];
	pp[VY]=ucon[2];
	pp[VZ]=ucon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }
    
#ifdef RADIATION
    ldouble urfcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};    
    conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom.gg,geom.GG);
    if(MYCOORDS!=CYLCOORDS) trans2_coco(geom.xxvec,urfcon,urfcon,MYCOORDS,BLCOORDS);
    if(urfcon[1]<0.) //inflow, resseting to atmosphere
      {
	//!! begin reset to floor
	//pp[EE0] = ERADATMMIN;
	//!! end reset to floor

	//atmosphere in radiation
	//set_radatmosphere(pp,xxvec,gg,GG,0);
	urfcon[1]=0.;
	if(MYCOORDS!=CYLCOORDS) trans2_coco(geomBL.xxvec,urfcon,urfcon,BLCOORDS,MYCOORDS);
	conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,geom.gg,geom.GG);
	pp[FX0]=urfcon[1];
	pp[FY0]=urfcon[2];
	pp[FZ0]=urfcon[3];//atmosphere in rho,uint and velocities and zero magn. field
      }
#endif

      }
#ifdef SPECIAL_BC_CHECK 
    else if(bc_check == -1 || bc_check == 1)
      {
        //printf("Using special BC. \n");
        if(disk_bc_check < 0)
        { 
          for(iv=0;iv<NV;iv++)
          {
            if( (ix+TOI) == STREAM_IX )
            {
  	      //r component
              #ifdef MAGNFIELD
	      if(iv==VX || iv==B1 || iv==FX0)
              #else
	      if(iv==VX || iv==FX0)
              #endif
              {
                iix = ix - 1;
                iiy = iy; // iy +/- 1
                iiz = iz;
	        pp[iv]=-get_u(p,iv,iix,iiy,iiz);
              }
              #ifdef MAGNFIELD
	      else if(iv==VY || iv==B2 || iv==FY0) //theta component
              #else
	      else if(iv==VY || iv==FY0)
              #endif
              {
                iix = ix;
                iiy = iy + bc_check; // iy +/- 1
                iiz = iz;
	        pp[iv]=-get_u(p,iv,iix,iiy,iiz);
              }
	      else //average other quantities between r/theta
              {
                iix = ix - 1;
                iiy = iy + bc_check; // iy +/- 1
                pp_x = get_u(p,iv,iix,iy,iz);
	        pp_y = get_u(p,iv,ix,iiy,iz);
                pp[iv] = (pp_x + pp_y)/2.;
              }
            }
            else if( (ix+TOI) == (STREAM_IX+1) ) //2nd gc, use reflecting BC
            {
  	      if(iv==VY || iv==B2 || iv==FY0)
              {
                iix = ix;
                iiy = iy + bc_check; // iy +/- 1
                iiz = iz;
	        pp[iv]=-get_u(p,iv,iix,iiy,iiz);
              }
	      else
              {
                iix = ix;
                iiy = iy + bc_check; // iy +/- 1
                iiz = iz;
	        pp[iv]=get_u(p,iv,iix,iiy,iiz);
              }
            }
          }
        }
      }
#endif
    p2u(pp,uu,&geom);

    /*
      print_primitives(pp);
      print_conserved(uu);getchar();
    */

    return 0;  
  }

 else if(BCtype==XBCLO) //outflow inside BH
   {
     iix=0;
     iiy=iy;
     iiz=iz;

     //linear extrapolation
      for(iv=0;iv<NV;iv++)
       {
	 pp[iv]=get_u(p,iv,0,iiy,iiz);
       }

      if(RMIN>rhorizonBL) //boundary outside BH
	{
	  //checking for the gas outflow
	  ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]},uconBL[4];    
	  conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
	  trans2_coco(geom.xxvec,ucon,uconBL,MYCOORDS,BLCOORDS);
	
	  if(uconBL[1]>0.) //outflow
	    {
	      uconBL[1]=0.;
	      trans2_coco(geomBL.xxvec,uconBL,ucon,BLCOORDS,MYCOORDS);
	      conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
	      pp[VX]=ucon[1];
	      pp[VY]=ucon[2];
	      pp[VZ]=ucon[3];
	    }
	}

     p2u(pp,uu,&geom);
     return 0;
   }

//reflections in theta 
if(BCtype==YBCLO) //axis
  {      
    #ifdef SPECIAL_BC_CHECK
    int ret_val = 0;
    ret_val = if_indisk_bc_check(ix+TOI,iy+TOJ,iz+TOK);
    if(ret_val != 0) iiy = iy + ret_val;
    else iiy = -iy-1; 

    /*
    if(ix == 63 && iy == 41)
    { 
      printf("In bc check for YBCLO.\n");
      print_Nvector(pp,NV);
      printf("\n");
    }
    */
    #else
    iiy=-iy-1;
    #endif
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
    for(iv=0;iv<NV;iv++)
      {
	//theta component
        #ifdef MAGNFIELD
	if(iv==VY || iv==B2 || iv==FY0)
        #else
	if(iv==VY || iv==FY0)
        #endif
	  pp[iv]=-get_u(p,iv,iix,iiy,iiz);
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
       }

    p2u(pp,uu,&geom);
    return 0;
  }
if(BCtype==YBCHI)//axis or eq.plane
  {
    #ifdef SPECIAL_BC_CHECK
    int ret_val = 0;
    ret_val = if_indisk_bc_check(ix+TOI,iy+TOJ,iz+TOK);
    if(ret_val != 0) iiy = iy + ret_val;
    else iiy=NY-(iy-NY)-1; 

    /*
    if(ix == 63 && iy == 22)
    { 
      printf("In bc check for YBCHI.\n");
      print_Nvector(pp,NV);
      printf("\n");
    }
    */
    #else
    iiy=NY-(iy-NY)-1;
    #endif
    iiz=iz;
    iix=ix;
    gdet_src=get_g(g,3,4,iix,iiy,iiz);  
    gdet_bc=get_g(g,3,4,ix,iy,iz);  
  	  
    for(iv=0;iv<NV;iv++)
      {
        #ifdef MAGNFIELD
	if(iv==VY || iv==B2 || iv==FY0)
        #else
	if(iv==VY || iv==FY0)
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
if(iz<0) iiz=iz+NZ;
if(iz>NZ-1) iiz=iz-NZ;

for(iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u,iv,iix,iiy,iiz);
    pp[iv]=get_u(p,iv,iix,iiy,iiz);      
  }

return 0;

