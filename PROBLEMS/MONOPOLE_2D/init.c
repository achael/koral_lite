ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

//geometries
get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
xx=xxvecBL[1];
yy=xxvecBL[2];
zz=xxvecBL[3];

ldouble gg[4][5],GG[4][5],eup[4][4],elo[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,KERRCOORDS);
calc_G_arb(xxvecBL,GGBL,KERRCOORDS);

//ldouble eupBL[4][4],eloBL[4][4];
//ldouble tupBL[4][4],tloBL[4][4];
//calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
//calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);

ldouble r=geomBL.xx;
ldouble th=geomBL.yy;
ldouble rhor = (1. + sqrt(1. - BHSPIN*BHSPIN)) ;

rho = RHOATMMIN + (r/10./rhor)/pow(r,4)/B2RHORATIOMAX;
uint = UINTATMMIN + (r/10./rhor)/pow(r,4)/B2RHORATIOMAX;

//4-velocity in BL transformed to MYCOORDS
//ldouble ucon[4]={0.,0.,0.,0.};
//ucon[0] = -1/sqrt(-GGBL[0][0])
//ucon[3] = GG[0][3]*ucon[0]
//conv_vels(ucon,ucon,VEL3,VELPRIM,ggBL,GGBL);
   
pp[0]=rho;
pp[1]=uint;
pp[2]=0.; //zero in co-rotating frame
pp[3]=0.;
pp[4]=0.;

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.;
    pp[B1]=1./(r*r);
#endif

    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,xxvecBL,&geomBL,&geom);
    
#ifdef MAGNFIELD //MYCOORDS vector potential to calculate B's
    ldouble Acov[4];
    Acov[0]=Acov[1]=Acov[2]=0.;
    Acov[3]=(1.-cos(th));///(r*sin(th));

//pp[B1]=0.;
//pp[B2]=0.;
//pp[B3]=Acov[3];
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
