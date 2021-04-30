//definitions
ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T;
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

//coordinates
ldouble r=geomBL.xx;
ldouble th=geomBL.yy;

//printf("ix iy iz = %i %i %i\n",ix,iy,iz);
//printf("g : rr thth phiphi : %f %f %f \n",geomBL.gg[1][1],geomBL.gg[2][2],geomBL.gg[3][3]);
//printf("Kr : rr thth phiphi : %f %f %f \n",get_gKr(1,1,1,ix,iy,iz),get_gKr(2,2,2,ix,iy,iz),get_gKr(3,3,3,ix,iy,iz));

//ambient
set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);

#ifdef EVOLVEPHOTONNUMBER
    pp[NF0]=calc_NFfromE(pp[EE0]);
    //pp[NF]=calc_NFfromT(calc_PEQ_Tfromurho(pp[UU],pp[RHO]));
    //pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
    //pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/ATMTRADINIT;
#endif


#endif

#ifdef MAGNFIELD //setting to zero everywhere. Does this prevent initial failure?
pp[B1]=pp[B2]=pp[B3]=0.;
#endif

//entropy
pp[5]=calc_Sfromu(pp[0],pp[1],ix,iy,iz);
//to conserved
p2u(pp,uu,&geom);

//if((ix+TOI) == (STREAM_IX-1) && (iy+TOJ) >= STREAM_IYT && (iy+TOJ) <= STREAM_IYB)
//{
//  printf("%e \n",pp[0]);
//}

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
