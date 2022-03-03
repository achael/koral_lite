ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble uu[NV],pp[NV],ppback[NV],T;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

ldouble r = geomBL.xx;
ldouble th = geomBL.yy;
ldouble rhor = (1. + sqrt(1. - BHSPIN*BHSPIN)) ;

rho = RHOATMMIN + (r/10./rhor)/pow(r,4)/B2RHORATIOMAXINIT;
uint = UINTATMMIN + (r/10./rhor)/pow(r,4)/B2UURATIOMAXINIT;
   
pp[0]=rho;
pp[1]=uint;
pp[2]=0.; //zero initial velocity in co-rotating frame
pp[3]=0.;
pp[4]=0.;

//entropy
pp[5]=calc_Sfromu(pp[0],pp[1],geom.ix,geom.iy,geom.iz);

#ifdef MAGNFIELD // set to zero at first
pp[B1]=pp[B2]=pp[B3]=0.;
#endif

//transform primitives from BL to MYCOORDS
//trans_pall_coco(pp, pp, BLCOORDS, MYCOORDS,xxvecBL,&geomBL,&geom);
    
#ifdef MAGNFIELD //calculate vector potential
ldouble r_mag,th_mag;

#ifdef INIT_MAGN_CORNERS    
// ANDREW define vector potential at corners not cell centers
// TODO: Issues going to grid edge! We assume we are well inside grid!
ldouble xxvec_c[4],xxvecBL_c[4];    
xxvec_c[0] = global_time;
xxvec_c[1] = get_xb(ix,0);
xxvec_c[2] = get_xb(iy,1);
xxvec_c[3] = get_xb(iz,2);
coco_N(xxvec_c,xxvecBL_c,MYCOORDS,BLCOORDS);

ldouble r_c=xxvecBL_c[1];
ldouble th_c=xxvecBL_c[2];

r_mag=r_c;
th_mag=th_c;
#else
r_mag=r;
th_mag=th;
#endif

ldouble Acov[4];
Acov[0]=Acov[1]=Acov[2]=0.;
Acov[3]=(1.-cos(th_mag));

// Vector potential A is temporarily saved in pp[B1], pp[B2], pp[B3].
// These will be replaced by B components immediately
pp[B1]=Acov[1];
pp[B2]=Acov[2];
pp[B3]=Acov[3];
#endif

//all primitives to conserved
p2u(pp,uu,&geom);

/***********************************************/
// set the final values
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
