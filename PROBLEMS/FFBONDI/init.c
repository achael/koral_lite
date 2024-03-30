ldouble rho,mx,my,mz,m,E,uint,ur;
ldouble uu[NV],pp[NV],ppback[NV],T;

struct tfun_params
{
  ldouble n,r,C4,C5,T_crit,Kadiab;
};

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BONDICOORDS);

// init
init_bondi_full(pp, ix, iy, iz);

/*
//overwrite magnetic field
#if defined(MANGFIELD) && defined(VECPOTGIVEN)
#ifdef INIT_MAGN_CORNERS    
// ANDREW define vector potential at corners not cell centers
// TODO: Issues going to grid edge! We assume we are well inside grid!
ldouble xxvec_c[4],xxvecBL_c[4];    
xxvec_c[0] = global_time;
xxvec_c[1] = get_xb(ix,0);
xxvec_c[2] = get_xb(iy,1);
xxvec_c[3] = get_xb(iz,2);
coco_N(xxvec_c,xxvecBL_c,MYCOORDS,BONDICOORDS);

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
Acov[3]=bfac*(1.-cos(th_mag));

// Vector potential A is temporarily saved in pp[B1], pp[B2], pp[B3].
// These will be replaced by B components immediately
pp[B1]=Acov[1];
pp[B2]=Acov[2];
pp[B3]=Acov[3];
#endif //MAGNFIELD and VECPOTGIVEN
*/

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

#ifdef FORCEFREE
set_u_scalar(ffinvarr, ix, iy, iz, 1.0);
set_cflag(FFINVFLAG, ix,iy,iz,1);
set_cflag(MHDINVFLAG,ix,iy,iz,0);
#endif
