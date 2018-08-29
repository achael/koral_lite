// Initial torus based on Fishbone & Moncrief (ApJ, 207, 962, 1976)


/////////////////////////////////////////////////////////////////////
int init_dsandvels_fishbone_moncrief(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell);

/////////////////////////////////////////////////////////////////////

ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T,uintorg;
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

ldouble r=geomBL.xx;
ldouble th=geomBL.yy;

// call function init_dsandvels_fishbone_moncrief, which computes density rho, internal energy uint and angular momentum ell at current location (r,th)
init_dsandvels_fishbone_moncrief(r, th, BHSPIN, &rho, &uint, &ell); 
uintorg=uint;

if(rho<0.) //outside donut
{
  //atmosphere: atmtype=0. Need to supply RHOATMMIN, UINTATMMIN at rout=2.
  set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
  set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif
}
else //inside donut
{
  //calculate atmosphere values as backup
  set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
  set_radatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#endif
  
  pgas = GAMMAM1 * uint;
 
  
  ldouble ult,ulph,ucov[4],ucon[4];
  
  // This is the original version
  //ell*=-1.;
  //ulph = sqrt(-1./(geomBL.GG[0][0]/ell/ell + 2./ell*geomBL.GG[0][3] + geomBL.GG[3][3]));
  //ult = ulph / ell;

  //ulph = ell;
  //ult = (-geomBL.GG[0][3]*ell - sqrt((geomBL.GG[0][3]*geomBL.GG[0][3]*ell*ell - geomBL.GG[0][0]*(1+geomBL.GG[3][3]*ell*ell))))/geomBL.GG[0][0];
  
  //ucov[0]=ult;
  //ucov[1]=0.;
  //ucov[2]=0.;
  //ucov[3]=ulph;
  
  //indices_12(ucov,ucon,geomBL.GG);
  
  // The following is based on the statement in Kozlowski et al. that ell is actually u_phi u^t
  ldouble rhosq = geomBL.gg[0][0] * geomBL.gg[3][3] - geomBL.gg[0][3] * geomBL.gg[0][3];
  ldouble utsq = (-geomBL.gg[3][3] - sqrt(geomBL.gg[3][3] * geomBL.gg[3][3] - 4. * rhosq * ell * ell)) / (2. * rhosq);
    
  ucon[0] = sqrt(utsq);
  ucon[1] = 0.;
  ucon[2] = 0.;
  ucon[3] = ell / (geomBL.gg[3][3] * ucon[0]) - geomBL.gg[0][3] * ucon[0] / geomBL.gg[3][3];
  
  indices_12(ucon,ucov,geomBL.gg);
  
  conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
  
  // set pp[0], pp[1], but use the atmosphere values as floors
  pp[0] = rho;
  pp[1] = uint;
  if (pp[0] < ppback[0] || pp[1] < ppback[1])
  {
    pp[0] = ppback[0];
    pp[1] = ppback[1];
  }
  pp[2]=ucon[1];
  pp[3]=ucon[2];
  pp[4]=ucon[3];
  
#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
  pp[B1]=pp[B2]=pp[B3]=0.;
#endif
  
#ifdef RADIATION
  //distributing pressure
  ldouble P,aaa,bbb;
  P=GAMMAM1*uint;
  //solving for T satisfying P=pgas+prad=bbb T + aaa T^4
  aaa=4.*SIGMA_RAD/3.;
  bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
  ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));
  ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;
  
  E=calc_LTE_EfromT(T4);
  Fx=Fy=Fz=0.;
  uint=calc_PEQ_ufromTrho(T4,rho,ix,iy,iz);
  
  pp[UU]=my_max(uint,ppback[1]);
  pp[EE0]=my_max(E,ppback[EE0]);
  
  
  
  pp[FX0]=Fx;
  pp[FY0]=Fy;
  pp[FZ0]=Fz;
  
  //transforming from BL lab radiative primitives to code non-ortonormal primitives
  prad_ff2lab(pp,pp,&geomBL);
  
#endif
  
#ifdef EVOLVEPHOTONNUMBER
  pp[NF0]=calc_NFfromE(pp[EE0]);
#endif
  
  //transforming primitives from BL to MYCOORDS
  trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
  
#ifdef MAGNFIELD
  //MYCOORDS vector potential to calculate B's
  ldouble Acov[4];
  ldouble r_mag,th_mag;
  Acov[0]=Acov[1]=Acov[2]=0.;

#ifdef INIT_MAGN_CORNERS    
    //ANDREW define vecpot at corners not cell centers
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
    //define vecpot at cell centers
    r_mag=r;
    th_mag=th;
#endif
      
#ifndef MULTILOOPS
  //standard single poloidal loop: Aphi = max[(rho/rhomax)-FM_Aphi_cut, 0]
  //recompute rho -- possibly at corner
  ldouble u_av_chop, u_av_mid;
  init_dsandvels_fishbone_moncrief(r_mag, th_mag, BHSPIN, &rho, &uint, &ell);
  if(fabs((rho-pp[RHO])/pp[RHO])>.5) rho=pp[RHO];
  //vector potential
  Acov[3]=my_max((rho/FM_rho0) - FM_Aphi_cut, 0.);
#else
  ldouble lambda = 0.5;
  ldouble anorm = 1.;
  ldouble rchop = 90.;
  ldouble u_av = uintorg;
  ldouble u_av_chop, u_av_mid;
  //midplane at r_mag
  init_dsandvels_fishbone_moncrief(r_mag, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
  //midplane at rchop
  init_dsandvels_fishbone_moncrief(rchop, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
  //vector potential follows contours of UU
  ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
  ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane
  
  //ldouble rin=KT_R0/2.;
  ldouble rin=FM_rin;
  ldouble STARTFIELD = 1.25*rin;
  ldouble q, fr, fr_start, vpot=0.;
  if (r_mag > STARTFIELD && r_mag < rchop) {
    q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th_mag), 3); // * pow(tanh(r_mag/rsmooth),2);
  } else q = 0;
    
    
    if(q > 0.) {
      fr = (pow(r_mag,0.6)/0.6  + 0.5/0.4*pow(r_mag,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }
  
  //    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
  Acov[3]=vpot*sin((M_PI/2.-th_mag));;
#endif
  
  // Vector potential A is temporarily saved in pp[B1], pp[B2], pp[B3]. These will be replaced by B components immediately
  pp[B1]=Acov[1];
  pp[B2]=Acov[2];
  pp[B3]=Acov[3];
#endif
  
}

//entropy
pp[5]=calc_Sfromu(pp[0],pp[1],ix,iy,iz);
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
update_entropy(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);
