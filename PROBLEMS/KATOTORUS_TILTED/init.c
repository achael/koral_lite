int init_dsandvels_katotorus(FTYPE r, FTYPE th, FTYPE phic, FTYPE drot, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell);

ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T,uintorg;
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;
ldouble xd,yd,zd,xp,yp,zp,lp,sinphi,cosphi,sinphip,cosphip,sina,cosa,sinth,costh,sinthp,costhp; //primed angles and coords
ldouble dthdt,dphidt,omega,vmag,ellphi,ellth,pref1,pref2;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

ldouble r=geomBL.xx;
ldouble th=geomBL.yy;
ldouble phi=geomBL.zz;

init_dsandvels_katotorus(r, th, phi, KT_DROT, BHSPIN, &rho, &uint, &ell); 
uintorg=uint;

/*
if( (ix+TOI) == (TNX-1) && (iy+TOJ) == (TNY-1) && (iz+TOK) == (TNZ-1)) 
{
  ldouble r1,r2,sinth1,sinth2; //Needed for spherical curl
  r1 = exp(get_xb(ix,0));
  r2 = exp(get_xb(ix+1,0));
  sinth1 = sin(M_PI*get_xb(iy,1));
  sinth2 = sin(M_PI*get_xb(iy+1,1));

  printf("r r1 r2 : %e %e %e\n",r,r1,r2);
  printf("th th1 th2 : %e %e %e \n",th,M_PI*get_xb(iy,1),M_PI*get_xb(iy+1,1));
  printf("sinth sinth1 sinth2 : %e %e %e \n",sin(th),sinth1,sinth2);
}
*/

if(rho<0.) //outside donut
  {
    //ambient
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif
  }
 else //inside donut
   {
    //ambient
    set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#endif

    ldouble xxvec[4],dxdx[4][4];
    xxvec[0]=0.;xxvec[1]=get_xb(ix,0);xxvec[2]=get_xb(iy,1);xxvec[3]=get_xb(iz,2);
    dxdx_MKS22KS(xxvec,dxdx); //dx^mu/dx'^mu

    pgas = GAMMAM1 * uint;
    //ell*=-1.;

    //cartesian based on spherical coords
    xd = r*sin(th)*cos(phi);
    yd = r*sin(th)*sin(phi);
    zd = r*cos(th);

    // x and y after rotating by drot
    xp = xd;
    yp = yd*cos(KT_DROT) - zd*sin(KT_DROT);
    zp = yd*sin(KT_DROT) + zd*cos(KT_DROT);
    lp = sqrt(xp*xp + yp*yp);

    //trig quantities necessary for defining angular momentum in unrotated (BH centered coords)
    sinphi = sin(phi);
    cosphi = cos(phi);
    sinphip = yp/lp;
    cosphip = xp/lp;
    sina = sin(KT_DROT);
    cosa = cos(KT_DROT);
    sinth = sin(th);
    costh = cos(th);
    sinthp = lp/r;
    costhp = zp/r;

    omega = ell/(r*r*sinthp*sinthp);
    if (fabs(KT_DROT) > 0)
    {
      pref1 = cosa - sinphi*sina*costh/sinth;
      pref2 = cosphi*sina;
    }
    else
    {
      pref1 = 1.;
      pref2 = 0.;
    }

    dphidt = omega*pref1; // d(phi)/dt
    dthdt = omega*pref2; // d(theta)/dt
    //if(xd < 0)
    //{
    //  dthdt *= -1.;
    //}

    //vmag = omega*sinthp;
    vmag = sqrt(sinth*sinth*dphidt*dphidt + dthdt*dthdt);

    //angular momentum components
    //pref1 = sinphi*cosphip*cosa - cosphi*sinphip;
    //pref2 = cosphip*sina;
    //ellphi = ell*(cosphi*cosphip*cosa + sinphi*sinphip);
    //ellth = ell*sqrt( pref1*pref1 + pref2*pref2 );//*my_sign(KT_DROT);

    //if(xd < 0)
    //{
    //  ellth *= -1.; //dth/dt < 0 in x < 0 region
    //}

    ldouble ult,ult2,ulph,ulth,ucov[4],ucon[4];

    //ulph = sqrt(-1./(geomBL.GG[0][0]/ell/ell + 2./ell*geomBL.GG[0][3] + geomBL.GG[3][3]));
    //ulph = sqrt(-1./(geomBL.GG[0][0]/ellphi/ellphi + 2./ellphi*geomBL.GG[0][3] + geomBL.GG[3][3]));
    //ulph = sqrt(-1./(geomBL.GG[0][0]/ellphi/ellphi + 2./ellphi*geomBL.GG[0][3] + geomBL.GG[3][3] + 2.*ellth/ellphi/ellphi*geomBL.GG[0][2] + ellth*ellth/ellphi/ellphi*geomBL.GG[2][2]));

    //ulth = 0.;
    //if (KT_DROT > 0)
    //{
      //ulth = sqrt(-1./(geomBL.GG[0][0]/ellth/ellth + 2./ellth*geomBL.GG[0][2] + geomBL.GG[2][2]));
      //ulth = sqrt(-1./(geomBL.GG[0][0]/ellth/ellth + 2./ellth*geomBL.GG[0][2] + geomBL.GG[2][2] + 2.*ellphi/ellth/ellth*geomBL.GG[0][3] + ellphi*ellphi/ellth/ellth*geomBL.GG[3][3]))*my_sign(KT_DROT);
      //if(xd < 0)
      //{
      //  ulth *= -1.;
      //}
    //}
    ult = sqrt(-1./(geomBL.gg[0][0] + dthdt*dthdt*geomBL.gg[2][2] + 2.*dphidt*geomBL.gg[0][3] + dphidt*dphidt*geomBL.gg[3][3]));
    ulph = ult*dphidt;
    ulth = ult*dthdt;

    //printf("%f %f",ult,ult2);

    ucon[0]=ult;
    ucon[1]=0.;
    ucon[2]=ulth;
    ucon[3]=ulph;
    
    //indices_12(ucov,ucon,geomBL.GG);
    indices_21(ucon,ucov,geomBL.gg);

    conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
   
  
    pp[0]=my_max(rho,ppback[0]); 
    pp[1]=my_max(uint,ppback[1]);
    pp[2]=ucon[1]; 
    pp[3]=ucon[2]; //derived in BL?
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

#ifdef EVOLVEPHOTONNUMBER
    pp[NF0]=calc_NFfromE(pp[EE0]);
#endif

    //transforming from BL lab radiative primitives to code non-ortonormal primitives
    prad_ff2lab(pp,pp,&geomBL);

#endif

    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
    
#ifdef MAGNFIELD 
    //MYCOORDS vector potential to calculate B's
    ldouble Acov[4],A0;
    ldouble dthdthp,dphidthp;//Spherical coordinate transforms for rotated to unrotated derivative
    ldouble aaa1,aaa2,aaa3; //for ease of writing transforms
    Acov[0]=0.;
    Acov[1]=0.;
    Acov[2]=0.;

    aaa1 = (sinthp*cosa + costhp*sinphip*sina);
    aaa2 = pow((costhp*cosa - sinthp*sinphip*sina),2.);
    aaa3 = pow( (1. - aaa2),0.5);
    dthdthp = aaa1/aaa3;
    
    aaa1 = -sina/cosphip/sinthp/sinthp;
    aaa2 = 1. + pow( ( sinphip*cosa/cosphip + sina*costhp/sinthp/cosphip), 2.);
    dphidthp = aaa1/aaa2;

    ldouble Bpref1,Bpref2,Bpref3;
    Bpref1 = sqrt(cosphi*cosphi + (sinphi*cosa - sina*costh/sinth));
    Bpref2 = sina*cosphi/sinth;
    Bpref3 = cosa - sina*sinphi*costh/sinth;  
 
#if(NTORUS==1)
    //standard single poloidal loop
    A0=my_max(rho-1.e-1*KT_RHO0,0.);

    //sending MKS2 (converting from BL) - dxdx=0 and doesn't work in MKS3...
    Acov[2] = -sinth*(1.+cosa*cosphi*cosphi)*dphidthp*A0;
    Acov[3] = sinth*dthdthp*A0;
#elif (NTORUS==2)
    //standard single poloidal loop
    ldouble rin=KT_R0;
    ldouble STARTFIELD = 100.*rin;
    ldouble rout = 5.e3; //cover full range of torus?
    ldouble lambda = 2.*(rout-rin);

    A0 = my_max( pow(r, RADIUS_POWER) * (pp[RHO] - RHO_CUT_FACTOR * KT_RHO0) * pow(sinthp, SIN_POWER_THETA) * pow(sin(r/lambda), SIN_POWER) , 0.);
    Acov[2] = -sinth*(1.+cosa*cosphi*cosphi)*dphidthp*A0;
    Acov[3] = sinth*dthdthp*A0;
#elif (NTORUS==3) // dipolar loop
    ldouble lambda = 35.;
    ldouble anorm = 1.;
    ldouble rchop = 3.6e3;
    ldouble u_av = uintorg;
    ldouble u_av_chop, u_av_mid;
    //midplane at r
    init_dsandvels_katotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_katotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    //vetor potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=KT_R0/2.;
    ldouble STARTFIELD = 1.25*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r > STARTFIELD && r < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sinthp, SIN_POWER_THETA); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    
    if(q > 0.) {
      fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }
     
    //    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
    A0=vpot;//*sin((M_PI/2.-geomBL.yy));;
    Acov[2] = -sinth*(1.+cosa*cosphi*cosphi)*dphidthp*A0;
    Acov[3] = sinth*dthdthp*A0;
#elif (NTORUS==4) //Quadrupolar loops
    //standard single poloidal loop
    A0=costhp*my_max(rho-1.e-1*KT_RHO0,0.);

    //sending MKS2 (converting from BL) - dxdx=0 and doesn't work in MKS3...
    Acov[2] = -sinth*(1.+cosa*cosphi*cosphi)*dphidthp*A0;
    Acov[3] = sinth*dthdthp*A0;
#endif

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
