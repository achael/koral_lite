int init_dsandvels_katotorus(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell);

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

init_dsandvels_katotorus(r, th, BHSPIN, &rho, &uint, &ell); 
uintorg=uint;

if(rho<0.) //outside donut
  {
    //ambient
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif

#ifdef MAGNFIELD
#ifdef BATMZ //put field in the atmosphere and torus (Better for MAD with small torus?)
    if(r*sin(th) < BATMZ_MAXR && r*sin(th) > BATMZ_MINR)
    {
      pp[B1]=0.;//BATMZ_B0*cos(th)/r/r; //will this get auto rescaled?
      pp[B2]=0.;//BATMZ_B0*sin(th)/r/r;  
      pp[B3]=0.;
    }
#endif
#endif
  }
 else //inside donut
   {
    //ambient
    set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#endif

    pgas = GAMMAM1 * uint;
    ell*=-1.;

    ldouble ult,ulph,ucov[4],ucon[4];
    ulph = sqrt(-1./(geomBL.GG[0][0]/ell/ell + 2./ell*geomBL.GG[0][3] + geomBL.GG[3][3]));
    ult = ulph / ell;

    ucov[0]=ult;
    ucov[1]=0.;
    ucov[2]=0.;
    ucov[3]=ulph;
    
    indices_12(ucov,ucon,geomBL.GG);

    conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
   
  
    pp[0]=my_max(rho,ppback[0]); 
    pp[1]=my_max(uint,ppback[1]);
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
    ldouble Acov[4];
    Acov[0]=Acov[1]=Acov[2]=0.;

#if(NTORUS==1)
    //standard single poloidal loop
    Acov[3]=my_max(pp[RHO]-1.e-1*KT_RHO0,0.);
#elif (NTORUS==2)
    //standard single poloidal loop
    ldouble rin=KT_R0;
    ldouble STARTFIELD = KT_RIN;//0.5*rin;//100.*rin;
    ldouble rout = KT_ROUT;//5.e3; //cover full range of torus?
    ldouble lambda = 2.*(rout-rin);

//    if(r > STARTFIELD)
//    {
      Acov[3] = my_max( pow(r, RADIUS_POWER) * (pp[RHO] - RHO_CUT_FACTOR * KT_RHO0) * pow(sin(th), 3) * pow(sin(r/lambda), SIN_POWER) , 0.);
//    }
//    else
//    {
//      Acov[3] = my_max( (pp[RHO] - RHO_CUT_FACTOR * KT_RHO0) * pow(sin(th), 3), 0.);
//    }
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
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    
    if(q > 0.) {
      fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }
     
    //    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
    Acov[3]=vpot;//*sin((M_PI/2.-geomBL.yy));;

#elif (NTORUS==4) //Quadrupolar loops
    ldouble k = -0.4;
    ldouble lambda = 0.4;
    ldouble anorm = 1.;
    ldouble rchop = 4.e3;
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
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    
    if(q > 0.) {
      fr = (r  + 0.6*0.5/0.4) * pow(r,k) / lambda /0.6;
      fr_start = (STARTFIELD  + 0.6*0.5/0.4) * pow(STARTFIELD,k) / lambda / 0.6;
      vpot += q * sin(fr - fr_start) ;
    }
     
    //    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
    Acov[3]=vpot*sin((M_PI/2.-geomBL.yy));;
#elif (NTORUS==5) //Modified Quadrupolar loops
    ldouble k = -0.4;
    ldouble lambda = 0.4;
    ldouble anorm = 1.;
    ldouble rchop = 3.e2;
    ldouble rbreak = 1.e2;
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
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3.); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    
    if(q > 0.) {
      fr = r * pow(r*pow(1. + r/rbreak,0.5) ,k)/lambda;
      fr_start = STARTFIELD * pow(STARTFIELD*pow(1. + STARTFIELD/rbreak,0.5) ,k)/lambda;
      vpot += q * sin(fr - fr_start) ;
    }
     
    //    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
    Acov[3]=vpot*sin((M_PI/2.-geomBL.yy));;
#endif

    pp[B1]=Acov[1];
    pp[B2]=Acov[2];
    pp[B3]=Acov[3];

#ifdef BATMZ //put field in the atmosphere and torus (Better for MAD with small torus?)
    if(r*sin(th) < BATMZ_MAXR && r*sin(th) > BATMZ_MINR)
    {
      pp[B1]=BATMZ_TORUSRESCALE*BATMZ_B0*cos(th)/r/r; //will this get auto rescaled?
      pp[B2]=BATMZ_TORUSRESCALE*BATMZ_B0*sin(th)/r/r;  
      pp[B3]=0.;
    }
#endif //overwrites any previous B field initialization

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
