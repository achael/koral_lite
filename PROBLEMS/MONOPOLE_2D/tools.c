int
donut_analytical_solution(ldouble *pp,ldouble *xxvecBL,ldouble ggBL[][5],ldouble GGBL[][5] )
{
  ldouble xx=xxvecBL[1];
 
  ldouble podpierd=-(GGBL[0][0]-2.*ELL*GGBL[0][3]+ELL*ELL*GGBL[3][3]);
  ldouble ut=-1./sqrt(podpierd);

  ut/=UTPOT; //rescales rin
  ldouble Vphi,Vr;
  ldouble D,W,eps,uT,uphi,uPhi,rho,ucon[4],uint,E,Fx,Fy,Fz;
  if(ut<-1 || podpierd<0. || xx<3. || NODONUT || INFLOWING)
    return -1; //outside donut

  ldouble h=-1./ut;
  eps=(h-1.)/GAMMA;
  rho=pow(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.));
  uint=rho*eps;
  uphi=-ELL*ut;
  uT=GGBL[0][0]*ut+GGBL[0][3]*uphi;
  uPhi=GGBL[3][3]*uphi+GGBL[0][3]*ut;
  Vphi=uPhi/uT;
  Vr=0.;

  //4-velocity in BL
  pp[2]=pp[3]=pp[4]=0.;
  pp[0]=rho;
  pp[1]=uint;

#ifdef RADIATION
  ldouble P,aaa,bbb;
  P=GAMMAM1*uint;
  //solving for T satisfying P=pgas+prad=bbb T + aaa T^4
  aaa=4.*SIGMA_RAD/3.;
  bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
  ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));
  ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;

  E=calc_LTE_EfromT(T4);
  Fx=Fy=Fz=0.;
  uint=calc_PEQ_ufromTrho(T4,rho);

  pp[1]=uint;
  pp[6]=E;
  pp[7]=Fx;
  pp[8]=Fy;
  pp[9]=Fz;
#endif

  return 0;
}
