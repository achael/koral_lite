
int init_dsandvels_katotorus(FTYPE r, FTYPE th, FTYPE phic, FTYPE drot, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell)
{
  ldouble rho,pre,phi,phi0,uint,phit0,phit,L,L0,cs0,uint0,rho0,n,R,xd,yd,zd,pf;
  
  if(r<KT_RMIN) {*rhoout=-1.;return 0;}

  xd = r*sin(th)*cos(phic);
  yd = r*sin(th)*sin(phic);
  zd = r*cos(th);

  // rotation of disk
  pf = yd*cos(drot) - zd*sin(drot); //prefactor

  R=sqrt( xd*xd + pf*pf ); //cylindrical radius 

  L0=sqrt(KT_R0*KT_R0*KT_R0)/(KT_R0-2.);
  L=L0*pow(R/KT_R0,KT_A);
  uint0=calc_PEQ_ufromTrho(KT_T0,KT_RHO0,0,0,0);
  rho0=KT_RHO0;
  cs0=sqrt(GAMMA*(GAMMAM1*uint0)/rho0);
  phi = -1./(r-2.);
  phi0 = -1./(KT_R0-2.);
  phit = phi + 1/2./(1.-KT_A)*(L/R)*(L/R);
  phit0 = phi0 + 1/2./(1.-KT_A)*(L0/KT_R0)*(L0/KT_R0);
  n = 3.; //but GAMMA = 5./3. !
  
  rho = rho0  * pow(1. - GAMMA/cs0/cs0*(phit-phit0)/(n+1.),n);
  pre = rho0 * cs0 * cs0 / GAMMA * pow(rho/rho0,1.+1./n);
  uint = pre / GAMMAM1;

  #ifdef KT_RHO_FLOOR //This keeps the same Bernoulli number but reduces mass in disk
  if (rho < KT_RHO_FLOOR) {*rhoout=-1.;return 0;}
  #endif

  *rhoout = rho;
  *uuout = uint;
  *ell = L;

  return 0;

}

