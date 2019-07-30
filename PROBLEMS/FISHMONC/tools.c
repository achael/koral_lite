// Fishbone & Moncrief (ApJ, 207, 962, 1976) torus with inner edge at rin and pressure maximum at rmax.

int init_dsandvels_fishbone_moncrief(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ellout)
{
  // FM_rin, FM_rmax are the equatorial radii of the torus inner endge and pressure maximum
  
  ldouble rin = FM_rin;
  ldouble rmax = FM_rmax;
  
  ldouble Delta = r*r - 2.*r + a*a;
  ldouble Sigma = r*r + a*a*cos(th)*cos(th);
  ldouble A = (r*r + a*a)*(r*r + a*a) - Delta*a*a*sin(th)*sin(th);

  // The constant l is evalulated at rmax using eq (3.8) in FM76
  
  ldouble ell = (rmax*rmax*rmax*rmax + rmax*rmax*a*a - 2.*rmax*a*a - a*sqrt(rmax)*(rmax*rmax - a*a)) / (rmax*rmax - 3.*rmax + 2.*a*sqrt(rmax)) / pow(rmax,1.5);
  
  //ln h(r,th) (eq (3.6) in FM76)
  
  ldouble hlog = 0.5*log((1. + sqrt(1. + 4.*ell*ell*Sigma*Sigma*Delta/(A*sin(th)*A*sin(th))))/(Sigma*Delta/A)) - 0.5*sqrt(1 + 4.*ell*ell*Sigma*Sigma*Delta/(A*sin(th)*A*sin(th))) - 2.*a*r*ell/A - 0.5*log((1. + sqrt(1. + 4.*ell*ell*rin*rin*(rin*rin - 2.*rin + a*a)*pow((rin*rin*rin + rin*a*a + 2.*a*a),-2)))/(rin*(rin*rin - 2.*rin + a*a)/(rin*rin*rin + rin*a*a + 2.*a*a))) + 0.5*sqrt(1. + 4.*ell*ell*rin*rin*(rin*rin - 2.*rin + a*a)/pow((rin*rin*rin + rin*a*a + 2.*a*a),2)) + ell*2.*a/(rin*rin*rin + rin*a*a + 2.*a*a);
  ldouble h = exp(hlog);
                         
  ldouble Delta0 = rmax*rmax - 2.*rmax + a*a;
  ldouble Sigma0 = rmax*rmax;
  ldouble A0 = (rmax*rmax + a*a)*(rmax*rmax + a*a) - Delta0*a*a;

  //ln h(rmax,pi/2) (eq (3.6) in FM76)
  
  ldouble hlog0 = 0.5*log((1. + sqrt(1. + 4.*ell*ell*Sigma0*Sigma0*Delta0/(A0*A0)))/(Sigma0*Delta0/A0)) - 0.5*sqrt(1 + 4.*ell*ell*Sigma0*Sigma0*Delta0/(A0*A0)) - 2.*a*rmax*ell/A0 - 0.5*log((1. + sqrt(1. + 4.*ell*ell*rin*rin*(rin*rin - 2.*rin + a*a)*pow((rin*rin*rin + rin*a*a + 2.*a*a),-2)))/(rin*(rin*rin - 2.*rin + a*a)/(rin*rin*rin + rin*a*a + 2.*a*a))) + 0.5*sqrt(1. + 4.*ell*ell*rin*rin*(rin*rin - 2.*rin + a*a)/pow((rin*rin*rin + rin*a*a + 2.*a*a),2)) + ell*2.*a/(rin*rin*rin + rin*a*a + 2.*a*a);
  ldouble h0 = exp(hlog0);
  
  ldouble n = 1. / (GAMMA - 1.);
  ldouble rho, uu, p0, p;
  if (h >= 1. && r >= rin)  // the point is inside the torus
  {
    // Set rho = FM_rho0 ((h-1)/(h0-1))^n  (note: FM_rho0 is the maximum density)
    //      p0 = FM_rho0 (h0-1)/(n+1)  (note: p0 is the maximum pressure)
    //       p = p0 ((h-1)/(h0-1))^(n+1)
    //       u = p/(Gamma-1)
    
    rho = FM_rho0 * pow(((h - 1.) / (h0 - 1.)),n);
    p0 = FM_rho0 * (h0 - 1.) / (n + 1.);
    p = p0 * pow(((h - 1.) / (h0 - 1.)),(n+1));
    uu = p / (GAMMA - 1.);
  }
  else  // the point is outside the torus
  {
    // Set rho = u = -1 to indicate to the calling program that this is outside the torus
    
    rho = -1.;
    uu = -1.;
  }
  
  // Copy and return results
  
  *rhoout = rho;
  *uuout = uu;
  *ellout = ell;
  
  //printf("%e %e %e %e %e %e\n", r, h-1., h0-1., rho, p, ell);
  
  return 0;
}

