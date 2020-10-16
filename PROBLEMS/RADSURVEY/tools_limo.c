void compute_gd( FTYPE r, FTYPE th, FTYPE a, FTYPE *gdtt, FTYPE *gdtp, FTYPE *gdpp ) {
   FTYPE Sigma, tmp;

   Sigma = (pow(a*cos(th),2) + pow(r,2));
   tmp = 2*r*pow(sin(th),2)/Sigma;

   //metric (expressions taken from limotorus4.nb):
   *gdtt = -1 + 2*r/Sigma;
   *gdtp = -a*tmp;
   *gdpp = (pow(r,2) + pow(a,2)*(1+tmp)) * pow(sin(th),2);
}

// Keplerian, equatorial angular momentum density: l = u_phi / u_t
FTYPE lK(FTYPE r, FTYPE a) {
   FTYPE curlyF, curlyG;

   curlyF = 1 - 2*a/pow(r, 1.5) + pow(a/r, 2);
   curlyG = 1 - 2/r + a/pow(r, 1.5);
   return( pow(r,0.5) * curlyF / curlyG );
}

FTYPE l3d(FTYPE r, FTYPE a, FTYPE rbreak1, FTYPE rbreak2, FTYPE xi) {
   return ( xi * lK( r<=rbreak1 ? rbreak1 : r>=rbreak2 ? rbreak2 : r , a) );
}

//FTYPE l3d(FTYPE lam, FTYPE a, FTYPE lambreak1, FTYPE lambreak2, FTYPE xi) {
//   return ( xi * lK( lam<=lambreak1 ? lambreak1 : lam>=lambreak2 ? lambreak2 : lam , a) );
//}

FTYPE rtbis(FTYPE (*func)(FTYPE,FTYPE*), FTYPE *parms, FTYPE x1, FTYPE x2, FTYPE xacc)
//Taken from HARM:nrutil.c
{
   int j;
   FTYPE dx,f,fmid,xmid,rtb;
   f=(*func)(x1, parms);
   fmid=(*func)(x2, parms);
   if (f*fmid >= 0.0) {
      printf( "f(%g)=%g f(%g)=%g\n", x1, f, x2, fmid );
      printf("Root must be bracketed for bisection in rtbis\n");
   }
   rtb = (f < 0.0) ? (dx=x2-x1,x1) : (dx=x1-x2,x2); //Orient the search so that f>0 lies at x+dx.
   for (j=1;j<=100;j++) {
      fmid=(*func)(xmid=rtb+(dx *= 0.5),parms); //Bisection loop.
      if (fmid <= 0.0) {
         rtb=xmid;
      }
      if (fabs(dx) < xacc || fmid == 0.0) {
         return rtb;
      }
   }
   printf("Too many bisections in rtbis\n");
   return 0.0; //Never get here.
}


// von Zeipel cylinder radius corresponding to equatorial radius
FTYPE Lam_of_rEq(FTYPE req, FTYPE a, FTYPE rbreak1, FTYPE rbreak2, FTYPE xi) {

   FTYPE gdtt,gdtp,gdpp,l,lam,lam2;
   compute_gd(req, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   l = l3d(req, a, rbreak1, rbreak2, xi);
   lam2 = -l*(l*gdtp+gdpp)/(l*gdtt+gdtp);
   lam = sqrt(lam2);

   return lam;

}
FTYPE rEqofLam_func(FTYPE req, FTYPE *parms) {
   FTYPE lamguess, lam, a, rbreak1, rbreak2, xi;

   lam       = parms[0];
   a         = parms[1];
   rbreak1 = parms[2];
   rbreak2 = parms[3];
   xi        = parms[4];

   lamguess = Lam_of_rEq(req, a, rbreak1, rbreak2, xi);
   
   return ( lam*lam - lamguess*lamguess);
}

// equatorial radius corresponding to von Zeipel cylinder radius
FTYPE rEq_of_Lam(FTYPE lam, FTYPE a, FTYPE rbreak1, FTYPE rbreak2, FTYPE xi) {

   FTYPE parms[5];

   //store params in an array before function call
   parms[0] = lam;
   parms[1] = a;
   parms[2] = rbreak1;
   parms[3] = rbreak2;
   parms[4] = xi;
   
   //solve for equatorial radius using bisection, specify large enough root search range, (1e-3, 1e3)
   //demand accuracy 5x machine prec.
   //in non-rel limit req = lambda , use 10x that as the upper limit:

   
   if(rEqofLam_func(.01*lam,parms)*rEqofLam_func(10*lam,parms)>0) 
      return -1.;
   //  printf("%e %e %e %e %e \n",lam,Lam_of_rEq(.01*lam, a, rbreak1, rbreak2, xi),Lam_of_rEq(10*lam, a, rbreak1, rbreak2, xi), rEqofLam_func(.01*lam,parms),rEqofLam_func(10*lam,parms));
   else return( rtbis( &rEqofLam_func, parms, .01*lam, 10*lam, 5.*DBL_EPSILON ) );
}



FTYPE lamBL_func(FTYPE lam, FTYPE *parms) {
   FTYPE r,req,th, gdtt, gdtp, gdpp, gdtpeq, a, rbreak1, rbreak2, xi, l;

   gdtt      = parms[0];
   gdtp      = parms[1];
   gdpp      = parms[2];
   a         = parms[3];
   rbreak1 = parms[4];
   rbreak2 = parms[5];
   xi        = parms[6];
   r         = parms[7];
   th        = parms[8];
   
   req = rEq_of_Lam(lam, a, rbreak1, rbreak2, xi);
   l = l3d(req, a, rbreak1, rbreak2, xi);
    
   return ( lam*lam + l * (l*gdtp + gdpp) / (l*gdtt + gdtp) );
}

// von Zeipel cylinder radius of an arbitrary r, theta
FTYPE lamBL(FTYPE r, FTYPE th, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp, FTYPE a, FTYPE rbreak1, FTYPE rbreak2, FTYPE xi) {
   FTYPE R = r*sin(th);//, used as initial guess for lamBL

   FTYPE parms[9];

   //store params in an array before function call
   parms[0] = gdtt;
   parms[1] = gdtp;
   parms[2] = gdpp;
   parms[3] = a;
   parms[4] = rbreak1;
   parms[5] = rbreak2;
   parms[6] = xi;
   parms[7] = r;
   parms[8] = th;
   
   //solve for lin using bisection, specify large enough root search range, (1e-3, 1e3)
   //demand accuracy 5x machine prec.
   //in non-rel limit rml = r , use 10x that as the upper limit:
   if(th==M_PI_2)
      return( Lam_of_rEq(r, a, rbreak1, rbreak2, xi) );
   else
      return( rtbis( &lamBL_func, parms, R, 10*R, 5.*DBL_EPSILON ) );
}

FTYPE omega3d( FTYPE l, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp ) {
   return( -(gdtt*l + gdtp)*pow(gdpp + gdtp*l,-1) );
}

FTYPE compute_Agrav( FTYPE om, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp ){
   return (sqrt(fabs(1./ ( gdtt + 2*om*gdtp + pow(om,2)*gdpp ) )));
}

/*
FTYPE rmidlam( FTYPE x, FTYPE *parms ) {
   FTYPE lamsq, gdtt, gdtp, gdpp, a, rbreak1, rbreak2, xi;
   FTYPE lam_x, ans;

   lamsq     = parms[0];   // square of target lambda
   gdtt      = parms[1];
   gdtp      = parms[2];
   gdpp      = parms[3];
   a         = parms[4];
   rbreak1 = parms[5];
   rbreak2 = parms[6];
   xi        = parms[7];

   compute_gd(x, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lam_x = lamBL(x, M_PI_2, gdtt, gdtp, gdpp, a, rbreak1, rbreak2, xi);   // lambda at current value of x

   ans = lamsq - lam_x*lam_x;

   return(ans);
}

FTYPE limotorus_findrml(FTYPE lam, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp, FTYPE a, FTYPE rbreak1, FTYPE rbreak2, FTYPE xi) {
   FTYPE parms[8];
   FTYPE rml;

   //store params in an array before function call
   parms[0] = lam*lam;
   parms[1] = gdtt;
   parms[2] = gdtp;
   parms[3] = gdpp;
   parms[4] = a;
   parms[5] = rbreak1;
   parms[6] = rbreak2;
   parms[7] = xi;

   //solve for lin using bisection, specify large enough root search range, (1e-3, 1e3)
   //demand accuracy 5x machine prec.
   //in non-rel limit rml = lam , use 10x that as the upper limit:
   rml = rtbis( &rmidlam, parms, 3, 10*lam, 5.*DBL_EPSILON );

   return( rml );
}
*/

struct lnf_int_params {
   FTYPE a, rbreak1, rbreak2, xi;
};

// Integrand in expression for lnf. This function is called
// by the GSL quadrature routine.
// all quantities are equatorial
FTYPE lnf_integrand(FTYPE r, void *params) {
   struct lnf_int_params *pars = (struct lnf_int_params *) params;
   FTYPE a = pars->a, rbreak1 = pars->rbreak1, rbreak2 = pars->rbreak2, xi = pars->xi;
   FTYPE rhalf, r2, r3, r4,a2, lam, lam2, lamroot;
   FTYPE gdtt, gdtp, gdpp;
   FTYPE l, om, dl_dr, dx_dr, dom_dr, integrand;
   FTYPE term1, term2, term3, term4, D, E;
   FTYPE oneplusx, dx_dlam, dlK_dlam, om_numerator, om_denominator;

   // all values below are midplane values
   rhalf = sqrt(r);
   r2 = r*r;
   r3 = r2*r;
   r4 = r3*r;
   
   a2 = a*a;
   //compute_gd(r, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   //lam = Lam_of_rEq(r, a, rbreak1, rbreak2, xi);
   //lam = lamBL(r, M_PI_2, gdtt, gdtp, gdpp, a, rbreak1, rbreak2, xi);
   //lam2 = lam*lam;

   l = l3d(r, a, rbreak1, rbreak2, xi);

   term1 = (r-2) * l;
   term2 = 2*a;
   term3 = r3 + a2*(r+2);
   term4 = term2 * l;

   om_numerator = term1 + term2;
   om_denominator = term3 - term4;
   om = om_numerator / om_denominator;

   //om = omega3d( l, gdtt, gdtp, gdpp )
     
   // derivatives
   if (r <= rbreak1 || r >= rbreak2) {
      dl_dr = 0;
   } else {
     //dl_dr_A = (a2*r2*rhalf + 4*a*r3 - 3*r3*rhalf +0.5*r4*rhalf);
     //dl_dr_B = (a*r - 2*r*rhalf + r2*rhalf);
     //dl_dr = xi * dl_dr_A / (dl_dr_B * dl_dr_B);

      oneplusx = 1 - 2/r + a/r/rhalf;
      dx_dr = 2/r2 - 1.5*a/r2/rhalf;
      dl_dr = (oneplusx*0.5/rhalf + (1-rhalf)*dx_dr) / (oneplusx*oneplusx);

      //old
      //oneplusx = 1 - 2/lam + a*pow(lam,-1.5);
      //dx_dlam = 2*pow(lam,-2) - 1.5*a*pow(lam,-2.5);
      //lamroot = sqrt(lam);
      //dlK_dlam = (oneplusx*0.5*pow(lamroot,-1) + (a-lamroot)*dx_dlam) / pow(oneplusx,2);

      //D = term3 - 2*term4 - lam2 * (r-2);
      //E = l * (3*r2 + a2 - lam2);
      //dl_dr = E / ( 2*lam*om_numerator / (xi*dlK_dlam) - D);
   }

   dom_dr = ( om_denominator * (l + (r-2)*dl_dr) - om_numerator * (3*r2+a2+2*a*dl_dr) )
            / pow(om_denominator, 2);

   integrand = -l/(1 - om*l) * dom_dr;

   return( integrand );
}

int init_dsandvels_limotorus(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell)
{
   FTYPE kappa;
   FTYPE hh, eps;
   FTYPE rho, u, ur, uh, up;
   int pl;
   FTYPE R, rbreak1, rbreak2, xi;
   //FTYPE lambreak1, lambreak2;
   FTYPE lam, rml, l, om, Agrav, f3d;
   FTYPE lamin, lin, omin, Agravin, f3din;
   FTYPE lnf3d, lnferr;
   FTYPE gdtt, gdtp, gdpp;

   // GSL quadrature stuff
   gsl_integration_workspace *intwork;
   int worksize = 1000;   // workspace size
   gsl_function integrand;
   int intstatus;
   R = r*sin(th);
   if (R < LT_RIN) {*rhoout = -1.; return(0);}

   ///
   /// Parameters
   ///

   kappa = LT_KAPPA;   // AKMARK: entropy constant that appears in EOS
   xi = LT_XI;   // AKMARK: omega is set to this fraction of Keplerian omega
   rbreak1 = LT_R1; // AKMARK: locations of breaks in torus angular momentum profile
   rbreak2 = LT_R2;

   ///
   /// Computations at break radii
   ///

   // AKMARK: lamBL can calculate lambda at an arbitrary location, but it needs lambreak1,2;
   // how to calculate lambreak1,2 themselves?
   // Solution: note that lambreak1,2 are only needed by "l3d" in order to determine which region of the torus we are in.
   // But lambreak1,2 can be considered to be in region 2, so just need a way to make "l3d" believe that we are in region 2.
   // This can be done by setting lambreak1,2 to very small and very large values respectively.

   //compute_gd(rbreak1, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   //lambreak1 = lamBL(rbreak1, gdtt, gdtp, gdpp, a, 0, 200000, xi);
   //compute_gd(rbreak2, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   //lambreak2 = lamBL(rbreak2, gdtt, gdtp, gdpp, a, lambreak1, 200000, xi);

   ///
   /// Computations at torus inner edge
   ///

   compute_gd(LT_RIN, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lamin = lamBL(LT_RIN, M_PI_2, gdtt, gdtp, gdpp, a, rbreak1, rbreak2, xi);
   lin = l3d(LT_RIN, a, rbreak1, rbreak2, xi);
   omin = omega3d(lin, gdtt, gdtp, gdpp);
   Agravin = compute_Agrav(omin, gdtt, gdtp, gdpp);
   f3din = 1.;   // the way f3d is defined, it equals 1 at r=rin

   ///
   /// Computations at current point: r, th
   ///

   compute_gd(r, th, a, &gdtt, &gdtp, &gdpp);
   lam = lamBL(r, th, gdtt, gdtp, gdpp, a, rbreak1, rbreak2, xi);
   rml =  rEq_of_Lam(lam, a, rbreak1, rbreak2, xi); 
   l = l3d(rml, a, rbreak1, rbreak2, xi);
   om = omega3d(l, gdtt, gdtp, gdpp);
   Agrav = compute_Agrav(om, gdtt, gdtp, gdpp);

   //rml = limotorus_findrml(lam, gdtt, gdtp, gdpp, a, rbreak1, rbreak2, xi );
   // f3d requires numerical evaluation of an integral

   // First, set up things that GSL quadrature routine requires
   struct lnf_int_params pars = {a, rbreak1, rbreak2, xi};
   integrand.function = &lnf_integrand;   // integrand: function
   integrand.params = &pars;   // other parameters to pass to integrand (besides integration variable)
   intwork = gsl_integration_workspace_alloc(worksize);   // integration workspace

   // then perform integration to obtain ln(f3d) ...
   intstatus = gsl_integration_qags(&integrand, LT_RIN, rml, 0., 1.e-8, worksize, intwork, &lnf3d, &lnferr);
   gsl_integration_workspace_free( intwork );
   if (intstatus != GSL_SUCCESS) {
      printf("GSL integration failed during limotorus setup at r=%21.15g, th=%21.15g; setting density to zero at this point", r, th);
      lnf3d = GSL_NEGINF;   // cause density to be 0 at the current point
   }

   // ... and finally, f3d
   f3d = exp(-lnf3d);

   hh = f3din*Agrav / (f3d*Agravin);
   eps = (-1 + hh)*pow(LT_GAMMA,-1);

   if (eps < 0) rho = -1; else rho = pow((-1 + LT_GAMMA)*eps*pow(kappa,-1),pow(-1 + LT_GAMMA,-1));

   *rhoout = rho;
   *ell=l;
   *uuout = kappa * pow(rho, LT_GAMMA) / (LT_GAMMA - 1.);

   //my_err("wfsg");
   //printf("%e %e %e\n",rho,kappa,*uuout);
   //exit(0);   getch();

   return(0);

}

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
