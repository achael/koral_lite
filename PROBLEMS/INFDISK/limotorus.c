#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#define FTYPE double

#define LT_KAPPA 2.e3
#define LT_XI 0.91
#define LT_R1 30.
#define LT_R2 200.
#define LT_GAMMA 4./3.
#define LT_RIN 15.

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

FTYPE l3d(FTYPE lam, FTYPE a, FTYPE lambreak1, FTYPE lambreak2, FTYPE xi) {
   return ( xi * lK( lam<=lambreak1 ? lambreak1 : lam>=lambreak2 ? lambreak2 : lam , a) );
}

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

FTYPE lamBL_func(FTYPE lam, FTYPE *parms) {
   FTYPE gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi, l;

   gdtt      = parms[0];
   gdtp      = parms[1];
   gdpp      = parms[2];
   a         = parms[3];
   lambreak1 = parms[4];
   lambreak2 = parms[5];
   xi        = parms[6];

   l = l3d(lam, a, lambreak1, lambreak2, xi);

   return ( lam*lam + l * (l*gdtp + gdpp) / (l*gdtt + gdtp) );
}

// von Zeipel cylinder radius (needs to be calculated iteratively for nonzero spin)
FTYPE lamBL(FTYPE R, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp, FTYPE a, FTYPE lambreak1, FTYPE lambreak2, FTYPE xi) {
   // R = r*sin(th), used as initial guess for lamBL

   FTYPE parms[7];

   //store params in an array before function call
   parms[0] = gdtt;
   parms[1] = gdtp;
   parms[2] = gdpp;
   parms[3] = a;
   parms[4] = lambreak1;
   parms[5] = lambreak2;
   parms[6] = xi;

   //solve for lin using bisection, specify large enough root search range, (1e-3, 1e3)
   //demand accuracy 5x machine prec.
   //in non-rel limit rml = r , use 10x that as the upper limit:
   return( rtbis( &lamBL_func, parms, R, 10*R, 5.*DBL_EPSILON ) );

}

FTYPE omega3d( FTYPE l, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp ) {
   return( -(gdtt*l + gdtp)*pow(gdpp + gdtp*l,-1) );
}

FTYPE compute_Agrav( FTYPE om, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp ){
   return (sqrt(fabs(1./ ( gdtt + 2*om*gdtp + pow(om,2)*gdpp ) )));
}

FTYPE rmidlam( FTYPE x, FTYPE *parms ) {
   FTYPE lamsq, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi;
   FTYPE lam_x, ans;

   lamsq     = parms[0];   // square of target lambda
   gdtt      = parms[1];
   gdtp      = parms[2];
   gdpp      = parms[3];
   a         = parms[4];
   lambreak1 = parms[5];
   lambreak2 = parms[6];
   xi        = parms[7];

   compute_gd(x, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lam_x = lamBL(x, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi);   // lambda at current value of x

   ans = lamsq - lam_x*lam_x;

   return(ans);
}

FTYPE limotorus_findrml(FTYPE lam, FTYPE gdtt, FTYPE gdtp, FTYPE gdpp, FTYPE a, FTYPE lambreak1, FTYPE lambreak2, FTYPE xi) {
   FTYPE parms[8];
   FTYPE rml;

   //store params in an array before function call
   parms[0] = lam*lam;
   parms[1] = gdtt;
   parms[2] = gdtp;
   parms[3] = gdpp;
   parms[4] = a;
   parms[5] = lambreak1;
   parms[6] = lambreak2;
   parms[7] = xi;

   //solve for lin using bisection, specify large enough root search range, (1e-3, 1e3)
   //demand accuracy 5x machine prec.
   //in non-rel limit rml = r , use 10x that as the upper limit:
   rml = rtbis( &rmidlam, parms, 6, 10*parms[0], 5.*DBL_EPSILON );

   return( rml );
}

struct lnf_int_params {
   FTYPE a, lambreak1, lambreak2, xi;
};

// Integrand in expression for lnf. This function is called
// by the GSL quadrature routine.
FTYPE lnf_integrand(FTYPE r, void *params) {
   struct lnf_int_params *pars = (struct lnf_int_params *) params;
   FTYPE a = pars->a, lambreak1 = pars->lambreak1, lambreak2 = pars->lambreak2, xi = pars->xi;
   FTYPE r2, r3, a2, lam, lam2, lamroot;
   FTYPE gdtt, gdtp, gdpp;
   FTYPE l, om, dl_dr, dom_dr, integrand;
   FTYPE term1, term2, term3, term4, D, E;
   FTYPE oneplusx, dx_dlam, dlK_dlam, om_numerator, om_denominator;

   // all values below are midplane values
   r2 = pow(r,2);
   r3 = pow(r,3);
   a2 = a*a;
   compute_gd(r, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lam = lamBL(r, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi);
   lam2 = lam*lam;

   l = l3d(lam, a, lambreak1, lambreak2, xi);

   term1 = (r-2) * l;
   term2 = 2*a;
   term3 = r3 + a2*(r+2);
   term4 = term2 * l;

   om_numerator = term1 + term2;
   om_denominator = term3 - term4;
   om = om_numerator / om_denominator;

   // derivatives
   if (lam <= lambreak1 || lam >= lambreak2) {
      dl_dr = 0;
   } else {
      oneplusx = 1 - 2/lam + a*pow(lam,-1.5);
      dx_dlam = 2*pow(lam,-2) - 1.5*a*pow(lam,-2.5);
      lamroot = sqrt(lam);
      dlK_dlam = (oneplusx*0.5*pow(lamroot,-1) + (a-lamroot)*dx_dlam) / pow(oneplusx,2);
      D = term3 - 2*term4 - lam2 * (r-2);
      E = l * (3*r2 + a2 - lam2);
      dl_dr = E / ( 2*lam*om_numerator / (xi*dlK_dlam) - D);
   }

   dom_dr = ( om_denominator * (l + (r-2)*dl_dr) - om_numerator * (3*r2+a2+2*a*dl_dr) )
            / pow(om_denominator, 2);

   integrand = -l/(1-om*l) * dom_dr;

   return( integrand );
}

int init_dsandvels_limotorus(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell)
{
   FTYPE kappa;
   FTYPE hh, eps;
   FTYPE rho, u, ur, uh, up;
   int pl;
   FTYPE R, rbreak1, rbreak2, xi;
   FTYPE lambreak1, lambreak2;
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
   if (R < LT_RIN) {*rhoout = 0.; return(0);}

   ///
   /// Parameters
   ///

   kappa = LT_KAPPA;   // AKMARK: entropy constant that appears in EOS
   xi = LT_XI;   // AKMARK: omega is set to this fraction of Keplerian omega
   rbreak1 = LT_R1; // AKMARK: locations of breaks in torus angular momentum profile
   rbreak2 = LT_R2;

   //a=-0.7
   //rbreak1 = 44.481;
   //rbreak2 = 1000.;


   ///
   /// Computations at break radii
   ///

   // AKMARK: lamBL can calculate lambda at an arbitrary location, but it needs lambreak1,2;
   // how to calculate lambreak1,2 themselves?
   // Solution: note that lambreak1,2 are only needed by "l3d" in order to determine which region of the torus we are in.
   // But lambreak1,2 can be considered to be in region 2, so just need a way to make "l3d" believe that we are in region 2.
   // This can be done by setting lambreak1,2 to very small and very large values respectively.
   compute_gd(rbreak1, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lambreak1 = lamBL(rbreak1, gdtt, gdtp, gdpp, a, 0, 200000, xi);
   compute_gd(rbreak2, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lambreak2 = lamBL(rbreak2, gdtt, gdtp, gdpp, a, lambreak1, 200000, xi);

   ///
   /// Computations at torus inner edge
   ///

   compute_gd(LT_RIN, M_PI_2, a, &gdtt, &gdtp, &gdpp);
   lamin = lamBL(LT_RIN, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi);
   lin = l3d(lamin, a, lambreak1, lambreak2, xi);
   omin = omega3d(lin, gdtt, gdtp, gdpp);
   Agravin = compute_Agrav(omin, gdtt, gdtp, gdpp);
   f3din = 1.;   // the way f3d is defined, it equals 1 at r=rin

   ///
   /// Computations at current point: r, th
   ///

   compute_gd(r, th, a, &gdtt, &gdtp, &gdpp);
   lam = lamBL(R, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi);
   l = l3d(lam, a, lambreak1, lambreak2, xi);
   om = omega3d(l, gdtt, gdtp, gdpp);
   Agrav = compute_Agrav(om, gdtt, gdtp, gdpp);

   rml = limotorus_findrml( lam, gdtt, gdtp, gdpp, a, lambreak1, lambreak2, xi );

   // f3d requires numerical evaluation of an integral

   // First, set up things that GSL quadrature routine requires
   struct lnf_int_params pars = {a, lambreak1, lambreak2, xi};
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

   if (eps < 0) rho = 0; else rho = pow((-1 + LT_GAMMA)*eps*pow(kappa,-1),pow(-1 + LT_GAMMA,-1));

   *rhoout = rho;
   *ell=l;
   *uuout = kappa * pow(rho, LT_GAMMA) / (LT_GAMMA - 1.);

   return(0);

}



int main() {

   int nr = 30, nth = 16;
   FTYPE Rin = 2., Rout = 500.;   // a = 0
   //FTYPE Rin = 1.364, Rout = 2000.;   // a = 0.9
   FTYPE th1 = 0., th2 = M_PI_2;
   FTYPE dr = (Rout - Rin) / nr, dth = (th2 - th1) / nth;
   FTYPE factor;
   FTYPE r, th, rho, uu, ell;
   int i, j;
   FILE *outfile, *radfile;

//    FTYPE R0 = 1.05, startx1 = -0.6736, dx1 = 0.0372985;

   FTYPE a = 0.9;
   factor = log(Rout/Rin);

   outfile = fopen("slice.dat", "w");
   radfile = fopen("radslice.dat","w");

   FTYPE Sigma;

   for (i=0; i<=nr; i++) 
     {
       r = Rin*exp((FTYPE)i/nr*factor);
       printf("r: %.2f\n",r);

       Sigma=0.;
       for (j=0; j<=nth; j++) 
	 {
	   th = th1 + j*dth;
	   if(th>M_PI) th-=0.01;

	   init_dsandvels_limotorus(r, th, a, &rho, &uu, &ell);

	   Sigma+=rho*r*dth;

	   fprintf(outfile, "%g\t%g\t%g\t%g\t%g\n", r*sin(th), r*cos(th), rho, uu,ell);
	 }
      fprintf(outfile, "\n");
      fprintf(radfile,"%g %g\n",r,Sigma);

     }

   fclose(outfile);
   fclose(radfile);

   return(0);

}

//    intstatus = gsl_integration_qag(&integrand, lamin, lam, 1.e-7, 1.e-7, worksize, GSL_INTEG_GAUSS61, intwork, &lnf3d, &lnferr);
