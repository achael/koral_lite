#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>


struct tfun_params
{
  ldouble n,r,C4,C5,T_crit,Kadiab;
};


ldouble tfunc(ldouble t, void *params)
{
  struct tfun_params *p
    = (struct tfun_params *) params;
  ldouble n = p->n;
  ldouble r = p->r;
  ldouble C4 = p->C4;
  ldouble C5 = p->C5;

  ldouble fun = pow(1 + (1+n)*t,2)*(1 - 2./r + C4*C4/(pow(r,4)*pow(t,2*n))) - C5;
  return fun;
}

ldouble tfunc_deriv(ldouble t, void * params)
{
  struct tfun_params *p
    = (struct tfun_params *) params;
  ldouble n = p->n;
  ldouble r = p->r;
  ldouble C4 = p->C4;
  ldouble C5 = p->C5;

  ldouble dfun = 2*(1 + t +n*t)*((1+n)*pow(r,3)*(r-2) - C4*C4*(n+t*(n*n-1))*pow(t,-1-2*n));
  dfun /= pow(r,4);
  
  return dfun;
}

void tfunc_fdf(ldouble t, void *params, ldouble *y, ldouble *dy)
{
  struct tfun_params *p
    = (struct tfun_params *) params;
  ldouble n = p->n;
  ldouble r = p->r;
  ldouble C4 = p->C4;
  ldouble C5 = p->C5;

  ldouble fun = pow(1 + (1+n)*t,2)*(1 - 2./r + C4*C4/(pow(r,4)*pow(t,2*n))) - C5;
  ldouble dfun = 2*(1 + t + n*t)*((1+n)*pow(r,3)*(r-2) - C4*C4*(n+t*(n*n-1))*pow(t,-1-2*n));
  dfun = dfun/pow(r,4);

  *y = fun;
  *dy = dfun;
  return;
}
  
int init_bondi_hydro(void* params, ldouble *rhoout, ldouble *uuout, ldouble *urout)
{
  struct tfun_params *p
    = (struct tfun_params *) params;

  ldouble n = p->n;
  ldouble r = p->r;
  ldouble C4 = p->C4;
  ldouble C5 = p->C5;
  ldouble T_crit = p->T_crit;
  ldouble Kadiab=p->Kadiab;
  ldouble gamma = (1+n)/n;
  
  // function parameters
  gsl_function_fdf FDF;
  FDF.f = &tfunc;
  FDF.df = &tfunc_deriv;
  FDF.fdf = &tfunc_fdf;
  FDF.params = params;

  
  // root finder
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double t0, t = T_crit;
  
  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, t);

  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      t0 = t;
      t = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (t, t0, 0, 1e-3);

      //if (status == GSL_SUCCESS) printf ("Converged:\n");
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  // compute other quantities
  ldouble rho = pow(t/Kadiab, n);
  ldouble pgas = rho*t;
  ldouble uint = pgas/(gamma-1.);
  ldouble ur = -C4/pow(t,n)/pow(r,2);

  *rhoout = rho;
  *uuout = uint;
  *urout = ur;
}


int init_bondi_full(ldouble* pp, int ix, int iy, int iz)
{
 
  //geometries
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  struct geometry geomBL;
  fill_geometry_arb(ix,iy,iz,&geomBL,BONDICOORDS);

  ldouble r = geomBL.xx;
  ldouble th = geomBL.yy;

  // constants
  ldouble mdot = MDOTINIT;
  ldouble rc = RSONIC;
  ldouble gamma = GAMMA;
  ldouble n = 1./(gamma - 1);

  // values at the critical point
  ldouble ur_crit  = -sqrt(0.5/rc);
  ldouble cs_crit  =  sqrt(ur_crit*ur_crit / (1 - 3*ur_crit*ur_crit));
  ldouble rho_crit =  -mdot/(4*M_PI*rc*rc*ur_crit);
  ldouble T_crit = 1./(gamma/(cs_crit*cs_crit) - gamma/(gamma-1.));
  ldouble Kadiab = T_crit / pow(rho_crit, gamma-1.);

  // constants
  ldouble C4 = pow(Kadiab, n)*mdot/(4*M_PI);
  ldouble C5 = pow(1. + (1.+n)*T_crit, 2)*(1 - 2./rc + ur_crit*ur_crit);

  // struct
  struct tfun_params params = {n,r,C4,C5,T_crit,Kadiab};

  // hydro values
  ldouble rho,uint,ur;
  init_bondi_hydro(&params, &rho, &uint, &ur); 

  // fill prims 
  pp[RHO]=rho;
  pp[UU]=uint;
  pp[VX]=ur; 
  pp[VY]=0.;
  pp[VZ]=0.;

  // These will be made consistent with B direction after init
  #ifdef FORCEFREE
  pp[UUFF]=uint;
  pp[VXFF]=ur;
  pp[VYFF]=0.;
  pp[VZFF]=0.;
  #endif

  //entropy
  pp[5]=calc_Sfromu(pp[0],pp[1],geom.ix,geom.iy,geom.iz);

  #ifdef MAGNFIELD // set to zero at first
  pp[B1]=pp[B2]=pp[B3]=0.;
  #endif

  // magn field from direct monopole
  #ifdef MAGNFIELD 
  ldouble sigma_crit = SIGMACRIT;
  ldouble bfac = sqrt(sigma_crit*rho_crit);

  pp[B1]=bfac*pow(rc,2)/pow(r,2);
  pp[B2]=0.;
  pp[B3]=0.;
  #endif

  //transform primitives from BL to MYCOORDS
  trans_pmhd_coco(pp, pp, BONDICOORDS, MYCOORDS, geomBL.xxvec,&geomBL,&geom);

  return 0;
}

