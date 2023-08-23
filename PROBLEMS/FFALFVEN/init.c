//sets the initial conditions
//called from a loop going over ix,iy,iz
// https://arxiv.org/pdf/astro-ph/0503420.pdf appendix

/***********************************************/
//structure of geometry
//holds coordinates and metric, see ko.h
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
ldouble pp[NV],uu[NV];

// state on the LHS
int i;

// definition of the state on the LHS
// all equations assume cartesian minkowski
#ifdef K07PROBLEM
ldouble Bcon_lhs[4] = {0.,10.,1.,0.};
ldouble ucon_lhs[4] = {1.,0.,0.,0.};

#else

//ldouble Bcon_lhs[4] = {0,3.,3.,0};
ldouble Bcon_lhs[4] = {0,7.,1.,0.};
ldouble ucon_lhs[4] = {1.,0.,0.,0.};


//ldouble Bcon_lhs[4] = {0.,0,2.3094,0.};
//ldouble ucon_lhs[4] = {0.,0.57735,0.,0.};
#endif


ucon_lhs[0] = sqrt(1 + dot3nr(ucon_lhs,ucon_lhs));
ldouble bcon_lhs[4],bcov_lhs[4];
ldouble bsq_lhs;
bcon_lhs[0] = dot3nr(Bcon_lhs, ucon_lhs);
bcov_lhs[0] = -bcon_lhs[0];
for(i=1;i<4;i++)
{
  bcon_lhs[i] = (Bcon_lhs[i] + ucon_lhs[i]*bcon_lhs[0])/ucon_lhs[0];
  bcov_lhs[i] = bcon_lhs[i];
}
bsq_lhs = dot(bcon_lhs, bcov_lhs);
if(PROCID==0 && ix==0) printf("bsq_lhs: %e \n", bsq_lhs);

// set density
// TODO: add density slope 
// bsq should be constant everywhere in this setup (??)
ldouble rho,ugas;
#ifndef SIGMAINIT
rho = RHOINIT;
ugas = UUINIT;
#else
rho=bsq_lhs/SIGMAINIT;
ugas=rho*THETAINIT/(GAMMA-1.)/MU_GAS;
#endif

// calculate alfven wavespeed
ldouble Ec;
Ec=rho + GAMMA*ugas + bsq_lhs;

ldouble mu, gammamu; // todo right sign in plusminus? 
mu = (bcon_lhs[0]*bcon_lhs[1] - Ec*ucon_lhs[0]*ucon_lhs[1]);
mu -= sqrt(Ec*pow(bcon_lhs[1]*ucon_lhs[0] - bcon_lhs[0]*ucon_lhs[1],2));
mu /= (bcon_lhs[0]*bcon_lhs[0] - Ec*ucon_lhs[0]*ucon_lhs[0]);
gammamu = 1/sqrt(1-mu*mu);

if(PROCID==0 && ix==0) printf("Alvfen wavespeed: %e %e \n", mu, gammamu);

// Lorentz boost backwards to the wave frame
ldouble ucon_wave_lhs[4],bcon_wave_lhs[4];

ucon_wave_lhs[0] = gammamu*ucon_lhs[0] - mu*gammamu*ucon_lhs[1];
ucon_wave_lhs[1] = gammamu*ucon_lhs[1] - mu*gammamu*ucon_lhs[0];
ucon_wave_lhs[2] = ucon_lhs[2];
ucon_wave_lhs[3] = ucon_lhs[3];

bcon_wave_lhs[0] = gammamu*bcon_lhs[0] - mu*gammamu*bcon_lhs[1];
bcon_wave_lhs[1] = gammamu*bcon_lhs[1] - mu*gammamu*bcon_lhs[0];
bcon_wave_lhs[2] = bcon_lhs[2];
bcon_wave_lhs[3] = bcon_lhs[3];

ldouble Bcon_wave_lhs[4];
for(i=0;i<4;i++)
{
  Bcon_wave_lhs[i] = bcon_wave_lhs[i]*ucon_wave_lhs[0] - bcon_wave_lhs[0]*ucon_wave_lhs[i];
}

// calculate constant parameters of the Alfven wave ellipse
ldouble chi,ay,az,c,d,dd,a11,a22,a33,a12,a13,a23,byc,bzc;
chi = ucon_wave_lhs[1]/bcon_wave_lhs[1];
ay = (ucon_wave_lhs[2]-chi*bcon_wave_lhs[2])/(ucon_wave_lhs[0]-chi*bcon_wave_lhs[0]);
az = (ucon_wave_lhs[3]-chi*bcon_wave_lhs[3])/(ucon_wave_lhs[0]-chi*bcon_wave_lhs[0]);
c = chi*bsq_lhs / (ucon_wave_lhs[0]-chi*bcon_wave_lhs[0]);
d = bsq_lhs - bcon_wave_lhs[1]*bcon_wave_lhs[1];
dd = (pow(bcon_wave_lhs[1],2) - bsq_lhs*pow(ucon_wave_lhs[1],2)) / pow(Bcon_wave_lhs[1],2);
a11 = 1 - ay*ay;
a22 = 1 - az*az;
a33 = -(c*c + d);
a12 = -ay*az;
a13 = -c*ay;
a23 = -c*az;
byc = c*ay/dd;
bzc = c*az/dd;

// calculate local rotation angle theta
// TODO -- offset from origin center
ldouble awidth = PULSEWIDTH;
ldouble x0 = -0.5*awidth;
ldouble x1 = 0.5*awidth;
ldouble theta_lhs, thetarot, theta;

#ifdef K07PROBLEM
theta_lhs = 0.;
thetarot = 2*M_PI;

if(geom.xx < x0)
{
  theta = theta_lhs;
}
else if(geom.xx > x1)
{
  theta = theta_lhs + thetarot;
}
else
{
  ldouble xi = (geom.xx - x0)/(x1-x0);
  theta = theta_lhs + thetarot*(3*xi*xi - 2*xi*xi*xi);
}

#else
theta_lhs = atan((bcon_wave_lhs[3]-bzc)/(bcon_wave_lhs[2]-byc));
thetarot = M_PI;
//thetarot = 0.3;
if(geom.xx < x0)
{
  theta = theta_lhs;
}
else if(geom.xx > x1)
{
 theta = theta_lhs + thetarot;
}
else
{
  theta = theta_lhs + thetarot*pow(sin(M_PI*(geom.xx + 0.5*awidth)/(2*awidth)),2);
 //theta = theta_lhs + thetarot*0.5*(1+sin(5*M_PI*geom.xx));//??????
}
#endif

// get the fields in the wave frame
ldouble cth,sth,byz;
ldouble bcon_wave[4], ucon_wave[4];
cth = cos(theta);
sth = sin(theta);
byz = sqrt(d + c*c/dd)/sqrt(a11*cth*cth + 2*a12*sth*cth + a22*sth*sth);

bcon_wave[1] = bcon_wave_lhs[1];
bcon_wave[2] = byc + byz*cth;
bcon_wave[3] = bzc + byz*sth;

ucon_wave[1] = ucon_wave_lhs[1];
ucon_wave[2] = ucon_wave_lhs[2] + chi*(bcon_wave[2] - bcon_wave_lhs[2]);
ucon_wave[3] = ucon_wave_lhs[3] + chi*(bcon_wave[3] - bcon_wave_lhs[3]);

ucon_wave[0] = sqrt(1 + dot3nr(ucon_wave,ucon_wave));
bcon_wave[0] = dot3nr(bcon_wave,ucon_wave)/ucon_wave[0];

// boost back to the lab frame
ldouble ucon[4],bcon[4],Bcon[4];

ucon[0] = gammamu*ucon_wave[0] + mu*gammamu*ucon_wave[1];
ucon[1] = gammamu*ucon_wave[1] + mu*gammamu*ucon_wave[0];
ucon[2] = ucon_wave[2];
ucon[3] = ucon_wave[3];

bcon[0] = gammamu*bcon_wave[0] + mu*gammamu*bcon_wave[1];
bcon[1] = gammamu*bcon_wave[1] + mu*gammamu*bcon_wave[0];
bcon[2] = bcon_wave[2];
bcon[3] = bcon_wave[3];

for(i=0;i<4;i++)
  Bcon[i] = bcon[i]*ucon[0] - bcon[0]*ucon[i];

// fill primitives
pp[RHO] = rho;
pp[UU] = ugas;
pp[B1] = Bcon[1];
pp[B2] = Bcon[2];
pp[B3] = Bcon[3];
pp[VX] = ucon[1];
pp[VY] = ucon[2];
pp[VZ] = ucon[3];

//printf("%d | %e\n",ix,pp[VX]*pp[B1]+pp[VY]*pp[B2]+pp[VZ]*pp[B3]);


////////////////////////////////////////////////////////////////////////////////////////
// modify rho for linear sigma slope or tanh
// currently working with ALFVEN problem
/*
#if defined(INIT_SIGMA_TANH) || defined(INIT_SIGMA_LIN)
//ldouble bsq=Bsq/(gamma*gamma); //this is true if gamma=gammaperp

ldouble sigma;
ldouble sigmaleft,sigmaright;
#ifdef INIT_SIGMA_HIGHLEFT
sigmaleft=SIGMAINIT;
sigmaright=SIGMAINITMIN;
#else
sigmaleft=SIGMAINITMIN;
sigmaright=SIGMAINIT;
#endif

#if defined(INIT_SIGMA_TANH)
ldouble tw=SIGMAINITWIDTH;
ldouble to=SIGMAINITOFFSET;
sigma = sigmaleft + (sigmaright-sigmaleft)*0.5*(tanh(tw*(geom.xx-to))+1);

#elif defined(INIT_SIGMA_LIN)

if(geom.xx<-0.5*LLL)
{
  sigma = sigmaleft;
}
else if(geom.xx>0.5*LLL)
{
  sigma = sigmaright;
}
else
{
  sigma = (sigmaright-sigmaleft)*(geom.xx-0.5*LLL)/(LLL) + sigmaright;
}

#endif
pp[RHO]=bsq/sigma;
pp[UU]=pp[RHO]*THETAINIT/(GAMMA-1.)/MU_GAS;

#endif // defined(INIT_SIGMA_TANH) || defined(INIT_SIGMA_LIN)
*/


// get entropy
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);

//convert primitives to conserved
p2u(pp,uu,&geom);

//save to memory
int iv;
for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

// set force-free flag

#ifdef FORCEFREE
set_u_scalar(ffinvarr, ix, iy, iz, 1.0);
set_cflag(FFINVFLAG, ix,iy,iz,1);
set_cflag(MHDINVFLAG,ix,iy,iz,0);
#endif
