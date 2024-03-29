//sets the initial conditions
//called from a loop going over ix,iy,iz

/***********************************************/
//structure of geometry
//holds coordinates and metric, see ko.h
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
ldouble pp[NV],uu[NV];

///////////////////////
// Komissarov 1999/2002 problems
//////////////////////


ldouble mu, gammamu;
ldouble Ex,Ey,Ez,Bx,By,Bz,Bsq,Esq,gamma;
ldouble xp, Epx, Epy, Epz, Bpx, Bpy, Bpz;

// Fast wave test, Komissarov 2002
#if(FFPROBLEM==FASTWAVE)
if(ix==0) printf("FASTWAVE problem\n");


mu=1;

Bx=1.;
Bz=0.;
if(geom.xx<-0.1)
  By=1.;
else if(geom.xx>=-0.1 && geom.xx<=0.1)
  By=1.-1.5*(geom.xx+0.1);
else if(geom.xx>0.1)
  By=0.7;

Ex = 0.;
Ey = 0.;
//Ez = 1.-mu*By; //McKinney06 pg 8 
Ez = -mu*By; 

// Alfven wave test komissarov 2004
#elif(FFPROBLEM==ALFVEN)
if(ix==0) printf("ALFVEN problem\n");


mu=-0.5; // wave speed (3-velocity)

gammamu=1./sqrt(1-mu*mu);
xp = geom.xx*gammamu; // Lorentz contraction 

// Wave frame fields
ldouble f=1+sin(5*M_PI*xp);

Bpx = 1.;
Bpy = 1.;
if(xp<-0.1)
  Bpz=1.;
else if(xp>=-0.1 && xp<=0.1)
  Bpz=1.0 + 0.15*f;
else if(xp>0.1)
  Bpz=1.3;

Epx = -Bpz;
Epy = 0.;
Epz = 1.0;

// Lorentz transformation to grid frame
Ex = Epx;
Ey = gammamu*(Epy + mu*Bpz);
Ez = gammamu*(Epz - mu*Bpy);

Bx = Bpx;
By = gammamu*(Bpy - mu*Epz);
Bz = gammamu*(Bpz + mu*Epy);

// Degenerate Alfven wave test, Komissarov 2002
#elif(FFPROBLEM==ALFVEN2)
if(ix==0) printf("ALFVEN2 problem\n");

mu=0.5; // wave speed (3-velocity)

gammamu=1./sqrt(1-mu*mu);
xp = geom.xx*gammamu; // Lorentz contraction 

// Wave frame fields
ldouble phi;
if(xp<-0.1)
  phi=0.;
else if(xp>=-0.1 && xp<=0.1)
  phi=2.5*M_PI*(xp+0.1);
else if(xp>0.1)
  phi=0.5*M_PI;

Epx = 0;
Epy = 0;
Epz = 0;
Bpx = 0;
Bpy = 2*cos(phi);
Bpz = 2*sin(phi);

// Lorentz transformation to grid frame
Ex = Epx;
Ey = gammamu*(Epy + mu*Bpz);
Ez = gammamu*(Epz - mu*Bpy);

Bx = Bpx;
By = gammamu*(Bpy - mu*Epz);
Bz = gammamu*(Bpz + mu*Epy);


// Three wave test, Komissarov 2002
#elif(FFPROBLEM==THREEWAVE)
if(ix==0) printf("THREEWAVE problem\n");


if(geom.xx<0)
{
  Bx=1.0;
  By=1.5;
  Bz=3.5;
  Ex=-1.0;
  Ey=-0.5;
  Ez=0.5;
}
else if(geom.xx==0)
{
  Bx=0.;
  By=0.;
  Bz=0.;
  Ex=0.;
  Ey=0.;
  Ez=0.;
}
else
{
  Bx=1.;
  By=3.0;
  Bz=3.0;
  Ex=-1.5;
  Ey=2.0;
  Ez=-1.5;
}

//Force-free breakdown test, Komissarov 2002
#elif(FFPROBLEM==BREAKDOWN)
if(geom.xx<0)
{
  Bx=1.0;
  By=1.0;
  Bz=1.0;
  Ex=0.0;
  Ey=0.5;
  Ez=-0.5;
}
else if(geom.xx>0.2)
{
  Bx=1.0;
  By=-1.0;
  Bz=-1.0;
  Ex=0.0;
  Ey=0.5;
  Ez=-0.5;
}
else
{
  Bx=1.0;
  By=1.0 - 10*geom.xx;
  Bz=1.0 - 10*geom.xx;
  Ex=0.0;
  Ey=0.5;
  Ez=-0.5;
}


#else
printf("FFPROBLEM NOT RECOGNIZED!\n");
exit(-1);
#endif

// special relativistic problem, alpha=1 and gdet=1
// relationship from McKinney 06 eq 17 for velocity
// no difference between contravariant and covariant quantities

Bsq = Bx*Bx + By*By + Bz*Bz;
Esq = Ex*Ex + Ey*Ey + Ez*Ez;
gamma = sqrt(Bsq/(Bsq-Esq));

if(!isfinite(gamma))
{
  printf("init gamma not finite!\n");
  exit(-1);
}

pp[VX] = (gamma/Bsq)*(Ey*Bz - Ez*By);
pp[VY] = (gamma/Bsq)*(Ez*Bx - Ex*Bz);
pp[VZ] = (gamma/Bsq)*(Ex*By - Ey*Bx);

pp[B1] = Bx;
pp[B2] = By;
pp[B3] = Bz;

////////////////////////////////////////////////////////////////////////////////////////
// initial density set by sigma or rho,u

//printf("%d | %e\n",ix,pp[VX]*pp[B1]+pp[VY]*pp[B2]+pp[VZ]*pp[B3]);x
#ifndef SIGMAINIT

pp[RHO] = RHOINIT;
pp[UU] = UUINIT;

#else
ldouble b0 = pp[B1]*pp[VX]+pp[B2]*pp[VY]+pp[B3]*pp[VZ];
ldouble bsq = (Bsq+b0*b0)/(gamma*gamma);

pp[RHO]=bsq/SIGMAINIT;
pp[UU]=pp[RHO]*THETAINIT/(GAMMA-1.)/MU_GAS;

////////////////////////////////////////////////////////////////////////////////////////
// modify rho for linear sigma slope or tanh

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

#endif // if defined(INIT_SIGMA_TANH) || defined(INIT_SIGMA_LIN)
#endif // if ndef SIGMAINIT


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
