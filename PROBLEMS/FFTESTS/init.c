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

#if(FFPROBLEM==LINEARALFVEN)
ldouble rho=RHOINIT;
ldouble uint=UUINIT;
ldouble Bsq=SIGMAINIT*rho;
ldouble va = sqrt(Bsq)/sqrt(Bsq + rho + GAMMA*uint);

ldouble k = 2*M_PI/LLL;
ldouble omega = k*va;
ldouble Delta=1.e-3;

pp[B1] = Sqrt(Bsq);
pp[B2] = -Sqrt(Bsq)*Delta*Sin(geom.xx*k);
pp[B3] = Sqrt(Bsq)*Delta*(1+Cos(geom.xx*k));


#ifdef LINEARALFVENSINGLEPERIOD
if(geom.xx>0.5*LLL || geom.xx < -0.5*LLL)
{
  pp[B2]=0;
  pp[B3]=0;
}
#endif

/*
// these expressions for E work to give v when d is small
// otherwise there is a nonzero vx
Ex = 0.;
Ey = Sqrt(Bsq)*Delta*va*(1+Cos(geom.xx*k));
Ez = Sqrt(Bsq)*Delta*va*Sin(geom.xx*k);
*/

pp[VX] = 0;
pp[VY] = -va*pp[B2]/Sqrt(Bsq);
pp[VZ] = -va*pp[B3]/Sqrt(Bsq);

pp[RHO] = rho;
pp[UU] = uint;

ldouble ucon[4],ucov[4],bcon[4],bcov[4],bsq;
calc_ucon_ucov_from_prims(pp, &geom, ucon, ucov);
calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, &geom, bcon, bcov, &bsq);

ldouble w = rho + GAMMA*uint;
ldouble va_v2 = (bcon[1] + ucon[1]*sqrt(w+bsq))/(bcon[0] + ucon[0]*sqrt(w+bsq));
ldouble va_ff = (bcon[1] + ucon[1]*sqrt(bsq))/(bcon[0] + ucon[0]*sqrt(bsq));;
if(geom.ix==0 && PROCID==0)
{
  printf("linearalfven: %e %e %e %e\n",k,omega,2*M_PI/omega,2*M_PI/omega/DTOUT1);
  printf("alfvenspeeds: %e %e %e\n",va,va_v2,va_ff);
  //printf("%e %e %e %e \n",ucon[0],ucon[1],ucon[2],ucon[3]);
  //printf("%e %e\n",Bsq,bsq);
  //exit(-1);
}

#else // Komissarov 1999/2002 problems

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

// Alfven wave test, De Villiers and Hawley 2002
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
// currently working with ALFVEN problem

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
#endif // ndef SIGMAINIT
#endif // FFPROBLEM

////////////////////////////////////////////////////////////////////////////////////////
// modify rho for linear sigma slope

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
