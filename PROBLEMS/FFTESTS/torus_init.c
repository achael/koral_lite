//torus parameters
#define MASS 10.
#define TORUSENTR 1.e2
#define GAMMA (5./3.) //gamma of gas itself
#define EFFGAMMA (4./3.) //of mixture with opt. thick radiation

//physics in cgs
#define GGG (6.674e-8)
#define CCC (2.998e10)

//GM/C2 for Msun in cm
#define MSUNCM 147700.
#define GMC2 (MSUNCM*MASS)
#define GMC3 (GMC2/CCC0)

//conversions to/from G=M=c=1
#define tempCGS2GU(x)    (x)
#define tempGU2CGS(x)    (x)
#define lenCGS2GU(x)    (x/MASSCM)
#define lenGU2CGS(x)    (x*MASSCM)
#define timeCGS2GU(x)   (x/MASSCM*CCC)
#define timeGU2CGS(x)   (x*MASSCM/CCC)
#define velCGS2GU(x)    (x/CCC)
#define velGU2CGS(x)    (x*CCC)
#define rhoCGS2GU(x)    (x*GGG/CCC/CCC*MASSCM*MASSCM)
#define rhoGU2CGS(x)    (x/GGG*CCC*CCC/MASSCM/MASSCM)
#define surfdensCGS2GU(x)    (x*GGG/CCC/CCC*MASSCM)
#define surfdensGU2CGS(x)    (x/GGG*CCC*CCC/MASSCM)
#define massCGS2GU(x)    (x*GGG/CCC/CCC/MASSCM)
#define massGU2CGS(x)    (x/GGG*CCC*CCC*MASSCM)
#define kappaCGS2GU(x)  (x/GGG*CCC*CCC/MASSCM)
#define kappaGU2CGS(x)  (x*GGG/CCC/CCC*MASSCM)
#define endenCGS2GU(x) (x*GGG*MASSCM*MASSCM/CCC/CCC/CCC/CCC)
#define endenGU2CGS(x) (x/GGG/MASSCM/MASSCM*CCC*CCC*CCC*CCC)
#define fluxCGS2GU(x) (x*GGG*MASSCM*MASSCM/CCC/CCC/CCC/CCC/CCC)
#define fluxGU2CGS(x) (x/GGG/MASSCM/MASSCM*CCC*CCC*CCC*CCC*CCC)

//constants
#define K_BOLTZ (1.3806488e-16 * GGG / CCC / CCC / CCC / CCC / MASSCM)
#define M_PROTON massCGS2GU(1.67262158e-24)
#define M_ELECTR massCGS2GU(9.11e-28)
#define SIGMA_RAD (5.67e-5 * GGG / CCC / CCC / CCC / CCC / CCC * MASSCM * MASSCM)
#define A_RAD (4.*SIGMA_RAD)
#define MU_GAS 1.
#define Z_RATIO (1.0)
#define Pi (3.141592654)     
#define KAPPA_ES_COEFF (kappaCGS2GU(0.4))

//the torus

//angular momentum of the torus:
ldouble L=3.8;

//geom.gg[0][0] = g_tt, etc.
ldouble W = (1./2.)*log(-(geom.gg[0][0]*geom.gg[3][3])/(geom.gg[3][3]+L*L*geom.gg[0][0]));

//specific enthalpy at the inner edge - determines the location of the inner edge
ldouble Win = -0.049;

//some more stuff
ldouble w=exp(-(W-Win));
ldouble epsilon = (w-1.)/EFFGAMMA;  
ldouble effgamma=EFFGAMMA;

if(epsilon>0.) //interior of the torus
  {
    ldouble kappa = TORUSENTR; 
    //density
    ldouble rho=powl((effgamma-1)*epsilon/kappa,1/(effgamma-1)); 
    //internal energy density, pressure=(gamma-1.)*uint
    ldouble uint = kappa * pow(rho, effgamma) / (effgamma - 1.); 
    //angular velocity
    ldouble omega = -L*(geom.gg[0][0]/geom.gg[3][3]);
    //u^phi/u^t
    ldouble OMEGA1=sqrt(-omega*omega/(geom.gg[0][0]+omega*omega*geom.gg[3][3]));
    
    //radiative part
    //taking original torus pressure and distributing it between gas and radiation to satisfy local thermal equilibrium:
    //solving for T satisfying P=pgas+prad=bbb T + aaa T^4

    ldouble P,aaa,bbb;
    P=(effgamma-1.)*uint;
    aaa=4.*SIGMA_RAD/3.;
    bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
    ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));
    //effective temperature of radiation
    ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;
    //radiation energy density from E=sigma T4^4
    ldouble E=calc_LTE_EfromT(T4);  
    //mulitply by g_rr
    E*=geomBL.gg[1][1];
    //fluid frame radiative fluxes
    ldouble Fx,Fy,Fz;
    Fx=Fy=Fz=0.;
    
    //then calculate new gas pressure from (rho,T4)
    uint=calc_PEQ_ufromTrho(T4,rho);

    //the end
  }
