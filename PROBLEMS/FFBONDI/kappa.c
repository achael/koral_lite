kappa=calc_opacities_01(rho, Te, Ti, Tgas, Trad, pp, ggg, opac);

/*
//absorption opacities

ldouble rhocgs=rhoGU2CGS(rho);
Trad=Tgas;

#ifdef NCOMPTONIZATION //the color temperature of radiation
Trad = calc_ncompt_Thatrad_full(pp,ggg);
#endif

ldouble zeta = Trad/Tgas;
//test
//zeta=1.;

ldouble mpcgs=1.67262158e-24;

#ifndef SKIPFANCYOPACITIES		   

//absorbtion mean
ldouble kappaff,kappabe;
kappaff=kappaCGS2GU((6.6e-24/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas)*log(1.+1.6*zeta))*rho*(1.+4.4e-10*Tgas);
kappabe=kappaCGS2GU((5.0e-15/mpcgs/mpcgs)*rhocgs*pow(Tgas,-1.7)/Tgas/Tgas/Tgas*log(1.+1.6*zeta))*rho;

*kappagasAbs=kappaff+kappabe;

#ifdef OPACSKIPBE
*kappagasAbs=kappaff;
#endif

#ifdef OPACSMOOTH
ldouble Tgas0=3.e7;
ldouble zeta0=Trad/Tgas0;
//test
//zeta0=1.;
ldouble kappaff0=kappaCGS2GU((6.6e-24/mpcgs/mpcgs)*rhocgs/Tgas0/Tgas0/Tgas0/sqrt(Tgas0)*log(1.+1.6*zeta0))*rho*(1.+4.4e-10*Tgas0);
ldouble kappabe0=kappaCGS2GU((5.0e-15/mpcgs/mpcgs)*rhocgs*pow(Tgas0,-1.7)/Tgas0/Tgas0/Tgas0*log(1.+1.6*zeta0))*rho;

*kappagasAbs=kappaff+kappabe+10.*(kappaff0+kappabe0)*Tgas0*Tgas0*Tgas0*Tgas0/Tgas/Tgas/Tgas/Tgas;
#endif

*kapparadAbs=*kappagasAbs/zeta/zeta/zeta;
	
//Roseland mean - not used at all							  
*kappagasRos=kappaCGS2GU((2.1e-25/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas)*(1.-exp(-6.94*zeta)))*rho*(1.+4.4e-10*Tgas);
*kapparadRos=*kappagasAbs/zeta/zeta/zeta;

//default
kappa=*kappagasAbs;

//test
//*kappagasAbs=kappaff;
//*kapparadAbs=kappabe;

#else
//the simplest
kappa=kappaCGS2GU((6.6e-24/mpcgs/mpcgs)*rhocgs/Tgas/Tgas/Tgas/sqrt(Tgas))*rho*(1.+4.4e-10*Tgas);



#endif



*/
