//absorption

ldouble rhocgs=rhoGU2CGS(rho);
ldouble Tcgs=tempGU2CGS(Tgas);
ldouble kappaffcgs=6.4e22*rhocgs/Tcgs/Tcgs/Tcgs/sqrt(Tcgs);

//test
//return 0.;
kappa= kappaCGS2GU(kappaffcgs)*rho;

//ldouble kappabfcgs=4.8e-24/1.67262158e-24/1.67262158e-24*rhocgs/Tcgs/Tcgs/Tcgs/sqrt(Tcgs)*ZZsun;
//return (kappaCGS2GU(kappaffcgs)*rho+kappaCGS2GU(kappabfcgs)*rho);

