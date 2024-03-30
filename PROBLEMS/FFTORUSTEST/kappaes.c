//scattering
ldouble OpacCorrection = 1;

#ifdef SKIPESOPACITY
return 0.;
#endif

#ifndef SCATTERING
return 0.;
#endif

#ifdef KLEINNISHINA
ldouble Tkn;
Tkn = sqrt(Trad*Trad + Te*Te);
OpacCorrection = 1/(1 + pow(Tkn/4.5e8,0.86) );
#endif

return  kappaCGS2GU(0.2*(1. + HFRAC)*OpacCorrection)*rho;

