//scattering
ldouble OpacCorrection = 1;


#ifdef KLEINNISHINA
ldouble Tkn;
Tkn = sqrt(Trad*Trad + Te*Te);
Tkn = T;
OpacCorrection = 1/(1 + pow(Tkn/4.5e8,0.86) );
#endif

return  0.;// kappaCGS2GU(0.2*(1. + HFRAC)*OpacCorrection)*rho;

