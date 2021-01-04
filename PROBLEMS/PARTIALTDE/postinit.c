ldouble Mdot = DISKVR*DISKSIGMA*2.*M_PI*ROUT;
ldouble mdotscale = rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.);
if(PROCID==0) printf("Mdot_in at the edge: %f Mdot_Edd\n",-Mdot*mdotscale/calc_mdotEdd());
