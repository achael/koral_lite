struct struct_of_state state;
get_state(pp,&geom,&state);


v1=state.Tgas;
v2=state.Trad;
v3=state.kappaes;
v4=state.kappa;
v5=state.cs;
v7=state.Ehat;     

ldouble radlum,totlum;
calc_lum(xxvecout[1],3,&radlum,&totlum);
v6 = 1.e-40+radlum/calc_lumEdd()*(rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)*velGU2CGS(1.)*velGU2CGS(1.));

v8=state.entr;
v9=state.radentr;
v10=state.Gi[0];
v11=state.Gic[0];
