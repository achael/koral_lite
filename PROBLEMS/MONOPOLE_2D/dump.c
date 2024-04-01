//v1 - 24
//..
//v7 - 30

ldouble bcond[4],bcovd[4],bsqd;


#ifdef MAGNFIELD
calc_bcon_prim(pp,bcond,&geom);
indices_21(bcond,bcovd,gg); 
bsqd = dot(bcond,bcovd);
v1=bsqd/2.;

v2=calc_divB(ix,iy,iz);
#endif
