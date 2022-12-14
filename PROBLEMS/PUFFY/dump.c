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

#ifdef RADIATION
ldouble Rtt,Ehat,ucon[4],prad;
calc_ff_Rtt(pp,&Rtt,ucon,&geom);
Ehat=-Rtt; 
prad=Ehat/3.;
v3=prad;

//tau in radius
ldouble tau[3];
calc_tautot(pp,xxvec,dx,tau);
v4=tau[1];
calc_tauabs(pp,xxvec,dx,tau);
v5=tau[1];
#endif

//a'la plasma beta : pmag/(pgas+prad)
v6=v1/(GAMMAM1*pp[UU] + v3);


