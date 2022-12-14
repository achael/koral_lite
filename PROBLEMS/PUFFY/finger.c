//drives the density towards initial state
int ix,iy,iz,ii,iv;
/**************************/

#pragma omp parallel for private(ix,iy,iz,iv) schedule (dynamic)
for(ii=0;ii<Nloop_0;ii++) //domain only
  {
    ldouble pp[NV],pptorus[NV],uu[NV];
    struct geometry geom;
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2];
    fill_geometry(ix,iy,iz,&geom);

    ldouble xxBL[4];
    coco_N(geom.xxvec,xxBL,MYCOORDS, BLCOORDS);
    
    //Keplerian orbital period
    ldouble Omk = 1./(BHSPIN+sqrt(xxBL[1]*xxBL[1]*xxBL[1]));
    ldouble Pk = 2.*M_PI/Omk;
    
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,ix,iy,iz);
	pptorus[iv]=get_u(ptemp1,iv,ix,iy,iz);
      }

    //density driven toward initial 
    //temperature fixed
    ldouble drho=0.;
    ldouble duint=0;
    ldouble dvr=0.;
    ldouble dvth=0.;
    ldouble dvph=0.;

    #ifdef RESTORETORUS

    ldouble effscale=1.; //how fast restores
    if(pptorus[ENTR]>0.) //inside initial torus
      {
	drho =  - effscale
	  * global_dt / Pk 
	  * (pp[RHO]-pptorus[RHO]);

	//not to overshoot 0.
	while (drho+pp[RHO]<0.)
	  {
	    drho/=2.;
	  }

	//drive towards torus parameters
	duint = (pp[UU]*pp[RHO]+pptorus[UU]*drho)/(pp[RHO]+drho) - pp[UU];
	dvr = (pp[VX]*pp[RHO]+pptorus[VX]*drho)/(pp[RHO]+drho) - pp[VX];
	dvth = (pp[VY]*pp[RHO]+ pptorus[VY]*drho)/(pp[RHO]+drho) - pp[VY];
	dvph = (pp[VZ]*pp[RHO] + pptorus[VZ]*drho)/(pp[RHO]+drho) - pp[VZ];
      }

    //if(iy==NY/2 && ix==NX/2)
    //printf("%e %e %e\n",drho,pp[RHO],pptorus[RHO]);

#endif

    //modifying only density, leaving temperature and velocities and magn. field intact
    pp[RHO]+=drho;
    pp[UU]+=duint;
    pp[VX]+=dvr;
    pp[VY]+=dvth;
    pp[VZ]+=dvph;

    p2u(pp,uu,&geom);
    
    for(iv=0;iv<NV;iv++)
      {
	set_u(u,iv,ix,iy,iz,uu[iv]);
	set_u(p,iv,ix,iy,iz,pp[iv]);
      }
}
