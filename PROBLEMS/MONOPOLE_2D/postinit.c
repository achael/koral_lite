//scales enthalpy so sigma = sigmaconst everywhere
#if defined(SIGMAWCONSTINIT) && defined(FORCEFREE)

int ix,iy,iz;

#ifdef MAGNFIELD
if(PROCID==0) {printf("Renormalizing magnetic field... ");fflush(stdout);}

#pragma omp parallel for private(ix,iy,iz) schedule (dynamic)
for(iz=0;iz<NZ;iz++)
  {
    for(iy=0;iy<NY;iy++)
      {
	for(ix=0;ix<NX;ix++)
	  {
	    /***********************************************/
	    ldouble pp[NV],uu[NV];
	    struct geometry geom;
	    fill_geometry(ix,iy,iz,&geom);

	    struct geometry geomBL;
	    fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
 
	    int iv;
	    PLOOP(iv)
	      pp[iv]=get_u(p,iv,ix,iy,iz);

	    ldouble bcon[4],bcov[4],bsq;
	    calc_bcon_prim(pp,bcon,&geom);
	    indices_21(bcon,bcov,geom.gg); 
	    bsq = dot(bcon,bcov);
	    ldouble enthalpy=pp[UU]*GAMMA + pp[RHO];
	    
	    ldouble fac = bsq/enthalpy/SIGMAWCONSTINIT;

	    pp[UU] *= fac;
	    pp[RHO] *= fac;
	       
	    p2u(pp,uu,&geom);

	    PLOOP(iv)
	    {
	      set_u(p,iv,ix,iy,iz,pp[iv]);
	      set_u(u,iv,ix,iy,iz,uu[iv]);
	    }
	  }
      }
  }


if(PROCID==0) {printf("done!\n");fflush(stdout);}
#endif
#endif
