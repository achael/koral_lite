//scales magnetic pressure so MAXBETA = pmag/pgas
int ix,iy,iz;

#ifdef MAGNFIELD
if(PROCID==0) {printf("Renormalizing magnetic field... ");fflush(stdout);}
ldouble maxbeta=0.;
ldouble maxpmag=0;
ldouble maxpgas=0.;

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

	    ldouble bcon[4],bcov[4],bsq,pmag;
	    calc_bcon_prim(pp,bcon,&geom);
	    indices_21(bcon,bcov,geom.gg); 
	    bsq = dot(bcon,bcov);
	    pmag = bsq/2.;
	    
	    ldouble pgas=GAMMAM1*pp[UU];

#ifdef BETANORMFULL
		//normalizing wrt everywhere
#pragma omp critical
		if(pgas/pmag>maxbeta) 
		  {
		    maxbeta=pgas/pmag;
		    //printf("%d %d > %e %e %e\n",ix,iy,pmag,pgas,maxbeta);
		  }
		if(pmag>maxpmag) 
		  {
		    maxpmag=pmag;
		    //printf("%d %d > %e %e %e\n",ix,iy,pmag,pgas,maxbeta);
		  }
		if(pgas>maxpgas) 
		  {
		    maxpgas=pgas;
		    //printf("%d %d > %e %e %e\n",ix,iy,pmag,pgas,maxbeta);
		  }
		    
#else //normalizing wrt to the equatorial plane
		if(geom.iy==NY/2)
		  {
#pragma omp critical
		    if(pgas/pmag>maxbeta)
		    {
		      maxbeta=pgas/pmag;
		      //printf("%d %d > %e %e %e\n",ix,iy,pmag,pgas,maxbeta);
                    }
		    if(pmag>maxpmag) 
		    {
		      maxpmag=pmag;
		      //printf("%d %d > %e %e %e\n",ix,iy,pmag,pgas,maxbeta);
		    }
		    if(pgas>maxpgas) 
		    {
		      maxpgas=pgas;
		      //printf("%d %d > %e %e %e\n",ix,iy,pmag,pgas,maxbeta);
		    }

		  }
#endif
	  }
      }
  }

  ldouble global_maxbeta=maxbeta;
  ldouble global_maxpgas=maxpgas;
  ldouble global_maxpmag=maxpmag;

//choosing maximal maxbeta from the whole domain
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Allreduce(&maxbeta, &global_maxbeta, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(&maxpgas, &global_maxpgas, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(&maxpmag, &global_maxpmag, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);

  //maxbeta=global_maxbeta;
#endif

maxbeta = global_maxpgas/global_maxpmag;
ldouble fac=sqrt(maxbeta/MAXBETA);

//manual normalization - verify beta!
#ifdef BETANORMFACTOR
fac=BETANORMFACTOR;
#endif
//printf("rescaling magn.fields by %e (%e)\n",fac,maxbeta);

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
	    int iv;

	    PLOOP(iv)
	      pp[iv]=get_u(p,iv,ix,iy,iz);
	    
	    pp[B1]*=fac;
	    pp[B2]*=fac;
	    pp[B3]*=fac;

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
	    
