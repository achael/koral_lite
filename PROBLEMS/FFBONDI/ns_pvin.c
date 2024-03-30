
	//take from the cell symmetrical wrt. the boundary!
	iix=-ix-1;

	//copy primatives from the symmetric cell
    //only want the scalar primatives from symmetric cell
	for(iv=0;iv<NV;iv++)
	  {
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }

    // get geometry for the current ghost cell
    struct geometry geomGC;
    fill_geometry(ix, iiy, iiz, &geomGC);

    struct geometry geomBLGC;
    fill_geometry_arb(ix, iiy, iiz, &geomBLGC, BLCOORDS);


    // get appropriate geometries for symmetric cell 
    struct geometry geomD, geomBLD; 
    fill_geometry(iix,iiy,iiz,&geomD);
    fill_geometry_arb(iix,iiy,iiz,&geomBLD, BLCOORDS);  

    // get radius and velocity at cell 0 by transformiong ppat0 to BLCOORDS
    trans_pall_coco(pp, pp, MYCOORDS, BLCOORDS, geomD.xxvec, 
            &geomD, &geomBLD);
    ldouble ratD = geomBLD.xx;
    ldouble vatD = pp[VX];
    #ifdef NS_OMEGA
    ldouble omegaatD = pp[VZ];
    #endif
    ldouble fatD = pp[FX];
    //printf("at 0, F = %f\n", fat0);
    // desired vel at radius at boundary
    ldouble ratb = RMIN;
    ldouble vatb = -V_IN;
    ldouble fatb = -F_IN;
    #ifdef NS_OMEGA
    ldouble omegaatb = NS_OMEGA;
    #endif
    
    //calculate slope
    ldouble m = (vatD - vatb)/(ratD - ratb);
    ldouble mf = (fatD - fatb)/(ratD - ratb);
    #ifdef NS_OMEGA
    ldouble momega = (omegaatD - omegaatb)/(ratD - ratb); 
    #endif
    //printf("m,mf = %f, %f\n", m,mf);

    // get r for the ghost cell
    ldouble rgc = geomBLGC.xx;
    
    // extrapolate velocity to the ghost using slope and y-intercept
    ldouble vgc = m*(rgc - ratb) + vatb;
    ldouble fgc = mf*(rgc - ratb) + fatb;
    
    #ifdef NS_OMEGA
    ldouble omegagc = momega*(rgc - ratb) + omegaatb;
    #endif
    //fgc*=pp[EE];// multiply by the symmetric Ehat
    //printf("F = %e\n", fgc);

	// transform to BL in cell 0 geometry
    //trans_pall_coco(pp, pp, MYCOORDS,BLCOORDS, geomD.xxvec,&geomD,&geomDBL);

	// set VX to the extrapolated value, and copy for the flux
    // maybe it would be good to extrapolate flux seperately
	pp[VX] = vgc;
    if(vatD>0.) pp[VX]=0.;
    
   
    // copy tangential velocities note they should already be in BLCOORDS
    pp[VY] = pp[VY];
    pp[VZ] = pp[VZ];
    #ifdef NS_OMEGA
    pp[VZ] = omegagc;;
    #endif

    #ifdef NONROTATING
    pp[VZ]=0.;
    #endif
    
    #ifdef RADIATION
        pp[FX] = fgc;
	if (fatD>0.) pp[FX]=0.;
        pp[FY] = pp[FY];
        pp[FZ] = pp[FZ];
        //pp[FX] = 0.;
        //pp[FY] = 0.;
        //pp[FZ] = 0.;
        #ifdef FLUXEQLVEL
	        pp[FX] = pp[VX];
	        pp[FY] = pp[VY];
	        pp[FZ] = pp[VZ];
        #endif
	// allow radiation to flow in with the same BL velocity
	// do not all radiation to flow out
	// conserve ehat gdet accross the boundary
	#ifdef RAD_INFLOW
        // get primatives form cell 0
        ldouble pp0[NV];
	for(iv=0;iv<NV;iv++) 
	{
		pp0[iv] = get_u(p,iv, 0, iiy, iiz);
	}	
	//geom0,  should already be declared and filled from bc.c
    struct geometry geomBL0;
    fill_geometry_arb(0, iiy, iiz, &geomBL0, BLCOORDS);
	
	// transform to BL coordinates in zero geometry
    trans_pall_coco(pp0, pp0, MYCOORDS, BLCOORDS, geom0.xxvec, 
            &geom0, &geomBL0);
	//  copy BL velocity to primatives
	pp[FX] = pp0[FX];
	pp[FY] = pp0[FY];
        pp[FZ] = pp0[FZ];
	//no outflow
	if (pp0[FX]>0.)
	{
		pp[FX] = 0.;
	}
	//conserve EE
	pp[EE] = pp0[EE]*geom0.gdet/geom.gdet;
	// convertiong back to code coordinates happens later	
        #endif
    #endif
	

	// transform pp from BL coordinates to code coordinates 
    trans_pall_coco(pp, pp,BLCOORDS, MYCOORDS, 
            geomBLGC.xxvec,&geomBLGC, &geomGC);

    #ifdef MAGNFIELD
    
        // ensure that gdet B1 is conserved accross r
        // ldouble br0 = pp[B1];
        ldouble br0 = get_u(p, B1, 0, iiy, iiz);
        //brm1 = gdet from cell 0 * br from cell 0 / geom from cell -1
        ldouble brm1 = geom0.gdet*br0/geom.gdet;
        pp[B1] = brm1;

    #endif
