
	//take from the cell symmetrical wrt. the boundary!
	iix=-ix-1;

    // geometry in the symmetric cell in the domain
    struct geometry geomD;
    fill_geometry(iix,iiy,iiz,&geomD);

    struct geometry geomDBL;
    fill_geometry_arb(iix,iiy,iiz,&geomDBL,BLCOORDS);



	//copy primatives
	for(iv=0;iv<NV;iv++)
	  {
	    pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }
	// transform to BL in cell 0 geometry
	trans_pall_coco(pp, pp, MYCOORDS,BLCOORDS, geomD.xxvec,&geomD,&geomDBL);
	// reflect VX
	pp[VX] = -pp[VX];
	//reflect FX
    #ifdef RADIATION
	    pp[FX] = -pp[FX];
    #endif
    
    #ifdef NONROTATING
    pp[VZ]=0.;
    #endif


	// transform back after reflecting in ghost cell geometry
    trans_pall_coco(pp, pp,BLCOORDS, MYCOORDS, geomBL.xxvec,&geomBL, &geom);

    #ifdef MAGNFIELD
    
        // ensure that gdet B1 is conserved accross r
        // ldouble br0 = pp[B1];
        ldouble br0 = get_u(p, B1, 0, iiy, iiz);
        //brm1 = gdet from cell 0 * br from cell 0 / geom from cell -1
        ldouble brm1 = geom0.gdet*br0/geom.gdet;
        pp[B1] = brm1;

    #endif
