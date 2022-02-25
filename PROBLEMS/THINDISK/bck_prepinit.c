int init_dsandvels_limotorus(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell);

//precalculates the initial state
//done here to be able to refer here even after restart
//B1-B3 in ptemp1 contains only the initial vector potential!
//ENTR positive if inside torus, negative if outside

if(PROCID==0) {fflush(stdout);printf("Precalculating torus and the atmosphere...\n");}
int ix,iy,iz;
#pragma omp parallel for private(ix,iy,iz) schedule (dynamic)
for(iz=0;iz<NZ;iz++)
  {
    for(iy=0;iy<NY;iy++)
      {
	for(ix=0;ix<NX;ix++)
	  {

	    ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
	    ldouble uu[NV], pp[NV],ppback[NV],T,uintorg;
	    ldouble Vphi,Vr;
	    ldouble D,W,eps,uT,uphi,uPhi;

	    //geometries
	    struct geometry geom;
	    fill_geometry(ix,iy,iz,&geom);

	    struct geometry geomBL;
	    fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

	    ldouble r=geomBL.xx;
	    ldouble th=geomBL.yy;

	    init_dsandvels_limotorus(r, th, BHSPIN, &rho, &uint, &ell);
	    uintorg=uint;

	    if(rho<0.) //outside donut
	      {
		//ambient
		set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
		set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif
		pp[ENTR]=-1.;
	      }
	    else //inside donut
	      {
		pp[ENTR]=1.;
		//ambient
		set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
		set_radatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#endif

		uint=LT_KAPPA * pow(rho, LT_GAMMA) / (LT_GAMMA - 1.);
		pgas = GAMMAM1 * uint;
		ell*=-1.;

		ldouble ult,ulph,ucov[4],ucon[4];
		ulph = sqrt(-1./(geomBL.GG[0][0]/ell/ell + 2./ell*geomBL.GG[0][3] + geomBL.GG[3][3]));
		ult = ulph / ell;

		ucov[0]=ult;
		ucov[1]=0.;
		ucov[2]=0.;
		ucov[3]=ulph;
    
		indices_12(ucov,ucon,geomBL.GG);

		conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
   
		pp[0]=my_max(rho,ppback[0]); 
		pp[1]=my_max(uint,ppback[1]);
		pp[2]=ucon[1]; 
		pp[3]=ucon[2];
		pp[4]=ucon[3];

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
		pp[B1]=pp[B2]=pp[B3]=0.; 
#endif


#ifdef RADIATION
		//distributing pressure
		ldouble P,aaa,bbb;
		P=GAMMAM1*uint;
		//solving for T satisfying P=pgas+prad=bbb T + aaa T^4
		aaa=4.*SIGMA_RAD/3.;
		bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
		ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));
		ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;

		E=calc_LTE_EfromT(T4);
		Fx=Fy=Fz=0.;
		uint=calc_PEQ_ufromTrho(T4,rho);

		pp[UU]=my_max(uint,ppback[1]);
		pp[EE0]=my_max(E,ppback[EE0]);



		pp[FX0]=Fx;
		pp[FY0]=Fy;
		pp[FZ0]=Fz;

		//transforming from BL lab radiative primitives to code non-ortonormal primitives
		prad_ff2lab(pp,pp,&geomBL);

#endif

#ifdef NCOMPTONIZATION
		pp[NF0]=calc_NFfromE(pp[EE0]);
#endif

		//transforming primitives from BL to MYCOORDS
		trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
    
#ifdef MAGNFIELD 
		//MYCOORDS vector potential to calculate B's
		ldouble Acov[4];
		Acov[0]=Acov[1]=Acov[2]=0.;

#if(NTORUS==3)
		//LIMOFIELD from a=0 SANE harm init.c
		ldouble lambda = 2.5;
		ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
		ldouble rchop = 550.; //outer boundary of field loops
		ldouble u_av = pp[UU];
		ldouble u_av_chop, u_av_mid;
		//midplane at r
		init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
		//midplane at rchop
		init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
		//vetor potential follows contours of UU
		ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
		ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

		ldouble rin=LT_RIN;
		ldouble STARTFIELD = 2.5*rin;
		ldouble q, fr, fr_start, vpot=0.;
		if (r > STARTFIELD && r < rchop) {
		  q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
		} else q = 0;

		if(q > 0.) {
		  fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
		  fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
		  vpot += q * sin(fr - fr_start) ;
		}

		Acov[3]=vpot;

#elif (NTORUS==4)
		//LIMOFIELD from a=0 SANE harm init.c + denser loops
		ldouble lambda = 1.;
		ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
		ldouble rchop = 550.; //outer boundary of field loops
		ldouble u_av = pp[UU];
		ldouble u_av_chop, u_av_mid;
		//midplane at r
		init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
		//midplane at rchop
		init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
		//vetor potential follows contours of UU
		ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
		ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

		ldouble rin=LT_RIN;
		ldouble STARTFIELD = 2.5*rin;
		ldouble q, fr, fr_start, vpot=0.;
		if (r > STARTFIELD && r < rchop) {
		  q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
		} else q = 0;

		if(q > 0.) {
		  fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
		  fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
		  vpot += q * sin(fr - fr_start) ;
		}

		Acov[3]=vpot;

#elif (NTORUS==5)

		//LIMOFIELD from a=0 SANE harm init.c for mimic_dynamo
		ldouble lambda = 2.5;
		ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
		ldouble rchop = 550.; //outer boundary of field loops
		ldouble u_av = pp[UU];
		ldouble u_av_chop, u_av_mid;
		//midplane at r
		init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
		//midplane at rchop
		init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
		//vetor potential follows contours of UU
		ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
		ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

		ldouble rin=LT_RIN;
		ldouble STARTFIELD = 2.5*rin;
		ldouble q, fr, fr_start, vpot=0.;
		if (r > STARTFIELD && r < rchop) {
		  q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
		} else q = 0;

    
		if(q > 0.) {
		  fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
		  fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
		  vpot += q * sin(fr - fr_start) ;
		}
        
		Acov[2]=vpot*sin((M_PI/2.-geomBL.yy));;

    
#elif (NTORUS==7 )

		//LIMOFIELD from a=0 SANE harm init.c with flipping polarity in theta
		ldouble lambda = 2.5;
		ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
		ldouble rchop = 350.; //outer boundary of field loops
		ldouble u_av = uintorg;
		ldouble u_av_chop, u_av_mid;
		//midplane at r
		init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
		//midplane at rchop
		init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
		//vetor potential follows contours of UU
		ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
		ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

		ldouble rin=LT_RIN;
		ldouble STARTFIELD = 1.25*rin;
		ldouble q, fr, fr_start, vpot=0.;
		if (r > STARTFIELD && r < rchop) {
		  q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
		} else q = 0;

    
		if(q > 0.) {
		  fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
		  fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
		  vpot += q * sin(fr - fr_start) ;
		}
     
		//    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
		Acov[3]=vpot*sin((M_PI/2.-geomBL.yy));



#elif (NTORUS==8)

		Acov[3]=my_max(pow(pp[RHO]*geomBL.xx/4.e-22,2.)-0.02,0.)*sqrt(1.e-23)*sin((M_PI/2.-geomBL.yy)/0.1);


#elif (NTORUS==6)

		//LIMOFIELD from a=0 SANE harm init.c with flipping polarity in theta
		ldouble lambda = 2.5;
		ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
		ldouble rchop = 350.; //outer boundary of field loops
		ldouble u_av = uintorg;
		ldouble u_av_chop, u_av_mid;
		//midplane at r
		init_dsandvels_limotorus(r, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
		//midplane at rchop
		init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
		//vetor potential follows contours of UU
		ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
		ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

		ldouble rin=LT_RIN;
		ldouble STARTFIELD = 1.25*rin;
		ldouble q, fr, fr_start, vpot=0.;
		if (r > STARTFIELD && r < rchop) {
		  q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th), 3); // * pow(tanh(r/rsmooth),2);
		} else q = 0;

    
		if(q > 0.) {
		  fr = (pow(r,0.6)/0.6  + 0.5/0.4*pow(r,-0.4)) / lambda;
		  fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
		  vpot += q * sin(fr - fr_start) ;
		}
     
		//    if(iy==NY/2) printf("%d %f %f > %e %e %e %e\n",iy,r,th,uchop,u_av_mid,u_av, u_av_chop);
		Acov[3]=vpot*sin((M_PI/2.-geomBL.yy));;

#else //standard single poloidal loop
		Acov[3]=my_max(pow(pp[RHO]*geomBL.xx*geomBL.xx/4.e-20,2.)-0.02,0.)*sqrt(1.e-23)*pow(sin(fabs(geomBL.yy)),4.);
#endif

		pp[B1]=Acov[1];
		pp[B2]=Acov[2];
		pp[B3]=Acov[3];
#endif

	      }

	    
	    //to conserved
	    p2u(pp,uu,&geom);



	    /***********************************************/

	    int iv;

	    for(iv=0;iv<NV;iv++)
	      {
		set_u(ptemp1,iv,ix,iy,iz,pp[iv]);
	      }

	  }
      }
  }



