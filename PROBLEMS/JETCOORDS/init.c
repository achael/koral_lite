//int init_dsandvels_limotorus(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell);
int init_dsandvels_fishbone_moncrief(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell);

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

// good faith estimate of initial gamma_gas
#ifdef CONSISTENTGAMMA
ldouble gamma=calc_gammaintfromTei(1.e10,1.e10); 
set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif

// calculate torus rho, uint, ell from tools.c function
init_dsandvels_fishbone_moncrief(r, th, BHSPIN, &rho, &uint, &ell); 
//init_dsandvels_limotorus(r, th, BHSPIN, &rho, &uint, &ell); 
uintorg=uint;

if(rho<0.) //we are outside the donut, in the atmosphere
{
    //atmosphere: atmtype=0. Need to define RHOATMMIN, UINTATMMIN defined at rout=2.
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,ATMTYPE);  //density, temperature
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0); //radiation
#ifdef EVOLVEPHOTONNUMBER
    pp[NF]=1./2.70118/K_BOLTZ * pp[EE] / ATMTRADINIT;   //atmosphere has constant init temp.
#endif
#endif
}
else // we are inside the donut
{
    //calculate hd atmosphere values as backup
    set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,ATMTYPE);
#ifdef RADIATION
    //the initial radiation remains as background
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0); 
#endif

    //Calculate torus 4-velocity from ell; 
    ldouble ult,ulph,ucov[4],ucon[4];

    // old method
    /*
    ell*=-1.;
    ulph = sqrt(-1./(geomBL.GG[0][0]/ell/ell + 2./ell*geomBL.GG[0][3] + geomBL.GG[3][3]));
    ult = ulph / ell;

    ucov[0]=ult;
    ucov[1]=0.;
    ucov[2]=0.;
    ucov[3]=ulph;
    
    indices_12(ucov,ucon,geomBL.GG);
    conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
    */

    // Calculate torus 4-velocity from ell; new method
    // Based on the statement in Kozlowski et al. that ell is actually u_phi u^t
    ldouble rhosq = geomBL.gg[0][0] * geomBL.gg[3][3] - geomBL.gg[0][3] * geomBL.gg[0][3];
    ldouble utsq = (-geomBL.gg[3][3] - sqrt(geomBL.gg[3][3] * geomBL.gg[3][3] - 4. * rhosq * ell * ell)) / (2. * rhosq);
    
    ucon[0] = sqrt(utsq);
    ucon[1] = 0.;
    ucon[2] = 0.;
    ucon[3] = ell / (geomBL.gg[3][3] * ucon[0]) - geomBL.gg[0][3] * ucon[0] / geomBL.gg[3][3];
  
    indices_21(ucon,ucov,geomBL.gg);
    conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);

    // set pp[0], pp[1], but use the atmosphere values as floors
    pp[0]=my_max(rho,ppback[0]); //choice of higher of two densities between background and torus
    pp[1]=my_max(uint,ppback[1]);

    // set initial velocities
    pp[2]=ucon[1]; 
    pp[3]=ucon[2];
    pp[4]=ucon[3];

#ifdef MAGNFIELD
    //initially set B zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.; 
#endif
    
#ifdef RADIATION 

    //make initial radiation atm in torus rotate with the same velocity     
    pp[FX]=pp[VX];
    pp[FY]=pp[VY];
    pp[FZ]=pp[VZ];

#ifdef EVOLVEPHOTONNUMBER
    pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/ATMTRADINIT; //initial rad atm has constant temp
#endif
#endif

    //transforming primitives from BL to MYCOORDS -- take advantage of precomputed MY <--> BL
#ifdef PRECOMPUTE_MY2OUT
    trans_pall_coco_out2my(pp,pp,&geomBL,&geom);
#else      
    trans_pall_coco(pp, pp, BLCOORDS, MYCOORDS, geomBL.xxvec,&geomBL,&geom);
#endif

    //Set vector potential to compute magnetic field
#ifdef MAGNFIELD 

    ldouble Acov[4];
    ldouble r_mag,th_mag;
    Acov[0]=Acov[1]=Acov[2]=0.;

#ifdef INIT_MAGN_CORNERS
    
    // ANDREW define at corners not cell centers
    // ANDREW -- shouldn't this actually be edges? 
    // What about boundary corners? -- we assume we are well inside grid so that A=0 out there!
    
    ldouble xxvec_c[4],xxvecBL_c[4];    
    xxvec_c[0] = global_time;
    xxvec_c[1] = get_xb(ix,0);
    xxvec_c[2] = get_xb(iy,1);
    xxvec_c[3] = get_xb(iz,2);
    coco_N(xxvec_c,xxvecBL_c,MYCOORDS,BLCOORDS);

    ldouble r_c=xxvecBL_c[1];
    ldouble th_c=xxvecBL_c[2];

    r_mag=r_c;
    th_mag=th_c;
#else
    r_mag=r;
    th_mag=th;
#endif
    
#if (NTORUS==0) // single loop v1
    //standard single poloidal loop: Aphi = max[(rho/rhomax) - FM_Aphi_cut, 0]
    //recompute rho -- possibly at corner

    init_dsandvels_fishbone_moncrief(r_mag, th_mag, BHSPIN, &rho, &uint, &ell);
    if(fabs((rho-pp[RHO])/pp[RHO])>.5)
      rho=pp[RHO];

    //vector potential
    Acov[3]=my_max((rho/FM_rho0) - FM_Aphi_cut, 0.);

#elif (NTORUS==1) // single loop v2 -- MAD code comparison
    //standard single poloidal loop: Aphi = max[(rho/rhomax) - FM_Aphi_cut, 0]
    //recompute rho -- possibly at corner

    init_dsandvels_fishbone_moncrief(r_mag, th_mag, BHSPIN, &rho, &uint, &ell);
    if(fabs((rho-pp[RHO])/pp[RHO])>.5)
      rho=pp[RHO];

    ldouble q = (rho/FM_rho0)*pow(r_mag/FM_rin, 3);
    q *= pow(sin(th_mag),3);
    q *= exp(-r_mag / FM_rcut);
    q -= FM_Aphi_cut;
      
    //vector potential
    Acov[3]=my_max(q, 0.);

#elif (NTORUS==2) // dipolar loops
    ldouble lambda = 0.5;
    ldouble anorm = 1.;
    ldouble rchop = 90.;
    ldouble u_av = uintorg;
    ldouble u_av_chop, u_av_mid;

    //midplane at r_mag
    init_dsandvels_fishbone_moncrief(r_mag, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_fishbone_moncrief(rchop, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //vector potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane
  
  
    ldouble rin=FM_rin;
    ldouble STARTFIELD = 1.25*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r_mag > STARTFIELD && r_mag < rchop)
    {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th_mag), 3); 
    }
    else
      q = 0;
    
    if(q > 0.)
    {
      fr = (pow(r_mag,0.6)/0.6  + 0.5/0.4*pow(r_mag,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }
  
    Acov[3]=vpot
    
#elif (NTORUS==3) // quadrupolar loops
    ldouble lambda = 0.5;
    ldouble anorm = 1.;
    ldouble rchop = 90.;
    ldouble u_av = uintorg;
    ldouble u_av_chop, u_av_mid;

    //midplane at r_mag
    init_dsandvels_fishbone_moncrief(r_mag, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_fishbone_moncrief(rchop, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //vector potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane
  
    ldouble rin=FM_rin;
    ldouble STARTFIELD = 1.25*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r_mag > STARTFIELD && r_mag < rchop)
    {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th_mag), 3); 
    }
    else
      q = 0;
   
    if(q > 0.)
    {
      fr = (pow(r_mag,0.6)/0.6  + 0.5/0.4*pow(r_mag,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }
  
    Acov[3]=vpot*sin((M_PI/2.-th_mag));;
#endif

    // Vector potential A is temporarily saved in pp[B1], pp[B2], pp[B3]. These will be replaced by B components immediately
    pp[B1]=Acov[1];
    pp[B2]=Acov[2];
    pp[B3]=Acov[3];
#endif

} //if rho<0
    
// Entropy
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);

// Thermodynamics 
#ifdef EVOLVEELECTRONS

ldouble rhogas=pp[RHO];
ldouble Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);

//slightly colder electrons initially
ldouble ue=1./20.*pp[UU];
ldouble ui=(1.-1./20.)*pp[UU];
ldouble Te,Ti;
pp[ENTRE]=calc_Sefromrhou(calc_thermal_ne(pp)*MU_E*M_PROTON,ue,ELECTRONS);
pp[ENTRI]=calc_Sefromrhou(rhogas,ui,IONS);

#ifdef CONSISTENTGAMMA
Ti=solve_Teifromnmu(pp[RHO]/MU_I/M_PROTON, M_PROTON, ui, IONS); 
Te=solve_Teifromnmu(pp[RHO]/MU_E/M_PROTON, M_ELECTR, ue, ELECTRONS);


gamma=calc_gammaintfromTei(Te,Ti); //recompute gamma_gas from actual electron/ion temps
set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif 

pp[ENTRE]=calc_SefromrhoT(rhogas,Te,ELECTRONS);
pp[ENTRI]=calc_SefromrhoT(rhogas,Ti,IONS);

#endif //EVOLVEELECTRONS

//all primitives to conserved
p2u(pp,uu,&geom);

/***********************************************/
// set the final values
/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
{
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
}

//entropy
update_entropy_cell(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);

/*
if (ix==3 && iy==12)
  {

    printf("Gas: %e  %e %e\n",pp[UU+1],pp[UU+2],pp[UU+3]);
    printf("Rad: %e  %e %e\n",pp[EE0+1],pp[EE0+2],pp[EE0+3]);

      ldouble qsq=0.;
      int ii,jj;
      for(ii=1;ii<4;ii++)
	for(jj=1;jj<4;jj++)
	  qsq+=pp[UU+ii]*pp[UU+jj]*geom.gg[ii][jj];
      ldouble gamma2=1.+qsq;
      printf("Lorentz gamma gas: %e\n\n",sqrt(gamma2));

      qsq=0.;
      for(ii=1;ii<4;ii++)
	for(jj=1;jj<4;jj++)
	  qsq+=pp[EE0+ii]*pp[EE0+jj]*geom.gg[ii][jj];
      gamma2=1.+qsq;
      printf("Lorentz gamma rad: %e\n\n",sqrt(gamma2));

      exit(-1);
  }
*/
