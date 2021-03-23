int init_dsandvels_limotorus(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell);

ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T,uintorg;
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);


#ifdef CONSISTENTGAMMA
ldouble gamma=calc_gammaintfromTei(1.e10,1.e10); //good faith estimate
set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif

ldouble r=geomBL.xx;
ldouble th=geomBL.yy;

init_dsandvels_limotorus(r, th, BHSPIN, &rho, &uint, &ell); //calculate torus tools.c
uintorg=uint;

if(rho<0.) //outside donut
  {
    //ambient
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0); //density, temperature
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0); //radiation
#ifdef EVOLVEPHOTONNUMBER
    //pp[NF0]=calc_NFfromE(pp[EE0]);
    //pp[NF]=calc_NFfromT(calc_PEQ_Tfromurho(pp[UU],pp[RHO]));
    //pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
    pp[NF]=1./2.70118/K_BOLTZ * pp[EE] / ATMTRADINIT;
#endif
#endif
  }
 else //inside donut
   {
    //ambient
    set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0); //background radiation
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
   
    pp[0]=my_max(rho,ppback[0]); //choice of higher of two densities between background and torus
    pp[1]=my_max(uint,ppback[1]);
    pp[2]=ucon[1]; 
    pp[3]=ucon[2];
    pp[4]=ucon[3];

#ifdef RADIATION //make radiation rotate with the same velocity as the torus
    pp[FX]=pp[VX];
    pp[FY]=pp[VY];
    pp[FZ]=pp[VZ];
#endif

    
#ifdef MAGNFIELD//initially set B zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.; 
#endif
    

#ifdef EVOLVEPHOTONNUMBER
    //pp[NF0]=calc_NFfromE(pp[EE0]);
    //pp[NF]=calc_NFfromT(calc_PEQ_Tfromurho(pp[UU],pp[RHO]));
    //pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
    pp[NF]=1./2.70118/K_BOLTZ * pp[EE]/ATMTRADINIT;
#endif

    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);

    //Set Magnetic Field from vector potential
#ifdef MAGNFIELD 
    //MYCOORDS vector potential to calculate B's
    ldouble Acov[4];
    ldouble r_mag,th_mag;
    Acov[0]=Acov[1]=Acov[2]=0.;

#ifdef INIT_MAGN_CORNERS    
    //ANDREW define at corners not cell centers
    //what about boundary corners? -- assume we are well inside grid so that A=0 out there? 
    ldouble xxvec_c[4],xxvecBL_c[4];    
    xxvec_c[0] = global_time;
    xxvec_c[1] = get_xb(ix,0);
    xxvec_c[2] = get_xb(iy,1);
    xxvec_c[3] = get_xb(iz,2);
    coco_N(xxvec_c,xxvecBL_c,MYCOORDS,BLCOORDS);

    ldouble r_c=xxvecBL_c[1];
    ldouble th_c=xxvecBL_c[2];

    //define bfield at corners not cell centers
    r_mag=r_c;
    th_mag=th_c;
#else
    r_mag=r;
    th_mag=th;
#endif
    
#if (NTORUS==4)
    //LIMOFIELD from a=0 SANE harm init.c + denser dipolar loops 
    ldouble lambda = 1.75; //ANDREW changed from 1
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 550.; //outer boundary of field loops
    ldouble u_av = pp[UU];
    ldouble u_av_chop, u_av_mid;
    //midplane at r_mag
    init_dsandvels_limotorus(r_mag, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vetor potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 2.*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r_mag > STARTFIELD && r_mag < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th_mag), 3); // * pow(tanh(r/rsmooth),2);
    } else q = 0;

    if(q > 0.) {
      fr = (pow(r_mag,0.6)/0.6  + 0.5/0.4*pow(r_mag,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }

    Acov[3] = vpot;

#elif (NTORUS==10)
    //LIMOFIELD from a=0 SANE harm init.c + denser quadrupolar loops 
    ldouble lambda = 1.;
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 550.; //outer boundary of field loops
    ldouble u_av = pp[UU];
    ldouble u_av_chop, u_av_mid;
    //midplane at r_mag
    init_dsandvels_limotorus(r_mag, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vector potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 2.*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r_mag > STARTFIELD && r_mag < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th_mag), 3); // * pow(tanh(r_mag/rsmooth),2);
    } else q = 0;

    if(q > 0.) {
      fr = (pow(r_mag,0.6)/0.6  + 0.5/0.4*pow(r_mag,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }

    Acov[3] = vpot;

    //quadrupolar loops
    Acov[3]=vpot*sin((M_PI/2.-th_mag));
    
#elif (NTORUS==14)
    //torus from a=0  harm init.c + dipolar oops
    ldouble lambda = 15.;
    ldouble anorm=1.; //BOBMARK: not used, letting HARM normalize the field
    ldouble rchop = 550.; //outer boundary of field loops
    ldouble u_av = pp[UU];
    ldouble u_av_chop, u_av_mid;
    //midplane at r_mag
    init_dsandvels_limotorus(r_mag, M_PI/2., BHSPIN, &rho, &u_av_mid, &ell);
    //midplane at rchop
    init_dsandvels_limotorus(rchop, M_PI/2., BHSPIN, &rho, &u_av_chop, &ell);
    
    //vector potential follows contours of UU
    ldouble uchop = u_av - u_av_chop; //vpot->zero on contour of radius r=rchop
    ldouble uchopmid = u_av_mid - u_av_chop; //vpot->zero away from midplane

    ldouble rin=LT_RIN;
    ldouble STARTFIELD = 2.*rin;
    ldouble q, fr, fr_start, vpot=0.;
    if (r_mag > STARTFIELD && r_mag < rchop) {
      q = anorm * (uchop - 0.2*uchopmid) / (0.8*uchopmid) * pow(sin(th_mag), 3); // * pow(tanh(r_mag/rsmooth),2);
    } else q = 0;

    if(q > 0.) {
      fr = (pow(r_mag,0.6)/0.6  + 0.5/0.4*pow(r_mag,-0.4)) / lambda;
      fr_start = (pow(STARTFIELD,0.6)/0.6  + 0.5/0.4*pow(STARTFIELD,-0.4)) / lambda;
      vpot += q * sin(fr - fr_start) ;
    }

    Acov[3]=vpot;

#else //standard single poloidal loop
    ldouble dmax = rhoCGS2GU(RHOMAX_VPOT_CGS)*RMAX_VPOT*RMAX_VPOT;
    init_dsandvels_limotorus(r_mag, th_mag, BHSPIN, &rho, &uint, &ell);
    if(fabs((rho-pp[RHO])/pp[RHO])>.5) rho=pp[RHO];
    Acov[3]=my_max(pow(rho*r_mag*r_mag/dmax , 2.) - RHO_VPOT_CUT, 0.) * pow(sin(fabs(th_mag)), 4.);
#endif
    pp[B1]=Acov[1];
    pp[B2]=Acov[2];
    pp[B3]=Acov[3];
#endif

} //if rho<0
    
//entropy
pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);

//thermodynamics and relativistic electrons
#ifdef EVOLVEELECTRONS
#ifdef  RELELECTRONS
int ie;
for (ie=0; ie < NRELBIN; ie++)
{
  if(relel_gammas[ie]<RELEL_INJ_MIN || relel_gammas[ie]>RELEL_INJ_MAX)
  pp[NEREL(ie)] = 0.0 ; //Always initialize with zero nonthermal
  else
  pp[NEREL(ie)] = pow(relel_gammas[ie],-RELEL_HEAT_INDEX);
}

ldouble unth_tmp = calc_relel_uint(pp);
if (unth_tmp > 0)
{
    ldouble scalefac = 1.e-5*pp[UU]/unth_tmp;
    for(ie=0;ie<NRELBIN;ie++)
    {
      pp[NEREL(ie)] *= scalefac;
    }
}
#endif //RELELECtRONS
ldouble rhogas=pp[RHO];
ldouble Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);

//slightly colder electrons initially
ldouble ue=1./100.*pp[UU];
ldouble ui=(1.-1./100.)*pp[UU];
ldouble Te,Ti;
pp[ENTRE]=calc_Sefromrhou(calc_thermal_ne(pp)*MU_E*M_PROTON,ue,ELECTRONS);
pp[ENTRI]=calc_Sefromrhou(rhogas,ui,IONS);

#ifdef CONSISTENTGAMMA
Ti=solve_Teifromnmu(pp[RHO]/MU_I/M_PROTON, M_PROTON, ui, IONS); //solves in parallel for gamma and temperature
Te=solve_Teifromnmu(pp[RHO]/MU_E/M_PROTON, M_ELECTR, ue, ELECTRONS); //solves in parallel for gamma and temperature
    
#ifdef RELELECTRONS
gamma = calc_gammaint_relel(pp, Te, Ti);
#else
gamma=calc_gammaintfromTei(Te,Ti); //good faith estimate
#endif 

set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif //CONSISTENTGAMMA
pp[ENTRE]=calc_SefromrhoT(rhogas,Te,ELECTRONS);
pp[ENTRI]=calc_SefromrhoT(rhogas,Ti,IONS);
#endif //EVOLVEELECTRONS

//to conserved
p2u(pp,uu,&geom);

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
