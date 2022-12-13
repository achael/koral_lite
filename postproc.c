/*! \file postproc.c
 \brief Routines for postprocessing, used both on the go and separately
 */

#include "ko.h"

/*********************************************/
//int calc_radialprofiles(ldouble profiles[][NX])
/* calculates radial profiles - L(r) etc. */
/* uses mostly primitives, but may use avg quantities as well */
/* however, this part would work only when postprocessing with */
/* ./avg, otherwise pavg hold non-normalized numbers */

//surface density (2) (column)
//rest mass flux (3)
//rho-weighted minus radial velocity (4)
//<u_phi>* (5)
//Keplerian u_phi (6)
//abs optical depth (7)
//tot optical depth (8)
//net accretion rate at given radius (9)
//inflow accretion rate at given radius (10)
//outflow accretion rate at given radius (11)
//luminosity at given radius opt. thin (12)
//location of the photosphere (13)
//total mhd energy flux (14)
//outflowin mhd energy flux (15)
//jet mhd energy flux (16)
//total rad energy flux (17)
//outflowin rad energy flux (18)
//jet rad energy flux (19)
//outflowin mass flux (20)
//jet mass flux (21)
//luminosity at given radius everywhere (22)
//surface density in the inflow (23)
//rho-weighted minus radial velocity in the inflow (24)
//opt thin mhd energy flux (25)
//opt. thin rad energy flux (26)
//opt. thin mass flux (27)
//rho-weighted qtheta (28)
//rho-weighted temperature (29)
//rho-weighted magn.field angle <sqrt(grr gphph)b^r b^ph> / <bsq> (30)
//scale-height (31)
//rho-weighted beta (32)
//rho-wighted prad/pgas (33)
//alpha (34)
//rad. viscosity energy flux (35)
//rho-weighted minus radial velocity in the outflow (36) -> rho-weighted eccentricity for problem 87

//conserved flux rho ur transformed to OUTCOORDS (37)
//conserved flux rho ur in MYCOORDS (38)
//conserved flux rho ur+Trt in MYCOORDS (39)
//conserved flux for Rrt int MYCOORDS (40)

//surface density of energy for Be<0. = int (Ehat + uint) (41)
//rho-weighted radial velocity in the jet (42)
//magnetic flux in the jet (43)
//kinetic + binding flux in the jet (44)
//radial velocity close to the axis (45)                                        
//rho-weighted Bernoulli (46)  
//rho-weighted qphi (47)          
//magnetic flux everywhere (48)
//kinetic flux everywhere (49)
//radial profiles of rho-weighted radiation temperature (50)
//<u_phi> in the outflow (51)
//integrated heating/cooling rate fluid frame G^t (52)
//energy-averaged radial transport velocity (53)
//energy-averaged vertical transport velocity (54)
//inclination angle of the integrated angular momentum vector (55)
//thermal energy flux everywhere (56)
//convective hydro luminosity (57)
//integrated T^r_phi (58)
//rho-weighted Bernoulli without thermal component (59)  
//average dissipation (60)
//jet CM X (61) (see code comparison document)
//jet CM Y (62)
//jet Ixx (63)
//jet Iyy (64)
//jet Ixy (65)
//jet Area (66)
/*********************************************/

int calc_radialprofiles(ldouble profiles[][NX])
{
  //adjust NRADPROFILES in problem.h
  #ifdef CONSISTENTGAMMA
  if(PROCID==0) printf("1 gammas not consistent in postprocessing when CONSISTENTGAMMA\n");
  #endif

  //calculates scale-height etc.
  calc_avgs_throughout();     
  int ix;
  
  //define extra quantities for problem 134
  ldouble normalize[NX], rho134[NX], pgas134[NX], bfield134[NX], uconphi134[NX], ptot134[NX], betainv134[NX],scaleheight134[NX],rholambda134[NX];
  ldouble thetamin_134=M_PI/3.;
  ldouble thetamax_134=2*M_PI/3.;
  ldouble sigmagdet[NX];
   
  //search for appropriate radial index
#pragma omp parallel for private(ix)
  for(ix=0;ix<NX;ix++)
  {
    int iy,iz,iv,i,j;
    ldouble x0[3],x0l[3],x0r[3],xm1[3],xp1[3];
    ldouble dx0, dxm2, dxm1, dxp1, dxp2;
    ldouble xx[4],xxBL[4],dx[3],dxph[3],dxcgs[3],mdot,rho,rhouconrcons,uint,temp,temprad,ucon[4],utcon[4],ucon3[4],TttBe,TttBenoth;
    ldouble rhoucont,enden,rhouconr,rhouconrucovphi,Tij[4][4],Tij22[4][4],Rij[4][4],Rviscij[4][4],Trt,Trphi,Fluxx[NV],Rrt,Rtt,Ttt,Rviscrt,bsq,bcon[4],bcov[4];
    ldouble Trtmagn,Trtkin,Trtuint,endensimple,endenuconr,endenuconth,Omega;
    ldouble ucov[4],pp[NV],gg[4][5],GG[4][5],ggBL[4][5],GGBL[4][5],Ehat;
    ldouble tautot,tautotloc,tauabs,tauabsloc;
    ldouble avgsums[NV+NAVGVARS][NX];
    ldouble Bangle1,Bangle2,brbphi;
    ldouble diffflux[NV];
    ldouble fd_u0[NV],fd_up1[NV],fd_up2[NV],fd_um1[NV],fd_um2[NV];
    ldouble fd_p0[NV],fd_pp1[NV],fd_pp2[NV],fd_pm1[NV],fd_pm2[NV],fd_pm3[NV],fd_pp3[NV];
    ldouble fd_pl[NV],fd_pr[NV],fd_plm1[NV],fd_prm1[NV],fd_plp1[NV],fd_prp1[NV];
    ldouble fd_ul[NV],fd_ur[NV],fd_ulm1[NV],fd_urm1[NV],fd_ulp1[NV],fd_urp1[NV];
    ldouble du[NV],dul[NV],dur[NV],aaa[24],ahd,arad,Jvec[3];
    ldouble Gi[4],Giff[4];
    
    //vertically integrated/averaged profiles
    
    for(iv=0;iv<NRADPROFILES;iv++)
      profiles[iv][ix]=0.;
    
    Jvec[0]=Jvec[1]=Jvec[2]=0.;
    
    //keep track of extra quantities for problem 134/139
    if (PROBLEM == 134 || PROBLEM == 139)
    {
      normalize[ix] = 0.;
      rho134[ix] = 0.;
      pgas134[ix] = 0.;
      bfield134[ix] = 0.;
      uconphi134[ix] = 0.;
      ptot134[ix] = 0.;
      betainv134[ix] = 0.;
      scaleheight134[ix]=0.;
      rholambda134[ix]=0.;
    }
    
    //outside horizon?
    struct geometry geomBLtemp;
    fill_geometry_arb(ix,0,0,&geomBLtemp,OUTCOORDS);
    if(geomBLtemp.xx<=1.1*rhorizonBL) continue; //to avoid working inside horizon
    
    Bangle1=Bangle2=0.;
    sigmagdet[ix]=0.;
    ldouble jetsigma=0.;
    
    for(iv=0;iv<NAVGVARS;iv++)
      avgsums[iv][ix]=0.;
    
    tautot=tauabs=0.;
    
    for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
      {
        //metric
        pick_g(ix,iy,iz,gg);
        pick_G(ix,iy,iz,GG);

	// cell center geometries
        struct geometry geom;
        fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);
        
        struct geometry geomm1;
        fill_geometry_arb(ix-1,iy,iz,&geomm1,MYCOORDS);
        
        struct geometry geomp1;
        fill_geometry_arb(ix+1,iy,iz,&geomp1,MYCOORDS);
        
        struct geometry geomBL;
        fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
        
        struct geometry geomBLm1;
        fill_geometry_arb(ix-1,iy,iz,&geomBLm1,OUTCOORDS);
        
        struct geometry geomBLp1;
        fill_geometry_arb(ix+1,iy,iz,&geomBLp1,OUTCOORDS);

	// face geometries
        struct geometry geoml;
        fill_geometry_face(ix,iy,iz,0,&geoml);
        
        struct geometry geomr;
        fill_geometry_face(ix+1,iy,iz,0,&geomr);
        
        struct geometry geomBLl;
        fill_geometry_face_arb(ix,iy,iz,0,&geomBLl,OUTCOORDS);
        
        struct geometry geomBLr;
        fill_geometry_face_arb(ix+1,iy,iz,0,&geomBLr,OUTCOORDS);
        
        ldouble gdetuBL=geomBL.gdet;
        
#ifdef RADOUTPUTWITHINDISK
        if(fabs(geomBL.yy-M_PI/2)>scaleth_otg[ix])
          continue;
#endif
        
        
#ifdef RADOUTPUTWITHINDTHETA  // this limits the theta integral
	if(fabs(geomBL.yy-M_PI/2)>RADOUTPUTWITHINDTHETA)
          continue;
#endif
        
#if (GDETIN==0) //gdet out of derivatives
        gdetuBL=1.;
#endif
        
        //coordinates
	get_xx(ix,iy,iz,xx);
        #ifdef PRECOMPUTE_MY2OUT
        get_xxout(ix, iy, iz, xxBL); 
        #else
        coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
        #endif

	//cell dimensions
	//ANDREW put cell size code in a function with precompute option
        get_cellsize_out(ix, iy, iz, dx);

	/*
	ldouble xx1[4],xx2[4];
        xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_x(iz,2);
        xx2[0]=0.;xx2[1]=get_xb(ix+1,0);xx2[2]=get_x(iy,1);xx2[3]=get_x(iz,2);
        coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
        coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
        dx[0]=fabs(xx2[1]-xx1[1]);
        xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_x(iz,2);
        xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_xb(iy+1,1);xx2[3]=get_x(iz,2);
        coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
        coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
        dx[1]=fabs(xx2[2]-xx1[2]);
        xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_xb(iz,2);
        xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_x(iy,1);xx2[3]=get_xb(iz+1,2);
        coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
        coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
        dx[2]=fabs(xx2[3]-xx1[3]);
        */
	
	if(NZ==1) 
        {
          dx[2]=2.*M_PI;
        }
        else
        {
          #ifdef PHIWEDGE
          dx[2]*=(2.*M_PI/PHIWEDGE);
          #endif
        }

        ldouble dxph[3];
        dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
        dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
        dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);
        
        //primitives at the cell - either averaged or original, in BL or MYCOORDS
        for(iv=0;iv<NV;iv++)
        {
          pp[iv]=get_u(p,iv,ix,iy,iz);
        }


	#ifdef CONSERVED_FLUX_OUTPUT
        // ANDREW TODO
	// These are all defined on the face, so we can't use precomputed
        // can take a long time and don't seem necessary in radial profiles

        //primitives at radial neighbours
        for(i=0;i<NV;i++)
        {
          fd_p0[i]=pp[i];
          fd_pp1[i]=get_u(p,i,ix+1,iy,iz);
          fd_pm1[i]=get_u(p,i,ix-1,iy,iz);
          fd_pm2[i]=get_u(p,i,ix-2,iy,iz);
          fd_pp2[i]=get_u(p,i,ix+2,iy,iz);
        }
        
        //internal coordinates
        x0[0]=get_x(ix,0);
        
        x0l[0]=get_xb(ix,0);
        xm1[0]=get_x(ix-1,0);
        x0l[1]=xm1[1]=get_x(iy,1);
        x0l[2]=xm1[2]=get_x(iz,2);
        
        x0r[0]=get_xb(ix+1,0);
        xp1[0]=get_x(ix+1,0);
        x0r[1]=xp1[1]=get_x(iy,1);
        x0r[2]=xp1[2]=get_x(iz,2);
        
        dx0=get_size_x(ix,0);
        dxm1=get_size_x(ix-1,0);
        dxp1=get_size_x(ix+1,0);
        
        //interpolation to get left/right biased interpolated primtives
        avg2point(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2,0,MINMOD_THETA);
        avg2point(fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pp2,fd_plp1,fd_prp1,dxm1,dx0,dxp1,dxp2,dxp2,0,MINMOD_THETA);
        avg2point(fd_pm2,fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_plm1,fd_prm1,dxm2,dxm2,dxm1,dx0,dxp1,0,MINMOD_THETA);

        trans_pall_coco(fd_pl,fd_pl,MYCOORDS,OUTCOORDS,xx,&geoml,&geomBLl);
        trans_pall_coco(fd_pr,fd_pr,MYCOORDS,OUTCOORDS,xx,&geomr,&geomBLr);
        trans_pall_coco(fd_plp1,fd_plp1,MYCOORDS,OUTCOORDS,xx,&geoml,&geomBLl);
        trans_pall_coco(fd_prm1,fd_prm1,MYCOORDS,OUTCOORDS,xx,&geomr,&geomBLr);
        #else
	for(i=0;i<NV;i++) Fluxx[i]=0.; // no conserved fluxes if this flag is off
	#endif
        
	//transforming interpolated primitives to BL
        #ifdef PRECOMPUTE_MY2OUT
        trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
	//trans_pall_coco_my2out(fd_pm1,fd_pm1,&geomm1,&geomBLm1); // don't need these anymore
        //trans_pall_coco_my2out(fd_pp1,fd_pp1,&geomp1,&geomBLp1);
        #else      
        trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,xx,&geom,&geomBL);
	//trans_pall_coco(fd_pm1,fd_pm1,MYCOORDS,OUTCOORDS,xx,&geomm1,&geomBLm1); // don't need these anymore
        //trans_pall_coco(fd_pp1,fd_pp1,MYCOORDS,OUTCOORDS,xx,&geomp1,&geomBLp1);
        #endif
	
        ldouble vischeating=0.;
        
        if(doingavg)
        {
          rho=get_uavg(pavg,RHO,ix,iy,iz);
          uint=get_uavg(pavg,UU,ix,iy,iz);
          temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
          bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
          bcon[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
          bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
          bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
          bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
          utcon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
          utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
          utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
          utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
          vischeating=get_uavg(pavg,AVGVISCHEATING,ix,iy,iz);
          
          ldouble uconpp[4]={utcon[0],utcon[1],utcon[2],utcon[3]};
          conv_vels(uconpp,uconpp,VEL4,VEL4,geomBL.gg,geomBL.GG);
          conv_vels(uconpp,uconpp,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
          pp[VX]=uconpp[1];
          pp[VY]=uconpp[2];
          pp[VZ]=uconpp[3];
          
          Omega=utcon[3]/utcon[0];
          rhouconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
          rhoucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz);
          
          endenuconr = get_uavg(pavg,AVGUUUCON(1),ix,iy,iz) * sqrt(geomBL.gg[1][1]);
          endenuconth = get_uavg(pavg,AVGUUUCON(2),ix,iy,iz) * sqrt(geomBL.gg[2][2]);
          
          rhouconrcons=get_uavg(pavg,AVGRHOURDIFF,ix,iy,iz);
          rhouconrucovphi=get_uavg(pavg,AVGRHOUCONUCOV(1,3),ix,iy,iz);
          
          TttBe=get_uavg(pavg,AVGRHOUCONUCOV(0,0),ix,iy,iz)
          + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(0,0),ix,iy,iz)
          + get_uavg(pavg,AVGBSQUCONUCOV(0,0),ix,iy,iz)
          - get_uavg(pavg,AVGBCONBCOV(0,0),ix,iy,iz);
          
          TttBenoth=get_uavg(pavg,AVGRHOUCONUCOV(0,0),ix,iy,iz);
          
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)              
              Tij[i][j]=get_uavg(pavg,AVGTIJ(i,j),ix,iy,iz);
          
          for(i=0;i<NV;i++)
            Fluxx[i]=get_uavg(pavg,AVGFLUXXL(i),ix,iy,iz);
          
          Trphi=Tij[1][3];
          Trt=Tij[1][0];
          Ttt=Tij[0][0];
          enden = Tij[0][0]+rhoucont;
          endensimple = uint;
          
          Trtmagn = get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
          - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz);
          
          Trtkin = get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz);
          
          Trtuint =  get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz);
          
          indices_2122(Tij,Tij22,geomBL.GG);
          
#ifdef RADIATION
          for(i=0;i<4;i++)
            for(j=0;j<4;j++)
              Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
          
          enden += Rij[0][0];
          Rrt = Rij[1][0];
          Rtt = Rij[0][0];
          Ehat = get_uavg(pavg,AVGEHAT,ix,iy,iz);
          temprad=calc_ncompt_Thatrad_fromEN(Ehat,get_uavg(pavg,AVGNFHAT,ix,iy,iz));
          
          endensimple+=Ehat;
          
          endenuconr += get_uavg(pavg,AVGEHATUCON(1),ix,iy,iz) * sqrt(geomBL.gg[1][1]);
          endenuconth += get_uavg(pavg,AVGEHATUCON(2),ix,iy,iz) * sqrt(geomBL.gg[2][2]);
          
          
          int derdir[3]={0,0,0};
          calc_Rij_visc(pp,&geomBL,Rviscij,derdir);
          
          Rviscrt = Rviscij[1][0];
          
          for(iv=0;iv<4;iv++)
            Giff[iv]=get_uavg(pavg,AVGGHAT(iv),ix,iy,iz);
#endif
          
          //no need of transforming interpolated primitives to BL, already there
          
        }  // end if(doingavg)
        else //on the go from the primitives
        {
          rho=pp[0];
          uint=pp[1];
          utcon[1]=pp[2];
          utcon[2]=pp[3];
          utcon[3]=pp[4];
          temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
          
#ifdef MAGNFIELD
          calc_bcon_prim(pp,bcon,&geomBL);
          indices_21(bcon,bcov,geomBL.gg);
          bsq = dotB(bcon,bcov);
#endif
          
          conv_vels_both(utcon,utcon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
          Omega=utcon[3]/utcon[0];
          rhouconr=rhouconrcons=rho*utcon[1];
          
          endenuconr = uint*utcon[1]* sqrt(geomBL.gg[1][1]);
          endenuconth = uint*utcon[2]* sqrt(geomBL.gg[2][2]);
          
          calc_Tij(pp,&geomBL,Tij22);
          indices_2221(Tij22,Tij,geomBL.gg);
          
          Trphi=Tij[1][3];
          Trt = Tij[1][0];
          Ttt = Tij[0][0];
          
          TttBe=Ttt-(GAMMAM1*uint+bsq/2.);
          TttBenoth= rho*utcon[0]*ucov[0];
          
          Trtmagn = bsq*utcon[1]*ucov[0] - bcon[1]*bcov[0];
          Trtkin = rho*utcon[1]*ucov[0];
          
          Trtuint = uint*utcon[1]*ucov[0];
          enden = Tij[0][0] + rho*utcon[0];
          endensimple = uint;
          
#ifdef RADIATION
          calc_Rij(pp,&geomBL,Rij);
          indices_2221(Rij,Rij,geomBL.gg);
          
          Rrt = Rij[1][0];
          Rtt = Rij[0][0];
          
          ldouble Rtthat,uconr[4];
          calc_ff_Rtt(&get_u(p,0,ix,iy,iz),&Rtthat,uconr,&geomBL);
          Ehat=-Rtthat;
          enden+=Rij[0][0];
          endensimple+=Ehat;
          endenuconr+=Ehat*ucon[1]*sqrt(geomBL.gg[1][1]);
          endenuconth+=Ehat*ucon[2]*sqrt(geomBL.gg[2][2]);
          
          int derdir[3]={0,0,0};
          calc_Rij_visc(pp,&geomBL,Rviscij,derdir);
          
          Rviscrt = Rviscij[1][0];
          
          //lab-frame four fource
          calc_Gi(pp,&geomBL,Gi,0.0, 1, 0);//ANDREW Lab frame=1
          boost2_lab2ff(Gi,Giff,pp,geomBL.gg,geomBL.GG);
#endif

	  #ifdef CONSERVED_FLUX_OUTPUT
          //estimating the diffusive flux
          //to conserved
          p2u(fd_pl,fd_ul,&geomBLl);
          p2u(fd_pr,fd_ur,&geomBLr);
          p2u(fd_plp1,fd_ulp1,&geomBLr);
          p2u(fd_prm1,fd_urm1,&geomBLl);
          
          //gradient of conserved
          PLOOP(iv)
          {
            //getting rid of gdetu in conserved - integrated with gdet lateron
            dul[iv]=fd_ul[iv]-fd_urm1[iv];
            dur[iv]=fd_ulp1[iv]-fd_ur[iv];
            du[iv]=.5*(dul[iv]+dur[iv]);
            du[iv]/=gdetuBL;
            
          }
          
          //wavespeeds
          calc_wavespeeds_lr_pure(pp,&geomBL,aaa);
          ahd=my_max(fabs(aaa[0]),fabs(aaa[1]));
          arad=my_max(fabs(aaa[6]),fabs(aaa[7]));
          
          //diffusive flux
          PLOOP(iv)
          {
            if(iv<NVMHD) diffflux[iv]=-0.5*ahd*du[iv];
            else diffflux[iv]=-0.5*arad*du[iv];
          }
          
          //adding up to the conserved fluxes
          Fluxx[RHO]=geomBL.gdet*(rhouconr+diffflux[RHO]);
          Fluxx[UU]=geomBL.gdet*(rhouconr+Trt+diffflux[UU]);
          
          #ifdef RADIATION
          Fluxx[EE0]=geomBL.gdet*(Rrt+diffflux[EE0]);
          #endif
          #endif //CONSERVED_FLUX_OUTPUT   
        }  // end on the go from primitives
        
        ldouble muBe,Be,Benoth,betagamma2;
	betagamma2=(-Trt/rhouconr)*(-Trt/rhouconr) - 1.;
        muBe=-(Trt+rhouconr)/rhouconr;
        Be=-(TttBe+rhoucont)/rhoucont;
        Benoth=-(TttBenoth+rhoucont)/rhoucont;
#ifdef RADIATION
        muBe+=-Rrt/rhouconr;
        Be+=-Rtt/rhoucont;
#endif
        
        
        int isjet;
	//old criterion based on Bernoulli
        //if(muBe>0.05 && (xxBL[2]<M_PI/4. || xxBL[2]>3.*M_PI/4.))
        //  isjet=1;
	//new criterion based on beta*gamma
	if(betagamma2>1. && xxBL[2]<M_PI/3.) // top jet only
	  isjet=1;
        else 
          isjet=0;
        
        ldouble pregas = GAMMAM1*uint;
        ldouble ptot = pregas;
#ifdef MAGNFIELD
        ptot+=bsq/2.;
#endif
#ifdef RADIATION
        ldouble prerad = Ehat/3.;
        ptot+=prerad;
#endif
        
        //alpha
        boost22_lab2ff(Tij22,Tij22,pp,geomBL.gg,geomBL.GG);
        ldouble alpha=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij22[1][3]/ptot;
        
        //angular velocity
        ldouble Omega=utcon[3]/utcon[0];
        
        //MRI resolution parameters
        ldouble qtheta,qphi;
        calc_Qthetaphi(ix,iy,iz,&qtheta,&qphi);
        
        //to calculate magn. field angle
        //bsq and brbphi taken from avg if neeeded
        ldouble bconfake[4],bcovfake[4],bsqfake;
        calc_angle_brbphibsq(ix,iy,iz,&brbphi,&bsqfake,bconfake,bcovfake);
        Bangle1+=rho*brbphi*dxph[1];
        Bangle2+=rho*bsq*dxph[1];
        
        //optical depths
        struct opacities opac;
        ldouble tauabsloc = utcon[0]*calc_kappa(pp,&geomBL,&opac);
        ldouble tautotloc = utcon[0]*calc_kappaes(pp,&geomBL);
        tautot+=tautotloc*dxph[1];
        tauabs+=tauabsloc*dxph[1];
        
        //alpha (34) (column)
        profiles[32][ix]+=alpha*rho*dxph[1];
        
	//rho-weighted beta (32)
        ldouble prermhd = GAMMAM1*uint;
#ifdef RADIATION
        prermhd+=Ehat/3.;
#endif
        ldouble ibeta=bsq/2./(prermhd+bsq/2.);
        profiles[30][ix]+=rho*ibeta*dxph[1];
        if (PROBLEM==93){
          profiles[57][ix]+=bsq/2.*dxph[1];
          profiles[58][ix]+=prermhd*dxph[1];
        }
        
        //rho-weighted prad/pgas (33)
#ifdef RADIATION
        profiles[31][ix]+=rho*prerad/(pregas+prerad+0.5*bsq)*dxph[1];
#else
        profiles[31][ix]+=0.;
#endif
        
        //surface density (2) (column)
        profiles[0][ix]+=rho*dxph[1];
        
        //surface energy density (41)
        if(muBe<0.)
          profiles[39][ix]+=endensimple*dxph[1];
        
        //numerator of scale height (31) (column)
#ifndef CALCHRONTHEGO
        profiles[29][ix]+=rho*dxph[1]*pow(tan(fabs(M_PI/2.-xxBL[2])),2.);
#endif
        
        //surface density in the inflow (23)
        //if(utcon[1]<0.)
        //  profiles[21][ix]+=rho*dxph[1];
        
        //rho-weighted q-theta (28)
        profiles[26][ix]+=rho*qtheta*dxph[1];
	       
        //rho-weighted q-phi (47)
        profiles[45][ix]+=rho*qphi*dxph[1];
	       
        //rho-weighted temperature (29)
        profiles[27][ix]+=rho*temp*dxph[1];
        
        //rho-weighted rad temperature (50)
        profiles[48][ix]+=rho*temprad*dxph[1];
        
        //rest mass flux (3)
        profiles[1][ix]+=-rhouconr*dx[1]*dx[2]*geomBL.gdet;
        
        //conserved flux (rhour) transformed to OUTCOORDS (may be imprecise) (37)
        profiles[35][ix]+=-rhouconrcons*dx[1]*dx[2];
        
        //conserved flux (gdet rhour) in MYCOORDS (38)
        profiles[36][ix]+=-Fluxx[RHO]*get_size_x(iy,1)*dx[2];
        
        //conserved flux for Trt (39) in MYCOORDS
        profiles[37][ix]+=-Fluxx[UU]*get_size_x(iy,1)*dx[2];
        
        //temporary surface density to normalize scale height
	sigmagdet[ix] +=rho*dx[1]*dx[2]*geomBL.gdet;

	//scale height using code comparison paper defn;
	scaleheight134[ix] += rho*fabs(0.5*M_PI - xxBL[2]) * dx[1]*dx[2]*geomBL.gdet;

	//density-weighted MRI wavelength using code comparison paper defn
	//ANDREW TODO -- should it be bsq/2?
	ldouble lambdaMRI = 2*M_PI * fabs(bcon[2]) / sqrt(bsq + rho + uint + pregas) / fabs(Omega);
	rholambda134[ix] += rho*lambdaMRI* dx[1]*dx[2]*geomBL.gdet; 
	  
        //total mhd energy flux (14)
        profiles[12][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;
        
        //convective mhd energy flux (57)
        profiles[55][ix]+=rhouconr/rho*(Ttt)*dx[1]*dx[2]*geomBL.gdet;
        
        //integrated Trphi (58)
        profiles[56][ix]+=(Trphi-rhouconrucovphi)*Omega*dx[1]*dx[2]*geomBL.gdet;
        
        //magnetic mhd energy flux in jet (43)
        if(isjet==1)
          profiles[41][ix]+=(-Trtmagn)*dx[1]*dx[2]*geomBL.gdet;
        
        //kinetic + binding mhd energy flux in jet (44)
        if(isjet==1)
          profiles[42][ix]+=(-Trtkin)*dx[1]*dx[2]*geomBL.gdet;
        
        //magnetic mhd energy flux (48)
        profiles[46][ix]+=(-Trtmagn)*dx[1]*dx[2]*geomBL.gdet;
        
        //kinetic  mhd energy flux (49)
        profiles[47][ix]+=(-Trtkin)*dx[1]*dx[2]*geomBL.gdet;
        
        //thermal energy flux (56)
        profiles[54][ix]+=(-Trtuint)*dx[1]*dx[2]*geomBL.gdet;
        
        //integrated cooling rate (52)
        profiles[50][ix]+=Giff[0]*dx[1]*dx[2]*geomBL.gdet;
        
        //opt thin mhd energy flux (25)
        if(tautot<1.)
          profiles[23][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;
        
        //outflowin mhd energy flux (15)
        if(utcon[1]>0.)
          profiles[13][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;
        
        //jet mhd energy flux (16)
        if(isjet==1)
          profiles[14][ix]+=(-Trt)*dx[1]*dx[2]*geomBL.gdet;
        
        //jet gas velocity (42)
        if(isjet==1)
        {
          profiles[40][ix]+=rhouconr*dx[1]*dx[2]*geomBL.gdet;
          jetsigma+=rho*dx[1]*dx[2]*geomBL.gdet;
        }
        
        //gas velocity near the axis (45)
        if((doingavg && (iy==NCCORRECTPOLAR+1 || iy==(NY-NCCORRECTPOLAR-2))))
          profiles[43][ix]+=0.5*utcon[1];
        if((!doingavg &&  (iy==NCCORRECTPOLAR+1)))
          profiles[43][ix]+=utcon[1];
        
        //rho-weighted Bernoulli (46)
        profiles[44][ix]+=rho*Be*dxph[1];
        
        //rho-weighted Bernoulli without thermal component (59)
	profiles[57][ix]+=rho*Benoth*dxph[1];
        
        //integrated dissipation
        profiles[58][ix]+=vischeating*dxph[1];

	//projected TOP jet CM and MOI tensor (according to code comparison paper)
	//assumes flat space
	if(xxBL[2]>0. && xxBL[2]<=thetamin_134)
	{
	  ldouble sigloc = bsq/rho;
	  ldouble psiloc=0.;
	  if(sigloc>1) psiloc=1.;
	  else if(sigloc>0.1) psiloc=sigloc;
	  
	  ldouble xloc=xxBL[1]*sin(xxBL[2])*cos(xxBL[3]);
	  ldouble yloc=xxBL[1]*sin(xxBL[2])*sin(xxBL[3]); 
	  ldouble dA=psiloc*xxBL[1]*xxBL[1]*sin(xxBL[2])*cos(xxBL[2])*dx[1]*dx[2];
	  
	  profiles[59][ix] += xloc*dA;
	  profiles[60][ix] += yloc*dA;
	  profiles[61][ix] += xloc*xloc*dA;
	  profiles[62][ix] += yloc*yloc*dA;
	  profiles[63][ix] += xloc*yloc*dA;
	  profiles[64][ix] += dA;

	}
#ifdef RADIATION
        //total rad energy flux (17)
	profiles[15][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;
        
        //conserved flux for Rrt in MYCOORDS (40)
        profiles[38][ix]+=-Fluxx[EE0]*get_size_x(iy,1)*dx[2];
        
        //rad viscosity energy flux (35)
        profiles[33][ix]+=(-Rviscrt)*dx[1]*dx[2]*geomBL.gdet;
        
        //opt. thin rad energy flux (26)
        if(tautot<1.)
          profiles[24][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;
        
        //outflowin rad energy flux (18)
        if(utcon[1]>0.)
          profiles[16][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;
        
        //jet rad energy flux (19)
        if(isjet==1)
          profiles[17][ix]+=(-Rrt)*dx[1]*dx[2]*geomBL.gdet;
#endif
        //outflowin mass flux (20)
        if(utcon[1]>0.)
          profiles[18][ix]+=(-rhouconr)*dx[1]*dx[2]*geomBL.gdet;
        
        //opt. thin mass flux (27)
        if(tautot<1.)
          profiles[25][ix]+=(-rhouconr)*dx[1]*dx[2]*geomBL.gdet;
        
        //jet mass flux (21)
        if(isjet==1)
          profiles[19][ix]+=(-rhouconr)*dx[1]*dx[2]*geomBL.gdet;
        
        //rho-weighted minus radial velocity (4)
        profiles[2][ix]+=-rhouconr*dxph[1];
        
        //energy-averaged radial transport velocity (53)
        if(Be<0.)
          profiles[51][ix]+=endenuconr*dxph[1];
        
        //energy-averaged vertical transport velocity (54)
        if(Be<0.)
          profiles[52][ix]+=my_sign(geomBL.xxvec[2]-M_PI/2.)*endenuconth*dxph[1];
        
        //rho-weighted minus radial velocity in the inflow (24)
        if(utcon[1]<0.)
          profiles[22][ix]+=-rhouconr*dxph[1];
        
        //rho-weighted minus radial velocity in the outflow (36)
        #if(PROBLEM==87)
        ldouble angmomR = ucov[3];
        ldouble epsR =-(1.+ucov[0]);
        ldouble eccR = 2.*angmomR*angmomR*epsR;
        if(eccR > -1.) 
        {
          profiles[34][ix]+= rho*dxph[1]*sqrt(1. - 2.*ucov[3]*ucov[3]*(1.+ucov[0]));
          profiles[21][ix]+=rho*dxph[1];
        }
        //if(eccN > -1.) profiles[34][ix]+= rho*dxph[1]*sqrt(1. + eccN); //eccentricity
        #else
        if(utcon[1]<0.)
          profiles[21][ix]+=rho*dxph[1];
        if(utcon[1]>0.)
          profiles[34][ix]+=rhouconr*dxph[1];
        #endif

        
        //abs optical depth (7)
        profiles[5][ix]=tauabs;
        
        //tot optical depth (8)
        profiles[6][ix]=tautot;
        
        for(iv=0;iv<NV+NAVGVARS;iv++)
          avgsums[iv][ix]+=get_uavg(pavg,iv,ix,iy,iz)*dxph[1];
        
        //summing up 3d angular momentum
        ldouble Jx,Jy,Jz,th,ph,r;
        ldouble xxveccar[4],velvec[3],posvec[3],jvec[3];
        r=geomBL.xx;
        th=geomBL.yy;
        ph=geomBL.zz;
        coco_N(geomBL.xxvec,xxveccar,BLCOORDS,MINKCOORDS);
        posvec[0]=xxveccar[1];
        posvec[1]=xxveccar[2];
        posvec[2]=xxveccar[3];
        velvec[0]=sin(th)*cos(ph)*rho*utcon[1]
        + cos(th)*cos(ph)*rho*utcon[2]
        - sin(ph)*rho*utcon[3];
        velvec[1] = sin(th)*sin(ph)*rho*utcon[1]
        + cos(th)*sin(ph)*rho*utcon[2]
        + cos(ph)*rho*utcon[3];
        velvec[1]*=r;
        velvec[2] = cos(th)*rho*utcon[1]
        - sin(th)*rho*utcon[2];
        velvec[2]*=r*sin(th);
        jvec[0]=posvec[1]*velvec[2]-posvec[2]*velvec[1];
        jvec[1]=-(posvec[0]*velvec[2]-posvec[2]*velvec[0]);
        jvec[2]=posvec[0]*velvec[1]-posvec[1]*velvec[0];
        
        Jvec[0]+=jvec[0];
        Jvec[1]+=jvec[1];
        Jvec[2]+=jvec[2];
        
        //<(rho+u+bsq/2)u^r><u_phi> (5)
        //profiles[3][ix]+=get_uavg(pavg,AVGWUCON(1),ix,iy,iz)*get_uavg(pavg,AVGUCOV(3),ix,iy,iz)*dxph[1];
        
        //u_phi everywhere (5)
        //if(utcon[1]<0.)
        {
          if(doingavg)
            profiles[3][ix]+=get_uavg(pavg,AVGRHOUCOV(3),ix,iy,iz)*dxph[1];
          else
            profiles[3][ix]+=rho*ucov[3]*dxph[1];
        }
        
        //int the outflow (51)
        if(utcon[1]>0.)
        {
          if(doingavg)
            profiles[49][ix]+=get_uavg(pavg,AVGRHOUCOV(3),ix,iy,iz)*dxph[1];
          else
            profiles[49][ix]+=rho*ucov[3]*dxph[1];
        }
        
        // special quantities for problem 134//139
        if (PROBLEM == 134 || PROBLEM == 139)
	{
	  if(xxBL[2]<thetamax_134 && xxBL[2]>thetamin_134)
          {
          normalize[ix] += dx[1] * dx[2] * geomBL.gdet;
          rho134[ix] += rho * dx[1] * dx[2] * geomBL.gdet;
          ldouble pgass = uint * (GAMMA - 1.);
          pgas134[ix] += pgass * dx[1] * dx[2] * geomBL.gdet;
          bfield134[ix] += sqrt(bsq) * dx[1] * dx[2] * geomBL.gdet;
          uconphi134[ix] += utcon[3] * dx[1] * dx[2] * geomBL.gdet;
          ptot134[ix] += (pgass + (bsq/2.)) * dx[1] * dx[2] * geomBL.gdet;
          betainv134[ix] += (0.5 * bsq / pgass) * dx[1] * dx[2] * geomBL.gdet;
          }
        }
      }  // end for(iy=0;iy<NY;iy++)
    }  // end for(iz=0;iz<NZ;iz++)
    
    //calculating angle between the total angular momentum vector at this radius and the axis
    ldouble angle = acos( (Jvec[2])/(sqrt(Jvec[0]*Jvec[0]+Jvec[1]*Jvec[1]+Jvec[2]*Jvec[2])));
    profiles[53][ix]=angle;
        
    //normalizing by sigma
    ldouble sigmaout=profiles[0][ix]-profiles[21][ix];
    ldouble sigmain=profiles[21][ix];
    profiles[2][ix]/=profiles[0][ix]; //total
    profiles[51][ix]/=profiles[39][ix];
    profiles[52][ix]/=profiles[39][ix];
    
    profiles[49][ix]/=sigmaout; //normalized by surface density in the outlfow
    profiles[22][ix]/=profiles[21][ix]; //sigma in
    #if(PROBLEM==87)
    profiles[34][ix]/=profiles[21][ix];
    #else
    profiles[34][ix]/=sigmaout;
    #endif
    profiles[26][ix]/=profiles[0][ix];
    profiles[49][ix]/=(profiles[0][ix]-profiles[21][ix]); //normalized by surface density in the outlfow
    profiles[45][ix]/=profiles[0][ix];
    profiles[27][ix]/=profiles[0][ix];
    profiles[48][ix]/=profiles[0][ix];
    profiles[44][ix]/=profiles[0][ix];
    profiles[57][ix]/=profiles[0][ix];
#ifndef CALCHRONTHEGO
    profiles[29][ix]/=profiles[0][ix];
    profiles[29][ix]=sqrt(profiles[29][ix]); //scale height
#else
    profiles[29][ix]=scaleth_otg[ix]; //scale height
#endif
    
    profiles[30][ix]/=profiles[0][ix];
    profiles[31][ix]/=profiles[0][ix];
    profiles[32][ix]/=profiles[0][ix];
    
    profiles[40][ix]/=jetsigma;
    
    Bangle1/=profiles[0][ix];
    Bangle2/=profiles[0][ix];
    
    //rho-weighted magn.field angle -<sqrt(grr gphph)b^r b^ph> / <bsq> (30)
    profiles[28][ix]=-Bangle1/Bangle2;
    
    //normalizing by <(rho+u+bsq/2)u^r>
    //profiles[3][ix]/=avgsums[AVGWUCON(1)][ix];
    profiles[3][ix]/=profiles[0][ix];
    
    //Keplerian u_phi (6)
    ldouble r=xxBL[1];
    profiles[4][ix]=((r*r-2.*BHSPIN*sqrt(r)+BHSPIN*BHSPIN)/(sqrt(r*(r*r-3.*r+2.*BHSPIN*sqrt(r)))));  
    
    //net accretion rate at given radius (9)
    profiles[7][ix]=fabs(calc_mdot(xxBL[1],0));
    //inflow accretion rate at given radius (10)
    profiles[8][ix]=fabs(calc_mdot(xxBL[1],1));
    //outflow accretion rate at given radius (11)
    profiles[9][ix]=fabs(calc_mdot(xxBL[1],2));
    
    //luminosity at given radius (12) -- outside photosphere
    ldouble radlum,totallum;
    calc_lum(xxBL[1],0,&radlum,&totallum);
    profiles[10][ix]=radlum;
    //luminosity at given radius (22) -- everywhere
    calc_lum(xxBL[1],1,&radlum,&totallum);
    profiles[20][ix]=radlum;
    //location of the photosphere (13)
    profiles[11][ix]=calc_photloc(ix);
    
    // special quantities for problem 134/139
    if (PROBLEM == 134 || PROBLEM == 139)
    {
      profiles[7][ix] = rho134[ix] / normalize[ix];
      profiles[8][ix] = pgas134[ix] / normalize[ix];
      profiles[12][ix] = bfield134[ix] / normalize[ix];
      profiles[13][ix] = uconphi134[ix] / normalize[ix];
      profiles[17][ix] = ptot134[ix] / normalize[ix];
      profiles[18][ix] = betainv134[ix] / normalize[ix];
      profiles[29][ix] = scaleheight134[ix] / sigmagdet[ix];
      profiles[27][ix] = rholambda134[ix] / sigmagdet[ix];
    }

    // TOP jet profiles
    ldouble Ajet = profiles[64][ix];
    if (Ajet==0.) Ajet=1.;
    
    profiles[59][ix] /= Ajet;
    profiles[60][ix] /= Ajet;
    profiles[61][ix] -= (profiles[59][ix]*profiles[59][ix]*Ajet);
    profiles[62][ix] -= (profiles[60][ix]*profiles[60][ix]*Ajet);
    profiles[63][ix] -= (profiles[59][ix]*profiles[60][ix]*Ajet);
    
  }  // end for(ix=0;ix<NX;ix++)

  return 0;
}


/*********************************************/
/* calculates theta profiles  */
//total energy flux (2) (column)
/*********************************************/

int calc_thetaprofiles(ldouble profiles[][NY])
{
  //adjust NTHPROFILES in problem.h
  ldouble rho,uint,bsq,bcon[4],bcov[4],utcon[4],ucov[4],rhouconr,rhoucont,rhouconrucovt;
  ldouble Tij[4][4],Tij22[4][4],Rij[4][4],Trt,Trtkin,Rrt,Ehat,Rviscij[4][4],Rviscrt,Rtt,Ttt;
  ldouble pp[NV];

  //choose radius where to extract from
  int ix,i,j,iv,iz;

  //search for appropriate radial index
  ldouble xx[4],xxBL[4];
  ldouble radius=1.e3;
  #ifdef THPROFRADIUS
  radius=THPROFRADIUS;
  #endif
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      #ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix, 0, 0, xxBL);
      #else
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
      #endif

      if(xxBL[1]>radius) break;
    }

  int iy;
#pragma omp parallel for private(iy)
  for(iy=0;iy<NY;iy++)
    {
      for(iv=0;iv<NTHPROFILES;iv++)
	profiles[iv][iy]=0.;

      if(NZ==1) //phi-symmetry
	{
	  iz=0;
	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	  struct geometry geom;
	  fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);
	      
	  //primitives at the cell - either averaged or original, in BL or MYCOORDS
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  //to BL, res-files and primitives in avg in MYCOORDS	  
          #ifdef PRECOMPUTE_MY2OUT
          trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
          #else      
          trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, xx,&geom,&geomBL);
          #endif

	  if(doingavg)
	    {
	      rho=get_uavg(pavg,RHO,ix,iy,iz);
	      uint=get_uavg(pavg,UU,ix,iy,iz);
	      bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
	      bcon[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
	      bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
	      bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
	      bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
	      utcon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      rhouconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
	      rhoucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz);
	      rhouconrucovt = get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz);

	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)

		  Tij[i][j]=get_uavg(pavg,AVGRHOUCONUCOV(i,j),ix,iy,iz)
		    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(i,j),ix,iy,iz)
		    + get_uavg(pavg,AVGBSQUCONUCOV(i,j),ix,iy,iz)
		    - get_uavg(pavg,AVGBCONBCOV(i,j),ix,iy,iz); 

	      Trt=Tij[1][0];
	      Ttt=Tij[0][0];

#ifdef RADIATION  
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 

	      Rrt = Rij[1][0];
	      Rtt = Rij[0][0];
	      Ehat = get_uavg(pavg,AVGEHAT,ix,iy,iz);

	      int derdir[3]={0,0,0};
	      calc_Rij_visc(pp,&geomBL,Rviscij,derdir);      
	      Rviscrt = Rviscij[1][0];
#endif
		  
	      //no need of transforming interpolated primitives to BL, already there
		 
	    }
	  else
	    { 
	      rho=pp[0];
	      uint=pp[1];
	      utcon[1]=pp[2];
	      utcon[2]=pp[3];
	      utcon[3]=pp[4];
		  
#ifdef MAGNFIELD
	      calc_bcon_prim(pp,bcon,&geomBL);
	      indices_21(bcon,bcov,geomBL.gg); 
	      bsq = dotB(bcon,bcov); 
#endif

	      conv_vels_both(utcon,utcon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
	      rhouconr=rho*utcon[1];
	      rhoucont=rho*utcon[0];
	      

	      calc_Tij(pp,&geomBL,Tij22);
	      indices_2221(Tij22,Tij,geomBL.gg);

	      Trt = Tij[1][0];
	      Ttt = Tij[0][0];
	      rhouconrucovt =  rho*utcon[1]*ucov[0];
		 
#ifdef RADIATION
	      calc_Rij(pp,&geomBL,Rij);
	      indices_2221(Rij,Rij,geomBL.gg);

	      Rrt = Rij[1][0];
	      Rtt=Rij[0][0];

	      ldouble Rtthat,uconr[4];
	      calc_ff_Rtt(&get_u(p,0,ix,iy,iz),&Rtthat,uconr,&geomBL);
	      Ehat=-Rtthat; 	

	      int derdir[3]={0,0,0}; 
	      calc_Rij_visc(pp,&geomBL,Rviscij,derdir);
      
	      Rviscrt = Rviscij[1][0];
#endif
	      
	    }
	  
	  ldouble muBe,Be;
	  muBe=-(Trt+rhouconr)/rhouconr;
	  Be=-(Ttt+rhoucont)/rhoucont;
#ifdef RADIATION
	  muBe+=-Rrt/rhouconr;
	  Be+=-Rtt/rhoucont;
#endif

	  int isconverged;
	  if(fabs(utcon[1])>xxBL[1]/(global_time/2.)) 
	    isconverged=1;
	  else
	    isconverged=0;
	  
	  ldouble fluxconv=fluxGU2CGS(1.); 

	  //total energy flux in cgs  (2)
	  profiles[0][iy]=-fluxconv*(Rrt+rhouconr+Trt);
	      
	  //radiative flux in cgs (3)
	  profiles[1][iy]=-fluxconv*(Rrt);

	  //kinetic flux in cgs (8)
	  //ldouble binding = sqrt(-1./geomBL.GG[0][0]);
	  //Trtkin=rhouconrucovt - rhouconr*binding;
	  Trtkin = rhouconr*(utcon[0]-1.);
	  profiles[6][iy]=fluxconv*(Trtkin);

	  //gas velocity (4)
	  profiles[2][iy]=utcon[1];

	  //converged? (5)
	  profiles[3][iy]=(ldouble)isconverged;

	  //optical depths
	  int iix;
	  ldouble tau1,tau2,taueff1,taueff2,Rphot,Rphoteff;
	  tau1=tau2=taueff1=taueff2=0.;
	  Rphot=Rphoteff=-1.;
	  struct opacities opac;
  

	  for(iix=NX-1;iix>=0;iix--)
	    {
	      struct geometry geomBL2;
	      fill_geometry_arb(iix,iy,iz,&geomBL2,OUTCOORDS);
	      ldouble grr=geomBL2.gg[1][1];

	      //cell dimensions
	      ldouble dxph[3],dx[0];
	      
    	      //ANDREW put cell size code in a function with precompute option
              get_cellsize_out(iix, iy, iz, dx);
	      /*
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(iix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(iix+1,0);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,OUTCOORDS); 
	      coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      */
		
	      dxph[0]=dx[0]*sqrt(geomBL2.gg[1][1]);

	      ldouble rho2=get_u(p,RHO,iix,iy,iz);
              ldouble uint2=get_u(p,UU,iix,iy,iz);
              ldouble Tgas2=calc_PEQ_Tfromurho(uint2,rho2,ix,iy,iz);
	      ldouble kabsloc = calc_kappa(&get_u(p,0,iix,iy,iz),&geomBL2,&opac);
	      ldouble kscaloc = calc_kappaes(&get_u(p,0,iix,iy,iz),&geomBL2);
	      if(doingavg)
                {
                  utcon[0]=get_uavg(pavg,AVGRHOUCON(0),iix,iy,iz)/get_uavg(pavg,RHO,iix,iy,iz);
                  utcon[1]=get_uavg(pavg,AVGRHOUCON(1),iix,iy,iz)/get_uavg(pavg,RHO,iix,iy,iz);
                  utcon[2]=get_uavg(pavg,AVGRHOUCON(2),iix,iy,iz)/get_uavg(pavg,RHO,iix,iy,iz);
                  utcon[3]=get_uavg(pavg,AVGRHOUCON(3),iix,iy,iz)/get_uavg(pavg,RHO,iix,iy,iz);
		  conv_vels_both(utcon,utcon,ucov,VELPRIM,VEL4,geomBL2.gg,geomBL2.GG);
		}
	      else
		{
		  //temporary
		  ucov[0]=-1.; ucov[1]=0.;
		}

	      if(geomBL2.xx>rhorizonBL)
		{
		  tau1+=-(kabsloc+kscaloc)*(ucov[0]+ucov[1])*sqrt(grr)*dxph[0];
		  tau2=-(kabsloc+kscaloc)*(ucov[0]+ucov[1])*geomBL2.xx*sqrt(grr);
		  taueff1+=-sqrt(kabsloc*(kabsloc+kscaloc))*(ucov[0]+ucov[1])*sqrt(grr)*dxph[0];
		  taueff2=-sqrt(kabsloc*(kabsloc+kscaloc))*(ucov[0]+ucov[1])*geomBL2.xx*sqrt(grr);
		  if(Rphot<0. && my_max(tau1,tau2)>2./3.)
		    Rphot=geomBL2.xx;
		  if(Rphoteff<0. && my_max(taueff1,taueff2)>2./3.)
		    Rphoteff=geomBL2.xx;
		}
	      
	    }

	  //location of photosphere (6)
	  profiles[4][iy]=Rphot;

	  //location of the effective photosphere (7)
	  profiles[5][iy]=Rphoteff;

	  //Be parameter (9)
	  profiles[7][iy]=Be;

	}

    }

  return 0;
}


/*********************************************/
/* 
 Calculates scalars - total mass, accretion rate etc. These are output in scalars.dat, where scalars[i] is written in column i+2
 
  scalars[0]: total mass in the domain
  scalars[1]: mdot in GU
  scalars[2]: some kind of luminosity
  scalars[3]: magnetic flux through horizon
  scalars[4]: mdot in Eddington units
  scalars[5]: MRI Q_theta at some radius
  scalars[6]: mean temperature at some radius
  scalars[7]: MAD parameter (dipole)
  scalars[8]: scale height at some radius
  scalars[9]: MAD parameter (quadrupole)
  scalars[10]: luminosity proxy for problem 134,139
  scalars[11]: Edot for problem 134,139,140
  scalars[12]: Ldot for problem 134,139,140
 
 */
/*********************************************/

int calc_scalars(ldouble *scalars,ldouble t)
{
  //adjust NSCALARS in problem.h

  //total mass inside the domain (2nd column)
  scalars[0]=calc_totalmass();

  //accretion rate through horizon (3)
  ldouble mdotscale = rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.);

 ldouble rmdot=rhorizonBL;
#if(PROBLEM==7) //BONDI
  rmdot = RMIN;
#endif
  ldouble mdot=calc_mdot(rmdot,0);
  scalars[1]=-mdot;

  //accretion rate through horizon in Edd. units (6)
  scalars[4]=-mdot*mdotscale/calc_mdotEdd();

  //luminosities 
  ldouble rlum=5000.;
  ldouble taulimit=2./3.;
#if(PROBLEM==139) //MADCC
  rlum=15.;
#endif
#if(PROBLEM==69) //INJDISK
  rlum=2./3.*DISKRCIR;
#endif
#if(PROBLEM==7) //BONDI
  rlum=RMAX;
#endif
#if(PROBLEM==7) //RADTORUS
  rlum=13.;
#endif
#if(PROBLEM==91 || PROBLEM==94)  //TDEMILIO
  rlum=0.99*ROUT;
#endif
  ldouble radlum,totallum;
  calc_lum(rlum,1,&radlum,&totallum);
  #if(PROBLEM==87)
  calc_lum_tausurface(taulimit,&radlum);
  #endif

  //radiative luminosity everywhere (4)
  scalars[2]=radlum*mdotscale*CCC0*CCC0/calc_lumEdd();

  if(PROBLEM==89 || PROBLEM==79) //RADTORUS or SOFTBALL
  {
    //luminosity exiting through radial and theta boundaries (3)
    scalars[1]=calc_exitlum();
  }

  //total energy flux at infinity (rho ur + Trt + Rrt) (12)
  scalars[10]=totallum;
  
  //mri resolution parameter Q_theta (7) at rmri

  ldouble xx[4],xxBL[4];
  get_xx(NX-1,0,0,xx);
  #ifdef PRECOMPUTE_MY2OUT
  get_xxout(NX-1, 0, 0, xxBL);
  #else
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
  #endif
  
  ldouble rmri=xxBL[1]/2.;  // middle radial cell
#if(PROBLEM==69 || PROBLEM==140) //INJDISK or RADSURVEY
  rmri=20.;
#endif
  scalars[5]=calc_resmri(rmri);

  //rho-weighted temperature at rtemp (8)
  ldouble rtemp=15.;  // radius at which temperature is calculated
#if(PROBLEM==69) //INJDISK
  rtemp=20.;
#endif
#if(PROBLEM==7) //BONDI
  rtemp=RMAX;
#endif
  scalars[6]=calc_meantemp(rtemp);

  //magnetic flux through horizon parameter (5)
  ldouble Bfluxquad;
  ldouble Bflux;
  calc_Bflux(rhorizonBL,0.,&Bflux,&Bfluxquad);
  scalars[3]=Bflux;

  //MAD parameter (9)
  scalars[7]=Bflux/sqrt(fabs(mdot))*sqrt(4.*M_PI)/2.;

  //MAD-quad parameter (11)
  scalars[9]=(Bfluxquad/sqrt(fabs(mdot))*sqrt(4.*M_PI)/2.);

  //scaleheight at rscale (10)
  ldouble rscale=15.;
  scalars[8]=calc_scaleheight(rscale);

  //brightness at R=THPROFRADIUS (13)
  //fixed polar index - for power spectrum calculation                                                                                        
  //search for appropriate radial index                                                                                                                       
  ldouble radius=5.e3;
  #ifdef THPROFRADIUS
  radius=THPROFRADIUS;
  #endif
  int ix;
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      #ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix, 0, 0, xxBL);
      #else
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
      #endif

      if(xxBL[1]>radius) break;
    }

  double totlum;
  calc_local_lum(ix,NCCORRECTPOLAR+1,0,&radlum,&totlum);
  scalars[11]=totlum;
  
#if(PROBLEM==134 || PROBLEM==139 || PROBLEM==140)  // FISHMONC for EHT code comparison tests
  ldouble Edot;
  Edot = calc_Edot(rhorizonBL);
  scalars[11] = Edot - mdot;
  
  ldouble Ldot;
  Ldot = calc_Ldot(rhorizonBL);
  scalars[12] = Ldot;

  ldouble luminosity_proxy, theta_min = M_PI / 3., theta_max = 2. * M_PI / 3.;
  luminosity_proxy = calc_lum_proxy(rhorizonBL, theta_min, theta_max);
  scalars[10] = luminosity_proxy;
#endif

  // accretion rates through the outer edge for TDEMILIO 
#if(PROBLEM==91 || PROBLEM==94 || PROBLEM==87)
  rmdot = 500.;
  
  //total outflow (12)
  mdot=calc_mdot(rmdot,3);
  scalars[10]=mdot*mdotscale/calc_mdotEdd();
  
  //unbound outflow (13)
  mdot=calc_mdot(rmdot,1);
  scalars[11]=-mdot*mdotscale/calc_mdotEdd();
  
#endif

  /*********************************************/
  //Tgas Trad Egas Erad for testing Comptonization
  /*********************************************/

#ifdef TESTCOMPTINSCALARS4FLAT
  ldouble pp[NV];
  int iv;
  PLOOP(iv) pp[iv]=get_u(p,iv,0,0,0);
  struct geometry geom;
  fill_geometry(0,0,0,&geom);
  ldouble ugas[4],Te,Ti,Tgas,Trad,Tradbb,Rtt,Ehatrad;
  int dominates;
  calc_ff_Rtt(pp,&Rtt,ugas,&geom);
  Ehatrad=-Rtt;
  calc_PEQ_Teifrompp(pp,&Te,&Ti,0,0,0);
  Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],0,0,0);
  #ifndef EVOLVEPHOTONNUMBER
  Trad=Tradbb=calc_LTE_TfromE(Ehatrad);
  #else
  Tradbb=calc_LTE_TfromE(Ehatrad);
  Trad=calc_ncompt_Thatrad(pp,&geom,Ehatrad);
  #endif

  scalars[0]=Tgas;
  scalars[1]=Trad;
  scalars[2]=Te; //(4th)
  scalars[3]=Ti;
  scalars[4]=Tradbb;
  scalars[5]=numdensGU2CGS(pp[NF0]);
  scalars[6]=rhoGU2CGS(pp[RHO]);
#endif

  //ANDREW should this be in not just one cell??  
#if (PROBLEM==107 || PROBLEM==114) //RELELTEST or RELELEXPAND
  ldouble pp[NV];
  int iv;
  PLOOP(iv) pp[iv]=get_u(p,iv,2,0,0);
  struct geometry geom;
  fill_geometry(0,0,0,&geom);
  ldouble ugas[4],Te,Ti,Tgas,Trad,Tradbb,Rtt,Ehatrad;
  int dominates;
  calc_ff_Rtt(pp,&Rtt,ugas,&geom);
  Ehatrad=-Rtt;
  calc_PEQ_Teifrompp(pp,&Te,&Ti,0,0,0);
  Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO],0,0,0);
  #ifndef EVOLVEPHOTONNUMBER
  Trad=Tradbb=calc_LTE_TfromE(Ehatrad);
  #else
  Tradbb=calc_LTE_TfromE(Ehatrad);
  Trad=calc_ncompt_Thatrad(pp,&geom,Ehatrad);
  #endif
  
  int is;
  ldouble uth,uthi,uur;
  for (is=0;is<NSCALARS;is++) scalars[is]=0;

  uur=calc_relel_uint(pp);
  ldouble nur=calc_relel_ne(pp);
  ldouble ni=pp[RHO]/M_PROTON/MU_I;
  ldouble neth=calc_thermal_ne(pp);
  ldouble rhoeth=neth*M_PROTON*MU_E;
  
  uth  = calc_ufromSerho(pp[ENTRE],rhoeth,ELECTRONS,0,0,0);
  uthi = calc_ufromSerho(pp[ENTRI],pp[RHO],IONS,0,0,0);

  scalars[0]=Te; //(2nd)
  scalars[1]=Ti;
  scalars[2]=Trad;
  scalars[3]=numdensGU2CGS(calc_thermal_ne(pp));
  scalars[4]=endenGU2CGS(Ehatrad);
  scalars[5]=endenGU2CGS(uur);
  scalars[6]=endenGU2CGS(uth);
  scalars[7]=endenGU2CGS(uthi);
  scalars[8]=endenGU2CGS(pp[UU]);
  scalars[9]=numdensGU2CGS(calc_relel_ne(pp));
  scalars[10]=pp[ENTRE];
  scalars[11]=pp[ENTRI];
  int end=12;
 #ifdef RELELECTRONS
  int ie;
  for (ie=0 ; (ie < NRELBIN) && (ie+7 < NSCALARS); ie++)
    scalars[end+ie]=numdensGU2CGS(pp[NEREL(ie)]);
 #endif 
#endif


#if(PROBLEM==106 || PROBLEM==118) //SEDRIVEN or SEDRIVENTURBRELEL

  //total integrated internal energy
  ldouble ugtot=0;
  ldouble vischeattot=0.,vischeatnegetot=0.,kinentot=0.;
  int i,j;
  for(i=0;i<NX;i++)
    { 
      for(j=0;j<NY;j++)
	{
	  struct geometry geomBL;
	  fill_geometry_arb(i,j,0,&geomBL,OUTCOORDS);
	  ldouble dV=get_size_x(i,0)*get_size_x(j,1);
	  ugtot+=get_u(p,UU,i,j,0)*dV;

	  vischeattot+=get_u_scalar(vischeating,i,j,0)*dV;
	  vischeatnegetot+=get_u_scalar(vischeatingnegebalance,i,j,0)*dV;
	}
    }

  scalars[7]=vischeatnegetot; //(9)
  scalars[8]=vischeattot; //(10)
  scalars[9]=ugtot; //(11)

  //rms velocity
  ldouble rmsvel=0.;
  for(i=0;i<NX;i++)
    { 
      for(j=0;j<NY;j++)
	{
	  struct geometry geomBL;
	  fill_geometry_arb(i,j,0,&geomBL,OUTCOORDS);
	  ldouble dV=get_size_x(i,0)*get_size_x(j,1);
	  rmsvel+=(get_u(p,VX,i,j,0)*get_u(p,VX,i,j,0)+get_u(p,VY,i,j,0)*get_u(p,VY,i,j,0))*dV;
	  kinentot+=get_u(p,RHO,i,j,0)*(get_u(p,VX,i,j,0)*get_u(p,VX,i,j,0)+get_u(p,VY,i,j,0)*get_u(p,VY,i,j,0))*dV;
	}
    }
  rmsvel=sqrt(rmsvel);
  scalars[10]=rmsvel; //(12)

  //electron integrated internal energy
  ldouble uetot=0;
  ldouble uitot=0;
  ldouble uurtot=0;

  //average temperatures
  ldouble Tgav,Teav,Tiav,netot;
  ldouble pspecies,ptot,uspecies,utot;
  Tgav=Teav=Tiav=netot=0.;
  pspecies=ptot=uspecies=utot=0.;

#ifdef RELELECTRONS
         int ie;
         int end=12;
         for (ie=0 ; (ie < NRELBIN) && (ie+end < NSCALARS); ie++)
         scalars[end+ie] = 0.0;
#endif

  for(i=0;i<NX;i++)
    { 
      for(j=0;j<NY;j++)
	{
	  struct geometry geomBL;
	  fill_geometry_arb(i,j,0,&geomBL,OUTCOORDS);
	  ldouble dV=get_size_x(i,0)*get_size_x(j,1);
	  ldouble Tg,Te,Ti;
	  ldouble *pp=&get_u(p,0,i,j,0);
	  Tg=calc_PEQ_Teifrompp(pp,&Te,&Ti,i,j,0); 


	  ldouble rho=get_u(p,RHO,i,j,0);
	  ldouble uint=get_u(p,UU,i,j,0);

	  ldouble gamma=GAMMA;
#ifdef CONSISTENTGAMMA
	  gamma=pick_gammagas(i,j,0);
#endif
	  ldouble gammam1=gamma-1.;

          //total pressure
	  ldouble pre=gammam1*uint;
	  /**************/
	  //electrons
	  /**************/
	  ldouble ne=calc_thermal_ne(pp); //thermal only
	  ldouble pe=K_BOLTZ*ne*Te;
	  ldouble gammae=GAMMAE;
#ifdef CONSISTENTGAMMA
#ifndef FIXEDGAMMASPECIES
	  gammae=calc_gammaintfromtemp(Te,ELECTRONS);
#endif
#endif	  
	  ldouble ue=pe/(gammae-1.);

	  /**************/
	  //ions
	  /**************/
	  ldouble ni=rho/MU_I/M_PROTON; //number density of photons and electrons
	  ldouble pi=K_BOLTZ*ni*Ti;
	  ldouble gammai=GAMMAI;
#ifdef CONSISTENTGAMMA
#ifndef FIXEDGAMMASPECIES
	  gammai=calc_gammaintfromtemp(Ti,IONS);
#endif
#endif
	  ldouble ui=pi/(gammai-1.);

	  netot+=ne*dV;
	  Tgav+=ne*Tg*dV;
	  Teav+=ne*Te*dV;
	  Tiav+=ne*Ti*dV;

	  uspecies+=(ue+ui)*dV;
	  pspecies+=(pe+pi)*dV;
	  utot+=uint*dV;
	  ptot+=pre*dV;
	  uetot+=ue*dV;
          uitot+=ui*dV;

#ifdef RELELECTRONS
	 uurtot+=dV*calc_relel_uint(&get_u(p,0,i,j,0));
         for (ie=0 ; (ie < NRELBIN) && (ie+end < NSCALARS); ie++)
	 scalars[end+ie] += dV*get_u(p,NEREL(ie),i,j,0);
#endif
	}
    }

  scalars[11]=uetot; //(13)

  Tgav/=netot;
  Teav/=netot;
  Tiav/=netot;

  scalars[0]=Tgav; //(2)
  scalars[1]=Teav; //(3)
  scalars[2]=Tiav; //(4)
  scalars[3] = uetot;
  scalars[4] = uitot;
  scalars[5] = uurtot;

  scalars[6]=calc_totalmass(); //(8)

#ifdef RELELECTRONS
  for (ie=0 ; (ie < NRELBIN) && (ie+end < NSCALARS); ie++)
  scalars[end+ie] /= netot;
#endif


#endif //SEDDRIVEN or SEDRIVENTURBRELEL


  /*********************************************/
  //L2 measure
  /*********************************************/

#if(PROBLEM==79) //SOFTBALL

  //L2 norm of density
  ldouble L2=0;
  int i,j;
  for(i=0;i<NX;i++)
    { 
      for(j=0;j<NY;j++)
	{
	  struct geometry geomBL;
	  fill_geometry_arb(i,j,0,&geomBL,OUTCOORDS);
	  ldouble dV=get_size_x(i,0)*get_size_x(j,1)*2.*M_PI*geomBL.gdet;
	  L2+=get_u(p,RHO,i,j,0)*get_u(p,RHO,i,j,0)*dV;
	}
    }
  scalars[9]=L2;///(ldouble)NX;

  //average quantities
  //beta = p_rad / p_gas
  //temp, angular momentum
  ldouble beta=0;
  ldouble rho=0;
  ldouble pp[NV];
  int iv;
  struct geometry geom,geomBL;
  ldouble dV,Ehat,Rtt,prad,ugas[4],pgas,betaloc,rho2;

   for(i=0;i<NX;i++)
    { 
      for(j=0;j<NY;j++)
	{
	  PLOOP(iv) pp[iv]=get_u(p,iv,i,j,0);
	  fill_geometry(i,j,0,&geom);
	  calc_ff_Rtt(pp,&Rtt,ugas,&geom);
	  Ehat=-Rtt;
	  prad=1./3.*Ehat;
	  pgas=GAMMAM1*pp[UU];
	  dV=get_size_x(i,0)*get_size_x(j,1)*2.*M_PI*geom.gdet;
	  betaloc=prad/pgas;
	  rho=pp[RHO];
	  beta+=rho*rho*dV*betaloc;
	  rho2+=rho*rho*dV;
	}
    }
   scalars[3]=beta/rho2;

#endif //SOFTBALL

  /*********************************************/
  //L1 ERRRORS for some problems
  /*********************************************/

#ifdef CALCL1_RMHDWAVE
  //temporarily here: L1 error for RMHDWAVE
  ldouble L1=0;
  int i,j;
  for(i=0;i<NX;i++)
    {
      calc_primitives(i,0,0,0,0);
      ldouble xx=get_x(i,0);
      ldouble dx=get_size_x(i,0);
      ldouble myrho=RHOZERO+DRRE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DRIM/DRRE*sin(OMRE*t-KK*xx));
      //L1 in rho:
      L1+=fabs(get_u(p,RHO,i,0,0)-myrho)*dx;
    }
  scalars[0]=L1;///(ldouble)NX;
#endif

#ifdef CALCL1_HDWAVE
  //temporarily here: L1 error for HDWAVE
  ldouble L1=0;
  int i;
  for(i=0;i<NX;i++)
    {
      calc_primitives(i,0,0,0,0);
      ldouble xx=get_x(i,0);
      ldouble om=1./CC*2.*Pi;
      ldouble myrho=RHOZERO*(1.+AAA*cos(KK*xx-om*t));
      ldouble myuint=UINT*(1.+GAMMA*AAA*cos(KK*xx-om*t));
      ldouble mycs=1./CC;
      ldouble myvx=AAA*cos(KK*xx-om*t)*mycs;
      //L1 in rho:
      L1+=fabs(get_u(p,0,i,0,0)-myrho);
    }
  scalars[0]=L1/(ldouble)NX;
#endif //CALCL1_HDWAVE

#ifdef CALCL1_HUBBLE
  //temporarily here: L1 error for HUBBLE
  ldouble L1=0;
  int i;
  for(i=0;i<NX;i++)
    {
      ldouble xx=get_x(i,0);      
      ldouble myrho=RHO0 / (1.+VPRIME*t);
      ldouble myuint=UINT0 / pow(1.+VPRIME*t,GAMMA);
      ldouble myvx=VPRIME*xx / (1.+VPRIME*t);
      //L1 in rho:
      L1+=fabs(get_u(p,0,i,0,0)-myrho);
    }
  scalars[0]=L1/(ldouble)NX;
#endif //CALCL1_HUBBLE

  return 0;
}



//**********************************************************************
//integrates mass in the domain
//**********************************************************************



ldouble
calc_totalmass()
{
  int ix, iy, iz;
  
  ldouble mass = 0.;
  
#pragma omp parallel for private(ix,iy,iz) reduction(+:mass)
  for(iz = 0; iz < NZ; iz++)
  {
    ldouble xx[4], dx[3], rho, gdet;

    for(iy = 0; iy < NY; iy++)
    {
      for(ix = 0; ix < NX; ix++)
      {
        if(doingavg)
        {
          struct geometry geomBL;
          fill_geometry_arb(ix, iy, iz, &geomBL, OUTCOORDS);
          
          rho = get_uavg(pavg, RHO, ix, iy, iz);
          gdet = geomBL.gdet;

	  //cell dimensions
    	  //ANDREW put cell size code in a function with precompute option
          get_cellsize_out(ix, iy, iz, dx);
	  /*
	  ldouble xx1[4], xx2[4];
          xx1[0] = 0.; xx1[1] = get_xb(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_x(iz, 2);
          xx2[0] = 0.; xx2[1] = get_xb(ix+1, 0);xx2[2] = get_x(iy, 1); xx2[3] = get_x(iz, 2);
          coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
          coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
          dx[0] = fabs(xx2[1] -xx1[1]);

          xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_xb(iy, 1); xx1[3] = get_x(iz, 2);
          xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_xb(iy+1, 1); xx2[3] = get_x(iz, 2);
          coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
          coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
          dx[1] = fabs(xx2[2] - xx1[2]);

	  xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_xb(iz, 2);
          xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_x(iy, 1); xx2[3] = get_xb(iz+1, 2);
          coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
          coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
          dx[2] = fabs(xx2[3] - xx1[3]);
          */
	  if(NZ==1) 
          {
            dx[2] = 2. * M_PI;
          }
          else
          {
            #ifdef PHIWEDGE
            dx[2] *= (2. * M_PI / PHIWEDGE);
            #else

	    #endif
          }
        }
        else 
        {
          get_xx(ix, iy, iz, xx);
          dx[0] = get_size_x(ix, 0);
          dx[1] = get_size_x(iy, 1);
          dx[2] = get_size_x(iz, 2);
          #ifdef PHIWEDGE
          if(NZ > 1) dx[2]*=(2. * M_PI / PHIWEDGE);
          #endif
          gdet = calc_gdet(xx);
          rho = get_u(p, 0, ix, iy, iz);
        }
        
        mass += rho * dx[0] * dx[1] * dx[2] * gdet;
      }
    }
  }

  return mass;
}


//**********************************************************************
//calculates the Eddington mass accretion rate
//**********************************************************************

ldouble
calc_mdotEdd()
{
#if (PROBLEM==7) //spherical bondi
  ldouble mcgs=2.23/16.*1e18*MASS; //g/s
#else
  ldouble mcgs=1.09649*2.23e18*MASS*(0.057/etaNT); //g/s \propto 1/etaNT(a)
#endif

  return mcgs;
}


//**********************************************************************
//calculates the Eddington luminosity
//*********************************************************************

ldouble
calc_lumEdd()
{
  ldouble Lcgs=1.25e38*MASS; //erg/s

  return Lcgs;
}


//**********************************************************************
//calculates local radial fluxes of energy
//normalized to total sphere, taken at radius radlum
//**********************************************************************

int
calc_local_lum(int ix,int iy,int iz,ldouble *radlum, ldouble *totallum)
{
  int iv,i,j;
  ldouble xx[4],dx[3],pp[NV],Rrt,rhour,Tij[4][4],Trt;
  ldouble Rij[4][4],Rtt,ehat,ucongas[4];
  ldouble tautot[3],tau=0.;
  ldouble gdet;
  double lum,jet;

  for(iv=0;iv<NV;iv++)
    pp[iv]=get_u(p,iv,ix,iy,iz);

  struct geometry geomBL;
  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  if(doingavg)
    {
      PLOOP(iv)
	pp[iv]=get_uavg(pavg,iv,ix,iy,iz);

      ldouble ucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
      ldouble uconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);		  
      rhour=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
      Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
	+ GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
	+ get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
	- get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

#ifdef RADIATION
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
      Rrt=Rij[1][0];
#else
      Rrt=0.;
#endif

      lum=-geomBL.gdet*Rrt*4.*M_PI;
      jet=geomBL.gdet*(Trt+rhour+Rrt)*4.*M_PI;
    }
  else
    {
	  
      ucongas[1]=pp[2];
      ucongas[2]=pp[3];
      ucongas[3]=pp[4];	      
      conv_vels(ucongas,ucongas,VELPRIM,VEL4,geom.gg,geom.GG);

      rhour = pp[RHO]*ucongas[1];
	  
      calc_Tij(pp,&geom,Tij);
      indices_2221(Tij,Tij,geom.gg);
      Trt=Tij[1][0];


#ifdef RADIATION	      
      calc_Rij(pp,&geom,Rij); 
      indices_2221(Rij,Rij,geom.gg);
      Rrt=Rij[1][0];
#endif
     

      lum=-geom.gdet*Rrt*4.*M_PI;
      jet=geom.gdet*(rhour+Trt+Rrt)*4.*M_PI;

    }

  *radlum=lum;
  *totallum=jet;

  return 0;
}


//**********************************************************************
//calculates luminosity by integrating positive flux from the axis up to tau=1 surface
//normalized to total sphere, taken at radius radius
//type==0: R^r_t outside photosphere
//type==1: sum of positive R^r_t everywhere
//type==2: R^r_t in the outflow region
//type==3: sum of R^r_t everywhere
//**********************************************************************

int
calc_lum(ldouble radius,int type,ldouble *radlum, ldouble *totallum)
{

  int ix,iy,iz;
  ldouble xxC[4],xxBL[4];
 
  //search for appropriate radial index
  for(ix=0;ix<NX-1;ix++)
    {
      #ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix, 0, 0, xxBL);
      #else
      get_xx(ix,0,0,xxC);
      coco_N(xxC,xxBL,MYCOORDS,OUTCOORDS);
      #endif

      if(xxBL[1]>radius) break;
    }
  if(ix==NX) 
    {
      ix=NX-2;
      #ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix, 0, 0, xxBL);
      #else
      get_xx(ix,0,0,xxC);      
      coco_N(xxC,xxBL,MYCOORDS,OUTCOORDS);
      #endif
    }
      

  if(NY==1 && NZ==1) //spherical symmetry
  {

      iz=0; 
      iy=0;
      ldouble dx[3],pp[NV],Rrt,rhour,uintur,Tij[4][4],Trt;
      ldouble Rij[4][4],Rtt,ehat,ucongas[4];
      ldouble tautot[3],tau=0.;
      ldouble gdet;
      ldouble lum,jet;

      int iv;
      for(iv=0;iv<NV;iv++)
	pp[iv]=get_u(p,iv,ix,iy,iz);

      struct geometry geomBL;
      fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);

      if(doingavg)
      {
	  PLOOP(iv)
	    pp[iv]=get_uavg(pavg,iv,ix,iy,iz);

	  ldouble ucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	  ldouble uconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);		  
	  rhour=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
	  uintur=get_uavg(pavg,AVGUUUCON(1),ix,iy,iz);
	  Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
	    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
	    + get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
	    - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

	  Rrt=0.;
#ifdef RADIATION

	  int i,j;
	  if(type==0) //R^r_t outside photosphere
	    {
	      Rrt=0.;
	    }
	  else if(type==1) //R^r_t everywhere
	    {
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
	      Rrt=Rij[1][0];// + ehat*uconr);
	      if(Rrt<0.) Rrt=0.;
	    }
	  else if(type==2) //R^r_t everywhere in outflow
	    {
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
	      Rrt=Rij[1][0];// + ehat*uconr);
	      if(uconr<0. || Rrt<0.) Rrt=0.;
	    }
	  else
	    Rrt=0.;

	  lum=-geom.gdet*Rrt*4.*M_PI;
#else //RADIATION
	  lum=-geom.gdet*uintur*4.*M_PI;    
#endif
	  jet=geomBL.gdet*(Trt+rhour+Rrt)*4.*M_PI;
      }
      else //doingavg
      {
	
	  ucongas[1]=pp[2];
	  ucongas[2]=pp[3];
	  ucongas[3]=pp[4];	      
	  conv_vels(ucongas,ucongas,VELPRIM,VEL4,geom.gg,geom.GG);

	  rhour = pp[RHO]*ucongas[1];
	  uintur = pp[UU]*ucongas[1];
	  
	  calc_Tij(pp,&geom,Tij);
	  indices_2221(Tij,Tij,geom.gg);
	  Trt=Tij[1][0];

	  Rrt=0.;

#ifdef RADIATION	      
	  if(type==0) //R^r_t outside photosphere
	    {
	      Rrt=0.;
	    }
	  else if(type==1) //sum of positive R^r_t everywhere
	    {
	      calc_Rij(pp,&geom,Rij); 
	      Rrt=Rij[1][0];
	      if(Rrt<0.)
		Rrt=0.;
	    }
	  else if(type==2) //R^r_t in the outflow region
	    {
	      calc_Rij(pp,&geom,Rij); 
	      Rrt=Rij[1][0];
	      if(Rrt<0. || ucongas[1]<0.)
		Rrt=0.;
	    }
	  else if(type==3) //sum of R^r_t everywhere
	    {
	      calc_Rij(pp,&geom,Rij); 
	      Rrt=Rij[1][0];
	    }
	  else
	    Rrt=0.;
	  lum=-geom.gdet*Rrt*4.*M_PI;
#else //RADIATION 
	  lum=-geom.gdet*uintur*4.*M_PI;    
#endif
	  jet=geom.gdet*(rhour+Trt+Rrt)*4.*M_PI;
      }

      *radlum=lum;
      *totallum=jet;
      return 0.;
  }
  else //non-sph symmetry
  {

      ldouble lum=0.;
      ldouble jet=0.;
      #pragma omp parallel for private(iy,iz) reduction(+:lum) reduction(+:jet)
      for(iz=0;iz<NZ;iz++)
      {
        for(iy=0;iy<NY;iy++)
	{
	  ldouble xx[4],dx[3],pp[NV],Rrt,rhour,uintur,Tij[4][4],Trt;
          ldouble Rij[4][4],Rtt,ehat,ucongas[4],ucovgas[4];
          ldouble tautot[3],tau=0.;
          //ldouble gdet;

	  int iv;
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);

	  //get_xx(ix,iy,iz,xx);
	  //dx[0]=get_size_x(ix,0);
	  //dx[1]=get_size_x(iy,1);
	  //dx[2]=2.*M_PI;
	  //gdet=geom.gdet;

	  ldouble dxph[3],dxBL[3];
	  
	  //cell dimensions
	  //ANDREW put cell size code in a function with precompute option
          get_cellsize_out(ix, iy, iz, dxBL);
	  /*
	  ldouble xx1[4],xx2[4];
	  xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_x(iz,2);
	  xx2[0]=0.;xx2[1]=get_xb(ix+1,0);xx2[2]=get_x(iy,1);xx2[3]=get_x(iz,2);
	  coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
	  coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
	  dxBL[0]=fabs(xx2[1]-xx1[1]);
	  xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_x(iz,2);
	  xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_xb(iy+1,1);xx2[3]=get_x(iz,2);
	  coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
	  coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
	  dxBL[1]=fabs(xx2[2]-xx1[2]);
	  xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_xb(iz,2);
	  xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_x(iy,1);xx2[3]=get_xb(iz+1,2);
	  coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
	  coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
	  dxBL[2]=fabs(xx2[3]-xx1[3]);
          */
	  if(NZ==1) 
          {
            dxBL[2]=2.*M_PI;
          }
          else
          {
            #ifdef PHIWEDGE
            dxBL[2] *= (2. * M_PI / PHIWEDGE);
            #endif
          }
	  dxph[0]=dxBL[0]*sqrt(geomBL.gg[1][1]);
	  dxph[1]=dxBL[1]*sqrt(geomBL.gg[2][2]);
	  dxph[2]=dxBL[2]*sqrt(geomBL.gg[3][3]);
	  
	  if(doingavg)
	  {
	      PLOOP(iv)
		pp[iv]=get_uavg(pavg,iv,ix,iy,iz);

	      calc_tautot(pp,&geomBL,dxph,tautot);

	      ldouble ucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      ldouble uconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);		  

	      rhour=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
	      uintur=get_uavg(pavg,AVGUUUCON(1),ix,iy,iz);
	      
	      Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
		+ GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
		+ get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		- get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

	      tau+=ucont*tautot[1];

	      Rrt=0.;
#ifdef RADIATION
	      int i,j;
	      if(type==0) //R^r_t outside photosphere
		{
		  if(tau>1.) break;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  Rrt=Rij[1][0];
		  if(Rrt<0.) Rrt=0.;
		}
	      else if(type==1) //positive R^r_t everywhere
		{
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
		  Rrt=-Rij[1][0];
		  if(Rrt<0.) Rrt=0.;

		  
		}
	      else if(type==2) //R^r_t everywhere in outflow
		{
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
		  Rrt=-Rij[1][0];
		  if(uconr<0. || Rrt<0.) Rrt=0.;
		}
	      else if(type==3) //any R^r_t everywhere
		{
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
		  Rrt=-Rij[1][0];
		  
		}
	      else
		Rrt=0.;
	      
	      lum+=geomBL.gdet*Rrt*dxBL[1]*dxBL[2];
#else //RADIATION 
	      lum+=geomBL.gdet*uintur*dxBL[1]*dxBL[2];
#endif

	      jet+=geomBL.gdet*(rhour+Rrt+Rrt)*dxBL[1]*dxBL[2];
	  }
	  else //snapshot
	  { 
	      
	      //to BL
	      #ifdef PRECOMPUTE_MY2OUT
              trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
              #else      
              trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomBL);
              #endif

	      //hydro part may be inconsistent!
	      calc_tautot(pp,&geomBL,dxph,tautot);

	      ucongas[1]=pp[2];
	      ucongas[2]=pp[3];
	      ucongas[3]=pp[4];	      
	      conv_vels(ucongas,ucongas,VELPRIM,VEL4,geomBL.gg,geomBL.GG);

	      indices_21(ucongas,ucovgas,geomBL.gg);

	      rhour = pp[RHO]*ucongas[1];
	      uintur = pp[UU]*ucongas[1];
	  
	      calc_Tij(pp,&geomBL,Tij);
	      indices_2221(Tij,Tij,geomBL.gg);
	      Trt=Tij[1][0];

	      tau+=ucongas[0]*tautot[1];

	      Rrt=0.;
#ifdef RADIATION
	      if(type==0) //R^r_t outside photosphere
		{
		  if(tau>1.) break;	  
		  calc_Rij(pp,&geomBL,Rij); 
		  indices_2221(Rij,Rij,geomBL.gg);
		  Rrt=-Rij[1][0];
		  if(Rrt<0.) Rrt=0.;
		}
	      else if(type==1) //R^r_t everywhere
		{
		  calc_Rij(pp,&geomBL,Rij); 
		  indices_2221(Rij,Rij,geomBL.gg);
		  Rrt=-Rij[1][0];
		  if(Rrt<0.) Rrt=0.;
		}
	      else if(type==2) //R^r_t in the outflow region
		{
		  calc_Rij(pp,&geomBL,Rij); 
		  indices_2221(Rij,Rij,geomBL.gg);
		  Rrt=-Rij[1][0];
		  if(Rrt<0. || ucongas[1]<0.)
		    Rrt=0.;
		}
	      else
		Rrt=0.;

	      lum+=geomBL.gdet*Rrt*dxBL[1]*dxBL[2];
#else //RADIATION
	      lum+=geomBL.gdet*uintur*dxBL[1]*dxBL[2];
#endif
	      jet+=geomBL.gdet*(rhour+Trt+Rrt)*dxBL[1]*dxBL[2];

	      //ANDREW -- what is this definition of jet luminosity?
	      //jet+=-geomBL.gdet*rhour*(ucovgas[0]+sqrt(-geomBL.gg[0][0]))*dxBL[1]*dxBL[2];  

	  } //snapshot
	} //iy
      } // iz
      
      *radlum=lum;
      *totallum=jet;
      return 0.;
  } // axisymmetric or not

    return -1;
}


int
calc_lum_tausurface(ldouble taumax,ldouble *radlum)
{

  int ix,iy,iz;
  ldouble xxC[4],xxBL[4];
      
  if(NY>1) //non-sph symmetry
  {

    ldouble lum=0.;
    #pragma omp parallel for private(iy,iz) reduction(+:lum)
    for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
      {
        ldouble tau=0.;
        //search for appropriate radial index
        for(ix=NX-1;ix>-1;ix--)
        {
          #ifdef PRECOMPUTE_MY2OUT
          get_xxout(ix, 0, 0, xxBL);
          #else
          get_xx(ix,0,0,xxC);
          coco_N(xxC,xxBL,MYCOORDS,OUTCOORDS);
          #endif

	  ldouble xx[4],dx[3],pp[NV],Rrt,rhour,uintur,Tij[4][4],Trt;
          ldouble Rij[4][4],Rtt,ehat,ucongas[4],ucovgas[4];
          //ldouble gdet;

	  int iv;
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);

	  ldouble dxph[3],dxBL[3];
	  
	  //cell dimensions
	  //ANDREW put cell size code in a function with precompute option
          get_cellsize_out(ix, iy, iz, dxBL);

	  if(NZ==1) 
          {
            dxBL[2]=2.*M_PI;
          }
          else
          {
            #ifdef PHIWEDGE
            dxBL[2] *= (2. * M_PI / PHIWEDGE);
            #endif
          }
	  dxph[0]=dxBL[0]*sqrt(geomBL.gg[1][1]);
	  dxph[1]=dxBL[1]*sqrt(geomBL.gg[2][2]);
	  dxph[2]=dxBL[2]*sqrt(geomBL.gg[3][3]);
	  
	  if(doingavg)
	  {
	      PLOOP(iv)
		pp[iv]=get_uavg(pavg,iv,ix,iy,iz);

	      ldouble ucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      ldouble uconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);		  

	      rhour=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
	      uintur=get_uavg(pavg,AVGUUUCON(1),ix,iy,iz);
	      
	      Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
		+ GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
		+ get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		- get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

	      tau+=(ucongas[0]-ucongas[1])*calc_kappaes(pp,&geomBL)*dxph[0];

	      Rrt=0.;
#ifdef RADIATION
	      int i,j;
	      if(tau >= taumax) //R^r_t outside photosphere
		{
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  Rrt=-Rij[1][0];
		  if(Rrt<0.) Rrt=0.;
	          lum+=geomBL.gdet*Rrt*dxBL[1]*dxBL[2];
                  break;
		}
	      else
		Rrt=0.;

              if(xxBL[1] < 2.) break;
#else
	      if(tau >= taumax) //R^r_t outside photosphere
		{
		  Rrt=uintur;
		  if(Rrt<0.) Rrt=0.;
	          lum+=geomBL.gdet*Rrt*dxBL[1]*dxBL[2];
                  break;
		}
	      else
		Rrt=0.;

              if(xxBL[1] < 2.) break;
#endif
	  }
	  else //snapshot
	  { 
	      
	      //to BL
	      #ifdef PRECOMPUTE_MY2OUT
              trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
              #else      
              trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomBL);
              #endif

	      ucongas[1]=pp[2];
	      ucongas[2]=pp[3];
	      ucongas[3]=pp[4];	      
	      conv_vels(ucongas,ucongas,VELPRIM,VEL4,geomBL.gg,geomBL.GG);

	      indices_21(ucongas,ucovgas,geomBL.gg);

	      rhour = pp[RHO]*ucongas[1];
	      uintur = pp[UU]*ucongas[1];
	  
	      calc_Tij(pp,&geomBL,Tij);
	      indices_2221(Tij,Tij,geomBL.gg);
	      Trt=Tij[1][0];

	      tau+=(ucongas[0]-ucongas[1])*calc_kappaes(pp,&geomBL)*dxph[0];

	      Rrt=0.;
#ifdef RADIATION
	      if(tau >= taumax) //R^r_t outside photosphere
		{
		  calc_Rij(pp,&geomBL,Rij); 
		  indices_2221(Rij,Rij,geomBL.gg);
		  Rrt=-Rij[1][0];
		  if(Rrt<0.) Rrt=0.;
    	          lum+=geomBL.gdet*Rrt*dxBL[1]*dxBL[2];
                  break;
		}
	      else
		Rrt=0.;

              if(xxBL[1] < 2.) break;
#else
	      if(tau >= taumax) //R^r_t outside photosphere
		{
		  Rrt=uintur;
		  if(Rrt<0.) Rrt=0.;
    	          lum+=geomBL.gdet*Rrt*dxBL[1]*dxBL[2];
                  break;
		}
	      else
		Rrt=0.;

              if(xxBL[1] < 2.) break;
#endif
	  } //snapshot
        } //ix
      } //iy
    } // iz
      
      *radlum=lum;
      return 0.;
  } // axisymmetric or not

    return -1;
}


//**********************************************************************
//calculates luminosity escaping through radial and polar boundaries
//**********************************************************************

ldouble
calc_exitlum()
{
#ifdef RADIATION

  int ix,iy,iz,iv,i,j;
  ldouble xx[4],xxBL[4],dx[3],pp[NV],Rrt;
  ldouble Rij[4][4],Rtt,Rtht;
  ldouble tautot[3],tau=0.;
  ldouble gdet;

  ldouble lum=0.,jet=0.;

  int distfromboundary=5;

  if(NZ==1) //phi-symmetry only
    {
      //inner radial
      iz=0; ix=0+distfromboundary;
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=2.*M_PI;
	  if(!doingavg)
	    {
	      calc_Rij(pp,&geom,Rij); 
	      indices_2221(Rij,Rij,geom.gg);
	      Rrt=-Rij[1][0];
	      if(Rrt>0.) Rrt=0.; //neglect inflowin luminosity
	      lum+=geom.gdet*fabs(Rrt)*dx[1]*dx[2];
	    }
	}
      //outer radial
      iz=0; ix=NX-1-distfromboundary;
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=2.*M_PI;
	  if(!doingavg)
	    {
	      calc_Rij(pp,&geom,Rij); 
	      indices_2221(Rij,Rij,geom.gg);
	      Rrt=-Rij[1][0];
	      if(Rrt<0.) Rrt=0.; //neglect inflowin luminosity
	      lum+=geom.gdet*fabs(Rrt)*dx[1]*dx[2];
	    }
	}

      //upper theta 
      iz=0; iy=0+distfromboundary;
      for(ix=0;ix<NX;ix++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=2.*M_PI;
	  if(!doingavg)
	    {
	      calc_Rij(pp,&geom,Rij); 
	      indices_2221(Rij,Rij,geom.gg);
	      Rtht=-Rij[2][0];
	      if(Rtht>0.) Rtht=0.; //neglect inflowin luminosity
	      lum+=geom.gdet*fabs(Rtht)*dx[0]*dx[2];
	    }
	}

      //lower theta 
      iz=0; iy=NY-1-distfromboundary;
      for(ix=0;ix<NX;ix++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=2.*M_PI;
	  if(!doingavg)
	    {
	      calc_Rij(pp,&geom,Rij); 
	      indices_2221(Rij,Rij,geom.gg);
	      Rtht=-Rij[2][0];
	      if(Rtht<0.) Rtht=0.; //neglect inflowin luminosity
	      lum+=geom.gdet*fabs(Rtht)*dx[0]*dx[2];
	    }
	}
    }
      
  return lum;
#else
  return -1.;
#endif
}


//**********************************************************************
//calculates MRI resolution parameter Q_theta ar rmri
//**********************************************************************

ldouble
calc_resmri(ldouble radius)
{
#ifndef BHDISK_PROBLEMTYPE
  return -1.; //no MRI
#endif

#ifdef MAGNFIELD

  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4];
  ldouble qtheta=0.,sigma=0.;
 
  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      #ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix, 0, 0, xxBL);
      #else
      get_xx(ix,0,0,xx);      
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
      #endif

      if(xxBL[1]>radius) break;
    }

#pragma omp parallel for private(iy,iz) reduction(+:sigma) reduction(+:qtheta)
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      {
	ldouble dx[3];
	dx[1]=get_size_x(iy,1);
	ldouble rho=get_u(p,RHO,ix,iy,iz);
	
	ldouble q1,q2;
	calc_Qthetaphi(ix,iy,iz,&q1,&q2);

        sigma+=rho*dx[1];
	qtheta+=rho*fabs(q1)*dx[1];
	//qphi+=rho*q2*dx[1];
      }

  return qtheta/sigma;
  

#endif
    return -1.;
}


//**********************************************************************
//calculates mean temperature at rmri
//**********************************************************************

ldouble
calc_meantemp(ldouble radius)
{
  
#ifndef BHDISK_PROBLEMTYPE
  return -1.; //no disk no cry
#endif

  int ix, iy, iz;
  ldouble xx[4], xxBL[4];
  ldouble mtemp = 0., sigma = 0.;
 
  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    #ifdef PRECOMPUTE_MY2OUT
    get_xxout(ix, 0, 0, xxBL);
    #else
    get_xx(ix,0,0,xx);      
    coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
    #endif

    if(xxBL[1] > radius) break;
  }

#pragma omp parallel for private(iy,iz) reduction(+:sigma) reduction(+:mtemp)
  for (iz = 0; iz < NZ; iz++)
  {
    for(iy = 5; iy < NY-5; iy++)
    {
      ldouble dx[3];
      dx[1] = get_size_x(iy, 1);
      ldouble rho = get_u(p, RHO, ix, iy, iz);
      ldouble ugas = get_u(p, UU, ix, iy, iz);
      ldouble temp = calc_PEQ_Tfromurho(ugas, rho, ix, iy, iz);
      sigma += rho * dx[1];
      mtemp += rho * temp*dx[1];
    }
  }
  
  return mtemp/sigma;
}


//**********************************************************************
//calculates scale height at radius
//**********************************************************************

ldouble
calc_scaleheight(ldouble radius)
{
#ifndef BHDISK_PROBLEMTYPE
  return -1.; //no disk no cry
#endif


  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3];
  ldouble mtemp=0.,sigma=0.,rho,ugas;
 
  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {

      #ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix, 0, 0, xxBL);
      #else
      get_xx(ix,0,0,xx);      
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
      #endif

      if(xxBL[1]>radius) break;
    }

  return scaleth_otg[ix];
}


//**********************************************************************
//calculates theta corresponding to integrated tau=1 from the axis
//**********************************************************************

ldouble
calc_photloc(int ix)
{
  if(MYCOORDS != OUTCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS
     && MYCOORDS != MKS2COORDS && MYCOORDS != MKS3COORDS && MYCOORDS != JETCOORDS)
    return -1.; //no BH

  ldouble tau=0.,pp[NV],xx[4],xxBL[4],dx[3];

  int iz=0; int iy,iv; 

  if(NZ==1)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  ldouble tautot[3];

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);

	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=get_size_x(iz,2);
	  dx[0]=dx[0]*sqrt(geom.gg[1][1]);
	  dx[1]=dx[1]*sqrt(geom.gg[2][2]);
	  dx[2]=2.*M_PI*sqrt(geom.gg[3][3]);

	  #ifdef PRECOMPUTE_MY2OUT
          get_xxout(ix, iy, iz, xxBL);
          #else
	  get_xx(ix,iy,iz,xx);
          coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
          #endif

	  calc_tautot(pp,&geom,dx,tautot);
	  tau+=tautot[1];
	  if(tau>1.) break;
	}
      return xxBL[2];
    }
  else
    return -1;
}

//**********************************************************************
//calculates rest mass flux through r=radius within range of thetas
//normalized to 2pi in phi
//type == 0 (net)
//type == 1 (inflow only)
//type == 2 (outflow only)
//type == 3 (Be>0 outflow only)
//**********************************************************************

ldouble
calc_mdot(ldouble radius, int type)
{
  int ix, iy, iz;
  ldouble xx[4], xxBL[4];

  //variables for Bernoulli computation on the fly
  ldouble rhoucont,Tij[4][4],Ttt;
  ldouble Rij[4][4],Rtt,Be;
  
  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    #ifdef PRECOMPUTE_MY2OUT
    get_xxout(ix, 0, 0, xxBL);
    #else
    get_xx(ix,0,0,xx);      
    coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
    #endif

    if(xxBL[1] > radius) break;
  }

  ldouble mdot=0.;
#pragma omp parallel for private(iy,iz) reduction(+:mdot)
  for(iz = 0; iz < NZ; iz++)
  {
    int iv;
    ldouble dx[3], gdet, rho, rhouconr, ucon[4], pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5];
    for(iy = 0; iy < NY; iy++)
    {
      struct geometry geom;
      fill_geometry_arb(ix, iy, iz, &geom, MYCOORDS);
      
      struct geometry geomBL;
      fill_geometry_arb(ix, iy, iz, &geomBL, OUTCOORDS);
      
      if(doingavg)
      {
        rhouconr = get_uavg(pavg, AVGRHOUCON(1), ix, iy, iz);
        gdet = geomBL.gdet;

	//cell dimensions
    	//ANDREW put cell size code in a function with precompute option
        get_cellsize_out(ix, iy, iz, dx);
	/*  
	ldouble xx1[4], xx2[4];
        xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_xb(iy, 1); xx1[3] = get_x(iz, 2);
        xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_xb(iy+1, 1); xx2[3] = get_x(iz, 2);
        coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        dx[1] = fabs(xx2[2] - xx1[2]);
        xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_xb(iz, 2);
        xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_x(iy, 1); xx2[3] = get_xb(iz+1, 2);
        coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        dx[2] = fabs(xx2[3] - xx1[3]);
        */
	if(NZ==1) 
        {
          dx[2] = 2. * M_PI;
        }
        else
        {
          #ifdef PHIWEDGE
          dx[2] *= (2. * M_PI / PHIWEDGE);
          #endif
        }
      }
      else
      {
        for(iv = 0; iv < NVMHD; iv++)
          pp[iv] = get_u(p, iv, ix, iy, iz);
        
        dx[0] = get_size_x(ix, 0);
        dx[1] = get_size_x(iy, 1);
        dx[2] = get_size_x(iz, 2);
        if(NZ==1)
        {
          dx[2] = 2. * M_PI;
        }
        else
        {
          #ifdef PHIWEDGE
          dx[2] *= (2. * M_PI / PHIWEDGE);
          #endif
        }
        pick_g(ix, iy, iz, gg);
        pick_G(ix, iy, iz, GG);

        rho = pp[0];
        ucon[1] = pp[2];
        ucon[2] = pp[3];
        ucon[3] = pp[4];
        
        conv_vels(ucon, ucon, VELPRIM, VEL4, geom.gg, geom.GG);
        rhouconr = rho * ucon[1];
        gdet = geom.gdet;

        //compute bernoulli 
        Be=0.;
        Ttt=0.;
        Rtt=0.;

	calc_Tij(pp,&geom,Tij);
	indices_2221(Tij,Tij,geom.gg);
	Ttt=Tij[0][0];

        rhoucont= rho*ucon[0];
        Be += -(Ttt+rhoucont)/rhoucont;
        
        #ifdef RADIATION
	calc_Rij(pp,&geom,Rij);
	indices_2221(Rij,Rij,geom.gg);
	Rtt=Rij[0][0];

        Be += -Rtt/rhoucont;
        #endif
        //printf("%e \n",Be);
      }
      
      if(NY==1)
      {
        dx[1] = 2.;
        dx[2] = 2. * M_PI;
      }
      
      if(type==0 || (type==1 && rhouconr<0.) || (type==2 && rhouconr>0.) || (type==3 && Be > 0.) )
        mdot += gdet * rhouconr * dx[1] * dx[2];
    }
  }

  return mdot;
}


//**********************************************************************
//calculates the luminosity proxy for the EHT code comparison project
//**********************************************************************

ldouble
calc_lum_proxy(ldouble radius, ldouble theta_min, ldouble theta_max)
{
  int ix, iy, iz;

  ldouble xx[4], xxBL[4];
  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    #ifdef PRECOMPUTE_MY2OUT
    get_xxout(ix, 0, 0, xxBL);
    #else
    get_xx(ix,0,0,xx);      
    coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
    #endif

    if(xxBL[1] > radius) break;
  }
  
  int ixmin = ix;
  ldouble luminosity = 0.;

#pragma omp parallel for private(ix,iy,iz) reduction(+: luminosity)
  for(iz = 0; iz < NZ; iz++)
  {
    ldouble Edot, dx[3], gdet, pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5];
    ldouble rho, uu, pressure, bsq, bfield, ucon[4], ucov[4], bcon[4], bcov[4];

    for(iy = 0; iy < NY; iy++)
    {
      for(ix = ixmin; ix < NX; ix++)
      {
	#ifdef PRECOMPUTE_MY2OUT
        get_xxout(ix, iy, iz, xxBL);
        #else
        get_xx(ix,iy,iz,xx);      
        coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
        #endif
      
        // Check that theta lies within the required range
        if (xxBL[2] >= theta_min && xxBL[2] <= theta_max)
        {
          struct geometry geom;
          fill_geometry_arb(ix, iy, iz, &geom, MYCOORDS);

	  int iv;
          for(iv = 0; iv < NVMHD; iv++)
            pp[iv] = get_u(p, iv, ix, iy, iz);
          
          dx[0] = get_size_x(ix, 0);
          dx[1] = get_size_x(iy, 1);
          dx[2] = get_size_x(iz, 2);
          if(NZ==1)
          {
            dx[2] = 2. * M_PI;
          }
          else
          {
            #ifdef PHIWEDGE
            dx[2] *= (2. * M_PI / PHIWEDGE);
            #endif
          }
          pick_g(ix, iy, iz, gg);
          pick_G(ix, iy, iz, GG);
          
          rho = pp[0];
          uu = pp[UU];
          pressure = uu * (GAMMA - 1.);
          
          ucon[1] = pp[2];
          ucon[2] = pp[3];
          ucon[3] = pp[4];
          conv_vels(ucon, ucon, VELPRIM, VEL4, geom.gg, geom.GG);
          indices_21(ucon, ucov, geom.gg);
          
          calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, &geom, bcon, bcov, &bsq);
          bfield = sqrt(bsq);
          
          gdet = geom.gdet;
          
          // This is the luminosity proxy
          luminosity += gdet * dx[0] * dx[1] * dx[2] * rho * rho * rho * exp(-0.2 * pow(rho * rho / (bfield * pressure * pressure), (1./3.))) / (pressure * pressure);
        }
      }
    }
  }
  
  return luminosity;
}


//**********************************************************************
//calculates energy flux through r=radius:
//
//   Edot = - Integral_over_theta_phi[(rho + u + p + bsq) u^r u_t - b^r b_t]
//
//**********************************************************************

ldouble
calc_Edot(ldouble radius)
{
  int ix, iy, iz;
  ldouble xx[4], xxBL[4];
  
  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    #ifdef PRECOMPUTE_MY2OUT
    get_xxout(ix, 0, 0, xxBL);
    #else
    get_xx(ix,0,0,xx);      
    coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
    #endif

    if(xxBL[1] > radius) break;
  }
  
  ldouble Edot = 0.;

#pragma omp parallel for private(iy,iz) reduction(+:Edot)
  for(iz = 0; iz < NZ; iz++)
  {
    ldouble dx[3], gdet, rhouconrucovt, uuuconrucovt, bsquconrucovt, bconrbcovt;
    ldouble pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5], ucon[4], ucov[4], bcon[4], bcov[4], rho, uu, bsq;
    for(iy = 0; iy < NY; iy++)
    {
      struct geometry geom;
      fill_geometry_arb(ix, iy, iz, &geom, MYCOORDS);
      
      struct geometry geomBL;
      fill_geometry_arb(ix, iy, iz, &geomBL, OUTCOORDS);
      
      if(doingavg)
      {
        rhouconrucovt = get_uavg(pavg, AVGRHOUCONUCOV(1,0), ix, iy, iz);
        uuuconrucovt = get_uavg(pavg, AVGUUUCONUCOV(1,0), ix, iy, iz);
        bsquconrucovt = get_uavg(pavg, AVGBSQUCONUCOV(1,0), ix, iy, iz);
        bconrbcovt = get_uavg(pavg, AVGBCONBCOV(1,0), ix, iy, iz);
        gdet = geomBL.gdet;

	//cell dimensions
    	//ANDREW put cell size code in a function with precompute option
        get_cellsize_out(ix, iy, iz, dx);
        /*
	ldouble xx1[4], xx2[4];
        xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_xb(iy, 1); xx1[3] = get_x(iz, 2);
        xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_xb(iy+1, 1); xx2[3] = get_x(iz, 2);
        coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        dx[1] = fabs(xx2[2] - xx1[2]);
        xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_xb(iz, 2);
        xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_x(iy, 1); xx2[3] = get_xb(iz+1, 2);
        coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        dx[2] = fabs(xx2[3] - xx1[3]);
        */
        if(NZ==1) 
        {  
          dx[2] = 2. * M_PI;
        }
        else
        {
          #ifdef PHIWEDGE
          dx[2] *= (2. * M_PI / PHIWEDGE);
          #endif
        }
      }
      else //snapshot
      {
	int iv;
        for(iv = 0; iv < NVMHD; iv++)
          pp[iv] = get_u(p, iv, ix, iy, iz);
        
        dx[0] = get_size_x(ix, 0);
        dx[1] = get_size_x(iy, 1);
        dx[2] = get_size_x(iz, 2);
        if(NZ==1)
        {
          dx[2] = 2. * M_PI;
        }
        else
        { 
          #ifdef PHIWEDGE
          dx[2] *= (2. * M_PI / PHIWEDGE);
          #endif
        }
        pick_g(ix, iy, iz, gg);
        pick_G(ix, iy, iz, GG);
        
        rho = pp[0];
        uu = pp[UU];
        
        ucon[1] = pp[2];
        ucon[2] = pp[3];
        ucon[3] = pp[4];
        conv_vels(ucon, ucon, VELPRIM, VEL4, geom.gg, geom.GG);
        indices_21(ucon, ucov, geom.gg);
        
        calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, &geom, bcon, bcov, &bsq);
        
        rhouconrucovt = rho * ucon[1] * ucov[0];
        uuuconrucovt = uu * ucon[1] * ucov[0];
        bsquconrucovt = bsq * ucon[1] * ucov[0];
        bconrbcovt = bcon[1] * bcov[0];
        
        gdet = geom.gdet;
      }
      
        Edot += -gdet * dx[1] * dx[2] * (rhouconrucovt + GAMMA * uuuconrucovt + bsquconrucovt - bconrbcovt);
    }
  }
  
  return Edot;
}


//**********************************************************************
//calculates angular momentum flux through r=radius:
//
//   Ldot = Integral_over_theta_phi[(rho + u + p + bsq) u^r u_phi - b^r b_phi]
//
//**********************************************************************

ldouble
calc_Ldot(ldouble radius)
{
  int ix, iy, iz;
  ldouble xx[4], xxBL[4];
  
  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    #ifdef PRECOMPUTE_MY2OUT
    get_xxout(ix, 0, 0, xxBL);
    #else
    get_xx(ix,0,0,xx);      
    coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
    #endif

    if(xxBL[1] > radius) break;
  }
  
  ldouble Ldot = 0.;

#pragma omp parallel for private(iy,iz) reduction(+:Ldot)
  for(iz = 0; iz < NZ; iz++)
  {
    ldouble dx[3], gdet, rhouconrucovphi, uuuconrucovphi, bsquconrucovphi, bconrbcovphi;
    ldouble pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5], ucon[4], ucov[4], bcon[4], bcov[4], rho, uu, bsq;

    for(iy = 0; iy < NY; iy++)
    {
      struct geometry geom;
      fill_geometry_arb(ix, iy, iz, &geom, MYCOORDS);
      
      struct geometry geomBL;
      fill_geometry_arb(ix, iy, iz, &geomBL, OUTCOORDS);
      
      if(doingavg)
      {
        rhouconrucovphi = get_uavg(pavg, AVGRHOUCONUCOV(1,3), ix, iy, iz);
        uuuconrucovphi = get_uavg(pavg, AVGUUUCONUCOV(1,3), ix, iy, iz);
        bsquconrucovphi = get_uavg(pavg, AVGBSQUCONUCOV(1,3), ix, iy, iz);
        bconrbcovphi = get_uavg(pavg, AVGBCONBCOV(1,3), ix, iy, iz);
        gdet = geomBL.gdet;
	
	//cell dimensions
    	//ANDREW put cell size code in a function with precompute option
        get_cellsize_out(ix, iy, iz, dx);
	/*
	ldouble xx1[4], xx2[4];
        xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_xb(iy, 1); xx1[3] = get_x(iz, 2);
        xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_xb(iy+1, 1); xx2[3] = get_x(iz, 2);
        coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        dx[1] = fabs(xx2[2] - xx1[2]);
        xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_xb(iz, 2);
        xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_x(iy, 1); xx2[3] = get_xb(iz+1, 2);
        coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        dx[2] = fabs(xx2[3] - xx1[3]);
        */
        if(NZ==1) 
        {
          dx[2] = 2. * M_PI;
        }
        else
        {
          #ifdef PHIWEDGE
          dx[2] *= (2. * M_PI / PHIWEDGE);
          #endif
        }
      }
      else
      {
	int iv;
        for(iv = 0; iv < NVMHD; iv++)
          pp[iv] = get_u(p, iv, ix, iy, iz);
        
        dx[0] = get_size_x(ix, 0);
        dx[1] = get_size_x(iy, 1);
        dx[2] = get_size_x(iz, 2);
        if(NZ==1)
        {
          dx[2] = 2. * M_PI;
        }
        else
        {
          #ifdef PHIWEDGE
          dx[2] *= (2. * M_PI / PHIWEDGE);
          #endif
        }
        pick_g(ix, iy, iz, gg);
        pick_G(ix, iy, iz, GG);
        
        rho = pp[0];
        uu = pp[UU];
        
        ucon[1] = pp[2];
        ucon[2] = pp[3];
        ucon[3] = pp[4];
        conv_vels(ucon, ucon, VELPRIM, VEL4, geom.gg, geom.GG);
        indices_21(ucon, ucov, geom.gg);
        
        calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, &geom, bcon, bcov, &bsq);
        
        rhouconrucovphi = rho * ucon[1] * ucov[3];
        uuuconrucovphi = uu * ucon[1] * ucov[3];
        bsquconrucovphi = bsq * ucon[1] * ucov[3];
        bconrbcovphi = bcon[1] * bcov[3];
        
        gdet = geom.gdet;
      }
      
      Ldot += gdet * dx[1] * dx[2] * (rhouconrucovphi + GAMMA * uuuconrucovphi + bsquconrucovphi - bconrbcovphi);
    }
  }
  
  return Ldot;
}


//**********************************************************************
//calculates magnetic flux through r=radius within range of thetas
//normalized to 2pi in phi
//type == 0 (default)
//**********************************************************************

int
calc_Bflux(ldouble radius, int type, ldouble *Bflux, ldouble* Bfluxquad)
{
  *Bflux = *Bfluxquad = 0.;
  
  if(MYCOORDS != OUTCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS &&
     MYCOORDS != MKS2COORDS && MYCOORDS != MKS3COORDS && MYCOORDS!=JETCOORDS)
  {
    return -1.; //no BH
  }

#ifdef MAGNFIELD

  int ix, iy, iz;
  ldouble xx[4], xxBL[4];
  ldouble Psi, Psiquad;

  //search for appropriate radial index ix corresponding to the required radius
  for(ix = 0; ix < NX; ix++)
  {
    #ifdef PRECOMPUTE_MY2OUT
    get_xxout(ix, 0, 0, xxBL);
    #else
    get_xx(ix,0,0,xx);      
    coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
    #endif

    if(xxBL[1] > radius) break;
  }

  Psi = 0.;  // dipole flux
  Psiquad = 0.;  // quadrupole flux

  // We have changed the following so that it can handle both NZ=1 and NZ>1
#pragma omp parallel for private(iy,iz) reduction(+:Psi) reduction(+:Psiquad)
  for (iz = 0; iz < NZ; iz++)
  {
    ldouble rho, ucon[4], pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5];
    ldouble dx[3];
    for(iy = 0; iy < NY; iy++)
    {
      struct geometry geom;
      fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);
      
      struct geometry geomBL;
      fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
      
      
      if(doingavg)  // working with averages
      {
        ldouble bcon[4] = {get_uavg(pavg, AVGBCON(0), ix, iy, iz),
			   get_uavg(pavg, AVGBCON(1), ix, iy, iz),
			   get_uavg(pavg, AVGBCON(2), ix, iy, iz),
			   get_uavg(pavg, AVGBCON(3), ix, iy, iz)};
        
        ldouble ucon[4], ucov[4], bcov[4];
        
        ucon[1] = get_uavg(pavg, AVGRHOUCON(1), ix, iy, iz)/get_uavg(pavg, RHO, ix, iy, iz);
        ucon[2] = get_uavg(pavg, AVGRHOUCON(2), ix, iy, iz)/get_uavg(pavg, RHO, ix, iy, iz);
        ucon[3] = get_uavg(pavg, AVGRHOUCON(3), ix, iy, iz)/get_uavg(pavg, RHO, ix, iy, iz);
        conv_vels(ucon, ucon, VEL4, VEL4, geomBL.gg, geomBL.GG);
        indices_21(ucon, ucov, geomBL.gg);
        
        //Normalize b^0 to be orthogonal with u^\mu
        bcon[0] = -dot3nr(bcon, ucov) / ucov[0];
        indices_21(bcon, bcov, geomBL.gg);

        //Normalize b^mu to be equal to B^2
        ldouble bsq = get_uavg(pavg, AVGBSQ, ix, iy, iz);
        ldouble alphanorm = bsq / dotB(bcon, bcov);
        if(alphanorm < 0.)
	  my_err("alpha.lt.0 in b0 norm !!\n");
        int i4;
	for(i4 = 0; i4 < 4; i4++)
        {
          bcon[i4] *= sqrt(alphanorm);
        }
        
        //Bcon[1]
        ldouble Br = bcon[1] * ucon[0] - bcon[0] * ucon[1];

	//cell dimensions
    	//ANDREW put cell size code in a function with precompute option
        get_cellsize_out(ix, iy, iz, dx);
	/*  
        ldouble xx1[4], xx2[4];
        xx1[0] = 0.; xx1[1] = get_xb(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_x(iz, 2);
        xx2[0] = 0.; xx2[1] = get_xb(ix+1, 0); xx2[2] = get_x(iy, 1); xx2[3] = get_x(iz, 2);
        coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        dx[0] = fabs(xx2[1] - xx1[1]);
        xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_xb(iy, 1); xx1[3] = get_x(iz, 2);
        xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_xb(iy+1, 1); xx2[3] = get_x(iz, 2);
        coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        dx[1] = fabs(xx2[2] - xx1[2]);
        xx1[0] = 0.; xx1[1] = get_x(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_xb(iz, 2);
        xx2[0] = 0.; xx2[1] = get_x(ix, 0); xx2[2] = get_x(iy, 1); xx2[3] = get_xb(iz+1, 2);
        coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        dx[2] = fabs(xx2[3] - xx1[3]);
        */
        if(NZ==1)
        {
          dx[2] = 2. * M_PI;
        }
        else
        {
          #ifdef PHIWEDGE
          dx[2] *= (2. * M_PI / PHIWEDGE);
          #endif
        }
        
        if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
          Psi += geomBL.gdet * fabs(Br) * dx[1] * dx[2];
        
        if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
          Psiquad += geomBL.gdet * Br * my_sign(geomBL.yy - M_PI/2.) * dx[1] * dx[2];
      }
      else  // working with snapshot
      {
	int iv;
        for(iv = 0; iv < NVMHD; iv++)
        {
          pp[iv] = get_u(p, iv, ix, iy, iz);
        }
        
        dx[0] = get_size_x(ix, 0);
        dx[1] = get_size_x(iy, 1);
        if(NZ == 1)
        {
          dx[2] = 2.*M_PI;
        }
        else
        {
          dx[2] = get_size_x(iz, 2);

          #ifdef PHIWEDGE
          dx[2] *= (2. * M_PI / PHIWEDGE);
          #endif
        }
        
        struct geometry geom;
        fill_geometry_arb(ix, iy, iz, &geom, MYCOORDS);
        
        struct geometry geomBL;
        fill_geometry_arb(ix, iy, iz, &geomBL, OUTCOORDS);
        
        ldouble Br = pp[B1];
        
        if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
          Psi += geom.gdet * fabs(Br) * dx[1] * dx[2];
        
        if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
          Psiquad += geom.gdet * Br * my_sign(geomBL.yy - M_PI/2.) * dx[1] * dx[2];
      }  // if(doingavg)
    }  // for(iy=0;iy<NY;iy++)
  }  // for (iz = 0; iz < NZ; iz++)

  *Bflux = Psi;  // dipole flux
  *Bfluxquad = Psiquad;  // quadrupole flux
  
  return 0;
  
#else  // no magnetic field
  
  return -1;
  
#endif
}
