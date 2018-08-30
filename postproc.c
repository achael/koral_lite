/*! \file postproc.c
 \brief Routines for postprocessing, used both on the go and separately
 */


#include "ko.h"


/*********************************************/
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
//rho-weighted minus radial velocity in the outflow (36)
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
  ldouble normalize[NX], rho134[NX], pgas134[NX], bfield134[NX], uconphi134[NX], ptot134[NX], betainv134[NX];

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
    ldouble Sigmagdet;
    ldouble diffflux[NV];
    ldouble fd_u0[NV],fd_up1[NV],fd_up2[NV],fd_um1[NV],fd_um2[NV];
    ldouble fd_p0[NV],fd_pp1[NV],fd_pp2[NV],fd_pm1[NV],fd_pm2[NV],fd_pm3[NV],fd_pp3[NV];
    ldouble fd_pl[NV],fd_pr[NV],fd_plm1[NV],fd_prm1[NV],fd_plp1[NV],fd_prp1[NV];
    ldouble fd_ul[NV],fd_ur[NV],fd_ulm1[NV],fd_urm1[NV],fd_ulp1[NV],fd_urp1[NV];
    ldouble du[NV],dul[NV],dur[NV],aaa[24],ahd,arad,Jvec[3];
    ldouble Gi[4],Giff[4];
    int injet;
    
    //vertically integrated/averaged profiles
    
    for(iv=0;iv<NRADPROFILES;iv++)
      profiles[iv][ix]=0.;
    
    Jvec[0]=Jvec[1]=Jvec[2]=0.;
    
    //keep track of extra quantities for problem 134
    if (PROBLEM == 134)
    {
      normalize[ix] = 0.;
      rho134[ix] = 0.;
      pgas134[ix] = 0.;
      bfield134[ix] = 0.;
      uconphi134[ix] = 0.;
      ptot134[ix] = 0.;
      betainv134[ix] = 0.;
    }
    
    //outside horizon?
    struct geometry geomBLtemp;
    fill_geometry_arb(ix,0,0,&geomBLtemp,OUTCOORDS);
    if(geomBLtemp.xx<=1.1*rhorizonBL) continue; //to avoid working inside horizon
    
    Bangle1=Bangle2=0.;
    Sigmagdet=0.;
    ldouble jetsigma=0.;
    
    for(iv=0;iv<NAVGVARS;iv++)
      avgsums[iv][ix]=0.;
    
    tautot=tauabs=0.;
    
    // #ifdef BHDISK_PROBLEMTYPE
    for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
      {
        //metric
        pick_g(ix,iy,iz,gg);
        pick_G(ix,iy,iz,GG);
        
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
        coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
        ldouble dxph[3];
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
        
        dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
        dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
        dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);
        
        //primitives at the cell - either averaged or original, in BL or MYCOORDS
        for(iv=0;iv<NV;iv++)
        {
          pp[iv]=get_u(p,iv,ix,iy,iz);
        }
        
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
        
        //to BL, res-files and primitives in avg in MYCOORDS
        trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,xx,&geom,&geomBL);
        
        //transforming interpolated primitives to BL
        trans_pall_coco(fd_pm1,fd_pm1,MYCOORDS,OUTCOORDS,xx,&geomm1,&geomBLm1);
        trans_pall_coco(fd_pp1,fd_pp1,MYCOORDS,OUTCOORDS,xx,&geomp1,&geomBLp1);
        trans_pall_coco(fd_pl,fd_pl,MYCOORDS,OUTCOORDS,xx,&geoml,&geomBLl);
        trans_pall_coco(fd_pr,fd_pr,MYCOORDS,OUTCOORDS,xx,&geomr,&geomBLr);
        trans_pall_coco(fd_plp1,fd_plp1,MYCOORDS,OUTCOORDS,xx,&geoml,&geomBLl);
        trans_pall_coco(fd_prm1,fd_prm1,MYCOORDS,OUTCOORDS,xx,&geomr,&geomBLr);
        
        ldouble vischeating=0.;
        
        if(doingavg)
        {
          rho=get_uavg(pavg,RHO,ix,iy,iz);
          uint=get_uavg(pavg,UU,ix,iy,iz);
          //temperature
          temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
          //temp=get_uavg(pavg,AVGTGAS,ix,iy,iz);
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
              
            /*
             Tij[i][j]=get_uavg(pavg,AVGRHOUCONUCOV(i,j),ix,iy,iz)
             + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(i,j),ix,iy,iz)
             + get_uavg(pavg,AVGBSQUCONUCOV(i,j),ix,iy,iz)
             + delta(i,j)*(GAMMAM1*get_uavg(pavg,UU,ix,iy,iz)+get_uavg(pavg,AVGBSQ,ix,iy,iz)/2.)
             - get_uavg(pavg,AVGBCONBCOV(i,j),ix,iy,iz);
             */
              
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
          //temperature
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
            
            //test - right face only
            /*
             du[iv]=dur[iv];
             du[iv]/=gdetuBL; //de facto substracting ahd (getd_r (rho ut_rR - rho ut_rL))
             */
          }
          
          //test - right face only
          /*
           double ff1[NV],ff2[NV];
           f_flux_prime(fd_pr,0,ix+1,iy,iz,ff1,0);
           f_flux_prime(fd_plp1,0,ix+1,iy,iz,ff2,0);
           rhouconr=.5*(ff1[0]+ff2[0])/gdetuBL; //de facto plotting gdet_r * rhour_r
           Trt=.5*(ff1[1]-ff1[0]+ff2[1]-ff2[0])/gdetuBL;
           //Trt=.5*(ff1[1]+ff2[1])/gdetuBL;
           */
          
          //wavespeeds
          calc_wavespeeds_lr_pure(pp,&geomBL,aaa);
          ahd=my_max(fabs(aaa[0]),fabs(aaa[1]));
          arad=my_max(fabs(aaa[6]),fabs(aaa[7]));
          
          //test
          /*
           calc_wavespeeds_lr_pure(fd_pm1,&geomBLm1,aaa);
           ahd=my_max(ahd,my_max(fabs(aaa[0]),fabs(aaa[1])));
           arad=my_max(arad,my_max(fabs(aaa[6]),fabs(aaa[7])));
           
           calc_wavespeeds_lr_pure(fd_pp1,&geomBLp1,aaa);
           ahd=my_max(ahd,my_max(fabs(aaa[0]),fabs(aaa[1])));
           arad=my_max(arad,my_max(fabs(aaa[6]),fabs(aaa[7])));
           */
          
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
        }  // end on the go from primitives
        
        ldouble muBe,Be,Benoth;
        muBe=-(Trt+rhouconr)/rhouconr;
        Be=-(TttBe+rhoucont)/rhoucont;
        Benoth=-(TttBenoth+rhoucont)/rhoucont;
#ifdef RADIATION
        muBe+=-Rrt/rhouconr;
        Be+=-Rtt/rhoucont;
#endif
        
        
        int isjet;
        if(muBe>0.05 && (xxBL[2]<M_PI/4. || xxBL[2]>3.*M_PI/4.))
          isjet=1;
        else isjet=0;
        
        
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
        //ldouble ibeta=bsq/2./(GAMMAM1*uint);
        profiles[30][ix]+=rho*ibeta*dxph[1];
        if (PROBLEM==93){
          profiles[57][ix]+=bsq/2.*dxph[1];
          profiles[58][ix]+=prermhd*dxph[1];
        }
        
        //rho-weighted prad/pgas (33)
#ifdef RADIATION
        profiles[31][ix]+=rho*prerad/pregas*dxph[1];
#else
        profiles[31][ix]+=0.;
#endif
        
        //surface density (2) (column)
        profiles[0][ix]+=rho*dxph[1];
        //temporarily total pressure gas+radiation:
        //profiles[0][ix]+=(prermhd)*dxph[1];
        
        //surface energy density (41)
        if(muBe<0.)
          profiles[39][ix]+=endensimple*dxph[1];
        //temporarily magnetic pressure:
        //profiles[39][ix]+=(bsq/2.)*dxph[1];
        
        //numerator of scale height (31) (column)
#ifndef CALCHRONTHEGO
        profiles[29][ix]+=rho*dxph[1]*pow(tan(fabs(M_PI/2.-xxBL[2])),2.);
#endif
        
        //surface density in the inflow (23)
        if(utcon[1]<0.)
          profiles[21][ix]+=rho*dxph[1];
        
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
        
        //temporary surface density to normalize what is above
        Sigmagdet+=rho*dx[1]*dx[2]*geomBL.gdet;
        
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
        
        //total rad energy flux (17)
#ifdef RADIATION

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
        //if(fabs(geomBL.xxvec[2]-M_PI/2)<scaleth_otg[ix])
        if(Be<0.)
          profiles[52][ix]+=my_sign(geomBL.xxvec[2]-M_PI/2.)*endenuconth*dxph[1];
        
        //rho-weighted minus radial velocity in the inflow (24)
        if(utcon[1]<0.)
          profiles[22][ix]+=-rhouconr*dxph[1];
        
        //rho-weighted minus radial velocity in the outflow (36)
        if(utcon[1]>0.)
          profiles[34][ix]+=rhouconr*dxph[1];
        
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
        
        //u_phi evywhree (5)
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
        
        // special quantities for problem 134
        if (PROBLEM == 134)
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
    profiles[34][ix]/=sigmaout;
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
    
    //to get velocities
    //profiles[2][ix]=profiles[1][ix]/Sigmagdet;
    //profiles[22][ix]=profiles[8][ix]/Sigmagdet;
    //profiles[34][ix]=profiles[9][ix]/Sigmagdet;
    
    
    //luminosity at given radius (12)
    ldouble radlum,totallum;
    calc_lum(xxBL[1],0,&radlum,&totallum);
    profiles[10][ix]=radlum;
    //luminosity at given radius (22)
    calc_lum(xxBL[1],1,&radlum,&totallum);
    profiles[20][ix]=radlum;
    //location of the photosphere (13)
    profiles[11][ix]=calc_photloc(ix);
    
    // special quantities for problem 134
    if (PROBLEM == 134)
    {
      profiles[7][ix] = rho134[ix] / normalize[ix];
      profiles[8][ix] = pgas134[ix] / normalize[ix];
      profiles[12][ix] = bfield134[ix] / normalize[ix];
      profiles[13][ix] = uconphi134[ix] / normalize[ix];
      profiles[17][ix] = ptot134[ix] / normalize[ix];
      profiles[18][ix] = betainv134[ix] / normalize[ix];

    }
    
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
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
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
	  trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,xx,&geom,&geomBL);

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
	      //coordinates
	      ldouble dxph[3],dx[0];
	      ldouble xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(iix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(iix+1,0);xx2[2]=get_xb(iy,1);xx2[3]=get_xb(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
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
  scalars[10]: luminosity proxy for problem 134
  scalars[11]: Edot for problem 134
  scalars[12]: Ldot for problem 134
 
 */
/*********************************************/

int calc_scalars(ldouble *scalars,ldouble t)
{
  //adjust NSCALARS in problem.h

  /*********************************************/
  //base for BHDISK problems
  /*********************************************/

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

  ldouble xx[4],xxBL[4];
  get_xx(NX-1,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);

  //luminosities 
  ldouble rlum=15.;
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

  //radiative luminosity everywhere (4)
  scalars[2]=radlum*mdotscale*CCC0*CCC0/calc_lumEdd();

  if(PROBLEM==89 || PROBLEM==79) //RADTORUS or SOFTBALL
  {
    //luminosity exiting through radial and theta boundaries (3)
    scalars[1]=calc_exitlum();
  }

  //total energy flux at infinity (rho ur + Trt + Rrt) (12)
  //calc_lum(rlum,1,&radlum,&totallum);
  scalars[10]=totallum;
  
  //mri resolution parameter Q_theta (7) at rmri
  ldouble rmri=xxBL[1]/2.;  // middle radial cell
#if(PROBLEM==69) //INJDISK
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
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
      if(xxBL[1]>radius) break;
    }

  double totlum;
  calc_local_lum(ix,NCCORRECTPOLAR+1,0,&radlum,&totlum);
  scalars[11]=totlum;
  
#if(PROBLEM==134)  // FISHMONC for EHT code comparison tests
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

  /* accretion rates through the outer edge for TDEMILIO */
  
#if(PROBLEM==91 || PROBLEM==94)
  rmdot = 300.;
  
  //inflow (12)
  mdot=calc_mdot(rmdot,1);
  scalars[10]=-mdot*mdotscale/calc_mdotEdd();
  //outflow (13)
  mdot=calc_mdot(rmdot,2);
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

  //ANDREW should be in not just one cell??  
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

  //uth = calc_thermal_ne(pp)*K_BOLTZ * Te /(calc_gammaintfromtemp(Te,ELECTRONS)-1.);
  //uthi = (K_BOLTZ*pp[RHO]/MU_I/M_PROTON) * Ti /(calc_gammaintfromtemp(Ti,IONS)-1.);
  // uth = calc_thermal_ne(pp)*K_BOLTZ * Te /(4./3.-1.);
  // uthi = (K_BOLTZ*pp[RHO]/MU_I/M_PROTON) * Ti /(5./3.-1.);

  uur=calc_relel_uint(pp);
  ldouble nur=calc_relel_ne(pp);
  ldouble ni=pp[RHO]/M_PROTON/MU_I;
  ldouble neth=calc_thermal_ne(pp);
  ldouble rhoeth=neth*M_PROTON*MU_E;
  
  //uth  = calc_ufromS4n(pp[ENTRE],neth,ELECTRONS,0,0,0);
  //uthi = calc_ufromS4n(pp[ENTRI],ni,IONS,0,0,0);
  uth  = calc_ufromSen(pp[ENTRE],rhoeth,ELECTRONS,0,0,0);
  uthi = calc_ufromSen(pp[ENTRI],pp[RHO],IONS,0,0,0);

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

	  //	  if(i==TNX/2 && j==TNY/2) printf("%d %d : %f %f %f : %f %f : %e %e\n",i,j,gamma,gammae,gammai,(ue+ui)/uint,(pe+pi)/pre,ue,uint);

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
  //printf("%e\n",uurtot/utot);
  //scalars[3]=uspecies/utot;  //(5)
  //scalars[4]=pspecies/ptot;  //(6)
  //scalars[5]=kinentot; //(7)
  scalars[3] = uetot;
  scalars[4] = uitot;
  scalars[5] = uurtot;

  scalars[6]=calc_totalmass(); //(8)

#ifdef RELELECTRONS
  for (ie=0 ; (ie < NRELBIN) && (ie+end < NSCALARS); ie++)
  scalars[end+ie] /= netot;
#endif


#endif


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

#endif

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
#endif

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
#endif

  return 0;
}


/*********************************************/
/* calculates box-related scalars  */
/*********************************************/

int calc_boxcorrscalars(ldouble *boxcorrscalars,ldouble t)
{
#ifdef CONSISTENTGAMMA
  if(PROCID==0) printf("2 gammas not consistent in postprocessing when CONSISTENTGAMMA\n");
#endif

#if(BOXCORROUTPUT==1)
  //adjust NBOXCORRSCALARS in choices.h
  
  int ix,iy,iz,ii,iv,i,j;
  int bix1,bix2,biy1,biy2,giy;
  ldouble xx[4],xxBL[4];
  ldouble pp[NV];

  //zero scalars by default
  for(ii=0;ii<NBOXCORRSCALARS;ii++)
    boxcorrscalars[ii]=0.;

  //search for appropriate indices
  

  //limits of this tile
  int rmin,rmax;
  get_xx(0,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
  rmin=xxBL[1];
  get_xx(NX-1,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
  rmax=xxBL[1];

  //global index of the theta limits
  int giy1,giy2;
  mpi_local2globalidx(0,0,0,&ix,&giy1,&iz);
  mpi_local2globalidx(0,NY-1,0,&ix,&giy2,&iz);

  //if tile completely outside the box
  int ifoutsidebox=0;

  //theta indices of the box
  int ITH1, ITH2;
#ifdef BOXITH //symmetrical wrt eq.plane
  ITH1=TNY/2-BOXITH;
  ITH2=TNY/2+BOXITH-1;
#endif
#ifdef BOXITH1
  ITH1=BOXITH1;
#endif
#ifdef BOXITH2
  ITH2=BOXITH2;
#endif
  

  if((rmax<BOXR1) || (rmin>BOXR2))
    ifoutsidebox=1;

  if((giy2 < ITH1) || (giy1 >ITH2))
    ifoutsidebox=1;
  
  if(!ifoutsidebox)  //do the integrals only if given tile covers some part the box
    {
      //radial limits first 
      bix1=bix2=-1;
      for(ix=0;ix<NX;ix++)
	{
	  get_xx(ix,0,0,xx);
	  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
	  if(xxBL[1]>BOXR1 & bix1<0) bix1=ix;
	  if(xxBL[1]>BOXR2 & bix2<0) bix2=ix;
	}
      if(bix1<0) bix1=NX;
      if(bix2<0) bix2=NX;

      //then polar angle

      biy1=biy2=-1;
      for(iy=0;iy<NY;iy++)
	{
	  mpi_local2globalidx(0,iy,0,&ix,&giy,&iz);
      
	  if(giy>ITH1 && biy1<0) biy1=iy;
	  if(giy>ITH2 && biy2<0) biy2=iy;
	}
      if(biy1<0) biy1=NY;
      if(biy2<0) biy2=NY;
  
   

      for(ix=bix1;ix<bix2;ix++)
	for(iy=biy1;iy<biy2;iy++)
	  for(iz=0;iz<NZ;iz++)
	    {
	      struct geometry geom;
	      fill_geometry(ix,iy,iz,&geom);
	
	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	
	      //coordinate
	      ldouble dx[3];
	      get_cell_size_arb(ix,iy,iz,dx,OUTCOORDS);

	      //primitives at the cell - either averaged or original, in BL or MYCOORDS
	      for(iv=0;iv<NV;iv++)
		pp[iv]=get_u(p,iv,ix,iy,iz);

	      //to BL, res-files and primitives in avg in MYCOORDS
	      trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);
	      
	      //iy+1 cell, to get gradients:
	      struct geometry geom2;
	      fill_geometry(ix,iy+1,iz,&geom2);
	
	      struct geometry geomBL2;
	      fill_geometry_arb(ix,iy+1,iz,&geomBL2,OUTCOORDS);
	      
	      //from now on - working in BL coords
        
	      //primitives and derivatives
	      ldouble rho=pp[RHO];
	      ldouble uint=pp[UU];
	      ldouble Ti,Te,trad;
	      ldouble temp=calc_PEQ_Teifrompp(pp,&Te,&Ti,ix,iy,iz);
	      ldouble bsq=0.;
	      ldouble ucon[4],ucov[4],bcon[4],bcov[4],vel[4];
	      ucon[1]=pp[VX];
	      ucon[2]=pp[VY];
	      ucon[3]=pp[VZ];

	      ldouble brbphi, bfake[4],bsqint,bangle;
#ifdef MAGNFIELD
	      calc_angle_brbphibsq(ix,iy,iz,&brbphi,&bsqint,bfake,bfake);
	      bangle=-brbphi/bsqint;
	      calc_bcon_prim(pp,bcon,&geomBL);
	      indices_21(bcon,bcov,geomBL.gg); 
	      bsq = dotB(bcon,bcov);
#endif

	      conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);

	      ldouble rhouconr=rho*ucon[1];
	      ldouble Tij[4][4],Tij22[4][4],Rij[4][4];
	      calc_Tij(pp,&geomBL,Tij22);
	      indices_2221(Tij22,Tij,geomBL.gg);
	      ldouble Trt = Tij[1][0],Rrt=0.;
	      ldouble Trtmagn = bsq*ucon[1]*ucov[0] - bcon[1]*bcov[0];
	      ldouble Trtkin =  rho*ucon[1]*ucov[0];
	      ldouble enden = Tij[0][0] + rho*ucon[0];
	      ldouble Ehat=0.;
	      ldouble Ehat2=0.;
	      ldouble Gi[4],Giff[4]={0.,0.,0.,0.};
	      ldouble Gic[4],Gicff[4]={0.,0.,0.,0.};

	      
	      ldouble uconr[4];
#ifdef RADIATION
	      calc_Rij(pp,&geomBL,Rij);
	      indices_2221(Rij,Rij,geomBL.gg);
	      Rrt = Rij[1][0];
	      ldouble kappaes=calc_kappaes(pp,&geomBL);
	      ldouble Rtt;
	      calc_ff_Rtt(pp,&Rtt,uconr,&geomBL);
	      Ehat=-Rtt; 	

	      //radiative velocity
	      uconr[1]=pp[FX];
	      uconr[2]=pp[FY];
	      uconr[3]=pp[FZ];
	      conv_vels(uconr,uconr,VELPRIMRAD,VEL4,geomBL.gg,geomBL.GG);
	      trad=calc_LTE_TfromE(Ehat);
#ifdef EVOLVEPHOTONNUMBER
	      trad=calc_ncompt_Thatrad(pp,&geomBL,Ehat);
#endif
#endif

	      ldouble pregas = GAMMAM1*uint;
	      ldouble premag = bsq/2.;
	      ldouble prerad = 0.;
#ifdef RADIATION
	      prerad = Ehat/3.;
#endif
	      ldouble pretot = pregas + premag + prerad;

	      //angular velocity
	      ldouble Omega=ucon[3]/ucon[0];
	      //MRI resolution parameters
	      ldouble qtheta,qphi;
	      calc_Qthetaphi(ix,iy,iz,&qtheta,&qphi);

	      ldouble gamma=GAMMA;
#ifdef CONSISTENTGAMMA
	      gamma=pick_gammagas(ix,iy,iz);
#endif

	      //saving local integrals 
	      boxcorrscalars[0]=rho; // (2nd column) - density
	      boxcorrscalars[1]=ucon[1]; // (3) - radial velocity
	      boxcorrscalars[2]=bsq; // (4) - magn. en density
	      boxcorrscalars[3]=temp; // (5) - tgas
	      boxcorrscalars[4]=Te; // (6) - Te
	      boxcorrscalars[5]=Ti; // (7) - Ti
	      boxcorrscalars[6]=gamma; // (8) - gamma
	      
 
	      
	      if(PROCID==0)
		{
		  //printing boxcorrscalars
		  fprintf(fout_boxcorrscalars,"%e ",t);
		  for(iv=0;iv<NBOXCORRSCALARS;iv++)
		    fprintf(fout_boxcorrscalars,"%e ",boxcorrscalars[iv]);
		  fprintf(fout_boxcorrscalars,"\n");
		  fflush(fout_boxcorrscalars);
		}
	    }

    
   
    } //if(!ifoutsidebox)


  


 
#endif //BOXCORROUTPUT==1
  return 0;
}


/*********************************************/
/* calculates box-related scalars  */
/*********************************************/

int calc_boxscalars(ldouble *boxscalars,ldouble t)
{
#ifdef CONSISTENTGAMMA
  if(PROCID==0) printf("2 gammas not consistent in postprocessing when CONSISTENTGAMMA\n");
  #endif

#if(BOXOUTPUT==1)
  //adjust NBOXSCALARS in problem.h
  
  int ix,iy,iz,ii,iv,i,j;
  int bix1,bix2,biy1,biy2,giy;
  ldouble xx[4],xxBL[4];

  //zero scalars by default
  for(ii=0;ii<NBOXSCALARS;ii++)
    boxscalars[ii]=0.;

  //search for appropriate indices
  

  //limits of this tile
  int rmin,rmax;
  get_xx(0,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
  rmin=xxBL[1];
  get_xx(NX-1,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
  rmax=xxBL[1];

  //global index of the theta limits
  int giy1,giy2;
  mpi_local2globalidx(0,0,0,&ix,&giy1,&iz);
  mpi_local2globalidx(0,NY-1,0,&ix,&giy2,&iz);

  //if tile completely outside the box
  int ifoutsidebox=0;

  //theta indices of the box
  int ITH1, ITH2;
#ifdef BOXITH //symmetrical wrt eq.plane
  ITH1=TNY/2-BOXITH;
  ITH2=TNY/2+BOXITH-1;
#endif
#ifdef BOXITH1
  ITH1=BOXITH1;
#endif
#ifdef BOXITH2
  ITH2=BOXITH2;
#endif
  

  if((rmax<BOXR1) || (rmin>BOXR2))
    ifoutsidebox=1;

  if((giy2 < ITH1) || (giy1 >ITH2))
    ifoutsidebox=1;
  
  if(!ifoutsidebox)  //do the integrals only if given tile covers some part the box
    {
      //radial limits first 
      bix1=bix2=-1;
      for(ix=0;ix<NX;ix++)
	{
	  get_xx(ix,0,0,xx);
	  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
	  if(xxBL[1]>BOXR1 & bix1<0) bix1=ix;
	  if(xxBL[1]>BOXR2 & bix2<0) bix2=ix;
	}
      if(bix1<0) bix1=NX;
      if(bix2<0) bix2=NX;

      //then polar angle

      biy1=biy2=-1;
      for(iy=0;iy<NY;iy++)
	{
	  mpi_local2globalidx(0,iy,0,&ix,&giy,&iz);
      
	  if(giy>ITH1 && biy1<0) biy1=iy;
	  if(giy>ITH2 && biy2<0) biy2=iy;
	}
      if(biy1<0) biy1=NY;
      if(biy2<0) biy2=NY;
  
      //printf("PROCID: %d > bix: %d - %d > biy: %d - %d -> giy: %d - %d\n",PROCID,bix1,bix2,biy1,biy2,biy1+TOJ,biy2+TOJ);

      //first integrals / averages within the box
      ldouble mass,pgasint,pradint,ptotint,Gtint,Gctint,volume,Fdiffth,Fadvth,Fradth,bangle,banglep,scale,Trphiint,QplusintTrphi,Qplusint;
      mass=pgasint=pradint=ptotint=Gtint=Gctint=Qplusint=QplusintTrphi=Trphiint=0.;
      bangle=banglep=0.;
      ldouble brsq,bthsq,bphsq,brsqav,bthsqav,bphsqav,br,brav;
      brsq=bthsq=bphsq=br=0.;
      brsqav=bthsqav=bphsqav=brav=0.;
      ldouble tempav,qthetaav,alphaav,uradrav,uradthav,tradav,urav,uthav,uphav,Fdiffthav,Fadvthav,Fradthav,bangleav,banglepav,scaleav;
      tempav=tradav=qthetaav=alphaav=urav=uthav=uphav=uradthav=uradrav=volume=Fdiffthav=Fadvthav=Fradthav=bangleav=banglepav=scaleav=0.;
      ldouble rhodev=0.,vphidev=0.;
      ldouble pp[NV];

       ldouble Rrttot[3]={0.,0.,0.};
      ldouble Trttot[3]={0.,0.,0.};
      ldouble Trtkintot[3]={0.,0.,0.};
      ldouble Trtmagntot[3]={0.,0.,0.};
      ldouble rhouconrtot[3]={0.,0.,0.};
      ldouble Ehatuconrtot[3]={0.,0.,0.};
      ldouble areas[3]={0.,0.,0.};

      #ifdef BOXSCALARJUSTONECELL
      bix2=bix1+1;
      biy2=biy1+1;
      #endif

      /*
      int ix,iy,iz;
      mpi_global2localidx(10,32,0,&ix,&iy,&iz);
      if(if_indomain(ix,iy,iz))
	print_primitives(&get_u(p,0,ix,iy,iz));
      */

      for(ix=bix1;ix<bix2;ix++)
	for(iy=biy1;iy<biy2;iy++)
	  for(iz=0;iz<NZ;iz++)
	    {
	      struct geometry geom;
	      fill_geometry(ix,iy,iz,&geom);
	
	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	
	      //coordinate
	      ldouble dx[3];
	      get_cell_size_arb(ix,iy,iz,dx,OUTCOORDS);

	      //primitives at the cell - either averaged or original, in BL or MYCOORDS
	      for(iv=0;iv<NV;iv++)
		pp[iv]=get_u(p,iv,ix,iy,iz);

	      //to BL, res-files and primitives in avg in MYCOORDS
	      trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);
	      
	      //iy+1 cell, to get gradients:
	      struct geometry geom2;
	      fill_geometry(ix,iy+1,iz,&geom2);
	
	      struct geometry geomBL2;
	      fill_geometry_arb(ix,iy+1,iz,&geomBL2,OUTCOORDS);
	      
	      //primitives at the cell - either averaged or original, in BL or MYCOORDS
	      ldouble pp2[NV];
	      for(iv=0;iv<NV;iv++)
		pp2[iv]=get_u(p,iv,ix,iy+1,iz);

	      //to BL, res-files and primitives in avg in MYCOORDS
	      trans_pall_coco(pp2,pp2,MYCOORDS,OUTCOORDS,geom2.xxvec,&geom2,&geomBL2);

	      //from now on - working in BL coords
        
	      //primitives and derivatives
	      ldouble rho=pp[RHO];
	      ldouble uint=pp[UU];
	      ldouble Ti,Te,trad;
	      ldouble temp=calc_PEQ_Teifrompp(pp,&Te,&Ti,ix,iy,iz);
	      ldouble bsq=0.;
	      ldouble ucon[4],ucov[4],bcon[4],bcov[4],vel[4];
	      ucon[1]=pp[VX];
	      ucon[2]=pp[VY];
	      ucon[3]=pp[VZ];
		  
#ifdef MAGNFIELD
	      ldouble brbphi, bfake[4],bsqint;
	      calc_angle_brbphibsq(ix,iy,iz,&brbphi,&bsqint,bfake,bfake);
	      bangle=-brbphi/bsqint;
	      calc_bcon_prim(pp,bcon,&geomBL);
	      indices_21(bcon,bcov,geomBL.gg); 
	      bsq = dotB(bcon,bcov);
	      banglep=fabs(bcon[1]*sqrt(geomBL.gg[1][1])/(bcon[2]*sqrt(geomBL.gg[2][2])));
	      brsq=bcon[1]*bcov[1];
	      br=bcon[1]*sqrt(geomBL.gg[1][1]);
	      bthsq=bcon[2]*bcov[2];
	      bphsq=bcon[3]*bcov[3];
#endif

	      conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
	      ldouble rhouconr=rho*ucon[1];

	      ldouble Tij[4][4],Tij22[4][4],Rij[4][4];
	      calc_Tij(pp,&geomBL,Tij22);
	      indices_2221(Tij22,Tij,geomBL.gg);
	      ldouble Trt = Tij[1][0],Rrt=0.;
	      ldouble Trtmagn = bsq*ucon[1]*ucov[0] - bcon[1]*bcov[0];
	      ldouble Trtkin =  rho*ucon[1]*ucov[0];
	      ldouble enden = Tij[0][0] + rho*ucon[0];
	      ldouble Ehat=0.;
	      ldouble Ehat2=0.;
	      ldouble Gi[4],Giff[4]={0.,0.,0.,0.};
	      ldouble Gic[4],Gicff[4]={0.,0.,0.,0.};

	      
	      ldouble uconr[4];
#ifdef RADIATION
	      calc_Rij(pp,&geomBL,Rij);
	      indices_2221(Rij,Rij,geomBL.gg);
	      Rrt = Rij[1][0];
	      ldouble kappaes=calc_kappaes(pp,&geomBL);
	      ldouble Rtt;
	      calc_ff_Rtt(pp,&Rtt,uconr,&geomBL);
	      Ehat=-Rtt; 	

	      //iy+1
	      calc_ff_Ehat(pp2,&Ehat2,uconr,&geomBL2);

	      //diffusive flux
	      Fdiffth=1./3./kappaes*(Ehat2-Ehat)/((geomBL2.yy-geomBL.yy)*geomBL.xx);

	      //adv flux
	      //TODO : additional 4/3 factor should be used when converting to advective flux
	      Fadvth = -Ehat*ucon[2]*geomBL.xx;

	      //total flux
	      Fradth = Rij[2][0]*geomBL.xx;

	      enden+=Rij[0][0];

	      //radiative velocity
	      uconr[1]=pp[FX];
	      uconr[2]=pp[FY];
	      uconr[3]=pp[FZ];
	      conv_vels(uconr,uconr,VELPRIMRAD,VEL4,geomBL.gg,geomBL.GG);

	   


	      //four fource
	      calc_Gi(pp,&geomBL,Giff,0.0, 0,0); //!AC urel_old=0
#if defined(COMPTONIZATION) || defined(EVOLVEPHOTONNUMBER) || defined(NCOMPTONIZATION)
	      ldouble uconff[4];
              uconff[1]=uconff[2]=uconff[3]=0.;
              uconff[0]=1.;
	      calc_Compt_Gi(pp,&geomBL,Gicff,Ehat,Te,kappaes,uconff);
#endif 

	   

	      trad=calc_LTE_TfromE(Ehat);
#ifdef EVOLVEPHOTONNUMBER
	      trad=calc_ncompt_Thatrad(pp,&geomBL,Ehat);
#endif
#endif

	      ldouble pregas = GAMMAM1*uint;
	      ldouble premag = bsq/2.;
	      ldouble prerad = 0.;
#ifdef RADIATION
	      prerad = Ehat/3.;
#endif
	      ldouble pretot = pregas + premag + prerad;

	      //alpha 
	      ldouble Tijff[4][4],Tij22ff[4][4];
	      //boost to Keplerian frame! - makes little difference                                                                             
              ldouble ppKepl[NV];
	      for(iv=0;iv<NV;iv++)
                ppKepl[iv]=pp[iv];
              ppKepl[VX]=ppKepl[VY]=0.;
              ldouble Omk=sqrt(1./geomBL.xx/geomBL.xx/geomBL.xx);
              ppKepl[VZ]=Omk;
              boost22_lab2ff(Tij22,Tij22ff,ppKepl,geomBL.gg,geomBL.GG);
              indices_2221(Tij22ff,Tijff,geomBL.gg);
              ldouble Trphi=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij22ff[1][3];
              ldouble dOmdr=-3./2./geomBL.xx*Omk;
              ldouble qplusTrphi=Trphi*geomBL.xx*dOmdr;
	      ldouble qplus=calc_ViscousHeating(ix,iy,iz);
	      ldouble alpha=Trphi/pretot;

	      //angular velocity
	      ldouble Omega=ucon[3]/ucon[0];
	      //MRI resolution parameters
	      ldouble qtheta,qphi;
	      calc_Qthetaphi(ix,iy,iz,&qtheta,&qphi);

	      //PLACE - overwrite with avg quantities if required - currently not complete and not maintained
	      if(doingavg)
		{
		  rho=get_uavg(pavg,RHO,ix,iy,iz);
		  uint=get_uavg(pavg,UU,ix,iy,iz);
		  bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
		  temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
		  //temp=get_uavg(pavg,AVGTGAS,ix,iy,iz);
		  ucon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  ucon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  ucon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  ucon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  ldouble uconpp[4]={ucon[0],ucon[1],ucon[2],ucon[3]};
		  conv_vels(uconpp,uconpp,VEL4,VEL4,geomBL.gg,geomBL.GG);
		  conv_vels(uconpp,uconpp,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
		  pp[VX]=uconpp[1];
		  pp[VY]=uconpp[2];
		  pp[VZ]=uconpp[3];


		  Omega=ucon[3]/ucon[0];
		  
		  for(i=0;i<4;i++)
		    {
		      for(j=0;j<4;j++)
			{
			  Tij[i][j]=get_uavg(pavg,AVGRHOUCONUCOV(i,j),ix,iy,iz)
			    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(i,j),ix,iy,iz)
			    + get_uavg(pavg,AVGBSQUCONUCOV(i,j),ix,iy,iz)
			    - get_uavg(pavg,AVGBCONBCOV(i,j),ix,iy,iz); 
#ifdef RADIATION
			  Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 
#endif
			}
		      Giff[i]=get_uavg(pavg,AVGGHAT(i),ix,iy,iz);
		      Gicff[i]=get_uavg(pavg,AVGGHATCOMPT(i),ix,iy,iz);
		    }
		  pregas = GAMMAM1*uint;
		  premag = bsq/2.;
		  pretot = pregas + premag;
		  prerad = Ehat = 0.;
#ifdef RADIATION
		  Ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  trad=calc_ncompt_Thatrad_fromEN(Ehat,get_uavg(pavg,AVGNFHAT,ix,iy,iz));
		  prerad = Ehat/3.;
		  pretot+=prerad;
		  //radiative velocity not overwritten yet
#endif
		  //alpha not averaged properly - use snapshots rather or not
		  indices_2122(Tij,Tij22,geomBL.GG);
		  boost22_lab2ff(Tij22,Tij22,pp,geomBL.gg,geomBL.GG);
		  alpha=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij22[1][3]/pretot;

		  //bangle
		  brbphi=get_uavg(pavg,AVGBCONBCOV(1,3),ix,iy,iz)*geomBL.GG[3][3]*sqrt(geomBL.gg[1][1]*geomBL.gg[3][3]);
		  bangle=-brbphi/bsq;

		  banglep=fabs(get_uavg(pavg,AVGBCON(1),ix,iy,iz)*sqrt(geomBL.gg[1][1])/(get_uavg(pavg,AVGBCON(2),ix,iy,iz)*sqrt(geomBL.gg[2][2])));

		  br=get_uavg(pavg,AVGBCON(1),ix,iy,iz)*geomBL.gg[1][1];
		  brsq=get_uavg(pavg,AVGBCONBCOV(1,1),ix,iy,iz);
		  bthsq=get_uavg(pavg,AVGBCONBCOV(2,2),ix,iy,iz);
		  bphsq=get_uavg(pavg,AVGBCONBCOV(3,3),ix,iy,iz);
		  //qplus not handled
		  //neither Qtheta - stays the same
		}

	      ldouble dvol=dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      volume+=dvol;

	      //to calculate scaleheight
	      scale=rho*dvol*(M_PI/2. - geomBL.yy)*(M_PI/2. - geomBL.yy);
	      scaleav+=scale;

	      //integrals
	      mass+=rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      pgasint+=pregas*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      pradint+=prerad*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      ptotint+=pretot*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Gtint+=(Giff[0])*dx[0]*dx[1]*dx[2]*geomBL.gdet;

	      Gctint+=(Gicff[0])*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Qplusint+=qplus*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      QplusintTrphi+=qplusTrphi*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Trphiint+=Trphi*dx[0]*dx[1]*dx[2]*geomBL.gdet;


	      brav+=br*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      brsqav+=brsq*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      bthsqav+=bthsq*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      bphsqav+=bphsq*dx[0]*dx[1]*dx[2]*geomBL.gdet;

	      //rho-averages
	      tempav+=temp*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      tradav+=trad*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      qthetaav+=qtheta*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      alphaav+=alpha*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      urav+=ucon[1]*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      uthav+=ucon[2]*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      uphav+=ucon[3]*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      uradrav+=uconr[1]*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      uradthav+=uconr[2]*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Fdiffthav+=Fdiffth*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Fadvthav+=Fadvth*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Fradthav+=Fradth*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      bangleav+=bangle*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      banglepav+=banglep*rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;

	     
	    }

      //now once again over volume to calculate deviations from mean
      //these not MPI-compatible!
      ldouble meanrho=mass/volume;
      ldouble meanvphi=uphav/mass;

      
      for(ix=bix1;ix<bix2;ix++)
	for(iy=biy1;iy<biy2;iy++)
	  for(iz=0;iz<NZ;iz++)
	    {
	      struct geometry geom;
	      fill_geometry(ix,iy,iz,&geom);
	
	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	
	      //coordinate
	      ldouble dx[3];
	      get_cell_size_arb(ix,iy,iz,dx,OUTCOORDS);

	      //primitives at the cell - either averaged or original, in BL or MYCOORDS
	      for(iv=0;iv<NV;iv++)
		pp[iv]=get_u(p,iv,ix,iy,iz);

	      //to BL, res-files and primitives in avg in MYCOORDS
	      trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);

	      ldouble ucon[4],ucov[4];
	      ucon[1]=pp[VX];
	      ucon[2]=pp[VY];
	      ucon[3]=pp[VZ];
	      conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
	      
	      //primitives and derivatives
	      ldouble rho=pp[RHO];

	      //PLACE - overwrite with avg quantities if required - currently not complete and not maintained
	      if(doingavg)
		{
		  rho=get_uavg(pavg,RHO,ix,iy,iz);
		  ucon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		}

	      //integrals
	      rhodev+=(rho-meanrho)*(rho-meanrho)*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      vphidev+=(ucon[3]-meanvphi)*(ucon[3]-meanvphi)*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	    }

 
      //now integrals of fluxes over the walls [left radial, right radial, top+bottom]
     

      //left radial wall
      if(BOXR1>rmin && BOXR1<=rmax) //within this tile
	{
	  ix=bix1;
	  for(iy=biy1;iy<=biy2;iy++)
	    for(iz=0;iz<NZ;iz++)
	      {
		struct geometry geom;
		fill_geometry(ix,iy,iz,&geom);
	
		struct geometry geomBL;
		fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	
		//coordinates
		ldouble dx[3];
		get_cell_size_arb(ix,iy,iz,dx,OUTCOORDS);
	    

		//primitives at the cell - either averaged or original, in BL or MYCOORDS
		for(iv=0;iv<NV;iv++)
		  pp[iv]=get_u(p,iv,ix,iy,iz);

		//to BL, res-files and primitives in avg in MYCOORDS
		trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);

		//from now on - working in BL coords
        
		//primitives and derivatives
		ldouble rho=pp[RHO];
		ldouble uint=pp[UU];
		ldouble temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
		ldouble bsq=0.;
		ldouble ucon[4],ucov[4],bcon[4],bcov[4],vel[4];
		ucon[1]=pp[VX];
		ucon[2]=pp[VY];
		ucon[3]=pp[VZ];
		  
#ifdef MAGNFIELD
		calc_bcon_prim(pp,bcon,&geomBL);
		indices_21(bcon,bcov,geomBL.gg); 
		bsq = dotB(bcon,bcov); 
#endif

		conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
		ldouble rhouconr=rho*ucon[1];


		ldouble Tij[4][4],Tij22[4][4],Rij[4][4];
		calc_Tij(pp,&geomBL,Tij22);
		indices_2221(Tij22,Tij,geomBL.gg);
		ldouble Trt = Tij[1][0],Rrt=0.;
		ldouble Trtmagn = bsq*ucon[1]*ucov[0] - bcon[1]*bcov[0];
		ldouble Trtkin =  rho*ucon[1]*ucov[0];

		ldouble Ehat=0.;
#ifdef RADIATION
		calc_Rij(pp,&geomBL,Rij);
		indices_2221(Rij,Rij,geomBL.gg);
		Rrt = Rij[1][0];
		ldouble Rtt,uconr[4];
		calc_ff_Rtt(pp,&Rtt,uconr,&geomBL);
		Ehat=-Rtt;
#endif
		ldouble Ehatuconr = Ehat*ucon[1];

		//PLACE - overwrite with avg quantities if required
		if(doingavg)
		  {
		    Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
		      + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
		      + get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		      - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 
		    Trtmagn= get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		      - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 
		    Trtkin = get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz);
		    rhouconr = get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
#ifdef RADIATION
		    Rrt=get_uavg(pavg,AVGRIJ(1,0),ix,iy,iz); 
		    Ehatuconr = get_uavg(pavg,AVGEHATUCON(1),ix,iy,iz);
#endif
		  }

		//integrals
		rhouconrtot[0]+=rhouconr*dx[1]*dx[2]*geomBL.gdet;
		Trttot[0]+=Trt*dx[1]*dx[2]*geomBL.gdet;
		Trtkintot[0]+=Trtkin*dx[1]*dx[2]*geomBL.gdet;
		Trtmagntot[0]+=Trtmagn*dx[1]*dx[2]*geomBL.gdet;
		Rrttot[0]+=Rrt*dx[1]*dx[2]*geomBL.gdet;
		Ehatuconrtot[0]+=Ehatuconr*dx[1]*dx[2]*geomBL.gdet;
		areas[0]+=dx[2]*dx[1]*sqrt(geomBL.gg[3][3])*sqrt(geomBL.gg[2][2]);
	      }

	}


      //right radial wall (with minus sign)
      if(BOXR2>rmin && BOXR2<=rmax) //within this tile
	{
	  ix=bix2;
	  for(iy=biy1;iy<=biy2;iy++)
	    for(iz=0;iz<NZ;iz++)
	      {
		struct geometry geom;
		fill_geometry(ix,iy,iz,&geom);
	
		struct geometry geomBL;
		fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	
		//coordinates
		ldouble dx[3];
		get_cell_size_arb(ix,iy,iz,dx,OUTCOORDS);
	    
		//primitives at the cell - either averaged or original, in BL or MYCOORDS
		for(iv=0;iv<NV;iv++)
		  pp[iv]=get_u(p,iv,ix,iy,iz);

		//to BL, res-files and primitives in avg in MYCOORDS
		trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);

		//from now on - working in BL coords
        
		//primitives and derivatives
		ldouble rho=pp[RHO];
		ldouble uint=pp[UU];
		ldouble temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
		ldouble bsq=0.;
		ldouble ucon[4],ucov[4],bcon[4],bcov[4],vel[4];
		ucon[1]=pp[VX];
		ucon[2]=pp[VY];
		ucon[3]=pp[VZ];
		  
#ifdef MAGNFIELD
		calc_bcon_prim(pp,bcon,&geomBL);
		indices_21(bcon,bcov,geomBL.gg); 
		bsq = dotB(bcon,bcov); 
#endif

		conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
		ldouble rhouconr=rho*ucon[1];

		ldouble Tij[4][4],Tij22[4][4],Rij[4][4];
		calc_Tij(pp,&geomBL,Tij22);
		indices_2221(Tij22,Tij,geomBL.gg);
		ldouble Trt = Tij[1][0],Rrt=0.;
		ldouble Trtmagn = bsq*ucon[1]*ucov[0] - bcon[1]*bcov[0];
		ldouble Trtkin =  rho*ucon[1]*ucov[0];

		ldouble Ehat=0.;
#ifdef RADIATION
		calc_Rij(pp,&geomBL,Rij);
		indices_2221(Rij,Rij,geomBL.gg);
		Rrt = Rij[1][0];
		ldouble Rtt,uconr[4];
		calc_ff_Rtt(pp,&Rtt,uconr,&geomBL);
		Ehat=-Rtt;
#endif
		ldouble Ehatuconr = Ehat*ucon[1];

		//PLACE - overwrite with avg quantities if required
		if(doingavg)
		  {
		    Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
		      + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
		      + get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		      - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 
		    Trtmagn= get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
		      - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 
		    Trtkin = get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz);
		    rhouconr = get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
#ifdef RADIATION
		    Rrt=get_uavg(pavg,AVGRIJ(1,0),ix,iy,iz); 
		    Ehatuconr = get_uavg(pavg,AVGEHATUCON(1),ix,iy,iz);
#endif
		  }

		//integrals
		rhouconrtot[1]-=rhouconr*dx[1]*dx[2]*geomBL.gdet;
		Trttot[1]-=Trt*dx[1]*dx[2]*geomBL.gdet;
		Trtkintot[1]-=Trtkin*dx[1]*dx[2]*geomBL.gdet;
		Trtmagntot[1]-=Trtmagn*dx[1]*dx[2]*geomBL.gdet;
		Rrttot[1]-=Rrt*dx[1]*dx[2]*geomBL.gdet;
		Ehatuconrtot[1]-=Ehatuconr*dx[1]*dx[2]*geomBL.gdet;
		areas[1]+=dx[2]*dx[1]*sqrt(geomBL.gg[3][3])*sqrt(geomBL.gg[2][2]);
	      }

	}


      //top face (with minus sign)
      if(ITH1>=giy1 && ITH1<=giy2)
	{
	  iy=biy1;
	  for(ix=bix1;ix<=bix2;ix++)
	    for(iz=0;iz<NZ;iz++)
	      {
		struct geometry geom;
		fill_geometry(ix,iy,iz,&geom);
	
		struct geometry geomBL;
		fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	
		//coordinates
		ldouble dx[3];
		get_cell_size_arb(ix,iy,iz,dx,OUTCOORDS);
	    
		//primitives at the cell - either averaged or original, in BL or MYCOORDS
		for(iv=0;iv<NV;iv++)
		  pp[iv]=get_u(p,iv,ix,iy,iz);

		//to BL, res-files and primitives in avg in MYCOORDS
		trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);

		//from now on - working in BL coords
        
		//primitives and derivatives
		ldouble rho=pp[RHO];
		ldouble uint=pp[UU];
		ldouble temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
		ldouble bsq=0.;
		ldouble ucon[4],ucov[4],bcon[4],bcov[4],vel[4];
		ucon[1]=pp[VX];
		ucon[2]=pp[VY];
		ucon[3]=pp[VZ];
		  
#ifdef MAGNFIELD
		calc_bcon_prim(pp,bcon,&geomBL);
		indices_21(bcon,bcov,geomBL.gg); 
		bsq = dotB(bcon,bcov); 
#endif

		conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
		ldouble rhouconth=rho*ucon[2];

		ldouble Tij[4][4],Tij22[4][4],Rij[4][4];
		calc_Tij(pp,&geomBL,Tij22);
		indices_2221(Tij22,Tij,geomBL.gg);
		ldouble Ttht = Tij[2][0],Rtht=0.;
		ldouble Tthtmagn = bsq*ucon[2]*ucov[0] - bcon[2]*bcov[0];
		ldouble Tthtkin =  rho*ucon[2]*ucov[0];

		ldouble Ehat=0.;
#ifdef RADIATION
		calc_Rij(pp,&geomBL,Rij);
		indices_2221(Rij,Rij,geomBL.gg);
		Rtht = Rij[2][0];
		ldouble Rtt,uconr[4];
		calc_ff_Rtt(pp,&Rtt,uconr,&geomBL);
		Ehat=-Rtt;
#endif	
		ldouble Ehatuconth = Ehat*ucon[2];

		//PLACE - overwrite with avg quantities if required
		if(doingavg)
		  {
		    Ttht=get_uavg(pavg,AVGRHOUCONUCOV(2,0),ix,iy,iz)
		      + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(2,0),ix,iy,iz)
		      + get_uavg(pavg,AVGBSQUCONUCOV(2,0),ix,iy,iz)
		      - get_uavg(pavg,AVGBCONBCOV(2,0),ix,iy,iz); 
		    Tthtmagn= get_uavg(pavg,AVGBSQUCONUCOV(2,0),ix,iy,iz)
		      - get_uavg(pavg,AVGBCONBCOV(2,0),ix,iy,iz); 
		    Tthtkin = get_uavg(pavg,AVGRHOUCONUCOV(2,0),ix,iy,iz);
		    rhouconth = get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz);
#ifdef RADIATION
		    Rtht=get_uavg(pavg,AVGRIJ(2,0),ix,iy,iz); 
		    Ehatuconth = get_uavg(pavg,AVGEHATUCON(2),ix,iy,iz);
#endif
		  }

		//integrals (watch the sign!)
		rhouconrtot[2]+=rhouconth*dx[2]*dx[0]*geomBL.gdet;
		Trttot[2]+=Ttht*dx[2]*dx[0]*geomBL.gdet;
		Trtkintot[2]+=Tthtkin*dx[2]*dx[0]*geomBL.gdet;
		Trtmagntot[2]+=Tthtmagn*dx[2]*dx[0]*geomBL.gdet;
		Rrttot[2]+=Rtht*dx[2]*dx[0]*geomBL.gdet;
		Ehatuconrtot[2]+=Ehatuconth*dx[2]*dx[0]*geomBL.gdet;

		areas[2]+=dx[2]*dx[0]*sqrt(geomBL.gg[1][1])*sqrt(geomBL.gg[3][3]);
	      }

	}

      //bottom face
      if((ITH2+1)>giy1 && (ITH2+1)<=giy2)
	{
	  iy=biy2;
	  for(ix=bix1;ix<=bix2;ix++)
	    for(iz=0;iz<NZ;iz++)
	      {
		struct geometry geom;
		fill_geometry(ix,iy,iz,&geom);
	
		struct geometry geomBL;
		fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	
		//coordinates
		ldouble dx[3];
		get_cell_size_arb(ix,iy,iz,dx,OUTCOORDS);
	    
		//primitives at the cell - either averaged or original, in BL or MYCOORDS
		for(iv=0;iv<NV;iv++)
		  pp[iv]=get_u(p,iv,ix,iy,iz);

		//to BL, res-files and primitives in avg in MYCOORDS
		trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);

		//from now on - working in BL coords
        
		//primitives and derivatives
		ldouble rho=pp[RHO];
		ldouble uint=pp[UU];
		ldouble temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
		ldouble bsq=0.;
		ldouble ucon[4],ucov[4],bcon[4],bcov[4],vel[4];
		ucon[1]=pp[VX];
		ucon[2]=pp[VY];
		ucon[3]=pp[VZ];
		  
#ifdef MAGNFIELD
		calc_bcon_prim(pp,bcon,&geomBL);
		indices_21(bcon,bcov,geomBL.gg); 
		bsq = dotB(bcon,bcov); 
#endif

		conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
		ldouble rhouconth=rho*ucon[2];

		ldouble Tij[4][4],Tij22[4][4],Rij[4][4];
		calc_Tij(pp,&geomBL,Tij22);
		indices_2221(Tij22,Tij,geomBL.gg);
		ldouble Ttht = Tij[2][0],Rtht=0.;
		ldouble Tthtmagn = bsq*ucon[2]*ucov[0] - bcon[2]*bcov[0];
		ldouble Tthtkin =  rho*ucon[2]*ucov[0];

		ldouble Ehat=0.;
#ifdef RADIATION
		calc_Rij(pp,&geomBL,Rij);
		indices_2221(Rij,Rij,geomBL.gg);
		Rtht = Rij[2][0];
		ldouble Rtt,uconr[4];
		calc_ff_Rtt(pp,&Rtt,uconr,&geomBL);
		Ehat=-Rtt;
#endif	
		ldouble Ehatuconth = Ehat*ucon[2];

		//PLACE - overwrite with avg quantities if required
		if(doingavg)
		  {
		    Ttht=get_uavg(pavg,AVGRHOUCONUCOV(2,0),ix,iy,iz)
		      + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(2,0),ix,iy,iz)
		      + get_uavg(pavg,AVGBSQUCONUCOV(2,0),ix,iy,iz)
		      - get_uavg(pavg,AVGBCONBCOV(2,0),ix,iy,iz); 
		    Tthtmagn= get_uavg(pavg,AVGBSQUCONUCOV(2,0),ix,iy,iz)
		      - get_uavg(pavg,AVGBCONBCOV(2,0),ix,iy,iz); 
		    Tthtkin = get_uavg(pavg,AVGRHOUCONUCOV(2,0),ix,iy,iz);
		    rhouconth = get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz);
#ifdef RADIATION
		    Rtht=get_uavg(pavg,AVGRIJ(2,0),ix,iy,iz); 
		    Ehatuconth = get_uavg(pavg,AVGEHATUCON(2),ix,iy,iz);
#endif
		    

		  }

		//integrals (watch the sign!)
		rhouconrtot[2]+=-rhouconth*dx[2]*dx[0]*geomBL.gdet;
		Trttot[2]+=-Ttht*dx[2]*dx[0]*geomBL.gdet;
		Trtkintot[2]+=-Tthtkin*dx[2]*dx[0]*geomBL.gdet;
		Trtmagntot[2]+=-Tthtmagn*dx[2]*dx[0]*geomBL.gdet;
		Rrttot[2]+=-Rtht*dx[2]*dx[0]*geomBL.gdet;
		Ehatuconrtot[2]+=-Ehatuconth*dx[2]*dx[0]*geomBL.gdet;

		areas[2]+=dx[2]*dx[0]*sqrt(geomBL.gg[1][1])*sqrt(geomBL.gg[3][3]);
	      }
	}

      
      
      //saving local integrals 
      boxscalars[0]=mass; // (2nd column) - total mass inside the box
      boxscalars[1]=ptotint; // (3) - total integrated pressure | beta = pmag/ptot = (($3-$4-$5)/$3)
      boxscalars[2]=pgasint; // (4) - gas integrated pressure
      boxscalars[3]=pradint; // (5) - rad integrated pressure 
      boxscalars[4]=Gtint; // (6) - integrated total heating (\hat G_abs^t + \hat G_compt^t)
      boxscalars[5]=Gctint; // (7) - integrated Compton heating (\hat G_compt^t)

      boxscalars[6]=tempav; // (8) - averaged temperature
      boxscalars[7]=tradav; // (9) - averaged temperature of radiation
      boxscalars[8]=qthetaav; // (10) - averaged Qtheta
      boxscalars[9]=alphaav; // (11) - averaged alpha
      boxscalars[10]=urav; // (12) - averaged radial vel
      boxscalars[11]=uthav; // (13) - averaged polar vel
      boxscalars[12]=uradrav; // (14) - averaged radial vel
      boxscalars[13]=uradthav; // (15) - averaged polar vel
      boxscalars[14]=Fdiffthav; // (16) - diffusive flux in theta
      boxscalars[15]=Fadvthav; // (17) - adv flux in theta
      boxscalars[16]=Fradthav; // (18) - total flux in theta
      boxscalars[17]=bangleav; // (19) - averaged B field angle brbphi/bsq (?)
      boxscalars[18]=scaleav; // (20) - scale height (in squared radians from eq.plane) -

      boxscalars[19]=volume; // (21) - volume

      int nia=19; //number of boxscalars slots filled so far

      boxscalars[nia+1]=areas[0]; // (22)
      boxscalars[nia+2]=areas[1];
      boxscalars[nia+3]=areas[2];

      int nia2=nia+3;

      boxscalars[nia2+1]=rhouconrtot[0]; // (25) - mass flux through inner face
      boxscalars[nia2+2]=rhouconrtot[1]; // (+1) - (-)mass flux through right face
      boxscalars[nia2+3]=rhouconrtot[2]; // (+2) - mass flux outflowing through top and bottom together
      boxscalars[nia2+4]=Trttot[0]; // (+3) - Trt flux 
      boxscalars[nia2+5]=Trttot[1]; // (+4) - -Trt flux
      boxscalars[nia2+6]=Trttot[2]; // (+5) - Ttht flux
      boxscalars[nia2+7]=Rrttot[0]; // (+6) - Rrt flux
      boxscalars[nia2+8]=Rrttot[1]; // (+7) - -Rrt flux
      boxscalars[nia2+9]=Rrttot[2]; // (+8) - Rtht flux
      boxscalars[nia2+10]=Ehatuconrtot[0]; // (+9) - Ehatuconr flux
      boxscalars[nia2+11]=Ehatuconrtot[1]; // (+10) - -Ehatuconr flux
      boxscalars[nia2+12]=Ehatuconrtot[2]; // (+11) - Ehatuconth flux

      int nia3=nia2+12;
      
     boxscalars[nia3+1]=QplusintTrphi; // (37) - integrated heating (Trphi*dOm/dr)
     boxscalars[nia3+2]=Trphiint; // (38) - integrated Trphi in the fluid frame
     boxscalars[nia3+3]=Qplusint; // (39) - integrated heating (from entropy change)

     boxscalars[nia3+4]=rhodev; // (40) - deviation of rho from mean rho
     boxscalars[nia3+5]=vphidev; // (41) - deviation of vphi from mean rho
     boxscalars[40]=brsqav; // (42) - brsq
     boxscalars[41]=bthsqav; // (43) - bthsq
     boxscalars[42]=bphsqav; // (44) - bphsq
     boxscalars[43]=brav; // (45) - br
   
    } //if(!ifoutsidebox)

  //TEST
  //printf("pre %d > %e\n",PROCID,boxscalars[4]);

  
  //aggregating over all tiles to the master who will then print out
#ifdef MPI
  ldouble bscsum[NBOXSCALARS];
  MPI_Reduce(boxscalars, bscsum, NBOXSCALARS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  for(iv=0;iv<NBOXSCALARS;iv++)
    boxscalars[iv]=bscsum[iv];
#endif

  // if(PROCID==0) printf("sum %d > %e\n",PROCID,boxscalars[4]);

  if(PROCID==0)
    { 
      //normalizing the rho-averages by the total mass
      ldouble mass=boxscalars[0];
      for(iv=6;iv<=18;iv++)
	boxscalars[iv]/=mass;
      
      //taking the root of scaleheight
      boxscalars[18]=sqrt(boxscalars[18]);

      //normalizing by volume to obtain standard deviation
      boxscalars[38]=sqrt(boxscalars[38]/boxscalars[19]);
      boxscalars[39]=sqrt(boxscalars[39]/boxscalars[19]);
    }

  


 
#endif //BOXOUTPUT==1
  return 0;
}


/*******************************************************/
/* calculates vertical profiles of quantities within a box  */
/*******************************************************/

int calc_boxvertscalars(ldouble boxvertscalars[TNY][NBOXVERTSCALARS],ldouble t)
{
#ifdef CONSISTENTGAMMA
  if(PROCID==0) printf("gammas not consistent in boxvertscalars postprocessing when CONSISTENTGAMMA\n");
#endif
#ifdef MPI
  if(PROCID==0) printf("boxvertscalars rather does not work with MPI\n");
#endif
  
  int sizey;

#if(BOXVERTOUTPUT==1)
  //adjust NBOXVERTSCALARS in problem.h
  int ix,iy,iz,ii,iv,i,j;
  int bix1,bix2,biy1,biy2,giy;
  ldouble xx[4],xxBL[4];

  

  //search for appropriate indices
  //limits of this tile
  int rmin,rmax;
  get_xx(0,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
  rmin=xxBL[1];
  get_xx(NX-1,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
  rmax=xxBL[1];

  //global index of the theta limits
  int giy1,giy2;
  mpi_local2globalidx(0,0,0,&ix,&giy1,&iz);
  mpi_local2globalidx(0,NY-1,0,&ix,&giy2,&iz);

  //if tile completely outside the box
  int ifoutsidebox=0;

  //theta indices of the box
  int ITH1, ITH2;
#ifdef BOXITH //symmetrical wrt eq.plane
  ITH1=TNY/2-BOXITH;
  ITH2=TNY/2+BOXITH-1;
#endif
#ifdef BOXITH1
  ITH1=BOXITH1;
#endif
#ifdef BOXITH2
  ITH2=BOXITH2;
#endif


  

  //radial limits first 
  bix1=bix2=-1;
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
      if(xxBL[1]>BOXR1 & bix1<0) bix1=ix;
      if(xxBL[1]>BOXR2 & bix2<0) bix2=ix;
    }
  if(bix1<0) bix1=NX;
  if(bix2<0) bix2=NX;

  //then polar angle

  biy1=biy2=-1;
  for(iy=0;iy<NY;iy++)
    {
      mpi_local2globalidx(0,iy,0,&ix,&giy,&iz);
      
      if(giy>ITH1 && biy1<0) biy1=iy;
      if(giy>ITH2 && biy2<0) biy2=iy;
    }
  if(biy1<0) biy1=NY;
  if(biy2<0) biy2=NY;

  //vertical number of cells - to be returned at the end
  sizey=biy2-biy1;

  
 

  //printf("PROCID: %d > bix: %d - %d > biy: %d - %d -> giy: %d - %d\n",PROCID,bix1,bix2,biy1,biy2,biy1+TOJ,biy2+TOJ);

  //first integrals / averages within the box at given theta slice
  ldouble mass,pgasint,pradint,ptotint,Gtint,Gctint,volume,Fdiffth,Fadvth,Fadvr,dFadvdr,Fadv2th,Fradth,bangle,banglep,scale,Trphiint,QplusintTrphi,Qplusint;
  ldouble brsq,bthsq,bphsq,brsqav,bthsqav,bphsqav,br,brav;
  ldouble tempav,qthetaav,alphaav,uradrav,uradthav,tradav,urav,uthav,uphav,Fdiffthav,Fadvthav,Fadvrav,dFadvrdrav,dFtotrdrav,Fradthav,Fadv2thav,bangleav,banglepav,scaleav;
  ldouble dpgasdz,dpmagdz,dpraddz;
  ldouble rhodev=0.,vphidev=0.;
  ldouble pp[NV];
  ldouble Rrttot[3]={0.,0.,0.};
  ldouble Trttot[3]={0.,0.,0.};
  ldouble Trtkintot[3]={0.,0.,0.};
  ldouble Trtmagntot[3]={0.,0.,0.};
  ldouble rhouconrtot[3]={0.,0.,0.};
  ldouble Ehatuconrtot[3]={0.,0.,0.};
  ldouble areas[3]={0.,0.,0.};
  ldouble theta;

  for(iy=biy1;iy<biy2;iy++)
    {
      //zero scalars to start with
      for(j=0;j<NBOXVERTSCALARS;j++)
	boxvertscalars[iy-biy1][j]=0.;

      dpgasdz=dpraddz=dpmagdz=0.;
      mass=pgasint=pradint=ptotint=Gtint=Gctint=Qplusint=QplusintTrphi=Trphiint=0.;
      bangle=banglep=0.;
      brsq=bthsq=bphsq=br=0.;
      brsqav=bthsqav=bphsqav=brav=0.;
      tempav=tradav=qthetaav=alphaav=urav=uthav=uphav=uradthav=uradrav=volume=Fdiffthav=Fadvrav=dFadvrdrav=dFtotrdrav=Fadvthav=Fadv2thav=Fradthav=bangleav=banglepav=scaleav=0.;
 
      for(ix=bix1;ix<bix2;ix++)
	{
	  for(iz=0;iz<NZ;iz++)
	    {
	      struct geometry geom;
	      fill_geometry(ix,iy,iz,&geom);
	
	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);

	      theta=geomBL.yy;
	      //printf("%d %d %d > %f\n",ix,iy,iz,theta);
	
	      //coordinate
	      ldouble dx[3];
	      get_cell_size_arb(ix,iy,iz,dx,OUTCOORDS);

	      //primitives at the cell - either averaged or original, in BL or MYCOORDS
	      for(iv=0;iv<NV;iv++)
		pp[iv]=get_u(p,iv,ix,iy,iz);

	      //to BL, res-files and primitives in avg in MYCOORDS
	      trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);
	      
	      //iy+1 cell, to get gradients:
	      struct geometry geomp1;
	      fill_geometry(ix,iy+1,iz,&geomp1);
	
	      struct geometry geomm1;
	      fill_geometry(ix,iy-1,iz,&geomm1);
	
	      struct geometry geomBLp1;
	      fill_geometry_arb(ix,iy+1,iz,&geomBLp1,OUTCOORDS);
	      
	      struct geometry geomBLm1;
	      fill_geometry_arb(ix,iy-1,iz,&geomBLm1,OUTCOORDS);
	      
	      //primitives at the cell - either averaged or original, in BL or MYCOORDS
	      ldouble ppp1[NV];
	      for(iv=0;iv<NV;iv++)
		ppp1[iv]=get_u(p,iv,ix,iy+1,iz);
	      ldouble ppm1[NV];
	      for(iv=0;iv<NV;iv++)
		ppm1[iv]=get_u(p,iv,ix,iy-1,iz);

	      //to BL, res-files and primitives in avg in MYCOORDS
	      trans_pall_coco(ppp1,ppp1,MYCOORDS,OUTCOORDS,geomp1.xxvec,&geomp1,&geomBLp1);
	      trans_pall_coco(ppm1,ppm1,MYCOORDS,OUTCOORDS,geomm1.xxvec,&geomm1,&geomBLm1);
 
	      //from now on - working in BL coords
        
	      //primitives and derivatives
	      ldouble rho=pp[RHO];
	      ldouble uint=pp[UU];
	      ldouble Ti,Te,trad;
	      ldouble temp=calc_PEQ_Teifrompp(pp,&Te,&Ti,ix,iy,iz);
	      ldouble bsq,bsqp1,bsqm1;
	      bsq=bsqp1=bsqm1=0.;
	      ldouble ucon[4],ucov[4],bcon[4],bcov[4],vel[4];
	      ucon[1]=pp[VX];
	      ucon[2]=pp[VY];
	      ucon[3]=pp[VZ];
		  
#ifdef MAGNFIELD
	      ldouble brbphi, bfake[4],bsqint;

	      calc_bcon_prim(ppp1,bcon,&geomBL);
	      indices_21(bcon,bcov,geomBL.gg); 
	      bsqp1 = dotB(bcon,bcov);

	      calc_bcon_prim(ppm1,bcon,&geomBL);
	      indices_21(bcon,bcov,geomBL.gg); 
	      bsqm1 = dotB(bcon,bcov);

	      calc_angle_brbphibsq(ix,iy,iz,&brbphi,&bsqint,bfake,bfake);
	      bangle=-brbphi/bsqint;
	      calc_bcon_prim(pp,bcon,&geomBL);
	      indices_21(bcon,bcov,geomBL.gg); 
	      bsq = dotB(bcon,bcov);
	      banglep=fabs(bcon[1]*sqrt(geomBL.gg[1][1])/(bcon[2]*sqrt(geomBL.gg[2][2])));
	      brsq=bcon[1]*bcov[1];
	      br=bcon[1]*sqrt(geomBL.gg[1][1]);
	      bthsq=bcon[2]*bcov[2];
	      bphsq=bcon[3]*bcov[3];

	      
#endif

	      conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
	      ldouble rhouconr=rho*ucon[1];

	      ldouble Tij[4][4],Tij22[4][4],Rij[4][4];
	      calc_Tij(pp,&geomBL,Tij22);
	      indices_2221(Tij22,Tij,geomBL.gg);
	      ldouble Trt = Tij[1][0],Rrt=0.;
	      ldouble Trtmagn = bsq*ucon[1]*ucov[0] - bcon[1]*bcov[0];
	      ldouble Trtkin =  rho*ucon[1]*ucov[0];
	      ldouble enden = Tij[0][0] + rho*ucon[0];
	      ldouble Ehat=0.;
	      ldouble Ehatp1=0.;
	      ldouble Ehatm1=0.;
	      ldouble Gi[4],Giff[4]={0.,0.,0.,0.};
	      ldouble Gic[4],Gicff[4]={0.,0.,0.,0.};

	      
	      ldouble uconr[4];
#ifdef RADIATION
	      calc_Rij(pp,&geomBL,Rij);
	      indices_2221(Rij,Rij,geomBL.gg);
	      Rrt = Rij[1][0];
	      ldouble kappaes=calc_kappaes(pp,&geomBL);
	      ldouble Rtt;
	      calc_ff_Rtt(pp,&Rtt,uconr,&geomBL);
	      Ehat=-Rtt; 	

	      //iy+1
	      calc_ff_Ehat(ppp1,&Ehatp1,uconr,&geomBLp1);

	      //iy-1
	      calc_ff_Ehat(ppm1,&Ehatm1,uconr,&geomBLm1);

	      //diffusive flux
	      Fdiffth=-1./3./kappaes*(Ehatp1-Ehatm1)/((geomBLp1.yy-geomBLm1.yy)*geomBL.xx);

	      //adv flux (includes buoyant)
	      Fadvth = 4./3.*Ehat*ucon[2]*geomBL.xx;

	      //radial adv flux (includes buoyant)
	      Fadvr = 4./3.*Ehat*ucon[1];

	      //adv flux (no buoyant)
	      Fadv2th = 4./3.*Ehat*ucon[2]*geomBL.xx;

	      //total flux
	      Fradth = -Rij[2][0]*geomBL.xx;

	      enden+=Rij[0][0];

	      //radiative velocity
	      uconr[1]=pp[FX];
	      uconr[2]=pp[FY];
	      uconr[3]=pp[FZ];
	      conv_vels(uconr,uconr,VELPRIMRAD,VEL4,geomBL.gg,geomBL.GG);

	      //four fource
	      calc_Gi(pp,&geomBL,Giff,0.0, 0, 0); //ANDREW fluid frame=0 
#if defined(COMPTONIZATION) || defined(EVOLVEPHOTONNUMBER) || defined(NCOMPTONIZATION)
	      ldouble uconff[4];
	      uconff[1]=uconff[2]=uconff[3]=0.;
	      uconff[0]=1.;
	      calc_Compt_Gi(pp,&geomBL,Gicff,Ehat,Te,kappaes,uconff);
#endif 

	      trad=calc_LTE_TfromE(Ehat);
#ifdef EVOLVEPHOTONNUMBER
	      trad=calc_ncompt_Thatrad(pp,&geomBL,Ehat);
#endif
#endif

	      ldouble pregas = GAMMAM1*uint;
	      ldouble pregasp1 = GAMMAM1*ppp1[UU];
	      ldouble pregasm1 = GAMMAM1*ppm1[UU];

	      ldouble premag = bsq/2.;
	      ldouble premagp1 = bsqp1/2.;
	      ldouble premagm1 = bsqm1/2.;
	      ldouble prerad,preradp1,preradm1;
	      prerad=preradp1=preradm1= 0.;
#ifdef RADIATION
	      prerad = Ehat/3.;
	      preradp1 = Ehatp1/3.;
	      preradm1 = Ehatm1/3.;
#endif
	      ldouble pretot = pregas + premag + prerad;

	      //alpha 
	      ldouble Tijff[4][4],Tij22ff[4][4];
	      //boost to Keplerian frame! - makes little difference                                                                             
	      ldouble ppKepl[NV];
	      for(iv=0;iv<NV;iv++)
		ppKepl[iv]=pp[iv];
	      ppKepl[VX]=ppKepl[VY]=0.;
	      ldouble Omk=sqrt(1./geomBL.xx/geomBL.xx/geomBL.xx);
	      ppKepl[VZ]=Omk;
	      boost22_lab2ff(Tij22,Tij22ff,ppKepl,geomBL.gg,geomBL.GG);
	      indices_2221(Tij22ff,Tijff,geomBL.gg);
	      ldouble Trphi=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij22ff[1][3];
	      ldouble dOmdr=-3./2./geomBL.xx*Omk;
	      ldouble qplusTrphi=Trphi*geomBL.xx*dOmdr;
	      ldouble qplus=calc_ViscousHeating(ix,iy,iz);
	      ldouble alpha=Trphi/pretot;

	      //angular velocity
	      ldouble Omega=ucon[3]/ucon[0];
	      //MRI resolution parameters
	      ldouble qtheta,qphi;
	      calc_Qthetaphi(ix,iy,iz,&qtheta,&qphi);

	      //PLACE - overwrite with avg quantities if required - currently not complete 
	      if(doingavg)
		{
		  rho=get_uavg(pavg,RHO,ix,iy,iz);
		  uint=get_uavg(pavg,UU,ix,iy,iz);
		  bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
		  //temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
		  temp=get_uavg(pavg,AVGTGAS,ix,iy,iz);
		  ucon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  ucon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  ucon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  ucon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  ldouble uconpp[4]={ucon[0],ucon[1],ucon[2],ucon[3]};
		  conv_vels(uconpp,uconpp,VEL4,VEL4,geomBL.gg,geomBL.GG);
		  conv_vels(uconpp,uconpp,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
		  pp[VX]=uconpp[1];
		  pp[VY]=uconpp[2];
		  pp[VZ]=uconpp[3];
		  Omega=ucon[3]/ucon[0];
		  
		  for(i=0;i<4;i++)
		    {
		      for(j=0;j<4;j++)
			{
			  Tij[i][j]=get_uavg(pavg,AVGRHOUCONUCOV(i,j),ix,iy,iz)
			    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(i,j),ix,iy,iz)
			    + get_uavg(pavg,AVGBSQUCONUCOV(i,j),ix,iy,iz)
			    - get_uavg(pavg,AVGBCONBCOV(i,j),ix,iy,iz); 
#ifdef RADIATION
			  Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 
#endif
			}
		      Giff[i]=get_uavg(pavg,AVGGHAT(i),ix,iy,iz);
		      Gicff[i]=get_uavg(pavg,AVGGHATCOMPT(i),ix,iy,iz);
		    }
		  pregas = GAMMAM1*uint;
		  premag = bsq/2.;

		  pregasp1=GAMMAM1*get_uavg(pavg,UU,ix,iy+1,iz);
		  pregasm1=GAMMAM1*get_uavg(pavg,UU,ix,iy-1,iz);
		  premagp1=get_uavg(pavg,AVGBSQ,ix,iy+1,iz)/2.;
		  premagm1=get_uavg(pavg,AVGBSQ,ix,iy-1,iz)/2.;
		  

		  pretot = pregas + premag;
		  prerad = preradp1 = preradm1 = Ehat = 0.;
#ifdef RADIATION
		  Ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  Ehatp1=get_uavg(pavg,AVGEHAT,ix,iy+1,iz);
		  Ehatm1=get_uavg(pavg,AVGEHAT,ix,iy-1,iz);
		  
		  trad=0.;//get_uavg(pavg,AVGTRAD,ix,iy,iz);
		  prerad = Ehat/3.;

		  preradp1 = Ehatp1/3.;
		  preradm1 = Ehatm1/3.;
	      
		  pretot+=prerad;
		  //radiative velocity not overwritten yet

		  //diffusive flux
		  Fdiffth=-1./3./kappaes*(get_uavg(pavg,AVGEHAT,ix,iy+1,iz)-get_uavg(pavg,AVGEHAT,ix,iy-1,iz))/((geomBLp1.yy-geomBLm1.yy)*geomBL.xx);

		  //adv flux (with buoyant)
		  Fadvth = 4./3.*get_uavg(pavg,AVGEHATUCON(2),ix,iy,iz)*geomBL.xx;
  
		  //adv flux (with buoyant)
		  Fadvr = 4./3.*get_uavg(pavg,AVGEHATUCON(1),ix,iy,iz);
  
		  //adv flux (no buoyant)
		  //Fadv2th = 4./3.*get_uavg(pavg,AVGEHAT,ix,iy,iz)*get_uavg(pavg,AVGUCON(2),ix,iy,iz)*geomBL.xx;
		  Fadv2th = 4./3.*get_uavg(pavg,AVGEHAT,ix,iy,iz)*ucon[2]*geomBL.xx;
  
		  //total flux
		  Fradth =  -get_uavg(pavg,AVGRIJ(2,0),ix,iy,iz)*geomBL.xx;

		  //rad. energy weighted radiative velocity
		  uconr[1]=get_uavg(pavg,AVGEHATUCON(1),ix,iy,iz)/get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  uconr[2]=get_uavg(pavg,AVGEHATUCON(2),ix,iy,iz)/get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  uconr[3]=get_uavg(pavg,AVGEHATUCON(3),ix,iy,iz)/get_uavg(pavg,AVGEHAT,ix,iy,iz);
#endif
		  //alpha not averaged properly - use snapshots rather or not
		  indices_2122(Tij,Tij22,geomBL.GG);
		  boost22_lab2ff(Tij22,Tij22,pp,geomBL.gg,geomBL.GG);
		  alpha=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij22[1][3]/pretot;

		  //bangle
		  brbphi=get_uavg(pavg,AVGBCONBCOV(1,3),ix,iy,iz)*geomBL.GG[3][3]*sqrt(geomBL.gg[1][1]*geomBL.gg[3][3]);
		  bangle=-brbphi/bsq;

		  banglep=fabs(get_uavg(pavg,AVGBCON(1),ix,iy,iz)*sqrt(geomBL.gg[1][1])/(get_uavg(pavg,AVGBCON(2),ix,iy,iz)*sqrt(geomBL.gg[2][2])));

		  br=get_uavg(pavg,AVGBCON(1),ix,iy,iz)*geomBL.gg[1][1];
		  brsq=get_uavg(pavg,AVGBCONBCOV(1,1),ix,iy,iz);
		  bthsq=get_uavg(pavg,AVGBCONBCOV(2,2),ix,iy,iz);
		  bphsq=get_uavg(pavg,AVGBCONBCOV(3,3),ix,iy,iz);
		  //qplus not handled
		  //neither Qtheta - stays the same

		}

	      ldouble dvol=dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      volume+=dvol;

	      //vertical support
	      dpgasdz+=(pregasp1-pregasm1)/((geomBLp1.yy-geomBLm1.yy)*geomBL.xx)/rho*dvol;
	      dpraddz+=(preradp1-preradm1)/((geomBLp1.yy-geomBLm1.yy)*geomBL.xx)/rho*dvol;
	      dpmagdz+=(premagp1-premagm1)/((geomBLp1.yy-geomBLm1.yy)*geomBL.xx)/rho*dvol;

	      //integrals
	      mass+=rho*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      pgasint+=pregas*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      pradint+=prerad*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      ptotint+=pretot*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Gtint+=(Giff[0])*dx[0]*dx[1]*dx[2]*geomBL.gdet;

	      Gctint+=(Gicff[0])*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Qplusint+=qplus*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      QplusintTrphi+=qplusTrphi*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Trphiint+=Trphi*dx[0]*dx[1]*dx[2]*geomBL.gdet;

	      //averages
	      brav+=br*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      brsqav+=brsq*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      bthsqav+=bthsq*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      bphsqav+=bphsq*dx[0]*dx[1]*dx[2]*geomBL.gdet;

	      tempav+=temp*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      tradav+=trad*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      qthetaav+=qtheta*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      alphaav+=alpha*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      urav+=ucon[1]*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      uthav+=ucon[2]*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      uphav+=ucon[3]*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      uradrav+=uconr[1]*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      uradthav+=uconr[2]*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Fdiffthav+=Fdiffth*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Fadvrav+=Fadvr*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Fadvthav+=Fadvth*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Fadv2thav+=Fadv2th*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      Fradthav+=Fradth*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      bangleav+=bangle*dx[0]*dx[1]*dx[2]*geomBL.gdet;
	      banglepav+=banglep*dx[0]*dx[1]*dx[2]*geomBL.gdet;


	     
	    } //iz-loop
	} //ix-loop


      //radial divergences
      ix=bix1;
      ldouble Fadvrix1=0.;
      ldouble Ftotrix1=0.;
      
      ldouble r1;
     
      struct geometry geomBL;
      fill_geometry_arb(ix,iy,0,&geomBL,OUTCOORDS);
      r1=geomBL.xx;

      Fadvrix1/=NZ;
      for(iz=0;iz<NZ;iz++)
	{
	  Fadvrix1+=4./3.*get_uavg(pavg,AVGEHATUCON(1),ix,iy,iz);
	  Ftotrix1+= -get_uavg(pavg,AVGRIJ(1,0),ix,iy,iz);
	}
      Fadvrix1/=NZ;
      Ftotrix1/=NZ;

      ix=bix2;
      ldouble Fadvrix2=0.;
      ldouble Ftotrix2=0.;
      ldouble r2;
      fill_geometry_arb(ix,iy,0,&geomBL,OUTCOORDS);
      r2=geomBL.xx;
      for(iz=0;iz<NZ;iz++)
	{
	  Fadvrix2+=4./3.*get_uavg(pavg,AVGEHATUCON(1),ix,iy,iz);
	  Ftotrix2+= -get_uavg(pavg,AVGRIJ(1,0),ix,iy,iz);
	}
      Fadvrix2/=NZ;
      Ftotrix2/=NZ;

      dFadvrdrav=(Fadvrix2*(r2/r1)*(r2/r1)-Fadvrix1); //radial divergence of the radial advective flux at given theta
      dFtotrdrav=(Ftotrix2*(r2/r1)*(r2/r1)-Ftotrix1); //radial divergence of the radial advective flux at given theta
      

      //saving scalars
      boxvertscalars[iy-biy1][0]=theta; // (1st column) - theta coordinate
      boxvertscalars[iy-biy1][1]=mass/volume; // (2nd column) - density
      boxvertscalars[iy-biy1][2]=ptotint/volume;; // (3) - total integrated pressure | beta = pmag/ptot = (($3-$4-$5)/$3)
      boxvertscalars[iy-biy1][3]=pgasint/volume;; // (4) - gas integrated pressure
      boxvertscalars[iy-biy1][4]=pradint/volume;; // (5) - rad integrated pressure 
      boxvertscalars[iy-biy1][5]=Gtint/volume;; // (6) - integrated total heating (\hat G_abs^t + \hat G_compt^t)
      boxvertscalars[iy-biy1][6]=Gctint/volume;; // (7) - integrated Compton heating (\hat G_compt^t)

      boxvertscalars[iy-biy1][7]=tempav/volume; // (8) - averaged temperature
      boxvertscalars[iy-biy1][8]=tradav/volume; // (9) - averaged temperature of radiation
      boxvertscalars[iy-biy1][9]=qthetaav/volume; // (10) - averaged Qtheta
      boxvertscalars[iy-biy1][10]=alphaav/volume; // (11) - averaged alpha
      boxvertscalars[iy-biy1][11]=urav/volume; // (12) - averaged radial vel
      boxvertscalars[iy-biy1][12]=uthav/volume; // (13) - averaged polar vel
      boxvertscalars[iy-biy1][13]=uradrav/volume; // (14) - averaged radial vel
      boxvertscalars[iy-biy1][14]=uradthav/volume; // (15) - averaged polar vel
      boxvertscalars[iy-biy1][15]=Fdiffthav/volume; // (16) - diffusive flux in theta
      boxvertscalars[iy-biy1][16]=Fadvthav/volume; // (17) - adv flux in theta
      boxvertscalars[iy-biy1][17]=Fradthav/volume; // (18) - total flux in theta
      boxvertscalars[iy-biy1][18]=bangleav/volume; // (19) - averaged B field angle brbphi/bsq (?)
      boxvertscalars[iy-biy1][19]=dpgasdz/volume; // (20) - 1/rho dpgas/dz 
      boxvertscalars[iy-biy1][20]=dpraddz/volume; // (21) - 1/rho dprad/dz -1. 
      boxvertscalars[iy-biy1][21]=dpmagdz/volume; // (22) - 1/rho dpmag/dz -1.
      boxvertscalars[iy-biy1][22]=Fadv2thav/volume; // (23) - <E><uth> - no buoyant
      boxvertscalars[iy-biy1][23]=Fadvrav/volume; // (24) - adv flux in radius
      boxvertscalars[iy-biy1][24]=dFadvrdrav; // (24) - divergence of radial advective flux
      boxvertscalars[iy-biy1][25]=dFtotrdrav; // (25) - divergence of radial total flux
     
      //printf("%d %f\n",iy-biy1,theta);
    } //iy-loop

#endif //BOXVERTOUTPUT==1

  return sizey;
}


/*********************************************/
/* calculates var-related scalars  */
/*********************************************/

int calc_varscalars(ldouble *varscalars,ldouble t)
{
#if(VAROUTPUT==1)


  //adjust NVARSCALARS in problem.h
  int ix,iy,iz,ii,iv;
  ldouble pp[NV]; ldouble xx[4],xxBL[4],xx1[4],xxBL1[4];
  ldouble varscalarsloc[NVARSCALARS];
  //zero scalars by default
  for(ii=0;ii<NVARSCALARS;ii++)
    varscalars[ii]=varscalarsloc[ii]=0.;


 //limits of this tile
  int rmin,rmax;
  get_xx(0,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
  rmin=xxBL[1];
  get_xx(NX-1,0,0,xx);
  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
  rmax=xxBL[1];

  //if tile completely outside the box
  int ifoutsidebox=0;

  if((rmax<VARRADIUS) || (rmin>VARRADIUS))
    ifoutsidebox=1;

  #ifdef MPI
  if(TK!=0) ifoutsidebox=1; //only slice through giz=0
#endif

  if(!ifoutsidebox)  //do the calculations only if VARRADIUS inside given tile
    {

      //search for appropriate radial index
     
      ldouble radius=VARRADIUS;

      for(ix=0;ix<NX;ix++)
	{
	  get_xx(ix,0,0,xx);
	  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
	  if(xxBL[1]>radius) break;
	}

 
      //ix fixed by VARRADIUS
      iz=0; //no azimuthal averaging - slice through iz=0, use ./phisli first, or run on the go
  
      //divide 0-2*MPI uniformly into NVARCUTS and find appropriate polar indices
      int iys[NVARCUTS],i;
      ldouble th;
      for(i=0;i<NVARCUTS;i++)
	{
	  iy=0;
	  th=M_PI*(ldouble)i/(ldouble)(NVARCUTS-1);
	  do
	    {
	      iy++;
	      get_xx(ix,iy,iz,xx);
	      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
	      get_xx(ix,iy-1,iz,xx1);
	      coco_N(xx1,xxBL1,MYCOORDS,OUTCOORDS);
	    }
	  while(!(th<=xxBL[2] && th>xxBL1[2]) && iy<=NY-1);

	  if(iy>=NY) iy=-1; //no cut within this tile

	  #ifdef MPI
	  if(TJ==0 && iy==-1 && i==0) iy=0;
	  if(TJ==NTY-1 && iy==-1 && i==NVARCUTS-1) iy=NY-1;
	  #else
	  if(iy==-1 && i==0) iy=0;
	  if(iy==-1 && i==NVARCUTS-1) iy=NY-1;
	  #endif

	  iys[i]=iy;

	  //if(PROCID==0) printf("%d > %d > %d %f %f\n",PROCID,i,iy,th,xxBL[2]);
	}

      int idx=0;
      ldouble diy;	  
      for(i=0;i<NVARCUTS;i++,idx+=NVARVARSPERCUT)
	{
	  iy=iys[i];
	  if(iy<0) continue;
	  //printf("%d\n",iy);

	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	
	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	
	  //primitives at the cell - either averaged or original, in BL or MYCOORDS
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);
		
	  //to BL, res-files and primitives in avg in MYCOORDS
	  trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);

	  //from now on - working in BL coords
        
	  //primitives and derivatives
	  /****************************/
	  /****************************/
	  /****************************/
	  ldouble rho=pp[RHO];
	  ldouble uint=pp[UU];
	  ldouble temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
	  ldouble bsq=0.;
	  ldouble ucon[4],ucov[4],bcon[4],bcov[4],vel[4];
	  ucon[1]=pp[VX];
	  ucon[2]=pp[VY];
	  ucon[3]=pp[VZ];
		  
#ifdef MAGNFIELD
	  calc_bcon_prim(pp,bcon,&geomBL);
	  indices_21(bcon,bcov,geomBL.gg); 
	  bsq = dotB(bcon,bcov); 
#endif

	  conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
	  ldouble rhouconr=rho*ucon[1];

	  ldouble Tij[4][4],Tij22[4][4],Rij[4][4];
	  calc_Tij(pp,&geomBL,Tij22);
	  indices_2221(Tij22,Tij,geomBL.gg);
	  ldouble Trt = Tij[1][0],Rrt=0.;
	  ldouble Trtmagn = bsq*ucon[1]*ucov[0] - bcon[1]*bcov[0];
	  ldouble Trtkin =  rho*ucon[1]*ucov[0];

#ifdef RADIATION
	  calc_Rij(pp,&geomBL,Rij);
	  indices_2221(Rij,Rij,geomBL.gg);
	  Rrt = Rij[1][0];
#endif

	  //PLACE - overwrite with avg quantities if required
	  if(doingavg)
	    {
	  
	    }

	  ldouble fluxconv=fluxGU2CGS(1.); 
	  varscalarsloc[idx+0]=-fluxconv*(Trt+Rrt+rhouconr); // (2nd + 2*i column) - total energy, i - # of slice
	  varscalarsloc[idx+1]=-fluxconv*Rrt; // (3rd + 2*i column) - radiative flux


	}
    }

  //aggregating over all tiles to the master who will then print out
#ifdef MPI
  MPI_Reduce(varscalarsloc, varscalars, NVARSCALARS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  for(iv=0;iv<NVARSCALARS;iv++)
    varscalars[iv]=varscalarsloc[iv];
#endif

#endif //VAROUTPUT==1
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//integrates mass in the domain

ldouble
calc_totalmass()
{
  int ix, iy, iz;
  ldouble xx[4], dx[3], mass, rho, gdet;
  
  mass = 0.;
  
  for(iz = 0; iz < NZ; iz++)
  {
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
//**********************************************************************
//**********************************************************************
//calculates the Eddington mass accretion rate

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
//**********************************************************************
//*********************************************************************
//calculates the Eddington luminosity

ldouble
calc_lumEdd()
{
  ldouble Lcgs=1.25e38*MASS; //erg/s

  return Lcgs;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates local radial fluxes of energy
//normalized to total sphere, taken at radius radius

int
calc_local_lum(int ix,int iy,int iz,ldouble *radlum, ldouble *totallum)
{
  int iv,i,j;
  ldouble xx[4],xxBL[4],dx[3],pp[NV],Rrt,rhour,Tij[4][4],Trt;
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

      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);

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
      Rrt=Rij[1][0];// + ehat*uconr);
      //	  if(Rrt<0.) Rrt=0.;
#else
      Rrt=0.;
#endif

      lum=-geomBL.gdet*Rrt*4.*M_PI;
      jet=geomBL.gdet*(Trt+rhour+Rrt)*4.*M_PI;
    }
  else
    {
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
	  
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
      Rrt=Rij[1][0];// + ehat*ucongas[1];
      //if(Rrt<0.)	  	    Rrt=0.;
#endif
     

      lum=-geom.gdet*Rrt*4.*M_PI;
      jet=geom.gdet*(rhour+Trt+Rrt)*4.*M_PI;

      //printf("%e %e %e %e\n",xxBL[1],Rrt,geom.gdet,lum);
    }

  *radlum=lum;
  *totallum=jet;

  return 0;
}


//**********************************************************************
//calculates luminosity by integrating positive flux from the axis up to tau=1 surface
//normalized to total sphere, taken at radius radius

int
calc_lum(ldouble radius,int type,ldouble *radlum, ldouble *totallum)
{


  int ix,iy,iz,iv,i,j;
  ldouble xx[4],xxBL[4],dx[3],pp[NV],Rrt,rhour,uintur,Tij[4][4],Trt;
  ldouble Rij[4][4],Rtt,ehat,ucongas[4];
  ldouble tautot[3],tau=0.;
  ldouble gdet;

 
  //search for appropriate radial index
  for(ix=0;ix<NX-1;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
      if(xxBL[1]>radius) break;
    }
  if(ix==NX) 
    {
      ix=NX-2;
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
    }
      
  

  ldouble lum=0.,jet=0.;

  if(NY==1 && NZ==1) //spherical symmetry
    {
      iz=0; 
      iy=0;

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

	  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);

	  ldouble ucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	  ldouble uconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);		  
	  rhour=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
	  uintur=get_uavg(pavg,AVGUUUCON(1),ix,iy,iz);
	  Trt=get_uavg(pavg,AVGRHOUCONUCOV(1,0),ix,iy,iz)
	    + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(1,0),ix,iy,iz)
	    + get_uavg(pavg,AVGBSQUCONUCOV(1,0),ix,iy,iz)
	    - get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz); 

#ifdef RADIATION

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
#else
	  Rrt=0.;
#endif

	
#ifdef RADIATION
	  lum=-geom.gdet*Rrt*4.*M_PI;
#else 
	  lum=-geom.gdet*uintur*4.*M_PI;    
#endif
	  jet=geomBL.gdet*(Trt+rhour+Rrt)*4.*M_PI;
	}
      else
	{
	  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
	  
	  ucongas[1]=pp[2];
	  ucongas[2]=pp[3];
	  ucongas[3]=pp[4];	      
	  conv_vels(ucongas,ucongas,VELPRIM,VEL4,geom.gg,geom.GG);

	  rhour = pp[RHO]*ucongas[1];
	  uintur = pp[UU]*ucongas[1];
	  
	  calc_Tij(pp,&geom,Tij);
	  indices_2221(Tij,Tij,geom.gg);
	  Trt=Tij[1][0];


#ifdef RADIATION	      
	  if(type==0) //R^r_t outside photosphere
	    {
	      Rrt=0.;
	    }
	  else if(type==1) //sum of positive R^r_t everywhere
	    {
	      //calc_ff_Rtt(pp,&Rtt,ucongas,&geom);
	      //ehat=-Rtt;
	      calc_Rij(pp,&geom,Rij); 
	      //indices_2221(Rij,Rij,geom.gg);
	      Rrt=Rij[1][0];// + ehat*ucongas[1];
	      if(Rrt<0.)
		Rrt=0.;
	    }
	  else if(type==2) //R^r_t in the outflow region
	    {
	      //calc_ff_Rtt(pp,&Rtt,ucongas,&geom);
	      //ehat=-Rtt;
	      calc_Rij(pp,&geom,Rij); 
	      //	      indices_2221(Rij,Rij,geom.gg);
	      Rrt=Rij[1][0];// + ehat*ucongas[1];
	      if(Rrt<0. || ucongas[1]<0.)
		Rrt=0.;
	    }
	  else if(type==3) //sum of R^r_t everywhere
	    {
	      calc_Rij(pp,&geom,Rij); 
	      //indices_2221(Rij,Rij,geom.gg);
	      Rrt=Rij[1][0];// + ehat*ucongas[1];	      
	    }
	  else
	    Rrt=0.;
#else
	  Rrt=0.;
#endif

#ifdef RADIATION
	  lum=-geom.gdet*Rrt*4.*M_PI;
#else 
	  lum=-geom.gdet*uintur*4.*M_PI;    
#endif

	  jet=geom.gdet*(rhour+Trt+Rrt)*4.*M_PI;

	  //printf("%e %e %e %e\n",xxBL[1],Rrt,geom.gdet,lum);
	}

      *radlum=lum;
      *totallum=jet;
      return 0.;
    }
  else //non-sph symmetry
    {
      for(iz=0;iz<NZ;iz++)
      for(iy=0;iy<NY;iy++)
	{
	  for(iv=0;iv<NV;iv++)
	    pp[iv]=get_u(p,iv,ix,iy,iz);

	  struct geometry geomBL;
	  fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  get_xx(ix,iy,iz,xx);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=2.*M_PI;
	  gdet=geom.gdet;
	  ldouble dxph[3],dxBL[3];
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

	      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
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

#ifdef RADIATION
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
		  
		  Rrt=-Rij[1][0];// + ehat*uconr);
		  if(Rrt<0.) Rrt=0.;

		  
		}
	      else if(type==2) //R^r_t everywhere in outflow
		{
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
		  Rrt=-Rij[1][0];// + ehat*uconr);
		  if(uconr<0. || Rrt<0.) Rrt=0.;
		}
	      else if(type==3) //any R^r_t everywhere
		{
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		  
		  Rrt=-Rij[1][0];// + ehat*uconr);
		  
		}
	      else
		Rrt=0.;
#else
	      Rrt=0.;
#endif

#ifdef RADIATION
	      lum+=geomBL.gdet*Rrt*dxBL[1]*dxBL[2];
	      #else 
	      lum+=geomBL.gdet*uintur*dxBL[1]*dxBL[2];
	      #endif

	      jet+=geomBL.gdet*(rhour+Rrt+Rrt)*dxBL[1]*dxBL[2];
	    }
	  else //snapshot
	    { 
	      
	      //to BL
	      trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);
	      //hydro part may be insonsistent

	      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
	      calc_tautot(pp,&geomBL,dxph,tautot);

	      ucongas[1]=pp[2];
	      ucongas[2]=pp[3];
	      ucongas[3]=pp[4];	      
	      conv_vels(ucongas,ucongas,VELPRIM,VEL4,geomBL.gg,geomBL.GG);

	      rhour = pp[RHO]*ucongas[1];
	      uintur = pp[UU]*ucongas[1];
	  
	      calc_Tij(pp,&geomBL,Tij);
	      indices_2221(Tij,Tij,geomBL.gg);
	      Trt=Tij[1][0];

	      tau+=ucongas[0]*tautot[1];
	      
#ifdef RADIATION
	      if(type==0) //R^r_t outside photosphere
		{
		  if(tau>1.) break;	  
		  //trans_prad_coco(pp,pp,MYCOORDS,KERRCOORDS,xx,&geom,&geomBL);
		  //prad_lab2on(pp,pp,&geomBL);
		  calc_Rij(pp,&geomBL,Rij); 
		  indices_2221(Rij,Rij,geomBL.gg);
		  Rrt=-Rij[1][0];
		  if(Rrt<0.) Rrt=0.;
		}
	      else if(type==1) //R^r_t everywhere
		{
		  //calc_ff_Rtt(pp,&Rtt,ucongas,&geomBL);
		  //ehat=-Rtt;
		  calc_Rij(pp,&geomBL,Rij); 
		  indices_2221(Rij,Rij,geomBL.gg);
		  Rrt=-Rij[1][0];// + ehat*ucongas[1];
		  if(Rrt<0.) Rrt=0.;
		}
	      else if(type==2) //R^r_t in the outflow region
		{
		  //calc_ff_Rtt(pp,&Rtt,ucongas,&geomBL);
		  //ehat=-Rtt;
		  calc_Rij(pp,&geomBL,Rij); 
		  indices_2221(Rij,Rij,geomBL.gg);
		  Rrt=-Rij[1][0];// + ehat*ucongas[1];
		  if(Rrt<0. || ucongas[1]<0.)
		    Rrt=0.;
		}
	      else
		Rrt=0.;
#else
	      Rrt=0.;
#endif

	      jet+=geomBL.gdet*(rhour+Trt+Rrt)*dxBL[1]*dxBL[2];

	      #ifdef RADIATION
	      lum+=geomBL.gdet*Rrt*dxBL[1]*dxBL[2];
	      #else 
	      lum+=geomBL.gdet*uintur*dxBL[1]*dxBL[2];
	      #endif
	    }

	  //#ifdef CGSOUTPUT
	  //never!
	  //Rrt=fluxGU2CGS(Rrt);
	  //dx[1]=lenGU2CGS(dx[1]);
	  //dx[2]=lenGU2CGS(dx[2]);
	  //#endif		  


	  //printf("%e %e %e -> %e\n",Rrt,gdet,dx[1],dx[2],lum);getch();
	}
      
      *radlum=lum;
      *totallum=jet;
      return 0.;
    }

    return -1;
}


//**********************************************************************
//calculates luminosity escaping through radial and polar boundaries

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
//**********************************************************************
//**********************************************************************
//calculates MRI resolution parameter Q_theta ar rmri

ldouble
calc_resmri(ldouble radius)
{
#ifndef BHDISK_PROBLEMTYPE
  return -1.; //no MRI
#endif

#ifdef MAGNFIELD

  int ix,iy,iz,iv;
  ldouble xx[4],xxBL[4],dx[3];
  ldouble qtheta=0.,qphi=0.,sigma=0.,rho;
 
  //search for appropriate radial index
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
      if(xxBL[1]>radius) break;
    }

  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      {
	dx[1]=get_size_x(iy,1);
	rho=get_u(p,RHO,ix,iy,iz);

	sigma+=rho*dx[1];
	ldouble q1,q2;
	calc_Qthetaphi(ix,iy,iz,&q1,&q2);
	qtheta+=rho*q1*dx[1];
	qphi+=rho*q2*dx[1];
      }

  return qtheta/sigma;
  

#endif
    return -1.;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates mean temperature at rmri

ldouble
calc_meantemp(ldouble radius)
{
#ifndef BHDISK_PROBLEMTYPE
  return -1.; //no disk no cry
#endif


  int ix, iy, iz, iv;
  ldouble xx[4], xxBL[4], dx[3];
  ldouble mtemp = 0., sigma = 0., rho, ugas;
 
  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    get_xx(ix, 0, 0, xx);
    coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
    if(xxBL[1] > radius) break;
  }

  for (iz = 0; iz < NZ; iz++)
  {
    for(iy = 5; iy < NY-5; iy++)
    {
      dx[1] = get_size_x(iy, 1);
      rho = get_u(p, RHO, ix, iy, iz);
      ugas = get_u(p, UU, ix, iy, iz);
      ldouble temp = calc_PEQ_Tfromurho(ugas, rho, ix, iy, iz);
      sigma += rho * dx[1];
      mtemp += rho * temp*dx[1];
    }
  }
  
  return mtemp/sigma;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates mean temperature at rmri

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
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
      if(xxBL[1]>radius) break;
    }

  return scaleth_otg[ix];
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates theta corresponding to integrated tau from the axis

ldouble
calc_photloc(int ix)
{
  if(MYCOORDS != OUTCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS && MYCOORDS != MKS2COORDS && MYCOORDS != MKS3COORDS)
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

	  get_xx(ix,iy,iz,xx);
	  dx[0]=get_size_x(ix,0);
	  dx[1]=get_size_x(iy,1);
	  dx[2]=get_size_x(iz,2);
	  dx[0]=dx[0]*sqrt(geom.gg[1][1]);
	  dx[1]=dx[1]*sqrt(geom.gg[2][2]);
	  dx[2]=2.*M_PI*sqrt(geom.gg[3][3]);

	  coco_N(xx,xxBL,MYCOORDS,OUTCOORDS);
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
//**********************************************************************
//**********************************************************************
//calculates rest mass flux through r=radius within range of thetas
//normalized to 2pi in phi
//type == 0 (net)
//type == 1 (inflow only)
//type == 2 (outflow only)

ldouble
calc_mdot(ldouble radius, int type)
{
  int ix, iy, iz, iv;
  ldouble xx[4], xxBL[4], dx[3], mdot, gdet, rho, rhouconr, ucon[4], pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5];

  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    get_xx(ix, 0, 0, xx);
    coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
    if(xxBL[1] > radius) break;
  }
  
  mdot = 0.;
  
  for(iz = 0; iz < NZ; iz++)
  {
    for(iy = 0; iy < NY; iy++)
    {
      struct geometry geom;
      fill_geometry_arb(ix, iy, iz, &geom, MYCOORDS);
      
      struct geometry geomBL;
      fill_geometry_arb(ix, iy, iz, &geomBL, OUTCOORDS);
      
      if(doingavg)
      {
        //rho = get_uavg(pavg, RHO, ix, iy, iz);
        //ucon[1] = get_uavg(pavg, AVGRHOUCON(1), ix, iy, iz)/get_uavg(pavg, RHO, ix, iy, iz);
        //ucon[2] = get_uavg(pavg, AVGRHOUCON(2), ix, iy, iz)/get_uavg(pavg, RHO, ix, iy, iz);
        //ucon[3] = get_uavg(pavg, AVGRHOUCON(3), ix, iy, iz)/get_uavg(pavg, RHO, ix, iy, iz);
        rhouconr = get_uavg(pavg, AVGRHOUCON(1), ix, iy, iz);
        gdet = geomBL.gdet;
        ldouble xx1[4], xx2[4];
        //xx1[0] = 0.; xx1[1] = get_xb(ix, 0); xx1[2] = get_x(iy, 1); xx1[3] = get_x(iz, 2);
        //xx2[0] = 0.; xx2[1] = get_xb(ix+1, 0);xx2[2] = get_x(iy, 1); xx2[3] = get_x(iz, 2);
        //coco_N(xx1, xx1, MYCOORDS, OUTCOORDS);
        //coco_N(xx2, xx2, MYCOORDS, OUTCOORDS);
        //dx[0] = fabs(xx2[1] - xx1[1]);
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
        
        get_xx(ix, iy, iz, xx);
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
        
        coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
        
        
        rho = pp[0];
        ucon[1] = pp[2];
        ucon[2] = pp[3];
        ucon[3] = pp[4];
        
        conv_vels(ucon, ucon, VELPRIM, VEL4, geom.gg, geom.GG);
        rhouconr = rho * ucon[1];
        gdet = geom.gdet;
      }
      
      if(NY==1)
      {
        dx[1] = 2.;
        dx[2] = 2. * M_PI;
      }
      
      if(type==0 || (type==1 && ucon[1]<0.) || (type==2 && ucon[1]>0.))
        mdot += gdet * rhouconr * dx[1] * dx[2];
    }
  }
  
  return mdot;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates the luminosity proxy for the EHT code comparison project
//
//

ldouble
calc_lum_proxy(ldouble radius, ldouble theta_min, ldouble theta_max)
{
  int ix, iy, iz, iv;
  ldouble Edot, xx[4], xxBL[4], dx[3], gdet, pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5], rho, uu, pressure, bsq, bfield, ucon[4], ucov[4], bcon[4], bcov[4], luminosity;
  
  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    get_xx(ix, 0, 0, xx);
    coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
    if(xxBL[1] > radius) break;
  }
  int ixmin = ix;

  luminosity = 0.;
  
  for(iz = 0; iz < NZ; iz++)
  {
    for(iy = 0; iy < NY; iy++)
    {
      for(ix = ixmin; ix < NX; ix++)
      {
        get_xx(ix, iy, iz, xx);
        coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
      
        // Check that theta lies within the required range
        if (xxBL[2] >= theta_min && xxBL[2] <= theta_max)
        {
          struct geometry geom;
          fill_geometry_arb(ix, iy, iz, &geom, MYCOORDS);
          
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
//**********************************************************************
//**********************************************************************
//calculates energy flux through r=radius:
//
//   Edot = - Integral_over_theta_phi[(rho + u + p + bsq) u^r u_t - b^r b_t]
//

ldouble
calc_Edot(ldouble radius)
{
  int ix, iy, iz, iv;
  ldouble Edot, xx[4], xxBL[4], dx[3], gdet, rhouconrucovt, uuuconrucovt, bsquconrucovt, bconrbcovt, pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5], ucon[4], ucov[4], bcon[4], bcov[4], rho, uu, bsq;
  
  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    get_xx(ix, 0, 0, xx);
    coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
    if(xxBL[1] > radius) break;
  }
  
  Edot = 0.;
  
  for(iz = 0; iz < NZ; iz++)
  {
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
        
        get_xx(ix, iy, iz, xx);
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
        
        coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
        
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
//**********************************************************************
//**********************************************************************
//calculates angular momentum flux through r=radius:
//
//   Ldot = Integral_over_theta_phi[(rho + u + p + bsq) u^r u_phi - b^r b_phi]
//

ldouble
calc_Ldot(ldouble radius)
{
  int ix, iy, iz, iv;
  ldouble Ldot, xx[4], xxBL[4], dx[3], gdet, rhouconrucovphi, uuuconrucovphi, bsquconrucovphi, bconrbcovphi, pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5], ucon[4], ucov[4], bcon[4], bcov[4], rho, uu, bsq;
  
  //search for appropriate radial index
  for(ix = 0; ix < NX; ix++)
  {
    get_xx(ix, 0, 0, xx);
    coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
    if(xxBL[1] > radius) break;
  }
  
  Ldot = 0.;
  
  for(iz = 0; iz < NZ; iz++)
  {
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
        
        get_xx(ix, iy, iz, xx);
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
        
        coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
        
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
//**********************************************************************
//**********************************************************************
//calculates magnetic flux through r=radius within range of thetas
//normalized to 2pi in phi
//type == 0 (default)

int
calc_Bflux(ldouble radius, int type, ldouble *Bflux, ldouble* Bfluxquad)
{
  *Bflux = *Bfluxquad = 0.;
  
  if(MYCOORDS != OUTCOORDS && MYCOORDS != KSCOORDS && MYCOORDS != MKS1COORDS && MYCOORDS != MKS2COORDS && MYCOORDS != MKS3COORDS)
  {
    return -1.; //no BH
  }

#ifdef MAGNFIELD

  int ix, iy, iz, iv, i4;
  ldouble xx[4], xxBL[4], dx[3], Psi, Psiquad, rho, ucon[4], pp[NV], gg[4][5], GG[4][5], ggBL[4][5], GGBL[4][5];

  //search for appropriate radial index ix corresponding to the required radius
  for(ix = 0; ix < NX; ix++)
  {
    get_xx(ix, 0, 0, xx);
    coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
    if(xxBL[1] > radius) break;
  }

  Psi = 0.;  // dipole flux
  Psiquad = 0.;  // quadrupole flux

//  We have changed the following so that it can handle both NZ=1 and NZ>1
  for (iz = 0; iz < NZ; iz++)
  {
    for(iy=0;iy<NY;iy++)
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
        if(alphanorm < 0.) my_err("alpha.lt.0 in b0 norm !!\n");
        for(i4 = 0; i4 < 4; i4++)
        {
          bcon[i4] *= sqrt(alphanorm);
        }
        
        //Bcon[1]
        ldouble Br = bcon[1] * ucon[0] - bcon[0] * ucon[1];
        
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
        
        get_xx(ix, iy, iz, xx);
        coco_N(xx, xxBL, MYCOORDS, OUTCOORDS);
        
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


/*********************************************/
//calculates quantites with respect to the averaged solution
//and projects them on the equatorial plane
/*********************************************/

int calc_anarelradialprofiles(ldouble profiles[][NX])
{
  //adjust NRADPROFILES in problem.h

  int ix,iy,iz,iv,i,j;
  ldouble rhoavg,uintavg,tempavg,uconavg[4],utconavg[4],bconavg[4],ppavg[NV],bsqavg;
  ldouble rho,uint,temp,ucon[4],ucov[4],bcon[4],bcov[4],bsq;
  ldouble Tij[4][4],alpha,pp[NV],Rij[4][4],Ehat,urcon[4];

  //loop over radius
  for(ix=0;ix<NX;ix++)
    {
      //vertically integrated/averaged profiles

      for(iv=0;iv<NANARELRADPROFILES;iv++)
	profiles[iv][ix]=0.;
      
 #ifdef BHDISK_PROBLEMTYPE
     if(NZ==1) //phi-symmetry
	{
	  iz=0;
	  for(iy=0;iy<NY;iy++)
	    {
	      struct geometry geom;
	      fill_geometry_arb(ix,iy,iz,&geom,MYCOORDS);

	      struct geometry geomBL;
	      fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
	      
	      ldouble dxph[3],dx[3];
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

	      dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
	      dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
	      dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);
	      
	      /*******************************/
	      //averaged quantities (always VEL4 BL!)
	      rhoavg=get_uavg(pavg,RHO,ix,iy,iz);
	      uintavg=get_uavg(pavg,UU,ix,iy,iz);
	      bsqavg=get_uavg(pavg,AVGBSQ,ix,iy,iz);
	      bconavg[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
	      bconavg[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
	      bconavg[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
	      bconavg[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
	      uconavg[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      uconavg[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      uconavg[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      uconavg[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	      conv_vels(uconavg,utconavg,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
	      ppavg[RHO]=rhoavg;
	      ppavg[UU]=uintavg;
	      ppavg[VX]=utconavg[1];
	      ppavg[VY]=utconavg[2];
	      ppavg[VZ]=utconavg[3];
	      //magnetic field etc. empty!

	      /*******************************/
	      //snapshot quantities (always VEL4 BL!)
	      //use old ones instead
	      for(iv=0;iv<NV;iv++)
		  pp[iv]=get_u(p,iv,ix,iy,iz);
	      //to BL     
	      trans_pall_coco(pp,pp,MYCOORDS,OUTCOORDS,geom.xxvec,&geom,&geomBL);
	      rho=pp[0];
	      uint=pp[1];
	      ucon[1]=pp[2];
	      ucon[2]=pp[3];
	      ucon[3]=pp[4];
	      conv_vels_both(ucon,ucon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
		  
              #ifdef MAGNFIELD
	      calc_bcon_prim(pp,bcon,&geomBL);
	      indices_21(bcon,bcov,geomBL.gg); 
	      bsq = dotB(bcon,bcov); 
              #endif

              #ifdef RADIATION
	      calc_Rij(pp,&geomBL,Rij);
	      indices_2221(Rij,Rij,geomBL.gg);
	      calc_ff_Rtt(pp,&Ehat,urcon,&geomBL);
	      Ehat*=-1.;
              #endif


	      /*******************************/
	      //calculating alpha
	      
	      //stress energy tensor in lab frame
	      calc_Tij(pp,&geomBL,Tij);

	      //boosting it from lab frame to the average comoving frame
	      //what if there are secular trends? 
	      //one should take avg from around snapshot file
	      //or do as Penna+12 did, i.e., take averaged phi velocity=VZ
	      boost22_lab2ff(Tij,Tij,pp,geomBL.gg,geomBL.GG);

	      //pressure
	      ldouble ptot = GAMMAM1*uint + 1./3.*Ehat;
	      alpha=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij[1][3]/ptot;
	      
	      /*******************************/
	      //vertical integrals
	      //within scale-height
	      if(fabs(geomBL.yy-M_PI/2.) < scaleth_otg[ix])
		{
		  //surface density (2nd column)
		  profiles[2-2][ix]+=rho*dxph[1];

		  //alpha numerator (3) (sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij[1][3]) 
		  profiles[3-2][ix]+=sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*Tij[1][3]*rho*dxph[1];
		  
		  //alpha denominator (4) (ptot)
		  profiles[4-2][ix]+=ptot*rho*dxph[1];
		  
		  //direct alpha (5) 
		  profiles[5-2][ix]+=alpha*rho*dxph[1];
		}

	    }

	  //normalizing by sigma
	  profiles[1][ix]/=profiles[0][ix];
	  profiles[2][ix]/=profiles[0][ix];
	  profiles[3][ix]/=profiles[0][ix];
	 
	}

#endif
    }

  return 0;
}
