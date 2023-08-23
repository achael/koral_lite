/*! \file silo.c
 \brief Routines for SILO output files
 */

#ifndef NOSILO

//KORAL - silo.c
//routines for writing a silo file with quadratic mesh
//used both on the go and separately

#include "ko.h"
#include <silo.h>
#include <string.h>

/*********************************************/
/* writes silo file in dumps */
/*********************************************/
int fprint_silofile(ldouble time, int num, char* folder, char* prefix)
{

  char bufor[50];
  sprintf(bufor,"%s/%s%04d.silo",folder,prefix,num);

  mpi_exchangedata();
  calc_avgs_throughout();

  DBfile *file = NULL;    // The Silo file pointer 
  char *coordnames[3];    // Names of the coordinates 
  ldouble *nodex;         // The coordinate arrays 
  ldouble *nodey;
  ldouble *nodez;
  ldouble *coordinates[3];// The array of coordinatearrays 
  int dimensions[3];      // The number of nodes 

  // Create the Silo file 
  file = DBCreate(bufor, DB_CLOBBER, DB_LOCAL, NULL,DB_PDB);

  // Name the coordinate axes ‘X’ and ‘Y’ 
  coordnames[0] = strdup("X");
  coordnames[1] = strdup("Y");
  coordnames[2] = strdup("Z");
  
  // Give the cartesian coordinates of the mesh 
  int ix,iy,iz,iv,imx,imy,imz;
  int i,j;
  ldouble pp[NV],uu[NV],xxvec[4],xxveccar[4],xxvecsph[4],xx1[4],xx2[4];
  
  //number of zones
  int nx=NX;
  int ny=NY;
  int nz=NZ;
  
  //number of nodes
  int nnx=NX+1;
  int nny=NY+1;
  int nnz=NZ+1;
  if(NY==1) nny=NY;
  if(NZ==1) nnz=NZ;
  #ifdef FULLPHI //printing one more cell in phi to close the sphere
  nz++;nnz++;
  #endif

  // Allocate arrays for saved quantities
  nodex=(ldouble *)malloc(nnx*nny*nnz*sizeof(ldouble));
  nodey=(ldouble *)malloc(nnx*nny*nnz*sizeof(ldouble));
  nodez=(ldouble *)malloc(nnx*nny*nnz*sizeof(ldouble));

  ldouble *rho = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *entr = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *uint = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  ldouble *temp = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Omega = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *muBe = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Be = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  int *entropyinv = (int*)malloc(nx*ny*nz*sizeof(int));
  ldouble *expansion = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *NHr = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *entrlnT = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *u0 = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *u1 = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *u2 = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *u3 = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *lorentz = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  ldouble *lorentz_perp = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *lorentz_par = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *ffinv = (ldouble*)malloc(nx*ny*nz*sizeof(double));
    
  ldouble *vx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vy = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vz = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  ldouble *Edotx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Edoty = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Edotz = (ldouble*)malloc(nx*ny*nz*sizeof(double));


  #ifdef MAGNFIELD
  ldouble *Bangle = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *bsq = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *B1comp = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *B2comp = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *B3comp = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Bx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *By = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Bz = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *phi = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *phi_accum = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *beta = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *sigma = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *sigmaw = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *betainv = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *divB = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *jdens = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *OmegaF = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Qtheta = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Qphi = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  
  #ifdef BATTERY
  ldouble *dBxdtbat = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *dBydtbat = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *dBzdtbat = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *phibat = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif

  #ifdef MIMICDYNAMO  
  if(doingpostproc) 
  {
    set_bc(time,0);
    mimic_dynamo(1.); 
  }
  ldouble *Bxdyn = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Bydyn = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Bzdyn = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *phidyn = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Avecdyn = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif
  #endif //MAGNFIELD

  #ifdef PRINTVISCHEATINGTOSILO
  ldouble *dtauarr= (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vischeat= (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *vischeatnege= (ldouble*)malloc(nx*ny*nz*sizeof(double));  
  ldouble *vischeatnegi= (ldouble*)malloc(nx*ny*nz*sizeof(double));  
  ldouble *deltae= (ldouble*)malloc(nx*ny*nz*sizeof(double));  
  #endif

  #ifdef PRINTCOULOMBTOSILO
  ldouble *coulomb = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif
 
  #ifdef PRINTGAMMATOSILO
  ldouble *gammag = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif
  
  #ifdef EVOLVEELECTRONS
  ldouble *tempe = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *tempi = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *ue = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *ui = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif 

  #ifdef RELELECTRONS
  ldouble *gammabrk = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *G0relel = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *G0icrelel = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *G0synrelel = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *urelel = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *uratio_tot = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *uratio_th = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *nrelel = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *neth = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  int ie;
  ldouble gammapbrk[NRELBIN]; //shared gammabreak array for nonthermal  electrons
  for(ie=0; ie<NRELBIN; ie++) gammapbrk[ie] = pow(relel_gammas[ie], RELEL_HEAT_INDEX + 0.5);
  #endif  

  #ifdef RADIATION
  ldouble *trad = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *tradlte = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  ldouble *taucoupling = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *tausca = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *tauabs = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *taueff = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Erad = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Ehat = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Fx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Fy = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Fz = (ldouble*)malloc(nx*ny*nz*sizeof(double));  
  ldouble *Nph = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *uradx = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *urady = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *uradz = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *Gtff =  (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *G0icth = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  ldouble *tauscar = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *tauabsr = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *taueffr = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  ldouble *tauscar2 = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *tauabsr2 = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  ldouble *taueffr2 = (ldouble*)malloc(nx*ny*nz*sizeof(double));

  ldouble *entrrad = (ldouble*)malloc(nx*ny*nz*sizeof(double));
  #endif //RADIATION

    
  //first fill coordinates on nodes
#pragma omp parallel for private(ix,iy,iz,iv,imx,imy,imz,i,j,pp,uu,xxvec,xxveccar,xxvecsph,xx1,xx2) schedule (static)

#ifdef PRINTZGC_RIGHT
  for(iz=NG;iz<=nz+NG;iz++)
#else
  for(iz=0;iz<=nz;iz++)
#endif
  {

#ifdef PRINTYGC_RIGHT
    for(iy=NG;iy<=ny+NG;iy++)
#else
    for(iy=0;iy<=ny;iy++)
#endif
    {
      
#if defined(PRINTXGC_RIGHT)
      for(ix=NG;ix<=nx+NG;ix++)
#elif defined(PRINTXGC_LEFT)
      for(ix=-NG;ix<=nx-NG;ix++)
#else
      for(ix=0;ix<=nx;ix++)
#endif
      {

        if(NZ==1 && iz>0) continue;
	if(NY==1 && iy>0) continue;

	int iix,iiy,iiz;
	iix=ix;
	iiy=iy;
	iiz=iz;

	if(NZ>1 && NY>1)
	{
	    xxvec[0]=0.;xxvec[1]=get_xb(iix,0);xxvec[2]=get_xb(iiy,1);xxvec[3]=get_xb(iiz,2);
	}
	else if(NZ==1 && NY==1)
	{
	    xxvec[0]=0.;xxvec[1]=get_xb(iix,0);xxvec[2]=get_x(iiy,1);xxvec[3]=get_x(iiz,2);
	}
	else if(NZ==1)
	{
	    xxvec[0]=0.;xxvec[1]=get_xb(iix,0);xxvec[2]=get_xb(iiy,1);xxvec[3]=get_x(iiz,2);
	}
	else if(NY==1)
	{
	    xxvec[0]=0.;xxvec[1]=get_xb(iix,0);xxvec[2]=get_x(iiy,1);xxvec[3]=get_xb(iiz,2);
	}

	coco_N(xxvec,xxvecsph,MYCOORDS,SPHCOORDS);
	coco_N(xxvec,xxveccar,MYCOORDS,MINKCOORDS);

	imz=iz;
	imy=iy;
	imx=ix;
        #ifdef PRINTZGC_RIGHT
	imz=iz-NG;
        #endif
        #ifdef PRINTXGC_RIGHT
	imx=ix-NG;
        #endif
        #ifdef PRINTXGC_LEFT
	imx=ix+NG;
        #endif
        #ifdef PRINTYGC_RIGHT
	imy=iy-NG;
        #endif

	int nodalindex=imz*(nny*nnx) + imy*nnx + imx;

	//coordinates

	ldouble coordscale=1.;
	#ifdef COORDSINPCINSILO
	coordscale=1./(PARSECCGS/MASSCM);
	#endif

#if(SILOCOORDS==SPHCOORDS)
	nodex[nodalindex]=xxvecsph[1]*coordscale;
	nodey[nodalindex]=xxvecsph[2]*coordscale;
	nodez[nodalindex]=xxvecsph[3]*coordscale;
#else
	nodex[nodalindex]=xxveccar[1]*coordscale;
	nodey[nodalindex]=xxveccar[2]*coordscale;
	nodez[nodalindex]=xxveccar[3]*coordscale;
#endif	      
      }
    }
  }


  //then fill the zones with values
#pragma omp parallel for private(ix,iy,iz,iv,imx,imy,imz,i,j,pp,uu,xxvec,xxveccar,xxvecsph,xx1,xx2) schedule (static)

#ifdef PRINTZGC_RIGHT
  for(iz=NG;iz<nz+NG;iz++)
#else
  for(iz=0;iz<nz;iz++)
#endif
  {

#ifdef PRINTYGC_RIGHT
    for(iy=NG;iy<ny+NG;iy++)
#else
    for(iy=0;iy<ny;iy++)
#endif
    {
      
#if defined(PRINTXGC_RIGHT)
    for(ix=NG;ix<nx+NG;ix++)
#elif defined(PRINTXGC_LEFT)
    for(ix=-NG;ix<nx-NG;ix++)
#else
    for(ix=0;ix<nx;ix++)
#endif
    {

      int iix,iiy,iiz;
      iix=ix;
      iiy=iy;
      iiz=iz;

      struct geometry geom;
      fill_geometry(iix,iiy,iiz,&geom);
      struct geometry geomout;
      fill_geometry_arb(iix,iiy,iiz,&geomout,OUTCOORDS);

      //cell dimensions
      ldouble dxph[3],dx[3];
      get_cellsize_out(ix, iy, iz, dx);

      dxph[0]=dx[0]*sqrt(geomout.gg[1][1]);
      dxph[1]=dx[1]*sqrt(geomout.gg[2][2]);
      dxph[2]=dx[2]*sqrt(geomout.gg[3][3]);

      // cell coordinates
      get_xx(iix,iiy,iiz,xxvec);
      coco_N(xxvec,xxvecsph,MYCOORDS,SPHCOORDS);
      coco_N(xxvec,xxveccar,MYCOORDS,MINKCOORDS);

      ldouble r=xxvecsph[1];
      ldouble th=xxvecsph[2];
      ldouble ph=xxvecsph[3];
      
      // metric determinant
      ldouble gdet,gdetu;
      gdet=geomout.gdet;
      gdetu=gdet;
      #if (GDETIN==0) //gdet out of derivatives
      gdetu=1.;
      #endif

      // coordinates and gdet of cells +/- 1 in radius
      ldouble gdet1,gdet2;
      #ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix-1, iiy, iiz, xx1);
      get_xxout(ix+1, iiy, iiz, xx2);	      
      #else
      get_xx(iix-1,iiy,iiz,xx1);
      get_xx(iix+1,iiy,iiz,xx2);
      coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
      coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
      #endif

      gdet1=calc_gdet_arb(xx1,OUTCOORDS);
      gdet2=calc_gdet_arb(xx2,OUTCOORDS);

      // compute zonal index
      imz=iz;
      imy=iy;
      imx=ix;
      #ifdef PRINTZGC_RIGHT
      imz=iz-NG;
      #endif
      #ifdef PRINTXGC_RIGHT
      imx=ix-NG;
      #endif
      #ifdef PRINTXGC_LEFT
      imx=ix+NG;
      #endif
      #ifdef PRINTYGC_RIGHT
      imy=iy-NG;
      #endif

      int zonalindex=imz*(ny*nx) + imy*nx + imx;

      // copy over primitives
      for(iv=0;iv<NV;iv++)
      {
	if(doingavg)
	  pp[iv]=get_uavg(pavg,iv,ix,iy,iz); 
	else
	  pp[iv]=get_u(p,iv,iix,iiy,iiz);
      }

      // get (or recompute) gamma
      ldouble gamma=GAMMA;
      #ifdef CONSISTENTGAMMA
      #ifdef EVOLVEELECTRONS
      gamma=calc_gammagas(pp,ix,iy,iz);
      set_u_scalar(gammagas,ix,iy,iz,gamma);
      #endif
      gamma=pick_gammagas(ix,iy,iz);
      #endif

      // gas expansion
      ldouble exploc=0.;
      int derdir[3] = {0,0,0};
      ldouble shear[4][4];
      calc_shear_lab(pp, &geom, shear, &exploc, MHD, derdir);

      //primitives to OUTCOORDS
      #ifdef PRECOMPUTE_MY2OUT
      trans_pall_coco_my2out(pp,pp,&geom,&geomout);
      #else      
      trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);
      #endif

      // entropy inversion flag
      entropyinv[zonalindex]=get_cflag(ENTROPYFLAG3,ix,iy,iz);
	    
      // other quantities 
      ldouble vel[4],vcov[4],vcon[4],velprim[4];
      ldouble Tit[4],Tij[4][4];
      ldouble rhouconr,rhoucont;
      ldouble tempeloc,tempiloc;
      ldouble pregas,premag,ueloc,uiloc;
      ldouble nethloc,nrelelloc,G0relelloc,G0ic_relel_loc,G0syn_relel_loc,G0ic_th_loc,urelelloc,gbrkloc;
      ldouble bcon[4],bcov[4];
      
      if(doingavg==0) //using snapshot data
      {
	
	rho[zonalindex]  = pp[RHO];
	entr[zonalindex] = pp[ENTR];
	uint[zonalindex] = pp[UU];
	pregas = (gamma-1.)*pp[UU];

	vel[0] = 0.;
	vel[1]=pp[VX];
	vel[2]=pp[VY];
	vel[3]=pp[VZ];

	conv_vels(vel,vel,VELPRIM,VEL4,geomout.gg,geomout.GG);
	for(i=0;i<4;i++) vcon[i]=vel[i];
	indices_21(vcon,vcov,geomout.gg);

	rhouconr = rho[zonalindex]*vcon[1];
	rhoucont = rho[zonalindex]*vcon[0];

	// angular velocity
	Omega[zonalindex]=vel[3]/vel[0];

	#ifdef MAGNFIELD
	B1comp[zonalindex] = pp[B1];
        B2comp[zonalindex] = pp[B2];
        B3comp[zonalindex] = pp[B3];

        calc_bcon_prim(pp,bcon,&geomout);
	indices_21(bcon,bcov,geomout.gg); 
	bsq[zonalindex] = dotB(bcon,bcov);

	premag = 0.5*bsq[zonalindex];
	beta[zonalindex] = pregas/premag;
	betainv[zonalindex] = 1./beta[zonalindex];
	sigma[zonalindex] = bsq[zonalindex]/(rho[zonalindex]);
	sigmaw[zonalindex] = bsq[zonalindex]/(rho[zonalindex] + gamma*uint[zonalindex]);

        #endif		  

	// stress energy
	calc_Tij(pp,&geomout,Tij);
	indices_2221(Tij,Tij,geomout.gg);

	Tit[1]=Tij[1][0];
	Tit[2]=Tij[2][0];
	Tit[3]=Tij[3][0];

	// Bernoulli flux & Bernoulli number
	muBe[zonalindex] = -(Tij[1][0] + rhouconr)/rhouconr;
	Be[zonalindex] = -(Tij[0][0] + rhoucont)/rhoucont;

#ifdef MAGNFIELD		  
	Qtheta[zonalindex]=2.*M_PI/Omega[zonalindex]/dx[1]*fabs(bcon[2])/sqrt(rho[zonalindex]);
	Qphi[zonalindex]=2.*M_PI/Omega[zonalindex]/dx[2]*fabs(bcon[3])/sqrt(rho[zonalindex]);

	#ifdef SHEARINGBOX
	Qtheta[zonalindex]=2.*M_PI/SHEAROM/dx[2]*fabs(bcon[3])/sqrt(rho[zonalindex]);
	Qphi[zonalindex]=2.*M_PI/SHEAROM/dx[1]*fabs(bcon[2])/sqrt(rho[zonalindex]);
	#endif

	//to calculate magn. field angle
	ldouble brbphi,bsq,bfake[4];
	#ifdef BHDISK_PROBLEMTYPE
	calc_angle_brbphibsq(ix,iy,iz,&brbphi,&bsq,bfake,bfake);
	Bangle[zonalindex]=-brbphi/bsq;
	#else
	Bangle[zonalindex]=-1.;
	#endif

	if(ix==0 || (NY>1 && iy==0) || (NZ>1 && iz==0)) //divB left-biased
	  divB[zonalindex]=0;
	else
	  divB[zonalindex]=calc_divB(ix,iy,iz);
#endif //MAGNFIELD

#ifdef RADIATION
	taucoupling[zonalindex]=estimate_gas_radiation_coupling(pp,&geomout);
#endif //RADIATION

        // electrons & ions
	ueloc=uiloc=0.;
#ifdef EVOLVEELECTRONS
	// temperatures
	calc_PEQ_Teifrompp(pp,&tempeloc,&tempiloc,ix,iy,iz);

	// electrons
	ldouble ne=calc_thermal_ne(pp); 
	ldouble pe=K_BOLTZ*ne*tempeloc;
	ldouble gammae=GAMMAE;
        #ifdef CONSISTENTGAMMA
        #ifndef FIXEDGAMMASPECIES
	gammae=calc_gammaintfromtemp(tempeloc,ELECTRONS);
        #endif
        #endif
	ueloc=pe/(gammae-1.);

	//ions
	ldouble ni=pp[RHO]/MU_I/M_PROTON; //number density of photons and electrons
	ldouble pi=K_BOLTZ*ni*tempiloc;
	ldouble gammai=GAMMAI;
        #ifdef CONSISTENTGAMMA
        #ifndef FIXEDGAMMASPECIES
	gammai=calc_gammaintfromtemp(tempiloc,IONS);
        #endif
        #endif
	uiloc=pi/(gammai-1.);
#endif //EVOLVEELECTRONS

#ifdef PRINTVISCHEATINGTOSILO
	deltae[zonalindex]=calc_ViscousElectronHeatingFraction(&get_u(p,0,ix,iy,iz),&geomout);
	vischeat[zonalindex]=get_u_scalar(vischeating,ix,iy,iz);
	vischeatnege[zonalindex]=get_u_scalar(vischeatingnegebalance,ix,iy,iz);; 
	vischeatnegi[zonalindex]=get_u_scalar(vischeatingnegibalance,ix,iy,iz);; 
	dtauarr[zonalindex]=-1.;
#endif //PRINTVISCHEATINGTOSILO
      }
      else //using averaged data
      {

	 rho[zonalindex] = get_uavg(pavg,RHO,ix,iy,iz);
	 entr[zonalindex] = get_uavg(pavg,ENTR,ix,iy,iz);
	 uint[zonalindex] = get_uavg(pavg,UU,ix,iy,iz);

	 pregas = get_uavg(pavg,AVGPGAS,ix,iy,iz);
	 gamma = 1. + pregas/uint[zonalindex];

	 // velocity
	 vel[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	 vel[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	 vel[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	 vel[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	 for(i=0;i<4;i++) vcon[i]=vel[i];
	 indices_21(vcon,vcov,geomout.gg);

	 rhouconr=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz);
	 rhoucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz);

	 // angular velocity
	 Omega[zonalindex]=vel[3]/vel[0];

         // ANDREW these are 4-velocites, not VELR, so I don't think this is right...
	 //updates pp[VI] to have rho-weighted velocities there
	 //pp[VX]=vel[1]; 
	 //pp[VY]=vel[2];
	 //pp[VZ]=vel[3];

	 #ifdef MAGNFIELD
   	 B1comp[zonalindex] = get_uavg(pavg,B1,ix,iy,iz);
         B2comp[zonalindex] = get_uavg(pavg,B2,ix,iy,iz);
         B3comp[zonalindex] = get_uavg(pavg,B3,ix,iy,iz);

	 bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
	 bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
	 bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
	 bsq[zonalindex]=get_uavg(pavg,AVGBSQ,ix,iy,iz);

	 premag=0.5*bsq[zonalindex];
	 beta[zonalindex]=pregas/0.5/bsq[zonalindex];
	 betainv[zonalindex]=1./beta[zonalindex];
	 sigma[zonalindex]=bsq[zonalindex]/(rho[zonalindex] + uint[zonalindex] + pregas);
	 sigmaw[zonalindex]=bsq[zonalindex]/(rho[zonalindex] + gamma*uint[zonalindex]);

     
         #endif		  
	 
	 // stress energy tensor
         /*
	 for(i=0;i<4;i++)
	   for(j=0;j<4;j++)
	     Tij[i][j]=get_uavg(pavg,AVGRHOUCONUCOV(i,j),ix,iy,iz)
	       + GAMMA*get_uavg(pavg,AVGUUUCONUCOV(i,j),ix,iy,iz)
	       + get_uavg(pavg,AVGBSQUCONUCOV(i,j),ix,iy,iz)
	       + delta(i,j)*(GAMMA*get_uavg(pavg,UU,ix,iy,iz) + 1./2.*get_uavg(pavg,AVGBSQ,ix,iy,iz))
	       - get_uavg(pavg,AVGBCONBCOV(i,j),ix,iy,iz); 
	 */
	 
	 // averaged stress energy tensor
	 for(i=0;i<4;i++)
	   for(j=0;j<4;j++)
	     Tij[i][j]=get_uavg(pavg,AVGTIJ(i,j),ix,iy,iz);


	 Tit[1]=Tij[1][0];
	 Tit[2]=Tij[2][0];
	 Tit[3]=Tij[3][0];

	 // Bernoulli flux / number
	 muBe[zonalindex]=-(Tij[1][0] + rhouconr)/rhouconr;
	 Be[zonalindex]=-(Tij[0][0] + rhoucont)/rhoucont;

#ifdef MAGNFIELD

	 Qtheta[zonalindex]=2.*M_PI/Omega[zonalindex]/dx[1]*fabs(bcon[2])/sqrt(rho[zonalindex]);
	 Qphi[zonalindex]=2.*M_PI/Omega[zonalindex]/dx[2]*fabs(bcon[3])/sqrt(rho[zonalindex]);
	 ldouble brbphi,bsq,bfake[4];

	 #ifdef BHDISK_PROBLEMTYPE
	 calc_angle_brbphibsq(ix,iy,iz,&brbphi,&bsq,bfake,bfake); //to calculate magn. field angle
	 Bangle[zonalindex]=-brbphi/bsq;
	 #else
	 Bangle[zonalindex]=-1.;
	 #endif

	 if(ix==0 || (NY>1 && iy==0) || (NZ>1 && iz==0)) //divB left-biased
	   divB[zonalindex]=0;
	 else
	   divB[zonalindex]=calc_divB(ix,iy,iz);
#endif // MAGNFIELD

#ifdef RADIATION
	 taucoupling[zonalindex]=estimate_gas_radiation_coupling(pp,&geomout);
#endif // RADIATION

#ifdef EVOLVEELECTRONS
	 ldouble pe,pi;
	 pe=get_uavg(pavg,AVGPE,ix,iy,iz);
	 pi=get_uavg(pavg,AVGPI,ix,iy,iz);
	 #ifndef CONSISTENTGAMMA
	 pi=pregas-pe;
	 #endif

	 //electrons
	 ne=calc_thermal_ne(pp); 
	 tempeloc=pe/K_BOLTZ/ne;
	 ldouble gammae=GAMMAE;
	 #ifdef CONSISTENTGAMMA
	 #ifndef FIXEDGAMMASPECIES
	 gammae=calc_gammaintfromtemp(tempeloc,ELECTRONS);
	 #endif
	 #endif
	 ueloc=pe/(gammae-1.);

	 //ions
	 ldouble ni=rho[zonalindex]/MU_I/M_PROTON; 
	 tempiloc=pi/K_BOLTZ/ni;
	 ldouble gammai=GAMMAI;
	 #ifdef CONSISTENTGAMMA
	 #ifndef FIXEDGAMMASPECIES
	 gammai=calc_gammaintfromtemp(tempiloc,IONS);
	 #endif
	 #endif
	 uiloc=pi/(gammai-1.);
#endif // EVOLVEELECTRONS


#ifdef PRINTVISCHEATINGTOSILO
	 vischeat[zonalindex]=get_uavg(pavg,AVGVISCHEATING,ix,iy,iz);
	 vischeatnege[zonalindex]=get_uavg(pavg,AVGVISCHEATINGNEGE,ix,iy,iz);
	 vischeatnegi[zonalindex]=get_uavg(pavg,AVGVISCHEATINGNEGI,ix,iy,iz);
	 deltae[zonalindex]=-1.;
	 dtauarr[zonalindex]=-1.;
#endif //PRINTVISCHEATINGTOSILO	 


	 
      } //doingavg

      
// Nonthermal electron quantities
#ifdef EVOLVEELECTRONS
#ifdef RELELECTRONS
      nethloc=ne;
      nrelelloc = calc_relel_ne(pp);
      urelelloc = calc_relel_uint(pp);
      G0relelloc = -1.*calc_relel_G0_fluidframe(pp,&geomout, 0.0, 0); 
      G0ic_relel_loc = -1*calc_relel_G0_fluidframe_direct(pp, &geomout, 3);
      G0syn_relel_loc = -1*calc_relel_G0_fluidframe_direct(pp, &geomout, 1);

      //calculate synchrotron break frequency
      //ANDREW speed this up by saving an array
      gbrkloc=RELEL_INJ_MIN;
      //absolute maximum of g^4*n for g > RELGAMMAMIN
      ldouble nbrk=pp[NEREL(0)]*gammapbrk[0];
      ldouble nbrk2;
      for(ie=1;ie<NRELBIN;ie++)
      {
	  if (relel_gammas[ie] < RELEL_INJ_MIN)
	  {
	     gbrkloc=RELEL_INJ_MIN;
	     nbrk =  pp[NEREL(ie)]*gammapbrk[ie];
	  }

	  else 
	  {
	     nbrk2 =  pp[NEREL(ie)]*gammapbrk[ie];
	     if(nbrk2 > nbrk)
	     {
	       nbrk=nbrk2;
	       gbrkloc=relel_gammas[ie];
	     }
	  }
      }
#endif //RELELECTRONS
#endif //EVOLVEELECTRONS

      // put quantities into arrays
      expansion[zonalindex]=exploc;
      lorentz[zonalindex]=fabs(vel[0])/sqrt(fabs(geomout.GG[0][0]));
      u0[zonalindex]=fabs(vel[0]);
      u1[zonalindex]=vel[1];
      u2[zonalindex]=vel[2];
      u3[zonalindex]=vel[3];

#ifdef MAGNFIELD
      //OmegaF[zonalindex] = Omega[zonalindex] - (u1[zonalindex]/u0[zonalindex])*(B3comp[zonalindex]/B1comp[zonalindex]);

      OmegaF[zonalindex]=Omega[zonalindex];
      if(fabs(u2[zonalindex]*B3comp[zonalindex])>1.e-8)
	OmegaF[zonalindex] -= (u2[zonalindex]*B3comp[zonalindex])/(u0[zonalindex]*B2comp[zonalindex]);

      ldouble velff[4];
      ldouble vpar;

      vpar = get_driftvel_velr(pp,velff,&geomout);
      conv_vels(velff,velff,VELR,VEL4,geomout.gg,geomout.GG);
      
      lorentz_perp[zonalindex]=fabs(velff[0])/sqrt(fabs(geomout.GG[0][0]));
      lorentz_par[zonalindex]=1./sqrt(1+pow(lorentz[zonalindex],-2)-pow(lorentz_perp[zonalindex],-2));

#endif

#ifdef MAGNFIELD
#ifdef FORCEFREE // must have VELPRIM=VELR

      //ffinv[zonalindex]=get_cflag(FFINVFLAG,ix,iy,iz);
      ffinv[zonalindex] = get_u_scalar(ffinvarr, geom.ix,geom.iy,geom.iz);
      
      int derdir2[3] = {2,2,2};
      ldouble jcon[4],jcov[4];
      calc_current(&geom,jcon,derdir2);
      indices_21(jcon,jcov,geom.gg); // not in lab frame

      //ldouble FF[4][4];
      //calc_faraday(FF,geom.ix,geom.iy,geom.iz,0);
      jdens[zonalindex] = sqrt(fabs(dot(jcon,jcov)));

#else
      jdens[zonalindex]=0.;      
      ffinv[zonalindex]=0.;
#endif
#endif
      
      #ifdef EVOLVEELECTRONS
      tempi[zonalindex]=tempiloc; //ion temperature
      tempe[zonalindex]=tempeloc;  //electron temperature
      ui[zonalindex]=uiloc; //ion energy density
      ue[zonalindex]=ueloc;  //electron energy density

      #ifdef RELELECTRONS
      gammabrk[zonalindex]=gbrkloc;
      urelel[zonalindex]=urelelloc;
      nrelel[zonalindex]=nrelelloc;
      neth[zonalindex]=nethloc;
      uratio_tot[zonalindex]=urelelloc/uint[zonalindex];
      uratio_th[zonalindex]=urelelloc/ueloc;
      G0relel[zonalindex]=G0relelloc;
      G0icrelel[zonalindex]=G0ic_relel_loc;
      G0synrelel[zonalindex]=G0syn_relel_loc;
      #endif
      #endif 

      ldouble temploc=calc_PEQ_Tfromurho(uint[zonalindex],rho[zonalindex],ix,iy,iz);
      temp[zonalindex]=temploc;
      entrlnT[zonalindex]=kB_over_mugas_mp*entr[zonalindex]/pp[RHO];

#ifdef CGSOUTPUT
      rho[zonalindex]=rhoGU2CGS(rho[zonalindex]);
      uint[zonalindex]=endenGU2CGS(uint[zonalindex]);

      #ifdef MAGNFIELD 
      bsq[zonalindex]=4.*M_PI*endenGU2CGS(bsq[zonalindex]);
      #endif

      #ifdef EVOLVEELECTRONS
      ui[zonalindex]=endenGU2CGS(ui[zonalindex]);
      ue[zonalindex]=endenGU2CGS(ue[zonalindex]);

      #ifdef RELELECTRONS
      urelel[zonalindex]=endenGU2CGS(urelel[zonalindex]);
      nrelel[zonalindex]=numdensGU2CGS(nrelel[zonalindex]);
      neth[zonalindex]=numdensGU2CGS(neth[zonalindex]);
      G0relel[zonalindex]=endenGU2CGS(G0relel[zonalindex])*timeCGS2GU(1.);
      G0icrelel[zonalindex]=endenGU2CGS(G0icrelel[zonalindex])*timeCGS2GU(1.);
      G0synrelel[zonalindex]=endenGU2CGS(G0synrelel[zonalindex])*timeCGS2GU(1.);
      #endif
      #endif

      #ifdef  PRINTVISCHEATINGTOSILO
      vischeat[zonalindex]    =endenGU2CGS(vischeat[zonalindex])*timeCGS2GU(1.);
      vischeatnege[zonalindex]=endenGU2CGS(vischeatnege[zonalindex])*timeCGS2GU(1.);
      vischeatnegi[zonalindex]=endenGU2CGS(vischeatnegi[zonalindex])*timeCGS2GU(1.);
      #endif	      
#endif //CGSOUTPUT

#ifdef PRINTCOULOMBTOSILO
      coulomb[zonalindex]=calc_CoulombCoupling(pp,&geomout);
#endif

#ifdef PRINTGAMMATOSILO
      gammag[zonalindex]=gamma;
#endif     	      

      // vector velocities and energy fluxes
      vx[zonalindex]=vel[1];
      vy[zonalindex]=vel[2];
      vz[zonalindex]=vel[3];

      Edotx[zonalindex]=Tit[1];
      Edoty[zonalindex]=Tit[2];
      Edotz[zonalindex]=Tit[3];

      // transform velocity to cartesian
      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS ||
	  MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS || MYCOORDS==MKS3COORDS || MYCOORDS==JETCOORDS ||
	  MYCOORDS==MSPH1COORDS || MYCOORDS==MKER1COORDS)
      {
	  vel[2]*=r;
	  vel[3]*=r*sin(th);

	  vx[zonalindex] = sin(th)*cos(ph)*vel[1] 
	    + cos(th)*cos(ph)*vel[2]
	    - sin(ph)*vel[3];

	  vy[zonalindex] = sin(th)*sin(ph)*vel[1] 
	    + cos(th)*sin(ph)*vel[2]
	    + cos(ph)*vel[3];

	  vz[zonalindex] = cos(th)*vel[1] 
	    - sin(th)*vel[2];

	  Tit[2]*=r;
	  Tit[3]*=r*sin(th);

	  Edotx[zonalindex] = sin(th)*cos(ph)*Tit[1] 
	    + cos(th)*cos(ph)*Tit[2]
	    - sin(ph)*Tit[3];

	  Edoty[zonalindex] = sin(th)*sin(ph)*Tit[1] 
	    + cos(th)*sin(ph)*Tit[2]
	    + cos(ph)*Tit[3];

	  Edotz[zonalindex] = cos(th)*Tit[1] 
	    - sin(th)*Tit[2];
      }

      // vector magnetic field quantities
#ifdef MAGNFIELD
      //magnetic field	      
      Bx[zonalindex]=bcon[1];
      By[zonalindex]=bcon[2];
      Bz[zonalindex]=bcon[3];

      #ifdef BATTERY // battery field
      ldouble dBdtbat[4];
      estimate_Bgrowth_battery(ix,iy,iz,dBdtbat);
      dBxdtbat[zonalindex]=dBdtbat[1];
      dBydtbat[zonalindex]=dBdtbat[2];
      dBzdtbat[zonalindex]=dBdtbat[3];
      #endif

      #ifdef MIMICDYNAMO // dynamo field
      ldouble ppdyn[NV]; int idyn;
      PLOOP(idyn) ppdyn[idyn]=get_u(p,idyn,ix,iy,iz);	      
      ppdyn[B1]=get_u(pvecpot,1,ix,iy,iz);
      ppdyn[B2]=get_u(pvecpot,2,ix,iy,iz);
      ppdyn[B3]=get_u(pvecpot,3,ix,iy,iz);

      #ifdef PRECOMPUTE_MY2OUT
      trans_pmhd_coco_my2out(ppdyn, ppdyn, &geom, &geomout);
      #else      
      trans_pmhd_coco(ppdyn, ppdyn, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);
      #endif

      ldouble bcondyn[4],bcovdyn[4];
      calc_bcon_prim(ppdyn,bcondyn,&geomout);
      //indices_21(bcondyn,bcovdyn,geomout.gg); 

      Bxdyn[zonalindex]=bcondyn[1];
      Bydyn[zonalindex]=bcondyn[2];
      Bzdyn[zonalindex]=bcondyn[3];
      #endif

      // calculate vector potential
      //ANDREW changed phi defn to match harmpi
      int iphimin,iphimax;
      iphimin=0;
      iphimax=ny-1;

      #if defined(CORRECT_POLARAXIS) || defined(CORRECT_POLARAXIS_3D)
      //iphimin=NCCORRECTPOLAR; 
      #ifndef HALFTHETA
      //iphimax=ny-NCCORRECTPOLAR-1;
      #endif
      #endif
      if(iy==iphimin)
      {
	 //phi[zonalindex]=geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1)*2.*M_PI;
	 phi_accum[zonalindex]=geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1);//*2.*M_PI;
	 phi[zonalindex] = phi_accum[zonalindex] - 0.5*geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1);
      }
      else if(iy>iphimin && iy<=iphimax)
      {
	 imz=iz;imy=iy;imx=ix;
	 #ifdef PRINTXGC_RIGHT
	 imx=ix-NG;
	 #endif
	 #ifdef PRINTXGC_LEFT
	 imx=ix-NG;
	 #endif
	 #ifdef PRINTYGC_RIGHT
	 imy=iy-NG;
	 #endif
	 int idx=imz*(ny*nx) + (imy-1)*nx + imx;
	 //phi[zonalindex]=phi[idx]+geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1)*2.*M_PI;
	 phi_accum[zonalindex]=phi_accum[idx] + geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1); //*2.*M_PI;
	 phi[zonalindex] = phi_accum[zonalindex] - 0.5*geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1);
      }

#ifdef MIMICDYNAMO
      Avecdyn[zonalindex]=get_u(ptemp1,B3,ix,iy,iz);
      if(iy==iphimin)
      {
	  phidyn[zonalindex]=geom.gdet*get_u(pvecpot,1,ix,iy,iz)*get_size_x(iy,1)*2.*M_PI;
      }
      else if(iy>iphimin && iy<=iphimax)
      {
	  imz=iz;imy=iy;imx=ix;
	  #ifdef PRINTXGC_RIGHT
	  imx=ix-NG;
	  #endif
	  #ifdef PRINTXGC_LEFT
	  imx=ix+NG;
	  #endif
	  #ifdef PRINTYGC_RIGHT
	  imy=iy-NG;
	  #endif
	  int idx=imz*(ny*nx) + (imy-1)*nx + imx;
	  phidyn[zonalindex]=phidyn[idx]+geom.gdet*get_u(pvecpot,1,ix,iy,iz)*get_size_x(iy,1)*2.*M_PI;
      }
#endif

#ifdef BATTERY
      if(iy==iphimin)
      {
	  phibat[zonalindex]=geom.gdet*dBdtbat[1]*get_size_x(iy,1)*2.*M_PI;
      }
      else if(iy>iphimin && iy<=iphimax)
      {
	  imz=iz;imy=iy;imx=ix;
	  #ifdef PRINTXGC_RIGHT
	  imx=ix-NG;
	  #endif
	  #ifdef PRINTXGC_LEFT
	  imx=ix+NG;
	  #endif
	  #ifdef PRINTYGC_RIGHT
	  imy=iy-NG;
	  #endif
	  int idx=imz*(ny*nx) + (imy-1)*nx + imx;
	  phibat[zonalindex]=phibat[idx]+geom.gdet*dBdtbat[1]*get_size_x(iy,1)*2.*M_PI;
      }
#endif

      //transform Bfield to cartesian
      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS ||
	  MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS || MYCOORDS==MKS3COORDS || MYCOORDS==JETCOORDS ||
	  MYCOORDS==MSPH1COORDS || MYCOORDS==MKER1COORDS)
      {
	  bcon[2]*=r;
	  bcon[3]*=r*sin(th);

	  Bx[zonalindex] = sin(th)*cos(ph)*bcon[1] 
	    + cos(th)*cos(ph)*bcon[2]
	    - sin(ph)*bcon[3];

	  By[zonalindex] = sin(th)*sin(ph)*bcon[1] 
	    + cos(th)*sin(ph)*bcon[2]
	    + cos(ph)*bcon[3];

	  Bz[zonalindex] = cos(th)*bcon[1] 
	    - sin(th)*bcon[2];

#ifdef BATTERY
	  dBdtbat[2]*=r;
	  dBdtbat[3]*=r*sin(th);

	  dBxdtbat[zonalindex] = sin(th)*cos(ph)*dBdtbat[1] 
	    + cos(th)*cos(ph)*dBdtbat[2]
	    - sin(ph)*dBdtbat[3];

	  dBydtbat[zonalindex] = sin(th)*sin(ph)*dBdtbat[1] 
	    + cos(th)*sin(ph)*dBdtbat[2]
	    + cos(ph)*dBdtbat[3];

	  dBzdtbat[zonalindex] = cos(th)*dBdtbat[1] 
	    - sin(th)*dBdtbat[2];

#endif

#ifdef MIMICDYNAMO
	  bcondyn[2]*=r;
	  bcondyn[3]*=r*sin(th);

	  Bxdyn[zonalindex] = sin(th)*cos(ph)*bcondyn[1] 
	    + cos(th)*cos(ph)*bcondyn[2]
	    - sin(ph)*bcondyn[3];

	  Bydyn[zonalindex] = sin(th)*sin(ph)*bcondyn[1] 
	    + cos(th)*sin(ph)*bcondyn[2]
	    + cos(ph)*bcondyn[3];

	  Bzdyn[zonalindex] = cos(th)*bcondyn[1] 
	    - sin(th)*bcondyn[2];
#endif
        }
#endif //MAGNFIELD

#ifdef RADIATION

	ldouble Rtt,ehat,ugas[4],urad[4],rvel[4],Rij[4][4],Rij22[4][4],Giff[4],tradloc,tradlteloc;
	struct opacities opac;  
	ldouble tauabsloc = vcon[0]*calc_kappa(pp,&geomout,&opac);
	ldouble tauscaloc = vcon[0]*calc_kappaes(pp,&geomout);
	ldouble taueffloc = sqrt(tauabsloc*(tauabsloc+tauscaloc));

	if(doingavg==0) //from snapshot
	{
	    calc_ff_Rtt(pp,&Rtt,ugas,&geomout);
	    ehat=-Rtt;  	      							  
	    calc_Rij(pp,&geomout,Rij); //calculates R^munu in OUTCOORDS

	    //four fource
	    calc_Gi(pp,&geomout,Giff,0.0,0,0); // fluid frame 

	    #ifdef RADFLUXFFINOUTPUT
	    boost22_lab2ff(Rij,Rij,pp,geomout.gg,geomout.GG);
	    #endif
	    indices_2221(Rij,Rij,geomout.gg);

	    //calculating radiation temperature
	    tradlteloc=calc_LTE_TfromE(ehat);
	    tradloc=calc_LTE_TfromE(ehat);

	    #ifdef EVOLVEPHOTONNUMBER
	    tradloc=calc_ncompt_Thatrad(pp,&geomout,ehat);
	    #endif
	}
	else //using avg
	{
	   ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
	   tradloc=calc_LTE_TfromE(ehat);
	   tradlteloc=calc_LTE_TfromE(ehat);
	   for(j=0;j<4;j++)
	     Giff[j]=get_uavg(pavg,AVGGHAT(j),ix,iy,iz);

	   //ANDREW if Gi0 accidentally saved in lab frame:
	   #ifdef SIMOUTPUT_GILAB2FF
	   boost2_lab2ff(Giff,Giff,pp,geomout.gg,geomout.GG); //ANDREW avg Gff already in OUTCOORDS
	   #endif
	   #ifdef EVOLVEPHOTONNUMBER		  
	   tradloc=calc_ncompt_Thatrad_fromEN(ehat,get_uavg(pavg,AVGNFHAT,ix,iy,iz));
	   #endif
	   for(i=0;i<4;i++)
	     for(j=0;j<4;j++)
	       Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
	}
	indices_2122(Rij,Rij22,geomout.GG);

	// Add radiation to Bernoulli number
	muBe[zonalindex]+=-Rij[1][0]/rhouconr;
	Be[zonalindex]+=-Rij[0][0]/rhoucont;

	//correcting rad-velocities basing on <R^t_mu>
	int radcorr;

	p2u(pp,uu,&geomout);

	uu[EE0]=gdetu*Rij[0][0];
	uu[FX0]=gdetu*Rij[0][1];
	uu[FY0]=gdetu*Rij[0][2];
	uu[FZ0]=gdetu*Rij[0][3];

	//print_conserved(uu);
	u2p_rad(uu,pp,&geomout,&radcorr);

	#ifdef COMPTONIZATION
	//Calculate thermal comptonization
	ldouble Gic_tmp[4];
	ldouble ucon_ff[4];
	ucon_ff[1]=ucon_ff[2]=ucon_ff[3]=0.;
	ucon_ff[0]=1.;
	ldouble kappaes=calc_kappaes(pp, &geomout);
	calc_Compt_Gi(pp, &geomout, Gic_tmp, ehat, tempeloc, kappaes, ucon_ff);
	G0ic_th_loc = fabs(Gic_tmp[0]); 

	G0icth[zonalindex]=G0ic_th_loc;
	#endif

	Ehat[zonalindex]=ehat;
	Gtff[zonalindex]=-Giff[0];
	Erad[zonalindex]=Rij22[0][0];

	#ifdef EVOLVEPHOTONNUMBER
	Nph[zonalindex]=pp[NF0];
	#endif

	#ifdef CGSOUTPUT
	Ehat[zonalindex]=endenGU2CGS(Ehat[zonalindex]);
	Gtff[zonalindex]=endenGU2CGS(Gtff[zonalindex])*timeCGS2GU(1.);
	Erad[zonalindex]=endenGU2CGS(Erad[zonalindex]);
	#ifdef EVOLVEPHOTONNUMBER
	Nph[zonalindex]=numdensGU2CGS(Nph[zonalindex]);
	#endif
	#ifdef COMPTONIZATION
	G0icth[zonalindex]=endenGU2CGS(G0icth[zonalindex]) * timeCGS2GU(1.);
	#endif
	#endif

	trad[zonalindex]=tradloc;
	tradlte[zonalindex]=tradlteloc;
	entrrad[zonalindex]=(4./3.)*A_RAD*tradloc*tradloc*tradloc/pp[RHO];

	Fx[zonalindex]=Rij22[1][0];
	Fy[zonalindex]=Rij22[2][0];
	Fz[zonalindex]=Rij22[3][0];

	urad[1]=pp[FX0];
	urad[2]=pp[FY0];
	urad[3]=pp[FZ0];
	conv_vels(urad,urad,VELPRIM,VEL4,geomout.gg,geomout.GG);

	uradx[zonalindex]=urad[1];
	urady[zonalindex]=urad[2];
	uradz[zonalindex]=urad[3];

	//Angular integration (tau_theta)
	if(iy==0)
	{
	    tausca[zonalindex]=tauscaloc*dxph[1];
	    tauabs[zonalindex]=tauabsloc*dxph[1];
	    taueff[zonalindex]=taueffloc*dxph[1];
	}
	else
	{
	    imz=iz;
	    imy=iy;
	    imx=ix;
	    #ifdef PRINTXGC_RIGHT
	    imx=ix-NG;
	    #endif
	    #ifdef PRINTXGC_LEFT
	    imx=ix+NG;
	    #endif
	    #ifdef PRINTYGC_RIGHT
	    imy=iy-NG;
	    #endif
	    int idx=imz*(ny*nx) + (imy-1)*nx + imx;
	    if(iy<=NY/2) //proper integration only in the upper half
	    {
		tausca[zonalindex]=tausca[idx]+tauscaloc*dxph[1];
		tauabs[zonalindex]=tauabs[idx]+tauabsloc*dxph[1];
		taueff[zonalindex]=taueff[idx]+taueffloc*dxph[1];
	    }
	    else
	    {
		idx=imz*(ny*nx) + (NY/2-1)*nx + imx;
		tausca[zonalindex]=tausca[idx];
		tauabs[zonalindex]=tauabs[idx];
		taueff[zonalindex]=taueff[idx];
	    }
	}

	//transform rad flux to cartesian
	if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS ||
	    MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS || MYCOORDS==MKS3COORDS || MYCOORDS==JETCOORDS ||
	    MYCOORDS==MSPH1COORDS || MYCOORDS==MKER1COORDS)
	{
	    Rij22[2][0]*=sqrt(geomout.gg[2][2]);
	    Rij22[3][0]*=sqrt(geomout.gg[3][3]);

	    Fx[zonalindex] = sin(th)*cos(ph)*Rij22[1][0] 
	      + cos(th)*cos(ph)*Rij22[2][0]
	      - sin(ph)*Rij22[3][0];

	    Fy[zonalindex] = sin(th)*sin(ph)*Rij22[1][0] 
	      + cos(th)*sin(ph)*Rij22[2][0]
	      + cos(ph)*Rij22[3][0];

	    Fz[zonalindex] = cos(th)*Rij22[1][0] 
	      - sin(th)*Rij22[2][0];

	    urad[2]*=r;
	    urad[3]*=r*sin(th);

	    uradx[zonalindex] = sin(th)*cos(ph)*urad[1]
	      + cos(th)*cos(ph)*urad[2]
	      - sin(ph)*urad[3];

	    urady[zonalindex] = sin(th)*sin(ph)*urad[1]
	      + cos(th)*sin(ph)*urad[2]
	      + cos(ph)*urad[3];

	    uradz[zonalindex] = cos(th)*urad[1]
	      - sin(th)*urad[2];
	}

#endif //RADIATION	      
      }
    }
  }
  

  //Loop for radially integrated optical depth

#ifdef RADIATION  

#pragma omp parallel for private(ix,iy,iz,iv,imx,imy,imz,i,j,pp,uu,xxvec,xxveccar,xxvecsph,xx1,xx2) schedule (static)

#ifdef PRINTZGC_RIGHT
  for(iz=NG;iz<nz+NG;iz++)
#else
  for(iz=0;iz<nz;iz++)
#endif
  {
#if defined(PRINTXGC_RIGHT)
    for(ix = nx+NG-1; ix > NG-1; ix--)
#elif defined(PRINTXGC_LEFT)
    for(ix = nx-NG-1;ix > -NG-1; ix--)
#else
    for(ix=(nx-1);ix>-1;ix--)
#endif
    {
#ifdef PRINTYGC_RIGHT
      for(iy=NG;iy<ny+NG;iy++)
#else
      for(iy=0;iy<ny;iy++)
#endif
      {  

	ldouble gamma=GAMMA;
	#ifdef CONSISTENTGAMMA
	gamma=pick_gammagas(ix,iy,iz);
	#endif

	int iix,iiy,iiz;
	iix=ix;
	iiy=iy;
	iiz=iz;

	struct geometry geom;
	fill_geometry(iix,iiy,iiz,&geom);
	struct geometry geomout;
	fill_geometry_arb(iix,iiy,iiz,&geomout,OUTCOORDS);

	ldouble dxph[3],dx[3];
	//cell dimensions
	get_cellsize_out(ix, iy, iz, dx);
	dxph[0]=dx[0]*sqrt(geomout.gg[1][1]);
	dxph[1]=dx[1]*sqrt(geomout.gg[2][2]);
	dxph[2]=dx[2]*sqrt(geomout.gg[3][3]);

	get_xx(iix,iiy,iiz,xxvec);
	coco_N(xxvec,xxvecsph,MYCOORDS,SPHCOORDS);
	coco_N(xxvec,xxveccar,MYCOORDS,MINKCOORDS);

	ldouble r=xxvecsph[1];
	ldouble th=xxvecsph[2];
	ldouble ph=xxvecsph[3];

	// metric determinant
	ldouble gdet,gdetu;
	gdet=geomout.gdet;
	gdetu=gdet;
	#if (GDETIN==0) //gdet out of derivatives
	gdetu=1.;
	#endif

	//gdet and coordinates of cells +- 1 in radius
	ldouble gdet1,gdet2;
	#ifdef PRECOMPUTE_MY2OUT
	get_xxout(ix-1, iiy, iiz, xx1);
	get_xxout(ix+1, iiy, iiz, xx2);	      
	#else
	get_xx(iix-1,iiy,iiz,xx1);
	get_xx(iix+1,iiy,iiz,xx2);
	coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
	coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
	#endif

	gdet1=calc_gdet_arb(xx1,OUTCOORDS);
	gdet2=calc_gdet_arb(xx2,OUTCOORDS);

	// calculate zonal index
	imz=iz;
	imy=iy;
	imx=ix;
	#ifdef PRINTZGC_RIGHT
	imz=iz-NG;
	#endif
	#ifdef PRINTXGC_RIGHT
	imx=ix-NG;
	#endif
	#ifdef PRINTXGC_LEFT
	imx=ix+NG;
	#endif
	#ifdef PRINTYGC_RIGHT
	imy=iy-NG;
	#endif

	int zonalindex=imz*(ny*nx) + imy*nx + imx;
	for(iv=0;iv<NV;iv++)
	  {
	    if(doingavg)
	      pp[iv]=get_uavg(pavg,iv,ix,iy,iz); //this should not be used later but it is //ANDREW ??
	    else
	      pp[iv]=get_u(p,iv,iix,iiy,iiz);
	  }

	//primitives to OUTCOORDS
	#ifdef PRECOMPUTE_MY2OUT
	trans_pall_coco_my2out(pp,pp,&geom,&geomout);
	#else      
	trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);
	#endif

	//velocities etc
	ldouble vel[4],vcov[4],vcon[4],velprim[4];

	if(doingavg==0) //using snapshot data
	{		  
	    vel[1]=pp[VX];
	    vel[2]=pp[VY];
	    vel[3]=pp[VZ];

	    conv_vels(vel,vel,VELPRIM,VEL4,geomout.gg,geomout.GG);
	    for(i=0;i<4;i++) vcon[i]=vel[i];
	    indices_21(vel,vcov,geomout.gg);
	}
	else //using averaged data
	{
	     vel[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	     vel[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	     vel[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	     vel[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	     for(i=0;i<4;i++) vcon[i]=vel[i];
	     indices_21(vel,vcov,geomout.gg);
	}

	#ifdef RADIATION
	struct opacities opac;  
	//ldouble tauabslocr = vcon[0]*(1.-fabs(vcon[1]))*calc_kappa(pp,&geomout,&opac);
	//ldouble tauscalocr = vcon[0]*(1.-fabs(vcon[1]))*calc_kappaes(pp,&geomout);
	ldouble tauabslocr = (vcon[0]-vcon[1])*calc_kappa(pp,&geomout,&opac);
	ldouble tauscalocr = (vcon[0]-vcon[1])*calc_kappaes(pp,&geomout);

	ldouble tauefflocr = sqrt(tauabslocr*(tauabslocr+tauscalocr));

	ldouble tauabslocr2 = -(vcov[0]+vcov[1])*calc_kappa(pp,&geomout,&opac);
	ldouble tauscalocr2 = -(vcov[0]+vcov[1])*calc_kappaes(pp,&geomout);
	ldouble tauefflocr2 = sqrt(tauabslocr2*(tauabslocr2+tauscalocr2));
	#endif

	//Radial integration (tau_theta)
	if(ix==NX-1)
	{
	    #ifdef RADIATION
	    tauscar[zonalindex]=tauscalocr*dxph[0];
	    tauabsr[zonalindex]=tauabslocr*dxph[0];
	    taueffr[zonalindex]=tauefflocr*dxph[0];

	    tauscar2[zonalindex]=tauscalocr2*dxph[0];
	    tauabsr2[zonalindex]=tauabslocr2*dxph[0];
	    taueffr2[zonalindex]=tauefflocr2*dxph[0];
	    #endif
	    //Column density (N_H) in CGS
	    NHr[zonalindex] = rho[zonalindex]*dxph[0]/(MU_GAS*M_PROTON)/lenGU2CGS(1)/lenGU2CGS(1);
	}
	else
	{
	    imz=iz;
	    imy=iy;
	    imx=ix;
	    #ifdef PRINTXGC_RIGHT
	    imx=ix-NG;
	    #endif
	    #ifdef PRINTXGC_LEFT
	    imx=ix+NG;
	    #endif
	    #ifdef PRINTYGC_RIGHT
	    imy=iy-NG;
	    #endif
	    int idx=imz*(ny*nx) + imy*nx + (imx+1);
	    #ifdef RADIATION
	    tauscar[zonalindex]=tauscar[idx]+tauscalocr*dxph[0];
	    tauabsr[zonalindex]=tauabsr[idx]+tauabslocr*dxph[0];
	    taueffr[zonalindex]=taueffr[idx]+tauefflocr*dxph[0];

	    tauscar2[zonalindex]=tauscar2[idx]+tauscalocr2*dxph[0];
	    tauabsr2[zonalindex]=tauabsr2[idx]+tauabslocr2*dxph[0];
	    taueffr2[zonalindex]=taueffr2[idx]+tauefflocr2*dxph[0];
	    #endif
	    //Column density (N_H) in CGS
	    NHr[zonalindex] = NHr[idx] + rho[zonalindex]*dxph[0]/(MU_GAS*M_PROTON)/lenGU2CGS(1)/lenGU2CGS(1);		  
	}
      }
    }
  }
#endif //RADIATION

// assign grid 
  int ndim;
  int dimensionsnode[3];
  if(ny==1 && nz==1) //1d
  {
      ndim=1;
      dimensions[0]=nx;
      dimensionsnode[0]=nnx;
      coordinates[0]=nodex;
  }
  else if(nz==1)     //2d
  {
      ndim=2;

      dimensions[0] = nx;
      dimensions[1] = ny;
      dimensionsnode[0]=nnx;
      dimensionsnode[1]=nny;

      coordinates[0] = nodex;
      coordinates[1] = nodey;
      #ifdef SILO2D_XZPLANE
      coordinates[1] = nodez;
      #endif
  }
  else if(ny==1)     //2d, switch order
  {
      ndim=2;
      
      // How many nodes in each direction? 
      dimensions[0] = nx;
      dimensions[1] = nz;
      dimensionsnode[0]=nnx;
      dimensionsnode[1]=nnz;

      // Assign coordinates to coordinates array 
      coordinates[0] = nodex;
      coordinates[1] = nodey;  //works for spherical-like coordinates

      
      #ifdef SHEARINGBOX
      coordinates[1] = nodez;
      #endif

  }
  else //3d
  {
      ndim=3;
      
      // How many nodes in each direction? 
      dimensions[0] = nx;
      dimensions[1] = ny;
      dimensions[2] = nz;

      dimensionsnode[0]=nnx;
      dimensionsnode[1]=nny;
      dimensionsnode[2]=nnz;

      // Assign coordinates to coordinates array 
      coordinates[0] = nodex;
      coordinates[1] = nodey;
      coordinates[2] = nodez;
  }      
     
  DBoptlist *optList = DBMakeOptlist(3);
  float ftime=(float)time;
  DBAddOption(optList, DBOPT_DTIME, (void*)&time);
  DBAddOption(optList, DBOPT_TIME, (void*)&ftime);
  DBAddOption(optList, DBOPT_CYCLE, (void*)&nstep);

  // Write out the mesh to the file 
  DBPutQuadmesh(file, "mesh1", coordnames, coordinates,
  		dimensionsnode, ndim, DB_DOUBLE, DB_NONCOLLINEAR, optList);

  // Write scalars to the file
  DBPutQuadvar1(file, "entr","mesh1", entr,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  
  DBPutQuadvar1(file, "entrlnT","mesh1", entrlnT,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  
  DBPutQuadvar1(file, "rho","mesh1", rho,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  
  DBPutQuadvar1(file, "entropyinv","mesh1", entropyinv,
  		dimensions, ndim, NULL, 0, 
		DB_INT, DB_ZONECENT, optList);
  
  DBPutQuadvar1(file, "uint","mesh1", uint,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "temp","mesh1", temp,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "omega","mesh1", Omega,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);


  DBPutQuadvar1(file, "u0","mesh1", u0,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "u1","mesh1", u1,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "u2","mesh1", u2,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "u3","mesh1", u3,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  
  DBPutQuadvar1(file, "ffinv","mesh1", ffinv,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  
  DBPutQuadvar1(file, "lorentz","mesh1", lorentz,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  
  DBPutQuadvar1(file, "lorentz_perp","mesh1", lorentz_perp,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "lorentz_par","mesh1", lorentz_par,
		
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);


  DBPutQuadvar1(file, "muBe","mesh1", muBe,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "Be","mesh1", Be,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "expansion","mesh1", expansion,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

#ifdef PRINTVISCHEATINGTOSILO
  DBPutQuadvar1(file, "dtau","mesh1", dtauarr,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "vischeating","mesh1", vischeat,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "deltae","mesh1", deltae,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "vischeatnege","mesh1", vischeatnege,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "vischeatnegi","mesh1", vischeatnegi,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
#endif

#ifdef PRINTCOULOMBTOSILO
  DBPutQuadvar1(file, "coulomb","mesh1", coulomb,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
#endif

#ifdef PRINTGAMMATOSILO
  DBPutQuadvar1(file, "gammagas","mesh1", gammag,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
#endif



#ifdef EVOLVEELECTRONS
  DBPutQuadvar1(file, "tempe","mesh1", tempe,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "tempi","mesh1", tempi,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "ue","mesh1", ue,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "ui","mesh1", ui,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  
#ifdef RELELECTRONS
  DBPutQuadvar1(file, "urelel","mesh1", urelel,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "gammabrk","mesh1", gammabrk,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "uratio_tot","mesh1", uratio_tot,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "uratio_th","mesh1", uratio_th,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "nerelel","mesh1", nrelel,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
 
  DBPutQuadvar1(file, "neth","mesh1", neth,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  
  DBPutQuadvar1(file, "G0relel","mesh1", G0relel,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "G0releliC","mesh1", G0icrelel,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "G0relelsyn","mesh1", G0synrelel,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

#endif
#endif
  
#ifdef RADIATION
  DBPutQuadvar1(file, "erad","mesh1", Erad,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "ehat","mesh1", Ehat,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "Gtff","mesh1", Gtff,
  		dimensions, ndim, NULL, 0, 
	      DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "taucoupling","mesh1", taucoupling,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "tausca_theta","mesh1", tausca,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "tauabs_theta","mesh1", tauabs,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "taueff_theta","mesh1", taueff,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "tausca_r","mesh1", tauscar,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "tauabs_r","mesh1", tauabsr,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "taueff_r","mesh1", taueffr,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "entrrad","mesh1", entrrad,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "NH","mesh1", NHr,
               dimensions, ndim, NULL, 0, 
               DB_DOUBLE, DB_ZONECENT, optList);

  #ifdef EVOLVEPHOTONNUMBER
  DBPutQuadvar1(file, "nph","mesh1", Nph,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  #endif
  DBPutQuadvar1(file, "trad","mesh1", trad,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "tradlte","mesh1", tradlte,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  #ifdef COMPTONIZATION
  DBPutQuadvar1(file, "G0thiC","mesh1", G0icth,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  #endif
#endif

#ifdef MAGNFIELD

  DBPutQuadvar1(file, "bsq","mesh1", bsq,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "B1","mesh1", B1comp,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "B2","mesh1", B2comp,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "B3","mesh1", B3comp,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "OmegaF","mesh1", OmegaF,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "phi","mesh1", phi,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "Qtheta","mesh1", Qtheta,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "Qphi","mesh1", Qphi,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "Bangle","mesh1", Bangle,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "beta","mesh1", beta,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "sigma","mesh1", sigma,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "betainv","mesh1", betainv,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "divB","mesh1", divB,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  DBPutQuadvar1(file, "jdens","mesh1", jdens,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  #ifdef MIMICDYNAMO
  DBPutQuadvar1(file, "phidyn","mesh1", phidyn,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  DBPutQuadvar1(file, "Avecdyn","mesh1", Avecdyn,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  #endif
#endif


  // Write vectors to file
  char *names[3];  
  ldouble *vels[3];

  //velocity
  vels[0]=vx;
  vels[1]=vy;
  vels[2]=vz;
#ifdef SILO2D_XZPLANE
#if (MYCOORDS!=CYLCOORDS)
  vels[1]=vz;
#endif
  /*
   DBPutQuadvar1(file, "velocity_z","mesh1", vy,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  */
#endif

  names[0] = strdup("vel1");
  names[1] = strdup("vel2");
  names[2] = strdup("vel3");
  DBPutQuadvar(file, "velocity","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  //en. flux
  vels[0]=Edotx;
  vels[1]=Edoty;
  vels[2]=Edotz;
#ifdef SILO2D_XZPLANE
#if (MYCOORDS!=CYLCOORDS)

  vels[1]=Edotz;
#endif
  DBPutQuadvar1(file, "en_flux_z","mesh1", Edoty,
  		dimensions, ndim, NULL, 0, 
	      DB_DOUBLE, DB_ZONECENT, optList);
#endif

  names[0] = strdup("Tit1");
  names[1] = strdup("Tit2");
  names[2] = strdup("Tit3");
  DBPutQuadvar(file, "en_flux","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

#ifdef MAGNFIELD
  //magn field
  vels[0]=Bx;
  vels[1]=By;
  vels[2]=Bz;
#ifdef SILO2D_XZPLANE
#if (MYCOORDS!=CYLCOORDS) 
  vels[1]=Bz;
#endif
  /*  DBPutQuadvar1(file, "magn_field_z","mesh1", By,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  */
#endif
  names[0] = strdup("B1");
  names[1] = strdup("B2");
  names[2] = strdup("B3");
  DBPutQuadvar(file, "magn_field","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);



#ifdef BATTERY
  //magn field growth due to battery

  DBPutQuadvar1(file, "phibat","mesh1", phibat,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  vels[0]=dBxdtbat;
  vels[1]=dBydtbat;
  vels[2]=dBzdtbat;
#ifdef SILO2D_XZPLANE
#if (MYCOORDS!=CYLCOORDS)
  vels[1]=dBzdtbat;
#endif
 DBPutQuadvar1(file, "battery_dBdt_z","mesh1", dBydtbat,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
#endif
  names[0] = strdup("battery_dB1dt");
  names[1] = strdup("battery_dB2dt");
  names[2] = strdup("battery_dB3dt");
  DBPutQuadvar(file, "battery_dBdt","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
#endif

#ifdef MIMICDYNAMO
  //magn field due to dynamo
  vels[0]=Bxdyn;
  vels[1]=Bydyn;
  vels[2]=Bzdyn;
#ifdef SILO2D_XZPLANE
#if (MYCOORDS!=CYLCOORDS) 
  vels[1]=Bzdyn;
#endif
  DBPutQuadvar1(file, "magn_field_dyn_z","mesh1", Bydyn,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
  #endif
  names[0] = strdup("B1dyn");
  names[1] = strdup("B2dyn");
  names[2] = strdup("B3dyn");
  DBPutQuadvar(file, "magn_field_dyn","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
	       DB_DOUBLE, DB_ZONECENT, optList);
  #endif


#endif //MAGNFIELD

  #ifdef RADIATION 
  //radiative flux
  vels[0]=Fx;
  vels[1]=Fy;
  vels[2]=Fz;
#ifdef SILO2D_XZPLANE
#if (MYCOORDS!=CYLCOORDS)
  vels[1]=Fz;
#endif
  DBPutQuadvar1(file, "rad_flux_z","mesh1", Fy,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
#endif
  names[0] = strdup("F1");
  names[1] = strdup("F2");
  names[2] = strdup("F3");
  DBPutQuadvar(file, "rad_flux","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);

  //rad rest frame velocity
  vels[0]=uradx;
  vels[1]=urady;
  vels[2]=uradz;
#ifdef SILO2D_XZPLANE
#if (MYCOORDS!=CYLCOORDS)
 
  vels[1]=uradz;
#endif
  DBPutQuadvar1(file, "rad_vel_z","mesh1", uradz,
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
#endif
  names[0] = strdup("urad1");
  names[1] = strdup("urad2");
  names[2] = strdup("urad3");
  DBPutQuadvar(file, "rad_vel","mesh1", 3, names, vels, 
  		dimensions, ndim, NULL, 0, 
		DB_DOUBLE, DB_ZONECENT, optList);
#endif
 
  // Close the Silo file and free memory
  DBClose(file);

  free(nodex);
  free(nodey);
  free(nodez);

  free(rho);
  free(uint);
  free(entr);
  free(temp);
  free(Omega);
  free(muBe);
  free(Be);
  free(expansion);
  free(NHr);
  free(entrlnT);
  free(entropyinv);
  
  free(u0);
  free(u1);
  free(u2);
  free(u3);
  free(lorentz);

  free(ffinv);
  free(lorentz_perp);
  free(lorentz_par);
  free(vx);
  free(vy);
  free(vz);
  free(Edotx);
  free(Edoty);
  free(Edotz);

#ifdef MAGNFIELD
  free(Bangle);
  free(bsq);
  free(B1comp);
  free(B2comp);
  free(B3comp);
  free(Bx);
  free(By);
  free(Bz);
  free(phi);
  free(phi_accum);
  free(beta); 
  free(sigma);
  free(sigmaw);
  free(betainv); 
  free(divB);
  free(jdens);
  free(OmegaF);
  free(Qtheta);
  free(Qphi);

  #ifdef BATTERY
  free(dBxdtbat);
  free(dBydtbat);
  free(dBzdtbat);
  free(phibat);
  #endif

  #ifdef MIMICDYNAMO
  free(Bxdyn);
  free(Bydyn);
  free(Bzdyn);
  free(phidyn);
  free(Avecdyn);
  #endif
#endif
 
#ifdef PRINTVISCHEATINGTOSILO
  free(deltae);
  free(dtauarr);
  free(vischeat);
  free(vischeatnege);
  free(vischeatnegi);
#endif

#ifdef PRINTCOULOMBTOSILO
  free(coulomb);
#endif

#ifdef PRINTGAMMATOSILO
  free(gammag);
#endif


#ifdef EVOLVEELECTRONS
  free(tempe);
  free(tempi);
  free(ue);
  free(ui);
#ifdef RELELECTRONS
  free(gammabrk);
  free(urelel);
  free(nrelel); 
  free(neth);
  free(uratio_tot);
  free(uratio_th);
  free(G0relel);
  free(G0icrelel);
  free(G0synrelel);
#endif
#endif


#ifdef RADIATION
  free(trad);
  free(tradlte);

  free(tausca);
  free(tauabs);
  free(taueff);
  free(taucoupling);
  free(Erad);
  free(Ehat);
  
  free(Fx);
  free(Fy);
  free(Fz);
  free(Nph);
  free(uradx);
  free(urady);
  free(uradz);
  free(Gtff);
  free(G0icth);

  free(tauscar);
  free(tauabsr);
  free(taueffr);

  free(tauscar2);
  free(tauabsr2);
  free(taueffr2);

  free(entrrad);
#endif

  return (0);
}

#endif //NOSILO
