/*! \file magn.c
 \brief Magnetic field related routines
*/

#include "ko.h"

//***********************************************************************
/* calculate both magnetic field four-vectors and bsq knowing gas four-velocity ucov */
//***********************************************************************

void calc_bcon_bcov_bsq_from_4vel(ldouble *pr, ldouble *ucon, ldouble *ucov, void* ggg,
				  ldouble *bcon, ldouble *bcov, ldouble *bsq)
{

  int j;
  struct geometry *geom
  = (struct geometry *) ggg;

  // First calculate bcon0
  bcon[0] = pr[B1]*ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3] ;
  
  // Then spatial components of bcon
  
#ifdef NONRELMHD
  for(j = 1; j < 4; j++)
    bcon[j] = pr[B1-1+j]; //b^i=B^i

#else  // relativistic case
  
  ldouble u0inv = 1. / ucon[0];

#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
  for(j=1;j<4;j++)
    bcon[j] = (pr[B1-1+j] + bcon[0] * ucon[j]) * u0inv ;
  
#endif //NONRELMHD
  
  // Convert to bcov and calculate bsq
  indices_21(bcon, bcov, geom->gg);
  *bsq = dotB(bcon, bcov);

  return ;
}

//***********************************************************************
/* calculate magnetic field four-vector knowing gas four-velocity ucov */
//***********************************************************************

void calc_bcon_4vel(double *pr, double *ucon, double *ucov, double *bcon)
{
  int j;
  
  bcon[0] = pr[B1]*ucov[1] + pr[B2]*ucov[2] + pr[B3]*ucov[3] ;

#ifdef NONRELMHD
  for(j=1;j<4;j++)
    bcon[j] = pr[B1-1+j]; //b^i=B^i

#else  // relativistic case

  ldouble u0inv = 1. / ucon[0];

#ifdef APPLY_OMP_SIMD
 //#pragma omp simd
#endif
  for(j=1;j<4;j++)
      bcon[j] = (pr[B1-1+j] + bcon[0]*ucon[j]) * u0inv ;
  
#endif //NONRELMHD

  return ;
}


//***********************************************************************
/* calculate B^i from bcon and gas four-velocity ucon */
//***********************************************************************

void calc_Bcon_4vel(double *pr, double *ucon, double *bcon, double *Bcon)
{
  int j;
  
  Bcon[0]=0.;
  
#ifdef NONRELMHD
  for(j=1;j<4;j++)
  {
    Bcon[j] = bcon[j];
  }
#else
  for(j=1;j<4;j++)
  {
    Bcon[j] = bcon[j]*ucon[0] - bcon[0]*ucon[j];
  }
#endif //NONRELMHD
  
  return ;
}

//***********************************************************************
/* calculate magnetic field four-vector from primitives */
//***********************************************************************

void calc_bcon_prim(double *pp, double *bcon, void* ggg)
{
  int j;
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble ucon[4],ucov[4];

  calc_ucon_ucov_from_prims(pp, geom, ucon, ucov);

  calc_bcon_4vel(pp,ucon,ucov,bcon);

  return;
}

//***********************************************************************
/* calculate B^i from bcon and primitives */
//***********************************************************************

void calc_Bcon_prim(double *pp, double *bcon, double *Bcon, void* ggg)
{
  int j;
  
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble ucon[4],ucov[4];

  calc_ucon_ucov_from_prims(pp, geom, ucon, ucov);
  

  Bcon[0]=0.;

#ifdef NONRELMHD
   for(j=1;j<4;j++)
     Bcon[j] = bcon[j];  
#else
  for(j=1;j<4;j++)
    {
      Bcon[j] = bcon[j]*ucon[0] - bcon[0]*ucon[j];
    }
#endif //NONRELMHD

  return;
}

/***********************************************************************************************/
/* wrappers for missing cells / dimensions */
/***********************************************************************************************/

int fl_x(int i)
{
  if(NX==1) return 0;
  return i;
}

int fl_y(int i)
{
  if(NY==1) return 0;
  return i;
}

int fl_z(int i)
{
  if(NZ==1) return 0;
  return i;
}

int flx_x(int i)
{
  return i;
}

int flx_y(int i)
{
  return fl_y(i);
}

int flx_z(int i)
{
  return fl_z(i);
}

int fly_x(int i)
{
  return fl_x(i);
}

int fly_y(int i)
{
  return i;
}

int fly_z(int i)
{
  return fl_z(i);
}

int flz_x(int i)
{
  return fl_x(i);
}

int flz_y(int i)
{
  return fl_y(i);
}

int flz_z(int i)
{
  return i;
}


/***********************************************************************************************/
/* flux constrained transport */
/***********************************************************************************************/

// B^i = \dF^{it}
// E_i = - [ijk] v^j B^k  , such that (\detg B^i),t = - (\detg(B^i v^j - B^j v^i)),j
//                                                  = - (\detg [ijk] E_k),j = ([ijk] emf[k]),j
      
// -> E_1 = v^3 B^2 - v^2 B^3
// -> E_2 = v^1 B^3 - v^3 B^1
// -> E_3 = v^2 B^1 - v^1 B^2

// emf[i] = - \detg E_i

// And notice that Fj[Bi] = \dF^{ij} = B^i v^j - B^j v^i , where j=dir

// so:
// emf_1 = B^3 v^2 - B^2 v^3 = F2[B3] or -F3[B2]
// emf_2 = B^1 v^3 - B^3 v^1 = F3[B1] or -F1[B3]
// emf_3 = B^2 v^1 - B^1 v^2 = F1[B2] or -F2[B1]

// Notice only 6 independent ways.  The diagonal terms vanish (e.g. Fi[Bi]=0).
/***********************************************************************************************/


int
flux_ct()
{
#ifdef MAGNFIELD

  //requires GDETIN = 1
  if(GDETIN==0)
    {
      my_err("MAGNFIELD requires GDETIN==1\n");
      exit(0);
    }

  //TOTH algorithm from HARM's fluxct.c
  ldouble coefemf[4];

  if((NY>1)&&(NZ>1)) coefemf[1]=0.25;
  else coefemf[1]=0.5; 
  if((NX>1)&&(NZ>1)) coefemf[2]=0.25;
  else coefemf[2]=0.5; 
  if((NX>1)&&(NY>1)) coefemf[3]=0.25;
  else coefemf[3]=0.5; 
  
  int ix,iy,iz,iv,ii;
  
#pragma omp parallel for private(ix,iy,iz,iv) schedule (static)

  //calculating EMF at corners
  for(ii=0;ii<Nloop_4;ii++) 
    {
      ix=loop_4[ii][0];
      iy=loop_4[ii][1];
      iz=loop_4[ii][2]; 
      
      ////////////////////
      // EMF1
      ////////////////////
      
#if((NY>1)||(NZ>1))
      set_emf(1,ix,iy,iz,
	      coefemf[1] * (
                            #if(NY>1)
			    + get_ub(flby,B3,fly_x(ix),fly_y(iy),fly_z(iz),1) 
			    + get_ub(flby,B3,fly_x(ix),fly_y(iy),fly_z(iz-1),1)
                            #endif
                            #if(NZ>1)
			    - get_ub(flbz,B2,flz_x(ix),flz_y(iy),flz_z(iz),2) 
			    - get_ub(flbz,B2,flz_x(ix),flz_y(iy-1),flz_z(iz),2)
                            #endif
			    ));
#else  
      set_emf(1,ix,iy,iz,0.); // not really 0, but differences in emf will be 0
#endif 
      
	////////////////////
	// EMF2
	////////////////////
#if((NX>1)||(NZ>1))
      set_emf(2,ix,iy,iz,
	      coefemf[2] * (
                            #if(NZ>1)
			    + get_ub(flbz,B1,flz_x(ix),flz_y(iy),flz_z(iz),2) 
			    + get_ub(flbz,B1,flz_x(ix-1),flz_y(iy),flz_z(iz),2)
                            #endif
                            #if(NX>1)
			    - get_ub(flbx,B3,flx_x(ix),flx_y(iy),flx_z(iz),0) 
			    - get_ub(flbx,B3,flx_x(ix),flx_y(iy),flx_z(iz-1),0)
                            #endif
			    ));
#else  
      set_emf(2,ix,iy,iz,0.);
#endif 

	////////////////////
	// EMF3
	////////////////////
#if((NX>1)||(NY>1))
      set_emf(3,ix,iy,iz,
	      coefemf[3] * (
                            #if(NX>1)
			    + get_ub(flbx,B2,flx_x(ix),flx_y(iy),flx_z(iz),0) 
			    + get_ub(flbx,B2,flx_x(ix),flx_y(iy-1),flx_z(iz),0)
                            #endif
                            #if(NY>1)
			    - get_ub(flby,B1,fly_x(ix),fly_y(iy),fly_z(iz),1) 
			    - get_ub(flby,B1,fly_x(ix-1),fly_y(iy),fly_z(iz),1)
                            #endif
			    ));
#else  
      set_emf(3,ix,iy,iz,0.);
#endif
    }
  
  //reset certain emfs at the boundaries to ensure stationarity
  adjust_fluxcttoth_emfs();

#pragma omp parallel for private(ix,iy,iz,iv) schedule (static)
  for(ii=0;ii<Nloop_4;ii++) // 0...NX
    {
      ix=loop_4[ii][0];
      iy=loop_4[ii][1];
      iz=loop_4[ii][2]; 

      /////////////////////////////////////
      // F1
      ////////////////////////////////////
#if(NX>1)
	if(iy<NY && iz<NZ) //no need to fill x-face fluxes for iy=NY etc., 
	  {
	    set_ubx(flbx,B1,ix,iy,iz,0.);
	    set_ubx(flbx,B2,ix,iy,iz,0.5 * (get_emf(3,ix,iy,iz) + get_emf(3,ix,iy+1,iz)));
	    set_ubx(flbx,B3,ix,iy,iz,-0.5 * (get_emf(2,ix,iy,iz) + get_emf(2,ix,iy,iz+1)));
	  }
#endif

      /////////////////////////////////////
      // F2
      ////////////////////////////////////
#if(NY>1)
      if(ix<NX && iz<NZ)	
	{
	  set_uby(flby,B1,ix,iy,iz,-0.5 * (get_emf(3,ix,iy,iz) + get_emf(3,ix+1,iy,iz)));
	  set_uby(flby,B2,ix,iy,iz,0.);
	  set_uby(flby,B3,ix,iy,iz,0.5 * (get_emf(1,ix,iy,iz) + get_emf(1,ix,iy,iz+1)));
	}
#endif
      
			    
      /////////////////////////////////////
      // F3
      ////////////////////////////////////
#if(NZ>1)
	if(ix<NX && iy<NY)	
	{
	  set_ubz(flbz,B1,ix,iy,iz,0.5 * (get_emf(2,ix,iy,iz) + get_emf(2,ix+1,iy,iz)));
	  set_ubz(flbz,B2,ix,iy,iz,-0.5 * (get_emf(1,ix,iy,iz) + get_emf(1,ix,iy+1,iz)));
	  set_ubz(flbz,B3,ix,iy,iz,0.);
	}
#endif
    }
  
#endif //MAGNFIELD
      return 0;
}


/***********************************************************************************************/
// resets emf's near the boundaries where corresponding velocities are zero
/***********************************************************************************************/

int adjust_fluxcttoth_emfs()
{

#ifdef CORRECT_POLARAXIS
  int ix,iz;
#ifdef MPI
  if(TJ==0) //upper axis
#endif
    {
      //over all corners at the polar edge
      for(ix=0;ix<=NX;ix++)
	for(iz=0;iz<=NZ;iz++)
	  {
	    set_emf(1,ix,0,iz,0.);
	    set_emf(3,ix,0,iz,0.);
	  }
    }
  
#ifdef MPI
  if(TJ==NTY-1) //lower axis
#endif
    {
      //over all corners at the polar edge
      for(ix=0;ix<=NX;ix++)
	for(iz=0;iz<=NZ;iz++)
	  {
	    set_emf(1,ix,NY,iz,0.);
	    set_emf(3,ix,NY,iz,0.);
	  }
    }

#endif //CORRECTPOLARAXIS

  return 0;
}


/***********************************************************************************************/
//calculates magnetic field from vector potential given in pinput[B1..B3]
/***********************************************************************************************/

int
calc_BfromA(ldouble* pinput, int ifoverwrite)
{
  #ifdef MAGNFIELD

  int ix,iy,iz,iv,ii;
  
  //ANDERW loop_4 seems to be the same as loop_0? 
#pragma omp parallel for private(ix,iy,iz,iv) schedule (static)
  for(ii=0;ii<Nloop_4;ii++) //all corners of the inner domain
    {      
      ix=loop_4[ii][0];
      iy=loop_4[ii][1];
      iz=loop_4[ii][2];
    
      if(NZ==1 && iz>0) continue;
      if(NY==1 && iy>0) continue;

      //calculating A_i on corners by averaging neighbouring cell centers
      ldouble A[3];

      for(iv=0;iv<3;iv++)
	{
#if defined(MONOPOLE_FIELD_CORNERS) // explicit monopole B field solution
	  ldouble xxvec[4],xxvecBL[4],r,th;
	  if (iv==2)
	  {
	    xxvec[0] = global_time;
            xxvec[1] = get_xb(ix,0);
	    xxvec[2] = get_xb(iy,1);
	    xxvec[3] = get_xb(iz,2);
            coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
            th=xxvecBL[2];
	    A[iv] = 1.-cos(th);
	  }
	  else A[iv]=0.;
#elif defined(INIT_MAGN_CORNERS) // define the vecpot on the corners in the init.c
          A[iv]=get_u(pinput,B1+iv,ix,iy,iz);
#else // linear interpolation to  corners
	  if(NY==1 && NZ==1)
	    A[iv]=1./2.*(get_u(pinput,B1+iv,ix,iy,iz) + get_u(pinput,B1+iv,ix-1,iy,iz));

	  if(NY>1 && NZ==1)
	    A[iv]=1./4.*(get_u(pinput,B1+iv,ix,iy,iz) + get_u(pinput,B1+iv,ix,iy-1,iz) + 
			 get_u(pinput,B1+iv,ix-1,iy,iz) + get_u(pinput,B1+iv,ix-1,iy-1,iz));

	  if(NZ>1 && NY==1)
	    A[iv]=1./4.*(get_u(pinput,B1+iv,ix,iy,iz) + get_u(pinput,B1+iv,ix,iy,iz-1) + 
			 get_u(pinput,B1+iv,ix-1,iy,iz) + get_u(pinput,B1+iv,ix-1,iy,iz-1));

	  if(NZ>1 && NY>1)
	    A[iv]=1./8.*(get_u(pinput,B1+iv,ix,iy,iz) + get_u(pinput,B1+iv,ix,iy-1,iz) + 
			 get_u(pinput,B1+iv,ix-1,iy,iz) + get_u(pinput,B1+iv,ix-1,iy-1,iz)
			 +get_u(pinput,B1+iv,ix,iy,iz-1) + get_u(pinput,B1+iv,ix,iy-1,iz-1) 
			 +get_u(pinput,B1+iv,ix-1,iy,iz-1) + get_u(pinput,B1+iv,ix-1,iy-1,iz-1));
#endif 
	  //saving to pvecpot
 	  set_u(pvecpot,B1+iv,ix,iy,iz,A[iv]);
	}            
    } //cell loop
  
  //calculating curl to get B
  //new components of B^i in pvecpot[1...3]
  calc_BfromA_core();
  
  //overwriting vector potential with magnetic fields (e.g., at init)  
  if(ifoverwrite)
    {
      for(ix=0-NGCX;ix<NX+NGCX;ix++)
	for(iy=0-NGCY;iy<NY+NGCY;iy++)
	  for(iz=0-NGCZ;iz<NZ+NGCZ;iz++)
	    {

	      struct geometry geom;
	      fill_geometry(ix,iy,iz,&geom);
      
	      ldouble pp[NV],uu[NV];
	      PLOOP(iv)
		pp[iv]=get_u(p,iv,ix,iy,iz);

	      pp[B1]=get_u(pvecpot,1,ix,iy,iz);
	      pp[B2]=get_u(pvecpot,2,ix,iy,iz);
	      pp[B3]=get_u(pvecpot,3,ix,iy,iz);

	      p2u(pp,uu,&geom);

	      set_u(p,B1,ix,iy,iz,pp[B1]);
	      set_u(p,B2,ix,iy,iz,pp[B2]);
	      set_u(p,B3,ix,iy,iz,pp[B3]);
	      set_u(u,B1,ix,iy,iz,uu[B1]);
	      set_u(u,B2,ix,iy,iz,uu[B2]);
	      set_u(u,B3,ix,iy,iz,uu[B3]);     
	}
    }

#endif //MAGNFIELD

  return 0;
}


/***********************************************************************************************/
//calculates B-field from A given on corners in B1-B3 primitives of pvecpot
//new components of B^i in pvecpot[1...3]
/***********************************************************************************************/

int
calc_BfromA_core()
{
  #ifdef MAGNFIELD

  int ix,iy,iz,iv,ii;
  
  if(NY==1 && NZ==1)
    {
      my_err("1D calc_BfromA_core() not implemented.\n"); exit(-1);
    }

  if(NY>1 && NZ>1)
    if(PROCID==0) printf("calc_BfromA_core(): warning: assumes d/dz A_i = 0. \n");

#pragma omp parallel for private(ix,iy,iz,iv) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain, need vecpot defined on domain + 1 layer
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 
  
      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);

      ldouble B[4];
      ldouble dA[4][4];

      dA[1][1]=dA[2][2]=dA[3][3]=0.;

      if(TNZ==1) /* flux-ct-compatible, axisymmetric! */
	{
	  
	  dA[1][2] = -(get_u(pvecpot,B1,ix,iy,iz) - get_u(pvecpot,B1,ix,iy+1,iz)
		       + get_u(pvecpot,B1,ix+1,iy,iz) - get_u(pvecpot,B1,ix+1,iy+1,iz))/(2.*get_size_x(iy,1)) ;
	  dA[1][1] = -(get_u(pvecpot,B1,ix,iy,iz) + get_u(pvecpot,B1,ix,iy+1,iz)
		       - get_u(pvecpot,B1,ix+1,iy,iz) - get_u(pvecpot,B1,ix+1,iy+1,iz))/(2.*get_size_x(ix,0)) ;
	  dA[1][3] = 0.;

	  dA[2][2] = -(get_u(pvecpot,B2,ix,iy,iz) - get_u(pvecpot,B2,ix,iy+1,iz)
		       + get_u(pvecpot,B2,ix+1,iy,iz) - get_u(pvecpot,B2,ix+1,iy+1,iz))/(2.*get_size_x(iy,1)) ;
	  dA[2][1] = -(get_u(pvecpot,B2,ix,iy,iz) + get_u(pvecpot,B2,ix,iy+1,iz)
		       - get_u(pvecpot,B2,ix+1,iy,iz) - get_u(pvecpot,B2,ix+1,iy+1,iz))/(2.*get_size_x(ix,0)) ;
	  dA[2][3] = 0.;

	  dA[3][2] = -(get_u(pvecpot,B3,ix,iy,iz) - get_u(pvecpot,B3,ix,iy+1,iz)
		       + get_u(pvecpot,B3,ix+1,iy,iz) - get_u(pvecpot,B3,ix+1,iy+1,iz))/(2.*get_size_x(iy,1)) ;
	  dA[3][1] = -(get_u(pvecpot,B3,ix,iy,iz) + get_u(pvecpot,B3,ix,iy+1,iz)
		       - get_u(pvecpot,B3,ix+1,iy,iz) - get_u(pvecpot,B3,ix+1,iy+1,iz))/(2.*get_size_x(ix,0)) ;
	  dA[3][3] = 0.;

	}
      else //generalized to three dimensions
	{

	  dA[1][1] = (get_u(pvecpot,B1,ix+1,iy,iz) - get_u(pvecpot,B1,ix,iy,iz) +
		      get_u(pvecpot,B1,ix+1,iy+1,iz) - get_u(pvecpot,B1,ix,iy+1,iz) +
		      get_u(pvecpot,B1,ix+1,iy,iz+1) - get_u(pvecpot,B1,ix,iy,iz+1) +
		      get_u(pvecpot,B1,ix+1,iy+1,iz+1) - get_u(pvecpot,B1,ix,iy+1,iz+1) )/(4.*get_size_x(ix,0)) ;
	  dA[1][2] = (get_u(pvecpot,B1,ix,iy+1,iz) - get_u(pvecpot,B1,ix,iy,iz) +
		      get_u(pvecpot,B1,ix+1,iy+1,iz) - get_u(pvecpot,B1,ix+1,iy,iz) +
		      get_u(pvecpot,B1,ix,iy+1,iz+1) - get_u(pvecpot,B1,ix,iy,iz+1) +
		      get_u(pvecpot,B1,ix+1,iy+1,iz+1) - get_u(pvecpot,B1,ix+1,iy,iz+1))/(4.*get_size_x(iy,1)) ;
	  dA[1][3] = (get_u(pvecpot,B1,ix,iy,iz+1) - get_u(pvecpot,B1,ix,iy,iz) +
		      get_u(pvecpot,B1,ix+1,iy,iz+1) - get_u(pvecpot,B1,ix+1,iy,iz) +
		      get_u(pvecpot,B1,ix,iy+1,iz+1) - get_u(pvecpot,B1,ix,iy+1,iz) +
		      get_u(pvecpot,B1,ix+1,iy+1,iz+1) - get_u(pvecpot,B1,ix+1,iy+1,iz))/(4.*get_size_x(iz,2)) ;

	  dA[2][1] = (get_u(pvecpot,B2,ix+1,iy,iz) - get_u(pvecpot,B2,ix,iy,iz) +
		      get_u(pvecpot,B2,ix+1,iy+1,iz) - get_u(pvecpot,B2,ix,iy+1,iz) +
		      get_u(pvecpot,B2,ix+1,iy,iz+1) - get_u(pvecpot,B2,ix,iy,iz+1) +
		      get_u(pvecpot,B2,ix+1,iy+1,iz+1) - get_u(pvecpot,B2,ix,iy+1,iz+1) )/(4.*get_size_x(ix,0)) ;
	  dA[2][2] = (get_u(pvecpot,B2,ix,iy+1,iz) - get_u(pvecpot,B2,ix,iy,iz) +
		      get_u(pvecpot,B2,ix+1,iy+1,iz) - get_u(pvecpot,B2,ix+1,iy,iz) +
		      get_u(pvecpot,B2,ix,iy+1,iz+1) - get_u(pvecpot,B2,ix,iy,iz+1) +
		      get_u(pvecpot,B2,ix+1,iy+1,iz+1) - get_u(pvecpot,B2,ix+1,iy,iz+1))/(4.*get_size_x(iy,1)) ;
	  dA[2][3] = (get_u(pvecpot,B2,ix,iy,iz+1) - get_u(pvecpot,B2,ix,iy,iz) +
		      get_u(pvecpot,B2,ix+1,iy,iz+1) - get_u(pvecpot,B2,ix+1,iy,iz) +
		      get_u(pvecpot,B2,ix,iy+1,iz+1) - get_u(pvecpot,B2,ix,iy+1,iz) +
		      get_u(pvecpot,B2,ix+1,iy+1,iz+1) - get_u(pvecpot,B2,ix+1,iy+1,iz))/(4.*get_size_x(iz,2)) ;

	   dA[3][1] = (get_u(pvecpot,B3,ix+1,iy,iz) - get_u(pvecpot,B3,ix,iy,iz) +
		      get_u(pvecpot,B3,ix+1,iy+1,iz) - get_u(pvecpot,B3,ix,iy+1,iz) +
		      get_u(pvecpot,B3,ix+1,iy,iz+1) - get_u(pvecpot,B3,ix,iy,iz+1) +
		      get_u(pvecpot,B3,ix+1,iy+1,iz+1) - get_u(pvecpot,B3,ix,iy+1,iz+1) )/(4.*get_size_x(ix,0)) ;
	  dA[3][2] = (get_u(pvecpot,B3,ix,iy+1,iz) - get_u(pvecpot,B3,ix,iy,iz) +
		      get_u(pvecpot,B3,ix+1,iy+1,iz) - get_u(pvecpot,B3,ix+1,iy,iz) +
		      get_u(pvecpot,B3,ix,iy+1,iz+1) - get_u(pvecpot,B3,ix,iy,iz+1) +
		      get_u(pvecpot,B3,ix+1,iy+1,iz+1) - get_u(pvecpot,B3,ix+1,iy,iz+1))/(4.*get_size_x(iy,1)) ;
	  dA[3][3] = (get_u(pvecpot,B3,ix,iy,iz+1) - get_u(pvecpot,B3,ix,iy,iz) +
		      get_u(pvecpot,B3,ix+1,iy,iz+1) - get_u(pvecpot,B3,ix+1,iy,iz) +
		      get_u(pvecpot,B3,ix,iy+1,iz+1) - get_u(pvecpot,B3,ix,iy+1,iz) +
		      get_u(pvecpot,B3,ix+1,iy+1,iz+1) - get_u(pvecpot,B3,ix+1,iy+1,iz))/(4.*get_size_x(iz,2)) ;
	}

      B[1]=(dA[2][3] - dA[3][2])/geom.gdet;
      B[2]=(dA[3][1] - dA[1][3])/geom.gdet;
      B[3]=(dA[1][2] - dA[2][1])/geom.gdet;

      set_u(pvecpot,1,ix,iy,iz,B[1]);
      set_u(pvecpot,2,ix,iy,iz,B[2]);
      set_u(pvecpot,3,ix,iy,iz,B[3]);
     
    } //cell loop
  
#endif //MAGNFIELD

  return 0;
}


/***********************************************************************************************/
/** calculates div B for given cell *****************************************************/
/***********************************************************************************************/

ldouble
calc_divB(int ix,int iy,int iz)
{
  if(!if_indomain(ix,iy,iz)) return 0.; //do not calculate in ghost cells
  
  ldouble divB;
  
  if(NZ==1)
    {
      //this is corner based, but uses cell centered values 
      divB = (pick_gdet(ix,iy,iz)*get_u(p,B1,ix,iy,iz) + pick_gdet(ix,iy-1,iz)*get_u(p,B1,ix,iy-1,iz) 
	      - pick_gdet(ix-1,iy,iz)*get_u(p,B1,ix-1,iy,iz) - pick_gdet(ix-1,iy-1,iz)*get_u(p,B1,ix-1,iy-1,iz))/(2.*(get_x(ix+1,0)-get_x(ix,0)))
	+ (pick_gdet(ix,iy,iz)*get_u(p,B2,ix,iy,iz) + pick_gdet(ix-1,iy,iz)*get_u(p,B2,ix-1,iy,iz) 
	   - pick_gdet(ix,iy-1,iz)*get_u(p,B2,ix,iy-1,iz) - pick_gdet(ix-1,iy-1,iz)*get_u(p,B2,ix-1,iy-1,iz))/(2.*(get_x(iy+1,1)-get_x(iy,1)));
    }

  if(NZ>1)
    {
      divB = (pick_gdet(ix,iy,iz)*get_u(p,B1,ix,iy,iz) + pick_gdet(ix,iy-1,iz)*get_u(p,B1,ix,iy-1,iz) 
	      - pick_gdet(ix-1,iy,iz)*get_u(p,B1,ix-1,iy,iz) - pick_gdet(ix-1,iy-1,iz)*get_u(p,B1,ix-1,iy-1,iz)
	      + pick_gdet(ix,iy,iz-1)*get_u(p,B1,ix,iy,iz-1) + pick_gdet(ix,iy-1,iz-1)*get_u(p,B1,ix,iy-1,iz-1) 
	      - pick_gdet(ix-1,iy,iz-1)*get_u(p,B1,ix-1,iy,iz-1) - pick_gdet(ix-1,iy-1,iz-1)*get_u(p,B1,ix-1,iy-1,iz-1))
	/(4.*(get_x(ix,0)-get_x(ix-1,0)))
	+(pick_gdet(ix,iy,iz)*get_u(p,B2,ix,iy,iz) + pick_gdet(ix-1,iy,iz)*get_u(p,B2,ix-1,iy,iz) 
	  - pick_gdet(ix,iy-1,iz)*get_u(p,B2,ix,iy-1,iz) - pick_gdet(ix-1,iy-1,iz)*get_u(p,B2,ix-1,iy-1,iz)
	  +pick_gdet(ix,iy,iz-1)*get_u(p,B2,ix,iy,iz-1) + pick_gdet(ix-1,iy,iz-1)*get_u(p,B2,ix-1,iy,iz-1) 
	  - pick_gdet(ix,iy-1,iz-1)*get_u(p,B2,ix,iy-1,iz-1) - pick_gdet(ix-1,iy-1,iz-1)*get_u(p,B2,ix-1,iy-1,iz-1))
	/(4.*(get_x(iy,1)-get_x(iy-1,1)))
	+(pick_gdet(ix,iy,iz)*get_u(p,B3,ix,iy,iz) + pick_gdet(ix-1,iy,iz)*get_u(p,B3,ix-1,iy,iz) 
	  - pick_gdet(ix,iy,iz-1)*get_u(p,B3,ix,iy,iz-1) - pick_gdet(ix-1,iy,iz-1)*get_u(p,B3,ix-1,iy,iz-1)
	  +pick_gdet(ix,iy-1,iz)*get_u(p,B3,ix,iy-1,iz) + pick_gdet(ix-1,iy-1,iz)*get_u(p,B3,ix-1,iy-1,iz) 
	  - pick_gdet(ix,iy-1,iz-1)*get_u(p,B3,ix,iy-1,iz-1) - pick_gdet(ix-1,iy-1,iz-1)*get_u(p,B3,ix-1,iy-1,iz-1))
	/(4.*(get_x(iz,2)-get_x(iz-1,2)));
    }
   
  divB/=pick_gdet(ix,iy,iz);

  return divB;  
}


/***********************************************************************************************/
//calculate Q_theta MRI parameter
/***********************************************************************************************/

int
calc_Qthetaphi(int ix, int iy, int iz,ldouble *Qtheta,ldouble *Qphi)
{
  if(!doingavg)
    {
      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);

      ldouble rho=get_u(p,RHO,ix,iy,iz);
      ldouble bcon[4];
      calc_bcon_prim(&get_u(p,0,ix,iy,iz),bcon,&geom);
      ldouble ucon[4];
      ucon[1]=get_u(p,VX,ix,iy,iz);
      ucon[2]=get_u(p,VY,ix,iy,iz);
      ucon[3]=get_u(p,VZ,ix,iy,iz);
      conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
      ldouble Omega = ucon[3]/ucon[0];
      if(Omega==0.) Omega=BIG;
      ldouble dxth=get_xb(iy+1,1)-get_xb(iy,1);
      ldouble dxph=get_xb(iz+1,2)-get_xb(iz,2);
 
      *Qtheta = 2.*M_PI/fabs(Omega)/dxth*fabs(bcon[2])/sqrt(rho);
      *Qphi = 2.*M_PI/fabs(Omega)/dxph*fabs(bcon[3])/sqrt(rho);
    }
  else //doingavg, in BL
    {
      struct geometry geomBL;
      fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

      ldouble bcon2 = get_uavg(pavg,AVGBCON(2),ix,iy,iz);
      ldouble bcon3 = get_uavg(pavg,AVGBCON(3),ix,iy,iz);
      ldouble Omega = get_uavg(pavg,AVGUCON(3),ix,iy,iz)/get_uavg(pavg,AVGUCON(0),ix,iy,iz);
      ldouble rho = get_uavg(pavg,RHO,ix,iy,iz);
      ldouble dxth,dxph;
      ldouble xx1[4],xx2[4],dx[3];

      // ANDREW TODO -- do we need need this to be in BL? 
      if(OUTCOORDS==BLCOORDS)
      {
        get_cellsize_out(ix, iy, iz, dx);
	dxth = dx[1];
	dxph = dx[2];
      }
      else
      {
      
        xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_x(iz,2);
        xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_xb(iy+1,1);xx2[3]=get_x(iz,2);
        coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
        coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
        dxth=fabs(xx2[2]-xx1[2]);
        xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_xb(iz,2);
        xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_x(iy,1);xx2[3]=get_xb(iz,2);
        coco_N(xx1,xx1,MYCOORDS,BLCOORDS);
        coco_N(xx2,xx2,MYCOORDS,BLCOORDS);
        dxph=fabs(xx2[3]-xx1[3]);
      }
      *Qtheta = 2.*M_PI/fabs(Omega)/dxth*fabs(bcon2)/sqrt(rho);
      *Qphi= 2.*M_PI/fabs(Omega)/dxph*fabs(bcon3)/sqrt(rho);
    }
  return 0;
}


/***********************************************************************************************/
//calculates sqrt(g_rr g_phph) , (b^r b^phi) , and b^2
/***********************************************************************************************/
int
calc_angle_brbphibsq(int ix, int iy, int iz, ldouble *brbphi, ldouble *bsq, ldouble *bcon,ldouble *bcov)
{
  int i;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  struct geometry geomBL;
  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

  ldouble pp[NV];

  if(doingavg)
    {
      *bsq = get_uavg(pavg,AVGBSQ,ix,iy,iz);
      ldouble bcon1bcon3 = 
	get_uavg(pavg,AVGBCONBCOV(1,0),ix,iy,iz)*geomBL.GG[3][0] + 
	get_uavg(pavg,AVGBCONBCOV(1,1),ix,iy,iz)*geomBL.GG[3][1] +
	get_uavg(pavg,AVGBCONBCOV(1,2),ix,iy,iz)*geomBL.GG[3][2] +
	get_uavg(pavg,AVGBCONBCOV(1,3),ix,iy,iz)*geomBL.GG[3][3];
      *brbphi = sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*bcon1bcon3;
    }
  else
    {
      PLOOP(i)
	pp[i]=get_u(p,i,ix,iy,iz);
      
      //to BL
#ifdef PRECOMPUTE_MY2OUT
      trans_pmhd_coco_my2out(pp,pp,&geom,&geomBL);
#else      
      trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,geom.xxvec,&geom,&geomBL);
#endif
      
      //b^mu
      calc_bcon_prim(pp, bcon, &geomBL);
      
      //b_mu
      indices_21(bcon,bcov,geomBL.gg); 

      *bsq = dotB(bcon,bcov);
      *brbphi = sqrt(geomBL.gg[1][1]*geomBL.gg[3][3])*bcon[1]*bcon[3];
    }

  return 0;
}


/***********************************************************************************************/
//calculates (b^p b^phi) and b^2
//here cannot work on averages (bpoloidal bphi is not averaged)
/***********************************************************************************************/

int
calc_angle_bpbphibsq(int ix, int iy, int iz, ldouble *bpbphi, ldouble *bsq, ldouble *bcon, ldouble *bcov)
{
  int i;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  struct geometry geomBL;
  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

  ldouble pp[NV];
  PLOOP(i)
    pp[i]=get_u(p,i,ix,iy,iz);
	      
  //to BL
#ifdef PRECOMPUTE_MY2OUT
  trans_pmhd_coco_my2out(pp, pp, &geom, &geomBL);
#else      
  trans_pmhd_coco(pp,pp,MYCOORDS,BLCOORDS,geom.xxvec,&geom,&geomBL);
#endif
  
  //b^mu
  calc_bcon_prim(pp, bcon, &geomBL);
  
  //b_mu
  indices_21(bcon,bcov,geomBL.gg); 

  *bsq = dotB(bcon,bcov);

  ldouble br = sqrt(geomBL.gg[1][1])*bcon[1];
  ldouble bth = sqrt(geomBL.gg[2][2])*bcon[2];
  ldouble bp = sqrt(br*br + bth*bth)*my_sign(bcon[1]);

  *bpbphi = bp*sqrt(geomBL.gg[3][3])*bcon[3];

  return 0;
}


/***********************************************************************************************/
//calculates curl of a 3-vector field from array ptu starting from index idx
//returns to curl[1..4]
//TODO: idx not used
/***********************************************************************************************/

int
calc_curl(ldouble *ptu, int ix, int iy, int iz, void* ggg, ldouble *curl)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble dA[4][4];

  dA[1][1]=dA[2][2]=dA[3][3]=0.;

  //d_1 A_2
  if(NX==1)
    dA[1][2]=0.;
  else
    {
      if(ix>-NG && ix<NX+NG-1)
	dA[1][2]=(get_u(ptu,B2,ix+1,iy,iz)-get_u(ptu,B2,ix-1,iy,iz))/(get_x(ix+1,0)-get_x(ix-1,0));
      if(ix==-NG)
	dA[1][2]=(-3.*get_u(ptu,B2,ix,iy,iz)+4.*get_u(ptu,B2,ix+1,iy,iz)-get_u(ptu,B2,ix,iy,iz+2))/(get_x(ix+2,0)-get_x(ix,0));
      if(ix==NX+NG-1)
	dA[1][2]=(3.*get_u(ptu,B2,ix,iy,iz)-4.*get_u(ptu,B2,ix-1,iy,iz)+get_u(ptu,B2,ix-2,iy,iz))/(get_x(ix,0)-get_x(ix-2,0));
    }

  //d_1 A_3
  if(NX==1)
    dA[1][3]=0.;
  else
    {
      if(ix>-NG && ix<NX+NG-1)
	dA[1][3]=(get_u(ptu,B3,ix+1,iy,iz)-get_u(ptu,B3,ix-1,iy,iz))/(get_x(ix+1,0)-get_x(ix-1,0));
      if(ix==-NG)
	dA[1][3]=(-3.*get_u(ptu,B3,ix,iy,iz)+4.*get_u(ptu,B3,ix+1,iy,iz)-get_u(ptu,B3,ix,iy,iz+2))/(get_x(ix+2,0)-get_x(ix,0));
      if(ix==NX+NG-1)
	dA[1][3]=(3.*get_u(ptu,B3,ix,iy,iz)-4.*get_u(ptu,B3,ix-1,iy,iz)+get_u(ptu,B3,ix-2,iy,iz))/(get_x(ix,0)-get_x(ix-2,0));
    }

  //d_2 A_3
  if(NY==1)
    dA[2][3]=0.;
  else
    {
      if(iy>-NG && iy<NY+NG-1)
	dA[2][3]=(get_u(ptu,B3,ix,iy+1,iz)-get_u(ptu,B3,ix,iy-1,iz))/(get_x(iy+1,1)-get_x(iy-1,1));
      if(iy==-NG)
	dA[2][3]=(-3.*get_u(ptu,B3,ix,iy,iz)+4.*get_u(ptu,B3,ix,iy+1,iz)-get_u(ptu,B3,ix,iy+2,iz))/(get_x(iy+2,1)-get_x(iy,1));
      if(iy==NY+NG-1)
	dA[2][3]=(3.*get_u(ptu,B3,ix,iy,iz)-4.*get_u(ptu,B3,ix,iy-1,iz)+get_u(ptu,B3,ix,iy-2,iz))/(get_x(iy,1)-get_x(iy-2,1));
    }

  //d_2 A_1
  if(NY==1)
    dA[2][1]=0.;
  else
    {
      if(iy>-NG && iy<NY+NG-1)
	dA[2][1]=(get_u(ptu,B1,ix,iy+1,iz)-get_u(ptu,B1,ix,iy-1,iz))/(get_x(iy+1,1)-get_x(iy-1,1));
      if(iy==-NG)
	dA[2][1]=(-3.*get_u(ptu,B1,ix,iy,iz)+4.*get_u(ptu,B1,ix,iy+1,iz)-get_u(ptu,B1,ix,iy+2,iz))/(get_x(iy+2,1)-get_x(iy,1));
      if(iy==NY+NG-1)
	dA[2][1]=(3.*get_u(ptu,B1,ix,iy,iz)-4.*get_u(ptu,B1,ix,iy-1,iz)+get_u(ptu,B1,ix,iy-2,iz))/(get_x(iy,1)-get_x(iy-2,1));
    }

  //d_3 A_2
  if(NZ==1)
    dA[3][2]=0.;
  else
    {
      if(iz>-NG && iz<NZ+NG-1)
	dA[3][2]=(get_u(ptu,B2,ix,iy,iz+1)-get_u(ptu,B2,ix,iy,iz-1))/(get_x(iz+1,2)-get_x(iz-1,2));
      if(iz==-NG)
	dA[3][2]=(-3.*get_u(ptu,B2,ix,iy,iz)+4.*get_u(ptu,B2,ix,iy,iz+1)-get_u(ptu,B2,ix,iy,iz+2))/(get_x(iz+2,2)-get_x(iz,2));
      if(iz==NZ+NG-1)
	dA[3][2]=(3.*get_u(ptu,B2,ix,iy,iz)-4.*get_u(ptu,B2,ix,iy,iz-1)+get_u(ptu,B2,ix,iy,iz-2))/(get_x(iz,2)-get_x(iz-2,2));
    }

  //d_3 A_1
  if(NZ==1)
    dA[3][1]=0.;
  else
    {
      if(iz>-NG && iz<NZ+NG-1)
	dA[3][1]=(get_u(ptu,B1,ix,iy,iz+1)-get_u(ptu,B1,ix,iy,iz-1))/(get_x(iz+1,2)-get_x(iz-1,2));
      if(iz==-NG)
	dA[3][1]=(-3.*get_u(ptu,B1,ix,iy,iz)+4.*get_u(ptu,B1,ix,iy,iz+1)-get_u(ptu,B1,ix,iy,iz+2))/(get_x(iz+2,2)-get_x(iz,2));
      if(iz==NZ+NG-1)
	dA[3][1]=(3.*get_u(ptu,B1,ix,iy,iz)-4.*get_u(ptu,B1,ix,iy,iz-1)+get_u(ptu,B1,ix,iy,iz-2))/(get_x(iz,2)-get_x(iz-2,2));
    }

  //gdet B^i = d_j A_k eps^ijk

  curl[1]=(dA[2][3] - dA[3][2])/geom->gdet;
  curl[2]=(dA[3][1] - dA[1][3])/geom->gdet;
  curl[3]=(dA[1][2] - dA[2][1])/geom->gdet;

  return 0;
}


/***********************************************************************************************/
// mimics alpha-dynamo in axisymmetric sims involving MRI ***************/
/***********************************************************************************************/

int
mimic_dynamo(ldouble dtin)
{
#ifdef MAGNFIELD
#ifdef MIMICDYNAMO
#ifdef BHDISK_PROBLEMTYPE

  int ii;

#pragma omp parallel for schedule (static)
  for(ii=0;ii<Nloop_6;ii++) //inner domain plus 1-cell layer including corners
    {      
      int ix,iy,iz,iv;
      ldouble dt;
      ldouble Aphi,Pk,Omk,Om;
      ldouble xxBL[4],bcon[4],bcov[4],Bp;
      ldouble ucon[4],ucov[4];
      ldouble pp[NV], bsqin, bsqout, ugasin, ugasout;
      ugasin=ugasout=0.;

      ix=loop_6[ii][0];
      iy=loop_6[ii][1];
      iz=loop_6[ii][2]; 
  
      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);

      for(iv=0;iv<NV;iv++)
      {
        pp[iv]=get_u(p,iv,ix,iy,iz);
      }
      
      calc_ucon_ucov_from_prims(pp, &geom, ucon, ucov);

      //let's use ptemp1 to hold the extra vector potential
      set_u(ptemp1,B1,ix,iy,iz,0.);
      set_u(ptemp1,B2,ix,iy,iz,0.);
      set_u(ptemp1,B3,ix,iy,iz,0.);

      int j;
      ldouble angle,bbphi,bsq;

      //magnetic field angle      
      calc_angle_brbphibsq(ix,iy,iz,&bbphi,&bsq,bcon,bcov);
      angle=-bbphi/bsq;

      //BL radius
#ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix, iy, iz, xxBL);
#else           
      coco_N(geom.xxvec,xxBL,MYCOORDS, BLCOORDS);
#endif
      
      //to avoid BH
      if(xxBL[1]<1.0001*rhorizonBL) continue;

      //timescale
      Omk = 1./(BHSPIN+sqrt(xxBL[1]*xxBL[1]*xxBL[1]));
      Pk = 2.*M_PI/Omk;

      //angle
      ldouble facangle=0.;
      if(isfinite(angle))
	{
	  if(angle<-1.) angle=-1.;
	  facangle = my_max(0., (THETAANGLE-angle)/THETAANGLE);
	}

      //radius
      ldouble facradius = step_function(xxBL[1]-1.*rISCOBL,.1*rISCOBL);

      #ifdef MAXRADIUS4DYNAMO
      if(xxBL[1]>2.*MAXRADIUS4DYNAMO) continue;
      facradius *= step_function(MAXRADIUS4DYNAMO-xxBL[1],.1*MAXRADIUS4DYNAMO);
      #endif

      //pre(gas+rad)
      ldouble gamma=GAMMA;
#ifdef CONSISTENTGAMMA
      gamma=pick_gammagas(ix,iy,iz);
#endif
      ldouble gammam1=gamma-1.;
      ldouble prermhd = gammam1*get_u(p,UU,ix,iy,iz);
#ifdef RADIATION
      ldouble Rtt,Ehat,uconr[4],prad;
      calc_ff_Rtt(&get_u(p,0,ix,iy,iz),&Rtt,uconr,&geom);
      Ehat=-Rtt; 
      prad=Ehat/3.;
      prermhd+=prad;		
#endif
      
      //magnetic beta
      ldouble beta = 0.5*bsq/prermhd;
      
      //bsq/rho 
      ldouble betarho = 0.5*bsq/get_u(p,RHO,ix,iy,iz);

      ldouble facmag1 = step_function(1.-beta,0.1);      
      ldouble facmag2 = step_function(.1-betarho,.01);
      
    
#ifdef CALCHRONTHEGO //calculate HR on the go
      int gix=ix+TOI; //global radial index
      if(gix<0) gix=0; 
      if(gix>=TNX) gix=TNX-1;
      ldouble HRDTHETA = scaleth_otg[gix];
      if(HRDTHETA > 0.9*M_PI/2.) HRDTHETA=0.9*M_PI/2.; //to avoid the axis
#else      
      ldouble HRDTHETA = EXPECTEDHR * M_PI/2.;
#endif

      ldouble zH = (M_PI/2. - xxBL[2])/HRDTHETA;
      ldouble zHpow = 1.;
      ldouble faczH = my_max(0.,pow(1. - zH*zH,zHpow));
      ldouble facmagnetization = faczH;
      //ldouble facmagnetization = my_min(facmag1,facmag2);		       
      //ldouble facmagnetization = my_min(faczH,my_min(facmag1,facmag2));

      //the extra vector potential
      ldouble effalpha=ALPHADYNAMO;

      //dynamo proportional to vertical gravity ~ z
      #ifdef ALPHAFLIPSSIGN
      effalpha = - (M_PI/2. - xxBL[2])/(HRDTHETA/2.) * ALPHADYNAMO;  //2 to get average alpha = alphadynamo
      #endif

      //timestep
      dt=dtin;

      ldouble Bphi=get_u(p,B3,ix,iy,iz);
      Aphi = effalpha * (HRDTHETA/(M_PI/2.)) / 0.4 
	* dt / Pk  * xxBL[1] * geom.gg[3][3] * Bphi
	* facradius 
	* facmagnetization 
	* facangle;

      //saving vector potential to ptemp1
      set_u(ptemp1,B3,ix,iy,iz,Aphi); 

      //damping azimuthal component of magnetic field if beta exceeds DAMPBETA
#ifdef DAMPBETA      
     
      ldouble dBphi = - ALPHABETA 
	* facradius
	* faczH
	* dt / Pk 
	* my_max(0.,beta - BETASATURATED) / BETASATURATED 
	* Bphi;
      
      if((dBphi+Bphi)*Bphi<0.) dBphi=-Bphi; //not to overshoot zero 
                                                   
      set_u(p,B3,ix,iy,iz,Bphi+dBphi);    
      
#endif    


    }

      //once the whole array is filled with cell centered A^phi we can 
      //calculate the extra magnetic field returned through pvecpot[1..3]
      calc_BfromA(ptemp1,0);  
   
  //and superimpose it on the original one
#pragma omp parallel for schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain only!
    {
      int ix,iy,iz,iv;
      ldouble pp[NV], ucon[4], ucov[4], bcon[4], bcov[4], bsqin, bsqout, ugasin, ugasout;
      ugasin=ugasout=0.;

      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2];

      struct geometry geom;
      ldouble B[4]; ldouble xxBL[4];

      fill_geometry(ix,iy,iz,&geom);

      //BL radius
      //coco_N(geom.xxvec,xxBL,MYCOORDS, BLCOORDS);

      B[1]=get_u(pvecpot,1,ix,iy,iz);
      B[2]=get_u(pvecpot,2,ix,iy,iz);
      B[3]=get_u(pvecpot,3,ix,iy,iz);


      set_u(p,B1,ix,iy,iz,get_u(p,B1,ix,iy,iz)+B[1]);
      set_u(p,B2,ix,iy,iz,get_u(p,B2,ix,iy,iz)+B[2]);
      set_u(p,B3,ix,iy,iz,get_u(p,B3,ix,iy,iz)+B[3]);

     
      ldouble uutemp[NV];
      p2u(&get_u(p,0,ix,iy,iz),uutemp,&geom);
      set_u(u,B1,ix,iy,iz,uutemp[B1]);
      set_u(u,B2,ix,iy,iz,uutemp[B2]);
      set_u(u,B3,ix,iy,iz,uutemp[B3]);
      
    }

#endif
#endif
#endif

return 0;
}
