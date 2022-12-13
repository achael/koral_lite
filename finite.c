/*! \file finite.c
 \brief Routines related to finite difference and grid
*/

#include "ko.h"

//**********************************************************************
/*! \fn int reduce_order_check(ldouble *pm2,ldouble *pm1,ldouble *p0,ldouble *pp1,ldouble *pp2,int ix,int iy,int iz)
 \brief reduces reconstruction order to first order for small time steps or within certain radius
 \param[in] pm2, pm1, p0, pp1, pp2 Pointers to primitives in cell and neighbors
 \param[in] ix,iy,iz Cell indices
*/
//**********************************************************************
int
reduce_order_check(ldouble *pm2,ldouble *pm1,ldouble *p0,ldouble *pp1,ldouble *pp2,int ix,int iy,int iz)
{
  int reconstrpar=0;
#ifdef REDUCEORDERTEMP
  ldouble t0,tp1,tm1,tmin;
  t0 =calc_PEQ_Tfromurho(p0[UU], p0[RHO]);
  tp1=calc_PEQ_Tfromurho(pp1[UU],pp1[RHO]);
  tm1=calc_PEQ_Tfromurho(pm1[UU],pm1[RHO]);
  tmin=my_min(t0,my_min(tp1,tm1));	  
  if(tmin<REDUCEORDERTEMP)
    reconstrpar=1;
#endif

#ifdef REDUCEORDERRADIUS
  ldouble xxBL[4];
  get_xx_arb(ix,iy,iz,xxBL,BLCOORDS);
  if(xxBL[1]<REDUCEORDERRADIUS)
    reconstrpar=1;
#endif
  
  return reconstrpar;
}

//**********************************************************************
/*! \fn ldouble reduce_minmod_theta(ldouble *pm2,ldouble *pm1,ldouble *p0,ldouble *pp1,ldouble *pp2,int ix,int iy,int iz)
 \brief reduces minmod theta value to 1 near axis and inner boundary
 \param[in] pm2, pm1, p0, pp1, pp2 Pointers to primitives in cell and neighbors
 \param[in] ix,iy,iz Cell indices
*/
//**********************************************************************
ldouble
reduce_minmod_theta(ldouble *pm2,ldouble *pm1,ldouble *p0,ldouble *pp1,ldouble *pp2,int ix,int iy,int iz)
{
  ldouble theta=MINMOD_THETA;

#if (REDUCEMINMODTHETA==1) //reduce near the axis
  int gix,giy,giz;
  giy=iy+TOJ;
  int limit = NCCORRECTPOLAR + 1;
  if(giy<limit || giy> (TNY-limit-1))
    theta=1.0;
#endif

#if (REDUCEMINMODTHETA==2) //reduce near the inner boundary
  int gix,giy,giz;
  gix=ix+TOI;
  int limit = NCELLSREDUCEMINMODTHETA;
  if(gix<limit)
    theta=1.0;
#endif
  
  return theta;
}
 

//**********************************************************************
/*! \fn int avg2point(ldouble *um2,ldouble *um1,ldouble *u0,ldouble *up1,ldouble *up2,ldouble *ul,ldouble *ur,ldouble dxm2,ldouble dxm1,ldouble dx0,ldouble dxp1,ldouble dxp2,int param,ldouble theta)
 \brief Interpolates primitives to the left and right walls of current cell i
 
 @param[in] um2, um1, u0, up1, up2 values of primitive in the i-2, i-1, i, i+1, i+2 cells
 @param[out] ul, ur interpolated primitives at the left and right walls of cell i
 @param[in] dxm2, dxm1, dx0, dxp1, dxp2 sizes of the five cells
 @param[in] param reconstrpar
 @param[in] minmod_theta MINMOD_THETA
 
 Several interpolation schemes are available.
 
 INT_ORDER=0: basic donor cell\n
 INT_ORDER=1: Minmod (FLUXLIMITER=0), Monotonized Central (FLUXLIMITER=1), Superbee (FLUXLIMITER=4)\n
 INT_ORDER=2: Piecewise Parabolic Method (PPM)\n
 
 */
//**********************************************************************
int
avg2point(ldouble *um2,ldouble *um1,ldouble *u0,ldouble *up1,ldouble *up2,
	  ldouble *ul,ldouble *ur,
	  ldouble dxm2,ldouble dxm1,ldouble dx0,ldouble dxp1,ldouble dxp2,
	  int param,ldouble theta)
{
  ldouble r0[NV],rm1[NV],rp1[NV];

  if(param!=0) //overrule the standard reconstruction
  {
    int i;
    if(param==1) //DONOR CELL
    {
      for(i=0;i<NV;i++)
      {
        ur[i]=u0[i];
        ul[i]=u0[i];
      }
    }
  }  // if(param!=0)
  
  else if(INT_ORDER==0) //DONOR CELL
  {
    int i;
    
    for(i=0;i<NV;i++)
    {
      ur[i]=u0[i];
      ul[i]=u0[i];
    }
  }  // else if(INT_ORDER==0)
  
  else if(INT_ORDER==1) //linear
  {
    ldouble diffpar=theta;
    int i;
    
    for(i=0;i<NV;i++)
    {
      // Slope limiter code rewritten by Ramesh. No function-calls, no divisions.
      ldouble slope;
      ldouble deltam = u0[i]-um1[i];
      ldouble deltap = up1[i]-u0[i];
      
      if (deltam * deltap <= 0.)
      {
        // We are at a local maximum or minimum. Use zero slope (i.e., donor cell)
        
        ur[i] = u0[i];
        ul[i] = u0[i];
      }
      else
      {
        if (deltam > 0.)
        {
          // Slopes are positive. Choose limiter appropriately
          if (FLUXLIMITER == 0) // MinMod
          {
            slope = my_min(my_min(theta*deltam, 0.5*(deltam+deltap)), theta*deltap); // theta=1 is MinMod, theta=2 is MC
          }
          else if (FLUXLIMITER == 1) // MC
          {
            slope = my_min(my_min(2*deltam, 0.5*(deltam+deltap)), 2*deltap);
          }
          else if (FLUXLIMITER == 2) // Osher -- discouraged since it is not symmetric
          {
            printf("Error: Osher slope limiter is discouraged since it is not symmetric\n");
            exit(-2);
          }
          else if (FLUXLIMITER == 3) // Koren -- discouraged since it is not symmetric
          {
            printf("Error: Koren slope limiter is discouraged since it is not symmetric\n");
            exit(-3);
          }
          else // Superbee
          {
            slope = my_max(my_min(2*deltam, deltap), my_min(deltam, 2*deltap));
          }
        }
        else
        {
          // Slopes are negative. Choose limiter appropriately
          
          if (FLUXLIMITER == 0) // MinMod
          {
            slope = my_max(my_max(theta*deltam, 0.5*(deltam+deltap)), theta*deltap); // theta=1 is MinMod, theta=2 is MC
          }
          else if (FLUXLIMITER == 1) // MC
          {
            slope = my_max(my_max(2*deltam, 0.5*(deltam+deltap)), 2*deltap);
          }
          else if (FLUXLIMITER == 2) // Osher -- discouraged since it is not symmetric
          {
            printf("Error: Osher slope limiter is discouraged since it is not symmetric\n");
            exit(-2);
          }
          else if (FLUXLIMITER == 3) // Koren -- discouraged since it is not symmetric
          {
            printf("Error: Koren slope limiter is discouraged since it is not symmetric\n");
            exit(-3);
          }
          else // Superbee
          {
            slope = my_min(my_max(2*deltam, deltap), my_max(deltam, 2*deltap));
          }
        }
        
        ur[i] = u0[i] + 0.5*slope;
        ul[i] = u0[i] - 0.5*slope;
      }
      
      if(isnan(ur[i]) || isnan (ul[i])) printf("%d %e %e %e %e %e\n",i,um2[i],um1[i],u0[i],up1[i],up2[i]);
      //u0 remains intact - in linear reconstruction cell averaged equals cell centered
    } // for(i=0;i<NV;i++)
  }  // else if(INT_ORDER==1)

  else if(INT_ORDER==2) //parabolic PPM
  {
    //The following is based on Colella & Woodward (J. Comp. Phys. 54, 174, 1984).
    //It uses five points: m2, m1, 0, p1, p2.
    //The code has been checked and verified by Ramesh: July 14, 2017
    
    // Define various quantities that apear in the formula
    ldouble dxp2_plus_dxp1 = dxp2 + dxp1;
    ldouble dxp2_plus_dxp1_inv = 1. / dxp2_plus_dxp1;
    ldouble dxp1_plus_dx0 = dxp1 + dx0;
    ldouble dxp1_plus_dx0_inv = 1. / dxp1_plus_dx0;
    ldouble dx0_plus_dxm1 = dx0 + dxm1;
    ldouble dx0_plus_dxm1_inv = 1. / dx0_plus_dxm1;
    ldouble dxm1_plus_dxm2 = dxm1 + dxm2;
    ldouble dxm1_plus_dxm2_inv = 1. / dxm1_plus_dxm2;
    
    ldouble dxm1_plus_dx0_plus_dxp1_inv = 1. / (dxm1+dx0+dxp1);
    ldouble dx0_plus_dxp1_plus_dxp2_inv = 1. / (dx0+dxp1+dxp2);
    ldouble dxm2_plus_dxm1_plus_dx0_inv = 1. / (dxm2+dxm1+dx0);
    
    ldouble dx0_plus_twodxm1 = dx0 + 2. * dxm1;
    ldouble dxp1_plus_twodx0 = dxp1 + 2. * dx0;
    ldouble dxm1_plus_twodxm2 = dxm1 + 2. * dxm2;
    ldouble twodxp1_plus_dx0 = 2. * dxp1 + dx0;
    ldouble twodxp2_plus_dxp1 = 2. * dxp2 + dxp1;
    ldouble twodx0_plus_dxm1 = 2. * dx0 + dxm1;
    
    //ldouble l,r;
    int iv;
    ldouble dri[NV],drim1[NV],drip1[NV];
    
    for(iv=0;iv<NV;iv++)
    {
      // dri, drip1, drim1 are the slopes delta a_j, delta a_{j+1}, delta a_{j-1} in eq (1.7) of C&W
      dri[iv] = dx0 * dxm1_plus_dx0_plus_dxp1_inv *
      (dx0_plus_twodxm1 * dxp1_plus_dx0_inv * (up1[iv]-u0[iv]) +
       twodxp1_plus_dx0 * dx0_plus_dxm1_inv * (u0[iv]-um1[iv]));
      drip1[iv] = dxp1 * dx0_plus_dxp1_plus_dxp2_inv *
      (dxp1_plus_twodx0 * dxp2_plus_dxp1_inv * (up2[iv]-up1[iv]) +
       twodxp2_plus_dxp1 * dxp1_plus_dx0_inv * (up1[iv]-u0[iv]));
      drim1[iv] = dxm1 * dxm2_plus_dxm1_plus_dx0_inv *
      (dxm1_plus_twodxm2 * dx0_plus_dxm1_inv * (u0[iv]-um1[iv]) +
       twodx0_plus_dxm1 * dxm1_plus_dxm2_inv * (um1[iv]-um2[iv]));
      
      // Limit the slopes to be monotonic. This is eq (1.8) in C&W. (Note a minor typo in C&W: one of their _{j-1} should be _{j+1})
      if( (up1[iv]-u0[iv]) * (u0[iv]-um1[iv]) > 0.)
      {
        dri[iv] = my_min(fabs(dri[iv]), my_min(2. * fabs(u0[iv] - um1[iv]), 2. * fabs(u0[iv] - up1[iv]))) * my_sign(dri[iv]);
      }
      else
      {
        dri[iv]=0.;
      }
        
      if( (up2[iv]-up1[iv]) * (up1[iv]-u0[iv]) > 0.)
      {
        drip1[iv] = my_min(fabs(drip1[iv]), my_min(2. * fabs(up1[iv] - u0[iv]), 2. * fabs(up1[iv] - up2[iv]))) * my_sign(drip1[iv]);
      }
      else
      {
        drip1[iv]=0.;
      }
      
      if( (u0[iv]-um1[iv]) * (um1[iv]-um2[iv]) > 0.)
      {
        drim1[iv] = my_min(fabs(drim1[iv]), my_min(2. * fabs(um1[iv] - um2[iv]), 2. * fabs(um1[iv] - u0[iv]))) * my_sign(drim1[iv]);
      }
      else
      {
        drim1[iv]=0.;
      }
    }
    
    // Work on the right face of cell j
    ldouble Z1, Z2, DX_inv;
    Z1 = dx0_plus_dxm1 / dxp1_plus_twodx0;
    Z2 = dxp2_plus_dxp1 / twodxp1_plus_dx0;
    DX_inv = 1. / (dxm1+dx0+dxp1+dxp2);
    
    for(iv=0;iv<NV;iv++)
    {
      // This is a_{j+1/2) in eq (1.6) of Colella & Woodward
      
      ur[iv] = u0[iv] + dx0 * dxp1_plus_dx0_inv * (up1[iv]-u0[iv]) +
            DX_inv * ((2.*dxp1*dx0) * dxp1_plus_dx0_inv * (Z1-Z2) * (up1[iv]-u0[iv]) -
            dx0 * Z1 * drip1[iv] + dxp1 * Z2 * dri[iv]);
    }
    
    // Next work on the left face of cell j
    Z1 = dxm1_plus_dxm2 / dx0_plus_twodxm1;
    Z2 = dxp1_plus_dx0 / twodx0_plus_dxm1;
    DX_inv = 1. / (dxm2+dxm1+dx0+dxp1);
    
    for(iv=0;iv<NV;iv++)
    {
      // This is a_{j-1/2} in eq (1.6) of Colella & Woodward
      ul[iv] = um1[iv] + dxm1 * dx0_plus_dxm1_inv * (u0[iv]-um1[iv]) +
            DX_inv * ((2.*dx0*dxm1) * dx0_plus_dxm1_inv * (Z1-Z2) * (u0[iv]-um1[iv]) -
            dxm1 * Z1 * dri[iv] + dx0 * Z2 * drim1[iv]);
    }
    
    // Make sure that the parabola remains monotonic. The following is equivalent to eq (1.10) in C&W, though it looks different  
    for(iv=0;iv<NV;iv++)
    {
      if((ur[iv]-u0[iv])*(u0[iv]-ul[iv])<=0.)
      {
        ul[iv] = u0[iv];
        ur[iv] = u0[iv];
      }
      if((ur[iv] - ul[iv]) * (ul[iv] - (3. * u0[iv] - 2. * ur[iv])) < 0.)
      {
        ul[iv] = 3. * u0[iv] - 2.*ur[iv];
      }
      if((ur[iv] - ul[iv]) * ((3. * u0[iv] - 2. * ul[iv]) - ur[iv]) < 0.)
      {
        ur[iv] = 3. * u0[iv] - 2. * ul[iv];
      }
      
      //pass up reconstructed value at center - only if reconstructing average -> center
      //To check consistency, make sure that u0[iv] in the LHS below is equal to u0[iv]
      //u0[iv]=ul[iv]+.5*(ur[iv]-ul[iv] + 6.*(u0[iv]-.5*(ul[iv]+ur[iv]))*(1.-.5));
    }
  }  // else if(INT_ORDER==2)

  return 0;
}

//**********************************************************************
/*! \fn int calc_wavespeeds()
 \brief Calculates and saves wavespeeds for cells within the domain plus ghost cells
 
 Wavespeeds are in code units for the basic uniformly spaced coordinates x used in KORAL\n
 For GRMHD problems, assumes that the wave dispersion relation in the fluid frame is isotropic, 
with wave speed given by \f$v^2 = c_s^2 + v_A^2 - c_s^2 v_A^2\f$. 
 This is then transformed to the code coordinate frame.
 
 Calculates left, right and maximum wavespeeds for hydro and radiation for each cell (ix, iy, iz). 
 These velocities are saved in array aaa.\n
 The values in aaa are then transferred to global arrays:\n
 -Hydro left/right: ahdxl, ahdxr, ahdyl, ahdyr, ahdzl, ahdzr\n
 -Hydro maximum absolute: ahdx, ahdy, ahdz\n
 -Radiation left/right: aradxl, aradxr, aradyl, aradyr, aradzl, aradzr\n
 -Radiation maximum absolute: aradx, arady, aradz\n
 */
//**********************************************************************
int
calc_wavespeeds()
{
  
  int ii;
#pragma omp parallel for
  for(ii=0;ii<Nloop_1;ii++) //domain plus lim (=1 usually) ghost cells
  {
    int ix,iy,iz;
    ix=loop_1[ii][0];
    iy=loop_1[ii][1];
    iz=loop_1[ii][2];
    
#ifndef MPI4CORNERS
    if(if_outsidegc(ix,iy,iz)==1) continue; //avoid corners
#endif
    
    ldouble aaa[24];
    
    if(!is_cell_active(ix,iy,iz))
      continue;
    calc_wavespeeds_lr(ix,iy,iz,aaa);
    
    save_wavespeeds(ix,iy,iz,aaa);
  }
  
  return 0;
}


//**********************************************************************
/*! \fn int save_wavespeeds(int ix, int iy, int iz, ldouble *aaa)
 \brief saves characteristic wavespeeds from aaa[] to the global arrays. Also compute timesteps. 

 \param[in] ix,iy,iz cell indices
 \param[in] aaa array with cell wavespeeds
*/
//**********************************************************************
int
save_wavespeeds(int ix,int iy,int iz, ldouble *aaa)
{
  ldouble aaaxhd,aaaxrad,aaayhd,aaayrad,aaazhd,aaazrad;
  ldouble aaaxrad2,aaayrad2,aaazrad2;

  //hydro wavespeeds
  set_u_scalar(ahdxl,ix,iy,iz,aaa[0]);
  set_u_scalar(ahdxr,ix,iy,iz,aaa[1]);
  set_u_scalar(ahdyl,ix,iy,iz,aaa[2]);
  set_u_scalar(ahdyr,ix,iy,iz,aaa[3]);
  set_u_scalar(ahdzl,ix,iy,iz,aaa[4]);
  set_u_scalar(ahdzr,ix,iy,iz,aaa[5]);
  
  aaaxhd=my_max(fabs(aaa[0]),fabs(aaa[1]));
  aaayhd=my_max(fabs(aaa[2]),fabs(aaa[3]));
  aaazhd=my_max(fabs(aaa[4]),fabs(aaa[5]));

  set_u_scalar(ahdx,ix,iy,iz,aaaxhd);
  set_u_scalar(ahdy,ix,iy,iz,aaayhd);
  set_u_scalar(ahdz,ix,iy,iz,aaazhd);

  //original radiative wavespeeds in slots aaa[6::11] are used later
  //speeds in slots[12::] are limited by the optical depth
  //used to calculate the fluxes
  set_u_scalar(aradxl,ix,iy,iz,aaa[6+6]);
  set_u_scalar(aradxr,ix,iy,iz,aaa[6+7]);
  set_u_scalar(aradyl,ix,iy,iz,aaa[6+8]);
  set_u_scalar(aradyr,ix,iy,iz,aaa[6+9]);
  set_u_scalar(aradzl,ix,iy,iz,aaa[6+10]);
  set_u_scalar(aradzr,ix,iy,iz,aaa[6+11]);

  aaaxrad=my_max(fabs(aaa[6+6]),fabs(aaa[6+7]));
  aaayrad=my_max(fabs(aaa[6+8]),fabs(aaa[6+9]));
  aaazrad=my_max(fabs(aaa[6+10]),fabs(aaa[6+11]));
  
  set_u_scalar(aradx,ix,iy,iz,aaaxrad);
  set_u_scalar(arady,ix,iy,iz,aaayrad);
  set_u_scalar(aradz,ix,iy,iz,aaazrad);

  //searching for the maximal unlimited wavespeed for setting the timestep
  // here radiative wavespeeds not limited by the optical depth
  aaaxrad=my_max(fabs(aaa[6]),fabs(aaa[7]));
  aaayrad=my_max(fabs(aaa[8]),fabs(aaa[9]));
  aaazrad=my_max(fabs(aaa[10]),fabs(aaa[11]));

  ldouble dx=get_size_x(ix,0);
  ldouble dy=get_size_x(iy,1);
  ldouble dz=get_size_x(iz,2);

  ldouble wsx,wsy,wsz;
  wsx=aaaxhd;
  wsy=aaayhd;
  wsz=aaazhd;

  #ifdef RADIATION
  #ifndef SKIPRADEVOLUTION
  wsx=my_max(aaaxhd,aaaxrad);
  wsy=my_max(aaayhd,aaayrad);
  wsz=my_max(aaazhd,aaazrad);
  #endif
  #endif

  //determine the time step
  //only domain cells
  if(if_indomain(ix,iy,iz)==1) 
    {      
      ldouble tstepden,ws_ph;
   
      if(NZ>1 && NY>1)
	tstepden=wsx/dx + wsy/dy + wsz/dz;
      else if(NZ==1 && NY>1)
	tstepden=wsx/dx + wsy/dy;
      else if(NY==1 && NZ>1)
	tstepden=wsx/dx + wsz/dz;
      else
	tstepden=wsx/dx;   

      tstepden/=TSTEPLIM;

      set_u_scalar(cell_tstepden,ix,iy,iz,tstepden);

      //global variables for maximum/minimum (inverse) cell timestep
      ////#pragma omp critical
      if(tstepden>tstepdenmax) tstepdenmax=tstepden;  
      if(tstepden<tstepdenmin) tstepdenmin=tstepden;  
    }
  return 0;
}

//**********************************************************************
/*! \fn int save_timesteps
 \brief saves individual timesteps to be used inside time stepping
 Does not work with MPI
*/
//**********************************************************************
int
save_timesteps()
{

  int ii;
  ldouble dtminloc = BIG;

  //#pragma omp parallel for reduction(min:dtminloc)  // "min" option in reduction is not available in older OpenMP implementations
  for(ii = 0; ii < Nloop_0; ii++)
  {
    int ix, iy, iz;
    ix = loop_0[ii][0];
    iy = loop_0[ii][1];
    iz = loop_0[ii][2];
    
    ldouble cell_dt_local = 1. / get_u_scalar(cell_tstepden, ix, iy, iz);

    // timestep is shortest of current cell and neighbors.
    // ANDREW: is this doing the right thing on the x boundaries?
#ifdef SHORTERTIMESTEP
    ldouble dtm1, dtp1, dt;
    if(ix > 0)
      dtm1 = 1. / get_u_scalar(cell_tstepden, ix-1, iy, iz);
    else
      dtm1 = BIG;
    
    if(ix < NX)
      dtp1 = 1. / get_u_scalar(cell_tstepden, ix+1, iy, iz);
    else
      dtp1 = BIG;
    
    dt = 1. / get_u_scalar(cell_tstepden, ix, iy, iz);
    
    cell_dt_local = my_min(my_min(dtm1, dtp1), dt);
#endif
    
    set_u_scalar(cell_dt,ix,iy,iz,cell_dt_local);
    
    //find the shortest

    //global variables for shortest cell timestep
    if(cell_dt_local < dtminloc)
      dtminloc = cell_dt_local;
  }
  
  return 0;
}


//**********************************************************************
/*! \fn int calc_u2p(int type, int setflags)
 \brief Calculates all primitives from global u
 \param[in] type, not currently used
 \param[in] setflags, should always=1 to set flags for cell fixups
*/
//**********************************************************************
int
calc_u2p(int type, int setflags)
{
  int ii;

  //timer start
  struct timespec temp_clock;
  my_clock_gettime(&temp_clock);
  start_u2ptime=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
  
  //calculates the primitives
#pragma omp parallel for schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain only
  {
    int ix,iy,iz;
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2];
    
    //skip if cell is passive
    if(!is_cell_active(ix,iy,iz))
      continue;
    
    calc_primitives(ix,iy,iz,type,setflags);
  }
  
  //fixup here hd and rad after inversions
  cell_fixup(FIXUP_U2PMHD);
#ifdef RADIATION
  cell_fixup(FIXUP_U2PRAD);
#endif

  //re-set boundary conditions
  set_bc(global_time,0);
  
  //timer stop
  my_clock_gettime(&temp_clock);
  end_u2ptime=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
  
  return 0;
} 


//**********************************************************************
/*! \fn int do_correct()
 \brief Corrects conserved quantities on polar axis and neutron star surface
*/
//**********************************************************************
int
do_correct()
{
  int ix,iy,iz,ii;
  
  //correct the polar axis 
#ifdef CORRECT_POLARAXIS
  correct_polaraxis();
#endif
#ifdef CORRECT_POLARAXIS_3D
  ///correct_polaraxis(); //first fill at fixed phi including magnetic field
  correct_polaraxis_3d(); //then overwrite velocities with phi averages
#endif
#ifdef SMOOTH_POLARAXIS
  smooth_polaraxis();
#endif

  //correct the neutron star surface
#ifdef CORRECT_NSSURFACE
  correct_nssurface();
#endif

  return 0;
}

//**********************************************************************
/*! \fn int op_explicit(ldouble t, ldouble dtin)
 \brief Explicit evolution from time t to t+dtin
 
 \param[in] t, dtin initial time and time step
 
 Saves initial primitives and conserveds in ppreexplict, upreexplicit\n
 Interpolates primitives to cell faces using the required interpolation scheme, 
 and calculates left- and right-biased conserveds and fluxes\n
 Calculates wave speeds\n
 Computes combined fluxes at the faces using the selected approximate Riemann solver\n
*/
//**********************************************************************

int
op_explicit(ldouble t, ldouble dtin) 
{
  int ix,iy,iz,iv,ii;
  ldouble dt;
  
  int perform_sweep = 1; //Brandon - Added for special conditions where a sweep should not be performed under SPECIAL_BC_CHECK
  #ifdef SPECIAL_BC_CHECK
  int giix,giiy,giiz;
  #endif

  // Save conserveds and primitives over domain + ghost (no corners)
  copyi_u(1.,u,upreexplicit); //conserved quantities before explicit update
  copyi_u(1.,p,ppreexplicit); //primitive quantities before explicit update

  // calculates H/R and velocity averages
  calc_avgs_throughout(); 

  // calculates wavespeeds over the domain and ghost cells
  calc_wavespeeds();

#ifndef SKIPEVOLUTION

  //**********************************************************************
  // Next interpolate to the cell walls and calculate left and right-biased fluxes
#pragma omp parallel for private(ii,ix,iy,iz,iv)  schedule (static)
  for(ii=0;ii<Nloop_1;ii++) //domain plus lim (=1 usually) ghost cells
  {
      ix=loop_1[ii][0];
      iy=loop_1[ii][1];
      iz=loop_1[ii][2];

      #ifdef SPECIAL_BC_CHECK
      giix = ix + TOI;
      giiy = iy + TOJ;
      giiz = iz + TOK;
      #endif

      #ifndef MPI4CORNERS
      if(if_outsidegc(ix,iy,iz)==1) continue; //avoid corners
      #endif
     
      #ifdef SPECIAL_BC_CHECK
      #include PR_BC_SPECIAL_LOOP
      if(ret_val == 4) 
      {
        continue; //Exclude 'corners' in stream region
      }
      #endif

      //create arrays for interpolating conserved quantities
      struct geometry geom;
      //ldouble x0[3],x0l[3],x0r[3],xm1[3],xp1[3];
      ldouble fd_r0[NV],fd_rm1[NV],fd_rp1[NV];
      ldouble fd_u0[NV],fd_up1[NV],fd_up2[NV],fd_um1[NV],fd_um2[NV];
      ldouble fd_p0[NV],fd_pp1[NV],fd_pp2[NV],fd_pm1[NV],fd_pm2[NV],fd_pm3[NV],fd_pp3[NV];
      ldouble fd_s0[NV],fd_sp1[NV],fd_sp2[NV],fd_sm1[NV],fd_sm2[NV];
      ldouble fd_uLr[NV],fd_uLl[NV],fd_uRl[NV],fd_uRr[NV];
      ldouble fd_pLr[NV],fd_pLl[NV],fd_pRl[NV],fd_pRr[NV];
      ldouble fd_pr[NV],fd_pl[NV];
      ldouble fd_sLr[NV],fd_sLl[NV],fd_sRl[NV],fd_sRr[NV];
      ldouble fd_fstarl[NV],fd_fstarr[NV],fd_dul[3*NV],fd_dur[3*NV],fd_pdiffl[NV],fd_pdiffr[NV];
      ldouble a0[2],am1[2],ap1[2],al,ar,amax,dx;  
      ldouble ffRl[NV],ffRr[NV],ffLl[NV],ffLr[NV];
      ldouble ffl[NV],ffr[NV]; 
      ldouble dx0, dxm2, dxm1, dxp1, dxp2;
      ldouble minmod_theta=MINMOD_THETA;
      int reconstrpar;
      int i,dol,dor;



      //**********************************************************************
      // x 'sweep'
      //**********************************************************************

      perform_sweep = 1;

#ifdef MPI4CORNERS
      if(NX>1 && iy>=-1 && iy<NY+1 && iz>=-1 && iz<NZ+1 && perform_sweep == 1) //needed to calculate face fluxes for flux-CT divB enforcement
#else
      if(NX>1 && iy>=0 && iy<NY && iz>=0 && iz<NZ && perform_sweep == 1)
#endif
      {
		dol=dor=1;
		if(ix<0) dol=0;
		if(ix>=NX) dor=0;

                #ifdef SPECIAL_BC_CHECK //Don't do l/r fluxes when at GC - Brandon
                #ifndef SEARCH_STREAM_BOUNDARY
                if(TNY>1 && TNZ==1)
                {
                  if(giiy >= STREAM_IYT && giiy <= STREAM_IYB)
                  {
                    if(giix == STREAM_IX || giix == (STREAM_IX+1) || giix == (STREAM_IX+2) || giix == (STREAM_IX+3)) dor=0;
                  }
                }
                else if(TNY==1 && TNZ>1)
                {
                  if(giiz >= STREAM_IZT && giiz <= STREAM_IZB)
                  {
                    if(giix == STREAM_IX || giix == (STREAM_IX+1) || giix == (STREAM_IX+2) || giix == (STREAM_IX+3)) dor=0;
                  }
                }
                else if(TNY>1 && TNZ>1)
                {
                  #ifndef STREAM_RING
                  if(giiy >= STREAM_IYT && giiy <= STREAM_IYB && giiz >= STREAM_IZT && giiz <= STREAM_IZB)
                  #else
                  if(giiy >= STREAM_IYT && giiy <= STREAM_IYB)
                  #endif
                  {
                    if(giix == STREAM_IX || giix == (STREAM_IX+1) || giix == (STREAM_IX+2) || giix == (STREAM_IX+3)) dor=0;
                  }
                }
                #endif
                #endif

                // is_cell_active is currently always 1
		// skip flux calculation if not needed
		if((ix==0 && is_cell_active(ix,iy,iz)==0) || (ix>0 && is_cell_active(ix,iy,iz)==0 && is_cell_active(ix-1,iy,iz)==0))
		  dol=0;
		if((ix==NX-1 && is_cell_active(ix,iy,iz)==0) || (ix<NX-1 && is_cell_active(ix,iy,iz)==0 && is_cell_active(ix+1,iy,iz)==0))
		  dor=0;
		
                // x0[0] is x of current cell        
		//x0[0]=get_x(ix,0);

                // xm1[0], xp1[0] are x of left and right cell centers. Are these quantities used anywhere?
		//xm1[0]=get_x(ix-1,0);
	        //xp1[0]=get_x(ix+1,0);
		
                // x0l[0,1,2] are x, y, z of left x-wall, x0r[0,1,2] are x, y, z of right x-wall,
		//x0l[0]=get_xb(ix,0);
		//x0l[1]=xm1[1]=get_x(iy,1); 
		//x0l[2]=xm1[2]=get_x(iz,2);

		//x0r[0]=get_xb(ix+1,0);
		//x0r[1]=xp1[1]=get_x(iy,1);
		//x0r[2]=xp1[2]=get_x(iz,2);
	      
                // dx0, dxm1, dxp1 are x-sizes (wall to wall) of cells ix, ix-1, ix+1, dxm2m, dxp2 are sizes of cells ix-2, ix+2		
		dx0=get_size_x(ix,0);    
		dxm1=get_size_x(ix-1,0);    
		dxp1=get_size_x(ix+1,0);    
	  
		if(INT_ORDER>1)
		{
		  dxm2=get_size_x(ix-2,0);    
		  dxp2=get_size_x(ix+2,0);    
		}
	  
		for(i=0;i<NV;i++)
		{
		  //fd_p0, fd_pp1, fd_pm1 are primitives at current, left and right cells, fd_pm2, fd_pp2 are for next two cells
		  fd_p0[i] =get_u(p,i,ix,iy,iz);
		  fd_pp1[i]=get_u(p,i,ix+1,iy,iz);
		  fd_pm1[i]=get_u(p,i,ix-1,iy,iz);
            
		  if(INT_ORDER>1)
		  {
		    fd_pm2[i]=get_u(p,i,ix-2,iy,iz);
		    fd_pp2[i]=get_u(p,i,ix+2,iy,iz);
		  }
		}

		reconstrpar=0;

#ifdef REDUCEORDERWHENNEEDED
		reconstrpar = reduce_order_check(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif
		minmod_theta=MINMOD_THETA;
#ifdef REDUCEMINMODTHETA  // reduce minmod_theta near axis or inner boundary
		minmod_theta = reduce_minmod_theta(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif

                // Interpolate primitives to the left and right walls of current cell: fd_pl, fd_pr
		avg2point(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2,reconstrpar,minmod_theta);   

		// if(ix>0)
		if(dol) //no need to calculate at left face of first GC if dol=0
		{
                  // Left wall of current cell: compute fluxes and save in array ffl[NV]
 		  fill_geometry_face(ix,iy,iz,0,&geom);
		  check_floors_mhd(fd_pl,VELPRIM,&geom);
		  f_flux_prime(fd_pl,0,ix,iy,iz,ffl,1);
		}

		// if(ix<NX)
		if(dor) //no need to calculate at right face of first GC if dor=0
		{
                  // Right wall of current cell: compute fluxes and save in array ffr[NV]
		  fill_geometry_face(ix+1,iy,iz,0,&geom);
		  check_floors_mhd(fd_pr,VELPRIM,&geom);
		  f_flux_prime(fd_pr,0,ix+1,iy,iz,ffr,0);
		}

		//save interpolated values to memory
                //Note that l and r of a given cell ix are the left and right wall of that cell
		//whereas L and R of given ix are quantities to the left and right of wall ix
		for(i=0;i<NV;i++)
		{
                  // Save fd_pl in array pbRx (Primitive_R) of wall ix
                  // Save fd_pr in array pbLx (Primitive_L) of wall ix+1
		  set_ubx(pbRx,i,ix,iy,iz,fd_pl[i]);
		  set_ubx(pbLx,i,ix+1,iy,iz,fd_pr[i]);
		    
		  if(dol)
                  // Save ffl in array flRx (F_R) of wall ix
		  set_ubx(flRx,i,ix,iy,iz,ffl[i]);
		  if(dor)
                  // Save ffr in array flLx (F_L) of wall ix+1 
		  set_ubx(flLx,i,ix+1,iy,iz,ffr[i]);
		} 
       }  // if(NX>1 && iy>=0 && iy<NY && iz>=0 && iz<NZ...)
      
      //**********************************************************************
      //y 'sweep'
      //**********************************************************************
  
      perform_sweep = 1;
#ifdef SPECIAL_BC_CHECK
      //if(giix > STREAM_IX)
      //if(giix > STREAM_IX && giix <= (STREAM_IX+1))
      if(giix > STREAM_IX && giix <= (STREAM_IX+3))
      {
#ifdef MPI4CORNERS
        if(TNY>1 && TNZ==1)
        {
          if(giiy > (STREAM_IYT-1) && giiy < (STREAM_IYB+1)) perform_sweep = 0;
        }
        else if(TNY>1 && TNZ>1)
        {
          #ifndef STREAM_RING
          if(giiz > (STREAM_IZT-1) && giiz < (STREAM_IZB-1))
          { 
            if(giiy > (STREAM_IYT-1) && giiy < (STREAM_IYB+1)) perform_sweep = 0;
          }
          #else
          if(giiy > (STREAM_IYT-1) && giiy < (STREAM_IYB+1)) perform_sweep = 0;
          #endif
        }
#else
        if(TNY>1 && TNZ==1)
        {
          if(giiy >= (STREAM_IYT-1) && giiy <= (STREAM_IYB+1)) perform_sweep = 0;
        }
        else if(TNY>1 && TNZ>1)
        {
          #ifndef STREAM_RING
          if(giiz >= (STREAM_IZT-1) && giiz <= (STREAM_IZB-1))
          { 
            if(giiy >= (STREAM_IYT-1) && giiy <= (STREAM_IYB+1)) perform_sweep = 0;
          }
          #else
          if(giiy >= (STREAM_IYT-1) && giiy <= (STREAM_IYB+1)) perform_sweep = 0;
          #endif
        }
#endif
      }
#endif

#ifdef MPI4CORNERS
      if(NY>1 && ix>=-1 && ix<NX+1 && iz>=-1 && iz<NZ+1 && perform_sweep == 1)
#else
      if(NY>1 && ix>=0 && ix<NX && iz>=0 && iz<NZ && perform_sweep == 1)
#endif
      {
		dol=dor=1;
		if(iy<0) dol=0;
		if(iy>=NY) dor=0;

                // is_cell_active is currently always 1
                // skip flux calculation if not needed
		if((iy==0 && is_cell_active(ix,iy,iz)==0) || (iy>0 && is_cell_active(ix,iy,iz)==0 && is_cell_active(ix,iy-1,iz)==0))
		  dol=0;
             
		if((iy==NY-1 && is_cell_active(ix,iy,iz)==0) || (iy<NY-1 && is_cell_active(ix,iy,iz)==0 && is_cell_active(ix,iy+1,iz)==0))
		  dor=0;
        
                // x0[1] is y of current cell        
	        //x0[1]=get_x(iy,1);

                // xm1[1], xp1[1] are y of left and right cell centers. Are these quantities used anywhere?
		//xm1[1]=get_x(iy-1,1);
                //xp1[1]=get_x(iy+1,1);
		
                // x0l[0,1,2] are x, y, z of left y-wall, x0r[0,1,2] are x, y, z of right y-wall,		
	 	//x0l[1]=get_xb(iy,1);
		//x0l[0]=xm1[0]=get_x(ix,0); 
		//x0l[2]=xm1[2]=get_x(iz,2);

		//x0r[1]=get_xb(iy+1,1);
		//x0r[0]=xp1[0]=get_x(ix,0);
		//x0r[2]=xp1[2]=get_x(iz,2);

                // dx0, dxm1, dxp1 are y-sizes (wall to wall) of cells iy, iy-1, iy+1, dxm2m, dxp2 are sizes of cells iy-2, iy+2
		dx0=get_size_x(iy,1);    
		dxm1=get_size_x(iy-1,1);    
		dxp1=get_size_x(iy+1,1);    
	
		if(INT_ORDER>1)
		{
		  dxm2=get_size_x(iy-2,1);  
		  dxp2=get_size_x(iy+2,1);
		}    
		  
		for(i=0;i<NV;i++)
		{
                  //fd_p0, fd_pp1, fd_pm1 are primitives at current, left and right cells, fd_pm2, fd_pp2 are for next two cells
		  fd_p0[i]=get_u(p,i,ix,iy,iz);
		  fd_pp1[i]=get_u(p,i,ix,iy+1,iz);
		  fd_pm1[i]=get_u(p,i,ix,iy-1,iz);
		  if(INT_ORDER>1)
		  {
		    fd_pm2[i]=get_u(p,i,ix,iy-2,iz);
		    fd_pp2[i]=get_u(p,i,ix,iy+2,iz);
		  }
		}
	  
		reconstrpar=0;
#ifdef REDUCEORDERWHENNEEDED
		reconstrpar = reduce_order_check(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif

		minmod_theta=MINMOD_THETA;
#ifdef REDUCEMINMODTHETA  // reduce minmod_theta near axis or inner boundary
		minmod_theta = reduce_minmod_theta(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif

                // Interpolate primitives to the left and right walls of current cell: fd_pl, fd_pr
                avg2point(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2,reconstrpar,minmod_theta);

		//iy>0
		if(dol) //no need to calculate at left face of first GC if dol=0
		{
                  // Left wall of current cell: compute fluxes and save in array ffl[NV]
		  fill_geometry_face(ix,iy,iz,1,&geom);
		  check_floors_mhd(fd_pl,VELPRIM,&geom);
		  f_flux_prime(fd_pl,1,ix,iy,iz,ffl,1);
		}

		//iy<NY
		if(dor) //no need to calculate at right face of first GC if dor=0
		{
                  // Right wall of current cell: compute fluxes and save in array ffr[NV]
		  fill_geometry_face(ix,iy+1,iz,1,&geom);
		  check_floors_mhd(fd_pr,VELPRIM,&geom);
		  f_flux_prime(fd_pr,1,ix,iy+1,iz,ffr,0);   	          
		}

                //save interpolated values to memory
                //Note that l and r of a given cell iy are the left and right wall of that cell,
		//whereas L and R of given iy are quantities to the left and right of wall iy
		for(i=0;i<NV;i++)
		{
                  // Save fd_pl in array pbRy (Primitive_R) of wall iy
                  // Save fd_pr in array pbLy (Primitive_L) of wall iy+1
		  set_uby(pbRy,i,ix,iy,iz,fd_pl[i]);
		  set_uby(pbLy,i,ix,iy+1,iz,fd_pr[i]);

		  if(dol)
                  // Save ffl in array flRy (F_R) of wall iy
		  set_uby(flRy,i,ix,iy,iz,ffl[i]);
		  if(dor)
                  // Save ffr in array flLy (F_R) of wall iy+1
		  set_uby(flLy,i,ix,iy+1,iz,ffr[i]);
		}
      }  // if(NY>1 && ix>=0 && ix<NX && iz>=0 && iz<NZ...)

      //**********************************************************************
      //z 'sweep'
      //**********************************************************************
	      
      perform_sweep = 1;

#ifdef SPECIAL_BC_CHECK
      if(giix > STREAM_IX && giix <= (STREAM_IX+3))
      {
#ifdef MPI4CORNERS
        if(TNY==1 && TNZ>1)
        {
          if(giiz > (STREAM_IZT-1) && giiz < (STREAM_IZB+1)) perform_sweep = 0;
        }
        else if(TNY>1 && TNZ>1)
        {
          #ifndef STREAM_RING
          if(giiz > (STREAM_IZT-1) && giiz < (STREAM_IZB-1))
          { 
            if(giiy > (STREAM_IYT-1) && giiy < (STREAM_IYB+1)) perform_sweep = 0;
          }
          #else
          if(giiy > (STREAM_IYT-1) && giiy < (STREAM_IYB+1)) perform_sweep = 0;
          #endif
        }
#else
        if(TNY==1 && TNZ>1)
        {
          if(giiz >= (STREAM_IZT-1) && giiz <= (STREAM_IZB+1)) perform_sweep = 0;
        }
        else if(TNY>1 && TNZ>1)
        {
          #ifndef STREAM_RING
          if(giiz >= (STREAM_IZT-1) && giiz <= (STREAM_IZB-1))
          { 
            if(giiy >= (STREAM_IYT-1) && giiy <= (STREAM_IYB+1)) perform_sweep = 0;
          }
          #else
          if(giiy >= (STREAM_IYT-1) && giiy <= (STREAM_IYB+1)) perform_sweep = 0;
          #endif
        }
#endif
      }
#endif

#ifdef MPI4CORNERS
      if(NZ>1 && ix>=-1 && ix<NX+1 && iy>=-1 && iy<NY+1 && perform_sweep == 1)
#else
      if(NZ>1 && ix>=0 && ix<NX && iy>=0 && iy<NY && perform_sweep == 1)
#endif
      {
	        dol=dor=1;
                if(iz<0) dol=0;
                if(iz>=NZ) dor=0;
         
                // is_cell_active is currently always 1
                // skip flux calculation if not needed
                if((iz==0 && is_cell_active(ix,iy,iz)==0) || (iz>0 && is_cell_active(ix,iy,iz)==0 && is_cell_active(ix,iy,iz-1)==0))
                  dol=0;
                if((iz==NZ-1 && is_cell_active(ix,iy,iz)==0) || (iz<NZ-1 && is_cell_active(ix,iy,iz)==0 && is_cell_active(ix,iy,iz+1)==0))
                  dor=0;
         
                // x0[2] is z of current cell        
	        //x0[2]=get_x(iz,2);

                // xm1[2], xp1[2] are z of left and right cell centers. Are these quantities used anywhere?
		//xm1[2]=get_x(iz-1,2);
                //xp1[2]=get_x(iz+1,2);

                // x0l[0,1,2] are x, y, z of left z-wall, x0r[0,1,2] are x, y, z of right z-wall,		
                //x0l[2]=get_xb(iz,2);
                //x0l[0]=xm1[0]=get_x(ix,0);
                //x0l[1]=xm1[1]=get_x(iy,1);
         
                //x0r[2]=get_xb(iz+1,2);
                //x0r[0]=xp1[0]=get_x(ix,0);
                //x0r[1]=xp1[1]=get_x(iy,1);
         
                // dx0, dxm1, dxp1 are z-sizes (wall to wall) of cells iz, iz-1, iz+1, dxm2m, dxp2 are sizes of cells iz-2, iz+2
                dx0=get_size_x(iz,2);
                dxm1=get_size_x(iz-1,2);
                dxp1=get_size_x(iz+1,2);
         
                if(INT_ORDER>1)
                {
                  dxm2=get_size_x(iz-2,2);
                  dxp2=get_size_x(iz+2,2);
                }
         
                for(i=0;i<NV;i++)
                {
                  //fd_p0, fd_pp1, fd_pm1 are primitives at current, left and right cells, fd_pm2, fd_pp2 are for next two cells
                  fd_p0[i]=get_u(p,i,ix,iy,iz);
                  fd_pp1[i]=get_u(p,i,ix,iy,iz+1);
                  fd_pm1[i]=get_u(p,i,ix,iy,iz-1);
           
                  if(INT_ORDER>1)
                  {
                    fd_pm2[i]=get_u(p,i,ix,iy,iz-2);
                    fd_pp2[i]=get_u(p,i,ix,iy,iz+2);
                  }
                }
         
                reconstrpar=0;
#ifdef REDUCEORDERWHENNEEDED
                reconstrpar = reduce_order_check(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif
         
                minmod_theta=MINMOD_THETA;
#ifdef REDUCEMINMODTHETA  // reduce minmod_theta near axis or inner boundary
                minmod_theta = reduce_minmod_theta(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif
         
                // Interpolate primitives to the left and right walls of current cell: fd_pl, fd_pr
                avg2point(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2,reconstrpar,minmod_theta);

		//iz>0
                if(dol) //no need to calculate at left face of first GC if dol=0
                {
                  // Left wall of current cell: compute fluxes and save in array ffl[NV]
                  fill_geometry_face(ix,iy,iz,2,&geom);
                  check_floors_mhd(fd_pl,VELPRIM,&geom);
                  f_flux_prime(fd_pl,2,ix,iy,iz,ffl,0);
                }

		//iz<NZ
                if(dor) //no need to calculate at right face of first GC if dor=0
                {
                  // Right wall of current cell: compute fluxes and save in array ffr[NV]
                  fill_geometry_face(ix,iy,iz+1,2,&geom);
                  check_floors_mhd(fd_pr,VELPRIM,&geom);
                  f_flux_prime(fd_pr,2,ix,iy,iz+1,ffr,1);
                }
         
                //save interpolated values to memory
                //Note that l and r of a given cell iy are the left and right wall of that cell,
                //whereas L and R of given iy are quantities to the left and right of wall iy
                for(i=0;i<NV;i++)
                {
                  // Save fd_pl in array pbRz (Primitive_R) of wall iz
                  // Save fd_pr in array pbLz (Primitive_L) of wall iz+1
                  set_ubz(pbRz,i,ix,iy,iz,fd_pl[i]);
                  set_ubz(pbLz,i,ix,iy,iz+1,fd_pr[i]);
           
                  if(dol)
                  // Save ffl in array flRz (F_R) of wall iz
                  set_ubz(flRz,i,ix,iy,iz,ffl[i]);
                  if(dor)
                  // Save ffr in array flLz (F_R) of wall iz+1
                  set_ubz(flLz,i,ix,iy,iz+1,ffr[i]);   
                }
      }  // if(NZ>1 && ix>=0 && ix<NX && iy>=0 && iy<NY...)
	     
  }  // for(ii=0;ii<Nloop_1;ii++)
	    
  //**********************************************************************
  // Compute fluxes at the six walls of all cells using the selected approximation of the Riemann problem
  
#pragma omp barrier
#pragma omp parallel for private(ii,iy,iz,ix)  schedule (static)
  for(ii=0;ii<Nloop_1;ii++) //domain plus lim (=1 usually) ghost cells
  {
    ix=loop_1[ii][0];
    iy=loop_1[ii][1];
    iz=loop_1[ii][2];
    f_calc_fluxes_at_faces(ix,iy,iz);
  }

  //**********************************************************************
  // Constrained transport to preserve div.B=0

#ifdef MAGNFIELD
  #pragma omp barrier
  flux_ct();
#endif

  //**********************************************************************
  // Evolve the conserved quantities
  
#pragma omp barrier  
#pragma omp parallel for private(ii,ix,iy,iz,iv) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
  {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      // Source term
      ldouble ms[NV],gs[NV],val,du;

      if(is_cell_active(ix,iy,iz)==0)
      {
        // Source terms applied only for active cells	
	PLOOP(iv) ms[iv]=0.; 
      }
      else
      {
        // Get metric source terms ms[iv] and any other source terms gs[iv], and save combined source terms in ms[iv]
	f_metric_source_term(ix,iy,iz,ms);
	f_general_source_term(ix,iy,iz,gs);

	PLOOP(iv) ms[iv]+=gs[iv];
      }
      
      // Get the cell size in the three directions
      ldouble dx=get_size_x(ix,0);
      ldouble dy=get_size_x(iy,1);
      ldouble dz=get_size_x(iz,2);

      int doxl,doxr,doyl,doyr,dozl,dozr;
      doxl=doxr=doyl=doyr=dozl=dozr=1;

      //timestep
      dt=dtin;  // dtin is an input parameter to op_explicit
      
      //update all conserved according to fluxes and source terms
      for(iv=0;iv<NV;iv++)
      {
	ldouble flxr,flyr,flzr,flxl,flyl,flzl;
	  
        // Get the fluxes on the six faces.
	// Recall that flbx, flby, flbz are the fluxes at the left walls of cell ix, iy, iz.
	// To get the right fluxes, we need flbx(ix+1,iy,iz), etc.
	flxl=get_ub(flbx,iv,ix,iy,iz,0);
	flxr=get_ub(flbx,iv,ix+1,iy,iz,0);
	flyl=get_ub(flby,iv,ix,iy,iz,1);
	flyr=get_ub(flby,iv,ix,iy+1,iz,1);
	flzl=get_ub(flbz,iv,ix,iy,iz,2);
	flzr=get_ub(flbz,iv,ix,iy,iz+1,2);
		  
	// Compute Delta U from the six fluxes
	du = -(flxr-flxl)*dt/dx - (flyr-flyl)*dt/dy - (flzr-flzl)*dt/dz;

	// Compute new conserved by adding Delta U and the source term
	val = get_u(u,iv,ix,iy,iz) + du + ms[iv]*dt;

	// Save the new conserved to memory
#ifdef SKIPHDEVOLUTION
	if(iv>=NVMHD)
#endif
#ifdef RADIATION
#ifdef SKIPRADEVOLUTION
#ifdef EVOLVEPHOTONNUMBER
	if(iv!=EE && iv!=FX && iv!=FY && iv!=FZ && iv!=NF)
#else
	if(iv!=EE && iv!=FX && iv!=FY && iv!=FZ)
#endif
#endif  // SKIPRADEVOLUTION
#endif  // RADIATION
#ifdef SKIPHDBUTENERGY
	if(iv>=NVMHD || iv==UU)
#endif
	  set_u(u,iv,ix,iy,iz,val);	 

      }  // for(iv=0;iv<NV;iv++)

  }  // for(ii=0;ii<Nloop_0;ii++)


   /************************************************************************/
   /********* explicit *** RADIATION COUPLING  *****************************/
   /************************************************************************/
  
#ifdef RADIATION
#ifndef SKIPRADSOURCE
#ifdef EXPLICIT_LAB_RAD_SOURCE
  #pragma omp barrier
  #pragma omp parallel for private(ii,ix,iy,iz,iv,dt) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
  {
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2]; 

    if(is_cell_active(ix,iy,iz)==0) 
      continue;

    //timestep
    dt=dtin;

    // Use primitives from *p, i.e., from the beginning of this timestep.
    // explicit_rad_source_term computes the radiation source term and adds it to the conserveds.
    explicit_rad_source_term(ix,iy,iz,dt);
  }
  
#endif //EXPLICIT_LAB_RAD_SOURCE
#endif //SKIPRADSOURCE
#endif //RADIATION
#endif //SKIPEVOLUTION


  //**********************************************************************  
  // Compute postexplicit primitives and count entropy inversions
  calc_u2p(1,1);
  //**********************************************************************
	    
  // Entropy Mixing
  
#ifdef MIXENTROPIESPROPERLY
  mix_entropies(dt);
#endif

  return GSL_SUCCESS;
}

//**********************************************************************
/*! \fn int op_intermediate(ldouble t, ldouble dtin)
 \brief  applies particle heating after explicit evolution
 \param[in] t, dtin current global time and time step
*/
//**********************************************************************
int
op_intermediate(ldouble t, ldouble dt)
{

   // Apply adiabatic expansion to nonthermal electrons
#ifdef EVOLVEELECTRONS
#ifdef RELELECTRONS
#ifndef SKIPRELELEXPANSION
  relel_adiab(dt);
#endif
#endif
  
  // Apply viscous heating to thermal & nonthermal electrons and ions  
#ifndef HEATELECTRONSATENDRK2   //here? or separate after RK2

  heat_electronions_with_state(dt); 
#endif
#endif
  
  return 0;
}

//**********************************************************************
/*! \fn int apply_dynamo(ldouble t, ldouble dtin)
 \brief  mimics alpha-dynamo in axisymmetric sims involving MRI 
 \param[in] t, dtin current global time and time step
*/
//**********************************************************************
int
apply_dynamo(ldouble t, ldouble dt)
{

#ifdef MIMICDYNAMO
  
  //correlates ghost cells
  mpi_exchangedata();
  calc_avgs_throughout();
  set_bc(t,0);
  
  //mimics dynamo
  mimic_dynamo(dt); 

  //update primitives
  calc_u2p(0,0);

#endif

  return GSL_SUCCESS;
}


//**********************************************************************
/*! \fn int op_implicit(ldouble t, ldouble dtin)
 \brief The radiative implicit source term operator
 
 \param[in] t, dtin current global time and time step
*/
//**********************************************************************
int
op_implicit(ldouble t, ldouble dtin) 
{

  int ii;


  // counter for the average number of iterations in the implicit solver
  for(ii=0;ii<NGLOBALINTSLOT;ii++)
    global_int_slot[ii]=0.;

  /************************************************************************/
  /******** implicit **** RADIATION ***************************************/
  /************************************************************************/

#ifdef RADIATION
#ifndef SKIPRADSOURCE
#ifdef IMPLICIT_LAB_RAD_SOURCE  // this is the default

#pragma omp parallel for private(ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
  {
    int ix,iy,iz;
    ldouble dt;
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2]; 

    if(is_cell_active(ix,iy,iz)==0 || is_cell_corrected_polaraxis(ix,iy,iz)==1 || skip_cell_implicit(ix,iy,iz)==1)
      continue;

    //timestep
    dt=dtin;

    //uses values already in *p as the initial guess
    implicit_lab_rad_source_term(ix,iy,iz,dt);
  } 

  //fixup here after implicit source term 
  cell_fixup(FIXUP_RADIMP);

#endif //IMPLICIT_LAB_RAD_SOURCE
#endif //SKIPRADSOURCE
#endif //RADIATION

  return 0;
}


//************************************************************************
/*! \fn ldouble f_calc_fluxes_at_faces(int ix,int iy,int iz)
 \brief Calculates fluxes at cell faces
 
 \param[in] ix, iy, iz Indices of the three cell faces under consideration
 
 Uses the selected approximate Riemann solver to estimate the fluxes on the six faces of the cell (ix, iy, iz)\n
 The fluxes are saved in global arrays flbx, flby, flbz\n
 Note that the fluxes calculated here correspond to the three left walls of cell ix, iy, iz
 
*/
//************************************************************************

ldouble f_calc_fluxes_at_faces(int ix,int iy,int iz)
{
  int i;
  struct geometry geom;
   
  ldouble a0[2],am1[2],ap1[2],ag,al,ar,amax,cmin,cmax,csLl[2],csLr[2],csRl[2],csRr[2];
  ldouble am1l[2],am1r[2],ap1l[2],ap1r[2];
  ldouble fd_u0[NV],fd_up1[NV],fd_up2[NV],fd_um1[NV],fd_um2[NV],fd_r0[NV],fd_rm1[NV],fd_rp1[NV];
  ldouble fd_pLl[NV], fd_pRl[NV], fd_uLr[NV],fd_uLl[NV],fd_uRl[NV],fd_uRr[NV];
  ldouble fd_fstarl[NV],fd_fstarr[NV],fd_dul[3*NV],fd_dur[3*NV],fd_pdiffl[NV],fd_pdiffr[NV];
  ldouble gdet,gg[4][5],GG[4][5];
  //ldouble eup[4][4],elo[4][4];

  // flbx[NV], flby[NV], flbz[NV] are the fluxes at the three walls under consideration
  for(i=0;i<NV;i++) 
  {
    set_ubx(flbx,i,ix,iy,iz,0.);
    set_uby(flby,i,ix,iy,iz,0.);
    set_ubz(flbz,i,ix,iy,iz,0.);
  }

  //**********************************************************************
  //Work on the x-face at ix, iy, iz, which lies in between cells ix-1,iy,iz and ix,iy,iz
 
#ifdef MPI4CORNERS
  if(NX>1 && ix>=0 && ix<=NX && iy>=-1 && iy<NY+1 && iz>=-1 && iz<NZ+1)
#else
  if(NX>1 && ix>=0 && ix<=NX && iy>=0 && iy<NY && iz>=0 && iz<NZ)
#endif
  {
	// Characteristic wave speeds in the two adjoining cells of the current face,
        // which are used for combining the left-biased and right-biased fluxes at the face
        // ap1, am1 correspond to ix and ix-1, i.e., speeds on the right and left of the current face;
	// l and r correspond to left-going and right-going waves; if neither l nor r, it is the maximum speed
        // [0], [1] correspond to hydro and radiation wave speeds      

        ap1l[0]=get_u_scalar(ahdxl,ix,iy,iz);
	ap1r[0]=get_u_scalar(ahdxr,ix,iy,iz);
	ap1l[1]=get_u_scalar(aradxl,ix,iy,iz);
	ap1r[1]=get_u_scalar(aradxr,ix,iy,iz);
	ap1[0]=get_u_scalar(ahdx,ix,iy,iz);
	ap1[1]=get_u_scalar(aradx,ix,iy,iz);
        
	am1l[0]=get_u_scalar(ahdxl,ix-1,iy,iz);
	am1r[0]=get_u_scalar(ahdxr,ix-1,iy,iz);
	am1l[1]=get_u_scalar(aradxl,ix-1,iy,iz);
	am1r[1]=get_u_scalar(aradxr,ix-1,iy,iz);
	am1[0]=get_u_scalar(ahdx,ix-1,iy,iz);
	am1[1]=get_u_scalar(aradx,ix-1,iy,iz);

	//primitives at the face
	for(i=0;i<NV;i++)
	{
          // fd_pLl, fd_pRl are the left-biased and right-biased primitives at the current cell face
          fd_pLl[i]=get_ub(pbLx,i,ix,iy,iz,0);
          fd_pRl[i]=get_ub(pbRx,i,ix,iy,iz,0);
	}

	fill_geometry_face(ix,iy,iz,0,&geom);

	// fd_uLl, fd_uRl are the left-biased and right-biased conserveds at the current cell face
        p2u(fd_pLl,fd_uLl,&geom);
        p2u(fd_pRl,fd_uRl,&geom);
        
#ifdef WAVESPEEDSATFACES //re-calculate the wavespeeds directly at the face
	ldouble aaa[24];
	//left biased wavespeeds
	calc_wavespeeds_lr_pure(fd_uLl,&geom,aaa);
	am1l[0]=aaa[0];
	am1r[0]=aaa[1];
	am1l[1]=aaa[6];
	am1r[1]=aaa[7];
	    
	am1[0]=my_max(fabs(aaa[0]),fabs(aaa[1]));
        am1[1]=my_max(fabs(aaa[6]),fabs(aaa[7]));

	//right biased wavespeeds
	calc_wavespeeds_lr_pure(fd_uRl,&geom,aaa);
	ap1l[0]=aaa[0];
	ap1r[0]=aaa[1];
	ap1l[1]=aaa[6];
	ap1r[1]=aaa[7];

	ap1[0]=my_max(fabs(aaa[0]),fabs(aaa[1]));
	ap1[1]=my_max(fabs(aaa[6]),fabs(aaa[7]));
#endif
   
        // Loop over variables and calculate flux using Lax-Friedrichs or HLL
        for(i=0;i<NV;i++)
        {
          // Choose the proper characteristic speeds: al (left-going wave), ar (right-going wave), ag (maximum wavespeed)
          // Hydro and radiation are treated as two separate systems
#ifdef RADIATION
          if(i<NVMHD) //hydro characteristic speed
          {
            ag=my_max(ap1[0],am1[0]);
            al=my_min(ap1l[0],am1l[0]);
            ar=my_max(ap1r[0],am1r[0]);
            
#ifdef BATTERY
            //when radiation battery on - magnetic fields are affected by radiation
#ifdef BATTERYRADWAVESPEEDS
#ifdef BATTERYRADWAVESPEEDSBONLY
            if(i>=B1 && i<=B3)
#endif
            {
              ag=my_max(ag,my_max(ap1[1],am1[1]));
              al=my_min(al,my_min(ap1l[1],am1l[1]));
              ar=my_max(ar,my_max(ap1r[1],am1r[1]));
            }
#endif
#endif
          }
          else //radiative characteristic speed
          {
            ag=my_max(ap1[1],am1[1]);
            al=my_min(ap1l[1],am1l[1]);
            ar=my_max(ap1r[1],am1r[1]);
          }
          
#else //no RADIATION -- hydro only
          ag=my_max(ap1[0],am1[0]);
          al=my_min(ap1l[0],am1l[0]);
          ar=my_max(ap1r[0],am1r[0]);
#endif
          
          if (FLUXMETHOD==LAXF_FLUX) //Lax-Friedrichs Flux
          {
	    
            //Lax-Friedrichs: Flux = 0.5 * [FR + FL - ag * (UR - UL)]
            fd_fstarl[i] = .5*(get_ub(flRx,i,ix,iy,iz,0) + get_ub(flLx,i,ix,iy,iz,0) - ag * (fd_uRl[i] - fd_uLl[i]));     
            set_ubx(flbx,i,ix,iy,iz,fd_fstarl[i]);
	    
          } 
          
          if (FLUXMETHOD==HLL_FLUX) //HLL Flux
          {
            if(al>0.)
            {
              fd_fstarl[i] = get_ub(flLx,i,ix,iy,iz,0);
            }
            else if(ar<0.)
            {
              fd_fstarl[i] = get_ub(flRx,i,ix,iy,iz,0);
            }
            else
            {
              //HLL: Flux = [ar * FL - al * FR + al * ar * (UR - UL)] / (ar - al)
              fd_fstarl[i] = (-al*get_ub(flRx,i,ix,iy,iz,0) + ar*get_ub(flLx,i,ix,iy,iz,0) + al*ar* (fd_uRl[i] - fd_uLl[i]))/(ar-al);
            }
            set_ubx(flbx,i,ix,iy,iz,fd_fstarl[i]);
          } 
        }  // for(i=0;i<NV;i++)
	
  }  // if(NX>1 && ix>=0 && ix<=NX && iy>=0 && iy<NY && iz>=0 && iz<NZ...)

  
  //**********************************************************************
  //Work on the y-face at ix, iy, iz, which lies in between cells ix,iy-1,iz and ix,iy,iz  
#ifdef MPI4CORNERS
  if(NY>1 && iy>=0 && iy<=NY && ix>=-1 && ix<NX+1 && iz>=-1 && iz<NZ+1)
#else
  if(NY>1 && iy>=0 && iy<=NY  && ix>=0 && ix<NX && iz>=0 && iz<NZ)
#endif
  {
        // Characteristic wave speeds in the two adjoining cells of the current face,
        // which are used for combining the left-biased and right-biased fluxes at the face
        // ap1, am1 correspond to ix and ix-1, i.e., speeds on the right and left of the current face;
        // l and r correspond to left-going and right-going waves; if neither l nor r, it is the maximum speed
        // [0], [1] correspond to hydro and radiation wave speeds
      
	ap1l[0]=get_u_scalar(ahdyl,ix,iy,iz);
        ap1r[0]=get_u_scalar(ahdyr,ix,iy,iz);
	ap1l[1]=get_u_scalar(aradyl,ix,iy,iz);
        ap1r[1]=get_u_scalar(aradyr,ix,iy,iz);
	ap1[0]=get_u_scalar(ahdy,ix,iy,iz);
	ap1[1]=get_u_scalar(arady,ix,iy,iz);

	am1l[0]=get_u_scalar(ahdyl,ix,iy-1,iz);
	am1r[0]=get_u_scalar(ahdyr,ix,iy-1,iz);
	am1l[1]=get_u_scalar(aradyl,ix,iy-1,iz);
	am1r[1]=get_u_scalar(aradyr,ix,iy-1,iz);
	am1[0]=get_u_scalar(ahdy,ix,iy-1,iz);
	am1[1]=get_u_scalar(arady,ix,iy-1,iz);

	for(i=0;i<NV;i++)
	{
          // fd_pLl, fd_pRl are the left-biased and right-biased primitives at the current cell face
          fd_pLl[i]=get_ub(pbLy,i,ix,iy,iz,1);
          fd_pRl[i]=get_ub(pbRy,i,ix,iy,iz,1);
	}

	fill_geometry_face(ix,iy,iz,1,&geom);

	// fd_uLl, fd_uRl are the left-biased and right-biased conserveds at the current cell face
        p2u(fd_pLl,fd_uLl,&geom);
        p2u(fd_pRl,fd_uRl,&geom);
      
#ifdef WAVESPEEDSATFACES // recompute wavespeeds directly at face
	ldouble aaa[24];
	//left-biased wavespeeds
	calc_wavespeeds_lr_pure(fd_uLl,&geom,aaa);
	am1l[0]=aaa[2];
	am1r[0]=aaa[3];
	am1l[1]=aaa[8];
	am1r[1]=aaa[9];
	am1[0]=my_max(fabs(aaa[2]),fabs(aaa[3]));
	am1[1]=my_max(fabs(aaa[8]),fabs(aaa[9]));

	//right-biased wavespeeds
	calc_wavespeeds_lr_pure(fd_uRl,&geom,aaa);
	ap1l[0]=aaa[2];
	ap1r[0]=aaa[3];
	ap1l[1]=aaa[8];
	ap1r[1]=aaa[9];
	ap1[0]=my_max(fabs(aaa[2]),fabs(aaa[3]));
	ap1[1]=my_max(fabs(aaa[8]),fabs(aaa[9]));
#endif
 	    
        // Loop over variables and calculate flux using Lax-Friedrichs or HLL, as required
	for(i=0;i<NV;i++)
	{
            // Choose the proper characteristic speeds: al (left-going wave), ar (right-going wave), ag (maximum wavespeed)
            // Hydro and radiation are treated as two separate systems
#ifdef RADIATION
	    if(i<NVMHD) //hydro characteristic speeds
            {
              ag=my_max(ap1[0],am1[0]);
              al=my_min(ap1l[0],am1l[0]);
              ar=my_max(ap1r[0],am1r[0]);
              
#ifdef BATTERY
              //when radiation battery on - magnetic fields are affected by radiation
#ifdef BATTERYRADWAVESPEEDS
#ifdef BATTERYRADWAVESPEEDSBONLY
              if(i>=B1 && i<=B3)
#endif
              {
                ag=my_max(ag,my_max(ap1[1],am1[1]));
                al=my_min(al,my_min(ap1l[1],am1l[1]));
                ar=my_max(ar,my_max(ap1r[1],am1r[1]));
              }
#endif
#endif
            }
            else //radiative characteristic speeds
            {
              ag=my_max(ap1[1],am1[1]);
              al=my_min(ap1l[1],am1l[1]);
              ar=my_max(ap1r[1],am1r[1]);
            }
            
#else //no RADIATION -- use hydro  wavespeeds
            ag=my_max(ap1[0],am1[0]);
            al=my_min(ap1l[0],am1l[0]);
            ar=my_max(ap1r[0],am1r[0]);
#endif
            
            if (FLUXMETHOD==LAXF_FLUX) //Lax-Friedrichs Flux
            {
	      
              //Lax-Friedrichs: Flux = 0.5 * [FR + FL - ag * (UR - UL)]
              fd_fstarl[i] = .5*(get_ub(flRy,i,ix,iy,iz,1) + get_ub(flLy,i,ix,iy,iz,1) - ag * (fd_uRl[i] - fd_uLl[i]));
              
              set_uby(flby,i,ix,iy,iz,fd_fstarl[i]);
            }
            
            if (FLUXMETHOD==HLL_FLUX) //HLL Flux
            {
              if(al>0.)
              {
                fd_fstarl[i] = get_ub(flLy,i,ix,iy,iz,1);
              }
              else if(ar<0.)
              {
                fd_fstarl[i] = get_ub(flRy,i,ix,iy,iz,1);
              }
              else
              {
                //HLL: Flux = [ar * FL - al * FR + al * ar * (UR - UL)] / (ar - al)
                fd_fstarl[i] = (-al*get_ub(flRy,i,ix,iy,iz,1) + ar*get_ub(flLy,i,ix,iy,iz,1) + al*ar* (fd_uRl[i] - fd_uLl[i]))/(ar-al);
              }
              
              set_uby(flby,i,ix,iy,iz,fd_fstarl[i]);
            } 
	}  // for(i=0;i<NV;i++)
  }  // if(NY>1 && iy>=0 && iy<=NY  && ix>=0 && ix<NX && iz>=0 && iz<NZ...)

  
  //**********************************************************************
  // Work on the z-face at ix, iy, iz, which lies in between cells ix,iy,iz-1 and ix,iy,iz  
#ifdef MPI4CORNERS
  if(NZ>1 && iz>=0 && iz<=NZ && ix>=-1 && ix<NX+1 && iy>=-1 && iy<NY+1)
#else
  if(NZ>1 && iz>=0 && iz<=NZ && ix>=0 && ix<NX && iy>=0 && iy<NY)
#endif
  {
        // Characteristic wave speeds in the two adjoining cells of the current face,
        // which are used for combining the left-biased and right-biased fluxes at the face
        // ap1, am1 correspond to ix and ix-1, i.e., speeds on the right and left of the current face;
        // l and r correspond to left-going and right-going waves; if neither l nor r, it is the maximum speed
        // [0], [1] correspond to hydro and radiation wave speeds

        ap1l[0]=get_u_scalar(ahdzl,ix,iy,iz);
	ap1r[0]=get_u_scalar(ahdzr,ix,iy,iz);
	ap1l[1]=get_u_scalar(aradzl,ix,iy,iz);
	ap1r[1]=get_u_scalar(aradzr,ix,iy,iz);
	ap1[0]=get_u_scalar(ahdz,ix,iy,iz);
	ap1[1]=get_u_scalar(aradz,ix,iy,iz);

	am1l[0]=get_u_scalar(ahdzl,ix,iy,iz-1);
	am1r[0]=get_u_scalar(ahdzr,ix,iy,iz-1);
	am1l[1]=get_u_scalar(aradzl,ix,iy,iz-1);
	am1r[1]=get_u_scalar(aradzr,ix,iy,iz-1);
	am1[0]=get_u_scalar(ahdz,ix,iy,iz-1);
	am1[1]=get_u_scalar(aradz,ix,iy,iz-1);
 
	for(i=0;i<NV;i++)
	{
          // fd_pLl, fd_pRl are the left-biased and right-biased primitives at the current cell face  
          fd_pLl[i]=get_ub(pbLz,i,ix,iy,iz,2);
          fd_pRl[i]=get_ub(pbRz,i,ix,iy,iz,2);
	}

	fill_geometry_face(ix,iy,iz,2,&geom);

	// fd_uLl, fd_uRl are the left-biased and right-biased conserveds at the current cell face
        p2u(fd_pLl,fd_uLl,&geom);
        p2u(fd_pRl,fd_uRl,&geom);
      
#ifdef WAVESPEEDSATFACES // recompute wavespeeds directly at face
	ldouble aaa[24];
	//left-biased wavespeeds
	calc_wavespeeds_lr_pure(fd_uLl,&geom,aaa);
	am1l[0]=aaa[4];
	am1r[0]=aaa[5];
	am1l[1]=aaa[10];
	am1r[1]=aaa[11];
	am1[0]=my_max(fabs(aaa[4]),fabs(aaa[5]));
	am1[1]=my_max(fabs(aaa[10]),fabs(aaa[11]));

	//right-biased wavespeeds
	calc_wavespeeds_lr_pure(fd_uRl,&geom,aaa);
	ap1l[0]=aaa[4];
	ap1r[0]=aaa[5];
	ap1l[1]=aaa[10];
	ap1r[1]=aaa[11];
	ap1[0]=my_max(fabs(aaa[4]),fabs(aaa[5]));
	ap1[1]=my_max(fabs(aaa[10]),fabs(aaa[11]));
#endif

        // Loop over variables and calculate flux using Lax-Friedrichs or HLL, as required
	for(i=0;i<NV;i++)
	{
            // Choose the proper characteristic speeds: al (left-going wave), ar (right-going wave), ag (maximum wavespeed)
            // Hydro and radiation are treated as two separate systems
#ifdef RADIATION
	    if(i<NVMHD) // hydro characteristic speeds
            {
              ag=my_max(ap1[0],am1[0]);
              al=my_min(ap1l[0],am1l[0]);
              ar=my_max(ap1r[0],am1r[0]);
              
#ifdef BATTERY
              //when radiation battery on - magnetic fields are affected by radiation
#ifdef BATTERYRADWAVESPEEDS
#ifdef BATTERYRADWAVESPEEDSBONLY
              if(i>=B1 && i<=B3)
#endif
              {
                ag=my_max(ag,my_max(ap1[1],am1[1]));
                al=my_min(al,my_min(ap1l[1],am1l[1]));
                ar=my_max(ar,my_max(ap1r[1],am1r[1]));
              }
#endif
#endif
            }
            else //radiative characteristic speeds
            {
              ag=my_max(ap1[1],am1[1]);
              al=my_min(ap1l[1],am1l[1]);
              ar=my_max(ap1r[1],am1r[1]);
            }
#else //no radiation -- hydro characteristic speeds
            ag=my_max(ap1[0],am1[0]);
            al=my_min(ap1l[0],am1l[0]);
            ar=my_max(ap1r[0],am1r[0]);
#endif
            
            if (FLUXMETHOD==LAXF_FLUX) //Lax-Friedrichs Flux
            {
	      
              //Lax-Friedrichs: Flux = 0.5 * [FR + FL - ag * (UR - UL)]
              fd_fstarl[i] = .5*(get_ub(flRz,i,ix,iy,iz,2) + get_ub(flLz,i,ix,iy,iz,2) - ag * (fd_uRl[i] - fd_uLl[i]));
              
              set_ubz(flbz,i,ix,iy,iz,fd_fstarl[i]);
            } 
            
            if (FLUXMETHOD==HLL_FLUX) //HLL Flux
            {
              if(al>0.)
              {
                fd_fstarl[i] = get_ub(flLz,i,ix,iy,iz,2);
              }
              else if(ar<0.)
              {
                fd_fstarl[i] = get_ub(flRz,i,ix,iy,iz,2);
              }
              else
              {
                //HLL: Flux = [ar * FL - al * FR + al * ar * (UR - UL)] / (ar - al)  
                fd_fstarl[i] = (-al*get_ub(flRz,i,ix,iy,iz,2) + ar*get_ub(flLz,i,ix,iy,iz,2) + al*ar* (fd_uRl[i] - fd_uLl[i]))/(ar-al);
              }
              
              set_ubz(flbz,i,ix,iy,iz,fd_fstarl[i]);
            } 
	}  // for(i=0;i<NV;i++)
  }  // if(NZ>1 && iz>=0 && iz<=NZ && ix>=0 && ix<NX && iy>=0 && iy<NY...)
	
  
  return 0;
}


//***********************************************************************
/*! \fn int set_grid(ldouble *mindx,ldouble *mindy, ldouble *mindz, ldouble *maxdtfac)
  \brief Sets up the grid, namely locations of cell centers and cell boundaries
 
 \param[out] mindx smallest interval in x
 \param[out] mindy smallest interval in y
 \param[out] mindz smallest interval in z
 \param[out] maxdtfac set to unity (unclear why it is needed here)
 
 The grid is set up as follows:
 
 Global array xb: This corresponds to cell wall locations in the current tile. It is a 1D array that stacks x, y, z wall locations in sequence. 
 The coordinates of the walls are distributed uniformly from MINX to MAXX in x, MINY to MAXY in y, MINZ to MAXZ in z. 
 The first wall is at ix = -NG, the next at ix = -NG+1, ..., NX+NG, followed by iy = -NG, etc., for a total of (NX+NY+NZ)+6*NG+3 entries.
 
 Global array x: This 1D array corresponds to cell centers in the current tile, which are computed by taking averages of wall locations. 
 The first cell is at ix = -NG, next at ix = -NG+1, ..., NX+NG-1, followed by iy = -NG, etc., for a total of (NX+NY+NZ)+6*NG entries.
 
 Note: For a given i, x[i] is the cell center and xb[i] is the LEFT wall of that cell
 
 Note: xb[i], x[i] values are uniformly spaced. Conversion to real coordinates are done through `coco' routines
 
*/
//***********************************************************************
int
set_grid(ldouble *mindx,ldouble *mindy, ldouble *mindz, ldouble *maxdtfac)
{
  int i1,i2,ix,iy,iz;
  ldouble mdx,mdy,mdz,dx,dy,dz,gloc[4][5],xx[4];
  mdx=mdy=mdz=-1;
  ldouble maxdt=-1;

  int ix1,ix2,iy1,iy2,iz1,iz2;
  
  ix1=-0;
  ix2=NX+0;
  iy1=-0;
  iy2=NY+0;
  iz1=-0;
  iz2=NZ+0;
  
  //x
  for(i1=ix1-NG;i1<=ix2+NG;i1++)
    {
      set_xb(i1,0,calc_xb(i1,0));  
      if(i1>-NG) set_x(i1-1,0,.5*(get_xb(i1,0)+get_xb(i1-1,0)));
     }
  //y
  for(i1=iy1-NG;i1<=iy2+NG;i1++)
    {
      set_xb(i1,1,calc_xb(i1,1));  
      if(i1>-NG) set_x(i1-1,1,.5*(get_xb(i1,1)+get_xb(i1-1,1)));
   }
  //z
  for(i1=iz1-NG;i1<=iz2+NG;i1++)
    {
      set_xb(i1,2,calc_xb(i1,2));  
      if(i1>-NG) set_x(i1-1,2,.5*(get_xb(i1,2)+get_xb(i1-1,2)));
    }

  //consistency check
#if(MYCOORDS==BLCOORDS)
  if(get_x(-1,0)<=rhorizonBL)
    {
      printf("ix %d > %f\n",-1,get_x(-1,0));
      my_err("-1 cell inside horizon\n");
    }
#endif

  // what is the minimal cell size
  for(ix=ix1;ix<ix2;ix++)
    {
      for(iy=iy1;iy<iy2;iy++)
	{
	  for(iz=iz1;iz<iz2;iz++)
	    {
	      dx=get_size_x(ix,0);
	      dy=get_size_x(iy,1);
	      dz=get_size_x(iz,2);

	      if((dx<mdx || mdx<0.)) mdx=dx;
	      if((dy<mdx || mdy<0.)) mdy=dy;
	      if((dz<mdx || mdz<0.)) mdz=dz;
	    }
	}
    }  
  
  *mindx=mdx;
  *mindy=mdy;
  *mindz=mdz;
  *maxdtfac=maxdt;

#ifdef PRECOMPUTE_MY2OUT
  set_grid_outcoords();
#endif
  return 0;
}

// sets the output grid coordinates at the cell centers and faces
// used with PRECOMPUTE_MY2OUT defined
int set_grid_outcoords()
{

  int ix,iy,iz,ii;
  #pragma omp parallel for private(ix,iy,iz,ii) 
  for(ix=-NGCX;ix<NX+NGCX;ix++)
    for(iy=-NGCY;iy<NY+NGCY;iy++)
      for(iz=-NGCZ;iz<NZ+NGCZ;iz++)
      {
	ldouble xx[4], xxout[4];
	
	//cell centers
	get_xx(ix,iy,iz,xx);
        coco_N(xx,xxout,MYCOORDS,OUTCOORDS);

	for(ii=0;ii<3;ii++)
	  set_xout(ii,ix,iy,iz,xxout[ii+1]);

	// TODO these are only needed when postproc==1 ? 
        //x-faces
	if(ix==-NG)
	{
          xx[0] = global_time; xx[1] = get_xb(ix, 0); xx[2] = get_x(iy, 1); xx[3] = get_x(iz, 2);
          coco_N(xx,xxout,MYCOORDS,OUTCOORDS);
	  for(ii=0;ii<3;ii++)
	    set_xbout_xface(ii,ix,iy,iz,xxout[ii+1]);
	}
	xx[0] = global_time; xx[1] = get_xb(ix+1,0); xx[2] = get_x(iy,1); xx[3] = get_x(iz,2);
	coco_N(xx,xxout,MYCOORDS,OUTCOORDS);
        for(ii=0;ii<3;ii++)
	  set_xbout_xface(ii,ix+1,iy,iz,xxout[ii+1]);
	  
	//y-faces
	if(iy==-NG)
	{
          xx[0] = global_time; xx[1] = get_x(ix,0); xx[2] = get_xb(iy,1); xx[3] = get_x(iz,2);
          coco_N(xx,xxout,MYCOORDS,OUTCOORDS);
	  for(ii=0;ii<3;ii++)
	    set_xbout_yface(ii,ix,iy,iz,xxout[ii+1]);
	}
	xx[0] = global_time; xx[1] = get_x(ix,0); xx[2] = get_xb(iy+1,1); xx[3] = get_x(iz,2);
	coco_N(xx,xxout,MYCOORDS,OUTCOORDS);
	for(ii=0;ii<3;ii++)
	  set_xbout_yface(ii,ix,iy+1,iz,xxout[ii+1]);

	//z-faces
	if(iz==-NG)
	{
          xx[0] = global_time; xx[1] = get_x(ix,0); xx[2] = get_x(iy,1); xx[3] = get_xb(iz,2);
          coco_N(xx,xxout,MYCOORDS,OUTCOORDS);
	  for(ii=0;ii<3;ii++)
	    set_xbout_zface(ii,ix,iy,iz,xxout[ii+1]);
	}
	xx[0] = global_time; xx[1] = get_x(ix,0); xx[2] = get_x(iy,1); xx[3] = get_xb(iz+1,2);
	coco_N(xx,xxout,MYCOORDS,OUTCOORDS);
        for(ii=0;ii<3;ii++)
	  set_xbout_zface(ii,ix,iy,iz+1,xxout[ii+1]);	  
     }

 return 0;
}

//**********************************************************************
/*! \fn int alloc_loops(int init,ldouble t,ldouble dt)
 \brief Sets up loops 0 to 6 to speed up parallel for loops
 
 loop_0[i][0,1,2] = ix, iy, iz of cells within the domain of the present tile\n
 Nloop_0 = number of cells in loop_0 = NX*NY*NZ
 
 loop_1[i][0,1,2] = ix, iy, iz of cells within domain plus ghost cells of width lim\n
 Nloop_1 = number of cells in loop_1 = (NX+2*lim)*(NY+2*lim)*(NZ+2*lim) in 3D
 
 loop_2[i][0,1,2] = ghost cells only, without corners\n
 Nloop_2 = number of ghost cells without corners
 
 loop_02[i][0,1,2] = domain plus ghost cells, but no corners\n
 Nloop_02
 
 loop_4[i][0,1,2] = unclear
 
 loop_5[i]0,1,2] = all cells: domain+ghost+corners. Is it different from loop_1?\n
 
 loop_6[i][0,1,2] = ix, iy, iz of cells within domain plus 1 cell all round (like loop_1 with lim=1)\n
 Nloop_6 = number of cells in loop_1 = (NX+2)*(NY+2)*(NZ+2) in 3D
 
 \todo Go through the definitions of the various loop_i and make sure that their respective ranges (domain, lim, ghost, corner) are consistent.
 */
//**********************************************************************

int
alloc_loops()
{
  int zone=-1;
  int ix,iy,iz,i,ii,jj ;
  int ix1,ix2,iy1,iy2,iz1,iz2,szix1,szix2;
  int toi,tsi,tnx,imaxx=-1;

  //by default cover the whole local tile
  ix1=0;
  iy1=0;
  iz1=0;

  ix2=NX;
  iy2=NY;
  iz2=NZ;  

  szix1=ix1;
  szix2=ix2;
    
    //**********************************************************************
    //inside domain only
    Nloop_0=0;

    if((loop_0=(int **)malloc(SX*SY*SZ*sizeof(int*)))==NULL) my_err("malloc err. - loops 1\n");
    for(i=0;i<SX*SY*SZ;i++) if((loop_0[i]=(int *)malloc(3*sizeof(int)))==NULL) my_err("malloc err. - loops 2\n");

    for(ix=ix1;ix<ix2;ix++)
      {
	for(iy=iy1;iy<iy2;iy++)
	  {
	    for(iz=iz1;iz<iz2;iz++)
	      {	
                #ifdef SPECIAL_BC_CHECK
                #include PR_BC_SPECIAL_LOOP
                if(ret_val != 0) continue;
                #endif

		loop_0[Nloop_0][0]=ix;
		loop_0[Nloop_0][1]=iy;
		loop_0[Nloop_0][2]=iz;

		Nloop_0++;
	      }
	  }
      }
 
    //shuffling:
    #if (SHUFFLELOOPS)
    shuffle_loop(loop_0,Nloop_0);
    #endif 
 
    //**********************************************************************
    //inside + ghost cells - number depending on the order of reconstruction
    //used to indicate where calculate fluxes
    //excluding corners (does it actually exclude corners??)
    int xlim1,xlim2,ylim1,ylim2,zlim1,zlim2;
    int xlim,ylim,zlim;
    int lim;

    if(INT_ORDER==1) lim=1;
    if(INT_ORDER==2) lim=1;
    if(INT_ORDER==4) lim=2;

    if(TNX>1) xlim1=xlim2=lim; else xlim1=xlim2=0;  
    if(TNY>1) ylim1=ylim2=lim; else ylim1=ylim2=0;
    if(TNZ>1) zlim1=zlim2=lim; else zlim1=zlim2=0;

    Nloop_1=0;
    if((loop_1=(int **)malloc(SX*SY*SZ*sizeof(int*)))==NULL) my_err("malloc err. - loops 3\n");
    for(i=0;i<SX*SY*SZ;i++) if((loop_1[i]=(int *)malloc(3*sizeof(int)))==NULL) my_err("malloc err. - loops 4\n");

    for(ix=-xlim1+ix1;ix<ix2+xlim2;ix++)
      {
	for(iy=-ylim1+iy1;iy<iy2+ylim2;iy++)
	  {
	    for(iz=-zlim1+iz1;iz<iz2+zlim2;iz++)
	      {	
                #ifdef SPECIAL_BC_CHECK
                #include PR_BC_SPECIAL_LOOP
                if(ret_val == 4) continue;
                #endif

		loop_1[Nloop_1][0]=ix;
		loop_1[Nloop_1][1]=iy;
		loop_1[Nloop_1][2]=iz;

		Nloop_1++;
	      
		//loop_1=(int **)realloc(loop_1,(Nloop_1+1)*sizeof(int*));
		//loop_1[Nloop_1]=(int *)malloc(3*sizeof(int));	      
	      }
	  }
      }

    //shuffling:
    #if (SHUFFLELOOPS)
    shuffle_loop(loop_1,Nloop_1);
    #endif

    //**********************************************************************
    //**********************************************************************
    //only ghost cells (with no corners)

    //reduction of size basing on the dimension
    //for constrained transport fluxes copied onto missing dimensions
    if(TNX>1) xlim1=xlim2=NG; else xlim1=xlim2=0;  
    if(TNY>1) ylim1=ylim2=NG; else ylim1=ylim2=0;
    if(TNZ>1) zlim1=zlim2=NG; else zlim1=zlim2=0;

    Nloop_2=0;
    if((loop_2=(int **)malloc(SX*SY*SZ*sizeof(int*)))==NULL) my_err("malloc err. - loops 5\n");
    for(i=0;i<SX*SY*SZ;i++) if((loop_2[i]=(int *)malloc(3*sizeof(int)))==NULL) my_err("malloc err. - loops 6\n");

    for(ix=-xlim1+ix1;ix<ix2+xlim2;ix++)
      {
	for(iy=-ylim1+iy1;iy<iy2+ylim2;iy++)
	  {
	    for(iz=-zlim1+iz1;iz<iz2+zlim2;iz++)
	      {	 
		//within domain:
                #ifdef SPECIAL_BC_CHECK
                #include PR_BC_SPECIAL_LOOP
                if(ret_val == 0) //"globally" not in corner or domain, but maybe not on tile - Brandon
                { 
                  if(if_indomain(ix,iy,iz)==1) continue; //in domain on tile
		  if(if_outsidegc(ix,iy,iz)==1) continue; //corner on tile
                }
                else if(ret_val == 4)
                {
                  continue;
                }
                #else
		if(if_indomain(ix,iy,iz)==1) continue;
		//but not at the corners
		if(if_outsidegc(ix,iy,iz)==1) continue;
                #endif

		loop_2[Nloop_2][0]=ix;
		loop_2[Nloop_2][1]=iy;
		loop_2[Nloop_2][2]=iz;

		Nloop_2++;
	      
		//loop_2=(int **)realloc(loop_2,(Nloop_2+1)*sizeof(int*));
		//loop_2[Nloop_2]=(int *)malloc(3*sizeof(int));
	      }
	  }
      }	
 
    //shuffling:
    #if (SHUFFLELOOPS)
    shuffle_loop(loop_2,Nloop_2);
    #endif

    //**********************************************************************
    //**********************************************************************
    //domain and all ghost cells (no corners)
    Nloop_02=Nloop_0+Nloop_2;
    
    if((loop_02=(int **)malloc(SX*SY*SZ*sizeof(int*)))==NULL) my_err("malloc err. - loops 7\n");
    for(i=0;i<SX*SY*SZ;i++) if((loop_02[i]=(int *)malloc(3*sizeof(int)))==NULL) my_err("malloc err. - loops 8\n");

    for(ix=0;ix<Nloop_0;ix++)
      {
	loop_02[ix][0]=loop_0[ix][0];
	loop_02[ix][1]=loop_0[ix][1];
	loop_02[ix][2]=loop_0[ix][2];
      }
    for(ix=0;ix<Nloop_2;ix++)
      {
	loop_02[ix+Nloop_0][0]=loop_2[ix][0];
	loop_02[ix+Nloop_0][1]=loop_2[ix][1];
	loop_02[ix+Nloop_0][2]=loop_2[ix][2];
      }
  
    //shuffling:
    #if (SHUFFLELOOPS)
    shuffle_loop(loop_02,Nloop_02);
    #endif

    //**********************************************************************
    //**********************************************************************
    //1-deep surfaces on corners only
    /*
    if(NX>1) xlim=NG; else xlim=0;  
    if(NY>1) ylim=NG; else ylim=0;
    if(NZ>1) zlim=NG; else zlim=0;

    Nloop_3=0;
    loop_3=(int **)malloc(SX*SY*SZ*sizeof(int*));
    for(i=0;i<SX*SY*SZ;i++) loop_3[i]=(int *)malloc(3*sizeof(int));
     //test
    for(ix=-0;ix<NX+0;ix++)
      {
	for(iy=-ylim;iy<NY+ylim;iy++)
	  {
	    for(iz=-zlim;iz<NZ+zlim;iz++)
	      {
		loop_3[Nloop_3][0]=ix;
		loop_3[Nloop_3][1]=iy;
		loop_3[Nloop_3][2]=iz;

		Nloop_3++;
	      
		loop_3=(int **)realloc(loop_3,(Nloop_3+1)*sizeof(int*));
		loop_3[Nloop_3]=(int *)malloc(3*sizeof(int));	      
	      }
	  }
      }

    //shuffling:
#if (SHUFFLELOOPS)
    shuffle_loop(loop_3,Nloop_3);
#endif
*/
    
    //**********************************************************************
    //**********************************************************************
    //all corners of the domain - like in staggered grid
    //ANDREW -- this doesn't seem right...
    Nloop_4=0;
    if((loop_4=(int **)malloc((SX+1)*(SY+1)*(SZ+1)*sizeof(int*)))==NULL) my_err("malloc err. - loops 9\n");
    for(i=0;i<(SX+1)*(SY+1)*(SZ+1);i++) if((loop_4[i]=(int *)malloc(3*sizeof(int)))==NULL) my_err("malloc err. - loops 10\n");

    xlim2=ix2;
    ylim2=iy2;
    zlim2=iz2;
    //if(TNY>1) ylim2=iy2; else ylim2=iy1;
    //if(TNZ>1) zlim2=iz2; else zlim2=iz1;

    for(ix=ix1;ix<=xlim2;ix++)
      {
	for(iy=iy1;iy<=ylim2;iy++)
	  {
	    for(iz=iz1;iz<=zlim2;iz++)
	      {	
		loop_4[Nloop_4][0]=ix;
		loop_4[Nloop_4][1]=iy;
		loop_4[Nloop_4][2]=iz;

		Nloop_4++;
	      
		//loop_4=(int **)realloc(loop_4,(Nloop_4+1)*sizeof(int*));
		//loop_4[Nloop_4]=(int *)malloc(3*sizeof(int));	      
	      }
	  }
      }

    //shuffling:
    #if (SHUFFLELOOPS)
    shuffle_loop(loop_4,Nloop_4);
    #endif
       
    //**********************************************************************
    //**********************************************************************
    //domain + ghost cells + corners  = total  
    Nloop_5=0;

    if((loop_5=(int **)malloc(SX*SY*SZ*sizeof(int*)))==NULL) my_err("malloc err. - loops 11\n");
    for(i=0;i<SX*SY*SZ;i++) if((loop_5[i]=(int *)malloc(3*sizeof(int)))==NULL) my_err("malloc err. - loops 12\n");
  
    if(TNX>1) xlim1=xlim2=NGCX; else xlim1=xlim2=0;  
    if(TNY>1) ylim1=ylim2=NGCY; else ylim1=ylim2=0;
    if(TNZ>1) zlim1=zlim2=NGCZ; else zlim1=zlim2=0;


    for(ix=-xlim1+ix1;ix<ix2+xlim2;ix++)
      {
	for(iy=-ylim1+iy1;iy<iy2+ylim2;iy++)
	  {
	    for(iz=-zlim1+iz1;iz<iz2+zlim2;iz++)
	      {	
		loop_5[Nloop_5][0]=ix;
		loop_5[Nloop_5][1]=iy;
		loop_5[Nloop_5][2]=iz;

		Nloop_5++;
	      
		//loop_5=(int **)realloc(loop_5,(Nloop_5+1)*sizeof(int*));
		//loop_5[Nloop_5]=(int *)malloc(3*sizeof(int));	      
	      }
	  }
      }

    //shuffling:
#if (SHUFFLELOOPS)
    shuffle_loop(loop_5,Nloop_5);
#endif
      
    //**********************************************************************
    //**********************************************************************
    //inner domain plus 1-cell layer including corners
    Nloop_6=0;
    if((loop_6=(int **)malloc(SX*SY*SZ*sizeof(int*)))==NULL) my_err("malloc err. - loops 13\n");
    for(i=0;i<SX*SY*SZ;i++) if((loop_6[i]=(int *)malloc(3*sizeof(int)))==NULL) my_err("malloc err. - loops 14\n");
  
    if(TNX>1) xlim1=xlim2=1; else xlim1=xlim2=0;  
    if(TNY>1) ylim1=ylim2=1; else ylim1=ylim2=0;
    if(TNZ>1) zlim1=zlim2=1; else zlim1=zlim2=0;

    for(ix=-xlim1+ix1;ix<ix2+xlim2;ix++)
      {
	for(iy=-ylim1+iy1;iy<iy2+ylim2;iy++)
	  {
	    for(iz=-zlim1+iz1;iz<iz2+zlim2;iz++)
	      {
		loop_6[Nloop_6][0]=ix;
		loop_6[Nloop_6][1]=iy;
		loop_6[Nloop_6][2]=iz;

		Nloop_6++;
	      
		//loop_6=(int **)realloc(loop_6,(Nloop_6+1)*sizeof(int*));
		//loop_6[Nloop_6]=(int *)malloc(3*sizeof(int));	      
	      }
	  }
      }

    //shuffling:
#if (SHUFFLELOOPS)
    shuffle_loop(loop_6,Nloop_6);
#endif

return zone;
}


//**********************************************************************
//* prints grid **********************************************************
//**********************************************************************
int
print_grid(ldouble min_dx, ldouble min_dy, ldouble min_dz)
{
  int i1;
  for(i1=-NG;i1<=NX+NG-1;i1++)
    printf("x: %6d %8.3f|%8.3f|%8.3f (%8.3f)\n",i1,get_xb(i1,0),get_x(i1,0),get_xb(i1+1,0),get_size_x(i1,0));
  printf("\n");
  for(i1=-NG;i1<=NY+NG-1;i1++)
    printf("y: %6d %8.3f|%8.3f|%8.3f (%8.3f)\n",i1,get_xb(i1,1),get_x(i1,1),get_xb(i1+1,1),get_size_x(i1,1));
  printf("\n");
  for(i1=-NG;i1<=NZ+NG-1;i1++)
    printf("z: %6d %8.3f|%8.3f|%8.3f (%8.3f)\n",i1,get_xb(i1,2),get_x(i1,2),get_xb(i1+1,2),get_size_x(i1,2));

  printf("\n min_dx = %8.3f\n min_dy = %8.3f\n min_dz = %8.3f\n",min_dx,min_dy,min_dz);
      
  //getchar(); 

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************

//routines to handle arrays - many of them replaced by defs in ko.h

//**********************************************************************
//**********************************************************************
//**********************************************************************

//*********************************************
//sets cell center location
//*********************************************
int set_x(int ic, int idim, ldouble val)
{  
  if(idim==0)
    x[ic+NG]=val;
  if(idim==1)
    x[ic+NG + NX+2*NG]=val;
  if(idim==2)
    x[ic+NG + NX+2*NG + NY+2*NG]=val;
  return 0;
}

//*********************************************
//* sets locations of cell boundaries *********
//*********************************************
int set_xb(int ic, int idim,ldouble val)
{  
  if(idim==0)
    xb[ic+NG]=val;
  if(idim==1)
    xb[ic+NG + NX+2*NG +1]=val;
  if(idim==2)
    xb[ic+NG + NX+2*NG + 1 + NY+2*NG + 1]=val;
  return 0;
}

//*********************************************
//* returns size of cell **********************
//*********************************************
ldouble get_size_x(int ic, int idim)
{
  return get_xb(ic+1,idim)-get_xb(ic,idim);
}

//***********************************************************
//returns four-vector of coordinates
//***********************************************************
int
get_xx(int ix,int iy, int iz, ldouble *xx)
{
  xx[0]=global_time;
  xx[1]=get_x(ix,0);
  xx[2]=get_x(iy,1);
  xx[3]=get_x(iz,2);
  return 0;
}

// in OUTCOORDS
// NOTE that the time coordinate is useless here!!
int
get_xxout(int ix, int iy, int iz, ldouble *xx)
{
  xx[0]=-1;
  xx[1]=get_xout(0,ix,iy,iz);
  xx[2]=get_xout(1,ix,iy,iz);
  xx[3]=get_xout(2,ix,iy,iz);
  return 0;
}

//***********************************************************
//returns four-vector of cell-centered coordinates in arbitrary coordinates
//***********************************************************
int 
get_xx_arb(int ix,int iy,int iz,ldouble *xx,int COORDSOUT)
{
  
#ifdef PRECOMPUTE_MY2OUT // use precomputed coordinates if COORDS == OUTCOORDS
  if(COORDSOUT == OUTCOORDS)
  {
    get_xxout(ix, iy, iz, xx); // time will be nonsense! seems ok everywhere this is used
  }
  else
  {
    ldouble xx0[4];
    get_xx(ix,iy,iz,xx0);
    coco_N(xx0,xx,MYCOORDS,COORDSOUT);
  }
#else
    ldouble xx0[4];
    get_xx(ix,iy,iz,xx0);
    coco_N(xx0,xx,MYCOORDS,COORDSOUT);
#endif  
  return 0;
}


//*********************************************
//* returns locations of cell boundaries ******
//* in equally spaced code coordinates ********
//*********************************************
ldouble
calc_xb(int i,int idim)
{
  ldouble c0,c1,dc,xb;
  if(idim==0)
  {
    c0=MINX;
    c1=MAXX;
    dc=(c1-c0)/(ldouble)TNX;
    xb=c0+(ldouble)(i+TOI)*dc;
  }
  if(idim==1)
  {
    c0=MINY;
    c1=MAXY;
    dc=(c1-c0)/(ldouble)TNY;
    xb=c0+(ldouble)(i+TOJ)*dc;
  }
  if(idim==2)
  {
    c0=MINZ;
    c1=MAXZ;
    dc=(c1-c0)/(ldouble)TNZ;
    xb=c0+(ldouble)(i+TOK)*dc;
  }
  
  return xb;
} 

///////////////////////////////////////////////////////////////
//deals with arrays [NX+NG x NY+NG x NZ+NG x gSIZE] - cell centers metric
/*
int set_g(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG)
    my_err("problem w/ set_g - index out of range");
  
  uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + \
               (iY(iy)+NGCY)*(SX)*gSIZE + \
               (iZMET(iz)+NGCZMET)*(SY)*(SX)*gSIZE \
      ] = value;
  return 0;
}
*/

///////////////////////////////////////////////////////////////
//deals with arrays [NX+NG x NY+NG x NZ+NG x 16] - 4x4 tensors
//metric specific 
/*
int set_T(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w set_T - index ouf of range");
  
  uarr[i*4+j + (ix+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZMET(iz)+NGCZMET)*(SY)*(SX)*16] = value;
  return 0;
}
*/

///////////////////////////////////////////////////////////////
//deals with arrays [NX+NG x NY+NG x NZ+NG x 16] - 4x4 tensors
/*
int set_Tfull(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w set_T - index ouf of range");
  
  uarr[i*4+j + (ix+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZ(iz)+NGCZ)*(SY)*(SX)*16] = value;
  return 0;
}
*/

//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x NV] - cell boundaries in idim

int set_ub(ldouble* uarr,int iv,int ix,int iy,int iz,ldouble value,int idim)
{
  if(idim==0)
    {
      uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX+1)*NV + (iZ(iz)+NGCZ)*(SY)*(SX+1)*NV] = value;
    }
  if(idim==1)
    {
      uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY+1)*(SX)*NV] = value;
    }
  if(idim==2)
    {
      uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY)*(SX)*NV] = value;
    }
  return 0;
}


///////////////////////////////////////////////////////////////
//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x gSIZE] - cell boundaries in idim metric

int set_gb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim)
{
  if(idim==0)
    {
      uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX+1)*gSIZE + (iZMET(iz)+NGCZMET)*(SY)*(SX+1)*gSIZE] = value;
    }
  if(idim==1)
    {
      uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX)*gSIZE + (iZMET(iz)+NGCZMET)*(SY+1)*(SX)*gSIZE] = value;
    }
  if(idim==2)
    {
      uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX)*gSIZE + (iZMET(iz)+NGCZMET)*(SY)*(SX)*gSIZE] = value;
    }
  return 0;
}


///////////////////////////////////////////////////////////////
//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x 16] - tensors at cell boundaries in idim metric

int set_Tb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim)
{
  if(idim==0)
    {
      uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX+1)*16 + (iZMET(iz)+NGCZMET)*(SY)*(SX+1)*16] = value;
    }
  if(idim==1)
    {
      uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZMET(iz)+NGCZMET)*(SY+1)*(SX)*16] = value;
    }
  if(idim==2)
    {
      uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZMET(iz)+NGCZMET)*(SY)*(SX)*16] = value;
    }
  return 0;
}

///////////////////////////////////////////////////////////////
//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x gSIZE] - cell boundaries in idim metric
/*
ldouble get_gb(ldouble* uarr,int i,int j,int ix,int iy,int iz,int idim)
{
  if(idim==0)
    {
      if(ix<-NG || ix>NX+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_gb x - index ouf of range");  
      return uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG+1)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*gSIZE];
    }
  if(idim==1)
    {
      if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_gb y - index ouf of range");  
      return uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*gSIZE];
    }
  if(idim==2)
    {
      if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ+NG) my_err("blont w get_gb z - index ouf of range");  
      return uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG)*gSIZE];
    }
  return 0;
}
*/

///////////////////////////////////////////////////////////////
//deals with arrays [NX+NG x NY+NG x NZ+NG] - cell centers 
/*
int set_u_scalar(ldouble* uarr,int ix,int iy,int iz,ldouble value)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w set_u_scalar - index ouf of range");

  uarr[ix+NG + (iy+NG)*(NX+2*NG) + (iz+NG)*(NY+2*NG)*(NX+2*NG)] = value;
  return 0;
}
*/

///////////////////////////////////////////////////////////////
//deals with arrays [NX+NG x NY+NG x NZ+NG x NV] - cell centers 
/*
ldouble get_u_scalar(ldouble* uarr,int ix,int iy,int iz)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_u_scalar - index ouf of range");
  
  //TODO something better, so far skipping calculating wave speeds at ghosts
  //as it is easier due to extrapolation of primitives quantities 

  //printf("%4d %4d %4d\n",ix,iy,iz); getchar();
  
  return uarr[ix+NG + (iy+NG)*(NX+2*NG) + (iz+NG)*(NY+2*NG)*(NX+2*NG)];
}
*/



//**********************************************************************
//* Variable copying, multiplication, addition
//**********************************************************************


///////////////////////////////////////////////////////////////
//array multiplication
//uu2=factor*uu1
int
copy_u_core(ldouble factor,ldouble *uu1,ldouble* uu2, long long N)
{
  long long i;
#ifdef APPLY_OMP_SIMD
  #pragma omp parallel for simd
#else
  #pragma omp parallel for private(i)
#endif
  for (i=0;i<N;i++)
    uu2[i]=uu1[i]*factor;
  return 0;
}


///////////////////////////////////////////////////////////////
int
copy_u(ldouble factor,ldouble *uu1,ldouble* uu2 )
{
  long long Ngrid=SX*SY*SZ;
  long long Nprim=Ngrid*NV;
  copy_u_core(factor,uu1,uu2,Nprim);
  return 0;
}


///////////////////////////////////////////////////////////////
int
copyi_u(ldouble factor,ldouble *uu1,ldouble* uu2)	
{
  int ii;
#ifdef APPLY_OMP_SIMD
  #pragma omp parallel for simd
#else
  #pragma omp parallel for private(ii)
#endif
  for(ii=0;ii<Nloop_02;ii++) //domain + ghost cells, but no corners
  {
    int ix,iy,iz,iv;

    ix=loop_02[ii][0];
    iy=loop_02[ii][1];
    iz=loop_02[ii][2];
    
    PLOOP(iv)
    set_u(uu2,iv,ix,iy,iz,factor*get_u(uu1,iv,ix,iy,iz));
  }

  return 0;
}

///////////////////////////////////////////////////////////////
//array multiplication plus addition
//uu3=f1*uu1+f2*uu2
int
add_u_core(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble *uu3, long long N)
{
  long long i;
#ifdef APPLY_OMP_SIMD
  #pragma omp parallel for simd private(i)
#else
  #pragma omp parallel for private(i)
#endif
  for (i=0;i<N;i++)
    uu3[i]=uu1[i]*f1+uu2[i]*f2;
  return 0;
}


///////////////////////////////////////////////////////////////

int
add_u(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble *uu3)
{
  long long Ngrid=SX*SY*SZ;
  long long Nprim=Ngrid*NV;
  add_u_core(f1,uu1,f2,uu2,uu3,Nprim);
  return 0;
}


///////////////////////////////////////////////////////////////
int
addi_u(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble *uu3)
{
  int ii;
#ifdef APPLY_OMP_SIMD
  #pragma omp parallel for simd
#else
  #pragma omp parallel for
#endif
  for(ii=0;ii<Nloop_0;ii++) //domain only
    {
      int ix,iy,iz,iv;
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2];
      PLOOP(iv)
	set_u(uu3,iv,ix,iy,iz,f1*get_u(uu1,iv,ix,iy,iz)+f2*get_u(uu2,iv,ix,iy,iz));
    }

  return 0;
}



///////////////////////////////////////////////////////////////
//array multiplication plus addition on 3 matrices
//uu3=f1*uu1+f2*uu2
int
add_u_core_3(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble f3, ldouble *uu3, ldouble *uu4,long long N)
{
  long long i;
#ifdef APPLY_OMP_SIMD
  #pragma omp parallel for simd private(i)
#else
  #pragma omp parallel for private(i)
#endif
  for (i=0;i<N;i++)
    uu4[i]=uu1[i]*f1+uu2[i]*f2+uu3[i]*f3;
  return 0;
}


///////////////////////////////////////////////////////////////
int
add_u_3(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble f3, ldouble *uu3, ldouble *uu4)
{
  long long Ngrid=SX*SY*SZ;
  long long Nprim=Ngrid*NV;
  add_u_core_3(f1,uu1,f2,uu2,f3,uu3,uu4,Nprim);
  return 0;
}


///////////////////////////////////////////////////////////////
int
addi_u_3(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble f3, ldouble *uu3, ldouble *uu4)
{
  int ii;
#ifdef APPLY_OMP_SIMD
  #pragma omp parallel for simd
#else
  #pragma omp parallel for
#endif
  for(ii=0;ii<Nloop_0;ii++) //domain only
    {
      int ix,iy,iz,iv;
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2];
      PLOOP(iv)
	set_u(uu4,iv,ix,iy,iz,f1*get_u(uu1,iv,ix,iy,iz)+f2*get_u(uu2,iv,ix,iy,iz)+f3*get_u(uu3,iv,ix,iy,iz));
    }

  return 0;
}


//**********************************************************************
// Domain checks
//**********************************************************************

//checks if cell is inside main domain
int
if_indomain(int ix,int iy,int iz)
{
  if(ix>=0 && ix<NX && iy>=0 && iy<NY && iz>=0 && iz<NZ) return 1;
  else return 0;
}


///////////////////////////////////////////////////////////////
//checks if cell outside both domain and ghostcells, i.e. if cells in corners

int
if_outsidegc(int ix,int iy,int iz)
{  
  if(((ix<0 || ix>=NX) && (iy>=0 && iy<NY) && (iz>=0 && iz<NZ)) ||
     ((ix>=0 && ix<NX) && (iy<0 || iy>=NY) && (iz>=0 && iz<NZ)) || 
     ((ix>=0 && ix<NX) && (iy>=0 && iy<NY) && (iz<0 || iz>=NZ)) ||
     (ix>=0 && ix<NX && iy>=0 && iy<NY && iz>=0 && iz<NZ))
    return 0;
  else
    return 1;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************



/*************************************************************/
/* returns problem specific BC defined in PROBLEMS/XXX/bc.c */
/*************************************************************/

int
calc_bc(int ix,int iy,int iz,ldouble t,ldouble *uu,ldouble *pp,int ifinit,int BCtype)
{
  
#include PR_BC
  
  return 0;
} 

// boundary conditions - sets conserved in a given ghost cell
int set_bc_core(int ix,int iy,int iz,double t,ldouble *uval,ldouble *pval,int ifinit,int BCtype)
{
  int iix,iiy,iiz,iv;
  iix=ix;
  iiy=iy;
  iiz=iz;
          
#ifdef SPECIFIC_BC  // BC specific for given problem
  calc_bc(ix,iy,iz,t,uval,pval,ifinit,BCtype);
#else   //standard BC  
 
  if(BCtype==XBCLO || BCtype==XBCHI)
    {       
#ifdef PERIODIC_XBC
      iix=ix;
      if(ix<0) iix=ix+NX;
      if(ix>NX-1) iix=ix-NX;
#endif
#ifdef COPY_XBC
      iix=ix;
      if(ix<0) iix=0;
      if(ix>NX-1) iix=NX-1;
#endif
    }

  if(BCtype==YBCLO || BCtype==YBCHI)
    {       
#ifdef PERIODIC_YBC
      iiy=iy;
      if(iy<0) iiy=iy+NY;
      if(iy>NY-1) iiy=iy-NY;
      if(NY<NG) iiy=0;
#endif
#ifdef COPY_YBC
      iiy=iy;
      if(iy<0) iiy=0;
      if(iy>NY-1) iiy=NY-1;
#endif
    }

  if(BCtype==ZBCLO || BCtype==ZBCHI)
    {       
#ifdef PERIODIC_ZBC
      iiz=iz;
      if(iz<0) iiz=iz+NZ;
      if(iz>NZ-1) iiz=iz-NZ;
      if(NZ<NG) iiz=0;
#endif
#ifdef COPY_ZBC
      iiz=iz;
      if(iz<0) iiz=0;
      if(iz>NZ-1) iiz=NZ-1;
#endif
    }


  for(iv=0;iv<NV;iv++)
    pval[iv]=get_u(p,iv,iix,iiy,iiz);
    
  #ifdef CONSISTENTGAMMA
  set_u_scalar(gammagas, ix, iy, ix, get_u_scalar(gammagas,iix,iiy,iiz));
  #endif

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  p2u(pval,uval,&geom);
     
#endif //SPECIFIC_BC   

  return 0;
}

// boundary conditions - sets conserved in all the ghost cells
int set_bc(ldouble t,int ifinit)
{
  int ii;
  
  //first fill the GC with no corners
#ifdef APPLY_OMP_SIMD
  #pragma omp parallel for simd schedule(static)
#else
  #pragma omp parallel for schedule(static)
#endif
  for(ii=0;ii<Nloop_2;ii++) //ghost cells only, no corners
    {
      int ix,iy,iz,iv;

      ix=loop_2[ii][0];
      iy=loop_2[ii][1];
      iz=loop_2[ii][2];

      //to avoid corners here - treated later if necessary
      if(if_outsidegc(ix,iy,iz)==1) continue;

      //type of BC
      int BCtype=-1;
      if(ix<0) BCtype=XBCLO;
      if(ix>=NX) BCtype=XBCHI;
      if(iy<0) BCtype=YBCLO;
      if(iy>=NY) BCtype=YBCHI;
      if(iz<0) BCtype=ZBCLO;
      if(iz>=NZ) BCtype=ZBCHI;

      #ifdef SPECIAL_BC_CHECK
      //need special check in mpi for if BCs are real BCs
      int is_disk_XBC = 0;
      int is_disk_YBC = 0;
      int is_disk_ZBC = 0;
      //if( (ix+TOI) == STREAM_IX || (ix+TOI) == (STREAM_IX+1))
      if( (ix+TOI) == STREAM_IX || (ix+TOI) == (STREAM_IX+1) || (ix+TOI) == (STREAM_IX+2) || (ix+TOI) == (STREAM_IX+3))
      {
        if(TNY>1 && TNZ==1) //2D r-theta
        { 
          if( (iy+TOJ) >= STREAM_IYT && (iy+TOJ) <= STREAM_IYB ) is_disk_XBC = 1;
        }
        else if(TNY==1 && TNZ>1) //2D r-phi
        { 
          if( (iz+TOK) >= STREAM_IZT && (iz+TOK) <= STREAM_IZB ) is_disk_XBC = 1;
        }
        else if(TNY>1 && TNZ>1) //2D r-phi
        { 
          #ifndef STREAM_RING
          if( (iy+TOJ) >= STREAM_IYT && (iy+TOJ) <= STREAM_IYB && (iz+TOK) >= STREAM_IZT && (iz+TOK) <= STREAM_IZB ) is_disk_XBC = 1;
          #else
          if( (iy+TOJ) >= STREAM_IYT && (iy+TOJ) <= STREAM_IYB ) is_disk_XBC = 1;
          #endif
        }
      }
 
      #include PR_BC_SPECIAL //special BC check - Brandon
      #endif

      if(BCtype==-1) my_err("wrong GC in loop_2\n");

      #ifdef SPECIAL_BC_CHECK
      if(mpi_isitBC(BCtype)==0 && is_disk_XBC == 0 && is_disk_YBC == 0 && is_disk_ZBC == 0) //this border exchanged through MPI
      #else
      if(mpi_isitBC(BCtype)==0) //this border exchanged through MPI
      #endif
	{
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	}
      else //need for real BC
	{
	  ldouble uval[NV],pval[NV];

	  set_bc_core(ix,iy,iz,t,uval,pval,ifinit,BCtype);

	  for(iv=0;iv<NV;iv++)
	    {
	      set_u(u,iv,ix,iy,iz,uval[iv]);
	      set_u(p,iv,ix,iy,iz,pval[iv]);	      
	    }
	}
    }

  //treating the corners of tiles

  int ix,iy,iz,iv;

#ifdef MPI4CORNERS

  /*****************************************************************/
  /* now calculate conserved in corners from exchanged primitives */
  /*****************************************************************/  
  struct geometry geom;

  if(TNY==1 && TNZ>1) //2D
    {
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCLO)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iz=-NG;iz<0;iz++)
	      for(iy=0;iy<NY;iy++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iz=-NG;iz<0;iz++)
	       for(iy=0;iy<NY;iy++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iz=NZ;iz<NZ+NG;iz++)
	       for(iy=0;iy<NY;iy++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iz=NZ;iz<NZ+NG;iz++)
	       for(iy=0;iy<NY;iy++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
    }
  
  if(TNZ==1 && TNY>1) //2D
    {
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iy=-NG;iy<0;iy++)
	      for(iz=0;iz<NZ;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iy=-NG;iy<0;iy++)
	      for(iz=0;iz<NZ;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iy=NY;iy<NY+NG;iy++)
	      for(iz=0;iz<NZ;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iy=NY;iy<NY+NG;iy++)
	      for(iz=0;iz<NZ;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
    }
  
  if(TNZ>1 && TNY>1) //full 3d
    {
      //elongated blocks first
      //along z
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iy=-NG;iy<0;iy++)
	      for(iz=0;iz<NZ;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iy=-NG;iy<0;iy++)
	      for(iz=0;iz<NZ;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iy=NY;iy<NY+NG;iy++)
	      for(iz=0;iz<NZ;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iy=NY;iy<NY+NG;iy++)
	      for(iz=0;iz<NZ;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      //along y
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCLO)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iz=-NG;iz<0;iz++)
	      for(iy=0;iy<NY;iy++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCLO)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iz=-NG;iz<0;iz++)
	      for(iy=0;iy<NY;iy++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCHI)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iz=NZ;iz<NZ+NG;iz++)
	      for(iy=0;iy<NY;iy++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCHI)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iz=NZ;iz<NZ+NG;iz++)
	      for(iy=0;iy<NY;iy++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      //along x
      if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0) 
	{
	  for(iy=-NG;iy<0;iy++)
	    for(iz=-NG;iz<0;iz++)
	      for(ix=0;ix<NX;ix++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0) 
	{
	  for(iy=NY;iy<NY+NG;iy++)
	    for(iz=-NG;iz<0;iz++)
	      for(ix=0;ix<NX;ix++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0) 
	{
	  for(iy=-NG;iy<0;iy++)
	    for(iz=NZ;iz<NZ+NG;iz++)
	      for(ix=0;ix<NX;ix++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0) 
	{
	  for(iy=NY;iy<NY+NG;iy++)
	    for(iz=NZ;iz<NZ+NG;iz++)
	      for(ix=0;ix<NX;ix++)
		{	    
		  fill_geometry(ix,iy,iz,&geom);
		  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      //now cubic corners corners
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iy=-NG;iy<0;iy++)
	      for(iz=-NG;iz<0;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom); p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iy=-NG;iy<0;iy++)
	      for(iz=-NG;iz<0;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom); p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iy=NY;iy<NY+NG;iy++)
	      for(iz=-NG;iz<0;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom); p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iy=NY;iy<NY+NG;iy++)
	      for(iz=-NG;iz<0;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom); p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iy=-NG;iy<0;iy++)
	      for(iz=NZ;iz<NZ+NG;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom); p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iy=-NG;iy<0;iy++)
	      for(iz=NZ;iz<NZ+NG;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom); p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0) 
	{
	  for(ix=-NG;ix<0;ix++)
	    for(iy=NY;iy<NY+NG;iy++)
	      for(iz=NZ;iz<NZ+NG;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom); p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0) 
	{
	  for(ix=NX;ix<NX+NG;ix++)
	    for(iy=NY;iy<NY+NG;iy++)
	      for(iz=NZ;iz<NZ+NG;iz++)
		{	    
		  fill_geometry(ix,iy,iz,&geom); p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
 


    }


  //corners of the whole domain are sometimes not real BC
  //so need to fill them with something  
  /**************************************/
  // Brandon - Note that here we use mpi_isitBC_forcorners since these values correspond to 
  // boundaries where we would NOT have passed values via mpi    
  // since here mpi_isitBC(XBC)=1 but mpi_isitBC(YBC)=0 
  // at upper and lower poles for transmitting boundary
  /**************************************/
  int xlim,ylim,zlim;
  int lim,i,j,k;
  int ix1,iy1,iz1;
  int ix2,iy2,iz2;

  if(TNZ==1 && TNY>1) //2d
    {
      iz=0;

      //total corners, filling one cell deep surfaces

      //bottom-left corner
      if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCLO)==1)
	{

	  for(i=0;i<NG-1;i++)
	    {
	      PLOOP(iv)
	      {
		set_u(p,iv,-NG+i,-1,iz,get_u(p,iv,-NG+i,0,iz));
		set_u(p,iv,-1,-NG+i,iz,get_u(p,iv,0,-NG+i,iz));
	      }

	      fill_geometry(-NG+i,-1,iz,&geom);
	      p2u(&get_u(p,0,-NG+i,-1,iz),&get_u(u,0,-NG+i,-1,iz),&geom);
	      fill_geometry(-1,-NG+i,iz,&geom);
	      p2u(&get_u(p,0,-1,-NG+i,iz),&get_u(u,0,-1,-NG+i,iz),&geom);
	    }
      
	  //averaging <(-1,0),(0,-1)> -> (-1,-1)
	  ix1=-1;iy1=0;ix2=0;iy2=-1;
          #if defined(PERIODIC_YBC) && !defined(MPI) //periodic boundary conditions more tricky in MPI, in this case use simpler total corners
	  ix1=ix2=-1;iy1=iy2=NY-1;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=NX-1;iy1=iy2=-1;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,-1,-1,iz,.5*(get_u(p,iv,ix1,iy1,iz)+get_u(p,iv,ix2,iy2,iz)));

	  fill_geometry(-1,-1,iz,&geom);
	  p2u(&get_u(p,0,-1,-1,iz),&get_u(u,0,-1,-1,iz),&geom);

	  //averaging <(-2,-1),(-1,-2)> -> (-2,-2)
	  ix1=-2;iy1=-1;ix2=-1;iy2=-2;
          #if defined(PERIODIC_YBC) && !defined(MPI)
	  ix1=ix2=-2;iy1=iy2=NY-2;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=NX-2;iy1=iy2=-2;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,-2,-2,iz,.5*(get_u(p,iv,ix1,iy1,iz)+get_u(p,iv,ix2,iy2,iz)));

	  fill_geometry(-2,-2,iz,&geom);
	  p2u(&get_u(p,0,-2,-2,iz),&get_u(u,0,-2,-2,iz),&geom);
      

	}
  
      //top-left corner
      if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCHI)==1)
	{
 
	  for(i=0;i<NG-1;i++)
	    {
	      PLOOP(iv)
	      {
		set_u(p,iv,-NG+i,NY,iz,get_u(p,iv,-NG+i,NY-1,iz));
		set_u(p,iv,-1,NY+i+1,iz,get_u(p,iv,0,NY+i+1,iz));
	      }

	      fill_geometry(-NG+i,NY,iz,&geom);
	      p2u(&get_u(p,0,-NG+i,NY,iz),&get_u(u,0,-NG+i,NY,iz),&geom);
	      fill_geometry(-1,NY+i+1,iz,&geom);
	      p2u(&get_u(p,0,-1,NY+i+1,iz),&get_u(u,0,-1,NY+i+1,iz),&geom);
	    }

	  //averaging <(-1,NY-1),(0,NY)> -> (-1,NY)
	  ix1=-1;iy1=NY-1;ix2=0;iy2=NY;
	  #if defined(PERIODIC_YBC) && !defined(MPI)
	  ix1=ix2=-1;iy1=iy2=0;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=NX-1;iy1=iy2=NY;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,-1,NY,iz,.5*(get_u(p,iv,ix1,iy1,iz)+get_u(p,iv,ix2,iy2,iz)));


	  fill_geometry(-1,NY,iz,&geom);
	  p2u(&get_u(p,0,-1,NY,iz),&get_u(u,0,-1,NY,iz),&geom);

	  //averaging <(-2,NY),(-1,NY+1)> -> (-2,NY+1)
	  ix1=-2;iy1=NY;ix2=-1;iy2=NY+1;
          #if defined(PERIODIC_YBC) && !defined(MPI)
	  ix1=ix2=-2;iy1=iy2=1;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=NX-2;iy1=iy2=NY+1;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,-2,NY+1,iz,.5*(get_u(p,iv,ix1,iy1,iz)+get_u(p,iv,ix2,iy2,iz)));

	  fill_geometry(-2,NY+1,iz,&geom);
	  p2u(&get_u(p,0,-2,NY+1,iz),&get_u(u,0,-2,NY+1,iz),&geom);
	}

  //bottom-right corner
      if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCLO)==1)
	{
	  for(i=0;i<NG-1;i++)
	    {
	      PLOOP(iv)
	      {
		set_u(p,iv,NX+i+1,-1,iz,get_u(p,iv,NX+i+1,0,iz));
		set_u(p,iv,NX,-NG+i,iz,get_u(p,iv,NX-1,-NG+i,iz));
	      }

	      fill_geometry(NX+i+1,-1,iz,&geom);
	      p2u(&get_u(p,0,NX+i+1,-1,iz),&get_u(u,0,NX+i+1,-1,iz),&geom);
	      fill_geometry(NX,-NG+i,iz,&geom);
	      p2u(&get_u(p,0,NX,-NG+i,iz),&get_u(u,0,NX,-NG+i,iz),&geom);
	    }

	  //averaging <(NX-1,-1),(NX,0)> -> (NX,-1)
	  ix1=NX-1;iy1=-1;ix2=NX;iy2=0;
          #if defined(PERIODIC_YBC) && !defined(MPI)
	  ix1=ix2=NX;iy1=iy2=NY-1;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=0;iy1=iy2=-1;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,NX,-1,iz,.5*(get_u(p,iv,ix1,iy1,iz)+get_u(p,iv,ix2,iy2,iz)));
	 
	  fill_geometry(NX,-1,iz,&geom);
	  p2u(&get_u(p,0,NX,-1,iz),&get_u(u,0,NX,-1,iz),&geom);

	  //averaging <(NX,-2),(NX+1,-1)> -> (NX+1,-2)
	  ix1=NX;iy1=-2;ix2=NX+1;iy2=-1;
	  #if defined(PERIODIC_YBC) && !defined(MPI)
	  ix1=ix2=NX+1;iy1=iy2=NY-2;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=1;iy1=iy2=-2;
          #endif

	   PLOOP(iv)
	    set_u(p,iv,NX+1,-2,iz,.5*(get_u(p,iv,ix1,iy1,iz)+get_u(p,iv,ix2,iy2,iz)));

	  fill_geometry(NX+1,-2,iz,&geom);
	  p2u(&get_u(p,0,NX+1,-2,iz),&get_u(u,0,NX+1,-2,iz),&geom);
	}

      //top-right corner
      if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCHI)==1)
	{
	  for(i=0;i<NG-1;i++)
	    {
	      PLOOP(iv)
	      {
		set_u(p,iv,NX+i+1,NY,iz,get_u(p,iv,NX+i+1,NY-1,iz));
		set_u(p,iv,NX,NY+i+1,iz,get_u(p,iv,NX-1,NY+i+1,iz));
	      }
	      fill_geometry(NX+i+1,NY,iz,&geom);
	      p2u(&get_u(p,0,NX+i+1,NY,iz),&get_u(u,0,NX+i+1,NY,iz),&geom);
	      fill_geometry(NX,NY+i+1,iz,&geom);
	      p2u(&get_u(p,0,NX,NY+i+1,iz),&get_u(u,0,NX,NY+i+1,iz),&geom);
	    }

	  //averaging <(NX-1,NY),(NX,NY-1)> -> (NX,NY)
	  ix1=NX-1;iy1=NY;ix2=NX;iy2=NY-1;
          #if defined(PERIODIC_YBC) && !defined(MPI)
	  ix1=ix2=NX;iy1=iy2=0;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=0;iy1=iy2=NY;
          #endif
	  
	  PLOOP(iv)
	    set_u(p,iv,NX,NY,iz,.5*(get_u(p,iv,ix1,iy1,iz)+get_u(p,iv,ix2,iy2,iz)));

	  fill_geometry(NX,NY,iz,&geom);
	  p2u(&get_u(p,0,NX,NY,iz),&get_u(u,0,NX,NY,iz),&geom);
	  
	  //averaging <(NX,NY+1),(NX+1,NY)> -> (NX+1,NY+1)
	  ix1=NX;iy1=NY+1;ix2=NX+1;iy2=NY;
          #if defined(PERIODIC_YBC) && !defined(MPI)
	  ix1=ix2=NX+1;iy1=iy2=1;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=1;iy1=iy2=NY+1;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,NX+1,NY+1,iz,.5*(get_u(p,iv,ix1,iy1,iz)+get_u(p,iv,ix2,iy2,iz)));

	  fill_geometry(NX+1,NY+1,iz,&geom);
	  p2u(&get_u(p,0,NX+1,NY+1,iz),&get_u(u,0,NX+1,NY+1,iz),&geom);
	}
    }

  if(TNY==1 && TNZ>1) //2d
    {
      iy=0;

      //total corners, filling one cell deep surfaces

      //bottom-left corner
      if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
#ifdef SHEARINGBOX //in shearing box copy into corners along z
	  PLOOP(iv)
	      {
		set_u(p,iv,-1,0,-1,get_u(p,iv,-1,0,0));
		set_u(p,iv,-2,0,-1,get_u(p,iv,-2,0,0));
		set_u(p,iv,-1,0,-2,get_u(p,iv,-1,0,0));
		set_u(p,iv,-2,0,-2,get_u(p,iv,-2,0,0));
	      }
	  fill_geometry(-1,0,-1,&geom);p2u(&get_u(p,0,-1,0,-1),&get_u(u,0,-1,0,-1),&geom);
	  fill_geometry(-2,0,-1,&geom);p2u(&get_u(p,0,-1,0,-1),&get_u(u,0,-1,0,-1),&geom);
	  fill_geometry(-1,0,-2,&geom);p2u(&get_u(p,0,-1,0,-2),&get_u(u,0,-1,0,-2),&geom);
	  fill_geometry(-2,0,-2,&geom);p2u(&get_u(p,0,-2,0,-2),&get_u(u,0,-2,0,-2),&geom);

#else

	  for(i=0;i<NG-1;i++)
	    {

	      PLOOP(iv)
	      {
		set_u(p,iv,-NG+i,iy,-1,get_u(p,iv,-NG+i,iy,0));
		set_u(p,iv,-1,iy,-NG+i,get_u(p,iv,0,iy,-NG+i));
	      }

	      fill_geometry(-NG+i,iy,-1,&geom);
	      p2u(&get_u(p,0,-NG+i,iy,-1),&get_u(u,0,-NG+i,iy,-1),&geom);
	      fill_geometry(-1,iy,-NG+i,&geom);
	      p2u(&get_u(p,0,-1,iy,-NG+i),&get_u(u,0,-1,iy,-NG+i),&geom);
	    }
      
	  //averaging <(-1,0),(0,-1)> -> (-1,-1)
	  ix1=-1;iz1=0;ix2=0;iz2=-1;
          #if defined(PERIODIC_ZBC) && !defined(MPI) //periodic boundary conditions more tricky in MPI, in this case use simpler total corners
	  ix1=ix2=-1;iz1=iz2=NZ-1;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=NX-1;iz1=iz2=-1;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,-1,iy,-1,.5*(get_u(p,iv,ix1,iy,iz1)+get_u(p,iv,ix2,iy,iz2)));

	  fill_geometry(-1,iy,-1,&geom);
	  p2u(&get_u(p,0,-1,iy,-1),&get_u(u,0,-1,iy,-1),&geom);

	  //averaging <(-2,-1),(-1,-2)> -> (-2,-2)
	  ix1=-2;iz1=-1;ix2=-1;iz2=-2;
          #if defined(PERIODIC_ZBC) && !defined(MPI)
	  ix1=ix2=-2;iz1=iz2=NZ-2;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=NX-2;iz1=iz2=-2;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,-2,iy,-2,.5*(get_u(p,iv,ix1,iy,iz1)+get_u(p,iv,ix2,iy,iz2)));

	  fill_geometry(-2,iy,-2,&geom);
	  p2u(&get_u(p,0,-2,iy,-2),&get_u(u,0,-2,iy,-2),&geom);
 #endif     

	}
  
      //top-left corner
      if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{

#ifdef SHEARINGBOX //in shearing box copy into corners along z
	  PLOOP(iv)
	      {
		set_u(p,iv,-1,0,NZ,get_u(p,iv,-1,0,NZ-1));
		set_u(p,iv,-2,0,NZ,get_u(p,iv,-2,0,NZ-1));
		set_u(p,iv,-1,0,NZ+1,get_u(p,iv,-1,0,NZ-1));
		set_u(p,iv,-2,0,NZ+1,get_u(p,iv,-2,0,NZ-1));
	      }
	  fill_geometry(-1,0,NZ,&geom);p2u(&get_u(p,0,-1,0,NZ),&get_u(u,0,-1,0,NZ),&geom);
	  fill_geometry(-2,0,NZ,&geom);p2u(&get_u(p,0,-1,0,NZ),&get_u(u,0,-1,0,NZ),&geom);
	  fill_geometry(-1,0,NZ+1,&geom);p2u(&get_u(p,0,-1,0,NZ+1),&get_u(u,0,-1,0,NZ+1),&geom);
	  fill_geometry(-2,0,NZ+1,&geom);p2u(&get_u(p,0,-2,0,NZ+1),&get_u(u,0,-2,0,NZ+1),&geom);
#else

 
	  for(i=0;i<NG-1;i++)
	    {
	      PLOOP(iv)
	      {
		set_u(p,iv,-NG+i,iy,NZ,get_u(p,iv,-NG+i,iy,NZ-1));
		set_u(p,iv,-1,iy,NZ+i+1,get_u(p,iv,0,iy,NZ+i+1));
	      }

	      fill_geometry(-NG+i,iy,NZ,&geom);
	      p2u(&get_u(p,0,-NG+i,iy,NZ),&get_u(u,0,-NG+i,iy,NZ),&geom);
	      fill_geometry(-1,iy,NZ+i+1,&geom);
	      p2u(&get_u(p,0,-1,iy,NZ+i+1),&get_u(u,0,-1,iy,NZ+i+1),&geom);
	    }

	  //averaging <(-1,NY-1),(0,NY)> -> (-1,NY)
	  ix1=-1;iz1=NZ-1;ix2=0;iz2=NZ;
	  #if defined(PERIODIC_ZBC) && !defined(MPI)
	  ix1=ix2=-1;iz1=iz2=0;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=NX-1;iz1=iz2=NZ;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,-1,iy,NZ,.5*(get_u(p,iv,ix1,iy,iz1)+get_u(p,iv,ix2,iy,iz2)));


	  fill_geometry(-1,iy,NZ,&geom);
	  p2u(&get_u(p,0,-1,iy,NZ),&get_u(u,0,-1,iy,NZ),&geom);

	  //averaging <(-2,NY),(-1,NY+1)> -> (-2,NY+1)
	  ix1=-2;iz1=NZ;ix2=-1;iz2=NZ+1;
          #if defined(PERIODIC_ZBC) && !defined(MPI)
	  ix1=ix2=-2;iz1=iz2=1;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=NX-2;iz1=iz2=NZ+1;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,-2,iy,NZ+1,.5*(get_u(p,iv,ix1,iy,iz1)+get_u(p,iv,ix2,iy,iz2)));

	  fill_geometry(-2,iy,NZ+1,&geom);
	  p2u(&get_u(p,0,-2,iy,NZ+1),&get_u(u,0,-2,iy,NZ+1),&geom);
	  #endif
	}

  //bottom-right corner
      if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
	  #ifdef SHEARINGBOX //in shearing box copy into corners along z
	  PLOOP(iv)
	      {
		set_u(p,iv,NX,0,-1,get_u(p,iv,NX,0,0));
		set_u(p,iv,NX+1,0,-1,get_u(p,iv,NX+1,0,0));
		set_u(p,iv,NX,0,-2,get_u(p,iv,NX,0,0));
		set_u(p,iv,NX+1,0,-2,get_u(p,iv,NX+1,0,0));
	      }
	  fill_geometry(NX,0,-1,&geom);p2u(&get_u(p,0,NX,0,-1),&get_u(u,0,NX,0,-1),&geom);
	  fill_geometry(NX+1,0,-1,&geom);p2u(&get_u(p,0,NX+1,0,-1),&get_u(u,0,NX+1,0,-1),&geom);
	  fill_geometry(NX,0,-2,&geom);p2u(&get_u(p,0,NX,0,-2),&get_u(u,0,NX,0,-2),&geom);
	  fill_geometry(NX+1,0,-2,&geom);p2u(&get_u(p,0,NX+1,0,-2),&get_u(u,0,NX+1,0,-2),&geom);
	      	
#else

	  for(i=0;i<NG-1;i++)
	    {
	      PLOOP(iv)
	      {
		set_u(p,iv,NX+i+1,iy,-1,get_u(p,iv,NX+i+1,iy,0));
		set_u(p,iv,NX,iy,-NG+i,get_u(p,iv,NX-1,iy,-NG+i));
	      }

	      fill_geometry(NX+i+1,iy,-1,&geom);
	      p2u(&get_u(p,0,NX+i+1,iy,-1),&get_u(u,0,NX+i+1,iy,-1),&geom);
	      fill_geometry(NX,iy,-NG+i,&geom);
	      p2u(&get_u(p,0,NX,iy,-NG+i),&get_u(u,0,NX,iy,-NG+i),&geom);
	    }

	  //averaging <(NX-1,-1),(NX,0)> -> (NX,-1)
	  ix1=NX-1;iz1=-1;ix2=NX;iz2=0;
          #if defined(PERIODIC_ZBC) && !defined(MPI)
	  ix1=ix2=NX;iz1=iz2=NZ-1;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=0;iz1=iz2=-1;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,NX,iy,-1,.5*(get_u(p,iv,ix1,iy,iz1)+get_u(p,iv,ix2,iy,iz2)));
	 
	  fill_geometry(NX,iy,-1,&geom);
	  p2u(&get_u(p,0,NX,iy,-1),&get_u(u,0,NX,iy,-1),&geom);

	  //averaging <(NX,-2),(NX+1,-1)> -> (NX+1,-2)
	  ix1=NX;iz1=-2;ix2=NX+1;iz2=-1;
	  #if defined(PERIODIC_ZBC) && !defined(MPI)
	  ix1=ix2=NX+1;iz1=iz2=NZ-2;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=1;iz1=iz2=-2;
          #endif

	   PLOOP(iv)
	     set_u(p,iv,NX+1,iy,-2,.5*(get_u(p,iv,ix1,iy,iz1)+get_u(p,iv,ix2,iy,iz2)));

	  fill_geometry(NX+1,iy,-2,&geom);
	  p2u(&get_u(p,0,NX+1,iy,-2),&get_u(u,0,NX+1,iy,-2),&geom);
#endif
	}

      //top-right corner
      if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{
	  #ifdef SHEARINGBOX //in shearing box copy into corners along z
	  PLOOP(iv)
	      {
		set_u(p,iv,NX,0,NZ,get_u(p,iv,NX,0,NZ-1));
		set_u(p,iv,NX+1,0,NZ,get_u(p,iv,NX+1,0,NZ-1));
		set_u(p,iv,NX,0,NZ+1,get_u(p,iv,NX,0,NZ-1));
		set_u(p,iv,NX+1,0,NZ+1,get_u(p,iv,NX+1,0,NZ-1));
	      }
	  fill_geometry(NX,0,NZ,&geom);p2u(&get_u(p,0,NX,0,NZ),&get_u(u,0,NX,0,NZ),&geom);
	  fill_geometry(NX+1,0,NZ,&geom);p2u(&get_u(p,0,NX+1,0,NZ),&get_u(u,0,NX+1,0,NZ),&geom);
	  fill_geometry(NX,0,NZ+1,&geom);p2u(&get_u(p,0,NX,0,NZ+1),&get_u(u,0,NX,0,NZ+1),&geom);
	  fill_geometry(NX+1,0,NZ+1,&geom);p2u(&get_u(p,0,NX+1,0,NZ+1),&get_u(u,0,NX+1,0,NZ+1),&geom);
#else

	  for(i=0;i<NG-1;i++)
	    {
	      PLOOP(iv)
	      {
		set_u(p,iv,NX+i+1,iy,NZ,get_u(p,iv,NX+i+1,iy,NZ-1));
		set_u(p,iv,NX,iy,NZ+i+1,get_u(p,iv,NX-1,iy,NZ+i+1));
	      }

	      fill_geometry(NX+i+1,iy,NZ,&geom);
	      p2u(&get_u(p,0,NX+i+1,iy,NZ),&get_u(u,0,NX+i+1,iy,NZ),&geom);
	      fill_geometry(NX,iy,NZ+i+1,&geom);
	      p2u(&get_u(p,0,NX,iy,NZ+i+1),&get_u(u,0,NX,iy,NZ+i+1),&geom);
	    }

	  //averaging <(NX-1,NY),(NX,NY-1)> -> (NX,NY)
	  ix1=NX-1;iz1=NZ;ix2=NX;iz2=NZ-1;
          #if defined(PERIODIC_ZBC) && !defined(MPI)
	  ix1=ix2=NX;iz1=iz2=0;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=0;iz1=iz2=NZ;
          #endif
	  
	  PLOOP(iv)
	    set_u(p,iv,NX,iy,NZ,.5*(get_u(p,iv,ix1,iy,iz1)+get_u(p,iv,ix2,iy,iz2)));

	  fill_geometry(NX,iy,NZ,&geom);
	  p2u(&get_u(p,0,NX,iy,NZ),&get_u(u,0,NX,iy,NZ),&geom);
	  
	  //averaging <(NX,NY+1),(NX+1,NY)> -> (NX+1,NY+1)
	  ix1=NX;iz1=NZ+1;ix2=NX+1;iz2=NZ;
          #if defined(PERIODIC_ZBC) && !defined(MPI)
	  ix1=ix2=NX+1;iz1=iz2=1;
          #endif
	  #if defined(PERIODIC_XBC) && !defined(MPI)
	  ix1=ix2=1;iz1=iz2=NZ+1;
          #endif

	  PLOOP(iv)
	    set_u(p,iv,NX+1,iy,NZ+1,.5*(get_u(p,iv,ix1,iy,iz1)+get_u(p,iv,ix2,iy,iz2)));

	  fill_geometry(NX+1,iy,NZ+1,&geom);
	  p2u(&get_u(p,0,NX+1,iy,NZ+1),&get_u(u,0,NX+1,iy,NZ+1),&geom);
	  #endif
	}
    }

  /**************************************/
  if(TNZ>1 && TNY>1) //full 3d
    {
      #ifdef SHEARINGBOX 
      //in shearing box copy into corners along z for top and bottom 
      //but first apply periodic on horizontal edges

      //elongated along z - periodic
      for(iz=0;iz<NZ;iz++) 
	for(iy=-NG;iy<0;iy++)
	  for(ix=-NG;ix<0;ix++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy+NY,iz));
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
      for(iz=0;iz<NZ;iz++) 
	for(iy=NY;iy<NY+NG;iy++)
	  for(ix=-NG;ix<0;ix++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy-NY,iz));
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
      for(iz=0;iz<NZ;iz++) 
	for(iy=-NG;iy<0;iy++)
	  for(ix=NX;ix<NX+NG;ix++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy+NY,iz));
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
      for(iz=0;iz<NZ;iz++) 
	for(iy=NY;iy<NY+NG;iy++)
	  for(ix=NX;ix<NX+NG;ix++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy-NY,iz));
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
      //elongated along x - copying along z
       for(ix=0;ix<NX;ix++) 
	for(iy=-NG;iy<0;iy++)
	  for(iz=-NG;iz<0;iz++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,0));
	      if(get_u(p,VZ,ix,iy,iz)>0.) set_u(p,VZ,ix,iy,iz,0.);
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
      for(ix=0;ix<NX;ix++) 
	for(iy=NY;iy<NY+NG;iy++)
	  for(iz=-NG;iz<0;iz++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,0));
	      if(get_u(p,VZ,ix,iy,iz)>0.) set_u(p,VZ,ix,iy,iz,0.);
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
       for(ix=0;ix<NX;ix++) 
	for(iy=-NG;iy<0;iy++)
	  for(iz=NZ;iz<NZ+NG;iz++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,NZ-1));
	      if(get_u(p,VZ,ix,iy,iz)<0.) set_u(p,VZ,ix,iy,iz,0.);
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
      for(ix=0;ix<NX;ix++) 
	for(iy=NY;iy<NY+NG;iy++)
	  for(iz=NZ;iz<NZ+NG;iz++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,NZ-1));
	      if(get_u(p,VZ,ix,iy,iz)<0.) set_u(p,VZ,ix,iy,iz,0.);
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
      //elongated along y - copying along z
       for(iy=0;iy<NY;iy++) 
	for(ix=-NG;ix<0;ix++)
	  for(iz=-NG;iz<0;iz++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,0));
	      if(get_u(p,VZ,ix,iy,iz)>0.) set_u(p,VZ,ix,iy,iz,0.);
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
       for(iy=0;iy<NY;iy++) 
	 for(ix=NX;ix<NX+NG;ix++)
	  for(iz=-NG;iz<0;iz++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,0));
	      if(get_u(p,VZ,ix,iy,iz)>0.) set_u(p,VZ,ix,iy,iz,0.);
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
       for(iy=0;iy<NY;iy++) 
	for(ix=-NG;ix<0;ix++)
	  for(iz=NZ;iz<NZ+NG;iz++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,NZ-1));
	      if(get_u(p,VZ,ix,iy,iz)<0.) set_u(p,VZ,ix,iy,iz,0.);
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
       for(iy=0;iy<NY;iy++) 
	for(ix=NX;ix<NX+NG;ix++)
	  for(iz=NZ;iz<NZ+NG;iz++)
	    {
	      PLOOP(iv)
		set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,NZ-1));
	      if(get_u(p,VZ,ix,iy,iz)<0.) set_u(p,VZ,ix,iy,iz,0.);
	      fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
       //total lower corners
       for(ix=-NG;ix<0;ix++) 
	 for(iy=-NG;iy<0;iy++)
	   for(iz=-NG;iz<0;iz++)
	     {
	       PLOOP(iv)
		 set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,0));
	      if(get_u(p,VZ,ix,iy,iz)>0.) set_u(p,VZ,ix,iy,iz,0.);
	       fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	     }
       for(ix=-NG;ix<0;ix++) 
	 for(iy=NY;iy<NY+NG;iy++)
	   for(iz=-NG;iz<0;iz++)
	     {
	       PLOOP(iv)
		 set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,0));
	      if(get_u(p,VZ,ix,iy,iz)>0.) set_u(p,VZ,ix,iy,iz,0.);
	       fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	     }
       for(ix=NX;ix<NX+NG;ix++) 
	 for(iy=-NG;iy<0;iy++)
	   for(iz=-NG;iz<0;iz++)
	     {
	       PLOOP(iv)
		 set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,0));
	      if(get_u(p,VZ,ix,iy,iz)>0.) set_u(p,VZ,ix,iy,iz,0.);
	       fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	     }
       for(ix=NX;ix<NX+NG;ix++) 
	 for(iy=NY;iy<NY+NG;iy++)
	   for(iz=-NG;iz<0;iz++)
	     {
	       PLOOP(iv)
		 set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,0));
	      if(get_u(p,VZ,ix,iy,iz)>0.) set_u(p,VZ,ix,iy,iz,0.);
	       fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	     }
       //total upper corners
       for(ix=-NG;ix<0;ix++) 
	 for(iy=-NG;iy<0;iy++)
	   for(iz=NZ;iz<NZ+NG;iz++)
	     {
	       PLOOP(iv)
		 set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,NZ-1));
	       if(get_u(p,VZ,ix,iy,iz)<0.) set_u(p,VZ,ix,iy,iz,0.);
	       fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	     }
       for(ix=-NG;ix<0;ix++) 
	 for(iy=NY;iy<NY+NG;iy++)
	   for(iz=NZ;iz<NZ+NG;iz++)
	     {
	       PLOOP(iv)
		 set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,NZ-1));
	       if(get_u(p,VZ,ix,iy,iz)<0.) set_u(p,VZ,ix,iy,iz,0.);
	       fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	     }
       for(ix=NX;ix<NX+NG;ix++) 
	 for(iy=-NG;iy<0;iy++)
	   for(iz=NZ;iz<NZ+NG;iz++)
	     {
	       PLOOP(iv)
		 set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,NZ-1));
	       if(get_u(p,VZ,ix,iy,iz)<0.) set_u(p,VZ,ix,iy,iz,0.);
	       fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	     }
       for(ix=NX;ix<NX+NG;ix++) 
	 for(iy=NY;iy<NY+NG;iy++)
	   for(iz=NZ;iz<NZ+NG;iz++)
	     {
	       PLOOP(iv)
		 set_u(p,iv,ix,iy,iz,get_u(p,iv,ix,iy,NZ-1));
	       if(get_u(p,VZ,ix,iy,iz)<0.) set_u(p,VZ,ix,iy,iz,0.);
	       fill_geometry(ix,iy,iz,&geom);p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	     }
	
#else //not SHEARINGBOX, regular

      //elongated corners along z
      //filling one cell deep surfaces, and averaging diagonally
      if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCLO)==1)
	{
	  for(iz=-NG;iz<NZ+NG;iz++) { //in the total total corners it fills crap but overwritten below!
	    for(i=0;i<NG-1;i++) {
	      PLOOP(iv) {
		set_u(p,iv,-NG+i,-1,iz,get_u(p,iv,-NG+i,0,iz));
		set_u(p,iv,-1,-NG+i,iz,get_u(p,iv,0,-NG+i,iz)); }
	      fill_geometry(-NG+i,-1,iz,&geom);  p2u(&get_u(p,0,-NG+i,-1,iz),&get_u(u,0,-NG+i,-1,iz),&geom);
	      fill_geometry(-1,-NG+i,iz,&geom);  p2u(&get_u(p,0,-1,-NG+i,iz),&get_u(u,0,-1,-NG+i,iz),&geom);
	    }
      
	    //averaging <(-1,0),(0,-1)> -> (-1,-1)
	    PLOOP(iv)
	      set_u(p,iv,-1,-1,iz,.5*(get_u(p,iv,-1,0,iz)+get_u(p,iv,0,-1,iz)));
	    fill_geometry(-1,-1,iz,&geom);  p2u(&get_u(p,0,-1,-1,iz),&get_u(u,0,-1,-1,iz),&geom);

	    //averaging <(-2,-1),(-1,-2)> -> (-2,-2)
	    PLOOP(iv)
	      set_u(p,iv,-2,-2,iz,.5*(get_u(p,iv,-2,-1,iz)+get_u(p,iv,-1,-2,iz)));
	    fill_geometry(-2,-2,iz,&geom);   p2u(&get_u(p,0,-2,-2,iz),&get_u(u,0,-2,-2,iz),&geom);
	  }
	}

      if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCHI)==1)
	{
	  for(iz=-NG;iz<NZ+NG;iz++) {
	    for(i=0;i<NG-1;i++) {
	      PLOOP(iv) {
		set_u(p,iv,-NG+i,NY,iz,get_u(p,iv,-NG+i,NY-1,iz));
		set_u(p,iv,-1,NY+i+1,iz,get_u(p,iv,0,NY+i+1,iz)); }
	      fill_geometry(-NG+i,NY,iz,&geom);  p2u(&get_u(p,0,-NG+i,NY,iz),&get_u(u,0,-NG+i,NY,iz),&geom);
	      fill_geometry(-1,NY+i+1,iz,&geom);  p2u(&get_u(p,0,-1,NY+i+1,iz),&get_u(u,0,-1,NY+i+1,iz),&geom);
	    }
	    
	    PLOOP(iv)
	      set_u(p,iv,-1,NY,iz,.5*(get_u(p,iv,-1,NY-1,iz)+get_u(p,iv,0,NY,iz)));
	    fill_geometry(-1,NY,iz,&geom);  p2u(&get_u(p,0,-1,NY,iz),&get_u(u,0,-1,NY,iz),&geom);

	    PLOOP(iv)
	      set_u(p,iv,-2,NY+1,iz,.5*(get_u(p,iv,-2,NY,iz)+get_u(p,iv,-1,NY+1,iz)));
	    fill_geometry(-2,NY+1,iz,&geom);   p2u(&get_u(p,0,-2,NY+1,iz),&get_u(u,0,-2,NY+1,iz),&geom);
	  }
	}

      if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCLO)==1)
	{
	  for(iz=-NG;iz<NZ+NG;iz++) {
	    for(i=0;i<NG-1;i++) {
	      PLOOP(iv) {
		set_u(p,iv,NX,-NG+i,iz,get_u(p,iv,NX-1,-NG+i,iz));
		set_u(p,iv,NX+i+1,-1,iz,get_u(p,iv,NX+i+1,0,iz)); }
	      fill_geometry(NX+i+1,-1,iz,&geom);    p2u(&get_u(p,0,NX+i+1,-1,iz),&get_u(u,0,NX+i+1,-1,iz),&geom);
	      fill_geometry(NX,-NG+i,iz,&geom);    p2u(&get_u(p,0,NX,-NG+i,iz),&get_u(u,0,NX,-NG+i,iz),&geom);
	    }

	    PLOOP(iv)
	      set_u(p,iv,NX,-1,iz,.5*(get_u(p,iv,NX-1,-1,iz)+get_u(p,iv,NX,0,iz)));
	    fill_geometry(NX,-1,iz,&geom); p2u(&get_u(p,0,NX,-1,iz),&get_u(u,0,NX,-1,iz),&geom);

	    PLOOP(iv)
	      set_u(p,iv,NX+1,-2,iz,.5*(get_u(p,iv,NX,-2,iz)+get_u(p,iv,NX+1,-1,iz)));
	    fill_geometry(NX+1,-2,iz,&geom);  p2u(&get_u(p,0,NX+1,-2,iz),&get_u(u,0,NX+1,-2,iz),&geom);
	  }
	}

       if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCHI)==1)
	{
	  for(iz=-NG;iz<NZ+NG;iz++) {
	    for(i=0;i<NG-1;i++)
	      {
		PLOOP(iv) {
		  set_u(p,iv,NX+i+1,NY,iz,get_u(p,iv,NX+i+1,NY-1,iz));
		  set_u(p,iv,NX,NY+i+1,iz,get_u(p,iv,NX-1,NY+i+1,iz));
		}
		fill_geometry(NX+i+1,NY,iz,&geom);		p2u(&get_u(p,0,NX+i+1,NY,iz),&get_u(u,0,NX+i+1,NY,iz),&geom);
		fill_geometry(NX,NY+i+1,iz,&geom);		p2u(&get_u(p,0,NX,NY+i+1,iz),&get_u(u,0,NX,NY+i+1,iz),&geom);
	      }

	    PLOOP(iv)
	      set_u(p,iv,NX,NY,iz,.5*(get_u(p,iv,NX-1,NY,iz)+get_u(p,iv,NX,NY-1,iz)));

	    fill_geometry(NX,NY,iz,&geom);	    p2u(&get_u(p,0,NX,NY,iz),&get_u(u,0,NX,NY,iz),&geom);

	    PLOOP(iv)
	      set_u(p,iv,NX+1,NY+1,iz,.5*(get_u(p,iv,NX,NY+1,iz)+get_u(p,iv,NX+1,NY,iz)));

	    fill_geometry(NX+1,NY+1,iz,&geom);	    p2u(&get_u(p,0,NX+1,NY+1,iz),&get_u(u,0,NX+1,NY+1,iz),&geom);
	  }
	}

       //elongated corners along y, filling one cell deep surfaces, and averaging diagonally
      if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
	  for(iy=-NG;iy<NY+NG;iy++) {
	    for(i=0;i<NG-1;i++) {
	      PLOOP(iv) {
		set_u(p,iv,-NG+i,iy,-1,get_u(p,iv,-NG+i,iy,0));
		set_u(p,iv,-1,iy,-NG+i,get_u(p,iv,0,iy,-NG+i)); }
	      fill_geometry(-NG+i,iy,-1,&geom);  p2u(&get_u(p,0,-NG+i,iy,-1),&get_u(u,0,-NG+i,iy,-1),&geom);
	      fill_geometry(-1,iy,-NG+i,&geom);  p2u(&get_u(p,0,-1,iy,-NG+i),&get_u(u,0,-1,iy,-NG+i),&geom);
	    }
      
	    PLOOP(iv)
	      set_u(p,iv,-1,iy,-1,.5*(get_u(p,iv,-1,iy,0)+get_u(p,iv,0,iy,-1)));
	    fill_geometry(-1,iy,-1,&geom);  p2u(&get_u(p,0,-1,iy,-1),&get_u(u,0,-1,iy,-1),&geom);

	    PLOOP(iv)
	      set_u(p,iv,-2,iy,-2,.5*(get_u(p,iv,-2,iy,-1)+get_u(p,iv,-1,iy,-2)));
	    fill_geometry(-2,iy,-2,&geom);   p2u(&get_u(p,0,-2,iy,-2),&get_u(u,0,-2,iy,-2),&geom);
	  }
	}

      if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{
	  for(iy=-NG;iy<NY+NG;iy++) {
	    for(i=0;i<NG-1;i++) {
	      PLOOP(iv) {
		set_u(p,iv,-NG+i,iy,NZ,get_u(p,iv,-NG+i,iy,NZ-1));
		set_u(p,iv,-1,iy,NZ+i+1,get_u(p,iv,0,iy,NZ+i+1)); }

	      fill_geometry(-NG+i,iy,NZ,&geom);  p2u(&get_u(p,0,-NG+i,iy,NZ),&get_u(u,0,-NG+i,iy,NZ),&geom);
	      fill_geometry(-1,iy,NZ+i+1,&geom);  p2u(&get_u(p,0,-1,iy,NZ+i+1),&get_u(u,0,-1,iy,NZ+i+1),&geom);
	    }
      
	    PLOOP(iv)
	      set_u(p,iv,-1,iy,NZ,.5*(get_u(p,iv,-1,iy,NZ-1)+get_u(p,iv,0,iy,NZ)));
	    fill_geometry(-1,iy,NZ,&geom);  p2u(&get_u(p,0,-1,iy,NZ),&get_u(u,0,-1,iy,NZ),&geom);

	    PLOOP(iv)
	      set_u(p,iv,-2,iy,NZ+1,.5*(get_u(p,iv,-2,iy,NZ)+get_u(p,iv,-1,iy,NZ+1)));
	    fill_geometry(-2,iy,NZ+1,&geom);   p2u(&get_u(p,0,-2,iy,NZ+1),&get_u(u,0,-2,iy,NZ+1),&geom);
	  }
	}

      if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
	  for(iy=-NG;iy<NY+NG;iy++) {
	    for(i=0;i<NG-1;i++) {
	      PLOOP(iv) {
		set_u(p,iv,NX,iy,-NG+i,get_u(p,iv,NX-1,iy,-NG+i));
		set_u(p,iv,NX+i+1,iy,-1,get_u(p,iv,NX+i+1,iy,0)); }
	      fill_geometry(NX+i+1,iy,-1,&geom);    p2u(&get_u(p,0,NX+i+1,iy,-1),&get_u(u,0,NX+i+1,iy,-1),&geom);
	      fill_geometry(NX,iy,-NG+i,&geom);    p2u(&get_u(p,0,NX,iy,-NG+i),&get_u(u,0,NX,iy,-NG+i),&geom);
	    }

	    PLOOP(iv)
	      set_u(p,iv,NX,iy,-1,.5*(get_u(p,iv,NX-1,iy,-1)+get_u(p,iv,NX,iy,0)));
	    fill_geometry(NX,iy,-1,&geom); p2u(&get_u(p,0,NX,iy,-1),&get_u(u,0,NX,iy,-1),&geom);

	    PLOOP(iv)
	      set_u(p,iv,NX+1,iy,-2,.5*(get_u(p,iv,NX,iy,-2)+get_u(p,iv,NX+1,iy,-1)));
	    fill_geometry(NX+1,iy,-2,&geom);  p2u(&get_u(p,0,NX+1,iy,-2),&get_u(u,0,NX+1,iy,-2),&geom);

	  }
	}

       if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{
	  for(iy=-NG;iy<NY+NG;iy++) {
	    for(i=0;i<NG-1;i++)
	      {
		PLOOP(iv) {
		  set_u(p,iv,NX+i+1,iy,NZ,get_u(p,iv,NX+i+1,iy,NZ-1));
		  set_u(p,iv,NX,iy,NZ+i+1,get_u(p,iv,NX-1,iy,NZ+i+1));
		}

		fill_geometry(NX+i+1,iy,NZ,&geom);		p2u(&get_u(p,0,NX+i+1,iy,NZ),&get_u(u,0,NX+i+1,iy,NZ),&geom);
		fill_geometry(NX,iy,NZ+i+1,&geom);		p2u(&get_u(p,0,NX,iy,NZ+i+1),&get_u(u,0,NX,iy,NZ+i+1),&geom);
	      }

	    PLOOP(iv)
	      set_u(p,iv,NX,iy,NZ,.5*(get_u(p,iv,NX-1,iy,NZ)+get_u(p,iv,NX,iy,NZ-1)));
	    fill_geometry(NX,iy,NZ,&geom);	    p2u(&get_u(p,0,NX,iy,NZ),&get_u(u,0,NX,iy,NZ),&geom);

	    PLOOP(iv)
	      set_u(p,iv,NX+1,iy,NZ+1,.5*(get_u(p,iv,NX,iy,NZ+1)+get_u(p,iv,NX+1,iy,NZ)));
	    fill_geometry(NX+1,iy,NZ+1,&geom);	    p2u(&get_u(p,0,NX+1,iy,NZ+1),&get_u(u,0,NX+1,iy,NZ+1),&geom);
	  }
	}

      //elongated corners along x, filling one cell deep surfaces, and averaging diagonally
      if(mpi_isitBC_forcorners(YBCLO)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
	  for(ix=-NG;ix<NX+NG;ix++) {
	    for(i=0;i<NG-1;i++) {
	      PLOOP(iv) {
		set_u(p,iv,ix,-NG+i,-1,get_u(p,iv,ix,-NG+i,0));
		set_u(p,iv,ix,-1,-NG+i,get_u(p,iv,ix,0,-NG+i)); }
	      fill_geometry(ix,-NG+i,-1,&geom);  p2u(&get_u(p,0,ix,-NG+i,-1),&get_u(u,0,ix,-NG+i,-1),&geom);
	      fill_geometry(ix,-1,-NG+i,&geom);  p2u(&get_u(p,0,ix,-1,-NG+i),&get_u(u,0,ix,-1,-NG+i),&geom);
	    }
      
	    PLOOP(iv)
	      set_u(p,iv,ix,-1,-1,.5*(get_u(p,iv,ix,-1,0)+get_u(p,iv,ix,0,-1)));
	    fill_geometry(ix,-1,-1,&geom);  p2u(&get_u(p,0,ix,-1,-1),&get_u(u,0,ix,-1,-1),&geom);

	    PLOOP(iv)
	      set_u(p,iv,ix,-2,-2,.5*(get_u(p,iv,ix,-2,-1)+get_u(p,iv,ix,-1,-2)));
	    fill_geometry(ix,-2,-2,&geom);   p2u(&get_u(p,0,ix,-2,-2),&get_u(u,0,ix,-2,-2),&geom);
	  }
	}

      if(mpi_isitBC_forcorners(YBCLO)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{
	  for(ix=-NG;ix<NX+NG;ix++) {
	    for(i=0;i<NG-1;i++) {
	      PLOOP(iv) {
		set_u(p,iv,ix,-NG+i,NZ,get_u(p,iv,ix,-NG+i,NZ-1));
		set_u(p,iv,ix,-1,NZ+i+1,get_u(p,iv,ix,0,NZ+i+1)); }
	      fill_geometry(ix,-NG+i,NZ,&geom);  p2u(&get_u(p,0,ix,-NG+i,NZ),&get_u(u,0,ix,-NG+i,NZ),&geom);
	      fill_geometry(ix,-1,NZ+i+1,&geom);  p2u(&get_u(p,0,ix,-1,NZ+i+1),&get_u(u,0,ix,-1,NZ+i+1),&geom);
	    }
      
	    PLOOP(iv)
	      set_u(p,iv,ix,-1,NZ,.5*(get_u(p,iv,ix,-1,NZ-1)+get_u(p,iv,ix,0,NZ)));
	    fill_geometry(ix,-1,NZ,&geom);  p2u(&get_u(p,0,ix,-1,NZ),&get_u(u,0,ix,-1,NZ),&geom);

	    PLOOP(iv)
	      set_u(p,iv,ix,-2,NZ+1,.5*(get_u(p,iv,ix,-2,NZ)+get_u(p,iv,ix,-1,NZ+1)));
	    fill_geometry(ix,-2,NZ+1,&geom);   p2u(&get_u(p,0,ix,-2,NZ+1),&get_u(u,0,ix,-2,NZ+1),&geom);
	  }
	}

      if(mpi_isitBC_forcorners(YBCHI)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
	  for(ix=-NG;ix<NX+NG;ix++) {
	    for(i=0;i<NG-1;i++) {
	      PLOOP(iv) {
		set_u(p,iv,ix,NY,-NG+i,get_u(p,iv,ix,NY-1,-NG+i));
		set_u(p,iv,ix,NY+i+1,-1,get_u(p,iv,ix,NY+i+1,0)); }

	      fill_geometry(ix,NY+i+1,-1,&geom);    p2u(&get_u(p,0,ix,NY+i+1,-1),&get_u(u,0,ix,NY+i+1,-1),&geom);
	      fill_geometry(ix,NY,-NG+i,&geom);    p2u(&get_u(p,0,ix,NY,-NG+i),&get_u(u,0,ix,NY,-NG+i),&geom);
	    }

	    PLOOP(iv)
	      set_u(p,iv,ix,NY,-1,.5*(get_u(p,iv,ix,NY-1,-1)+get_u(p,iv,ix,NY,0)));
	    fill_geometry(ix,NY,-1,&geom); p2u(&get_u(p,0,ix,NY,-1),&get_u(u,0,ix,NY,-1),&geom);

	    PLOOP(iv)
	      set_u(p,iv,ix,NY+1,-2,.5*(get_u(p,iv,ix,NY,-2)+get_u(p,iv,ix,NY+1,-1)));
	    fill_geometry(ix,NY+1,-2,&geom);  p2u(&get_u(p,0,ix,NY+1,-2),&get_u(u,0,ix,NY+1,-2),&geom);
	  }
	}

       if(mpi_isitBC_forcorners(YBCHI)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{
	  for(ix=-NG;ix<NX+NG;ix++) {
	    for(i=0;i<NG-1;i++)
	      {
		PLOOP(iv) {
		  set_u(p,iv,ix,NY+i+1,NZ,get_u(p,iv,ix,NY+i+1,NZ-1));
		  set_u(p,iv,ix,NY,NZ+i+1,get_u(p,iv,ix,NY-1,NZ+i+1));
		}
		fill_geometry(ix,NY+i+1,NZ,&geom);		p2u(&get_u(p,0,ix,NY+i+1,NZ),&get_u(u,0,ix,NY+i+1,NZ),&geom);
		fill_geometry(ix,NY,NZ+i+1,&geom);		p2u(&get_u(p,0,ix,NY,NZ+i+1),&get_u(u,0,ix,NY,NZ+i+1),&geom);
	      }

	    PLOOP(iv)
	      set_u(p,iv,ix,NY,NZ,.5*(get_u(p,iv,ix,NY-1,NZ)+get_u(p,iv,ix,NY,NZ-1)));
	    fill_geometry(ix,NY,NZ,&geom);	    p2u(&get_u(p,0,ix,NY,NZ),&get_u(u,0,ix,NY,NZ),&geom);

	    PLOOP(iv)
	      set_u(p,iv,ix,NY+1,NZ+1,.5*(get_u(p,iv,ix,NY,NZ+1)+get_u(p,iv,ix,NY+1,NZ)));
	    fill_geometry(ix,NY+1,NZ+1,&geom);	    p2u(&get_u(p,0,ix,NY+1,NZ+1),&get_u(u,0,ix,NY+1,NZ+1),&geom);
	  }
	}

       //total total corners
       //TODO - so far very simplified!!!
       if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCLO)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
	  for(ix=-NG;ix<0;ix++) 
	    for(iy=-NG;iy<0;iy++) 
	      for(iz=-NG;iz<0;iz++) 
		{
		  PLOOP(iv) 
		    set_u(p,iv,ix,iy,iz,get_u(p,iv,0,0,0));
		  fill_geometry(ix,iy,iz,&geom);  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}

       if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCLO)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
	  for(ix=NX;ix<NX+NG;ix++) 
	    for(iy=-NG;iy<0;iy++) 
	      for(iz=-NG;iz<0;iz++) 
		{
		  PLOOP(iv) 
		    set_u(p,iv,ix,iy,iz,get_u(p,iv,NX-1,0,0));
		  fill_geometry(ix,iy,iz,&geom);  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
       
       if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCHI)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
	  for(ix=-NG;ix<0;ix++) 
	    for(iy=NY;iy<NY+NG;iy++) 
	      for(iz=-NG;iz<0;iz++) 
		{
		  PLOOP(iv) 
		    set_u(p,iv,ix,iy,iz,get_u(p,iv,0,NY-1,0));
		  fill_geometry(ix,iy,iz,&geom);  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}

       if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCHI)==1 && mpi_isitBC_forcorners(ZBCLO)==1)
	{
	  for(ix=NX;ix<NX+NG;ix++) 
	    for(iy=NY;iy<NY+NG;iy++) 
	      for(iz=-NG;iz<0;iz++) 
		{
		  PLOOP(iv) 
		    set_u(p,iv,ix,iy,iz,get_u(p,iv,NX-1,NY-1,0));
		  fill_geometry(ix,iy,iz,&geom);  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}

       if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCLO)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{
	  for(ix=-NG;ix<0;ix++) 
	    for(iy=-NG;iy<0;iy++) 
	      for(iz=NZ;iz<NZ+NG;iz++) 
		{
		  PLOOP(iv) 
		    set_u(p,iv,ix,iy,iz,get_u(p,iv,0,0,NZ-1));
		  fill_geometry(ix,iy,iz,&geom);  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}

       if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCLO)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{
	  for(ix=NX;ix<NX+NG;ix++) 
	    for(iy=-NG;iy<0;iy++) 
	      for(iz=NZ;iz<NZ+NG;iz++) 
		{
		  PLOOP(iv) 
		    set_u(p,iv,ix,iy,iz,get_u(p,iv,NX-1,0,NZ-1));
		  fill_geometry(ix,iy,iz,&geom);  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}
       
       if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCHI)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{
	  for(ix=-NG;ix<0;ix++) 
	    for(iy=NY;iy<NY+NG;iy++) 
	      for(iz=NZ;iz<NZ+NG;iz++) 
		{
		  PLOOP(iv) 
		    set_u(p,iv,ix,iy,iz,get_u(p,iv,0,NY-1,NZ-1));
		  fill_geometry(ix,iy,iz,&geom);  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}

       if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCHI)==1 && mpi_isitBC_forcorners(ZBCHI)==1)
	{
	  for(ix=NX;ix<NX+NG;ix++) 
	    for(iy=NY;iy<NY+NG;iy++) 
	      for(iz=NZ;iz<NZ+NG;iz++) 
		{
		  PLOOP(iv) 
		    set_u(p,iv,ix,iy,iz,get_u(p,iv,NX-1,NY-1,NZ-1));
		  fill_geometry(ix,iy,iz,&geom);  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
		}
	}

    }



 ldouble uval[NV],pval[NV];

 /*****************************************************************/
 //corners in the middle - apply boundary condition on what is already in ghost cells
 /*****************************************************************/

 if(TNY>1 && TNZ==1)
   {
     iz=0;
     if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCLO)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=-NG;j<0;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCLO);
	       PLOOP(iv)
	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC_forcorners(XBCLO)==1 && mpi_isitBC_forcorners(YBCHI)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCLO);
	       PLOOP(iv)
	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCLO)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCHI);
	       PLOOP(iv)
	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC_forcorners(XBCHI)==1 && mpi_isitBC_forcorners(YBCHI)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCHI);
	       PLOOP(iv)
	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC_forcorners(YBCLO)==1 && mpi_isitBC_forcorners(XBCLO)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=-NG;j<0;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCLO);
	       PLOOP(iv)
	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC_forcorners(YBCLO)==1 && mpi_isitBC_forcorners(XBCHI)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCLO);
	       PLOOP(iv)
	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC_forcorners(YBCHI)==1 && mpi_isitBC_forcorners(XBCLO)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCHI);
	       PLOOP(iv)
	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC_forcorners(YBCHI)==1 && mpi_isitBC_forcorners(XBCHI)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCHI);
	       PLOOP(iv)
	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
   }

 /****************************/
 // Brandon - here we need to call mpi_isitBC since this is asking the question, 
 // is this border exchanged through MPI? 
 // If the answer is no, it will use regular bc, which we don't want for transmitting boundary
 /****************************/
 if(TNY>1 && TNZ>1) //full 3d
   {
     //elongated along z
     if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCLO)==0)
       {
	 for(iz=0;iz<NZ;iz++)
	   for(i=-NG;i<0;i++)
	     for(j=-NG;j<0;j++)
	       {
		 set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCLO);
		 PLOOP(iv)		 {
		   set_u(u,iv,i,j,iz,uval[iv]);
		   set_u(p,iv,i,j,iz,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCHI)==0)
       {
	 for(iz=0;iz<NZ;iz++)
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCLO);
	       PLOOP(iv)	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCLO)==0)
       {
	 for(iz=0;iz<NZ;iz++)
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCHI);
	       PLOOP(iv)	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCHI)==0)
       {
	 for(iz=0;iz<NZ;iz++)
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCHI);
	       PLOOP(iv)	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC(YBCLO)==1 && mpi_isitBC(XBCLO)==0)
       {
	 for(iz=0;iz<NZ;iz++)
	 for(i=-NG;i<0;i++)
	   for(j=-NG;j<0;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCLO);
	       PLOOP(iv)	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC(YBCLO)==1 && mpi_isitBC(XBCHI)==0)
       {
	 for(iz=0;iz<NZ;iz++)
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCLO);
	       PLOOP(iv)	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC(YBCHI)==1 && mpi_isitBC(XBCLO)==0)
       {
	 for(iz=0;iz<NZ;iz++)
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCHI);
	       PLOOP(iv)	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }
     if(mpi_isitBC(YBCHI)==1 && mpi_isitBC(XBCHI)==0)
       {
	 for(iz=0;iz<NZ;iz++)
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     {
	       set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCHI);
	       PLOOP(iv)	       {
		 set_u(u,iv,i,j,iz,uval[iv]);
		 set_u(p,iv,i,j,iz,pval[iv]);	      
	       }
	     }
       }

     //elongated along y
     if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(ZBCLO)==0)
       {
	 for(iy=0;iy<NY;iy++)
	   for(i=-NG;i<0;i++)
	     for(j=-NG;j<0;j++)
	       {
		 set_bc_core(i,iy,j,t,uval,pval,ifinit,XBCLO);
		 PLOOP(iv)		 {
		   set_u(u,iv,i,iy,j,uval[iv]);
		   set_u(p,iv,i,iy,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(ZBCHI)==0)
       {
	 for(iy=0;iy<NY;iy++)
	   for(i=-NG;i<0;i++)
	     for(j=NZ;j<NZ+NG;j++)
	       {
		 set_bc_core(i,iy,j,t,uval,pval,ifinit,XBCLO);
		 PLOOP(iv)	       {
		   set_u(u,iv,i,iy,j,uval[iv]);
		   set_u(p,iv,i,iy,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(ZBCLO)==0)
       {
	 for(iy=0;iy<NY;iy++)
	   for(i=NX;i<NX+NG;i++)
	     for(j=-NG;j<0;j++)
	       {
		 set_bc_core(i,iy,j,t,uval,pval,ifinit,XBCHI);
		 PLOOP(iv)	       {
		   set_u(u,iv,i,iy,j,uval[iv]);
		   set_u(p,iv,i,iy,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(ZBCHI)==0)
       {
	 for(iy=0;iy<NY;iy++)
	   for(i=NX;i<NX+NG;i++)
	     for(j=NZ;j<NZ+NG;j++)
	       {
		 set_bc_core(i,iy,j,t,uval,pval,ifinit,XBCHI);
		 PLOOP(iv)	       {
		   set_u(u,iv,i,iy,j,uval[iv]);
		   set_u(p,iv,i,iy,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(ZBCLO)==1 && mpi_isitBC(XBCLO)==0)
       {
	 for(iy=0;iy<NY;iy++)
	   for(i=-NG;i<0;i++)
	     for(j=-NG;j<0;j++)
	       {
		 set_bc_core(i,iy,j,t,uval,pval,ifinit,ZBCLO);
		 PLOOP(iv)	       {
		   set_u(u,iv,i,iy,j,uval[iv]);
		   set_u(p,iv,i,iy,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(ZBCLO)==1 && mpi_isitBC(XBCHI)==0)
       {
	 for(iy=0;iy<NY;iy++)
	   for(i=NX;i<NX+NG;i++)
	     for(j=-NG;j<0;j++)
	       {
		 set_bc_core(i,iy,j,t,uval,pval,ifinit,ZBCLO);
		 PLOOP(iv)	       {
		   set_u(u,iv,i,iy,j,uval[iv]);
		   set_u(p,iv,i,iy,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(ZBCHI)==1 && mpi_isitBC(XBCLO)==0)
       {
	 for(iy=0;iy<NY;iy++)
	   for(i=-NG;i<0;i++)
	     for(j=NZ;j<NZ+NG;j++)
	       {
		 set_bc_core(i,iy,j,t,uval,pval,ifinit,ZBCHI);
		 PLOOP(iv)	       {
		   set_u(u,iv,i,iy,j,uval[iv]);
		   set_u(p,iv,i,iy,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(ZBCHI)==1 && mpi_isitBC(XBCHI)==0)
       {
	 for(iy=0;iy<NY;iy++)
	   for(i=NX;i<NX+NG;i++)
	     for(j=NZ;j<NZ+NG;j++)
	       {
		 set_bc_core(i,iy,j,t,uval,pval,ifinit,ZBCHI);
		 PLOOP(iv)	       {
		   set_u(u,iv,i,iy,j,uval[iv]);
		   set_u(p,iv,i,iy,j,pval[iv]);	      
		 }
	       }
       }

     //elongated along x
     if(mpi_isitBC(YBCLO)==1 && mpi_isitBC(ZBCLO)==0)
       {
	 for(ix=0;ix<NX;ix++)
	   for(i=-NG;i<0;i++)
	     for(j=-NG;j<0;j++)
	       {
		 set_bc_core(ix,i,j,t,uval,pval,ifinit,YBCLO);
		 PLOOP(iv)		 {
		   set_u(u,iv,ix,i,j,uval[iv]);
		   set_u(p,iv,ix,i,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(YBCLO)==1 && mpi_isitBC(ZBCHI)==0)
       {
	 for(ix=0;ix<NX;ix++)
	   for(i=-NG;i<0;i++)
	     for(j=NZ;j<NZ+NG;j++)
	       {
		 set_bc_core(ix,i,j,t,uval,pval,ifinit,YBCLO);
		 PLOOP(iv)	       {
		   set_u(u,iv,ix,i,j,uval[iv]);
		   set_u(p,iv,ix,i,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(YBCHI)==1 && mpi_isitBC(ZBCLO)==0)
       {
	 for(ix=0;ix<NX;ix++)
	   for(i=NY;i<NY+NG;i++)
	     for(j=-NG;j<0;j++)
	       {
		 set_bc_core(ix,i,j,t,uval,pval,ifinit,YBCHI);
		 PLOOP(iv)	       {
		   set_u(u,iv,ix,i,j,uval[iv]);
		   set_u(p,iv,ix,i,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(YBCHI)==1 && mpi_isitBC(ZBCHI)==0)
       {
	 for(ix=0;ix<NX;ix++)
	   for(i=NY;i<NY+NG;i++)
	     for(j=NZ;j<NZ+NG;j++)
	       {
		 set_bc_core(ix,i,j,t,uval,pval,ifinit,YBCHI);
		 PLOOP(iv)	       {
		   set_u(u,iv,ix,i,j,uval[iv]);
		   set_u(p,iv,ix,i,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(ZBCLO)==1 && mpi_isitBC(YBCLO)==0)
       {
	 for(ix=0;ix<NX;ix++)
	   for(i=-NG;i<0;i++)
	     for(j=-NG;j<0;j++)
	       {
		 set_bc_core(ix,i,j,t,uval,pval,ifinit,ZBCLO);
		 PLOOP(iv)	       {
		   set_u(u,iv,ix,i,j,uval[iv]);
		   set_u(p,iv,ix,i,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(ZBCLO)==1 && mpi_isitBC(YBCHI)==0)
       {
	 for(ix=0;ix<NX;ix++)
	   for(i=NY;i<NY+NG;i++)
	     for(j=-NG;j<0;j++)
	       {
		 set_bc_core(ix,i,j,t,uval,pval,ifinit,ZBCLO);
		 PLOOP(iv)	       {
		   set_u(u,iv,ix,i,j,uval[iv]);
		   set_u(p,iv,ix,i,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(ZBCHI)==1 && mpi_isitBC(YBCLO)==0)
       {
	 for(ix=0;ix<NX;ix++)
	   for(i=-NG;i<0;i++)
	     for(j=NZ;j<NZ+NG;j++)
	       {
		 set_bc_core(ix,i,j,t,uval,pval,ifinit,ZBCHI);
		 PLOOP(iv)	       {
		   set_u(u,iv,ix,i,j,uval[iv]);
		   set_u(p,iv,ix,i,j,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(ZBCHI)==1 && mpi_isitBC(YBCHI)==0)
       {
	 for(ix=0;ix<NX;ix++)
	   for(i=NY;i<NY+NG;i++)
	     for(j=NZ;j<NZ+NG;j++)
	       {
		 set_bc_core(ix,i,j,t,uval,pval,ifinit,ZBCHI);
		 PLOOP(iv)	       {
		   set_u(u,iv,ix,i,j,uval[iv]);
		   set_u(p,iv,ix,i,j,pval[iv]);	      
		 }
	       }
       }

     //corners corners but withing the domain
     //protruding only in x
     if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=-NG;j<0;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,XBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=-NG;j<0;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,XBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,XBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,XBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
      if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,XBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,XBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,XBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,XBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     //protruding only in y
     if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==1 && mpi_isitBC(ZBCLO)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=-NG;j<0;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,YBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==1 && mpi_isitBC(ZBCHI)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=-NG;j<0;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,YBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==1 && mpi_isitBC(ZBCLO)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,YBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==1 && mpi_isitBC(ZBCHI)==0)
       {
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,YBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==1 && mpi_isitBC(ZBCLO)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,YBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==1 && mpi_isitBC(ZBCHI)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,YBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==1 && mpi_isitBC(ZBCLO)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,YBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==1 && mpi_isitBC(ZBCHI)==0)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,YBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }

     //protruding only in z
     if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==1)
       {
	 for(i=-NG;i<0;i++)
	   for(j=-NG;j<0;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,ZBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==1)
       {
	 for(i=-NG;i<0;i++)
	   for(j=-NG;j<0;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,ZBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==1)
       {
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,ZBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==1)
       {
	 for(i=-NG;i<0;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,ZBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
      if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==1)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,YBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==1)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=-NG;j<0;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,ZBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==1)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=-NG;k<0;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,ZBCLO);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }
     if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==1)
       {
	 for(i=NX;i<NX+NG;i++)
	   for(j=NY;j<NY+NG;j++)
	     for(k=NZ;k<NZ+NG;k++)
	       {
		 set_bc_core(i,j,k,t,uval,pval,ifinit,ZBCHI);
		 PLOOP(iv)     	 {
		   set_u(u,iv,i,j,k,uval[iv]);
		   set_u(p,iv,i,j,k,pval[iv]);	      
		 }
	       }
       }

#endif //SHEARINGBOX

   }
  
#endif

 #ifdef CONSISTENTGAMMA
 //update adiabatic indices
 set_gammagas(0);
 #endif

  return 0;
} //set_bc()


//*********************************************************
//fixing after failed MHD main inversion
//*********************************************************

int
cell_fixup(int type)
{
  if(DOFIXUPS==0)
    return 0;

  if(DOU2PMHDFIXUPS==0 && type==FIXUP_U2PMHD)
    return 0;

  if(DOU2PRADFIXUPS==0 && type==FIXUP_U2PRAD)
    return 0;

  if(DORADIMPFIXUPS==0 && type==FIXUP_RADIMP)
    return 0;

  int ix,iy,iz,iv;
  int in,ii,iii;
  int verbose=1;

  copyi_u(1.,u,u_bak_fixup);
  copyi_u(1.,p,p_bak_fixup);

  //gets the neighboring the primitives
#pragma omp parallel for private(ix,iy,iz,iv,iii,in) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain only
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      //do not correct if overwritten later on
      if(is_cell_corrected_polaraxis(ix,iy,iz)) continue;


      if(((get_cflag(HDFIXUPFLAG,ix,iy,iz)!=0 && type==FIXUP_U2PMHD) ||
	  (get_cflag(RADFIXUPFLAG,ix,iy,iz)!=0 && type==FIXUP_U2PRAD) ||
	  (get_cflag(RADIMPFIXUPFLAG,ix,iy,iz)!=0 && type==FIXUP_RADIMP)) && is_cell_active(ix,iy,iz)
	)
	{
	  //ANDREW - I think maybe this should be here, to set the flag back to its default? 
	  //this should not be here, should it?
	  /*
	  if(type==FIXUP_U2PMHD) set_cflag(HDFIXUPFLAG,ix,iy,iz,0); //try only once
	  if(type==FIXUP_U2PRAD) set_cflag(RADFIXUPFLAG,ix,iy,iz,0); //try only once
	  if(type==FIXUP_RADIMP) set_cflag(RADIMPFIXUPFLAG,ix,iy,iz,0); //try only once
	  */

	  //looking around for good neighbors
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);

	  ldouble ppn[6][NV],pp[NV],uu[NV];

	  //should care about global but at the stage where it is called knowns not about the boundaries
	  //so fixups idividually in tiles but later exchanged 

	  in=0; //number of successful neighbors
		  
	  if(ix-1>=0 &&  
	     ((get_cflag(HDFIXUPFLAG,ix-1,iy,iz)==0 && type==FIXUP_U2PMHD) ||
	      (get_cflag(RADFIXUPFLAG,ix-1,iy,iz)==0 && type==FIXUP_U2PRAD) ||
	      (get_cflag(RADIMPFIXUPFLAG,ix-1,iy,iz)==0 && type==FIXUP_RADIMP)))	      
	    {
              #ifdef SPECIAL_BC_CHECK //make sure that ix-1 is not a stream cell
              if(TNY>1 && TNZ==1)
              {
                if((iy+TOJ) >= STREAM_IYT && (iy+TOJ) <= STREAM_IYB && (ix+TOI) > (STREAM_IX+1))
                {
                  if((ix-1+TOI) > (STREAM_IX+3)) //not a stream cell in ix-1
                  {
                    in++;
                    for(iv=0;iv<NV;iv++)
                      ppn[in-1][iv]=get_u(p,iv,ix-1,iy,iz);
                  }
                }
                else
                {
                  in++;
                  for(iv=0;iv<NV;iv++)
                    ppn[in-1][iv]=get_u(p,iv,ix-1,iy,iz);
                }
              }
              #else
	      in++;
	      for(iv=0;iv<NV;iv++)
    		ppn[in-1][iv]=get_u(p,iv,ix-1,iy,iz);
              #endif
	    }

	  if(ix+1<NX &&
	     ((get_cflag(HDFIXUPFLAG,ix+1,iy,iz)==0 && type==FIXUP_U2PMHD) ||
	      (get_cflag(RADFIXUPFLAG,ix+1,iy,iz)==0 && type==FIXUP_U2PRAD) ||
	      (get_cflag(RADIMPFIXUPFLAG,ix+1,iy,iz)==0 && type==FIXUP_RADIMP)))
	    {
              #ifdef SPECIAL_BC_CHECK //make sure that ix-1 is not a stream ghost cell
              if(TNY>1 && TNZ==1)
              {
                if((iy+TOJ) >= STREAM_IYT && (iy+TOJ) <= STREAM_IYB && (ix+TOI) < STREAM_IX)
                {
                  if((ix+1+TOI) < STREAM_IX) //not a stream cell in ix+1
                  {
                    in++;
                    for(iv=0;iv<NV;iv++)
                      ppn[in-1][iv]=get_u(p,iv,ix+1,iy,iz);
                  }
                }
                else //not in danger of using a stream ghost to do a fixup
                {
                  in++;
                  for(iv=0;iv<NV;iv++)
                    ppn[in-1][iv]=get_u(p,iv,ix+1,iy,iz);
                }
              }
              #else
	      in++;
	      for(iv=0;iv<NV;iv++)
		ppn[in-1][iv]=get_u(p,iv,ix+1,iy,iz);
              #endif
	    }

	  if(iy-1>=0 &&
	     ((get_cflag(HDFIXUPFLAG,ix,iy-1,iz)==0 && type==FIXUP_U2PMHD) ||
	      (get_cflag(RADFIXUPFLAG,ix,iy-1,iz)==0 && type==FIXUP_U2PRAD) ||
	      (get_cflag(RADIMPFIXUPFLAG,ix,iy-1,iz)==0 && type==FIXUP_RADIMP)))
	    {
              #ifdef SPECIAL_BC_CHECK //make sure that ix-1 is not a stream ghost cell
              if(TNY>1 && TNZ==1)
              {
                if((ix+TOI) >= STREAM_IX && (ix+TOI) <= (STREAM_IX+3) && (iy+TOJ) > STREAM_IYB)
                {
                  if((iy-1+TOJ) > STREAM_IYB) //not a stream cell in iy-1
                  {
                    in++;
                    for(iv=0;iv<NV;iv++)
                      ppn[in-1][iv]=get_u(p,iv,ix,iy-1,iz);
                  }
                }
                else //not in danger of using a stream ghost to do a fixup
                {
                  in++;
                  for(iv=0;iv<NV;iv++)
                    ppn[in-1][iv]=get_u(p,iv,ix,iy-1,iz);
                }
              }
              #else
	      in++;
	      for(iv=0;iv<NV;iv++)
		ppn[in-1][iv]=get_u(p,iv,ix,iy-1,iz);
              #endif
	    }

	  if(iy+1<NY &&
	     ((get_cflag(HDFIXUPFLAG,ix,iy+1,iz)==0 && type==FIXUP_U2PMHD) ||
	      (get_cflag(RADFIXUPFLAG,ix,iy+1,iz)==0 && type==FIXUP_U2PRAD) ||
	      (get_cflag(RADIMPFIXUPFLAG,ix,iy+1,iz)==0 && type==FIXUP_RADIMP)))
	    {
              #ifdef SPECIAL_BC_CHECK //make sure that iy+1 is not a stream ghost cell
              if(TNY>1 && TNZ==1)
              {
                if((ix+TOI) >= STREAM_IX && (ix+TOI) <= (STREAM_IX+3) && (iy+TOJ) < STREAM_IYT)
                {
                  if((iy+1+TOJ) < STREAM_IYT) //not a stream cell in iy+1
                  {
                    in++;
                    for(iv=0;iv<NV;iv++)
                      ppn[in-1][iv]=get_u(p,iv,ix,iy+1,iz);
                  }
                }
                else //not in danger of using a stream ghost to do a fixup
                {
                  in++;
                  for(iv=0;iv<NV;iv++)
                    ppn[in-1][iv]=get_u(p,iv,ix,iy+1,iz);
                }
              }
              #else
	      in++;
	      for(iv=0;iv<NV;iv++)
		ppn[in-1][iv]=get_u(p,iv,ix,iy+1,iz);
              #endif
	    }

	  if(iz-1>=0 &&
	     ((get_cflag(HDFIXUPFLAG,ix,iy,iz-1)==0 && type==FIXUP_U2PMHD) ||
	      (get_cflag(RADFIXUPFLAG,ix,iy,iz-1)==0 && type==FIXUP_U2PRAD) ||
	      (get_cflag(RADIMPFIXUPFLAG,ix,iy,iz-1)==0 && type==FIXUP_RADIMP)))
	    {
	      in++;
	      for(iv=0;iv<NV;iv++)
		ppn[in-1][iv]=get_u(p,iv,ix,iy,iz-1);
	    }

	  if(iz+1<NZ  &&
	     ((get_cflag(HDFIXUPFLAG,ix,iy,iz+1)==0 && type==FIXUP_U2PMHD) ||
	      (get_cflag(RADFIXUPFLAG,ix,iy,iz+1)==0 && type==FIXUP_U2PRAD) ||
	      (get_cflag(RADIMPFIXUPFLAG,ix,iy,iz+1)==0 && type==FIXUP_RADIMP)))
	    {
	      in++;
	      for(iv=0;iv<NV;iv++)
		ppn[in-1][iv]=get_u(p,iv,ix,iy,iz+1);
	    }

	  if((NZ==1 && NY==1 && in>=1) ||
	     (NZ==1 && in>=2) ||
	     (NY==1 && in>=2) ||
	     in>=3) //are there sufficient number of neighbors to do fixup
	    {
	      for(iv=0;iv<NV;iv++)
		{
		  int fixthis=0;

		  if(type==FIXUP_U2PMHD &&  iv!=RHO && iv<B1) //skip correcting magnetic field not to disrupt div B
		    fixthis=1;

		  if(type==FIXUP_U2PRAD &&  iv!=RHO && iv>=EE && iv<=FZ) //fix only radiative quantites
		    fixthis=1;

		  if(type==FIXUP_RADIMP &&  iv!=RHO && (iv<B1 || (iv>=EE && iv<=FZ))) //fix both mhd and rad but not magn field
		    fixthis=1;

                  #ifdef EVOLVEPHOTONNUMBER
		  if(type==FIXUP_RADIMP &&  iv==NF) //fix number of photons
		    fixthis=1;
                  #endif

                  #ifdef EVOLVEELECTRONS
		  if(type==FIXUP_RADIMP &&  (iv==ENTRE || iv==ENTRI)) //fix electron/ion entropy
		    fixthis=1;
                  #endif

		  if(fixthis==1)
		    {
		      pp[iv]=0.;
		      for(iii=0;iii<in;iii++)
			pp[iv]+=ppn[iii][iv];
		      pp[iv]/=(ldouble)in;  
		    }
		  else //leave as was
		    {
		      pp[iv]=get_u(p,iv,ix,iy,iz); 
		    }
		}
	      
	      p2u(pp,uu,&geom);

	      if(verbose>1) 
		{
		  if(type==FIXUP_U2PMHD) printf("%4d > %4d %4d %4d > U2PMHD > fixing up with %d neighbors\n",PROCID,ix+TOI,iy+TOJ,iz+TOK,in);
		  if(type==FIXUP_U2PRAD) printf("%4d > %4d %4d %4d > U2PRAD > fixing up with %d neighbors\n",PROCID,ix+TOI,iy+TOJ,iz+TOK,in);
		  if(type==FIXUP_RADIMP) printf("%4d > %4d %4d %4d > RADIMP > fixing up with %d neighbors\n",PROCID,ix+TOI,iy+TOJ,iz+TOK,in);
		}
	      
	      //save to updated arrays memory
	      for(iv=0;iv<NV;iv++)
		{
		  set_u(u_bak_fixup,iv,ix,iy,iz,uu[iv]);
		  set_u(p_bak_fixup,iv,ix,iy,iz,pp[iv]);
		}
	     }
	     else
	     {
	       if(type==FIXUP_U2PMHD) printf("%4d > %4d %4d %4d > U2PMHD > didn't manage to hd fixup\n",PROCID,ix+TOI,iy+TOJ,iz+TOK);
	       if(type==FIXUP_U2PRAD) printf("%4d > %4d %4d %4d > U2PRAD > didn't manage to hd fixup\n",PROCID,ix+TOI,iy+TOJ,iz+TOK);
	       if(type==FIXUP_RADIMP) printf("%4d > %4d %4d %4d > RADIMP > didn't manage to hd fixup\n",PROCID,ix+TOI,iy+TOJ,iz+TOK);
	     }
	}
    }

  //restoring to memory
  copyi_u(1.,u_bak_fixup,u);
  copyi_u(1.,p_bak_fixup,p);

  return 0;
}


//**********************************************************************
//averages the corrected cells over azimuth not touching magn. field
//**********************************************************************

int
smooth_polaraxis()
{
  int nc=NCCORRECTPOLAR; //correct velocity in nc most polar cells;

  //spherical like coords
  if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS ||
      MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS || MYCOORDS==MKS3COORDS || MYCOORDS==TKS3COORDS ||
      MYCOORDS==MSPH1COORDS || MYCOORDS==JETCOORDS)
    {
      int ix;
//#pragma omp parallel for private(ix)
      for(ix=0;ix<NX;ix++)
	{
      int iy,iz,iv,ic,iysrc,ixsrc;
      ldouble pp[NV],uu[NV];
      struct geometry geom;
      
	  //upper axis
#ifdef MPI
	  if(TJ==0)
#endif
	    {
	      for(ic=0;ic<nc;ic++)
		{
		  iy=ic;
		  PLOOP(iv)
		    pp[iv]=0.;
		  for(iz=0;iz<NZ;iz++)
		    {
		      PLOOP(iv)
			pp[iv]+=get_u(p,iv,ix,iy,iz);
		    }
		  PLOOP(iv)
		    pp[iv]/=NZ;
		  for(iz=0;iz<NZ;iz++)
		    {
		      fill_geometry(ix,iy,iz,&geom);
		      //recover magnetic field from this cell
#ifdef MAGNFIELD
		      pp[B1]=get_u(p,B1,ix,iy,iz);
		      pp[B2]=get_u(p,B2,ix,iy,iz);
		      pp[B3]=get_u(p,B3,ix,iy,iz);
#endif		      
		      p2u(pp,uu,&geom);
		      PLOOP(iv)
		      {
			if(iv<B1 || iv>B3) { //skip magnetic field
			  set_u(p,iv,ix,iy,iz,pp[iv]);  
			  set_u(u,iv,ix,iy,iz,uu[iv]); }
		      }
		    }
		}
	    }
	  //lower axis
#ifndef HALFTHETA
#ifdef MPI
	  if(TJ==NTY-1)
#endif
	    {
	      for(ic=0;ic<nc;ic++)
		{
		  iy=NY-1-ic;
		  PLOOP(iv)
		    pp[iv]=0.;
		  for(iz=0;iz<NZ;iz++)
		    {
		      PLOOP(iv)
			pp[iv]+=get_u(p,iv,ix,iy,iz);
		    }
		  PLOOP(iv)
		    pp[iv]/=NZ;
		  for(iz=0;iz<NZ;iz++)
		    {
		      fill_geometry(ix,iy,iz,&geom);
		      //recover magnetic field from this cell
#ifdef MAGNFIELD
		      pp[B1]=get_u(p,B1,ix,iy,iz);
		      pp[B2]=get_u(p,B2,ix,iy,iz);
		      pp[B3]=get_u(p,B3,ix,iy,iz);
#endif
		      p2u(pp,uu,&geom);
		      PLOOP(iv)
		      {
			if(iv<B1 || iv>B3) { //skip magnetic field
			  set_u(p,iv,ix,iy,iz,pp[iv]);  
			  set_u(u,iv,ix,iy,iz,uu[iv]); }
		      }
		    }
		}
	    }
#endif
	}
    }

return 0;
}


//**********************************************************************
// correct cells near the surface of the neutron star
// by interpolating the velocities to v_in
//**********************************************************************

#ifdef CORRECT_NSSURFACE
int
correct_nssurface()
{   
    #ifdef MPI
    if(TI != 0) return 0; // only first X cells in the first radial tile do fixup
    #endif

    // number of cells to correct
    int nc = NCCORRECT_NSSURFACE;
    
    // variables to hold cell index
    int ix, iy, iz, iv, ic;
    
    // iterate over all y and z and ix up to ic
    #pragma omp parallel for private(ic,ix,iy,iz,iv) schedule (static)
    for(ix=0; ix<nc; ix++) // interating up to cell ix=nc-1
    {
        for(iy=0; iy<NY; iy++)
        {
            for(iz=0; iz<NZ; iz++)
            {   
                // arrays of primatives and conserved quantities for this cell
                ldouble pp[NV], uu[NV];

                // get the geometries at cell nc to extrapolate from(ix=nc)
                struct geometry geomNC;
                fill_geometry(nc, iy, iz, &geomNC);
                struct geometry geomBLNC;
                fill_geometry_arb(nc, iy, iz, &geomBLNC,BLCOORDS);
                
                // fill the primatives at the cell to extrapolate from
                ldouble ppnc[NV];
                for(iv=0; iv<NV; iv++) ppnc[iv]= get_u(p, iv, nc, iy, iz);
                
                //transform to BL coordinates at the cell to extrapolate from
#ifdef PRECOMPUTE_MY2OUT
		trans_pall_coco_my2out(ppnc, ppnc, &geomNC, &geomBLNC);
#else
                trans_pall_coco(ppnc, ppnc, MYCOORDS, BLCOORDS, geomNC.xxvec, &geomNC, &geomBLNC);
#endif

	    
                // get the value for velocity at cell nc to extrapolate
                ldouble vnc;
                vnc = ppnc[VX];
                // get the value of r at cel nc
                ldouble rnc;
                rnc = geomBLNC.xx; // 0 should correspond to x
                
                // get mdot from r, v, rho at cel lnc
                ldouble rhonc, mdotnc;
                rhonc = ppnc[RHO];
                mdotnc = rhonc*vnc*rnc*rnc;
                         
                // boundary values from define.h
                ldouble vb = -V_IN;
                ldouble rb = RNS;

                // calculate slope for v
                ldouble m = (vb - vnc)/(rb-rnc);
                
                // get the geometries at the current cell to correct
                struct geometry geom;
                fill_geometry(ix, iy, iz, &geom);

                struct geometry geomBL;
                fill_geometry_arb(ix, iy, iz, &geomBL, BLCOORDS);

                // get r of the current cell to correct
                ldouble r;
                r = geomBL.xx; // 0 should corresopnd to x

                // linear extrapolation to the current cell
                ldouble v;
                v = m*(r-rb) + vb;
                
                // set rho to give consitant mdot
                ldouble rho;
                rho = mdotnc/v/r/r; 

                // get the primatives at the cell to correct
                for (iv=0; iv<NV; iv++) pp[iv]=get_u(p, iv, ix, iy, iz);

                //transform to BL coordinates at cell to correct
#ifdef PRECOMPUTE_MY2OUT
                trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
#else
                trans_pall_coco(pp,pp,MYCOORDS, BLCOORDS, geom.xxvec, &geom, &geomBL);
#endif
                //write the new v to pp
                pp[VX] = v;
                
                // write the new rho
                pp[RHO] = rho;

                // transfrom back to code coordinates
#ifdef PRECOMPUTE_MY2OUT
                trans_pall_coco_out2my(pp,pp,&geomBL,&geom,);
#else		
                trans_pall_coco(pp,pp,BLCOORDS, MYCOORDS,geomBL.xxvec, &geomBL, &geom);
#endif                
                // compute conserved quantities
                p2u(pp,uu,&geom);
                //save to memory
                for(iv=0;iv<NV; iv++)
                {
                    set_u(p,iv,ix,iy,iz,pp[iv]);
                    set_u(u,iv,ix,iy,iz,uu[iv]);
                }


            }   
        }
    }



    return 0;
}
#endif


//**********************************************************************
//treats the most polar cells is a special way, correcting them, not evolving them
//**********************************************************************

int
correct_polaraxis()
{
      int nc=NCCORRECTPOLAR; //correct velocity in nc most polar cells;

      int ix,iy,iz,iv,ic,iysrc,ixsrc;
      
      //spherical like coords
      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS ||
	  MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS || MYCOORDS==MKS3COORDS || MYCOORDS==TKS3COORDS ||
	  MYCOORDS==MSPH1COORDS || MYCOORDS==JETCOORDS)
	{
	  #pragma omp parallel for private(ic,ix,iy,iz,iv,iysrc) schedule (static)
	  for(ix=0;ix<NX;ix++)
	    {
	      for(iz=0;iz<NZ;iz++)
		{
		  ldouble th,thsrc,thaxis;
		  ldouble pp[NV],uu[NV];
		  struct geometry geom;

		  //upper
#ifdef MPI
		  if(TJ==0) //tile number
#endif
		    {
		      thaxis=get_xb(0,1);
		      for(ic=0;ic<nc;ic++)
			{
			  iy=ic;iysrc=nc;
			  th=get_x(iy,1);
			  thsrc=get_x(iysrc,1);	      
	      	  
			  fill_geometry(ix,iy,iz,&geom);
	  
			  PLOOP(iv)
			    pp[iv]=get_u(p,iv,ix,iy,iz);
			  
			  //gas densities
			  pp[RHO]=get_u(p,RHO,ix,iysrc,iz);
			  pp[UU]=get_u(p,UU,ix,iysrc,iz);
			  pp[ENTR]=get_u(p,ENTR,ix,iysrc,iz);		  
#ifdef EVOLVEELECTRONS
			  pp[ENTRE]=get_u(p,ENTRE,ix,iysrc,iz);		  
			  pp[ENTRI]=get_u(p,ENTRI,ix,iysrc,iz);		  
#ifdef RELELECTRONS
                          int ie;
                          for(ie=0;ie<NRELBIN;ie++) pp[NEREL(ie)]=get_u(p,NEREL(ie),ix,iysrc,iz);
#endif
#endif
		  		  		  
			  //gas velocities
			  pp[VX]=get_u(p,VX,ix,iysrc,iz);
			  pp[VZ]=get_u(p,VZ,ix,iysrc,iz);
			  pp[VY]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,VY,ix,iysrc,iz);
			  
#ifdef MAGNFIELD
#ifdef CORRECTMAGNFIELD
			  //do overwrite magnetic field, div B!=0 does not propagate in
			  pp[B1]=get_u(p,B1,ix,iysrc,iz);
			  pp[B3]=get_u(p,B3,ix,iysrc,iz);
			  pp[B2]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,B2,ix,iysrc,iz);
#endif
#endif

#ifdef RADIATION
			  //rad density
			  pp[EE0]=get_u(p,EE0,ix,iysrc,iz);

#ifdef EVOLVEPHOTONNUMBER
			  //no. of photons
			  pp[NF0]=get_u(p,NF0,ix,iysrc,iz);
#endif

			  //rad velocities
			  pp[FX0]=get_u(p,FX0,ix,iysrc,iz);
			  pp[FZ0]=get_u(p,FZ0,ix,iysrc,iz);
			  pp[FY0]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,FY0,ix,iysrc,iz);
       

#endif 
			  
			  //update gamma in the corrected cell
#ifdef CONSISTENTGAMMA
			  ldouble newgamma=calc_gammagas(pp,ix,iy,iz);
			  set_u_scalar(gammagas,ix,iy,iz,newgamma);
#endif

			  p2u(pp,uu,&geom);

			  PLOOP(iv)
			  {
			    set_u(p,iv,ix,iy,iz,pp[iv]);  
			    set_u(u,iv,ix,iy,iz,uu[iv]); 
			  }
			}
		    }

		  //lower
#ifndef HALFTHETA
#ifdef MPI
		  if(TJ==NTY-1)
#endif
		    {
		      thaxis=get_xb(NY,1);
		      for(ic=0;ic<nc;ic++)
			{
			  iy=NY-1-ic;iysrc=NY-1-nc;
			  th=get_x(iy,1);
			  thsrc=get_x(iysrc,1);	      
	      	  
			  fill_geometry(ix,iy,iz,&geom);
	  
			  PLOOP(iv)
			    pp[iv]=get_u(p,iv,ix,iy,iz);
  
			  //gas densities
			  pp[RHO]=get_u(p,RHO,ix,iysrc,iz);
			  pp[UU]=get_u(p,UU,ix,iysrc,iz);
			  pp[ENTR]=get_u(p,ENTR,ix,iysrc,iz);	

#ifdef EVOLVEELECTRONS
			  pp[ENTRE]=get_u(p,ENTRE,ix,iysrc,iz);		  
			  pp[ENTRI]=get_u(p,ENTRI,ix,iysrc,iz);
#ifdef RELELECTRONS

                          int ie;
                          for(ie=0;ie<NRELBIN;ie++) pp[NEREL(ie)]=get_u(p,NEREL(ie),ix,iysrc,iz);
#endif		  
#endif	  

			  //gas velocities
			  pp[VX]=get_u(p,VX,ix,iysrc,iz);
			  pp[VZ]=get_u(p,VZ,ix,iysrc,iz);
			  pp[VY]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,VY,ix,iysrc,iz);

			  
#ifdef MAGNFIELD
#ifdef CORRECTMAGNFIELD
			  pp[B1]=get_u(p,B1,ix,iysrc,iz);
			  pp[B3]=get_u(p,B3,ix,iysrc,iz);
			  pp[B2]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,B2,ix,iysrc,iz);
#endif
#endif


#ifdef RADIATION
			  //rad density
			  pp[EE0]=get_u(p,EE0,ix,iysrc,iz);

#ifdef EVOLVEPHOTONNUMBER
			  //no. of photons
			  pp[NF0]=get_u(p,NF0,ix,iysrc,iz);
#endif

			  //rad velocities
			  pp[FX0]=get_u(p,FX0,ix,iysrc,iz);
			  pp[FZ0]=get_u(p,FZ0,ix,iysrc,iz);
			  pp[FY0]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,FY0,ix,iysrc,iz);

#endif 

			  //update gamma in the corrected cell
#ifdef CONSISTENTGAMMA
			  ldouble newgamma=calc_gammagas(pp,ix,iy,iz);
			  set_u_scalar(gammagas,ix,iy,iz,newgamma);
#endif


			  p2u(pp,uu,&geom);

			  PLOOP(iv)
			  {
			    set_u(p,iv,ix,iy,iz,pp[iv]);  
			    set_u(u,iv,ix,iy,iz,uu[iv]); 
			  }
			}
		    }
#endif
		}
	    }
	}

      //cylindrical like coords
      if (MYCOORDS==CYLCOORDS || MYCOORDS==MCYL1COORDS)
	{
	  //#pragma omp parallel for private(ic,ix,iy,iz,iv,ixsrc) schedule (static)
	  for(iy=0;iy<NY;iy++)
	    {
	      for(iz=0;iz<NZ;iz++)
		{
		  ldouble R,Rsrc,Raxis;
		  ldouble pp[NV],uu[NV];
		  struct geometry geom;

		  //upper
		  Raxis=get_xb(0,0);
		  for(ic=0;ic<nc;ic++)
		    {
		      ix=ic;ixsrc=nc;
		      R=get_x(ix,0);
		      Rsrc=get_x(ixsrc,0);	      
	      	  
		      fill_geometry(ix,iy,iz,&geom);
	  
		      PLOOP(iv)
			pp[iv]=get_u(p,iv,ix,iy,iz);

		      //gas densities
		      pp[RHO]=get_u(p,RHO,ixsrc,iy,iz);
		      pp[UU]=get_u(p,UU,ixsrc,iy,iz);
		      pp[ENTR]=get_u(p,ENTR,ixsrc,iy,iz);		  
#ifdef EVOLVEELECTRONS
		      pp[ENTRE]=get_u(p,ENTRE,ixsrc,iy,iz);		  
		      pp[ENTRI]=get_u(p,ENTRI,ixsrc,iy,iz);		  
#ifdef RELELECTRONS
                      int ie;
                      for(ie=0;ie<NRELBIN;ie++) pp[NEREL(ie)]=get_u(p,NEREL(ie),ixsrc,iy,iz);
#endif
#endif
		      //gas velocities
		      pp[VY]=get_u(p,VY,ixsrc,iy,iz);
		      pp[VZ]=get_u(p,VZ,ixsrc,iy,iz);
		      pp[VX]=fabs((R-Raxis)/(Rsrc-Raxis))*get_u(p,VX,ixsrc,iy,iz);

		      
#ifdef MAGNFIELD
#ifdef CORRECTMAGNFIELD
		      
		      pp[B2]=get_u(p,B2,ixsrc,iy,iz);
		      pp[B3]=get_u(p,B3,ixsrc,iy,iz);
		      pp[B1]=fabs((R-Raxis)/(Rsrc-Raxis))*get_u(p,B1,ixsrc,iy,iz);
		      
#endif
#endif
	

#ifdef RADIATION
		      //rad density
		      pp[EE0]=get_u(p,EE0,ixsrc,iy,iz);

#ifdef EVOLVEPHOTONNUMBER
		      //no. of photons
		      pp[NF0]=get_u(p,NF0,ixsrc,iy,iz);
#endif

		      //rad velocities
		      pp[FY0]=get_u(p,FY0,ixsrc,iy,iz);
		      pp[FZ0]=get_u(p,FZ0,ixsrc,iy,iz);
		      pp[FX0]=fabs((R-Raxis)/(Rsrc-Raxis))*get_u(p,FX0,ixsrc,iy,iz);

#endif 

		      //update gamma in the corrected cell
#ifdef CONSISTENTGAMMA
		      ldouble newgamma=calc_gammagas(pp,ix,iy,iz);
		      set_u_scalar(gammagas,ix,iy,iz,newgamma);
#endif

		      p2u(pp,uu,&geom);

		      PLOOP(iv)
		      {

			set_u(p,iv,ix,iy,iz,pp[iv]);  
			set_u(u,iv,ix,iy,iz,uu[iv]); 
		      }
		    }
		}
	    }
	}
      

    
  return 0; 
}

//**********************************************************************
//treats the most polar cells is a special way, correcting them, not evolving them
//**********************************************************************

int
correct_polaraxis_3d()
{
      int nc=NCCORRECTPOLAR; //correct velocity in nc most polar cells;

      //spherical like coords
      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS ||
	  MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS|| MYCOORDS==MKS3COORDS|| MYCOORDS==TKS3COORDS ||
	  MYCOORDS==MSPH1COORDS || MYCOORDS==JETCOORDS)
	{
          int ix;
//#pragma omp parallel for private(ix)
	  for(ix=0;ix<NX;ix++)
	    {
          int iy,iz,iv,ic,iysrc,ixsrc,gix;
          ldouble ucon[4];
          struct geometry geom,geomBL;
	      fill_geometry_arb(ix,0,0,&geomBL,BLCOORDS);

	      //to avoid amibous VEL4 after 
	      //if(geomBL.xx < 1.*rhorizonBL )      continue;

	      gix=ix+TOI;
	      
	      //overwriting
	      for(iz=0;iz<NZ;iz++)
		{
		  ldouble r,th;
		  ldouble pp[NV],uu[NV];

		  //upper axis
#ifdef MPI
		  if(TJ==0)
#endif
		    {

		      iysrc=nc;
		      for(ic=0;ic<nc;ic++)
			{
			  iy=ic;
	      	 
			  fill_geometry(ix,iy,iz,&geom);
			  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

			  PLOOP(iv) pp[iv]=get_u(p,iv,ix,iysrc,iz);

#ifdef POLARAXISAVGIN3D

			  ldouble r=geomBL.xx;
			  ldouble th=geomBL.yy;
			  ldouble ph=geomBL.zz;


			  ldouble vr,vth,vph,vx,vy,vz;
			  ldouble cosph,sinth,costh,sinph;
			  sinth=sin(th);
			  costh=cos(th);
			  sinph=sin(ph);
			  cosph=cos(ph);
	  
			  //gas
			  pp[RHO]=axis1_primplus[RHO][gix];
			  pp[UU]=axis1_primplus[UU][gix];
			  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);

			  //if(geomBL.xx > 1.*rhorizonBL ) 
			  {
			    //gas velocities
			    vx=axis1_primplus[VX][gix];
			    vy=axis1_primplus[VY][gix];
			    vz=axis1_primplus[VZ][gix];
			    vr =-((-(cosph*sinth*vx) - sinph*sinth*vy - Power(cosph,2)*costh*vz - 
				   costh*Power(sinph,2)*vz)/
				  ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
			    vth = -((-(cosph*costh*vx) - costh*sinph*vy + Power(cosph,2)*sinth*vz + 
				     Power(sinph,2)*sinth*vz)/
				    ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
			    vph = -((sinph*vx - cosph*vy)/(Power(cosph,2) + Power(sinph,2)));

			    //vth /= r;
			    //vph /= r*sinth;
			    vr/=sqrt(geom.gg[1][1]);
			    vth/=sqrt(geom.gg[2][2]);
			    vph/=sqrt(geom.gg[3][3]);

			    ucon[1]=vr; ucon[2]=vth; ucon[3]=vph;
			    /*
			      ldouble xxvec[4],xxvecBL[4];
			      get_xx(ix,iy,iz,xxvec);
			      coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
			      printf("%d > %e %e %e\n",ix,r,th,ph);
			      print_4vector(xxvec);
			      print_4vector(xxvecBL);
			      print_metric(geomBL.gg);
			      print_metric(geom.gg);
			      print_4vector(ucon);
			    */

			    //conv_vels(ucon,ucon,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
			    //trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
			    //conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);
			    //			  print_4vector(ucon);
			    //getch();

			    pp[VX]=ucon[1];
			    pp[VY]=ucon[2];
			    pp[VZ]=ucon[3];
			    //add average rotation
			    pp[VZ]+=axis1_primplus[NV][gix];
			  }
			  //print_primitives(pp);getch();
		     
#ifdef MAGNFIELD
			  //do not overwrite magnetic field, not to break div B=0 there
#endif

#ifdef RADIATION
			  //rad density
			  pp[EE]=axis1_primplus[EE][gix];

#ifdef EVOLVEPHOTONNUMBER
			  //no. of photons
			  pp[NF]=axis1_primplus[NF][gix];
#endif
			  //if(geomBL.xx > 1.*rhorizonBL )
			  {
			    //rad velocities
			    vx=axis1_primplus[FX][gix];
			    vy=axis1_primplus[FY][gix];
			    vz=axis1_primplus[FZ][gix];
			    vr =-((-(cosph*sinth*vx) - sinph*sinth*vy - Power(cosph,2)*costh*vz - 
				   costh*Power(sinph,2)*vz)/
				  ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
			    vth = -((-(cosph*costh*vx) - costh*sinph*vy + Power(cosph,2)*sinth*vz + 
				     Power(sinph,2)*sinth*vz)/
				    ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
			    vph = -((sinph*vx - cosph*vy)/(Power(cosph,2) + Power(sinph,2)));

			    //vth /= r;
			    //vph /= r*sinth;
			    vr/=sqrt(geom.gg[1][1]);
			    vth/=sqrt(geom.gg[2][2]);
			    vph/=sqrt(geom.gg[3][3]);
			    //if(ic==0 && ix==10)
			    //  printf("1 %d > %e %e | %e %e %e\n",iy,vth,sqrt(geom.gg[2][2]),vx,vy,vz);
			    ucon[1]=vr; ucon[2]=vth; ucon[3]=vph;
			    //conv_vels(ucon,ucon,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
			    //trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
			    //conv_vels(ucon,ucon,VEL4,VELPRIMRAD,geom.gg,geom.GG);
			  
			    pp[FX]=ucon[1];
			    pp[FY]=ucon[2];
			    pp[FZ]=ucon[3];
			    
			    //add average rotation
			    pp[FZ]+=axis1_primplus[NV+1][gix];
			  }

#endif //RADIATION
#endif //POLARAXISAVGIN3D

			  //update gamma in the corrected cell
#ifdef CONSISTENTGAMMA
			  ldouble newgamma=calc_gammagas(pp,ix,iy,iz);
			  set_u_scalar(gammagas,ix,iy,iz,newgamma);
#endif

			  p2u(pp,uu,&geom);
			  
			  PLOOP(iv)
			  {
			    if(iv<B1 || iv>B3)
			      {
				set_u(p,iv,ix,iy,iz,pp[iv]);  
				set_u(u,iv,ix,iy,iz,uu[iv]);
			      }
			  }
			}
		    }

		  //bottom axis
#ifdef MPI
		  if(TJ==NTY-1)
#endif
		    {
		      iysrc=NY-1-nc; 
		      for(ic=0;ic<nc;ic++)
			{
			  iy=NY-1-ic;

			  fill_geometry(ix,iy,iz,&geom);
			  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
			  
			  PLOOP(iv) pp[iv]=get_u(p,iv,ix,iysrc,iz);
	      	 
#ifdef POLARAXISAVGIN3D
			  ldouble r=geomBL.xx;
			  ldouble th=geomBL.yy;
			  ldouble ph=geomBL.zz;
			  ldouble vr,vth,vph,vx,vy,vz;
			  ldouble cosph,sinth,costh,sinph;
			  sinth=sin(th);
			  costh=cos(th);
			  sinph=sin(ph);
			  cosph=cos(ph);
	  
			  //gas
			  pp[RHO]=axis2_primplus[RHO][gix];
			  pp[UU]=axis2_primplus[UU][gix];
			  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);

			  //if(geomBL.xx > 1.*rhorizonBL ) 	 
			  {  
			    //gas velocities
			    vx=axis2_primplus[VX][gix];
			    vy=axis2_primplus[VY][gix];
			    vz=axis2_primplus[VZ][gix];
			    vr =-((-(cosph*sinth*vx) - sinph*sinth*vy - Power(cosph,2)*costh*vz - 
				   costh*Power(sinph,2)*vz)/
				  ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
			    vth = -((-(cosph*costh*vx) - costh*sinph*vy + Power(cosph,2)*sinth*vz + 
				     Power(sinph,2)*sinth*vz)/
				    ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
			    vph = -((sinph*vx - cosph*vy)/(Power(cosph,2) + Power(sinph,2)));

			    //vth /= r;
			    //vph /= r*sinth;
			    vr/=sqrt(geom.gg[1][1]);
			    vth/=sqrt(geom.gg[2][2]);
			    //if(ic==0 && ix==10)
			    // printf("2 %d > %e %e | %e %e %e\n",iy,vth,sqrt(geom.gg[2][2]),vx,vy,vz);
			    vph/=sqrt(geom.gg[3][3]);

			    ucon[1]=vr; ucon[2]=vth; ucon[3]=vph;
			    //conv_vels(ucon,ucon,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
			    //trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
			    //conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);

			    pp[VX]=ucon[1];
			    pp[VY]=ucon[2];
			    pp[VZ]=ucon[3];
			    
			    //add average rotation
			    pp[VZ]+=axis2_primplus[NV][gix];
			  }
#ifdef MAGNFIELD
			  //do not overwrite magnetic field, not to break div B=0 there
#endif

#ifdef RADIATION
			  //rad density
			  pp[EE]=axis2_primplus[EE][gix];

#ifdef EVOLVEPHOTONNUMBER
			  //no. of photons
			  pp[NF]=axis2_primplus[NF][gix];
#endif
			  //if(geomBL.xx > 1.*rhorizonBL ) 	
			  {  
			    //rad velocities
			    vx=axis2_primplus[FX][gix];
			    vy=axis2_primplus[FY][gix];
			    vz=axis2_primplus[FZ][gix];
			    vr =-((-(cosph*sinth*vx) - sinph*sinth*vy - Power(cosph,2)*costh*vz - 
				   costh*Power(sinph,2)*vz)/
				  ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
			    vth = -((-(cosph*costh*vx) - costh*sinph*vy + Power(cosph,2)*sinth*vz + 
				     Power(sinph,2)*sinth*vz)/
				    ((Power(cosph,2) + Power(sinph,2))*(Power(costh,2) + Power(sinth,2))));
			    vph = -((sinph*vx - cosph*vy)/(Power(cosph,2) + Power(sinph,2)));

			    //vth /= r;
			    //vph /= r*sinth;
			    vr/=sqrt(geom.gg[1][1]);
			    vth/=sqrt(geom.gg[2][2]);
			    vph/=sqrt(geom.gg[3][3]);

			    ucon[1]=vr; ucon[2]=vth; ucon[3]=vph;
			    //conv_vels(ucon,ucon,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
			    //trans2_coco(geomBL.xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
			    //conv_vels(ucon,ucon,VEL4,VELPRIMRAD,geom.gg,geom.GG);

			    pp[FX]=ucon[1];
			    pp[FY]=ucon[2];
			    pp[FZ]=ucon[3];
			    
			    //add average rotation
			    pp[FZ]+=axis2_primplus[NV+1][gix];
			  }
#endif //RADIATION
#endif //POLARAXISAVGIN3D

			  //update gamma in the corrected cell
#ifdef CONSISTENTGAMMA
			  ldouble newgamma=calc_gammagas(pp,ix,iy,iz);
			  set_u_scalar(gammagas,ix,iy,iz,newgamma);
#endif

			  p2u(pp,uu,&geom);

			  PLOOP(iv)
			  {
			    if(iv<B1 || iv>B3)
			      {
				set_u(p,iv,ix,iy,iz,pp[iv]);  
				set_u(u,iv,ix,iy,iz,uu[iv]);
			      }
			  }
			}
		    }

		}
	     
	    }
	}

  return 0; 
}

//**********************************************************************
// say if given cell is within NCCORRECTPOLAR from axis */
//**********************************************************************

int
is_cell_corrected_polaraxis(int ix, int iy, int iz)
{

#if defined(CORRECT_POLARAXIS) || defined(CORRECT_POLARAXIS_3D)
#ifdef MPI
  if(TJ==0) //tile
#endif
    if(iy<NCCORRECTPOLAR) 
      return 1;
#ifndef HALFTHETA
#ifdef MPI
  if(TJ==NTY-1) //tile
#endif   
    if(iy>(NY-NCCORRECTPOLAR-1))
      return 1;
#endif
#endif
  
  return 0;
}

//**********************************************************************
/* do we skip implicit for some reason*/
//**********************************************************************

int
skip_cell_implicit(int ix, int iy, int iz)
{

  int skip=0;

#if defined(SKIPIMPLICIT_ATBH) || defined(SKIPIMPLICIT_HIGHSIGMA)
  struct geometry geom;
  
  fill_geometry(ix,iy,iz,&geom);
  
#ifdef SKIPIMPLICIT_ATBH
  ldouble xxBL[4];

  #ifdef PRECOMPUTE_MY2OUT
  get_xxout(geom.ix, geom.iy, geom.iz, xxBL);
  #else
  coco_N(geom.xxvec,xxBL,MYCOORDS,BLCOORDS);
  #endif

  if(xxBL[1]<rhorizonBL) return 1;
#endif

#ifdef SKIPIMPLICIT_HIGHSIGMA
  int iv;
  ldouble pp[NV];
 
  PLOOP(iv) pp[iv]=get_u(p, iv, ix, iy, iz);

  ldouble ucon[4],ucov[4];
  ldouble bcon[4],bcov[4];
  ldouble bsq;
  calc_ucon_ucov_from_prims(pp, &geom, ucon, ucov);
  calc_bcon_bcov_bsq_from_4vel(pp, ucon, ucov, &geom, bcon, bcov, &bsq);

  ldouble sigma=bsq/pp[RHO];

  if(sigma>SKIPIMPLICIT_HIGHSIGMA) return 1;
#endif
#endif
  
  return skip;
}

//**********************************************************************
/* say if given cell is evolved or rather corrected */
//**********************************************************************

int
is_cell_active(int ix, int iy, int iz)
{
  //by default ALWAYS active
  //ANDREW TODO: use for AMR? 
  return 1;
}
//**********************************************************************
//evolves entropies using fluxes at faces to calculate mixing rates
//**********************************************************************

int
get_factors_entropies_following_gas(int ix,int iy,int iz,ldouble *f0,
				    ldouble *fxl, ldouble *fxr,
				    ldouble* fyl, ldouble *fyr,
				    ldouble *fzl, ldouble *fzr,
				    ldouble dt, int iv)
{
      ldouble dx=get_size_x(ix,0);
      ldouble dy=get_size_x(iy,1);
      ldouble dz=get_size_x(iz,2);
      
      ldouble flxr,flyr,flzr,flxl,flyl,flzl;
      ldouble uinit,ufinal,duxl,duxr,duyl,duyr,duzl,duzr;
      
      flxl=get_ub(flbx,iv,ix,iy,iz,0);
      flxr=get_ub(flbx,iv,ix+1,iy,iz,0);
      flyl=get_ub(flby,iv,ix,iy,iz,1);
      flyr=get_ub(flby,iv,ix,iy+1,iz,1);
      flzl=get_ub(flbz,iv,ix,iy,iz,2);
      flzr=get_ub(flbz,iv,ix,iy,iz+1,2);

      //unsplit scheme
      //du=-(flxr*mxr-flxl*mxl)*dt/dx - (flyr*myr-flyl*myl)*dt/dy - (flzr*mzr-flzl*mzl)*dt/dz;

      duxl=flxl*dt/dx;
      duxr=-flxr*dt/dx;
      duyl=flyl*dt/dy;
      duyr=-flyr*dt/dy;
      duzl=flzl*dt/dz;
      duzr=-flzr*dt/dz;
      
      uinit=get_u(upreexplicit,iv,ix,iy,iz);
      ufinal = uinit + duxl + duxr + duyl + duyr + duzl + duzr;
  
      *f0 = uinit / ufinal;
      *fxl = duxl / ufinal;
      *fxr = duxr / ufinal;
      *fyl = duyl / ufinal;
      *fyr = duyr / ufinal;
      *fzl = duzl / ufinal;
      *fzr = duzr / ufinal;

      /*
      //negative contributions - lost mass
      *f0 = (uinit + my_min(0., duxl) + my_min(0., duxr)
	    + my_min(0., duyl) + my_min(0., duyr)
	    + my_min(0., duzl) + my_min(0., duzr)) / ufinal;

      //positive contributions - gained mass
      *fxl = my_max(0., duxl) / ufinal;
      *fxr = my_max(0., duxr) / ufinal;
      *fyl = my_max(0., duyl) / ufinal;
      *fyr = my_max(0., duyr) / ufinal;
      *fzl = my_max(0., duzl) / ufinal;
      *fzr = my_max(0., duzr) / ufinal;
      */

      return 0;
}


//**********************************************************************
// Mixing entropy adjustment
// TODO: still A LOT of work here
// ANDREW - made compatible with S4 entropy (02/15/16)
//**********************************************************************

int
mix_entropies(ldouble dt)
{     
  int ix,iy,iz,ii;

 #pragma omp parallel for private(ii,ix,iy,iz) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2];

      ldouble gamma=GAMMA;
      #ifdef CONSISTENTGAMMA
      gamma=pick_gammagas(ix,iy,iz);
      #endif
   
      ldouble gammai,gammae,ginv,gm1;
      ldouble sinit,kappainit,scalefac;      
      ldouble f0,fxl,fxr,fyl,fyr,fzl,fzr;
      ldouble sfinal,Sfinal,s0,sxl,sxr,syl,syr,szl,szr;
      ldouble rsxl,rsxr,rsyl,rsyr,rszl,rszr;
      ldouble u0,uxl,uxr,uyl,uyr,uzl,uzr;
      ldouble p0,pxl,pxr,pyl,pyr,pzl,pzr;
      ldouble ufinal;
      ldouble Sbak,Sbake,Sbaki;
      ldouble gfac, entsum, Tefinal, Tifinal,Tgasfinal;
      ldouble s_scale;
      ldouble pintfinalpre ,pintfinalxl,pintfinalxr;
      ldouble pintfinalyl,pintfinalyr ,pintfinalzl,pintfinalzr;

      ldouble uintfinalpre ,uintfinalxl,uintfinalxr;
      ldouble uintfinalyl,uintfinalyr ,uintfinalzl,uintfinalzr;
      
      ldouble rhoinit = get_u(ppreexplicit, RHO, ix,iy,iz);
      ldouble ninit = one_over_mugas_mp * rhoinit;
      ldouble neinit = calc_thermal_ne(&get_u(ppreexplicit,0,ix,iy,iz));
      ldouble niinit = rhoinit / (MU_I*M_PROTON);
    
      ldouble rhofinal = get_u(p,RHO,ix,iy,iz);
      ldouble nfinal = rhofinal * one_over_mugas_mp;
      ldouble nefinal = calc_thermal_ne(&get_u(p,0,ix,iy,iz));
      ldouble nifinal = rhofinal/(MU_I*M_PROTON);
    
      //**************************************
      //get fractions of the entropies coming from each cell
      //**************************************
      
      get_factors_entropies_following_gas(ix,iy,iz,&f0,&fxl,&fxr,&fyl,&fyr,&fzl,&fzr,dt,RHO);     

      //***********************************
      //******** total gas entropy ********
      //***********************************
      
#ifndef DONOTMIXGASENTROPY //if entropy ~ p/rho^gamma, then no problem with mixing
      
      //initial state temperature and gamma
#ifdef MIXENTROPIES_CONST_PRESSURE
      ldouble Tgasinit=calc_TfromS(get_u(ppreexplicit,ENTR,ix,iy,iz), rhoinit, ix, iy, iz);
      sinit = get_u(ppreexplicit,ENTR,ix,iy,iz) / rhoinit; //entropy per particle of initial state
      kappainit = (ninit*K_BOLTZ*Tgasinit) / pow(rhoinit, gamma);

      //s_scale= (MU_GAS * M_PROTON);
      ginv = 1./gamma;
      gm1 = gamma-1.;
      gfac = gm1;      
    
      #ifdef NOLOGINS
      scalefac = kappainit / (sinit*gfac); 
      #else
      scalefac = kappainit / (exp(sinit*gfac));
      #endif
#endif //MIXENTROPIES_CONST_PRESSURE
      

      //calculates entropy per particle
      //assumes entropy flux is (n gdet s u^i) = (rho / mu / M_PROTON gdet s u^i) \propto (rho gdet s u^i)
      //what is below is actually  mu / M_PROTON s but rho will return later 
      s0=get_u(upreexplicit,ENTR,ix,iy,iz)/get_u(upreexplicit,RHO,ix,iy,iz);
      
      //based on ratio of fluxes (which may misbehave)
      sxl=get_ub(flbx,ENTR,ix,iy,iz,0)/get_ub(flbx,RHO,ix,iy,iz,0);
      sxr=get_ub(flbx,ENTR,ix+1,iy,iz,0)/get_ub(flbx,RHO,ix+1,iy,iz,0);
      syl=get_ub(flby,ENTR,ix,iy,iz,1)/get_ub(flby,RHO,ix,iy,iz,1);
      syr=get_ub(flby,ENTR,ix,iy+1,iz,1)/get_ub(flby,RHO,ix,iy+1,iz,1);
      szl=get_ub(flbz,ENTR,ix,iy,iz,2)/get_ub(flbz,RHO,ix,iy,iz,2);
      szr=get_ub(flbz,ENTR,ix,iy,iz+1,2)/get_ub(flbz,RHO,ix,iy,iz+1,2);
    
      
#if defined(UPWINDENTROPYMIXING) || defined(PARTLYUPWINDENTROPYMIXING)
      //upwind scheme
      ldouble rs_c=get_u(ppreexplicit,ENTR,ix,iy,iz)/get_u(ppreexplicit,RHO,ix,iy,iz);
      if(fxl<0) rsxl=rs_c;
      else rsxl=get_u(ppreexplicit,ENTR,ix-1,iy,iz)/get_u(ppreexplicit,RHO,ix-1,iy,iz) ;
      if(fxr<0) rsxr=rs_c;
      else rsxr=get_u(ppreexplicit,ENTR,ix+1,iy,iz)/get_u(ppreexplicit,RHO,ix+1,iy,iz) ;
      if(fyl<0) rsyl=rs_c;
      else rsyl=get_u(ppreexplicit,ENTR,ix,iy-1,iz)/get_u(ppreexplicit,RHO,ix,iy-1,iz) ;
      if(fyr<0) rsyr=rs_c;
      else rsyr=get_u(ppreexplicit,ENTR,ix,iy+1,iz)/get_u(ppreexplicit,RHO,ix,iy+1,iz) ;
      if(fzl<0) rszl=rs_c;
      else rszl=get_u(ppreexplicit,ENTR,ix,iy,iz-1)/get_u(ppreexplicit,RHO,ix,iy,iz-1) ;
      if(fzr<0) rszr=rs_c;
      else rszr=get_u(ppreexplicit,ENTR,ix,iy,iz+1)/get_u(ppreexplicit,RHO,ix,iy,iz+1) ;

      #ifdef UPWINDENTROPYMIXING
      sxl=rsxl;sxr=rsxr;syl=rsyl;syr=rsyr;szl=rszl;szr=rszr;
      #endif
      
      #ifdef PARTLYUPWINDENTROPYMIXING //corrects to cell center when outflowing
      if(fxl<0) sxl=rsxl;
      if(fxr<0) sxr=rsxr;
      if(fyl<0) syl=rsyl;
      if(fyr<0) syr=rsyr;
      if(fzl<0) szl=rszl;
      if(fzr<0) szr=rszr;
      #endif

#elif defined(PARTLYRECONSTRUCTEDUPWINDENTROPYMIXING)
      //upwind scheme
      if(fxl<0) rsxl=get_ub(pbRx,ENTR,ix,iy,iz,0)/get_ub(pbRx,RHO,ix,iy,iz,0);
      else rsxl=get_ub(pbLx,ENTR,ix,iy,iz,0)/get_ub(pbLx,RHO,ix,iy,iz,0);
      if(fxr<0) rsxr=get_ub(pbLx,ENTR,ix+1,iy,iz,0)/get_ub(pbLx,RHO,ix+1,iy,iz,0);
      else rsxr=get_ub(pbRx,ENTR,ix+1,iy,iz,0)/get_ub(pbRx,RHO,ix+1,iy,iz,0);
      if(fyl<0) rsyl=get_ub(pbRx,ENTR,ix,iy,iz,1)/get_ub(pbRx,RHO,ix,iy,iz,1);
      else rsyl=get_ub(pbLx,ENTR,ix,iy,iz,1)/get_ub(pbLx,RHO,ix,iy,iz,1);
      if(fyr<0) rsyr=get_ub(pbLx,ENTR,ix,iy+1,iz,1)/get_ub(pbLx,RHO,ix,iy+1,iz,1);
      else rsyr=get_ub(pbRx,ENTR,ix,iy+1,iz,1)/get_ub(pbRx,RHO,ix,iy+1,iz,1);
      if(fzl<0) rszl=get_ub(pbRx,ENTR,ix,iy,iz,2)/get_ub(pbRx,RHO,ix,iy,iz,2);
      else rszl=get_ub(pbLx,ENTR,ix,iy,iz,2)/get_ub(pbLx,RHO,ix,iy,iz,2);
      if(fzr<0) rszr=get_ub(pbLx,ENTR,ix,iy,iz+1,2)/get_ub(pbLx,RHO,ix,iy,iz+1,2);
      else rszr=get_ub(pbRx,ENTR,ix,iy,iz+1,2)/get_ub(pbRx,RHO,ix,iy,iz+1,2);
      
      sxl=rsxl;sxr=rsxr;syl=rsyl;syr=rsyr;szl=rszl;szr=rszr;
#endif
      

      if(!isfinite(sxl)) sxl=s0;
      if(!isfinite(sxr)) sxr=s0;
      if(!isfinite(syl)) syl=s0;
      if(!isfinite(syr)) syr=s0;
      if(!isfinite(szl)) szl=s0;
      if(!isfinite(szr)) szr=s0;

      //reduces to regular Godunov explicit     
      Sbak =  rhofinal*(f0*s0 + fxl*sxl + fxr*sxr + fyl*syl + fyr*syr + fzl*szl + fzr*szr); 
            
#if defined(OLDENTROPYMIXING) //does not correct for entropy of mixing
      Sfinal = Sbak;
      
#elif defined(MIXENTROPIES_CONST_PRESSURE)       //ANDREW constant pressure mixing     
      
      ldouble rhofinalg = pow(rhofinal, gamma);
      #ifdef NOLOGINS
      p0  = s0 *gfac*scalefac*rhofinalg;
      pxl = sxl*gfac*scalefac*rhofinalg;
      pxr = sxr*gfac*scalefac*rhofinalg;
      pyl = syl*gfac*scalefac*rhofinalg;
      pyr = syr*gfac*scalefac*rhofinalg;
      pzl = szl*gfac*scalefac*rhofinalg;
      pzr = szr*gfac*scalefac*rhofinalg;
      #else
      p0  = exp(s0 *gfac)*scalefac*rhofinalg;
      pxl = exp(sxl*gfac)*scalefac*rhofinalg;
      pxr = exp(sxr*gfac)*scalefac*rhofinalg;
      pyl = exp(syl*gfac)*scalefac*rhofinalg;
      pyr = exp(syr*gfac)*scalefac*rhofinalg;
      pzl = exp(szl*gfac)*scalefac*rhofinalg;
      pzr = exp(szr*gfac)*scalefac*rhofinalg;
      #endif

      #ifndef DONOTLIMITENTRINMIXING
      pintfinalpre = gm1*get_u(ppreexplicit,UU,ix,iy,iz);
      pintfinalxl =  gm1* get_u(ppreexplicit,UU,ix-1,iy,iz);
      pintfinalxr =  gm1* get_u(ppreexplicit,UU,ix+1,iy,iz);
      pintfinalyl =  gm1* get_u(ppreexplicit,UU,ix,iy-1,iz);
      pintfinalyr =  gm1* get_u(ppreexplicit,UU,ix,iy+1,iz);
      pintfinalzl =  gm1* get_u(ppreexplicit,UU,ix,iy,iz-1);
      pintfinalzr =  gm1* get_u(ppreexplicit,UU,ix,iy,iz+1);

      p0  = my_min(p0, LIMITFACTORINMIXING*pintfinalpre);
      pxl = my_min(pxl, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalxl));
      pxr = my_min(pxr, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalxr));
      pyl = my_min(pyl, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalyl));
      pyr = my_min(pyr, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalyr));
      pzl = my_min(pzl, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalzl));
      pzr = my_min(pzr, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalzr));
      #endif

      entsum = f0*pow(p0,ginv) + fxl*pow(pxl,ginv) + fxr*pow(pxr,ginv) + fyl*pow(pyl,ginv) + fyr*pow(pyr,ginv) + fzl*pow(pzl,ginv) + fzr*pow(pzr,ginv);
      Tgasfinal = pow(entsum, gamma) / (nfinal * K_BOLTZ);
      Sfinal = calc_SfromT(rhofinal, Tgasfinal, ix, iy, iz);
      
#else //constant density
      //accounts for adiabatic expansion/compression independently for each component, only then mixes

      #ifdef DONOTLIMITENTRINMIXING
      u0=calc_ufromS(s0*rhofinal,rhofinal,ix,iy,iz);
      uxl=calc_ufromS(sxl*rhofinal,rhofinal,ix,iy,iz);
      uxr=calc_ufromS(sxr*rhofinal,rhofinal,ix,iy,iz);
      uyl=calc_ufromS(syl*rhofinal,rhofinal,ix,iy,iz);
      uyr=calc_ufromS(syr*rhofinal,rhofinal,ix,iy,iz);
      uzl=calc_ufromS(szl*rhofinal,rhofinal,ix,iy,iz);
      uzr=calc_ufromS(szr*rhofinal,rhofinal,ix,iy,iz);

      #else
      uintfinalpre = get_u(ppreexplicit,UU,ix,iy,iz);
      uintfinalxl =  get_u(ppreexplicit,UU,ix-1,iy,iz);
      uintfinalxr =  get_u(ppreexplicit,UU,ix+1,iy,iz);
      uintfinalyl =  get_u(ppreexplicit,UU,ix,iy-1,iz);
      uintfinalyr =  get_u(ppreexplicit,UU,ix,iy+1,iz);
      uintfinalzl =  get_u(ppreexplicit,UU,ix,iy,iz-1);
      uintfinalzr =  get_u(ppreexplicit,UU,ix,iy,iz+1);

      u0 =my_min(calc_ufromS(s0*rhofinal,rhofinal,ix,iy,iz),LIMITFACTORINMIXING*uintfinalpre);
      uxl=my_min(calc_ufromS(sxl*rhofinal,rhofinal,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalxl));
      uxr=my_min(calc_ufromS(sxr*rhofinal,rhofinal,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalxr));
      uyl=my_min(calc_ufromS(syl*rhofinal,rhofinal,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalyl));
      uyr=my_min(calc_ufromS(syr*rhofinal,rhofinal,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalyr));
      uzl=my_min(calc_ufromS(szl*rhofinal,rhofinal,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalzl));
      uzr=my_min(calc_ufromS(szr*rhofinal,rhofinal,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalzr));
      #endif
      
      ufinal=f0*u0+fxl*uxl+fxr*uxr+fyl*uyl+fyr*uyr+fzl*uzl+fzr*uzr;
      Sfinal=calc_Sfromu(rhofinal,ufinal,ix,iy,iz);
      
#endif
       
       if(!isfinite(Sfinal))
	 {
	   //Godunov backup
	   Sfinal=Sbak;
	 }
       set_u(p,ENTR,ix,iy,iz,Sfinal);

#endif //DONOTMIXGASENTROPY
     
#ifdef EVOLVEELECTRONS
#ifndef DONOTMIXSPECIESENTROPY

      //*****************************************      
      //******** electrons **********************
      //*****************************************

      //initial state temperature and gamma
#ifdef MIXENTROPIES_CONST_PRESSURE
      ldouble rhoeinit=MU_E*M_PROTON*neinit;      
      ldouble Teinit=calc_TfromSerho(get_u(ppreexplicit,ENTRE,ix,iy,iz), rhoeinit, ELECTRONS, ix, iy, iz);
      sinit = get_u(ppreexplicit,ENTRE,ix,iy,iz) / neinit; //entropy per particle of initial state

      //assume gammae constant over initial state
      #if defined(CONSISTENTGAMMA) && !defined(FIXEDGAMMASPECIES)
      gammae = calc_gammaintfromtemp(Teinit, ELECTRONS);
      #else
      gammae = GAMMAE;
      #endif

      //find the scale factor that relates kappa to exp(s)
      // ANDREW ONLY WORKS FOR S2,S3,S4 (not original S1)
      ginv = 1./gammae;
      gm1 = gammae-1.;
      kappainit = (neinit*K_BOLTZ*Teinit) / pow(rhoeinit, gammae);
      
      #if !defined(CONSISTENTGAMMA) && defined(NOLOGINS2)
      gfac = (gammae-1.);      
      scalefac = kappainit / (sinit*gfac); 
      #else
      gfac = (gammae-1.)/K_BOLTZ;      
      scalefac = kappainit / (exp(sinit*gfac));
      #endif

#endif //MIXENTROPIES_CONST_PRESSURE
      
      //entropy flux is (n_thermal gdet s_per_particle u^i)
      s0 =get_u(upreexplicit,ENTRE,ix,iy,iz) / calc_thermal_ne(&get_u(upreexplicit,0,ix,iy,iz));
      sxl=get_ub(flbx,ENTRE,ix,iy,iz,0)   / calc_thermal_ne((get_ub_ptr(flbx,0,ix,iy,iz,0))); 
      sxr=get_ub(flbx,ENTRE,ix+1,iy,iz,0) / calc_thermal_ne((get_ub_ptr(flbx,0,ix+1,iy,iz,0)));
      syl=get_ub(flby,ENTRE,ix,iy,iz,1)   / calc_thermal_ne((get_ub_ptr(flby,0,ix,iy,iz,1)));
      syr=get_ub(flby,ENTRE,ix,iy+1,iz,1) / calc_thermal_ne((get_ub_ptr(flby,0,ix,iy+1,iz,1)));
      szl=get_ub(flbz,ENTRE,ix,iy,iz,2)   / calc_thermal_ne((get_ub_ptr(flbz,0,ix,iy,iz,2)));
      szr=get_ub(flbz,ENTRE,ix,iy,iz+1,2) / calc_thermal_ne((get_ub_ptr(flbz,0,ix,iy,iz+1,2)));
      
      //upwind schemes
#if defined(UPWINDENTROPYMIXING) || defined(PARTLYUPWINDENTROPYMIXING)
      ldouble rs_0=get_u(ppreexplicit,ENTRE,ix,iy,iz)/calc_thermal_ne(&get_u(ppreexplicit,0,ix,iy,iz));
      if(fxl<0) rsxl=rs_0;
      else rsxl=get_u(ppreexplicit,ENTRE,ix-1,iy,iz)/calc_thermal_ne(&get_u(ppreexplicit,0,ix-1,iy,iz)) ;
      if(fxr<0) rsxr=rs_0;
      else rsxr=get_u(ppreexplicit,ENTRE,ix+1,iy,iz)/calc_thermal_ne(&get_u(ppreexplicit,0,ix+1,iy,iz)) ;
      if(fyl<0) rsyl=rs_0;
      else rsyl=get_u(ppreexplicit,ENTRE,ix,iy-1,iz)/calc_thermal_ne(&get_u(ppreexplicit,0,ix,iy-1,iz)) ;
      if(fyr<0) rsyr=rs_0;
      else rsyr=get_u(ppreexplicit,ENTRE,ix,iy+1,iz)/calc_thermal_ne(&get_u(ppreexplicit,0,ix,iy+1,iz)) ;
      if(fzl<0) rszl=rs_0;
      else rszl=get_u(ppreexplicit,ENTRE,ix,iy,iz-1)/calc_thermal_ne(&get_u(ppreexplicit,0,ix,iy,iz-1)) ;
      if(fzr<0) rszr=rs_0;
      else rszr=get_u(ppreexplicit,ENTRE,ix,iy,iz+1)/calc_thermal_ne(&get_u(ppreexplicit,0,ix,iy,iz+1)) ;

      #ifdef UPWINDENTROPYMIXING //always corrects to cell center      
      sxl=rsxl;
      sxr=rsxr;
      syl=rsyl;
      syr=rsyr;
      szl=rszl;
      szr=rszr;
      #endif

      #ifdef PARTLYUPWINDENTROPYMIXING //only corrects to cell center when outflowing
      if(fxl<0) sxl=rsxl;
      if(fxr<0) sxr=rsxr;
      if(fyl<0) syl=rsyl;
      if(fyr<0) syr=rsyr;
      if(fzl<0) szl=rszl;
      if(fzr<0) szr=rszr;
      #endif
 
#elif defined(PARTLYRECONSTRUCTEDUPWINDENTROPYMIXING)
      if(fxl<0) rsxl=get_ub(pbRx,ENTRE,ix,iy,iz,0) / calc_thermal_ne(get_ub_ptr(pbRx,0,ix,iy,iz,0));
      else rsxl=get_ub(pbLx,ENTRE,ix,iy,iz,0) / calc_thermal_ne(get_ub_ptr(pbLx,0,ix,iy,iz,0));
      if(fxr<0) rsxr=get_ub(pbLx,ENTRE,ix+1,iy,iz,0) / calc_thermal_ne(get_ub_ptr(pbLx,0,ix+1,iy,iz,0));
      else rsxr=get_ub(pbRx,ENTRE,ix+1,iy,iz,0) / calc_thermal_ne(get_ub_ptr(pbRx,0,ix+1,iy,iz,0));
      if(fyl<0) rsyl=get_ub(pbRx,ENTRE,ix,iy,iz,1) / calc_thermal_ne(get_ub_ptr(pbRx,0,ix,iy,iz,1));
      else rsyl=get_ub(pbLx,ENTRE,ix,iy,iz,1) / calc_thermal_ne(get_ub_ptr(pbLx,0,ix,iy,iz,1));
      if(fyr<0) rsyr=get_ub(pbLx,ENTRE,ix,iy+1,iz,1) / calc_thermal_ne(get_ub_ptr(pbLx,0,ix,iy+1,iz,1));
      else rsyr=get_ub(pbRx,ENTRE,ix,iy+1,iz,1)/calc_thermal_ne(get_ub_ptr(pbRx,0,ix,iy+1,iz,1));
      if(fzl<0) rszl=get_ub(pbRx,ENTRE,ix,iy,iz,2)/calc_thermal_ne(get_ub_ptr(pbRx,0,ix,iy,iz,2));
      else rszl=get_ub(pbLx,ENTRE,ix,iy,iz,2)/calc_thermal_ne(get_ub_ptr(pbLx,0,ix,iy,iz,2));
      if(fzr<0) rszr=get_ub(pbLx,ENTRE,ix,iy,iz+1,2)/calc_thermal_ne(get_ub_ptr(pbLx,0,ix,iy,iz+1,2));
      else rszr=get_ub(pbRx,ENTRE,ix,iy,iz+1,2)/calc_thermal_ne(get_ub_ptr(pbRx,0,ix,iy,iz+1,2));
      
      sxl=rsxl;sxr=rsxr;syl=rsyl;syr=rsyr;szl=rszl;szr=rszr;
#endif

      //floors
      if(!isfinite(sxl)) sxl=s0;
      if(!isfinite(sxr)) sxr=s0;
      if(!isfinite(syl)) syl=s0;
      if(!isfinite(syr)) syr=s0;
      if(!isfinite(szl)) szl=s0;
      if(!isfinite(szr)) szr=s0;

      //old entropy mixing -- reproduces Finite Volume (with possible upwind corrections)      
      Sbake = nefinal*(f0*s0 + fxl*sxl + fxr*sxr + fyl*syl + fyr*syr + fzl*szl + fzr*szr);
      
#if defined(OLDENTROPYMIXING)
      Sfinal = Sbake;
      
#elif defined(MIXENTROPIES_CONST_PRESSURE)    //ANDREW constant pressure mixing    
      
    
      ldouble rhoefinal = MU_E*M_PROTON*nefinal;
      ldouble rhoefinalge = pow(rhoefinal, gammae);
      #if !defined(CONSISTENTGAMMA) && defined(NOLOGINS2)      
      p0  = s0 *gfac*scalefac*rhoefinalge;
      pxl = sxl*gfac*scalefac*rhoefinalge;
      pxr = sxr*gfac*scalefac*rhoefinalge;
      pyl = syl*gfac*scalefac*rhoefinalge;
      pyr = syr*gfac*scalefac*rhoefinalge;
      pzl = szl*gfac*scalefac*rhoefinalge;
      pzr = szr*gfac*scalefac*rhoefinalge;
      #else
      p0  = exp(s0 *gfac)*scalefac*rhoefinalge;
      pxl = exp(sxl*gfac)*scalefac*rhoefinalge;
      pxr = exp(sxr*gfac)*scalefac*rhoefinalge;
      pyl = exp(syl*gfac)*scalefac*rhoefinalge;
      pyr = exp(syr*gfac)*scalefac*rhoefinalge;
      pzl = exp(szl*gfac)*scalefac*rhoefinalge;
      pzr = exp(szr*gfac)*scalefac*rhoefinalge;
      #endif
      
      #ifndef DONOTLIMITENTRINMIXING
      //uses full rho, but since it's a limiter only it should be ok
      
      pintfinalpre = gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy,iz),  get_u(ppreexplicit,RHO,ix,iy,iz),ELECTRONS,ix,iy,iz);
      pintfinalxl =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix-1,iy,iz),get_u(ppreexplicit,RHO,ix-1,iy,iz),ELECTRONS,ix-1,iy,iz);
      pintfinalxr =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix+1,iy,iz),get_u(ppreexplicit,RHO,ix+1,iy,iz),ELECTRONS,ix+1,iy,iz);
      pintfinalyl =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy-1,iz),get_u(ppreexplicit,RHO,ix,iy-1,iz),ELECTRONS,ix,iy-1,iz);
      pintfinalyr =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy+1,iz),get_u(ppreexplicit,RHO,ix,iy+1,iz),ELECTRONS,ix,iy+1,iz);
      pintfinalzl =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy,iz-1),get_u(ppreexplicit,RHO,ix,iy,iz-1),ELECTRONS,ix,iy,iz-1);
      pintfinalzr =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy,iz+1),get_u(ppreexplicit,RHO,ix,iy,iz+1),ELECTRONS,ix,iy,iz+1);
      
      /*
      pintfinalpre =  gm1*get_u(ppreexplicit,UU,ix,iy,iz);
      pintfinalxl =  gm1* get_u(ppreexplicit,UU,ix-1,iy,iz);
      pintfinalxr =  gm1* get_u(ppreexplicit,UU,ix+1,iy,iz);
      pintfinalyl =  gm1* get_u(ppreexplicit,UU,ix,iy-1,iz);
      pintfinalyr =  gm1* get_u(ppreexplicit,UU,ix,iy+1,iz);
      pintfinalzl =  gm1* get_u(ppreexplicit,UU,ix,iy,iz-1);
      pintfinalzr =  gm1*get_u(ppreexplicit,UU,ix,iy,iz+1);
      */
      
      p0  = my_min(p0,  LIMITFACTORINMIXING*pintfinalpre);
      pxl = my_min(pxl, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalxl));
      pxr = my_min(pxr, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalxr));
      pyl = my_min(pyl, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalyl));
      pyr = my_min(pyr, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalyr));
      pzl = my_min(pzl, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalzl));
      pzr = my_min(pzr, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalzr));
      #endif

      entsum = f0*pow(p0,ginv)+fxl*pow(pxl,ginv)+fxr*pow(pxr,ginv)+fyl*pow(pyl,ginv)+fyr*pow(pyr,ginv)+fzl*pow(pzl,ginv)+fzr*pow(pzr,ginv);
      Tefinal = pow(entsum, gammae) / (nefinal * K_BOLTZ);
      Sfinal=calc_SefromrhoT(rhoefinal, Tefinal, ELECTRONS);
      
#else //constant density
      //accounts for adiabatic expansion/compression independently for each component, only then mixes
    
      ldouble rhoefinal = MU_E*M_PROTON*nefinal;
      #ifdef DONOTLIMITENTRINMIXING
      u0 =calc_ufromSerho(s0 *nefinal,rhoefinal,ELECTRONS,ix,iy,iz);
      uxl=calc_ufromSerho(sxl*nefinal,rhoefinal,ELECTRONS,ix,iy,iz);
      uxr=calc_ufromSerho(sxr*nefinal,rhoefinal,ELECTRONS,ix,iy,iz);
      uyl=calc_ufromSerho(syl*nefinal,rhoefinal,ELECTRONS,ix,iy,iz);
      uyr=calc_ufromSerho(syr*nefinal,rhoefinal,ELECTRONS,ix,iy,iz);
      uzl=calc_ufromSerho(szl*nefinal,rhoefinal,ELECTRONS,ix,iy,iz);
      uzr=calc_ufromSerho(szr*nefinal,rhoefinal,ELECTRONS,ix,iy,iz);
      #else

      //uses full rho, but since it's a limiter only it should be ok
      uintfinalpre = get_u(ppreexplicit,UU,ix,iy,iz);
      uintfinalxl =  get_u(ppreexplicit,UU,ix-1,iy,iz);
      uintfinalxr =  get_u(ppreexplicit,UU,ix+1,iy,iz);
      uintfinalyl =  get_u(ppreexplicit,UU,ix,iy-1,iz);
      uintfinalyr =  get_u(ppreexplicit,UU,ix,iy+1,iz);
      uintfinalzl =  get_u(ppreexplicit,UU,ix,iy,iz-1);
      uintfinalzr =  get_u(ppreexplicit,UU,ix,iy,iz+1);
      /*
      uintfinalpre = calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy,iz),  get_u(ppreexplicit,RHO,ix,iy,iz),ELECTRONS,ix,iy,iz);
      uintfinalxl =  calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix-1,iy,iz),get_u(ppreexplicit,RHO,ix-1,iy,iz),ELECTRONS,ix-1,iy,iz);
      uintfinalxr =  calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix+1,iy,iz),get_u(ppreexplicit,RHO,ix+1,iy,iz),ELECTRONS,ix+1,iy,iz);
      uintfinalyl =  calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy-1,iz),get_u(ppreexplicit,RHO,ix,iy-1,iz),ELECTRONS,ix,iy-1,iz);
      uintfinalyr =  calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy+1,iz),get_u(ppreexplicit,RHO,ix,iy+1,iz),ELECTRONS,ix,iy+1,iz);
      uintfinalzl =  calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy,iz-1),get_u(ppreexplicit,RHO,ix,iy,iz-1),ELECTRONS,ix,iy,iz-1);
      uintfinalzr =  calc_ufromSerho(get_u(ppreexplicit,ENTRE,ix,iy,iz+1),get_u(ppreexplicit,RHO,ix,iy,iz+1),ELECTRONS,ix,iy,iz+1);
      */
      u0 =my_min(calc_ufromSerho(s0*nefinal,rhoefinal,ELECTRONS,ix,iy,iz),LIMITFACTORINMIXING*uintfinalpre);
      uxl=my_min(calc_ufromSerho(sxl*nefinal,rhoefinal,ELECTRONS,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalxl));
      uxr=my_min(calc_ufromSerho(sxr*nefinal,rhoefinal,ELECTRONS,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalxr));
      uyl=my_min(calc_ufromSerho(syl*nefinal,rhoefinal,ELECTRONS,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalyl));
      uyr=my_min(calc_ufromSerho(syr*nefinal,rhoefinal,ELECTRONS,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalyr));
      uzl=my_min(calc_ufromSerho(szl*nefinal,rhoefinal,ELECTRONS,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalzl));
      uzr=my_min(calc_ufromSerho(szr*nefinal,rhoefinal,ELECTRONS,ix,iy,iz),LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalzr));
      #endif 

      ufinal=f0*u0+fxl*uxl+fxr*uxr+fyl*uyl+fyr*uyr+fzl*uzl+fzr*uzr;
      Sfinal=calc_Sefromrhou(rhoefinal,ufinal,ELECTRONS);
      
#endif //OLDENTROPYMIXING
      
      if(!isfinite(Sfinal))
	{
	  Sfinal=Sbake;
	}

      set_u(p,ENTRE,ix,iy,iz,Sfinal);

      //*****************************************
      //******** ions ***************************
      //*****************************************
      
#ifdef MIXENTROPIES_CONST_PRESSURE
      //initial state temperature and gamma
      ldouble Tiinit=calc_TfromSerho(get_u(ppreexplicit,ENTRI,ix,iy,iz), rhoinit, IONS, ix, iy, iz);
      sinit = get_u(ppreexplicit,ENTRI,ix,iy,iz) / niinit; //scale by rho below to get ents per particle in this case


      #if defined(CONSISTENTGAMMA) && !defined(FIXEDGAMMASPECIES)
      gammai = calc_gammaintfromtemp(Tiinit, ELECTRONS);
      #else
      gammai = GAMMAI;
      #endif
      gm1 = gammai-1.;
      ginv = 1./gammai;
    
      kappainit = (niinit*K_BOLTZ*Tiinit) / pow(rhoinit, gammai);
      s_scale= (MU_I * M_PROTON);
      
      //find the scale factor that relates kappa to exp(s)
      #if !defined(CONSISTENTGAMMA) && defined(NOLOGINS2)
      gfac = (gammai-1.);
      scalefac = kappainit / (sinit*gfac); 
      #else
      gfac = (gammai-1.)/K_BOLTZ;      
      scalefac = kappainit / (exp(sinit*gfac));
      #endif

#endif
      
      //ANDREW -- OK to scale by rho, since n_i = rho/mu_i/m_proton everywhere      
      s0=get_u(upreexplicit,ENTRI,ix,iy,iz)/get_u(upreexplicit,RHO,ix,iy,iz);
      sxl=get_ub(flbx,ENTRI,ix,iy,iz,0)/get_ub(flbx,RHO,ix,iy,iz,0);
      sxr=get_ub(flbx,ENTRI,ix+1,iy,iz,0)/get_ub(flbx,RHO,ix+1,iy,iz,0);
      syl=get_ub(flby,ENTRI,ix,iy,iz,1)/get_ub(flby,RHO,ix,iy,iz,1);
      syr=get_ub(flby,ENTRI,ix,iy+1,iz,1)/get_ub(flby,RHO,ix,iy+1,iz,1);
      szl=get_ub(flbz,ENTRI,ix,iy,iz,2)/get_ub(flbz,RHO,ix,iy,iz,2);
      szr=get_ub(flbz,ENTRI,ix,iy,iz+1,2)/get_ub(flbz,RHO,ix,iy,iz+1,2);

      //upwind scheme
#if defined(UPWINDENTROPYMIXING) || defined(PARTLYUPWINDENTROPYMIXING)
      rs_0 = get_u(ppreexplicit,ENTRI,ix,iy,iz)/get_u(ppreexplicit,RHO,ix,iy,iz);
      if(fxl<0) rsxl=rs_0;
      else rsxl=get_u(ppreexplicit,ENTRI,ix-1,iy,iz)/get_u(ppreexplicit,RHO,ix-1,iy,iz) ;
      if(fxr<0) rsxr=rs_0;
      else rsxr=get_u(ppreexplicit,ENTRI,ix+1,iy,iz)/get_u(ppreexplicit,RHO,ix+1,iy,iz) ;
      if(fyl<0) rsyl=rs_0;
      else rsyl=get_u(ppreexplicit,ENTRI,ix,iy-1,iz)/get_u(ppreexplicit,RHO,ix,iy-1,iz) ;
      if(fyr<0) rsyr=rs_0;
      else rsyr=get_u(ppreexplicit,ENTRI,ix,iy+1,iz)/get_u(ppreexplicit,RHO,ix,iy+1,iz) ;
      if(fzl<0) rszl=rs_0;
      else rszl=get_u(ppreexplicit,ENTRI,ix,iy,iz-1)/get_u(ppreexplicit,RHO,ix,iy,iz-1) ;
      if(fzr<0) rszr=rs_0;
      else rszr=get_u(ppreexplicit,ENTRI,ix,iy,iz+1)/get_u(ppreexplicit,RHO,ix,iy,iz+1) ;
     
      #ifdef UPWINDENTROPYMIXING //corrects to cell center everywhere
      sxl=rsxl;sxr=rsxr;syl=rsyl;syr=rsyr;szl=rszl;szr=rszr;
      #endif
     
      #ifdef PARTLYUPWINDENTROPYMIXING //corrects to cell center when outflowing
      if(fxl<0) sxl=rsxl;
      if(fxr<0) sxr=rsxr;
      if(fyl<0) syl=rsyl;
      if(fyr<0) syr=rsyr;
      if(fzl<0) szl=rszl;
      if(fzr<0) szr=rszr;
      #endif

#elif defined(PARTLYRECONSTRUCTEDUPWINDENTROPYMIXING)
      if(fxl<0) rsxl=get_ub(pbRx,ENTRI,ix,iy,iz,0)/get_ub(pbRx,RHO,ix,iy,iz,0);
      else rsxl=get_ub(pbLx,ENTRI,ix,iy,iz,0)/get_ub(pbLx,RHO,ix,iy,iz,0);
      if(fxr<0) rsxr=get_ub(pbLx,ENTRI,ix+1,iy,iz,0)/get_ub(pbLx,RHO,ix+1,iy,iz,0);
      else rsxr=get_ub(pbRx,ENTRI,ix+1,iy,iz,0)/get_ub(pbRx,RHO,ix+1,iy,iz,0);
      if(fyl<0) rsyl=get_ub(pbRx,ENTRI,ix,iy,iz,1)/get_ub(pbRx,RHO,ix,iy,iz,1);
      else rsyl=get_ub(pbLx,ENTRI,ix,iy,iz,1)/get_ub(pbLx,RHO,ix,iy,iz,1);
      if(fyr<0) rsyr=get_ub(pbLx,ENTRI,ix,iy+1,iz,1)/get_ub(pbLx,RHO,ix,iy+1,iz,1);
      else rsyr=get_ub(pbRx,ENTRI,ix,iy+1,iz,1)/get_ub(pbRx,RHO,ix,iy+1,iz,1);
      if(fzl<0) rszl=get_ub(pbRx,ENTRI,ix,iy,iz,2)/get_ub(pbRx,RHO,ix,iy,iz,2);
      else rszl=get_ub(pbLx,ENTRI,ix,iy,iz,2)/get_ub(pbLx,RHO,ix,iy,iz,2);
      if(fzr<0) rszr=get_ub(pbLx,ENTRI,ix,iy,iz+1,2)/get_ub(pbLx,RHO,ix,iy,iz+1,2);
      else rszr=get_ub(pbRx,ENTRI,ix,iy,iz+1,2)/get_ub(pbRx,RHO,ix,iy,iz+1,2);
      
      sxl=rsxl;sxr=rsxr;syl=rsyl;syr=rsyr;szl=rszl;szr=rszr;
#endif
      
      //floors
      if(!isfinite(sxl)) sxl=s0; 
      if(!isfinite(sxr)) sxr=s0;
      if(!isfinite(syl)) syl=s0;
      if(!isfinite(syr)) syr=s0;
      if(!isfinite(szl)) szl=s0;
      if(!isfinite(szr)) szr=s0;

      //ANDREW -- ok to scale by rho, since n_i = rho/mu_i/m_proton everywhere
      Sbaki = rhofinal*(f0*s0 + fxl*sxl + fxr*sxr + fyl*syl + fyr*syr + fzl*szl + fzr*szr);

#if defined(OLDENTROPYMIXING)
      Sfinal = Sbaki;

#elif defined(MIXENTROPIES_CONST_PRESSURE)      

      ldouble rhofinalgi = pow(rhofinal, gammai);

      //ANDREW this should be right....but the answer looks wrong!
      gfac = gfac*s_scale;
      
      #if !defined(CONSISTENTGAMMA) && defined(NOLOGINS2)      
      p0  = s0 *gfac*scalefac*rhofinalgi;
      pxl = sxl*gfac*scalefac*rhofinalgi;
      pxr = sxr*gfac*scalefac*rhofinalgi;
      pyl = syl*gfac*scalefac*rhofinalgi;
      pyr = syr*gfac*scalefac*rhofinalgi;
      pzl = szl*gfac*scalefac*rhofinalgi;
      pzr = szr*gfac*scalefac*rhofinalgi;
      #else
      p0  = exp(s0 *gfac)*scalefac*rhofinalgi;
      pxl = exp(sxl*gfac)*scalefac*rhofinalgi;
      pxr = exp(sxr*gfac)*scalefac*rhofinalgi;
      pyl = exp(syl*gfac)*scalefac*rhofinalgi;
      pyr = exp(syr*gfac)*scalefac*rhofinalgi;
      pzl = exp(szl*gfac)*scalefac*rhofinalgi;
      pzr = exp(szr*gfac)*scalefac*rhofinalgi;
      #endif
      
      #ifndef DONOTLIMITENTRINMIXING
      
      pintfinalpre = gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy,iz),  get_u(ppreexplicit,RHO,ix,iy,iz),  IONS,ix,iy,iz);
      pintfinalxl =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix-1,iy,iz),get_u(ppreexplicit,RHO,ix-1,iy,iz),IONS,ix-1,iy,iz);
      pintfinalxr =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix+1,iy,iz),get_u(ppreexplicit,RHO,ix+1,iy,iz),IONS,ix+1,iy,iz);
      pintfinalyl =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy-1,iz),get_u(ppreexplicit,RHO,ix,iy-1,iz),IONS,ix,iy-1,iz);
      pintfinalyr =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy+1,iz),get_u(ppreexplicit,RHO,ix,iy+1,iz),IONS,ix,iy+1,iz);
      pintfinalzl =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy,iz-1),get_u(ppreexplicit,RHO,ix,iy,iz-1),IONS,ix,iy,iz-1);
      pintfinalzr =  gm1*calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy,iz+1),get_u(ppreexplicit,RHO,ix,iy,iz+1),IONS,ix,iy,iz+1);
      
      /*
      pintfinalpre = gm1* get_u(ppreexplicit,UU,ix,iy,iz);
      pintfinalxl =  gm1* get_u(ppreexplicit,UU,ix-1,iy,iz);
      pintfinalxr =  gm1* get_u(ppreexplicit,UU,ix+1,iy,iz);
      pintfinalyl =  gm1* get_u(ppreexplicit,UU,ix,iy-1,iz);
      pintfinalyr =  gm1* get_u(ppreexplicit,UU,ix,iy+1,iz);
      pintfinalzl =  gm1* get_u(ppreexplicit,UU,ix,iy,iz-1);
      pintfinalzr =  gm1* get_u(ppreexplicit,UU,ix,iy,iz+1);
      */
      
      p0  = my_min(p0, LIMITFACTORINMIXING*pintfinalpre);
      pxl = my_min(pxl, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalxl));
      pxr = my_min(pxr, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalxr));
      pyl = my_min(pyl, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalyl));
      pyr = my_min(pyr, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalyr));
      pzl = my_min(pzl, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalzl));
      pzr = my_min(pzr, LIMITFACTORINMIXING*my_max(pintfinalpre, pintfinalzr));
      #endif

      entsum = f0*pow(p0,ginv)+fxl*pow(pxl,ginv)+fxr*pow(pxr,ginv)+fyl*pow(pyl,ginv)+fyr*pow(pyr,ginv)+fzl*pow(pzl,ginv)+fzr*pow(pzr,ginv);
      Tifinal = pow(entsum, gammai) / (nifinal * K_BOLTZ);
      Sfinal = calc_SefromrhoT(rhofinal,Tifinal,IONS);


#else
      //ANDREW mix entropies constant density
      //ANDREW not as consistent as mixing at constant pressure

      #ifdef DONOTLIMITENTRINMIXING 
      u0 =calc_ufromSerho(s0*rhofinal, rhofinal,IONS,ix,iy,iz);
      uxl=calc_ufromSerho(sxl*rhofinal,rhofinal,IONS,ix,iy,iz);
      uxr=calc_ufromSerho(sxr*rhofinal,rhofinal,IONS,ix,iy,iz);
      uyl=calc_ufromSerho(syl*rhofinal,rhofinal,IONS,ix,iy,iz);
      uyr=calc_ufromSerho(syr*rhofinal,rhofinal,IONS,ix,iy,iz);
      uzl=calc_ufromSerho(szl*rhofinal,rhofinal,IONS,ix,iy,iz);
      uzr=calc_ufromSerho(szr*rhofinal,rhofinal,IONS,ix,iy,iz);
      #else
      
      /*
      uintfinalpre = get_u(ppreexplicit,UU,ix,iy,iz);
      uintfinalxl =  get_u(ppreexplicit,UU,ix-1,iy,iz);
      uintfinalxr =  get_u(ppreexplicit,UU,ix+1,iy,iz);
      uintfinalyl =  get_u(ppreexplicit,UU,ix,iy-1,iz);
      uintfinalyr =  get_u(ppreexplicit,UU,ix,iy+1,iz);
      uintfinalzl =  get_u(ppreexplicit,UU,ix,iy,iz-1);
      uintfinalzr =  get_u(ppreexplicit,UU,ix,iy,iz+1);
      */
      
      uintfinalpre = calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy,iz),get_u(ppreexplicit,RHO,ix,iy,iz),IONS,ix,iy,iz);
      uintfinalxl =  calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix-1,iy,iz),get_u(ppreexplicit,RHO,ix-1,iy,iz),IONS,ix-1,iy,iz);
      uintfinalxr =  calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix+1,iy,iz),get_u(ppreexplicit,RHO,ix+1,iy,iz),IONS,ix+1,iy,iz);
      uintfinalyl =  calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy-1,iz),get_u(ppreexplicit,RHO,ix,iy-1,iz),IONS,ix,iy-1,iz);
      uintfinalyr =  calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy+1,iz),get_u(ppreexplicit,RHO,ix,iy+1,iz),IONS,ix,iy+1,iz);
      uintfinalzl =  calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy,iz-1),get_u(ppreexplicit,RHO,ix,iy,iz-1),IONS,ix,iy,iz-1);
      uintfinalzr =  calc_ufromSerho(get_u(ppreexplicit,ENTRI,ix,iy,iz+1),get_u(ppreexplicit,RHO,ix,iy,iz+1),IONS,ix,iy,iz+1);
            
      u0 =my_min(calc_ufromSerho(s0*rhofinal, rhofinal,IONS,ix,iy,iz), LIMITFACTORINMIXING*uintfinalpre);
      uxl=my_min(calc_ufromSerho(sxl*rhofinal,rhofinal,IONS,ix,iy,iz), LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalxl));
      uxr=my_min(calc_ufromSerho(sxr*rhofinal,rhofinal,IONS,ix,iy,iz), LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalxr));
      uyl=my_min(calc_ufromSerho(syl*rhofinal,rhofinal,IONS,ix,iy,iz), LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalyl));
      uyr=my_min(calc_ufromSerho(syr*rhofinal,rhofinal,IONS,ix,iy,iz), LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalyr));
      uzl=my_min(calc_ufromSerho(szl*rhofinal,rhofinal,IONS,ix,iy,iz), LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalzl));
      uzr=my_min(calc_ufromSerho(szr*rhofinal,rhofinal,IONS,ix,iy,iz), LIMITFACTORINMIXING*my_max(uintfinalpre,uintfinalzr));
      #endif
      
      ufinal=f0*u0+fxl*uxl+fxr*uxr+fyl*uyl+fyr*uyr+fzl*uzl+fzr*uzr;
      Sfinal=calc_Sefromrhou(rhofinal,ufinal,IONS);
#endif
      
      if(!isfinite(Sfinal))
	{
	  Sfinal=Sbaki;
	}

      //Sfinal=Sbaki;
      set_u(p,ENTRI,ix,iy,iz,Sfinal);
  
#endif //EVOLVEELECTRONS
#endif // DONOTMIXSPECIESENTROPY

      int sanitycheck=0;
      if(!isfinite(get_u(p,ENTR,ix,iy,iz))) sanitycheck=-1;
#ifdef EVOLVEELECTRONS
      if(!isfinite(get_u(p,ENTRE,ix,iy,iz))) sanitycheck=-1;
      if(!isfinite(get_u(p,ENTRI,ix,iy,iz))) sanitycheck=-1;
#endif

      if(sanitycheck < 0) //something failed, let's use Godunov evolution
	{
	  set_u(p,ENTR,ix,iy,iz,Sbak);
#ifdef EVOLVEELECTRONS
#ifndef DONOTMIXSPECIESENTROPY
	  set_u(p,ENTRE,ix,iy,iz,Sbake);
	  set_u(p,ENTRI,ix,iy,iz,Sbaki);
#endif
#endif	  
	}
	  
      //recalculating conserved quantities
      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);
      p2u_mhd(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);

    }
  
      return 0;
}

