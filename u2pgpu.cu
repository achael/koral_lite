/*! \file u2pgpu.cu
 \brief MHD Conserved to primitives conversion
 */

extern "C" {

#include "ko.h"

}

#include "kogpu.h"

// reminder of default settings
//#define U2PCONV 1.e-10
//#define U2P_EQS U2P_EQS_NOBLE
//#define U2P_SOLVER U2P_SOLVER_W


// todo deleted type
__global__ void calc_primitives_kernel(int Nloop_0, int setflags,
				       int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
				       ldouble *u_arr, ldouble *p_arr,
				       ldouble *x_arr, ldouble *g_arr, ldouble *G_arr)
{

  int verbose=0;
  
  // get index for this thread
  // Nloop_0 is number of cells to update;
  // usually Nloop_0=NX*NY*NZ, but sometimes there are weird bcs inside domain 
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_0) return;
    
  // get indices from 1D arrays
  int ix=loop_0_ix[ii];
  int iy=loop_0_iy[ii];
  int iz=loop_0_iz[ii]; 
  
  //skip if cell is passive
  if(!is_cell_active_device(ix,iy,iz))
    return;
  
  struct geometry geom;
  //fill_geometry(ix,iy,iz,&geom);
  fill_geometry_device(ix,iy,iz, x_arr,&geom,g_arr, G_arr);

   
  int u2pret,u2pretav;
  ldouble uu[NV],pp[NV];
  
  int corrected[3]={0,0,0}, fixups[2]={0,0};
  for(int iv=0;iv<NV;iv++)
  {
    uu[iv]=get_u(u_arr,iv,ix,iy,iz);
    pp[iv]=get_u(p_arr,iv,ix,iy,iz);
  }

  //TODO -- put in flags
  /*
  if(setflags)
  {
    set_cflag(ENTROPYFLAG,ix,iy,iz,0);
    set_cflag(ENTROPYFLAG2,ix,iy,iz,0);
  }
  */

  //TODO -- put in check for corrected_polaraxis
  //u to p inversion is done here
  //if(is_cell_corrected_polaraxis(ix,iy,iz))
  //{
  //  u2p_solver_Bonly(uu,pp,&geom); // invert only the magnetic field, the rest will be overwritten
  //}
  //else
  //{
  //u2p_device(uu,pp,&geom,corrected,fixups); // regular inversion
  //}

    //TODO
  /*
  //set flags for entropy solver
  if(corrected[0]==1 && setflags) //hd correction - entropy solver
  {
    set_cflag(ENTROPYFLAG,ix,iy,iz,1);
  }
  
  if(corrected[2]==1 && setflags) //borrowing energy from radiation didn't work
  {  
    set_cflag(ENTROPYFLAG2,ix,iy,iz,1);
  }
  
  //check hd floors
  int floorret=0;
  
  if(is_cell_active(ix,iy,iz) && !is_cell_corrected_polaraxis(ix,iy,iz))
  {
    floorret=check_floors_mhd(pp,VELPRIM,&geom);
  }
  
  if(floorret<0.)
  {
    corrected[0]=1;
  }

  /*
  //check rad floors
#ifdef RADIATION
  floorret=0;
  if(is_cell_active(ix,iy,iz) &&  !is_cell_corrected_polaraxis(ix,iy,iz))
  {
    floorret=check_floors_rad(pp,VELPRIMRAD,&geom);
  }
  
  if(floorret<0.)
  {
    corrected[1]=1;
  }
#endif
  */
    
  //set new primitives and conserved
  for(int iv=0;iv<NV;iv++)
  { 
    set_u(p_arr,iv,ix,iy,iz,pp[iv]);
  }

  //TODO
  /*
  //set flags for fixups of unsuccessful cells
  if(setflags)
  {
    if(fixups[0]>0)
    {
      set_cflag(HDFIXUPFLAG,ix,iy,iz,1);
      global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]++;
    }
    else
      set_cflag(HDFIXUPFLAG,ix,iy,iz,0);
    
    if(fixups[1]>0)
    {
      set_cflag(RADFIXUPFLAG,ix,iy,iz,-1);
      global_int_slot[GLOBALINTSLOT_NTOTALRADFIXUPS]++;
    }
    else
      set_cflag(RADFIXUPFLAG,ix,iy,iz,0); 
  }
    */

} 




//**********************************************************************
//high-level u2p solver
// type: not used
//**********************************************************************
/*
__device__ __host__ int u2p_device(ldouble *uu0, ldouble *pp, void *ggg, int corrected[3], int fixups[2])
{
  struct geometry *geom
  = (struct geometry *) ggg;
  
  ldouble uu[NV];
  int iv;
  PLOOP(iv) uu[iv]=uu0[iv];
  
  ldouble (*gg)[5], (*GG)[5], gdet, gdetu, gdetu_inv;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif
  gdetu_inv = 1. / gdetu;
  
  int verbose=0;
  int hdcorr=0;
  int radcor=0;
  corrected[0]=corrected[1]=0;
  fixups[0]=fixups[1]=0;
  
  int u2pret,u2pentrret,ret;
  ldouble ppbak[NV];
  for(u2pret=0;u2pret<NV;u2pret++)
    ppbak[u2pret]=pp[u2pret];

  //************************************
  //magneto-hydro part
  //************************************
  ldouble u0=pp[1];
  
  //************************************
  //hot hydro - conserving energy
  ret=0;
  u2pret=-1;
  
  //test
  ldouble ppold[NV];
  ldouble ppentr[NV];
  PLOOP(iv)
  {
    ppold[iv]=pp[iv];
    ppentr[iv]=-1.; // negative value indicates that entropy inversion not yet calculated
  }
  
  //negative uu[0] = rho u^t
  if(uu[0] * gdetu_inv < 0.)
  {
    int gix,giy,giz;
    mpi_local2globalidx(geom->ix,geom->iy,geom->iz,&gix,&giy,&giz);
    if(verbose) printf("%4d > %4d %4d %4d > NEGUU  > neg uu[0] - requesting fixup\n",PROCID,gix,giy,giz);
    pp[0]=RHOFLOOR; //used when not fixing up
    uu[0]=RHOFLOOR*gdetu;
    ret=-2;    //to request fixup
               //ANDREW -- but ret=-1 if energy inversion failes but entropy inversion does not!
               //ANDREW -- do we always want a fixup if we have negative uu[0] ? 
    u2pret=-1; // indicates that inversion is needed
    
#ifndef SWAPPAPC
    global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]++;  //but count as fixup
#endif
  }
    
  if(u2pret!=0)  // u2pret=-1 at this stage, so this is always satisfied
  {
#ifdef ENFORCEENTROPY  
    u2pret=-1;  //skip hot energy-conserving inversion and go to entropy inversion
#else
    u2pret = u2p_solver(uu,pp,ggg,U2P_HOT,0);  // invert using the hot energy equation    
#endif //ENFORCEENTROPY
  }
  
  if(ALLOWENTROPYU2P)  // Inversion with entropy equation -- on by default (see choices.h)
  {
    if(u2pret<0)  // true if energy equation failed, or if energy equation was not required (because ENFORCEENTROPY is defined)
    {
      ret=-1;
      
      if(verbose>2 )
      {
        printf("u2p_entr     >>> %d %d <<< %d >>> %e > %e\n",geom->ix + TOI, geom->iy + TOJ,u2pret,u0,pp[1]);
      }
      
      //************************************
      //entropy solver - conserving entropy
      if(ppentr[RHO]<0.) //if not yet calculated
      {
        u2pret=u2p_solver(uu,pp,ggg,U2P_ENTROPY,0);  // invert using entropy equation
      }
      
      if(u2pret<0)
      {
        ret=-2;
        
        if(verbose>1)
        {
          printf("u2p_entr err No. %4d > %e %e %e > %e %e > %4d %4d %4d\n",u2pret,uu[0],uu[1],uu[5],pp[0],pp[1],geom->ix,geom->iy,geom->iz);
        }
	
      } // if(u2pret<0) // second time -- entropy eqn
    } // if(u2pret<0) // first time -- energy eqn
  }  // if(ALLOWENTROPYU2P)
  
  if(u2pret<0)  // entropy equation also failed
  {
 
    //leaving primitives unchanged - should not happen
    if(verbose>1 || 1)
    {
      printf("%4d > %4d %4d %4d > MHDU2PFAIL > u2p prim. unchanged > %d \n",PROCID,geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK,u2pret);
    }
    ret=-3;
    for(u2pret=0;u2pret<NV;u2pret++)
      pp[u2pret]=ppbak[u2pret];
  }
  
  if(ret<0) //to update conserved
    hdcorr=1;  
  if(ret<-1) //request fixup when entropy failed
    fixups[0]=1;
  else
    fixups[0]=0;
  
  //************************************
  //radiation part
  //************************************
  
  corrected[2]=0;
  
#ifdef RADIATION  
#ifdef BALANCEENTROPYWITHRADIATION
  
  //trying to balance gain of energy because of entropy inversion
  //by borrowing from the radiation field
  if(ret==-1) //entropy u2p was used in MHD part
  {
    ldouble uunew[NV],ppnew[NV];
    PLOOP(iv) { uunew[iv]=uu[iv]; ppnew[iv]=pp[iv]; }
    p2u_mhd(pp,uunew,geom);
    ldouble dugas = uunew[UU] - uu[UU];  //this much energy was introduced
    if(fabs(dugas)<0.1*fabs(uunew[EE0])) //correction relatively small - is this general enough?
    {
      uunew[EE0]-=dugas; //balancing with radiation
      u2p_rad(uunew,ppnew,geom,&radcor);
    }
    else
      radcor=1;
    
    if(radcor==0) //there was enough energy to borrow from and uunew inverts with hot
    {
      PLOOP(iv)
      uu[iv]=uunew[iv];
      //printf("entropy correction worked at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
    }
    else
    {
      corrected[2]=1; //entropy correction didn't work
      //printf("entropy correction didn't work at %d %d\n",geom->ix+TOI,geom->iy+TOJ);
    }
  }
#endif //BALANCEENTROPYWITHRADIATION

  //Do the radiative inversion from u2p_rad.c
  u2p_rad(uu,pp,geom,&radcor);

#endif // RADIATION

  //************************************  
  //output
  //************************************
  
  //rad fixups only for critical failure in implicit
  if(radcor>0)     
    fixups[1]=1;
  else
    fixups[1]=0;
    
  if(hdcorr>0) corrected[0]=1;
  if(radcor>0) corrected[1]=1;
  
  return ret;
} 
*/


int calc_u2p_gpu(int type, int setflags)
{

  ldouble *d_u_arr, *d_p_arr;

  cudaError_t err = cudaSuccess;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // NOTE: size of xb is copied from initial malloc in misc.c
  // these need to be long long if the grid is on one tile and large (~256^3)
  long long Nprim  = (SX)*(SY)*(SZ)*NV;

  // allocate and sync prims and cons to device 
  err = cudaMalloc(&d_p_arr,    sizeof(ldouble)*Nprim);  
  err = cudaMalloc(&d_u_arr,    sizeof(ldouble)*Nprim);

  err = cudaMemcpy(d_p_arr, p, sizeof(ldouble)*Nprim, cudaMemcpyHostToDevice);  // is this used to seed?
  err = cudaMemcpy(d_u_arr, u, sizeof(ldouble)*Nprim, cudaMemcpyHostToDevice);
  

  // launch calc_primitives_kernel

  int threadblocks = (Nloop_0 / TB_SIZE) + ((Nloop_0 % TB_SIZE)? 1:0);

  cudaEventRecord(start);
  calc_primitives_kernel<<<threadblocks, TB_SIZE>>>(Nloop_0, setflags, 
                                                    d_loop0_ix, d_loop0_iy, d_loop0_iz,
                                                    d_u_arr, d_p_arr,
                                                    d_x, d_gcov, d_gcon);

  cudaEventRecord(stop);
  err = cudaPeekAtLastError();
  cudaDeviceSynchronize(); //TODO: do we need this, does cudaMemcpy synchronize?

  cudaEventSynchronize(stop);
  float tms = 0.;
  cudaEventElapsedTime(&tms, start,stop);
  printf("gpu update time: %0.2f \n",tms);

  // ======= TODO
  // Free Device Memory
  cudaFree(d_u_arr);
  cudaFree(d_p_arr);

  return 0;
}



