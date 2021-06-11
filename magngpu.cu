/*! \file magn.c
 \brief Magnetic field related routines
*/

extern "C" {

#include "ko.h"

}

#include "kogpu.h"

static __device__ int fl_x(int i);
static __device__ int fl_y(int i);
static __device__ int fl_z(int i);
static __device__ int flx_x(int i);
static __device__ int flx_y(int i);
static __device__ int flx_z(int i);
static __device__ int fly_x(int i);
static __device__ int fly_y(int i);
static __device__ int fly_z(int i);
static __device__ int flz_x(int i);
static __device__ int flz_y(int i);
static __device__ int flz_z(int i);

//***********************************************************************
/* calculate both magnetic field four-vectors and bsq knowing gas four-velocity ucov */
//***********************************************************************

__device__ __host__ void calc_bcon_bcov_bsq_from_4vel_device(ldouble *pr,
							     ldouble *ucon, ldouble *ucov, void* ggg,
		                        		     ldouble *bcon, ldouble *bcov, ldouble *bsq)
{

  struct geometry *geom
    = (struct geometry *) ggg;

  // First calculate bcon0
  bcon[0] = pr[B1]*ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3] ;
  
  // Then spatial components of bcon
#ifdef NONRELMHD
  for(int j = 1; j < 4; j++)
  {
    bcon[j] = pr[B1-1+j]; //b^i=B^i
  }
#else  // relativistic case
  
  ldouble u0inv = 1. / ucon[0];
  for(int j=1;j<4;j++)
  {
    bcon[j] = (pr[B1-1+j] + bcon[0] * ucon[j]) * u0inv ;
  }
#endif //NONRELMHD
  
  // Convert to bcov and calculate bsq
  indices_21_device(bcon, bcov, geom->gg);
  *bsq = dotB(bcon, bcov); //NOTE: preprocessor macro, ok

  return ;
}

/***********************************************************************************************/
/* wrappers for missing cells / dimensions */
/***********************************************************************************************/

static __device__ int fl_x(int i)
{
  if(NX==1) return 0;
  return i;
}

static __device__ int fl_y(int i)
{
  if(NY==1) return 0;
  return i;
}

static __device__ int fl_z(int i)
{
  if(NZ==1) return 0;
  return i;
}

static __device__ int flx_x(int i)
{
  return i;
}

static __device__ int flx_y(int i)
{
  return fl_y(i);
}

static __device__ int flx_z(int i)
{
  return fl_z(i);
}

static __device__ int fly_x(int i)
{
  return fl_x(i);
}

static __device__ int fly_y(int i)
{
  return i;
}

static __device__ int fly_z(int i)
{
  return fl_z(i);
}

static __device__ int flz_x(int i)
{
  return fl_x(i);
}

static __device__ int flz_y(int i)
{
  return fl_y(i);
}

static __device__ int flz_z(int i)
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


__global__ void flux_ct_setemf_kernel(int Nloop_4,
			              int* loop_4_ix, int* loop_4_iy, int* loop_4_iz,
			              ldouble* emf_arr,
			              ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr)
{
#ifdef MAGNFIELD

  //TODO err
  /*
  //requires GDETIN = 1
  if(GDETIN==0)
    {
      my_err("MAGNFIELD requires GDETIN==1\n");
      exit(0);
    }
  */

  // get index for this thread
  // Nloop_4 is number of cell corners to update;
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_4) return;

  // cell indices
  int ix=loop_4_ix[ii];
  int iy=loop_4_iy[ii];
  int iz=loop_4_iz[ii]; 
  
  //TOTH algorithm from HARM's fluxct.c
  ldouble coefemf[4];
  if((NY>1)&&(NZ>1)) coefemf[1]=0.25;
  else coefemf[1]=0.5; 
  if((NX>1)&&(NZ>1)) coefemf[2]=0.25;
  else coefemf[2]=0.5; 
  if((NX>1)&&(NY>1)) coefemf[3]=0.25;
  else coefemf[3]=0.5; 

  ////////////////////
  // EMF1
  ////////////////////
#if((NY>1)||(NZ>1))
  set_emf(emf_arr,1,ix,iy,iz,
	  coefemf[1] * (
                            #if(NY>1)
			    + get_ub(flby_arr,B3,fly_x(ix),fly_y(iy),fly_z(iz),1) 
			    + get_ub(flby_arr,B3,fly_x(ix),fly_y(iy),fly_z(iz-1),1)
                            #endif
                            #if(NZ>1)
			    - get_ub(flbz_arr,B2,flz_x(ix),flz_y(iy),flz_z(iz),2) 
			    - get_ub(flbz_arr,B2,flz_x(ix),flz_y(iy-1),flz_z(iz),2)
                            #endif
			));
#else  
  set_emf(emf_arr,1,ix,iy,iz,0.); // not really 0, but differences in emf will be 0
#endif 
      
  ////////////////////
  // EMF2
  ////////////////////
#if((NX>1)||(NZ>1))
  set_emf(emf_arr,2,ix,iy,iz,
	  coefemf[2] * (
                            #if(NZ>1)
			    + get_ub(flbz_arr,B1,flz_x(ix),flz_y(iy),flz_z(iz),2) 
			    + get_ub(flbz_arr,B1,flz_x(ix-1),flz_y(iy),flz_z(iz),2)
                            #endif
                            #if(NX>1)
			    - get_ub(flbx_arr,B3,flx_x(ix),flx_y(iy),flx_z(iz),0) 
			    - get_ub(flbx_arr,B3,flx_x(ix),flx_y(iy),flx_z(iz-1),0)
                            #endif
			));
#else  
  set_emf(emf_arr,2,ix,iy,iz,0.);
#endif 

  ////////////////////
  // EMF3
  ////////////////////
#if((NX>1)||(NY>1))
  set_emf(emf_arr,3,ix,iy,iz,
	  coefemf[3] * (
                            #if(NX>1)
			    + get_ub(flbx_arr,B2,flx_x(ix),flx_y(iy),flx_z(iz),0) 
			    + get_ub(flbx_arr,B2,flx_x(ix),flx_y(iy-1),flx_z(iz),0)
                            #endif
                            #if(NY>1)
			    - get_ub(flby_arr,B1,fly_x(ix),fly_y(iy),fly_z(iz),1) 
			    - get_ub(flby_arr,B1,fly_x(ix-1),fly_y(iy),fly_z(iz),1)
                            #endif
			));
#else  
  set_emf(emf_arr,3,ix,iy,iz,0.);
#endif

  ///////////////////
  // Correct at corners
  // this is adjust_fluxcttoth_emfs() in magn.c
  //////////////////
#ifdef CORRECTPOLARAXIS
#ifdef MPI
  //if(TJ==0) //TODO MPI
#endif
  {
    // upper axis
    if(iy==0 && ix>=0 && ix<=NX && iz>=0 && iz<=NZ)
    {
      set_emf(emf_arr,1,ix,iy,iz,0.);
      set_emf(emf_arr,3,ix,iy,iz,0.);
    }
  }
#ifdef MPI
  // if(TJ==NTY-1) // TODO MPI
#endif
  {
    // lower axis
    if(iy==NY && ix>=0 && ix<=NX && iz>=0 && iz<=NZ)
    {
      set_emf(emf_arr,1,ix,iy,iz,0.);
      set_emf(emf_arr,3,ix,iy,iz,0.);
    }
  }
#
#endif //CORRECTPOLARAXIS
#endif //MAGNFIELD
}

__global__ void flux_ct_getemf_kernel(int Nloop_4,
			              int* loop_4_ix, int* loop_4_iy, int* loop_4_iz,
			              ldouble* emf_arr,
			              ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr)
{

#ifdef MAGNFIELD
  // get index for this thread
  // Nloop_4 is number of cell corners to update;
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_4) return;

  // cell indices
  int ix=loop_4_ix[ii];
  int iy=loop_4_iy[ii];
  int iz=loop_4_iz[ii]; 

  /////////////////////////////////////
  // F1
  ////////////////////////////////////
#if(NX>1)
   if(iy<NY && iz<NZ) //no need to fill x-face fluxes for iy=NY etc., 
   {
     set_ubx(flbx_arr,B1,ix,iy,iz,0.);
     set_ubx(flbx_arr,B2,ix,iy,iz,0.5  * (get_emf(emf_arr,3,ix,iy,iz) + get_emf(emf_arr,3,ix,iy+1,iz)));
     set_ubx(flbx_arr,B3,ix,iy,iz,-0.5 * (get_emf(emf_arr,2,ix,iy,iz) + get_emf(emf_arr,2,ix,iy,iz+1)));
   }
#endif

  /////////////////////////////////////
  // F2
  ////////////////////////////////////
#if(NY>1)
  if(ix<NX && iz<NZ)	
  {
    set_uby(flby_arr,B1,ix,iy,iz,-0.5 * (get_emf(emf_arr,3,ix,iy,iz) + get_emf(emf_arr,3,ix+1,iy,iz)));
    set_uby(flby_arr,B2,ix,iy,iz,0.);
    set_uby(flby_arr,B3,ix,iy,iz,0.5  * (get_emf(emf_arr,1,ix,iy,iz) + get_emf(emf_arr,1,ix,iy,iz+1)));
  }
#endif
      
			    
  /////////////////////////////////////
  // F3
  ////////////////////////////////////
#if(NZ>1)
  if(ix<NX && iy<NY)	
  {
    set_ubz(flbz_arr,B1,ix,iy,iz,0.5  * (get_emf(emf_arr,2,ix,iy,iz) + get_emf(emf_arr,2,ix+1,iy,iz)));
    set_ubz(flbz_arr,B2,ix,iy,iz,-0.5 * (get_emf(emf_arr,1,ix,iy,iz) + get_emf(emf_arr,1,ix,iy+1,iz)));
    set_ubz(flbz_arr,B3,ix,iy,iz,0.);
  }
#endif
	
#endif //MAGNFIELD

}


ldouble flux_ct_gpu()
{
  cudaError_t err = cudaSuccess;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  int threadblocks = (Nloop_4 / TB_SIZE) + ((Nloop_4 % TB_SIZE)? 1:0);f
  //printf("\nTest %d\n", threadblocks); fflush(stdout);

  cudaEventRecord(start);

  // set emf kernel
  flux_ct_setemf_kernel<<threadblocks, TB_SIZE>>(Nloop_4,
			                         d_loop4_ix, d_loop4_iy, d_loop4_iz,
			                         d_emf_arr,
			                         d_flbx_arr, d_flby_arr, d_flbz_arr);

  // synchronize
  err = cudaPeekAtLastError();
  // printf("ERROR-Kernel (error code %s)!\n", cudaGetErrorString(err));
  cudaDeviceSynchronize();

  // get emf kernel
  flux_ct_getemf_kernel<<threadblocks, TB_SIZE>>(Nloop_4,
			                         d_loop4_ix, d_loop4_iy, d_loop4_iz,
			                         d_emf_arr,
			                         d_flbx_arr, d_flby_arr, d_flbz_arr);

      
  // synchronize
  err = cudaPeekAtLastError();
  // printf("ERROR-Kernel (error code %s)!\n", cudaGetErrorString(err));
  cudaDeviceSynchronize();
  
  cudaEventRecord(stop);
  
  // timing information // TODO include Synchronize calls
  cudaEventSynchronize(stop);
  float tms = 0.;
  cudaEventElapsedTime(&tms, start,stop);
  printf("gpu flux_ct time: %0.2f \n",tms);
 
#ifdef CPUKO 
  ldouble* flbx_tmp;
  long long NfluxX = (SX+1)*(SY)*(SZ)*NV;
  if((flbx_tmp=(ldouble*)malloc(NfluxX*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  err = cudaMemcpy(flbx_tmp, d_flbx_arr, NfluxX*sizeof(ldouble), cudaMemcpyDeviceToHost);
  if(err != cudaSuccess) printf("failed cudaMemcpy of d_flbx_arr to flbx_tmp\n");
  printf("gpu flux_ct flbx[NV]: ");
  for(int iv=0;iv<NV;iv++)
    printf("%e ", get_ub(flbx_tmp, iv, ixTEST, iyTEST, izTEST,0));
  printf("\n");
  free(flbx_tmp);
#endif

  return (ldouble)tms;
}

