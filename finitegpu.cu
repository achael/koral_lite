extern "C" {

#include "ko.h"

}

#include "kogpu.h"

// persistent arrays, extern'd in kogpu.h
ldouble *d_p_arr, *d_u_arr;
ldouble *d_flbx_arr, *d_flby_arr, *d_flbz_arr;
ldouble *d_pbLx_arr, *d_pbLy_arr, *d_pbLz_arr;
ldouble *d_pbRx_arr, *d_pbRy_arr, *d_pbRz_arr;
ldouble *d_flLx_arr, *d_flLy_arr, *d_flLz_arr;
ldouble *d_flRx_arr, *d_flRy_arr, *d_flRz_arr;
ldouble *d_ahdxl_arr, *d_ahdyl_arr, *d_ahdzl_arr;
ldouble *d_ahdxr_arr, *d_ahdyr_arr, *d_ahdzr_arr;
ldouble *d_ahdx_arr, *d_ahdy_arr, *d_ahdz_arr;
ldouble *d_aradxl_arr, *d_aradyl_arr, *d_aradzl_arr;
ldouble *d_aradxr_arr, *d_aradyr_arr, *d_aradzr_arr;
ldouble *d_aradx_arr, *d_arady_arr, *d_aradz_arr;
ldouble *d_emf_arr;
ldouble *d_cell_tstepden_arr, *d_cell_dt_arr;
int *d_cellflag_arr, *d_int_slot_arr;

//taken from https://stackoverflow.com/questions/55140908/can-anybody-help-me-with-atomicmin-function-syntax-for-cuda
// TODO check if it works
// TODO is this the best solution for tstepdenmin/max? 
__device__ double atomicMin_double(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}
__device__ double atomicMax_double(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(fmax(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}


int prealloc_arrays_gpu()
{
  cudaError_t err = cudaSuccess;

  long long Ngrid  = (SX)*(SY)*(SZ);
  long long Nprim  = (SX)*(SY)*(SZ)*NV;
  long long NfluxX = (SX+1)*(SY)*(SZ)*NV;
  long long NfluxY = (SX)*(SY+1)*(SZ)*NV;
  long long NfluxZ = (SX)*(SY)*(SZ+1)*NV;
  long long Nemf = (NX+1)*(NY+1)*(NZ+1)*3;
  long long Ncellflag = (SX)*(SY)*(SZ)*NFLAGS;
  
  err = cudaMalloc(&d_p_arr,    sizeof(ldouble)*Nprim);
  err = cudaMalloc(&d_u_arr,    sizeof(ldouble)*Nprim);

  err = cudaMalloc(&d_pbLx_arr, sizeof(ldouble)*NfluxX);
  err = cudaMalloc(&d_pbLy_arr, sizeof(ldouble)*NfluxY);
  err = cudaMalloc(&d_pbLz_arr, sizeof(ldouble)*NfluxZ);

  err = cudaMalloc(&d_pbRx_arr, sizeof(ldouble)*NfluxX);
  err = cudaMalloc(&d_pbRy_arr, sizeof(ldouble)*NfluxY);
  err = cudaMalloc(&d_pbRz_arr, sizeof(ldouble)*NfluxZ);

  err = cudaMalloc(&d_flLx_arr, sizeof(ldouble)*NfluxX);
  err = cudaMalloc(&d_flLy_arr, sizeof(ldouble)*NfluxY);
  err = cudaMalloc(&d_flLz_arr, sizeof(ldouble)*NfluxZ);

  err = cudaMalloc(&d_flRx_arr, sizeof(ldouble)*NfluxX);
  err = cudaMalloc(&d_flRy_arr, sizeof(ldouble)*NfluxY);
  err = cudaMalloc(&d_flRz_arr, sizeof(ldouble)*NfluxZ);
  
  err = cudaMalloc(&d_flbx_arr, sizeof(ldouble)*NfluxX);
  err = cudaMalloc(&d_flby_arr, sizeof(ldouble)*NfluxY);
  err = cudaMalloc(&d_flbz_arr, sizeof(ldouble)*NfluxZ);

  err = cudaMalloc(&d_ahdxl_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_ahdyl_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_ahdzl_arr, sizeof(ldouble)*Ngrid);

  err = cudaMalloc(&d_ahdxr_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_ahdyr_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_ahdzr_arr, sizeof(ldouble)*Ngrid);

  err = cudaMalloc(&d_ahdx_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_ahdy_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_ahdz_arr, sizeof(ldouble)*Ngrid);

  err = cudaMalloc(&d_aradxl_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_aradyl_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_aradzl_arr, sizeof(ldouble)*Ngrid);

  err = cudaMalloc(&d_aradxr_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_aradyr_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_aradzr_arr, sizeof(ldouble)*Ngrid);

  err = cudaMalloc(&d_aradx_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_arady_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_aradz_arr, sizeof(ldouble)*Ngrid);
  
  err = cudaMalloc(&d_emf_arr,  sizeof(ldouble)*Nemf);

  err = cudaMalloc(&d_cell_tstepden_arr, sizeof(ldouble)*Ngrid);
  err = cudaMalloc(&d_cell_dt_arr, sizeof(ldouble)*Ngrid);
    
  err = cudaMalloc(&d_cellflag_arr, sizeof(int)*Ncellflag);
  err = cudaMalloc(&d_int_slot_arr, sizeof(int)*NGLOBALINTSLOT);

  // TODO: add error checks
  return 1;
}

int free_arrays_gpu()
{
  cudaFree(d_p_arr);
  cudaFree(d_u_arr);

  cudaFree(d_pbLx_arr);
  cudaFree(d_pbLy_arr);
  cudaFree(d_pbLz_arr);
  
  cudaFree(d_pbRx_arr);
  cudaFree(d_pbRy_arr);
  cudaFree(d_pbRz_arr);

  cudaFree(d_flLx_arr);
  cudaFree(d_flLy_arr);
  cudaFree(d_flLz_arr);
  
  cudaFree(d_flRx_arr);
  cudaFree(d_flRy_arr);
  cudaFree(d_flRz_arr);
  
  cudaFree(d_flbx_arr);
  cudaFree(d_flby_arr);
  cudaFree(d_flbz_arr);

  cudaFree(d_ahdxl_arr);
  cudaFree(d_ahdyl_arr);
  cudaFree(d_ahdzl_arr);

  cudaFree(d_ahdxr_arr);
  cudaFree(d_ahdyr_arr);
  cudaFree(d_ahdzr_arr);

  cudaFree(d_ahdx_arr);
  cudaFree(d_ahdy_arr);
  cudaFree(d_ahdz_arr);

  cudaFree(d_aradxl_arr);
  cudaFree(d_aradyl_arr);
  cudaFree(d_aradzl_arr);

  cudaFree(d_aradxr_arr);
  cudaFree(d_aradyr_arr);
  cudaFree(d_aradzr_arr);

  cudaFree(d_aradx_arr);
  cudaFree(d_arady_arr);
  cudaFree(d_aradz_arr);
  
  cudaFree(d_emf_arr);

  cudaFree(d_cell_tstepden_arr);
  cudaFree(d_cell_dt_arr);
  
  cudaFree(d_cellflag_arr);
  cudaFree(d_int_slot_arr);
  
  return 1;
}

int push_pu_gpu()
{
  // TODO: probably don't want to do it this way...

  cudaError_t err = cudaSuccess;
  
  if(doTEST==1) printf("H u: %e \n", get_u(u,ivTEST,ixTEST,iyTEST,izTEST));

  // copy prims, cons from host to device
  long long Ngrid  = (SX)*(SY)*(SZ);
  long long Nprim  = (SX)*(SY)*(SZ)*NV;
  err = cudaMemcpy(d_u_arr, u, sizeof(ldouble)*Nprim, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_p_arr, p, sizeof(ldouble)*Nprim, cudaMemcpyHostToDevice);

  // copy tstep arrays
  // TODO: do these need to be initialized or will they be entirely internal to GPU eventually?
  err = cudaMemcpy(d_cell_tstepden_arr, cell_tstepden, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_cell_dt_arr, cell_dt, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  
  // copy fluxes and wavespeeds from host to device
  // TODO: in the future, this will be entirely internal to the GPU
  
  long long NfluxX = (SX+1)*(SY)*(SZ)*NV;
  long long NfluxY = (SX)*(SY+1)*(SZ)*NV;
  long long NfluxZ = (SX)*(SY)*(SZ+1)*NV;

  err =  cudaMemcpy(d_pbLx_arr, pbLx, sizeof(ldouble)*NfluxX, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_pbLy_arr, pbLy, sizeof(ldouble)*NfluxY, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_pbLz_arr, pbLz, sizeof(ldouble)*NfluxZ, cudaMemcpyHostToDevice);

  err =  cudaMemcpy(d_pbRx_arr, pbRx, sizeof(ldouble)*NfluxX, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_pbRy_arr, pbRy, sizeof(ldouble)*NfluxY, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_pbRz_arr, pbRz, sizeof(ldouble)*NfluxZ, cudaMemcpyHostToDevice);

  err =  cudaMemcpy(d_flLx_arr, flLx, sizeof(ldouble)*NfluxX, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_flLy_arr, flLy, sizeof(ldouble)*NfluxY, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_flLz_arr, flLz, sizeof(ldouble)*NfluxZ, cudaMemcpyHostToDevice);

  err =  cudaMemcpy(d_flRx_arr, flRx, sizeof(ldouble)*NfluxX, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_flRy_arr, flRy, sizeof(ldouble)*NfluxY, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_flRz_arr, flRz, sizeof(ldouble)*NfluxZ, cudaMemcpyHostToDevice);
  
  err =  cudaMemcpy(d_flbx_arr, flbx, sizeof(ldouble)*NfluxX, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_flby_arr, flby, sizeof(ldouble)*NfluxY, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_flbz_arr, flbz, sizeof(ldouble)*NfluxZ, cudaMemcpyHostToDevice);

  err = cudaMemcpy(d_ahdxl_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdyl_arr, ahdyl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdzl_arr, ahdzl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_ahdxr_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdyr_arr, ahdyr, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdzr_arr, ahdzr, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_ahdx_arr, ahdx, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdy_arr, ahdy, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdz_arr, ahdz, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_aradxl_arr, aradxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradyl_arr, aradyl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradzl_arr, aradzl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_aradxr_arr, aradxr, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradyr_arr, aradyr, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradzr_arr, aradzr, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_aradx_arr, aradx, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_arady_arr, arady, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradz_arr, aradz, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  
  
  // TODO: add error checks
  return 1;
}

int pull_pu_gpu()
{
  // TODO: probably only want p, maybe rename
  cudaError_t err = cudaSuccess;

  // copy prims, cons back from device to host 
  long long Nprim  = (SX)*(SY)*(SZ)*NV;
  err = cudaMemcpy(p, d_p_arr, sizeof(ldouble)*Nprim, cudaMemcpyDeviceToHost); 
  err = cudaMemcpy(u, d_u_arr, sizeof(ldouble)*Nprim, cudaMemcpyDeviceToHost); 
  
  return 1;
}

// TODO this is here
void print_double_array_at(FILE *fp, double *array, int ix, int iy, int iz)
{
  int iv = 0;
  fprintf(fp, "[ [%d, %d, %d], [%g", ix, iy, iz, get_u(array, iv, ix, iy, iz));

  for (iv=1; iv<NV; iv++) {
    fprintf(fp, ", %g", get_u(array, iv, ix, iy, iz)); 
  } 

  fprintf(fp, "] ]");
}

int output_state_debug(const char *fname, const char *header, const char *ctimes, const char *gtimes)
{ 
  // writes a diagnostic output file in json format

  FILE *fp = fopen(fname, "w");
  fprintf(fp, "{\n");

  // write header and times
  fprintf(fp, "%s\n", header);
  fprintf(fp, "\"cpu_timing\": %s,\n", ctimes);
  fprintf(fp, "\"gpu_timing\": %s,\n", gtimes);  

  // TODO loop over zones
  fprintf(fp, "\"cpu_prims\":[\n");
  print_double_array_at(fp, p, ixTEST, iyTEST, izTEST);
  fprintf(fp, "\n],\n\"cpu_cons\":[\n");
  print_double_array_at(fp, u, ixTEST, iyTEST, izTEST);
  fprintf(fp, "\n],\n");

  // TODO, make this more modular
  {
  long long Nprim = (SX)*(SY)*(SZ)*NV;
  ldouble* p_tmp, *u_tmp;

  if((p_tmp=(ldouble*)malloc(sizeof(ldouble)*Nprim))==NULL) my_err("malloc err.\n");
  if((u_tmp=(ldouble*)malloc(sizeof(ldouble)*Nprim))==NULL) my_err("malloc err.\n");

  cudaError_t err = cudaSuccess;

  err = cudaMemcpy(p_tmp, d_p_arr, sizeof(ldouble)*Nprim, cudaMemcpyDeviceToHost);
  if(err != cudaSuccess) printf("failed cudaMemcpy of d_p_arr to p_tmp\n");

  err = cudaMemcpy(u_tmp, d_u_arr, sizeof(ldouble)*Nprim, cudaMemcpyDeviceToHost);
  if(err != cudaSuccess) printf("failed cudaMemcpy of d_u_arr to u_tmp\n");

  fprintf(fp, "\"gpu_prims\":[\n");
  print_double_array_at(fp, p_tmp, ixTEST, iyTEST, izTEST);
  fprintf(fp, "\n],\n\"gpu_cons\":[\n");
  print_double_array_at(fp, u_tmp, ixTEST, iyTEST, izTEST);
  fprintf(fp, "\n]\n");

  free(u_tmp);
  free(p_tmp);
  }

  fprintf(fp, "}");
  fclose(fp);

  // TODO error checking...
  return 1;
}

__device__ __host__ int is_cell_active_device(int ix, int iy, int iz)
{
  //NOTE: by default ALWAYS active -- this may change
  return 1;
}


__device__ __host__ int is_cell_corrected_polaraxis_device(int ix, int iy, int iz)
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

//checks if cell is inside main domain
__device__ __host__ int is_cell_indomain_device(int ix, int iy, int iz)
{
  if(ix>=0 && ix<NX && iy>=0 && iy<NY && iz>=0 && iz<NZ)
    return 1;
  else
    return 0;
}


// TODO replace get_x, get_xb everywhere

// get grid coordinate at the cell center indexed ic in dimeinsion idim
// copied from get_x macro in ko.h
__device__ __host__ ldouble get_x_device(ldouble* x_arr, int ic, int idim)
{
  ldouble x_out;
  x_out = (idim==0 ? x_arr[ic+NG] :		     
          (idim==1 ? x_arr[ic+NG + NX+2*NG] :  
	  (idim==2 ? x_arr[ic+NG + NX+2*NG + NY+2*NG ] : 0.)));

  return x_out;
}

// get grid coordinate on the cell wall indexed ic in dimension idim
// copied from get_xb macro in ko.h
__device__ __host__ ldouble get_xb_device(ldouble* xb_arr, int ic, int idim)
{
  ldouble xb_out;
  xb_out = (idim==0 ? xb_arr[ic+NG] :		     
           (idim==1 ? xb_arr[ic+NG + NX+2*NG + 1] :  
	   (idim==2 ? xb_arr[ic+NG + NX+2*NG +1 + NY+2*NG +1 ] : 0.)));

  return xb_out;
}

// get size of cell indexed ic in dimension idim
// copied from get_size_x in finite.c
__device__ __host__ ldouble get_size_x_device(ldouble* xb_arr, int ic, int idim)
{
  ldouble dx;
  dx = get_xb_device(xb_arr, ic+1,idim) - get_xb_device(xb_arr, ic, idim);
  return dx;
}


//**********************************************************************//
//Interpolates primitives to the left and right walls of current cell i
//um2, um1, u0, up1, up2 values of primitive in the i-2, i-1, i, i+1, i+2 cells
//ul, ur interpolated primitives at the left and right walls of cell i
//dxm2, dxm1, dx0, dxp1, dxp2 sizes of the five cells
//reconstrpar -- overrides standard reconstruction to reduce to donor cell
//theta MINMOD_THETA

//Several interpolation schemes are available to define at compilation
 
//INT_ORDER=0: basic donor cell
//INT_ORDER=1: Minmod (FLUXLIMITER=0), Monotonized Central (FLUXLIMITER=1), Superbee (FLUXLIMITER=4)
//INT_ORDER=2: Piecewise Parabolic Method (PPM)
 
//**********************************************************************//
__device__ __host__ int avg2point_device(ldouble *um2,ldouble *um1,ldouble *u0,ldouble *up1,ldouble *up2,
	                                 ldouble *ul,ldouble *ur,
	                                 ldouble dxm2,ldouble dxm1,ldouble dx0,ldouble dxp1,ldouble dxp2,
	                                 int reconstrpar,ldouble theta)
{
  
  if(INT_ORDER==0 || reconstrpar==1) // donor cell, no interpolation
  { 
    for(int iv=0;iv<NV;iv++)
    {
      ur[iv]=u0[iv];
      ul[iv]=u0[iv];
    }
  }  
  
  else if(INT_ORDER==1) // linear interpolation 
  {
    
    for(int iv=0;iv<NV;iv++)
    {
      // Slope limiter code rewritten by Ramesh. No function-calls, no divisions.
      ldouble slope;
      ldouble deltam = u0[iv]-um1[iv];
      ldouble deltap = up1[iv]-u0[iv];
      
      if (deltam * deltap <= 0.)
      {
        // We are at a local maximum or minimum. Use zero slope (i.e., donor cell)
        ur[iv] = u0[iv];
        ul[iv] = u0[iv];
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
            return -2;
	    //exit(-2); // TODO
          }
          else if (FLUXLIMITER == 3) // Koren -- discouraged since it is not symmetric
          {
            printf("Error: Koren slope limiter is discouraged since it is not symmetric\n");
            return -3;
	    //exit(-3); // TODO
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
	    return -2;
            //exit(-2); // TODO 
          }
          else if (FLUXLIMITER == 3) // Koren -- discouraged since it is not symmetric
          {
            printf("Error: Koren slope limiter is discouraged since it is not symmetric\n");
	    return -3;
            //exit(-3); // TODO
          }
          else // Superbee
          {
            slope = my_min(my_max(2*deltam, deltap), my_max(deltam, 2*deltap));
          }
        }
        
        ur[iv] = u0[iv] + 0.5*slope;
        ul[iv] = u0[iv] - 0.5*slope;
      }
      
      if(isnan(ur[iv]) || isnan (ul[iv]))
	printf("%d %e %e %e %e %e\n",iv,um2[iv],um1[iv],u0[iv],up1[iv],up2[iv]);
    }  // for(int iv=0;iv<NV;iv++)
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
    ldouble dri[NV],drim1[NV],drip1[NV];
    
    for(int iv=0;iv<NV;iv++)
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
    
    for(int iv=0;iv<NV;iv++)
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
    
    for(int iv=0;iv<NV;iv++)
    {
      // This is a_{j-1/2} in eq (1.6) of Colella & Woodward
      ul[iv] = um1[iv] + dxm1 * dx0_plus_dxm1_inv * (u0[iv]-um1[iv]) +
            DX_inv * ((2.*dx0*dxm1) * dx0_plus_dxm1_inv * (Z1-Z2) * (u0[iv]-um1[iv]) -
            dxm1 * Z1 * dri[iv] + dx0 * Z2 * drim1[iv]);
    }
    
    // Make sure that the parabola remains monotonic.
    // The following is equivalent to eq (1.10) in C&W, though it looks different  
    for(int iv=0;iv<NV;iv++)
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
    }
    
  }  // else if(INT_ORDER==2)

  return 0;
}


//**********************************************************************
// kernels
//**********************************************************************

// TODO wrap up xyz wavespeeds in their own array
__global__ void calc_wavespeeds_kernel(int Nloop_1,
				       int* loop_1_ix, int* loop_1_iy, int* loop_1_iz,
				       ldouble* x_arr, ldouble* xb_arr,
				       ldouble* g_arr, ldouble* G_arr,
				       ldouble* p_arr,
				       ldouble* ahdxl_arr, ldouble* ahdyl_arr, ldouble* ahdzl_arr,
				       ldouble* ahdxr_arr, ldouble* ahdyr_arr, ldouble* ahdzr_arr,
				       ldouble* ahdx_arr,  ldouble* ahdy_arr,  ldouble* ahdz_arr,	   
				       ldouble* aradxl_arr, ldouble* aradyl_arr, ldouble* aradzl_arr,
				       ldouble* aradxr_arr, ldouble* aradyr_arr, ldouble* aradzr_arr,
				       ldouble* aradx_arr,  ldouble* arady_arr,  ldouble* aradz_arr,
				       ldouble* cell_tstepden_arr, ldouble* tstepdenmin_ptr, ldouble* tstepdenmax_ptr)

{
  
  // get index for this thread
  // Nloop_1 is number of cells to compute bcs for
  // domain plus lim (=1 usually) ghost cells
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_1) return;
    
  // get indices from 1D arrays
  int ix=loop_1_ix[ii];
  int iy=loop_1_iy[ii];
  int iz=loop_1_iz[ii]; 

  // TODO MPI
#ifndef MPI4CORNERS
  if(if_outsidegc(ix,iy,iz)==1) return; //avoid corners
#endif
    
  // currently is_cell_active always returns 1
  if(!is_cell_active_device(ix,iy,iz)) return;

  // get geometry
  struct geometry geom;
  fill_geometry_device(ix,iy,iz,&geom, x_arr,g_arr,G_arr);
  
  // get primitives 
  ldouble pp[NV];
  for(int iv=0;iv<NV;iv++)
    pp[iv]=get_u(p_arr,iv,ix,iy,iz);

  ldouble aaa[18];
  calc_wavespeeds_lr_pure_device(pp,&geom,aaa);

  // save the wavespeeds
  // NOTE: this part was in save_wavespeeds(), which is not used elsewhere

  // hydro wavespeeds
  set_u_scalar(ahdxl_arr,ix,iy,iz,aaa[0]);
  set_u_scalar(ahdxr_arr,ix,iy,iz,aaa[1]);
  set_u_scalar(ahdyl_arr,ix,iy,iz,aaa[2]);
  set_u_scalar(ahdyr_arr,ix,iy,iz,aaa[3]);
  set_u_scalar(ahdzl_arr,ix,iy,iz,aaa[4]);
  set_u_scalar(ahdzr_arr,ix,iy,iz,aaa[5]);
  
  ldouble aaaxhd=my_max(fabs(aaa[0]),fabs(aaa[1]));
  ldouble aaayhd=my_max(fabs(aaa[2]),fabs(aaa[3]));
  ldouble aaazhd=my_max(fabs(aaa[4]),fabs(aaa[5]));

  set_u_scalar(ahdx_arr,ix,iy,iz,aaaxhd);
  set_u_scalar(ahdy_arr,ix,iy,iz,aaayhd);
  set_u_scalar(ahdz_arr,ix,iy,iz,aaazhd);

  // TODO -- put behind RADIATION gate? 
  // radiative wavespeeds in slots[12::] are limited by the optical depth
  // used to calculate the fluxes
  set_u_scalar(aradxl_arr,ix,iy,iz,aaa[12]);
  set_u_scalar(aradxr_arr,ix,iy,iz,aaa[13]);
  set_u_scalar(aradyl_arr,ix,iy,iz,aaa[14]);
  set_u_scalar(aradyr_arr,ix,iy,iz,aaa[15]);
  set_u_scalar(aradzl_arr,ix,iy,iz,aaa[16]);
  set_u_scalar(aradzr_arr,ix,iy,iz,aaa[17]);

  ldouble aaaxrad=my_max(fabs(aaa[12]),fabs(aaa[13]));
  ldouble aaayrad=my_max(fabs(aaa[14]),fabs(aaa[15]));
  ldouble aaazrad=my_max(fabs(aaa[16]),fabs(aaa[17]));
  
  set_u_scalar(aradx_arr,ix,iy,iz,aaaxrad);
  set_u_scalar(arady_arr,ix,iy,iz,aaayrad);
  set_u_scalar(aradz_arr,ix,iy,iz,aaazrad);

  // search for the maximal unlimited wavespeed to set the timestep
  // only domain cells
  if(is_cell_indomain_device(ix,iy,iz)==1) 
  {      

    ldouble wsx=aaaxhd;
    ldouble wsy=aaayhd;
    ldouble wsz=aaazhd;

    #ifdef RADIATION
    #ifndef SKIPRADEVOLUTION
    // here use radiative wavespeeds not limited by the optical depth
    // in slots 6-11
    aaaxrad=my_max(fabs(aaa[6]),fabs(aaa[7]));
    aaayrad=my_max(fabs(aaa[8]),fabs(aaa[9]));
    aaazrad=my_max(fabs(aaa[10]),fabs(aaa[11]));

    wsx=my_max(aaaxhd,aaaxrad);
    wsy=my_max(aaayhd,aaayrad);
    wsz=my_max(aaazhd,aaazrad);
    #endif
    #endif

    ldouble dx=get_size_x_device(xb_arr,ix,0);
    ldouble dy=get_size_x_device(xb_arr,iy,1);
    ldouble dz=get_size_x_device(xb_arr,iz,2);

    ldouble tstepden;
    if(NZ>1 && NY>1)
      tstepden=wsx/dx + wsy/dy + wsz/dz;
    else if(NZ==1 && NY>1)
      tstepden=wsx/dx + wsy/dy;
    else if(NY==1 && NZ>1)
      tstepden=wsx/dx + wsz/dz;
    else
      tstepden=wsx/dx;   

    // apply limiter
    tstepden/=TSTEPLIM;

    // set_global_arr
    set_u_scalar(cell_tstepden_arr,ix,iy,iz,tstepden);

    // TODO TODO -- do these atomics work?
    // TODO is there a faster way? 
    //global variables for maximum/minimum (inverse) cell timestep
    atomicMin_double(tstepdenmin_ptr, tstepden);
    atomicMax_double(tstepdenmax_ptr, tstepden);
    //if(tstepden>tstepdenmax) tstepdenmax=tstepden;  
    //if(tstepden<tstepdenmin) tstepdenmin=tstepden;  
  }
  
}


__global__ void calc_interp_kernel(int Nloop_1,
				   int* loop_1_ix, int* loop_1_iy, int* loop_1_iz,
				   ldouble* x_arr, ldouble* xb_arr,
				   ldouble* gbx_arr, ldouble* gby_arr, ldouble* gbz_arr,
				   ldouble* Gbx_arr, ldouble* Gby_arr, ldouble* Gbz_arr,				 
				   ldouble* pbLx_arr, ldouble* pbLy_arr, ldouble* pbLz_arr,
				   ldouble* pbRx_arr, ldouble* pbRy_arr, ldouble* pbRz_arr,
				   ldouble* flLx_arr, ldouble* flLy_arr, ldouble* flLz_arr,
				   ldouble* flRx_arr, ldouble* flRy_arr, ldouble* flRz_arr,				   
				   ldouble* p_arr)
{

  // get index for this thread
  // Nloop_1 is number of cells to compute bcs for
  // domain plus lim (=1 usually) ghost cells
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_1) return;
    
  // get indices from 1D arrays
  int ix=loop_1_ix[ii];
  int iy=loop_1_iy[ii];
  int iz=loop_1_iz[ii]; 

  #ifdef SPECIAL_BC_CHECK
  int giix,giiy,giiz;
  #endif
  int perform_sweep = 1; //Brandon - Added for special conditions where a sweep should not be performed under SPECIAL_BC_CHECK

  // TODO MPI indices
  #ifdef SPECIAL_BC_CHECK
  giix = ix + TOI;
  giiy = iy + TOJ;
  giiz = iz + TOK;
  #endif

  // TODO MPI
  #ifndef MPI4CORNERS
  if(if_outsidegc(ix,iy,iz)==1) return; //avoid corners
  #endif

  // TODO -- what is ret_val? defined in PR_BC_SPECIAL_LOOP?
  #ifdef SPECIAL_BC_CHECK
  #include PR_BC_SPECIAL_LOOP
  if(ret_val == 4) 
  {
    return; //Exclude 'corners' in stream region
  }
  #endif

  //create arrays for interpolating conserved quantities
  struct geometry geom;
  ldouble fd_p0[NV],fd_pp1[NV],fd_pp2[NV],fd_pm1[NV],fd_pm2[NV];
  ldouble fd_pr[NV],fd_pl[NV];
  ldouble ffl[NV],ffr[NV]; 
  ldouble dx0, dxm2, dxm1, dxp1, dxp2;
  ldouble minmod_theta=MINMOD_THETA;
  int reconstrpar;
  int dol,dor;

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

    // TODO - check these !
#ifdef SPECIAL_BC_CHECK //Don't do l/r fluxes when at GC - Brandon
#ifndef SEARCH_STREAM_BOUNDARY
    if(TNY>1 && TNZ==1)
    {
      if(giiy >= STREAM_IYT && giiy <= STREAM_IYB)
      {
        if(giix == STREAM_IX || giix == (STREAM_IX+1) || giix == (STREAM_IX+2) || giix == (STREAM_IX+3))
	  dor=0;
      }
    }
    else if(TNY==1 && TNZ>1)
    {
      if(giiz >= STREAM_IZT && giiz <= STREAM_IZB)
      {
        if(giix == STREAM_IX || giix == (STREAM_IX+1) || giix == (STREAM_IX+2) || giix == (STREAM_IX+3))
	  dor=0;
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
#endif //SPECIAL_BC_CHECK

    // skip flux calculation if not needed
    // is_cell_active is currently always 1
    if((ix==0 && is_cell_active_device(ix,iy,iz)==0) || (ix>0 && is_cell_active_device(ix,iy,iz)==0 && is_cell_active_device(ix-1,iy,iz)==0))
      dol=0;
    if((ix==NX-1 && is_cell_active_device(ix,iy,iz)==0) || (ix<NX-1 && is_cell_active_device(ix,iy,iz)==0 && is_cell_active_device(ix+1,iy,iz)==0))
      dor=0;
			      
    // dx0, dxm1, dxp1 are x-sizes (wall to wall) of cells ix, ix-1, ix+1, dxm2m, dxp2 are sizes of cells ix-2, ix+2		
    dx0  = get_size_x_device(xb_arr,ix,  0);    
    dxm1 = get_size_x_device(xb_arr,ix-1,0);    
    dxp1 = get_size_x_device(xb_arr,ix+1,0);    
	  
    if(INT_ORDER>1)
    {
      dxm2 = get_size_x_device(xb_arr,ix-2,0);    
      dxp2 = get_size_x_device(xb_arr,ix+2,0);    
    }

    //fd_p0, fd_pp1, fd_pm1 are primitives at current, left and right cells, fd_pm2, fd_pp2 are for next two cells   
    for(int iv=0;iv<NV;iv++)
    {
      fd_p0[iv]  = get_u(p_arr,iv,ix,  iy,iz);
      fd_pp1[iv] = get_u(p_arr,iv,ix+1,iy,iz);
      fd_pm1[iv] = get_u(p_arr,iv,ix-1,iy,iz);
            
      if(INT_ORDER>1)
      {
        fd_pm2[iv] = get_u(p_arr,iv,ix-2,iy,iz);
        fd_pp2[iv] = get_u(p_arr,iv,ix+2,iy,iz);
      }
    }

    reconstrpar=0;
    minmod_theta=MINMOD_THETA;
 
    // TODO -- reduce order
#ifdef REDUCEORDERWHENNEEDED
    reconstrpar = reduce_order_check(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif

    // TODOD -- reduce minmod
#ifdef REDUCEMINMODTHETA  // reduce minmod_theta near axis or inner boundary
    minmod_theta = reduce_minmod_theta(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif

    // Interpolate primitives to the left and right walls of current cell: fd_pl, fd_pr
    avg2point_device(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2,reconstrpar,minmod_theta);   

    // if(ix>0)
    if(dol) //no need to calculate at left face of first GC if dol=0
    {
      // Left wall of current cell: compute fluxes and save in array ffl[NV]
      fill_geometry_face_device(ix,iy,iz,0,&geom,
                                x_arr,xb_arr,gbx_arr,gby_arr,gbz_arr,Gbx_arr,Gby_arr,Gbz_arr);
      check_floors_mhd_device(fd_pl,VELPRIM,&geom);
     
      f_flux_prime_device(fd_pl,0,ffl,&geom);
      //            if(ix==ixTEST && iy==iyTEST && iz==izTEST) printf("gdet ix: %e\n",geom.gdet);
    }

    // if(ix<NX)
    if(dor) //no need to calculate at right face of first GC if dor=0
    {
      // Right wall of current cell: compute fluxes and save in array ffr[NV]
      fill_geometry_face_device(ix+1,iy,iz,0,&geom,
                                x_arr,xb_arr,gbx_arr,gby_arr,gbz_arr,Gbx_arr,Gby_arr,Gbz_arr);
      check_floors_mhd_device(fd_pr,VELPRIM,&geom);
     
      f_flux_prime_device(fd_pr,0,ffr,&geom);
      //            if(ix==ixTEST && iy==iyTEST && iz==izTEST) printf("gdet ix+1: %e\n",geom.gdet);
    }

    //save interpolated values to memory
    //Note that l and r of a given cell ix are the left and right wall of that cell
    //whereas L and R of given ix are quantities to the left and right of wall ix
    for(int iv=0;iv<NV;iv++)
    {
      // Save fd_pl in array pbRx (Primitive_R) of wall ix
      // Save fd_pr in array pbLx (Primitive_L) of wall ix+1
      set_ubx(pbRx_arr,iv,ix,  iy,iz,fd_pl[iv]);
      set_ubx(pbLx_arr,iv,ix+1,iy,iz,fd_pr[iv]);

      // Save ffl in array flRx (F_R) of wall ix     
      // Save ffr in array flLx (F_L) of wall ix+1 
      if(dol)
        set_ubx(flRx_arr,iv,ix,  iy,iz,ffl[iv]);
      if(dor)
        set_ubx(flLx_arr,iv,ix+1,iy,iz,ffr[iv]);
    } 
  }  // if(NX>1 && iy>=0 && iy<NY && iz>=0 && iz<NZ...)
      
  //**********************************************************************
  //y 'sweep'
  //**********************************************************************
  
  perform_sweep = 1;

  // TODO - check these !      
#ifdef SPECIAL_BC_CHECK
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

    // skip flux calculation if not needed
    // is_cell_active is currently always 1
    if((iy==0 && is_cell_active_device(ix,iy,iz)==0) || (iy>0 && is_cell_active_device(ix,iy,iz)==0 && is_cell_active_device(ix,iy-1,iz)==0))
      dol=0;
             
    if((iy==NY-1 && is_cell_active_device(ix,iy,iz)==0) || (iy<NY-1 && is_cell_active_device(ix,iy,iz)==0 && is_cell_active_device(ix,iy+1,iz)==0))
      dor=0;

    // dx0, dxm1, dxp1 are y-sizes (wall to wall) of cells iy, iy-1, iy+1, dxm2m, dxp2 are sizes of cells iy-2, iy+2
    dx0  = get_size_x_device(xb_arr,iy,  1);    
    dxm1 = get_size_x_device(xb_arr,iy-1,1);    
    dxp1 = get_size_x_device(xb_arr,iy+1,1);    
	
    if(INT_ORDER>1)
    {
      dxm2 = get_size_x_device(xb_arr,iy-2,1);  
      dxp2 = get_size_x_device(xb_arr,iy+2,1);
    }

    //fd_p0, fd_pp1, fd_pm1 are primitives at current, left and right cells, fd_pm2, fd_pp2 are for next two cells		  
    for(int iv=0;iv<NV;iv++)
    {
      fd_p0[iv] =  get_u(p_arr,iv,ix,iy,  iz);
      fd_pp1[iv] = get_u(p_arr,iv,ix,iy+1,iz);
      fd_pm1[iv] = get_u(p_arr,iv,ix,iy-1,iz);
      if(INT_ORDER>1)
      {
        fd_pm2[iv]=get_u(p_arr,iv,ix,iy-2,iz);
        fd_pp2[iv]=get_u(p_arr,iv,ix,iy+2,iz);
      }
    }
	  
    reconstrpar=0;
    minmod_theta=MINMOD_THETA;
		
    // TODO -- reduce order
#ifdef REDUCEORDERWHENNEEDED
   reconstrpar = reduce_order_check(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif

   // TODO -- minmod
#ifdef REDUCEMINMODTHETA  // reduce minmod_theta near axis or inner boundary
    minmod_theta = reduce_minmod_theta(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif

    // Interpolate primitives to the left and right walls of current cell: fd_pl, fd_pr
    avg2point_device(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2,reconstrpar,minmod_theta);

    //iy>0
    if(dol) //no need to calculate at left face of first GC if dol=0
    {
      // Left wall of current cell: compute fluxes and save in array ffl[NV]
      fill_geometry_face_device(ix,iy,iz,1,&geom,
 			       x_arr,xb_arr,gbx_arr,gby_arr,gbz_arr,Gbx_arr,Gby_arr,Gbz_arr);
      check_floors_mhd_device(fd_pl,VELPRIM,&geom);
     
     f_flux_prime_device(fd_pl,1,ffl,&geom);
     //           if(ix==ixTEST && iy==iyTEST && iz==izTEST) printf("gdet iy: %e\n",geom.gdet);
   }


    //iy<NY
    if(dor) //no need to calculate at right face of first GC if dor=0
    {

      // Right wall of current cell: compute fluxes and save in array ffr[NV]
      fill_geometry_face_device(ix,iy+1,iz,1,&geom,
 			       x_arr,xb_arr,gbx_arr,gby_arr,gbz_arr,Gbx_arr,Gby_arr,Gbz_arr);
      check_floors_mhd_device(fd_pr,VELPRIM,&geom);
     
      f_flux_prime_device(fd_pr,1,ffr,&geom);
     //           if(ix==ixTEST && iy==iyTEST && iz==izTEST) printf("gdet iy+1: %e\n",geom.gdet);
    }

    //save interpolated values to memory
    //Note that l and r of a given cell iy are the left and right wall of that cell,
    //whereas L and R of given iy are quantities to the left and right of wall iy
    for(int iv=0;iv<NV;iv++)
    {
      // Save fd_pl in array pbRy (Primitive_R) of wall iy
      // Save fd_pr in array pbLy (Primitive_L) of wall iy+1
      set_uby(pbRy_arr,iv,ix,iy,  iz,fd_pl[iv]);
      set_uby(pbLy_arr,iv,ix,iy+1,iz,fd_pr[iv]);
 
      // Save ffl in array flRy (F_R) of wall iy
      // Save ffr in array flLy (F_R) of wall iy+1
      if(dol)              
        set_uby(flRy_arr,iv,ix,iy,  iz,ffl[iv]);
      if(dor)
        set_uby(flLy_arr,iv,ix,iy+1,iz,ffr[iv]);
    } 
 }  // if(NY>1 && ix>=0 && ix<NX && iz>=0 && iz<NZ...)


 //**********************************************************************
 //z 'sweep'
 //**********************************************************************
	      
  perform_sweep = 1;

//TODO -- check these! 
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

    // skip flux calculation if not needed    
    // is_cell_active is currently always 1
    if((iz==0 && is_cell_active_device(ix,iy,iz)==0) || (iz>0 && is_cell_active_device(ix,iy,iz)==0 && is_cell_active_device(ix,iy,iz-1)==0))
      dol=0;
    if((iz==NZ-1 && is_cell_active_device(ix,iy,iz)==0) || (iz<NZ-1 && is_cell_active_device(ix,iy,iz)==0 && is_cell_active_device(ix,iy,iz+1)==0))
      dor=0;

    // dx0, dxm1, dxp1 are z-sizes (wall to wall) of cells iz, iz-1, iz+1, dxm2m, dxp2 are sizes of cells iz-2, iz+2
    dx0  = get_size_x_device(xb_arr,iz, 2);
    dxm1 = get_size_x_device(xb_arr,iz-1,2);
    dxp1 = get_size_x_device(xb_arr,iz+1,2);
         
    if(INT_ORDER>1)
    {
      dxm2 = get_size_x_device(xb_arr,iz-2,2);
      dxp2 = get_size_x_device(xb_arr,iz+2,2);
    }

    //fd_p0, fd_pp1, fd_pm1 are primitives at current, left and right cells, fd_pm2, fd_pp2 are for next two cells
    for(int iv=0;iv<NV;iv++)
    {
      fd_p0[iv]  = get_u(p_arr,iv,ix,iy,iz);
      fd_pp1[iv] = get_u(p_arr,iv,ix,iy,iz+1);
      fd_pm1[iv] = get_u(p_arr,iv,ix,iy,iz-1);
           
      if(INT_ORDER>1)
      {
        fd_pm2[iv]=get_u(p_arr,iv,ix,iy,iz-2);
        fd_pp2[iv]=get_u(p_arr,iv,ix,iy,iz+2);
      }
    }
         
    reconstrpar=0;
    minmod_theta=MINMOD_THETA;
		
    // TODO -- reduce order
#ifdef REDUCEORDERWHENNEEDED
    reconstrpar = reduce_order_check(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif
         
    // TODO -- minmod
#ifdef REDUCEMINMODTHETA  // reduce minmod_theta near axis or inner boundary
    minmod_theta = reduce_minmod_theta(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,ix,iy,iz);
#endif
         
    // Interpolate primitives to the left and right walls of current cell: fd_pl, fd_pr
    avg2point_device(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2,reconstrpar,minmod_theta);

    //iz>0
    if(dol) //no need to calculate at left face of first GC if dol=0
    {
      // Left wall of current cell: compute fluxes and save in array ffl[NV]
      fill_geometry_face_device(ix,iy,iz,2,&geom,
			        x_arr,xb_arr,gbx_arr,gby_arr,gbz_arr,Gbx_arr,Gby_arr,Gbz_arr);
      check_floors_mhd_device(fd_pl,VELPRIM,&geom);   
      f_flux_prime_device(fd_pl,2,ffl,&geom);
      //      if(ix==ixTEST && iy==iyTEST && iz==izTEST) printf("gdet iz: %e\n",geom.gdet);
    }

    //iz<NZ
    if(dor) //no need to calculate at right face of first GC if dor=0
    {
      // Right wall of current cell: compute fluxes and save in array ffr[NV]
      fill_geometry_face_device(ix,iy,iz+1,2,&geom,
			        x_arr,xb_arr,gbx_arr,gby_arr,gbz_arr,Gbx_arr,Gby_arr,Gbz_arr);
      check_floors_mhd_device(fd_pr,VELPRIM,&geom);   
      f_flux_prime_device(fd_pr,2,ffr,&geom);
      //      if(ix==ixTEST && iy==iyTEST && iz==izTEST) printf("gdet iz+1: %e\n",geom.gdet);
    }
         
    //save interpolated values to memory
    //Note that l and r of a given cell iy are the left and right wall of that cell,
    //whereas L and R of given iy are quantities to the left and right of wall iy
    for(int iv=0;iv<NV;iv++)
    {
      // Save fd_pl in array pbRz (Primitive_R) of wall iz
      // Save fd_pr in array pbLz (Primitive_L) of wall iz+1
      set_ubz(pbRz_arr,iv,ix,iy,iz,  fd_pl[iv]);
      set_ubz(pbLz_arr,iv,ix,iy,iz+1,fd_pr[iv]);

      // Save ffl in array flRz (F_R) of wall iz
      // Save ffr in array flLz (F_R) of wall iz+1
      if(dol)
        set_ubz(flRz_arr,iv,ix,iy,iz,  ffl[iv]);
      if(dor)
        set_ubz(flLz_arr,iv,ix,iy,iz+1,ffr[iv]);

      if(ix==ixTEST && iy==iyTEST && iz==izTEST)
	{
	  //printf("%d %d %e %e %e %e\n",dol, dor, fd_pl[iv],fd_pr[iv],ffl[iv],ffr[iv]);
	}
    }
  }  // if(NZ>1 && ix>=0 && ix<NX && iy>=0 && iy<NY...)
	     
}


//************************************************************************//
/*! \fn __global__ void calc_fluxes_kernel()
 \brief Calculates fluxes at cell faces
 
 Uses the selected approximate Riemann solver to estimate the fluxes on the six faces of the cell (ix, iy, iz)\n
 The fluxes are saved in global arrays flbx, flby, flbz\n
 Note that the fluxes calculated here correspond to the three left walls of cell ix, iy, iz
 
*/
//************************************************************************//
// TODO -- wrap up xyz wavespeeds and xyz boundary metric and xyz primitives and xyz LR fluxes in single arrays
__global__ void calc_fluxes_kernel(int Nloop_1,
                                   int* loop_1_ix, int* loop_1_iy, int* loop_1_iz,
		       	           ldouble* x_arr, ldouble* xb_arr,
				   ldouble* gbx_arr, ldouble* gby_arr, ldouble* gbz_arr,
				   ldouble* Gbx_arr, ldouble* Gby_arr, ldouble* Gbz_arr,
				   ldouble* pbLx_arr, ldouble* pbLy_arr, ldouble* pbLz_arr,
				   ldouble* pbRx_arr, ldouble* pbRy_arr, ldouble* pbRz_arr,
				   ldouble* flLx_arr, ldouble* flLy_arr, ldouble* flLz_arr,
				   ldouble* flRx_arr, ldouble* flRy_arr, ldouble* flRz_arr,
				   ldouble* ahdxl_arr, ldouble* ahdyl_arr, ldouble* ahdzl_arr,
				   ldouble* ahdxr_arr, ldouble* ahdyr_arr, ldouble* ahdzr_arr,
				   ldouble* ahdx_arr,  ldouble* ahdy_arr,  ldouble* ahdz_arr,	   
				   ldouble* aradxl_arr, ldouble* aradyl_arr, ldouble* aradzl_arr,
				   ldouble* aradxr_arr, ldouble* aradyr_arr, ldouble* aradzr_arr,
				   ldouble* aradx_arr,  ldouble* arady_arr,  ldouble* aradz_arr,	 
				   ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr)
{

  // get index for this thread
  // Nloop_1 is number of cells to compute bcs for
  // domain plus lim (=1 usually) ghost cells
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_1) return;
    
  // get indices from 1D arrays
  int ix=loop_1_ix[ii];
  int iy=loop_1_iy[ii];
  int iz=loop_1_iz[ii]; 

  struct geometry geom; 
  ldouble ag,al,ar;
  ldouble am1l[2],am1r[2],am1[2];
  ldouble ap1l[2],ap1r[2],ap1[2];
  ldouble fd_pL[NV], fd_pR[NV], fd_uL[NV],fd_uR[NV];
  

  // flbx[NV], flby[NV], flbz[NV] are the fluxes at the three walls under consideration
  // set all to 0 to start
  for(int iv=0;iv<NV;iv++) 
  {
    set_ubx(flbx_arr,iv,ix,iy,iz,0.);
    set_uby(flby_arr,iv,ix,iy,iz,0.);
    set_ubz(flbz_arr,iv,ix,iy,iz,0.);
  }

  //**********************************************************************//
  // Work on the x-face at ix, iy, iz, which lies in between cells (ix-1,iy,iz) and (ix,iy,iz)
#ifdef MPI4CORNERS
  if(NX>1 && ix>=0 && ix<=NX && iy>=-1 && iy<NY+1 && iz>=-1 && iz<NZ+1)
#else
  if(NX>1 && ix>=0 && ix<=NX && iy>=0 && iy<NY && iz>=0 && iz<NZ)
#endif
  {

    // fd_pL, fd_pR are the left-biased and right-biased primitives at the current cell face
    for(int iv=0;iv<NV;iv++)
    {      
      fd_pL[iv] = get_ub(pbLx_arr,iv,ix,iy,iz,0);
      fd_pR[iv] = get_ub(pbRx_arr,iv,ix,iy,iz,0);
    }

    // fd_uL, fd_uR are the left-biased and right-biased conserveds at the current cell face
    fill_geometry_face_device(ix,iy,iz,0,&geom,x_arr,xb_arr,
			      gbx_arr,gby_arr,gbz_arr,Gbx_arr,Gby_arr,Gbz_arr);
    p2u_device(fd_pL,fd_uL,&geom);
    p2u_device(fd_pR,fd_uR,&geom);
    
    // Characteristic wave speeds in the two adjoining cells of the current face,
    // ap1, am1 correspond to ix and ix-1, i.e., speeds on the right and left of the current face;
    // l and r correspond to left-going and right-going waves; if neither l nor r, it is the maximum speed
    // [0], [1] correspond to hydro and radiation wave speeds      

    ap1l[0] = get_u_scalar(ahdxl_arr,ix,iy,iz);
    ap1r[0] = get_u_scalar(ahdxr_arr,ix,iy,iz);
    ap1[0]  = get_u_scalar(ahdx_arr, ix,iy,iz);
   
    am1l[0] = get_u_scalar(ahdxl_arr,ix-1,iy,iz);
    am1r[0] = get_u_scalar(ahdxr_arr,ix-1,iy,iz);
    am1[0]  = get_u_scalar(ahdx_arr, ix-1,iy,iz);

    ap1l[1] = get_u_scalar(aradxl_arr,ix,iy,iz);
    ap1r[1] = get_u_scalar(aradxr_arr,ix,iy,iz);
    ap1[1]  = get_u_scalar(aradx_arr, ix,iy,iz);
    
    am1l[1] = get_u_scalar(aradxl_arr,ix-1,iy,iz);
    am1r[1] = get_u_scalar(aradxr_arr,ix-1,iy,iz);
    am1[1]  = get_u_scalar(aradx_arr, ix-1,iy,iz);
	
    /* //TODO
#ifdef WAVESPEEDSATFACES //re-calculate the wavespeeds directly at the face
    ldouble aaa[18];
    //left biased wavespeeds
    calc_wavespeeds_lr_pure_device(fd_uL,&geom,aaa);
    am1l[0]=aaa[0];
    am1r[0]=aaa[1];
    am1l[1]=aaa[6];
    am1r[1]=aaa[7];
	    
    am1[0]=my_max(fabs(aaa[0]),fabs(aaa[1]));
    am1[1]=my_max(fabs(aaa[6]),fabs(aaa[7]));

    //right biased wavespeeds
    calc_wavespeeds_lr_pure_device(fd_uR,&geom,aaa);
    ap1l[0]=aaa[0];
    ap1r[0]=aaa[1];
    ap1l[1]=aaa[6];
    ap1r[1]=aaa[7];

    ap1[0]=my_max(fabs(aaa[0]),fabs(aaa[1]));
    ap1[1]=my_max(fabs(aaa[6]),fabs(aaa[7]));
#endif
    */
	
    // Loop over variables and calculate flux using Lax-Friedrichs or HLL
    for(int iv=0;iv<NV;iv++)
    {
      // Choose the proper characteristic speeds: al (left-going wave), ar (right-going wave), ag (maximum wavespeed)
      // Hydro and radiation are treated as two separate systems
#ifdef RADIATION
      if(iv<NVMHD) //hydro characteristic speed
      {
        ag=my_max(ap1[0], am1[0]);
        al=my_min(ap1l[0],am1l[0]);
        ar=my_max(ap1r[0],am1r[0]);
            
#ifdef BATTERY
        //when radiation battery on - magnetic fields are affected by radiation
#ifdef BATTERYRADWAVESPEEDS
#ifdef BATTERYRADWAVESPEEDSBONLY
        if(iv>=B1 && iv<=B3)
#endif
        {
          ag=my_max(ag,my_max(ap1[1], am1[1]));
          al=my_min(al,my_min(ap1l[1],am1l[1]));
          ar=my_max(ar,my_max(ap1r[1],am1r[1]));
        }
#endif
#endif
      }
      else //radiative characteristic speed
      {
        ag=my_max(ap1[1], am1[1]);
        al=my_min(ap1l[1],am1l[1]);
        ar=my_max(ap1r[1],am1r[1]);
      }
          
#else //no RADIATION -- hydro only
      ag=my_max(ap1[0], am1[0]);
      al=my_min(ap1l[0],am1l[0]);
      ar=my_max(ap1r[0],am1r[0]);
#endif

      ldouble fd_fstar;
      if(FLUXMETHOD==LAXF_FLUX) //Lax-Friedrichs Flux
      {
        //Lax-Friedrichs: Flux = 0.5 * [FR + FL - ag * (UR - UL)]
        fd_fstar = 0.5*(get_ub(flRx_arr,iv,ix,iy,iz,0) + get_ub(flLx_arr,iv,ix,iy,iz,0) - ag * (fd_uR[iv] - fd_uL[iv]));     
        set_ubx(flbx_arr,iv,ix,iy,iz,fd_fstar);	    
      } 
      else if(FLUXMETHOD==HLL_FLUX) //HLL Flux
      {
        if(al>0.)
        {
          fd_fstar = get_ub(flLx_arr,iv,ix,iy,iz,0);
        }
        else if(ar<0.)
        {
          fd_fstar = get_ub(flRx_arr,iv,ix,iy,iz,0);
        }
        else
        {
          //HLL: Flux = [ar * FL - al * FR + al * ar * (UR - UL)] / (ar - al)25
          fd_fstar = (-al*get_ub(flRx_arr,iv,ix,iy,iz,0) + ar*get_ub(flLx_arr,iv,ix,iy,iz,0) + al*ar* (fd_uR[iv] - fd_uL[iv]))/(ar-al);
        }
        set_ubx(flbx_arr,iv,ix,iy,iz,fd_fstar);
      } 
    } // for(iv=0;i<NV;i++)
  }  // if(NX>1 && ix>=0 && ix<=NX && iy>=0 && iy<NY && iz>=0 && iz<NZ...)

  
  //**********************************************************************//
  //Work on the y-face at ix, iy, iz, which lies in between cells ix,iy-1,iz and ix,iy,iz  
#ifdef MPI4CORNERS
  if(NY>1 && iy>=0 && iy<=NY && ix>=-1 && ix<NX+1 && iz>=-1 && iz<NZ+1)
#else
  if(NY>1 && iy>=0 && iy<=NY  && ix>=0 && ix<NX && iz>=0 && iz<NZ)
#endif
  {
    // fd_pL, fd_pR are the left-biased and right-biased primitives at the current cell face
    for(int iv=0;iv<NV;iv++)
    {
      fd_pL[iv]=get_ub(pbLy_arr,iv,ix,iy,iz,1);
      fd_pR[iv]=get_ub(pbRy_arr,iv,ix,iy,iz,1);
    }
    
    // fd_uL, fd_uR are the left-biased and right-biased conserveds at the current cell face
    fill_geometry_face_device(ix,iy,iz,1,&geom,x_arr,xb_arr,
			      gbx_arr,gby_arr,gbz_arr,Gbx_arr,Gby_arr,Gbz_arr);
    p2u_device(fd_pL,fd_uL,&geom);
    p2u_device(fd_pR,fd_uR,&geom);

    // Characteristic wave speeds in the two adjoining cells of the current face,
    // ap1, am1 correspond to ix and ix-1, i.e., speeds on the right and left of the current face;
    // l and r correspond to left-going and right-going waves; if neither l nor r, it is the maximum speed
    // [0], [1] correspond to hydro and radiation wave speeds    
    ap1l[0]=get_u_scalar(ahdyl_arr,ix,iy,iz);
    ap1r[0]=get_u_scalar(ahdyr_arr,ix,iy,iz);
    ap1[0] =get_u_scalar(ahdy_arr, ix,iy,iz);
    
    am1l[0]=get_u_scalar(ahdyl_arr,ix,iy-1,iz);
    am1r[0]=get_u_scalar(ahdyr_arr,ix,iy-1,iz);
    am1[0] =get_u_scalar(ahdy_arr,ix,iy-1,iz);    
    
    ap1l[1]=get_u_scalar(aradyl_arr,ix,iy,iz);
    ap1r[1]=get_u_scalar(aradyr_arr,ix,iy,iz);
    ap1[1] =get_u_scalar(arady_arr, ix,iy,iz);
    
    am1l[1]=get_u_scalar(aradyl_arr,ix,iy-1,iz);
    am1r[1]=get_u_scalar(aradyr_arr,ix,iy-1,iz);
    am1[1] =get_u_scalar(arady_arr, ix,iy-1,iz);

    /* //TODO
#ifdef WAVESPEEDSATFACES // recompute wavespeeds directly at face
    ldouble aaa[18];
    //left-biased wavespeeds
    calc_wavespeeds_lr_pure(fd_uL,&geom,aaa);
    am1l[0]=aaa[2];
    am1r[0]=aaa[3];
    am1l[1]=aaa[8];
    am1r[1]=aaa[9];
    am1[0]=my_max(fabs(aaa[2]),fabs(aaa[3]));
    am1[1]=my_max(fabs(aaa[8]),fabs(aaa[9]));

    //right-biased wavespeeds
    calc_wavespeeds_lr_pure(fd_uR,&geom,aaa);
    ap1l[0]=aaa[2];
    ap1r[0]=aaa[3];
    ap1l[1]=aaa[8];
    ap1r[1]=aaa[9];
    ap1[0]=my_max(fabs(aaa[2]),fabs(aaa[3]));
    ap1[1]=my_max(fabs(aaa[8]),fabs(aaa[9]));
#endif
    */	    

    // Loop over variables and calculate flux using Lax-Friedrichs or HLL, as required
    for(int iv=0;iv<NV;iv++)
    {
      // Choose the proper characteristic speeds: al (left-going wave), ar (right-going wave), ag (maximum wavespeed)
      // Hydro and radiation are treated as two separate systems
#ifdef RADIATION
      if(iv<NVMHD) //hydro characteristic speeds
      {
        ag=my_max(ap1[0],  am1[0]);
        al=my_min(ap1l[0],am1l[0]);
        ar=my_max(ap1r[0],am1r[0]);
              
#ifdef BATTERY
      //when radiation battery on - magnetic fields are affected by radiation
#ifdef BATTERYRADWAVESPEEDS
#ifdef BATTERYRADWAVESPEEDSBONLY
       if(vi>=B1 && iv<=B3)
#endif
       {
         ag=my_max(ag,my_max(ap1[1], am1[1]));
         al=my_min(al,my_min(ap1l[1],am1l[1]));
         ar=my_max(ar,my_max(ap1r[1],am1r[1]));
       }
#endif
#endif
     }
     else //radiative characteristic speeds
     {
       ag=my_max(ap1[1], am1[1]);
       al=my_min(ap1l[1],am1l[1]);
       ar=my_max(ap1r[1],am1r[1]);
     }
            
#else //no RADIATION -- use hydro  wavespeeds
     ag=my_max(ap1[0], am1[0]);
     al=my_min(ap1l[0],am1l[0]);
     ar=my_max(ap1r[0],am1r[0]);
#endif

     ldouble fd_fstar;
     if(FLUXMETHOD==LAXF_FLUX) //Lax-Friedrichs Flux
     {
       //Lax-Friedrichs: Flux = 0.5 * [FR + FL - ag * (UR - UL)]
       fd_fstar = 0.5*(get_ub(flRy_arr,iv,ix,iy,iz,1) + get_ub(flLy_arr,iv,ix,iy,iz,1) - ag * (fd_uR[iv] - fd_uL[iv])); 
       set_uby(flby_arr,iv,ix,iy,iz,fd_fstar);
     }
     else if(FLUXMETHOD==HLL_FLUX) //HLL Flux
     {
       if(al>0.)
       {
         fd_fstar = get_ub(flLy_arr,iv,ix,iy,iz,1);
       }
       else if(ar<0.)
       {
         fd_fstar = get_ub(flRy_arr,iv,ix,iy,iz,1);
       }
       else
       {
         //HLL: Flux = [ar * FL - al * FR + al * ar * (UR - UL)] / (ar - al)
         fd_fstar = (-al*get_ub(flRy_arr,iv,ix,iy,iz,1) + ar*get_ub(flLy_arr,iv,ix,iy,iz,1) + al*ar* (fd_uR[iv] - fd_uL[iv]))/(ar-al);
       }
       set_uby(flby_arr,iv,ix,iy,iz,fd_fstar);
     } 
   }  // for(iv=0;i<NV;i++)
 }  // if(NY>1 && iy>=0 && iy<=NY  && ix>=0 && ix<NX && iz>=0 && iz<NZ...)


  //**********************************************************************//
  // Work on the z-face at ix, iy, iz, which lies in between cells ix,iy,iz-1 and ix,iy,iz  
#ifdef MPI4CORNERS
  if(NZ>1 && iz>=0 && iz<=NZ && ix>=-1 && ix<NX+1 && iy>=-1 && iy<NY+1)
#else
  if(NZ>1 && iz>=0 && iz<=NZ && ix>=0 && ix<NX && iy>=0 && iy<NY)
#endif
  {
    // fd_pL, fd_pR are the left-biased and right-biased primitives at the current cell face  
    for(int iv=0;iv<NV;iv++)
    {
      fd_pL[iv]=get_ub(pbLz_arr,iv,ix,iy,iz,2);
      fd_pR[iv]=get_ub(pbRz_arr,iv,ix,iy,iz,2);
    }

    // fd_uL, fd_uR are the left-biased and right-biased conserveds at the current cell face
    fill_geometry_face_device(ix,iy,iz,2,&geom,x_arr,xb_arr,
			      gbx_arr,gby_arr,gbz_arr,Gbx_arr,Gby_arr,Gbz_arr);
    p2u_device(fd_pL,fd_uL,&geom);
    p2u_device(fd_pR,fd_uR,&geom);

    // Characteristic wave speeds in the two adjoining cells of the current face,
    // ap1, am1 correspond to ix and ix-1, i.e., speeds on the right and left of the current face;
    // l and r correspond to left-going and right-going waves; if neither l nor r, it is the maximum speed
    // [0], [1] correspond to hydro and radiation wave speeds

    ap1l[0]=get_u_scalar(ahdzl_arr,ix,iy,iz);
    ap1r[0]=get_u_scalar(ahdzr_arr,ix,iy,iz);
    ap1[0] =get_u_scalar(ahdz_arr, ix,iy,iz);

    am1l[0]=get_u_scalar(ahdzl_arr,ix,iy,iz-1);
    am1r[0]=get_u_scalar(ahdzr_arr,ix,iy,iz-1);
    am1[0] =get_u_scalar(ahdz_arr, ix,iy,iz-1);
	
    ap1l[1]=get_u_scalar(aradzl_arr,ix,iy,iz);
    ap1r[1]=get_u_scalar(aradzr_arr,ix,iy,iz);
    ap1[1] =get_u_scalar(aradz_arr, ix,iy,iz);

    am1l[1]=get_u_scalar(aradzl_arr,ix,iy,iz-1);
    am1r[1]=get_u_scalar(aradzr_arr,ix,iy,iz-1);
    am1[1] =get_u_scalar(aradz_arr, ix,iy,iz-1);

	/* //TODO
#ifdef WAVESPEEDSATFACES // recompute wavespeeds directly at face
    ldouble aaa[18];
    //left-biased wavespeeds
    calc_wavespeeds_lr_pure(fd_uL,&geom,aaa);
    am1l[0]=aaa[4];
    am1r[0]=aaa[5];
    am1l[1]=aaa[10];
    am1r[1]=aaa[11];
    am1[0]=my_max(fabs(aaa[4]),fabs(aaa[5]));
    am1[1]=my_max(fabs(aaa[10]),fabs(aaa[11]));

    //right-biased wavespeeds
    calc_wavespeeds_lr_pure(fd_uR,&geom,aaa);
    ap1l[0]=aaa[4];
    ap1r[0]=aaa[5];
    ap1l[1]=aaa[10];
    ap1r[1]=aaa[11];
    ap1[0]=my_max(fabs(aaa[4]),fabs(aaa[5]));
    ap1[1]=my_max(fabs(aaa[10]),fabs(aaa[11]));
#endif
	*/
	
    // Loop over variables and calculate flux using Lax-Friedrichs or HLL, as required
    for(int iv=0;iv<NV;iv++)
    {
      // Choose the proper characteristic speeds: al (left-going wave), ar (right-going wave), ag (maximum wavespeed)
      // Hydro and radiation are treated as two separate systems
#ifdef RADIATION
      if(iv<NVMHD) // hydro characteristic speeds
      {
        ag=my_max(ap1[0], am1[0]);
        al=my_min(ap1l[0],am1l[0]);
        ar=my_max(ap1r[0],am1r[0]);
              
#ifdef BATTERY
        //when radiation battery on - magnetic fields are affected by radiation
#ifdef BATTERYRADWAVESPEEDS
#ifdef BATTERYRADWAVESPEEDSBONLY
        if(iv>=B1 && iv<=B3)
#endif
        {
          ag=my_max(ag,my_max(ap1[1], am1[1]));
          al=my_min(al,my_min(ap1l[1],am1l[1]));
          ar=my_max(ar,my_max(ap1r[1],am1r[1]));
        }
#endif
#endif
      }
      else //radiative characteristic speeds
      {
        ag=my_max(ap1[1], am1[1]);
        al=my_min(ap1l[1],am1l[1]);
        ar=my_max(ap1r[1],am1r[1]);
      }
#else //no radiation -- hydro characteristic speeds
      ag=my_max(ap1[0], am1[0]);
      al=my_min(ap1l[0],am1l[0]);
      ar=my_max(ap1r[0],am1r[0]);
#endif

      ldouble fd_fstar;
      if(FLUXMETHOD==LAXF_FLUX) //Lax-Friedrichs Flux
      {      
        //Lax-Friedrichs: Flux = 0.5 * [FR + FL - ag * (UR - UL)]
        fd_fstar = 0.5*(get_ub(flRz_arr,iv,ix,iy,iz,2) + get_ub(flLz_arr,iv,ix,iy,iz,2) - ag * (fd_uR[iv] - fd_uL[iv]));      
        set_ubz(flbz_arr,iv,ix,iy,iz,fd_fstar);
      } 
      else if(FLUXMETHOD==HLL_FLUX) //HLL Flux
      {
        if(al>0.)
        {
          fd_fstar = get_ub(flLz_arr,iv,ix,iy,iz,2);
        }
        else if(ar<0.)
        {
          fd_fstar = get_ub(flRz_arr,iv,ix,iy,iz,2);
        }
        else
        {
          //HLL: Flux = [ar * FL - al * FR + al * ar * (UR - UL)] / (ar - al)  
          fd_fstar = (-al*get_ub(flRz_arr,iv,ix,iy,iz,2) + ar*get_ub(flLz_arr,iv,ix,iy,iz,2) + al*ar* (fd_uR[iv] - fd_uL[iv]))/(ar-al);
        }
        set_ubz(flbz_arr,iv,ix,iy,iz,fd_fstar);
      } 
    }  // for(iv=0;iv<NV;iv++)
  }  // if(NZ>1 && iz>=0 && iz<=NZ && ix>=0 && ix<NX && iy>=0 && iy<NY...)
}


__global__ void calc_update_kernel(int Nloop_0, 
                                   int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
		       	           ldouble* x_arr, ldouble* xb_arr,
                                   ldouble* g_arr, ldouble* G_arr, ldouble* gKr_arr,
				   ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr,
				   ldouble* u_arr, ldouble* p_arr, ldouble dtin)
{
  // get index for this thread
  // Nloop_0 is number of cells to update;
  // usually Nloop_0=NX*NY*NZ, but sometimes there are weird bcs inside domain 
  int ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_0) return;
    
  // get indices from 1D arrays
  int ix=loop_0_ix[ii];
  int iy=loop_0_iy[ii];
  int iz=loop_0_iz[ii]; 

  // Source term
  ldouble ms[NV];
  //ldouble gs[NV]; //NOTE gs[NV] is for artifical sources, rarely used
#ifdef NOSOURCES
  for(int iv=0;iv<NV;iv++) ms[iv]=0.;
#else
  if(is_cell_active_device(ix,iy,iz)==0) // NOTE: is_cell_active currently always returns 1 
  {
     // Source terms applied only for active cells	  
     for(int iv=0;iv<NV;iv++) ms[iv]=0.; 
  }
  else
  {
     // Get metric source terms ms[iv]
     // and any other source terms gs[iv]
     struct geometry geom;
     fill_geometry_device(ix,iy,iz,&geom, x_arr,g_arr,G_arr);

     ldouble *pp = &get_u(p_arr,0,ix,iy,iz); 
     f_metric_source_term_device(pp, ms, &geom, gKr_arr);

     //f_metric_source_term_device(ix,iy,iz,ms,p_arr, x_arr,g_arr,G_arr,gKr_arr); // OLD
     //f_general_source_term(ix,iy,iz,gs); //NOTE TODO: *very* rarely used, ignore for now
     //for(int iv=0;iv<NV;iv++) ms[iv]+=gs[iv];
  }
#endif


 if(doTEST==1 && ix==ixTEST && iy==iyTEST && iz==izTEST)
   printf("D ms[NV]: %e %e %e %e %e %e %e %e %e\n", ms[0],ms[1],ms[2],ms[3],ms[4],ms[5],ms[6],ms[7],ms[8]);
  
  // Get the cell size in the three directions
  ldouble dx = get_size_x_device(xb_arr,ix,0); 
  ldouble dy = get_size_x_device(xb_arr,iy,1); 
  ldouble dz = get_size_x_device(xb_arr,iz,2); 
  
  //update all conserved according to fluxes and source terms      
  for(int iv=0;iv<NV;iv++)
  {	

    // Get the initial value of the conserved quantity
    ldouble val = get_u(u_arr,iv,ix,iy,iz);
    
    if(doTEST==1 && ix==ixTEST && iy==iyTEST && iz==izTEST && iv==ivTEST)
      printf("D u: %e\n", val);
    
    // Get the fluxes on the six faces.
    // flbx, flby, flbz are the fluxes at the LEFT walls of cell ix, iy, iz.
    // To get the RIGHT fluxes, we need flbx(ix+1,iy,iz), etc.
    ldouble flxl=get_ub(flbx_arr,iv,ix,iy,iz,0);
    ldouble flxr=get_ub(flbx_arr,iv,ix+1,iy,iz,0);
    ldouble flyl=get_ub(flby_arr,iv,ix,iy,iz,1);
    ldouble flyr=get_ub(flby_arr,iv,ix,iy+1,iz,1);
    ldouble flzl=get_ub(flbz_arr,iv,ix,iy,iz,2);
    ldouble flzr=get_ub(flbz_arr,iv,ix,iy,iz+1,2);

    // Compute Delta U from the six fluxes
    ldouble du = -(flxr-flxl)*dtin/dx - (flyr-flyl)*dtin/dy - (flzr-flzl)*dtin/dz;

    // Compute the new conserved by adding Delta U and the source term
    val += (du + ms[iv]*dtin);

    // Save the new conserved to memory

    //TODO
//#ifdef SKIPHDEVOLUTION
//  if(iv>=NVMHD)
//#endif
//#ifdef RADIATION
//#ifdef SKIPRADEVOLUTION
//#ifdef EVOLVEPHOTONNUMBER
//  if(iv!=EE && iv!=FX && iv!=FY && iv!=FZ && iv!=NF)
//#else
//  if(iv!=EE && iv!=FX && iv!=FY && iv!=FZ)
//#endif
//#endif  
//#endif  
//#ifdef SKIPHDBUTENERGY
//  if(iv>=NVMHD || iv==UU)
//#endif
	
    set_u(u_arr,iv,ix,iy,iz,val);	 

  }  
}


ldouble calc_wavespeeds_gpu()
{
  cudaError_t err = cudaSuccess;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // Determine number of threadblocks
  int threadblocks = (Nloop_1 / TB_SIZE) + ((Nloop_1 % TB_SIZE)? 1:0);
  //printf("\nTest %d\n", threadblocks); fflush(stdout);

  // Launch kernel
  cudaEventRecord(start);

  ldouble* d_tstepdenmin;
  ldoulbe* d_tstepdenmax;
  err = cudaMalloc(&d_tstepdenmin, sizeof(ldouble));
  err = cudaMalloc(&d_tstepdenmax, sizeof(ldouble));
  err = cudaMemcpy(d_tstepdenmin, &tstepdenmin, sizeof(ldouble), cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_tstepdenmax, &tstepdenmax, sizeof(ldouble), cudaMemcpyHostToDevice);
  
  calc_wavespeeds_kernel<<<threadblocks, TB_SIZE>>>(Nloop_1,
                                                    d_loop1_ix, d_loop1_iy, d_loop1_iz,
		                                    d_x, d_xb,
						    d_gcov, d_gcon,
						    d_p_arr,
						    d_ahdxl_arr, d_ahdyl_arr, d_ahdzl_arr,
		                                    d_ahdxr_arr, d_ahdyr_arr, d_ahdzr_arr,
                                                    d_ahdx_arr,  d_ahdy_arr,  d_ahdz_arr,
		                                    d_aradxl_arr, d_aradyl_arr, d_aradzl_arr,
		                                    d_aradxr_arr, d_aradyr_arr, d_aradzr_arr,
		                                    d_aradx_arr,  d_arady_arr,  d_aradz_arr,
		                                    d_cell_tstepden_arr, d_tstepdenmin, d_tstepdenmax
						    );
  cudaEventRecord(stop);
  err = cudaPeekAtLastError();
  // printf("ERROR-Kernel (error code %s)!\n", cudaGetErrorString(err));

  // synchronize
  cudaDeviceSynchronize();
    
  // timing information
  cudaEventSynchronize(stop);
  float tms = 0.;
  cudaEventElapsedTime(&tms, start,stop);
  printf("gpu calc_wavespeeds time: %0.2f \n",tms);

  // TODO replace tstepdenmin_tmp,tstepdenmax_tmp with actual tstepdenmin/max when confident reduction works
  ldouble tstepdenmin_tmp, tstepdenmax_tmp;
  err = cudaMemcpy(&tstepdenmin_tmp, d_tstepdemin, sizeof(ldouble), cudaMemcpyDeviceToHost);
  err = cudaMemcpy(&tstepdenmax_tmp, d_tstepdenmax, sizeof(ldouble), cudaMemcpyDeviceToHost);
  cudaFree(d_tstepdenmin);
  cudaFree(d_tstepdenmax);

  printf("gpu calc wavespeeds tstepdensmin/max: %e %e\n",tstepdenmin_tmp,tstepdenmax_tmp);
  return (ldouble)tms;
}

ldouble calc_interp_gpu()
{
  cudaError_t err = cudaSuccess;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // Determine number of threadblocks
  int threadblocks = (Nloop_1 / TB_SIZE) + ((Nloop_1 % TB_SIZE)? 1:0);
  //printf("\nTest %d\n", threadblocks); fflush(stdout);

  // Launch kernel
  cudaEventRecord(start);
  calc_interp_kernel<<<threadblocks, TB_SIZE>>>(Nloop_1,
                                                d_loop1_ix, d_loop1_iy, d_loop1_iz,
		                                d_x, d_xb,
		                                d_gbx, d_gby, d_gbz,
		                                d_Gbx, d_Gby, d_Gbz,
		                                d_pbLx_arr, d_pbLy_arr, d_pbLz_arr,
		                                d_pbRx_arr, d_pbRy_arr, d_pbRz_arr,
		                                d_flLx_arr, d_flLy_arr, d_flLz_arr,
		                                d_flRx_arr, d_flRy_arr, d_flRz_arr,
						d_p_arr);
  cudaEventRecord(stop);
  err = cudaPeekAtLastError();
  // printf("ERROR-Kernel (error code %s)!\n", cudaGetErrorString(err));

  // synchronize
  cudaDeviceSynchronize();
  
  // timing information
  cudaEventSynchronize(stop);
  float tms = 0.;
  cudaEventElapsedTime(&tms, start,stop);
  printf("gpu calc_interp time: %0.2f \n",tms);
 
#ifdef CPUKO 
  ldouble* f_tmp;
  long long NfluxX = (SX+1)*(SY)*(SZ)*NV;
  long long NfluxY = (SX)*(SY+1)*(SZ)*NV;
  long long NfluxZ = (SX)*(SY)*(SZ+1)*NV;
  
  if((f_tmp=(ldouble*)malloc(NfluxZ*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  err = cudaMemcpy(f_tmp, d_flLz_arr, NfluxZ*sizeof(ldouble), cudaMemcpyDeviceToHost);

  
  if(err != cudaSuccess) printf("failed cudaMemcpy of d_pbLx_arr to pbLx_tmp\n");
  printf("gpu calc_interp flLz[NV]: ");
  for(int iv=0;iv<NV;iv++)
    printf("%e ", get_ub(f_tmp, iv, ixTEST, iyTEST, izTEST,2));
  printf("\n");
  free(f_tmp);
#endif

  return (ldouble)tms;
}

ldouble calc_fluxes_gpu()
{
  cudaError_t err = cudaSuccess;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // Determine number of threadblocks
  int threadblocks = (Nloop_1 / TB_SIZE) + ((Nloop_1 % TB_SIZE)? 1:0);
  //printf("\nTest %d\n", threadblocks); fflush(stdout);

  // Launch kernel
  cudaEventRecord(start);
  calc_fluxes_kernel<<<threadblocks, TB_SIZE>>>(Nloop_1,
                                                d_loop1_ix, d_loop1_iy, d_loop1_iz,
		                                d_x, d_xb,
		                                d_gbx, d_gby, d_gbz,
		                                d_Gbx, d_Gby, d_Gbz,
		                                d_pbLx_arr, d_pbLy_arr, d_pbLz_arr,
		                                d_pbRx_arr, d_pbRy_arr, d_pbRz_arr,
		                                d_flLx_arr, d_flLy_arr, d_flLz_arr,
		                                d_flRx_arr, d_flRy_arr, d_flRz_arr,
		                                d_ahdxl_arr, d_ahdyl_arr, d_ahdzl_arr,
		                                d_ahdxr_arr, d_ahdyr_arr, d_ahdzr_arr,
                                                d_ahdx_arr,  d_ahdy_arr,  d_ahdz_arr,
		                                d_aradxl_arr, d_aradyl_arr, d_aradzl_arr,
		                                d_aradxr_arr, d_aradyr_arr, d_aradzr_arr,
		                                d_aradx_arr,  d_arady_arr,  d_aradz_arr,
		                                d_flbx_arr, d_flby_arr, d_flbz_arr);
  
  cudaEventRecord(stop);
  err = cudaPeekAtLastError();
  // printf("ERROR-Kernel (error code %s)!\n", cudaGetErrorString(err));

  // synchronize
  cudaDeviceSynchronize();
  
  // timing information
  cudaEventSynchronize(stop);
  float tms = 0.;
  cudaEventElapsedTime(&tms, start,stop);
  printf("gpu calc_fluxes time: %0.2f \n",tms);
 
#ifdef CPUKO 
  ldouble* flbx_tmp;
  long long NfluxX = (SX+1)*(SY)*(SZ)*NV;
  if((flbx_tmp=(ldouble*)malloc(NfluxX*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  err = cudaMemcpy(flbx_tmp, d_flbx_arr, NfluxX*sizeof(ldouble), cudaMemcpyDeviceToHost);
  if(err != cudaSuccess) printf("failed cudaMemcpy of d_flbx_arr to flbx_tmp\n");
  printf("gpu calc_fluxes flbx[NV]: ");
  for(int iv=0;iv<NV;iv++)
    printf("%e ", get_ub(flbx_tmp, iv, ixTEST, iyTEST, izTEST,0));
  printf("\n");
  free(flbx_tmp);
#endif


  return (ldouble)tms;
}

ldouble calc_update_gpu(ldouble dtin)
{
  cudaError_t err = cudaSuccess;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  // Determine number of threadblocks
  int threadblocks = (Nloop_0 / TB_SIZE) + ((Nloop_0 % TB_SIZE)? 1:0);
  //printf("\nTest %d\n", threadblocks); fflush(stdout);

  // Launch kernel
  cudaEventRecord(start);
  calc_update_kernel<<<threadblocks, TB_SIZE>>>(Nloop_0, 
						d_loop0_ix, d_loop0_iy, d_loop0_iz,
						d_x, d_xb,d_gcov, d_gcon, d_Kris,
						d_flbx_arr, d_flby_arr, d_flbz_arr,
						d_u_arr, d_p_arr, dtin);
  
  cudaEventRecord(stop);
  err = cudaPeekAtLastError();
  // printf("ERROR-Kernel (error code %s)!\n", cudaGetErrorString(err));

  // synchronize
  cudaDeviceSynchronize();
  
  // timing information
  cudaEventSynchronize(stop);
  float tms = 0.;
  cudaEventElapsedTime(&tms, start,stop);
  printf("gpu update time: %0.2f \n",tms);

  // output for comparison
#ifdef CPUKO 
  ldouble* u_tmp;
  long long Nprim  = (SX)*(SY)*(SZ)*NV;
  if((u_tmp=(ldouble*)malloc(Nprim*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  err = cudaMemcpy(u_tmp, d_u_arr, Nprim*sizeof(ldouble), cudaMemcpyDeviceToHost);
  if(err != cudaSuccess) printf("failed cudaMemcpy of d_p_arr to p_tmp\n");
  printf("gpu update uu[NV]: ");
  for(int iv=0;iv<NV;iv++)
    printf("%e ", get_u(u_tmp, iv, ixTEST, iyTEST, izTEST));
  printf("\n");
  free(u_tmp);
#endif

  // set global timestep dt
  dt = dtin;

  return (ldouble)tms;
}
