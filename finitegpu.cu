extern "C" {

#include "ko.h"

}

// get grid location on boundary indexed ic in dimension idim
// copied from get_xb macro in ko.h
__device__ ldouble get_xb_gpu_kernel(ldouble* xb_arr, int ic, int idim)
{
  ldouble xb_out;
  xb_out = (idim==0 ? xb_arr[ic+NG] :		     \
           (idim==1 ? xb_arr[ic+NG + NX+2*NG + 1] :  \
	   (idim==2 ? xb_arr[ic+NG + NX+2*NG +1 + NY+2*NG +1 ] : 0.)));

  return xb_out;
}

// get size of cell indexed ic in dimension idim
// copied from get_size_x in finite.c
__device__ ldouble get_size_x_gpu_kernel(ldouble* xb_arr, int ic, int idim)
{
  ldouble dx;
  dx = get_xb_gpu_kernel(xb_arr, ic+1,idim) - get_xb_gpu_kernel(xb_arr, ic, idim);
  return dx;
}


__global__ void calc_update_gpu_kernel(ldouble dtin, int Nloop_0, int* d_array, 
                                       int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
				       ldouble* xb_arr)
{

  int ii;
  int ix,iy,iz,iv;
  ldouble dx,dy,dz;
  ldouble val,du;
  ldouble ms[NV];
  //ldouble gs[NV]; //NOTE gs[NV] is for artifical sources, rarely used

  // get index for this thread
  ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_0) return;
  
  atomicAdd(d_array, 1); // NOTE TODO: placeholder test
  
  // get indices from 1D arrays
  ix=loop_0_ix[ii];
  iy=loop_0_iy[ii];
  iz=loop_0_iz[ii]; 

  if(ii==22222){
    printf("D   : %d %d %d %d\n",ii, ix,iy,iz);
  }

  // Source term
  // check if cell is active
  // NOTE: is_cell_active always returns 1 -- a placeholder function put in long ago
  
  if(0) //if(is_cell_active(ix,iy,iz)==0)
  {
    // Source terms applied only for active cells	  
    for(iv=0;iv<NV;iv++) ms[iv]=0.; 
  }
  else
  {
     // Get metric source terms ms[iv]
     // and any other source terms gs[iv] 

     //f_metric_source_term(ix,iy,iz,ms);  //TODO: somewhat complicated
     //f_general_source_term(ix,iy,iz,gs); //NOTE: *very* rarely used,ignore for now
     for(iv=0;iv<NV;iv++)
     {
       ms[iv] = 0; // TODO: placeholder metric term of 0
       //ms[iv]+=gs[iv];
     }
  }

  
      
  // Get the cell size in the three directions
  dx = get_size_x_gpu_kernel(xb_arr,ix,0); //dx=get_size_x(ix,0);
  dy = get_size_x_gpu_kernel(xb_arr,iy,1); //dy=get_size_x(iy,1);
  dz = get_size_x_gpu_kernel(xb_arr,iz,2); //dz=get_size_x(iz,2);

  if(ii==0)
  {
    printf("D size_x 0 %e \n", get_size_x_gpu_kernel(xb_arr,11,0));
    printf("D size_x 1 %e \n", get_size_x_gpu_kernel(xb_arr,13,0));
    printf("D size_x 2 %e \n", get_size_x_gpu_kernel(xb_arr,5,0));
  }
  /*
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
	
//#ifdef SKIPHDEVOLUTION
//	if(iv>=NVMHD)
//#endif
//#ifdef RADIATION
//#ifdef SKIPRADEVOLUTION
//#ifdef EVOLVEPHOTONNUMBER
//	if(iv!=EE && iv!=FX && iv!=FY && iv!=FZ && iv!=NF)
//#else
//	if(iv!=EE && iv!=FX && iv!=FY && iv!=FZ)
//#endif
//#endif  // SKIPRADEVOLUTION
//#endif  // RADIATION
//#ifdef SKIPHDBUTENERGY
//	if(iv>=NVMHD || iv==UU)
//#endif
	
	
	set_u(u,iv,ix,iy,iz,val);	 

  }  // for(iv=0;iv<NV;iv++)
 */

}

int calc_update_gpu(ldouble dtin)
{

  int TB_SIZE = 64;   
  int *d_temp, h_temp=0;
  int *d_loop0_ix,*d_loop0_iy,*d_loop0_iz;
  int *h_loop0_ix,*h_loop0_iy,*h_loop0_iz;
  ldouble *d_xb_arr;
  
  cudaError_t err = cudaSuccess;

  // Allocate device arrays 
  
  // printf("ERROR (error code %s)!\n", cudaGetErrorString(err));
  err = cudaMalloc(&d_temp, sizeof(int));

  err = cudaMalloc(&d_loop0_ix, sizeof(int)*Nloop_0);
  err = cudaMalloc(&d_loop0_iy, sizeof(int)*Nloop_0);
  err = cudaMalloc(&d_loop0_iz, sizeof(int)*Nloop_0);

  err = cudaMalloc(&d_xb_arr, sizeof(ldouble)*(NX+1+NY+1+NZ+1+6*NG));
  
  // Copy data to device arrays
  
  // NOTE: when we add more functions to device, most of these should only be copied once
  // Make 1D arrays of ix,iy,iz indicies and copy to device
  h_loop0_ix = (int*)malloc(sizeof(int)*Nloop_0);
  h_loop0_iy = (int*)malloc(sizeof(int)*Nloop_0);
  h_loop0_iz = (int*)malloc(sizeof(int)*Nloop_0);

  for(int ii=0; ii<Nloop_0; ii++){
    h_loop0_ix[ii] = loop_0[ii][0];     
    h_loop0_iy[ii] = loop_0[ii][1];     
    h_loop0_iz[ii] = loop_0[ii][2];
    if (ii==22222) printf("H   :  %d %d %d %d\n",ii,h_loop0_0[ii],h_loop0_1[ii],h_loop0_2[ii]) ;
  }

  err =  cudaMemcpy(d_loop0_ix, h_loop0_ix, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_loop0_iy, h_loop0_iy, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_loop0_iz, h_loop0_iz, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);

  free(h_loop0_ix);
  free(h_loop0_iy);
  free(h_loop0_iz);

  // copy grid boundary data xb (global array) to device
  // NOTE: size of xb is copied from initial malloc in misc.c 
  printf("H size_x 0 %e \n", get_size_x(11,0));
  printf("H size_x 1 %e \n", get_size_x(13,0));
  printf("H size_x 2 %e \n", get_size_x(5,0));
  
  err =  cudaMemcpy(d_xb_arr, xb, sizeof(ldouble)*(NX+1+NY+1+NZ+1+6*NG), cudaMemcpyHostToDevice);
  
  // Launch calc_update_gpu_kernel

  int threadblocks = (Nloop_0 / TB_SIZE) + ((Nloop_0 % TB_SIZE)? 1:0);
  printf("Test %d\n", threadblocks); fflush(stdout);

  err = cudaMemset(d_temp,0,sizeof(int));
  // printf("ERRORMEMESET (error code %s)!\n", cudaGetErrorString(err));

  calc_update_gpu_kernel<<<threadblocks, TB_SIZE>>>(dtin, Nloop_0, d_temp,
						    d_loop0_ix, d_loop0_iy, d_loop0_iz);
  err = cudaPeekAtLastError();
  cudaDeviceSynchronize();
  // printf("ERROR-Kernel (error code %s)!\n", cudaGetErrorString(err));

  err =  cudaMemcpy(&h_temp, d_temp, sizeof(int), cudaMemcpyDeviceToHost);
  // printf("ERROR-Memcpy (error code %s)!\n", cudaGetErrorString(err));

  printf("back from device %d\n\n",h_temp);

  // Free Device Memory
  cudaFree(d_loop0_ix);
  cudaFree(d_loop0_iy);
  cudaFree(d_loop0_iz);
  cudaFree(d_xb_arr);
  
  cudaFree(d_temp);

  // set global timestep dt
  dt = dtin;

  return 0;
}
