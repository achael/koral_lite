extern "C" {

#include "ko.h"

}

#define ixTEST 13
#define iyTEST 21
#define izTEST 8
#define iiTEST 22222
#define ivTEST 0

// get data value from array u_arr of the quantity indexed iv
// at the cell center indexed ix,iy,iz
// copied from get_u macro in ko.h
// iX(ix), iY(iy), iZ(iz) are macros defined in ko.h, that return either the index or 0
// depending on the problem dimension
__device__ ldouble get_u_device(ldouble* u_arr,int iv,int ix,int iy,int iz)
{
  ldouble u_out;
  u_out = u_arr[iv + (iX(ix)+(NGCX))*NV + \
		     (iY(iy)+(NGCY))*(SX)*NV + \
		     (iZ(iz)+(NGCZ))*(SY)*(SX)*NV];
  return u_out;
}

// get data value from array ub_arr of quantity indexed iv
// on the left wall of cell indexed ix,iy,iz in dimension idim
// copied from get_ub macro in ko.h
__device__ ldouble get_ub_device(ldouble* ub_arr, int iv, int ix, int iy, int iz, int idim)
{
  ldouble ub_out;
  ub_out = (idim==0 ? ub_arr[iv + (iX(ix)+(NGCX))*NV + \
				  (iY(iy)+(NGCY))*(SX+1)*NV + \
				  (iZ(iz)+(NGCZ))*(SY)*(SX+1)*NV] : \
	   (idim==1 ? ub_arr[iv + (iX(ix)+(NGCX))*NV + \
			          (iY(iy)+(NGCY))*(SX)*NV + \
		                  (iZ(iz)+(NGCZ))*(SY+1)*(SX)*NV] : \
	   (idim==2 ? ub_arr[iv + (iX(ix)+(NGCX))*NV + \
			          (iY(iy)+(NGCY))*(SX)*NV + \
			          (iZ(iz)+(NGCZ))*(SY)*(SX)*NV] : 0.)));
  return ub_out;
}

// get grid coordinate on the cell wall indexed ic in dimension idim
// copied from get_xb macro in ko.h
__device__ ldouble get_xb_device(ldouble* xb_arr, int ic, int idim)
{
  ldouble xb_out;
  xb_out = (idim==0 ? xb_arr[ic+NG] :		     \
           (idim==1 ? xb_arr[ic+NG + NX+2*NG + 1] :  \
	   (idim==2 ? xb_arr[ic+NG + NX+2*NG +1 + NY+2*NG +1 ] : 0.)));

  return xb_out;
}

// get size of cell indexed ic in dimension idim
// copied from get_size_x in finite.c
__device__ ldouble get_size_x_device(ldouble* xb_arr, int ic, int idim)
{
  ldouble dx;
  dx = get_xb_device(xb_arr, ic+1,idim) - get_xb_device(xb_arr, ic, idim);
  return dx;
}


__global__ void calc_update_gpu_kernel(ldouble dtin, int Nloop_0, int* d_array, 
                                       int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
				       ldouble* xb_arr,
				       ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr,
				       ldouble* u_arr)
{

  int ii;
  int ix,iy,iz,iv;
  ldouble dx,dy,dz;
  ldouble flxl,flxr,flyl,flyr,flzl,flzr;
  ldouble val,du;
  ldouble ms[NV];
  //ldouble gs[NV]; //NOTE gs[NV] is for artifical sources, rarely used

  // get index for this thread
  // Nloop_0 is number of cells to update;
  // usually Nloop_0=NX*NY*NZ, but sometimes there are weird bcs inside domain 
  ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_0) return;
  
  atomicAdd(d_array, 1); // NOTE TODO: placeholder test
  
  // get indices from 1D arrays
  ix=loop_0_ix[ii];
  iy=loop_0_iy[ii];
  iz=loop_0_iz[ii]; 

  if(ii==iiTEST){
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
  dx = get_size_x_device(xb_arr,ix,0); //dx=get_size_x(ix,0);
  dy = get_size_x_device(xb_arr,iy,1); //dy=get_size_x(iy,1);
  dz = get_size_x_device(xb_arr,iz,2); //dz=get_size_x(iz,2);

  // test sizes 
  if(ii==0)
  {
    printf("D size_x 0 %e \n", get_size_x_device(xb_arr,ixTEST,0));
    printf("D size_x 1 %e \n", get_size_x_device(xb_arr,iyTEST,1));
    printf("D size_x 2 %e \n", get_size_x_device(xb_arr,izTEST,2));
  }
  
  //update all conserved according to fluxes and source terms      
  for(iv=0;iv<NV;iv++)
  {	

    // Get the initial value of the conserved quantity
    val = get_u_device(u_arr,iv,ix,iy,iz);
    if(ix==ixTEST && iy==iyTEST && iz==izTEST && iv==ivTEST)
      printf("D u: %e\n", val);
    
    // Get the fluxes on the six faces.
    // flbx, flby, flbz are the fluxes at the LEFT walls of cell ix, iy, iz.
    // To get the RIGHT fluxes, we need flbx(ix+1,iy,iz), etc.
    flxl=get_ub_device(flbx_arr,iv,ix,iy,iz,0);
    flxr=get_ub_device(flbx_arr,iv,ix+1,iy,iz,0);
    flyl=get_ub_device(flby_arr,iv,ix,iy,iz,1);
    flyr=get_ub_device(flby_arr,iv,ix,iy+1,iz,1);
    flzl=get_ub_device(flbz_arr,iv,ix,iy,iz,2);
    flzr=get_ub_device(flbz_arr,iv,ix,iy,iz+1,2);
	   
    if(ix==ixTEST && iy==iyTEST && iz==izTEST && iv==ivTEST)
      printf("D fluxes: %e %e %e %e %e %e\n", flxl,flxr,flyl,flyr,flzl,flzr);

    // Compute Delta U from the six fluxes
    du = -(flxr-flxl)*dtin/dx - (flyr-flyl)*dtin/dy - (flzr-flzl)*dtin/dz;

    // Compute the new conserved by adding Delta U and the source term
    val += (du + ms[iv]*dtin);

    // Save the new conserved to memory
    
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
	
    u_arr[iv] = val;
    //set_u(u,iv,ix,iy,iz,val);	 

  }  
}

int calc_update_gpu(ldouble dtin)
{

  int TB_SIZE = 64;   
  int *d_temp, h_temp=0;
  int *d_loop0_ix,*d_loop0_iy,*d_loop0_iz;
  int *h_loop0_ix,*h_loop0_iy,*h_loop0_iz;
  ldouble *d_xb_arr;
  ldouble *d_u_arr;
  ldouble *d_flbx_arr,*d_flby_arr,*d_flbz_arr;
  
  cudaError_t err = cudaSuccess;

  // Allocate device arrays 
  
  // printf("ERROR (error code %s)!\n", cudaGetErrorString(err));
  err = cudaMalloc(&d_temp, sizeof(int));

  err = cudaMalloc(&d_loop0_ix, sizeof(int)*Nloop_0);
  err = cudaMalloc(&d_loop0_iy, sizeof(int)*Nloop_0);
  err = cudaMalloc(&d_loop0_iz, sizeof(int)*Nloop_0);

  // NOTE: size of xb,flbx,flby,flbz is copied from initial malloc in misc.c
  // these need to be long long if the grid is on one tile and large (~256^3)
  long long Nxb    = (NX+1+NY+1+NZ+1+6*NG);
  long long Nprim  = (SX)*(SY)*(SZ)*NV;
  long long NfluxX = (SX+1)*(SY)*(SZ)*NV;
  long long NfluxY = (SX)*(SY+1)*(SZ)*NV;
  long long NfluxZ = (SX)*(SY)*(SZ+1)*NV;
  
  err = cudaMalloc(&d_xb_arr,   sizeof(ldouble)*Nxb);
  err = cudaMalloc(&d_u_arr,    sizeof(ldouble)*Nprim);
  err = cudaMalloc(&d_flbx_arr, sizeof(ldouble)*NfluxX);
  err = cudaMalloc(&d_flby_arr, sizeof(ldouble)*NfluxY);
  err = cudaMalloc(&d_flbz_arr, sizeof(ldouble)*NfluxZ);
  
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
    if (ii==iiTEST) printf("H   :  %d %d %d %d\n",ii,h_loop0_ix[ii],h_loop0_iy[ii],h_loop0_iz[ii]) ;
  }

  err =  cudaMemcpy(d_loop0_ix, h_loop0_ix, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_loop0_iy, h_loop0_iy, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_loop0_iz, h_loop0_iz, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);

  free(h_loop0_ix);
  free(h_loop0_iy);
  free(h_loop0_iz);

  // copy grid boundary data from xb (global array) to device
  printf("H size_x 0 %e \n", get_size_x(ixTEST,0));
  printf("H size_x 1 %e \n", get_size_x(iyTEST,1));
  printf("H size_x 2 %e \n", get_size_x(izTEST,2));
  err =  cudaMemcpy(d_xb_arr, xb, sizeof(ldouble)*Nxb, cudaMemcpyHostToDevice);

  // copy conserved quantities from u (global array) to device
  printf("H u: %e \n", get_u(u,ivTEST,ixTEST,iyTEST,izTEST));
  err = cudaMemcpy(d_u_arr, u, sizeof(ldouble)*Nprim, cudaMemcpyHostToDevice);
  
  // copy fluxes data from flbx,flby,flbz (global arrays) to device
  printf("H fluxes: %e %e %e %e %e %e\n",
	 get_ub(flbx,ivTEST,ixTEST,iyTEST,izTEST,0),
	 get_ub(flbx,ivTEST,ixTEST+1,iyTEST,izTEST,0),
         get_ub(flby,ivTEST,ixTEST,iyTEST,izTEST,1),
	 get_ub(flby,ivTEST,ixTEST,iyTEST+1,izTEST,1),
	 get_ub(flbz,ivTEST,ixTEST,iyTEST,izTEST,2),
	 get_ub(flbz,ivTEST,ixTEST,iyTEST,izTEST+1,2));
  err =  cudaMemcpy(d_flbx_arr, flbx, sizeof(ldouble)*NfluxX, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_flby_arr, flby, sizeof(ldouble)*NfluxY, cudaMemcpyHostToDevice);
  err =  cudaMemcpy(d_flbz_arr, flbz, sizeof(ldouble)*NfluxZ, cudaMemcpyHostToDevice);

  // Launch calc_update_gpu_kernel

  int threadblocks = (Nloop_0 / TB_SIZE) + ((Nloop_0 % TB_SIZE)? 1:0);
  printf("\nTest %d\n", threadblocks); fflush(stdout);

  err = cudaMemset(d_temp,0,sizeof(int));
  // printf("ERRORMEMESET (error code %s)!\n", cudaGetErrorString(err));

  calc_update_gpu_kernel<<<threadblocks, TB_SIZE>>>(dtin, Nloop_0, d_temp,
						    d_loop0_ix, d_loop0_iy, d_loop0_iz,
						    d_xb_arr,
						    d_flbx_arr, d_flby_arr, d_flbz_arr,
						    d_u_arr);
  err = cudaPeekAtLastError();
  cudaDeviceSynchronize();
  // printf("ERROR-Kernel (error code %s)!\n", cudaGetErrorString(err));

  err =  cudaMemcpy(&h_temp, d_temp, sizeof(int), cudaMemcpyDeviceToHost);
  // printf("ERROR-Memcpy (error code %s)!\n", cudaGetErrorString(err));

  printf("back from device %d\n\n",h_temp);

  // TODO Copy updated u back from device?
  //err = cudaMemcpy(&u, d_u_arr, sizeof(ldouble)*Nprim, cudaMemcpyDeviceToHost);
  
  // Free Device Memory
  cudaFree(d_loop0_ix);
  cudaFree(d_loop0_iy);
  cudaFree(d_loop0_iz);
  
  cudaFree(d_xb_arr);
  cudaFree(d_flbx_arr);
  cudaFree(d_flby_arr);
  cudaFree(d_flbz_arr);
  cudaFree(d_u_arr);
  
  cudaFree(d_temp);

  // set global timestep dt
  dt = dtin;

  return 0;
}
