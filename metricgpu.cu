extern "C" {

#include "ko.h"

}

#include "kogpu.h"


// persistent arrays
int *d_loop0_ix, *d_loop0_iy, *d_loop0_iz;
ldouble *d_x;    //[(NX+NY+NZ+6*NG)*sizeof(ldouble)]
ldouble *d_xb;   //[(NX+1+NY+1+NZ+1+6*NG)*sizeof(ldouble)]
ldouble *d_gcov; //[SX*SY*SZMET*sizeof(ldouble)]
ldouble *d_gcon; //[SX*SY*SZMET*sizeof(ldouble)]
ldouble *d_Kris; //[(SX)*(SY)*(SZMET)*64*sizeof(ldouble)];

int GEOMETRY_HAS_BEEN_ALLOCED = 0;

int push_geometry_gpu()
{
 
  if (GEOMETRY_HAS_BEEN_ALLOCED == 1) {
    fprintf(stderr, "geom already alloced\n");
    return 0;
  }
  fprintf(stderr, "geom being alloced\n");

  cudaError_t err = cudaSuccess;

  // array sizes
  long long Nx     = (NX+NY+NZ+6*NG);
  long long Nxb    = (NX+1+NY+1+NZ+1+6*NG);
  long long Nmet   = (SX)*(SY)*(SZMET)*gSIZE;
  long long Nkris  = (SX)*(SY)*(SZMET)*64;

  // allocate device arrays
  err = cudaMalloc(&d_loop0_ix, sizeof(int)*Nloop_0);
  err = cudaMalloc(&d_loop0_iy, sizeof(int)*Nloop_0);
  err = cudaMalloc(&d_loop0_iz, sizeof(int)*Nloop_0);

  err = cudaMalloc(&d_x,        sizeof(ldouble)*Nx);
  err = cudaMalloc(&d_xb,       sizeof(ldouble)*Nxb);

  err = cudaMalloc(&d_gcov,     sizeof(ldouble)*Nmet);
  err = cudaMalloc(&d_gcon,     sizeof(ldouble)*Nmet);
  err = cudaMalloc(&d_Kris,     sizeof(ldouble)*Nkris);

  // Make 1D arrays of ix,iy,iz indicies for easier copy to device
  int *h_loop0_ix = (int*)malloc(sizeof(int)*Nloop_0);
  int *h_loop0_iy = (int*)malloc(sizeof(int)*Nloop_0);
  int *h_loop0_iz = (int*)malloc(sizeof(int)*Nloop_0);

  for(int ii=0; ii<Nloop_0; ii++){
    h_loop0_ix[ii] = loop_0[ii][0];
    h_loop0_iy[ii] = loop_0[ii][1];
    h_loop0_iz[ii] = loop_0[ii][2];
  }

  // copy index arrays
  // this reads from loop_0 above, which is set in finite.c:alloc_loops()
  err =  cudaMemcpy(d_loop0_ix, h_loop0_ix, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop0_ix to device failed.\n");
  err =  cudaMemcpy(d_loop0_iy, h_loop0_iy, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop0_iy to device failed.\n");
  err =  cudaMemcpy(d_loop0_iz, h_loop0_iz, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop0_iz to device failed.\n");

  free(h_loop0_ix);
  free(h_loop0_iy);
  free(h_loop0_iz);

  // copy coordinate arrays, which are set in finite.c:set_grid(...) & finite.c:set_x[b](...)
  err =  cudaMemcpy(d_x, x, sizeof(ldouble)*Nx, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing x to device failed.\n");
  err =  cudaMemcpy(d_xb, xb, sizeof(ldouble)*Nxb, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing xb to device failed.\n");

  // copy metric/Christoffel arrays, which are set in metric.c:calc_metric()
  err = cudaMemcpy(d_gcov, g, sizeof(double)*Nmet, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing g to device failed.\n");
  err = cudaMemcpy(d_gcon, G, sizeof(double)*Nmet, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing G to device failed.\n");
  err = cudaMemcpy(d_Kris, gKr, sizeof(double)*Nkris, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing gKr to device failed.\n");

  GEOMETRY_HAS_BEEN_ALLOCED = 1;

  return 0;
}


int free_geometry_gpu()
{
  if (GEOMETRY_HAS_BEEN_ALLOCED == 0) {
    fprintf(stderr, "can't free what we don't have\n");
    return 0;
  }
  fprintf(stderr, "freeing geometry\n");

  cudaFree(d_loop0_ix);
  cudaFree(d_loop0_iy);
  cudaFree(d_loop0_iz);

  cudaFree(d_x);
  cudaFree(d_xb);

  cudaFree(d_gcov);
  cudaFree(d_gcon);
  cudaFree(d_Kris);

  GEOMETRY_HAS_BEEN_ALLOCED = 0;

  return 0;
}


