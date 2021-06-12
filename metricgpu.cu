extern "C" {

#include "ko.h"

}

#include "kogpu.h"

// persistent arrays
int *d_loop0_ix, *d_loop0_iy, *d_loop0_iz;
int *d_loop1_ix, *d_loop1_iy, *d_loop1_iz;
int *d_loop4_ix, *d_loop4_iy, *d_loop4_iz;
ldouble *d_x;    //[(NX+NY+NZ+6*NG)*sizeof(ldouble)]
ldouble *d_xb;   //[(NX+1+NY+1+NZ+1+6*NG)*sizeof(ldouble)]
ldouble *d_gcov; //[SX*SY*SZMET*sizeof(ldouble)]
ldouble *d_gcon; //[SX*SY*SZMET*sizeof(ldouble)]
ldouble *d_gbx, *d_gby, *d_gbz;
ldouble *d_Gbx, *d_Gby, *d_Gbz;
ldouble *d_Kris; //[(SX)*(SY)*(SZMET)*64*sizeof(ldouble)];

int GEOMETRY_HAS_BEEN_ALLOCED = 0;


//**********************************************************************
//fills geometry structure for cell ix,iy,iz
//**********************************************************************

__device__ __host__ int fill_geometry_device(int ix,int iy,int iz, void* geom,
					     ldouble* x_arr,ldouble* g_arr, ldouble* G_arr)
{

  struct geometry *ggg 
    = (struct geometry *) geom;

  ggg->par=-1;
  ggg->ifacedim = -1;
  ggg->coords=MYCOORDS;

  // coordinates
  ggg->ix=ix;  ggg->iy=iy;  ggg->iz=iz;
  ggg->xxvec[0]=0.;
  ggg->xxvec[1]=get_x_device(x_arr, ix, 0);
  ggg->xxvec[2]=get_x_device(x_arr, iy, 1);
  ggg->xxvec[3]=get_x_device(x_arr, iz, 2);
  ggg->xx=ggg->xxvec[1];
  ggg->yy=ggg->xxvec[2];
  ggg->zz=ggg->xxvec[3];

  // geometry
  //pick_g(ix,iy,iz,ggg->gg);
  //pick_G(ix,iy,iz,ggg->GG);
  for(int i=0;i<4;i++)
  {
    for(int j=0;j<5;j++)
    {
      ggg->gg[i][j]=get_g(g_arr,i,j,ix,iy,iz);
      ggg->GG[i][j]=get_g(G_arr,i,j,ix,iy,iz);
    }
  }
  
  ggg->alpha=sqrt(-1./ggg->GG[0][0]);
  ggg->gdet=ggg->gg[3][4];
  ggg->gttpert=ggg->GG[3][4];
  
  return 0;  
}


//**********************************************************************
//fills geometry structure for cell face ix,iy,iz in idim
//**********************************************************************

__host__ __device__ int fill_geometry_face_device(int ix,int iy,int iz,int idim, void *geom,
						  ldouble* x_arr, ldouble* xb_arr,
						  ldouble* gbx_arr,ldouble* gby_arr,ldouble* gbz_arr,
						  ldouble* Gbx_arr,ldouble* Gby_arr,ldouble* Gbz_arr)
{
  //TODO
  /*
  if(doingpostproc) //not precalculated
    {
      fill_geometry_face_arb_device(ix,iy,iz,idim,geom,MYCOORDS);
      return 0;
    }
  */
  
  struct geometry *ggg 
    = (struct geometry *) geom;

  ggg->par=-1;
  ggg->ifacedim = idim;
  ggg->coords=MYCOORDS;

  // coordinates
  ggg->ix=ix;  ggg->iy=iy;  ggg->iz=iz;  
  ggg->xxvec[0]=0.;
  if(idim==0) //x-face
  {
    ggg->xxvec[1]=get_xb_device(xb_arr,ix,0);
    ggg->xxvec[2]=get_x_device (x_arr, iy,1);
    ggg->xxvec[3]=get_x_device (x_arr, iz,2);
  }
  if(idim==1) //y-face
  {
    ggg->xxvec[1]=get_x_device (x_arr, ix,0);
    ggg->xxvec[2]=get_xb_device(xb_arr,iy,1);
    ggg->xxvec[3]=get_x_device (x_arr, iz,2);
  }
  if(idim==2) //z-face
  {
    ggg->xxvec[1]=get_x_device (x_arr, ix,0);
    ggg->xxvec[2]=get_x_device (x_arr, iy,1);
    ggg->xxvec[3]=get_xb_device(xb_arr,iz,2);
  }
  ggg->xx=ggg->xxvec[1];
  ggg->yy=ggg->xxvec[2];
  ggg->zz=ggg->xxvec[3];

  // geometry
  //pick_gb(ix,iy,iz,idim,ggg->gg);
  //pick_Gb(ix,iy,iz,idim,ggg->GG);
  for(int i=0;i<4;i++)
  {
    for(int j=0;j<5;j++)
    {
      if(idim==0)
      {
        ggg->gg[i][j]=get_gb(gbx_arr,i,j,ix,iy,iz,0);
	ggg->GG[i][j]=get_gb(Gbx_arr,i,j,ix,iy,iz,0);
      }
      else if(idim==1)
      {
        ggg->gg[i][j]=get_gb(gby_arr,i,j,ix,iy,iz,1);
	ggg->GG[i][j]=get_gb(Gby_arr,i,j,ix,iy,iz,1);
      }
      else if(idim==2)
      {
        ggg->gg[i][j]=get_gb(gbz_arr,i,j,ix,iy,iz,2);
	ggg->GG[i][j]=get_gb(Gbz_arr,i,j,ix,iy,iz,2);
      }      
    }
  }

  ggg->alpha=sqrt(-1./ggg->GG[0][0]);
  ggg->gdet=ggg->gg[3][4];
  ggg->gttpert=ggg->GG[3][4];

  return 0;
}


//**********************************************************************
//pushes geometry data to the gpu
//**********************************************************************

int push_geometry_gpu()
{
  if (GEOMETRY_HAS_BEEN_ALLOCED == 1) {
    return 0;
  }

  cudaError_t err = cudaSuccess;

  // array sizes
  long long Nx     = (NX+NY+NZ+6*NG);
  long long Nxb    = (NX+1+NY+1+NZ+1+6*NG);
  long long Nmet   = (SX)*(SY)*(SZMET)*gSIZE;
  long long Nkris  = (SX)*(SY)*(SZMET)*64;

  long long NmetX = (SX+1)*(SY)*(SZMET)*gSIZE;
  long long NmetY = (SX)*(SY+1)*(SZMET)*gSIZE;
  long long NmetZ = (SX)*(SY)*(SZMET+1)*gSIZE;
	       
  // allocate device arrays
  err = cudaMalloc(&d_loop0_ix, sizeof(int)*Nloop_0);
  err = cudaMalloc(&d_loop0_iy, sizeof(int)*Nloop_0);
  err = cudaMalloc(&d_loop0_iz, sizeof(int)*Nloop_0);
  err = cudaMalloc(&d_loop1_ix, sizeof(int)*Nloop_1);
  err = cudaMalloc(&d_loop1_iy, sizeof(int)*Nloop_1);
  err = cudaMalloc(&d_loop1_iz, sizeof(int)*Nloop_1);
  err = cudaMalloc(&d_loop4_ix, sizeof(int)*Nloop_4);
  err = cudaMalloc(&d_loop4_iy, sizeof(int)*Nloop_4);
  err = cudaMalloc(&d_loop4_iz, sizeof(int)*Nloop_4);

  err = cudaMalloc(&d_x,        sizeof(ldouble)*Nx);
  err = cudaMalloc(&d_xb,       sizeof(ldouble)*Nxb);

  err = cudaMalloc(&d_gcov,     sizeof(ldouble)*Nmet);
  err = cudaMalloc(&d_gcon,     sizeof(ldouble)*Nmet);
  err = cudaMalloc(&d_Kris,     sizeof(ldouble)*Nkris);

  err = cudaMalloc(&d_gbx,     sizeof(ldouble)*NmetX);
  err = cudaMalloc(&d_Gby,     sizeof(ldouble)*NmetX);
  err = cudaMalloc(&d_gbz,     sizeof(ldouble)*NmetY);
  err = cudaMalloc(&d_Gbx,     sizeof(ldouble)*NmetY);  
  err = cudaMalloc(&d_gby,     sizeof(ldouble)*NmetZ);
  err = cudaMalloc(&d_Gbz,     sizeof(ldouble)*NmetZ);
  
  // Make 1D arrays of ix,iy,iz indicies for easier copy to device
  int *h_loop0_ix = (int*)malloc(sizeof(int)*Nloop_0);
  int *h_loop0_iy = (int*)malloc(sizeof(int)*Nloop_0);
  int *h_loop0_iz = (int*)malloc(sizeof(int)*Nloop_0);
  for(int ii=0; ii<Nloop_0; ii++)
  {
    h_loop0_ix[ii] = loop_0[ii][0];
    h_loop0_iy[ii] = loop_0[ii][1];
    h_loop0_iz[ii] = loop_0[ii][2];
  }

  int *h_loop1_ix = (int*)malloc(sizeof(int)*Nloop_1);
  int *h_loop1_iy = (int*)malloc(sizeof(int)*Nloop_1);
  int *h_loop1_iz = (int*)malloc(sizeof(int)*Nloop_1);
  for(int ii=0; ii<Nloop_1; ii++)
  {
    h_loop1_ix[ii] = loop_1[ii][0];
    h_loop1_iy[ii] = loop_1[ii][1];
    h_loop1_iz[ii] = loop_1[ii][2];
  }

  int *h_loop4_ix = (int*)malloc(sizeof(int)*Nloop_4);
  int *h_loop4_iy = (int*)malloc(sizeof(int)*Nloop_4);
  int *h_loop4_iz = (int*)malloc(sizeof(int)*Nloop_4);
  for(int ii=0; ii<Nloop_4; ii++){
    h_loop4_ix[ii] = loop_4[ii][0];
    h_loop4_iy[ii] = loop_4[ii][1];
    h_loop4_iz[ii] = loop_4[ii][2];
  }

  // copy index arrays
  // this reads from loop_0 above, which is set in finite.c:alloc_loops()
  err =  cudaMemcpy(d_loop0_ix, h_loop0_ix, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop0_ix to device failed.\n");
  err =  cudaMemcpy(d_loop0_iy, h_loop0_iy, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop0_iy to device failed.\n");
  err =  cudaMemcpy(d_loop0_iz, h_loop0_iz, sizeof(int)*Nloop_0, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop0_iz to device failed.\n");

  err =  cudaMemcpy(d_loop1_ix, h_loop1_ix, sizeof(int)*Nloop_1, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop1_ix to device failed.\n");
  err =  cudaMemcpy(d_loop1_iy, h_loop1_iy, sizeof(int)*Nloop_1, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop1_iy to device failed.\n");
  err =  cudaMemcpy(d_loop1_iz, h_loop1_iz, sizeof(int)*Nloop_1, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop1_iz to device failed.\n");
  
  err =  cudaMemcpy(d_loop4_ix, h_loop4_ix, sizeof(int)*Nloop_4, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop4_ix to device failed.\n");
  err =  cudaMemcpy(d_loop4_iy, h_loop4_iy, sizeof(int)*Nloop_4, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop4_iy to device failed.\n");
  err =  cudaMemcpy(d_loop4_iz, h_loop4_iz, sizeof(int)*Nloop_4, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing d_loop4_iz to device failed.\n");

  free(h_loop0_ix);
  free(h_loop0_iy);
  free(h_loop0_iz);
  
  free(h_loop1_ix);
  free(h_loop1_iy);
  free(h_loop1_iz);
  
  free(h_loop4_ix);
  free(h_loop4_iy);
  free(h_loop4_iz);

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

  // copy boundary metric arrays
  err = cudaMemcpy(d_gbx, gbx, sizeof(double)*NmetX, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing gbx to device failed.\n");
  err = cudaMemcpy(d_gby, gby, sizeof(double)*NmetY, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing gby to device failed.\n");
  err = cudaMemcpy(d_gbz, gbz, sizeof(double)*NmetZ, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing gbz to device failed.\n");
  
  err = cudaMemcpy(d_Gbx, Gbx, sizeof(double)*NmetX, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing Gbx to device failed.\n");
  err = cudaMemcpy(d_Gby, Gby, sizeof(double)*NmetY, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing Gby to device failed.\n");
  err = cudaMemcpy(d_Gbz, Gbz, sizeof(double)*NmetZ, cudaMemcpyHostToDevice);
  if(err != cudaSuccess) printf("Passing Gbz to device failed.\n");

  // indicate we allocated geometry
  GEOMETRY_HAS_BEEN_ALLOCED = 1;

  return 0;
}


//**********************************************************************
//frees gpu geometry data 
//**********************************************************************

int free_geometry_gpu()
{
  if (GEOMETRY_HAS_BEEN_ALLOCED == 0) {
    return 0;
  }

  cudaFree(d_loop0_ix);
  cudaFree(d_loop0_iy);
  cudaFree(d_loop0_iz);
  cudaFree(d_loop1_ix);
  cudaFree(d_loop1_iy);
  cudaFree(d_loop1_iz);
  cudaFree(d_loop4_ix);
  cudaFree(d_loop4_iy);
  cudaFree(d_loop4_iz);

  cudaFree(d_x);
  cudaFree(d_xb);

  cudaFree(d_gcov);
  cudaFree(d_gcon);
  cudaFree(d_Kris);

  cudaFree(d_gbx);
  cudaFree(d_gby);
  cudaFree(d_gbz);

  cudaFree(d_Gbx);
  cudaFree(d_Gby);
  cudaFree(d_Gbz);
  
  GEOMETRY_HAS_BEEN_ALLOCED = 0;

  return 0;
}


