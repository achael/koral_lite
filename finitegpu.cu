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

int *d_cellflag_arr, *d_int_slot_arr;

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
  long long Nprim  = (SX)*(SY)*(SZ)*NV;
  err = cudaMemcpy(d_u_arr, u, sizeof(ldouble)*Nprim, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_p_arr, p, sizeof(ldouble)*Nprim, cudaMemcpyHostToDevice);

  // copy fluxes and wavespeeds from host to device
  // TODO: in the future, this will be entirely internal to the GPU
  long long Ngrid  = (SX)*(SY)*(SZ);
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
  err = cudaMemcpy(d_ahdyl_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdzl_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_ahdxr_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdyr_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdzr_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_ahdx_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdy_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_ahdz_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_aradxl_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradyl_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradzl_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_aradxr_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradyr_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradzr_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  

  err = cudaMemcpy(d_aradx_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_arady_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);
  err = cudaMemcpy(d_aradz_arr, ahdxl, sizeof(ldouble)*Ngrid, cudaMemcpyHostToDevice);  
  
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


// TODO replace get_x, get_xb, and get_gKr   everywherex

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
/*
__device__ __host__ ldouble get_gKr_device(ldouble* gKr_arr, int i,int j, int k,
				  int ix, int iy, int iz)
{
  ldouble gKr_out = gKr_arr[i*4*4+j*4+k + (iX(ix)+(NGCX))*64 + \
				          (iY(iy)+(NGCY))*(SX)*64 + \
			                  (iZMET(iz)+(NGCZMET))*(SY)*(SX)*64];
  return gKr_out;
}
*/
// get size of cell indexed ic in dimension idim
// copied from get_size_x in finite.c
__device__ __host__ ldouble get_size_x_device(ldouble* xb_arr, int ic, int idim)
{
  ldouble dx;
  dx = get_xb_device(xb_arr, ic+1,idim) - get_xb_device(xb_arr, ic, idim);
  return dx;
}


// Metric source term
// TODO: deleted RADIATION and SHEARINGBOX parts
__device__ __host__ int f_metric_source_term_device(int ix, int iy, int iz, ldouble* ss,
		      	                            ldouble* p_arr, ldouble* x_arr,
			                            ldouble* g_arr, ldouble* G_arr, ldouble* gKr_arr)
{

     
  struct geometry geom;
  fill_geometry_device(ix,iy,iz,&geom,x_arr,g_arr,G_arr);

      
  ldouble (*gg)[5],(*GG)[5],gdetu;
  ldouble *pp = &get_u(p_arr,0,ix,iy,iz);
  gg=geom.gg;
  GG=geom.GG;

  #if (GDETIN==0) //no metric determinant inside derivatives
  gdetu=1.;
  #else
  gdetu=geom.gdet;
  #endif
  
  ldouble T[4][4];
  //calculating stress energy tensor components
  calc_Tij_device(pp,&geom,T); 
  indices_2221_device(T,T,gg);

  for(int i=0;i<4;i++)
  {
    for(int j=0;j<4;j++)
    {
	if(isnan(T[i][j])) 
	{
	    printf("%d %d %e\n",i,j,T[i][j]);
	    printf("nan in metric_source_terms\n");
	    //my_err("nan in metric_source_terms\n");//TODO
	}
    }
  }
  
  // zero out all source terms initially
  for(int iv=0;iv<NV;iv++)
    ss[iv]=0.;  

  //terms with Christoffels
  for(int k=0;k<4;k++)
  {
    for(int l=0;l<4;l++)
    {
      ss[1]+=gdetu*T[k][l]*get_gKr(gKr_arr,l,0,k,ix,iy,iz);
      ss[2]+=gdetu*T[k][l]*get_gKr(gKr_arr,l,1,k,ix,iy,iz);
      ss[3]+=gdetu*T[k][l]*get_gKr(gKr_arr,l,2,k,ix,iy,iz);
      ss[4]+=gdetu*T[k][l]*get_gKr(gKr_arr,l,3,k,ix,iy,iz);       
    }
  }


#if (GDETIN==0)
  //gdet derivatives
  ldouble dlgdet[3];
  dlgdet[0]=gg[0][4]; //D[gdet,x1]/gdet
  dlgdet[1]=gg[1][4]; //D[gdet,x2]/gdet
  dlgdet[2]=gg[2][4]; //D[gdet,x3]/gdet

  //get 4-velocity
  ldouble vcon[4],ucon[4];
  vcon[1]=pp[2];
  vcon[2]=pp[3];
  vcon[3]=pp[4];  
  conv_vels_device(vcon,ucon,VELPRIM,VEL4,gg,GG); 

  //terms with dloggdet  
  for(int l=1;l<4;l++)
  {
    ss[0]+=-dlgdet[l-1]*pp[RHO]*ucon[l];
    ss[1]+=-dlgdet[l-1]*(T[l][0]+pp[RHO]*ucon[l]);
    ss[2]+=-dlgdet[l-1]*(T[l][1]);
    ss[3]+=-dlgdet[l-1]*(T[l][2]);
    ss[4]+=-dlgdet[l-1]*(T[l][3]);
    ss[5]+=-dlgdet[l-1]*pp[ENTR]*ucon[l];
  }   
#endif
  
  return 0;
}

//**********************************************************************
// calculate stress energy tensor
//**********************************************************************
__device__ __host__ int calc_Tij_device(ldouble *pp, void* ggg, ldouble T[][4])
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble utcon[4],ucon[4],ucov[4];  
  ldouble bcon[4],bcov[4],bsq=0.;
  
  //converts to 4-velocity
  utcon[0]=0.;
  for(int iv=1;iv<4;iv++)
    utcon[iv]=pp[1+iv];
  conv_vels_both_device(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

#ifdef NONRELMHD
  ucon[0]=1.;
  ucov[0]=-1.;
#endif

#ifdef MAGNFIELD
  calc_bcon_bcov_bsq_from_4vel_device(pp, ucon, ucov, geom, bcon, bcov, &bsq); 
#else
  bcon[0]=bcon[1]=bcon[2]=bcon[3]=0.;
  bsq=0.;
#endif
  
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  //gamma=pick_gammagas(geom->ix,geom->iy,geom->iz); //TODO
  #endif
  ldouble gammam1=gamma-1.;

  ldouble rho=pp[RHO];
  ldouble uu=pp[UU];  
  ldouble p=(gamma-1.)*uu; 
  ldouble w=rho+uu+p;
  ldouble eta=w+bsq;
  ldouble ptot=p+0.5*bsq;

#ifndef NONRELMHD  
  for(int i=0;i<4;i++)
    for(int j=0;j<4;j++)
      T[i][j]=eta*ucon[i]*ucon[j] + ptot*GG[i][j] - bcon[i]*bcon[j];
#else
  
  ldouble v2=dot3nr(ucon,ucov); //TODO
  for(int i=1;i<4;i++)
    for(int j=1;j<4;j++)
      T[i][j]=(rho)*ucon[i]*ucon[j] + ptot*GG[i][j] - bcon[i]*bcon[j];

  T[0][0]=uu + bsq/2. + rho*v2/2.;
  for(int i=1;i<4;i++)
    T[0][i]=T[i][0]=(T[0][0] + ptot) *ucon[i]*ucon[0] + ptot*GG[i][0] - bcon[i]*bcon[0];

#endif  // ifndef NONRELMHD

  return 0;
}


//**********************************************************************
// calculate total gas entropy from density & energy density
//**********************************************************************
__device__ __host__ ldouble calc_Sfromu_device(ldouble rho,ldouble u,int ix,int iy,int iz)
{
  ldouble gamma=GAMMA;
  #ifdef CONSISTENTGAMMA
  //gamma=pick_gammagas(ix,iy,iz); //TODO
  #endif
  ldouble gammam1=gamma-1.;
  ldouble indexn=1.0/gammam1;
  ldouble pre=gammam1*u;
  #ifdef NOLOGINS
  ldouble ret = rho*u / pow(rho,gamma);
  #else
  ldouble ret = rho*log(pow(pre,indexn)/pow(rho,indexn+1.));
  #endif

  return ret;
}


//**********************************************************************
// kernels
//**********************************************************************

__global__ void calc_update_kernel(int Nloop_0, 
                                   int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
		       	           ldouble* x_arr, ldouble* xb_arr,
                                   ldouble* gcov_arr, ldouble* gcon_arr, ldouble* gKr_arr,
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
     f_metric_source_term_device(ix,iy,iz,ms,p_arr, x_arr,gcov_arr,gcon_arr,gKr_arr);
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
  ldouble ag,al,ar,amax;
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
    ldouble aaa[24];
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
    ldouble aaa[24];
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
    ldouble aaa[24];
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
    for(int iv=0;i<NV;i++)
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
  calc_fluxes_kernel(Nloop_1,
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

