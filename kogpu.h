#pragma once

#define TB_SIZE 64

///////////////////////////////////////////////////////////////
// finitegpu.cu ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//Declaring constant geometrical arrays for device			                   
extern int *d_loop0_ix, *d_loop0_iy, *d_loop0_iz; 
extern ldouble *d_x;    //[(NX+NY+NZ+6*NG)*sizeof(ldouble)]
extern ldouble *d_xb;   //[(NX+1+NY+1+NZ+1+6*NG)*sizeof(ldouble)]
extern ldouble *d_gcov; //[SX*SY*SZMET*sizeof(ldouble)]
extern ldouble *d_gcon; //[SX*SY*SZMET*sizeof(ldouble)]
extern ldouble *d_Kris; //[(SX)*(SY)*(SZMET)*64*sizeof(ldouble)]

__device__ __host__ int is_cell_active_device (int ix, int iy, int iz);
__device__ __host__ ldouble get_xb_device(ldouble* xb_arr, int ic, int idim);
__device__ __host__ ldouble get_gKr_device(ldouble* gKr_arr, int i,int j, int k,
				           int ix, int iy, int iz);

__device__ __host__ ldouble get_size_x_device(ldouble* xb_arr, int ic, int idim);


__device__ __host__ int fill_geometry_device(int ix,int iy,int iz, ldouble* x_arr, void* geom,
				             ldouble* g_arr, ldouble* G_arr);
__device__ __host__ int f_metric_source_term_device(int ix, int iy, int iz, ldouble* ss,
			                            ldouble* p_arr,
			                            ldouble* g_arr, ldouble* G_arr, ldouble* l_arr);
			                   


// kernels below

__global__ void calc_update_gpu_kernel(ldouble dtin, int Nloop_0, 
                                       int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
				       ldouble* x_arr, ldouble* xb_arr,
                                       ldouble* gcov_arr, ldouble* gcon_arr, ldouble* gKr_arr,
				       ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr,
				       ldouble* u_arr, ldouble* p_arr);


///////////////////////////////////////////////////////////////
// relelegpu.cu ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

__device__ __host__ int indices_21_device(ldouble A1[4],ldouble A2[4],ldouble gg[][5]);
__device__ __host__ int indices_2211_device(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5]);
__device__ __host__ int indices_2221_device(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5]);

__device__ __host__ int calc_ucon_ucov_from_prims_device(ldouble *pr, void *ggg, ldouble *ucon, ldouble *ucov);
__device__ __host__ int conv_vels_device(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
__device__ __host__ int conv_vels_ut_device(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
__device__ __host__ int conv_vels_both_device(ldouble *u1,ldouble *u2con,ldouble *u2cov,int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
__device__ __host__ int conv_vels_core_device(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5],int);
__device__ __host__ ldouble calc_alpgam_device(ldouble *u1, ldouble gg[][5], ldouble GG[][5]);

__device__ __host__ int fill_utinvel3_device(ldouble *u1,double gg[][5],ldouble GG[][5]);
__device__ __host__ int fill_utinucon_device(ldouble *u1,double gg[][5],ldouble GG[][5]);

__device__ __host__ int calc_Tij_device(ldouble *pp, void* ggg, ldouble T[][4]);

///////////////////////////////////////////////////////////////
// magn .cu  (TODO ???) ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

__device__ __host__ void calc_bcon_bcov_bsq_from_4vel_device(ldouble *pr, ldouble *ucon, ldouble *ucov, void* ggg,
		                        		     ldouble *bcon, ldouble *bcov, ldouble *bsq);

///////////////////////////////////////////////////////////////
// u2pgpu.cu  /////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
__device__ __host__ int u2p_device(ldouble *uu0, ldouble *pp, void *ggg, int corrected[3], int fixups[2]);

__device__ __host__ int u2p_solver_W_device(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose)
// kernel
__global__ void calc_primitives_kernel(int Nloop_0, int setflags,
				       int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
				       ldouble *u_arr, ldouble *p_arr,
				       ldouble *x_arr, ldouble *g_arr, ldouble *G_arr);
