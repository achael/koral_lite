#pragma once

#define TB_SIZE 64


///////////////////////////////////////////////////////////////
// metricgpu.cu ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//Declaring constant geometrical arrays for device
extern int *d_loop0_ix, *d_loop0_iy, *d_loop0_iz;
extern int *d_loop1_ix, *d_loop1_iy, *d_loop1_iz;
extern int *d_loop4_ix, *d_loop4_iy, *d_loop4_iz;
extern ldouble *d_x;    //[(NX+NY+NZ+6*NG)*sizeof(ldouble)]
extern ldouble *d_xb;   //[(NX+1+NY+1+NZ+1+6*NG)*sizeof(ldouble)]
extern ldouble *d_gcov; //[SX*SY*SZMET*sizeof(ldouble)]
extern ldouble *d_gcon; //[SX*SY*SZMET*sizeof(ldouble)]
extern ldouble *d_gbx, *d_gby, *d_gbz;
extern ldouble *d_Gbx, *d_Gby, *d_Gbz;
extern ldouble *d_Kris; //[(SX)*(SY)*(SZMET)*64*sizeof(ldouble)]

__device__ __host__ int fill_geometry_device(int ix,int iy,int iz, void* geom,
				             ldouble* x_arr, ldouble* g_arr, ldouble* G_arr);
__host__ __device__ int fill_geometry_face_device(int ix,int iy,int iz,int idim, void *geom,
						  ldouble* x_arr, ldouble* xb_arr,
						  ldouble* gbx_arr,ldouble* gby_arr,ldouble* gbz_arr,
						  ldouble* Gbx_arr,ldouble* Gby_arr,ldouble* Gbz_arr);
///////////////////////////////////////////////////////////////
// finitegpu.cu ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

extern ldouble *d_u_arr;
extern ldouble *d_p_arr;

extern ldouble *d_flbx_arr,*d_flby_arr, *d_flbz_arr;
extern ldouble *d_pbLx_arr, *d_pbLy_arr, *d_pbLz_arr;
extern ldouble *d_pbRx_arr, *d_pbRy_arr, *d_pbRz_arr;
extern ldouble *d_flLx_arr, *d_flLy_arr, *d_flLz_arr;
extern ldouble *d_flRx_arr, *d_flRy_arr, *d_flRz_arr;
extern ldouble *d_ahdxl_arr, *d_ahdyl_arr, *d_ahdzl_arr;
extern ldouble *d_ahdxr_arr, *d_ahdyr_arr, *d_ahdzr_arr;
extern ldouble *d_ahdx_arr, *d_ahdy_arr, *d_ahdz_arr;
extern ldouble *d_aradl_arr, *d_aradl_arr, *d_aradl_arr;
extern ldouble *d_aradr_arr, *d_aradr_arr, *d_aradr_arr;
extern ldouble *d_aradx_arr, *d_arady_arr, *d_aradz_arr;

extern ldouble *d_emf_arr;

extern int *d_cellflag_arr;
extern int *d_int_slot_arr;

__device__ __host__ int is_cell_active_device(int ix, int iy, int iz);
__device__ __host__ int is_cell_corrected_polaraxis_device(int ix, int iy, int iz);
__device__ __host__ ldouble get_x_device(ldouble* x_arr, int ic, int idim);
__device__ __host__ ldouble get_xb_device(ldouble* xb_arr, int ic, int idim);
__device__ __host__ ldouble get_gKr_device(ldouble* gKr_arr, int i,int j, int k,
				           int ix, int iy, int iz);

__device__ __host__ ldouble get_size_x_device(ldouble* xb_arr, int ic, int idim);

// kernels

__global__ void calc_update_kernel(int Nloop_0, 
                                   int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
		        	   ldouble* x_arr, ldouble* xb_arr,
                                   ldouble* gcov_arr, ldouble* gcon_arr, ldouble* gKr_arr,
				   ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr,
				   ldouble* u_arr, ldouble* p_arr, ldouble dtin);
__global__ void calc_fluxes_kernel(int Nloop_1,
                                   int* loop_1_ix, int* loop_0_iy, int* loop_1_iz,
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
				   ldouble* aradx_arr,  ldouble* arady_arr,  ldouble* aradz_arr,	   				    ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr);


///////////////////////////////////////////////////////////////
// relelegpu.cu ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

__device__ __host__ int indices_12_device(ldouble A1[4],ldouble A2[4],ldouble GG[][5]);
__device__ __host__ int indices_21_device(ldouble A1[4],ldouble A2[4],ldouble gg[][5]);
__device__ __host__ int indices_2211_device(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5]);
__device__ __host__ int indices_2221_device(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5]);

__device__ __host__ int calc_ucon_ucov_from_prims_device(ldouble *pr, void *ggg, ldouble *ucon, ldouble *ucov);
__device__ __host__ int conv_vels_device(ldouble *u1,ldouble *u2,
					 int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
__device__ __host__ int conv_vels_ut_device(ldouble *u1,ldouble *u2,
					    int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
__device__ __host__ int conv_vels_both_device(ldouble *u1,ldouble *u2con,ldouble *u2cov,
					      int which1,int which2,ldouble gg[][5],ldouble GG[][5]);
__device__ __host__ int conv_vels_core_device(ldouble *u1,ldouble *u2,
					      int which1,int which2,ldouble gg[][5],ldouble GG[][5],int);

__device__ __host__ ldouble calc_alpgam_device(ldouble *u1, ldouble gg[][5], ldouble GG[][5]);
__device__ __host__ int fill_utinvel3_device(ldouble *u1,double gg[][5],ldouble GG[][5]);
__device__ __host__ int fill_utinucon_device(ldouble *u1,double gg[][5],ldouble GG[][5]);
__device__ __host__ int calc_normalobs_ncon_device(ldouble GG[][5], ldouble alpha, ldouble *ncon);

///////////////////////////////////////////////////////////////
// physics.cu  (TODO ???) /////////////////////////////////////
///////////////////////////////////////////////////////////////

__device__ __host__ int calc_Tij_device(ldouble *pp, void* ggg, ldouble T[][4]);
__device__ __host__ ldouble calc_Sfromu_device(ldouble rho,ldouble u,int ix,int iy,int iz);
__device__ __host__ int f_metric_source_term_device(int ix, int iy, int iz, ldouble* ss,
			                            ldouble* p_arr, ldouble* x_arr,
			                            ldouble* g_arr, ldouble* G_arr, ldouble* l_arr);
  
///////////////////////////////////////////////////////////////
// magngpu.cu   ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

__device__ __host__ void calc_bcon_bcov_bsq_from_4vel_device(ldouble *pr, ldouble *ucon,
							     ldouble *ucov, void* ggg,
		                        		     ldouble *bcon, ldouble *bcov, ldouble *bsq);

// kernels
__global__ void flux_ct_setemf_kernel(int Nloop_4,
			              int* loop_4_ix, int* loop_4_iy, int* loop_4_iz,
			              ldouble* emf_arr,
			              ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr);

__global__ void flux_ct_getemf_kernel(int Nloop_4,
			              int* loop_4_ix, int* loop_4_iy, int* loop_4_iz,
			              ldouble* emf_arr,
			              ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr);

///////////////////////////////////////////////////////////////
// u2pgpu.cu  /////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
__device__ __host__ int set_cflag_device(int *cellflag_arr, int iflag,int ix,int iy,int iz, int val);

__device__ __host__ int p2u_device(ldouble *p, ldouble *u, void *ggg);
__device__ __host__ int p2u_mhd_device(ldouble *p, ldouble *u, void *ggg);
__device__ __host__ ldouble calc_utp1_device(ldouble *vcon, ldouble *ucon, void *ggg);
  
__device__ int u2p_device(ldouble *uu0, ldouble *pp, void *ggg,
			  int corrected[3], int fixups[2], int int_slot_arr[NGLOBALINTSLOT]);

__device__ __host__ int u2p_solver_Bonly_device(ldouble *uu, ldouble *pp, void *ggg);
__device__ __host__ int u2p_solver_device(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose);
__device__ __host__ int u2p_solver_W_device(ldouble *uu, ldouble *pp, void *ggg,int Etype,int verbose);
__device__ __host__ int check_floors_mhd_device(ldouble *pp, int whichvel,void *ggg);
  
// kernel
__global__ void calc_primitives_kernel(int Nloop_0,
				       int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
                                       ldouble *x_arr, ldouble *g_arr, ldouble *G_arr,
				       ldouble *u_arr, ldouble *p_arr,
				       int setflags, int* cellflag_arr, int int_slot_arr[NGLOBALINTSLOT]);
				       
