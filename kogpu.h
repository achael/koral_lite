#pragma once

///////////////////////////////////////////////////////////////
// finitegpu.cu ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

__global__ ldouble get_xb_device(ldouble* xb_arr, int ic, int idim);
__global__ ldouble get_gKr_device(ldouble* gKr_arr, int i,int j, int k,
				  int ix, int iy, int iz);

__device__ ldouble get_size_x_device(ldouble* xb_arr, int ic, int idim);

__global__ int indices_2211_device(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5]);

__device__ int fill_geometry_device(int ix,int iy,int iz,void* geom,
				    ldouble* g_arr, ldouble* G_arr);
__device__ int f_metric_source_term_device(int ix, int iy, int iz, ldouble* ss,
			                   ldouble* p_arr,
			                   ldouble* g_arr, ldouble* G_arr, ldouble* l_arr);
__global__ void calc_update_gpu_kernel(ldouble dtin, int Nloop_0, 
                                       int* loop_0_ix, int* loop_0_iy, int* loop_0_iz,
				       ldouble* xb_arr,
				       ldouble* flbx_arr, ldouble* flby_arr, ldouble* flbz_arr,
				       ldouble* u_arr);


