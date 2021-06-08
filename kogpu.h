#pragma once

///////////////////////////////////////////////////////////////
// finitegpu.cu ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

__device__ ldouble get_xb_device(ldouble* xb_arr, int ic, int idim);
__device__ ldouble get_size_x_device(ldouble* xb_arr, int ic, int idim);
__device__ int fill_geometry_device(int ix,int iy,int iz,void* geom,
				    ldouble* g_arr, ldouble* G_arr);
__device__ int f_metric_source_term_device(int ix, int iy, int iz, ldouble* ss,
			                   ldouble* p_arr,
			                   ldouble* g_arr, ldouble* G_arr, ldouble* l_arr);
			                   
//Defining metric and Krystoffels for device			                   
ldouble *d_gcov;//[SX*SY*SZMET*sizeof(ldouble)];
ldouble *d_gcon;//[SX*SY*SZMET*sizeof(ldouble)];
ldouble *d_Kris;//[(SX)*(SY)*(SZMET)*64*sizeof(ldouble)];
