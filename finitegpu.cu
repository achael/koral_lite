extern "C" {

#include "ko.h"

}

__global__ void calc_update_gpu_kernel(ldouble dtin, int Nloop_0)
{

  int ii;
  ii = blockIdx.x * blockDim.x + threadIdx.x;
  if(ii >= Nloop_0) return;
  
  /*
  int ix,iy,iz,iv;
  ix=loop_0[ii][0];
  iy=loop_0[ii][1];
  iz=loop_0[ii][2]; 

  // Source term
  //ldouble ms[NV],gs[NV],val,du;
      
  
  if(is_cell_active(ix,iy,iz)==0)
  {
    // Source terms applied only for active cells	  
    for(iv=0;iv<NV;iv++) ms[iv]=0.; 
  }
  else
  {
     // Get metric source terms ms[iv] and any other source terms gs[iv], and save combined source terms in ms[iv]
     f_metric_source_term(ix,iy,iz,ms);
     f_general_source_term(ix,iy,iz,gs);
     for(iv=0;iv<NV;iv++) ms[iv]+=gs[iv];
  }
      
      
  // Get the cell size in the three directions
  ldouble dx=get_size_x(ix,0);
  ldouble dy=get_size_x(iy,1);
  ldouble dz=get_size_x(iz,2);

  int doxl,doxr,doyl,doyr,dozl,dozr;
  doxl=doxr=doyl=doyr=dozl=dozr=1;

  //timestep
  //dt=dtin;  //TODO global?  // dtin is an input parameter to op_explicit
      
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

  printf("Test\n");fflush(stdout);

  int TB_SIZE = 64;   
  int threadblocks = (Nloop_0 / TB_SIZE) + ((Nloop_0 % TB_SIZE)? 1:0);
  // calc_update_gpu_kernel<<<threadblocks, TB_SIZE>>>(dtin, Nloop_0);
 
  return 0;
}
