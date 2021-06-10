/*! \file relelegpu.cu
 \brief Some relativistic routines
*/

extern "C" {

#include "ko.h"

}

#include "kogpu.h"


__device__ __host__ int indices_12_device(ldouble A1[4],ldouble A2[4],ldouble GG[][5])
{

  for(int i = 0; i < 4; i++)
  {
    A2[i] = 0.;
  }

  for(int i = 0; i < 4; i++)
  {
    for(int k = 0; k < 4; k++)
    {
      A2[i] += A1[k] * GG[i][k];
    }
  }

  return 0;
}

__device__ __host__ int indices_21_device(ldouble A1[4],ldouble A2[4],ldouble gg[][5])
{
  
  ldouble At[4];

  for(int i=0;i<4;i++)
  {
    At[i]=0.;
  }

  for(int i=0;i<4;i++)
  {
    for(int j = 0; j < 4; j++)
    {
      At[i] += A1[j] * gg[i][j];
    }
  }

  for(int i=0; i<4; i++)
  {
    A2[i] = At[i];
  }
  
  return 0;
}

__device__ __host__ int indices_2211_device(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5])
{

  ldouble Tt[4][4];

  for(int i=0;i<4;i++)
  {
     for(int j=0;j<4;j++)
     {
       Tt[i][j]=0.;
       for(int k=0;k<4;k++)
       {
	 for(int l=0;l<4;l++)
	 {
	   Tt[i][j]+=T1[k][l]*gg[i][k]*gg[j][l];
	 }	  
       }
    }
  }

  for(int i=0;i<4;i++)
  {
    for(int j=0;j<4;j++)
    {
      T2[i][j]=Tt[i][j];
    }
  }

  return 0;
}


__device__ __host__ int indices_2221_device(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5])
{
 
  ldouble Tt[4][4];

  for(int i=0;i<4;i++)
  {
    for(int j=0;j<4;j++)
    {
      Tt[i][j]=0.;
    }
  }

  for(int i=0;i<4;i++)
  {
    for(int j=0;j<4;j++)
    {
      for(int k=0;k<4;k++)
      {
        Tt[i][j]+=T1[i][k]*gg[k][j];
      }
    }
  }

  for(int i=0;i<4;i++)
  {
    for(int j=0;j<4;j++)
    {
      T2[i][j]=Tt[i][j];
    }
  }

  return 0;
}

//*********************************************************************
//Takes primitives and computes ucon, ucov in VEL4 frame
//**********************************************************************

__device__ __host__ int calc_ucon_ucov_from_prims_device(ldouble *pr, void *ggg, ldouble *ucon, ldouble *ucov)
{
  struct geometry *geom
  = (struct geometry *) ggg;

  ucon[0] = 0.;
  ucon[1] = pr[VX];
  ucon[2] = pr[VY];
  ucon[3] = pr[VZ];
  
#ifdef NONRELMHD //only three-velocity used;
  fill_utinucon_device(ucon,geom->gg,geom->GG); 
  indices_21_device(ucon,ucov,geom->gg);
  
#else
  
  conv_vels_both_device(ucon,ucon,ucov,VELPRIM,VEL4,geom->gg,geom->GG);

#endif
  
  return 0;
}


//*********************************************************************
//computes ut and then calculates ucon
//**********************************************************************

__device__ __host__ int conv_vels_device(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5])
{
  
#ifdef NONRELMHD //only three-velocity used;
  u2[1]=u1[1];u2[2]=u1[2];u2[3]=u1[3];
  fill_utinucon_device(u2,gg,GG); 

#else

  conv_vels_core_device(u1,u2,which1,which2,gg,GG,0);  // 0 means u^t is not yet known

#endif
  return 0;
}


//**********************************************************************
//calculates ucon, assuming ut is known
//**********************************************************************

__device__ __host__ int conv_vels_ut_device(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5])
{
  
#ifdef NONRELMHD //only three-velocity used;
  u2[1]=u1[1];u2[2]=u1[2];u2[3]=u1[3];
  fill_utinucon_device(u2,gg,GG); 

#else
  
  conv_vels_core_device(u1,u2,which1,which2,gg,GG,1);  // 1 means u^t is known

#endif
  
  return 0;
}


//**********************************************************************
//calculates both ucon and ucov, assuming ut is unknown 
//**********************************************************************

__device__ __host__ int conv_vels_both_device (ldouble *u1,ldouble *u2con,ldouble *u2cov,int which1,int which2,ldouble gg[][5],ldouble GG[][5])
{
  
#ifdef NONRELMHD //only three-velocity used;
  u2con[1]=u1[1];u2con[2]=u1[2];u2con[3]=u1[3];
  fill_utinucon_device(u2con,gg,GG); 
  indices_21_device(u2con,u2cov,gg);

#else

  if(which2!=VEL4)
  {
    printf("conv_vels_both only works with which2==VEL4: %d -> %d\n",which1,which2);
    return -1;
  }
  
  conv_vels_core_device(u1,u2con,which1,which2,gg,GG,0); //0 means u^t is not yet known
  indices_21_device(u2con,u2cov,gg);

#endif
  return 0;
}


//**********************************************************************
//converts contravariant velocities u1 to contravariant u2con and covariant u2cov (if which2==VEL4)
// July 7, 17, Ramesh: Large reorganization
// sub-calculations done in fill_utinucon, fill_utinvel3. This version has been tested with test_con_vel.c.
//**********************************************************************

__device__ __host__ int conv_vels_core_device(ldouble *u1,ldouble *u2conout,int which1,int which2,
					      ldouble gg[][5],ldouble GG[][5],int utknown)
{
  
  ldouble u2con[4];

  int verbose=0;
  if(verbose)
  {
    printf("conv_vels: %d -> %d\n",which1,which2);
    //print_4vector(u1); //TODO
  }

  /*************** VEL3 -> VEL3 ***************/
  if(which1==VEL3 && which2==VEL3)
  {
    for(int i=0;i<4;i++)
    {
      u2con[i]=u1[i];
    }
  }
  
  /*************** VEL4 -> VEL4 ***************/
  else if(which1==VEL4 && which2==VEL4)
  {
    if(utknown==0)  // u^t is not known
    {
      fill_utinucon_device(u1, gg, GG);
    }
    
    for(int i=0;i<4;i++)
    {
      u2con[i]=u1[i];
    }
  }
  
  /*************** VELR -> VELR ***************/
  else if(which1==VELR && which2==VELR)
  {
    for(int i=0;i<4;i++)
    {
      u2con[i]=u1[i];
    }
  }
  
  /*************** VEL4 -> VEL3 ***************/
  else if(which1==VEL4 && which2==VEL3)
  {
    if(utknown==0)  // u^t is not known
    {
      fill_utinucon_device(u1, gg, GG);
    }
    
    for(int i=0;i<4;i++)
    {
      u2con[i]=u1[i]/u1[0];
    }
  }
  
  /*************** VEL3 -> VEL4 ***************/
  else if(which1==VEL3 && which2==VEL4)
  {
    fill_utinvel3_device(u1, gg, GG);
    u2con[0] = u1[0];
    
    if(u2con[0] < 1. || isnan(u2con[0]))
    {
      printf("ut.nan in conv_vels(%d,%d) VEL3->VEL4 - returning error\n",which1,which2);
      return -1;  
    }
    
    u2con[1] = u1[1] * u2con[0];
    u2con[2] = u1[2] * u2con[0];
    u2con[3] = u1[3] * u2con[0];
  }
  
  /*************** VEL3 -> VELR ***************/
  else if(which1==VEL3 && which2==VELR)
  {
    fill_utinvel3_device(u1, gg, GG);
    u2con[0] = u1[0];
    
    if(u2con[0] < 1. || isnan(u2con[0]))
    {
      printf("ut.nan in conv_vels(%d,%d) VEL3->VELR - returning error\n",which1,which2);
      return -1;
    }
    
    //to 4-velocity
    u2con[1] = u1[1] * u2con[0];
    u2con[2] = u1[2] * u2con[0];
    u2con[3] = u1[3] * u2con[0];
    
    //to relative velocity
    for(int i=1; i<4; i++)
    {
      u2con[i] = u2con[i] - u2con[0] * GG[0][i] / GG[0][0];
    }
  }
  
  /*************** VEL4 -> VELR ***************/
  else if(which1==VEL4 && which2==VELR)
  {
    if(utknown==0)  // u^t is not known
    {
      fill_utinucon_device(u1, gg, GG);
    }
    u2con[0] = u1[0];

    // to relative velocity
    for(int i=1; i<4; i++)
    {
      u2con[i] = u1[i] - u2con[0] * GG[0][i] / GG[0][0];
    }
  }

  /*************** VELR -> VEL4 ***************/
  else if(which1==VELR && which2==VEL4)
  {
    ldouble alpgam = calc_alpgam_device(u1, gg, GG);
    
    u2con[0]=-alpgam*GG[0][0];
    if(u2con[0]<0) // TODO warning? 
      u2con[0] = fabs(u2con[0]);
          
    u2con[1]=u1[1]-alpgam*GG[0][1];
    u2con[2]=u1[2]-alpgam*GG[0][2];
    u2con[3]=u1[3]-alpgam*GG[0][3];
  }
  
  /*************** VELR -> VEL3 ***************/
  else if(which1==VELR && which2==VEL3)
  {
    ldouble alpgam = calc_alpgam_device(u1, gg, GG);

    u2con[0]=-alpgam*GG[0][0];
    if(u2con[0]<0)
      u2con[0] = fabs(u2con[0]);

    u2con[1]=u1[1]/u2con[0] + GG[0][1]/GG[0][0];
    u2con[2]=u1[2]/u2con[0] + GG[0][2]/GG[0][0];
    u2con[3]=u1[3]/u2con[0] + GG[0][3]/GG[0][0];

  }

  /*************** not supported  ***************/
  else
  {
    //my_err("velocity conversion not supported.\n");
    return -1;
  }

  for(int i=0; i<4; i++)
  {
    u2conout[i] = u2con[i];
  }

  if(verbose)
  {
    //print_4vector(u2con); //TODO 
    printf("conv_vels done %d %d\n",which1,which2);
  }
  
  return 0;
}


//**********************************************************************
// July 9, 17, Ramesh: This is Andrew's version of alpgam, which ensures a positive quantity
//**********************************************************************

__device__ __host__ ldouble calc_alpgam_device(ldouble *u1, ldouble gg[][5], ldouble GG[][5])
{

  ldouble qsq=0.;
  
  for(int i=1;i<4;i++)
  {
    for(int j=1;j<4;j++)
    {
      qsq+=u1[i]*u1[j]*gg[i][j];
    }
  }
  
  ldouble gamma2=(1. + qsq);
  ldouble alpha2=(-1./GG[0][0]);
  ldouble alpgam2=alpha2*gamma2;
  if(alpgam2<0.)
  {
    //printf("alpgam2.lt.0 in VELR->VEL4\n");
    return 1.;
  }

  ldouble alpgam=sqrt(alpgam2);
  
  return alpgam;
}

//**********************************************************************
// July 7, 17, Ramesh: Calculates u^t from spatial components of three-velocity VEL3
// We solve: ut^2 * (a + 2*b + c) = -1
//   where a = g00, b = g0i*ui, c = gij*ui*uj
//   solution: ut = sqrt(-1/(a + 2*b + c))
//**********************************************************************

__device__ __host__ int fill_utinvel3_device(ldouble *u1,double gg[][5],ldouble GG[][5])
{

  ldouble a, b, c;
  a = gg[0][0];
  b = c = 0.;
  
  for(int i=1; i<4; i++)
  {
    b += u1[i] * gg[0][i];
    
    for(int j=1;j<4;j++)
    {
      c += u1[i] * u1[j] * gg[i][j];
    }
  }
  
  u1[0]=sqrt(-1. / (a + 2. * b + c));
  
  return 0;
}


//**********************************************************************
// Calculates u^t from spatial components of four-velocity u^mu
// July 7, 17, Ramesh: streamlined the code to improve efficiency
// We solve quadratic: a*ut^2 + 2*b*ut + c = 0
//   where a = g00, b = g0i*ui, c = 1 + gij*ui*uj
//   solution: ut = (-b +/- sqrt(b^2-a*c))/a
//**********************************************************************
__device__ __host__ int fill_utinucon_device(ldouble *u1,double gg[][5],ldouble GG[][5])
{

  ldouble a, b, c, delta;
  
  a = gg[0][0];
  b = 0.;
  c = 1.;
  
  for(int i=1; i<4; i++)
  {
    b += u1[i] * gg[0][i];
    
    for(int j=1; j<4; j++)
    {
      c += u1[i] * u1[j] * gg[i][j];
    }
  }
  
  delta = b * b - a * c;  // Note: b here is half the usual value
  if(delta < 0.)
  {
    printf("delta.lt.0 in fill_utinucon\n");
    //my_err("delta.lt.0 in fill_utinucon\n");
  }
  
  if(a < 0.)
  {
    u1[0] = (-b - sqrt(delta)) / a;
  }
  else //this is in ergoregion
  {
    //u1[0] = (-b + sqrt(delta)) / a; //ANDREW THIS WAS WRONG, should be minus sign everywhere
    u1[0] = (-b - sqrt(delta)) / a;
  }
  
  return 0;
}

__device__ __host__ void calc_bcon_bcov_bsq_from_4vel_device(ldouble *pr, ldouble *ucon, ldouble *ucov, void* ggg,
		                        		     ldouble *bcon, ldouble *bcov, ldouble *bsq)
{

  struct geometry *geom
    = (struct geometry *) ggg;

  // First calculate bcon0
  bcon[0] = pr[B1]*ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3] ;
  
  // Then spatial components of bcon
#ifdef NONRELMHD
  for(int j = 1; j < 4; j++)
  {
    bcon[j] = pr[B1-1+j]; //b^i=B^i
  }
#else  // relativistic case
  
  ldouble u0inv = 1. / ucon[0];
  for(int j=1;j<4;j++)
  {
    bcon[j] = (pr[B1-1+j] + bcon[0] * ucon[j]) * u0inv ;
  }
#endif //NONRELMHD
  
  // Convert to bcov and calculate bsq
  indices_21_device(bcon, bcov, geom->gg);
  *bsq = dotB(bcon, bcov); //NOTE: preprocessor macro, ok

  return ;
}

//**********************************************************************
//returns contravariant four-velocity of a normal observer, given the value of alpha
//n_mu=(-alp,0,0,0);  nu^mu = nu_0 * GG[mu][0]
//**********************************************************************

__device__ __host__ int calc_normalobs_ncon_device(ldouble GG[][5], ldouble alpha, ldouble *ncon)
{

  for(int i=0; i<4; i++)
  {
    ncon[i] = -alpha * GG[i][0];
  }
  
  return 0.;
}
