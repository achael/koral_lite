/*! \file frames.c
 \brief Boosting, moving indices etc.
 */


#include "ko.h"


/*****************************************************************/
/********** all primitives between coordinates  ***************/
/*****************************************************************/

int
trans_pall_coco(ldouble *pp1, ldouble *pp2, int CO1,int CO2, ldouble *xxvec, void* ggg1, void* ggg2)
{
  trans_pmhd_coco(pp1, pp2, CO1,CO2,xxvec,ggg1,ggg2);
#ifdef RADIATION
  trans_prad_coco(pp1, pp2, CO1,CO2,xxvec,ggg1,ggg2);
#endif
  return 0;
}

/*****************************************************************/
/********** hydro primitives (E,F^i) between coordinates  ********/
/********** does not touch radiative primitives ******************/
/*****************************************************************/

int
trans_pmhd_coco(ldouble *ppin, ldouble *ppout, int CO1,int CO2, ldouble *xxvec, void* ggg1,void* ggg2)
{
  struct geometry *geom1
  = (struct geometry *) ggg1;
  struct geometry *geom2
  = (struct geometry *) ggg2;
  
  int i,iv;
  ldouble pp1[NV],pp2[NV];
  for(i=0;i<NV;i++)
  {
    ppout[i]=ppin[i];
    pp1[i]=ppin[i];
    pp2[i]=ppout[i];
  }
  
  if(CO1==CO2)
  {    
    for (i = 0; i < 5; i++)
    {
      pp2[i] = pp1[i];
    }
  }
  else
  {
    pp2[0]=pp1[0];
    pp2[1]=pp1[1];
    
    //bcon in CO1
    ldouble ucon[4], ucov[4],uconback[4];    
    calc_ucon_ucov_from_prims(pp1, geom1, ucon, ucov);
    

#ifdef MAGNFIELD
    ldouble bcon[4],Bcon[4];
    
    //magnetic field 4-vector
    calc_bcon_4vel(pp1, ucon, ucov, bcon);
#endif

    for(iv=0;iv<4;iv++)
      uconback[iv]=ucon[iv];
    
    //convert ucon (and bcon) to CO2
    trans2_coco(xxvec,ucon,ucon,CO1,CO2);
    
#ifdef MAGNFIELD
    trans2_coco(xxvec,bcon,bcon,CO1,CO2);
#endif

    
    //to VELPRIM
    conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geom2->gg,geom2->GG);
    
    pp2[2]=ucon[1];
    pp2[3]=ucon[2];
    pp2[4]=ucon[3];
    
#ifdef MAGNFIELD

    //coming back to primitive B^i
    calc_Bcon_prim(pp2,bcon,Bcon,geom2);
        
    for (i = 0; i < 3; i++)
    {
      pp2[B1+i] = Bcon[1+i];
    }
#endif
  }
  
  for(i=0;i<NVMHD;i++)
  {
    ppout[i]=pp2[i];
  }
  
  return 0;
}


/*****************************************************************/
/********** radiative primitives (E,F^i) between coordinates  ****/
/********** does not touch hydro primitives **********************/
/*****************************************************************/

int
trans_prad_coco(ldouble *ppin, ldouble *ppout, int CO1,int CO2, ldouble *xxvec, void* ggg1, void* ggg2)
{
  
  struct geometry *geom1
    = (struct geometry *) ggg1;
  struct geometry *geom2
    = (struct geometry *) ggg2;

  int i;
  
  ldouble pp1[NV],pp2[NV];
  for(i=0;i<NV;i++) 
    {
      ppout[i]=ppin[i];
      pp1[i]=ppin[i];
      pp2[i]=ppout[i];
    }      
#ifdef RADIATION 
  if(CO1==CO2)
    {
      for(i=0;i<4;i++)
	pp2[EE0+i]=pp1[EE0+i];
     }
  else
    {
      //Erf unchanged
      pp2[EE0]=pp1[EE0];

      //velocity in CO1
      ldouble ucon[4];
      ucon[0]=0;
      ucon[1]=pp1[FX0];
      ucon[2]=pp1[FY0];
      ucon[3]=pp1[FZ0];

      conv_vels(ucon,ucon,VELPRIMRAD,VEL4,geom1->gg,geom1->GG);

      //converting to CO2
      trans2_coco(xxvec,ucon,ucon,CO1,CO2);

      //to VELPRIM
      conv_vels_ut(ucon,ucon,VEL4,VELPRIMRAD,geom2->gg,geom2->GG);

      pp2[FX0]=ucon[1]; 
      pp2[FY0]=ucon[2];
      pp2[FZ0]=ucon[3];

    }
#endif //RADIATION
  for(i=NVMHD;i<NV;i++)     
    {
      ppout[i]=pp2[i];
    }      
 
  return 0;
}

/*****************************************************************/
//calculates Lorenz matrix for lab -> ff
/*****************************************************************/

int
calc_Lorentz_lab2ff(ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble L[][4])
{
  int i,j,k;
  int verbose=0;

  //calculating the four-velocity of fluid in lab frame
  ldouble utcon[4],ucon[4],ucov[4],vpr[3];
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

  if(verbose>0) print_4vector(ucon);

  //four velocity of the lab frame
  ldouble alpha=sqrt(-1./GG[0][0]);
  ldouble wcon[4],wcov[4]={-alpha,0.,0.,0.};
  indices_12(wcov,wcon,GG);

  if(verbose>0) print_4vector(wcon);

  //temporary Om matrix
  ldouble Om[4][4];

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Om[i][j]=ucon[i]*wcov[j]-wcon[i]*ucov[j];
  
  //Lorentz factor = -w^mu u_mu
  ldouble gam=-dot(wcon,ucov);

  ldouble Omsum;
  
  //Lorentz matrix components
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Omsum=0.;
	for(k=0;k<4;k++)
	  Omsum+=Om[i][k]*Om[k][j];
	
	L[i][j]=kron(i,j)+1./(1.+gam)*Omsum+Om[i][j];
      }
  return 0;
}


int
calc_Lorentz_lab2ff_4vel(ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble L[][4], ldouble ucon[4], ldouble ucov[4])
{
  int i,j,k;
  int verbose=0;
    
  if(verbose>0) print_4vector(ucon);
  
  //four velocity of the lab frame
  ldouble alpha=sqrt(-1./GG[0][0]);
  ldouble wcon[4],wcov[4]={-alpha,0.,0.,0.};
  indices_12(wcov,wcon,GG);
  
  if(verbose>0) print_4vector(wcon);
  
  //temporary Om matrix
  ldouble Om[4][4];
  
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Om[i][j]=ucon[i]*wcov[j]-wcon[i]*ucov[j];
  
  //Lorentz factor = -w^mu u_mu
  ldouble gam=-dot(wcon,ucov);
  
  ldouble Omsum;

  //Lorentz matrix components
  ldouble one_over_one_plus_gam = 1./(1.+gam);
  for(i = 0;i < 4; i++)
  {
    for(j = 0; j < 4; j++)
    {
      Omsum = 0.;
      for(k = 0;k < 4; k++)
      {
        Omsum += Om[i][k] * Om[k][j];
      }
      
      L[i][j] = kron(i,j) + one_over_one_plus_gam * Omsum + Om[i][j];
    }
  }
  
  return 0;
}


/*****************************************************************/
//calculates Lorenz matrix for ff -> lab
/*****************************************************************/

int
calc_Lorentz_ff2lab(ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble L[][4])
{
  int i,j,k;
  int verbose=0;

  //calculating the four-velocity of fluid in lab frame
  ldouble wcon[4],wtcon[4],wcov[4];
  wtcon[1]=pp[2];
  wtcon[2]=pp[3];
  wtcon[3]=pp[4];
  conv_vels_both(wtcon,wcon,wcov,VELPRIM,VEL4,gg,GG);
 
  if(verbose>0) print_4vector(wcon);

  //four velocity of the lab frame
  ldouble alpha=sqrt(-1./GG[0][0]);
  ldouble ucon[4],ucov[4]={-alpha,0.,0.,0.};
  indices_12(ucov,ucon,GG);

  if(verbose>0) print_4vector(ucon);

  //temporary Om matrix
  ldouble Om[4][4];

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Om[i][j]=ucon[i]*wcov[j]-wcon[i]*ucov[j];
  
  //Lorentz factor = -w^mu u_mu
  ldouble gam=-dot(wcon,ucov);

  ldouble Omsum;

  //Lorentz matrix components
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Omsum=0.;
	for(k=0;k<4;k++)
	  Omsum+=Om[i][k]*Om[k][j];
	
	L[i][j]=kron(i,j)+1./(1.+gam)*Omsum+Om[i][j];
      }
  return 0;
}

/*****************************************************************/
//T^ij Lorentz boost from lab to fluid frame
/*****************************************************************/

int
boost22_lab2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5])
{ 
  int i,j,k,l;
  ldouble Tt[4][4];

  int verbose=0;

  if(verbose>0) print_tensor(T1);

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_lab2ff(pp,gg,GG,L);


  //copying 
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    }

  
  //correcting for ortonormality
  ldouble alpha=sqrt(-1./GG[0][0]);
  for(i=0;i<4;i++)
    {
      T2[i][0]*=alpha;
      T2[0][i]*=alpha;
    }
      

  if(verbose>0) print_tensor(T2);

  if(verbose>0) getchar();

  return 0;
}


/*****************************************************************/
//T^ij Lorentz boost from fluid frame to lab
/*****************************************************************/

int
boost22_ff2lab(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5])
{ 
  int i,j,k,l;
  ldouble Tt[4][4];

  int verbose=0;

  if(verbose>0) print_tensor(T1);

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_ff2lab(pp,gg,GG,L);

  //copying
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
  //correcting for ortonormality
  ldouble alpha=sqrt(-1./GG[0][0]);
  for(i=0;i<4;i++)
    {
      T1[i][0]/=alpha;
      T1[0][i]/=alpha;
    }
 
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    } 


  if(verbose>0) print_tensor(T2);

  if(verbose>0) getchar();

  return 0;
}


/*****************************************************************/
//T^ij Lorentz boost from radiation rest frame to lab
/*****************************************************************/

int
boost22_rf2lab(ldouble T1[][4],ldouble T2[][4],ldouble *pp0,ldouble gg[][5],ldouble GG[][5])
{ 
#ifdef LABRADFLUXES
  my_err("boost22_rf2lab() not working for LABRADFLUXES\n");
#endif

  int i,j,k,l;
  ldouble Tt[4][4];
  ldouble pp[NV];
  for(i=0;i<NV;i++)
    pp[i]=pp0[i];

  int verbose=0;

  if(verbose>0) print_tensor(T1);

  //artificial and temporary substitution
  ldouble urf[4]={0.,pp[FX],pp[FY],pp[FZ]};
  conv_vels(urf,urf,VELPRIMRAD,VELPRIM,gg,GG);
  pp[2]=urf[1];
  pp[3]=urf[2];
  pp[4]=urf[3];

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_ff2lab(pp,gg,GG,L);

  //copying
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
 
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    } 

  //dividing by lapse to express T2 in no-frame
  ldouble alpha=sqrt(-1./GG[0][0]);
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=T2[i][j]/alpha;
	}
    }

  if(verbose>0) print_tensor(T2);

  if(verbose>0) getchar();

  return 0;
}


/*****************************************************************/
//T^ij Lorentz boost from lab frame to radiation rest frame
/*****************************************************************/

int
boost22_lab2rf(ldouble T1[][4],ldouble T2[][4],ldouble *pp0,ldouble gg[][5],ldouble GG[][5])
{ 
#ifdef LABRADFLUXES
  my_err("boost22_lab2rf() not working for LABRADFLUXES\n");
#endif


  int i,j,k,l;
  ldouble Tt[4][4];
  ldouble pp[NV];
  for(i=0;i<NV;i++)
    pp[i]=pp0[i];

  int verbose=0;

  if(verbose>0) print_tensor(T1);

  //artificial and temporary substitution
  ldouble urf[4]={0.,pp[FX],pp[FY],pp[FZ]};
  conv_vels(urf,urf,VELPRIMRAD,VELPRIM,gg,GG);
  pp[2]=urf[1];
  pp[3]=urf[2];
  pp[4]=urf[3];

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_lab2ff(pp,gg,GG,L);

  //multiplying by lapse to express T1 in ZAMO
  ldouble alpha=sqrt(-1./GG[0][0]);
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j]*alpha;
	}
    }

  //copying
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
   
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    } 

  if(verbose>0) print_tensor(T2);

  if(verbose>0) getchar();

  return 0;
}


/*****************************************************************/
//A^i Lorentz boost from lab to fluid frame
/*****************************************************************/

int
boost2_lab2ff(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5])
{ 
  int i,j,k,l;
  ldouble At[4]   ;

  int verbose=0;

  if(verbose>0) print_4vector(A1);

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_lab2ff(pp,gg,GG,L);

  //copying
  for(i=0;i<4;i++)
    {
      At[i]=A1[i];
    }

  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      A2[i]=0.;
      for(k=0;k<4;k++)
	{
	  A2[i]+=L[i][k]*At[k];
	}
    }

  //laps to make it ortonormal
  ldouble alpha=sqrt(-1./GG[0][0]);
  for(i=0;i<4;i++)
    {
      A2[i]*=alpha;
    }

   if(verbose>0) print_4vector(A2);
   
  if(verbose>0) getchar();

  return 0; 
}


int
boost2_lab2ff_4vel(ldouble A1[4], ldouble A2[4], ldouble *pp, ldouble gg[][5], ldouble GG[][5], ldouble ucon[4], ldouble ucov[4])
{
  int i,j,k,l;
  ldouble At[4]   ;
  
  int verbose=0;
  
  if(verbose>0) print_4vector(A1);
  
  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_lab2ff_4vel(pp, gg, GG, L, ucon, ucov);
  
  //copying
  for(i=0;i<4;i++)
  {
    At[i]=A1[i];
  }
  
  if(verbose>0) print_tensor(L);
  
  //boosting
  for(i=0;i<4;i++)
  {
    A2[i]=0.;
    for(k=0;k<4;k++)
    {
      A2[i]+=L[i][k]*At[k];
    }
  }
  
  //laps to make it ortonormal
  ldouble alpha=sqrt(-1./GG[0][0]);
  for(i=0;i<4;i++)
  {
    A2[i]*=alpha;
  }
  
  if(verbose>0) print_4vector(A2);
  
  if(verbose>0) getchar();
  
  return 0;
}


/*****************************************************************/
//A^i Lorentz boost from lab to radiation rest frame
/*****************************************************************/

int
boost2_lab2rf(ldouble A1[4],ldouble A2[4],ldouble *pp0,ldouble gg[][5],ldouble GG[][5])
{ 
#ifdef LABRADFLUXES
  my_err("boost2_lab2rf() not working for LABRADFLUXES\n");
#endif
  int i,j,k,l;
  ldouble At[4]   ;
  ldouble pp[NV];
  for(i=0;i<NV;i++)
    pp[i]=pp0[i];

  int verbose=0;

  //artificial and temporary substitution
  ldouble urf[4]={0.,pp[FX],pp[FY],pp[FZ]};
  conv_vels(urf,urf,VELPRIMRAD,VELPRIM,gg,GG);
  pp[2]=urf[1];
  pp[3]=urf[2];
  pp[4]=urf[3];

  if(verbose>0) print_4vector(A1);

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_lab2ff(pp,gg,GG,L);

  //copying and multiplying by lapse to express A1 in ZAMO
  ldouble alpha=sqrt(-1./GG[0][0]);
  for(i=0;i<4;i++)
    {
      At[i]=A1[i]*alpha;
    }
  
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      A2[i]=0.;
      for(k=0;k<4;k++)
	{
	  A2[i]+=L[i][k]*At[k];
	}
    }

  if(verbose>0) print_4vector(A2);

  if(verbose>0) getchar();

  return 0; 
}


/*****************************************************************/
//A^i Lorentz boost from fluid to lab frame
/*****************************************************************/

int
boost2_ff2lab(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5])
{ 
  int i,j,k,l;
  ldouble At[4]   ;

  int verbose=0;

  if(verbose>0) print_4vector(A1);

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_ff2lab(pp,gg,GG,L);

  //copying and ortonormality
  ldouble alpha=sqrt(-1./GG[0][0]);
  for(i=0;i<4;i++)
    {
      At[i]=A1[i]/alpha;
    }
  

  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      A2[i]=0.;
      for(k=0;k<4;k++)
	{
	  A2[i]+=L[i][k]*At[k];
	}
    }


  if(verbose>0) print_4vector(A2);

  if(verbose>0) getchar();

  return 0; 
}


/*****************************************************************/
//multiplies 22 tensor T1 by 21 tensor A
//T2^ij = A^i_k A^j_l T1^kl
/*****************************************************************/

int
multiply22(ldouble T1[][4],ldouble T2[][4],ldouble A[][4])
{
  int i,j,k,l;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=A[i][k]*A[j][l]*Tt[k][l];
		}
	    }
	}
    }
  return 0;
}


/*****************************************************************/
//multiplies 11 tensor T1 by 21 tensor A
//T2_ij = A^k_i A^l_j T1_kl
/*****************************************************************/

int
multiply11(ldouble T1[][4],ldouble T2[][4],ldouble A[][4])
{
  int i,j,k,l;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=A[k][i]*A[l][j]*Tt[k][l];
		}
	    }
	}
    }
  return 0;
}

/*****************************************************************/
//multiplies 2 vector u1 by 21 tensor A
//u2^i = A^i_j u1^j
/*****************************************************************/
int
multiply2(ldouble *u1,ldouble *u2,ldouble A[][4])
{
  int i,j;
  ldouble ut[4];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=A[i][j]*ut[j];
	}
    }

  return 0;
}


/*****************************************************************/
// transforms spatial 3-vectors between coordinates
/*****************************************************************/

int
coco_3vector(ldouble A1[3],ldouble A2[3],int CO1,int CO2,void* ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble xxx[4]={0.,geom->xx,geom->yy,geom->zz};

  int i1,i2;
  if(CO1==CO2) 
    {
      for(i1=0;i1<3;i1++)
	A2[i1]=A1[i1];
      return 0;
    }
  else if(CO1==CYLCOORDS && CO2==MINKCOORDS)
    {
      ldouble ph=geom->zz;
      ldouble At[3];

      At[0]=A1[0]*cos(ph) - A1[2]*sin(ph);  //x
      At[2]=A1[0]*sin(ph) + A1[2]*cos(ph);  //y
      At[1]=A1[1];                          //z-component
      for(i1=0;i1<3;i1++)
	A2[i1]=At[i1];

      return 0;
    }
  else if(CO1==MCYL1COORDS && CO2==MINKCOORDS)
    {
      //MCYL1 modifies only radius and A1[] is in ortonormal basis so the same approach should work
      ldouble ph=geom->zz;
      ldouble At[3];

      At[0]=A1[0]*cos(ph) - A1[2]*sin(ph);  //x
      At[2]=A1[0]*sin(ph) + A1[2]*cos(ph);  //y
      At[1]=A1[1];                          //z-component
      for(i1=0;i1<3;i1++)
	A2[i1]=At[i1];

      return 0;
    }
  else
    my_err("transformation not implemented in coco_3vector\n");
  
  return -1;
}

/*****************************************************************/
//u^i transfromation between coordinates
/*****************************************************************/
		   
int
trans2_coco(ldouble *xx,ldouble *u1,ldouble *u2,int CO1, int CO2)
{
  ldouble dxdx[4][4];
  ldouble xx2[4];
  if(CO1==CO2)
    {
      u2[0]=u1[0];
      u2[1]=u1[1];
      u2[2]=u1[2];
      u2[3]=u1[3];
    }
  else if(((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==SPHCOORDS) ||
	  ((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==SPHCOORDS))
    {
      u2[0]=u1[0];
      u2[1]=u1[1];
      u2[2]=u1[2];
      u2[3]=u1[3];
    }
  else if(CO1==KSCOORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS))
    {
      dxdx_KS2BL(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==KSCOORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==MKS1COORDS && CO2==KSCOORDS)
    {
      dxdx_MKS12KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==MKS2COORDS && CO2==KSCOORDS)
    {
      dxdx_MKS22KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==MKS3COORDS && CO2==KSCOORDS)
    {
      dxdx_MKS32KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==JETCOORDS && CO2==KSCOORDS)
    {
      dxdx_JET2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==TKS3COORDS && CO2==KSCOORDS)
    {
      dxdx_TKS32KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==KSCOORDS && CO2==MKS1COORDS)
    {
      dxdx_KS2MKS1(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==KSCOORDS && CO2==MKS2COORDS)
    {
      dxdx_KS2MKS2(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==KSCOORDS && CO2==MKS3COORDS)
    {
      dxdx_KS2MKS3(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==KSCOORDS && CO2==JETCOORDS)
    {
      dxdx_KS2JET(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==KSCOORDS && CO2==TKS3COORDS)
    {
      dxdx_KS2TKS3(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==MCYL1COORDS && CO2==CYLCOORDS)
    {
      dxdx_MCYL12CYL(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==CYLCOORDS && CO2==MCYL1COORDS)
    {
      dxdx_CYL2MCYL1(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==MSPH1COORDS && (CO2==SPHCOORDS || CO2==SCHWCOORDS || CO2==KERRCOORDS))
    {
      dxdx_MSPH12SPH(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MSPH1COORDS)
    {
      dxdx_SPH2MSPH1(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==MKER1COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      dxdx_MKER12KER(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MKER1COORDS)
    {
      dxdx_KER2MKER1(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if (CO1==MKS1COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      dxdx_MKS12KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2BL(xx2,dxdx);
      multiply2(u2,u2,dxdx);
    }
  else if (CO1==MKS2COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      dxdx_MKS22KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2BL(xx2,dxdx);
      multiply2(u2,u2,dxdx);
    }
  else if (CO1==MKS3COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      dxdx_MKS32KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2BL(xx2,dxdx);
      multiply2(u2,u2,dxdx);
    }
  else if (CO1==JETCOORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      dxdx_JET2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2BL(xx2,dxdx);
      multiply2(u2,u2,dxdx);
    }
  else if (CO1==TKS3COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      dxdx_TKS32KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2BL(xx2,dxdx);
      multiply2(u2,u2,dxdx);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MKS1COORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2MKS1(xx2,dxdx);
      multiply2(u2,u2,dxdx);  
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MKS2COORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2MKS2(xx2,dxdx);
      multiply2(u2,u2,dxdx);  
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MKS3COORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2MKS3(xx2,dxdx);
      multiply2(u2,u2,dxdx);  
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==JETCOORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2JET(xx2,dxdx);
      multiply2(u2,u2,dxdx);  
    }  
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==TKS3COORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2TKS3(xx2,dxdx);
      multiply2(u2,u2,dxdx);  
    }
  else if ((CO1==BLCOORDS || CO1==SPHCOORDS) && CO2==CYLCOORDS)
    {
      dxdx_SPH2CYL(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if ((CO2==BLCOORDS || CO2==SPHCOORDS) && CO1==CYLCOORDS)
    {
      dxdx_CYL2SPH(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }

  else
    {
      printf("transformation not implemented in trans2_coco(): %d -> %d\n",CO1,CO2);
      getch();
    }

  return 0;
}


/*****************************************************************/
//T^ij transformation between coordinates
/*****************************************************************/

int
trans22_coco(ldouble *xx,ldouble T1[][4],ldouble T2[][4],int CO1, int CO2)
{
  ldouble dxdx[4][4];
  ldouble xx2[4];
  if(CO1==CO2)
    {
      int i,j;
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  T2[i][j]=T1[i][j];
    }
  else if(((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==SPHCOORDS) ||
	  ((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==SPHCOORDS))
    {
      int i,j;
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  T2[i][j]=T1[i][j];
    }
  else if(CO1==KSCOORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS))
    {
      dxdx_KS2BL(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==KSCOORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==MKS1COORDS && CO2==KSCOORDS)
    {
      dxdx_MKS12KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==MKS2COORDS && CO2==KSCOORDS)
    {
      dxdx_MKS22KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==KSCOORDS && CO2==MKS1COORDS)
    {
      dxdx_KS2MKS1(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==KSCOORDS && CO2==MKS2COORDS)
    {
      dxdx_KS2MKS2(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==MCYL1COORDS && CO2==CYLCOORDS)
    {
      dxdx_MCYL12CYL(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==CYLCOORDS && CO2==MCYL1COORDS)
    {
      dxdx_CYL2MCYL1(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==MSPH1COORDS && (CO2==SPHCOORDS || CO2==SCHWCOORDS || CO2==KERRCOORDS))
    {
      dxdx_MSPH12SPH(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MSPH1COORDS)
    {
      dxdx_SPH2MSPH1(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==MKER1COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      dxdx_MKER12KER(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MKER1COORDS)
    {
      dxdx_KER2MKER1(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if (CO1==MKS1COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      dxdx_MKS12KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2BL(xx2,dxdx);
      multiply22(T2,T2,dxdx);
    }
  else if (CO1==MKS2COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      dxdx_MKS22KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2BL(xx2,dxdx);
      multiply22(T2,T2,dxdx);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MKS1COORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2MKS1(xx2,dxdx);
      multiply22(T2,T2,dxdx);  
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MKS2COORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2MKS2(xx2,dxdx);
      multiply22(T2,T2,dxdx);  
    }
  else
    {
    printf("transformation not implemented in trans22_coco() %d -> %d\n",CO1,CO2);
    getch();
    }

  return 0;
}


/*****************************************************************/
// T_ij -> T^ij
/*****************************************************************/

int
indices_1122(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5])
{
  int i,j,k,l;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  Tt[i][j]+=T1[k][l]*GG[i][k]*GG[j][l];
		}	  
	    }
	}
    }

   for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=Tt[i][j];
	}
    }

  return 0;
}


/*****************************************************************/
// T^ij -> T_ij
/*****************************************************************/

int
indices_2211(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5])
{
  int i,j,k,l;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  Tt[i][j]+=T1[k][l]*gg[i][k]*gg[j][l];
		}	  
	    }
	}
    }

   for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=Tt[i][j];
	}
    }

  return 0;
}


/*****************************************************************/
// T^i_j -> T^ij
/*****************************************************************/

int
indices_2122(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5])
{
  int i,j,k;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      Tt[i][j]+=T1[i][k]*GG[k][j];
	    }	  
	}
    }

   for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=Tt[i][j];
	}
    }

  return 0;
}

/*****************************************************************/
// T_ij -> T^i_j
/*****************************************************************/

int
indices_1121(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5])
{
  int i,j,k;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      Tt[i][j]+=T1[k][j]*GG[j][i];
	    }	  
	}
    }

   for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=Tt[i][j];
	}
    }

  return 0;
}


/*****************************************************************/
// T^ij -> T^i_j
/*****************************************************************/

int
indices_2221(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5])
{
  int i;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
  {
    int j;
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
    for(j=0;j<4;j++)
    {
      Tt[i][j]=0.;
    }
  }

  for(i=0;i<4;i++)
  {
    int j;
    for(j=0;j<4;j++)
    {
      int k;
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
      for(k=0;k<4;k++)
      {
        Tt[i][j]+=T1[i][k]*gg[k][j];
      }
    }
  }

  for(i=0;i<4;i++)
  {
    int j;
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
    for(j=0;j<4;j++)
    {
      T2[i][j]=Tt[i][j];
    }
  }

  return 0;
}

/*****************************************************************/
// A_i -> A^j
/*****************************************************************/

int
indices_12(ldouble A1[4],ldouble A2[4],ldouble GG[][5])
{
  int i;

#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
  for(i = 0; i < 4; i++)
  {
    A2[i] = 0.;
  }

#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
  for(i = 0; i < 4; i++)
  {
    int k;
    for(k = 0; k < 4; k++)
    {
      A2[i] += A1[k] * GG[i][k];
    }
  }

  return 0;
}


/*****************************************************************/
// A^i -> A^_j
/*****************************************************************/
// 7/8/17, Ramesh: The code breaks if we try to calculate A2 directly, as in indices_12. It is necessary to compute At first and then copy it to A2.

int
indices_21(ldouble A1[4],ldouble A2[4],ldouble gg[][5])
{
  int i;
  ldouble At[4];

#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
  for(i=0;i<4;i++)
  {
    At[i]=0.;
  }

  for(i=0;i<4;i++)
  {
    int j;
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
    for(j = 0; j < 4; j++)
    {
      At[i] += A1[j] * gg[i][j];
    }
  }

#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
  for (i=0; i<4; i++)
  {
    A2[i] = At[i];
  }
  
  return 0;
}

/*****************************************************************/
// precomputed MYCOORDS -> OUTCOORDS transformation
// ONLY WORKS at the cell center
// which = 0 : MYCOORDS -> OUTCOORDS
// which = 1 : OUTCOORDS -> MYCOORDS
/*****************************************************************/
int
trans_pall_coco_my2out(ldouble *pp1, ldouble *pp2, void* ggg1, void* ggg2) 
{
  trans_pmhd_coco_my2out(pp1, pp2, ggg1,ggg2);
#ifdef RADIATION
  trans_prad_coco_my2out(pp1, pp2, ggg1,ggg2);
#endif
  return 0;
}

int
trans_pall_coco_out2my(ldouble *pp1, ldouble *pp2, void* ggg1, void* ggg2) 
{
  trans_pmhd_coco_out2my(pp1, pp2, ggg1,ggg2);
#ifdef RADIATION
  trans_prad_coco_out2my(pp1, pp2, ggg1,ggg2);
#endif
  return 0;
}

int
trans_pmhd_coco_my2out(ldouble *ppin, ldouble *ppout, void* ggg1, void* ggg2)
{
  int out = trans_pmhd_coco_precompute(ppin, ppout, ggg1, ggg2, 0);
  return out;
}

int
trans_pmhd_coco_out2my(ldouble *ppin, ldouble *ppout, void* ggg1, void* ggg2)
{
  int out = trans_pmhd_coco_precompute(ppin, ppout, ggg1, ggg2, 1);
  return out;
}

int
trans_pmhd_coco_precompute(ldouble *ppin, ldouble *ppout, void* ggg1,void* ggg2, int which)
{
  struct geometry *geom1
  = (struct geometry *) ggg1;
  struct geometry *geom2
  = (struct geometry *) ggg2;
  
  int i,iv;
  ldouble pp1[NV],pp2[NV];
  for(i=0;i<NV;i++)
  {
    ppout[i]=ppin[i];
    pp1[i]=ppin[i];
    pp2[i]=ppout[i];
  }
  
  if(OUTCOORDS==MYCOORDS)
  {    
    for (i = 0; i < 5; i++)
    {
      pp2[i] = pp1[i];
    }
  }
  else
  {
    pp2[0]=pp1[0];
    pp2[1]=pp1[1];

    // four velocity ucon
    ldouble ucon[4], ucov[4],uconback[4];    
    calc_ucon_ucov_from_prims(pp1, geom1, ucon, ucov);

    
#ifdef MAGNFIELD
    ldouble bcon[4],Bcon[4];
    
    //magnetic field 4-vector bcon
    calc_bcon_4vel(pp1, ucon, ucov, bcon);
#endif

    for(iv=0;iv<4;iv++)
      uconback[iv]=ucon[iv];
    
    //convert ucon (and bcon) to OUTCOORDS using precomputed matrices
    trans2_coco_precompute(ucon, ucon, geom1->ix, geom1->iy, geom1->iz, which);
        
#ifdef MAGNFIELD
    trans2_coco_precompute(bcon, bcon, geom1->ix, geom1->iy, geom1->iz, which);
#endif
    
    //to VELPRIM
    conv_vels_ut(ucon, ucon,VEL4,VELPRIM, geom2->gg, geom2->GG);
    
    pp2[2]=ucon[1];
    pp2[3]=ucon[2];
    pp2[4]=ucon[3];
    
#ifdef MAGNFIELD

    // back to primitive B^i
    calc_Bcon_prim(pp2,bcon,Bcon,geom2);
        
    for (i = 0; i < 3; i++)
    {
      pp2[B1+i] = Bcon[1+i];
    }
#endif
  }
  
  for(i=0;i<NVMHD;i++)
  {
    ppout[i]=pp2[i];
  }
  
  return 0;
}

int
trans_prad_coco_my2out(ldouble *ppin, ldouble *ppout, void* ggg1, void* ggg2)
{
  int out = trans_prad_coco_precompute(ppin, ppout, ggg1, ggg2, 0);
  return out;
}

int
trans_prad_coco_out2my(ldouble *ppin, ldouble *ppout, void* ggg1, void* ggg2)
{
  int out = trans_prad_coco_precompute(ppin, ppout, ggg1, ggg2, 1);
  return out;
}

int
trans_prad_coco_precompute(ldouble *ppin, ldouble *ppout, void* ggg1, void* ggg2, int which)
{
  
  struct geometry *geom1
    = (struct geometry *) ggg1;
  struct geometry *geom2
    = (struct geometry *) ggg2;

  int i;
  
  ldouble pp1[NV],pp2[NV];
  for(i=0;i<NV;i++) 
    {
      ppout[i]=ppin[i];
      pp1[i]=ppin[i];
      pp2[i]=ppout[i];
    }      
#ifdef RADIATION 
  if(OUTCOORDS==MYCOORDS)
    {
      for(i=0;i<4;i++)
	pp2[EE0+i]=pp1[EE0+i];
     }
  else
    {
      //Erf unchanged
      pp2[EE0]=pp1[EE0];

      //velocity in CO1
      ldouble ucon[4];
      ucon[0]=0;
      ucon[1]=pp1[FX0];
      ucon[2]=pp1[FY0];
      ucon[3]=pp1[FZ0];

      conv_vels(ucon,ucon,VELPRIMRAD,VEL4,geom1->gg,geom1->GG);
      
      //converting to CO2
      trans2_coco_precompute(ucon, ucon, geom1->ix, geom1->iy, geom1->iz, which);

      //to VELPRIM
      conv_vels_ut(ucon,ucon,VEL4,VELPRIMRAD,geom2->gg,geom2->GG);

      pp2[FX0]=ucon[1]; 
      pp2[FY0]=ucon[2];
      pp2[FZ0]=ucon[3];

    }
#endif //RADIATION
  for(i=NVMHD;i<NV;i++)     
    {
      ppout[i]=pp2[i];
    }      
 
  return 0;
}

int
trans2_coco_my2out(ldouble *u1, ldouble *u2, int ix, int iy, int iz)
{
  int out = trans2_coco_precompute(u1, u2, ix, iy, iz, 0);
  return out;

}

int
trans2_coco_out2my(ldouble *u1, ldouble *u2, int ix, int iy, int iz)
{
  int out = trans2_coco_precompute(u1, u2, ix, iy, iz, 1);
  return out;

}

int
trans2_coco_precompute(ldouble *u1, ldouble *u2, int ix, int iy, int iz, int which)
{
  ldouble dxdx[4][4];
  int i,j;

  if(OUTCOORDS==MYCOORDS)
  {
    for(i=0;i<4;i++)
      u2[i] = u1[i];
  }
  else
  {
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
      {
	if(which==0)
	{
          dxdx[i][j] = get_dxdx(dxdx_my2out,i,j,ix,iy,iz);

	}
	else if(which==1)
	{
          dxdx[i][j] = get_dxdx(dxdx_out2my,i,j,ix,iy,iz);
        }
	else
	{
          printf("In trans2_coco_precompute, which flag must be 0 or 1!\n");
	  getch();
	}
      }
    multiply2(u1,u2,dxdx);
  }
  
  return 0;
}

int
trans22_coco_my2out(ldouble T1[][4], ldouble T2[][4], int ix, int iy, int iz)
{
  int out = trans22_coco_precompute(T1, T2, ix, iy, iz, 0);
  return out;
}

int
trans22_coco_out2my(ldouble T1[][4], ldouble T2[][4], int ix, int iy, int iz)
{
  int out = trans22_coco_precompute(T1, T2, ix, iy, iz, 1);
  return out;
}

int
trans22_coco_precompute(ldouble T1[][4], ldouble T2[][4], int ix, int iy, int iz, int which)
{
  ldouble dxdx[4][4];
  int i,j;

  if(OUTCOORDS==MYCOORDS)
  {
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        T2[i][j] = T1[i][j];
  }
  else
  {
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
      {
	if(which==0)
	{
          dxdx[i][j] = get_dxdx(dxdx_my2out,i,j,ix,iy,iz);
        }
	else if(which==1)
	{
          dxdx[i][j] = get_dxdx(dxdx_out2my,i,j,ix,iy,iz);
        }
	else
	{
          printf("In trans22_coco_precompute, which flag must be 0 or 1!\n");
	  getch();
	}
      }
    
    multiply22(T1,T2,dxdx);
  }
  
  return 0;
}


/*********************************************************************************************************************/
/****** radiative ff primitives (\ehat,\hat urf^i) -> primitives in lab frame  ***************************************/
/****** has changed! previously took (Ehat, F^i), no longer! that was not consistent and did not tranform well *******/
/*********************************************************************************************************************/

int prad_ff2lab(ldouble *pp1, ldouble *pp2, void* ggg)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  int i,j;

  //print_primitives(pp1);

  ldouble (*gg)[5],(*GG)[5],gdetu;
  gg=geom->gg;
  GG=geom->GG;
  ldouble tlo[4][4];
  //calc_tetrades(geom->gg,tup,tlo,MYCOORDS);
  //approximate:
  DLOOP(i,j) tlo[i][j]=0.;
  DLOOPA(i) tlo[i][i]=1./sqrt((gg[i][i]));
  tlo[0][0]=1.;

  gdetu=geom->gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  ldouble Rij[4][4];

  int verbose=0;
 
  calc_Rij(pp1,geom,Rij);  
  boost22_ff2lab_with_alpha(Rij, Rij, pp1, gg, GG, geom->alpha);

  indices_2221(Rij,Rij,gg);  

  for(i=0;i<NVMHD;i++)
    pp2[i]=pp1[i];

  //temporarily store conserved in uu[]
 
  ldouble uu[NV];

  uu[EE0]=gdetu*Rij[0][0];
  uu[FX0]=gdetu*Rij[0][1];
  uu[FY0]=gdetu*Rij[0][2];
  uu[FZ0]=gdetu*Rij[0][3];
 
  #ifdef EVOLVEPHOTONNUMBER
  ldouble nphff=pp1[NF0];
  #endif 

  //convert to real primitives
  int corrected;
  u2p_rad(uu,pp2,geom,&corrected);

  #ifdef EVOLVEPHOTONNUMBER
  //velocities of the frames
  ldouble ut[4];ut[1]=pp2[VX];ut[2]=pp2[VY];ut[3]=pp2[VZ];
  ldouble uffcov[4],uffcon[4];
  conv_vels_both(ut,uffcon,uffcov,VELPRIM,VEL4,gg,GG);
  ldouble urfcov[4],urfcon[4];
  ut[1]=pp2[FX0];ut[2]=pp2[FY0];ut[3]=pp2[FZ0];
  conv_vels_both(ut,urfcon,urfcov,VELPRIMRAD,VEL4,gg,GG);

  ldouble relgamma = urfcon[0]*uffcov[0] + urfcon[1]*uffcov[1] +urfcon[2]*uffcov[2] +urfcon[3]*uffcov[3]; 
  ldouble nphrf = -nphff/relgamma;

  pp2[NF0]=nphrf;
  #endif

  //print_primitives(pp2);getch();

  return 0;
} 


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost from fluid frame to lab

int
boost22_ff2lab_with_alpha(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5], ldouble alpha)
{
  int i,j,k,l;
  ldouble Tt[4][4];
  
  int verbose=0;
  
  if(verbose>0) print_tensor(T1);
  
  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_ff2lab(pp,gg,GG,L);
  
  //copying
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
  //correcting for ortonormality
  //ldouble alpha=sqrt(-1./GG[0][0]);
  for(i=0;i<4;i++)
    {
      T1[i][0]/=alpha;
      T1[0][i]/=alpha;
    }
  
  if(verbose>0) print_tensor(L);
  
  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    }
  
  /*
   //dividing by lapse to express T2 in no-frame
   ldouble alpha=sqrt(-1./GG[0][0]);
   for(i=0;i<4;i++)
   {
   for(j=0;j<4;j++)
   {
     T2[i][j]=T2[i][j]/alpha;
   }
   }
  */
  
  if(verbose>0) print_tensor(T2);
  
  if(verbose>0) getchar();
  
  return 0;
}
