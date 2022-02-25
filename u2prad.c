//KORAL - u2p.c
//radiative routines related to u2p inversion

#include "ko.h"



// get's gamma assuming fixed E rather than using original R^t_t that we assume is flawed near floor regions.  We want to preserve R^t_i (i.e conserved momentum)
static int get_m1closure_gammarel2_cold(int verbose,
					void *ggg, 
					FTYPE *Avcon, 
					FTYPE *Avcov, 
					FTYPE *gammarel2return, 
					FTYPE *deltareturn, 
					FTYPE *numeratorreturn, 
					FTYPE *divisorreturn, 
					FTYPE *Erfreturn, 
					FTYPE *urfconrel)
{
   struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  FTYPE gamma2,gammarel2,delta;
  FTYPE Erf;
  FTYPE alpha=geom->alpha;
  int jj;

  //static
  FTYPE gctt, gn11, gn12,  gn13,  gn14,  gn22,  gn23,  gn24,  gn33,  gn34,  gn44,  Rtt,  Rtx,  Rty,  Rtz,  Rdtt,  Rdtx,  Rdty,  Rdtz;
  gn11=GG[0][0];
  gn12=GG[0][1];
  gn13=GG[0][2];
  gn14=GG[0][3];
  gn22=GG[1][1];
  gn23=GG[1][2];
  gn24=GG[1][3];
  gn33=GG[2][2];
  gn34=GG[2][3];
  gn44=GG[3][3]; 

  Rtt=Avcon[0];
  Rtx=Avcon[1];
  Rty=Avcon[2];
  Rtz=Avcon[3];

  Rdtt=Avcov[0];
  Rdtx=Avcov[1];
  Rdty=Avcov[2];
  Rdtz=Avcov[3];

  // choose gamma
  if(gammarel2return==NULL){
    FTYPE gammamax=GAMMAMAXRAD;
    gammarel2=gammamax*gammamax;
  }
  else gammarel2=*gammarel2return; // feed in desired gammarel2

  FTYPE utsq=gammarel2/(alpha*alpha);

  FTYPE Avcovorig[4],Avconorig[4];
  DLOOPA(jj) Avcovorig[jj]=Avcov[jj];
  DLOOPA(jj) Avconorig[jj]=Avcon[jj];

  // get new Avcov[0]=R^t_t

  // NOTEMARK: Note that Sqrt() is only ever negative when gammarel2<0, so never has to be concern.
  Avcov[0]=(0.25*(-4.*(gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*utsq*(gn11 + utsq) + 
                        Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
                              1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                                       gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2))))/(gn11*utsq*(gn11 + utsq));

  Erf=(0.75*Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
                           1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                                    gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2)))/(utsq*(gn11 + utsq)*(gn11 + 4.*utsq));

  if(verbose) printf("NOR SOL: Avcov0new=%g Avcov0old=%g Erf=%g :: %g %g %g\n",Avcov[0],Avcovorig[0],Erf,Rtx,Rty,Rtz);

  //modify Avcon
  indices_12(Avcov,Avcon,GG);
  if(verbose) DLOOPA(jj) printf("jj=%d Avconorig=%g Avconnew=%g\n",jj,Avconorig[jj],Avcon[jj]);

 
  delta=0; // not yet

  *gammarel2return=gammarel2;
  *deltareturn=delta;

  // get good relative velocity
  FTYPE gammarel=sqrt(gammarel2);

  // get relative 4-velocity
  if(Erf>0.0) SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*GG[0][jj]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
  else SLOOPA(jj) urfconrel[jj] = 0.0;

  *Erfreturn=Erf; // pass back new Erf to pointer

  return(0);
}




// get's gamma^2 for lab-frame gamma  using Rd and gcon
static int get_m1closure_gammarel2(int verbose,void *ggg, ldouble *Avcon, ldouble *Avcov, ldouble *gammarel2return,ldouble *deltareturn, ldouble *numeratorreturn, ldouble *divisorreturn)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble gamma2,gammarel2,delta,numerator,divisor;
  ldouble gamma2a,gamma2b;

  // mathematica solution that avoids catastrophic cancellation when Rtt very small (otherwise above gives gamma2=1/2 oddly when gamma2=1) -- otherwise same as above
  // well, then had problems for R~1E-14 for some reason when near BH.  Couldn't quickly figure out, so use no replacement of gv11.
  // see u2p_inversion.nb
  //static
  ldouble gctt, gn11, gn12,  gn13,  gn14,  gn22,  gn23,  gn24,  gn33,  gn34,  gn44,  Rtt,  Rtx,  Rty,  Rtz,  Rdtt,  Rdtx,  Rdty,  Rdtz;
  gn11=GG[0][0];
  gn12=GG[0][1];
  gn13=GG[0][2];
  gn14=GG[0][3];
  gn22=GG[1][1];
  gn23=GG[1][2];
  gn24=GG[1][3];
  gn33=GG[2][2];
  gn34=GG[2][3];
  gn44=GG[3][3];

  Rtt=Avcon[0];
  Rtx=Avcon[1];
  Rty=Avcon[2];
  Rtz=Avcon[3];

  Rdtt=Avcov[0];
  Rdtx=Avcov[1];
  Rdty=Avcov[2];
  Rdtz=Avcov[3];

  gamma2a=(-0.25*(2.*Power(gn11,2)*Power(Rdtt,2) + (gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*
        (gn12*Rdtx + gn13*Rdty + gn14*Rdtz + Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
            gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
               8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2)))) + 
       gn11*(4.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 
          2.*gn34*Rdty*Rdtz + gn44*Power(Rdtz,2) + Rdtt*
           (4.*gn13*Rdty + 4.*gn14*Rdtz + Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
               gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                  8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2)))))))/
   (gn11*Power(Rdtt,2) + 2.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn13*Rdtt*Rdty + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 
    2.*(gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2));

  if(verbose>1)     printf("gamma2a: %e\n",gamma2a);

  if( gamma2a<GAMMASMALLLIMIT || !isfinite(gamma2a) ){
    gamma2b=(0.25*(-2.*Power(gn11,2)*Power(Rdtt,2) - 1.*gn11*(4.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 
                                                              Rdty*(4.*gn13*Rdtt + 2.*gn23*Rdtx + gn33*Rdty) + 2.*(2.*gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2)) + 
                   gn11*Rdtt*Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
                                  gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                                        8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2))) + 
                   (gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*(-1.*gn12*Rdtx - 1.*gn13*Rdty - 1.*gn14*Rdtz + 
                                                        Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
                                                             gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                                                                   8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2))))))/
      (gn11*Power(Rdtt,2) + 2.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn13*Rdtt*Rdty + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 
       2.*(gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2));
    gamma2=gamma2b;

    if(verbose)     printf("gamma2b: %e\n",gamma2b);
  }
  else{
    // choose
    gamma2=gamma2a;
  }

  //  if(isnan(gamma2a) && !isnan(gamma2b))
  //    gamma2=gamma2b;

  //  if(!isnan(gamma2a) && isnan(gamma2b))
  //    gamma2=gamma2a;

  ////////////////////////
  //
  //cap on u^t
  //
  ///////////////////////
  ldouble alpha=geom->alpha;


  // get relative 4-velocity, that is always >=1 even in GR
  gammarel2 = gamma2*alpha*alpha;

  //to consider it a failure - leads to instability
  //if(isnan(gammarel2)) return -1;
  
  // check for machine error away from 1.0 that happens sometimes
  //if(isnan(gammarel2) || (gammarel2>GAMMASMALLLIMIT && gammarel2<1.0)) - WTF?
  
  if((gammarel2>GAMMASMALLLIMIT && gammarel2<1.0))
  {
    if(verbose) printf("Hit machine error of gammarel2=%27.20g fixed to be 1.0\n",gammarel2);

    gammarel2=1.0;
  }

  /*
  if(isnan(gammarel2))
    {
      printf("nan in get_m1closure_gammarel2\n");
      printf("%e %e %e\n",gamma2,gamma2a,gamma2b);
      print_4vector(Avcon);
      print_4vector(Avcov);
      //      getchar();
    }
  */

  *gammarel2return=gammarel2; 
  *deltareturn=delta=0;
  *numeratorreturn=numerator=0;
  *divisorreturn=divisor=0;
  return(0);
}


// get Erf
static int get_m1closure_Erf(void *ggg, ldouble *Avcon, ldouble gammarel2, ldouble *Erfreturn)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble alpha=geom->alpha;

  ////////////
  //
  // get initial attempt for Erf
  // If delta<0, then gammarel2=nan and Erf<RADLIMIT check below will fail as good.
  //
  ////////////
  *Erfreturn = 3.*Avcon[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

  /*
  if(isnan(*Erfreturn)) 
    {
      printf("nan in get_m1closure_Erf\n %d %d %d %e %e\n",geom->ix,geom->iy,geom->iz,*Erfreturn,gammarel2);
      print_4vector(Avcon);
      getchar();
    }
  */

  return(0);
}

// get contravariant relative 4-velocity in lab frame
static int get_m1closure_urfconrel(int verbose, 
				   void *ggg, 
				   ldouble *pp, 
				   ldouble *Avcon, 
				   ldouble *Avcov, 
				   ldouble gammarel2, 
				   ldouble delta, 
				   ldouble numerator, 
				   ldouble divisor, 
				   ldouble *Erfreturn,
				   ldouble *urfconrel, 
				   int *corflag)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble alpha=geom->alpha;

  ldouble Erf=*Erfreturn; // get initial Erf
  ldouble Erfslow,Erffast;
  ldouble gammamax=GAMMAMAXRAD; 
  int ii,jj,kk,usingfast;

  //////////////////////
  //
  // Fix-up inversion if problem with gamma (i.e. velocity) or energy density in radiation rest-frame (i.e. Erf)
  //
  //////////////////////

  // NOTE: gammarel2 just below 1.0 already fixed to be =1.0
  int nonfailure=gammarel2>=1.0 && Erf>ERADFLOOR && gammarel2<=gammamax*gammamax/GAMMASMALLLIMIT/GAMMASMALLLIMIT;
  // falilure1 : gammarel2 normal, but already Erf<ERADFLOOR (note for M1 that gammarel2>=1/4 for any reasonable chance for correct non-zero Erf
  int failure1=Avcon[0]<0.0 || (gammarel2>0.0 && gammarel2<=0.25 && delta>=0.0 && divisor!=0.0) || numerator==0.0 || (gammarel2>=1.0 && delta>=0.0 && divisor!=0.0 && Erf<ERADFLOOR);
  //int failure1=Avcon[0]<0.0 || gammarel2>0.0 && gammarel2<=0.25 && delta>=0.0 && divisor!=0.0 || numerator==0.0 || gammarel2>=1.0 && delta>=0.0 && divisor!=0.0 && Erf<ERADFLOOR;
  // gamma probably around 1
  int failure2=gammarel2<1.0 && gammarel2>0.0 && delta>=0.0;
  // i.e. all else, so not really used below.
  int failure3=(gammarel2>gammamax*gammamax && Erf>=ERADFLOOR) || gammarel2<0.0 || delta<0.  || (divisor==0.0 && numerator==0.0) || (divisor==0.0 && numerator!=0.0);
  //int failure3=gammarel2>gammamax*gammamax && Erf>=ERADFLOOR || gammarel2<0.0 || delta<0.  || divisor==0.0 && numerator==0.0 || divisor==0.0 && numerator!=0.0;

  // any failure
  int failure=!nonfailure || isinf(gammarel2) || isinf(Erf);

  if(failure && (failure1==0 && failure2==0 && failure3==0)){
    printf("Undetected failure, now considered\n");
  }


  if(nonfailure){
    // get good relative velocity
    ldouble gammarel=sqrt(gammarel2);
 
    for(ii=0;ii<4;ii++)
      {	  
	urfconrel[ii] = alpha * (Avcon[ii] + 1./3.*Erf*GG[0][ii]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
      }

    *Erfreturn=Erf; // pass back new Erf to pointer
    *corflag=0;
    return(0);
  }
  else{
    if(verbose) 
      {
	printf("failure: %d %d %d %d %d\n",nonfailure,failure1,failure2,failure3,failure);
	printf("%e %e %e\n",Erf,gammarel2,Avcov[0]);
	print_4vector(Avcon);
	print_4vector(Avcov);
	if(isnan(Erf))	getchar();
      }
    // get \gammarel=1 case
    ldouble gammarel2slow=pow(1.0+10.0*NUMEPSILON,2.0);
    ldouble Avconslow[4],Avcovslow[4],urfconrelslow[4];
    for(jj=0;jj<4;jj++)
      {
	Avconslow[jj]=Avcon[jj];
	Avcovslow[jj]=Avcov[jj];
      }
    Erfslow=Erf;
    get_m1closure_gammarel2_cold(verbose,ggg,Avconslow,Avcovslow,&gammarel2slow,&delta,&numerator,&divisor,&Erfslow,urfconrelslow);

    // get \gammarel=gammamax case
    ldouble gammarel2fast=gammamax*gammamax;
    ldouble Avconfast[4],Avcovfast[4],urfconrelfast[4];
    for(jj=0;jj<4;jj++)
      {
	Avconfast[jj]=Avcon[jj];
	Avcovfast[jj]=Avcov[jj];
      }
    Erffast=Erf;
    get_m1closure_gammarel2_cold(verbose,ggg,Avconfast,Avcovfast,&gammarel2fast,&delta,&numerator,&divisor,&Erffast,urfconrelfast);

    usingfast=1;
    // choose by which Avcov[0] is closest to original&&
    //    if( fabs(Avcovslow[0]-Avcov[0])>fabs(Avcovfast[0]-Avcov[0])){
    if( fabs(Avcovslow[0]-Avcov[0])>fabs(Avcovfast[0]-Avcov[0])){
      usingfast=1;
      for(jj=0;jj<4;jj++)
	{
	  Avcon[jj]=Avconfast[jj];
	  Avcov[jj]=Avcovfast[jj];
	  urfconrel[jj]=urfconrelfast[jj];
	}
      gammarel2=gammarel2fast;
      Erf=Erffast;
    }
    else{
      usingfast=0;
      for(jj=0;jj<4;jj++)
	{
	  Avcon[jj]=Avconslow[jj];
	  Avcov[jj]=Avcovslow[jj];
	  urfconrel[jj]=urfconrelslow[jj];
	}
      gammarel2=gammarel2slow;
      Erf=Erfslow;
    }    

    // report
    *corflag=1;
    if(verbose) printf("CASEGEN: gammarel>gammamax (cold, usingfast=%d): gammarel2=%g Erf=%g : i=%d j=%d k=%d\n",usingfast,gammarel2,Erf,geom->ix,geom->iy,geom->iz);
  }


  *Erfreturn=Erf; // pass back new Erf to pointer

  if(verbose) 
      {
	printf("end: %e %e %e %e\n",Erf,*Erfreturn,Erfslow,Erffast);
      }

  if(!isfinite(Erf) || !isfinite(gammarel2) || !isfinite(urfconrel[1])|| !isfinite(urfconrel[2])|| !isfinite(urfconrel[3]))
    {
      if(verbose)      
	{
	  int gix,giy,giz;
	  mpi_local2globalidx(geom->ix,geom->iy,geom->iz,&gix,&giy,&giz);
	  printf("------------\nJONNAN: %e %e %e %e\n",Erf,*Erfreturn,Erfslow,Erffast);
	  if(usingfast) printf("JONNAN: usingfast==1\n"); else printf("JONNAN: usingfast==0\n");
	  
	  printf("JONNAN: ijk=%d %d %d :  %g %g : %g %g %g %g : %d %d %d %d : %g %g %g %g\n",gix,giy,giz,Erf,gammarel2,urfconrel[0],urfconrel[1],urfconrel[2],urfconrel[3],failure1,failure2,failure3,failure,Avcon[0],Avcon[1],Avcon[2],Avcon[3]);
	}
    

      return -1;
  }

  return(0);
}



//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver for radiation
//uses M1 closure in arbitrary frame/metric
//radiative primitives: (E,\tilde u^i)
//  E - radiative energy density in the rad.rest frame
//  u^i - relative velocity of the rad.rest frame
//takes conserved R^t_mu in uu
//**********************************************************************
//**********************************************************************
int
u2p_rad_urf(ldouble *uu, ldouble *pp,void* ggg, int *corrected)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int verbose=0;
  
  //whether primitives corrected for caps, floors etc. - if so, conserved will be updated 
  *corrected=0;

  ldouble Rij[4][4];
  ldouble urfcon[4],urfcov[4],Erf;
  //conserved - R^t_mu
  ldouble Avcov[4]={uu[EE]/gdetu,uu[FX]/gdetu,uu[FY]/gdetu,uu[FZ]/gdetu};
  ldouble Avcon[4];
  //indices up - R^tmu
  indices_12(Avcov,Avcon,GG);

  ldouble gammarel2,delta,numerator,divisor;

  // get \gamma^2 for relative 4-velocity
  if(get_m1closure_gammarel2(verbose,ggg,Avcon,Avcov,&gammarel2,&delta,&numerator,&divisor)<0)
    {
      //printf("get_m1closure_gammarel2 failed at %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK);
      return -1;
    }

  // get E in radiation frame
  get_m1closure_Erf(ggg,Avcon,gammarel2,&Erf);

  // get relative 4-velocity
  if(get_m1closure_urfconrel(verbose,ggg,pp,Avcon,Avcov,gammarel2,delta,numerator,divisor,&Erf,urfcon,corrected)<0)
    {
      //printf("get_m1closure_urfconrel failed at %d %d %d\n",geom->ix+TOI,geom->iy+TOJ,geom->iz+TOK);
      return -1;
    }

  conv_vels(urfcon,urfcon,VELR,VELPRIMRAD,gg,GG);

  //new primitives
  pp[EE]=Erf;
  pp[FX]=urfcon[1];
  pp[FY]=urfcon[2];
  pp[FZ]=urfcon[3];


  return 0;
}


int
u2p_rad(ldouble *uu, ldouble *pp, void *ggg, int *corrected)
{
  //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
  *corrected=0;

  int u2pret=u2p_rad_urf(uu,pp,ggg,corrected);

  #ifdef EVOLVEPHOTONNUMBER
  struct geometry *geom
   = (struct geometry *) ggg;
  ldouble gdetu=geom->gdet;
  #if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
  #endif
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX];
  urfcon[2]=pp[FY];
  urfcon[3]=pp[FZ];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom->gg,geom->GG);
 
  pp[NF]=uu[NF]/urfcon[0]/gdetu;
  #endif

  //M1
  return u2pret;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical conserved to primitives solver for radiation
//works in ortonormal fluid frame
//used e.g. for not-frame-invariant  Eddington apr. 
//solves in 4dimensions using frame boosts etc.
int f_u2prad_num(ldouble *uu,ldouble *pp, void* ggg,ldouble *f)
{
  /*struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

  ldouble Rij[4][4];

  calc_Rij_M1_ff(pp,Rij);  
  trans22_on2cc(Rij,Rij,tlo);  
  boost22_ff2lab(Rij,Rij,pp,gg,GG); 
  indices_2221(Rij,Rij,gg);  

  f[0]=-Rij[0][0]+uu[EE0];
  f[1]=-Rij[0][1]+uu[FX0];
  f[2]=-Rij[0][2]+uu[FY0];
  f[3]=-Rij[0][3]+uu[FZ0];
*/
  return 0;
} 

int
print_state_u2prad_num (int iter, ldouble *x, ldouble *f)
{
  printf ("iter = %3d x = % .3e % .3e % .3e % .3e "
	  "f(x) = % .3e % .3e % .3e % .3e\n",
	  iter,
	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
  return 0;
}

int
u2p_rad_onff(ldouble *uu, ldouble *pp, void* ggg, int *corrected)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble pp0[NV],pporg[NV];
  ldouble J[4][4],iJ[4][4];
  ldouble x[4],f1[4],f2[4],f3[4];
  int i,j,k,iter=0;

  ldouble EPS = 1.e-6;
  ldouble CONV = 1.e-6;

  int verbose=0;

  for(i=EE0;i<NV;i++)
    {
      pporg[i]=pp[i];
    }

  //converting radiative primitives to fluid frame ortonormal
  //ad_lab2ff(pp,pp,geom);

  if(verbose!=0)   print_Nvector(uu,NV);
  do
    {
      iter++;
      for(i=EE0;i<NV;i++)
	{
	  pp0[i]=pp[i];
	}

      //valueas at zero state
      f_u2prad_num(uu,pp,geom,f1);
 
      //calculating approximate Jacobian
      for(j=0;j<4;j++)
	{
	  pp[j+EE0]=pp[j+EE0]+EPS*pp[EE0];
	    
	  if(verbose>0)    print_Nvector(pp,NV);
 	  f_u2prad_num(uu,pp,geom,f2);
	  if(verbose>0)    print_state_u2prad_num (iter,x,f2); 
     
	  for(i=0;i<4;i++)
	    {
	      J[i][j]=(f2[i] - f1[i])/(EPS*pp[EE0]);
	    }

	  pp[j+EE0]=pp0[j+EE0];
	}
	

      //inversion
      inverse_44matrix(J,iJ);

      //updating x
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+EE0];
	}

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      x[i]-=iJ[i][j]*f1[j];
	    }
	}
      if(verbose>0)    print_state_u2prad_num (iter,x,f1); 

      for(i=0;i<4;i++)
	{
	  pp[i+EE0]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+EE0]-pp0[i+EE0]);
	  f3[i]=fabs(f3[i]/pp0[EE0]);
	}

      if(f3[0]<CONV && f3[1]<CONV && f3[2]<CONV && f3[3]<CONV)
	break;

      if(iter>20)
	{
	  printf("iter exceeded in u2prad_num() %d %d %d\n",geom->ix,geom->iy,geom->iz);getchar();

	  pp[EE0]=pporg[EE0];
	  pp[FX0]=pp[FY0]=pp[FZ0]=0.;
	  
	  *corrected=1;
	  return -1;

	  break;
	}
     
    }
  while(1);
  
  if(pp[EE0]<EEFLOOR) 
    {
      printf("enegative u2prad() %d %d %d\n",geom->ix,geom->iy,geom->iz); getchar();
      pp[EE0]=EEFLOOR;
      pp[7]=pp[8]=pp[9]=0.;
      *corrected=1;
      return -1;
    }
  
  //converting to lab primitives
  //no - in EDD_APR I use fluid frame fluxes as primitives
  //prad_ff2lab(pp,pp,geom);
  
  if(verbose!=0)   {print_Nvector(pp,NV);}
  if(verbose>0)   {printf("----\n");}//getchar();}

  *corrected=0;
  return 0;

}

/// Debora - taken from koral_lite to count in photon number
int
check_floors_rad(ldouble *pp, int whichvel,void *ggg)
{
  //skip floors for some time
  //return 0;

  int verbose=0;
  int ret=0;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;


#ifdef RADIATION  
  ldouble pp2[NV];
  int iv;

  //ff rad energy density
  ldouble Rtt,Ehat,ugas[4],Eratio;
  calc_ff_Rtt(pp,&Rtt,ugas,geom);
  Ehat=-Rtt;
  Eratio=pp[EE0]/Ehat; //ratio of energy densities in two frames

  //absolute rad-frame EE:
  if(pp[EE0]<ERADFLOOR) 
    {
      if(verbose) printf("rhd_floors CASE R0 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[0],pp[EE0]);
      pp[EE0]=ERADFLOOR;
      ret=-1;
     }

#ifndef SKIPRADSOURCE

  //Ehat/rho ratios 
  if(Ehat<EERHORATIOMIN*pp[0]) 
    {
      if(verbose) printf("hd_floors CASE R2 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[0],pp[EE0]);
      pp[EE0]=Eratio*EERHORATIOMIN*pp[0];
      ret=-1;
    }

  if(Ehat>EERHORATIOMAX*pp[0]) 
    {
      if(verbose) printf("hd_floors CASE R3 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[0],Ehat);
      pp[0]=1./EERHORATIOMAX*Ehat;
      ret=-1;
    }
  
  //Ehat/uint ratios  
  if(Ehat<EEUURATIOMIN*pp[1]) 
    {
      if(verbose) printf("hd_floors CASE R4 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[1],Ehat);
      pp[EE0]=Eratio*EEUURATIOMIN*pp[1];
      ret=-1;
    }

  if(Ehat>EEUURATIOMAX*pp[1]) 
    {
      if(verbose) printf("hd_floors CASE R5 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[1],Ehat);
      pp[1]=1./EEUURATIOMAX*Ehat;
      ret=-1;
    }

  //ANDREW why is this SKIP_MAGNFIELD and not just ifndef MAGNFIELD
#ifdef SKIP_MAGNFIELD

  ldouble ucond[4],ucovd[4];
  ldouble bcond[4],bcovd[4],bsq,magpre;
  ldouble etacon[4],etarel[4];
  
  calc_ucon_ucov_from_prims(pp, geom, ucond, ucovd);
  calc_bcon_bcov_bsq_from_4vel(pp, ucond, ucovd, geom, bcond, bcovd, &bsq);
  magpre = 0.5 * bsq;

  //Ehat/uint ratios   
  if(magpre>B2EERATIOMAX*Ehat) 
    {
      if(verbose) printf("rad_floors CASE MR4 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,magpre,Ehat);
      pp[EE0]=Eratio*magpre/B2EERATIOMAX;
      ret=-1;
    }
  
  /* don't check this one
  if(Ehat>EEB2RATIOMAX*magpre) 
    {
    }
  */

#endif //SKIP_MAGNFIELD
#endif //SKIPRADSOURCE

#ifdef EVOLVEPHOTONNUMBER
  
  //velocities of the frames
  ldouble ut[4];ut[1]=pp[VX];ut[2]=pp[VY];ut[3]=pp[VZ];
  ldouble uffcov[4],uffcon[4];
  conv_vels_both(ut,uffcon,uffcov,VELPRIM,VEL4,gg,GG);
  ldouble urfcov[4],urfcon[4];
  ut[1]=pp[FX0];ut[2]=pp[FY0];ut[3]=pp[FZ0];
  conv_vels_both(ut,urfcon,urfcov,VELPRIMRAD,VEL4,gg,GG);

  //radiative stress tensor in the lab frame
  ldouble Rij[4][4];
  calc_Rij(pp,ggg,Rij);

  //Erad in fluid frame
  ldouble Ehatrad = Ehat;

  //radiation temperature
  ldouble Thatrad = calc_ncompt_Thatrad_4vel(pp,ggg,Ehatrad,urfcon,uffcov);
  ldouble Thatrad0 = Thatrad;
  ldouble ThatradBB=calc_LTE_TfromE(Ehatrad);

#ifdef MAXDIFFTRADS

  ldouble maxfac=MAXDIFFTRADS;

  //ANDREW: Trad limiter is stronger near BH
  #ifdef MAXDIFFTRADSNEARBH
  ldouble xxBL[4];
  
  #ifdef PRECOMPUTE_MY2OUT
  get_xxout(geom->ix, geom->iy, geom->iz, xxBL);
  #else
  coco_N(geom->xxvec,xxBL,MYCOORDS,BLCOORDS);
  #endif
  
  ldouble fac=step_function(xxBL[1]-2.5*rhorizonBL,.5*rhorizonBL);
  maxfac = MAXDIFFTRADSNEARBH + fac*(MAXDIFFTRADS-MAXDIFFTRADSNEARBH);
  #endif
  
  //Ceiling: Ramesh: we adjust the photon number when changing the temperature
  if (Thatrad > maxfac * ThatradBB)
  {
    pp[NF] *= (Thatrad / (maxfac * ThatradBB));
    Thatrad = maxfac * ThatradBB;
  }

  // Floor: Ramesh: different floors, depending on whether or not we have synchrotron. Again, adjust photon number
  #ifndef SYNCHROTRON
  if (Thatrad < ThatradBB)
  {
    pp[NF] *= (Thatrad / ThatradBB);
    Thatrad = ThatradBB;
  }
  #else
  if (Thatrad < ThatradBB / maxfac)
  {
    pp[NF] *= (Thatrad * maxfac / ThatradBB);
    Thatrad = ThatradBB / maxfac;
  }
  #endif //SYNCHROTRON
#endif //MAXDIFFTRADS

#ifdef PUTNFFLOOR
  //floor on number of photons - does not work in general
  if(pp[NF]<SMALL)
    {
      //printf("NF floor imposed at  (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[NF],SMALL);
      pp[NF]=SMALL;
    }
#endif
  
#endif //EVOLVEPHOTONNUMBER
#endif //RADIATION
 

  return ret;

}



// Debora - old function  
//**********************************************************************
//**********************************************************************
//**********************************************************************
//checks if rad primitives make sense
/*int
check_floors_rad(ldouble *pp, int whichvel,void *ggg)
{
  //skip floors for some time
  //return 0;

  int verbose=0;
  int ret=0;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;


#ifdef RADIATION  
  ldouble pp2[NV];
  int iv;

  //ff rad energy density
  ldouble Rtt,Ehat,ugas[4],Eratio;
  calc_ff_Rtt(pp,&Rtt,ugas,geom);
  Ehat=-Rtt;
  Eratio=pp[EE0]/Ehat; //ratio of energy densities in two frames

  //absolute rad-frame EE:
  if(pp[EE0]<ERADFLOOR) 
    {
      if(verbose) printf("rhd_floors CASE R0 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[0],pp[EE0]);
      pp[EE0]=ERADFLOOR;
      ret=-1;
     }

  #ifndef SKIPRADSOURCE

  //Ehat/rho ratios 
  if(Ehat<EERHORATIOMIN*pp[0]) 
    {
      if(verbose) printf("hd_floors CASE R2 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[0],pp[EE0]);
      pp[EE0]=Eratio*EERHORATIOMIN*pp[0];
      ret=-1;
    }

  if(Ehat>EERHORATIOMAX*pp[0]) 
    {
      if(verbose) printf("hd_floors CASE R3 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[0],Ehat);
      pp[0]=1./EERHORATIOMAX*Ehat;
      ret=-1;
    }
  
  //Ehat/uint ratios  
  if(Ehat<EEUURATIOMIN*pp[1]) 
    {
      if(verbose) printf("hd_floors CASE R4 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[1],Ehat);
      pp[EE0]=Eratio*EEUURATIOMIN*pp[1];
      ret=-1;
    }

  if(Ehat>EEUURATIOMAX*pp[1]) 
    {
      if(verbose) printf("hd_floors CASE R5 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[1],Ehat);
      pp[1]=1./EEUURATIOMAX*Ehat;
      ret=-1;
    }


#ifdef SKIP_MAGNFIELD

  ldouble ucond[4],ucovd[4];
  ldouble bcond[4],bcovd[4],magpre;  
  ldouble etacon[4],etarel[4];
  for(iv=1;iv<4;iv++)
    ucond[iv]=pp[1+iv];
  conv_vels(ucond,ucond,VELPRIM,VEL4,gg,GG);
 
  indices_21(ucond,ucovd,gg);
  calc_bcon_4vel(pp,ucond,ucovd,bcond);
  indices_21(bcond,bcovd,gg); 
  magpre = dotB(bcond,bcovd)/2.;

  //Ehat/uint ratios 
  
  if(magpre>B2EERATIOMAX*Ehat) 
    {
      if(verbose) printf("rad_floors CASE MR4 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,magpre,Ehat);
      pp[EE0]=Eratio*magpre/B2EERATIOMAX;
      ret=-1;
    }
  


#endif

#endif

#ifdef PUTNFFLOOR
  //floor on number of photons - does not work in general
  if(pp[NF]<SMALL)
    {
      //printf("NF floor imposed at  (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp[NF],SMALL);
      pp[NF]=SMALL;
    }
#endif
#endif
 

  return ret;

}
*/
