/*! \file metric.c
 \brief Metric-related routines
 */


#include "ko.h"


//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns metric determinant sqrt(-g)

ldouble
calc_gdet(ldouble *xx)
{
  return calc_gdet_arb(xx,MYCOORDS);
}
 

///////////////////////////////////////////////////////////////

ldouble
calc_gdet_arb(ldouble *xx,int COORDS)
{
#ifdef METRICNUMERIC
  if(COORDS==MKS1COORDS  || COORDS==MKS2COORDS || COORDS==MKS3COORDS || COORDS==TKS3COORDS || COORDS==MSPH1COORDS || COORDS==TFLATCOORDS)
    return calc_gdet_arb_num(xx,COORDS);
  else
    return calc_gdet_arb_ana(xx,COORDS);
#else
  return calc_gdet_arb_ana(xx,COORDS);
#endif
}


///////////////////////////////////////////////////////////////

ldouble
calc_gdet_arb_ana(ldouble *xx,int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
 
if(coords==SPHCOORDS) {
  return sqrt(Power(x1,4)*Power(Sin(x2),2));
 }

if(coords==CYLCOORDS) {
  return x1;
 }

if(coords==MINKCOORDS) {
  return 1.;
 }

if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;

  return Sqrt(Power(Power(a,2) + 2*Power(x1,2) + 
       Power(a,2)*Cos(2*x2),2)*Power(Sin(x2),2))
    /2.;
 }

if(coords==KSCOORDS) {
  ldouble a=BHSPIN;
return Sqrt(Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)*
	    Power(Sin(x2),2));
 }

if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
  return Sqrt(Power(exp(1.0),2*x1)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)*Power(Sin(x2),2))/2.;
#endif
 } 

if(coords==MKS2COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
  ldouble H0=0.;
#if(MYCOORDS==MKS2COORDS)
  R0=MKSR0;
  H0=MKSH0;
  return (Power(Pi,2)*Sqrt((Power(exp(1.0),2*x1)*Power(H0,2)*Power(Cot(1.5707963267948966*H0),2)*Power(Csc(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2),4))/((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/2.;
#endif
 } 


if(coords==MCYL1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKSR0;
#endif
  return Sqrt(Power(exp(1.0),2*x1)*Power(exp(x1) + R0,2));
 } 

if(coords==MSPH1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MSPH1COORDS)
  R0=MKSR0;
#endif
  return Sqrt(Power(exp(1.0),2*x1)*Power(exp(x1) + R0,4)*Power(Sin(x2),2));
 } 

if(coords==MKER1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKER1COORDS)
  R0=MKSR0;
#endif
  return Sqrt(Power(exp(1.0),2*x1)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)*Power(Sin(x2),2))/2.;
 } 

 return 0.;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns D[gdet,x^idim]/gdet

ldouble
calc_dlgdet(ldouble *xx, int idim)
{
  return calc_dlgdet_arb(xx,idim,MYCOORDS);
}


///////////////////////////////////////////////////////////////

ldouble
calc_dlgdet_arb(ldouble *xx, int idim,int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
 
if(coords==SPHCOORDS) {
;if(idim==0) return  2/x1
;if(idim==1) return  Cot(x2)
;if(idim==2) return  0
;
}

if(coords==CYLCOORDS) {
  if(idim==0) return 1./x1;
  if(idim==1) return 0.;
  if(idim==2) return 0.;
}

if(coords==MINKCOORDS) {
  return 0.;
}

if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;
;if(idim==0) return  (4*x1)/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))
;if(idim==1) return  ((-Power(a,2) + 2*Power(x1,2) + 3*Power(a,2)*Cos(2*x2))*Cot(x2))/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))
;if(idim==2) return  0
;
}

if(coords==KSCOORDS) {
 ldouble a=BHSPIN;
;if(idim==0) return  (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;if(idim==1) return  ((-Power(a,2) + 2*Power(x1,2) + 3*Power(a,2)*Cos(2*x2))*Cot(x2))/(2.*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;if(idim==2) return  0
;
}

if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
;if(idim==0) return  (Power(a,2) + 6*Power(exp(1.0),2*x1) + 8*exp(x1)*R0 + 2*Power(R0,2) + Power(a,2)*Cos(2*x2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;if(idim==1) return  ((-Power(a,2) + 2*Power(exp(x1) + R0,2) + 3*Power(a,2)*Cos(2*x2))*Cot(x2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;if(idim==2) return  0
;
#endif
 }

if(coords==MKS2COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
  ldouble H0=0.;
#if(MYCOORDS==MKS2COORDS)
  R0=MKSR0;
  H0=MKSH0;
;if(idim==0) return  (-(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 8*exp(x1)*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))
;if(idim==1) return  (H0*Power(Sec(H0*Pi*(-0.5 + x2)),2)*(-(Power(Pi,2)*Cot(1.5707963267948966*H0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + 4*Power(a,2)*Power(Pi,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Cot(1.5707963267948966*H0)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (12.566370614359172*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan(0. + 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2))))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))/(1. + 1.*Power(Tan(H0*Pi*(-0.5 + x2)),2))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))
;if(idim==2) return  0
;
#endif

 }


if(coords==MCYL1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKSR0;
#endif
;if(idim==0) return  (2*exp(x1) + R0)/(exp(x1) + R0)
;if(idim==1) return  0
;if(idim==2) return  0
;
 }

if(coords==MSPH1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MSPH1COORDS)
  R0=MKSR0;
#endif
;if(idim==0) return  (3*exp(x1) + R0)/(exp(x1) + R0)
;if(idim==1) return  Cot(x2)
;if(idim==2) return  0
;

 }

if(coords==MKER1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKER1COORDS)
  R0=MKSR0;
#endif
;if(idim==0) return  (Power(a,2) + 6*Power(exp(1.0),2*x1) + 8*exp(x1)*R0 + 2*Power(R0,2) + Power(a,2)*Cos(2*x2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;if(idim==1) return  ((-Power(a,2) + 2*Power(exp(x1) + R0,2) + 3*Power(a,2)*Cos(2*x2))*Cot(x2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;if(idim==2) return  0
;
 }

 return 0;

}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//g_ij

int
calc_g(ldouble *xx, ldouble g[][5])
{
  calc_g_arb(xx,g,MYCOORDS);
  return 0;
}
  

///////////////////////////////////////////////////////////////
//g_ij

int
calc_g_arb(ldouble *xx, ldouble g[][5],int COORDS)
{
  #ifdef METRICNUMERIC
  if(COORDS==MKS1COORDS  || COORDS==MKS2COORDS || COORDS==MKS3COORDS || COORDS==TKS3COORDS || COORDS==MSPH1COORDS || COORDS==TFLATCOORDS)
    calc_g_arb_num(xx,g,COORDS); 
  else
    calc_g_arb_ana(xx,g,COORDS);
  #else
  calc_g_arb_ana(xx,g,COORDS);
  #endif

  return 0;
}


///////////////////////////////////////////////////////////////

int
calc_g_arb_ana(ldouble *xx, ldouble g[][5],int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

  //g[3][4]=calc_gdet_arb(xx,coords);

if(coords==MKER1COORDS) {
#if(MYCOORDS==MKER1COORDS)
  ldouble a=BHSPIN;
  ldouble R0;
  R0=MKSR0;

;g[0][0]= -((Power(a,2) + 2*(Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0) + Power(a,2)*Cos(2*x2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2)))
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= (-4*a*(exp(x1) + R0)*Power(Sin(x2),2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[1][0]= 0
;g[1][1]= (Power(exp(1.0),2*x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))/(Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)
;g[2][3]= 0
;g[3][0]= (-4*a*(exp(x1) + R0)*Power(Sin(x2),2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= ((Power(a,4) + 2*Power(exp(x1) + R0,4) + Power(a,2)*(3*Power(exp(1.0),2*x1) + R0*(2 + 3*R0) + exp(x1)*(2 + 6*R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos(2*x2))*Power(Sin(x2),2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;
#endif
}

if(coords==MCYL1COORDS) {
#if(MYCOORDS==MCYL1COORDS)
  ldouble R0;
  R0=MKSR0;

;g[0][0]= -1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= Power(exp(1.0),2*x1)
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= 1
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= Power(exp(x1) + R0,2)
;
#endif
}

if(coords==MSPH1COORDS) {
#if(MYCOORDS==MSPH1COORDS)
  ldouble R0;
  R0=MKSR0;
;g[0][0]= -1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= Power(exp(1.0),2*x1)
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(exp(x1) + R0,2)
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= Power(exp(x1) + R0,2)*Power(Sin(x2),2)
;
#endif
}

if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0;
#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
;g[0][0]= -((Power(a,2) + 2*(Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0) + Power(a,2)*Cos(2*x2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2)))
;g[0][1]= (4*exp(x1)*(exp(x1) + R0))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[0][2]= 0
;g[0][3]= (-4*a*(exp(x1) + R0)*Power(Sin(x2),2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[1][0]= (4*exp(x1)*(exp(x1) + R0))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[1][1]= (4*Power(exp(1.0),2*x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;g[1][2]= 0
;g[1][3]= (-4*a*exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)
;g[2][3]= 0
;g[3][0]= (-4*a*(exp(x1) + R0)*Power(Sin(x2),2))/(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[3][1]= (-4*a*exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;g[3][2]= 0
;g[3][3]= (4*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;
#endif
}


if(coords==MKS2COORDS) {
  ldouble a=BHSPIN;
  ldouble R0,H0;
  R0=H0=0.;
#if(MYCOORDS==MKS2COORDS)
  R0=MKSR0;
  H0=MKSH0;
;g[0][0]= (-(Power(a,2)/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2))) + ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))/((Power(a,2)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) + (((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))
;g[0][1]= (Power(H0,2)*Power(Pi,4)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))*Power(Cot(1.5707963267948966*H0),2)*Power(Csc(1.5707963267948966 - 1.*ArcTan(0.6366197723675814*Tan(1.5707963267948966*H0)*(-1.5707963267948966 + 1.5707963267948966*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))))),4)*((-8*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0.6366197723675814*Tan(1.5707963267948966*H0)*(-1.5707963267948966 + 1.5707963267948966*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))))),4)*Power(Tan(1.5707963267948966*H0),2))/(Power(H0,2)*Power(Pi,4)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) - (8*R0*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0.6366197723675814*Tan(1.5707963267948966*H0)*(-1.5707963267948966 + 1.5707963267948966*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))))),4)*Power(Tan(1.5707963267948966*H0),2))/(exp(x1)*Power(H0,2)*Power(Pi,4)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3))))/(4.*((Power(a,2)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) + (((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))))
;g[0][2]= 0
;g[0][3]= (2*a*(exp(x1) + R0))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)*((Power(a,2)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) + (((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))))
;g[1][0]= (-2*(exp(x1) + R0)*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(exp(x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)*((Power(a,2)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) + (((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))))
;g[1][1]= -((((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)*((Power(a,2)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) + (((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))))
;g[1][2]= 0
;g[1][3]= (a*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(exp(x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)*((Power(a,2)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) + (((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))))
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= (Power(H0,2)*Power(Pi,4)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))*Power(Cot(1.5707963267948966*H0),2)*Power(Csc(1.5707963267948966 - 1.*ArcTan(0.6366197723675814*Tan(1.5707963267948966*H0)*(-1.5707963267948966 + 1.5707963267948966*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))))),4))/4.
;g[2][3]= 0
;g[3][0]= (2*a*(exp(x1) + R0))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)*((Power(a,2)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) + (((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))))
;g[3][1]= (a*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(exp(x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)*((Power(a,2)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) + (((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))))
;g[3][2]= 0
;g[3][3]= ((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))/((Power(a,2)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),3)) + (((-4*Power(exp(x1) + R0,2))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)) - ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))/(Power(exp(1.0),2*x1)*Power(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2),2)))*Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))
;

#endif


}


if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;
;g[0][0]= -1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= (-2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[1][0]= 0
;g[1][1]= (Power(x1,2) + Power(a,2)*Power(Cos(x2),2))/(Power(a,2) + (-2 + x1)*x1)
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(x1,2) + Power(a,2)*Power(Cos(x2),2)
;g[2][3]= 0
;g[3][0]= (-2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= Power(Sin(x2),2)*(Power(a,2) + Power(x1,2) + (2*Power(a,2)*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;
}


if(coords==KSCOORDS) {
 ldouble a=BHSPIN;
;g[0][0]= -1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[0][1]= (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[0][2]= 0
;g[0][3]= (-2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[1][0]= (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[1][1]= 1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[1][2]= 0
;g[1][3]= -(a*(1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))*Power(Sin(x2),2))
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(x1,2) + Power(a,2)*Power(Cos(x2),2)
;g[2][3]= 0
;g[3][0]= (-2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[3][1]= -(a*(1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))*Power(Sin(x2),2))
;g[3][2]= 0
;g[3][3]= Power(Sin(x2),2)*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2) + Power(a,2)*(1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))*Power(Sin(x2),2))
;
}

if(coords==SPHCOORDS) {
;g[0][0]= -1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= 1
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(x1,2)
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
   ;g[3][3]= Power(x1,2)*Power(Sin(x2),2)
;
}


if(coords==CYLCOORDS) {
;g[0][0]= -1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= 1
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= 1
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= x1*x1 
;
}

if(coords==MINKCOORDS) {
;g[0][0]= -1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= 1
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= 1
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= 1
;
}

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//g^ij

int
calc_G(ldouble *xx, ldouble G[][5])
{
  calc_G_arb(xx,G,MYCOORDS);
  return 0;
}
  

///////////////////////////////////////////////////////////////

int
calc_G_arb(ldouble *xx, ldouble G[][5],int COORDS)
{
  #ifdef METRICNUMERIC
  if(COORDS==MKS1COORDS  || COORDS==MKS2COORDS || COORDS==MKS3COORDS || COORDS==TKS3COORDS || COORDS==MSPH1COORDS || COORDS==TFLATCOORDS)
    calc_G_arb_num(xx,G,COORDS); 
  else
    calc_G_arb_ana(xx,G,COORDS); 
  #else
  calc_G_arb_ana(xx,G,COORDS); 
  #endif

  return 0;
}


///////////////////////////////////////////////////////////////

int
calc_G_arb_ana(ldouble *xx, ldouble G[][5],int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

  if(coords==MCYL1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKSR0;
#endif
;G[0][0]= -1
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= Power(exp(1.0),-2*x1)
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= Power(exp(x1) + R0,-2)
;

  }

  if(coords==MSPH1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MSPH1COORDS)
  R0=MKSR0;
#endif
;G[0][0]= -1
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= Power(exp(1.0),-2*x1)
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= Power(exp(x1) + R0,-2)
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= Power(Csc(x2),2)/Power(exp(x1) + R0,2)
;



  }

  if(coords==MKER1COORDS) {
    ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKER1COORDS)
  R0=MKSR0;
;G[0][0]= -((Power(a,4) + 2*Power(exp(x1) + R0,4) + Power(a,2)*(exp(x1) + R0)*(2 + 3*(exp(x1) + R0)) + Power(a,2)*(Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*Cos(2*x2))/((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))))
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= (-4*a*(exp(x1) + R0))/((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2)))
;G[1][0]= 0
;G[1][1]= (Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))/(Power(exp(1.0),2*x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][3]= 0
;G[3][0]= (-4*a*(exp(x1) + R0))/((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2)))
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= (2*((-2 + exp(x1) + R0)*(exp(x1) + R0) + Power(a,2)*Power(Cos(x2),2))*Power(Csc(x2),2))/((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2)))
;
#endif
  }

  if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
;G[0][0]= -(((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos(x2),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[0][1]= (2*(exp(x1) + R0))/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= (2*(exp(x1) + R0))/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[1][1]= (Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))/(Power(exp(1.0),2*x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[1][2]= 0
;G[1][3]= a/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= a/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[3][2]= 0
;G[3][3]= Power(Csc(x2),2)/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;
#endif
  }

  if(coords==MKS2COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
  ldouble H0=0.;
#if(MYCOORDS==MKS2COORDS)
  R0=MKSR0;
  H0=MKSH0;
;G[0][0]= -(((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))
;G[0][1]= (2*(exp(x1) + R0))/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= (2*(exp(x1) + R0))/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))
;G[1][1]= (Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))/(Power(exp(1.0),2*x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))
;G[1][2]= 0
;G[1][3]= a/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= (4*Power(Cos(H0*Pi*(-0.5 + (0.5*H0 + 0.3183098861837907*ArcTan(0.3183098861837907*Tan(1.5707963267948966*H0)*(-3.141592653589793 + 3.141592653589793*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))))/H0)),4)*Power(Tan(1.5707963267948966*H0),2))/(Power(H0,2)*Power(Pi,4)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= a/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)))
;G[3][2]= 0
;G[3][3]= Power(Csc((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2)/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos((Pi*(1 + Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2))))/2.),2))
;
#endif



  }

if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;
;G[0][0]= -((Power(a,4) + 2*Power(x1,4) + Power(a,2)*x1*(2 + 3*x1) + Power(a,2)*(Power(a,2) + (-2 + x1)*x1)*Cos(2*x2))/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))))
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= (-4*a*x1)/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;G[1][0]= 0
;G[1][1]= (Power(a,2) + (-2 + x1)*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][3]= 0
;G[3][0]= (-4*a*x1)/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= (2*((-2 + x1)*x1 + Power(a,2)*Power(Cos(x2),2))*Power(Csc(x2),2))/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;
}

/*
if(coords==SCHWCOORDS) {
;G[0][0]= x1/(2 - x1)
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= (-2 + x1)/x1
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= Power(x1,-2)
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= Power(Csc(x2),2)/Power(x1,2)
;

}
*/

if(coords==KSCOORDS) {
 ldouble a=BHSPIN;
;G[0][0]= -((x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;G[0][1]= (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[1][1]= (Power(a,2) + (-2 + x1)*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[1][2]= 0
;G[1][3]= a/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= a/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[3][2]= 0
;G[3][3]= Power(Csc(x2),2)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;
}


if(coords==SPHCOORDS) {
;G[0][0]= -1
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= 1
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= Power(x1,-2)
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= Power(Csc(x2),2)/Power(x1,2)
;
}

if(coords==CYLCOORDS) {
;G[0][0]= -1
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= 1
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= Power(x1,-2)
;
}

if(coords==MINKCOORDS) {
;G[0][0]= -1
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= 1
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= 1
;
}

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//Christofells: \Gamma^i_jk

int
calc_Krzysie(ldouble *xx, ldouble Krzys[][4][4])
{
  calc_Krzysie_arb(xx,Krzys,MYCOORDS);
  return 0;
}


///////////////////////////////////////////////////////////////

int
calc_Krzysie_arb(ldouble *xx, ldouble Krzys[][4][4],int COORDS)
{
  #ifdef METRICNUMERIC
  if(COORDS==MKS1COORDS  || COORDS==MKS2COORDS || COORDS==MKS3COORDS || COORDS==TKS3COORDS || COORDS==MSPH1COORDS || COORDS==TFLATCOORDS)
    calc_Krzysie_arb_num(xx,Krzys,COORDS);
  else
    calc_Krzysie_arb_ana(xx,Krzys,COORDS);
  #else
  calc_Krzysie_arb_ana(xx,Krzys,COORDS);
  #endif
  return 0;
}


///////////////////////////////////////////////////////////////

int
calc_Krzysie_at_center(int ix,int iy,int iz, ldouble Krzys[][4][4])
{
  ldouble xx[4];

  xx[0]=global_time;
  xx[1]=get_x(ix,0);
  xx[2]=get_x(iy,1);
  xx[3]=get_x(iz,2);

#if(MODYFIKUJKRZYSIE==0)
  calc_Krzysie(xx,Krzys);
#else

  ldouble Krzys_org[4][4][4];
  //analytical at center
  calc_Krzysie(xx,Krzys_org);
  calc_Krzysie(xx,Krzys);

  //modifying \Gamma ^mu_mu_k

  int kappa,mu;
  ldouble Sk,Wk[4],dS[4],gdet[3],D[4],dxk,Ck;
  for(kappa=1;kappa<=3;kappa++)
    {
      Ck=Krzys_org[0][kappa][0]
	+ Krzys_org[1][kappa][1]
	+ Krzys_org[2][kappa][2]
	+ Krzys_org[3][kappa][3];
 
      Sk=1.e-300 + fabs(Krzys_org[0][kappa][0])
	+ fabs(Krzys_org[1][kappa][1])
	+ fabs(Krzys_org[2][kappa][2])
	+ fabs(Krzys_org[3][kappa][3]);

      for(mu=0;mu<4;mu++)
	Wk[mu]=fabs(Krzys_org[mu][kappa][mu])/Sk;

     
      //center
      xx[0]=global_time;
      xx[1]=get_x(ix,0);
      xx[2]=get_x(iy,1);
      xx[3]=get_x(iz,2);
      gdet[0]=calc_gdet_arb(xx,MYCOORDS);

      //upper face
      xx[0]=global_time;
      xx[1]=get_x(ix,0);
      xx[2]=get_x(iy,1);
      xx[3]=get_x(iz,2);
      if(kappa==1) xx[1]=get_xb(ix+1,0);
      if(kappa==2) xx[2]=get_xb(iy+1,1);
      if(kappa==3) xx[3]=get_xb(iz+1,2);
      gdet[1]=calc_gdet_arb(xx,MYCOORDS);

      //lower face
      xx[0]=global_time;
      xx[1]=get_x(ix,0);
      xx[2]=get_x(iy,1);
      xx[3]=get_x(iz,2);
      if(kappa==1) xx[1]=get_xb(ix,0);
      if(kappa==2) xx[2]=get_xb(iy,1);
      if(kappa==3) xx[3]=get_xb(iz,2);
      gdet[2]=calc_gdet_arb(xx,MYCOORDS);

      //numerical differencing
      if(kappa==1) dxk=get_size_x(ix,0);
      if(kappa==2) dxk=get_size_x(iy,1);
      if(kappa==3) dxk=get_size_x(iz,2);
      D[kappa]=(gdet[1]-gdet[2])/(dxk*gdet[0]);

      //correcting Krzysie
      for(mu=0;mu<4;mu++)
	{
	  Krzys[mu][kappa][mu]+=(D[kappa]-Ck)*Wk[mu];
	  Krzys[mu][mu][kappa]=Krzys[mu][kappa][mu];
	}
    }

#endif
  return 0;
}


///////////////////////////////////////////////////////////////

int
calc_Krzysie_arb_ana(ldouble *xx, ldouble Krzys[][4][4],int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];


  if(coords==MCYL1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKSR0;
#endif
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= 0
;Krzys[0][0][2]= 0
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= 0
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= 0
;Krzys[0][2][0]= 0
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= 0
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= 0
;Krzys[0][3][2]= 0
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= 0
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= 0
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= 1
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= 0
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -1 - R0/exp(x1)
;Krzys[2][0][0]= 0
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= 0
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= 0
;Krzys[2][1][2]= 0
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= 0
;Krzys[2][2][2]= 0
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= 0
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= 0
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= 0
;Krzys[3][0][2]= 0
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= 0
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= exp(x1)/(exp(x1) + R0)
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= 0
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= exp(x1)/(exp(x1) + R0)
;Krzys[3][3][2]= 0
;Krzys[3][3][3]= 0
;
  }


  if(coords==MSPH1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MSPH1COORDS)
  R0=MKSR0;
#endif
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= 0
;Krzys[0][0][2]= 0
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= 0
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= 0
;Krzys[0][2][0]= 0
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= 0
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= 0
;Krzys[0][3][2]= 0
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= 0
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= 0
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= 1
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= -1 - R0/exp(x1)
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(((exp(x1) + R0)*Power(Sin(x2),2))/exp(x1))
;Krzys[2][0][0]= 0
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= 0
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= 0
;Krzys[2][1][2]= exp(x1)/(exp(x1) + R0)
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= exp(x1)/(exp(x1) + R0)
;Krzys[2][2][2]= 0
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= 0
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= -(Cos(x2)*Sin(x2))
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= 0
;Krzys[3][0][2]= 0
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= 0
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= exp(x1)/(exp(x1) + R0)
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= Cot(x2)
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= exp(x1)/(exp(x1) + R0)
;Krzys[3][3][2]= Cot(x2)
;Krzys[3][3][3]= 0
;


  }

  if(coords==MKER1COORDS) {
ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKER1COORDS)
  R0=MKSR0;
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= (2*exp(x1)*(Power(a,2) + Power(exp(x1) + R0,2))*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/((Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][0][2]= (-4*Power(a,2)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= (2*exp(x1)*(Power(a,2) + Power(exp(x1) + R0,2))*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/((Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= (-2*a*exp(x1)*(-Power(a,4) + 3*Power(a,2)*Power(exp(x1) + R0,2) + 6*Power(exp(x1) + R0,4) + Power(a,2)*(-Power(a,2) + Power(exp(x1) + R0,2))*Cos(2*x2))*Power(Sin(x2),2))/((Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][2][0]= (-4*Power(a,2)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= (8*Power(a,3)*(exp(x1) + R0)*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= (-2*a*exp(x1)*(-Power(a,4) + 3*Power(a,2)*Power(exp(x1) + R0,2) + 6*Power(exp(x1) + R0,4) + Power(a,2)*(-Power(a,2) + Power(exp(x1) + R0,2))*Cos(2*x2))*Power(Sin(x2),2))/((Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][3][2]= (8*Power(a,3)*(exp(x1) + R0)*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= (4*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/(exp(x1)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= (4*a*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/(exp(x1)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= ((exp(x1) + R0)*(Power(exp(1.0),3*x1) + 3*Power(exp(1.0),2*x1)*(-1 + R0) + (-2 + R0)*Power(R0,2) + Power(a,2)*(2*exp(x1) + R0) + exp(x1)*R0*(-5 + 3*R0)) + Power(a,2)*(Power(a,2) + exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Cos(x2),2))/((Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][2]= -(((exp(x1) + R0)*(Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0)))/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))))
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= (4*a*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/(exp(x1)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(-((Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*(4*Power(exp(x1) + R0,3) + Power(a,2)*(1 + 3*exp(x1) + 3*R0) + Power(a,2)*(-1 + exp(x1) + R0)*Cos(2*x2))) + 2*(exp(x1) + R0)*(Power(a,4) + 2*Power(exp(x1) + R0,4) + Power(a,2)*(3*Power(exp(1.0),2*x1) + R0*(2 + 3*R0) + exp(x1)*(2 + 6*R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos(2*x2)))*Power(Sin(x2),2))/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[2][0][0]= (-8*Power(a,2)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= (8*a*(exp(x1) + R0)*(Power(a,2) + Power(exp(x1) + R0,2))*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= (Power(a,2)*Power(exp(1.0),2*x1)*Cos(x2)*Sin(x2))/((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][1][2]= (exp(x1)*(exp(x1) + R0))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= (exp(x1)*(exp(x1) + R0))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][2][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= (8*a*(exp(x1) + R0)*(Power(a,2) + Power(exp(x1) + R0,2))*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= (Sin(x2)*(-(Cos(x2)*(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*(Power(a,4) + 2*Power(exp(x1) + R0,4) + Power(a,2)*(3*Power(exp(1.0),2*x1) + R0*(2 + 3*R0) + exp(x1)*(2 + 6*R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos(2*x2))) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Sin(x2)*Sin(2*x2) - Power(a,2)*(Power(a,4) + 2*Power(exp(x1) + R0,4) + Power(a,2)*(3*Power(exp(1.0),2*x1) + R0*(2 + 3*R0) + exp(x1)*(2 + 6*R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos(2*x2))*Sin(x2)*Sin(2*x2)))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= (4*a*exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/((Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[3][0][2]= (-8*a*(exp(x1) + R0)*Cot(x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= (4*a*exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/((Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= (2*exp(x1)*(((-2 + exp(x1) + R0)*(exp(x1) + R0) + Power(a,2)*Power(Cos(x2),2))*((Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*(4*Power(exp(x1) + R0,3) + Power(a,2)*(1 + 3*exp(x1) + 3*R0) + Power(a,2)*(-1 + exp(x1) + R0)*Cos(2*x2)) - 2*(exp(x1) + R0)*(Power(a,4) + 2*Power(exp(x1) + R0,4) + Power(a,2)*(3*Power(exp(1.0),2*x1) + R0*(2 + 3*R0) + exp(x1)*(2 + 6*R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos(2*x2))) + 4*Power(a,2)*(exp(x1) + R0)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2)))/((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[3][2][0]= (-8*a*(exp(x1) + R0)*Cot(x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= ((3*Power(a,4) + 8*Power(a,2)*exp(x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(exp(1.0),4*x1) + 8*Power(a,2)*R0 + 16*Power(a,2)*exp(x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 32*exp(x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*(Power(exp(1.0),2*x1) + (-1 + R0)*R0 + exp(x1)*(-1 + 2*R0)))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= (2*exp(x1)*(((-2 + exp(x1) + R0)*(exp(x1) + R0) + Power(a,2)*Power(Cos(x2),2))*((Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*(4*Power(exp(x1) + R0,3) + Power(a,2)*(1 + 3*exp(x1) + 3*R0) + Power(a,2)*(-1 + exp(x1) + R0)*Cos(2*x2)) - 2*(exp(x1) + R0)*(Power(a,4) + 2*Power(exp(x1) + R0,4) + Power(a,2)*(3*Power(exp(1.0),2*x1) + R0*(2 + 3*R0) + exp(x1)*(2 + 6*R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos(2*x2))) + 4*Power(a,2)*(exp(x1) + R0)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2)))/((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[3][3][2]= ((3*Power(a,4) + 8*Power(a,2)*exp(x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(exp(1.0),4*x1) + 8*Power(a,2)*R0 + 16*Power(a,2)*exp(x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 32*exp(x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*(Power(exp(1.0),2*x1) + (-1 + R0)*R0 + exp(x1)*(-1 + 2*R0)))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][3]= 0
;
#endif

  }

  if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0;
#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
;Krzys[0][0][0]= (8*(exp(x1) + R0)*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][0][1]= (4*exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][0][2]= (-4*Power(a,2)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][0][3]= (8*a*(exp(x1) + R0)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][1][0]= (4*exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][1][1]= (2*Power(exp(1.0),2*x1)*(-3*Power(a,4) - 4*Power(a,2)*exp(x1) + 8*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),4*x1) - 4*Power(a,2)*R0 + 24*Power(exp(1.0),2*x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 24*exp(x1)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 8*Power(R0,3) + 32*exp(x1)*Power(R0,3) + 8*Power(R0,4) - 4*Power(a,2)*(Power(a,2) + exp(x1) + R0)*Cos(2*x2) - Power(a,4)*Cos(4*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][1][2]= (-4*Power(a,2)*exp(x1)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][1][3]= (-4*a*exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][2][0]= (-4*Power(a,2)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][2][1]= (-4*Power(a,2)*exp(x1)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][2][2]= (-2*Power(exp(x1) + R0,2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[0][2][3]= (8*Power(a,3)*(exp(x1) + R0)*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][0]= (8*a*(exp(x1) + R0)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][3][1]= (-4*a*exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][3][2]= (8*Power(a,3)*(exp(x1) + R0)*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][3]= (-2*(exp(x1) + R0)*(Power(a,4) + 3*Power(a,4)*exp(x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*exp(x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*exp(x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*exp(x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(exp(x1) + R0)*(Power(a,2) + exp(x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*exp(x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + exp(x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][0][0]= (4*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/(exp(x1)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][0][1]= (-4*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2))*(2*(exp(x1) + R0) - Power(a,2)*Power(Sin(x2),2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= (4*a*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/(exp(x1)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][1][0]= (-4*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2))*(2*(exp(x1) + R0) - Power(a,2)*Power(Sin(x2),2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][1][1]= (10*Power(a,6) + 8*Power(a,4)*exp(x1) + 32*Power(a,2)*Power(exp(1.0),2*x1) + 36*Power(a,4)*Power(exp(1.0),2*x1) + 16*Power(a,2)*Power(exp(1.0),3*x1) - 64*Power(exp(1.0),4*x1) + 48*Power(a,2)*Power(exp(1.0),4*x1) - 32*Power(exp(1.0),5*x1) + 32*Power(exp(1.0),6*x1) + 32*Power(a,2)*exp(x1)*R0 + 72*Power(a,4)*exp(x1)*R0 + 32*Power(a,2)*Power(exp(1.0),2*x1)*R0 - 192*Power(exp(1.0),3*x1)*R0 + 192*Power(a,2)*Power(exp(1.0),3*x1)*R0 - 128*Power(exp(1.0),4*x1)*R0 + 192*Power(exp(1.0),5*x1)*R0 + 36*Power(a,4)*Power(R0,2) + 16*Power(a,2)*exp(x1)*Power(R0,2) - 192*Power(exp(1.0),2*x1)*Power(R0,2) + 288*Power(a,2)*Power(exp(1.0),2*x1)*Power(R0,2) - 192*Power(exp(1.0),3*x1)*Power(R0,2) + 480*Power(exp(1.0),4*x1)*Power(R0,2) - 64*exp(x1)*Power(R0,3) + 192*Power(a,2)*exp(x1)*Power(R0,3) - 128*Power(exp(1.0),2*x1)*Power(R0,3) + 640*Power(exp(1.0),3*x1)*Power(R0,3) + 48*Power(a,2)*Power(R0,4) - 32*exp(x1)*Power(R0,4) + 480*Power(exp(1.0),2*x1)*Power(R0,4) + 192*exp(x1)*Power(R0,5) + 32*Power(R0,6) + Power(a,2)*(15*Power(a,4) + 16*Power(a,2)*(3*Power(exp(1.0),2*x1) + 3*Power(R0,2) + exp(x1)*(1 + 6*R0)) + 16*(3*Power(exp(1.0),4*x1) + 3*Power(R0,4) + Power(exp(1.0),3*x1)*(-1 + 12*R0) + 2*Power(exp(1.0),2*x1)*(1 - R0 + 9*Power(R0,2)) + exp(x1)*R0*(2 - R0 + 12*Power(R0,2))))*Cos(2*x2) + 2*Power(a,4)*(3*Power(a,2) + 6*Power(exp(1.0),2*x1) + 6*Power(R0,2) + 4*exp(x1)*(1 + 3*R0))*Cos(4*x2) + Power(a,6)*Cos(6*x2))/(4.*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][1][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][3]= (a*(Power(a,4) - 8*Power(a,2)*exp(x1) + 3*Power(a,4)*exp(x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 16*Power(exp(1.0),3*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) - 8*Power(a,2)*R0 + 3*Power(a,4)*R0 - 8*Power(a,2)*exp(x1)*R0 + 48*Power(exp(1.0),2*x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 48*exp(x1)*Power(R0,2) + 24*Power(a,2)*exp(x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 16*Power(R0,3) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*exp(x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(exp(x1) + R0)*(-2 + Power(a,2) + exp(x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*exp(x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + exp(x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][2]= -(((exp(x1) + R0)*(Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0)))/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))))
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= (4*a*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/(exp(x1)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][3][1]= (a*(Power(a,4) - 8*Power(a,2)*exp(x1) + 3*Power(a,4)*exp(x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 16*Power(exp(1.0),3*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) - 8*Power(a,2)*R0 + 3*Power(a,4)*R0 - 8*Power(a,2)*exp(x1)*R0 + 48*Power(exp(1.0),2*x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 48*exp(x1)*Power(R0,2) + 24*Power(a,2)*exp(x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 16*Power(R0,3) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*exp(x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(exp(x1) + R0)*(-2 + Power(a,2) + exp(x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*exp(x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + exp(x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(((Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(a,4) + 3*Power(a,4)*exp(x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*exp(x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*exp(x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*exp(x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(exp(x1) + R0)*(Power(a,2) + exp(x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*exp(x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + exp(x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/(exp(x1)*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)))
;Krzys[2][0][0]= (-8*Power(a,2)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][0][1]= (-8*Power(a,2)*exp(x1)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= (8*a*(exp(x1) + R0)*(Power(a,2) + Power(exp(x1) + R0,2))*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][1][0]= (-8*Power(a,2)*exp(x1)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][1][1]= (-8*Power(a,2)*Power(exp(1.0),2*x1)*(exp(x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][1][2]= (exp(x1)*(exp(x1) + R0))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][1][3]= (a*exp(x1)*(3*Power(a,4) + 16*Power(a,2)*exp(x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 16*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),4*x1) + 16*Power(a,2)*R0 + 16*Power(a,2)*exp(x1)*R0 + 48*Power(exp(1.0),2*x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*exp(x1)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 16*Power(R0,3) + 32*exp(x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*Power(exp(x1) + R0,2))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Sin(2*x2))/(2.*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= (exp(x1)*(exp(x1) + R0))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][2][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= (8*a*(exp(x1) + R0)*(Power(a,2) + Power(exp(x1) + R0,2))*Sin(2*x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][3][1]= (a*exp(x1)*(3*Power(a,4) + 16*Power(a,2)*exp(x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 16*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),4*x1) + 16*Power(a,2)*R0 + 16*Power(a,2)*exp(x1)*R0 + 48*Power(exp(1.0),2*x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*exp(x1)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 16*Power(R0,3) + 32*exp(x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*Power(exp(x1) + R0,2))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Sin(2*x2))/(2.*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= (4*Sin(x2)*(-(Cos(x2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Cos(x2),2))*(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos(x2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2) + Power(a,2)*Cos(x2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Cos(x2),2))*(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2) - 2*Power(a,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Cos(x2),2))*Sin(x2)*Sin(2*x2)))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[3][0][0]= (-4*a*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][0][1]= (4*a*exp(x1)*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][0][2]= (-8*a*(exp(x1) + R0)*Cot(x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][0][3]= (4*Power(a,2)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][1][0]= (4*a*exp(x1)*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][1][1]= (4*a*Power(exp(1.0),2*x1)*(-Power(a,2) + 2*Power(exp(x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][1][2]= (-4*a*exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*Cot(x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][1][3]= (exp(x1)*(Power(a,4) + 3*Power(a,4)*exp(x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*exp(x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*exp(x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*exp(x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(exp(x1) + R0)*(Power(a,2) + exp(x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*exp(x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + exp(x1) + R0)*Cos(4*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][2][0]= (-8*a*(exp(x1) + R0)*Cot(x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][2][1]= (-4*a*exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*Cot(x2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][2][2]= -((a*(exp(x1) + R0))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[3][2][3]= ((3*Power(a,4) + 8*Power(a,2)*exp(x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(exp(1.0),4*x1) + 8*Power(a,2)*R0 + 16*Power(a,2)*exp(x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 32*exp(x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*(Power(exp(1.0),2*x1) + (-1 + R0)*R0 + exp(x1)*(-1 + 2*R0)))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][0]= (4*Power(a,2)*(Power(a,2) - 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][3][1]= (exp(x1)*(Power(a,4) + 3*Power(a,4)*exp(x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*exp(x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*exp(x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*exp(x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(exp(x1) + R0)*(Power(a,2) + exp(x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*exp(x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + exp(x1) + R0)*Cos(4*x2)))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][3][2]= ((3*Power(a,4) + 8*Power(a,2)*exp(x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(exp(1.0),4*x1) + 8*Power(a,2)*R0 + 16*Power(a,2)*exp(x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 32*exp(x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*(Power(exp(1.0),2*x1) + (-1 + R0)*R0 + exp(x1)*(-1 + 2*R0)))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][3]= -((a*(Power(a,4) + 3*Power(a,4)*exp(x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*exp(x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*exp(x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*exp(x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(exp(x1) + R0)*(Power(a,2) + exp(x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*exp(x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + exp(x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(exp(x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;
#endif
}



  if(coords==MKS2COORDS) {
  ldouble a=BHSPIN;
  ldouble R0,H0;
#if(MYCOORDS==MKS2COORDS)
  R0=MKSR0;
  H0=MKSH0;
;Krzys[0][0][0]= ((exp(x1) + R0)*(-((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*(exp(x1) + R0)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(-1 + exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[0][0][1]= -(exp(x1)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(exp(x1) + R0)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(-1 + exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[0][0][2]= (H0*Power(Pi,2)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*(4*Power(exp(x1) + R0,2)*(-(Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) - ((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[0][0][3]= (2*a*(exp(x1) + R0)*(-((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[0][1][0]= -(exp(x1)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(exp(x1) + R0)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(-1 + exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[0][1][1]= (exp(x1)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(-2*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(-(exp(x1)*(exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*exp(x1)*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + (exp(x1) + R0)*(-(exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*exp(x1)*(exp(x1) + R0)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*exp(x1)*(1 + exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[0][1][2]= (Power(a,2)*exp(x1)*H0*Power(Pi,2)*(exp(x1) + R0)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))/((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))
;Krzys[0][1][3]= -((a*exp(x1)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2)))
;Krzys[0][2][0]= (H0*Power(Pi,2)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*(4*Power(exp(x1) + R0,2)*(-(Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) - ((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[0][2][1]= (Power(a,2)*exp(x1)*H0*Power(Pi,2)*(exp(x1) + R0)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))/((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))
;Krzys[0][2][2]= -(Power(H0,2)*Power(Pi,4)*Power(exp(x1) + R0,2)*Power(Cot(1.5707963267948966*H0),2)*Power(Csc(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))
;Krzys[0][2][3]= -((Power(a,3)*H0*Power(Pi,2)*(exp(x1) + R0)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))/((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))
;Krzys[0][3][0]= (2*a*(exp(x1) + R0)*(-((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[0][3][1]= -((a*exp(x1)*((exp(x1) + R0)*(2 + exp(x1) + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2)))
;Krzys[0][3][2]= -((Power(a,3)*H0*Power(Pi,2)*(exp(x1) + R0)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))/((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))
;Krzys[0][3][3]= ((exp(x1) + R0)*((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(2*Power(exp(x1) + R0,3) + Power(a,2)*(1 + exp(x1) + R0) + Power(a,2)*(-1 + exp(x1) + R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][0][0]= ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(-((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*(exp(x1) + R0)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(-1 + exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(2.*exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][0][1]= (Power(a,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + (exp(x1) + R0)*((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(exp(x1) + R0)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(-1 + exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= (a*(Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(-((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][1][0]= (Power(a,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + (exp(x1) + R0)*((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(exp(x1) + R0)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(-1 + exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][1][1]= (8*(exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(-(exp(x1)*(exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*exp(x1)*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*Power(a,2)*(exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*exp(x1)*(exp(x1) + R0)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*exp(x1)*(1 + exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + (Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(-(exp(x1)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*exp(x1)*(exp(x1) + R0)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*exp(x1)*(1 + exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][1][2]= (H0*Power(Pi,2)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*(Power(a,2)*((Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 4*Power(exp(x1) + R0,2)*(-(Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + (Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(-(Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(a,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][1][3]= (a*(-((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(2*Power(exp(x1) + R0,3) + Power(a,2)*(1 + exp(x1) + R0) + Power(a,2)*(-1 + exp(x1) + R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 4*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= (H0*Power(Pi,2)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*(Power(a,2)*((Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 4*Power(exp(x1) + R0,2)*(-(Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + (Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(-(Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(a,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][2][2]= -(Power(H0,2)*Power(Pi,4)*(exp(x1) + R0)*(Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*Power(Cot(1.5707963267948966*H0),2)*Power(Csc(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4))/(4.*exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= (a*(Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*(-((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][3][1]= (a*(-((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(2*Power(exp(x1) + R0,3) + Power(a,2)*(1 + exp(x1) + R0) + Power(a,2)*(-1 + exp(x1) + R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 4*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= ((Power(a,2) + (-2 + exp(x1) + R0)*(exp(x1) + R0))*((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(2*Power(exp(x1) + R0,3) + Power(a,2)*(1 + exp(x1) + R0) + Power(a,2)*(-1 + exp(x1) + R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(2.*exp(x1)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[2][0][0]= (2*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Tan(1.5707963267948966*H0)*(-((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(H0*Power(Pi,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[2][0][1]= (4*exp(x1)*(exp(x1) + R0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Tan(1.5707963267948966*H0)*(Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(H0*Power(Pi,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= (4*a*(exp(x1) + R0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Tan(1.5707963267948966*H0)*(-((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(H0*Power(Pi,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[2][1][0]= (4*exp(x1)*(exp(x1) + R0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Tan(1.5707963267948966*H0)*(Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(H0*Power(Pi,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[2][1][1]= (2*Power(exp(1.0),2*x1)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Tan(1.5707963267948966*H0)*(Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - Power(a,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(H0*Power(Pi,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[2][1][2]= (exp(x1)*(exp(x1) + R0))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))
;Krzys[2][1][3]= (2*a*exp(x1)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Tan(1.5707963267948966*H0)*(-((Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(H0*Power(Pi,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= (exp(x1)*(exp(x1) + R0))/(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))
;Krzys[2][2][2]= (0.020531964509368672*H0*Power(Sec(H0*Pi*(-0.5 + x2)),2)*(Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Cot(1.5707963267948966*H0)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(240.34729839382604 + 240.34729839382604*Power(Tan(H0*Pi*(-0.5 + x2)),2)) + (306.0196847852814*Power(exp(1.0),2*x1) + 612.0393695705628*exp(x1)*R0 + 306.0196847852814*Power(R0,2) + 306.0196847852814*Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Tan(0. + 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2))))))/((Power(exp(1.0),2*x1) + 2.*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(1. + 1.*Power(Tan(H0*Pi*(-0.5 + x2)),2)))
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= (4*a*(exp(x1) + R0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Tan(1.5707963267948966*H0)*(-((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(H0*Power(Pi,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[2][3][1]= (2*a*exp(x1)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Tan(1.5707963267948966*H0)*(-((Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(H0*Power(Pi,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= (2*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Power(Sin(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4)*Tan(1.5707963267948966*H0)*((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(H0*Power(Pi,2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][0][0]= (a*(-((-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*(exp(x1) + R0)*(-Power(a,2) + (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(-1 + exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][0][1]= (a*exp(x1)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][0][2]= (a*H0*Power(Pi,2)*(exp(x1) + R0)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + Power(a,2)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))
;Krzys[3][0][3]= (Power(a,2)*(-((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][1][0]= (a*exp(x1)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][1][1]= (a*Power(exp(1.0),2*x1)*(Power(exp(x1) + R0,3)*(Power(exp(x1) + R0,3) - Power(a,2)*Power(2 + exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),4) + Power(a,2)*Power(exp(x1) + R0,2)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*((exp(x1) + R0)*Power(2 + exp(x1) + R0,2) + (Power(exp(1.0),3*x1) - 2*Power(a,2)*(2 + exp(x1) + R0) + Power(exp(1.0),2*x1)*(5 + 3*R0) + R0*(4 + 5*R0 + Power(R0,2)) + exp(x1)*(4 + 10*R0 + 3*Power(R0,2)))*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)) + Power(a,4)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(exp(x1) + R0 + (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)) + (exp(x1) + R0)*(2*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0)) + (-Power(a,2) + 2*Power(exp(1.0),2*x1) + R0*(3 + 2*R0) + exp(x1)*(3 + 4*R0))*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][1][2]= (a*exp(x1)*H0*Power(Pi,2)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,3)*(2 + exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + Power(a,2)*(2*(exp(x1) + Power(exp(1.0),2*x1) + R0 + 2*exp(x1)*R0 + Power(R0,2)) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))
;Krzys[3][1][3]= (exp(x1)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(-((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(2*Power(exp(x1) + R0,3) + Power(a,2)*(1 + exp(x1) + R0) + Power(a,2)*(-1 + exp(x1) + R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][2][0]= (a*H0*Power(Pi,2)*(exp(x1) + R0)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + Power(a,2)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))
;Krzys[3][2][1]= (a*exp(x1)*H0*Power(Pi,2)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,3)*(2 + exp(x1) + R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + Power(a,2)*(2*(exp(x1) + Power(exp(1.0),2*x1) + R0 + 2*exp(x1)*R0 + Power(R0,2)) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))
;Krzys[3][2][2]= -(a*Power(H0,2)*Power(Pi,4)*(exp(x1) + R0)*Power(Cot(1.5707963267948966*H0),2)*Power(Csc(1.5707963267948966 - 1.*ArcTan(0. + 1.*Tan(H0*Pi*(-0.5 + x2)))),4))/(4.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))
;Krzys[3][2][3]= (H0*Power(Pi,2)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*(Power(a,2)*((Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(-((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][3][0]= (Power(a,2)*(-((exp(x1) + R0)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*Power(exp(x1) + R0,2)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + (Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][3][1]= (exp(x1)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(-((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + 2*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(2*Power(exp(x1) + R0,3) + Power(a,2)*(1 + exp(x1) + R0) + Power(a,2)*(-1 + exp(x1) + R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + 2*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][3][2]= (H0*Power(Pi,2)*Cot(1.5707963267948966*H0)*Power(Sec(H0*Pi*(-0.5 + x2)),2)*(Power(a,2)*((Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(1.0),2*x1) + 2*exp(x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))) + Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*(-((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.) - Power(a,2)*((Power(a,2)*Sin(Pi*Power(Csc(1.5707963267948966*H0),2)*Csc(6.283185307179586*H0*(-0.5 + x2))*Sin(3.141592653589793*H0)*Power(Sin(H0*Pi*(-0.5 + x2)),2)))/2. - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)*Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)))) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) + Power(a,2)*Cos((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;Krzys[3][3][3]= (a*((Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(a,2) + 3*Power(exp(x1) + R0,2))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) + (Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - 2*Power(a,2)*(1 + exp(x1) + R0 - (-1 + exp(x1) + R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*(2*Power(exp(x1) + R0,3) + Power(a,2)*(1 + exp(x1) + R0) + Power(a,2)*(-1 + exp(x1) + R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))) - 2*(exp(x1) + R0)*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)))))/(2.*(Power(exp(x1) + R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2))*Power((exp(x1) + R0)*(Power(exp(x1) + R0,3) + Power(a,2)*(2 + exp(x1) + R0))*Power(Sec((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - Power(a,2)*(2*exp(x1) + Power(exp(1.0),2*x1) + 2*R0 + 2*exp(x1)*R0 + Power(R0,2) + Power(a,2)*Power(Sin((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2) - (Power(a,2) + Power(exp(1.0),2*x1) + 2*exp(x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Tan((Pi*Cot(1.5707963267948966*H0)*Tan(H0*Pi*(-0.5 + x2)))/2.),2)),2))
;
#endif


}
  
if(coords==KSCOORDS) {
 ldouble a=BHSPIN;
;Krzys[0][0][0]= (2*x1*(Power(x1,2) - Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][0][1]= ((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][0][2]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][0][3]= (2*a*x1*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][1][0]= ((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][1][1]= (2*(Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(x1 + Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][1][2]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][1][3]= (a*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][2][0]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][2][1]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][2][2]= (-2*Power(x1,2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[0][2][3]= (2*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][3][0]= (2*a*x1*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][3][1]= (a*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][3][2]= (2*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][3][3]= (-2*x1*Power(Sin(x2),2)*(Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][0][0]= -(((Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][0][1]= -(((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(2*x1 - Power(a,2)*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= (a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][1][0]= -(((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(2*x1 - Power(a,2)*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][1][1]= -(((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Cos(2*x2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][1][2]= -((Power(a,2)*Sin(2*x2))/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;Krzys[1][1][3]= (a*Power(Sin(x2),2)*(Power(a,4)*x1*Power(Cos(x2),4) + Power(x1,2)*(2*x1 + Power(x1,3) - Power(a,2)*Power(Sin(x2),2)) + Power(a,2)*Power(Cos(x2),2)*(2*x1*(-1 + Power(x1,2)) + Power(a,2)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= -((Power(a,2)*Sin(2*x2))/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;Krzys[1][2][2]= -((x1*(Power(a,2) + (-2 + x1)*x1))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= (a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][3][1]= (a*Power(Sin(x2),2)*(Power(a,4)*x1*Power(Cos(x2),4) + Power(x1,2)*(2*x1 + Power(x1,3) - Power(a,2)*Power(Sin(x2),2)) + Power(a,2)*Power(Cos(x2),2)*(2*x1*(-1 + Power(x1,2)) + Power(a,2)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(((Power(a,2) + (-2 + x1)*x1)*Power(Sin(x2),2)*(Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[2][0][0]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][0][1]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= (a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][1][0]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][1][1]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][1][2]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][1][3]= (a*Cos(x2)*Sin(x2)*(Power(x1,3)*(2 + x1) + 2*Power(a,2)*x1*(1 + x1)*Power(Cos(x2),2) + Power(a,4)*Power(Cos(x2),4) + 2*Power(a,2)*x1*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][2][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= (a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][3][1]= (a*Cos(x2)*Sin(x2)*(Power(x1,3)*(2 + x1) + 2*Power(a,2)*x1*(1 + x1)*Power(Cos(x2),2) + Power(a,4)*Power(Cos(x2),4) + 2*Power(a,2)*x1*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= -((Cos(x2)*Sin(x2)*(Power(a,6)*Power(Cos(x2),6) + Power(Cos(x2),4)*(3*Power(a,4)*Power(x1,2) + Power(a,6)*Power(Sin(x2),2)) + Power(Cos(x2),2)*(3*Power(a,2)*Power(x1,4) + 2*Power(a,4)*Power(x1,2)*Power(Sin(x2),2)) + x1*(Power(x1,5) + Power(a,2)*Power(x1,2)*(4 + x1)*Power(Sin(x2),2) + 2*Power(a,4)*Power(Sin(x2),4) + Power(a,4)*Power(Sin(2*x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[3][0][0]= (a*Power(x1,2) - Power(a,3)*Power(Cos(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][0][1]= (a*(Power(x1,2) - Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][0][2]= (-2*a*x1*Cot(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[3][0][3]= (Power(a,2)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][1][0]= (a*(Power(x1,2) - Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][1][1]= (a*Power(x1,2) - Power(a,3)*Power(Cos(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][1][2]= -((a*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))*Cot(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2))
;Krzys[3][1][3]= (Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][2][0]= (-2*a*x1*Cot(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[3][2][1]= -((a*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))*Cot(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2))
;Krzys[3][2][2]= -((a*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[3][2][3]= ((Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)*Cot(x2))/4. + Power(a,2)*x1*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[3][3][0]= (Power(a,2)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][3][1]= (Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][3][2]= ((Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)*Cot(x2))/4. + Power(a,2)*x1*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[3][3][3]= -((a*Power(Sin(x2),2)*(Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;
}

if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= (-2*(Power(a,2) + Power(x1,2))*(Power(a,2) - 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][0][2]= (-4*Power(a,2)*x1*Sin(2*x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= (-2*(Power(a,2) + Power(x1,2))*(Power(a,2) - 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= (2*a*(Power(a,4) - 3*Power(a,2)*Power(x1,2) - 6*Power(x1,4) + Power(a,2)*(Power(a,2) - Power(x1,2))*Cos(2*x2))*Power(Sin(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][2][0]= (-4*Power(a,2)*x1*Sin(2*x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= (8*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= (2*a*(Power(a,4) - 3*Power(a,2)*Power(x1,2) - 6*Power(x1,4) + Power(a,2)*(Power(a,2) - Power(x1,2))*Cos(2*x2))*Power(Sin(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][3][2]= (8*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= -(((Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= (a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= ((Power(a,2) - x1)*x1 - Power(a,2)*(-1 + x1)*Power(Cos(x2),2))/((Power(a,2) + (-2 + x1)*x1)*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][2]= -((x1*(Power(a,2) + (-2 + x1)*x1))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= (a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(((Power(a,2) + (-2 + x1)*x1)*Power(Sin(x2),2)*(Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[2][0][0]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= (a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= (Power(a,2)*Cos(x2)*Sin(x2))/((Power(a,2) + (-2 + x1)*x1)*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][1][2]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][2][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= (a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= -((Cos(x2)*Sin(x2)*(2*Power(a,2)*Power(x1,2)*(Power(a,2) + Power(x1,2))*Power(Cos(x2),2) + Power(a,4)*(Power(a,2) + Power(x1,2))*Power(Cos(x2),4) + x1*(Power(a,2)*Power(x1,3) + Power(x1,5) + 4*Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + 2*Power(a,4)*Power(Sin(x2),4) + Power(a,4)*Power(Sin(2*x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= (4*a*Power(x1,2) - 4*Power(a,3)*Power(Cos(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][0][2]= (-8*a*x1*Cot(x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= (4*a*Power(x1,2) - 4*Power(a,3)*Power(Cos(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= (4*((-2 + x1)*Power(x1,4) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(a,2)*Power(Cos(x2),2)*(2*(-1 + x1)*Power(x1,2) + Power(a,2)*Power(Sin(x2),2))))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][2][0]= (-8*a*x1*Cot(x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= ((3*Power(a,4) + 8*Power(a,2)*x1 + 8*Power(a,2)*Power(x1,2) + 8*Power(x1,4) + 4*Power(a,2)*(Power(a,2) + 2*(-1 + x1)*x1)*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= (4*((-2 + x1)*Power(x1,4) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(a,2)*Power(Cos(x2),2)*(2*(-1 + x1)*Power(x1,2) + Power(a,2)*Power(Sin(x2),2))))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][2]= ((3*Power(a,4) + 8*Power(a,2)*x1 + 8*Power(a,2)*Power(x1,2) + 8*Power(x1,4) + 4*Power(a,2)*(Power(a,2) + 2*(-1 + x1)*x1)*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][3]= 0
;
}


if(coords==SPHCOORDS) {
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= 0
;Krzys[0][0][2]= 0
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= 0
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= 0
;Krzys[0][2][0]= 0
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= 0
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= 0
;Krzys[0][3][2]= 0
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= 0
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= 0
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= 0
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= -x1
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(x1*Power(Sin(x2),2))
;Krzys[2][0][0]= 0
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= 0
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= 0
;Krzys[2][1][2]= 1/x1
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= 1/x1
;Krzys[2][2][2]= 0
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= 0
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= -(Cos(x2)*Sin(x2))
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= 0
;Krzys[3][0][2]= 0
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= 0
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= 1/x1
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= Cot(x2)
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= 1/x1
;Krzys[3][3][2]= Cot(x2)
;Krzys[3][3][3]= 0
;
}

if(coords==CYLCOORDS) {
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= 0
;Krzys[0][0][2]= 0
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= 0
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= 0
;Krzys[0][2][0]= 0
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= 0
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= 0
;Krzys[0][3][2]= 0
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= 0
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= 0
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= 0
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= 0
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -x1
;Krzys[2][0][0]= 0
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= 0
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= 0
;Krzys[2][1][2]= 0
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= 0
;Krzys[2][2][2]= 0
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= 0
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= 0
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= 0
;Krzys[3][0][2]= 0
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= 0
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= 1/x1
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= 0
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= 1/x1
;Krzys[3][3][2]= 0
;Krzys[3][3][3]= 0
;
}

if(coords==MINKCOORDS) {
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= 0
;Krzys[0][0][2]= 0
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= 0
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= 0
;Krzys[0][2][0]= 0
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= 0
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= 0
;Krzys[0][3][2]= 0
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= 0
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= 0
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= 0
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= 0
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= 0
;Krzys[2][0][0]= 0
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= 0
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= 0
;Krzys[2][1][2]= 0
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= 0
;Krzys[2][2][2]= 0
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= 0
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= 0
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= 0
;Krzys[3][0][2]= 0
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= 0
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= 0
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= 0
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= 0
;Krzys[3][3][2]= 0
;Krzys[3][3][3]= 0
;
}

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//fills geometry structure for cell ix,iy,iz

int
fill_geometry(int ix,int iy,int iz,void *geom)
{
  struct geometry *ggg 
    = (struct geometry *) geom;

  ggg->par=-1;
  ggg->ifacedim = -1;
  pick_g(ix,iy,iz,ggg->gg);
  pick_G(ix,iy,iz,ggg->GG);
  /*
  pick_T(tmuup,ix,iy,iz,ggg->tup);
  pick_T(tmulo,ix,iy,iz,ggg->tlo);
  pick_T(emuup,ix,iy,iz,ggg->eup);
  pick_T(emulo,ix,iy,iz,ggg->elo);
  */
  ggg->alpha=sqrt(-1./ggg->GG[0][0]);
  ggg->ix=ix;  ggg->iy=iy;  ggg->iz=iz;
  ggg->xxvec[0]=0.;
  ggg->xxvec[1]=get_x(ix,0);
  ggg->xxvec[2]=get_x(iy,1);
  ggg->xxvec[3]=get_x(iz,2);
  ggg->xx=ggg->xxvec[1];
  ggg->yy=ggg->xxvec[2];
  ggg->zz=ggg->xxvec[3];
  ggg->gdet=ggg->gg[3][4];
  ggg->gttpert=ggg->GG[3][4];
  ggg->coords=MYCOORDS;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//fills geometry structure for cell face ix,iy,iz in idim

int
fill_geometry_face(int ix,int iy,int iz,int idim, void *geom)
{
  if(doingpostproc) //not precalculated
    {
      fill_geometry_face_arb(ix,iy,iz,idim,geom,MYCOORDS);
      return 0;
    }
    

  struct geometry *ggg 
    = (struct geometry *) geom;

  pick_gb(ix,iy,iz,idim,ggg->gg);
  pick_Gb(ix,iy,iz,idim,ggg->GG);

  ggg->par=-1;
  ggg->ifacedim = idim;
  ggg->coords=MYCOORDS;

  if(idim==0)
    {
      /*
      //if these tetrads required at some point then they should be calculated on the go here
	
      pick_Tb(tmuupbx,ix,iy,iz,idim,ggg->tup);
      pick_Tb(tmulobx,ix,iy,iz,idim,ggg->tlo);
      pick_Tb(emuupbx,ix,iy,iz,idim,ggg->eup);
      pick_Tb(emulobx,ix,iy,iz,idim,ggg->elo);
      */
      ggg->xxvec[1]=get_xb(ix,0);
      ggg->xxvec[2]=get_x(iy,1);
      ggg->xxvec[3]=get_x(iz,2);
    }
  if(idim==1)
    {
      /*
      pick_Tb(tmuupby,ix,iy,iz,idim,ggg->tup);
      pick_Tb(tmuloby,ix,iy,iz,idim,ggg->tlo);
      pick_Tb(emuupby,ix,iy,iz,idim,ggg->eup);
      pick_Tb(emuloby,ix,iy,iz,idim,ggg->elo);
      */
      ggg->xxvec[1]=get_x(ix,0);
      ggg->xxvec[2]=get_xb(iy,1);
      ggg->xxvec[3]=get_x(iz,2);
    }
  if(idim==2)
    {
      /*
      pick_Tb(tmuupbz,ix,iy,iz,idim,ggg->tup);
      pick_Tb(tmulobz,ix,iy,iz,idim,ggg->tlo);
      pick_Tb(emuupbz,ix,iy,iz,idim,ggg->eup);
      pick_Tb(emulobz,ix,iy,iz,idim,ggg->elo);
      */
      ggg->xxvec[1]=get_x(ix,0);
      ggg->xxvec[2]=get_x(iy,1);
      ggg->xxvec[3]=get_xb(iz,2);
    }
  ggg->alpha=sqrt(-1./ggg->GG[0][0]);
  ggg->ix=ix;  ggg->iy=iy;  ggg->iz=iz;

  ggg->xxvec[0]=0.;
  ggg->xx=ggg->xxvec[1];
  ggg->yy=ggg->xxvec[2];
  ggg->zz=ggg->xxvec[3];

  ggg->gdet=ggg->gg[3][4];
  ggg->gttpert=ggg->GG[3][4];

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//fills geometry structure for cell ix,iy,iz in arbitrary metric

int
fill_geometry_arb(int ix,int iy,int iz,void *geom,int COORDS)
{
  struct geometry *ggg 
    = (struct geometry *) geom;

  ldouble xxvec[4],xxvecC[4];

  get_xx(ix,iy,iz,xxvec);
  coco_N(xxvec,xxvecC,MYCOORDS,COORDS);

  calc_g_arb(xxvecC,ggg->gg,COORDS);
  calc_G_arb(xxvecC,ggg->GG,COORDS);

  calc_tetrades(ggg->gg,ggg->tup,ggg->tlo,COORDS);
  calc_ZAMOes(ggg->gg,ggg->eup,ggg->elo,COORDS);

  ggg->alpha=sqrt(-1./ggg->GG[0][0]);
  ggg->ix=ix;  ggg->iy=iy;  ggg->iz=iz;
  ggg->ifacedim=-1;

  ggg->xxvec[0]=0.;
  ggg->xxvec[1]=xxvecC[1];
  ggg->xxvec[2]=xxvecC[2];
  ggg->xxvec[3]=xxvecC[3];  

  ggg->xx=xxvecC[1];
  ggg->yy=xxvecC[2];
  ggg->zz=xxvecC[3];
  
  ggg->gdet=calc_gdet_arb(xxvecC,COORDS);
  ggg->gttpert=calc_gttpert_arb(xxvecC,COORDS);
  
  ggg->coords=COORDS;


  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//fills geometry structure for face ix,iy,iz in idim in arbitrary metric

int
fill_geometry_face_arb(int ix,int iy,int iz,int idim, void *geom,int COORDS)
{
  struct geometry *ggg 
    = (struct geometry *) geom;

  ldouble xxvec[4],xxvecC[4];

  get_xx(ix,iy,iz,xxvec);
  
  if(idim==0) //x-face
    xxvec[1]=get_xb(ix,0);
  if(idim==1) //y-face
    xxvec[2]=get_xb(iy,1);
  if(idim==2) //z-face
    xxvec[3]=get_xb(iz,2);

  coco_N(xxvec,xxvecC,MYCOORDS,COORDS);

  calc_g_arb(xxvecC,ggg->gg,COORDS);
  calc_G_arb(xxvecC,ggg->GG,COORDS);

  calc_tetrades(ggg->gg,ggg->tup,ggg->tlo,COORDS);
  calc_ZAMOes(ggg->gg,ggg->eup,ggg->elo,COORDS);

  ggg->par=-1;
  ggg->alpha=sqrt(-1./ggg->GG[0][0]);
  ggg->ix=ix;  ggg->iy=iy;  ggg->iz=iz; ggg->ifacedim=idim;

  ggg->xxvec[0]=0.;
  ggg->xxvec[1]=xxvecC[1];
  ggg->xxvec[2]=xxvecC[2];
  ggg->xxvec[3]=xxvecC[3];  

  ggg->xx=xxvecC[1];
  ggg->yy=xxvecC[2];
  ggg->zz=xxvecC[3];
 
  ggg->gdet=calc_gdet_arb(xxvecC,COORDS);
  ggg->gttpert=calc_gttpert_arb(xxvecC,COORDS);
  ggg->coords=COORDS;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates orthonormal tetrad
//so far limited to grt.neq.0

int
calc_tetrades(ldouble g[][5], ldouble tmuup[][4], ldouble tmulo[][4],int coords)
{
  ldouble blob1,blob1sq;
  
  int method=2;
  int verbose=0;

  //numerical search for eigenvectors
  if(method==1)
    {

      double metric[]={g[0][0],g[0][1],g[0][2],g[0][3],
		       g[1][0],g[1][1],g[1][2],g[1][3],
		       g[2][0],g[2][1],g[2][2],g[2][3],
		       g[3][0],g[3][1],g[3][2],g[3][3]};		       

      gsl_matrix_view m = gsl_matrix_view_array (metric, 4, 4);     
      gsl_vector *eval = gsl_vector_alloc (4);
      gsl_matrix *evec = gsl_matrix_alloc (4, 4);     
      gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (4);       
      gsl_eigen_symmv (&m.matrix, eval, evec, w);     
      gsl_eigen_symmv_free (w);
       
      int i,j;
     
      for (i = 0; i < 4; i++)
	{
	  double eval_i 
	    = gsl_vector_get (eval, i);
	  gsl_vector_view evec_i 
	    = gsl_matrix_column (evec, i);
     
	  if(verbose) 
	    {
	      printf ("eigenvalue = %g\n", eval_i);
	      printf ("eigenvector = \n");
	      gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
	    }

	  for(j=0;j<4;j++) tmulo[j][i]=gsl_matrix_get(evec,i,j)/sqrt(fabs(eval_i));

	}          


      //sorting to get maximal values on the diagonal
      int neworder[4];
      ldouble max;
      for(i=0;i<4;i++)
	{
	  max=-1.;
	  for(j=0;j<4;j++)
	    if(fabs(tmulo[i][j])>max)
	      {
		neworder[i]=j;
		max=fabs(tmulo[i][j]);
	      }
	}
      
      if(verbose)
	{
	  print_tensor(tmulo);
	  printf("no: %d %d %d %d\n",neworder[0],neworder[1],neworder[2],neworder[3]);
	}
      
      ldouble temp[4][4];
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  temp[neworder[i]][j]=tmulo[i][j];
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  tmulo[i][j]=temp[i][j];
      //changing orientation to positive
      for(i=0;i<4;i++) 
	if(tmulo[i][i]<0.) 
	  for(j=0;j<4;j++) tmulo[i][j]*=-1.;

      if(verbose) print_tensor(tmulo);

      //ortonormal -> lab
      for(i=0;i<4;i++)
	{     
	  for(j=0;j<4;j++)
	    tmuup[i][j]=tmulo[i][j];
	}

      //make them covariant
      for(i=0;i<4;i++)
	{	  
	  indices_21(&tmuup[i][0],&tmuup[i][0],g);
	  if(tmuup[i][i]<0.) 
	  for(j=0;j<4;j++) tmuup[i][j]*=-1.;
	}

      if(verbose) print_tensor(tmuup);

      gsl_vector_free (eval);
      gsl_matrix_free (evec);
    }

  //algorithm from HARM's tetrad.c from Mathematica assuming only grt non-zero
  if(method==2)
    {
      ldouble gtt,grr,grt,ghh,gpp;
      int i,j;

      gtt=g[0][0];
      grr=g[1][1];
      grt=g[0][1];
      ghh=g[2][2]; //gthth
      gpp=g[3][3]; //gphph

      //lab -> ortonormal
      blob1=grr - sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt;
      blob1sq=blob1*blob1;
      tmulo[0][0]= (sqrt(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt)/
		    pow((4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt)))*blob1sq,0.25));
      tmulo[0][1] =  (2.0*grt)/(sqrt(4.0*((grt)*(grt)) + (grr - gtt)*(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt))*
				sqrt(fabs(grr - sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt)));
      tmulo[0][2] = 0.0;
      tmulo[0][3] = 0.0;
      tmulo[1][0] = (2.0*grt)/sqrt((4.0*((grt)*(grt)) + (grr - gtt)*(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt))*
				   (grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt));
      tmulo[1][1] = sqrt((grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt)/
			 (sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt)))*(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt)));
      tmulo[1][2] = 0.0;
      tmulo[1][3] = 0.0;
      tmulo[2][0] = 0.0;
      tmulo[2][1] = 0.0;
      tmulo[2][2] = 1.0/fabs(sqrt(ghh));
      tmulo[2][3] = 0.0;
      tmulo[3][0] = 0.0;
      tmulo[3][1] = 0.0;
      tmulo[3][2] = 0.0;
      tmulo[3][3] = 1.0/fabs(sqrt(gpp));

      //ortonormal -> lab
      for(i=0;i<4;i++)
	{     
	  for(j=0;j<4;j++)
	    tmuup[i][j]=tmulo[i][j];
	}

      //make them covariant
      for(i=0;i<4;i++)
	{	  
	  indices_21(tmuup[i],tmuup[i],g);
	  if(tmuup[i][i]<0.) 
	  for(j=0;j<4;j++) tmuup[i][j]*=-1.;
	}
    }

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates base vectors and 1-forms of LNRF to transform lab <--> LNRF

int
calc_ZAMOes(ldouble g[][5], ldouble emuup[][4], ldouble emulo[][4], int coords)
{
  ldouble e2nu,e2psi,e2mu1,e2mu2,omega;
  ldouble gtt,gtph,gphph,grr,gthth;
  int i,j;

  //TODO: to it in a general way
  //the following do not work but are not supposed to be used
  //so the error message suppressed
  //if(coords==KSCOORDS || coords==MKSCOORDS)
  //my_err("ZAMO tetrad calculation not implemented yet.\n");

  gtt=g[0][0];
  gtph=g[0][3];
  gphph=g[3][3];
  grr=g[1][1];
  gthth=g[2][2];

  //Bardeen's 72 coefficients:
  e2nu=-gtt+gtph*gtph/gphph;
  e2psi=gphph;
  e2mu1=grr;
  e2mu2=gthth;
  omega=-gtph/gphph;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	emuup[i][j]=0.;
	emulo[i][j]=0.;
      }

  emuup[0][0]=sqrt(fabs(e2nu));
  emuup[1][1]=sqrt(fabs(e2mu1));
  emuup[2][2]=sqrt(fabs(e2mu2));
  emuup[0][3]=-omega*sqrt(fabs(e2psi));
  emuup[3][3]=sqrt(fabs(e2psi));

  emulo[3][0]=omega*1./sqrt(fabs(e2nu));
  emulo[0][0]=1./sqrt(fabs(e2nu));
  emulo[1][1]=1./sqrt(fabs(e2mu1));
  emulo[2][2]=1./sqrt(fabs(e2mu2));
  emulo[3][3]=1./sqrt(fabs(e2psi));

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for BL -> KS

int
coco_BL2KS(ldouble *xBL, ldouble *xKS)
{
  ldouble r=xBL[1];
  ldouble a=BHSPIN;
  ldouble delta=r*r-2.*r+a*a;
  ldouble sqrta=sqrt(1.-a*a);

  //TODO: these transformations wrong

  //t
  xKS[0]=xBL[0]+2./sqrta*atanh(sqrta/(1.-r))+log(delta);
  //r
  xKS[1]=xBL[1];
  //theta
  xKS[2]=xBL[2];
  //phi
  xKS[3]=xBL[3];
  //?
  //xKS[3]=xBL[3]+a/sqrta*atanh(sqrta/(1.-r));

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for KS -> BL

int
coco_KS2BL(ldouble *xKS, ldouble *xBL)
{
  ldouble r=xKS[1];
  ldouble a=BHSPIN;
  ldouble delta=r*r-2.*r+a*a;
  ldouble sqrta=sqrt(1.-a*a);

  //t
  xBL[0]=xKS[0]-2./sqrta*atanh(sqrta/(1.-r))-log(delta);
  //r
  xBL[1]=xKS[1];
  //theta
  xBL[2]=xKS[2];
  //phi  
  xBL[3]=xKS[3];
  //?
  //xBL[3]=xKS[3]-a/sqrta*atanh(sqrta/(1.-r));

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for KS -> MKS1

int
coco_KS2MKS1(ldouble *xKS, ldouble *xMKS1)
{
  ldouble KSx0=xKS[0];
  ldouble KSx1=xKS[1];
  ldouble KSx2=xKS[2];
  ldouble KSx3=xKS[3];
  ldouble R0=0.;

#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
#endif

  xMKS1[0]
    = KSx0
    ;
  xMKS1[1]
    = log(KSx1-R0)
    ;
  xMKS1[2]
    = KSx2
    ;
  xMKS1[3]
    = KSx3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for KS -> MKS1

int
coco_KS2MKS2(ldouble *xKS, ldouble *xMKS1)
{
  ldouble KSx0=xKS[0];
  ldouble KSx1=xKS[1];
  ldouble KSx2=xKS[2];
  ldouble KSx3=xKS[3];
  ldouble R0=0.,H0=0.;

#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
#endif
#if(MYCOORDS==MKS2COORDS)
  R0=MKSR0;
  H0=MKSH0;
#endif

  xMKS1[0]
    = KSx0
    ;
  xMKS1[1]
    = log(KSx1-R0)
    ;
  xMKS1[2]
    = (0.5*H0+0.3183098861837907*
       ArcTan(0.3183098861837907*(-3.141592653589793
				  +2.*KSx2)*Tan(1.5707963267948966*H0)))/H0
    ;
  xMKS1[3]
    = KSx3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for KS -> MKS3

int
coco_KS2MKS3(ldouble *xKS, ldouble *xMKS)
{
  ldouble KSx0=xKS[0];
  ldouble KSx1=xKS[1];
  ldouble KSx2=xKS[2];
  ldouble KSx3=xKS[3];
  ldouble x0,x1,x2,x3;
  ldouble R0,H0,MY1,MY2,MP0;
  R0=H0=MY1=MY2=0.;

  
#if(MYCOORDS==MKS3COORDS)
  R0=MKSR0;
  H0=MKSH0;
  MY1=MKSMY1;
  MY2=MKSMY2;
  MP0=MKSMP0;
#endif

  x0
= KSx0
;
x1
= Log(KSx1 - R0)
;
x2
= (-(H0*Power(KSx1,MP0)*Pi) - Power(2,1 + MP0)*H0*MY1*Pi + 2*H0*Power(KSx1,MP0)*MY1*Pi + Power(2,1 + MP0)*H0*MY2*Pi + 2*Power(KSx1,MP0)*ArcTan(((-2*KSx2 + Pi)*Tan((H0*Pi)/2.))/Pi))/(2.*H0*(-Power(KSx1,MP0) - Power(2,1 + MP0)*MY1 + 2*Power(KSx1,MP0)*MY1 + Power(2,1 + MP0)*MY2)*Pi)
;
x3
= KSx3
;

  xMKS[0]=x0;
  xMKS[1]=x1;
  xMKS[2]=x2;
  xMKS[3]=x3;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for KS -> TKS3

int
coco_KS2TKS3(ldouble *xKS, ldouble *xMKS)
{
  ldouble KSx0=xKS[0];
  ldouble KSx1=xKS[1];
  ldouble KSx2=xKS[2];
  ldouble KSx3=xKS[3];
  ldouble x0,x1,x2,x3;
  ldouble R0,H0,MY1,MY2,MP0,T0;
  R0=H0=MY1=MY2=MP0=T0=0.;

  
#if(MYCOORDS==TKS3COORDS)
  T0=TKST0;
  R0=MKSR0;
  H0=MKSH0;
  MY1=MKSMY1;
  MY2=MKSMY2;
  MP0=MKSMP0;
#endif

  x1
= Log(KSx1 - R0)
;
x0
= (Power(2,T0)*KSx0)/Power(KSx1,T0)
;
x2
= (-(H0*Power(KSx1,MP0)*Pi) - Power(2,1 + MP0)*H0*MY1*Pi + 2*H0*Power(KSx1,MP0)*MY1*Pi + Power(2,1 + MP0)*H0*MY2*Pi + 2*Power(KSx1,MP0)*ArcTan(((-2*KSx2 + Pi)*Tan((H0*Pi)/2.))/Pi))/(2.*H0*(-Power(KSx1,MP0) - Power(2,1 + MP0)*MY1 + 2*Power(KSx1,MP0)*MY1 + Power(2,1 + MP0)*MY2)*Pi)
;
x3
= KSx3
;

  xMKS[0]=x0;
  xMKS[1]=x1;
  xMKS[2]=x2;
  xMKS[3]=x3;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MINK -> TFLAT

int
coco_MINK2TFLAT(ldouble *xKS, ldouble *xMKS)
{
  ldouble KSx0=xKS[0];
  ldouble KSx1=xKS[1];
  ldouble KSx2=xKS[2];
  ldouble KSx3=xKS[3];


  ldouble T0,x0,x1,x2,x3;
  T0=0.;
#if(MYCOORDS==TFLATCOORDS)
  T0=TFLATT0;
#endif

x1
= KSx1
;
x0
= KSx0/Power(KSx1,T0)
;
x2
= KSx2
;
x3
= KSx3
;

  xMKS[0]=x0;
  xMKS[1]=x1;
  xMKS[2]=x2;
  xMKS[3]=x3;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MKS1 -> KS

int
coco_MKS12KS(ldouble *xMKS1, ldouble *xKS)
{
  ldouble x0=xMKS1[0];
  ldouble x1=xMKS1[1];
  ldouble x2=xMKS1[2];
  ldouble x3=xMKS1[3];
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
#endif

  xKS[0]
    = x0
    ;
  xKS[1]
    = exp(x1) + R0
    ;
  xKS[2]
    = x2
    ;
  xKS[3]
    = x3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MKS2 -> KS

int
coco_MKS22KS(ldouble *xMKS1, ldouble *xKS)
{
  ldouble x0=xMKS1[0];
  ldouble x1=xMKS1[1];
  ldouble x2=xMKS1[2];
  ldouble x3=xMKS1[3];
  ldouble R0,H0;
  R0=H0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
#endif
#if(MYCOORDS==MKS2COORDS)
  R0=MKSR0;
  H0=MKSH0;
#endif

  xKS[0]
    = x0
    ;
  xKS[1]
    = exp(x1) + R0
    ;
  xKS[2]
    = 0.5*M_PI* (1.+Cot(M_PI/2.*H0)* Tan(H0 *M_PI* (-0.5+x2)))
    ;
  xKS[3]
    = x3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MKS3 -> KS

int
coco_MKS32KS(ldouble *xMKS, ldouble *xKS)
{
  ldouble x0=xMKS[0];
  ldouble x1=xMKS[1];
  ldouble x2=xMKS[2];
  ldouble x3=xMKS[3];
  ldouble R0,H0,MY1,MY2,MP0;
  ldouble KSx0,KSx1,KSx2,KSx3;
  R0=H0=MY1=MY2=0.;

#if(MYCOORDS==MKS3COORDS)
  R0=MKSR0;
  H0=MKSH0;
  MY1=MKSMY1;
  MY2=MKSMY2;
  MP0=MKSMP0;
#endif

  KSx0
= x0
;
KSx1
= exp(x1) + R0
;
KSx2
= (Pi*(1 + Cot((H0*Pi)/2.)*Tan(H0*Pi*(-0.5 + (MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0))*(1 - 2*x2) + x2))))/2.
;
KSx3
= x3
;

  xKS[0]=KSx0;
  xKS[1]=KSx1;
  xKS[2]=KSx2;
  xKS[3]=KSx3;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for TKS3 -> KS

int
coco_TKS32KS(ldouble *xMKS, ldouble *xKS)
{
  ldouble x0=xMKS[0];
  ldouble x1=xMKS[1];
  ldouble x2=xMKS[2];
  ldouble x3=xMKS[3];
  ldouble R0,H0,MY1,MY2,MP0,T0;
  ldouble KSx0,KSx1,KSx2,KSx3;
  R0=H0=MY1=MY2=MP0=T0=0.;

#if(MYCOORDS==TKS3COORDS)
  T0=TKST0;
  R0=MKSR0;
  H0=MKSH0;
  MY1=MKSMY1;
  MY2=MKSMY2;
  MP0=MKSMP0;
#endif

  KSx0
= (Power(exp(x1) + R0,T0)*x0)/Power(2,T0)
;
KSx1
= exp(x1) + R0
;
KSx2
= (Pi*(1 + Cot((H0*Pi)/2.)*Tan(H0*Pi*(-0.5 + (MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0))*(1 - 2*x2) + x2))))/2.
;
KSx3
= x3
;

  xKS[0]=KSx0;
  xKS[1]=KSx1;
  xKS[2]=KSx2;
  xKS[3]=KSx3;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for TFLAT -> MINK

int
coco_TFLAT2MINK(ldouble *xMKS, ldouble *xKS)
{
  ldouble x0=xMKS[0];
  ldouble x1=xMKS[1];
  ldouble x2=xMKS[2];
  ldouble x3=xMKS[3];


  ldouble T0=0.,KSx0,KSx1,KSx2,KSx3;
#if(MYCOORDS==TFLATCOORDS)
  T0=TFLATT0;
#endif

KSx0
= x0*Power(x1,T0)
;
KSx1
= x1
;
KSx2
= x2
;
KSx3
= x3
;


  xKS[0]=KSx0;
  xKS[1]=KSx1;
  xKS[2]=KSx2;
  xKS[3]=KSx3;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MCYL1 -> CYL

int
coco_MCYL12CYL(ldouble *xMCYL1, ldouble *xCYL)
{
  ldouble x0=xMCYL1[0];
  ldouble x1=xMCYL1[1];
  ldouble x2=xMCYL1[2];
  ldouble x3=xMCYL1[3];
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKSR0;
#endif

  xCYL[0]
    = x0
    ;
  xCYL[1]
    = exp(x1) + R0
    ;
  xCYL[2]
    = x2
    ;
  xCYL[3]
    = x3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for CYL -> MCYL1

int
coco_CYL2MCYL1(ldouble *xCYL, ldouble *xMCYL1)
{
  ldouble CYLx0=xCYL[0];
  ldouble CYLx1=xCYL[1];
  ldouble CYLx2=xCYL[2];
  ldouble CYLx3=xCYL[3];
  ldouble R0=0.;

#if(MYCOORDS==MCYL1COORDS)
  R0=MKSR0;
#endif

  xMCYL1[0]
    = CYLx0
    ;
  xMCYL1[1]
    = log(CYLx1-R0)
    ;
  xMCYL1[2]
    = CYLx2
    ;
  xMCYL1[3]
    = CYLx3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for CYL -> SPH

int
coco_CYL2SPH(ldouble *xCYL, ldouble *xSPH)
{
  ldouble t=xCYL[0];
  ldouble R=xCYL[1];
  ldouble z=xCYL[2];
  ldouble phi=xCYL[3];
  ldouble r,th;

  r=sqrt(R*R+z*z);
  th=atan2(R,z);

  xSPH[0]
    = t;
  ;
  xSPH[1]
    = r;
  ;
  xSPH[2]
    = th;
  ;
  xSPH[3]
    = phi;
  ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for SPH -> CYL

int
coco_SPH2CYL(ldouble *xCYL, ldouble *xSPH)
{
  ldouble t=xSPH[0];
  ldouble r=xSPH[1];
  ldouble th=xSPH[2];
  ldouble phi=xSPH[3];
  ldouble R,z;

  R=r*sin(th);
  z=r*cos(th);

  xCYL[0]
    = t;
  ;
  xCYL[1]
    = R;
  ;
  xCYL[2]
    = z;
  ;
  xCYL[3]
    = phi;
  ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MSPH1 -> SPH

int
coco_MSPH12SPH(ldouble *xMSPH1, ldouble *xSPH)
{
  ldouble x0=xMSPH1[0];
  ldouble x1=xMSPH1[1];
  ldouble x2=xMSPH1[2];
  ldouble x3=xMSPH1[3];
  ldouble R0=0.;
#if(MYCOORDS==MSPH1COORDS)
  R0=MKSR0;
#endif
  xSPH[0]
    = x0
    ;
  xSPH[1]
    = exp(x1) + R0
    ;
  xSPH[2]
    = x2
    ;
  xSPH[3]
    = x3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for SPH -> MSPH1

int
coco_SPH2MSPH1(ldouble *xSPH, ldouble *xMSPH1)
{
  ldouble SPHx0=xSPH[0];
  ldouble SPHx1=xSPH[1];
  ldouble SPHx2=xSPH[2];
  ldouble SPHx3=xSPH[3];
  ldouble R0=0.;

#if(MYCOORDS==MSPH1COORDS)
  R0=MKSR0;
#endif

  xMSPH1[0]
    = SPHx0
    ;
  xMSPH1[1]
    = log(-R0 + SPHx1)
    ;
  xMSPH1[2]
    = SPHx2
    ;
  xMSPH1[3]
    = SPHx3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MKER1 -> KER

int
coco_MKER12KER(ldouble *xMKER1, ldouble *xKER)
{
  ldouble x0=xMKER1[0];
  ldouble x1=xMKER1[1];
  ldouble x2=xMKER1[2];
  ldouble x3=xMKER1[3];
  ldouble R0=0.;
#if(MYCOORDS==MKER1COORDS)
  R0=MKSR0;
#endif

  xKER[0]
    = x0
    ;
  xKER[1]
    = exp(x1) + R0
    ;
  xKER[2]
    = x2
    ;
  xKER[3]
    = x3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for KER -> MKER1

int
coco_KER2MKER1(ldouble *xKER, ldouble *xMKER1)
{
  ldouble KERx0=xKER[0];
  ldouble KERx1=xKER[1];
  ldouble KERx2=xKER[2];
  ldouble KERx3=xKER[3];
  ldouble R0=0.;

#if(MYCOORDS==MKER1COORDS)
  R0=MKSR0;
#endif

  xMKER1[0]
    = KERx0
    ;
  xMKER1[1]
    = log(KERx1-R0)
    ;
  xMKER1[2]
    = KERx2
    ;
  xMKER1[3]
    = KERx3
    ;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for SPH -> CARTESIAN

int
coco_SPH2MINK(ldouble *xSPH, ldouble *xMINK)
{
  ldouble r=xSPH[1];
  ldouble th=xSPH[2];
  ldouble ph=xSPH[3];
  
  xMINK[0]=xSPH[0];
  xMINK[1]=r*sin(th)*cos(ph);
  xMINK[2]=r*sin(th)*sin(ph);
  xMINK[3]=r*cos(th);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MINK -> SPH

int
coco_MINK2SPH(ldouble *xMINK, ldouble *xSPH)
{
  ldouble x=xMINK[1];
  ldouble y=xMINK[2];
  ldouble z=xMINK[3];
  
  xSPH[0]=xMINK[0];
  xSPH[1]=sqrt(x*x+y*y+z*z);
  xSPH[2]=acos(z/xSPH[1]);
  xSPH[3]=my_atan2(y,x);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//wrapper to convert coordinates

int
coco_N(ldouble *x1, ldouble *x2,int CO1, int CO2)
{
  if(CO1==CO2)
    {
      x2[0]=x1[0];
      x2[1]=x1[1];
      x2[2]=x1[2];
      x2[3]=x1[3];
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==CYLCOORDS)
    {
      ldouble r,th,ph;
      r=x1[1];
      th=x1[2];
      ph=x1[3];
            
      x2[0]=x1[0];
      x2[3]=ph;
      x2[2]=r*cos(th);
      x2[1]=r*sin(th);
    }
  else if((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==MCYL1COORDS)
    {
      //to CYL
      coco_MCYL12CYL(x1,x1); 

      //to BL/SPH
      ldouble R,z,ph;
      R=x1[1];
      z=x1[2];
      ph=x1[3];
            
      x2[0]=x1[0];
      x2[3]=ph;
      x2[1]=sqrt(R*R+z*z);
      x2[2]=asin(R/x2[1]);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==MCYL1COORDS)
    {
      ldouble r,th,ph;

      r=x1[1];
      th=x1[2];
      ph=x1[3];

      //to CYL
      x2[0]=x1[0];
      x2[3]=ph;
      x2[2]=r*cos(th);
      x2[1]=r*sin(th);
      
      //to MCYL1
      coco_CYL2MCYL1(x2,x2);  
    }
  else if((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==CYLCOORDS)
    {
      ldouble R,z,ph;
      R=x1[1];
      z=x1[2];
      ph=x1[3];
            
      x2[0]=x1[0];
      x2[3]=ph;
      x2[1]=sqrt(R*R+z*z);
      x2[2]=atan2(z,R)+M_PI/2.;//asin(R/x2[1]);
    }
  
  else if(((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==SPHCOORDS) ||
	  ((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==SPHCOORDS))
    {
      x2[0]=x1[0];
      x2[1]=x1[1];
      x2[2]=x1[2];
      x2[3]=x1[3];
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==KSCOORDS)
    coco_BL2KS(x1,x2);
  else if (CO1==KSCOORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    coco_KS2BL(x1,x2);
  else if (CO1==KSCOORDS && CO2==MKS1COORDS)
    coco_KS2MKS1(x1,x2);
  else if (CO1==KSCOORDS && CO2==MKS2COORDS)
    coco_KS2MKS2(x1,x2);
  else if (CO1==KSCOORDS && CO2==MKS3COORDS)
    coco_KS2MKS3(x1,x2);
  else if (CO1==KSCOORDS && CO2==TKS3COORDS)
    coco_KS2TKS3(x1,x2);
  else if (CO1==MKS1COORDS && CO2==KSCOORDS)
    coco_MKS12KS(x1,x2);
  else if (CO1==MKS2COORDS && CO2==KSCOORDS)
    coco_MKS22KS(x1,x2);
  else if (CO1==MKS3COORDS && CO2==KSCOORDS)
    coco_MKS32KS(x1,x2);
  else if (CO1==TKS3COORDS && CO2==KSCOORDS)
    coco_TKS32KS(x1,x2);
  else if (CO1==MCYL1COORDS && CO2==CYLCOORDS)
    coco_MCYL12CYL(x1,x2);
  else if (CO1==CYLCOORDS && CO2==MCYL1COORDS)
    coco_CYL2MCYL1(x1,x2);  
  else if (CO1==MSPH1COORDS && (CO2==SPHCOORDS || CO2==SCHWCOORDS || CO2==KERRCOORDS))
    coco_MSPH12SPH(x1,x2);
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MSPH1COORDS)
    coco_SPH2MSPH1(x1,x2);  
  else if (CO1==MKER1COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    coco_MKER12KER(x1,x2);
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MKER1COORDS)
    coco_KER2MKER1(x1,x2);  
  else if (CO1==MKS1COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      coco_MKS12KS(x1,x2);
      coco_KS2BL(x2,x2);
    }
  else if (CO1==MKS2COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      coco_MKS22KS(x1,x2);
      coco_KS2BL(x2,x2);
    }
 else if (CO1==MKS3COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      coco_MKS32KS(x1,x2);
      coco_KS2BL(x2,x2);
    }
 else if (CO1==TKS3COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      coco_TKS32KS(x1,x2);
      coco_KS2BL(x2,x2);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==MKS1COORDS)
    {
      coco_BL2KS(x1,x2);
      coco_KS2MKS1(x2,x2);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==MKS2COORDS)
    {
      coco_BL2KS(x1,x2);
      coco_KS2MKS2(x2,x2);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==MKS3COORDS)
    {
      coco_BL2KS(x1,x2);
      coco_KS2MKS3(x2,x2);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==TKS3COORDS)
    {
      coco_BL2KS(x1,x2);
      coco_KS2TKS3(x2,x2);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MINKCOORDS)
    {
      coco_SPH2MINK(x1,x2);
    }
   else if ((CO1==TFLATCOORDS) && CO2==MINKCOORDS)
    {
      coco_TFLAT2MINK(x1,x2);
    }
  else if (CO1==MSPH1COORDS && CO2==MINKCOORDS)
    {
      coco_MSPH12SPH(x1,x2);
      coco_SPH2MINK(x2,x2);
    }
  else if (CO1==MKER1COORDS && CO2==MINKCOORDS)
    {
      coco_MKER12KER(x1,x2);
      coco_SPH2MINK(x2,x2);
    }
  else if (CO1==MINKCOORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      coco_MINK2SPH(x1,x2);
    }
  else if (CO1==MINKCOORDS && (CO2==TFLATCOORDS))
    {
      coco_MINK2TFLAT(x1,x2);
    }
  else if (CO1==MINKCOORDS && CO2==MKS1COORDS)
    {
      coco_MINK2SPH(x1,x2);
      coco_BL2KS(x2,x2);
      coco_KS2MKS1(x2,x2);
    }
  else if (CO1==MINKCOORDS && CO2==MKS2COORDS)
    {
      coco_MINK2SPH(x1,x2);
      coco_BL2KS(x2,x2);
      coco_KS2MKS2(x2,x2);
    }
 else if (CO1==MINKCOORDS && CO2==MKS3COORDS)
    {
      coco_MINK2SPH(x1,x2);
      coco_BL2KS(x2,x2);
      coco_KS2MKS3(x2,x2);
    }
else if (CO1==MINKCOORDS && CO2==TKS3COORDS)
    {
      coco_MINK2SPH(x1,x2);
      coco_BL2KS(x2,x2);
      coco_KS2TKS3(x2,x2);
    }
 else if (CO1==KSCOORDS && CO2==MINKCOORDS)
    {      
      coco_KS2BL(x1,x2);
      coco_SPH2MINK(x2,x2);
    }
  else if (CO1==MKS1COORDS && CO2==MINKCOORDS)
    {
      coco_MKS12KS(x1,x2);
      coco_KS2BL(x2,x2);
      coco_SPH2MINK(x2,x2);
    }
  else if (CO1==MKS2COORDS && CO2==MINKCOORDS)
    {
      coco_MKS22KS(x1,x2);
      coco_KS2BL(x2,x2);
      coco_SPH2MINK(x2,x2);
    }
  else if (CO1==MKS3COORDS && CO2==MINKCOORDS)
    {
      coco_MKS32KS(x1,x2);
      coco_KS2BL(x2,x2);
      coco_SPH2MINK(x2,x2);
    }
   else if (CO1==TKS3COORDS && CO2==MINKCOORDS)
    {
      coco_TKS32KS(x1,x2);
      coco_KS2BL(x2,x2);
      coco_SPH2MINK(x2,x2);
    }
   else if (CO1==CYLCOORDS && CO2==MINKCOORDS)
    {
      coco_CYL2SPH(x1,x2);
      coco_SPH2MINK(x2,x2);
    }
   else if (CO1==CYLCOORDS && (CO2==SPHCOORDS || CO2==BLCOORDS))
    {
      coco_CYL2SPH(x1,x2);
    }
   else if ((CO1==SPHCOORDS || CO1==BLCOORDS) && CO2==CYLCOORDS)
    {
      coco_SPH2CYL(x1,x2);
    }
   else if (CO1==MINKCOORDS && CO2==CYLCOORDS)
     {
       coco_MINK2SPH(x1,x2);
       coco_SPH2CYL(x2,x2);
     }

  else
    {
      printf("coco: %d -> %d\n",CO1,CO2);
      my_err("coco coordinate conversion not implemented\n");
    }
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for BL -> KS

int
dxdx_BL2KS(ldouble *xx, ldouble dxdx[][4])
{
  ldouble t=xx[0];
  ldouble r=xx[1];
  ldouble th=xx[2];
  ldouble ph=xx[3];

  ldouble a=BHSPIN;
  ldouble delta=r*r-2.*r+a*a;

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[0][1]=2.*r/delta;
  dxdx[3][1]=a/delta;    

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for KS -> BL

int
dxdx_KS2BL(ldouble *xx, ldouble dxdx[][4])
{
  ldouble t=xx[0];
  ldouble r=xx[1];
  ldouble th=xx[2];
  ldouble ph=xx[3];

  ldouble a=BHSPIN;
  ldouble delta=r*r-2.*r+a*a;

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[0][1]=-2.*r/delta;
  dxdx[3][1]=-a/delta;    

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for KS -> MKS1

int
dxdx_KS2MKS1(ldouble *xx, ldouble dxdx[][4])
{
  ldouble KSx0=xx[0];
  ldouble KSx1=xx[1];
  ldouble KSx2=xx[2];
  ldouble KSx3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=1./(KSx1-R0);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for KS -> MKS3

int
dxdx_KS2MKS3(ldouble *xx, ldouble dxdx[][4])
{
  ldouble KSx0=xx[0];
  ldouble KSx1=xx[1];
  ldouble KSx2=xx[2];
  ldouble KSx3=xx[3];

  ldouble R0,H0,MY1,MY2,MP0;
  R0=H0=MY1=MY2=0.;

#if(MYCOORDS==MKS3COORDS)
  R0=MKSR0;
  H0=MKSH0;
  MY1=MKSMY1;
  MY2=MKSMY2;
  MP0=MKSMP0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
;dxdx[0][0]= 1
;dxdx[0][1]= 0
;dxdx[0][2]= 0
;dxdx[0][3]= 0
;dxdx[1][0]= 0
;dxdx[1][1]= 1/(KSx1 - R0)
;dxdx[1][2]= 0
;dxdx[1][3]= 0
;dxdx[2][0]= 0
;dxdx[2][1]= -((Power(2,1 + MP0)*Power(KSx1,-1 + MP0)*MP0*(MY1 - MY2)*ArcTan(((-2*KSx2 + Pi)*Tan((H0*Pi)/2.))/Pi))/(H0*Power(Power(KSx1,MP0)*(1 - 2*MY1) + Power(2,1 + MP0)*(MY1 - MY2),2)*Pi))
;dxdx[2][2]= (-2*Power(KSx1,MP0)*Tan((H0*Pi)/2.))/(H0*(Power(KSx1,MP0)*(-1 + 2*MY1) + Power(2,1 + MP0)*(-MY1 + MY2))*Power(Pi,2)*(1 + (Power(-2*KSx2 + Pi,2)*Power(Tan((H0*Pi)/2.),2))/Power(Pi,2)))
;dxdx[2][3]= 0
;dxdx[3][0]= 0
;dxdx[3][1]= 0
;dxdx[3][2]= 0
;dxdx[3][3]= 1
;



  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for KS -> TKS3

int
dxdx_KS2TKS3(ldouble *xx, ldouble dxdx[][4])
{
  ldouble KSx0=xx[0];
  ldouble KSx1=xx[1];
  ldouble KSx2=xx[2];
  ldouble KSx3=xx[3];

  ldouble R0,H0,MY1,MY2,MP0,T0;
  R0=H0=MY1=MY2=MP0=T0=0.;

#if(MYCOORDS==TKS3COORDS)			
  T0=TKST0;
  R0=MKSR0;
  H0=MKSH0;
  MY1=MKSMY1;
  MY2=MKSMY2;
  MP0=MKSMP0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
;dxdx[0][0]= Power(2,T0)/Power(KSx1,T0)
;dxdx[0][1]= -(Power(2,T0)*KSx0*Power(KSx1,-1 - T0)*T0)
;dxdx[0][2]= 0
;dxdx[0][3]= 0
;dxdx[1][0]= 0
;dxdx[1][1]= 1/(KSx1 - R0)
;dxdx[1][2]= 0
;dxdx[1][3]= 0
;dxdx[2][0]= 0
;dxdx[2][1]= -((Power(2,1 + MP0)*Power(KSx1,-1 + MP0)*MP0*(MY1 - MY2)*ArcTan(((-2*KSx2 + Pi)*Tan((H0*Pi)/2.))/Pi))/(H0*Power(Power(KSx1,MP0)*(1 - 2*MY1) + Power(2,1 + MP0)*(MY1 - MY2),2)*Pi))
;dxdx[2][2]= (-2*Power(KSx1,MP0)*Tan((H0*Pi)/2.))/(H0*(Power(KSx1,MP0)*(-1 + 2*MY1) + Power(2,1 + MP0)*(-MY1 + MY2))*Power(Pi,2)*(1 + (Power(-2*KSx2 + Pi,2)*Power(Tan((H0*Pi)/2.),2))/Power(Pi,2)))
;dxdx[2][3]= 0
;dxdx[3][0]= 0
;dxdx[3][1]= 0
;dxdx[3][2]= 0
;dxdx[3][3]= 1
;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for MINK -> TFLAT

int
dxdx_MINK2TFLAT(ldouble *xx, ldouble dxdx[][4])
{
  ldouble KSx0=xx[0];
  ldouble KSx1=xx[1];
  ldouble KSx2=xx[2];
  ldouble KSx3=xx[3];

  ldouble T0=0.;
#if(MYCOORDS==TFLATCOORDS)
  T0=TFLATT0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
;dxdx[0][0]= Power(KSx1,-T0)
;dxdx[0][1]= -(KSx0*Power(KSx1,-1 - T0)*T0)
;dxdx[0][2]= 0
;dxdx[0][3]= 0
;dxdx[1][0]= 0
;dxdx[1][1]= 1
;dxdx[1][2]= 0
;dxdx[1][3]= 0
;dxdx[2][0]= 0
;dxdx[2][1]= 0
;dxdx[2][2]= 1
;dxdx[2][3]= 0
;dxdx[3][0]= 0
;dxdx[3][1]= 0
;dxdx[3][2]= 0
;dxdx[3][3]= 1
;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for KS -> MKS1

int
dxdx_KS2MKS2(ldouble *xx, ldouble dxdx[][4])
{
  ldouble KSx0=xx[0];
  ldouble KSx1=xx[1];
  ldouble KSx2=xx[2];
  ldouble KSx3=xx[3];
  ldouble R0=0.;
  ldouble H0=0.;
#if(MYCOORDS==MKS2COORDS)
  R0=MKSR0;
  H0=MKSH0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=1./(KSx1-R0);
  dxdx[2][2]=(2*Power(Sin(1.5707963267948966 - 1.*ArcTan(0.6366197723675814*(-1.5707963267948966 + 1.*KSx2)*
							 Tan(1.5707963267948966*H0))),2)*Tan(1.5707963267948966*H0))/(H0*Power(Pi,2));

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for MKS1 -> KS

int
dxdx_MKS12KS(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKSR0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=exp(x1);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for MKS2 -> KS

int
dxdx_MKS22KS(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
  ldouble R0=0.;
  ldouble H0=0.;
#if(MYCOORDS==MKS2COORDS)
  R0=MKSR0;
  H0=MKSH0;
#endif

  int i;
#ifdef APPLY_OMP_SIMD
  //#pragma omp simd
#endif
  for(i=0;i<4;i++)
  {
    int j;
    for(j=0;j<4;j++)
    {
      dxdx[i][j]=delta(i,j);
    }
  }
  
  dxdx[1][1]=exp(x1);
  dxdx[2][2]=
    (H0*Power(Pi,2)*Cot(1.5707963267948966*H0)*Power(Sec(1.*ArcTan(1.*Tan(H0*Pi*(-0.5 + x2)))),2))/2.;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for MKS3 -> KS

int
dxdx_MKS32KS(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

  ldouble R0,H0,MY1,MY2,MP0;
  R0=H0=MY1=MY2=0.;

#if(MYCOORDS==MKS3COORDS)
  R0=MKSR0;
  H0=MKSH0;
  MY1=MKSMY1;
  MY2=MKSMY2;
  MP0=MKSMP0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
;dxdx[0][0]= 1
;dxdx[0][1]= 0
;dxdx[0][2]= 0
;dxdx[0][3]= 0
;dxdx[1][0]= 0
;dxdx[1][1]= exp(x1)
;dxdx[1][2]= 0
;dxdx[1][3]= 0
;dxdx[2][0]= 0
;dxdx[2][1]= -(Power(2,-1 + MP0)*exp(x1)*H0*MP0*(MY1 - MY2)*Power(Pi,2)*Power(exp(x1) + R0,-1 - MP0)*(-1 + 2*x2)*Cot((H0*Pi)/2.)*Power(Sec(H0*Pi*(-0.5 + (MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0))*(1 - 2*x2) + x2)),2))
;dxdx[2][2]= (H0*Power(Pi,2)*(1 - 2*(MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0)))*Cot((H0*Pi)/2.)*Power(Sec(H0*Pi*(-0.5 + (MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0))*(1 - 2*x2) + x2)),2))/2.
;dxdx[2][3]= 0
;dxdx[3][0]= 0
;dxdx[3][1]= 0
;dxdx[3][2]= 0
;dxdx[3][3]= 1
;



  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for TKS3 -> KS

int
dxdx_TKS32KS(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

  ldouble R0,H0,MY1,MY2,MP0,T0;
  R0=H0=MY1=MY2=MP0=T0=0.;

#if(MYCOORDS==TKS3COORDS)
  T0=TKST0;
  R0=MKSR0;
  H0=MKSH0;
  MY1=MKSMY1;
  MY2=MKSMY2;
  MP0=MKSMP0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
;dxdx[0][0]= Power(exp(x1) + R0,T0)/Power(2,T0)
;dxdx[0][1]= (exp(x1)*Power(exp(x1) + R0,-1 + T0)*T0*x0)/Power(2,T0)
;dxdx[0][2]= 0
;dxdx[0][3]= 0
;dxdx[1][0]= 0
;dxdx[1][1]= exp(x1)
;dxdx[1][2]= 0
;dxdx[1][3]= 0
;dxdx[2][0]= 0
;dxdx[2][1]= -(Power(2,-1 + MP0)*exp(x1)*H0*MP0*(MY1 - MY2)*Power(Pi,2)*Power(exp(x1) + R0,-1 - MP0)*(-1 + 2*x2)*Cot((H0*Pi)/2.)*Power(Sec(H0*Pi*(-0.5 + (MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0))*(1 - 2*x2) + x2)),2))
;dxdx[2][2]= (H0*Power(Pi,2)*(1 - 2*(MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0)))*Cot((H0*Pi)/2.)*Power(Sec(H0*Pi*(-0.5 + (MY1 + (Power(2,MP0)*(-MY1 + MY2))/Power(exp(x1) + R0,MP0))*(1 - 2*x2) + x2)),2))/2.
;dxdx[2][3]= 0
;dxdx[3][0]= 0
;dxdx[3][1]= 0
;dxdx[3][2]= 0
;dxdx[3][3]= 1
;

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for TFLAT -> MINK

int
dxdx_TFLAT2MINK(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
  
  ldouble T0=0.;
#if(MYCOORDS==TFLATCOORDS)
  T0=TFLATT0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);

;dxdx[0][0]= Power(x1,T0)
;dxdx[0][1]= T0*x0*Power(x1,-1 + T0)
;dxdx[0][2]= 0
;dxdx[0][3]= 0
;dxdx[1][0]= 0
;dxdx[1][1]= 1
;dxdx[1][2]= 0
;dxdx[1][3]= 0
;dxdx[2][0]= 0
;dxdx[2][1]= 0
;dxdx[2][2]= 1
;dxdx[2][3]= 0
;dxdx[3][0]= 0
;dxdx[3][1]= 0
;dxdx[3][2]= 0
;dxdx[3][3]= 1
;


  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for CYL -> MCYL1

int
dxdx_CYL2MCYL1(ldouble *xx, ldouble dxdx[][4])
{
  ldouble CYLx0=xx[0];
  ldouble CYLx1=xx[1];
  ldouble CYLx2=xx[2];
  ldouble CYLx3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKSR0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=1./(CYLx1-R0);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for MCYL1 -> CYL

int
dxdx_MCYL12CYL(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKSR0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=exp(x1);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for CYL -> SPH

int
dxdx_CYL2SPH(ldouble *xx, ldouble dxdx[][4])
{
  ldouble CYLx0=xx[0];
  ldouble CYLx1=xx[1];
  ldouble CYLx2=xx[2];
  ldouble CYLx3=xx[3];
  ldouble r=sqrt(CYLx1*CYLx1+CYLx2*CYLx2);

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=CYLx1/r;
  dxdx[1][2]=CYLx2/r;
  dxdx[2][1]=CYLx2/r/r;
  dxdx[2][2]=-CYLx1*CYLx2/(CYLx2*r*r);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for SPH -> CYL

int
dxdx_SPH2CYL(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=sin(x2);
  dxdx[1][2]=x1*cos(x2);
  dxdx[2][1]=cos(x2);
  dxdx[2][2]=-x1*sin(x2);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for SPH -> MSPH1

int
dxdx_SPH2MSPH1(ldouble *xx, ldouble dxdx[][4])
{
  ldouble SPHx0=xx[0];
  ldouble SPHx1=xx[1];
  ldouble SPHx2=xx[2];
  ldouble SPHx3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MSPH1COORDS)
  R0=MKSR0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=1./(SPHx1-R0);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for MSPH1 -> SPH

int
dxdx_MSPH12SPH(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MSPH1COORDS)
  R0=MKSR0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=exp(x1);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for KER -> MKER1

int
dxdx_KER2MKER1(ldouble *xx, ldouble dxdx[][4])
{
  ldouble KERx0=xx[0];
  ldouble KERx1=xx[1];
  ldouble KERx2=xx[2];
  ldouble KERx3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MKER1COORDS)
  R0=MKSR0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=1./(KERx1-R0);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for MKER1 -> KER

int
dxdx_MKER12KER(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MKER1COORDS)
  R0=MKSR0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=exp(x1);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************

int
print_Krzysie(ldouble g[][4][4])
{
  int i,j,k;
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  for(k=0;k<4;k++)
	    {
	      printf("%10f ",g[i][j][k]);
	    }
	  printf("\n");
	}
      printf("\n");
    }
  printf("\n");
  return -1;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************

int
print_g(ldouble g[][5])
{
  int i,j,k;
  for(j=0;j<4;j++)
    {
      for(k=0;k<4;k++)
	{
	  printf("%10f ",g[j][k]);
	}
      printf("\n");
    }
  printf("\n");
  return -1;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns location of the horizone in BL

ldouble
r_horizon_BL(ldouble a)
{
  return 1.+sqrt(1-a*a);
}


///////////////////////////////////////////////////////////////
//returns location of ISCO in BL

ldouble
r_ISCO_BL(ldouble ac)
{
  double Z1,Z2;
  Z1=1.+pow(1.-ac*ac,1./3.)*(pow(1.+ac,1./3.)+pow(1.-ac,1./3.));
  Z2=pow(3.*ac*ac+Z1*Z1,1./2.);
  return (3.+Z2-pow((3.-Z1)*(3.+Z1+2.*Z2),1./2.));
}


///////////////////////////////////////////////////////////////
//returns location of the co-rotating marginally bound orbit in BL

ldouble
r_mbound_BL(ldouble a)
{
  return 2.*(1.-a/2.+sqrt(1.-a));
}


///////////////////////////////////////////////////////////////
//returns location of the photon orbit in BL

ldouble
r_photon_BL(ldouble a)
{
  return 2.*(1.-cosl(2./3.*acosl(-a)));
}


//**********************************************************************
/*! \fn int calc_metric()
 \brief Calculates and saves metric and Christoffels at all cell centers and cell faces within the domain as well as ghost cells and corners
 
 The following global arrays are created: g, gbx, gby, gbz, G, Gbx, Gby, Gbz, gKr
 
 Metric arrays g_munu, G^munu correspond to the basic uniform computational grid, x, xb. They are obtained by first computing the metric in physical coordinates as defined by MYCOORDS, and then transforming via the Jacobian dx/dx to x, xb coordinates. The metric array is 4x5 in size, of which 4x4 are the usual components, and g[3][4]=gdet and G[3][4]=gttpert (which seems to be some measure of difference between the numerical gmunu and analytical gmunu)
 
 Kristoffels are 4x4x4 matrices at each cell center and cell wall
 
 For GRMHD, with most of our usual coordinates, the metric and Kristoffels are computed numerically.
 
 \todo Trace through the many functions where metric and Christoffel elements are computed, either analytically or numerically, and check that all the coordinate transformations ('coco' routines) are okay
 */
////////////////////////////////////////////////////////////////////////

int
calc_metric()
{
  int ix,iy,iz,ii;

  #ifdef METRICNUMERIC
  if(GDETIN==0)
    {
      my_err("METRICNUMERIC requires GDETIN==1!\n");
      exit(-1);
    }
  #endif

  //it is not loop_5 because NGCZ != NGCZMET
  // Note: NGCX=NG, NGCY=NG if NY>1 else 0, NGCZMET=NG if NZ>1 and the metric is non-axisymmetric, else NGCZMET=0
#pragma omp parallel for private(ix,iy,iz,ii) 
    for(ix=-NGCX;ix<NX+NGCX;ix++)
      for(iy=-NGCY;iy<NY+NGCY;iy++)
        for(iz=-NGCZMET;iz<NZ+NGCZMET;iz++)
          {


#ifdef METRICAXISYMMETRIC
	if(iz!=0) continue;
#endif

	ldouble gloc[4][5];
	ldouble Kr[4][4][4];
	ldouble eup[4][4],elo[4][4];
	ldouble tup[4][4],tlo[4][4];
	ldouble xx[4];
	int i,j,k;
	//cell centers
	xx[0]=global_time;
	xx[1]=get_x(ix,0);
	xx[2]=get_x(iy,1);
	xx[3]=get_x(iz,2);

	calc_g(xx,gloc);
	for(i=0;i<4;i++)
	  for(j=0;j<4;j++)
	    set_g(g,i,j,ix,iy,iz,gloc[i][j]);
	for(j=0;j<3;j++)
	  set_g(g,j,4,ix,iy,iz,calc_dlgdet(xx,j));
            
	set_g(g,3,4,ix,iy,iz,calc_gdet(xx));

	calc_G(xx,gloc);
	for(i=0;i<4;i++)
	  for(j=0;j<4;j++)
	    set_g(G,i,j,ix,iy,iz,gloc[i][j]);

	set_g(G,3,4,ix,iy,iz,calc_gttpert(xx));
            
	calc_Krzysie_at_center(ix,iy,iz,Kr);
	for(i=0;i<4;i++)
	  for(j=0;j<4;j++)
	    for(k=0;k<4;k++)
	      set_gKr(i,j,k,ix,iy,iz,Kr[i][j][k]);
	      	      
	if(doingpostproc==0) //metric at faces needed only for time evolution
	  {
	
	    //x-faces
	    if(ix==-NG)
	      {
		xx[0]=global_time;
		xx[1]=get_xb(ix,0);
		xx[2]=get_x(iy,1);
		xx[3]=get_x(iz,2);
		calc_g(xx,gloc);
		for(i=0;i<4;i++)
		  for(j=0;j<4;j++)
		    set_gb(gbx,i,j,ix,iy,iz,gloc[i][j],0);

		calc_G(xx,gloc);
		for(i=0;i<4;i++)
		  for(j=0;j<4;j++)
		    set_gb(Gbx,i,j,ix,iy,iz,gloc[i][j],0);
		for(j=0;j<3;j++)
		  set_gb(gbx,j,4,ix,iy,iz,calc_dlgdet(xx,j),0);
		set_gb(gbx,3,4,ix,iy,iz,calc_gdet(xx),0);

		set_gb(Gbx,3,4,ix,iy,iz,calc_gttpert(xx),0);
	      }
	    xx[0]=global_time;
	    xx[1]=get_xb(ix+1,0);
	    xx[2]=get_x(iy,1);
	    xx[3]=get_x(iz,2);
	    calc_g(xx,gloc);
	    for(i=0;i<4;i++)
	      for(j=0;j<4;j++)
		set_gb(gbx,i,j,ix+1,iy,iz,gloc[i][j],0);

	    calc_G(xx,gloc);
	    for(i=0;i<4;i++)
	      for(j=0;j<4;j++)
		set_gb(Gbx,i,j,ix+1,iy,iz,gloc[i][j],0);

	    for(j=0;j<3;j++)
	      set_gb(gbx,j,4,ix+1,iy,iz,calc_dlgdet(xx,j),0);
	    set_gb(gbx,3,4,ix+1,iy,iz,calc_gdet(xx),0);

	    set_gb(Gbx,3,4,ix+1,iy,iz,calc_gttpert(xx),0);

	    //y-faces
	    if(iy==-NG)
	      {
		xx[0]=global_time;
		xx[1]=get_x(ix,0);
		xx[2]=get_xb(iy,1);
		xx[3]=get_x(iz,2);
		calc_g(xx,gloc);
		for(i=0;i<4;i++)
		  for(j=0;j<4;j++)
		    set_gb(gby,i,j,ix,iy,iz,gloc[i][j],1);

		calc_G(xx,gloc);
		for(i=0;i<4;i++)
		  for(j=0;j<4;j++)
		    set_gb(Gby,i,j,ix,iy,iz,gloc[i][j],1);
		for(j=0;j<3;j++)
		  set_gb(gby,j,4,ix,iy,iz,calc_dlgdet(xx,j),1);
		set_gb(gby,3,4,ix,iy,iz,calc_gdet(xx),1);

		set_gb(Gby,3,4,ix,iy,iz,calc_gttpert(xx),1);
	      }


	    xx[0]=global_time;
	    xx[1]=get_x(ix,0);
	    xx[2]=get_xb(iy+1,1);
	    xx[3]=get_x(iz,2);
	    calc_g(xx,gloc);
	    for(i=0;i<4;i++)
	      for(j=0;j<4;j++)
		set_gb(gby,i,j,ix,iy+1,iz,gloc[i][j],1);

	    calc_G(xx,gloc);
	    for(i=0;i<4;i++)
	      for(j=0;j<4;j++)
		set_gb(Gby,i,j,ix,iy+1,iz,gloc[i][j],1);
	    for(j=0;j<3;j++)
	      set_gb(gby,j,4,ix,iy+1,iz,calc_dlgdet(xx,j),1);
	    set_gb(gby,3,4,ix,iy+1,iz,calc_gdet(xx),1);

	    set_gb(Gby,3,4,ix,iy+1,iz,calc_gttpert(xx),1);
		
		  
	    //z-faces
	    if(iz==-NG)
	      {
		xx[0]=global_time;
		xx[1]=get_x(ix,0);
		xx[2]=get_x(iy,1);
		xx[3]=get_xb(iz,2);
		calc_g(xx,gloc);
		for(i=0;i<4;i++)
		  for(j=0;j<4;j++)
		    set_gb(gbz,i,j,ix,iy,iz,gloc[i][j],2);

		calc_G(xx,gloc);
		for(i=0;i<4;i++)
		  for(j=0;j<4;j++)
		    set_gb(Gbz,i,j,ix,iy,iz,gloc[i][j],2);
		for(j=0;j<3;j++)
		  set_gb(gbz,j,4,ix,iy,iz,calc_dlgdet(xx,j),2);
		set_gb(gbz,3,4,ix,iy,iz,calc_gdet(xx),2);

		set_gb(Gbz,3,4,ix,iy,iz,calc_gttpert(xx),2);


	      }
	    xx[0]=global_time;
	    xx[1]=get_x(ix,0);
	    xx[2]=get_x(iy,1);
	    xx[3]=get_xb(iz+1,2);
	    calc_g(xx,gloc);
	    for(i=0;i<4;i++)
	      for(j=0;j<4;j++)
		set_gb(gbz,i,j,ix,iy,iz+1,gloc[i][j],2);	  

	    calc_G(xx,gloc);
	    for(i=0;i<4;i++)
	      for(j=0;j<4;j++)
		set_gb(Gbz,i,j,ix,iy,iz+1,gloc[i][j],2);	  
	    for(j=0;j<3;j++)
	      set_gb(gbz,j,4,ix,iy,iz+1,calc_dlgdet(xx,j),2);
	    set_gb(gbz,3,4,ix,iy,iz+1,calc_gdet(xx),2);
	    
	    set_gb(Gbz,3,4,ix,iy,iz+1,calc_gttpert(xx),2);

	  }
	  }

    //precalculating characteristic radii and parameters
    //works for all metrics but makes sense only for BH problems
    rhorizonBL = r_horizon_BL(BHSPIN);
    rISCOBL = r_ISCO_BL(BHSPIN);
    rmboundBL = r_mbound_BL(BHSPIN);
    rphotonBL = r_photon_BL(BHSPIN);
    etaNT = 1.-sqrt(1.-2./3./r_ISCO_BL(BHSPIN));

  return 0;
}


/******************************************************/
//calculates metric in COORDS using other coordinates
//and appropriate dxdx as base
/******************************************************/

int
calc_g_arb_num(ldouble *xx, ldouble gout[][5],int COORDS)
{
  //base system of coordinates
  ldouble xxb[4],Gb[4][4],dxdx[4][4],G[4][4],g[4][4],gtemp[4][5];
  int BASECOORDS,i,j;

  if(COORDS==MKS1COORDS  || COORDS==MKS2COORDS || COORDS==MKS3COORDS || COORDS==TKS3COORDS)
    BASECOORDS=KSCOORDS;
  else if(COORDS==MSPH1COORDS)
    BASECOORDS=SPHCOORDS;   
  else if(COORDS==MKER1COORDS)
    BASECOORDS=BLCOORDS;
  else if(COORDS==TFLATCOORDS)
    BASECOORDS=MINKCOORDS;
  else
    my_err("calc_g_arb_num() called with unsupported COORDS\n");

  //let's transform coordinates first
  coco_N(xx,xxb,COORDS,BASECOORDS);

  //analytical metric in BASECOORDS
  calc_G_arb_ana(xxb,gtemp,BASECOORDS);
  DLOOP(i,j) Gb[i][j]=gtemp[i][j];

  //transformation matrix
  if(COORDS==MKS1COORDS) dxdx_KS2MKS1(xxb,dxdx);
  if(COORDS==MKS2COORDS) dxdx_KS2MKS2(xxb,dxdx);
  if(COORDS==MKS3COORDS) dxdx_KS2MKS3(xxb,dxdx);
  if(COORDS==MKER1COORDS) dxdx_KER2MKER1(xxb,dxdx);
  if(COORDS==TKS3COORDS) dxdx_KS2TKS3(xxb,dxdx);
  if(COORDS==TFLATCOORDS) dxdx_MINK2TFLAT(xxb,dxdx);
  if(COORDS==MSPH1COORDS) dxdx_SPH2MSPH1(xxb,dxdx);

  //G = dxdx dxdx Gb
  multiply22(Gb,G,dxdx);

  //g = G^-1
  inverse_44matrix(G,g);
  DLOOP(i,j) gout[i][j]=g[i][j];

  //voila!

  return 0;
}
   

/******************************************************/
//calculates metric in COORDS using other coordinates
//and appropriate dxdx as base
/******************************************************/

int
calc_G_arb_num(ldouble *xx, ldouble Gout[][5],int COORDS)
{
  //base system of coordinates
  ldouble xxb[4],Gb[4][4],dxdx[4][4],G[4][4],g[4][4],gtemp[4][5];
  int BASECOORDS,i,j;

  if(COORDS==MKS1COORDS  || COORDS==MKS2COORDS || COORDS==MKS3COORDS || COORDS==TKS3COORDS)
    BASECOORDS=KSCOORDS;
  else if(COORDS==MSPH1COORDS)
    BASECOORDS=SPHCOORDS;
  else if(COORDS==MKER1COORDS)
    BASECOORDS=BLCOORDS;
  else if(COORDS==TFLATCOORDS)
    BASECOORDS=MINKCOORDS;
  else
    my_err("calc_G_arb_num() called with unsupported COORDS\n");

  //let's transform coordinates first
  coco_N(xx,xxb,COORDS,BASECOORDS);

  //analytical metric in BASECOORDS
  calc_G_arb_ana(xxb,gtemp,BASECOORDS);
  DLOOP(i,j) Gb[i][j]=gtemp[i][j];

  //transformation matrix
  if(COORDS==MKS1COORDS) dxdx_KS2MKS1(xxb,dxdx);
  if(COORDS==MKS2COORDS) dxdx_KS2MKS2(xxb,dxdx);
  if(COORDS==MKS3COORDS) dxdx_KS2MKS3(xxb,dxdx);
  if(COORDS==MKER1COORDS) dxdx_KER2MKER1(xxb,dxdx);
  if(COORDS==TKS3COORDS) dxdx_KS2TKS3(xxb,dxdx);
  if(COORDS==TFLATCOORDS) dxdx_MINK2TFLAT(xxb,dxdx);
  if(COORDS==MSPH1COORDS) dxdx_SPH2MSPH1(xxb,dxdx);

  //G = dxdx dxdx Gb
  multiply22(Gb,G,dxdx);
  DLOOP(i,j) Gout[i][j]=G[i][j];

  //voila!

  return 0;
}
   

/******************************************************/
//calculates Kristoffels in COORDS by numerical
//differentiation of metric calculated through calc_g_arb()
/******************************************************/    

struct fg_params
{
  int i,j,k,COORDS;
  ldouble xx[4];
};


///////////////////////////////////////////////////////////////

ldouble fg (double x, void * params)
{
  struct fg_params *par
    = (struct fg_params *) params;
  ldouble xx[4]={par->xx[0],par->xx[1],par->xx[2],par->xx[3]};
  xx[par->k]=x;
  ldouble g[4][5];
  calc_g_arb(xx,g,par->COORDS);

  return g[par->i][par->j];
}


///////////////////////////////////////////////////////////////

int
calc_Krzysie_arb_num(ldouble *xx, ldouble Krzys[][4][4],int COORDS)
{
  //derivative of metric against coordinates
  ldouble dgdx[4][4][4];
  gsl_function F;
  double result, abserr;
  int i,j,k,l;
  struct fg_params par;
  F.function = &fg;
  F.params = &par;
  par.COORDS=COORDS;
  DLOOPA(i)
    par.xx[i]=xx[i];
  DLOOPB(i,j,k)
    {
      //dg_ij , k
      par.i=i; par.j=j; par.k=k;    
      gsl_deriv_central (&F, xx[k], 1.e-6, &result, &abserr);
      dgdx[i][j][k]=result;
    }

  //Kristoffels with all lower indices
  ldouble gamma[4][4][4];
  DLOOPB(i,j,k)
    {
      gamma[i][j][k]=0.5*(dgdx[i][j][k]+dgdx[i][k][j]-dgdx[j][k][i]);
    }

  //\Gamma^i_jk
  ldouble G[4][5];
  calc_G_arb(xx,G,COORDS);
  DLOOPB(i,j,k)
    {
      Krzys[i][j][k]=0.;
      DLOOPA(l)
	Krzys[i][j][k]+=G[i][l]*gamma[l][j][k];
    }

  return 0;
}


/******************************************************/
//calculates gdet by numerical determinant of calc_g_arb()
/******************************************************/    

ldouble
calc_gdet_arb_num(ldouble *xx,int COORDS)
{
  int i,j;
  ldouble g[4][5],a[4][4];
  calc_g_arb(xx,g,COORDS);
  DLOOP(i,j)
    a[i][j]=g[i][j];
  
  return sqrt(-determinant_44matrix(a));

}


/**********************/
//tests numerical calculation of metric
/**********************/

int
test_metric()
{
  //test perturbation to g_tt
  struct geometry geom,geomBL;
  fill_geometry_arb(NX-4,NY/2,0,&geomBL,BLCOORDS);
  fill_geometry(NX-4,NY/2,0,&geom);
  print_metric(geomBL.gg);
  printf("g_tt+1 = %.16e vs %.16e\n",geomBL.gttpert,geomBL.gg[0][0]+1.);
  print_metric(geom.gg);
  printf("g_tt+1 = %.16e vs %.16e\n",geom.gttpert,geomBL.gg[0][0]+1.);
  exit(1);

  //test time evolution of a TKS1 metric
  /*
  struct geometry geom,geomBL;
  for(global_time=0.;global_time<100.;global_time+=1.)
    {
      calc_metric();
      fill_geometry(NX/4,NY/2,0,&geom);
      fill_geometry_arb(NX/4,NY/2,0,&geomBL,MINKCOORDS);
      printf("%f %f %e %e %e\n",global_time,geomBL.xx,geom.gdet,geom.gg[0][1],geom.gg[1][1]);
      int i,j,k;
      ldouble g[4][5],a[4][4],xx[4];
      get_xx(NX/4,NY/2,0,xx);
      calc_g(xx,g);
      DLOOP(i,j)
	a[i][j]=g[i][j];
      ldouble gdet=sqrt(-determinant_44matrix(a));
      print_metric(g);
      print_metric(geom.GG);
      printf("gdet=%e\n",gdet);
      DLOOPB(i,j,k)
	printf("%d %d %d > %e\n",i,j,k,get_gKr(i,j,k,NX/4,NY/2,0));
      getch();
    }
  exit(1);
  */

  //test numerical metric computation
  /*
  struct geometry geom;
  fill_geometry(NX/2,NY/2,0,&geom);
  ldouble xxBL[4];
  coco_N(geom.xxvec,xxBL,MYCOORDS,BLCOORDS);
  print_4vector(geom.xxvec);
  print_4vector(xxBL);
  ldouble gnum[4][5],gana[4][5];
  calc_g_arb_num(geom.xxvec,gnum,MYCOORDS);
  calc_g_arb_ana(geom.xxvec,gana,MYCOORDS);
  print_metric(gnum);
  print_metric(gana);
  ldouble Knum[4][4][4];
  ldouble Kana[4][4][4];
  calc_Krzysie_arb_num(geom.xxvec, Knum, MYCOORDS);
  calc_Krzysie_arb_ana(geom.xxvec, Kana, MYCOORDS);

  int i,j,k;
  printf("gdets: %e vs %e\n",calc_gdet_arb_num(geom.xxvec,MYCOORDS),calc_gdet_arb_ana(geom.xxvec,MYCOORDS));
  DLOOPB(i,j,k)
    {
      //printf("%d %d %d > %e vs %e\n",i,j,k,Knum[i][j][k],Kana[i][j][k]);
    }

  fflush(stdout);
  exit(-1);
  */
  return 0;
}


///////////////////////////////////////////////////////////////
//calculates the perturbed part of g_tt

ldouble
calc_gttpert(ldouble *xx)
{
  return calc_gttpert_arb(xx,MYCOORDS);
}


///////////////////////////////////////////////////////////////
//assumes input (xx,g,G) in COORDS!

ldouble
calc_gttpert_arb(double *xx, int COORDS)
{
  int BASECOORDS,i,j;
  
  if(COORDS==MKS1COORDS  || COORDS==MKS2COORDS || COORDS==MKS3COORDS || COORDS==TKS3COORDS || COORDS==KSCOORDS)
    BASECOORDS=KSCOORDS;
  else if(COORDS==BLCOORDS || COORDS==MKER1COORDS)
    BASECOORDS=BLCOORDS;
  else
    BASECOORDS=-1;

  if(BASECOORDS==-1) //by default - flat spacetime
    return 0.;

  ldouble xxb[4],gttpert;
  //let's transform coordinates first
  coco_N(xx,xxb,COORDS,BASECOORDS);


  ldouble gpert[4][5],Gpert[4][5];
  ldouble gbase[4][5],Gbase[4][5];
  ldouble dxdxp[4][4];

  //analytical metric in BASECOORDS
  calc_g_arb_ana(xxb,gbase,BASECOORDS);
      
  DLOOP(i,j) gpert[i][j]=Gpert[i][j]=0.;

  if(BASECOORDS==KSCOORDS)
    {
      ldouble r,Sigma,costh;
      r=xxb[1];
      costh=cos(xxb[2]);
      Sigma=r*r + BHSPIN*BHSPIN*costh*costh;
      gpert[0][0]=gpert[1][1]=2.*r/Sigma;
      gpert[2][2]=gbase[2][2]-1.0;
      gpert[3][3]=gbase[3][3]-1.0;
     
      //transformation matrix
      if(COORDS==MKS1COORDS) dxdx_MKS12KS(xx,dxdxp);
      if(COORDS==MKS2COORDS) dxdx_MKS22KS(xx,dxdxp);
      if(COORDS==MKS3COORDS) dxdx_MKS32KS(xx,dxdxp);
      if(COORDS==TKS3COORDS) dxdx_TKS32KS(xx,dxdxp);
      if(COORDS==KSCOORDS) {DLOOP(i,j) if(i==j) dxdxp[i][i]=1.0; else dxdxp[i][j]=0.;}
    }


  if(BASECOORDS==BLCOORDS)
    {
      ldouble r,Sigma,costh,mup,DDp;
      r=xxb[1];
      costh=cos(xxb[2]);
      Sigma=r*r + BHSPIN*BHSPIN*costh*costh;
      mup=1.+BHSPIN*BHSPIN*costh*costh;
      DDp=1.-2./r+BHSPIN*BHSPIN/r/r;
      gpert[0][0]=2.*r/Sigma;
      gpert[1][1]=(mup-DDp)/(1.+DDp); //sign different than in HARM but this not used at all
      gpert[2][2]=gbase[2][2]-1.0;
      gpert[3][3]=gbase[3][3]-1.0;
     
      //transformation matrix
      if(COORDS==BLCOORDS) {DLOOP(i,j) if(i==j) dxdxp[i][i]=1.0; else dxdxp[i][j]=0.;}
      if(COORDS==MKER1COORDS) dxdx_MKER12KER(xx,dxdxp);
    }
      
  //from HARM, only _tt component
  int q=0,m,l;
  ldouble ftemp1=(-1.0 + dxdxp[q][q] * dxdxp[q][q]);
  ftemp1 *=-1.0;
  gttpert = (gpert[q][q] * dxdxp[q][q] * dxdxp[q][q]) + ftemp1;
  // now add 15 other terms
  ldouble ftemp2 = 0.;
  for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
      if((l!=q)&&(m!=q)) ftemp2+= gbase[l][m] * dxdxp[l][q] * dxdxp[m][q];
    }
  // add other 15 terms to answer for total of 16 terms
  gttpert+=ftemp2;

  return gttpert;
}


///////////////////////////////////////////////////////////////
//if coordinates spherical like

int
if_coords_sphericallike(int coords)
{
  if(coords==SPHCOORDS ||
     coords==KSCOORDS ||
     coords==BLCOORDS ||
     coords==MKS1COORDS ||
     coords==MKS2COORDS ||
     coords==MKS3COORDS ||
     coords==TKS3COORDS ||
     coords==MSPH1COORDS ||
     coords==MKER1COORDS)
    return 1;

  else
    return 0;
}


///////////////////////////////////////////////////////////////
//if coordinates cylindrical like

int
if_coords_cylindricallike(int coords)
{
  if(coords==CYLCOORDS ||
     coords==MCYL1COORDS)
    return 1;

  else
    return 0;
}
