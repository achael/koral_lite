/*! \file fileop.c
 \brief File operations
 */


#include "ko.h"


/*********************************************/
/*  adds up current quantities to the pavg array */
/*********************************************/

int
save_avg(ldouble dtin)
{
  int ix,iy,iz,iv,ii;

#pragma omp parallel for private(ix,iy,iz,iv) 
  for(ix=0;ix<NX;ix++) 
    {
      for(iy=0;iy<NY;iy++)
	{
	  ldouble avgz[NV+NAVGVARS];
	   for(iv=0;iv<NV+NAVGVARS;iv++)
	     avgz[iv]=0.;
	   for(iz=0;iz<NZ;iz++)
	    {
	      ldouble avg[NV+NAVGVARS];
	      p2avg(ix,iy,iz,avg);

	      //timestep
	      ldouble dt=dtin;

#ifdef RADIATION //if implicit failed, do not take this step into account at all for failed cells
	      if(get_cflag(RADIMPFIXUPFLAG,ix,iy,iz)==0)
#endif
		{

#if (AVGOUTPUT==2) //phi-averaged
		  for(iv=0;iv<NV+NAVGVARS;iv++)
		    avgz[iv]+=avg[iv];
#else //regular, without phi-averaging
		  set_u_scalar(avgselftime,ix,iy,iz,get_u_scalar(avgselftime,ix,iy,iz)+dt);
		  for(iv=0;iv<NV+NAVGVARS;iv++)
		    {
		      set_uavg(pavg,iv,ix,iy,iz, get_uavg(pavg,iv,ix,iy,iz)+avg[iv]*dt);
		    }
#endif
		}
	    }

#if (AVGOUTPUT==2) //phi-averaged
	  for(iv=0;iv<NV+NAVGVARS;iv++)
	    avgz[iv]/=NZ;
	  set_u_scalar(avgselftime,ix,iy,0,get_u_scalar(avgselftime,ix,iy,0)+dt);
	  for(iv=0;iv<NV+NAVGVARS;iv++)
	    {
	      set_uavg(pavg,iv,ix,iy,0, get_uavg(pavg,iv,ix,iy,0)+avgz[iv]*dt);
	    }
#endif
	}
    }

  avgtime+=dtin;

  
  return 0;
}


/*********************************************/
/* opens files etc. */
/*********************************************/

int
fprint_openfiles(char* folder)
{
  char bufor[100];

#ifdef NORESTART
  if(PROCID==0)
    {
      sprintf(bufor,"rm %s/*",folder);
      int i=system(bufor);
    }
  nfout1=0;
  nfout2=0;
#endif

#if(BOXOUTPUT==1) //this one MPI-safe
  if(PROCID==0)
    {
      sprintf(bufor,"%s/boxscalars.dat",folder);
      fout_boxscalars=fopen(bufor,"a");
    }
#endif

#if(VAROUTPUT==1) //this one (to be) MPI-safe
  if(PROCID==0)
    {
      sprintf(bufor,"%s/varscalars.dat",folder);
      fout_varscalars=fopen(bufor,"a");
    }
#endif

#ifndef MPI
  sprintf(bufor,"%s/scalars.dat",folder);
  fout_scalars=fopen(bufor,"a");

  sprintf(bufor,"%s/failures.dat",folder);
  fout_fail=fopen(bufor,"a");
#endif

  return 0;
}


/*********************************************/
/* closes file handles */
/*********************************************/

int
fprint_closefiles()
{
#if(BOXOUTPUT==1)
  if(PROCID==0)
    fclose(fout_boxscalars);
#endif

#if(VAROUTPUT==1)
  if(PROCID==0)
    fclose(fout_varscalars);
#endif

#ifndef MPI
  fclose(fout_scalars);
  fclose(fout_fail);  
#endif

  return 0;
}


/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_gridfile(char* folder)
{
  FILE* out;
  char bufor[50];
  sprintf(bufor,"%s/grid.dat",folder);
  out=fopen(bufor,"w");

  int ix,iy,iz,iv;
  ldouble pp[NV],uu[NV];
  for(iz=0;iz<NZ;iz++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(ix=0;ix<NX;ix++)
	    {
	      struct geometry geom;
	      fill_geometry(ix,iy,iz,&geom);

	      ldouble xxcar[4],xxsph[4];

	      coco_N(geom.xxvec,xxcar,MYCOORDS,MINKCOORDS); 
	      coco_N(geom.xxvec,xxsph,MYCOORDS,KERRCOORDS); 

	      fprintf(out,"%d %d %d %f %f %f %f %f %f %f %f %f\n",
		      ix,iy,iz, geom.xxvec[1],geom.xxvec[2],geom.xxvec[3], xxcar[1],xxcar[2],xxcar[3],xxsph[1],xxsph[2],xxsph[3]);

	    }
	}
    }

  fclose(out);

  return 0;
}


/*********************************************/
/*********************************************/
/*********************************************/
/* prints scalar quantities to scalars.dat */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_scalars(ldouble t, ldouble *scalars, int nscalars)
{
  #ifndef MPI
  calc_scalars(scalars,t);
  int iv;
  //printing scalars
#ifdef TIMEINSCALARSINSEC
      t=timeGU2CGS(t);
#endif
      
  fprintf(fout_scalars,"%e ",t);
  for(iv=0;iv<nscalars;iv++)
    fprintf(fout_scalars,"%.20e ",scalars[iv]);
  fprintf(fout_scalars,"\n");
  fflush(fout_scalars);
  #endif

  return 0;
}


/*********************************************/
/*********************************************/
/*********************************************/
/* prints scalars to boxcorrscalars.dat */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_boxcorrscalars(ldouble t)
{
  ldouble boxcorrscalars[NBOXCORRSCALARS];
  int iv;
  
  calc_boxcorrscalars(boxcorrscalars,t);

  //printing inside calc_boxcorrscalars()
  

  return 0;
}


/*********************************************/
/*********************************************/
/*********************************************/
/* prints scalars to boxscalars.dat */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_boxscalars(ldouble t)
{
  ldouble boxscalars[NBOXSCALARS];
  int iv;
  
  calc_boxscalars(boxscalars,t);
  if(PROCID==0)
    {
      //printing boxscalars
      fprintf(fout_boxscalars,"%e ",t);
      for(iv=0;iv<NBOXSCALARS;iv++)
	fprintf(fout_boxscalars,"%e ",boxscalars[iv]);
      fprintf(fout_boxscalars,"\n");
      fflush(fout_boxscalars);
    }
  

  return 0;
}


/*********************************************/
/*********************************************/
/*********************************************/
/* prints variability to varscalars.dat */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_varscalars(ldouble t)
{
  ldouble varscalars[NVARSCALARS];
  int iv;

  calc_varscalars(varscalars,t);
  if(PROCID==0)
    {
      //printing varscalars
      fprintf(fout_varscalars,"%e ",t);
      for(iv=0;iv<NVARSCALARS;iv++)
	fprintf(fout_varscalars,"%e ",varscalars[iv]);
      fprintf(fout_varscalars,"\n");
      fflush(fout_varscalars);
    }
  

  return 0;
}


/*********************************************/
/*********************************************/
/*********************************************/
/* prints radial profiles to radNNNN.dat */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_radprofiles(ldouble t, int nfile, char* folder, char* prefix)
{
  //#ifdef BHDISK_PROBLEMTYPE 
      char bufor[50],bufor2[50];
      sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

      fout_radprofiles=fopen(bufor,"w");

      ldouble mdotscale = (rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.))/calc_mdotEdd();
      ldouble lumscale = (fluxGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.))/calc_lumEdd();

      fprintf(fout_radprofiles,"# mdotGU2Edd: %e lumGU2Edd: %e\n",mdotscale,lumscale);

      int ix,iv;
      //calculating radial profiles
      ldouble profiles[NRADPROFILES][NX];
      calc_radialprofiles(profiles);
      //printing radial profiles  
      for(ix=0;ix<NX;ix++)
	{
	  ldouble xx[4],xxout[4];
	  get_xx(ix,0,0,xx);
	  coco_N(xx,xxout,MYCOORDS,OUTCOORDS); 
	  if(xxout[1]<rhorizonBL) continue;
	  fprintf(fout_radprofiles,"%e ",xxout[1]);
	  for(iv=0;iv<NRADPROFILES;iv++)
	    fprintf(fout_radprofiles,"%e ",profiles[iv][ix]);
	  fprintf(fout_radprofiles,"\n");
	}
      fclose(fout_radprofiles);
      //#endif
  
  return 0;
}
 

/*********************************************/
/*********************************************/
/*********************************************/
/* prints vertical profiles in a box to boxvertNNNN.dat */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_boxvert(ldouble t, int nfile, char* folder, char* prefix)
{
  //opening the file
  char bufor[50];
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);
  fout_boxvertscalars=fopen(bufor,"w");

  //calculating the numbers to be printed out
  ldouble boxvertscalars[TNY][NBOXVERTSCALARS];
  int iv,iy;
  int sizey=calc_boxvertscalars(boxvertscalars,t); //returns number of cells vertically in the box

  //printing 
  if(PROCID==0)
    {
       for(iy=0;iy<sizey;iy++)
	 {
	   //first scalar is the theta coordinate
	   for(iv=0;iv<NBOXVERTSCALARS;iv++)
	     fprintf(fout_boxvertscalars,"%e ",boxvertscalars[iy][iv]); 
	   fprintf(fout_boxvertscalars,"\n");
	 }
      fflush(fout_boxvertscalars);
    }
  
  //closing the file
  fclose(fout_boxvertscalars);
  
  return 0;
}
 

/*********************************************/
/*********************************************/
/*********************************************/
/* prints theta profiles to thNNNN.dat */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_thprofiles(ldouble t, int nfile, char* folder, char* prefix)
{
  //#ifdef BHDISK_PROBLEMTYPE 
  char bufor[50],bufor2[50];
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

  FILE *fout_thprofiles=fopen(bufor,"w");

  int ix,iy,iv;

  //search for appropriate radial index
  ldouble xx[4],xxBL[4];
  ldouble radius=1.e3;
  #ifdef THPROFRADIUS
  radius=THPROFRADIUS;
  #endif
  for(ix=0;ix<NX;ix++)
    {
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      if(xxBL[1]>radius) break;
    }

  //calculating theta profiles
  ldouble profiles[NTHPROFILES][NY];
  calc_thetaprofiles(profiles);
  //printing th profiles  
  for(iy=0;iy<NY;iy++)
    {
      get_xx(ix,iy,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS); 
      
      fprintf(fout_thprofiles,"%e ",xxBL[2]);
      for(iv=0;iv<NTHPROFILES;iv++)
	fprintf(fout_thprofiles,"%e ",profiles[iv][iy]);
      fprintf(fout_thprofiles,"\n");
    }
  fclose(fout_thprofiles);
  
  return 0;
}
 

/*********************************************/
/*********************************************/
/*********************************************/
/* prints radial profiles to anarelradNNNN.dat   */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_anarelradprofiles(ldouble t, int nfile, char* folder, char* prefix, ldouble profiles[][NX])
{
#ifdef BHDISK_PROBLEMTYPE 
      char bufor[50],bufor2[50];
      sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

      fout_radprofiles=fopen(bufor,"w");

      int ix,iv;
      //printing radial profiles  
      for(ix=0;ix<NX;ix++)
	{
	  ldouble xx[4],xxout[4];
	  get_xx(ix,0,0,xx);
	  coco_N(xx,xxout,MYCOORDS,BLCOORDS); 
	  if(xxout[1]<rhorizonBL) continue;
	  fprintf(fout_radprofiles,"%e ",xxout[1]);
	  for(iv=0;iv<NANARELRADPROFILES;iv++)
	    fprintf(fout_radprofiles,"%e ",profiles[iv][ix]);
	  fprintf(fout_radprofiles,"\n");
	}
      fclose(fout_radprofiles);
#endif
  
  return 0;
}
 

/*********************************************/
/*********************************************/
/*********************************************/
/* prints dumps to files outNNNN.dat and calls gnuplot */
/* codeprim == 1 - prints out code primitives, only coordinates converted to OUTCOORDS */
/* codeprim == 0 - prints ZAMO frame etc primitives - post processing, called by ana.c */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_outfile(ldouble t, int nfile, int codeprim, char* folder, char *prefix)
{
  char bufor[50],bufor2[50];
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);
  fout1=fopen(bufor,"w");
  
 
  //header
  //## nout time problem NX NY NZ
  fprintf(fout1,"## %d %e %d %d %d %d\n",nfout1,t,PROBLEM,NX,NY,NZ);

  sprintf(bufor2,"%s/%s%04d.%s",folder,prefix,nfile,IMAGETYPE);  

  int ix,iy,iz,iv;
  int gclx,gcrx,gcly,gcry,gclz,gcrz;

  //whether print ghost cells or not - default values
  gclx=gcly=gclz=0;
  gcrx=gcry=gcrz=0;
#ifdef PRINTGC_LEFT
  gclx=1;
#endif
#ifdef PRINTGC_RIGHT
  gcrx=1;
#endif
#ifdef PRINTXGC_LEFT
  gclx=1;
#endif
#ifdef PRINTXGC_RIGHT
  gcrx=1;
#endif
#ifdef PRINTYGC_LEFT
  gcly=1;
#endif
#ifdef PRINTYGC_RIGHT
  gcry=1;
#endif
#ifdef PRINTZGC_LEFT
  gclz=1;
#endif
#ifdef PRINTZGC_RIGHT
  gcrz=1;
#endif
#ifdef PRINTZONEMORE
  gcrz=1;
#endif
  
  
  /**************************/  
  /** writing order *********/  
  /**************************/  
 
#ifdef YZXDUMP
  for(iy=0;iy<NY;iy++)
    {
      for(iz=-gclz*NG;iz<NZ+gcrz*NG;iz++)
	{
	  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
	    {
#elif defined(ZXYDUMP)
	      for(iz=0;iz<NZ;iz++)
		{
		  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
		    {
		      for(iy=-gcly*NG;iy<NY+gcry*NG;iy++)
			{
#elif defined(YSLICE)
			  for(iy=YSLICE;iy<YSLICE+1;iy++)
			    {
			      for(iz=0;iz<NZ+gcrz*NG;iz++)
				{
				  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
				    {
#elif defined(ZSLICE)
				      for(iz=ZSLICE;iz<ZSLICE+1;iz++)
					{
					  for(iy=0;iy<NY;iy++)
					    {
					      for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
						{
#elif defined(YZSLICE)
						  for(iy=NY/2;iy<NY/2+1;iy++)
						    {
						      for(iz=NZ/2;iz<NZ/2+1;iz++)
							{
							  for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
							    {
#else
							      for(iz=0;iz<NZ;iz++)
								{
								  for(iy=-gcly*NG;iy<NY+gcry*NG;iy++)
								    {
								      for(ix=-gclx*NG;ix<NX+gcrx*NG;ix++)
									{
#endif

									  /**************************/  
									  /**************************/  
									  /**************************/  


									  //within domain:
									  //if(if_indomain(ix,iy,iz)==0 && if_outsidegc(ix,iy,iz)==1) continue;
									  struct geometry geom;
									  fill_geometry(ix,iy,iz,&geom);

									  ldouble mx,my,mz,E,e,xx,yy,zz,phipot,xxx[4],dx[3],vv[10],a0,a1,a2,v1,v2,dphidx,v3,Tgas,Trad,v4,v5,v6,v7,v8,v9,v10,v11,v12,Fx,Fy,Fz;
									  ldouble gg[4][5],GG[4][5];
									  ldouble pp[NV],uu[NV];
									  int i,j;

									  v1=v2=v3=v4=v5=v6=v7=v8=v9=v10=v11=v12=0.;

									  ldouble xxvec[4],xxvecout[4];

									  //transforming code coordinates to output coordinates
									  get_xx(ix,iy,iz,xxvec);
						      
									  coco_N(xxvec,xxvecout,MYCOORDS,OUTCOORDS);
									  
									  xx=xxvecout[1];
									  yy=xxvecout[2];
									  zz=xxvecout[3];

									  if(codeprim==0)
									    {
#ifndef PRINTINSIDEBH						  
									      if((OUTCOORDS==KERRCOORDS || OUTCOORDS==BLCOORDS) && xx<1.*rhorizonBL) continue;
#endif
									    }



									  xxx[0]=t;
									  xxx[1]=xx;
									  xxx[2]=yy;
									  xxx[3]=zz;

									  pick_g(ix,iy,iz,gg);
									  pick_G(ix,iy,iz,GG);
									  ldouble gdet=gg[3][4];

									  dx[0]=get_size_x(ix,0)*sqrt(gg[1][1]);
									  dx[1]=get_size_x(iy,1)*sqrt(gg[2][2]);
									  dx[2]=get_size_x(iz,2)*sqrt(gg[3][3]);   


									  ldouble pporg[NV];
									  for(iv=0;iv<NV;iv++)
									    {
									      uu[iv]=get_u(u,iv,ix,iy,iz);
									      pp[iv]=get_u(p,iv,ix,iy,iz);
									      pporg[iv]=get_u(p,iv,ix,iy,iz);
									    }	 



						  
									  //ldouble tup[4][4],tlo[4][4];    
									  //ldouble eup[4][4],elo[4][4];

									  //to transform primitives between coordinates if necessary
									  ldouble ggout[4][5],GGout[4][5];
									  struct geometry geomout;
									  calc_g_arb(xxvecout,ggout,OUTCOORDS);
									  calc_G_arb(xxvecout,GGout,OUTCOORDS);
									  fill_geometry_arb(ix,iy,iz,&geomout,OUTCOORDS);



									  if(MYCOORDS!=OUTCOORDS && codeprim==0)
									    {

#ifdef RADIATION
									      trans_prad_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);
#endif
									      trans_pmhd_coco(pp, pp, MYCOORDS,OUTCOORDS, xxvec,&geom,&geomout);

									      //from now on geom,gg, GG, tup, etc. defined in OUTCOORDS!
									      fill_geometry_arb(ix,iy,iz,&geom,OUTCOORDS);
									      for(i=0;i<4;i++)
										for(j=0;j<5;j++)
										  { gg[i][j]=ggout[i][j]; GG[i][j]=GGout[i][j]; }

									      //calc_tetrades(gg,tup,tlo,OUTCOORDS);
									      //calc_ZAMOes(gg,eup,elo,OUTCOORDS);
									    }



									  ldouble rho=pp[0];
									  ldouble uint=pp[1];
									  ldouble S=pp[5];
									  ldouble pre=(GAMMA-1.)*uint;
									  gdet=gg[3][4];
									  ldouble ut=uu[0]/gdet/rho;
									  Tgas=pre*MU_GAS*M_PROTON/K_BOLTZ/rho;

									  ldouble vx=pp[2];
									  ldouble vy=pp[3];
									  ldouble vz=pp[4];
									  ldouble vrel[4]={0,vx,vy,vz};	

									  if(codeprim==0)
									    {
									      conv_vels(vrel,vrel,VELPRIM,VEL4,gg,GG);						  
									      //TODO
									      //tetrads sie zesraly
									      //trans2_cc2on(vrel,vrel,tup);
									      #ifdef BHDISK_PROBLEMTYPE
									      vrel[2]*=xx;
									      vrel[3]*=xx*sin(yy);
									      #endif

									      //outvel - ortonormal VEL4
									      vx=vrel[1];
									      vy=vrel[2];
									      vz=vrel[3];
									    }
						  
						  


#ifdef RADIATION

									  if(codeprim==0)
									    {

                                                                              //#ifdef RADOUTPUTINFF
									      //prad_lab2ff(pp,pp,&geom);
                                                                              //#endif

                                                                              //#ifdef RADOUTPUTINZAMO 
									      //prad_lab2on(pp,pp,&geom);
                                                                              //#endif

#ifdef RADOUTPUTVELS
									      ldouble vrelrad[4]={0,pp[7],pp[8],pp[9]};
									      conv_vels(vrelrad,vrelrad,VELPRIMRAD,VEL4,gg,GG);						  
									      //trans2_cc2on(vrelrad,vrelrad,tup);
									      //rad outvel - ortonormal VEL4
									      pp[7]=vrelrad[1];
									      pp[8]=vrelrad[2];
									      pp[9]=vrelrad[3];	  
#endif

									    }
						  


									  //**********************************************************************
									  //**********************************************************************
									  //**********************************************************************

									  //summing up multifluids
									  int irf;
									  E=Fx=Fy=Fz=0.;
									  
									    {

									      //						      irf=0;
									      E+=pp[EE];
									      Fx+=pp[FX];
									      Fy+=pp[FY];
									      Fz+=pp[FZ];
									      //						      break;
									    }

									  #ifndef EVOLVEPHOTONNUMBER
									  Trad=calc_LTE_TfromE(E);
									  #else
									  Trad=calc_ncompt_Thatrad_full(pp,&geomout);
									  #endif
									  if(E<EEFLOOR || isnan(E)) E=EEFLOOR;
									  
#endif //RADIATION
									  /******************/
									  /* extra lines to calculate v1...v4 from PROBLEMS/XXX/dump.c */

#ifdef PR_DUMP
#include PR_DUMP
#endif
									  /******************/


									  fprintf(fout1,"%.4e %.4e %.4e "
										  "%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e "
										  "%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e "
										  "%.9e %.9e ",
										  xx,     //1
										  yy,     //2
										  zz,     //3		      
										  uu[0],  //4
										  uu[1],  //5
										  uu[2],  //6
										  uu[3],  //7
										  uu[4],  //8
										  uu[5],  //9
#ifdef RADIATION										  
										  uu[EE],  //10
										  uu[FX],  //11
										  uu[FY],  //12
										  uu[FZ],  //13
#else
										  0.,
										  0.,
										  0.,
										  0.,
#endif										  
										  rho,    //14
										  uint, 
										  vx,     //16
										  vy,     //17
										  vz,     //18
										  S,      //19
#ifdef RADIATION
										  E,
										  Fx,
										  Fy,
										  Fz
#elif defined(MAGNFIELD)
										  pp[B1],
										  pp[B2],
										  pp[B3],
										  0.
#else
										  0.,
										  0.,
										  0.,
										  0.
#endif
										  );
									  
								

									  fprintf(fout1,"%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
										  v1,     //24
										  v2,     //25
										  v3,     //26 
										  v4,     //27
										  v5,     //28
										  v6,     //29 
										  v7,     //30 
										  v8,     //31 
										  v9,     //32 
										  v10,     //33
										  v11,     //34 
										  v12     //35
										  );

									}
								      fprintf(fout1,"\n");
								    }
								  fprintf(fout1,"\n\n");
								}
							      fflush(fout1);
							      fclose(fout1);

							      //calling gnuplot to produce gifs
#ifdef YZSLICE
							      convert_out2gif_1d(bufor,bufor2,nfout1,t);
#else
							      if(NY>3 || NZ>3)
								convert_out2gif_2d(bufor,bufor2,nfout1,t);
							      else
								convert_out2gif_1d(bufor,bufor2,nfout1,t);
#endif
  
							      
							      return 0;

							    }
                              
                              
/*********************************************/
/*********************************************/
/*********************************************/
/* prints restart files */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_restartfile(ldouble t, char* folder)
{
  #ifdef MPI
  fprint_restartfile_mpi(t,folder);

  #else 

  fprint_restartfile_bin(t,folder); 

  #endif
  
  return 0;
}


/*********************************************/
/*********************************************/

int //parallel output to a single file
fprint_restartfile_mpi(ldouble t, char* folder)
{
  #ifdef MPI
  char bufor[250];

  //header
  if(PROCID==0)
    {
       sprintf(bufor,"%s/res%04d.head",folder,nfout1);
       fout1=fopen(bufor,"w"); 
       //## nout time problem NX NY NZ
       sprintf(bufor,"## %5d %5d %10.6e %5d %5d %5d %5d\n",nfout1,nfout2,t,PROBLEM,TNX,TNY,TNZ);
       fprintf(fout1,"%s",bufor);
       fclose(fout1);
    }
  
  //maybe not needed
  MPI_Barrier(MPI_COMM_WORLD);

  //body
  sprintf(bufor,"%s/res%04d.dat",folder,nfout1);

  MPI_File cFile;
  MPI_Status status;
  MPI_Request req;

  int rc = MPI_File_open( MPI_COMM_WORLD, bufor, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &cFile );
  if (rc) {
    printf( "Unable to open/create file %s\n", bufor );fflush(stdout); exit(-1);
    }

  /***** first write all the indices ******/

  //  int indices[NX*NY*NZ*3];
  int *indices;
  if((indices = (int *)malloc(NX*NY*NZ*3*sizeof(int)))==NULL) my_err("malloc err. - fileop 0\n");
  
  int ix,iy,iz,iv;
  int gix,giy,giz;

  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  indices[ix*NY*NZ*3+iy*NZ*3+iz*3+0]=gix;
	  indices[ix*NY*NZ*3+iy*NZ*3+iz*3+1]=giy;
	  indices[ix*NY*NZ*3+iy*NZ*3+iz*3+2]=giz;
	}

  //set the initial location at each process for indices
  MPI_Offset pos;
  pos=PROCID*NX*NY*NZ*(3*sizeof(int));  
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 

  //write all indices
  MPI_File_write( cFile, indices, NX*NY*NZ*3, MPI_INT, &status );
  
  /***** then primitives in the same order ******/

  //now let's try manually
  pos=TNX*TNY*TNZ*(3*sizeof(int)) + PROCID*NX*NY*NZ*(NV*sizeof(ldouble)); 
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 

  //ldouble pout[NX*NY*NZ*NV];
  ldouble *pout;
  if((pout = (ldouble *)malloc(NX*NY*NZ*NV*sizeof(ldouble)))==NULL) my_err("malloc err. - fileop 1\n");
    
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	for(iv=0;iv<NV;iv++)
	  {
	    pout[ix*NY*NZ*NV+iy*NZ*NV+iz*NV+iv]=get_u(p,iv,ix,iy,iz);

#ifdef OVERWRITEENTROPYINRESTARTFILESWITHNEGEHEATING
	    if(iv==ENTR)
	      pout[ix*NY*NZ*NV+iy*NZ*NV+iz*NV+iv]=get_u_scalar(vischeatingnegebalance,ix,iy,iz);
#endif

	  }

  MPI_File_write( cFile, pout, NX*NY*NZ*NV, MPI_LDOUBLE, &status );  
  MPI_File_close( &cFile );

  //maybe not needed
  MPI_Barrier(MPI_COMM_WORLD);

  if(PROCID==0)
    {
      //sprintf(bufor,"cp %s/res%04d.dat %s/reslast.dat",folder,nfout1,folder);
      sprintf(bufor,"rm %s/reslast.dat",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.dat %s/reslast.dat",nfout1,folder);
      iv=system(bufor);
      //sprintf(bufor,"cp %s/res%04d.head %s/reslast.head",folder,nfout1,folder);
      sprintf(bufor,"rm %s/reslast.head",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.head %s/reslast.head",nfout1,folder);
      iv=system(bufor);
    }

  free (indices);
  free (pout);
#endif
  return 0;
}


/*********************************************/
/*********************************************/

int //serial binary output
fprint_restartfile_bin(ldouble t, char* folder)
{
  char bufor[250];
  
  //header
  if(PROCID==0)
    {
       sprintf(bufor,"%s/res%04d.head",folder,nfout1);
       fout1=fopen(bufor,"w"); 
       //## nout time problem NX NY NZ
       sprintf(bufor,"## %5d %5d %10.6e %5d %5d %5d %5d\n",nfout1,nfout2,t,PROBLEM,TNX,TNY,TNZ);
       fprintf(fout1,"%s",bufor);
       fclose(fout1);
    }

  //body
  sprintf(bufor,"%s/res%04d.dat",folder,nfout1);
  fout1=fopen(bufor,"wb"); 

  int ix,iy,iz,iv;
  int gix,giy,giz;
  ldouble pp[NV];
  //indices first
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  fwrite(&gix,sizeof(int),1,fout1);
	  fwrite(&giy,sizeof(int),1,fout1);
	  fwrite(&giz,sizeof(int),1,fout1);
	}

  //then, in the same order, primitives
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  ldouble ppout[NV];
	  PLOOP(iv)
	    ppout[iv]=get_u(p,iv,ix,iy,iz);
#ifdef RESTARTOUTPUTINBL
	  struct geometry geom,geomBL;
	  fill_geometry(ix,iy,iz,&geom);
	  fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	  trans_pall_coco(ppout, ppout, MYCOORDS,BLCOORDS, geom.xxvec,&geom,&geomBL);
#endif

#ifdef OVERWRITEENTROPYINRESTARTFILESWITHNEGEHEATING
	  ppout[ENTR]=get_u_scalar(vischeatingnegebalance,ix,iy,iz);
#endif

	  fwrite(ppout,sizeof(ldouble),NV,fout1);
	}

  fclose(fout1);

  if(PROCID==0)
    {
      //sprintf(bufor,"cp %s/res%04d.dat %s/reslast.dat",folder,nfout1,folder)
      sprintf(bufor,"rm %s/reslast.dat",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.dat %s/reslast.dat",nfout1,folder);
      iv=system(bufor);
      //sprintf(bufor,"cp %s/res%04d.head %s/reslast.head",folder,nfout1,folder);
      sprintf(bufor,"rm %s/reslast.head",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.head %s/reslast.head",nfout1,folder);
      iv=system(bufor);    
    }

  return 0;
}

							  
/*********************************************/
/*********************************************/
/*********************************************/
/* reads dump file */
/* puts conserved into the memory */
/* converts them to primitives */
/*********************************************/
/*********************************************/
/*********************************************/

int
fread_restartfile(int nout1, char* folder,ldouble *t)
{


  int ret;
   
  
  #ifdef MPI
  ret=fread_restartfile_mpi(nout1,folder,t);

  #else //no MPI 
  ret=fread_restartfile_bin(nout1,folder,t);  
  #endif
  
  return ret;
}


/*********************************************/
/*********************************************/

int 
fread_restartfile_bin(int nout1, char *folder, ldouble *t)
{
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz;
  char fname[400],fnamehead[400];
  if(nout1>=0)
    {
      sprintf(fname,"%s/res%04d.dat",folder,nout1);
      #ifdef MPI
      sprintf(fnamehead,"%s/../0/res%04d.head",folder,nout1);
      #else
      sprintf(fnamehead,"%s/res%04d.head",folder,nout1);
      #endif
    }
  else
    {
      sprintf(fname,"%s/reslast.dat",folder);
      #ifdef MPI
      sprintf(fnamehead,"%s/../0/reslast.head",folder);
      #else
      sprintf(fnamehead,"%s/reslast.head",folder);
      #endif
    }

  FILE *fdump;

  /***********/
  //header file
  fdump=fopen(fnamehead,"r");
  if(fdump==NULL) return 1; //request start from scratch

  //reading parameters, mostly time
  int intpar[6];
  ret=fscanf(fdump,"## %d %d %lf %d %d %d %d\n",&intpar[0],&intpar[1],t,&intpar[2],&intpar[3],&intpar[4],&intpar[5]);
  if(PROCID==0)
    printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
	 fname,intpar[0],*t,intpar[2],intpar[3],intpar[4],intpar[5]); 
  nfout1=intpar[0]+1; //global file no.
  nfout2=intpar[1]; //global file no. for avg
  fclose(fdump);

  /***********/
  //body file
  fdump=fopen(fname,"rb");
 
  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;
  //int indices[NX*NY*NZ][3];
  int **indices;
  if((indices = (int **)malloc(NX*NY*NZ*sizeof(int*)))==NULL) my_err("malloc err. - fileop 2\n");
  for(i=0;i<NX*NY*NZ;i++)
    if((indices[i]=(int *)malloc(3*sizeof(int)))==NULL) my_err("malloc err. - fileop 3\n");

  //first indices
  for(ic=0;ic<NX*NY*NZ;ic++)
    {
      ret=fread(&gix,sizeof(int),1,fdump);
      ret=fread(&giy,sizeof(int),1,fdump);
      ret=fread(&giz,sizeof(int),1,fdump);

      mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);

      indices[ic][0]=ix;
      indices[ic][1]=iy;
      indices[ic][2]=iz;
    }

  //then primitives
   for(ic=0;ic<NX*NY*NZ;ic++)
    {
      /*
      printf("NV %d\n", NV);
      printf("NRELBIN %d\n",NRELBIN);
      printf("NVHD %d\n",NVHD);
      printf("NV-NRELBIN %d\n", NV-NRELBIN);
      printf("NE %d\n", NEREL(0));
      getch();
      */
#ifdef RESTARTFROMNORELEL
      //if (indices[ic][0]==0 && indices[ic][1]==53) print_Nvector(pp,NV);
      int nvold=NV-NRELBIN;

#ifdef RESTARTFROMNORELEL_NOCOMPT
      nvold += 1;
#endif
      
      ldouble ppold[nvold]; 
      ret=fread(ppold,sizeof(ldouble),nvold,fdump);
      //if (indices[ic][0]==0 && indices[ic][1]==53) print_Nvector(ppold,nvold);

      int ie;
      for (ie=0; ie<NVHD-NRELBIN; ie++) pp[ie] = ppold[ie];
      for (ie=0; ie<NRELBIN; ie++) pp[NVHD-NRELBIN+ie] = 0.0; 
      for (ie=0; ie<NV-NVHD; ie++) pp[NVHD+ie] = ppold[NVHD-NRELBIN+ie];

      //if (indices[ic][0]==0 && indices[ic][1]==53) {print_Nvector(pp,NV); getch();}

#else
      ret=fread(pp,sizeof(ldouble),NV,fdump);
#endif

      ix=indices[ic][0];
      iy=indices[ic][1];
      iz=indices[ic][2];

      fill_geometry(ix,iy,iz,&geom);

#ifdef RESTARTOUTPUTINBL
      struct geometry geomBL;
      fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
      trans_pall_coco(pp, pp, BLCOORDS,MYCOORDS, geomBL.xxvec,&geomBL,&geom);
#endif

      /*
#ifdef CONVERTLOGTONOLOGWHENRESTARTING
      ldouble gamma=GAMMA;
      ldouble gammam1=gamma-1.;
      ldouble ugas = pow((pow(pp[RHO],1./(gammam1)+1.)*exp(pp[ENTR]/pp[RHO])),gammam1)/(gamma-1.);
      pp[ENTR]=calc_Sfromu(pp[RHO],ugas,ix,iy,iz);

      ldouble ne=pp[RHO]/MU_E/M_PROTON;
      gamma=GAMMAE;
      //ldouble pe= pow((pow(pp[RHO],1./(gamma-1.)+1.)*exp(pp[ENTRE]/pp[RHO]/K_BOLTZ)),gamma-1.); 
      ldouble pe= pow((pow(pp[RHO],1./(gamma-1.)+1.)*exp(pp[ENTRE]/pp[RHO])),gamma-1.);
      ldouble Te=pe/K_BOLTZ/ne; 
      pp[ENTRE]=calc_S2fromrhoT(pp[RHO],Te,ELECTRONS);

      ldouble ni=pp[RHO]/MU_I/M_PROTON;
      gamma=GAMMAI;
      //ldouble pi= pow((pow(pp[RHO],1./(gamma-1.)+1.)*exp(pp[ENTRI]/pp[RHO]/K_BOLTZ)),gamma-1.); 
      ldouble pi= pow((pow(pp[RHO],1./(gamma-1.)+1.)*exp(pp[ENTRI]/pp[RHO])),gamma-1.);
      ldouble Ti=pi/K_BOLTZ/ni; 
      pp[ENTRI]=calc_S2fromrhoT(pp[RHO],Ti,IONS);

#endif
      */
      
      #ifdef CONSISTENTGAMMA
      ldouble Te,Ti;
      calc_PEQ_Teifrompp(pp,&Te,&Ti,ix,iy,iz);
      ldouble gamma=calc_gammaintfromTei(Te,Ti);
      set_u_scalar(gammagas,ix,iy,iz,gamma);
      #endif
      
      p2u(pp,uu,&geom);

#ifdef OVERWRITEENTROPYINRESTARTFILESWITHNEGEHEATING //recover entropy
      if(!doingpostproc)
	{
	  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);
	  uu[ENTR]=pp[ENTR]*(uu[RHO]/pp[RHO]);
	}
#endif


      //saving primitives
      for(iv=0;iv<NV;iv++)    
	{
	  set_u(u,iv,ix,iy,iz,uu[iv]);
	  set_u(p,iv,ix,iy,iz,pp[iv]);
	}
    }

  for(i=0;i<NX*NY*NZ;i++)
    free(indices[i]);
  free(indices);
  
  fclose(fdump);

  return 0;
}


/*********************************************/
/*********************************************/

int 
fread_restartfile_mpi(int nout1, char *folder, ldouble *t)
{
  #ifdef MPI
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz,tix,tiy,tiz;
  char fname[400],fnamehead[400];

  if(nout1>=0)
    {
      sprintf(fname,"%s/res%04d.dat",folder,nout1);
      sprintf(fnamehead,"%s/res%04d.head",folder,nout1);
    }
  else
    {
      sprintf(fname,"%s/reslast.dat",folder);
      sprintf(fnamehead,"%s/reslast.head",folder);
    }

  FILE *fdump;

  /***********/
  //header file
 
  fdump=fopen(fnamehead,"r");
  if(fdump==NULL) 
    {
      return 1; //request start from scratch
    }
  //reading parameters, mostly time
  int intpar[6];
  ret=fscanf(fdump,"## %d %d %lf %d %d %d %d\n",&intpar[0],&intpar[1],t,&intpar[2],&intpar[3],&intpar[4],&intpar[5]);
  if(PROCID==0)
  printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
	 fname,intpar[0],*t,intpar[2],intpar[3],intpar[4],intpar[5]); 
  nfout1=intpar[0]+1; //global file no.
  nfout2=intpar[1]; //global file no. for avg
  fclose(fdump);
    
  //maybe not needed
  MPI_Barrier(MPI_COMM_WORLD);

  /***********/
  //body file
  struct geometry geom;
  ldouble uu[NV],pp[NV],ftemp;

  MPI_File cFile;
  MPI_Status status;
  MPI_Request req;

  int rc = MPI_File_open( MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &cFile );
  if (rc) {
    printf( "Unable to open/create file %s\n", fname );fflush(stdout); exit(-1);
    }

  /***** first read all the indices ******/

 int nvold;
#ifdef RESTARTFROMNORELEL
  nvold=NV-NRELBIN;
#ifdef RESTARTFROMNORELEL_NOCOMPT
      nvold += 1;
#endif
#else
  nvold=NV;
#endif

  //first read the indices pretending to be a single process
  int *indices;
  if((indices = (int *)malloc(NX*NY*NZ*3*sizeof(int)))==NULL) my_err("malloc err. - fileop 5\n");
  int len=NX*NY*NZ;

  ldouble *pout;
  if((pout=(ldouble *)malloc(NX*NY*NZ*nvold*sizeof(ldouble)))==NULL) my_err("malloc err. - fileop 7\n");

  //set the initial location
  int procid=PROCID;
  MPI_Offset pos;

#ifdef RESTARTGENERALINDICES
  for(procid=0;procid<NTX*NTY*NTZ;procid++)
#endif
    {
      pos=procid*NX*NY*NZ*(3*sizeof(int));  

      MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
      //read them
      MPI_File_read( cFile, indices, 3*len, MPI_INT, &status );

      //convert to local
      for(ic=0;ic<len;ic++)
	{
	  gix=indices[ic*3+0];
	  giy=indices[ic*3+1];
	  giz=indices[ic*3+2];
	  mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);
	  indices[ic*3+0]=ix;
	  indices[ic*3+1]=iy;
	  indices[ic*3+2]=iz;
	}

      /***** then primitives in the same order ******/

      pos=TNX*TNY*TNZ*(3*sizeof(int)) + procid*NX*NY*NZ*(nvold*sizeof(ldouble)); 
      MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
      MPI_File_read( cFile, pout, len*nvold, MPI_LDOUBLE, &status );
 
      //rewriting to p
      int ppos;
      for(ic=0;ic<len;ic++)
	{
	  ix=indices[ic*3+0];
	  iy=indices[ic*3+1];
	  iz=indices[ic*3+2];

	  ppos=ic*nvold;

	  if(if_indomain(ix,iy,iz))
	    {
	      fill_geometry(ix,iy,iz,&geom);

#ifdef RESTARTFROMNORELEL
          int ie;
          for (ie=0; ie<8; ie++) set_u(p,ie,ix,iy,iz,pout[ppos+ie]);
	  for (ie=0; ie<(NV-NVHD); ie++) set_u(p,NVHD+ie,ix,iy,iz,pout[ppos+8+ie]);
          for (ie=0; ie<NRELBIN; ie++) set_u(p,8+ie,ix,iy,iz, 0.0); //set relel bins to zero 
#else
	  PLOOP(iv) set_u(p,iv,ix,iy,iz,pout[ppos+iv]);
#endif

	      
#ifdef CONSISTENTGAMMA
	      ldouble Te,Ti;
	      calc_PEQ_Teifrompp(&get_u(p,0,ix,iy,iz),&Te,&Ti,ix,iy,iz);
	      ldouble gamma=calc_gammaintfromTei(Te,Ti);
	      set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif


	      p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);


	      
#ifdef OVERWRITEENTROPYINRESTARTFILESWITHNEGEHEATING //recover entropy
	      if(!doingpostproc)
		{
		  ldouble pentre=calc_Sfromu(pp[RHO],pp[UU],ix,iy,iz);
		  ldouble uentre=pentre*(get_u(u,RHO,ix,iy,iz)/get_u(p,RHO,ix,iy,iz));
		  set_u(p,ENTR,ix,iy,iz,pentre);
		  set_u(u,ENTR,ix,iy,iz,uentre);
		}
#endif
	      

	    }
	}
    }


  MPI_File_close( &cFile );
  MPI_Barrier(MPI_COMM_WORLD);
  free(indices);
  free(pout);
#endif
  return 0;
}
							    

///////////////////////////////////////////////////////////////

int
fread_restartfile_mpi_org(int nout1, char *folder, ldouble *t)
{
  #ifdef MPI
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz,tix,tiy,tiz;
  char fname[400],fnamehead[400];

  if(nout1>=0)
    {
      sprintf(fname,"%s/res%04d.dat",folder,nout1);
      sprintf(fnamehead,"%s/res%04d.head",folder,nout1);
    }
  else
    {
      sprintf(fname,"%s/reslast.dat",folder);
      sprintf(fnamehead,"%s/reslast.head",folder);
    }


  FILE *fdump;

  /***********/
  //header file
 
  fdump=fopen(fnamehead,"r");
  if(fdump==NULL) 
    {
      return 1; //request start from scratch
    }
  //reading parameters, mostly time
  int intpar[6];
  ret=fscanf(fdump,"## %d %d %lf %d %d %d %d\n",&intpar[0],&intpar[1],t,&intpar[2],&intpar[3],&intpar[4],&intpar[5]);
  if(PROCID==0)
  printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
	 fname,intpar[0],*t,intpar[2],intpar[3],intpar[4],intpar[5]); 
  nfout1=intpar[0]+1; //global file no.
  nfout2=intpar[1]; //global file no. for avg
  fclose(fdump);
    
  //maybe not needed
  MPI_Barrier(MPI_COMM_WORLD);

  /***********/
  //body file
  struct geometry geom;
  ldouble uu[NV],pp[NV],ftemp;

  MPI_File cFile;
  MPI_Status status;
  MPI_Request req;

  int rc = MPI_File_open( MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &cFile );
  if (rc) {
    printf( "Unable to open/create file %s\n", fname );fflush(stdout); exit(-1);
    }

  /***** first read all the indices ******/


  //first read the indices
  #ifdef RESTARTGENERALINDICES
  //int indices[TNX*TNY*TNZ*3];
  int *indices;
  if((indices=(int *)malloc(TNX*TNY*TNZ*3*sizeof(int)))==NULL) my_err("malloc err. - fileop 4\n");
  int len=TNX*TNY*TNZ;
  #else
  //int indices[NX*NY*NZ*3];
  int *indices;
  if((indices = (int *)malloc(NX*NY*NZ*3*sizeof(int)))==NULL) my_err("malloc err. - fileop 5\n");
  int len=NX*NY*NZ;
  #endif

  //set the initial location
  MPI_Offset pos;
  #ifdef RESTARTGENERALINDICES
  pos=0;
  #else
  pos=PROCID*NX*NY*NZ*(3*sizeof(int));  
  #endif
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
  //read them
  MPI_File_read( cFile, indices, 3*len, MPI_INT, &status );

  //convert to local
  for(ic=0;ic<len;ic++)
    {
      gix=indices[ic*3+0];
      giy=indices[ic*3+1];
      giz=indices[ic*3+2];
      mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);
      indices[ic*3+0]=ix;
      indices[ic*3+1]=iy;
      indices[ic*3+2]=iz;
    }

  /***** then primitives in the same order ******/

 int nvold;
#ifdef RESTARTFROMNORELEL
  nvold=NV-NRELBIN;
#ifdef RESTARTFROMNORELEL_NOCOMPT
      nvold += 1;
#endif
#else
  nvold=NV;
#endif

  //new location in the second block
  #ifdef RESTARTGENERALINDICES
  pos=TNX*TNY*TNZ*(3*sizeof(int));
  //  ldouble pout[TNX*TNY*TNZ*NV];
  ldouble *pout;
  if((pout=(ldouble *)malloc(TNX*TNY*TNZ*nvold*sizeof(ldouble)))==NULL) my_err("malloc err. - fileop 6\n");
  #else
  pos=TNX*TNY*TNZ*(3*sizeof(int)) + PROCID*NX*NY*NZ*(nvold*sizeof(ldouble)); 
  //ldouble pout[NX*NY*NZ*NV];
  ldouble *pout;
  if((pout=(ldouble *)malloc(NX*NY*NZ*nvold*sizeof(ldouble)))==NULL) my_err("malloc err. - fileop 7\n");
  #endif
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
  //so far manually
  //ANDREW
  MPI_File_read( cFile, pout, len*nvold, MPI_LDOUBLE, &status );
 
  //rewriting to p
  int ppos;
  for(ic=0;ic<len;ic++)
    {
      ix=indices[ic*3+0];
      iy=indices[ic*3+1];
      iz=indices[ic*3+2];

      ppos=ic*nvold;

      if(if_indomain(ix,iy,iz))
	{
	  fill_geometry(ix,iy,iz,&geom);

	  //ANDREW
          //printf("%d %d %d \n", NV, NVHD, NRELBIN);
#ifdef RESTARTFROMNORELEL
          //ANDREW CHECK RESTART
          int ie;
          for (ie=0; ie<8; ie++) set_u(p,ie,ix,iy,iz,pout[ppos+ie]);
	  for (ie=0; ie<(NV-NVHD); ie++) set_u(p,NVHD+ie,ix,iy,iz,pout[ppos+8+ie]);
          for (ie=0; ie<NRELBIN; ie++) set_u(p,8+ie,ix,iy,iz, 0.0); //set relel bins to zero 
#else
	  PLOOP(iv) set_u(p,iv,ix,iy,iz,pout[ppos+iv]);
#endif
	  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	}
    }

  MPI_File_close( &cFile );
  MPI_Barrier(MPI_COMM_WORLD);
  free(indices);
  free(pout);
#endif
  return 0;
}
							    

/*********************************************/
/*********************************************/
/*********************************************/
/* prints avg files */
/*********************************************/
/*********************************************/
/*********************************************/

int
fprint_avgfile(ldouble t, char* folder,char* prefix)
{
  #ifdef MPI

  fprint_avgfile_mpi(t,folder,prefix);

  #else

  fprint_avgfile_bin(t,folder,prefix); 

  #endif
  
  return 0;
}


/*********************************************/
/*********************************************/

int //parallel output to a single file
fprint_avgfile_mpi(ldouble t, char* folder, char* prefix)
{
  #ifdef MPI
  char bufor[250];
  //header
  if(PROCID==0)
    {
      sprintf(bufor,"%s/%s%04d.head",folder,prefix,nfout2);
      fout1=fopen(bufor,"w"); 
      sprintf(bufor,"## %5d %10.6e %10.6e %10.6e\n",nfout2,t-avgtime,t,avgtime);
      fprintf(fout1,"%s",bufor);
      fclose(fout1);
    }

  //body
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfout2);

  MPI_File cFile;
  MPI_Status status;
  MPI_Request req;

 
  int rc = MPI_File_open( MPI_COMM_WORLD, bufor, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &cFile );
  if (rc) {
    printf( "Unable to open/create file %s\n", bufor );fflush(stdout); exit(-1);
    }

  /***** first write all the indices ******/

  int nz=NZ;
  int tnz=TNZ;
#ifdef AVGOUTPUT
  if(AVGOUTPUT==2)
    {
      nz=1;
      tnz=1;
    }
#endif

  int *indices;
  if((indices = (int *)malloc(NX*NY*nz*3*sizeof(int)))==NULL) my_err("malloc err. - fileop 8\n");
  
  int ix,iy,iz,iv;
  int gix,giy,giz;

  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<nz;iz++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  indices[ix*NY*nz*3+iy*nz*3+iz*3+0]=gix;
	  indices[ix*NY*nz*3+iy*nz*3+iz*3+1]=giy;
	  indices[ix*NY*nz*3+iy*nz*3+iz*3+2]=giz;
	}

  //set the initial location at each process for indices
  MPI_Offset pos;
  pos=PROCID*NX*NY*nz*(3*sizeof(int));  
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 

  //write all indices
  MPI_File_write( cFile, indices, NX*NY*nz*3, MPI_INT, &status );
  
  /***** then primitives in the same order ******/

  //now let's try manually
  pos=TNX*TNY*tnz*(3*sizeof(int)) + PROCID*NX*NY*nz*((NV+NAVGVARS)*sizeof(ldouble)); 
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 

  ldouble *pout;
  if((pout=(ldouble *) malloc(NX*NY*nz*(NV+NAVGVARS)*sizeof(ldouble)))==NULL) my_err("malloc err. - fileop 9\n");
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<nz;iz++)
	for(iv=0;iv<(NV+NAVGVARS);iv++)
	  pout[ix*NY*nz*(NV+NAVGVARS)+iy*nz*(NV+NAVGVARS)+iz*(NV+NAVGVARS)+iv]=get_uavg(pavg,iv,ix,iy,iz);

  MPI_File_write( cFile, pout, NX*NY*nz*(NV+NAVGVARS), MPI_LDOUBLE, &status );

  free(pout);
  free(indices);

  MPI_File_close( &cFile );

#endif
  return 0;
}


/*********************************************/
/*********************************************/

int //serial binary output
fprint_avgfile_bin(ldouble t, char* folder,char *prefix)
{
  char bufor[250];
  
  //header
  if(PROCID==0)
    {
      sprintf(bufor,"%s/%s%04d.head",folder,prefix,nfout2);
      fout1=fopen(bufor,"w"); 
      sprintf(bufor,"## %5d %10.6e %10.6e %10.6e\n",nfout2,t-avgtime,t,avgtime);
      fprintf(fout1,"%s",bufor);
      fclose(fout1);
    }

  //body
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfout2);
  fout1=fopen(bufor,"wb"); 



  int ix,iy,iz,iv;
  int gix,giy,giz;
  ldouble pp[NV];
  //indices first
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  fwrite(&gix,sizeof(int),1,fout1);
	  fwrite(&giy,sizeof(int),1,fout1);
	  fwrite(&giz,sizeof(int),1,fout1);
	}

  //then, in the same order, primitives
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	{
	  fwrite(&get_uavg(pavg,0,ix,iy,iz),sizeof(ldouble),NV+NAVGVARS,fout1);
	}

  fclose(fout1);


  return 0;
}


///////////////////////////////////////////////////////////////

int //serial binary output
fprint_avgfile_bin_old(ldouble t, char* folder,char *prefix)
{
  char bufor[250];
  
  //header
  if(PROCID==0)
    {
      sprintf(bufor,"%s/%s%04d.head",folder,prefix,nfout2);
      fout1=fopen(bufor,"w"); 
      sprintf(bufor,"## %5d %10.6e %10.6e %10.6e\n",nfout2,t-avgtime,t,avgtime);
      fprintf(fout1,"%s",bufor);
      fclose(fout1);
    }

  //body
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfout2);
  fout1=fopen(bufor,"wb"); 

  int ix,iy,iz,iv;
  int gix,giy,giz;
  ldouble pp[NV];
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  fwrite(&gix,sizeof(int),1,fout1);
	  fwrite(&giy,sizeof(int),1,fout1);
	  fwrite(&giz,sizeof(int),1,fout1);
	  fwrite(&get_uavg(pavg,0,ix,iy,iz),sizeof(ldouble),NV+NAVGVARS,fout1);
	}

  fclose(fout1);

  return 0;
}

							  
/*********************************************/
/*********************************************/
/*********************************************/
/* reads dump file */
/* puts conserved into the memory */
/* converts them to primitives */
/*********************************************/
/*********************************************/
/*********************************************/

int
fread_avgfile(int nout1, char* base,ldouble *pavg, ldouble *dt,ldouble *t)
{
  char bufor[250];

  #ifdef MPI
  my_err("fread_avgfile should not be used with MPI\n");
  exit(1);
  
  //fread_avgfile_mpi(nout1,base,pavg,dt,t);
  
  #else //no MPI

  fread_avgfile_bin(nout1,base,pavg,dt,t);

  #endif
  
  return 0;
}

                              
/*********************************************/
/*********************************************/

/*********************************************/
/*********************************************/

int 
fread_avgfile_bin(int nout1, char *base, ldouble *pavg, ldouble *dt, ldouble *t)
{
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz;
  char fname[40],fnamehead[40];

  /*
  sprintf(fname,"%s/avg%04d.dat",folder,nout1);
#ifdef MPI
  sprintf(fnamehead,"%s/../0/avg%04d.head",folder,nout1);
#else
  sprintf(fnamehead,"%s/avg%04d.head",folder,nout1);
#endif
  */

  printf("%s%04d.dat\n",base,nout1);
  printf("%s%04d.head\n",base,nout1);
  
  sprintf(fname,"%s%04d.dat",base,nout1);
  sprintf(fnamehead,"%s%04d.head",base,nout1);

  FILE *fdump;

  /***********/
  //header file
  fdump=fopen(fnamehead,"r");

  //reading parameters, mostly time
  int intpar[5];
  ldouble ldpar[5];
  ret=fscanf(fdump,"## %d %lf %lf %lf\n",&intpar[0],&ldpar[0],&ldpar[1],&ldpar[2]);
  if(PROCID==0) printf("avg file (%s) read no. %d at times: %.6e to %.6e (dt=%.6e)\n",
	 fname,intpar[0],ldpar[0],ldpar[1],ldpar[2]); 
  
  *t=.5*(ldpar[0]+ldpar[1]);
  *dt=ldpar[2];
  fclose(fdump);
 
  /***********/
  //body file

  fdump=fopen(fname,"rb");

  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;
  //int indices[NX*NY*NZ][3];
  int **indices;
  if((indices = (int **)malloc(NX*NY*NZ*sizeof(int*)))==NULL) my_err("malloc err. - fileop 10\n");
  for(i=0;i<NX*NY*NZ;i++)
    if((indices[i]=(int *)malloc(3*sizeof(int)))==NULL) my_err("malloc err. - fileop 11\n");

  //to mark unfilled slots
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	set_uavg(pavg,RHO,ix,iy,iz,-1.);

  //first indices
  for(ic=0;ic<NX*NY*NZ;ic++)
    {
      ret=fread(&gix,sizeof(int),1,fdump);
      ret=fread(&giy,sizeof(int),1,fdump);
      ret=fread(&giz,sizeof(int),1,fdump);

      mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);

      if(ix<0 || ix>=NX) {ix=0; printf("bad idx in avg: %d %d | %d %d %d\n",ic,NX*NY*NZ,ix,iy,iz);}
      if(iy<0 || iy>=NY) iy=0;
      if(iz<0 || iz>=NZ) iz=0;

      indices[ic][0]=ix;
      indices[ic][1]=iy;
      indices[ic][2]=iz;
    }

  //then averages
   for(ic=0;ic<NX*NY*NZ;ic++)
    {
      ix=indices[ic][0];
      iy=indices[ic][1];
      iz=indices[ic][2];

      ret=fread(&get_uavg(pavg,0,ix,iy,iz),sizeof(ldouble),NV+NAVGVARS,fdump);

#ifdef CONSISTENTGAMMA
      ldouble gamma = 1. + get_uavg(pavg,AVGPGAS,ix,iy,iz)/get_uavg(pavg,UU,ix,iy,iz);
      set_u_scalar(gammagas,ix,iy,iz,gamma);
#endif
    }

  for(i=0;i<NX*NY*NZ;i++)
    free(indices[i]);
  free(indices);

  fclose(fdump);

  return 0;
}

                              
/*********************************************/
/*********************************************/

int 
fread_avgfile_mpi(int nout1, char *folder,ldouble *pavg, ldouble *dt,ldouble *t)
{
#ifdef MPI
  my_err("fread_avgfile should not be used with MPI\n");
  exit(1);

  /*
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz;
  char fname[40],fnamehead[40];

  sprintf(fname,"%s/avg%04d.dat",folder,nout1);
  sprintf(fnamehead,"%s/avg%04d.head",folder,nout1);

  FILE *fdump;

  //header file
  fdump=fopen(fnamehead,"r");

  //reading parameters, mostly time
  int intpar[5];
  ldouble ldpar[5];
  ret=fscanf(fdump,"## %d %lf %lf %lf\n",&intpar[0],&ldpar[0],&ldpar[1],&ldpar[2]);
  if(PROCID==0) printf("avg file (%s) read no. %d at times: %.6e to %.6e (dt=%.6e)\n",
	 fname,intpar[0],ldpar[0],ldpar[1],ldpar[2]); 
  *t=.5*(ldpar[0]+ldpar[1]);
  *dt=ldpar[2];
  fclose(fdump);

  //body file
  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;

  MPI_File cFile;
  MPI_Status status;
  MPI_Request req;

  int rc = MPI_File_open( MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &cFile );
  if (rc) {
    printf( "Unable to open/create file %s\n", fname );fflush(stdout); exit(-1);
    }

  //set the initial location
  MPI_Offset pos;
  if(PROCID==0) pos=0;
  else
    pos=PROCID*NX*NY*NZ*(3*sizeof(int)+(NV+NAVGVARS)*sizeof(ldouble));
  MPI_File_seek( cFile, pos, MPI_SEEK_SET ); 
  
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	{
	  mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
	  
	  MPI_File_read( cFile, &gix, 1, MPI_INT, &status );
	  MPI_File_read( cFile, &giy, 1, MPI_INT, &status );
	  MPI_File_read( cFile, &giz, 1, MPI_INT, &status );

	  mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);

	  MPI_File_read( cFile, &get_uavg(pavg,0,ix,iy,iz), NV+NAVGVARS, MPI_LDOUBLE, &status );

	}

  MPI_File_close( &cFile );
  */
#endif

  return 0;
}

                              
/*********************************************/
/*********************************************/
/*********************************************/
/* wrapper for coordinate output */
/*********************************************/
/*********************************************/
/*********************************************/

int fprint_coordfile(char* folder,char* prefix)
{
#if (COORDOUTPUT==1)
  fprint_coordBL(folder,prefix);
#endif
#if (COORDOUTPUT==2)
  fprint_coordBL_shell(folder,prefix);
#endif

  return 0;
}


/*********************************************/
/*********************************************/
/*********************************************/
/* prints BL coordinates,  */
/*********************************************/
/*********************************************/
/*********************************************/
                              
int fprint_coordBL(char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%sBL.dat",folder,prefix);
   FILE* fout1=fopen(bufor,"w");

   int ix,iy,iz,iv;
   ldouble pp[NV];
   for(iz=0;iz<NZ;iz++)
     {
       for(iy=0;iy<NY;iy++)
	 {
	   for(ix=0;ix<NX;ix++)
	     {
	       struct geometry geom,geomBL;
	       fill_geometry(ix,iy,iz,&geom);
	       fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	       ldouble r=geomBL.xx;
	       ldouble th=geomBL.yy;
	       ldouble ph=geomBL.zz;
	     
	       fprintf(fout1,"%d %d %d ",ix,iy,iz);

	       fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph);

	       fprintf(fout1,"\n");
	     }
	 }
     }

   fflush(fout1);
   fclose(fout1);

   return 0;
 }

                              
/*********************************************/
/*********************************************/
/*********************************************/
/* prints BL polar/azimuthal coordinates on a shell  */
/*********************************************/
/*********************************************/
/*********************************************/
                              
int fprint_coordBL_shell(char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%sBL.dat",folder,prefix);
   FILE* fout1=fopen(bufor,"w");

   int ix,iy,iz,iv;
   ldouble pp[NV];

   ix=NX-1;

   for(iz=0;iz<NZ;iz++)
     {
       for(iy=0;iy<NY;iy++)
	 {
	   struct geometry geom,geomBL;
	   fill_geometry(ix,iy,iz,&geom);
	   fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	   ldouble r=geomBL.xx;
	   ldouble th=geomBL.yy;
	   ldouble ph=geomBL.zz;
	     
	   fprintf(fout1,"%d %d %d ",ix,iy,iz);

	   fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph);

	   fprintf(fout1,"\n");
	 }
     }

   fflush(fout1);
   fclose(fout1);

   return 0;
 }

                              
/*********************************************/
/*********************************************/
/*********************************************/
/* wrapper for simple output */
                              
int fprint_simplefile(ldouble t, int nfile, char* folder,char* prefix)
{
#if (SIMOUTPUT==1)
  fprint_simplecart(t,nfile,folder,prefix);
#endif

#if (SIMOUTPUT==2)
  fprint_simplesph(t,nfile,folder,prefix);
//#ifdef SIMOUTPUT_GU
//  fprint_simplesph_GU(t,nfile,folder,prefix);
//#endif
#endif

  //#if (SIMOUTPUT==3)
  //fprint_simplebondi(t,nfile,folder,prefix);
  //#endif

#if (SIMOUTPUT==4)
  char prefixext[100];
  sprintf(prefixext,"%sext",prefix);
  fprint_simpleextended(t,nfile,folder,prefixext);
#endif

  return 0;
}

                              
/*********************************************/
/*********************************************/
/*********************************************/
/* prints in ASCII indices, cart coordinates,*/
/* primitives, velocities in cartesian       */
/*********************************************/
/*********************************************/
/*********************************************/
                              
int fprint_simplecart(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);
   fout1=fopen(bufor,"w");
  
   //header
   //## nout time problem NX NY NZ
   fprintf(fout1,"## %d %e %d %d %d %d\n",nfout1,t,PROBLEM,NX,NY,NZ);

   /***********************************/  
   /** writing order is fixed  ********/  
   /***********************************/  
 
   int ix,iy,iz,iv;
   ldouble pp[NV];
   int nz=NZ;
#if (PROBLEM==120)
   nz=1;
#endif
   for(iz=0;iz<nz;iz++)
     {
       for(iy=0;iy<NY;iy++)
	 {
	   for(ix=0;ix<NX+2;ix++)
	     {
	       for(iv=0;iv<NV;iv++)
		 {
                   pp[iv]=get_u(p,iv,ix,iy,iz);
		 }

	       struct geometry geom,geomcart,geomout,geomsph;
	       fill_geometry(ix,iy,iz,&geom);
	       fill_geometry_arb(ix,iy,iz,&geomcart,MINKCOORDS);
	       fill_geometry_arb(ix,iy,iz,&geomout,OUTCOORDS);
	       fill_geometry_arb(ix,iy,iz,&geomsph,SPHCOORDS);

	       ldouble dx[3];
	       dx[0]=get_size_x(ix,0);
	       dx[1]=get_size_x(iy,1);
	       dx[2]=get_size_x(iz,2);
	       ldouble gdet=geom.gdet;
	       ldouble volume=dx[0]*dx[1]*dx[2]*gdet;
	       trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomout);
	       ldouble rho=rhoGU2CGS(pp[RHO]);
	       #ifdef SIMOUTPUTINTERNAL
	       rho=pp[RHO];
	       #endif
	       ldouble temp=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);
	       //ldouble tracer=pp[TRA];
	       ldouble vel[4]={0,pp[VX],pp[VY],pp[VZ]};	
	       ldouble vx,vy,vz;
	       conv_vels(vel,vel,VELPRIM,VEL4,geomout.gg,geomout.GG);

	       //TODO
	       //tetrads sie zesraly
	       //trans2_cc2on(vel,vel,geomout.tup); //ANDREW replaced with below 

	       vx=vel[1];
	       vy=vel[2];
	       vz=vel[3];
	       
	       //transform to cartesian
	       if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS   || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS ||
		   MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS || MYCOORDS==MKS3COORDS || MYCOORDS==MSPH1COORDS || MYCOORDS==MKER1COORDS)
		 {
		   ldouble r=geomsph.xx;
		   ldouble th=geomsph.yy;
		   ldouble ph=geomsph.zz;

		   vel[2]*=r;
		   vel[3]*=r*sin(th);
		    
		   vx = sin(th)*cos(ph)*vel[1] 
		     + cos(th)*cos(ph)*vel[2]
		     - sin(ph)*vel[3];

		   vy = sin(th)*sin(ph)*vel[1] 
		     + cos(th)*sin(ph)*vel[2]
		     + cos(ph)*vel[3];

		   vz = cos(th)*vel[1] 
		     - sin(th)*vel[2];
		 }
	     
	       fprintf(fout1,"%d %d %d ",ix,iy,iz); //1-3

	       fprintf(fout1,"%.5e %.5e %.5e ",geomcart.xx,geomcart.yy,geomcart.zz);//4-6

	       fprintf(fout1,"%.5e %.5e ",rho,temp);//7-8

	       fprintf(fout1,"%.5e %.5e %.5e ",vx,vy,vz);//9-11

	       fprintf(fout1,"%.5e ",volume);//12

	       #ifdef RADIATION
	       ldouble Rtt,ehat,Rij[4][4];
	       ldouble ugas[4],Fx,Fy,Fz;
	       if(doingavg==0)
		{
		  calc_ff_Rtt(pp,&Rtt,ugas,&geomout);
		  ehat=-Rtt;  
		  calc_Rij(pp,&geomout,Rij); //calculates R^munu in OUTCOORDS
		  indices_2221(Rij,Rij,geomout.gg);	      							  
		}
	      else
		{
		  ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  int i,j;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz);
		}

	       //transform to cartesian
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
		{
		  ldouble r=geomsph.xx;
		  ldouble th=geomsph.yy;
		  ldouble ph=geomsph.zz;

		  Rij[2][0]*=r;
		  Rij[3][0]*=r*sin(th);

		  Fx = sin(th)*cos(ph)*Rij[1][0] 
		    + cos(th)*cos(ph)*Rij[2][0]
		    - sin(ph)*Rij[3][0];

		  Fy = sin(th)*sin(ph)*Rij[1][0] 
		    + cos(th)*sin(ph)*Rij[2][0]
		    + cos(ph)*Rij[3][0];

		  Fz = cos(th)*Rij[1][0] 
		    - sin(th)*Rij[2][0];
		}
	       
	      fprintf(fout1,"%.5e %.5e %.5e %.5e ",endenGU2CGS(ehat),fluxGU2CGS(Fx),fluxGU2CGS(Fy),fluxGU2CGS(Fz));//13-16
#endif

#if (PROBLEM==115 || PROBLEM==135) //SHOCKELECTRONTEST
	      ldouble uugas;
	      ldouble Tg,Te,Ti;
	      uugas=pp[UU];
	      Tg=calc_PEQ_Teifrompp(pp,&Te,&Ti,ix,iy,iz); //temperatures after explicit
	      
	      /**************/
	      //electrons
	      /**************/
	      ldouble ne=rho/MU_E/M_PROTON; //number density of photons and electrons
	      ldouble pe=K_BOLTZ*ne*Te;
	      ldouble gammae=GAMMAE;
#ifdef CONSISTENTGAMMA
#ifndef FIXEDGAMMASPECIES
	      gammae=calc_gammaintfromtemp(Te,ELECTRONS);
#endif
#endif
	      ldouble ue=pe/(gammae-1.);

	      /**************/
	      //ions
	      /**************/
	      ldouble ni=rho/MU_I/M_PROTON; //number density of photons and electrons
	      ldouble pi=K_BOLTZ*ni*Ti;
	      ldouble gammai=GAMMAI;
#ifdef CONSISTENTGAMMA
#ifndef FIXEDGAMMASPECIES
	      gammai=calc_gammaintfromtemp(Ti,IONS);
#endif
#endif
	      ldouble ui=pi/(gammai-1.);

	      ue = calc_ufromSerho(pp[ENTRE], rho, ELECTRONS,ix,iy,iz);
	      ui = calc_ufromSerho(pp[ENTRI], rho, IONS,ix,iy,iz);

	    if(ix==550)
	    {
	    
	      printf("energy frac diff ix=550 in output: %e \n", (pp[UU] - ue - ui)/pp[UU]); 
	    }

	    ldouble gammagas=calc_gammagas(pp, ix, iy, iz);
	    gammagas=pick_gammagas(ix,iy,iz);
	    fprintf(fout1,"%e %e %e %.5e %.5e %.5e %.5e %.5e %.5e %.5e",uugas,ui,ue,get_u_scalar(vischeating,ix,iy,iz),get_u_scalar(vischeatingnegebalance,ix,iy,iz),gammagas,Te,Ti,gammae,gammai);//17-21 with rad, 13-17 without
#endif //PROBLEM==115

	      fprintf(fout1,"\n");
	     }
	 }
     }

   fflush(fout1);
   fclose(fout1);

   return 0;
 }

                              
/*********************************************/
/*********************************************/
/*********************************************/
/* prints in ASCII & BL coordinates,  */
/*********************************************/
/*********************************************/
/*********************************************/
                              
int fprint_simplesph(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   #ifdef GRTRANSSIMOUTPUT
   sprintf(bufor,"%s/%s%04d_grtnew.dat",folder,prefix,nfile);
   #else
   sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);
   #endif
   fout1=fopen(bufor,"w");
  
   /***********************************/  
   /** writing order is fixed  ********/  
   /***********************************/  
 
   int ix,iy,iz,iv;
   int iix;
   ldouble pp[NV],phi,tausca,tauscar,lorentz,vel[4],vcon[4],tauscarloc;
   int nz=NZ;
   struct geometry geom,geomBL;
 
#if (PROBLEM==120)
   nz=1;
#endif

   int xmin=-2;
  
   //ANDREW header for grtrans
#ifdef GRTRANSSIMOUTPUT
   for(iix=-2;iix<NX;iix++)
     {
        fill_geometry_arb(iix,NY/2,NZ/2,&geomBL,OUTCOORDS);
	if(geomBL.xx>=rhorizonBL)
	  {
	    xmin=iix;
	    break;
	  }
     }
       
     //Can't have NaN in grtrans files!
   if(NZ==1)
     fprintf(fout1,"%.5e %5d %5d %.5e %.5e ",t,NX+2,NY,BHSPIN,MASS);
   else
     fprintf(fout1,"%.5e %5d %5d %5d %.5e %.5e ",t,NX+2,NY,NZ,BHSPIN,MASS);     
   fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e\n",MKSR0,MKSH0,MKSMY1,MKSMY2,MKSMP0);
#ifdef RELELECTRONS
   fprintf(fout1,"%5d %.5e %.5e\n",NRELBIN, RELGAMMAMIN, RELGAMMAMAX);
#endif
#endif
   
#ifdef RELELECTRONS //ANDREW array for finding gamma break
  int ie;
  ldouble gammapbrk[NRELBIN];
  for(ie=0; ie<NRELBIN; ie++) gammapbrk[ie] = pow(relel_gammas[ie], RELEL_HEAT_INDEX + 0.5);
#endif  

   for(iz=0;iz<nz;iz++)
     {
       #ifndef RAD_INTEGRATION
       for(iix=-2;iix<NX;iix++)
	 {
	   phi=0.;
	   tausca=0.;
       #else //Start from outermost radial cell
       for(iy=0;iy<NY;iy++)
	 {
	   tauscar=0.;
       #endif
           #ifndef RAD_INTEGRATION
	   for(iy=0;iy<NY;iy++)
	   {
           #else //Start from outermost radial cell
           for(iix=NX-1;iix>-3;iix--)
	   {
           #endif

               ix=iix;
 
	       fill_geometry(ix,iy,iz,&geom);
	       fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);

#ifdef SIMOUTPUTWITHINDTHETA 
	       if(fabs(geomBL.yy-M_PI/2)>SIMOUTPUTWITHINDTHETA)
		 continue;
#endif

	       ldouble r=geomBL.xx;
	       ldouble th=geomBL.yy;
	       ldouble ph=geomBL.zz;

#ifdef GRTRANSSIMOUTPUT //ANDREW 2D ONLY
	       ldouble x1=geom.xx;
	       ldouble x2=geom.yy;
	       ldouble x3=geom.zz;
	       //if(NZ==1)
	       //{
               //fprintf(fout1,"%d %d ",ix,iy); //(1-2)
	       //fprintf(fout1,"%.5e %.5e ",x1,x2); //(3-4)
	       //fprintf(fout1,"%.5e %.5e ",r,th); //(5-6)
               //}
	       //else
	       //{
	       fprintf(fout1,"%d %d %d ",ix,iy,iz); //(1-3)
	       fprintf(fout1,"%.5e %.5e %.5e ",x1,x2,x3); //(4-6)
	       fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph); //(7-9)
               //}
	       
               //ANDREW fill values below horizon with values right above horizon
               //if(ix<xmin)
	       if(r<rhorizonBL)
	       {
                 ix=xmin;
                 fill_geometry(ix,iy,iz,&geom);
	         fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);
               }
#else	     
	       fprintf(fout1,"%d %d %d ",ix,iy,iz); //(1-3)
	       fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph); //(4-6)
#endif

	      for(iv=0;iv<NV;iv++)
	      {
		  if(doingavg)
		    pp[iv]=get_uavg(pavg,iv,ix,iy,iz);
		  else
		    pp[iv]=get_u(p,iv,ix,iy,iz);
	      }

	      ldouble dxph[3],dx[3],xx1[4],xx2[4];
	      xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_x(iz,2);
	      xx2[0]=0.;xx2[1]=get_xb(ix+1,0);xx2[2]=get_x(iy,1);xx2[3]=get_x(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
	      dx[0]=fabs(xx2[1]-xx1[1]);
	      xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_x(iz,2);
	      xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_xb(iy+1,1);xx2[3]=get_x(iz,2);
	      coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
	      dx[1]=fabs(xx2[2]-xx1[2]);
	      xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_xb(iz,2);
	      xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_x(iy,1);xx2[3]=get_xb(iz+1,2);
	      coco_N(xx1,xx1,MYCOORDS,OUTCOORDS);
	      coco_N(xx2,xx2,MYCOORDS,OUTCOORDS);
	      dx[2]=fabs(xx2[3]-xx1[3]);

	      dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
	      dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
	      dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);

	       ldouble gdet=geom.gdet;
	       ldouble volume=gdet*get_size_x(ix,0)*get_size_x(iy,1)*get_size_x(iz,2);
	       
               //ANDREW volume is just gdet for grtrans output 
#ifdef GRTRANSSIMOUTPUT
               volume=gdet;
#endif
	       ldouble rho,uint,pgas,temp,bsq,bcon[4],bcov[4],utcon[4],utcov[4],ucon[4],ucov[4],Tij[4][4],Tij22[4][4];
	       ldouble Ti,Te;

	       int i,j;

	       ldouble gamma=GAMMA;
#ifdef CONSISTENTGAMMA
	       gamma=pick_gammagas(ix,iy,iz);
#endif
	     
	       if(doingavg)
		 {
                   //ANDREW need pp for some relel computations below
		   trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomBL);

		   rho=get_uavg(pavg,RHO,ix,iy,iz);
		   uint=get_uavg(pavg,UU,ix,iy,iz);
		   pgas=get_uavg(pavg,AVGPGAS,ix,iy,iz);
		   //temperature
		   temp=calc_PEQ_Tfromprho(pgas,rho,ix,iy,iz);

		   vel[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   vel[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   vel[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   vel[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
                   for(i=0;i<4;i++) vcon[i]=vel[i];
                   lorentz = fabs(vel[0])/sqrt(fabs(geomBL.GG[0][0]));

		   //temp=get_uavg(pavg,AVGTGAS,ix,iy,iz);
		   utcon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);

                   int ii,jj;
                   for(ii=0;ii<4;ii++)
		     for(jj=0;jj<4;jj++)
		       Tij[ii][jj]=get_uavg(pavg,AVGTIJ(ii,jj),ix,iy,iz);                 

		   indices_2122(Tij,Tij22,geomBL.gg);  
                 
                   //ANDREW NORMALIZE u^0
#ifdef GRTRANSSIMOUTPUT
                   fill_utinucon(utcon,geomBL.gg,geomBL.GG);
		   indices_21(utcon,utcov,geomBL.gg); 
#endif
                   pp[RHO]=rho;
		   pp[UU]=uint;
#ifdef MAGNFIELD
		   bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
		   bcon[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
		   bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
		   bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
		   bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);

#ifdef GRTRANSSIMOUTPUT
                  //ANDREW NORMALIZE b^0 to be orthogonal with u^\mu
		  bcon[0]=-dot3nr(bcon,utcov)/utcov[0];
		  indices_21(bcon,bcov,geomBL.gg);

                  //ANDREW NORMALIZE b^mu to be equal to B^2
		  ldouble alphanorm = bsq/dotB(bcon,bcov);
		  if(alphanorm<0.) my_err("alpha.lt.0 in b0 norm !!\n");
                  for(i=0;i<4;i++)
		  {
		   bcon[i]*=sqrt(alphanorm);
		  }
#endif
#endif

		  
#ifdef EVOLVEELECTRONS
		  ldouble pe,pi;
		  pe=get_uavg(pavg,AVGPE,ix,iy,iz);
		  pi=get_uavg(pavg,AVGPI,ix,iy,iz);
		  //electrons
		  ldouble ne=get_uavg(pavg,RHO,ix,iy,iz)/MU_E/M_PROTON; 
		  //ions
		  ldouble ni=get_uavg(pavg,RHO,ix,iy,iz)/MU_I/M_PROTON; 

                  #ifdef RELELECTRONS
                  #ifndef NORELELAVGS
                  ne=get_uavg(pavg,AVGNETH,ix,iy,iz);            
                  #endif
                  #endif 
		  Te=pe/K_BOLTZ/ne;
		  Ti=pi/K_BOLTZ/ni;

		  //write these temperatures into the primitives as corresponding entropies
                  //#ifdef RELELENTROPY
                  //pp[ENTRE]=calc_S4fromnT(ne,Te,ELECTRONS);
                  //pp[ENTRI]=calc_S4fromnT(ni,Ti,IONS);
                  //#else
		  ldouble rhoeth=ne*MU_E*M_PROTON;
		  pp[ENTRE]=calc_SefromrhoT(rhoeth,Te,ELECTRONS);
		  pp[ENTRI]=calc_SefromrhoT(rho,Ti,IONS);
                  //#endif
		 
#endif                
		  //ldouble Tet,Tit;
		  //temp=calc_PEQ_Teifrompp(pp,&Tet,&Tit,geomBL.ix,geomBL.iy,geomBL.iz);

		 }
	       else //on the go from the primitives
		 { 
		   trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomBL);

		   rho=pp[0];
		   uint=pp[1];
		   pgas=(gamma-1.)*uint;
                   ldouble vel[4],vcon[4],tauscarloc;
                   
                   //obtain 4 velocity
	           vel[1]=pp[VX];
	           vel[2]=pp[VY];
	           vel[3]=pp[VZ];
	           conv_vels(vel,vel,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
                   //printf("%f %f %f %f %f\n",lorentz,vel[0],vel[1],vel[2],vel[3]);
                   for(i=0;i<4;i++) vcon[i]=vel[i];
                   lorentz = fabs(vel[0])/sqrt(fabs(geomBL.GG[0][0]));

                   calc_Tij(pp,&geomBL,Tij22);
		   indices_2221(Tij22,Tij,geomBL.gg);
		   //utcon[1]=pp[2];
		   //utcon[2]=pp[3];
		   //utcon[3]=pp[4];
           //conv_vels_both(utcon,utcon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
           calc_ucon_ucov_from_prims(pp, &geomBL, utcon, ucov);
           temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
           temp=calc_PEQ_Teifrompp(pp,&Te,&Ti,geomBL.ix,geomBL.iy,geomBL.iz);
#ifdef MAGNFIELD
		   //calc_bcon_prim(pp,bcon,&geomBL);
		   //indices_21(bcon,bcov,geomBL.gg);
		   //bsq = dotB(bcon,bcov);
           calc_bcon_bcov_bsq_from_4vel(pp, utcon, ucov, &geomBL, bcon, bcov, &bsq);
#endif

		 }

	       rho=rhoGU2CGS(pp[RHO]);
#ifdef OUTPUTINGU
	       rho=pp[RHO];
#endif
	     
#ifdef RHOLABINSIM
	       rho*=utcon[0];
#endif

	       fprintf(fout1,"%.5e %.5e ",rho,temp); //(7-8) //or (10-11) for grtrans 
	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",utcon[0],utcon[1],utcon[2],utcon[3]); //(9-12) //or (12-15) for grtrans

	       //fprintf(fout1,"%.5e ",volume);// (13)
	       //ANDREW AA CHANGED THIS
	       fprintf(fout1,"%.5e ", volume);// (13) //(16) for grtran

	       //ldouble dBdt[4];
	       //estimate_Bgrowth_battery(ix,iy,iz,dBdt);
	       //fprintf(fout1,"%.5e ",dBdt[3]);// (13)

	       ldouble ehat=0.;
#ifdef RADIATION
	       ldouble Rtt,Rij[4][4],Rij22[4][4],vel[4];
	       ldouble ugas[4],Fx,Fy,Fz;
	       ldouble Gi[4],Giff[4]={0.,0.,0.,0.};
	       ldouble Gic[4],Gicff[4]={0.,0.,0.,0.};
	       ldouble Trad;

	       ldouble CoulombCoupling=0.;
	       if(doingavg==0) 
		{
		  calc_ff_Rtt(pp,&Rtt,ugas,&geomBL);
		  ehat=-Rtt;  
		  calc_Rij(pp,&geomBL,Rij22); //calculates R^munu in OUTCOORDS
		  indices_2221(Rij22,Rij,geomBL.gg);

		  //four fource
		  calc_Gi(pp,&geomBL,Giff,0.0,0,0); //ANDREW 0 for fluid frame
		  //boost2_lab2ff(Gi,Giff,pp,geomBL.gg,geomBL.GG);
#if defined(COMPTONIZATION) || defined(EVOLVEPHOTONNUMBER)
		  ldouble kappaes=calc_kappaes(pp,&geomBL);
		  //test - directly in ff
		  vel[1]=vel[2]=vel[3]=0.; vel[0]=1.;
		  calc_Compt_Gi(pp,&geomBL,Gicff,ehat,Te,kappaes,vel);
		  //boost2_lab2ff(Gic,Gicff,pp,geomBL.gg,geomBL.GG);
#endif 

#ifdef EVOLVEPHOTONNUMBER //the color temperature of radiation
  Trad = calc_ncompt_Thatrad_full(pp,&geomBL);
#endif
	
#ifdef EVOLVEELECTRONS	 
#ifndef  SKIPCOULOMBCOUPLING
  CoulombCoupling=calc_CoulombCoupling(pp,&geomBL); 
#endif
#endif
		}
	      else
		{
		  ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  int i,j;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 
		  for(j=0;j<4;j++)
		    Giff[j]=get_uavg(pavg,AVGGHAT(j),ix,iy,iz);
		    
		  indices_2122(Rij,Rij22,geomBL.gg);  
                  //ANDREW recompute if Giff in avg accidentally saved as lab frame
                  #ifdef SIMOUTPUT_GILAB2FF
                  //boost2_lab2ff(Giff,Giff,pp,geomBL.gg,geomBL.GG); //ANDREW avg Gff already in OUTCOORDS
		  calc_Gi(pp,&geomBL,Giff,0.0,0,0); //ANDREW 0 for fluid frame, 2 for fluid frame thermal only
                  #endif

#if defined(COMPTONIZATION) || defined(EVOLVEPHOTONNUMBER)
		  for(j=0;j<4;j++)
		    Gicff[j]=get_uavg(pavg,AVGGHATCOMPT(j),ix,iy,iz);
#endif			  
		  Trad=calc_ncompt_Thatrad_fromEN(ehat,get_uavg(pavg,AVGNFHAT,ix,iy,iz));
		}
	       
	       //flux
	       Fx=Rij[1][0];
	       Fy=Rij[2][0];
	       Fz=Rij[3][0];
#ifdef OUTPUTINGU
	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",ehat,Fx,Fy,Fz); //(14) - (17) 
	       fprintf(fout1,"%.5e %.5e ",Giff[0],Gicff[0]); //(18)-(19)
               fprintf(fout1,"%.5e %.5e ",ehat, Trad); //(20), (21)         
#else
	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",endenGU2CGS(ehat),fluxGU2CGS(Fx),fluxGU2CGS(Fy),fluxGU2CGS(Fz)); //(14) - (17)  //17-20 for grtrans
	       ldouble conv=kappaGU2CGS(1.)*rhoGU2CGS(1.)*endenGU2CGS(1.)*CCC; //because (cE-4piB) in non-geom
	       fprintf(fout1,"%.5e %.5e ",Giff[0]*conv,Gicff[0]*conv); //(18)-(19)  //21-22 for grtrans
               fprintf(fout1,"%.5e %.5e ",CoulombCoupling*conv, Trad); //(20), (21) //23-24 for grtrans      
#endif

#endif //RADIATION

	       ldouble gammam1=gamma-1.;

	       ldouble betarad=ehat/3./(pgas);
	       //ldouble betarad=ehat/3./(pgas + bsq/2.); - Brandon: Should this be the new calculation?

	       //magn. field components
#ifdef MAGNFIELD
	       if(doingavg==0) 
		{
		  //calc_bcon_prim(pp,bcon,&geomBL);
		  //indices_21(bcon,bcov,geomBL.gg);
		  //bsq = dotB(bcon,bcov);
          calc_ucon_ucov_from_prims(pp, &geomBL, utcon, ucov);
          calc_bcon_bcov_bsq_from_4vel(pp, utcon, ucov, &geomBL, bcon, bcov, &bsq);
		}
	      else
		{
		  bsq=get_uavg(pavg,AVGBSQ,ix,iy, iz);
		  bcon[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
		  bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
		  bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
		  bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);

                  #ifdef GRTRANSSIMOUTPUT
                  //ANDREW NORMALIZE b^0 to be orthogonal with u^\mu
		  bcon[0]=-dot3nr(bcon,utcov)/utcov[0];
		  indices_21(bcon,bcov,geomBL.gg);

                  //ANDREW NORMALIZE b^mu to be equal to B^2
		  ldouble alphanorm = bsq/dotB(bcon,bcov);
		  if(alphanorm<0.) my_err("alpha.lt.0 in b0 norm !!\n");
                  for(i=0;i<4;i++)
		  {
		   bcon[i]*=sqrt(alphanorm);
		  }
                  #endif
		  
		}	       

	       ldouble betamag = bsq/2./(pgas + ehat/3.	+ bsq/2.);
	       //ldouble betamag = bsq/2./(pgas + ehat/3.);	- Brandon: Should this be the calculation?

	       //to CGS!
	       #ifndef OUTPUTINGU
	       bsq=endenGU2CGS(bsq);
	       ldouble scaling=endenGU2CGS(1.);
	       for(i=0;i<4;i++)
		 {
                   //ANDREW made scaling the same for i=0
                   /*
		   if(i==0) bcon[i]*=scaling;
		   else*/
		   bcon[i]*=sqrt(scaling);
		 }
               #endif

	       //magnetic flux parameter
	       int iphimin,iphimax;
	       iphimin=0;
	       iphimax=TNY-1;

#if defined(CORRECT_POLARAXIS) || defined(CORRECT_POLARAXIS_3D)
	       iphimin=NCCORRECTPOLAR; 
#ifndef HALFTHETA
	       iphimax=TNY-NCCORRECTPOLAR-1;
#endif
#endif
	       //
	       //if(iy>=iphimin && iy<=iphimax)
		// {
		//   phi+=geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1)*2.*M_PI;
		//   //phi+=bcon[1]*dxph[1]*2.*M_PI;
		//   tausca+=rho*0.34*lenGU2CGS(dxph[1]); //rho already converted to cgs
		// }
	       if(iy>=iphimin && iy<=iphimax)
		 {
                   #ifndef RAD_INTEGRATION
		   phi+=geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1)*2.*M_PI;
		   //phi+=bcon[1]*dxph[1]*2.*M_PI;
		   tausca+=rho*0.34*lenGU2CGS(dxph[1]); //rho already converted to cgs
                   #endif
		 }
               #ifdef RAD_INTEGRATION
               #ifdef RADIATION
	       if(ix>=0)
		 {
		   tauscarloc = vcon[0]*(1.-abs(vcon[1]))*calc_kappaes(pp,&geomBL); //rho already converted to cgs
                   if(ix==NX-1)
		   {
		       tauscar=tauscarloc*dxph[0];
      	           }
                   else
                   {
                       tauscar += tauscarloc*dxph[0];
                   }
		 }
               #endif
               #endif

               // ANDREW no b^2, just b[0] in output for grtrans
#ifdef GRTRANSSIMOUTPUT	       
	       fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e %.5e ",bcon[0],bcon[1],bcon[2],bcon[3],bsq,phi); //(14) - (19) or (22) - (27) if radiation included (+3 with grtrans 3d)
	       fprintf(fout1,"%.5e %.5e ",betamag,betarad); // (28) - (29) when rad and magn field on, (20) - (21) with no radiation (+3 with grtrans 3d)

#else    
	       fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e %.5e ",bsq,bcon[1],bcon[2],bcon[3],phi,betamag); //(14) - (19) or (22) - (27) if radiation included (+3 with  grtrans 3d)
               #ifndef RAD_INTEGRATION
	       fprintf(fout1,"%.5e %.5e ",betarad,tausca); // (28) - (29) when rad and magn field on, (20) - (21) with no radiation (+3 with grtrans 3d)
               #else
	       fprintf(fout1,"%.5e %.5e ",betarad,tauscar); // (28) - (29) when rad and magn field on, (20) - (21) with no radiation (+3 with grtrans 3d)
               #endif

#endif
#endif //MAGNFIELD
	       

#ifdef EVOLVEELECTRONS
	       
	       fprintf(fout1,"%.5e %.5e %.5e ",Te, Ti, gamma); // (30) - (32) when rad and magn field on, (22) - (24) with no radiation (+3 with grtrans 3d)
	       
	       ldouble vischeat,pe,ue,gammae,ne,tempeloc;
	       gammae=GAMMAE;
	       if(doingavg)
		 {
		   vischeat=get_uavg(pavg,AVGVISCHEATING,ix,iy,iz);
		   pe=get_uavg(pavg,AVGPE,ix,iy,iz);
		   ne=calc_thermal_ne(pp); 
		   tempeloc=pe/K_BOLTZ/ne;
	           gammae=GAMMAE;
                   #ifdef CONSISTENTGAMMA
		   #ifndef FIXEDGAMMASPECIES
		   gammae=calc_gammaintfromtemp(tempeloc,ELECTRONS);
                   #endif
		   #endif
		   ue=pe/(gammae-1.);
		 }
	       else
		 {
		   vischeat=get_u_scalar(vischeating,ix,iy,iz);
                   //#ifdef RELELENTROPY
                   //ue=calc_ufrromS4n(pp[ENTRE],calc_thermal_ne(pp),ELECTRONS,ix,iy,iz); 
                   //#else
		   ldouble rhoeth=calc_thermal_ne(pp)*M_PROTON*MU_E;
		   ue=calc_ufromSerho(pp[ENTRE],rhoeth,ELECTRONS,ix,iy,iz); 
                   //#endif
		 }


	       //ANDREW
	       //in  avg vischeat was averaged as du, not du/dtau
	       //recompute dt and use that as an estimate
	       //copied from problem.c              
               #ifdef DIVIDEVISCHEATBYDT
               dt=get_u_scalar(cell_dt,ix,iy,iz); //individual time step
	       ldouble dtau = dt/vel[0];
	       vischeat/=dtau;
               #endif
               

	       //skip for now
	       /*
	       if(doingavg) 
		 {
		   //convert erg/cm3 to erg/cm3/s 
		   //in dt: time over which avg integrated - but that only if ACCUMULATEVISCOUS HEATING!!!
		   vischeat/=dt;
		   //and then rescale per orbital time to obtain again erg/cm3 - visc heating per dynamical time
		   ldouble Pk=2.*M_PI*sqrt(geomBL.xx*geomBL.xx*geomBL.xx);
		   vischeat*=Pk;
		 }
	       */

	       ldouble meandeltae=get_uavg(pavg,AVGVISCHEATINGTIMESDELTAE,ix,iy,iz)/get_uavg(pavg,AVGVISCHEATING,ix,iy,iz);

	       #ifndef OUTPUTINGU
               ue = endenGU2CGS(ue);
               vischeat=endenGU2CGS(vischeat)*timeCGS2GU(1.);
               #endif
	       
	       fprintf(fout1,"%.5e %.5e %.5e ",meandeltae,vischeat,ue); //(33) - (35) if rad and magn on, (25) - (27) if not (+3 with grtrans)

	       //ANDREW rel electron quantities
	       //(36) -- if rad and magn on, (28) -- if not (+3 with grtrans included)
#ifdef RELELECTRONS
   ldouble nrelel, urelel, G0relel, gammabrk;
               if(doingavg==0)
	       {
		  urelel=calc_relel_uint(pp);
		  nrelel=calc_relel_ne(pp);
               }
	       else
	       {
	       #ifndef NORELELAVGS
	          nrelel=get_uavg(pavg,AVGNRELEL,ix,iy,iz);
                  urelel=get_uavg(pavg,AVGURELEL,ix,iy,iz);
               #else
		  urelel=calc_relel_uint(pp);
		  nrelel=calc_relel_ne(pp);
               #endif
	       }

               G0relel = -1.*calc_relel_G0_fluidframe(pp, &geomBL, 0.0, 0); //ANDREW - fluid frame

	       gammabrk=RELEL_INJ_MIN;
	       //absolute maximum of g^4*n for g > RELGAMMAMIN
	       ldouble nbrk=pp[NEREL(0)]*gammapbrk[0];
	       ldouble nbrk2;
	       for(ie=1;ie<NRELBIN;ie++)
	       {
		 if (relel_gammas[ie] < RELEL_INJ_MIN)
		 {
		   gammabrk=RELEL_INJ_MIN;
		   nbrk =  pp[NEREL(ie)]*gammapbrk[ie];
		 }

	         else 
	         {
               	   nbrk2 =  pp[NEREL(ie)]*gammapbrk[ie];
		   if(nbrk2 > nbrk)
	           {
		     nbrk=nbrk2;
	             gammabrk=relel_gammas[ie];
         	   }
	         }
		}
	       
	       	 #ifndef OUTPUTINGU
                 nrelel = numdensGU2CGS(nrelel);
                 urelel = endenGU2CGS(urelel);
		 G0relel = G0relel*kappaGU2CGS(1.)*rhoGU2CGS(1.)*endenGU2CGS(1.)*CCC; //because (cE-4piB) in non-geom
		 //G0relel = endenGU2CGS(G0relel)*timeCGS2GU(1.);
		 //printf("CGS: %e \n",nrelel);
                 #endif

		 fprintf(fout1,"%.5e %.5e %.5e %.5e",urelel,nrelel,G0relel,gammabrk); //(36) - (39) if rad and magn on , (38)-(41) with grtrans

	       ldouble nbin;
	       for(ie=0; ie<NRELBIN; ie++)
	       {
                 if(doingavg)
                   nbin=get_uavg(pavg,NEREL(ie),ix,iy,iz);
		 else
		   nbin=pp[NEREL(ie)];
		 #ifndef OUTPUTINGU
                 nbin = numdensGU2CGS(nbin);
		 #endif
	         fprintf(fout1," %.5e ",nbin); //(40)-- if rad and magn on, (42)-- with grtrans 3d
	       }	 	       
#endif //RELELECTRONS
#endif //EVOLVEELECTRONS

               //Full output - Leave off for output for HEROIC
               #ifndef GRTRANSSIMOUTPUT
               #ifdef FULLOUTPUT
               fprintf(fout1," %.5e ",lorentz); //(30) if rad and magn on and no relele, (40) with relele
	
               fprintf(fout1," %.5e %.5e %.5e %.5e ",ucon[0],ucon[1],ucon[2],ucon[3]); // 31-34  
               fprintf(fout1," %.5e %.5e %.5e %.5e ",ucov[0],ucov[1],ucov[2],ucov[3]); // 35-38
               #ifdef MAGNFIELD 
               fprintf(fout1," %.5e %.5e %.5e %.5e ",bcon[0],bcon[1],bcon[2],bcon[3]); // 39-42   
               fprintf(fout1," %.5e %.5e %.5e %.5e ",bcov[0],bcov[1],bcov[2],bcov[3]); // 43-46   
               #endif

               int iii,jjj;
               //output T^munu (columns 46-62)
               for(iii=0;iii<4;iii++)
               {
                 for(jjj=0;jjj<4;jjj++)
                 {
                   fprintf(fout1," %.5e ",Tij22[iii][jjj]);
                 }
               }
               //output T^mu_nu (columns 63-78)
               for(iii=0;iii<4;iii++)
               {
                 for(jjj=0;jjj<4;jjj++)
                 {
                   fprintf(fout1," %.5e ",Tij[iii][jjj]);
                 }
               }
               #ifdef RADIATION
               //output R^munu (columns 79-94)
               for(iii=0;iii<4;iii++)
               {
                 for(jjj=0;jjj<4;jjj++)
                 {
                   fprintf(fout1," %.5e ",Rij22[iii][jjj]);
                 }
               }
               //output R^mu_nu (columns 95-110)
               for(iii=0;iii<4;iii++)
               {
                 for(jjj=0;jjj<4;jjj++)
                 {
                   fprintf(fout1," %.5e ",Rij[iii][jjj]);
                 }
               }
               #endif
               #endif
               #endif

	       fprintf(fout1,"\n");
	     }
	 }
     }

   fflush(fout1);
   fclose(fout1);

   return 0;
 }


/*********************************************/
/*********************************************/
/*********************************************/
/* prints in ASCII & BL coordinates,  */
/*********************************************/
/*********************************************/
/*********************************************/
                              
int fprint_simpleextended(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

   fout1=fopen(bufor,"w");
  
   /***********************************/  
   /** writing order is fixed  ********/  
   /***********************************/  
 
   int ix,iy,iz,iv;
   ldouble pp[NV];
   for(iz=0;iz<NZ;iz++)
     {
       for(iy=0;iy<NY;iy++)
	 {
	   for(ix=-NG;ix<NX;ix++)
	     {
	       for(iv=0;iv<NV;iv++)
		 {
		   pp[iv]=get_u(p,iv,ix,iy,iz);
		 }
	       struct geometry geom,geomBL;
	       fill_geometry(ix,iy,iz,&geom);
	       fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);

	       ldouble dx[3];
	       dx[0]=get_size_x(ix,0);
	       dx[1]=get_size_x(iy,1);
	       dx[2]=get_size_x(iz,2);
	       ldouble gdet=geom.gdet;
	       ldouble volume=dx[0]*dx[1]*dx[2]*gdet;

	       
	       ldouble rho,uint,temp,bsq,bcon[4],bcov[4],utcon[4],ucon[4],ucov[4],Tij[4][4],Tij22[4][4],Ehatucon[4];
	       ldouble Rij[4][4],Ehat,temprad,Giff[4],Gi[4],Titmagn[4],Titkin[4],Tituint[4];
	       int i,j;
	      if(doingavg)
		{
		  rho=get_uavg(pavg,RHO,ix,iy,iz);
		  uint=get_uavg(pavg,UU,ix,iy,iz);
		  //temperature
		  
		  

		  ldouble gamma=GAMMA;
		  //#ifdef SKIP_CONSISTENTGAMMA
		  //gamma=calc_gammafromtemp(temp,GAS);  //this is not properly averaged, I should dump GAMMA*stuff to avg
		  //#endif
		  ldouble gammam1=gamma-1.;

		  temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
		  temp=get_uavg(pavg,AVGTGAS,ix,iy,iz);
#ifdef MAGNFIELD
		  bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
		  bcon[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
		  bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
		  bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
		  bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);
#endif
		  utcon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		  pp[RHO]=rho;
		  pp[UU]=uint;
		  
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Tij[i][j]=get_uavg(pavg,AVGRHOUCONUCOV(i,j),ix,iy,iz)
			+ gamma*get_uavg(pavg,AVGUUUCONUCOV(i,j),ix,iy,iz)
			+ get_uavg(pavg,AVGBSQUCONUCOV(i,j),ix,iy,iz)
			- get_uavg(pavg,AVGBCONBCOV(i,j),ix,iy,iz); 

		  for(i=0;i<4;i++)
		    {
		      Titmagn[i] = get_uavg(pavg,AVGBSQUCONUCOV(i,0),ix,iy,iz)
			- get_uavg(pavg,AVGBCONBCOV(i,0),ix,iy,iz); 

		      Titkin[i]  = get_uavg(pavg,AVGRHOUCONUCOV(i,0),ix,iy,iz);

		      Tituint[i]  =  get_uavg(pavg,AVGUUUCONUCOV(i,0),ix,iy,iz);
		    }

#ifdef RADIATION  
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 

		  Ehat = get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  for(iv=0;iv<4;iv++)
		    Giff[iv]=get_uavg(pavg,AVGGHAT(iv),ix,iy,iz);
		  for(iv=1;iv<4;iv++)
		    Ehatucon[iv]=get_uavg(pavg,AVGEHATUCON(iv),ix,iy,iz); 
		    
		  conv_vels(utcon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
		  pp[VX]=ucon[1];
		  pp[VY]=ucon[2];
		  pp[VZ]=ucon[3];
		  boost2_ff2lab(Giff,Gi,pp,geomBL.gg,geomBL.GG);
		  indices_21(Gi,Gi,geomBL.gg);
#endif
		}
	      else //on the go from the primitives
		{ 
		  trans_pall_coco(pp, pp, MYCOORDS,BLCOORDS, geom.xxvec,&geom,&geomBL);

		  rho=pp[0];
		  uint=pp[1];
		  //utcon[1]=pp[2];
		  //utcon[2]=pp[3];
		  //utcon[3]=pp[4];
          //conv_vels_both(utcon,utcon,ucov,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
          calc_ucon_ucov_from_prims(pp, &geomBL, utcon, ucov);
		  temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
#ifdef MAGNFIELD
		  //calc_bcon_prim(pp,bcon,&geomBL);
		  //indices_21(bcon,bcov,geomBL.gg);
		  //bsq = dotB(bcon,bcov);
          calc_bcon_bcov_bsq_from_4vel(pp, utcon, ucov, &geomBL, bcon, bcov, &bsq);
#endif
		  calc_Tij(pp,&geomBL,Tij);
		  indices_2221(Tij,Tij,geomBL.gg);

		  for(i=0;i<4;i++)
		    {
		      Titmagn[i] = bsq*utcon[i]*ucov[0] - bcon[i]*bcov[0];
		      Titkin[i] = rho*utcon[i]*ucov[0];
		      Tituint[i] = uint*utcon[i]*ucov[0];
		    }

#ifdef RADIATION
		  calc_Rij(pp,&geomBL,Rij);

		  #ifdef RADFLUXFFINOUTPUT
		  boost22_lab2ff(Rij,Rij,pp,geomBL.gg,geomBL.GG);
		  #endif
		  
		  indices_2221(Rij,Rij,geomBL.gg);
		  ldouble Rtthat,uconr[4];
		  calc_ff_Rtt(&get_u(p,0,ix,iy,iz),&Rtthat,uconr,&geomBL);
		  Ehat=-Rtthat; 
		  for(iv=1;iv<4;iv++)
		    Ehatucon[iv]=Ehat*utcon[iv];
		  //lab-frame four fource
		  calc_Gi(pp,&geomBL,Gi,0.0,1,0);//ANDREW 1 for lab frame
		  indices_21(Gi,Gi,geomBL.gg);
#endif
		}

	       //rho=rhoGU2CGS(rho);
	       //temp=tempGU2CGS(temp);

	       ldouble r=geomBL.xx;
	       ldouble th=geomBL.yy;
	       ldouble ph=geomBL.zz;
	     
	       //gas prims
	       fprintf(fout1,"%d %d %d ",ix,iy,iz); //(1-3)

	       fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph); //(4-6)

	       fprintf(fout1,"%.5e %.5e ",rho,uint); //(7-8)

	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",utcon[0],utcon[1],utcon[2],utcon[3]); //(9-12)

	       fprintf(fout1,"%.5e ",volume);// (13)

	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",bsq,bcon[1],bcon[2],bcon[3]); //(14) - (17)
    
	       //stress-energy tensor T^mu_nu
	       for(i=0;i<4;i++)
		 for(j=0;j<4;j++)
		   {
		     fprintf(fout1,"%.5e ",Tij[i][j]); //(18-33)
		   }
	       
	       #ifdef RADIATION
	       //radiative energy denisity in the fluid frame
	       fprintf(fout1,"%.5e ",Ehat); //(34)
	       
	       //radiative stress-energy tensor R^mu_nu
	       for(i=0;i<4;i++)
		 for(j=0;j<4;j++)
		   {
		     fprintf(fout1,"%.5e ",Rij[i][j]); //(35-50)
		   }

	       //radiative four-force G_nu
	       for(i=0;i<4;i++)
		   {
		     fprintf(fout1,"%.5e ",Gi[i]); //(51-54)
		   }
	       ldouble Trad;
	       Trad = calc_ncompt_Thatrad(pp, &geomBL, Ehat);
	       fprintf(fout1, "%.5e ", Trad); //55
	       fprintf(fout1, "%.5e ", pp[NF]); //56

	       //advective flux of radiative energy
	       for(i=1;i<4;i++)
		   {
		     fprintf(fout1,"%.5e ",Ehatucon[i]); //(57-59)
		   }

	       //flux of magnetic energy
	       for(i=0;i<4;i++)
		 fprintf(fout1,"%.5e ",Titmagn[i]); //(60-63)

	       //flux of ~kinetic energy
	       for(i=0;i<4;i++)
		 fprintf(fout1,"%.5e ",Titkin[i]);//(64-67)

	       //flux of thermal energy
	       for(i=0;i<4;i++)
		 fprintf(fout1,"%.5e ",Tituint[i]);//(68-71)

	       #endif

	       fprintf(fout1,"\n");
	       
	     }
	 }
     }

   fflush(fout1);
   fclose(fout1);

   return 0;
 }

                              
/*********************************************/
/*********************************************/
/*********************************************/
/* prints quantities in first slice  */
/*********************************************/
/*********************************************/
/*********************************************/
                              
int fprint_slice(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

   fout1=fopen(bufor,"w");
  
   /***********************************/  
   /** writing order is fixed  ********/  
   /***********************************/  
 
   int ix,iy,iz,iv;
   ldouble pp[NV],uu[NV],uucent[NV];
   
   iy=0;
   iz=0;
#ifdef PRINTXGC
   for(ix=-NGCX;ix<NX+NGCX;ix++)
#else
   for(ix=0;ix<NX;ix++)
#endif

     {
       struct geometry geom;
       fill_geometry(ix,iy,iz,&geom);
       fprintf(fout1,"%.5e ",geom.xx); //(1)
       for(iv=0;iv<NV;iv++) {
	 pp[iv]=get_u(p,iv,ix,iy,iz);
	 fprintf(fout1,"%.5e ",pp[iv]); //(2,3,4,5,6,7) for pure hydro
       }
       for(iv=0;iv<NV;iv++) {
	 uu[iv]=get_u(u,iv,ix,iy,iz);
	 fprintf(fout1,"%.5e ",uu[iv]); //(8-13)
       }

#if (PROBLEM==97 || PROBLEM == 98) //HUBBLE
       ldouble rhotrue =  AMBRHO/(1.+VPRIM*t);
       ldouble utrue  = AMBU/pow(1.+VPRIM*t,GAMMA);
       fprintf(fout1,"%.5e ", t);
       fprintf(fout1,"%.5e ", (pp[RHO]-rhotrue)/rhotrue); //error in rho
       fprintf(fout1,"%.5e ", (pp[UU]-utrue)/utrue); //error in uint
       fprintf(fout1,"\n");
#endif
     }

   fflush(fout1);
   fclose(fout1);

   return 0;
 }
