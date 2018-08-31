#include <stdio.h>
#include <stdlib.h>

#include "ko.h"

int
main
(int argc, char **argv)
{
  #ifdef MPI
  printf("phiavg works on shared memory only, do not use MPI, please\n");
  exit(-1);
  #endif
  
  //which files to read
  int no1,no2,nostep,ifavg,ifoutavg;
  
  int ifile,itot=0;
  int i,j,k,iv,ix,iy,iz,ic,gix,giy,giz;
  int nx,ny,nz,nv,ret,problem;
  double time;

  if(argc!=6 && argc!=5 && argc!=4)
    {
      printf("Not enough input arguments. Asks for ./phiavg no1 no2 nostep [ifavg=0 ifoutavg=0]\n");
      return -1;
    }
  if(argc==6)
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);
      ifavg=atof(argv[4]);
      ifoutavg=atof(argv[5]);
    }
  if(argc==5)
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);
      ifavg=atof(argv[4]);
      ifoutavg=0;
    }
  if(argc==4)
    {      
      no1=atof(argv[1]);
      no2=atof(argv[2]);
      nostep=atof(argv[3]);
      ifavg=0;
      ifoutavg=0;
    }

  printf("phiavg: %d %d %d %d %d\n",no1,no2,nostep,ifavg,ifoutavg);
  //if(ifavg) printf("for avg set the proper TNZ!\n");

  char folder[100],folderout[100],bufer[100],base[10],fname[100],fnamehead[100],fnameout[100],fnameheadout[100],bufor[500];
  sprintf(folder,"%s","dumps");
  sprintf(folderout,"%s","dumps_phiavg");

  nv=NV;

  if (ifoutavg)
  {
    if(ifavg)
      sprintf(base,"avgavg%04d-",no1);   
    else
      sprintf(base,"avgres%04d-",no1);
    no1=no2;
    nv+=NAVGVARS;    
  }

  else
  {
    if(ifavg)
    {
      sprintf(base,"%s","avg");
      nv+=NAVGVARS;
    }
    else
      sprintf(base,"%s","res");
  }

  //allocate memory
  nx=NX;ny=NY;
  double ***prims,*pp;
  pp=(double*)malloc(nv*sizeof(double));
  prims=(double***)malloc(nx*sizeof(double**));
  int **ngood=(int**)malloc(nx*sizeof(int*));
  for(i=0;i<nx;i++)
    {
      ngood[i]=(int*)malloc(ny*sizeof(int));
      prims[i]=(double**)malloc(ny*sizeof(double*));
      for(j=0;j<ny;j++)
	prims[i][j]=(double*)malloc(nv*sizeof(double));
    }
  int **indices;
      
  for(ifile=no1;ifile<=no2;ifile+=nostep)
    {
      itot++;

      //reading the files
      sprintf(fname,"%s/%s%04d.dat",folder,base,ifile);
      sprintf(fnamehead,"%s/%s%04d.head",folder,base,ifile);
      sprintf(fnameheadout,"%s/%s%04d.head",folderout,base,ifile);
      sprintf(fnameout,"%s/%s%04d.dat",folderout,base,ifile);

      FILE *fdump,*fout;

      /***********/
      //header file
      fdump=fopen(fnamehead,"r");

      //reading damp file parameters
      int intpar[6];
      
   
      if(ifavg || ifoutavg) //avg follows old writing sequence - to correct!!!
	{
	  sprintf(bufor,"cp %s %s\n",fnamehead,fnameheadout);
	  system(bufor);
	  nx=TNX;
	  ny=TNY;
	  nz=TNZ;
	  problem=PROBLEM;

	  fout=fopen(fnameout,"w");
	  
	  //body file
	  fdump=fopen(fname,"rb");
 
	  printf("avg file (%s) read no. %d\n",fname,itot);
	}
      else
	{
	  ret=fscanf(fdump,"## %d %d %lf %d %d %d %d\n",&intpar[0],&intpar[1],&time,&intpar[2],&intpar[3],&intpar[4],&intpar[5]);
	  problem=intpar[2];
	  nx=intpar[3];
	  ny=intpar[4];
	  nz=intpar[5];
	  fclose(fdump);
      
	  fout=fopen(fnameheadout,"w");
	  sprintf(bufor,"## %5d %5d %10.6e %5d %5d %5d %5d\n",intpar[0],intpar[1],time,intpar[2],intpar[3],intpar[4],intpar[5]);
	  fprintf(fout,"%s",bufor);
	  fprintf(fout,"phi-averaged to NZ=1\n");
	  fclose(fout);
       
	  /***********/
	  fout=fopen(fnameout,"w");
      
	  if(nx!=TNX || ny!=TNY || problem!=PROBLEM)
	    {
	      printf("PROBLEM files inconsistent with the dump file\n");
	      exit (-1);
	    }

	  printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n",
		 fname,intpar[0],time,intpar[2],intpar[3],intpar[4],intpar[5]); 


	  printf("time: %e resolution: %d x %d x %d NV: %d\n",time,nx,ny,nz,nv);

	  //body file
	  fdump=fopen(fname,"rb");
	}
 
      for(j=0;j<ny;j++)
	for(i=0;i<nx;i++)
	  {
	    ngood[i][j]=0;
	    for(iv=0;iv<nv;iv++)
	      prims[i][j][iv]=0.;
	  }

      if(itot==1)
	{
      indices = (int **)malloc(nx*ny*nz*sizeof(int*));
      for(i=0;i<nx*ny*nz;i++)
	indices[i]=(int *)malloc(3*sizeof(int));
	}

      //first indices
      for(ic=0;ic<nx*ny*nz;ic++)
	{
	  ret=fread(&ix,sizeof(int),1,fdump);
	  ret=fread(&iy,sizeof(int),1,fdump);
	  ret=fread(&iz,sizeof(int),1,fdump);

	  //phiavg only serial
	  //mpi_global2localidx(gix,giy,giz,&ix,&iy,&iz);

	  indices[ic][0]=ix;
	  indices[ic][1]=iy;
	  indices[ic][2]=iz;
	}

      //then primitives
      for(ic=0;ic<nx*ny*nz;ic++)
	{
	  ret=fread(pp,sizeof(double),nv,fdump);

	  gix=indices[ic][0];
	  giy=indices[ic][1];
	  giz=indices[ic][2];
	  
	  if(gix<0 || gix>=nx ||giy<0 || giy>=ny)
	    printf("blont: %d %d %d | %d %d %d\n",gix,giy,giz,nx,ny,nz);
	  else
	    {
	      ngood[gix][giy]++;
	      for(iv=0;iv<nv;iv++)
		prims[gix][giy][iv]+=pp[iv];
	    }
	}

      for(j=0;j<ny;j++)
	for(i=0;i<nx;i++)
	  for(iv=0;iv<nv;iv++)
	    prims[i][j][iv]/=(double)ngood[i][j];

      //print to a file
      
      //indices first
      iz=0;
      for(ix=0;ix<nx;ix++)
	for(iy=0;iy<ny;iy++)
	  {
	    fwrite(&ix,sizeof(int),1,fout);
	    fwrite(&iy,sizeof(int),1,fout);
	    fwrite(&iz,sizeof(int),1,fout);
	  }

      //then, in the same order, primitives
      for(ix=0;ix<nx;ix++)
	for(iy=0;iy<ny;iy++)
	  {
	    fwrite(&prims[ix][iy][0],sizeof(ldouble),nv,fout);
	  }

    

      fclose(fdump);
      fclose(fout);
    }

  for(i=0;i<nx*ny*nz;i++)
	free(indices[i]);
      free(indices);
  return 0;
}
