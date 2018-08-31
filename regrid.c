//regrids 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


int 
main(int argc, char **argv)
{  
  int i1,i2,i3,iv;

  if (argc!=17 && argc!=12)
    {
      printf("usage: ./regrid NX1 NY1 NZ1 coord1.dat res0001.dat NX2 NY2 NZ2 coord2.dat res0002.dat NV [regridinx=1 regridiny=1 regridinz=1 periodicphi=1 pert=0.0]\n");
      exit(1);
    }

  int NX1,NY1,NZ1,NX2,NY2,NZ2,NV;
  char file_coord1[100],file_coord2[100],file_res1[100],file_res2[100];
  int periodicphi,regridx,regridy,regridz;
  double pert=0.0;

  //random number gen. initialization
  srand ( time(NULL) );

  NX1=atoi(argv[1]);
  NY1=atoi(argv[2]);
  NZ1=atoi(argv[3]);

  sprintf(file_coord1,"%s",argv[4]);
  sprintf(file_res1,"%s",argv[5]);

  NX2=atoi(argv[6]);
  NY2=atoi(argv[7]);
  NZ2=atoi(argv[8]);

  sprintf(file_coord2,"%s",argv[9]);
  sprintf(file_res2,"%s",argv[10]);

  NV=atoi(argv[11]);

  if(argc==17)
    {
      regridx=atoi(argv[12]);
      regridy=atoi(argv[13]);
      regridz=atoi(argv[14]);
      
      periodicphi=atoi(argv[15]);
      pert=atof(argv[16]);
    }
  else
    {
      regridx=1;
      regridy=1;
      regridz=1;
      periodicphi=1;
      pert=0.0;
    }

  printf("regrid - projects snap shot file onto a new grid\n");
  printf("note: works only on rectangular-like (not MKS3) grids so far\n");
  printf("input: %d x %d x %d | %s | %s\n",NX1,NY1,NZ1,file_coord1,file_res1);
  printf("output: %d x %d x %d | %s | %s\n",NX2,NY2,NZ2,file_coord2,file_res2);
  if(regridx) printf("regriding in x\n");
  if(regridy) printf("regriding in y\n");
  if(regridz) printf("regriding in z\n");
  if(periodicphi) printf("periodic in phi\n");
  if(pert!=0.0) printf("pert=%f\n",pert);

  //allocating memory
  double *x1,*y1,*z1,*x2,*y2,*z2;
  double ****datain,****dataout;
  int *projx21,*projy21,*projz21;

  x1=(double*)malloc(NX1*sizeof(double));
  y1=(double*)malloc(NY1*sizeof(double));
  z1=(double*)malloc(NZ1*sizeof(double));
  x2=(double*)malloc(NX2*sizeof(double));
  y2=(double*)malloc(NY2*sizeof(double));
  z2=(double*)malloc(NZ2*sizeof(double));
  projx21=(int*)malloc(NX2*sizeof(int));
  projy21=(int*)malloc(NY2*sizeof(int));
  projz21=(int*)malloc(NZ2*sizeof(int));
  datain=(double****)malloc(NX1*sizeof(double***));
  for(i1=0;i1<NX1;i1++)
    {
      datain[i1]=(double***)malloc(NY1*sizeof(double**));
      for(i2=0;i2<NY1;i2++)
	{
	  datain[i1][i2]=(double**)malloc(NZ1*sizeof(double*));
	  for(i3=0;i3<NZ1;i3++)
	    {    
	      datain[i1][i2][i3]=(double*)malloc(NV*sizeof(double));
	   }
	}
    }
  dataout=(double****)malloc(NX2*sizeof(double***));
  for(i1=0;i1<NX2;i1++)
    {
      dataout[i1]=(double***)malloc(NY2*sizeof(double**));
      for(i2=0;i2<NY2;i2++)
	{
	  dataout[i1][i2]=(double**)malloc(NZ2*sizeof(double*));
	  for(i3=0;i3<NZ2;i3++)
	    {    
	      dataout[i1][i2][i3]=(double*)malloc(NV*sizeof(double));
	   }
	}
    }
   
  
  //reading in coordinates
  FILE *fcoord1,*fcoord2;

  fcoord1=fopen(file_coord1,"r");
  if(fcoord1==NULL)
    {
      printf("missing %s file.\n",file_coord1); exit(-1);
    }

  fcoord2 = fopen(file_coord2,"r");
  if(fcoord2==NULL)
    {
      printf("missing %s file.\n",file_coord2); exit(-1);
    }

  int ix,iy,iz;
  double v1,v2,v3;

  //input coordinates
  while(fscanf(fcoord1,"%d %d %d %le %le %le\n",&ix,&iy,&iz,&v1,&v2,&v3)!=EOF)
    {
      //serious redundancy here - but not harmful
      x1[ix]=v1;
      y1[iy]=v2;
      z1[iz]=v3;
    }

  /*
  for(i1=0;i1<NX1;i1++)
    printf("x1: %d %f\n",i1,x1[i1]);
  for(i1=0;i1<NY1;i1++)
    printf("y1: %d %f\n",i1,y1[i1]);
  for(i1=0;i1<NZ1;i1++)
    printf("z1: %d %f\n",i1,z1[i1]);
  */

  //output coordinates
  while(fscanf(fcoord2,"%d %d %d %le %le %le\n",&ix,&iy,&iz,&v1,&v2,&v3)!=EOF)
    {
      //serious redundancy here - but not harmful
      x2[ix]=v1;
      y2[iy]=v2;
      z2[iz]=v3;
    }

  /*
  for(i1=0;i1<NX2;i1++)
    printf("x2: %d %f\n",i1,x2[i1]);
  for(i1=0;i1<NY2;i1++)
    printf("y2: %d %f\n",i1,y2[i1]);
  for(i1=0;i1<NZ2;i1++)
    printf("z2: %d %f\n",i1,z2[i1]);
  */

  fclose(fcoord1);
  fclose(fcoord2);

  //searching for projection
  double dist,mindist;
  int besti;

  //in x
  for(i2=0;i2<NX2;i2++)
    {
      mindist=1.e50;
      besti=-1;
      for(i1=0;i1<NX1;i1++)
	{
	  dist=fabs(x2[i2]-x1[i1]);
	  if(dist<mindist || besti<0)
	    {
	      mindist=dist;
	      besti=i1;
	    }
	}
      projx21[i2]=besti;
      //printf("projz21[%d]=%d\n",i2,besti);
    }

  //in y
  for(i2=0;i2<NY2;i2++)
    {
      mindist=1.e50;
      besti=-1;
      for(i1=0;i1<NY1;i1++)
	{
	  dist=fabs(y2[i2]-y1[i1]);
	  if(dist<mindist || besti<0)
	    {
	      mindist=dist;
	      besti=i1;
	    }
	}
      projy21[i2]=besti;
      printf("projy21[%d]=%d\n",i2,besti);
    }

  //in z   
  for(i2=0;i2<NZ2;i2++)
    {
      mindist=1.e50;
      besti=-1;
      for(i1=0;i1<NZ1;i1++)
	{
	  dist=fabs(z2[i2]-z1[i1]);
	  if(periodicphi && (2.*M_PI-dist)<dist)
	    dist=2.*M_PI-dist;
	  if(dist<mindist || besti<0)
	    {
	      mindist=dist;
	      besti=i1;
	    }
	}
      projz21[i2]=besti;
      //printf("projz21[%d]=%d\n",i2,besti);
    }

  //reading in the data file
  FILE *fin = fopen(file_res1,"rb");
  if(fin==NULL) 
    {
      printf("missing %s file.\n",file_res1); exit(-1);
    }
  int **indices,ic;
  double *pp;
  pp=(double*)malloc(NV*sizeof(double));
  indices = (int **)malloc(NX1*NY1*NZ1*sizeof(int*));
  for(ic=0;ic<NX1*NY1*NZ1;ic++)
    indices[ic]=(int *)malloc(3*sizeof(int));

  //first indices
  for(ic=0;ic<NX1*NY1*NZ1;ic++)
    {
      fread(&ix,sizeof(int),1,fin);
      fread(&iy,sizeof(int),1,fin);
      fread(&iz,sizeof(int),1,fin);
      indices[ic][0]=ix;
      indices[ic][1]=iy;
      indices[ic][2]=iz;
    }

  //then primitives
   for(ic=0;ic<NX1*NY1*NZ1;ic++)
    {
      fread(pp,sizeof(double),NV,fin);
      ix=indices[ic][0];
      iy=indices[ic][1];
      iz=indices[ic][2];
      for(iv=0;iv<NV;iv++)    
	{
	  datain[ix][iy][iz][iv]=pp[iv];
	}
    }
   fclose(fin);

   //constructing the output data
   for(i1=0;i1<NX2;i1++)
     for(i2=0;i2<NY2;i2++)
       for(i3=0;i3<NZ2;i3++)
	 {
	   int tix,tiy,tiz;
	   if(regridx) 
	     tix=projx21[i1];
	   else
	     tix=i1;

	   if(regridy) 
	     tiy=projy21[i2];
	   else
	     tiy=i2;

	   if(regridz) 
	     tiz=projz21[i3];
	   else
	     tiz=i3;


	   for(iv=0;iv<NV;iv++)    
	     dataout[i1][i2][i3][iv]=datain[tix][tiy][tiz][iv];
	   if(pert>0.0) dataout[i1][i2][i3][5]*=(1.+((double)rand()/(double)RAND_MAX-0.5)*2.*pert);
	 }

   //writing to the new file
   FILE *fout = fopen(file_res2,"wb");
  
   //indices first
   for(ix=0;ix<NX2;ix++)
    for(iy=0;iy<NY2;iy++)
      for(iz=0;iz<NZ2;iz++)
	{
	  fwrite(&ix,sizeof(int),1,fout);
	  fwrite(&iy,sizeof(int),1,fout);
	  fwrite(&iz,sizeof(int),1,fout);
	}

  //then, in the same order, primitives
  for(ix=0;ix<NX2;ix++)
    for(iy=0;iy<NY2;iy++)
      for(iz=0;iz<NZ2;iz++)
	{
	  fwrite(dataout[ix][iy][iz],sizeof(double),NV,fout);
	}

  fclose(fout);

  printf("done! remember to modify res????.head manually if resolution has changed.\n");

  return 0;
}
