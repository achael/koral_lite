//regrids, keeps Bfield divergence 0 in 2D

//New usage example
// "./regrid old_NX old_NY old_NZ ./old_dumps/old_coord.BL ./old_dumps/old_res.dat new_NX new_NY new_NZ ./new_dumps/new_coord.BL ./new_dumps/new_res.dat NV regrid_x? regrid_y? regrid_z? periodicphi? perterbation interpolate?

// NV setting depending on simmulation
//    NV=6 for pure hydrodynamics
//    NV=9 for pure magnetohydrodynamics
//    NV=8 for HD + EVOLVEELECTRONS
//    NV=11 for MHD + EVOLVEELECTRONS
//    NV=13 for radiation MHD without EVOLVEPHOTONS
//    NV=14 for radiation MHD with EVOLVEPHOTONS
//    NV>13 for RELELECTRONS

//regrids 
#include "ko.h"

int 
main(int argc, char **argv)
{  

  #ifdef MPI
  mpi_myinit(argc,argv);
  if(PROCID==0) {
      printf("regrid works on shared memory only, do not use MPI\n");
      exit(-1);
  }
  #endif
  
  #ifdef OMP
  omp_myinit();  
  #endif

  if (argc!=18 && argc!=12)
    {
      printf("usage: ./regrid NX1 NY1 NZ1 coord1.dat res0001.dat NX2 NY2 NZ2 coord2.dat res0002.dat NV [regridinx=1 regridiny=1 regridinz=1 periodicphi=1 pert=0.0 interpolate=0]\n");
      exit(1);
    }

  //initialize pointers to functions
  init_pointers();

  //initialize constants
  initialize_constants();

  //no gsl error messages
  gsl_set_error_handler_off();
  
  //random number gen. initialization
  srand ( time(NULL) );

  //preparing arrays
  initialize_arrays();

  //sets the grid
  set_grid(&min_dx,&min_dy,&min_dz,&max_dt);
  alloc_loops(1,0.,0.);

  //precalculates metric etc.
  calc_metric();

  int i1,i2,i3,iv,j1,j2;
  int NX1,NY1,NZ1,NX2,NY2,NZ2;//,NV;
  char file_coord1[100],file_coord2[100],file_res1[100],file_res2[100];
  int periodicphi,regridx,regridy,regridz,interpolate;
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

  //NV=atoi(argv[11]);

  if(argc==18)
    {
      regridx=atoi(argv[12]);
      regridy=atoi(argv[13]);
      regridz=atoi(argv[14]);
      
      periodicphi=atoi(argv[15]);
      pert=atof(argv[16]);

      interpolate=atoi(argv[17]);
    }
  else
    {
      regridx=1;
      regridy=1;
      regridz=1;
      periodicphi=1;
      pert=0.0;
      interpolate=0;
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
  if(interpolate) printf("interpolating\n");

  //allocating memory
  double *x1,**y1,*z1,*x2,**y2,*z2;
  double ****datain,****dataout;
  int *projx21,**projy21,*projz21;

  x1=(double*)malloc(NX1*sizeof(double));
  y1=(double**)malloc(NX1*sizeof(double*)); //2D in case of radius dependent y-coords
  for(i1=0;i1<NX1;i1++) 
  {
    y1[i1]=(double*)malloc(NY1*sizeof(double));
  }
  z1=(double*)malloc(NZ1*sizeof(double));

  x2=(double*)malloc(NX2*sizeof(double));
  y2=(double**)malloc(NX2*sizeof(double*)); //2D in case of radius dependent y-coords
  for(i2=0;i2<NX2;i2++) 
  {
    y2[i2]=(double*)malloc(NY2*sizeof(double));
  }
  z2=(double*)malloc(NZ2*sizeof(double));

  projx21=(int*)malloc(NX2*sizeof(int));
  projy21=(int**)malloc(NX2*sizeof(int*));
  for(i2=0;i2<NX2;i2++) 
  {
    projy21[i2]=(int*)malloc(NY2*sizeof(int));
  }
  projz21=(int*)malloc(NZ2*sizeof(int));

  int *int_x1,**int_y1,*int_z1,*int_x2,**int_y2,*int_z2,*regridcell_x,**regridcell_y,*regridcell_z;

  // lower and upper indices for interpolation (stored as arrays)
  int_x1=(int*)malloc(NX2*sizeof(int));
  int_x2=(int*)malloc(NX2*sizeof(int));
  int_y1=(int**)malloc(NX2*sizeof(int*));
  for(i2=0;i2<NX2;i2++)
  {
    int_y1[i2]=(int*)malloc(NY2*sizeof(int));
  }
  int_y2=(int**)malloc(NX2*sizeof(int*));
  for(i2=0;i2<NX2;i2++)
  {
    int_y2[i2]=(int*)malloc(NY2*sizeof(int));
  }
  int_z1=(int*)malloc(NZ2*sizeof(int));
  int_z2=(int*)malloc(NZ2*sizeof(int));
  //defined so we can avoid regridding within the horizon and at the boundaries in general
  regridcell_x=(int*)malloc(NX2*sizeof(int)); 
  regridcell_y=(int**)malloc(NX2*sizeof(int*));
  for(i2=0;i2<NX2;i2++)
  {
    regridcell_y[i2]=(int*)malloc(NY2*sizeof(int));
  }
  regridcell_z=(int*)malloc(NZ2*sizeof(int));

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
      y1[ix][iy]=v2;
      z1[iz]=v3;
    }

  /*
  for(i1=0;i1<NX1;i1++)
    printf("x1: %d %f\n",i1,x1[i1]);
  for(i1=0;i1<NY1;i1++)
    printf("y1: %d %f\n",i1,y1[0][i1]);
  for(i1=0;i1<NZ1;i1++)
    printf("z1: %d %f\n",i1,z1[i1]);
  */

  //output coordinates
  while(fscanf(fcoord2,"%d %d %d %le %le %le\n",&ix,&iy,&iz,&v1,&v2,&v3)!=EOF)
    {
      //serious redundancy here - but not harmful
      x2[ix]=v1;
      y2[ix][iy]=v2;
      z2[iz]=v3;
    }

  /*
  for(i1=0;i1<NX2;i1++)
    printf("x2: %d %f\n",i1,x2[i1]);
  for(i1=0;i1<NY2;i1++)
    printf("y2: %d %f\n",i1,y2[0][i1]);
  for(i1=0;i1<NZ2;i1++)
    printf("z2: %d %f\n",i1,z2[i1]);
  */

  fclose(fcoord1);
  fclose(fcoord2);

  //searching for projection
  double dist,mindist;
  int besti,bestj;

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
      if( i2 == 0 || i2 == (NX2-1) || regridx == 0 )
      {
         regridcell_x[i2] = 0;
      }
      else
      {
        if( x2[i2] <= x1[0] ) // x coord in new grid extends beyond old grid limits, don't interpolate
        {
          regridcell_x[i2] = 0;
          besti = 0;
        }
        else if( x2[i2] >= x1[NX1-1] ) // x coord in new grid extends beyond old grid limits, don't interpolate
        {
          regridcell_x[i2] = 0;
          besti = NX1-1;
        }
        else
        {
          regridcell_x[i2] = 1;
          if( x2[i2] >= x1[besti] )
          { 
            int_x1[i2] = besti;
            int_x2[i2] = besti + 1;
          }
          else if( x2[i2] < x1[besti] )
          { 
            int_x1[i2] = besti - 1;
            int_x2[i2] = besti;
          }
        } 
      }
      projx21[i2]=besti;
      //printf("projz21[%d]=%d\n",i2,besti);
    }

  //in y
  for(i2=0;i2<NX2;i2++)
  {
    //use best radial index (found above) in case of non-uniform polar coords (theta(r))
    besti=projx21[i2];

    for(j2=0;j2<NY2;j2++)
    {
      mindist=1.e50;
      bestj=-1;
      for(j1=0;j1<NY1;j1++)
	{
	  dist=fabs(y2[i2][j2]-y1[besti][j1]);
	  if(dist<mindist || bestj<0)
	    {
	      mindist=dist;
	      bestj=j1;
	    }
	}
      if(NY2 >= 1) // Don't interpolate if in 1D/2D with x-z
      {
        if( j2 == 0 || j2 == (NY2-1) || regridy == 0 )
        {
           regridcell_y[i2][j2] = 0;
        }
        else
        {
          if( y2[i2][j2] <= y1[besti][0] ) // y coord in new grid extends beyond old grid limits, don't interpolate
          {
            regridcell_y[i2][j2] = 0;
            bestj = 0;
          }
          else if( y2[i2][j2] >= y1[besti][NY1-1] ) // y coord in new grid extends beyond old grid limits, don't interpolate
          {
            regridcell_y[i2][j2] = 0;
            bestj = NY1-1;
          }
          else
          {
            regridcell_y[i2][j2] = 1;
            if( y2[i2][j2] >= y1[besti][bestj] )
            { 
              int_y1[i2][j2] = bestj;
              int_y2[i2][j2] = bestj + 1;
            }
            else if( y2[i2][j2] < y1[besti][bestj] )
            { 
              int_y1[i2][j2] = bestj - 1;
              int_y2[i2][j2] = bestj;
            }
          }
        }
      }
      projy21[i2][j2]=bestj;
      //printf("x1=%f, x2=%f, y1=%f, y2=%f, projy21[%d][%d]=%d\n",x1[besti],x2[i2],y1[besti][bestj],y2[i2][j2],i2,j2,bestj);
    }
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
      if(NZ2 >= 1) // Don't interpolate if in 1D/2D with x-y
      {
        if( i2 == 0 || i2 == (NZ2-1) || regridz == 0)
        {
           regridcell_z[i2] = 0;
        }
        else
        {
          if( z2[i2] <= z1[0] ) // z coord in new grid extends beyond old grid limits, don't interpolate
          {
            regridcell_z[i2] = 0;
            besti = 0;
          }
          else if( z2[i2] >= z1[NZ1-1] ) // z coord in new grid extends beyond old grid limits, don't interpolate
          {
            regridcell_z[i2] = 0;
            besti = NZ1-1;
          }
          else
          {
            regridcell_z[i2] = 1;
            if( z2[i2] >= z1[besti] )
            { 
              int_z1[i2] = besti;
              int_z2[i2] = besti + 1;
            } 
            else if( z2[i2] < z1[besti] )
            { 
              int_z1[i2] = besti - 1;
              int_z2[i2] = besti;
            }
          }
        }
      }
      projz21[i2]=besti;
      //printf("projz21[%d]=%d\n",i2,besti);
    }

  //printf("Found nearest neighbors. \n");

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

  //printf("Read data. \n");
  //printf("%e %e \n",datain[0][0][0][0],dataout[0][0][0][0]);

   double dx1,dx2,dy1,dy2,dz1,dz2,DX,DY,DZ;
   double fx1y1z1,fx2y1z1,fx1y2z1,fx2y2z1,fx1y1z2,fx2y1z2,fx1y2z2,fx2y2z2;
   double fxy1z1,fxy2z1,fxy1z2,fxy2z2,fxyz1,fxyz2,fxyz; // 3D interp
   double fxy1,fxy2,fxy; //2D interp
   double fxz1,fxz2,fxz;
   double fyz1,fyz2,fyz;
   double fx,fy,fz;
   //variables for calculation of Eb
   double Bxx,Byy,Bzz,dxx,dyy,dzz,old_Eb,new_Eb,alpha;
   double pratio;
   //need to convert from MYCOORDS to OUTCOORDS as we do in silo? Coord file is coordBL
   double xxvec[4];
   struct geometry geomin; //just used to output r, th values
   struct geometry geomBL; //just used to output r, th values
   struct geometry geomout; //just used to output r, th values

   //initialize Eb to 0
   old_Eb = 0.;
   new_Eb = 0.;

   //constructing the output data
   for(i1=0;i1<NX2;i1++)
     for(i2=0;i2<NY2;i2++)
       for(i3=0;i3<NZ2;i3++)
	 {
           //printf("%i %i %i %i %i %i \n",i1,i2,i3,NX2,NY2,NZ2);
           int tix,tiy,tiz,ix,iy,iz;
	   if(regridx) 
	     tix=projx21[i1];
	   else
	     tix=i1;

	   if(regridy) 
	     tiy=projy21[i1][i2];
	   else
	     tiy=i2;

	   if(regridz) 
	     tiz=projz21[i3];
	   else
	     tiz=i3;

           //variables for interpolation
           ix = i1;
           iy = i2;
           iz = i3;
           //printf("Made it through setting integers.\n");            

           if(regridx)
           {
             if(regridcell_x[ix])
             {
               dx2 = x1[int_x2[ix]] - x2[ix];           
               dx1 = -x1[int_x1[ix]] + x2[ix];          
               DX = x1[int_x2[ix]] - x1[int_x1[ix]];
             }
             else
             {
               dx1=0;
               dx2=0;
               DX=1.;
             }
           }

           //printf("%e %e %e \n",dx1,dx2,DX);            

           if(regridy)
           {
             if(regridcell_y[ix][iy])
             {
               dy2 = y1[tix][int_y2[ix][iy]] - y2[ix][iy];
               dy1 = -y1[tix][int_y1[ix][iy]] + y2[ix][iy];
               DY = y1[tix][int_y2[ix][iy]] - y1[tix][int_y1[ix][iy]];
             }
             else
             {
               dy1=0;
               dy2=0;
               DY=1.;
             }
           }

           //printf("%e %e %e \n",dy1,dy2,DY);            

           //printf("%e %e \n",z1[0],z2[0]);            

           if(regridz)
           {
             if(regridcell_z[iz])
             {
               dz2 = z1[int_z2[iz]] - z2[iz];
               dz1 = -z1[int_z1[iz]] + z2[iz];
               DZ = z1[int_z2[iz]] - z1[int_z1[iz]];
             }
             else
             {
               dz1=0;
               dz2=0;
               DZ=1.;
             }
           }

           //compute magnetic energy using old field values
           #ifdef MAGNFIELD
           dxx = get_size_x(ix, 0);
           dyy = get_size_x(iy, 1);
           dzz = get_size_x(iz, 2);
           if(NZ2 == 1)
           {
             dzz *= 2.*M_PI;
           }
           #ifdef PHIWEDGE
           if(NZ2 > 1) dzz*=(2. * M_PI / PHIWEDGE);
           #endif
           Bxx = datain[projx21[ix]][projy21[ix][iy]][projz21[iz]][B1];
           Byy = datain[projx21[ix]][projy21[ix][iy]][projz21[iz]][B2];
           Bzz = datain[projx21[ix]][projy21[ix][iy]][projz21[iz]][B3];
           old_Eb += (Bxx*Bxx + Byy*Byy + Bzz*Bzz)*dxx*dyy*dzz*pick_gdet(ix,iy,iz);
           #endif

           //printf("%i %i %i \n",regridcell_x[0],regridcell_y[0],regridcell_z[0]);

           //printf("Made it through defining steps.\n");

           //printf("%i %i %i ",ix,iy,iz);

	   for(iv=0;iv<NV;iv++)    
           {
            if(interpolate)
            {
             if(regridcell_x[ix] && regridcell_y[ix][iy] && regridcell_z[iz])
             {
               fx1y1z1 = datain[int_x1[ix]][int_y1[ix][iy]][int_z1[iz]][iv];
               fx2y1z1 = datain[int_x2[ix]][int_y1[ix][iy]][int_z1[iz]][iv];
               fx1y2z1 = datain[int_x1[ix]][int_y2[ix][iy]][int_z1[iz]][iv];
               fx2y2z1 = datain[int_x2[ix]][int_y2[ix][iy]][int_z1[iz]][iv];
               fx1y1z2 = datain[int_x1[ix]][int_y1[ix][iy]][int_z2[iz]][iv];
               fx2y1z2 = datain[int_x2[ix]][int_y1[ix][iy]][int_z2[iz]][iv];
               fx1y2z2 = datain[int_x1[ix]][int_y2[ix][iy]][int_z2[iz]][iv];
               fx2y2z2 = datain[int_x2[ix]][int_y2[ix][iy]][int_z2[iz]][iv];

               fxy1z1 = (dx2/DX)*fx1y1z1 + (dx1/DX)*fx2y1z1;
               fxy2z1 = (dx2/DX)*fx1y2z1 + (dx1/DX)*fx2y2z1;
               fxy1z2 = (dx2/DX)*fx1y1z2 + (dx1/DX)*fx2y1z2;
               fxy2z2 = (dx2/DX)*fx1y2z2 + (dx1/DX)*fx2y2z2;

               fxyz1 = (dy2/DY)*fxy1z1 + (dy1/DY)*fxy2z1;
               fxyz2 = (dy2/DY)*fxy1z2 + (dy1/DY)*fxy2z2;

               fxyz = (dz2/DZ)*fxyz1 + (dz1/DZ)*fxyz2;
               dataout[i1][i2][i3][iv] = fxyz;
             }
             else if(regridcell_x[ix] && regridcell_y[ix][iy] && regridcell_z[iz] == 0)
             {
               fx1y1z1 = datain[int_x1[ix]][int_y1[ix][iy]][tiz][iv]; 
               fx2y1z1 = datain[int_x2[ix]][int_y1[ix][iy]][tiz][iv]; 
               fx1y2z1 = datain[int_x1[ix]][int_y2[ix][iy]][tiz][iv]; 
               fx2y2z1 = datain[int_x2[ix]][int_y2[ix][iy]][tiz][iv];

               fxy1 = (dx2/DX)*fx1y1z1 + (dx1/DX)*fx2y1z1;
               fxy2 = (dx2/DX)*fx1y2z1 + (dx1/DX)*fx2y2z1;

               fxy = (dy2/DY)*fxy1 + (dy1/DY)*fxy2;
               dataout[i1][i2][i3][iv] = fxy;
             }
             else if(regridcell_x[ix] && regridcell_y[ix][iy] == 0 && regridcell_z[iz])
             {
               fx1y1z1 = datain[int_x1[ix]][tiy][int_z1[iz]][iv]; 
               fx2y1z1 = datain[int_x2[ix]][tiy][int_z1[iz]][iv]; 
               fx1y1z2 = datain[int_x1[ix]][tiy][int_z2[iz]][iv]; 
               fx2y1z2 = datain[int_x2[ix]][tiy][int_z2[iz]][iv]; 

               fxz1 = (dx2/DX)*fx1y1z1 + (dx1/DX)*fx2y1z1;
               fxz2 = (dx2/DX)*fx1y1z2 + (dx1/DX)*fx2y1z2;

               fxz = (dz2/DZ)*fxz1 + (dz1/DZ)*fxz2;
               dataout[i1][i2][i3][iv] = fxz;
             }
             else if(regridcell_x[ix] == 0 && regridcell_y[ix][iy] && regridcell_z[iz])
             {
               fx1y1z1 = datain[tix][int_y1[ix][iy]][int_z1[iz]][iv]; 
               fx1y2z1 = datain[tix][int_y2[ix][iy]][int_z1[iz]][iv]; 
               fx1y1z2 = datain[tix][int_y1[ix][iy]][int_z2[iz]][iv]; 
               fx1y2z2 = datain[tix][int_y2[ix][iy]][int_z2[iz]][iv]; 

               fyz1 = (dy2/DY)*fx1y1z1 + (dy1/DY)*fx1y2z1;
               fyz2 = (dy2/DY)*fx1y1z2 + (dy1/DY)*fx1y2z2;

               fyz = (dz2/DZ)*fyz1 + (dz1/DZ)*fyz2;
               dataout[i1][i2][i3][iv] = fyz;
             }
             else if(regridcell_x[ix] && regridcell_y[ix][iy] == 0 && regridcell_z[iz] == 0)
             {
               fx1y1z1 = datain[int_x1[ix]][tiy][tiz][iv];
               fx2y1z1 = datain[int_x2[ix]][tiy][tiz][iv];

               fx = (dx2/DX)*fx1y1z1 + (dx1/DX)*fx2y1z1;
               dataout[i1][i2][i3][iv] = fx;
             }
             else if(regridcell_x[ix] == 0 && regridcell_y[ix][iy] && regridcell_z[iz] == 0)
             {
               fx1y1z1 = datain[tix][int_y1[ix][iy]][tiz][iv]; 
               fx1y2z1 = datain[tix][int_y2[ix][iy]][tiz][iv]; 

               fy = (dy2/DY)*fx1y1z1 + (dy1/DY)*fx1y2z1;
               dataout[i1][i2][i3][iv] = fy;
             }
             else if(regridcell_x[ix] == 0 && regridcell_y[ix][iy] == 0 && regridcell_z[iz])
             {
               fx1y1z1 = datain[tix][tiy][int_z1[iz]][iv];
               fx1y1z2 = datain[tix][tiy][int_z2[iz]][iv]; 

               fz = (dz2/DZ)*fx1y1z1 + (dz1/DZ)*fx1y1z2;
               dataout[i1][i2][i3][iv] = fz;
             }
             else if(regridcell_x[ix] == 0 && regridcell_y[ix][iy] == 0 && regridcell_z[iz] == 0)
             {
               dataout[i1][i2][i3][iv]=datain[tix][tiy][tiz][iv];
             }
            }
            else
            {
               dataout[i1][i2][i3][iv]=datain[tix][tiy][tiz][iv];
            }

            //convert MKS2 -> JETCOORDS
	    //fill_geometry_arb(i1,i2,i3,&geomin,MKS2COORDS);

            //get_xx(i1,i2,i3,xxvec);
            //xxvec[0]=0.;xxvec[1]=get_x(ix,0);xxvec[2]=get_x(iy,1);xxvec[3]=get_x(iz,2);
            //#ifdef PRECOMPUTE_MY2OUT
            //trans_pall_coco_my2out(dataout[i1][i2][i3],dataout[i1][i2][i3],&geomin,&geomout);
            //#else      
            //trans_pall_coco(dataout[i1][i2][i3], dataout[i1][i2][i3], MKS2COORDS,JETCOORDS, xxvec,&geomin,&geomout);
            //#endif

            //if(x2[ix] <= 2.) // damp B inside horizon?
            //{
            //  if(iv == B1 || iv == B2 || iv == B3) 
            //  { 
            //    printf("%i %i %i old pp[%i] = %e ",ix,iy,iz,iv,dataout[ix][iy][iz][iv]);
            //    dataout[ix][iy][iz][iv] *= 0.1;
            //    printf(", new pp[%i] = %e \n",iv,dataout[ix][iy][iz][iv]);
            //  }
            //}

            //pratio = fabs(dataout[i1][i2][i3][iv])/fabs(datain[tix][tiy][tiz][iv]);
            //printf(" %f ",pratio);
            //if(iv==(NV-1)) printf("\n");
            //if(pratio > 10.) //|| pratio < 0.1)
            //{
            //  printf("ix iy iz iv pratio : %i %i %i %i %f \n",ix,iy,iz,iv,pratio);
            //  printf("dx1 dx2 DX dy1 dy2 DY : %f %f %f %f %f %f \n",dx1,dx2,DX,dy1,dy2,DY);
            //  printf("px1y1 px2y1 px1y2 px2y2 : %e %e %e %e \n \n",datain[int_x1[ix]][int_y1[ix][iy]][tiz][iv],datain[int_x2[ix]][int_y1[ix][iy]][tiz][iv],datain[int_x1[ix]][int_y2[ix][iy]][tiz][iv],datain[int_x2[ix]][int_y2[ix][iy]][tiz][iv]);
            //}
           }

           if(i1==0)// && i2==0)
           { 
             printf("Filling xxvec...\n");
             get_xx(ix,iy,iz,xxvec);
             printf("%f %f %f %f\n",xxvec[0],xxvec[1],xxvec[2],xxvec[3]);

             printf("Filling geomout\n");
	     fill_geometry(ix,iy,iz,&geomout);
             printf("Filling geomin\n");
	     fill_geometry_arb(ix,iy,iz,&geomin,BLCOORDS);//MKS2COORDS);
             printf("%f %f \n",geomin.gdet,geomout.gdet);
             printf("%f %f %f \n",geomin.xx,geomin.yy,geomin.zz);
           }

	   if(pert>0.0) dataout[i1][i2][i3][VZ]*=(1.+((double)rand()/(double)RAND_MAX-0.5)*2.*pert);
	 }

   //printf("Interpolated/copied data. \n");

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

  double Br11,Br12,Br21,Br22,Bth11,Bth12,Bth21,Bth22;
  double divB,a1,a2,a3; //prefactors for solving for new B-field

  //printf("%e %e %e \n",pick_gdet(0,0,0),get_x(0,0),get_x(0,1));

  //perform divergence cleaning (only in 2D) if simulation has B-field
  // (this is non-relativistic, may have to expand code to handle this properly
  #ifdef MAGNFIELD //Only use for MHD
  for(ix=0;ix<NX2;ix++)
    for(iy=0;iy<NY2;iy++)
      for(iz=0;iz<NZ2;iz++)
	{
            if( NX2 > 1 && NY2 > 1 && NZ2 == 1 )
            {
              if( ix > 0 && iy > 0 )
              {
                Br11 = dataout[ix-1][iy-1][iz][B1];
                Br12 = dataout[ix-1][iy][iz][B1];
                Br21 = dataout[ix][iy-1][iz][B1];
                Br22 = dataout[ix][iy][iz][B1];

                Bth11 = dataout[ix-1][iy-1][iz][B2];
                Bth12 = dataout[ix-1][iy][iz][B2];
                Bth21 = dataout[ix][iy-1][iz][B2];
                Bth22 = dataout[ix][iy][iz][B2];

                //divergence in spherical coords
                //divB = (0.5/dx)*(Bx22 + Bx21 - Bx12 - Bx11) + (0.5/dy)*(Bx22 + Bx21 - Bx12 - Bx11);
                divB = (pick_gdet(ix,iy,iz)*Br22 + pick_gdet(ix,iy-1,iz)*Br21 
	      - pick_gdet(ix-1,iy,iz)*Br12 - pick_gdet(ix-1,iy-1,iz)*Br11)/(2.*(get_x(ix+1,0)-get_x(ix,0)))
	+ (pick_gdet(ix,iy,iz)*Bth22 + pick_gdet(ix-1,iy,iz)*Bth12
	   - pick_gdet(ix,iy-1,iz)*Bth21 - pick_gdet(ix-1,iy-1,iz)*Bth11)/(2.*(get_x(iy+1,1)-get_x(iy,1)));
                //printf("ix iy iz Br Bth divB : %i %i %i %e %e %e %e %e \n",ix,iy,iz,Br22,Bth22,Bx22,By22,divB);
                //printf("ix iy iz g B1 B2 x y divB: %i %i %i %e %e %e %e %e %e \n",ix,iy,iz,pick_gdet(ix,iy,iz),Br22,Bth22,get_x(ix+1,0),get_x(iy+1,1),fabs(divB));

                // Set divergence to zero if within domain
                if(ix > 0 && iy > 0)
                {
                    if(fabs(divB) > 1.e-20)
                    {
                      a1 = (2.*(get_x(ix+1,0)-get_x(ix,0)))*(pick_gdet(ix,iy,iz)*Bth22 + pick_gdet(ix-1,iy,iz)*Bth12
	   - pick_gdet(ix,iy-1,iz)*Bth21 - pick_gdet(ix-1,iy-1,iz)*Bth11)/(2.*(get_x(iy+1,1)-get_x(iy,1)));
                      a2 = (pick_gdet(ix,iy-1,iz)*Br21 
	      - pick_gdet(ix-1,iy,iz)*Br12 - pick_gdet(ix-1,iy-1,iz)*Br11);
                      a3 = pick_gdet(ix,iy,iz);
                      Br22 = -(a1+a2)/a3;
                      //Calc new divB
                      divB = (pick_gdet(ix,iy,iz)*Br22 + pick_gdet(ix,iy-1,iz)*Br21 
	      - pick_gdet(ix-1,iy,iz)*Br12 - pick_gdet(ix-1,iy-1,iz)*Br11)/(2.*(get_x(ix+1,0)-get_x(ix,0)))
	+ (pick_gdet(ix,iy,iz)*Bth22 + pick_gdet(ix-1,iy,iz)*Bth12
	   - pick_gdet(ix,iy-1,iz)*Bth21 - pick_gdet(ix-1,iy-1,iz)*Bth11)/(2.*(get_x(iy+1,1)-get_x(iy,1)));

                      //fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
                      //printf("ix iy iz x y r th g B1 B2 divB: %i %i %i %e %e %e %e %e %e %e %e \n\n",ix,iy,iz,get_x(ix+1,0),get_x(iy+1,1),geomBL.xx,geomBL.yy,pick_gdet(ix,iy,iz),Br22,Bth22,divB);

                      dataout[ix][iy][iz][B1] = Br22;
                    }
                }
              }
            }
	}
   #endif

   //calc new Eb
   #ifdef MAGNFIELD
   for(ix=0;ix<NX2;ix++)
    for(iy=0;iy<NY2;iy++)
      for(iz=0;iz<NZ2;iz++)
	{
          //compute magnetic energy using new field values
          dxx = get_size_x(ix, 0);
          dyy = get_size_x(iy, 1);
          dzz = get_size_x(iz, 2);
          if(NZ2 == 1)
          {
            dzz *= 2.*M_PI;
          }
          #ifdef PHIWEDGE
          if(NZ2 > 1) dzz*=(2. * M_PI / PHIWEDGE);
          #endif
          Bxx = dataout[ix][iy][iz][B1];
          Byy = dataout[ix][iy][iz][B2];
          Bzz = dataout[ix][iy][iz][B3];
          new_Eb += (Bxx*Bxx + Byy*Byy + Bzz*Bzz)*dxx*dyy*dzz*pick_gdet(ix,iy,iz);
	}
   #endif

  //compute factor for rescaling B-field to maintain same total magnetic energy in domain
  alpha = 1.;
  #ifdef MAGNFIELD
  alpha = sqrt(old_Eb/new_Eb);
  #endif

  double newnew_Eb;
  newnew_Eb = 0.;

  for(ix=0;ix<NX2;ix++)
    for(iy=0;iy<NY2;iy++)
      for(iz=0;iz<NZ2;iz++)
	{ 
          #ifdef MAGNFIELD
          //rescale magnetic field
          dataout[ix][iy][iz][B1] *= alpha;
          dataout[ix][iy][iz][B2] *= alpha;
          dataout[ix][iy][iz][B3] *= alpha;

          //compute magnetic energy using new field values
          dxx = get_size_x(ix, 0);
          dyy = get_size_x(iy, 1);
          dzz = get_size_x(iz, 2);
          if(NZ2 == 1)
          {
            dzz *= 2.*M_PI;
          }
          #ifdef PHIWEDGE
          if(NZ2 > 1) dzz*=(2. * M_PI / PHIWEDGE);
          #endif
          Bxx = dataout[ix][iy][iz][B1];
          Byy = dataout[ix][iy][iz][B2];
          Bzz = dataout[ix][iy][iz][B3];
          newnew_Eb += (Bxx*Bxx + Byy*Byy + Bzz*Bzz)*dxx*dyy*dzz*pick_gdet(ix,iy,iz);
          #endif

	  fwrite(dataout[ix][iy][iz],sizeof(double),NV,fout);
	}

  #ifdef MAGNFIELD
  //printf("old_Eb new_Eb newnew_Eb a a^2 : %e %e %e %e %e \n",old_Eb,new_Eb,newnew_Eb,alpha,alpha*alpha);
  #endif

  fclose(fout);

  printf("done! remember to modify res????.head manually if resolution has changed.\n");

  return 0;
}
