/*! \file fileop.c
 \brief File operations
*/

#include "ko.h"
#include <string.h>
/*************************************************/
/*  adds up current quantities to the pavg array */
/*************************************************/

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

#ifndef MPI
  sprintf(bufor,"%s/scalars.dat",folder);
  fout_scalars=fopen(bufor,"a");
  printf("fopen %s in function openfiles\n", bufor);

  sprintf(bufor,"%s/failures.dat",folder);
  fout_fail=fopen(bufor,"a");
  printf("fopen %s in function openfiles\n", bufor);

#endif

  return 0;
}


/*********************************************/
/* closes file handles */
/*********************************************/

int
fprint_closefiles()
{

#ifndef MPI
  fclose(fout_scalars);
  fclose(fout_fail);  
#endif

  return 0;
}


/*********************************************/
// coordinate grid file *********************//
/*********************************************/

int
fprint_gridfile(char* folder)
{
  FILE* out;
  char bufor[50];
  sprintf(bufor,"%s/grid.dat",folder);
  out=fopen(bufor,"w");
  printf("fopen %s in function fprint_gridfile\n", bufor);

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
		      ix,iy,iz,
		      geom.xxvec[1],geom.xxvec[2],geom.xxvec[3],
		      xxcar[1],xxcar[2],xxcar[3],
		      xxsph[1],xxsph[2],xxsph[3]);

	    }
	}
    }

  fclose(out);

  return 0;
}


/*********************************************/
/* prints scalar quantities to scalars.dat */
/*********************************************/

int
fprint_scalars(ldouble t, ldouble *scalars, int nscalars)
{

#ifndef MPI
  calc_scalars(scalars,t);
  int iv;

#ifdef TIMEINSCALARSINSEC
      t=timeGU2CGS(t);
#endif

  fprintf(fout_scalars,"%e ",t);
  for(iv=0;iv<nscalars;iv++)
    fprintf(fout_scalars,"%.20e ",scalars[iv]);
  fprintf(fout_scalars,"\n");
  fflush(fout_scalars);
#endif //ifndef MPI

  return 0;
}


/*********************************************/
/* prints radial profiles to radNNNN.dat */
/*********************************************/
int
fprint_radprofiles(ldouble t, int nfile, char* folder, char* prefix)
{

      char bufor[50],bufor2[50];
      sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

      fout_radprofiles=fopen(bufor,"w");
      printf("fopen %s in function fprint_radprofiles\n", bufor);

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

          #ifdef PRECOMPUTE_MY2OUT
          get_xxout(ix, 0, 0, xxout);
          #else
	  get_xx(ix,0,0,xx);
          coco_N(xx,xxout,MYCOORDS,OUTCOORDS);
          #endif
	  
	  if(xxout[1]<rhorizonBL) continue;
	  fprintf(fout_radprofiles,"%e ",xxout[1]);
	  for(iv=0;iv<NRADPROFILES;iv++)
	    fprintf(fout_radprofiles,"%e ",profiles[iv][ix]);
	  fprintf(fout_radprofiles,"\n");
	}
      fclose(fout_radprofiles);
  
  return 0;
}
 
 
/*********************************************/
/* prints theta profiles to thNNNN.dat */
/*********************************************/

int
fprint_thprofiles(ldouble t, int nfile, char* folder, char* prefix)
{

  char bufor[50],bufor2[50];
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);

  FILE *fout_thprofiles=fopen(bufor,"w");
  printf("fopen %s in function fprint_thprofiles\n", bufor);

  int ix,iy,iv;

  //search for appropriate radial index
  ldouble xx[4],xxBL[4];
  ldouble radius=1.e3;
  #ifdef THPROFRADIUS
  radius=THPROFRADIUS;
  #endif
  for(ix=0;ix<NX;ix++)
    {
      #ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix, 0, 0, xxBL);
      #else
      get_xx(ix,0,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      #endif

      if(xxBL[1]>radius) break;
    }

  //calculating theta profiles
  ldouble profiles[NTHPROFILES][NY];
  calc_thetaprofiles(profiles);
  //printing th profiles  
  for(iy=0;iy<NY;iy++)
    {
      #ifdef PRECOMPUTE_MY2OUT
      get_xxout(ix, iy, 0, xxBL);
      #else
      get_xx(ix,iy,0,xx);
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      #endif
      
      fprintf(fout_thprofiles,"%e ",xxBL[2]);
      for(iv=0;iv<NTHPROFILES;iv++)
	fprintf(fout_thprofiles,"%e ",profiles[iv][iy]);
      fprintf(fout_thprofiles,"\n");
    }
  fclose(fout_thprofiles);
  
  return 0;
}

/*********************************************/
/* prints restart files */
/*********************************************/

int
fprint_restartfile(ldouble t, char* folder)
{
#ifdef DUMPS_WRITE_HDF5

  #ifdef MPI
  fprint_restartfile_mpi_hdf5(t, folder);
  #else
  fprint_restartfile_serial_hdf5(t, folder);
  #endif

#else
  
  #ifdef MPI
  fprint_restartfile_mpi(t,folder);
  #else 
  fprint_restartfile_bin(t,folder); 
  #endif

#endif
  
  return 0;
}


/*********************************************/
//parallel output to a single restart file
/*********************************************/

int 
fprint_restartfile_mpi(ldouble t, char* folder)
{
  #ifdef MPI

  if (PROCID == 0)
  {
    printf("Entering fprint_restartfile_mpi\n");
  }

  char bufor[250];

  //header
  if(PROCID==0)
    {
       sprintf(bufor,"%s/res%04d.head",folder,nfout1);
       fout1=fopen(bufor,"w"); 
       printf("fopen %s in function fprint_restartfile_mpi, PROCID = %d\n", bufor, PROCID);
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
  //printf("MPI_File_open in fprint_restartfile_mpi, PROCID = %d\n", PROCID);
  if (rc)
  {
    printf( "Unable to open/create file %s\n", bufor );fflush(stdout); exit(-1);
  }

  /***** first write all the indices ******/

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

  ldouble *pout;
  if((pout = (ldouble *)malloc(NX*NY*NZ*NV*sizeof(ldouble)))==NULL) my_err("malloc err. - fileop 1\n");
    
  for(ix=0;ix<NX;ix++)
    for(iy=0;iy<NY;iy++)
      for(iz=0;iz<NZ;iz++)
	for(iv=0;iv<NV;iv++)
	  {
	    pout[ix*NY*NZ*NV+iy*NZ*NV+iz*NV+iv]=get_u(p,iv,ix,iy,iz);


	  }

  MPI_File_write( cFile, pout, NX*NY*NZ*NV, MPI_LDOUBLE, &status );  
  MPI_File_close( &cFile );

  //maybe not needed
  MPI_Barrier(MPI_COMM_WORLD);

  //move reslast.head and reslast.dat symbolic links
  if(PROCID==0)
    {

      sprintf(bufor,"rm %s/reslast.dat",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.dat %s/reslast.dat",nfout1,folder);
      iv=system(bufor);

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
//serial binary restart file output
/*********************************************/

int 
fprint_restartfile_bin(ldouble t, char* folder)
{
  printf("Entering fprint_restartfile_bin\n");
  char bufor[250];
  
  //header
  if(PROCID==0)
    {
       sprintf(bufor,"%s/res%04d.head",folder,nfout1);
       fout1=fopen(bufor,"w"); 
       printf("fopen %s in function fprint_restartfile_bin\n", bufor);
       sprintf(bufor,"## %5d %5d %10.6e %5d %5d %5d %5d\n",nfout1,nfout2,t,PROBLEM,TNX,TNY,TNZ);
       fprintf(fout1,"%s",bufor);
       fclose(fout1);
    }

  //body
  sprintf(bufor,"%s/res%04d.dat",folder,nfout1);
  fout1=fopen(bufor,"wb"); 
  printf("fopen %s in function fprint_restartfile_bin\n", bufor);

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

	  fwrite(ppout,sizeof(ldouble),NV,fout1);
	}

  fclose(fout1);

  //move reslast.head and reslast.dat symbolic links  
  if(PROCID==0)
    {
      sprintf(bufor,"rm %s/reslast.dat",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.dat %s/reslast.dat",nfout1,folder);
      iv=system(bufor);

      sprintf(bufor,"rm %s/reslast.head",folder);
      iv=system(bufor);
      sprintf(bufor,"ln -s res%04d.head %s/reslast.head",nfout1,folder);
      iv=system(bufor);    
    }

  return 0;
}


int //parallel MPI hdf5 output
fprint_restartfile_mpi_hdf5(ldouble t, char* folder)
{
#if defined DUMPS_WRITE_HDF5 && defined MPI

  #ifndef FOLDER_HDF5
  #define FOLDER_HDF5 "./dumps"
  #endif

  if (PROCID == 0)
  {
    printf("Entering fprint_restartfile_mpi_hdf5:  FOLDER_HDF5 = %s\n", FOLDER_HDF5);
  }

  MPI_Comm comm  = MPI_COMM_WORLD;
  MPI_Info info  = MPI_INFO_NULL;

  char bufor[250];
  
  // Write out header information in group HEADER in HDF5 file
  // Only PROCID=0 does this, but all processes open the file, group and dataspace

  hid_t dumps_file_id, dumps_group_id, dumps_dataspace_scalar, dumps_dataspace_array, dataspace_array, dumps_memspace_array, dumps_dataset_int, dumps_dataset_double, dumps_dataset_array, dumps_attribute_id, plist_id;
  hsize_t dims_h5[3], chunk_dims[3], offset_3d[3], stride_3d[3], count_3d[3], block_3d[3];
  herr_t status;
    
  int file_number = nfout1, file_avg = nfout2, problem_number = PROBLEM, nxx = TNX, nyy = TNY, nzz = TNZ, nprimitives = NV;
    
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
    
  char fname_h5[256];
  sprintf(fname_h5, "%s/res%04d.h5", FOLDER_HDF5, nfout1);
  dumps_file_id = H5Fcreate(fname_h5, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  dumps_group_id = H5Gcreate2(dumps_file_id, "/HEADER", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
  dumps_dataspace_scalar = H5Screate(H5S_SCALAR);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/FILE_NUMBER", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (PROCID == 0) status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_number);
  status = H5Dclose(dumps_dataset_int);
    
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/FILE_AVG", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (PROCID == 0) status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_avg);
  status = H5Dclose(dumps_dataset_int);
    
  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/HEADER/TIME", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (PROCID == 0) status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &t);
  status = H5Dclose(dumps_dataset_double);
    
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/PROBLEM_NUMBER", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (PROCID == 0) status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &problem_number);
  status = H5Dclose(dumps_dataset_int);
    
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/NX", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (PROCID == 0) status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nxx);
  status = H5Dclose(dumps_dataset_int);
    
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/NY", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (PROCID == 0) status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nyy);
  status = H5Dclose(dumps_dataset_int);
    
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/NZ", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (PROCID == 0) status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nzz);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/NPRIM", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (PROCID == 0) status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nprimitives);
  status = H5Dclose(dumps_dataset_int);
  
  status = H5Pclose(plist_id);
  status = H5Gclose(dumps_group_id);

  // Now work on the body of the dumps file. First set up the array and chunking information

  int ret, ix, iy, iz, iv, i, ic, gix, giy, giz;
  int ixx[NX][NY][NZ], iyy[NX][NY][NZ], izz[NX][NY][NZ];
  
  for(ix = 0; ix < NX; ix++)
    for(iy = 0; iy < NY; iy++)
      for(iz = 0; iz < NZ; iz++)
      {
        mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);

        ixx[ix][iy][iz] = gix;
        iyy[ix][iy][iz] = giy;
        izz[ix][iy][iz] = giz;
      }

  dims_h5[0] = TNX;
  dims_h5[1] = TNY;
  dims_h5[2] = TNZ;

  chunk_dims[0] = TNX / NTX;
  chunk_dims[1] = TNY / NTY;
  chunk_dims[2] = TNZ / NTZ;

  count_3d[0] = 1;
  count_3d[1] = 1;
  count_3d[2] = 1;

  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;

  block_3d[0] = chunk_dims[0];
  block_3d[1] = chunk_dims[1];
  block_3d[2] = chunk_dims[2];

  offset_3d[0] = ixx[0][0][0];
  offset_3d[1] = iyy[0][0][0];
  offset_3d[2] = izz[0][0][0];

  /*
  // Save indices in HDF5 file. Is this needed? If not, get rid of ixx, iyy, izz arrays
    
  for (iv = 0; iv < 3; iv++)
  {
    // This sequence of H5 commands, which is taken from an example code, looks very complicated, but it works. Probably could be streamlined a lot.

    dumps_dataspace_array = H5Screate_simple(3, dims_h5, NULL);
    dumps_memspace_array = H5Screate_simple(3, chunk_dims, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 3, chunk_dims);

    if (iv == 0)
    {
      dumps_dataset_array = H5Dcreate2(dumps_file_id, "/GIX", H5T_STD_I32BE, dumps_dataspace_array, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    }
    else if (iv == 1)
    { 
      dumps_dataset_array = H5Dcreate2(dumps_file_id, "/GIY", H5T_STD_I32BE, dumps_dataspace_array, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    } 
    else
    { 
      dumps_dataset_array = H5Dcreate2(dumps_file_id, "/GIZ", H5T_STD_I32BE, dumps_dataspace_array, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    } 

    H5Pclose(plist_id);
    H5Sclose(dumps_dataspace_array);

    dumps_dataspace_array = H5Dget_space(dumps_dataset_array);
    status = H5Sselect_hyperslab(dumps_dataspace_array, H5S_SELECT_SET, offset_3d, stride_3d, count_3d, block_3d);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    if (iv == 0)
    {
      status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_INT, dumps_memspace_array, dumps_dataspace_array, plist_id, ixx);
    }
    else if (iv == 1)
    { 
      status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_INT, dumps_memspace_array, dumps_dataspace_array, plist_id, iyy);
    } 
    else
    { 
      status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_INT, dumps_memspace_array, dumps_dataspace_array, plist_id, izz);
    } 

    status = H5Dclose(dumps_dataset_array);
  }
   */

  // Next work on the primitives one by one

  double primitive[NX][NY][NZ];
  int iattribute = 0;
  char prim_name[256], attribute_name[256];

  //dumps_dataspace_array = H5Screate_simple(3, dims_h5, NULL);
  //dumps_memspace_array = H5Screate_simple(3, chunk_dims, NULL);

  for(iv = 0; iv < NV; iv++)
  {

  // First copy the primitive to the 3D array

  for(ix = 0; ix < NX; ix++)
    for(iy = 0; iy < NY; iy++)
      for(iz = 0; iz < NZ; iz++)
      {
        primitive[ix][iy][iz] = get_u(p, iv, ix, iy, iz);
      }
    
  // Then find the primitive name and save in the hdf5 file

    get_prim_name(prim_name, iv);
    //printf("  prim_name: %s\n", prim_name);

    dumps_dataspace_array = H5Screate_simple(3, dims_h5, NULL);
    dumps_memspace_array = H5Screate_simple(3, chunk_dims, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 3, chunk_dims);

    dumps_dataset_array = H5Dcreate2(dumps_file_id, prim_name, H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(dumps_dataspace_array);

    dataspace_array = H5Dget_space(dumps_dataset_array);
    status = H5Sselect_hyperslab(dataspace_array, H5S_SELECT_SET, offset_3d, stride_3d, count_3d, block_3d);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, dumps_memspace_array, dataspace_array, plist_id, primitive);
    status = H5Dclose(dumps_dataset_array);
  }  

  status = H5Pclose(plist_id);
  status = H5Sclose(dumps_dataspace_scalar);
  status = H5Sclose(dumps_memspace_array);
  status = H5Fclose (dumps_file_id);

#ifdef RESTARTOUTPUTINBL
  printf("\n RESTARTOUTPUTINBL not allowed in hdf5 output!\n";
  exit(1);
#endif

#ifdef OVERWRITEENTROPYINRESTARTFILESWITHNEGEHEATING
  printf("\n OVERWRITEENTROPYINRESTARTFILESWITHNEGEHEATING not allowed in hdf5 output!\n";
  exit(1);
#endif

  // Finally, redefine restart files to the current dump

  if (PROCID == 0)
  {
    sprintf(bufor, "rm %s/reslast.h5", FOLDER_HDF5);
    iv=system(bufor);
    sprintf(bufor,"ln -s res%04d.h5 %s/reslast.h5", nfout1, FOLDER_HDF5);
    iv=system(bufor);
  }

#endif  // DUMPS_WRITE_HDF5

  return 0;
}

							  
/*********************************************/
/*********************************************/

int //serial hdf5 output
fprint_restartfile_serial_hdf5(ldouble t, char* folder)
{
#ifdef DUMPS_WRITE_HDF5

  char bufor[250];
  
  #ifndef FOLDER_HDF5
  #define FOLDER_HDF5 "./dumps"
  #endif

  printf("Entering fprint_restartfile_serial_hdf5: FOLDER_HDF5 = %s\n", FOLDER_HDF5);

  // Write out header information in group HEADER in HDF5 file

  hid_t dumps_file_id, dumps_group_id, dumps_dataspace_scalar, dumps_dataspace_array, dumps_dataset_int, dumps_dataset_double, dumps_dataset_array, dumps_attribute_id;
  hsize_t dims_h5[3];
  herr_t status;
    
  int file_number = nfout1, file_avg = nfout2, problem_number = PROBLEM, nxx = TNX, nyy = TNY, nzz = TNZ, nprimitives = NV;
    
  dims_h5[0] = NX;
  dims_h5[1] = NY;
  dims_h5[2] = NZ;
    
  char fname_h5[256];
  sprintf(fname_h5, "%s/res%04d.h5", FOLDER_HDF5, nfout1);
  dumps_file_id = H5Fcreate (fname_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  dumps_group_id = H5Gcreate2(dumps_file_id, "/HEADER", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
  dumps_dataspace_scalar = H5Screate(H5S_SCALAR);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/FILE_NUMBER", H5T_STD_I32BE, dumps_dataspace_scalar,
				 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &file_number);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/FILE_AVG", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &file_avg);
  status = H5Dclose(dumps_dataset_int);
 
  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/HEADER/TIME", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &t);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/PROBLEM_NUMBER", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &problem_number);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/NX", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &nxx);
  status = H5Dclose(dumps_dataset_int);
  
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/NY", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &nyy);
  status = H5Dclose(dumps_dataset_int);
  
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/NZ", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &nzz);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/NPRIM", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &nprimitives);
  status = H5Dclose(dumps_dataset_int);
  
  status = H5Gclose(dumps_group_id);

  // Now work on the body of the dumps file

  // Create arrays of indices

  int ret, ix, iy, iz, iv, i, ic, gix, giy, giz;
  int ixx[NX][NY][NZ], iyy[NX][NY][NZ], izz[NX][NY][NZ];
  
  for(ix = 0; ix < NX; ix++)
    for(iy = 0; iy < NY; iy++)
      for(iz = 0; iz < NZ; iz++)
      {
        mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);

        ixx[ix][iy][iz] = gix;
        iyy[ix][iy][iz] = giy;
        izz[ix][iy][iz] = giz;
      }

  // Save indices in HDF5 file. Is this needed? If not, get rid of ixx, iyy, izz arrays
    
  dumps_dataspace_array = H5Screate_simple(3, dims_h5, NULL);
  

  /*  
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/GIX", H5T_STD_I32BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    ixx);
  status = H5Dclose(dumps_dataset_array);
    
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/GIY", H5T_STD_I32BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    iyy);
  status = H5Dclose(dumps_dataset_array);
    
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/GIZ", H5T_STD_I32BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    izz);
  status = H5Dclose(dumps_dataset_array);
   */
    
  // Next work on the primitives one by one

  double primitive[NX][NY][NZ];
  int iattribute = 0;
  char prim_name[256], attribute_name[256];

  for(iv = 0; iv < NV; iv++)
  {

  // First copy the primitive to the 3D array

    for(ix = 0; ix < NX; ix++)
      for(iy = 0; iy < NY; iy++)
	for(iz = 0; iz < NZ; iz++)
	{
	  primitive[ix][iy][iz] = get_u(p, iv, ix, iy, iz);
	}
    
   // Then find the primitive name and save in the hdf5 file

    get_prim_name(prim_name, iv);
    //printf("  prim_name: %s\n", prim_name);
    
    dumps_dataset_array = H5Dcreate2(dumps_file_id, prim_name, H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT, &(primitive[0][0][0]));
    status = H5Dclose(dumps_dataset_array);
  }  

  status = H5Sclose(dumps_dataspace_scalar);
  status = H5Sclose(dumps_dataspace_array);
  status = H5Fclose (dumps_file_id);

#ifdef RESTARTOUTPUTINBL
  printf("\n RESTARTOUTPUTINBL not allowed in hdf5 output!\n";
  exit(1);
#endif

#ifdef OVERWRITEENTROPYINRESTARTFILESWITHNEGEHEATING
  printf("\n OVERWRITEENTROPYINRESTARTFILESWITHNEGEHEATING not allowed in hdf5 output!\n";
  exit(1);
#endif

  // Finally, redefine restart files to the current dump

  sprintf(bufor, "rm %s/reslast.h5", FOLDER_HDF5);
  iv=system(bufor);
  sprintf(bufor,"ln -s res%04d.h5 %s/reslast.h5", nfout1, FOLDER_HDF5);
  iv=system(bufor);

#endif  // DUMPS_WRITE_HDF5

  return 0;
}


/*********************************************/
/* reads dump file */
/* puts conserved into the memory */
/* converts them to primitives */
/*********************************************/

int
fread_restartfile(int nout1, char* folder,ldouble *t)
{


  int ret;
     
#ifdef DUMPS_READ_HDF5
  
  #ifdef MPI
  ret = fread_restartfile_mpi_hdf5(nout1, folder, t);
  #else
  ret = fread_restartfile_serial_hdf5(nout1, folder, t);
  #endif

#else

  #ifdef MPI
  ret=fread_restartfile_mpi(nout1,folder,t);
  #else //no MPI 
  ret=fread_restartfile_bin(nout1,folder,t);  
  #endif

#endif
  
  return ret;
}


/*********************************************/
//Read binary restart file to single thread
/*********************************************/

int 
fread_restartfile_bin(int nout1, char *folder, ldouble *t)
{
  printf("Entering fread_restartfile_bin\n");
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
  printf("fopen %s in function fread_restartfile_bin\n", fnamehead);
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
  printf("fopen %s in function fread_restartfile_bin\n", fname);

  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;

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


      // restart options
      // TODO currently can do EITHER RESTARTFROMMHD (for radiation and/or electrons)
      // OR RESTARTFROMNORELEL (FROM radiation+electrons to relelectrons)
      // SHOULD add more options
#if defined(RESTARTFROMMHD) //restart from mhd only
#ifdef RADIATION

      int nvold = 9;
      ldouble ppold[nvold]; 
      ret=fread(ppold,sizeof(ldouble),nvold,fdump);

      // copy over mhd primitives
      int ie;
      for (ie=0; ie<nvold; ie++) pp[ie] = ppold[ie]; 

      // initialize radiation primitives
      //ldouble Erad=INITERAD;
      ldouble Erad=pp[UU]*INITURADFRAC;   //initial radiation energy density
      pp[EE]=Erad;
      pp[FX]=pp[VX]; //initial radiation velocity is same as fluid     
      pp[FY]=pp[VY];
      pp[FZ]=pp[VZ];

#ifdef EVOLVEPHOTONNUMBER
      pp[NF]=calc_NFfromE(Erad); //initially blackbody
#endif

#ifdef EVOLVEELECTRONS
      // initialize electron primitives
      ldouble rhogas=pp[RHO];
      ldouble ue=(INITUEFRAC)*pp[UU];
      ldouble ui=(1.-INITUEFRAC)*pp[UU];

      pp[ENTRE]=calc_Sefromrhou(rhogas,ue,ELECTRONS); //TODO -- non hydrogen?
      pp[ENTRI]=calc_Sefromrhou(rhogas,ui,IONS);
#endif
#endif //RADIATION


#elif defined(RESTARTFROMNORELEL)
#ifdef RELELECTRONS
      int nvold=NV-NRELBIN;

#ifdef RESTARTFROMNORELEL_NOCOMPT
      nvold += 1;
#endif

      ldouble ppold[nvold]; 
      ret=fread(ppold,sizeof(ldouble),nvold,fdump);

      int ie;
      for (ie=0; ie<NVHD-NRELBIN; ie++) pp[ie] = ppold[ie];
      for (ie=0; ie<NRELBIN; ie++) pp[NVHD-NRELBIN+ie] = 0.0; 
      for (ie=0; ie<NV-NVHD; ie++) pp[NVHD+ie] = ppold[NVHD-NRELBIN+ie];

#endif
      
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

      
      #ifdef CONSISTENTGAMMA
      ldouble Te,Ti;
      calc_PEQ_Teifrompp(pp,&Te,&Ti,ix,iy,iz);
      ldouble gamma=calc_gammaintfromTei(Te,Ti);
      set_u_scalar(gammagas,ix,iy,iz,gamma);
      #endif

#ifdef RESTARTFROMMHD
#ifdef EVOLVEELECTRONS
      pp[ENTRE]=calc_SefromrhoT(rhogas,Te,ELECTRONS);
      pp[ENTRI]=calc_SefromrhoT(rhogas,Ti,IONS);
#endif
#endif

      p2u(pp,uu,&geom);


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
//read restart file to individual MPI tile
/*********************************************/

int 
fread_restartfile_mpi(int nout1, char *folder, ldouble *t)
{
  #ifdef MPI

  if (PROCID == 0)
  {
    printf("Entering fread_restartfile_mpi\n");
  }

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

  /***********/
  //Read header file from process 0 if it exists. If not return for start from scratch
 
  int exists = 0;

  if (PROCID == 0)
  {
    FILE *fdump = fopen(fnamehead,"r");
    printf("fopen %s in function fread_restartfile_mpi, PROCID = %d\n", fnamehead, PROCID);
    exists = (fdump != NULL); // cast file open to an integer for broadcasting
    if(exists == 1)
    {
      //reading parameters, mostly time
      int intpar[6];
      ret=fscanf(fdump,"## %d %d %lf %d %d %d %d\n",&intpar[0],&intpar[1],t,&intpar[2],&intpar[3],&intpar[4],&intpar[5]);
      if(PROCID==0)
      printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n", fname,intpar[0],*t,intpar[2],intpar[3],intpar[4],intpar[5]); 
      nfout1=intpar[0]+1; //global file no.
      nfout2=intpar[1]; //global file no. for avg
      fclose(fdump);
    }
  }
 
  // Broadcast `exists` from process 0 to all other processors using a blocking function         

  //printf("Before MPI_Bcast: PROCID, exists, nfout1, nfout2, t = %d %d %d %d %e\n", PROCID, exists, nfout1, nfout2, *t);
  int intexchange[3];
  intexchange[0] = exists;
  intexchange[1] = nfout1;
  intexchange[2] = nfout2;
  MPI_Bcast(intexchange, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(t, 1, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
  exists = intexchange[0];
  nfout1 = intexchange[1];
  nfout2 = intexchange[2];
  //printf("After MPI_Bcast: PROCID, exists, nfout1, nfout2, t = %d %d %d %d %e\n", PROCID, exists, nfout1, nfout2, *t);

  if (exists == 0)
  {
    return 1; //request start from scratch
  }
  
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
  //printf("MPI_File_open in fread_restartfile_mpi, PROCID = %d\n", PROCID);
  if (rc)
  {
    printf( "Unable to open/create file %s\n", fname );fflush(stdout); exit(-1);
  }

  /***** first read ALL the indices ******/

 int nvold;
#if defined(RESTARTFROMMHD)
  nvold=9;
#elif defined(RESTARTFROMNORELEL) //TODO can only do EITHRE RESTARTFROMMHD OR RESTARTFROMNORELEL
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

      /***** then read primitives in the same order ******/

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

#if defined(RESTARTFROMMHD) //restart from mhd only
#ifdef RELELECTRONS
              my_err("RESTARTFROMMHD does not work with RELELECTRONS currently!\n");
#endif
              // copy over mhd primitives
              int ie;
              for(ie=0; ie<6; ie++)
                set_u(p,ie,ix,iy,iz,pout[ppos+ie]); //density, energy, velocities, entropy
	      for(ie=0; ie<3; ie++)
		set_u(p,NVHD+ie,ix,iy,iz,pout[ppos+6+ie]); //magnetic field

#ifdef RADIATION		       
	      // initialize radiation primitives
	      //ldouble Erad=INITERAD;
	      ldouble Erad=pout[ppos+UU]*INITURADFRAC;   //initial radiation energy density

	      set_u(p,EE,ix,iy,iz,Erad);
	      set_u(p,FX,ix,iy,iz,pout[ppos+VX]); //initial radiation velocity is same as fluid     
	      set_u(p,FY,ix,iy,iz,pout[ppos+VY]);
	      set_u(p,FZ,ix,iy,iz,pout[ppos+VZ]);

#ifdef EVOLVEPHOTONNUMBER
	      ldouble Nrad=calc_NFfromE(Erad); //initially blackbody
	      set_u(p,NF,ix,iy,iz,Nrad);
#endif

#ifdef EVOLVEELECTRONS
	      // initialize electron primitives
	      ldouble rhogas=pout[ppos+RHO];
	      ldouble ue=(INITUEFRAC)*pout[ppos+UU];
	      ldouble ui=(1.-INITUEFRAC)*pout[ppos+UU];
	      ldouble Se=calc_Sefromrhou(rhogas,ue,ELECTRONS); //TODO -- non hydrogen?
	      ldouble Si=calc_Sefromrhou(rhogas,ui,IONS);
	      set_u(p,ENTRE,ix,iy,iz,Se);
	      set_u(p,ENTRI,ix,iy,iz,Si);
#endif
#endif //RADIATION
	  
#elif defined(RESTARTFROMNORELEL) // TODO only an option if NOT RESTARTFROMMHD
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

#ifdef RESTARTFROMMHD
#ifdef EVOLVEELECTRONS
              Se=calc_SefromrhoT(rhogas,Te,ELECTRONS);
	      Si=calc_SefromrhoT(rhogas,Ti,IONS);
	      set_u(p,ENTRE,ix,iy,iz,Se);
	      set_u(p,ENTRI,ix,iy,iz,Si);
#endif
#endif

	      p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	      

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

	 
/*********************************************/
/*********************************************/

int //parallel MPI hdf5 input
fread_restartfile_mpi_hdf5(int nout1, char *folder, ldouble *t)
{
#if defined DUMPS_READ_HDF5 && defined MPI

  #ifndef FOLDER_HDF5
  #define FOLDER_HDF5 "./dumps"
  #endif

  if (PROCID == 0)
  {
    printf("Entering fread_restartfile_mpi_hdf5:  FOLDER_HDF5 = %s\n", FOLDER_HDF5);
  }

  MPI_Comm comm  = MPI_COMM_WORLD;
  MPI_Info info  = MPI_INFO_NULL;

  hid_t dumps_file_id, dumps_group_id, dumps_dataspace_scalar, dumps_dataspace_array, dataspace_array, dumps_memspace_array, dumps_dataset_int, dumps_dataset_double, dumps_dataset_array, dumps_attribute_id, plist_id;
  hsize_t dims_h5[3], chunk_dims[3], offset_3d[3], stride_3d[3], count_3d[3], block_3d[3];
  herr_t status;
    
  int file_number, file_avg, problem_number, nxx, nyy, nzz, nprimitives;
    
  char fname_h5[400];

  if(nout1>=0)
  {
    sprintf(fname_h5, "%s/res%04d.h5", FOLDER_HDF5, nout1);
  }
  else
  {
    sprintf(fname_h5, "%s/reslast.h5", FOLDER_HDF5);
  }

  int exists;

  // Check if hdf5 file exists from process 0
  if (PROCID == 0) {
    FILE *fhdf5 = fopen(fname_h5, "r");
    printf("fopen %s in fread_restartfile_mpi_hdf5: PROCID = %d\n", fname_h5, PROCID);
    exists = (fhdf5 != NULL); // cast to an integer so it can be safely broadcasted
    if (exists == 1) fclose(fhdf5);
  }

  // Broadcast `exists` from process 0 to all other processors using a blocking function
  MPI_Bcast(&exists, 1, MPI_INT, 0, comm);

  // Request start from scratch if hdf5 file does not exist
  if(!exists)
    {
      return 1;
    }

  /***********/
  // Open hdf5 file and read header information from group HEADER

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, comm, info);
  dumps_file_id = H5Fopen(fname_h5, H5F_ACC_RDONLY, plist_id);
  //printf("H5Fopen in function fread_restartfile_mpi_hdf5, PROCID = %d\n", PROCID);
  H5Pclose(plist_id);

  if(&dumps_file_id==NULL)
  {
    printf("Input hdf5 file not available! PROCID = %d\n", PROCID);
    exit(1);
  }

  dumps_group_id = H5Gopen2(dumps_file_id, "/HEADER", H5P_DEFAULT);
  dumps_dataspace_scalar = H5Screate(H5S_SCALAR);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "FILE_NUMBER", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_number);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "FILE_AVG", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_avg);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_double = H5Dopen2(dumps_group_id, "TIME", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "PROBLEM_NUMBER", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &problem_number);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "NX", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nxx);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "NY", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nyy);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "NZ", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nzz);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "NPRIM", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nprimitives);
  status = H5Dclose(dumps_dataset_int);
  
  status = H5Pclose(plist_id);
  status = H5Gclose(dumps_group_id);
  status = H5Sclose(dumps_dataspace_scalar);

  if (PROCID == 0) printf("PROCID: %d  restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d and NV: %d\n", PROCID, fname_h5, file_number, *t, problem_number, nxx, nyy, nzz, nprimitives);

  nfout1 = file_number + 1; //global file no.
  nfout2 = file_avg; //global file no. for avg

  if (nxx != TNX || nyy != TNY || nzz != TNZ)
  {
    printf("Array sizes do not match!!\nnxx nyy nzz TNX TNY TNZ: %d %d %d %d %d %d\n", nxx, nyy, nzz, TNX, TNY, TNZ);
    exit(1);
  }

  if (nprimitives != NV)
  {
    printf("Number of primitives does not match!!\nprimitves NV: %d %d\n", nprimitives, NV);
    exit(1);
  }
  
  /***********/
  // Now read the primitives one by one and save in array p

  int ret, ix, iy, iz, iv, i, ic, gix, giy, giz;

  mpi_local2globalidx(0, 0, 0, &gix, &giy, &giz);

  dims_h5[0] = TNX;
  dims_h5[1] = TNY;
  dims_h5[2] = TNZ;

  chunk_dims[0] = TNX / NTX;
  chunk_dims[1] = TNY / NTY;
  chunk_dims[2] = TNZ / NTZ;

  count_3d[0] = 1;
  count_3d[1] = 1;
  count_3d[2] = 1;

  stride_3d[0] = 1;
  stride_3d[1] = 1;
  stride_3d[2] = 1;

  block_3d[0] = chunk_dims[0];
  block_3d[1] = chunk_dims[1];
  block_3d[2] = chunk_dims[2];

  offset_3d[0] = gix;
  offset_3d[1] = giy;
  offset_3d[2] = giz;
    
  double primitive[NX][NY][NZ];
  char prim_name[256];

  for(iv = 0; iv < NV; iv++)
  {
    // Find the primitive name and read primitive values from the hdf5 file

    get_prim_name(prim_name, iv);
    //printf("  prim_name: %s\n", prim_name);

    dumps_dataspace_array = H5Screate_simple(3, dims_h5, NULL);
    dumps_memspace_array = H5Screate_simple(3, chunk_dims, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 3, chunk_dims);

    dumps_dataset_array = H5Dopen2(dumps_file_id, prim_name, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(dumps_dataspace_array);

    dataspace_array = H5Dget_space(dumps_dataset_array);
    status = H5Sselect_hyperslab(dataspace_array, H5S_SELECT_SET, offset_3d, stride_3d, count_3d, block_3d);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dread(dumps_dataset_array, H5T_NATIVE_DOUBLE, dumps_memspace_array, dataspace_array, plist_id, primitive);
    status = H5Dclose(dumps_dataset_array);

    // Copy the primitives to p array

    for(ix = 0; ix < NX; ix++)
    {
      for(iy = 0; iy < NY; iy++)
      {
	for(iz = 0; iz < NZ; iz++)
	{
          set_u(p, iv, ix, iy, iz, primitive[ix][iy][iz]);
	}
      }
    }  
  }

  status = H5Pclose(plist_id);
  status = H5Sclose(dumps_memspace_array);
  status = H5Fclose (dumps_file_id);

  //printf("Done reading hdf5 file\n");

#ifdef RESTARTFROMNORELEL
  printf("RESTARTFROMNORELEL not yet available with HDF5!\n");
  exit(1);
#endif

  // Loop over cells, calculate geometry, modify primitives as needed, compute conserveds and save

  struct geometry geom;
  ldouble xxvec[4], xxvecout[4];
  ldouble uu[NV], pp[NV];

  for(ix = 0; ix < NX; ix++)
  {
    for(iy = 0; iy < NY; iy++)
    {
      for(iz = 0; iz < NZ; iz++)
      {
	for (iv = 0; iv < NV; iv++)
	{
	  pp[iv] = get_u(p, iv, ix, iy, iz);
	}

	fill_geometry(ix,iy,iz,&geom);

#ifdef RESTARTOUTPUTINBL
	struct geometry geomBL;
	fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	trans_pall_coco(pp, pp, BLCOORDS,MYCOORDS, geomBL.xxvec,&geomBL,&geom);
#endif

#ifdef CONVERTLOGTONOLOGWHENRESTARTING
	printf("CONVERTLOGTONOLOGWHENRESTARTING not yet available with HDF5!\n");
	exit(1);
#endif
      
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

	//save primitives and conserveds
	for(iv=0;iv<NV;iv++)    
	{
	  set_u(u,iv,ix,iy,iz,uu[iv]);
	  set_u(p,iv,ix,iy,iz,pp[iv]);
	}
      }
    }
  }

#endif  // DUMPS_READ_HDF5

  return 0;
}


/*********************************************/
/*********************************************/

int //serial hdf5 input
fread_restartfile_serial_hdf5(int nout1, char *folder, ldouble *t)
{
  #ifdef DUMPS_READ_HDF5

#ifndef FOLDER_HDF5
#define FOLDER_HDF5 "./dumps"
#endif

  printf("Entering fread_restartfile_serial_hdf5:  FOLDER_HDF5 = %s\n", FOLDER_HDF5);

  hid_t dumps_file_id, dumps_group_id, dumps_dataspace_scalar, dumps_dataspace_array, dumps_dataset_int, dumps_dataset_double, dumps_dataset_array, dumps_attribute_id;
  hsize_t dims_h5[3];
  herr_t status;
    
  int file_number, file_avg, problem_number, nxx, nyy, nzz, nprimitives;
    
  dims_h5[0] = NX;
  dims_h5[1] = NY;
  dims_h5[2] = NZ;
    
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz;
  char fname_h5[400];

  if(nout1>=0)
  {
    sprintf(fname_h5, "%s/res%04d.h5", FOLDER_HDF5, nout1);
  }
  else
  {
    sprintf(fname_h5, "%s/reslast.h5", FOLDER_HDF5);
  }

  // Check if hdf5 file exists
  FILE *fhdf5;
  fhdf5=fopen(fname_h5,"r");
  printf("fopen %s in function fread_restartfile_mpi_hdf5, PROCID = %d\n", fname_h5, PROCID);
  if(fhdf5==NULL) 
  {
    return 1; //request start from scratch
  }

  /***********/
  // Open hdf5 file and read header information from group HEADER

  dumps_file_id = H5Fopen (fname_h5, H5F_ACC_RDWR, H5P_DEFAULT);
  //printf("H5Fopen in function fread_restartfile_serial_hdf5, PROCID = %d\n", PROCID);
  if(&dumps_file_id==NULL)
  {
    printf("Input hdf5 file not available!  PROCID = %d\n", PROCID);
    exit(1);
  }

  dumps_group_id = H5Gopen2(dumps_file_id, "/HEADER", H5P_DEFAULT);
  dumps_dataspace_scalar = H5Screate(H5S_SCALAR);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "FILE_NUMBER", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_number);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "FILE_AVG", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &file_avg);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_double = H5Dopen2(dumps_group_id, "TIME", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "PROBLEM_NUMBER", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &problem_number);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "NX", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nxx);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "NY", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nyy);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "NZ", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nzz);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dopen2(dumps_group_id, "NPRIM", H5P_DEFAULT);
  status = H5Dread(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nprimitives);
  status = H5Dclose(dumps_dataset_int);
  
  status = H5Gclose(dumps_group_id);
  status = H5Sclose(dumps_dataspace_scalar);

  printf("restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d and NV: %d\n", fname_h5, file_number, *t, problem_number, nxx, nyy, nzz, nprimitives);
  nfout1 = file_number + 1; //global file no.
  nfout2 = file_avg; //global file no. for avg

  if (nxx != NX || nyy != NY || nzz != NZ)
  {
    printf("Array sizes do not match!!\nnxx nyy nzz NX NY NZ: %d %d %d %d %d %d\n", nxx, nyy, nzz, NX, NY, NZ);
    exit(1);
  }

  if (nprimitives != NV)
  {
    printf("Number of primitives does not match!!\nnxx nprimitives NV: %d %d\n", nprimitives, NV);
    exit(1);
  }
  
  /***********/
  // Now read the primitives one by one and save in array p

  double primitive[NX][NY][NZ];
  char prim_name[256];

  dumps_dataspace_array = H5Screate_simple(3, dims_h5, NULL);

  for(iv = 0; iv < NV; iv++)
  {
    // Find the primitive name and read primitive values from the hdf5 file

    get_prim_name(prim_name, iv);
    //printf("  prim_name: %s\n", prim_name);

    dumps_dataset_array = H5Dopen2(dumps_file_id, prim_name, H5P_DEFAULT);
    status = H5Dread(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, primitive);
    status = H5Dclose(dumps_dataset_array);

    // Copy the primitives to p array

    for(ix = 0; ix < NX; ix++)
    {
      for(iy = 0; iy < NY; iy++)
      {
	for(iz = 0; iz < NZ; iz++)
	{
          set_u(p, iv, ix, iy, iz, primitive[ix][iy][iz]);
	}
      }
    }  
  }

  status = H5Sclose(dumps_dataspace_array);
  status = H5Fclose (dumps_file_id);

  //printf("Done reading hdf5 file\n");

#ifdef RESTARTFROMNORELEL
  printf("RESTARTFROMNORELEL not yet available with HDF5!\n");
  exit(1);
#endif

  // Loop over cells, calculate geometry, modify primitives as needed, compute conserveds and save

  struct geometry geom;
  ldouble xxvec[4], xxvecout[4];
  ldouble uu[NV], pp[NV];

  for(ix = 0; ix < NX; ix++)
  {
    for(iy = 0; iy < NY; iy++)
    {
      for(iz = 0; iz < NZ; iz++)
      {
	for (iv = 0; iv < NV; iv++)
	{
	  pp[iv] = get_u(p, iv, ix, iy, iz);
	}

	fill_geometry(ix,iy,iz,&geom);

#ifdef RESTARTOUTPUTINBL
	struct geometry geomBL;
	fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
	trans_pall_coco(pp, pp, BLCOORDS,MYCOORDS, geomBL.xxvec,&geomBL,&geom);
#endif

#ifdef CONVERTLOGTONOLOGWHENRESTARTING
	printf("CONVERTLOGTONOLOGWHENRESTARTING not yet available with HDF5!\n");
	exit(1);
#endif
      
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

	//save primitives and conserveds
	for(iv=0;iv<NV;iv++)    
	{
	  set_u(u,iv,ix,iy,iz,uu[iv]);
	  set_u(p,iv,ix,iy,iz,pp[iv]);
	}
      }
    }
  }

  #endif  // DUMPS_READ_HDF5

  return 0;
}

/*********************************************/
/* prints avg files */
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
//parallel output to a single avgfile
/*********************************************/

int
fprint_avgfile_mpi(ldouble t, char* folder, char* prefix)
{
  #ifdef MPI
  char bufor[250];
  
  //header
  if(PROCID==0)
    {
      sprintf(bufor,"%s/%s%04d.head",folder,prefix,nfout2);
      fout1=fopen(bufor,"w"); 
      printf("fopen %s in function fprint_avgfile_mpi, PROCID = %d\n", bufor, PROCID);
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
  //printf("MPI_File_open in fprint_avgfile_mpi, PROCID = %d\n", PROCID);
  if (rc)
  {
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
//serial binary output to avg file
/*********************************************/

int 
fprint_avgfile_bin(ldouble t, char* folder,char *prefix)
{
  char bufor[250];
  
  //header
  if(PROCID==0)
    {
      sprintf(bufor,"%s/%s%04d.head",folder,prefix,nfout2);
      fout1=fopen(bufor,"w"); 
      printf("fopen %s in function fprint_avgfile_bin, PROCID = %d\n", bufor, PROCID);
      sprintf(bufor,"## %5d %10.6e %10.6e %10.6e\n",nfout2,t-avgtime,t,avgtime);
      fprintf(fout1,"%s",bufor);
      fclose(fout1);
    }

  //body
  sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfout2);
  fout1=fopen(bufor,"wb"); 
  printf("fopen %s in function fprint_avgfile_bin, PROCID = %d\n", bufor, PROCID);

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

							  
/*********************************************/
/* reads avg dump file */
/* puts conserved into the memory */
/* converts them to primitives */
/*********************************************/
int
fread_avgfile(int nout1, char* base,ldouble *pavg, ldouble *dt,ldouble *t)
{
  char bufor[250];

  #ifdef MPI
  my_err("fread_avgfile should not be used with MPI\n");
  exit(1);
  
  #else //no MPI

  fread_avgfile_bin(nout1,base,pavg,dt,t);

  #endif
  
  return 0;
}

int 
fread_avgfile_bin(int nout1, char *base, ldouble *pavg, ldouble *dt, ldouble *t)
{
  int ret, ix,iy,iz,iv,i,ic,gix,giy,giz;
  char fname[40],fnamehead[40];


  printf("%s%04d.dat\n",base,nout1);
  printf("%s%04d.head\n",base,nout1);
  
  sprintf(fname,"%s%04d.dat",base,nout1);
  sprintf(fnamehead,"%s%04d.head",base,nout1);

  FILE *fdump;

  /***********/
  //header file
  fdump=fopen(fnamehead,"r");
  printf("fopen %s in function fread_avgfile_bin, PROCID = %d\n", fnamehead, PROCID);

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
  printf("fopen %s in function fread_avgfile_bin, PROCID = %d\n", fname, PROCID);

  struct geometry geom;
  ldouble xxvec[4],xxvecout[4];
  ldouble uu[NV],pp[NV],ftemp;
  char c;

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
/* wrapper for coordinate output */
/*********************************************/

int fprint_coordfile(char* folder,char* prefix)
{
#if (COORDOUTPUT>0)
  fprint_coordBL(folder,prefix);
#endif

  return 0;
}

/*********************************************/
/* prints BL coordinates,  */
/*********************************************/
                              
int fprint_coordBL(char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%sBL.dat",folder,prefix);
   FILE* fout1=fopen(bufor,"w");
   printf("fopen %s in function fprint_coordBL\n", bufor);

#ifdef COORDOUTPUT_HDF5
   hid_t coordBL_file_id, coordBL_dataspace_id, r_dataset_id, theta_dataset_id, phi_dataset_id;
   hsize_t dims_h5[3];
   herr_t status;
   
   char bufor_h5[256];
   sprintf(bufor_h5, "./%s/%sBL.h5", folder, prefix);
   printf("\n%s\n", bufor_h5);
   
   coordBL_file_id = H5Fcreate (bufor_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   
   dims_h5[0] = NX;
   dims_h5[1] = NY;
   dims_h5[2] = NZ;
   coordBL_dataspace_id = H5Screate_simple(3, dims_h5, NULL);
   
   double r_h5[NX][NY][NZ], theta_h5[NX][NY][NZ], phi_h5[NX][NY][NZ];
   
   r_dataset_id = H5Dcreate2(coordBL_file_id, "/r", H5T_IEEE_F32BE, coordBL_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   theta_dataset_id = H5Dcreate2(coordBL_file_id, "/theta", H5T_IEEE_F32BE, coordBL_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   phi_dataset_id = H5Dcreate2(coordBL_file_id, "/phi", H5T_IEEE_F32BE, coordBL_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
   
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
#ifdef COORDOUTPUT_HDF5
           r_h5[ix][iy][iz] = r;
           theta_h5[ix][iy][iz] = th;
           phi_h5[ix][iy][iz] = ph;
#endif
	       
	     }
	 }
     }

   fflush(fout1);
   fclose(fout1);

 #ifdef COORDOUTPUT_HDF5
   status = H5Dwrite(r_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     r_h5);
   status = H5Dwrite(theta_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     theta_h5);
   status = H5Dwrite(phi_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     phi_h5);
   
   status = H5Dclose(r_dataset_id);
   status = H5Dclose(theta_dataset_id);
   status = H5Dclose(phi_dataset_id);
   
   status = H5Sclose(coordBL_dataspace_id);
   status = H5Fclose (coordBL_file_id);
#endif
  
   return 0;
 }

                              
/*********************************************/
/* wrapper for simple output */
/*********************************************/
                              
int fprint_simplefile(ldouble t, int nfile, char* folder,char* prefix)
{
  //Cartesian output
#if (SIMOUTPUT==1)
  fprint_simplecart(t,nfile,folder,prefix);
#endif

  //Spherical output
#if (SIMOUTPUT==2)
  fprint_simplesph(t,nfile,folder,prefix);
#endif

#ifdef SIMOUTPUT_PHIAVG
  fprint_simple_phiavg(t, nfile, folder, prefix);
#endif

#ifdef SIMOUTPUT_PHICORR
  fprint_simple_phicorr(t, nfile, folder, prefix);
#endif
  return 0;
}

/*********************************************/
/* prints in ASCII indices, cart coordinates,*/
/* primitives, velocities in cartesian       */
/*********************************************/
//ANDREW -- TODO: update

int fprint_simplecart(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);
   fout1=fopen(bufor,"w");
   printf("fopen %s in function fprint_simplecart\n", bufor);

   //header
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

               #ifdef PRECOMPUTE_MY2OUT
               trans_pall_coco_my2out(pp,pp,&geom,&geomout);
               #else      
               trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomout);
               #endif
	       
	       ldouble rho=rhoGU2CGS(pp[RHO]);
	       #ifdef SIMOUTPUTINTERNAL
	       rho=pp[RHO];
	       #endif
	       ldouble temp=calc_PEQ_Tfromurho(pp[UU],pp[RHO],ix,iy,iz);
	       ldouble vel[4]={0,pp[VX],pp[VY],pp[VZ]};	
	       ldouble vx,vy,vz;
	       conv_vels(vel,vel,VELPRIM,VEL4,geomout.gg,geomout.GG);

	       vx=vel[1];
	       vy=vel[2];
	       vz=vel[3];
	       
	       //transform to cartesian
	       if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS   || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS ||
		   MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS || MYCOORDS==MKS3COORDS || MYCOORDS==JETCOORDS ||
		   MYCOORDS==MSPH1COORDS || MYCOORDS==MKER1COORDS)
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
	      if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS ||
		  MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
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

	      ldouble gammagas=calc_gammagas(pp, ix, iy, iz);
	      gammagas=pick_gammagas(ix,iy,iz);
	      fprintf(fout1,"%e %e %e %.5e %.5e %.5e %.5e %.5e %.5e %.5e",
	      	      uugas,ui,ue,
		      get_u_scalar(vischeating,ix,iy,iz),
		      get_u_scalar(vischeatingnegebalance,ix,iy,iz),
		      gammagas,Te,Ti,gammae,gammai);//17-21 with rad, 13-17 without
#endif //PROBLEM==115 || PROBLEM==135

	      fprintf(fout1,"\n");
	     }
	 }
     }

   fflush(fout1);
   fclose(fout1);

   return 0;
 }

                              
/*********************************************/
/* prints in ASCII & BL coordinates,  */
/*********************************************/
                              
int fprint_simplesph(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   #if defined(GRTRANSSIMOUTPUT)
   sprintf(bufor,"%s/%s%04d_grtnew.dat",folder,prefix,nfile);
   #elif defined(GRTRANSSIMOUTPUT_2)
   sprintf(bufor,"%s/%s%04d_simcgs.dat",folder,prefix,nfile);
   #else
   sprintf(bufor,"%s/%s%04d.dat",folder,prefix,nfile);
   #endif
   fout1=fopen(bufor,"w");
   printf("fopen %s in function fprint_simplesph\n", bufor);

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
  
   //ANDREW -- header for grtrans
#if defined (GRTRANSSIMOUTPUT) || defined(GRTRANSSIMOUTPUT_2)
   for(iix=-2;iix<NX;iix++)
     {
        fill_geometry_arb(iix,NY/2,NZ/2,&geomBL,OUTCOORDS);
	if(geomBL.xx>=rhorizonBL)
	  {
	    xmin=iix;
	    break;
	  }
     }
       
   if(NZ==1)
     fprintf(fout1,"%.5e %5d %5d %.5e %.5e ",t,NX+2,NY,BHSPIN,MASS);
   else
     fprintf(fout1,"%.5e %5d %5d %5d %.5e %.5e ",t,NX+2,NY,NZ,BHSPIN,MASS);
#if(MYCOORDS==MKS3COORDS)
   fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e\n",MKSR0,MKSH0,MKSMY1,MKSMY2,MKSMP0);
#endif
#if(MYCOORDS==MKS2COORDS)
   fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e\n",MKSR0,MKSH0,-1.,-1.,-1.);
#endif
#ifdef RELELECTRONS
   fprintf(fout1,"%5d %.5e %.5e\n",NRELBIN, RELGAMMAMIN, RELGAMMAMAX);
#endif
#endif //GRTRANSSIMOUTPUT
   
#ifdef RELELECTRONS //ANDREW array for finding nonthermal gamma break
  int ie;
  ldouble gammapbrk[NRELBIN];
  for(ie=0; ie<NRELBIN; ie++) gammapbrk[ie] = pow(relel_gammas[ie], RELEL_HEAT_INDEX + 0.5);
#endif  

   // loop over all cells  
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

#if defined(GRTRANSSIMOUTPUT) || defined(GRTRANSSIMOUTPUT_2)
	       ldouble x1=geom.xx;
	       ldouble x2=geom.yy;
	       ldouble x3=geom.zz;
	       fprintf(fout1,"%d %d %d ",ix,iy,iz); //(1-3)
	       fprintf(fout1,"%.5e %.5e %.5e ",x1,x2,x3); //(4-6)
	       fprintf(fout1,"%.5e %.5e %.5e ",r,th,ph); //(7-9)
	       
               //ANDREW in grtrans, fill values below horizon with values right above horizon (??)
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

	      //cell dimensions
     	      //ANDREW put cell size code in a function with precompute option
              ldouble dxph[3],dx[3];
              get_cellsize_out(ix, iy, iz, dx);

	      /*
	      ldouble xx1[4],xx2[4];
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
              */
	      dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
	      dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
	      dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);

	      ldouble gdet=geom.gdet;
	      ldouble volume=gdet*get_size_x(ix,0)*get_size_x(iy,1)*get_size_x(iz,2);
	       
               
#if defined(GRTRANSSIMOUTPUT) || defined(GRTRANSSIMOUTPUT_2)
	       volume=gdet; // For grtrans output, replace volume with gdet
#endif
	       ldouble rho,rhoucont,uint,pgas,temp,bsq,bcon[4],bcov[4];
	       ldouble utcon[4],utcov[4],ucon[4],ucov[4],Tij[4][4],Tij22[4][4];
	       ldouble Ti,Te;
	       ldouble gamma=GAMMA;

	       int i,j;

#ifdef CONSISTENTGAMMA
	       gamma=pick_gammagas(ix,iy,iz);
#endif
	     
	       if(doingavg)
		 {
                   //ANDREW we need pp for some relel computations below
                   #ifdef PRECOMPUTE_MY2OUT
                   trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
                   #else      
                   trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomBL);
                   #endif
		   
		   rhoucont=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz);

		   rho=get_uavg(pavg,RHO,ix,iy,iz);
		   uint=get_uavg(pavg,UU,ix,iy,iz);
		   pgas=get_uavg(pavg,AVGPGAS,ix,iy,iz);
		   temp=calc_PEQ_Tfromprho(pgas,rho,ix,iy,iz);

		   vel[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   vel[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   vel[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   vel[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
                   for(i=0;i<4;i++) vcon[i]=vel[i];
                   lorentz = fabs(vel[0])/sqrt(fabs(geomBL.GG[0][0]));

		   utcon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);

                   int ii,jj;
                   for(ii=0;ii<4;ii++)
		     for(jj=0;jj<4;jj++)
		       Tij[ii][jj]=get_uavg(pavg,AVGTIJ(ii,jj),ix,iy,iz);                 

		   indices_2122(Tij,Tij22,geomBL.gg);  
                 
#if defined(GRTRANSSIMOUTPUT) || defined(GRTRANSSIMOUTPUT_2)
                   //ANDREW NORMALIZE u^0 for grtrans
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

#if defined(GRTRANSSIMOUTPUT) || defined(GRTRANSSIMOUTPUT_2)
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
		  ldouble rhoeth=ne*MU_E*M_PROTON;
		  pp[ENTRE]=calc_SefromrhoT(rhoeth,Te,ELECTRONS);
		  pp[ENTRI]=calc_SefromrhoT(rho,Ti,IONS);
		 
#endif                
		 }
	       else //not doingavg; on the go from the primitives
		 { 
                   #ifdef PRECOMPUTE_MY2OUT
                   trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
                   #else      
                   trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomBL);
                   #endif

		   rho=pp[0];
		   uint=pp[1];
		   pgas=(gamma-1.)*uint;
                   ldouble vel[4],vcon[4],tauscarloc;
                   
                   //obtain 4 velocity
	           vel[1]=pp[VX];
	           vel[2]=pp[VY];
	           vel[3]=pp[VZ];
	           conv_vels(vel,vel,VELPRIM,VEL4,geomBL.gg,geomBL.GG);
                   for(i=0;i<4;i++) vcon[i]=vel[i];
                   lorentz = fabs(vel[0])/sqrt(fabs(geomBL.GG[0][0]));
		   rhoucont=pp[RHO]*vcon[0];

                   calc_Tij(pp,&geomBL,Tij22);
		   indices_2221(Tij22,Tij,geomBL.gg);
                   calc_ucon_ucov_from_prims(pp, &geomBL, utcon, ucov);
                   temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
                   temp=calc_PEQ_Teifrompp(pp,&Te,&Ti,geomBL.ix,geomBL.iy,geomBL.iz);
#ifdef MAGNFIELD
                   calc_bcon_bcov_bsq_from_4vel(pp, utcon, ucov, &geomBL, bcon, bcov, &bsq);
#endif
		 }

#ifdef OUTPUTINGU
	       rho=pp[RHO];
#else
	       rho=rhoGU2CGS(pp[RHO]);
#endif
	     
#ifdef RHOLABINSIM
	       rho*=utcon[0];
#endif

	       //!! replaced temp with pgas
	       //fprintf(fout1,"%.5e %.5e ",rho,pgas); //(7-8) , (10-11) for grtrans 

	       fprintf(fout1,"%.5e %.5e ",rho,temp); //(7-8) , (10-11) for grtrans 
	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",utcon[0],utcon[1],utcon[2],utcon[3]); //(9-12) , (12-15) for grtrans
#ifndef GRTRANSSIMOUTPUT_2
	       fprintf(fout1,"%.5e ", volume); // (13) , (16) for grtrans
#endif
	       ldouble ehat=0.;
#ifdef RADIATION
	       ldouble Rtt,Rij[4][4],Rij22[4][4],vel[4];
	       ldouble ugas[4],Fx,Fy,Fz;
	       ldouble Gi[4],Giff[4]={0.,0.,0.,0.};
	       ldouble Gic[4],Gicff[4]={0.,0.,0.,0.};
	       ldouble Trad;

	       ldouble CoulombCoupling=0.;
	       if(doingavg==0) // from primitives
		{
		  calc_ff_Rtt(pp,&Rtt,ugas,&geomBL);
		  ehat=-Rtt;  
		  calc_Rij(pp,&geomBL,Rij22); //calculates R^munu in OUTCOORDS
		  indices_2221(Rij22,Rij,geomBL.gg);

		  //four fource
		  calc_Gi(pp,&geomBL,Giff,0.0,0,0); //ANDREW 0 for fluid frame
		  
#if defined(COMPTONIZATION)
		  ldouble kappaes=calc_kappaes(pp,&geomBL);
		  //directly in ff
		  vel[1]=vel[2]=vel[3]=0.; vel[0]=1.;
		  calc_Compt_Gi(pp,&geomBL,Gicff,ehat,Te,kappaes,vel);
#endif 
                  Trad = calc_LTE_TfromE(ehat);
#ifdef EVOLVEPHOTONNUMBER //the color temperature of radiation
                  Trad = calc_ncompt_Thatrad_full(pp,&geomBL);
#endif
	
#ifdef EVOLVEELECTRONS	 
#ifndef  SKIPCOULOMBCOUPLING
  CoulombCoupling=calc_CoulombCoupling(pp,&geomBL); 
#endif
#endif
		}
	       else //from avg
		{
		  ehat=get_uavg(pavg,AVGEHAT,ix,iy,iz);
		  int i,j;
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      Rij[i][j]=get_uavg(pavg,AVGRIJ(i,j),ix,iy,iz); 
		  for(j=0;j<4;j++)
		    Giff[j]=get_uavg(pavg,AVGGHAT(j),ix,iy,iz);
		    
		  indices_2122(Rij,Rij22,geomBL.gg);  

#if defined(COMPTONIZATION)
		  for(j=0;j<4;j++)
		    Gicff[j]=get_uavg(pavg,AVGGHATCOMPT(j),ix,iy,iz);
#endif

		  Trad = calc_LTE_TfromE(ehat);
#ifdef EVOLVEPHOTONNUMBER
		  Trad=calc_ncompt_Thatrad_fromEN(ehat,get_uavg(pavg,AVGNFHAT,ix,iy,iz));
#endif
		}
	       
	       //flux
	       Fx=Rij[1][0];
	       Fy=Rij[2][0];
	       Fz=Rij[3][0];
#ifndef GRTRANSSIMOUTPUT_2
#ifdef OUTPUTINGU
	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",ehat,Fx,Fy,Fz); //(14) - (17) 
	       fprintf(fout1,"%.5e %.5e ",Giff[0],Gicff[0]); //(18)-(19)
               fprintf(fout1,"%.5e %.5e ",ehat, Trad); //(20), (21)         
#else
	       fprintf(fout1,"%.5e %.5e %.5e %.5e ",endenGU2CGS(ehat),fluxGU2CGS(Fx),fluxGU2CGS(Fy),fluxGU2CGS(Fz)); //(14) - (17), (17)-(20) for grtrans
	       ldouble conv=kappaGU2CGS(1.)*rhoGU2CGS(1.)*endenGU2CGS(1.)*CCC; //because (cE-4piB) in non-geom
	       fprintf(fout1,"%.5e %.5e ",Giff[0]*conv,Gicff[0]*conv); //(18)-(19) , (21)-(22) for grtrans
               fprintf(fout1,"%.5e %.5e ",CoulombCoupling*conv, Trad); //(20), (21) , (23)-(24) for grtrans      
#endif
#endif
#endif //RADIATION

	       ldouble gammam1=gamma-1.;
	       ldouble betarad=ehat/3./(pgas);
#ifdef RADIATION
	       //ldouble muBe = (-Tij[1][0]-Rij[1][0] - rhouconr)/rhouconr;
	       ldouble bernoulli=(-Tij[0][0] -Rij[0][0] - rhoucont)/rhoucont;
#else
	       //ldouble muBe = (-Tij[1][0] - rhouconr)/rhouconr;
	       ldouble bernoulli=(-Tij[0][0] - rhoucont)/rhoucont;
#endif
	       
	       //magn. field components
#ifdef MAGNFIELD
	       if(doingavg==0) 
	        {
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

                  #if defined(GRTRANSSIMOUTPUT) || defined(GRTRANSSIMOUTPUT_2)
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
	       
	       //to CGS!
	       #ifndef OUTPUTINGU
	       bsq=endenGU2CGS(bsq);
	       ldouble scaling=endenGU2CGS(1.);
	       for(i=0;i<4;i++)
		 {
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
	      
	       if(iy>=iphimin && iy<=iphimax)
		 {
                   #ifndef RAD_INTEGRATION
		   phi+=geom.gdet*get_u(p,B1,ix,iy,iz)*get_size_x(iy,1)*2.*M_PI;
		   tausca+=rho*0.34*lenGU2CGS(dxph[1]); //rho already converted to cgs
                   #endif
		 }
               #ifdef RAD_INTEGRATION
               #ifdef RADIATION
	       if(ix>=0)
		 {
		   //tauscarloc = vcon[0]*(1.-abs(vcon[1]))*calc_kappaes(pp,&geomBL); //rho already converted to cgs
		   tauscarloc = (vcon[0]-vcon[1])*calc_kappaes(pp,&geomBL); //rho already converted to cgs

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

	       
	       //(14) - (19) or (22) - (27) if radiation included (+3 with grtrans 3d)
#if defined(GRTRANSSIMOUTPUT)
	       fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e %.5e ",bcon[0],bcon[1],bcon[2],bcon[3],bsq,phi);
#elif defined(GRTRANSSIMOUTPUT_2)
	       fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e ", bcon[0],bcon[1],bcon[2],bcon[3],bsq); 
#else    
	       fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e %.5e ",bsq,bcon[1],bcon[2],bcon[3],phi,betamag);
#endif //GRTRANSSIMOUTPUT

	       // (28) - (29) when rad and magn field on, (20) - (21) with no radiation (+3 with grtrans 3d)
#ifndef GRTRANSSIMOUTPUT_2
#ifndef RAD_INTEGRATION
               fprintf(fout1,"%.5e %.5e ",betarad,bernoulli); 	       
	       //fprintf(fout1,"%.5e %.5e ",betarad,tausca); 
	       //fprintf(fout1,"%.5e %.5e ",betarad,muBe); 
#else
	       fprintf(fout1,"%.5e %.5e ",betarad,bernoulli); 	       
	       //fprintf(fout1,"%.5e %.5e ",betarad,tauscar); 
	       //fprintf(fout1,"%.5e %.5e ",betarad,muBe); 
#endif //RAD_INTEGRATION
#endif //GRTRANSSIMOUTPUT_2
#endif //MAGNFIELD
	       
	       

#ifdef EVOLVEELECTRONS

	       // (30) - (32) when rad and magn field on, (22) - (24) with no radiation (+3 with grtrans 3d)
	       fprintf(fout1,"%.5e %.5e %.5e ",Te, Ti, gamma); 
	       
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
		   ldouble rhoeth=calc_thermal_ne(pp)*M_PROTON*MU_E;
		   ue=calc_ufromSerho(pp[ENTRE],rhoeth,ELECTRONS,ix,iy,iz); 
		 }


	       //ANDREW
	       //in avg, vischeat was averaged as du, not du/dtau
	       //recompute dt and use that as an estimate
               #ifdef DIVIDEVISCHEATBYDT
               dt=get_u_scalar(cell_dt,ix,iy,iz); //individual time step
	       ldouble dtau = dt/vel[0];
	       vischeat/=dtau;
               #endif
               
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
/* prints phi-averaged ASCII */
/* only code-comparison quantities for now */
/*********************************************/
                              
int fprint_simple_phiavg(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%s%04d_simphiavg.dat",folder,prefix,nfile);
   fout1=fopen(bufor,"w");
   printf("fopen %s in function fprint_simple_phiavg\n", bufor);

   /***********************************/  
   /** writing order is fixed  ********/  
   /***********************************/  
 
   int ix,iy,iz,iv;
   int iix;
   ldouble pp[NV],lorentz;
   int nz=NZ;
   struct geometry geom,geomBL;

   int xmin=0;

   //HEADER
   int cgsout=0;
   #ifdef CGSOUTPUT
   cgsout=1;
   #endif
   fprintf(fout1,"t=%.5e M=%.5e a=%.5e MYCOORDS=%d OUTCOORDS=%d\n CGSOUT=%d\n",t,BHSPIN,MASS,MYCOORDS,OUTCOORDS,cgsout);
   fprintf(fout1,"--------------------------------------------------------------------\n");
   
   // loop over all cells  
   for(iix=xmin;iix<NX;iix++)
   {
     ix=iix;
     for(iy=0;iy<NY;iy++)
     {

       ldouble rhoavg=0.;
       ldouble uintavg=0.;
       ldouble pgasavg=0.;
       ldouble bsqavg=0.;
       ldouble betaavg=0.;
       ldouble sigmaavg=0.;
       ldouble betainvavg=0.;
       ldouble gdetavg=0.;
       ldouble uconavg[4], Bavg[3];
       uconavg[0]=0.; uconavg[1]=0.; uconavg[2]=0.;uconavg[3]=0.;
       Bavg[0]=0.;Bavg[1]=0.;Bavg[2]=0.;

       ldouble r,th,ph;
       ldouble x1,x2,x3;
       for(iz=0;iz<nz;iz++)
       {
           
	      fill_geometry(ix,iy,iz,&geom);
	      fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS);

	      r=geomBL.xx; th=geomBL.yy; ph=geomBL.zz;
	      x1=geom.xx; x2=geom.yy;x3=geom.zz;

	      // fill in pp
	      for(iv=0;iv<NV;iv++)
	      {
	        if(doingavg)
	          pp[iv]=get_uavg(pavg,iv,ix,iy,iz);
	        else
	          pp[iv]=get_u(p,iv,ix,iy,iz);
	      }

	      //cell dimensions
      	      //ANDREW put cell size code in a function with precompute option
              ldouble dxph[3],dx[3];
	      get_cellsize_out(ix, iy, iz, dx);
	      /*
	      ldouble xx1[4],xx2[4];
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
              */
		
	      dxph[0]=dx[0]*sqrt(geomBL.gg[1][1]);
	      dxph[1]=dx[1]*sqrt(geomBL.gg[2][2]);
	      dxph[2]=dx[2]*sqrt(geomBL.gg[3][3]);

	      ldouble gdet=geom.gdet;
	      ldouble volume=gdet*get_size_x(ix,0)*get_size_x(iy,1)*get_size_x(iz,2);
	       
              ldouble rho,uint,pgas,temp,bsq,bcon[4],bcov[4],Bprim[3];
	      ldouble utcon[4],utcov[4],ucon[4],ucov[4],Tij[4][4],Tij22[4][4];
	      ldouble Ti,Te;
	      ldouble gamma=GAMMA; 
	      int i,j;

              #ifdef CONSISTENTGAMMA
	      gamma=pick_gammagas(ix,iy,iz);
              #endif
	     
	      if(doingavg)
		 {
                   #ifdef PRECOMPUTE_MY2OUT
                   trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
                   #else      
                   trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomBL);
                   #endif
		   
		   rho=get_uavg(pavg,RHO,ix,iy,iz);
		   uint=get_uavg(pavg,UU,ix,iy,iz);
		   pgas=get_uavg(pavg,AVGPGAS,ix,iy,iz);
		   temp=calc_PEQ_Tfromprho(pgas,rho,ix,iy,iz);

		   utcon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
		   utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);

		   lorentz = fabs(utcon[0])/sqrt(fabs(geomBL.GG[0][0]));

                   //ANDREW NORMALIZE u^0 for grtrans
                   fill_utinucon(utcon,geomBL.gg,geomBL.GG);
		   indices_21(utcon,utcov,geomBL.gg); 

		   /*
                   int ii,jj;
                   for(ii=0;ii<4;ii++)
		     for(jj=0;jj<4;jj++)
		       Tij[ii][jj]=get_uavg(pavg,AVGTIJ(ii,jj),ix,iy,iz);                 
		   indices_2122(Tij,Tij22,geomBL.gg);  
                   */

                   pp[RHO]=rho;
		   pp[UU]=uint;
#ifdef MAGNFIELD
		   bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
		   bcon[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
		   bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
		   bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
		   bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);

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
	       else //not doingavg; on the go from the primitives
	       { 
                   #ifdef PRECOMPUTE_MY2OUT
                   trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
                   #else      
                   trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS, geom.xxvec,&geom,&geomBL);
                   #endif

		   rho=pp[0];
		   uint=pp[1];
		   pgas=(gamma-1.)*uint;

		   /*
                   calc_Tij(pp,&geomBL,Tij22);
		   indices_2221(Tij22,Tij,geomBL.gg);
                   */
		   
                   calc_ucon_ucov_from_prims(pp, &geomBL, utcon, ucov);                   
                   lorentz = fabs(utcon[0])/sqrt(fabs(geomBL.GG[0][0]));
                   temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
                   //temp=calc_PEQ_Teifrompp(pp,&Te,&Ti,geomBL.ix,geomBL.iy,geomBL.iz);
#ifdef MAGNFIELD
                   calc_bcon_bcov_bsq_from_4vel(pp, utcon, ucov, &geomBL, bcon, bcov, &bsq);
#endif
	       }
  	     
	       ldouble betamag = -1;
	       ldouble betamaginv = 0.;
	       ldouble sigmamag = 0.;
	       Bprim[0]=0.;Bprim[1]=0.;Bprim[2]=0.;
#ifdef MAGNFIELD
	       betamag = 2.*(pgas)/bsq;
	       betamaginv = 1./betamag;
	       sigmamag = bsq/rho;
	       Bprim[0] = pp[B1];
	       Bprim[1] = pp[B2];
	       Bprim[2] = pp[B3];
#endif
	       //average
	       gdetavg += gdet*dx[2];
	       rhoavg += rho*gdet*dx[2];      
	       uintavg += uint*gdet*dx[2];
	       pgasavg += pgas*gdet*dx[2];
	       uconavg[0] += utcon[0]*gdet*dx[2];
	       uconavg[1] += utcon[1]*gdet*dx[2];      
	       uconavg[2] += utcon[2]*gdet*dx[2];      
	       uconavg[3] += utcon[3]*gdet*dx[2];      
	       Bavg[0] += Bprim[0]*gdet*dx[2];
	       Bavg[1] += Bprim[1]*gdet*dx[2];
	       Bavg[2] += Bprim[2]*gdet*dx[2];
	       bsqavg += bsq*gdet*dx[2];
	       betaavg += betamag*gdet*dx[2];
	       betainvavg += betamaginv*gdet*dx[2];      	       
	       sigmaavg += sigmamag*gdet*dx[2];      
	 
       } // end loop over z

       rhoavg /= gdetavg;
       uintavg /= gdetavg;
       pgasavg /= gdetavg;
       bsqavg /= gdetavg;
       betaavg /= gdetavg;
       betainvavg /= gdetavg;
       sigmaavg /= gdetavg;
       uconavg[0] /= gdetavg;
       uconavg[1] /= gdetavg;
       uconavg[2] /= gdetavg;
       uconavg[3] /= gdetavg;
       Bavg[0] /= gdetavg;
       Bavg[1] /= gdetavg;
       Bavg[2] /= gdetavg;

       //convert to cgs
       #ifdef CGSOUTPUT
       rhoavg = rhoGU2CGS(rhoavg);
       uintavg = endenGU2CGS(uintavg);
       pgasavg = endenGU2CGS(pgasavg);
       #endif
       
       // print profile
       fprintf(fout1,"%d %d ",ix,iy); //(1-2)
       fprintf(fout1,"%.5e %.5e ",x1,x2); //(3-4)
       fprintf(fout1,"%.5e %.5e ",r,th); //(5-6)
       fprintf(fout1,"%.5e %.5e %.5e ",rhoavg,uintavg,pgasavg); //(7-9)
       fprintf(fout1,"%.5e %.5e %.5e %.5e ",uconavg[0],uconavg[1],uconavg[2],uconavg[3]); //(10-13)
       fprintf(fout1,"%.5e %.5e %.5e ",Bavg[0],Bavg[1],Bavg[2]); //(14-16)
       fprintf(fout1,"%.5e %.5e %.5e %.5e ",bsqavg,betaavg,betainvavg,sigmaavg); //(17-20)
       fprintf(fout1,"\n");

     } // iy loop
   } // ix loop

   fflush(fout1);
   fclose(fout1);

   return 0;
}

/*********************************************/
/* print radius and phi dependent correlation functions in rho and betainv ASCII */
/* see code comparison paper */
/*********************************************/
                              
int fprint_simple_phicorr(ldouble t, int nfile, char* folder,char* prefix)
 {
   char bufor[50];
   sprintf(bufor,"%s/%s%04d_phicorr.dat",folder,prefix,nfile);
   fout1=fopen(bufor,"w");
   printf("fopen %s in function fprint_simple_phicorr\n", bufor);

   /***********************************/  
   /** writing order is fixed  ********/  
   /***********************************/  
 
   int ix,iy,iz,izz,iv;
   ldouble pp[NV];
   int nz=NZ;
   struct geometry geom,geomBL;

   int xmin=0;
   
   int iymin=NY/2 -1;
   int iymax=NY/2;

   // loop over all cells  
   for(ix=xmin;ix<NX;ix++)
   {

     ldouble r,th,ph;
     ldouble x1,x2,x3;

     //loop over phi values for output print
     for(izz=0;izz<nz;izz++)
     {

       ldouble rhocorr=0.;
       ldouble betainvcorr=0.;
       ldouble rhoavg=0.;
       ldouble betainvavg=0.;
       ldouble volavg=0.;
       
       for(iy=iymin;iy<(iymax+1);iy++)
       {
         ldouble rho,uint,pgas,bsq,betamaginv;
         ldouble rhotmp[nz], betainvtmp[nz], dxtmp[nz];

	 fill_geometry(ix,iy,izz,&geom);
	 fill_geometry_arb(ix,iy,izz,&geomBL,OUTCOORDS);
	 r=geomBL.xx; th=geomBL.yy; ph=geomBL.zz;
	 x1=geom.xx; x2=geom.yy;x3=geom.zz;

         // first loop and make a temporary array to store necessary data
         for(iz=0;iz<nz;iz++)
         {
	 
	   ldouble dx[3];
	   //cell dimensions
    	   //ANDREW put cell size code in a function with precompute option
           get_cellsize_out(ix, iy, iz, dx);
	   /*
	      ldouble x1[4],xx2[4];
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
	   */
	    if(doingavg)
	    {
	     rho=get_uavg(pavg,RHO,ix,iy,iz);
	     pgas=get_uavg(pavg,AVGPGAS,ix,iy,iz);
	     bsq=0.;
	     #ifdef MAGNFIELD
	     bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
	     #endif
	    }
            else
	    {
	     for(iv=0;iv<NV;iv++)
	       pp[iv]=get_u(p,iv,ix,iy,iz);
	     rho=pp[0];
	     uint=pp[1];
	     ldouble gamma=GAMMA; 
             #ifdef CONSISTENTGAMMA
	     gamma=pick_gammagas(ix,iy,iz);
             #endif
	     pgas=(gamma-1.)*uint;
	     bsq=0.;
	     ldouble utcon[4],utcov[4],bcon[4],bcov[4];
	     #ifdef MAGNFIELD
             calc_bcon_bcov_bsq_from_4vel(pp, utcon, utcov, &geomBL, bcon, bcov, &bsq);
	     #endif
	    }
	    betamaginv = bsq/(2*pgas);
	  
	    //fill correlations temporary arrays
	    betainvtmp[iz] = betamaginv;
	    rhotmp[iz] = rho;
	    dxtmp[iz] = dx[2]*dx[1];
         } //end first loop filling temporary arrays

	 //now integrate over z
         for(iz=0;iz<nz;iz++)
	 {
	   int iz_loop = iz+izz;
	   if(iz_loop>=NZ)
	     iz_loop -= NZ;

	   rhoavg += rhotmp[iz]*dxtmp[iz];
	   betainvavg += betainvtmp[iz]*dxtmp[iz];
	   volavg += dxtmp[iz];

	   rhocorr += rhotmp[iz]*rhotmp[iz_loop]*dxtmp[iz];
	   betainvcorr += betainvtmp[iz]*betainvtmp[iz_loop]*dxtmp[iz];
	 } // end integrating loop over z
       } // end integrating loop over y

       rhoavg /= volavg;
       betainvavg /= volavg;
       rhocorr -= volavg*rhoavg*rhoavg;
       betainvcorr -= volavg*betainvavg*betainvavg;
       
       // print profile
       fprintf(fout1,"%d %d ",ix,iz); //(1-2)
       fprintf(fout1,"%.5e %.5e ",r,ph); //(3-4)
       fprintf(fout1,"%.5e %.5e ",rhoavg,betainvavg); //(5-6)
       fprintf(fout1,"%.5e %.5e ",rhocorr,betainvcorr); //(7-8)
       fprintf(fout1,"\n");

       } //iz loop print
   } // ix loop

   fflush(fout1);
   fclose(fout1);

   return 0;
}

/////////////////////////////////////////////////       
void get_prim_name
  (
   char* prim_name,
  int iv
   )
{
  if (iv == RHO)
    {
      sprintf(prim_name,"/RHO");
      return;
    }
  else if (iv == UU)
    {
      sprintf(prim_name,"/UU");
      return;
    }
  else if (iv == VX)
    {
      sprintf(prim_name,"/VX");
      return;
    }
  else if (iv == VY)
    {
      sprintf(prim_name,"/VY");
      return;
    }
  else if (iv == VZ)
    {
      sprintf(prim_name,"/VZ");
      return;
    }
  else if (iv == ENTR)
    {
      sprintf(prim_name,"/ENTR");
      return;
    }
  
#ifdef MAGNFIELD
  if (iv == B1)
    {
      sprintf(prim_name,"/B1");
      return;
    }
  else if (iv == B2)
    {
      sprintf(prim_name,"/B2");
      return;
    }
  else if (iv == B3)
    {
      sprintf(prim_name,"/B3");
      return;
    }
#endif
  
#ifdef EVOLVEELECTRONS
  if (iv == ENTRE)
    {
      sprintf(prim_name,"/ENTRE");
      return;
    }
  else if (iv == ENTRI)
    {
      sprintf(prim_name,"/ENTRI");
      return;
    }
#endif
  
#ifdef RADIATION
  if (iv == EE)
    {
      sprintf(prim_name,"/EE");
      return;
    }
  else if (iv == FX)
    {
      sprintf(prim_name,"/FX");
      return;
    }
  else if (iv == FY)
    {
      sprintf(prim_name,"/FY");
      return;
    }
  else if (iv == FZ)
    {
      sprintf(prim_name,"/FZ");
      return;
    }
  
#ifdef EVOLVEPHOTONNUMBER
  if (iv == NF)
    {
      sprintf(prim_name,"/NF");
      return;
    }
#endif
#endif

#ifdef RELELECTRONS
  if (iv >= 8 && iv < NVHD)
    {
      int relbin = iv - 8;
      sprintf(prim_name,"/RELBIN%d", relbin);
      return;
    }
#endif
  
  sprintf(prim_name,"/UNKNOWN");
  return;
}

//serial hdf5 output for postprocessing
//derived from illinois standard
//https://github.com/AFD-Illinois/docs/wiki/GRMHD-Output-Format

int 
fprint_anaout_hdf5(ldouble t, char* folder, char* prefix)
{
#ifdef ANAOUT_HDF5


  // Write out header information in group HEADER in HDF5 file

  hid_t dumps_file_id, dumps_group_id, dumps_group_id2, dumps_group_id3,dumps_dataspace_scalar, dumps_dataspace_str, dumps_dataspace_array;
  hid_t dumps_dataset_int, dumps_dataset_double, dumps_dataset_str, dumps_dataset_array, dumps_attribute_id,strtype;
  hsize_t dims_h5[3];
  herr_t status;
    
  dims_h5[0] = NX;
  dims_h5[1] = NY;
  dims_h5[2] = NZ;

  // Make HDF5 file
  char fname_h5[256];
  //sprintf(fname_h5, "%s/ana%04d.h5", folder, nfout1);
  sprintf(fname_h5,"%s/%s%04d.h5",folder,prefix,nfout1);
  dumps_file_id = H5Fcreate(fname_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Dataspaces
  dumps_dataspace_scalar = H5Screate(H5S_SCALAR);
  dumps_dataspace_array = H5Screate_simple(3, dims_h5, NULL);
 
  // Dump time at top level (following illinois convention)
  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/t", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &t);
  status = H5Dclose(dumps_dataset_double);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Make Header
  int file_number = nfout1, problem_number = PROBLEM, nxx = TNX, nyy = TNY, nzz = TNZ, nprimitives = NV;
  int has_radiation=0, has_electrons=0;
  int ndim;
  char version[40];
  char metric_run[40];
  char metric_out[40];
  ldouble gam=GAMMA;
  ldouble bhspin=BHSPIN;
  
  if(nxx>1 && nyy>1 && nzz>1)
    ndim=3;
  else if((nxx>1 && nyy>1) || (nyy>1 && nzz>1) || (nxx>1 && nzz>1))
    ndim=2;
  else
    ndim=3;

  if(MYCOORDS==JETCOORDS)
    sprintf(metric_run,"%s","JETCOORDS");
  else if(MYCOORDS==MKS3COORDS)
    sprintf(metric_run,"%s","MKS3");
  else if(MYCOORDS==MKS2COORDS)
    sprintf(metric_run,"%s","MKS2");
  else if(MYCOORDS==KSCOORDS)
    sprintf(metric_run,"%s","KS");
  else
  {
    printf("COORDS can only be JETCOORDS,MKS3,MKS2,KS for HDF5 output!\n");
    exit(-1);
  }

#ifndef OUTCOORDS2
#define OUTCOORDS2 OUTCOORDS
#endif
  
  if(OUTCOORDS2==JETCOORDS)
    sprintf(metric_out,"%s","JETCOORDS");
  else if(OUTCOORDS2==MKS3COORDS)
    sprintf(metric_out,"%s","MKS3");
  else if(OUTCOORDS2==MKS2COORDS)
    sprintf(metric_out,"%s","MKS2");
  else if(OUTCOORDS2==KSCOORDS)
    sprintf(metric_out,"%s","KS");
  else if(OUTCOORDS2==BLCOORDS)
    sprintf(metric_out,"%s","BL");
  else
    {
      printf("OUTCOORDS2 can only be JETCOORDS,MKS3,MKS2,KS,BL for HDF5 output!\n");
      exit(-1);
    }
  
#ifdef RADIATION
  has_radiation=1;
#endif
#ifdef EVOLVEELECTRONS
  has_electrons=1;
#endif
  
  dumps_group_id = H5Gcreate2(dumps_file_id, "/header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  #ifdef ANAOUT_HDF5_V1
  sprintf(version,"%s","KORALv2");
  #else
  sprintf(version,"%s","KORALv2");
  #endif
  strtype=H5Tcopy(H5T_C_S1);
  status=H5Tset_size(strtype,strlen(metric_run));		     
  dumps_dataset_str = H5Dcreate2(dumps_file_id, "/header/version", strtype, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_str, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &version); 
  status = H5Dclose(dumps_dataset_str);

  strtype=H5Tcopy(H5T_C_S1);
  status=H5Tset_size(strtype,strlen(metric_run));		     
  dumps_dataset_str = H5Dcreate2(dumps_file_id, "/header/metric_run", strtype, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_str, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &metric_run); 
  status = H5Dclose(dumps_dataset_str);

  strtype=H5Tcopy(H5T_C_S1);
  status=H5Tset_size(strtype,strlen(metric_out));		     
  dumps_dataset_str = H5Dcreate2(dumps_file_id, "/header/metric_out", strtype, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); 
  status = H5Dwrite(dumps_dataset_str, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &metric_out); 
  status = H5Dclose(dumps_dataset_str);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/ndim", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &ndim);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/n1", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &nxx);
  status = H5Dclose(dumps_dataset_int);
  
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/n2", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &nyy);
  status = H5Dclose(dumps_dataset_int);
  
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/n3", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &nzz);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/n_prim", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    &nprimitives);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/gam", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &gam);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/bhspin", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &bhspin);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/file_number", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &file_number);
  status = H5Dclose(dumps_dataset_int);
 
  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/problem_number", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &problem_number);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/has_radiation", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &has_radiation);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/has_electrons", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &has_electrons);
  status = H5Dclose(dumps_dataset_int);
  

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////  
  // Units Data
  dumps_group_id2 = H5Gcreate2(dumps_file_id, "/header/units", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ldouble M_bh=MASS;
  ldouble L_unit=lenGU2CGS(1.);
  ldouble T_unit=timeGU2CGS(1.);
  ldouble M_unit=rhoGU2CGS(1.);               
  ldouble U_unit=endenGU2CGS(1.);
  ldouble B_unit=sqrt(4*M_PI*endenGU2CGS(1.));

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/units/M_bh", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &M_bh);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/units/L_unit", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &L_unit);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/units/T_unit", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &T_unit);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/units/M_unit", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &M_unit);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/units/U_unit", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &U_unit);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/units/B_unit", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &B_unit);
  status = H5Dclose(dumps_dataset_double);

  status = H5Gclose(dumps_group_id2); //close /header/units

  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////  
  // Geometry Data
  dumps_group_id2 = H5Gcreate2(dumps_file_id, "/header/geom", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ldouble startx1=MINX;
  ldouble startx2=MINY;
  ldouble startx3=MINZ;
  ldouble dx1=(MAXX-MINX)/(ldouble)TNX;
  ldouble dx2=(MAXY-MINY)/(ldouble)TNY;
  ldouble dx3=(MAXZ-MINZ)/(ldouble)TNZ;
  
  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/startx1", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &startx1);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/startx2", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &startx2);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/startx3", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &startx3);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/dx1", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &dx1);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/dx2", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &dx2);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/dx3", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &dx3);
  status = H5Dclose(dumps_dataset_double);

  
  // METRIC data
#if(MYCOORDS==JETCOORDS)
  dumps_group_id3 = H5Gcreate2(dumps_file_id, "/header/geom/jetcoords", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ldouble mksr0=MKSR0;
  ldouble rmin=RMIN;
  ldouble rmax=RMAX;
  ldouble rbrk=HYPRBRK;
  ldouble fjet=FJET;
  ldouble fdisk=FDISK;
  ldouble runi=RUNI;
  ldouble rcoll_jet=RCOLL_JET;
  ldouble rdecoll_jet=RDECOLL_JET;
  ldouble rcoll_disk=RCOLL_DISK;
  ldouble rdecoll_disk=RDECOLL_DISK;
  ldouble alpha1=ALPHA_1;
  ldouble alpha2=ALPHA_2;
#ifdef CYLINDRIFY
  int cylindrify=1;
  ldouble rcyl=RCYL;
  ldouble ncyl=NCYL;
#else
  int cylindrify=0;
  ldouble rcyl=-1.;
  ldouble ncyl=-1.;
#endif
  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/bhspin", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &bhspin);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/mksr0", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &mksr0);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/rmin", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &rmin);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/rmax", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &rmax);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/rbrk", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &rbrk);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/fjet", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &fjet);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/fdisk", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &fdisk);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/runi", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &runi);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/rcoll_jet", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &rcoll_jet);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/rdecoll_jet", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &rdecoll_jet);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/rcoll_disk", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &rcoll_disk);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/rdecoll_disk", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &rdecoll_disk);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/alpha1", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &alpha1);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/alpha2", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &alpha2);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_int = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/cylindrify", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &cylindrify);
  status = H5Dclose(dumps_dataset_int);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/rcyl", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &rcyl);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/jetcoords/ncyl", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &ncyl);
  status = H5Dclose(dumps_dataset_double);
  
  status = H5Gclose(dumps_group_id3);
#endif //JETCOORDS
  
#if(MYCOORDS==MKS3COORDS)
  dumps_group_id3 = H5Gcreate2(dumps_file_id, "/header/geom/mks3", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ldouble mksr0=MKSR0;
  ldouble mksh0=MKSH0;
  ldouble mksmy1=MKSMY1;
  ldouble mksmy2=MKSMY2;
  ldouble mksmp0=MKSMP0;

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/mks3/bhspin", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &bhspin);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/mks3/mksr0", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &mksr0);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/mks3/mksh0", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &mksh0);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/mks3/mksmy1", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &mksmy1);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/mks3/mksmy2", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &mksmy2);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/mks3/mksmp0", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &mksmp0);
  status = H5Dclose(dumps_dataset_double);

  status = H5Gclose(dumps_group_id3);
#endif //MKS3COORDS

#if(MYCOORDS==MKS2COORDS)
  dumps_group_id3 = H5Gcreate2(dumps_file_id, "/header/geom/mks2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  ldouble mksr0=MKSR0;
  ldouble mksh0=MKSH0;

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/mks2/bhspin", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &bhspin);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/mks2/mksr0", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &mksr0);
  status = H5Dclose(dumps_dataset_double);

  dumps_dataset_double = H5Dcreate2(dumps_file_id, "/header/geom/mks2/mksh0", H5T_IEEE_F64BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_double, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &mksh0);
  status = H5Dclose(dumps_dataset_double);

  status = H5Gclose(dumps_group_id3);
#endif //MKS3COORDS

  status = H5Gclose(dumps_group_id2); // close /header/geom
  status = H5Gclose(dumps_group_id);  // close /header

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////  
  // Body Data
  // loop over all cells and fill in the arrays with cell data

  ldouble x1_arr[NX][NY][NZ];
  ldouble x2_arr[NX][NY][NZ];
  ldouble x3_arr[NX][NY][NZ];
  ldouble gdet_arr[NX][NY][NZ];
  
  ldouble r_arr[NX][NY][NZ];
  ldouble th_arr[NX][NY][NZ];
  ldouble ph_arr[NX][NY][NZ];

  ldouble rho_arr[NX][NY][NZ];
  ldouble pgas_arr[NX][NY][NZ];
  ldouble uint_arr[NX][NY][NZ];
  
  ldouble u0_arr[NX][NY][NZ];
  ldouble u1_arr[NX][NY][NZ];
  ldouble u2_arr[NX][NY][NZ];
  ldouble u3_arr[NX][NY][NZ];

  ldouble U1_arr[NX][NY][NZ];
  ldouble U2_arr[NX][NY][NZ];
  ldouble U3_arr[NX][NY][NZ];

  #ifdef MAGNFIELD
  ldouble b0_arr[NX][NY][NZ];
  ldouble b1_arr[NX][NY][NZ];
  ldouble b2_arr[NX][NY][NZ];
  ldouble b3_arr[NX][NY][NZ];

  ldouble B1_arr[NX][NY][NZ];
  ldouble B2_arr[NX][NY][NZ];
  ldouble B3_arr[NX][NY][NZ];
  #endif
  
  #ifdef EVOLVEELECTRONS
  ldouble te_arr[NX][NY][NZ];
  ldouble ti_arr[NX][NY][NZ];
  ldouble gam_arr[NX][NY][NZ];
  #endif
  
  int iz,iy,ix,iv,i;
#pragma omp parallel for private(iz,iy,ix,iv,i)
  for(iz=0;iz<NZ;iz++)
  {
    for(iy=0;iy<NY;iy++)
    {
      for(ix=0;ix<NX;ix++)
      {

	// coordinates
        struct geometry geom,geomBL,geomBL0;
	fill_geometry(ix,iy,iz,&geom);
	fill_geometry_arb(ix,iy,iz,&geomBL0,OUTCOORDS);
	fill_geometry_arb(ix,iy,iz,&geomBL,OUTCOORDS2);

	ldouble r=geomBL.xx;
	ldouble th=geomBL.yy;
	ldouble ph=geomBL.zz;
        ldouble x1=geom.xx;
	ldouble x2=geom.yy;
	ldouble x3=geom.zz;
        ldouble gdet=geom.gdet;
	
        // primitives
        ldouble pp[NV];
	for(iv=0;iv<NV;iv++)
	{
	  if(doingavg)
	    pp[iv]=get_uavg(pavg,iv,ix,iy,iz);
	  else
	    pp[iv]=get_u(p,iv,ix,iy,iz);
	}

	
	// get output quantities in OUTCOORDS2 frame
	ldouble rho,uint,pgas,temp,Te,Ti,utcon[4],utcov[4],urel[4],bcon[4],bcov[4],Bcon[4],bsq;
	ldouble gamma=GAMMA;
	#ifdef CONSISTENTGAMMA
	gamma=pick_gammagas(ix,iy,iz);
        #endif

	if(doingavg)
	{
	  if(OUTCOORDS2!=OUTCOORDS)
	    my_err("to use avg in hdf5 output, OUTCOORDS2 must equal OUTCOORDS used at run!!\n");
	      		   
	  rho=get_uavg(pavg,RHO,ix,iy,iz);
	  uint=get_uavg(pavg,UU,ix,iy,iz);
	  pgas=get_uavg(pavg,AVGPGAS,ix,iy,iz);
	  temp=calc_PEQ_Tfromprho(pgas,rho,ix,iy,iz);                   

	  utcon[0]=get_uavg(pavg,AVGRHOUCON(0),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	  utcon[1]=get_uavg(pavg,AVGRHOUCON(1),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	  utcon[2]=get_uavg(pavg,AVGRHOUCON(2),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
	  utcon[3]=get_uavg(pavg,AVGRHOUCON(3),ix,iy,iz)/get_uavg(pavg,RHO,ix,iy,iz);
                 
          //NORMALIZE u^0 to be  consistent with u1,u2,u3
          fill_utinucon(utcon,geomBL.gg,geomBL.GG);
	  indices_21(utcon,utcov,geomBL.gg);

	  // conv vels to primitive VELR (but in OUTCOORDS)
	  conv_vels_ut(utcon, urel, VEL4, VELR, geomBL.gg, geomBL.GG);

          #ifdef MAGNFIELD
	  bsq=get_uavg(pavg,AVGBSQ,ix,iy,iz);
	  bcon[0]=get_uavg(pavg,AVGBCON(0),ix,iy,iz);
	  bcon[1]=get_uavg(pavg,AVGBCON(1),ix,iy,iz);
	  bcon[2]=get_uavg(pavg,AVGBCON(2),ix,iy,iz);
	  bcon[3]=get_uavg(pavg,AVGBCON(3),ix,iy,iz);

          //NORMALIZE b^0 to be orthogonal with u^\mu
	  bcon[0]=-dot3nr(bcon,utcov)/utcov[0];
	  indices_21(bcon,bcov,geomBL.gg);

          //NORMALIZE b^mu to be equal to B^2
	  ldouble alphanorm = bsq/dotB(bcon,bcov);
	  if(alphanorm<0.) my_err("alpha.lt.0 in b0 norm !!\n");
          for(i=0;i<4;i++) bcon[i]*=sqrt(alphanorm);

	  // back to primitive B^i (but in OUTCOORDS)
	  int j;
	  for(j=1;j<4;j++)
          {
             Bcon[j] = bcon[j]*utcon[0] - bcon[0]*utcon[j];
          }

          #endif //MAGNFIELD
		  
          #ifdef EVOLVEELECTRONS
	  
	  ldouble pe=get_uavg(pavg,AVGPE,ix,iy,iz);
	  ldouble pi=get_uavg(pavg,AVGPI,ix,iy,iz);
	  ldouble ne=get_uavg(pavg,RHO,ix,iy,iz)/MU_E/M_PROTON; 
	  ldouble ni=get_uavg(pavg,RHO,ix,iy,iz)/MU_I/M_PROTON; 

	  Te=pe/K_BOLTZ/ne;
	  Ti=pi/K_BOLTZ/ni;
      	 
          #endif                
        } //doingavg==1
	else 
	{
	  // transform primitives to OUTCOORDS
          #if defined(PRECOMPUTE_MY2OUT) && (OUTCOORDS==OUTCOORDS2)
          trans_pall_coco_my2out(pp,pp,&geom,&geomBL);
          #elif defined(PRECOMPUTE_MY2OUT)
          trans_pall_coco_my2out(pp,pp,&geom,&geomBL0); //transform to OUTCOORDS using precompute, then to OUTCOORDS2
          trans_pall_coco(pp,pp,OUTCOORDS,OUTCOORDS2,geomBL0.xxvec,&geomBL0,&geomBL);
          #else
          trans_pall_coco(pp, pp, MYCOORDS,OUTCOORDS2, geom.xxvec,&geom,&geomBL);
          #endif

	  // get scalars
	  rho=pp[0];
	  uint=pp[1];
	  pgas=(gamma-1.)*uint;
          temp=calc_PEQ_Tfromurho(uint,rho,ix,iy,iz);
          
          // get 4-velocity
          calc_ucon_ucov_from_prims(pp, &geomBL, utcon, utcov);

	  // conv vels to primitive VELR (but in OUTCOORDS)
	  conv_vels_ut(utcon, urel, VEL4, VELR, geomBL.gg, geomBL.GG);

	  
          #ifdef MAGNFIELD
          // calculate bcon 4 vector
	  calc_bcon_bcov_bsq_from_4vel(pp, utcon, utcov, &geomBL, bcon, bcov, &bsq);
        
          // back to primitive B^i (but in OUTCOORDS -- this should just be pp[B1],pp[B2],pp[B3])
          int j;
	  for(j=1;j<4;j++)
          {
             Bcon[j] = bcon[j]*utcon[0] - bcon[0]*utcon[j];
          }
          #endif

          #ifdef EVOLVEELECTRONS
	  temp=calc_PEQ_Teifrompp(pp,&Te,&Ti,geomBL.ix,geomBL.iy,geomBL.iz);
	  #endif
	} //doingavg==0

        x1_arr[ix][iy][iz] = x1;
        x2_arr[ix][iy][iz] = x2;
        x3_arr[ix][iy][iz] = x3;
        gdet_arr[ix][iy][iz] = gdet;
	
        r_arr[ix][iy][iz] = r;
        th_arr[ix][iy][iz] = th;
        ph_arr[ix][iy][iz] = ph;

        rho_arr[ix][iy][iz] = rho;
        pgas_arr[ix][iy][iz] = pgas;
        uint_arr[ix][iy][iz] = uint;
	
        u0_arr[ix][iy][iz] = utcon[0];
        u1_arr[ix][iy][iz] = utcon[1];
        u2_arr[ix][iy][iz] = utcon[2];
        u3_arr[ix][iy][iz] = utcon[3];

        U1_arr[ix][iy][iz] = urel[1];
        U2_arr[ix][iy][iz] = urel[2];
        U3_arr[ix][iy][iz] = urel[3];

        #ifdef MAGNFIELD
        b0_arr[ix][iy][iz] = bcon[0];
        b1_arr[ix][iy][iz] = bcon[1];
        b2_arr[ix][iy][iz] = bcon[2];
        b3_arr[ix][iy][iz] = bcon[3];

        B1_arr[ix][iy][iz] = Bcon[1];
        B2_arr[ix][iy][iz] = Bcon[2];
        B3_arr[ix][iy][iz] = Bcon[3];
        #endif
  
        #ifdef EVOLVEELECTRONS
        te_arr[ix][iy][iz] = Te;
        ti_arr[ix][iy][iz] = Ti;
	gam_arr[ix][iy][iz] = gamma;
        #endif

      }
    }
  }

  
  // Write runtime grid
  dumps_group_id = H5Gcreate2(dumps_file_id, "/grid_run", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/grid_run/gdet", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(gdet_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

 // !AC -- we don't need coordinates with the left corner defined in /header/geom
  /*
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/grid_run/x1", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(x1_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/grid_run/x2", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(x2_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/grid_run/x3", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(x3_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  status = H5Gclose(dumps_group_id);  // close /grid_run
  */
  
  // Write output grid
  dumps_group_id = H5Gcreate2(dumps_file_id, "/grid_out", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/grid_out/r", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(r_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/grid_out/th", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(th_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/grid_out/ph", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(ph_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  status = H5Gclose(dumps_group_id);  // close /grid_out

  // Write output quantities
  dumps_group_id = H5Gcreate2(dumps_file_id, "/quants", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/rho", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(rho_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

#ifdef ANAOUT_HDF5_V1
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/pgas", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(pgas_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/u0", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(u0_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/u1", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(u1_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/u2", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(u2_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/u3", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(u3_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  #ifdef MAGNFIELD
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/b0", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(b0_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/b1", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(b1_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/b2", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(b2_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/b3", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(b3_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);
  #endif
  
#else //ifdef ANAOUT_HDF5_V1
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/uint", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(uint_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/U1", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(U1_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/U2", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(U2_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/U3", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(U3_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  #ifdef MAGNFIELD

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/B1", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(B1_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/B2", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(B2_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/B3", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(B3_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);
  #endif
  
#endif  //ifdef ANAOUT_HDF5_V1
    
  #ifdef EVOLVEELECTRONS
  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/te", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(te_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);

  dumps_dataset_array = H5Dcreate2(dumps_file_id, "/quants/ti", H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, dumps_dataspace_array, H5P_DEFAULT,
		    &(ti_arr[0][0][0]));
  status = H5Dclose(dumps_dataset_array);
  #endif

  status = H5Gclose(dumps_group_id);  // close /quants

  // close file
  status = H5Sclose(dumps_dataspace_scalar);
  status = H5Sclose(dumps_dataspace_array);
  status = H5Fclose(dumps_file_id);

#endif  // ANAOUT_HDF5

  return 0;
}
