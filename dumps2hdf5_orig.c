/*KORAL - dumps2hdf5.c
 
 Converts binary dump files to hdf5 format and saves as .h5 files in specified folder.
 
 By default the .h5 files will be created in the same directory "./dumps" where the original binary dumps are located. If the .h5 files should go into a different directory, add in define.h a command like the following:
 
      #define FOLDER_HDF5 "./hdf5"
 
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "ko.h"
#include "hdf5.h"

#if defined(PR_WRITE_OUTPUT) && !defined(INCLUDE_WRITE_OUTPUT)
#include PR_WRITE_OUTPUT
#define INCLUDE_WRITE_OUTPUT
#endif


int main(int argc, char **argv)
{  
#ifdef MPI
  printf("dumps2hdf5 does not work with MPI. Please undefine MPI and rerun\n");
  exit(1);
#endif
  
#ifdef DUMPS_READ_HDF5
#undef DUMPS_READ_HDF5
#endif
  
#ifndef DUMPS_WRITE_HDF5
#define DUMPS_WRITE_HDF5
#endif

#ifndef FOLDER_HDF5
#define FOLDER_HDF5 "./dumps"
#endif

  //initialize pointers to functions
  init_pointers();

  //initialize constants
  initialize_constants();

  //which files to read
  int no1, no2, nostep;
  if(argc != 4)
  {
    printf("Wrong number of input arguments. Command is ./dumps2hdf5 init final step\n");
    exit(1);
  }
  else
  {
    no1=atof(argv[1]);
    no2=atof(argv[2]);
    nostep=atof(argv[3]);
  }

  char folder[100],bufer[100];
  sprintf(folder,"%s","./dumps");

  //no gsl error messages
  gsl_set_error_handler_off();
  
  //preparing arrays
  initialize_arrays();

  //sets the grid
  set_grid(&min_dx,&min_dy,&min_dz,&max_dt);
  alloc_loops(1, 0., 0.);

  //precalculates metric etc.
  calc_metric();

#ifdef RELELECTRONS
  set_relel_gammas();
#endif

  //precalculating problem related numbers
  int j;
#ifdef PR_PREPINIT
#include PR_PREPINIT
#endif
  
  int ifile, itot=0, readret;
  int ret, ix, iy, iz, iv, i, ic, gix, giy, giz;
  int indices[NX*NY*NZ][3];
  //int ixx[NX][NY][NZ], iyy[NX][NY][NZ], izz[NX][NY][NZ];
  
  double t, ttot;
  double primitive[NX][NY][NZ];
  ttot=0.;

  printf("Working on files #%04d to #%04d with %d step \n", no1, no2, nostep);

  ldouble pp[NV];

  for(ifile = no1; ifile <= no2; ifile += nostep)
  {
    itot++;
    
    //Set up file names to be read
    
    char fnamehead[400], fname[400];
    sprintf(fnamehead, "%s/res%04d.head", folder, ifile);
    sprintf(fname, "%s/res%04d.dat", folder, ifile);
    FILE *fdump;
    
    //Read data from header file
    
    fdump = fopen(fnamehead,"r");
    if (fdump == NULL)
    {
      printf("Restart file #%d not available\n", ifile);
      exit(1);
    }
    int intpar[6];
    ret = fscanf(fdump, "## %d %d %lf %d %d %d %d\n", &intpar[0], &intpar[1], &t, &intpar[2], &intpar[3], &intpar[4], &intpar[5]);
    printf("Restart file (%s) read no. %d at time: %f of PROBLEM: %d with NXYZ: %d %d %d\n", fname, intpar[0], t, intpar[2], intpar[3], intpar[4], intpar[5]);
    
    fclose(fdump);

    // Write out header information in group HEADER in HDF5 file
    
    hid_t dumps_file_id, dumps_group_id, dumps_dataspace_scalar, dumps_dataspace_array, dumps_dataset_int, dumps_dataset_double, dumps_dataset_array, dumps_attribute_id;
    hsize_t dims_h5[3];
    herr_t status;
    
    int file_number = intpar[0], file_avg = intpar[1], problem_number = intpar[2], nxx = intpar[3], nyy = intpar[4], nzz = intpar[5];
    
    if (nxx != TNX || nyy != TNY || nzz != TNZ)
    {
      printf("Array size does not match\nNX, NY, NZ = %d %d %d\nnxx, nyy, nzz = %d %d %d\n", TNX, TNY, TNZ, nxx, nyy, nzz);
      exit(1);
    }
    
    dims_h5[0] = NX;
    dims_h5[1] = NY;
    dims_h5[2] = NZ;
    
    char fname_h5[256];
    sprintf(fname_h5, "%s/res%04d.h5", FOLDER_HDF5, ifile);
    dumps_file_id = H5Fcreate (fname_h5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dumps_group_id = H5Gcreate2(dumps_file_id, "/HEADER", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    dumps_dataspace_scalar = H5Screate(H5S_SCALAR);

    dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/FILE_NUMBER", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

    int nprimitives = NV;
    dumps_dataset_int = H5Dcreate2(dumps_file_id, "/HEADER/NPRIM", H5T_STD_I32BE, dumps_dataspace_scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dumps_dataset_int, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      &nprimitives);
    status = H5Dclose(dumps_dataset_int);
    
    status = H5Gclose(dumps_group_id);
    
    // Read indices from dumps datafile
    
    fdump=fopen(fname,"rb");
    
    //ldouble uu[NV],pp[NV],ftemp;
    char prim_name[256], attribute_name[256];

    for(ic = 0; ic < NX*NY*NZ; ic++)
    {
      ret = fread(&gix, sizeof(int), 1, fdump);
      ret = fread(&giy, sizeof(int), 1, fdump);
      ret = fread(&giz, sizeof(int), 1, fdump);
      
      mpi_global2localidx(gix, giy, giz, &ix, &iy, &iz);
      
      //ixx[ix][iy][iz] = gix;
      //iyy[ix][iy][iz] = giy;
      //izz[ix][iy][iz] = giz;
      
      indices[ic][0] = ix;
      indices[ic][1] = iy;
      indices[ic][2] = iz;
    }

    dumps_dataspace_array = H5Screate_simple(3, dims_h5, NULL);
    
    /*    
    // Save indices in HDF5 file. Is this needed? If not, get rid of ixx, iyy, izz arrays
    
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
    
    // Read primitives and save in array p
    
    for(ic = 0; ic < NX * NY * NZ; ic++)
    {
      ret=fread(pp, sizeof(ldouble), NV, fdump);
      
      ix = indices[ic][0];
      iy = indices[ic][1];
      iz = indices[ic][2];
      
      for(iv = 0; iv < NV; iv++)
      {
        set_u(p, iv, ix, iy, iz, pp[iv]);
      }
    }
    
    // Now save primitives one by one in HDF5 file
    
    int iattribute = 0;
    for(iv = 0; iv < NV; iv++)
    {
      for(ic = 0; ic < NX * NY * NZ; ic++)
      {
        ix = indices[ic][0];
        iy = indices[ic][1];
        iz = indices[ic][2];

        primitive[ix][iy][iz] = get_u(p, iv, ix, iy, iz);
      }
      
      //sprintf(prim_name, "/PRIM%d", iv);
      get_prim_name(prim_name, iv);
      dumps_dataset_array = H5Dcreate2(dumps_file_id, prim_name, H5T_IEEE_F64BE, dumps_dataspace_array, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(dumps_dataset_array, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, primitive);
      status = H5Dclose(dumps_dataset_array);
    }
    
    status = H5Sclose(dumps_dataspace_scalar);
    status = H5Sclose(dumps_dataspace_array);
    status = H5Fclose (dumps_file_id);
  }

  return 0;
}



