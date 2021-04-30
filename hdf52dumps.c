/*KORAL - hdf52dumps.c
 
 Converts .h5 files in ./dumps folder to header and binary dump files in the same folder
 
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "ko.h"
#include "hdf5.h"

#ifdef FOLDER_HDF5
#undef FOLDER_HDF5
#endif
#define FOLDER_HDF5 "./dumps"


int main(int argc, char **argv)
{  
#ifdef MPI
  printf("hdf52dumps does not work with MPI. Please undefine MPI and rerun\n");
  exit(1);
#endif
  
#if !defined DUMPS_READ_HDF5 || defined DUMPS_WRITE_HDF5
  printf("Please correct define.h\nYou must define DUMPS_READ_HDF5 and undefine DUMPS_WRITE_HDF5\n");
  exit(1);
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
  //calc_metric();

#ifdef RELELECTRONS
  set_relel_gammas();
#endif

  //precalculating problem related numbers
  int j;
#ifdef PR_PREPINIT
#include PR_PREPINIT
#endif
  
  int ifile;
  int ret;
  
  ldouble t;

  printf("Working on files #%04d to #%04d with %d step \n", no1, no2, nostep);

  for(ifile = no1; ifile <= no2; ifile += nostep)
  {
    // Read in dump data from .h5 file

    ret = fread_restartfile_serial_hdf5(ifile, folder, &t);
    printf("Finished reading file %d  ret = %d\n", ifile, ret);

    // Write out dump file into header and binary dump file

    nfout1 = ifile;
    ret = fprint_restartfile_bin(t, folder);
    printf("Finished writing file %d  ret = %d\n", ifile, ret);
  }
    
  return 0;
}



