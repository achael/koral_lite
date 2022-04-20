/*! \file ko.c
 */

/*! \mainpage KORAL - ko.c
  //radiative hydrodynamical code
*/

#include "ko.h"

int
main(int argc, char **argv)
{
  
  int i,j;
  int ix,iy,iz;
  
  //check input arguments
  if(argc!=NUM_INPUTARG+1)
  {
    my_err("Not enough input arguments.\n");
    return -1;
  }
  else
  {      
    for (i=0;i<NUM_INPUTARG;i++)
    {
      inputarg[i]=atof(argv[i+1]);
      printf("%d: %f\n",i,inputarg[i]);
    }
  }


  //**************
  //initialization
  //**************
  
  // Choose MPI, OMP or single processor
  #ifdef MPI
  mpi_myinit(argc,argv);
  #else
  omp_myinit();
  #endif

  // check if definitions agree with each other and there are no conflicts
  am_i_sane();

  // initialize pointers to entropy-related functions: fixed GAMMA or consistent GAMMA
  init_pointers();

  // folder to dump to
  char folder[100];
  sprintf(folder,"%s","dumps");

  //this is not a postprocessing script
  doingavg=0;
  doingpostproc=0;
  doingpostproc_avg=0;
  
  // initialize global time
  global_time=0.;

  //gsl errors off
  gsl_set_error_handler_off();
  
  //random number gen. initialization
  #ifdef SRANDSEED
  srand(SRANDSEED);
  #else
  srand(time(NULL));
  #endif

  //allocate memory for global arrays
  initialize_arrays();
  
  //initialize constants
  initialize_constants();
  
  //set gamma bins for relativistic electrons
  #ifdef RELELECTRONS
  set_relel_gammas();
  #endif

  //Set up the grid arrays x[] and xb[]
  set_grid(&min_dx,&min_dy,&min_dz,&max_dt);
  
  //allocate loops: Nloop_0 to Nloop_6, and arrays loop_0 to loop_6
  alloc_loops();

  //fill array of initial gammagas values
  fill_arrays_at_init();

  //save grid file
  #if(GRIDOUTPUT==1)
  fprint_gridfile(folder);
  #endif

  //precalculates metric, Christoffels, etc.
  calc_metric();
  calc_cells_under_horiz();
 
  //save coordinate file
  #ifdef COORDOUTPUT
  fprint_coordfile(folder,"coord");
  #endif

#ifdef RADIATION
  
  //prepare arrays for accelerating radiative viscosity
  #if(RADVISCOSITY==SHEARVISCOSITY)
  reset_radviscaccel();
  #endif

#endif 

  //print scalings GU->CGS
  if(PROCID==0) print_scalings();

  //**************
  //tests
  //**************

  //test_gammagas(); exit(1);
  //test_epsilon(); exit(1);
  //test_metric(); exit(1);
  //test_inversion(); exit(1);
  //test_inversion_nonrel(); exit(1);
  //test_solve_implicit_lab_file(); exit(1);  
  //test_solve_implicit_lab(); exit(1);
  //test_Gi();  exit(-1);
  //test_opacities(); exit(1);
  //test_Ccoupling(); exit(1);
  //test_calcgamma(); exit(1);
  //test_Tsynch(); exit(1);
  //test_heatfit(); exit(1);

  //**************
  // Precalculate problem related numbers
  //**************

  #ifdef PR_PREPINIT
  #include PR_PREPINIT
  #endif

  //*************************************
  // Startup from restart file
  //*************************************
  //printf("NV: %d NX: %d NY: %d NZ: %d \n",NV,NX,NY,NZ);
  int ifinit=1;
  ldouble tstart;  

#ifndef NORESTART
  
  struct timespec temp_clock;
  my_clock_gettime(&temp_clock);    
  start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

  // attempt to read restart file
  
  ifinit=fread_restartfile(RESTARTNUM,folder,&tstart);

  my_clock_gettime(&temp_clock);    
  end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
  if(PROCID==0)
    printf( "-- I/O took %.3f seconds\n", end_time-start_time);
  
  if(!ifinit) global_time=tstart;
  
  if(!ifinit) //succesfully read initial state from file
  {
      if(PROCID==0)
      {
        printf("Sending initial data... ");
	fflush(stdout);
      }

#ifdef FORCEFREE
      fill_ffprims(); // make force-free primitives consistent
#endif

      //exchange initial state
      mpi_exchangedata();  
      calc_avgs_throughout();
      set_bc(tstart,1);
      
      #ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif

      if(PROCID==0)
      {
	printf("done!\n");
	fflush(stdout);
      }

  }
  else
  {
      if(PROCID==0)
        printf("Failed to find restart file \n");      
  }

#endif  // NORESTART
  
  //********************************************
  // NORESTART defined  or no restart file found.
  // Initialize problem via set_initial_profile
  //********************************************
  
  if(ifinit==1)
  {
      struct timespec temp_clock;
      my_clock_gettime(&temp_clock);    
      start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

      // Initialize new problem via init.c in the appropriate PROBLEMS folder
      set_initial_profile();
      tstart=0.;
      
      //exchange initial state
      if(PROCID==0)
      {
        printf("Sending initial data... ");
	fflush(stdout);
      }
      mpi_exchangedata();
      calc_avgs_throughout();
      set_bc(tstart,1);
      
      #ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      
      if(PROCID==0)
      {
	printf("done!\n");
	fflush(stdout);
      }

      //compute magnetic field for initial configuration
#ifdef MAGNFIELD
#ifdef VECPOTGIVEN
      if(PROCID==0)
      {
	printf("Calculating magn. field... \n");
	fflush(stdout);
      }
      calc_BfromA(p,1);

#ifdef FORCEFREE
      fill_ffprims(); // make force-free primitives consistent
#endif
      //exchange magn. field calculated in domain
      mpi_exchangedata();
      calc_avgs_throughout();
      set_bc(tstart, 1);
      
      #ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      if(PROCID==0)
      {
	printf("done!\n");
	fflush(stdout);
      }
#endif  // VECPOTGIVEN
#endif  // MAGNFIELD

      // other  post-init code goes here
      #ifdef PR_POSTINIT
      #include PR_POSTINIT
      #endif

      my_clock_gettime(&temp_clock);
      end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

      if(PROCID==0)
        printf( "-- initialization took %.3f seconds\n", end_time-start_time);

  }  // if(ifinit==1)

  //*****************************
  // Rescalings and Perturbations
  //*****************************
  int iix,iiy,iiz;
  struct geometry geom;
  ldouble PERTURBATION;

  
  // reset Te so that ue/ugas is constant
#ifdef RESETELECTRONTEMPERATURETOUEUGASRATIO

  if(PROCID==0) printf("RRESETELECTRONTEMPERATURETOUEUGASRATIO: %f\n",RESETELECTRONTEMPERATURETOUEUGASRATIO);
 
  for(i=0;i<Nloop_0;i++) //domain                                                                                                                               
  {
      iix=loop_0[i][0];
      iiy=loop_0[i][1];
      iiz=loop_0[i][2];

      #ifdef EVOLVEELECTRONS
      ldouble uint=get_u(p,UU,iix,iiy,iiz);
      ldouble rho=get_u(p,RHO,iix,iiy,iiz);
      ldouble ue,ui,se,si;
      ue=RESETELECTRONTEMPERATURETOUEUGASRATIO*uint;
      ui=(1.-RESETELECTRONTEMPERATURETOUEUGASRATIO)*uint;
      ldouble rhoeth=MU_E*M_PROTON*calc_thermal_ne(&get_u(p,0,iix,iiy,iiz));
      se=calc_Sefromrhou(rhoeth,ue,ELECTRONS);
      si=calc_Sefromrhou(rho,ui,IONS);
      
      set_u(p,ENTRE,iix,iiy,iiz,se);
      set_u(p,ENTRI,iix,iiy,iiz,si);
      #endif  // EVOLVEELECTRONS
   

      fill_geometry(iix,iiy,iiz,&geom);
      p2u(&get_u(p,0,iix,iiy,iiz),&get_u(u,0,iix,iiy,iiz),&geom);
  } 
#endif  // RESETELECTRONTEMPERATURETOUEUGASRATIO

  // Rescale density to adjust  accretion rate
#ifdef RESCALEDENSITY

  if(PROCID==0) printf("RESCALEDENSITY: %e\n",RESCALEDENSITY);
 
  for(i=0;i<Nloop_0;i++)                                                                                                                                
  {
      iix=loop_0[i][0];
      iiy=loop_0[i][1];
      iiz=loop_0[i][2];

      #ifdef EVOLVEELECTRONS
      ldouble Te, Ti; //calculate electron and ion temperatures to keep them fixed when rescaling
      calc_PEQ_Teifrompp(&get_u(p,0,iix,iiy,iiz),&Te, &Ti, iix,iiy,iiz);
      #endif

      set_u(p,RHO,iix,iiy,iiz,get_u(p,RHO,iix,iiy,iiz)*RESCALEDENSITY);
      set_u(p,UU,iix,iiy,iiz,get_u(p,UU,iix,iiy,iiz)*RESCALEDENSITY);
      set_u(p,B1,iix,iiy,iiz,get_u(p,B1,iix,iiy,iiz)*sqrt(RESCALEDENSITY));
      set_u(p,B2,iix,iiy,iiz,get_u(p,B2,iix,iiy,iiz)*sqrt(RESCALEDENSITY));
      set_u(p,B3,iix,iiy,iiz,get_u(p,B3,iix,iiy,iiz)*sqrt(RESCALEDENSITY));
      
      #ifdef EVOLVEELECTRONS
      ldouble rhoeth=MU_E*M_PROTON*calc_thermal_ne(&get_u(p,0,iix,iiy,iiz));
      set_u(p,ENTRE,iix,iiy,iiz,calc_S3fromrhoT(rhoeth,Te,ELECTRONS));
      set_u(p,ENTRI,iix,iiy,iiz,calc_S3fromrhoT(get_u(p,RHO,iix,iiy,iiz),Ti,IONS));
      #ifdef RELELECTRONS
      int ie;    
      for (ie=0; ie<NRELBIN; ie++) 
	set_u(p,NEREL(ie),iix,iiy,iiz,get_u(p,NEREL(ie),iix,iiy,iiz)*RESCALEDENSITY);
      #endif //RELELECTRONS
      #endif //EVOLVEELECTRONS

      #ifdef RADIATION
      set_u(p,EE,iix,iiy,iiz,get_u(p,EE,iix,iiy,iiz)*RESCALEDENSITY);
      #ifdef EVOLVEPHOTONNUMBER
      set_u(p,NF,iix,iiy,iiz,get_u(p,NF,iix,iiy,iiz)*RESCALEDENSITY);
      #endif //EVOLVEPHOTONNUMBER
      #endif //RADIATION
      
      fill_geometry(iix,iiy,iiz,&geom);
      p2u(&get_u(p,0,iix,iiy,iiz),&get_u(u,0,iix,iiy,iiz),&geom);
  } 
#endif  // RESCALEDENSITY

  //Perturb azimuthal velocity to trigger MRI
#ifdef PERTURBVZ

  if(PROCID==0) printf("AZIMUTHAL VELOCITY PERTURBATION: %f\n",PERTURBVZ);

  for(i=0;i<Nloop_0;i++) 
  {
        
      iix=loop_0[i][0];
      iiy=loop_0[i][1];
      iiz=loop_0[i][2];
     
      PERTURBATION =  2.*PERTURBVZ*(((double)rand()/(double)(RAND_MAX)) - 0.5);
        
      set_u(p,VZ,iix,iiy,iiz,get_u(p,VZ,iix,iiy,iiz)*(1.+PERTURBATION));
        
      fill_geometry(iix,iiy,iiz,&geom);
      p2u(&get_u(p,0,iix,iiy,iiz),&get_u(u,0,iix,iiy,iiz),&geom);
  }
#endif  // PERTURBVZ

  //Perturb density to trigger MRI
#ifdef PERTURBRHO

  if(PROCID==0) printf("DENSITY PERTURBATION: %f\n",PERTURBRHO);

  for(i=0;i<Nloop_0;i++) 
  {
        
      iix=loop_0[i][0];
      iiy=loop_0[i][1];
      iiz=loop_0[i][2];
      
      PERTURBATION =  2.*PERTURBRHO*(((double)rand()/(double)(RAND_MAX))-0.5 );
        
      set_u(p,RHO,iix,iiy,iiz,get_u(p,RHO,iix,iiy,iiz)*(1.+PERTURBATION));
        
      fill_geometry(iix,iiy,iiz,&geom);
      p2u(&get_u(p,0,iix,iiy,iiz),&get_u(u,0,iix,iiy,iiz),&geom);
  }
#endif  // PERTURBRHO
      
  //Perturb gas internal energy to trigger MRI    
#ifdef PERTURBUINT
  
  if(PROCID==0) printf("internal energy PERTURBATION: %f\n",PERTURBUINT);
    
  for(i=0;i<Nloop_0;i++)
  {
        
      iix=loop_0[i][0];
      iiy=loop_0[i][1];
      iiz=loop_0[i][2];
      
      PERTURBATION =  2.*PERTURBUINT*(((double)rand()/(double)(RAND_MAX) ) - 0.5 );
        
      set_u(p,UU,iix,iiy,iiz,get_u(p,UU,iix,iiy,iiz)*(1.+PERTURBATION));
        
      fill_geometry(iix,iiy,iiz,&geom);
      p2u(&get_u(p,0,iix,iiy,iiz),&get_u(u,0,iix,iiy,iiz),&geom);
  }
#endif  // ifdef PERTURBUINT
      
  
  //********************
  // Prepare initial arrays
  //********************

  fprint_openfiles(folder);

  //copies initial primitives to pinit
  copy_u(1.,p,pinit);

  //zero the avg array
  long long Ngrid=SX*SY*SZ;
  long long Navg=Ngrid*(NV+NAVGVARS);
  copy_u_core(0.,pavg,pavg,Navg);	  
  avgtime=0.;

  //zero the self time array
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
	set_u_scalar(avgselftime,ix,iy,iz,0.);

  //********************
  // Save initial files if starting from scratch
  //********************
  if(ifinit==1)
  {
    fprint_restartfile(tstart,folder);

#ifndef MPI // dump on-the-fly only for shared memory
    #if(SCAOUTPUT==1) //scalar dumps
    fprint_scalars(tstart,scalars,NSCALARS);
    #endif
      
    #if(RADOUTPUT==1) //radial profiles
    fprint_radprofiles(tstart,nfout1,folder,"rad");
    #endif

    #if(THOUTPUT==1) //theta profiles
    fprint_thprofiles(t,nfout1,folder,"th");
    #endif
	  
    #if(SILOOUTPUT==1) //silo files
    #ifndef NOSILO
    fprint_silofile(tstart,nfout1,folder,"sil");
    #endif
    #endif 
      
    #if(SIMOUTPUT!=0) //simple files
    fprint_simplefile(tstart,nfout1,folder,"sim");
    #endif

    #if(PRIMOUTPUT!=0) // primitive file
    fprint_primitive_file(tstart,nfout1,folder,"prim");
    #endif

    #if(RELELSPECTRUMOUTPUT==1) //nonthermal spectrum
    fprint_relel_spectrum(tstart,NTH_SPEC_IX,NTH_SPEC_IY,NTH_SPEC_IZ,nfout1,folder,"spe",0);
    #endif

#endif  // ifndef MPI

    nfout1++; //total number of output
  }  // if(ifinit==1)


  //******************************************************
  // Main evolution is all here
  //******************************************************
  //#include <gperftools/profiler.h>
  //ProfilerStart("profile.log");
  solve_the_problem(tstart, folder);
  //ProfilerStop();
  //******************************************************
  // Evolution done. Clean up and finish.
  //******************************************************
  
  free_arrays();
  fprint_closefiles();
  mpi_myfinalize();

  return 0;
} 



