//KORAL - ana.c
//dump files postprocessing

#include "ko.h"

//problem-specific output
#if defined(PR_WRITE_OUTPUT) && !defined(INCLUDE_WRITE_OUTPUT)
#include PR_WRITE_OUTPUT
#define INCLUDE_WRITE_OUTPUT
#endif

int 
main(int argc, char **argv)
{  

  #ifdef MPI
  mpi_myinit(argc,argv);
  if(PROCID==0)
      printf("ana works mostly on shared memory only, do not use MPI\n");
  #endif
  
  #ifdef OMP
  omp_myinit();  
  #endif

  //initialize pointers to functions
  init_pointers();

  //initialize constants
  initialize_constants();

  //which files to read
  int no1,no2,nostep,ifphiavg,ifrunonthego;
  if(argc<4 || argc>6)
  {
    printf("Wrong number of input arguments.\n");
    printf("./ana no1 no2 nostep [[0-snap (default)/1-phiavg/2-phisli/3-thsli] ifrunonthego=0]\n");
    return -1;
  }
  else
  {      
    no1=atof(argv[1]); //start file
    no2=atof(argv[2]); //stop file
    nostep=atof(argv[3]); //step through files

    ifphiavg=0; //ph averaged or not
    ifrunonthego=0; //do evolution from each output file to get evolution-dependent quantities. 

    if(argc==5)
    {
      ifphiavg=atof(argv[4]);
    }

    if(argc==6)
    {
      ifphiavg=atof(argv[4]);
      ifrunonthego=atof(argv[5]);
    }
  }

  if (ifphiavg>0 && ifphiavg<3 && NZ!=1)
  {
    printf("ifphiavg=%d -- make sure code is compiled with NZ=1!\n",ifphiavg);
    exit(-1);
  }

  if (ifphiavg==3 && ifphiavg<3 && NY!=1)
  {
    printf("ifphiavg=3 -- make sure code is compiled with NY=1!\n");
    exit(-1);
  }

  doingavg=0; //ana never runs on avg files
  doingpostproc=1; //we are post-processing
  if(ifrunonthego)
  {
    doingpostproc=0; //if one wants to run the evolution then all arrays must be allocated
    if(PROCID==0)
    {
      printf("------\n"
	     "Runonthego requested: after each restart file read in the code will \n"
	     "proceed with the evolution. Define NSTEPSTOP to specify the number\n"
	     "of steps required and make sure DTOUTs are large enough to prevent\n"
	     "any outputs. \n"
	     "-------\n");
      #ifdef RESCALEDENSITYPOSTPROC
      printf("Unclear how runonthego interacts with  RESCALEDENSITYPOSTPROC! -- avoid for now\n");
      exit(-1);
      #endif
    }
  }

  //folder to write to
  char folder[100],bufer[100];
  if(ifphiavg==0)
    sprintf(folder,"%s","dumps");
  else if(ifphiavg==1)
    sprintf(folder,"%s","dumps_phiavg");
  else if(ifphiavg==2)
    sprintf(folder,"%s","dumps_phisli");
  else if(ifphiavg==3)
    sprintf(folder,"%s","dumps_thsli");
  
  //no gsl error messages
  gsl_set_error_handler_off();
  
  //random number gen. initialization
  srand ( time(NULL) );

  //preparing arrays
  initialize_arrays();

  //sets the grid
  set_grid(&min_dx,&min_dy,&min_dz,&max_dt);
  alloc_loops();

  //precalculates metric etc.
  calc_metric();

  #ifdef RELELECTRONS
  set_relel_gammas();
  #endif

  //outputs coordinate files
  #ifdef COORDOUTPUT
  fprint_coordfile("analysis","coord");
  #endif

  //precalculating problem related numbers
  int i,j;
  #ifdef PR_PREPINIT
  #include PR_PREPINIT
  #endif
  
  //open the scalar files	
  char scalarsname[100];
  #if(SCAOUTPUT==1)					
  sprintf(scalarsname,"analysis/scalars_%04d_%04d_%04d.dat",no1,no2,nostep);  
  fout_scalars=fopen(scalarsname,"w");
  #endif

  //loop through files
  int ifile,itot=0,readret;
  ldouble pp[NV],uu[NV]; 
  ldouble t,ttot; 
  ttot=0.;

  if(PROCID==0)
  {
    printf("working on files #%04d to #%04d with %d step \n",no1,no2,nostep);
  }

  for(ifile=no1;ifile<=no2;ifile+=nostep)
  {
    itot++;

    //read restart file
    readret=fread_restartfile(ifile,folder,&t);
    nfout1=ifile;
    global_time=t;

    //exchange initial state
    mpi_exchangedata();
    
    //calculates scaleheight etc.
    calc_avgs_throughout();

    //set bc
    set_bc(t,1);

    //run the evolution for a bit if requested
    if(ifrunonthego)
    {
      //copies initial primitives to pinit
      copy_u(1.,p,pinit);

      //evolves
      solve_the_problem(t, "dummy");

      //rewrite *ppostimplicit to *p and use below for outputs
      copyi_u(1.,ppostimplicit,p);
    }
         
#ifdef RESCALEDENSITYPOSTPROC  // Rescale density in post processing
    //ANDREW -- what about ghost cells?
    if(PROCID==0) printf("RESCALEDENSITYPOSTPROC: %f\n",RESCALEDENSITYPOSTPROC); 
    int iix,iiy,iiz;
    for(iiz=0;iiz<NZ;iiz++)
      for(iiy=0;iiy<NY;iiy++)
	for(iix=0;iix<NX;iix++)
	{
	    //calculate electron and ion temperatures to keep them fixed when rescaling
            #ifdef EVOLVEELECTRONS 
	    ldouble Te, Ti; 
	    calc_PEQ_Teifrompp(&get_u(p,0,iix,iiy,iiz),&Te, &Ti, iix,iiy,iiz);
            #endif

	    //rescale density, energy density, and b field
	    set_u(p,RHO,iix,iiy,iiz,get_u(p,RHO,iix,iiy,iiz)*RESCALEDENSITYPOSTPROC);
	    set_u(p,UU,iix,iiy,iiz,get_u(p,UU,iix,iiy,iiz)*RESCALEDENSITYPOSTPROC);
	    set_u(p,B1,iix,iiy,iiz,get_u(p,B1,iix,iiy,iiz)*sqrt(RESCALEDENSITYPOSTPROC));
	    set_u(p,B2,iix,iiy,iiz,get_u(p,B2,iix,iiy,iiz)*sqrt(RESCALEDENSITYPOSTPROC));
	    set_u(p,B3,iix,iiy,iiz,get_u(p,B3,iix,iiy,iiz)*sqrt(RESCALEDENSITYPOSTPROC));

	    //restore electron and ion  temps
            #ifdef EVOLVEELECTRONS
	    ldouble rhoeth=MU_E*M_PROTON*calc_thermal_ne(&get_u(p,0,iix,iiy,iiz));
	    set_u(p,ENTRE,iix,iiy,iiz,calc_S3fromrhoT(rhoeth,Te,ELECTRONS));
	    set_u(p,ENTRI,iix,iiy,iiz,calc_S3fromrhoT(get_u(p,RHO,iix,iiy,iiz),Ti,IONS));

	    //rescale relativistic electron density
            #ifdef RELELECTRONS
	    int ie;    
	    for (ie=0; ie<NRELBIN; ie++) 
	    {
		set_u(p,NEREL(ie),iix,iiy,iiz,get_u(p,NEREL(ie),iix,iiy,iiz)*RESCALEDENSITYPOSTPROC);
	    }
            #endif //RELELECTRONS
            #endif //EVOLVEELECTRONS 

	    //Scale radiation energy and photon density
            #ifdef RADIATION
	    set_u(p,EE,iix,iiy,iiz,get_u(p,EE,iix,iiy,iiz)*RESCALEDENSITYPOSTPROC);
            #ifdef EVOLVEPHOTONNUMBER
	    set_u(p,NF,iix,iiy,iiz,get_u(p,NF,iix,iiy,iiz)*RESCALEDENSITYPOSTPROC);
            #endif
            #endif
		
	    //p2u scaled primitives
	    struct geometry geom;
	    fill_geometry(iix,iiy,iiz,&geom);
	    p2u(&get_u(p,0,iix,iiy,iiz),&get_u(u,0,iix,iiy,iiz),&geom);
	}  
#endif //RESCALEDENSITYPOSTPROC

    //calculate scalars
    #if(SCAOUTPUT==1)
    ldouble scalars[NSCALARS];
    calc_scalars(scalars,t);
    #endif

    //suffix and prefix for saved files depending on phiavg
    char prefix[40];
    char suffix[10];
    sprintf(suffix,"");

    //sprintf(suffix,"");
    if(ifphiavg==1)
      sprintf(suffix,"%sphiavg",suffix);
    if(ifphiavg==2)
      sprintf(suffix,"%sphisli",suffix);
    if(ifphiavg==3)
      sprintf(suffix,"%sthsli",suffix);

    //th-sliced - these save files only make sense for th-slices    
    if(ifphiavg==3) 
    {
      //silo output
#if(SILOOUTPUT==1)
#ifndef NOSILO
      sprintf(prefix,"sil%s",suffix);  
      fprint_silofile(t,nfout1,"analysis",prefix);
#endif
#endif

      //sim output
#if(SIMOUTPUT!=0)
      sprintf(prefix,"sim%s",suffix);
      fprint_simplefile(t,nfout1,"analysis",prefix);
#endif
    }
    
    //th-sliced - these save files only make sense for phi-slices        
    else if(ifphiavg==2) //phisliced - only these below make sense for phi-slices
    {

      //silo output
#if(SILOOUTPUT==1)
#ifndef NOSILO
      sprintf(prefix,"sil%s",suffix);  
      fprint_silofile(t,nfout1,"analysis",prefix);
#endif
#endif

      //sim output
#if(SIMOUTPUT!=0)	  
      sprintf(prefix,"sim%s",suffix);  
      fprint_simplefile(t,nfout1,"analysis",prefix);
#endif
    }
    
    //regular dump files
    else
    {

      //scalar output
#if(SCAOUTPUT==1)
      fprint_scalars(t,scalars,NSCALARS);
#endif

      //radial profiles
#if(RADOUTPUT==1)
      sprintf(prefix,"rad%s",suffix);  
      fprint_radprofiles(t,nfout1,"analysis",prefix);
#endif

      //theta profiles
#if(THOUTPUT==1)
      sprintf(prefix,"th%s",suffix);  
      fprint_thprofiles(t,nfout1,"analysis",prefix);
#endif

      //silo output
#if(SILOOUTPUT==1)
#ifndef NOSILO
      sprintf(prefix,"sil%s",suffix);  
      fprint_silofile(t,nfout1,"analysis",prefix);
#endif
#endif
      
      //sim output
#if(SIMOUTPUT!=0)	  
      sprintf(prefix,"sim%s",suffix);  
      fprint_simplefile(t,nfout1,"analysis",prefix);
#endif

      //hdf5 analysis output (replaces sim)
#ifdef ANAOUT_HDF5
      #ifdef ANAOUT_HDF5_V1
      sprintf(prefix,"ana%s",suffix); //old version for grtrans
      #else
      sprintf(prefix,"ipole%s",suffix); //new version for ipole
      #endif
      fprint_anaout_hdf5(t, "analysis",prefix);
#endif
      
      //relativistic electron spectrum
#if(RELELSPECTRUMOUTPUT==1)
      int ixx, iyy, izz;
      ixx=50; iyy=127; izz=0;   //hardcoded -- equator at approx 10M
      fprint_relel_spectrum(t, ixx, iyy, izz, nfout1, folder, "spe",0);
      fprint_relel_avg_spectrum(t, ixx, iyy, izz, nfout1, folder, "spe_avg",0);
#endif
      
    } //if(phiavg==3)
  }  //for(ifile=no1;ifile<=no2;ifile+=nostep)

  //Done with all res files  
  //close scalar files
#if(SCAOUTPUT==1)
  fclose(fout_scalars);
  char cpcommand[200];
  if(PROCID==0)
    {
      sprintf(cpcommand, "cp %s analysis/scalars.dat",scalarsname);
      system(cpcommand);
    }
#endif

  return 0;
}

