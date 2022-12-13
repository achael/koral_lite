//KORAL - avg.c
//avg files postprocessing

//guide to ifavg & ifavgavg combinations: 
//ifavg=1 + ifoutavg=1 : run on avgavg files
//ifavg=0 + ifoutavg=1 : run on avgres files
//ifavg=0 + ifoutavg=0 : run on res files (really should be using ana instead!)
//ifavg=1 + ifoutavg=0 : run on avg files


#include "ko.h"

int 
main(int argc, char **argv)
{  
  #ifdef MPI
  printf("avg works on shared memory only, do not use MPI, please\n");
  exit(-1);
  #endif

  #ifdef AVGOUTPUT
  if(AVGOUTPUT==2)
    printf("AVGOUTPUT==2 - avg files were dumped in 2d, copy to dumps_phiavg, and use as phiavg\n");
  #endif

  #ifdef AVGAVGOUTPUT
  printf("AVGAVGOUTPUT is deprecated -- run the separate script outavg.c instead!");
  exit(-1);
  #endif
  
  #ifdef OMP
  omp_myinit();  
  #endif

  //initialize pointers to functions
  init_pointers();

  //initialize constants
  initialize_constants();
  
  //which files to read
  int no1,no2,nostep,procotg,ifphiavg,ifavg,ifoutavg;
  if(argc!=8 && argc!=4 && argc!=7)
  {
    printf("Not enough input arguments.\n");
    printf("Asks for ./avg no1 no2 nostep [ifavg=1 procotg=0 ifphiavg=0 ifoutavg=0]\n");
    printf("NOTE: proctog is DEPRECATED and does nothing\n");
    return -1;
  }
  else
  {
    no1=atof(argv[1]); //the first file
    no2=atof(argv[2]); //the second file
    nostep=atof(argv[3]); //step size through files
    if(argc==8)
    {
      ifavg=atoi(argv[4]); //run on avg files or average res files
      procotg=atoi(argv[5]); //DEPRECATED -- print scalars
      ifphiavg=atoi(argv[6]); //run on phiavg files or not
      ifoutavg=atoi(argv[7]); //run on average of avg (or res) files or not
    }
    else if(argc==7)
    {
      ifavg=atoi(argv[4]);
      procotg=atoi(argv[5]);
      ifphiavg=atoi(argv[6]);
      ifoutavg=0;
    }
    else
    {
      ifavg=1;
      procotg=0;
      ifphiavg=0;
      ifoutavg=0;
    }
  }

  if (ifphiavg>0 && NZ!=1)
  {
    printf("ifphiavg>0 -- make sure code is compiled with NZ=1!\n");
    exit(-1);
  }
 
  //folder to write to
  char folder[100];
  sprintf(folder,"analysis");
  
  //folder to write from
  char folderin[100],bufor[100],base[100]; 
  if(ifphiavg==0)
    sprintf(folderin,"%s","dumps");
  else if(ifphiavg==1)
    sprintf(folderin,"%s","dumps_phiavg");
  else if(ifphiavg==2)
    my_err("using phisli(ced) data with ./avg makes no sense.\n");

  //file base names
  int no1_loop,no2_loop;
  if (ifoutavg)
  {
    if(ifavg)
      sprintf(base,"%s/avgavg%04d-",folderin,no1);   
    else
      sprintf(base,"%s/avgres%04d-",folderin,no1);

    no1_loop=no2;
    no2_loop=no2;
    doingavg=1; //we always run off of avg files ifoutavg==1

    printf("working on single file %s%04d.dat \n",base,no2);
  }
  else
  {
    if(ifavg)
    {
      sprintf(base,"%s/%s",folderin,"avg");
      doingavg=1; //run off normal avg files
    }
    else
    {
      sprintf(base,"%s/%s",folderin,"res");
      doingavg=0; //run off res files
    }
    no1_loop=no1;
    no2_loop=no2;

    printf("working on files %s%04d.dat to %s%04d.dat with %d step \n",base,no1,base,no2,nostep);
  }
  
  doingpostproc=1; 
  if (doingavg==0) doingpostproc_avg=1; //need flux arrays if running off res files
  else doingpostproc_avg=0;

  //currently gsl is not used
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

  //print coordinate file
  #ifdef COORDOUTPUT
  fprint_coordfile("analysis","coord");
  #endif

  #if(SCAOUTPUT==1)
  //opens the scalar file
  //if(procotg)
  //{
  sprintf(bufor,"%s/avgscalars.dat",folder);
  fout_scalars=fopen(bufor,"w");
  //}
  #endif

  //initialize avg array to 0
  long long i;  
  long long ncells = (SX)*(SY)*(SZ);
  long long avgarrsize=ncells*(NV+NAVGVARS);
  long long avgarrsizeMalloc=avgarrsize*sizeof(ldouble);

  ldouble *pavgtot=(ldouble*)malloc(avgarrsizeMalloc);
  for(i=0;i<avgarrsize;i++)
    pavg[i]=pavgtot[i]=0.;

  //load avg array from file
  int ifile,readret;
  int ix,iy,iz,iv;
  ldouble pp[NV],uu[NV],scalars[NSCALARS];
  ldouble ttot=0.;
  ldouble t=global_time;
  for(ifile=no1_loop; ifile<=no2_loop; ifile+=nostep)
  {
      if(doingavg) //reading avg file
      {
          #ifdef RESCALEDENSITYPOSTPROC
	  printf("Cannot RESCALEDENSITYPOSTPROC on avg file input!\n");
          printf("-- must run on res files (ifoutavg=0, ifavg=0)\n");
	  printf("Proceeding with avg with NO rescaling\n");
          #endif
	  
	  //read avg file
          readret=fread_avgfile(ifile,base,pavg,&dt,&t);	  

	  //flag errors in the avg file
	  int flag[NX][NY][NZ];
	  ldouble val;
	  for(iz=0;iz<NZ;iz++)
	    for(iy=0;iy<NY;iy++)
	      for(ix=0;ix<NX;ix++)
	      {
		flag[ix][iy][iz]=0;
		for(iv=0;iv<NV+NAVGVARS;iv++)
		{
		  val=get_uavg(pavg,iv,ix,iy,iz);
		  if(!isfinite(val))
		  {
		    //printf("%d  %d %d %d not finite!\n",iv,ix,iy,iz);
		    flag[ix][iy][iz]=-1;
		    //ANDREW encountered problems in these fields in simulations with no relel
                    #ifndef RELELECTRONS
		    if(iv>=AVGNRELEL && iv<=AVGNETH) flag[ix][iy][iz]=0;
		    #endif
		  }
		}
		if(get_uavg(pavg,RHO,ix,iy,iz)<SMALL || get_uavg(pavg,RHO,ix,iy,iz)>BIG)
		    flag[ix][iy][iz]=-1;
	      }

	  //correct for flagged errors
	  for(ix=0;ix<NX;ix++) 
	    for(iz=0;iz<NZ;iz++)
	      for(iy=0;iy<NY;iy++)
		if(flag[ix][iy][iz]<0) 
		{
		  if(TNZ==1) //corecting in 2d in theta
		  {
                    int cell_corrected = 0;
		    int iiy=iy;
		    if(iy>TNY/2) //lower half -- search toward equatorial plane
		    {
		      do
			iiy--;
		      while(flag[ix][iiy][iz]<0 && iiy>=0);
		    }
		    else //upper half -- search toward equatorial plane
		    {
		      do
		        iiy++;
		      while(flag[ix][iiy][iz]<0 && iiy<TNY);
		    }

		    if(iiy>=0 && iiy<TNY)
		    {
		      //printf("correcting2 %d %d with %d\n",ix,iy,iiy);
		      for(iv=0;iv<NV+NAVGVARS;iv++)
		      {
			set_uavg(pavg,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iiy,iz));
		      }
                      cell_corrected=1;
		    }
		  
                    //BRANDON - Correcting in theta only does not always work within the horizon.
		    //Correcting in R fixes this problem
                    #ifdef CORRECT_IN_R 
                    if(cell_corrected < 1)
                    {
		      int iix=ix;
		      if(ix>TNX/2) //outer half -- search toward interior
		      {
			do
			  iix--;
			while(flag[iix][iy][iz]<0 && iix>=0);
		      }
		      else //inner half -- search outward
		      {
			do
			 iix++;
			while(flag[iix][iy][iz]<0 && iix<TNX);
		      }
		      
		      if(iix>=0 && iix<TNX)
		      {
			//printf("correcting2 %d %d with %d\n",ix,iy,iiy);
			for(iv=0;iv<NV+NAVGVARS;iv++)
			{
			  set_uavg(pavg,iv,ix,iy,iz,get_uavg(pavg,iv,iix,iy,iz));
			} 
                        cell_corrected=1;
		      }
	            
                    } //if(cell_corrected<1)
                    #endif //CORRECT_IN_R
		    if(cell_corrected<1)
		    {
		      printf("was not able to find good neighbor for %d %d .\n",ix,iy);
		    }

		  } //TNZ==1
		  else //in 3d - in phi
		  {
		    int iiz = iz;
		    do
		      iiz--;
		    while(flag[ix][iy][iiz]<0 && iiz>0);
		    if(iiz<=0)
		    {
		      iiz=iz;
                      do
		       iiz++;
		      while(flag[ix][iy][iiz]<0 && iiz<TNZ);
		    }
		    if(iiz<TNZ && iiz>0)
		    {
		      printf("correcting3D %d %d %d with %d/%d\n",ix,iy,iz,iiz,TNZ);
		      for(iv=0;iv<NV+NAVGVARS;iv++)
		      {
		        set_uavg(pavg,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iy,iiz));
		      }
		    }		    
		    //printf("%d %d %d uncorrected : correcting in 3d not implemented yet.\n",ix,iy,iz); 
		  }
		}//if(flag[ix][iy][iz]<0) 

      } //if(doingavg)
      else //read res file
      {

	  readret=fread_restartfile(ifile,folderin,&t);
	 
	  if(ifile==no1)
            dt=1.;  
	  else
            dt=t-global_time;
         
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

	  //write primitives to pavg array
	  for(iz=0;iz<NZ;iz++)
	    for(iy=0;iy<NY;iy++)
	      for(ix=0;ix<NX;ix++)
		{
		  p2avg(ix,iy,iz,&get_uavg(pavg,0,ix,iy,iz));
		}
	  
      } //if(doingavg)

      //save scalars for single file
      global_time=t;

      #if(SCAOUTPUT==1)
      //if(procotg)
      //{
	//calculates scaleheight etc.
	calc_avgs_throughout();
	  
	//sets bc
	set_bc(t,0);
	  
	//calculate and save scalars
	calc_scalars(scalars,t);

	fprint_scalars(t,scalars,NSCALARS);

      //}
      #endif
      
      //Total average of all files
      add_u_core(1.,pavgtot,dt,pavg,pavgtot,avgarrsize);

      //Advance file loop
      ttot+=dt;
      
  }  //for(ifile=no1_loop; ifile<=no2_loop; ifile+=nostep)

  //Done reading files
  printf("ttot: %f\n",ttot);

  //copy total averaged quantities to pavg
  copy_u_core(1./ttot,pavgtot,pavg,avgarrsize);

  //rewrite primitives to p and run p2u
  for(iz=0;iz<NZ;iz++)
    for(iy=0;iy<NY;iy++)
      for(ix=0;ix<NX;ix++)
      {
	for(iv=0;iv<NV;iv++)
	  set_u(p,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iy,iz));
	struct geometry geom; //ANDREW -- what if pavg was written in BL coordinates?
	fill_geometry(ix,iy,iz,&geom);
	p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
      }
  
  //calculates scaleheight etc.
  calc_avgs_throughout();

  //projects on ghost cells
  set_bc(t,0);

  //ANDREW
  //in  avg vischeat was originally averaged as du, not du/dtau
  //recompute single timestep dt and use that as an estimate of dtau
  #ifdef DIVIDEVISCHEATBYDT
  set_gammagas(0);
  calc_wavespeeds();
  save_timesteps();
  #endif
  
  //Prefix and Suffix for output files
  char prefix[40];
  char suffix[10];
  sprintf(suffix,"");

  if(ifavg)
    sprintf(suffix,"");
  if(ifavg==0)
    sprintf(suffix,"res");
  if(ifphiavg)
    sprintf(suffix,"%sphiavg",suffix);

  //Radial Profiles Output
#if(RADOUTPUT==1)
  sprintf(prefix,"radavg%s%04d-",suffix,no1);
  fprint_radprofiles(t,no2,"analysis",prefix);
#endif

  //Theta Profiles Output
#if(THOUTPUT==1)
  sprintf(prefix,"thavg%s%04d-",suffix,no1);
  fprint_thprofiles(t,no2,"analysis",prefix);
#endif
  
  //Silo file output
#if(SILOOUTPUT==1)
#ifndef NOSILO
  sprintf(prefix,"silavg%s%04d-",suffix,no1);
  fprint_silofile(t,no2,"analysis",prefix);
#endif
#endif

  //Sim file output
#if(SIMOUTPUT!=0)	  
  sprintf(prefix,"simavg%s%04d-",suffix,no1);
  fprint_simplefile(t,no2,"analysis",prefix);
#endif

  //Relativisitc electron spectrum output
#if(RELELSPECTRUMOUTPUT==1)
  int nfout1=no2;
  //ANDREW -- hardcoded cell output -- At equator approx 10M
  int ixx=50;
  int iyy=132;
  int izz=0;  
  fprint_relel_spectrum(t, ixx, iyy, izz, nfout1, folder, "avg_spe",1);
  fprint_relel_avg_spectrum(t, ixx, iyy, izz, nfout1, folder, "avg_spe_cellavg",1);
#endif

#if(SCAOUTPUT==1)
  //if(procotg)
  //{
    fclose(fout_scalars);
  //}
#endif
  return 0;
}

