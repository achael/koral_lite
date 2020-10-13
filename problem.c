/*! \file problem.c
 \brief problem-related but not problem-specific routines
 */

 #include "ko.h"

/******************************************************/
/***************** time integration ********************/
/******************************************************/
/*! \fn int solve_the_problem(ldouble tstart, char* folder)
 \brief Outer routine that contains the main loop over time
 
 \param[in] tstart Start time of the simulation
 \param[in] folder Name of folder for output files ("dumps")
*/

int
solve_the_problem(ldouble tstart, char* folder)
{
  ldouble t = tstart, t1 = TMAX;
  ldouble totalmass=0.;

#ifdef DTOUT_LOG
  ldouble dtout;
  if(tstart==0.)
    dtout = pow(10,DTOUT1_LOG_INIT);
  else
    dtout = pow(10,floor(log10(tstart)+1));
  printf("dtout 1 %e\n",dtout);
#else  
  ldouble dtout = DTOUT1;
#endif
  ldouble dtoutavg = DTOUT2;
  
  ldouble dtsaveavg;
  #ifndef DTAVG
  dtsaveavg = 1;
  #else
  dtsaveavg = DTAVG;
  #endif
  
  ldouble dtsource, taim;
  ldouble fprintf_time = 0.;
  ldouble lasttout_floor=floor(t/dtout); 
  ldouble lasttoutavg_floor=floor(t/dtoutavg);
  ldouble lastsaveavg_floor=floor(t/dtsaveavg);
  
  int i1=0.,i2=0.;
  int fprintf_nstep=0;  
  int i,j,ii;
  int ix,iy,iz,iv;
  int loopsallociter;
  int spitoutput,lastzone;
  int nentr[8],nentr2[8];
  
  struct timespec temp_clock;

  // Initial timestep guess
  dt=-1.;
  max_ws[0]=max_ws[1]=max_ws[2]=10000.;
  if(NZ>1)
    tstepdenmax=max_ws[0]/min_dx + max_ws[1]/min_dy + max_ws[2]/min_dz;
  else if(NY>1)
    tstepdenmax=max_ws[0]/min_dx + max_ws[1]/min_dy;
  else
    tstepdenmax=max_ws[0]/min_dx;
  tstepdenmax/=TSTEPLIM;
  tstepdenmin=tstepdenmax;

  // Calculate and set consistent gamma over domain + ghost cells + corners
  set_gammagas(0);

  // Set initial timestep info in all cells
  for(ii=0;ii<Nloop_0;ii++) 
  {
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2];
    
    set_u_scalar(cell_tstepden,ix,iy,iz,tstepdenmax);
    set_u_scalar(cell_dt,ix,iy,iz,1./tstepdenmax);
  }

  // Choose the smallest timestep 
  mpi_synchtiming(&t);
  
  /***********************************************************************/  
  /***********************************************************************/
  // START OF MAIN SIMULATION TIME LOOP 
  /***********************************************************************/
  /***********************************************************************/
  
  nstep=0;
  while (t < t1 && nfout1<=NOUTSTOP && nstep<NSTEPSTOP)
  {   
      
      global_int_slot[GLOBALINTSLOT_NTOTALRADIMPFIXUPS]=0; //counting number of critical failures
      global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS]=0;    //counting mhd fixups
      global_int_slot[GLOBALINTSLOT_NTOTALRADFIXUPS]=0;    //counding rad fixups
 
      spitoutput=0;
      global_time=t;
      nstep++;
      
      //choose the smallest timestep etc.
      mpi_synchtiming(&t);

      //initial time mark
      my_clock_gettime(&temp_clock);
      start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
      ldouble tstepden;

      // dt is based on the estimate from the last timestep
      dt=1./tstepdenmax;
      global_dt=dt;
      
      if(t+dt>t1)
      {
	dt=t1-t;
      }
   
      // Resetting wavespeeds
      tstepdenmax=-1.;
      tstepdenmin=BIG;

      // Calculate viscosity tensor
      #if (RADVISCOSITY==SHEARVISCOSITY)
      calc_Rij_visc_total();
      #endif

      //**********************************************************************
      // Take a step in time according to the specified TIMESTEPPING SCHEME
      //**********************************************************************

      if(TIMESTEPPING==-100) //skip evolution completely
      {
	    save_timesteps(); 
	    t+=dt;
      }

      /************************** RK2IMEX **********************************/      
      else if(TIMESTEPPING==RK2IMEX)
      {
	    
	    // User defined finger
	    my_finger(global_time);

	    ldouble gamma=1.-1./sqrt(2.);
	    ldouble dtcell;

	    save_timesteps(); 

            // Calculate and set consistent gamma over domain + ghost cells + corners
	    set_gammagas(0);

	    dtcell=dt;

	    /******* 1st implicit **********/
            // Set ut0 = u over domain + ghost cells
	    copyi_u(1.,u,ut0);
        
            // Set ptm1 = p over domain
            copy_u(1.,p,ptm1);
        
            // Implicit evolution of radiation terms 
	    op_implicit(t, dt*gamma); //U(0) in *ut0;  U(1) in *u
        
	    global_impdt=dt*gamma;

            // Set ppostimplicit = p over domain + ghost cells
            copyi_u(1.,p,ppostimplicit);

	    // Calculate 1st implicit deriv
	    // R(U(1)) in *drt1;`
            #pragma omp parallel for private(ii,ix,iy,iz,iv) firstprivate(dtcell)
	    for(ii=0;ii<Nloop_0;ii++)
            {
              ix=loop_0[ii][0];
              iy=loop_0[ii][1];
              iz=loop_0[ii][2];

	      PLOOP(iv) set_u(drt1,iv,ix,iy,iz,(1./(dtcell*gamma))*get_u(u,iv,ix,iy,iz)+(-1./(dtcell*gamma))*get_u(ut0,iv,ix,iy,iz)); 
	    } 



	    /******* 1st explicit **********/

	    // Set ut1 = u over domain + ghost cells
	    copyi_u(1.,u,ut1);

	    // ANDREW is this excessive? should be consistent after implicit? 
	    calc_u2p(0,1);
        
            // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[0],&nentr2[0]);
        
            // Special treatment near axis (or inner surface)
	    do_correct();
	    
	    // Exchange MPI data for boundaries
	    mpi_exchangedata();
        
            // Set boundary conditions on conserveds in the ghost cells
	    set_bc(t,0);
	    
            // Explicit evolution (advection plus source terms) from t to t+dt
	    op_explicit (t, dt);  //U(1) in *ut1;

	    
            // Artifical dynamo (ifdef MIMICDYNAMO)
	    apply_dynamo(t,dt);
        
            // Intermediate step between explicit and implicit for relativistic electrons
	    op_intermediate (t,dt);
	    
	    global_expdt=dt;
 
            // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[1],&nentr2[1]);
	    copy_entropycount();

	    // Calculate 1st explicit deriv
	    //F(U(1)) in *dut1;
            #pragma omp parallel for private(ii,ix,iy,iz,iv) firstprivate(dtcell)
	    for(ii=0;ii<Nloop_0;ii++)
            {
              ix=loop_0[ii][0];
              iy=loop_0[ii][1];
              iz=loop_0[ii][2];

	      PLOOP(iv) set_u(dut1,iv,ix,iy,iz,(1./(dtcell))*get_u(u,iv,ix,iy,iz)+(-1./(dtcell))*get_u(ut1,iv,ix,iy,iz)); 
	    }
	   
	    /******* 1st together **********/
	    
            //(U(0) + dt F(U(1)) + dt (1-2gamma) R(U(1))) in *u
            #pragma omp parallel for private(ii,ix,iy,iz,iv) firstprivate(dtcell)
	    for(ii=0;ii<Nloop_0;ii++)
            {
              ix=loop_0[ii][0];
              iy=loop_0[ii][1];
              iz=loop_0[ii][2];

	      PLOOP(iv) set_u(u,iv,ix,iy,iz,get_u(ut0,iv,ix,iy,iz)+(dtcell)*get_u(dut1,iv,ix,iy,iz)+(dtcell*(1.-2.*gamma))*get_u(drt1,iv,ix,iy,iz)); 
	    }
	    
	    /******* 2nd implicit **********/
            // Set uforget = u over domain + ghost cells
	    copyi_u(1.,u,uforget);

	    // Invert to primitives
	    calc_u2p(0,1);
        
            // Set ptm1 = p over domain
            copy_u(1.,p,ptm1); 

#pragma omp barrier
            // Special treatment near axis (or inner surface)
	    do_correct();
#pragma omp barrier

            // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[2],&nentr2[2]);
	    
            // Implicit evolution of radiation terms
	    op_implicit (t,gamma*dt); //U(2) in *u

	    global_impdt=gamma*dt;
        
            // Set ppostimplicit = p over domain + ghost cells
	    copyi_u(1.,p,ppostimplicit);
        
            // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[3],&nentr2[3]);

	    #if(AVGOUTPUT>0) // Save to avg arrays
            #ifdef DTAVG //Dont save every step
            if(lastsaveavg_floor!=floor(t/dtsaveavg))
	    {
	      if(nstep>1)
	      {
		save_avg(dt);
                lastsaveavg_floor=floor(t/dtsaveavg);	 
              }
	    }
	    #else
	    if(nstep>1) save_avg(dt);
            #endif 
            #endif

	    // Calculate 2nd implicit deriv
	    // R(U(2)) in *drt2;
            #pragma omp parallel for private(ii,ix,iy,iz,iv) firstprivate(dtcell)
	    for(ii=0;ii<Nloop_0;ii++)
            {
              ix=loop_0[ii][0];
              iy=loop_0[ii][1];
              iz=loop_0[ii][2];

	      PLOOP(iv) set_u(drt2,iv,ix,iy,iz,(1./(dtcell*gamma))*get_u(u,iv,ix,iy,iz)+(-1./(dtcell*gamma))*get_u(uforget,iv,ix,iy,iz)); 
	    }
	  
	    /******* 2nd explicit **********/
            // Set ut2 = u over domain + ghost cells
	    copyi_u(1.,u,ut2);

	    // Invert to primitives
	    // ANDREW  is this excessive? should still be consistent after implicit!
	    calc_u2p(0,1);

            // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
            count_entropy(&nentr[4],&nentr2[4]);

            // Special treatment near axis (or inner surface)
	    do_correct();

	    // Exchange MPI data for boundaries
	    mpi_exchangedata();
        
            // Set boundary conditions on conserveds in the ghost cells
	    set_bc(t,0);

            // Explicit evolution (advection plus source terms) from t to t+dt
	    op_explicit (t,dt); //U(2) in *ut2;

            // Artifical dynamo (ifdef MIMICDYNAMO)
	    apply_dynamo(t,dt);
	    
            // Intermediate step between explicit and implicit for relativistic electrons
	    op_intermediate (t,dt);
        
	    global_expdt=dt;
	    
            // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[5],&nentr2[5]);

	    // Calculate 2nd explicit deriv
	    // F(U(2)) in *dut2;
            #pragma omp parallel for private(ii,ix,iy,iz,iv) firstprivate(dtcell)
	    for(ii=0;ii<Nloop_0;ii++)
            {
              ix=loop_0[ii][0];
              iy=loop_0[ii][1];
              iz=loop_0[ii][2];

	      PLOOP(iv) set_u(dut2,iv,ix,iy,iz,(1./(dtcell))*get_u(u,iv,ix,iy,iz)+(-1./(dtcell))*get_u(ut2,iv,ix,iy,iz)); 
	    }
	 
            /******* explicit together **********/
            //U(0) + dt/2 (F(U(1)) + F(U(2))) in *u
            #pragma omp parallel for private(ii,ix,iy,iz,iv) firstprivate(dtcell)
	    for(ii=0;ii<Nloop_0;ii++)
            {
              ix=loop_0[ii][0];
              iy=loop_0[ii][1];
              iz=loop_0[ii][2];

	      PLOOP(iv) set_u(u,iv,ix,iy,iz,get_u(ut0,iv,ix,iy,iz)+(dtcell/2.)*get_u(dut1,iv,ix,iy,iz)+(dtcell/2.)*get_u(dut2,iv,ix,iy,iz)); 
	    }
	    
	    /******* implicit together ***********/
	    //ANDREW why should this be a separate loop?
	    //u += dt/2 (R(U(1)) + R(U(2))) in *u
            #pragma omp parallel for private(ii,ix,iy,iz,iv) firstprivate(dtcell)
	    for(ii=0;ii<Nloop_0;ii++)
            {
              ix=loop_0[ii][0];
              iy=loop_0[ii][1];
              iz=loop_0[ii][2];

	      PLOOP(iv) set_u(u,iv,ix,iy,iz,get_u(u,iv,ix,iy,iz)+(dtcell/2.)*get_u(drt1,iv,ix,iy,iz)+(dtcell/2.)*get_u(drt2,iv,ix,iy,iz));
	    }
	   
	    // Final inversion
	    calc_u2p(0,1);
	    
	    // Heat species at end
            #ifdef HEATELECTRONSATENDRK2
	    heat_electronions_with_state(dt);
            #endif

            // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[6],&nentr2[6]);

	    // ANDREW should the next 3 steps be moved to start of loop?  
            // Special treatment near axis (or inner surface)
	    do_correct();

	    // Exchange MPI data for boundary cells
	    mpi_exchangedata();
        
            // Set boundary conditions on conserveds in the ghost cells
	    set_bc(t,0);
	
	    t+=dt;	 

            // Compute new entropy from rho and uint over the domain
	    update_entropy();

      }  // else if(TIMESTEPPING==RK2IMEX)

      /******************************* RK2 Heun **********************************/      
      else if(TIMESTEPPING==RK2HEUN)
      {
	    //User defined finger
	    my_finger(global_time);
 
	    save_timesteps();

            /********** 1st **********/	    	    
	    // Set ut0 = u over domain + ghost cells
	    copyi_u(1.,u,ut0);
        
	    // Set ptm1 = p over domain
	    copy_u(1.,p,ptm1);

	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[0],&nentr2[0]);
 
	    // Special treatment near axis (or inner surface)
	    do_correct();

	    // Exchange mpi boundary cells
	    mpi_exchangedata();
        
	    // Set boundary conditions on conserveds in the ghost cells
	    set_bc(t,0);
        
	    // Calculate and set consistent gamma over domain + ghost cells + corners
	    set_gammagas(0);

	    // Explicit evolution (advection plus source terms) from t to t+dt
	    global_expdt=dt;
	    op_explicit (t, dt);

	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[1],&nentr2[1]);
	    copy_entropycount();

	    // Artifical dynamo (ifdef MIMICDYNAMO)
	    apply_dynamo(t, dt);
        
	    // Intermediate step between explicit and implicit for relativistic electrons
	    op_intermediate(t, dt);

#ifdef RADIATION
	    // Special treatment near axis (or inner surface)
	    do_correct();
#endif
	    
	    // Implicit evolution of radiation terms
	    global_impdt=dt;
	    op_implicit (t, dt);
        
	    // Set ppostimplicit = p over domain + ghost cells
	    copyi_u(1.,p,ppostimplicit);
        
	    addi_u(1.,u,-1.,ut0,ut2); // dt*R(U(1)) 

            /*********** 2nd  **********/

	    // Set ut1 = u over domain + ghost cells
	    copyi_u(1.,u,ut1);

	    // ANDREW is this excessive?  Should be consistent after implicit!
	    calc_u2p(0,1);
        
	    // Set ptm1 = p over domain
	    copy_u(1.,p,ptm1); 
	    
	    // Special treatment near axis (or inner surface)
	    do_correct();
        
	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[2],&nentr2[2]);

	    // Exchange mpi data for boundaries
	    mpi_exchangedata();
        
	    // Set boundary conditions on conserveds in the ghost cells
	    set_bc(t,0);
	    
	    // Explicit evolution (advection plus source terms) from t to t+dt
	    global_expdt=dt;
	    op_explicit (t,dt);
        
	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[3],&nentr2[3]);
        
	    // Artifical dynamo (ifdef MIMICDYNAMO)
	    apply_dynamo(t,dt);
        
	    // Intermediate step between explicit and implicit for relativistic electrons
	    op_intermediate(t, dt);

#ifdef RADIATION	    
	    // Special treatment near axis (or inner surface)
	    do_correct();
#endif
	    
	    // Implicit evolution of radiation terms
	    op_implicit (t,dt);        
	    global_impdt=dt;
        
	    // Set ppostimplicit = p over domain + ghost cells
	    copyi_u(1.,p,ppostimplicit);

	    #if(AVGOUTPUT>0) // Save to avg arrays
            #ifdef DTAVG // Don't average every step
            if(lastsaveavg_floor!=floor(t/dtsaveavg))
	    {
	      if(nstep>1)
	      {
		save_avg(dt);
                lastsaveavg_floor=floor(t/dtsaveavg);	 
              }
	    }
	    #else
	    if(nstep>1) save_avg(dt);
            #endif 
            #endif

	    // Together      
	    addi_u(1.,u,-1.,ut1,ut3);  // dt*R(U(2))
	    addi_u_3(1.,ut0,1./2.,ut2,1./2.,ut3,u); //u = U(0) + dt/2 (R(U(1)) + R(U(2))) in *u

	    // Calculate primitves
	    calc_u2p(0,1); //do not calculate visc. heating, do count entropy inversions
	    
	    // Heat species at end
            #ifdef HEATELECTRONSATENDRK2
	    heat_electronions_with_state(dt);
            #endif

	    // Update entropy from rho and uint over the domain
	    update_entropy();
	    t+=dt;
      }  // else if(TIMESTEPPING==RK2HEUN)

      /************************** RK2 **********************************/      
      else if(TIMESTEPPING==RK2)
      { 
	    // User defined  finger
	    my_finger(global_time);
 
	    save_timesteps();
	 
	    /************ 1st ************/	 
	    // Set ut0 = u over domain + ghost cells
	    copyi_u(1.,u,ut0);
        
	    // Set ptm1 = p over domain
	    copy_u(1.,p,ptm1); 

	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[0],&nentr2[0]);
 
	    // Special treatment near axis (or inner surface)
	    do_correct();

	    // Exchange MPI data for boundaries
	    mpi_exchangedata();
        
	    // Set boundary conditions on conserveds in the ghost cells
	    set_bc(t,0);
        
	    // Calculate and set consistent gamma over domain + ghost cells + corners
	    set_gammagas(0);

	    // Explicit evolution (advection plus source terms) from t to t + 0.5*dt
	    global_expdt=.5*dt;
	    op_explicit (t,.5*dt);

	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[1],&nentr2[1]);
	    copy_entropycount();

	    // Artifical dynamo (ifdef MIMICDYNAMO)
	    apply_dynamo(t,.5*dt);
        
	    // Intermediate step between explicit and implicit for relativistic electrons
	    op_intermediate(t, .5*dt);

#ifdef RADIATION
	    // Special treatment near axis (or inner surface)
	    do_correct();
#endif

	    // Implicit evolution of radiation terms
	    global_impdt=.5*dt;
	    op_implicit (t,.5*dt);
        
	    // Set ppostimplicit = p over domain + ghost cells
	    copyi_u(1.,p,ppostimplicit);
        
	    addi_u(1.,u,-1.,ut0,ut2); // dt*R(U(1))
	                             

	    /************ 2nd ************/	 
	    // Set ut1 = u over domain + ghost cells
	    copyi_u(1.,u,ut1);

	    // ANDREW is this excessive? Should be consistent after implicit!
	    calc_u2p(0,1);
        
	    // Set ptm1 = p over domain
	    copy_u(1.,p,ptm1);

	    // Special treatment near axis (or inner surface)
	    do_correct();
        
	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[2],&nentr2[2]);

	    // Exchange MPI data for boundaries
	    mpi_exchangedata();
        
	    // Set boundary conditions on conserveds in the ghost cells
	    set_bc(t,0);

	    // Explicit evolution (advection plus source terms) from t to t+dt
	    global_expdt=dt;
	    op_explicit (t,dt);
        
	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[3],&nentr2[3]);
        
	    // Artifical dynamo (ifdef MIMICDYNAMO)
	    apply_dynamo(t,dt);
        
	    // Intermediate step between explicit and implicit for relativistic electrons
	    op_intermediate(t,dt);

#ifdef RADIATION	    
	    // Special treatment near axis (or inner surface)
	    do_correct();
#endif
	    
	    // Implicit evolution of radiation terms
	    op_implicit (t,dt);
	    global_impdt=dt;
        
	    // Set ppostimplicit = p over domain + ghost cells
	    copyi_u(1.,p,ppostimplicit);

	    #if(AVGOUTPUT>0) // Save to avg arrays
            #ifdef DTAVG //Don't average every step
            if(lastsaveavg_floor!=floor(t/dtsaveavg))
	    {
	      if(nstep>1)
	      {
		save_avg(dt);
                lastsaveavg_floor=floor(t/dtsaveavg);	 
              }
	    }
	    #else
	    if(nstep>1) save_avg(dt);
            #endif 
            #endif
	    
	    // Together      
	    addi_u(1.,u,-1.,ut1,ut3);  // dt * (R(U(2)))
	    addi_u_3(1.,ut0,0.,ut2,1.,ut3,u); //U(0) + dt R(U(2)) in *u

	    // Calculate primitves
	    calc_u2p(0,1); //do not calculate visc. heating, do count entropy inversions
	    
	    // Heat species at end
            #ifdef HEATELECTRONSATENDRK2
	    heat_electronions_with_state(dt);
            #endif

	    // Compute entropy from rho and uint over the domain
	    update_entropy();

	    t+=dt;
      }  // else if(TIMESTEPPING==RK2)

      /************************** RK1 **********************************/      
      else if(TIMESTEPPING==RK1) //only for testing!!!
      { 
            // User defined finger
	    my_finger(global_time);
	    save_timesteps();
	 
	    // Set ut0 = u over domain + ghost cells
	    copyi_u(1.,u,ut0);
        
	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[0],&nentr2[0]);
 
	    // Special treatment near axis (or inner surface)
	    do_correct();

	    // Exchange MPI boundary data
	    mpi_exchangedata();
        
	    // Set boundary conditions on conserveds in the ghost cells
	    set_bc(t,0);
        
	    // Calculate and set consistent gamma over domain + ghost cells + corners
	    set_gammagas(0);
	    
	    // Explicit evolution (advection plus source terms) from t to t+dt
	    global_expdt=dt;
	    op_explicit (t,1.*dt); 

	    // Count number of entropy inversions: ENTROPYFLAG, ENTROPYFLAG2
	    count_entropy(&nentr[1],&nentr2[1]); copy_entropycount();

	    // Artifical dynamo (ifdef MIMICDYNAMO)
	    apply_dynamo(t,dt);
        
	    // Intermediate step between explicit and implicit for relativistic electrons
	    op_intermediate(t, dt);

#ifdef RADIATION
	    // Special treatment near axis (or inner surface)
	    do_correct();
#endif
        
	    // Implicit evolution of radiation terms
	    global_impdt=dt;
	    op_implicit (t,1.*dt);
        
	    // Set ppostimplicit = p over domain + ghost cells
	    copyi_u(1.,p,ppostimplicit);

            #if(AVGOUTPUT>0) // Save to avg arrays
            #ifdef DTAVG //Don't average every step
            if(lastsaveavg_floor!=floor(t/dtsaveavg))
	    {
	      if(nstep>1)
	      {
		save_avg(dt);
                lastsaveavg_floor=floor(t/dtsaveavg);
	      }
	    }
	    #else
	    if(nstep>1) save_avg(dt);
            #endif 
            #endif
	    
	    // Calculate primitves
	    // ANDREW is this excessive? Should be consistent after implicit!
	    calc_u2p(0,1); //do not calculate visc. heating, do count entropy inversions
	    
	    // Heat species at end
            #ifdef HEATELECTRONSATENDRK2
	    heat_electronions_with_state(dt);
            #endif

	    // Update to new time: t+dt
	    t+=dt;
      }  // else if(TIMESTEPPING==RK1)
      
      else  // TIMESTEPPING ERROR
      {
	  my_err("wrong time stepping specified\n");
      }
       
      //**********************************************************************
      //************************* outputs ************************************
      //**********************************************************************

      //for outputs - use what came out of 2nd implicit: //ANDREW why??
      
      // Set uforget = p and p = ppostimplicit over domain
      copy_u(1.,p,uforget); //backup current primitives
      copy_u(1.,ppostimplicit,p);
   
      //counting faiures and average parameters of the implicit solver
      int nfailures[3],nfailuresloc[3]={global_int_slot[GLOBALINTSLOT_NTOTALRADIMPFAILURES],
					global_int_slot[GLOBALINTSLOT_NTOTALMHDFIXUPS],
					global_int_slot[GLOBALINTSLOT_NTOTALRADFIXUPS]};
      
#ifdef MPI
      MPI_Allreduce(nfailuresloc, nfailures, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);      
#else
      for(i=0;i<3;i++) nfailures[i]=nfailuresloc[i];
#endif
      
      //quit if we have exceeded maxiumum number of failures
      int maxfailures=TNX*TNY*TNZ/1;
      if(nfailures[0]>maxfailures || nfailures[1]>maxfailures || nfailures[2]>maxfailures)
      {
	  printf("exceeded # of failures (%d %d %d) - exiting.\n",
		 nfailures[0],nfailures[1],nfailures[2]);
	  exit(-1);
      }
      
#ifdef RADIATION
      //get average number of iterations in the implicit solver
      ldouble avimpitloc[5],avimpit[5];
      int impnumsloc[7],impnums[7];
      
      avimpitloc[0]=global_int_slot[GLOBALINTSLOT_NIMPENERMHD]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENERMHD]
	/((ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERMHD]+(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERMHDFF]);
      avimpitloc[1]=global_int_slot[GLOBALINTSLOT_NIMPENERRAD]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENERRAD]
	/((ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERRAD]+(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENERRADFF]);
      avimpitloc[2]=global_int_slot[GLOBALINTSLOT_NIMPENTRMHD]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENTRMHD]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENTRMHD];
      avimpitloc[3]=global_int_slot[GLOBALINTSLOT_NIMPENTRRAD]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPENTRRAD]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPENTRRAD];
      avimpitloc[4]=global_int_slot[GLOBALINTSLOT_NIMPLTE]==0 ? 0. : 
	(ldouble)global_int_slot[GLOBALINTSLOT_ITERIMPLTE]/(ldouble)global_int_slot[GLOBALINTSLOT_NIMPLTE];
            
      impnumsloc[0]=global_int_slot[GLOBALINTSLOT_NIMPENERMHD];
      impnumsloc[1]=global_int_slot[GLOBALINTSLOT_NIMPENERRAD];
      impnumsloc[2]=global_int_slot[GLOBALINTSLOT_NIMPENERMHDFF];
      impnumsloc[3]=global_int_slot[GLOBALINTSLOT_NIMPENERRADFF];
      impnumsloc[4]=global_int_slot[GLOBALINTSLOT_NIMPENTRMHD];
      impnumsloc[5]=global_int_slot[GLOBALINTSLOT_NIMPENTRRAD];
      impnumsloc[6]=0;

#ifdef MPI
      MPI_Allreduce(impnumsloc, impnums, 7, MPI_INT, MPI_SUM,
		    MPI_COMM_WORLD);  
      MPI_Allreduce(avimpitloc, avimpit, 5, MPI_LDOUBLE, MPI_MAX,
                   MPI_COMM_WORLD);  
#else
      for(i=0;i<5;i++) avimpit[i]=avimpitloc[i];
      for(i=0;i<7;i++) impnums[i]=impnumsloc[i];
#endif 
#endif  
      
      //time mark
      my_clock_gettime(&temp_clock);    
      end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

      //performance
      ldouble znps=TNX*TNY*TNZ/(end_time-start_time);

      //GM/c3 per day
      ldouble tgpd = dt/(end_time-start_time)*3600.*24.;

      //save avg files	  
      //avg goes first so that what is later can use it
#if(AVGOUTPUT>0) 
      if(lasttoutavg_floor!=floor(t/dtoutavg))
      {
	  my_clock_gettime(&temp_clock);    
	  start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
	  if(PROCID==0)
	    printf("%d > avg file no #%6d dumped ",PROCID,nfout2);

	  //save avg array to file
	  for(iz=0;iz<NZ;iz++)
	    for(iy=0;iy<NY;iy++)
	      for(ix=0;ix<NX;ix++)
	      {
		  for(iv=0;iv<(NV+NAVGVARS);iv++)
		    set_uavg(pavg,iv,ix,iy,iz,get_uavg(pavg,iv,ix,iy,iz) / get_u_scalar(avgselftime,ix,iy,iz));
		  set_u_scalar(avgselftime,ix,iy,iz,0.);
	      }              
	  
	  //save avg file
	  fprint_avgfile(t,folder,"avg");
      
	  //zero out avg values over domain

    long long Navg = (long long) SX*SY*SZ*(NV+NAVGVARS); // RN: Apr 2, 2019
	  copy_u_core(0., pavg, pavg, Navg);
	  avgtime=0.;
	  
	  nfout2++;

	  lasttoutavg_floor=floor(t/dtoutavg);	 

	  my_clock_gettime(&temp_clock);    
	  end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
	  if(PROCID==0)
	    printf( "-- I/O took %.3f seconds\n", end_time-start_time);

      }  // if(lasttoutavg_floor!=floor(t/dtoutavg))
#endif  // if(AVGOUTPUT>0)


      //save snapshot files
#ifdef DTOUT_LOG
      if(lasttout_floor!=floor(t/dtout))
#else
      if(lasttout_floor!=floor(t/dtout) || ALLSTEPSOUTPUT || t>.9999999*t1 || spitoutput==1)
#endif
      {
	  my_clock_gettime(&temp_clock);    
	  start_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
	  if(PROCID==0)
	    printf("%d > snap file no #%6d dumped at t=%12.5e\n",PROCID,nfout1,t);
	  
	  // Set boundary conditions on conserveds in the ghost cells
	  set_bc(t,0);

	  //print restart file
	  fprint_restartfile(t,folder);

	  //dump on-the-go dump files if not MPI
#ifndef MPI

          #if(SCAOUTPUT==1) // scalar  dumpfiles
	  fprint_scalars(t,scalars,NSCALARS);
          #endif
      
          #if(RADOUTPUT==1) //radial  profiles
	  fprint_radprofiles(t,nfout1,folder,"rad");
          #endif
	  
          #if(THOUTPUT==1) //theta profiles
	  fprint_thprofiles(t,nfout1,folder,"th");
          #endif
      
          #if(SILOOUTPUT==1) //silo files
          #ifndef NOSILO
	  fprint_silofile(t,nfout1,folder,"sil");
          #endif
          #endif 
      
          #if(SIMOUTPUT!=0) //sim files
	  fprint_simplefile(tstart,nfout1,folder,"sim");
          #endif
      
	  #if(RELELSPECTRUMOUTPUT==1) //nonthermal spectrum
          fprint_relel_spectrum(t,NTH_SPEC_IX,NTH_SPEC_IY,NTH_SPEC_IZ,nfout1,folder,"spe",0);
          #endif

#endif  // ifndef MPI


	  nfout1++;

	  
	  #ifdef DTOUT_LOG
          dtout = pow(10, DTOUT1_LOG_INIT + DTOUT_LOG*(nfout1-1));
	  printf("dtout %d  %e\n",nfout1, dtout);
          #endif
	  
	  lasttout_floor=floor(t/dtout);	 

	  my_clock_gettime(&temp_clock);    
	  end_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;
	  if(PROCID==0)
	    printf( "-- I/O took %.3f seconds\n", end_time-start_time);

      }  // if(lasttout_floor!=floor(t/dtout) || ALLSTEPSOUTPUT || t>.9999999*t1 || spitoutput==1)
      
      //print performance to screen only every second
      int printeacht=0;
      #ifdef PRINTEACHT
      printeacht = 1;
      #endif
      if((end_time-fprintf_time>1.  || printeacht) && PROCID==0)
      {
	  
	  printf("st #%6d t=%12.5e dt=%.2e mpi=%3.1f znps=%.0f tgpd=%.2e fail# %1d %1d %1d "
		 ,nstep,t,dt,2.*maxmp_time/(end_time-start_time),znps,tgpd,
		  nfailures[0],nfailures[1],nfailures[2]);
#ifdef RADIATION
#ifndef SKIPRADSOURCE
	  printf("imp# %d %d %d %d %d %d %d | %.1f %.1f %.1f %.1f %.1f ",
		  impnums[0],impnums[1],impnums[2],impnums[3],impnums[4],impnums[5],impnums[6],
		  avimpit[0],avimpit[1],avimpit[2],avimpit[3],avimpit[4]);
#endif
#endif

	  printf("\n");
	  fflush(stdout);

	  fprintf_time = end_time;
	  fprintf_nstep = nstep;
      }  // if((end_time-fprintf_time>1.  || printeacht) && PROCID==0)
      
      // restore original final primitives after outputs
      // Set p = uforget over domain
      copy_u(1.,uforget,p);

  } // while (t < t1 && nfout1<=NOUTSTOP && nstep<NSTEPSTOP)

  /***********************************************************************/  
  /***********************************************************************/
  // END OF MAIN SIMULATION TIME LOOP 
  /***********************************************************************/
  /***********************************************************************/
  return 0;
}  

/*********************************************/
/*! \fn int my_finger(ldouble t)
 \brief "stirs" the simulation as defined in input file in PROBLEMS folder
 
 \param[in] t Global time
 */
/*********************************************/

int my_finger(ldouble t)
{

#ifdef PR_FINGER
#include PR_FINGER
#endif

  return 0;
} 

/*********************************************/
/* suplementary routines */
/* uses PROBLEMS/XXX/tools.c */   
/*********************************************/
#ifdef PR_TOOLS
#include PR_TOOLS
#endif


