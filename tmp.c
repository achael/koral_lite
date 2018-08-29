

int
initialize_arrays()
{
  long long i,j,k;

  /********************* Basic arrays ***************************/

  long long Ngrid=(SX)*(SY)*(SZ);
  long long GridSize=Ngrid*sizeof(ldouble);
  long long Nprim=Ngrid*NV;
  long long PrimSize=Nprim*sizeof(ldouble);
  long long Navg=Ngrid*(NV+NAVGVARS);
  long long AvgSize=Navg*sizeof(ldouble);
    
  long long Nmet=(SX)*(SY)*(SZMET)*gSIZE;
  long long MetSize=Nmet*sizeof(ldouble);
  long long Nkris=(SX)*(SY)*(SZMET)*64;
  long long KrisSize=Nkris*sizeof(ldouble);
  
  long long Ntensor=Ngrid*16;
  long long TensorSize=Ntensor*sizeof(ldouble);
  
  //grid
  if((x=(ldouble*)malloc((NX+NY+NZ+6*NG)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
  if((xb=(ldouble*)malloc((NX+1+NY+1+NZ+1+6*NG)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");

  //primitives at cell centers
  if((p=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

  //quantities to average in time
  if((pavg=(ldouble*)malloc(AvgSize))==NULL) my_err("malloc err.\n");
  
  //conserved averages
  if((u=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

  //flags at cell centers
  if((cellflag=(int*)malloc(Ngrid*NFLAGS*sizeof(int)))==NULL) my_err("malloc err.\n");
 
  //metric at cell centers
  if((g=(ldouble*)malloc(MetSize))==NULL) my_err("malloc err.\n");
  if((G=(ldouble*)malloc(MetSize))==NULL) my_err("malloc err.\n");

  //Kristofels at cell centers
  if((gKr=(ldouble*)malloc(KrisSize))==NULL) my_err("malloc err.\n");

  //primitives at cell centers at initial state - used for fixed boundary conditions
  if((pinit=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
 
  //primitives at cell centers at the beginning and end of explicit operator
  if((upreexplicit=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
  if((ppreexplicit=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
   
  //primitives at cell centers at the end of implicit operator
  if((ppostimplicit=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

  //arrays for temporary use (e.g., vector potential, mimic_dynamo)
  if((pproblem1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
  if((ptemp1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
  if((pvecpot=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

  //arrays for avg time
  if((avgselftime=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");  

  //arrays for radiation tensor
#ifdef RADIATION
#if (RADVISCOSITY==SHEARVISCOSITY)
  if((Rijviscprev=(ldouble*)malloc(TensorSize))==NULL) my_err("malloc err.\n");
  if((Rijviscglobal=(ldouble*)malloc(TensorSize))==NULL) my_err("malloc err.\n");
  if((radvisclasttime=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
#endif
#endif

  //arrays for viscous heating
  if((vischeating=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
  if((vischeatingnegebalance=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
  if((vischeatingnegibalance=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
  if((vischeatingtimesdeltae=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");

  //these will aggregate over time so should start with zeros
  for(i=0;i<Ngrid;i++)
    vischeatingnegebalance[i]=vischeatingnegibalance[i]=0.;

  //gamma of gas at the beginnning of timestep
  if((gammagas=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
 
  /****************** extra arrays, used only for time evolution **********************/
  //we might need some of these arrays in postproc (if doingpostproc_avg==1
  #ifdef DIVIDEVISCHEATBYDT
  if(1.)
  #else
  if(doingpostproc==0 || doingpostproc_avg==1)
  #endif
  {     
     //buffer for sending/receiving messages
#ifdef MPI
     if((msgbufs=(ldouble**)malloc(MPIMSGBUFSIZE*sizeof(ldouble*)))==NULL) my_err("malloc err.\n");
     for(i=0;i<MPIMSGBUFSIZE;i++)
       if((msgbufs[i]=(ldouble*)malloc(my_max3(NX*NY*NV*NG,NY*NZ*NV*NG,NZ*NX*NV*NG)*sizeof(ldouble)))==NULL) my_err("malloc err.\n");
#endif

     long long NMetX = (SX+1)*(SY)*(SZMET)*gSIZE;
     long long MetXSize=NMetX*sizeof(ldouble);
     long long NMetY = (SX)*(SY+1)*(SZMET)*gSIZE;
     long long MetYSize=NMetY*sizeof(ldouble);
     long long NMetZ = (SX)*(SY)*(SZMET+1)*gSIZE;
     long long MetZSize=NMetZ*sizeof(ldouble);     
     long long NMetVec = (SX)*(SY)*(SZMET)*16;
     long long MetVecSize=NMetVec*sizeof(ldouble);
            
     //metric at cell x-faces
     if((gbx=(ldouble*)malloc(MetXSize))==NULL) my_err("malloc err.\n");
     if((Gbx=(ldouble*)malloc(MetXSize))==NULL) my_err("malloc err.\n");

     //metric at cell y-faces
     if((gby=(ldouble*)malloc(MetYSize))==NULL) my_err("malloc err.\n");
     if((Gby=(ldouble*)malloc(MetYSize))==NULL) my_err("malloc err.\n");

     //metric at cell z-faces
     if((gbz=(ldouble*)malloc(MetZSize))==NULL) my_err("malloc err.\n");
     if((Gbz=(ldouble*)malloc(MetZSize))==NULL) my_err("malloc err.\n");
      
     //LNRF basis one-forms and vectors
     if((emuup=(ldouble*)malloc(MetVecSize))==NULL) my_err("malloc err.\n");
     if((emulo=(ldouble*)malloc(MetVecSize))==NULL) my_err("malloc err.\n");

     //tetrad one-forms and vectors
     if((tmuup=(ldouble*)malloc(MetVecSize))==NULL) my_err("malloc err.\n");
     if((tmulo=(ldouble*)malloc(MetVecSize))==NULL) my_err("malloc err.\n");

     //Fluxes and wavespeeds
     long long NfluxX = (SX+1)*(SY)*(SZ)*NV;
     long long fluxXSize = NfluxX*sizeof(double);
     long long NfluxY = (SX)*(SY+1)*(SZ)*NV;
     long long fluxYSize = NfluxY*sizeof(double);
     long long NfluxZ = (SX)*(SY)*(SZ+1)*NV;
     long long fluxZSize = NfluxZ*sizeof(double);

     //wavespeeds hd and rad - max(al,ar)
     if((ahdx=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdy=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdz=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradx=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((arady=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradz=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     
     //wavespeeds hd and rad - leftgoing
     if((ahdxl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdyl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdzl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradxl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradyl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradzl=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
  
     //wavespeeds hd and rad - rightgoing
     if((ahdxr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdyr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((ahdzr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradxr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradyr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((aradzr=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
       
     //left-interpolated primitives at cell faces
     if((pbLx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((pbLy=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((pbLz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //right-interpolated primitives at cell faces
     if((pbRx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((pbRy=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((pbRz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //corrected flux at faces
     if((flbx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((flby=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((flbz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //flux based on left-interpolated conserved at cell faces
     if((flLx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((flLy=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((flLz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //flux based on right-interpolated conserved at cell faces
     if((flRx=(ldouble*)malloc(fluxXSize))==NULL) my_err("malloc err.\n");
     if((flRy=(ldouble*)malloc(fluxYSize))==NULL) my_err("malloc err.\n");
     if((flRz=(ldouble*)malloc(fluxZSize))==NULL) my_err("malloc err.\n");

     //auxiliary primitive arrays
     if((pproblem2=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ptm1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ut0=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ut1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ut2=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((ut3=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((dut0=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((dut1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((dut2=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((drt0=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((drt1=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((drt2=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((uforget=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((u_bak_fixup=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");
     if((p_bak_fixup=(ldouble*)malloc(PrimSize))==NULL) my_err("malloc err.\n");

     //timesteps required by each cell
     if((cell_tstepden=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");
     if((cell_dt=(ldouble*)malloc(GridSize))==NULL) my_err("malloc err.\n");

#ifdef MAGNFIELD
     //electromotive force at corners
     long long Nemf = (NX+1)*(NY+1)*(NZ+1)*3;
     long long EmfSize = Nemf*sizeof(ldouble);
     if((emf=(ldouble*)malloc(EmfSize))==NULL) my_err("malloc err.\n");
#endif
  
  }

  init_all_kappa_table();

  return 0;
}

