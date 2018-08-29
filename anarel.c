//KORAL - anarel.c
//loop over avg files and do ??? -- ANDREW: unclear what is saved here
//loop over res files and save averaged profiles if ANARELRADOUTPUT==1

#include "ko.h"

int 
main(int argc, char **argv)
{  
  #ifdef MPI
  printf("anarel works on one core only, do not use MPI, please\n");
  exit(-1);
  #endif

  long long i,j;

  //which files to read
  int no1res,no2res,nostepres;
  int no1avg,no2avg,nostepavg;
  if(argc!=7)
  {
    printf("Not enough input arguments.\n");
    printf("./anarel (no1 no2 nostep)_avg (no1 no2 nostep)_res\n");
    return -1;
  }
  else
  {      
    no1avg=atof(argv[1]);
    no2avg=atof(argv[2]);
    nostepavg=atof(argv[3]);
    no1res=atof(argv[4]);
    no2res=atof(argv[5]);
    nostepres=atof(argv[6]);
  }

  //currently gsl is not used
  gsl_set_error_handler_off();
  
  //random number gen. initialization
  srand ( time(NULL) );

  //preparing arrays
  initialize_arrays();

  //sets the grid
  set_grid(&min_dx,&min_dy,&min_dz,&max_dt);

  //precalculates metric etc.
  calc_metric();

  //print coordinate file
  #ifdef COORDOUTPUT
  fprint_coordfile("analysis","coord");
  #endif

  //folder to write in
  char folder[100],bufor[100];
  sprintf(folder,"analysis");

  /*********************/
  //loop over avg files first
  /*********************/
  long long ncells = (SX)*(SY)*(SZ);
  long long avgarrsize=ncells*(NV+NAVGVARS);
  long long avgarrsizeMalloc=avgarrsize*sizeof(ldouble);

  ldouble *pavgtot=(ldouble*)malloc(avgarrsizeMalloc);
  for(i=0;i<avgarrsize;i++)
    pavg[i]=pavgtot[i]=0.;

  printf("using avg files #%04d to #%04d with %d step \n",no1avg,no2avg,nostepavg);
  int ifile,readret;
  int ix,iy,iz,iv;
  ldouble pp[NV],uu[NV];
  ldouble ttot=0.;
  ldouble t=global_time;

  //loop over avg files and average them
  doingpostproc=1;  
  doingavg=1;
  for(ifile=no1avg;ifile<=no2avg;ifile+=nostepavg)
  {
      readret=fread_avgfile(ifile,"dumps/avg",pavg,&dt,&t);
      add_u_core(1.,pavgtot,dt,pavg,pavgtot,avgarrsize);
      ttot+=dt;
  }

  //copy average array to pavgtot
  copy_u_core(1./ttot,pavgtot,pavg,avgarrsize);
 
  /*********************/
  //at this point pavg holds averaged quantities
  //loop over snap files now
  /*********************/

  //initialize profiles
  ldouble profiles[NANARELRADPROFILES][NX];      
  ldouble profilesavg[NANARELRADPROFILES][NX];   
  ldouble tprev=-1.;

  for(i=0;i<NANARELRADPROFILES;i++)
    for(j=0;j<NX;j++)
      profilesavg[i][j]=0.;

  //loop over res files
  doingavg=0;
  ttot=0.;
  for(ifile=no1res;ifile<=no2res;ifile+=nostepres)
  {
    //reading res file
    readret=fread_restartfile(ifile,"dumps",&t);
    nfout1=ifile;

    //put previous profile into the average
    if(tprev>0.)
    {
      dt=t-tprev;
      for(i=0;i<NANARELRADPROFILES;i++)
	for(j=0;j<NX;j++)
	  profilesavg[i][j]+=profiles[i][j]*dt;

      ttot+=dt;	  
    }
      
    //sets bc
    set_bc(t,0); 
						
    //calculates scale-height
    calc_avgs_throughout();      
      
    //calculate and dump profiles
#if(ANARELRADOUTPUT==1)
    calc_anarelradialprofiles(profiles);
    fprint_anarelradprofiles(t,nfout1,"analysis","relrad",profiles);
#endif

    tprev=t;
  }
  //the last one ignored

  //divide profiles by ttot
  for(i=0;i<NANARELRADPROFILES;i++)
    for(j=0;j<NX;j++)
      profilesavg[i][j]/=ttot;

  //save averaged profiles
  char prefix[40];
  char suffix[10];
  sprintf(suffix,"");
#if(ANARELRADOUTPUT==1)
  sprintf(prefix,"relradavg%s%04d-",suffix,no1res);
  fprint_anarelradprofiles(t,no2res,"analysis",prefix,profilesavg);
#endif

  return 0;
}
