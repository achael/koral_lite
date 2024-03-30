
//perturbed randomly velocity
int ix,iy,iz,ii,iv;
/**************************/

//perturbing according to given power spectrum
#ifdef PERT_SPECTRUM
/*
if(global_double1<0. || (global_time-global_double1)>0.001*CS0_ZERO)
  {
    global_double1=global_time;
*/

int verbose=0;
fftw_complex *delvx, *delvy,*Aamp;
fftw_plan plan;
//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TNX * TNY);
delvx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TNX * TNY);
delvy =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TNX * TNY);
Aamp =  (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * TNX * TNY);


//initializing random numbers
const gsl_rng_type * Tran;
gsl_rng * ran;
gsl_rng_env_setup();
Tran = gsl_rng_default;
unsigned long int s = rand();

//TEST - not to change perturbation pattern at every timestep
//s=(unsigned long int)global_time*10./(1./CS0_INIT);

//not change at all
//s=10;

ran = gsl_rng_alloc (Tran);
gsl_rng_set(ran,s);


//allocating for the FFT
//ldouble Aamp[TNX][TNY][2]; //amplitudes of vector potential, z -componen
//ldouble vampx[TNX][TNY][2]; //amplitudes of velocity perturbation, x - component
//ldouble vampy[TNX][TNY][2]; //amplitudes of velocity perturbation, y - component
ldouble kx[TNX];
ldouble ky[TNY];
ldouble Lx=(MAXX-MINX);
ldouble Ly=(MAXY-MINY);
ldouble kpeakx=2.*M_PI/Lx*PERT_SPECTRUM_MODE;
ldouble kpeaky=2.*M_PI/Ly*PERT_SPECTRUM_MODE;

kx[0]=0.;
for(ix=1;ix<=TNX/2;ix++)
  {
    kx[ix]=2.*M_PI/(Lx/(double)ix);
    kx[TNX-ix]=-kx[ix];
  }
ky[0]=0.;
for(iy=1;iy<=TNY/2;iy++)
  {
    ky[iy]=2.*M_PI/(Ly/(double)iy);
    ky[TNY-iy]=-ky[iy];
  }

//original
/*
kx[0]=0.;
for(ix=1;ix<TNX;ix++)
  {
    kx[ix]=2.*M_PI/(Lx/(double)ix);
  }
ky[0]=0.;
for(iy=1;iy<TNY;iy++)
  {
    ky[iy]=2.*M_PI/(Ly/(double)iy);
  }
*/

for(ix=0;ix<TNX;ix++)
  for(iy=0;iy<TNY;iy++)
  {
    ldouble k=sqrt(kx[ix]*kx[ix]+ky[iy]*ky[iy]);
    ldouble Alocx = kx[ix]*kx[ix]*exp(-4.*kx[ix]/kpeakx);
    ldouble Alocy = ky[ix]*ky[ix]*exp(-4.*ky[ix]/kpeaky);
    ldouble Alocz = k*k*exp(-4.*k/kpeakx);
    //ldouble Alocz = Alocx;//+Alocy;
    ldouble r;
    r=gsl_ran_gaussian(ran,1.0);
    //if(ix==0 ) printf("%e\n",r);
    Aamp[ix*TNY+iy][0]=Alocz*r; //Re
    r=gsl_ran_gaussian(ran,1.0);
    Aamp[ix*TNY+iy][1]=Alocz*r; //Im
  }

//  overwrite: 

/*
for(ix=0;ix<TNX;ix++)
  for(iy=0;iy<TNY;iy++)
  {
    Aamp[ix*TNY+iy][0]=Aamp[ix*TNY+iy][1]=0.;
  }

ldouble r;
r=gsl_ran_gaussian(ran,1.0);
Aamp[1*TNY+3][0]=r;
r=gsl_ran_gaussian(ran,1.0);
Aamp[0*TNY+1][0]=-r;
r=gsl_ran_gaussian(ran,1.0);
Aamp[2*TNY+4][1]=r;
r=gsl_ran_gaussian(ran,1.0);
Aamp[4*TNY+4][0]=r;
r=gsl_ran_gaussian(ran,1.0);
Aamp[1*TNY+1][0]=r;
r=gsl_ran_gaussian(ran,1.0);
Aamp[7*TNY+4][0]=r;
r=gsl_ran_gaussian(ran,1.0);
Aamp[7*TNY+3][1]=-r;
r=gsl_ran_gaussian(ran,1.0);
Aamp[6*TNY+2][0]=-r;


for(ix=0;ix<TNX;ix++)
  for(iy=0;iy<=TNY/2;iy++)
  {
    r=gsl_ran_gaussian(ran,1.0);
    Aamp[ix*TNY+iy][0]=Aamp[ix*TNY+iy][1]=r;
  }
*/

for(ix=0;ix<TNX;ix++)
  for(iy=0;iy<TNY;iy++)
    {
      delvx[ix*TNY+iy][0]= delvx[ix*TNY+iy][1]= delvy[ix*TNY+iy][0]= delvy[ix*TNY+iy][1]=0.;
    }

ldouble norm=DELTAV_NORM;

//|delta \vec v| = \vec k \cross \vec A
for(ix=0;ix<TNX;ix++)
  for(iy=0;iy<=TNY/2;iy++)
  {
    delvx[ix*TNY+iy][0]=norm*ky[iy]*Aamp[ix*TNY+iy][0];
    delvx[ix*TNY+iy][1]=norm*ky[iy]*Aamp[ix*TNY+iy][1];
    delvy[ix*TNY+iy][0]=-norm*kx[ix]*Aamp[ix*TNY+iy][0];
    delvy[ix*TNY+iy][1]=-norm*kx[ix]*Aamp[ix*TNY+iy][1];
  }

//overwrite
/*
for(ix=0;ix<=TNX/2;ix++)
  for(iy=0;iy<=TNY/2;iy++)
  {
    vampx[ix][iy][0]=vampx[ix][iy][1]=0.;
    vampy[ix][iy][0]=vampy[ix][iy][1]=0.;
  }


vampx[4][0][0]=1.5;    
vampx[2][0][0]=1.;    
vampx[0][4][0]=1.5;    
vampx[0][2][0]=1.;    
*/


int ixt,iyt;
//make it Hermittean

iy=0;
  for(ix=TNX/2+1;ix<TNX;ix++)
  {
    ixt=TNX/2-(ix-TNX/2);
    iyt=iy;
    if(ixt>=TNX) ixt-=TNX;
    if(iyt>=TNY) iyt-=TNY;
    if(ixt<0) ixt+=TNX;
    if(iyt<0) iyt+=TNY;
    delvx[ix*TNY+iy][0]=delvx[ixt*TNY+iyt][0];
    delvx[ix*TNY+iy][1]=-delvx[ixt*TNY+iyt][1];
    delvy[ix*TNY+iy][0]=delvy[ixt*TNY+iyt][0];
    delvy[ix*TNY+iy][1]=-delvy[ixt*TNY+iyt][1];
  }

iy=TNY/2;
  for(ix=TNX/2+1;ix<TNX;ix++)
  {
    ixt=TNX/2-(ix-TNX/2);
    iyt=iy;
    if(ixt>=TNX) ixt-=TNX;
    if(iyt>=TNY) iyt-=TNY;
    if(ixt<0) ixt+=TNX;
    if(iyt<0) iyt+=TNY;
    delvx[ix*TNY+iy][0]=delvx[ixt*TNY+iyt][0];
    delvx[ix*TNY+iy][1]=-delvx[ixt*TNY+iyt][1];
    delvy[ix*TNY+iy][0]=delvy[ixt*TNY+iyt][0];
    delvy[ix*TNY+iy][1]=-delvy[ixt*TNY+iyt][1];
  }


for(iy=TNY/2+1;iy<TNY;iy++)
  for(ix=0;ix<TNX;ix++)
  {
    ixt=TNX/2-(ix-TNX/2);
    iyt=TNY/2-(iy-TNY/2);
    if(ixt>=TNX) ixt-=TNX;
    if(iyt>=TNY) iyt-=TNY;
    if(ixt<0) ixt+=TNX;
    if(iyt<0) iyt+=TNY;
    delvx[ix*TNY+iy][0]=delvx[ixt*TNY+iyt][0];
    delvx[ix*TNY+iy][1]=-delvx[ixt*TNY+iyt][1];
    delvy[ix*TNY+iy][0]=delvy[ixt*TNY+iyt][0];
    delvy[ix*TNY+iy][1]=-delvy[ixt*TNY+iyt][1];
  }

delvx[0*TNY+0][1]=0.;
delvy[0*TNY+0][1]=0.;
delvx[TNX/2*TNY+0][1]=0.;
delvy[TNX/2*TNY+0][1]=0.;
delvx[TNX/2*TNY+TNY/2][1]=0.;
delvy[TNX/2*TNY+TNY/2][1]=0.;
delvx[0*TNY+TNY/2][1]=0.;
delvy[0*TNY+TNY/2][1]=0.;


if(verbose==2)
  {
for(ix=0;ix<TNX;ix++)
  for(iy=0;iy<TNY;iy++)
    {
      printf("%d %d %e %e\n",ix,iy,delvx[ix*TNY+iy][0],delvx[ix*TNY+iy][1]);
    }

printf("\n*******\n\n");
  }

//inverse FFT producing real data

//x-component 
/*
for(ix=0;ix<TNX;ix++)
  for(iy=0;iy<TNY;iy++)
    {
      in[ix*TNY+iy][0]=vampx[ix][iy][0];
      in[ix*TNY+iy][1]=vampx[ix][iy][1];
    }
*/
plan =  fftw_plan_dft_2d ( TNX, TNY, delvx, delvx, FFTW_BACKWARD, FFTW_ESTIMATE );
fftw_execute(plan);

//y-component 
/*
for(ix=0;ix<TNX;ix++)
  for(iy=0;iy<TNY;iy++)
    {
      in[ix*TNY+iy][0]=vampy[ix][iy][0];
      in[ix*TNY+iy][1]=vampy[ix][iy][1];
    }
*/
plan =  fftw_plan_dft_2d ( TNX, TNY, delvy, delvy, FFTW_BACKWARD, FFTW_ESTIMATE );
fftw_execute(plan);

if(verbose==2)
  {
for(ix=0;ix<TNX;ix++)
  for(iy=0;iy<TNY;iy++)
    {
    printf("%d %d %e %e %e %e\n",ix,iy,delvx[ix*TNY+iy][0],delvx[ix*TNY+iy][1],delvy[ix*TNY+iy][0],delvy[ix*TNY+iy][1]);
    }

getch();//exit(1);
  }

fftw_destroy_plan(plan);

//freeing random numbers
gsl_rng_free(ran);
#endif


ldouble dmx=0,dmy=0.;

#pragma omp parallel for private(ix,iy,iz,iv) schedule (dynamic)
for(ii=0;ii<Nloop_0;ii++) //domain only
  {
    ldouble pp[NV],uu[NV];
    struct geometry geom;
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2];
    fill_geometry(ix,iy,iz,&geom);

    ldouble xxBL[4];
    coco_N(geom.xxvec,xxBL,MYCOORDS, BLCOORDS);
    
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,ix,iy,iz);
      }
    ldouble drho=0.;
    ldouble duint=0;
    ldouble dvx=0.;
    ldouble dvy=0.;
    ldouble dvz=0.;

#ifdef PERT_SPECTRUM
    dvx=delvx[ix*TNY+iy][0];
    dvy=delvy[ix*TNY+iy][0];

    //TEST - to show the additional component only
    //dvx=-pp[VX]+delvx[ix*TNY+iy][0];
    //dvy=-pp[VY]+delvy[ix*TNY+iy][0];
    
    
    if(fabs(delvx[ix*TNY+iy][1]/delvx[ix*TNY+iy][0])>1.e-4)
      printf("non Real x-matrix > %d %d %e %e\n",ix,iy,delvx[ix*TNY+iy][0],delvx[ix*TNY+iy][1]);
    if(fabs(delvy[ix*TNY+iy][1]/delvy[ix*TNY+iy][0])>1.e-4)
      printf("non Real y-matrix > %d %d %e %e\n",iy,iy,delvy[ix*TNY+iy][0],delvy[ix*TNY+iy][1]);
#endif



#ifdef PERT_RANDOM

    dvx=VPERT*((double)rand()/(double)RAND_MAX-0.5);
    dvy=VPERT*((double)rand()/(double)RAND_MAX-0.5);
    
    
#endif



    pp[RHO]+=drho;
    pp[UU]+=duint;
    pp[VX]+=dvx;
    pp[VY]+=dvy;
    pp[VZ]+=dvz;


    //aggregiating momenta
    dmx+=dvx*pp[RHO];
    dmy+=dvy*pp[RHO];

   p2u(pp,uu,&geom);
    
    for(iv=0;iv<NV;iv++)
      {
	set_u(u,iv,ix,iy,iz,uu[iv]);
	set_u(p,iv,ix,iy,iz,pp[iv]);
      }

}

//substracting net momentum change

#pragma omp parallel for private(ix,iy,iz,iv) schedule (dynamic)
for(ii=0;ii<Nloop_0;ii++) //domain only
  {
    ldouble pp[NV],uu[NV];
    struct geometry geom;
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2];
    fill_geometry(ix,iy,iz,&geom);

    ldouble xxBL[4];
    coco_N(geom.xxvec,xxBL,MYCOORDS, BLCOORDS);
    
   
    for(iv=0;iv<NV;iv++)
      {
	pp[iv]=get_u(p,iv,ix,iy,iz);
      }

    ldouble dvx=-dmx/pp[RHO]/TNX/TNY; 
    ldouble dvy=-dmy/pp[RHO]/TNX/TNY;

    pp[VX]+=dvx;
    pp[VY]+=dvy;

    p2u(pp,uu,&geom);
    
    for(iv=0;iv<NV;iv++)
      {
	set_u(u,iv,ix,iy,iz,uu[iv]);
	set_u(p,iv,ix,iy,iz,pp[iv]);
      }
  }

#ifdef PERT_SPECTRUM
fftw_free(delvx);
fftw_free(delvy);
#endif
