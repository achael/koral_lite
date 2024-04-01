
//resets to initial state
int ix,iy,iz,ii,iv;
/**************************/
/*
#pragma omp parallel for private(ix,iy,iz,iv) schedule (dynamic)
for(ii=0;ii<Nloop_0;ii++) //domain only
  {
    ldouble pp[NV],uu[NV];
    struct geometry geom;
    ix=loop_0[ii][0];
    iy=loop_0[ii][1];
    iz=loop_0[ii][2]; 
    
    for(iv=0;iv<NV;iv++)
      pp[iv]=get_u(pinit,iv,ix,iy,iz);

    p2u(pp,uu,&geom);
    
    for(iv=0;iv<NV;iv++)
      {
	set_u(u,iv,ix,iy,iz,uu[iv]);
	set_u(p,iv,ix,iy,iz,pp[iv]);
      }
}
*/
