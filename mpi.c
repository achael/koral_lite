/*! \file mpi.c
 \brief MPI-related routines
*/

#include "ko.h"

int
mpi_exchangedata()
{
  //time mark
  struct timespec temp_clock;
  my_clock_gettime(&temp_clock);    
  mid1_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

#ifdef MPI  
  MPI_Request reqs[MPIMSGBUFSIZE];
  int nreqs=0;
  mpi_senddata(reqs,&nreqs);
  mpi_recvdata(reqs,&nreqs);

  MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
  mpi_savedata();  

#endif

  my_clock_gettime(&temp_clock);    
  mid2_time=(ldouble)temp_clock.tv_sec+(ldouble)temp_clock.tv_nsec/1.e9;

  return 0;
}

#ifdef MPI //comments out send and recv so that non-mpi compiler can compile without MPI_Request

int
mpi_senddata(MPI_Request *reqs, int *nreqs)
{  
  int i,j,k,iv;
  int tx,ty,tz;
  int verbose=0;
  ldouble temp;

  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif

      for(i=0;i<NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[0][(i)*NY*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[0], NY*NZ*NV*NG, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d sent MPI_MSG_XLO to %d\n",PROCID,mpi_tile2procid(tx,ty,tz));
    }

  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[1][(i-NX+NG)*NY*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[1],NG*NY*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d sent MPI_MSG_XHI to %d\n",PROCID,mpi_tile2procid(tx,ty,tz));
    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0) //send YLO domain to YLO GC pi away at pole
      {
        tx=TI;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        for(i=0;i<NX;i++)
        {
	  for(j=0;j<NG;j++)
	  {  
            for(k=0;k<NZ;k++)
            {
	      for(iv=0;iv<NV;iv++)
              {
	        msgbufs[2][i*NG*NZ*NV + (j)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
              }
            }
          }
        } 
        MPI_Isend(msgbufs[2], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ-1;
        tz=TK;
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        for(i=0;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
	        msgbufs[2][i*NG*NZ*NV + (j)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
        MPI_Isend(msgbufs[2], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1) //send YHI domain to YHI GC pi away at pole
      {
        tx=TI;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=0;i<NX;i++)
   	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
	        msgbufs[3][i*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
        MPI_Isend(msgbufs[3], NX*NG*NZ*NV, MPI_DOUBLE,
		  mpi_tile2procid(tx,ty,tz), MPI_MSG_YHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ+1;
        tz=TK;
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        for(i=0;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
	        msgbufs[3][i*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
        MPI_Isend(msgbufs[3], NX*NG*NZ*NV, MPI_DOUBLE,
		  mpi_tile2procid(tx,ty,tz), MPI_MSG_YHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
        #ifdef TRANSMITTING_YBC
        }
        #endif
    }

  //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      tx=TI;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NG;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[4][i*NY*NG*NV + j*NG*NV + (k)*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[4], NX*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_ZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      tx=TI;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      msgbufs[5][i*NY*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
      MPI_Isend(msgbufs[5], NX*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_ZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }

#ifdef MPI4CORNERS

  /***************************/
  //elongated corners - along z
  //lower x lower y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0) //send YLO domain to YLO GC pi away at pole
      {
        tx=TI-1;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif

        for(i=0;i<NG;i++)
  	  for(j=0;j<NG;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
                {
	  	  msgbufs[12][i*NG*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
                }  
	      }
        MPI_Isend(msgbufs[12], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ-1;
        tz=TK;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif

        for(i=0;i<NG;i++)
  	  for(j=0;j<NG;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
	  	  msgbufs[12][i*NG*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[12], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

  //lower x higher y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1) //send YHI domain to YHI GC pi away at pole
      {
        tx=TI-1;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
 
        for(i=0;i<NG;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[13][i*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[13], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ+1;
        tz=TK;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif

        for(i=0;i<NG;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[13][i*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[13], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }
  
  //higher x lower y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0) //send YLO domain to YLO GC pi away at pole
      {
        tx=TI+1;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif

        for(i=NX-NG;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[14][(i-NX+NG)*NG*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[14], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ-1;
        tz=TK;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif

        for(i=NX-NG;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[14][(i-NX+NG)*NG*NZ*NV + j*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[14], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

  //higher x higher y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1) //send YHI domain to YHI GC pi away at pole
      {
        tx=TI+1;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
 
        for(i=NX-NG;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[15][(i-NX+NG)*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[15], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ+1;
        tz=TK;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif

        for(i=NX-NG;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[15][(i-NX+NG)*NG*NZ*NV + (j-NY+NG)*NZ*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[15], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

  /***************************/
  //elongated corners - along y
  //lower x lower z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=0;i<NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[16][i*NY*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[16], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  //lower x higher z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=0;i<NG;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[17][i*NY*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[17], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  
  //higher x lower z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[18][(i-NX+NG)*NY*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[18], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

  //higher x higher z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif

      for(i=NX-NG;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ-NG;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		msgbufs[19][(i-NX+NG)*NY*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	    }
      MPI_Isend(msgbufs[19], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }  

   /***************************/
  //elongated corners - along x
  //lower y lower z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0) //send YLO domain to YLO GC pi away at pole
      {
        tx=TI;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        for(i=0;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[20][i*NG*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[20], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ-1;
        tz=TK-1;
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif

        for(i=0;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[20][i*NG*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[20], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

  //lower y higher z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)  //send YLO domain to YLO GC pi away at pole
      {
        tx=TI;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        for(i=0;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[21][i*NG*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[21], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ-1;
        tz=TK+1;
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif

        for(i=0;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[21][i*NG*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[21], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }
  
  //higher y lower z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1) //send YHI domain to YHI GC pi away at pole
      {
        tx=TI;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=0;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[22][i*NG*NG*NV + (j-NY+NG)*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[22], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ+1;
        tz=TK-1;
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif

        for(i=0;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[22][i*NG*NG*NV + (j-NY+NG)*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[22], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

  //higher y higher z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1) //send YHI domain to YHI GC pi away at pole
      {
        tx=TI;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=0;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[23][i*NG*NG*NV + (j-NY+NG)*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[23], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ+1;
        tz=TK+1;
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif

        for(i=0;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[23][i*NG*NG*NV + (j-NY+NG)*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[23], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

  /***************************/
  //corners corners 
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0) //send YLO domain to YLO GC pi away at pole
      {
        tx=TI-1;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=0;i<NG;i++)
  	  for(j=0;j<NG;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[24][i*NG*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[24], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ-1;
        tz=TK-1;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif

        for(i=0;i<NG;i++)
  	  for(j=0;j<NG;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[24][i*NG*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[24], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0) //send YLO domain to YLO GC pi away at pole
      {
        tx=TI-1;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=0;i<NG;i++)
	  for(j=0;j<NG;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[25][i*NG*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[25], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ-1;
        tz=TK+1;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif

        for(i=0;i<NG;i++)
	  for(j=0;j<NG;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[25][i*NG*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[25], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

 if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == (NTY-1)) //send YHI domain to YHI GC pi away at pole
      {
        tx=TI-1;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=0;i<NG;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[26][i*NG*NG*NV + (j-NY+NG)*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[26], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ+1;
        tz=TK-1;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif

        for(i=0;i<NG;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[26][i*NG*NG*NV + (j-NY+NG)*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[26], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

 if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == (NTY-1)) //send YHI domain to YHI GC pi away at pole
      {
        tx=TI-1;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=0;i<NG;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[27][i*NG*NG*NV + (j-NY+NG)*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[27], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ+1;
        tz=TK+1;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif

        for(i=0;i<NG;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[27][i*NG*NG*NV + (j-NY+NG)*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[27], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0) //send YLO domain to YLO GC pi away at pole
      {
        tx=TI+1;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=NX-NG;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[28][(i-NX+NG)*NG*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[28], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ-1;
        tz=TK-1;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif

        for(i=NX-NG;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[28][(i-NX+NG)*NG*NG*NV + j*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[28], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0) //send YLO domain to YLO GC pi away at pole
      {
        tx=TI+1;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=NX-NG;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[29][(i-NX+NG)*NG*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[29], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ-1;
        tz=TK+1;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif

        for(i=NX-NG;i<NX;i++)
	  for(j=0;j<NG;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[29][(i-NX+NG)*NG*NG*NV + j*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[29], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

 if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == (NTY-1)) //send YHI domain to YHI GC pi away at pole
      {
        tx=TI+1;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=NX-NG;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[30][(i-NX+NG)*NG*NG*NV + (j-NY+NG)*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[30], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ+1;
        tz=TK-1;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif

        for(i=NX-NG;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=0;k<NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[30][(i-NX+NG)*NG*NG*NV + (j-NY+NG)*NG*NV + k*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[30], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  

 if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == (NTY-1)) //send YHI domain to YHI GC pi away at pole
      {
        tx=TI+1;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
 
        for(i=NX-NG;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[31][(i-NX+NG)*NG*NG*NV + (j-NY+NG)*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[31], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ+1;
        tz=TK+1;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif

        for(i=NX-NG;i<NX;i++)
	  for(j=NY-NG;j<NY;j++)
	    for(k=NZ-NG;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  msgbufs[31][(i-NX+NG)*NG*NG*NV + (j-NY+NG)*NG*NV + (k-NZ+NG)*NV + iv]=get_u(p,iv,i,j,k);
	      }
        MPI_Isend(msgbufs[31], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }  
#endif

  return 0;
}

int
mpi_recvdata(MPI_Request *reqs, int *nreqs)
{
  int i,j,k,iv;
  MPI_Status status;
  double temp;
  int verbose=0;
  int tx,ty,tz;

  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      MPI_Irecv(msgbufs[6], NY*NZ*NV*NG, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_XLO from %d\n",PROCID,mpi_tile2procid(TI+1,TJ,TK));
   }
  
  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      MPI_Irecv(msgbufs[7], NG*NY*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_XHI from %d\n",PROCID,mpi_tile2procid(TI-1,TJ,TK));
    } 

  //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YHI GC from YHI domain pi away at pole
      if(TJ == NTY-1)
      {
        tx=TI;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
        
        MPI_Irecv(msgbufs[8], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
        if(verbose) printf("%d received MPI_MSG_YHI from %d\n",PROCID,mpi_tile2procid(TI,TJ+1,TK));
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ+1;
        tz=TK;
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        MPI_Irecv(msgbufs[8], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
        if(verbose) printf("%d received MPI_MSG_YLO from %d\n",PROCID,mpi_tile2procid(TI,TJ+1,TK));
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YLO GC from YLO domain pi away at pole
      if(TJ==0)
      {
        tx=TI;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[9], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLO, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
        if(verbose) printf("%d received MPI_MSG_YLO from %d\n",PROCID,mpi_tile2procid(TI,TJ-1,TK));
         
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ-1;
        tz=TK;
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        MPI_Irecv(msgbufs[9], NX*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHI, MPI_COMM_WORLD,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
        if(verbose) printf("%d received MPI_MSG_YHI from %d\n",PROCID,mpi_tile2procid(TI,TJ-1,TK));
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {							
      tx=TI;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[10], NX*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_ZLO, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
    }
  //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      tx=TI;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[11], NX*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_ZHI, MPI_COMM_WORLD,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
      if(verbose) printf("%d received MPI_MSG_ZHI from %d\n",PROCID,mpi_tile2procid(TI,TJ,TK-1));
    }

  //corners
#ifdef MPI4CORNERS
  //elongated along z
  //upper x upper y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YHI GC from YHI domain pi away at pole
      if(TJ==(NTY-1))
      {
        tx=TI+1;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif

        MPI_Irecv(msgbufs[32], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ+1;
        tz=TK;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        MPI_Irecv(msgbufs[32], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }
  //upper x lower y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YLO GC from YLO domain pi away at pole
      if(TJ==0)
      {
        tx=TI+1;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif

        MPI_Irecv(msgbufs[33], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ-1;
        tz=TK;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        MPI_Irecv(msgbufs[33], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }
  //lower x upper y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YHI GC from YHI domain pi away at pole
      if(TJ==(NTY-1))
      {
        tx=TI-1;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif

        MPI_Irecv(msgbufs[34], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ+1;
        tz=TK;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        MPI_Irecv(msgbufs[34], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }
  //lower x lower y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YLO GC from YLO domain pi away at pole
      if(TJ==0)
      {
        tx=TI-1;
        ty=TJ;
        tz=TK+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif

        MPI_Irecv(msgbufs[35], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ-1;
        tz=TK;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        MPI_Irecv(msgbufs[35], NG*NG*NZ*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

  //elongated along y
  //upper x upper z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[36], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //upper x lower z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI+1;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx>=NTX) tx-=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[37], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //lower x upper z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK+1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz>=NTZ) tz-=NTZ;
      #endif
      MPI_Irecv(msgbufs[38], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }
  //lower x lower z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      tx=TI-1;
      ty=TJ;
      tz=TK-1;
      #ifdef PERIODIC_XBC
      if(tx<0) tx+=NTX;
      #endif
      #ifdef PERIODIC_ZBC
      if(tz<0) tz+=NTZ;
      #endif
      MPI_Irecv(msgbufs[39], NG*NY*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
      *nreqs=*nreqs+1;
   }

   //elongated along x
  //upper y upper z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YHI GC from YHI domain pi away at pole
      if(TJ==(NTY-1))
      {
        tx=TI;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[40], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ+1;
        tz=TK+1;
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif
        MPI_Irecv(msgbufs[40], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }
  //upper y lower z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YHI GC from YHI domain pi away at pole
      if(TJ==(NTY-1))
      {
        tx=TI;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[41], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ+1;
        tz=TK-1;
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif
        MPI_Irecv(msgbufs[41], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }
  //lower y upper z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YLO GC from YLO domain pi away at pole
      if(TJ==0)
      {
        tx=TI;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[42], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ-1;
        tz=TK+1;
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif
        MPI_Irecv(msgbufs[42], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }
  //lower x lower z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YLO GC from YLO domain pi away at pole
      if(TJ==0)
      {
        tx=TI;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[43], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI;
        ty=TJ-1;
        tz=TK-1;
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif
        MPI_Irecv(msgbufs[43], NX*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_YHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

  /********** corners corners ************/
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YHI GC from YHI domain pi away at pole
      if(TJ==(NTY-1))
      {
        tx=TI+1;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[44], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ+1;
        tz=TK+1;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif
        MPI_Irecv(msgbufs[44], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YHI GC from YHI domain pi away at pole
      if(TJ==(NTY-1))
      {
        tx=TI+1;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[45], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ+1;
        tz=TK-1;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif
        MPI_Irecv(msgbufs[45], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

 if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YLO GC from YLO domain pi away at pole
      if(TJ==0)
      {
        tx=TI+1;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[46], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ-1;
        tz=TK+1;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif
        MPI_Irecv(msgbufs[46], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

 if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YLO GC from YLO domain pi away at pole
      if(TJ==0)
      {
        tx=TI+1;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[47], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI+1;
        ty=TJ-1;
        tz=TK-1;
        #ifdef PERIODIC_XBC
        if(tx>=NTX) tx-=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif
        MPI_Irecv(msgbufs[47], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XLOYHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YHI GC from YHI domain pi away at pole
      if(TJ==(NTY-1))
      {
        tx=TI-1;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[48], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ+1;
        tz=TK+1;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif
        MPI_Irecv(msgbufs[48], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YHI GC from YHI domain pi away at pole
      if(TJ==(NTY-1))
      {
        tx=TI-1;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[49], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ+1;
        tz=TK-1;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty>=NTY) ty-=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif
        MPI_Irecv(msgbufs[49], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

 if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YLO GC from YLO domain pi away at pole
      if(TJ==0)
      {
        tx=TI-1;
        ty=TJ;
        tz=TK+1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[50], NG*NG*NG*NV, MPI_DOUBLE,
 		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ-1;
        tz=TK+1;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz>=NTZ) tz-=NTZ;
        #endif
        MPI_Irecv(msgbufs[50], NG*NG*NG*NV, MPI_DOUBLE,
 		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZLO, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

 if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC //recieve YLO GC from YLO domain pi away at pole
      if(TJ==0)
      {
        tx=TI-1;
        ty=TJ;
        tz=TK-1+(NTZ/2);
        if(tz > NTZ-1) tz -= NTZ;

        MPI_Irecv(msgbufs[51], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYLOZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      }
      else
      {
      #endif
        tx=TI-1;
        ty=TJ-1;
        tz=TK-1;
        #ifdef PERIODIC_XBC
        if(tx<0) tx+=NTX;
        #endif
        #ifdef PERIODIC_YBC
        if(ty<0) ty+=NTY;
        #endif
        #ifdef PERIODIC_ZBC
        if(tz<0) tz+=NTZ;
        #endif
        MPI_Irecv(msgbufs[51], NG*NG*NG*NV, MPI_DOUBLE,
		mpi_tile2procid(tx,ty,tz), MPI_MSG_XHIYHIZHI, MPI_COMM_WORLD ,&reqs[*nreqs]);
        *nreqs=*nreqs+1;
      #ifdef TRANSMITTING_YBC
      }
      #endif
   }

#endif  

  return 0;
}

int
mpi_savedata()
{
  int i,j,k,iv;
  int verbose=0;

  //upper x
  if(mpi_isitBC(XBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[6][(i-NX)*NY*NZ*NV + j*NZ*NV + k*NV + iv]);
	    }
    }
  
  //lower x
  if(mpi_isitBC(XBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=0;j<NY;j++)
	  for(k=0;k<NZ;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[7][(i+NG)*NY*NZ*NV + j*NZ*NV + k*NV + iv]);
    } 

 //upper y
  if(mpi_isitBC(YBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1)
      {
        int jp = 0; //Needed since the ghost cell index is inverted relative to what is sent?

        for(i=0;i<NX;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j + (2*NY + NG - 1);
                #ifdef MAGNFIELD //need to flip sign of fluxes and velocities across pole, theta and phi
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[8][i*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[8][i*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=0;i<NX;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
	        set_u(p,iv,i,j,k,msgbufs[8][i*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  //lower y
  if(mpi_isitBC(YBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)
      {
        int jp = 0; //Needed since the ghost cell index is inverted relative to what is sent?

        for(i=0;i<NX;i++)
	  for(j=-NG;j<0;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j - (NG+1);
                #ifdef MAGNFIELD //need to flip sign of fluxes and velocities across pole, theta and phi
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[9][i*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[9][i*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=0;i<NX;i++)
	  for(j=-NG;j<0;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
	        set_u(p,iv,i,j,k,msgbufs[9][i*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  //upper z
  if(mpi_isitBC(ZBCHI)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[10][i*NY*NG*NV + j*NG*NV + (k-NZ)*NV + iv]);
    }
   //lower z
  if(mpi_isitBC(ZBCLO)==0)
    {
      for(i=0;i<NX;i++)
	for(j=0;j<NY;j++)
	  for(k=-NG;k<0;k++)
	    for(iv=0;iv<NV;iv++)
	      set_u(p,iv,i,j,k,msgbufs[11][i*NY*NG*NV + j*NG*NV + (k+NG)*NV + iv]);
    }

  //corners
#ifdef MPI4CORNERS
  //elongated along z
  //upper x upper y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1)
      {
        int jp = 0; //Needed since the ghost cell index is inverted relative to what is sent?

        for(i=NX;i<NX+NG;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j + (2*NY + NG - 1);
                #ifdef MAGNFIELD //need to flip sign of fluxes and velocities across pole, theta and phi
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[32][(i-NX)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[32][(i-NX)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=NX;i<NX+NG;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[32][(i-NX)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }
  //upper x lower y
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)
      {
        int jp = 0;
        for(i=NX;i<NX+NG;i++)
	  for(j=-NG;j<0;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j - (NG+1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[33][(i-NX)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[33][(i-NX)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=NX;i<NX+NG;i++)
	  for(j=-NG;j<0;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[33][(i-NX)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }
  //lower x upper y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1)
      {
        int jp = 0; //Needed since the ghost cell index is inverted relative to what is sent?

        for(i=-NG;i<0;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j + (2*NY + NG - 1);
                #ifdef MAGNFIELD //need to flip sign of fluxes and velocities across pole, theta and phi
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[34][(i+NG)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[34][(i+NG)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=-NG;i<0;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[34][(i+NG)*NG*NZ*NV + (j-NY)*NZ*NV + k*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }
  //lower x lower y
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)
      {
        int jp = 0;
        for(i=-NG;i<0;i++)
	  for(j=-NG;j<0;j++)
	    for(k=0;k<NZ;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j - (NG+1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[35][(i+NG)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[35][(i+NG)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=-NG;i<0;i++)
	  for(j=-NG;j<0;j++)
	    for(k=0;k<NZ;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[35][(i+NG)*NG*NZ*NV + (j+NG)*NZ*NV + k*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  //elongated along y
  //upper x upper z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[36][(i-NX)*NY*NG*NV + j*NG*NV + (k-NZ)*NV + iv]);
	    }
    }
  //upper x lower z
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=0;j<NY;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[37][(i-NX)*NY*NG*NV + j*NG*NV + (k+NG)*NV + iv]);
	    }
    }
  //lower x upper z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=0;j<NY;j++)
	  for(k=NZ;k<NZ+NG;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[38][(i+NG)*NY*NG*NV + j*NG*NV + (k-NZ)*NV + iv]);
	    }
    }
  //lower x lower z
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=0;j<NY;j++)
	  for(k=-NG;k<0;k++)
	    {
	      for(iv=0;iv<NV;iv++)
		set_u(p,iv,i,j,k,msgbufs[39][(i+NG)*NY*NG*NV + j*NG*NV + (k+NG)*NV + iv]);
	    }
    }

  //elongated along x
  //upper y upper z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1)
      {
        int jp = 0; //Needed since the ghost cell index is inverted relative to what is sent?

        for(i=0;i<NX;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j + (2*NY + NG - 1);
                #ifdef MAGNFIELD //need to flip sign of fluxes and velocities across pole, theta and phi
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[40][i*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[40][i*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=0;i<NX;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[40][i*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }
  //upper y lower z
  if(mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == NTY-1)
      {
        int jp = 0; //Needed since the ghost cell index is inverted relative to what is sent?

        for(i=0;i<NX;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=-NG;k<0;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j + (2*NY + NG - 1);
                #ifdef MAGNFIELD //need to flip sign of fluxes and velocities across pole, theta and phi
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[41][i*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[41][i*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=0;i<NX;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=-NG;k<0;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[41][i*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }
  //lower y upper z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)
      {
        int jp = 0;
        for(i=0;i<NX;i++)
	  for(j=-NG;j<0;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j - (NG+1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[42][i*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[42][i*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=0;i<NX;i++)
	  for(j=-NG;j<0;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[42][i*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }
  //lower y lower z
  if(mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)
      {
        int jp = 0;
        int index = 0;
        int gix,giy,giz;

        for(i=0;i<NX;i++)
	  for(j=-NG;j<0;j++)
	    for(k=-NG;k<0;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j - (NG+1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[43][i*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[43][i*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=0;i<NX;i++)
	  for(j=-NG;j<0;j++)
	    for(k=-NG;k<0;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[43][i*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  //corners corners
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == (NTY-1))
      {
        int jp = 0;
        for(i=NX;i<NX+NG;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j + (2*NY + NG - 1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[44][(i-NX)*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[44][(i-NX)*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=NX;i<NX+NG;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[44][(i-NX)*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == (NTY-1))
      {
        int jp = 0;
        for(i=NX;i<NX+NG;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=-NG;k<0;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j + (2*NY + NG - 1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[45][(i-NX)*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[45][(i-NX)*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=NX;i<NX+NG;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=-NG;k<0;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[45][(i-NX)*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }


  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)
      {
        int jp = 0;
        for(i=NX;i<NX+NG;i++)
	  for(j=-NG;j<0;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j - (NG+1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[46][(i-NX)*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[46][(i-NX)*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=NX;i<NX+NG;i++)
	  for(j=-NG;j<0;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[46][(i-NX)*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)
      {
        int jp = 0;
        for(i=NX;i<NX+NG;i++)
	  for(j=-NG;j<0;j++)
	    for(k=-NG;k<0;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j - (NG+1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[47][(i-NX)*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[47][(i-NX)*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=NX;i<NX+NG;i++)
	  for(j=-NG;j<0;j++)
	    for(k=-NG;k<0;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[47][(i-NX)*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == (NTY-1))
      {
        int jp = 0;
        for(i=-NG;i<0;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j + (2*NY + NG - 1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[48][(i+NG)*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[48][(i+NG)*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=-NG;i<0;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[48][(i+NG)*NG*NG*NV + (j-NY)*NG*NV + (k-NZ)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == (NTY-1))
      {
        int jp = 0;
        for(i=-NG;i<0;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=-NG;k<0;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j + (2*NY + NG - 1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[49][(i+NG)*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[49][(i+NG)*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=-NG;i<0;i++)
	  for(j=NY;j<NY+NG;j++)
	    for(k=-NG;k<0;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[49][(i+NG)*NG*NG*NV + (j-NY)*NG*NV + (k+NG)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }


  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCHI)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)
      {
        int jp = 0;
        for(i=-NG;i<0;i++)
	  for(j=-NG;j<0;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j - (NG+1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[50][(i+NG)*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[50][(i+NG)*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=-NG;i<0;i++)
	  for(j=-NG;j<0;j++)
	    for(k=NZ;k<NZ+NG;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[50][(i+NG)*NG*NG*NV + (j+NG)*NG*NV + (k-NZ)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }

  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0 && mpi_isitBC(ZBCLO)==0)
    {
      #ifdef TRANSMITTING_YBC
      if(TJ == 0)
      {
        int jp = 0;
        for(i=-NG;i<0;i++)
	  for(j=-NG;j<0;j++)
	    for(k=-NG;k<0;k++)
	      for(iv=0;iv<NV;iv++)
              {
                jp = -j - (NG+1);
                #ifdef MAGNFIELD
                if(iv==VY || iv==B2 || iv==FY0 || iv==VZ || iv==B3 || iv==FZ0)
                #else
                if(iv==VY || iv==FY0 || iv==VZ || iv==FZ0)
                #endif
                {
	          set_u(p,iv,i,jp,k,(-1)*msgbufs[51][(i+NG)*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
                }
                else
                {
	          set_u(p,iv,i,jp,k,msgbufs[51][(i+NG)*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
                }
              }
      }
      else
      {
      #endif
        for(i=-NG;i<0;i++)
	  for(j=-NG;j<0;j++)
	    for(k=-NG;k<0;k++)
	      {
	        for(iv=0;iv<NV;iv++)
		  set_u(p,iv,i,j,k,msgbufs[51][(i+NG)*NG*NG*NV + (j+NG)*NG*NV + (k+NG)*NV + iv]);
	      }
      #ifdef TRANSMITTING_YBC
      }
      #endif
    }


#endif


  return 0; 
}
#endif //MPI


///////////////////////////////////////////////////////////////
//verify if we are on an outer MPI tile
int
mpi_hasBC()
{
#ifndef MPI
  return 1;
#else
  if(TI==0 || TI==NTX-1 || TJ==0 || TJ==NTY-1 || TK==0 || TK==NTZ-1)
    return 1; //this cell has some real BC
  else
    return 0;
#endif
}

///////////////////////////////////////////////////////////////
//verify if our current tile has a particular BC type
int
mpi_isitBC(int BCtype)
{
#ifndef MPI
  return 1;
#else //check here if real BC
  int perx,pery,perz; //is periodic in x,y,z?
  perx=pery=perz=0; 
  #ifdef PERIODIC_XBC
  perx=1;
  #endif
  #ifdef PERIODIC_YBC
  pery=1;
  #endif
  #ifdef PERIODIC_ZBC
  perz=1;
  #endif
  #ifdef TRANSMITTING_YBC //Special pery setting to tell mpi code to pass information for transmitting boundary
  pery=2;
  #endif

  if(BCtype==XBCLO && TI==0 && perx==0)
    return 1;
  if(BCtype==XBCHI && TI==NTX-1 && perx==0)
    return 1;
  if(BCtype==YBCLO && TJ==0 && pery==0)
    return 1;
  if(BCtype==YBCHI && TJ==NTY-1 && pery==0)
    return 1;
  if(BCtype==ZBCLO && TK==0 && perz==0)
    return 1;
  if(BCtype==ZBCHI && TK==NTZ-1 && perz==0)
    return 1;

  if(BCtype==YBCLO && TJ==0 && pery==2)
    return 0; //Set GC at YBCLO using tile pi away
  if(BCtype==YBCHI && TJ==NTY-1 && pery==2)
    return 0; //Set GC at YBCLO using tile pi away

  return 0; 
#endif
}

///////////////////////////////////////////////////////////////
//verify if given cell from set_bc() falls into a real BC, only used in setbc since transmitting boundary still needs to have corners filled when at real corner (i.e. XBCHI & YBCLO)

int
mpi_isitBC_forcorners(int BCtype)
{
#ifndef MPI
  return 1;
#else //check here if real BC
  int perx,pery,perz; //is periodic in x,y,z?
  perx=pery=perz=0; 
  #ifdef PERIODIC_XBC
  perx=1;
  #endif
  #ifdef PERIODIC_YBC
  pery=1;
  #endif
  #ifdef PERIODIC_ZBC
  perz=1;
  #endif
  #ifdef TRANSMITTING_YBC //Special pery setting to tell mpi code to pass information for transmitting boundary
  pery=2;
  #endif

  if(BCtype==XBCLO && TI==0 && perx==0)
    return 1;
  if(BCtype==XBCHI && TI==NTX-1 && perx==0)
    return 1;
  if(BCtype==YBCLO && TJ==0 && pery==0)
    return 1;
  if(BCtype==YBCHI && TJ==NTY-1 && pery==0)
    return 1;
  if(BCtype==ZBCLO && TK==0 && perz==0)
    return 1;
  if(BCtype==ZBCHI && TK==NTZ-1 && perz==0)
    return 1;

  if(BCtype==YBCLO && TJ==0 && pery==2)
    return 1; //Return 1 so that setbc code will fill corners (see setbc in finite.c)
  if(BCtype==YBCHI && TJ==NTY-1 && pery==2)
    return 1; //Return 1 so that setbc code will fill corners (see setbc in finite.c)

  return 0; 
#endif
}

///////////////////////////////////////////////////////////////
void
mpi_synchtiming(ldouble *time)
{
#ifdef MPI

  MPI_Barrier(MPI_COMM_WORLD);
  
  global_time=*time;
  
  MPI_Allreduce(&tstepdenmax, &global_tstepdenmax, 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);

  //maximal time taken by information exchange
  ldouble localmp_time = mid2_time-mid1_time;
  
  //only PROCID==0 writes to stdout
  MPI_Reduce(&localmp_time, &maxmp_time, 1, MPI_DOUBLE, MPI_MAX, 0,
	     MPI_COMM_WORLD);   

  //total operation time
  ldouble local_u2ptime=end_u2ptime-start_u2ptime;
  struct {
    double time;
    int ti;
  } in,outmin,outmax;
  in.time=local_u2ptime;
  in.ti=TI;

  MPI_Reduce(&in, &outmax, 1, MPI_DOUBLE_INT, MPI_MAXLOC,0,
		 MPI_COMM_WORLD);  
  MPI_Reduce(&in, &outmin, 1, MPI_DOUBLE_INT, MPI_MINLOC,0,
		 MPI_COMM_WORLD);  

  max_u2ptime=outmax.time;
  max_u2ptime_loc=outmax.ti;
  min_u2ptime=outmin.time;
  min_u2ptime_loc=outmin.ti;  

  MPI_Bcast(&global_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  MPI_Barrier(MPI_COMM_WORLD);

  *time=global_time;
  tstepdenmax=global_tstepdenmax;
  tstepdenmin=global_tstepdenmin;
#endif
}

///////////////////////////////////////////////////////////////
void
mpi_myinit(int argc, char *argv[])
{
#ifdef MPI
  int i,j,k;

  // Check if MPI tiling is okay  
  if (TNX % NTX != 0)
  {
    printf("\nERROR!!  NTX = %d  TNX = %d  not divisible!\n\n", NTX, TNX);
    exit(1);
  }
  if (TNY % NTY != 0)
  {
    printf("\nERROR!!  NTY = %d  TNY = %d  not divisible!\n\n", NTY, TNY);
    exit(1);
  }
  if (TNZ % NTZ != 0)
  {
    printf("\nERROR!!  NTZ = %d  TNZ = %d  not divisible!\n\n", NTZ, TNZ);
    exit(1);
  }
  
  //check for conflicts in declarations
  #ifndef OUTPUTPERCORE
  #ifdef RESOUTPUT_ASCII
  my_err("RESOUTPUT_ASCII requires MPI_OUTPUTPERCORE\n");exit(-1);
  #endif
  #ifdef AVGOUTPUT_ASCII
  my_err("RESOUTPUT_ASCII requires MPI_OUTPUTPERCORE\n");exit(-1);
  #endif
  #endif
  #ifdef CALCHRONTHEGO
  if(TNZ>1) //3D tiling
    {
      if(PROCID==0) printf("CALCHRONTHEGO not implemented for 3D.\n");
      exit(-1);
    }
  #endif

  //initialize MPI
  int provided, required = MPI_THREAD_FUNNELED;
  MPI_Init_thread(&argc, &argv, required, &provided);
  
  if (provided >= required)
  {
    if(PROCID==0)
    {
      //printf("\nThread Support: required %d  provided %d\n\n", required, provided);
      //printf("Number of OMP Threads: %d\n\n", omp_get_num_threads());
    }
  }
  else
  {
    // Insufficient support, degrade to 1 thread and warn the user
    if (PROCID==0)
    {
      printf("\nWarning: This MPI implementation provides insufficient threading support\n");
      printf("Thread Support: required %d  provided %d\n\n", required, provided);
    }
    exit(1);
    //omp_set_num_threads(1);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &PROCID);
  MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);

  if(NPROCS!=NTX*NTY*NTZ)
    {
      if(PROCID==0) printf("Wrong number of processes. Problem set up for: %d x %d x %d = %d processes. \n But NPROCS = %d\n",NTX,NTY,NTZ,NTX*NTY*NTZ,NPROCS);
      exit(-1);
    }
  
  mpi_procid2tile(PROCID,&TI,&TJ,&TK);
  mpi_tileorigin(TI,TJ,TK,&TOI,&TOJ,&TOK);

  // useful info on how processors are assigned
  //printf("PROCID TI TJ TK TOI TOJ TOK: %d %d %d %d %d %d %d\n", PROCID, TI, TJ,TK, TOI, TOJ, TOK);

  if(PROCID==0) printf("pid: %d/%d; tot.res: %dx%dx%d; tile.res:  %dx%dx%d\n"
	 "tile: %d,%d,%d; tile orig.: %d,%d,%d\n",PROCID,NPROCS,TNX,TNY,TNZ,NX,NY,NZ,TI,TJ,TK,TOI,TOJ,TOK);

#else //no MPI
  TI=TJ=TK=0;
  TOI=TOJ=TOK=0;
  PROCID=0;
  NPROCS=1;
#endif
}


///////////////////////////////////////////////////////////////
void
mpi_myfinalize()
{
#ifdef MPI
  MPI_Finalize();
#endif
}


///////////////////////////////////////////////////////////////
void
mpi_tileorigin(int ti, int tj, int tk, int* toi, int* toj, int* tok)
{
  *toi = ti * (TNX/NTX);
  *toj = tj * (TNY/NTY);
  *tok = tk * (TNZ/NTZ);
}


///////////////////////////////////////////////////////////////
void
mpi_global2localidx(int gix,int giy, int giz, int *lix, int *liy, int *liz)
{
  #ifdef MPI
  int tilei,tilej,tilek;
  int toi,toj,tok;
  mpi_procid2tile(PROCID,&tilei,&tilej,&tilek);
  mpi_tileorigin(tilei,tilej,tilek,&toi,&toj,&tok);
  *lix = gix - toi;
  *liy = giy - toj;
  *liz = giz - tok;
  #else
  *lix = gix;
  *liy = giy;
  *liz = giz;
  #endif
}


///////////////////////////////////////////////////////////////
void
mpi_local2globalidx(int lix,int liy, int liz, int *gix, int *giy, int *giz)
{
  #ifdef MPI
  int tilei,tilej,tilek;
  int toi,toj,tok;
  mpi_procid2tile(PROCID,&tilei,&tilej,&tilek);
  mpi_tileorigin(tilei,tilej,tilek,&toi,&toj,&tok);
  *gix = lix + toi;
  *giy = liy + toj;
  *giz = liz + tok;
  #else
  *gix = lix;
  *giy = liy;
  *giz = liz;
  #endif
}


///////////////////////////////////////////////////////////////
void
mpi_procid2tile(int procid, int* tilei, int* tilej, int* tilek)
{
  *tilek = floor(procid / (NTX * NTY));
  *tilej = floor((procid - (*tilek) * NTX * NTY) / NTX);
  *tilei = procid - NTX * NTY * (*tilek) - NTX * (*tilej);
}


///////////////////////////////////////////////////////////////
int
mpi_tile2procid(int tilei, int tilej, int tilek)
{
  return tilek * NTX * NTY + tilej * NTX + tilei;
}


///////////////////////////////////////////////////////////////
//parasitic openMP
int
omp_myinit()
{
#ifdef OMP
#ifdef MPI
  printf("MPI does not work with OMP.\n"); exit(-1);
#endif

  NPROCS=omp_get_num_threads();
  PROCID=0;
  TI=TJ=TK=0;
  TOI=TOJ=TOK=0;
  PROCID=0;  

#else
  NPROCS=1;
  TI=TJ=TK=0;
  TOI=TOJ=TOK=0;
  PROCID=0;  
#endif

  return 0;
}

///////////////////////////////////////////////////////////////
//get averages for unique mpi tiles
int
calc_avgs_throughout()
{

  /***************************/
  // scale height at each radius
#ifdef CALCHRONTHEGO 
  
  int ix,iy,iz,gix,giy,giz;
  struct geometry geom, geomBL;
  ldouble sigma,scaleth,xxBL[4];
  
#pragma omp parallel for
  for(gix=0;gix<TNX;gix++)
  {
    sigma_otg_temp[gix] = 0.;
    scaleth_otg_temp[gix] = 0.;
  }

#pragma omp parallel for private(iy, iz, gix, sigma, scaleth, geom, xxBL)
  for(ix=0;ix<NX;ix++)
  {
    sigma=scaleth=0.;
    for(iy=0;iy<NY;iy++)
    {
      for(iz=0;iz<NZ;iz++) // todo: does not work with multiple z tiles yet!
      {
        fill_geometry(ix,iy,iz,&geom);
	
#ifdef PRECOMPUTE_MY2OUT
        get_xxout(ix, iy, iz, xxBL);
#else
        coco_N(geom.xxvec,xxBL,MYCOORDS,BLCOORDS);
#endif
	sigma+=get_u(p,RHO,ix,iy,iz)*geom.gdet;
        scaleth+=get_u(p,RHO,ix,iy,iz)*geom.gdet*(M_PI/2. - xxBL[2])*(M_PI/2. - xxBL[2]);
      }
    }
    gix=ix+TOI;
    
    sigma_otg_temp[gix]=sigma;
    scaleth_otg_temp[gix]=scaleth;
  }

#ifdef MPI
  MPI_Allreduce(sigma_otg_temp, sigma_otg, TNX, MPI_LDOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
  MPI_Allreduce(scaleth_otg_temp, scaleth_otg, TNX, MPI_LDOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
#else
#pragma omp parallel for private(gix)
  for(ix=0;ix<NX;ix++)
  {
    gix=ix+TOI;
    sigma_otg[gix]=sigma_otg_temp[gix];
    scaleth_otg[gix]=scaleth_otg_temp[gix];
  }
#endif
  
#pragma omp parallel for private(gix)
  for(ix=-NGCX;ix<NX+NGCX;ix++)
  {
    gix=ix+TOI;
    if(gix>0 && gix<TNX)
      scaleth_otg[gix]=sqrt(scaleth_otg[gix]/sigma_otg[gix]);
    
  }
#endif  // ifdef CALCHRONTHEGO
  /***************************/

  /***************************/
  //average velocities 
#ifdef CORRECT_POLARAXIS_3D
#ifdef POLARAXISAVGIN3D

#ifdef PHIWEDGE
 if(PHIWEDGE<0.999*2.*M_PI)
   my_warning("POLARAXISAVGIN3D requires full 2 pi in azimuth!\n");
#endif
 
 int ix,iy,iz,gix,giy,giz,iv; 
 struct geometry geom, geomBL;
 ldouble v[4],ucon[4],r,th,ph;

#pragma omp parallel for private(iv)  // not tested
 for(gix=0;gix<TNX;gix++)
 {
   for(iv=0;iv<NV+2;iv++)
   {
     axis1_primplus_temp[iv][gix] = 0.;
     axis1_primplus[iv][gix] = 0.;
     axis2_primplus_temp[iv][gix] = 0.;
     axis2_primplus[iv][gix] = 0.;
   }
 }

#pragma omp parallel for private(gix, iz, iv, geom, geomBL, v)  // not tested
 for(ix=0;ix<NX;ix++)
 {
   gix=ix+TOI;
   for(iz=0;iz<NZ;iz++)
   {
#ifdef MPI
     if(TJ==0) //topmost tile
#endif
     {
       fill_geometry(ix,NCCORRECTPOLAR,iz,&geom);
       fill_geometry_arb(ix,NCCORRECTPOLAR,iz,&geomBL,BLCOORDS);
       
       axis1_primplus_temp[RHO][gix]+=get_u(p,RHO,ix,NCCORRECTPOLAR,iz);
       axis1_primplus_temp[UU][gix]+=get_u(p,UU,ix,NCCORRECTPOLAR,iz);

       //cartesian velocities in VX..VZ slots
       decompose_vels(&get_u(p,0,ix,NCCORRECTPOLAR,iz),VX,v,&geom,&geomBL);
       axis1_primplus_temp[VX][gix]+=v[1];
       axis1_primplus_temp[VY][gix]+=v[2];
       axis1_primplus_temp[VZ][gix]+=v[3];

       //angular velocity in NV slot!!!
       axis1_primplus_temp[NV][gix]+=get_u(p,VZ,ix,NCCORRECTPOLAR,iz);
       
#ifdef RADIATION
       axis1_primplus_temp[EE][gix]+=get_u(p,EE,ix,NCCORRECTPOLAR,iz);
#ifdef EVOLVEPHOTONNUMBER
       axis1_primplus_temp[NF][gix]+=get_u(p,NF,ix,NCCORRECTPOLAR,iz);
#endif
       //cartesian rad velocities in FX..FZ slots
       decompose_vels(&get_u(p,0,ix,NCCORRECTPOLAR,iz),FX,v,&geom,&geomBL);
       axis1_primplus_temp[FX][gix]+=v[1];
       axis1_primplus_temp[FY][gix]+=v[2];
       axis1_primplus_temp[FZ][gix]+=v[3];

       //rad angular velocity in NV+1 slot!!!
       axis1_primplus_temp[NV+1][gix]+=get_u(p,FZ,ix,NCCORRECTPOLAR,iz);
#endif //RADIATION
     }
     
#ifdef MPI
     if(TJ==NTY-1) //bottommost tile
#endif
     {
       fill_geometry(ix,NY-NCCORRECTPOLAR-1,iz,&geom);
       fill_geometry_arb(ix,NY-NCCORRECTPOLAR-1,iz,&geomBL,BLCOORDS);
       
       axis2_primplus_temp[RHO][gix]+=get_u(p,RHO,ix,NY-NCCORRECTPOLAR-1,iz);
       axis2_primplus_temp[UU][gix]+=get_u(p,UU,ix,NY-NCCORRECTPOLAR-1,iz);
       
       //cartesian velocities in VX..VZ slots
       decompose_vels(&get_u(p,0,ix,NY-NCCORRECTPOLAR-1,iz),VX,v,&geom,&geomBL);
       axis2_primplus_temp[VX][gix]+=v[1];
       axis2_primplus_temp[VY][gix]+=v[2];
       axis2_primplus_temp[VZ][gix]+=v[3];

       //angular velocity in NV slot!!!
       axis2_primplus_temp[NV][gix]+=get_u(p,VZ,ix,NY-NCCORRECTPOLAR-1,iz);
       
#ifdef RADIATION
       axis2_primplus_temp[EE][gix]+=get_u(p,EE,ix,NY-NCCORRECTPOLAR-1,iz);
#ifdef EVOLVEPHOTONNUMBER
       axis2_primplus_temp[NF][gix]+=get_u(p,NF,ix,NY-NCCORRECTPOLAR-1,iz);
#endif 

       //cartesian rad velocities in FX..FZ slots
       decompose_vels(&get_u(p,0,ix,NY-NCCORRECTPOLAR-1,iz),FX,v,&geom,&geomBL);
       axis2_primplus_temp[FX][gix]+=v[1];
       axis2_primplus_temp[FY][gix]+=v[2];
       axis2_primplus_temp[FZ][gix]+=v[3];

       //rad angular velocity in NV+1 slot!!!
       axis2_primplus_temp[NV+1][gix]+=get_u(p,FZ,ix,NY-NCCORRECTPOLAR-1,iz);
#endif //RADIATION
     }
   }
   for(iv=0;iv<NV+2;iv++)
   {
     axis1_primplus_temp[iv][gix]/=NZ;
     axis2_primplus_temp[iv][gix]/=NZ;
   }
 }

#ifdef MPI
 MPI_Allreduce(&axis1_primplus_temp[0][0], &axis1_primplus[0][0], TNX*(NV+2), MPI_LDOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
 MPI_Allreduce(&axis2_primplus_temp[0][0], &axis2_primplus[0][0], TNX*(NV+2), MPI_LDOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
#else 
#pragma omp parallel for private(gix, iv)  // not tested
  for(ix=0;ix<NX;ix++)
  {
    gix=ix+TOI;
    for(iv=0;iv<NV+2;iv++)
    {
      axis1_primplus[iv][gix]=axis1_primplus_temp[iv][gix];
      axis2_primplus[iv][gix]=axis2_primplus_temp[iv][gix];
    }
  }
#endif
  
#endif //POLARAXISAVGIN3D
#endif  // CORRECT_POLARAXIS_3D

  
  /***************************/

return 0;
}
