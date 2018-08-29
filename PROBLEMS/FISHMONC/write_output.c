// Special purpose routine to write out some relevant quantities in gravitational units, as required for the EHT GRMHD code comparison project

int write_initial()
{
  char bufor[50];
  fout1 = fopen("./PROBLEMS/FISHMONC/initial_equator.dat","w");

  // The following lets the user choose the coordinates in which the output quantities are calculated, but it is not working properly
  
  //int OUTPUTCOORDS = OUTCOORDS;  // use this for the original OUTCOORDS
  int OUTPUTCOORDS = KSCOORDS;  // use this for Kerr-Schild coordinates
  //int OUTPUTCOORDS = BLCOORDS;  // use this for BL coordinates
  
  int ix, iy, iz, iv, iix;
  ldouble pp[NV], phi, tausca;
  
  iy = NY / 2;
  iz = 0;
  for (iix = -2; iix < NX; iix++)
  {
    ix = iix;
    
    struct geometry geom, geomOUTPUT,geomBL;
    fill_geometry(ix, iy, iz, &geom);
    fill_geometry_arb(ix, iy, iz, &geomOUTPUT, OUTPUTCOORDS);
    fill_geometry_arb(ix,iy,iz,&geomBL, BLCOORDS);
    
    ldouble r=geomOUTPUT.xx;
    ldouble th=geomOUTPUT.yy;
    ldouble ph=geomOUTPUT.zz;

    fprintf(fout1,"%d %d %d ", ix, iy, iz);  //(1-3)
    fprintf(fout1,"%.5e %.5e %.5e ", r, th, ph);  //(4-6)
    
    for(iv = 0; iv < NV; iv++)
    {
      pp[iv] = get_u(p, iv, ix, iy, iz);
    }

    trans_pall_coco(pp, pp, MYCOORDS, OUTPUTCOORDS, geom.xxvec, &geom, &geomOUTPUT);
    
    ldouble dxph[3],dx[3],xx1[4],xx2[4];
    xx1[0]=0.;xx1[1]=get_xb(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_x(iz,2);
    xx2[0]=0.;xx2[1]=get_xb(ix+1,0);xx2[2]=get_x(iy,1);xx2[3]=get_x(iz,2);
    coco_N(xx1,xx1,MYCOORDS,OUTPUTCOORDS);
    coco_N(xx2,xx2,MYCOORDS,OUTPUTCOORDS);
    dx[0]=fabs(xx2[1]-xx1[1]);
    xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_xb(iy,1);xx1[3]=get_x(iz,2);
    xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_xb(iy+1,1);xx2[3]=get_x(iz,2);
    coco_N(xx1,xx1,MYCOORDS,OUTPUTCOORDS);
    coco_N(xx2,xx2,MYCOORDS,OUTPUTCOORDS);
    dx[1]=fabs(xx2[2]-xx1[2]);
    xx1[0]=0.;xx1[1]=get_x(ix,0);xx1[2]=get_x(iy,1);xx1[3]=get_xb(iz,2);
    xx2[0]=0.;xx2[1]=get_x(ix,0);xx2[2]=get_x(iy,1);xx2[3]=get_xb(iz+1,2);
    coco_N(xx1,xx1,MYCOORDS,OUTPUTCOORDS);
    coco_N(xx2,xx2,MYCOORDS,OUTPUTCOORDS);
    dx[2]=fabs(xx2[3]-xx1[3]);
    
    dxph[0]=dx[0]*sqrt(geomOUTPUT.gg[1][1]);
    dxph[1]=dx[1]*sqrt(geomOUTPUT.gg[2][2]);
    dxph[2]=dx[2]*sqrt(geomOUTPUT.gg[3][3]);
    
    ldouble gdet=geom.gdet;
    ldouble volume=gdet*get_size_x(ix,0)*get_size_x(iy,1)*get_size_x(iz,2);
    
    fprintf(fout1,"%.5e %.5e ", gdet, volume);  //(7-8)
    
    ldouble rho, uint, pgas, pmag, ptot, betainv, bsq, bcon[4], bcov[4], ucon[4], ucov[4];
    ldouble gamma = GAMMA;
    
    rho = pp[0];
    uint = pp[1];
    
    calc_ucon_ucov_from_prims(pp, &geomOUTPUT, ucon, ucov);
    ldouble lorentz = fabs(ucon[0]) / sqrt(fabs(geomOUTPUT.GG[0][0]));
    //ldouble lorentz2 = fabs(ucon[0]) / sqrt(fabs(geomBL.GG[0][0]));
  
    calc_bcon_prim(pp,bcon,&geomOUTPUT);
    indices_21(bcon,bcov,geomOUTPUT.gg);
    bsq = dot(bcon,bcov);
    
    pgas = (gamma - 1.) * uint;
    pmag = bsq / 2.;
    ptot = pgas + pmag;
    betainv = (pmag + 1.e-90) / pgas;

    fprintf(fout1,"%.5e %.5e ", rho, pgas);  //(9-10)
    fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e ", ucon[0], ucon[1], ucon[2], ucon[3], lorentz); //(11-15)
    fprintf(fout1,"%.5e %.5e %.5e %.5e %.5e ", bsq, bcon[0], bcon[1], bcon[2], bcon[3]); //(16-20)
    fprintf(fout1,"%.5e %.5e %.5e %.5e ", pgas, pmag, ptot, betainv); //(21-24)
    
    fprintf(fout1,"\n");
  }

  fflush(fout1);
  fclose(fout1);
  
  return 0;
}
