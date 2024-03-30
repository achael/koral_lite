//int
//convert_out2gif_1d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

/****************************************/
/****************************************/
/****************************************/
   fprintf(fgnu,
	  "set style line 1 lw 2 lc 1 lt 1\n"
	  "set style line 2 lw 2 lc 2 lt 1 pt 7 ps .5 \n"
	  "set style line 3 lw 2 lc 3 lt 1 pt 7 ps .5 \n"
	  "set style line 4 lw 2 lc 4 lt 1 pt 7 ps .5 \n"
	  "set style line 10 lw 2 lc 3 lt 3\n"
	  "set style line 11 lw 2 lc 4 lt 3\n"
	  "set style line 12 lw 2 lc 5 lt 3\n"
	  "set style line 13 lw 2 lc 0 lt 3\n"
	  "set style line 14 lw 2 lc 7 lt 3\n"
	  "set style line 100 lw 2 lc -1 lt 3\n"
	  "set style line 101 lw 1 lc -1 lt 3\n"
	  "set term gif large size 1200,600\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set autoscale\n" 
	   "set log\n"
#ifdef RADIATION
	  "set label \"t=%.2e (%.2e s) (Prad/Pgas)_out=%.2e (Mdot/Mdot_Edd)_out=%.2e Rmin=%.0f Rmax=%.0f RBondi_init=%.0f\" at screen .13, .98\n"
#else
	  "set label \"t=%.2e (%.2e s) (Prad/Pgas)_out=0*%.2e (Mdot/Mdot_Edd)_out=%.2e Rmin=%.0f Rmax=%.0f RBondi_init=%.0f\" at screen .13, .98\n"
#endif

	  "set lmargin at screen 0.07\n"
	  "set rmargin at screen 0.33\n"
	  "set bmargin at screen .07\n"
	  "set tmargin at screen .42\n"
	  "set autoscale\n"
	  "set xrange [%f:%f]\n"
	  //	  "unset log y\n"
	  "set format x \"%%.1f\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
#ifdef RADIATION
	  "plot \"%s\" u 1:(($29+1.e-20)**2)**.5 w lp ls 4 ti \"luminosity\"\n"
#else
	  "plot \"%s\" u 1:($14) w lp ls 4 ti \"nothing\"\n"
#endif


	  "set lmargin at screen 0.40\n"
	  "set rmargin at screen 0.66\n"
	  "set bmargin at screen .07\n"
	  "set tmargin at screen .42\n"
	  //	  "unset log y\n"
	  "set format x \"%%.1f\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  //	  "plot \"%s\" u 1:27 w lp ls 2 pt 7 ps .5  ti \"tau_abs\", \"%s\" u 1:26 w lp ls 3 pt 7 ps .5  ti \"tau_tot\"\n"
#if defined(RADIATION) && defined(NCOMPTONIZATION)
	   "plot \"%s\" u 1:($32) w lp ls 2 ti \"nph\"\n"
#else
	  "plot \"%s\" u 1:(1) w lp ls 2 pt 7 ps .5  ti \"nothing\"\n"
#endif

	  "set lmargin at screen 0.73\n"
	  "set rmargin at screen 0.99\n"
	  "set bmargin at screen .07\n"
	  "set tmargin at screen .42\n"
	  //	  "unset log y\n"
	  "set format x \"%%.1f\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	   "set log y\n"
	   "plot \"%s\" u 1:(-4.*3.1416*$14*$1*$1*$16*%e) w lp ls 3   ti \"Mdot/Mdot_init\","
	   " \"%s\" u 1:(-$16/$28) w lp ls 4 ti \"Mach number\"\n"
	   "set log y\n"

	  "set lmargin at screen 0.07\n"
	  "set rmargin at screen 0.33\n"
	  "set bmargin at screen .5\n"
	  "set tmargin at screen .9\n"
	  //	  "unset log y\n"
	  "set format x \"\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  //	  "set log y\n"
#ifdef RADIATION
	  "plot \"%s\" u 1:($14) w lp ls 4 pt 7 ps .5 ti \"rho\","
" \"%s\" u 1:15 w lp ls 3 ti \"u_int\", "
" \"%s\" u 1:20 w lp ls 2 ti \"E_rad\", "
"\"dumps/out0000.dat\" u 1:($14) w l ls 100 ti \"rho hydro Bondi\","
"\"dumps/out0000.dat\" u 1:($15) w l ls 101 lt 1 ti \"u_int hydro Bondi\"\n"
#else
	  "plot \"%s\" u 1:($14) w lp ls 4 pt 7 ps .5 ti \"rho\","
" \"%s\" u 1:15 w lp ls 3 ti \"u_int\", "
" \"%s\" u 1:20 w lp ls 2 ti \"u_int\", "
"\"dumps/out0000.dat\" u 1:($14) w l ls 100 ti \"rho hydro Bondi\","
"\"dumps/out0000.dat\" u 1:($15) w l ls 101 lt 1 ti \"u_int hydro Bondi\"\n"
#endif
	  //	  "unset log y\n"

	  "set lmargin at screen 0.40\n"
	  "set rmargin at screen 0.66\n"
	  "set bmargin at screen .5\n"
	  "set tmargin at screen .9\n"
	  "set format x \"\"\n"
	  "set format y \"%%.1e\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
#ifdef RADIATION
	  "plot \"%s\" u 1:24 w lp ls 4 pt 7 ti \"Tgas\", \"%s\" u 1:25 w lp ls 2 ti \"Trad\", "
"\"dumps/out0000.dat\" u 1:($24) w l ls 100 ti \"T_gas hydro Bondi\"\n"
#else
	  "plot \"%s\" u 1:24 w lp ls 4 pt 7 ti \"Tgas\", \"%s\" u 1:24 w lp ls 4 ti \"Tgas\", "
"\"dumps/out0000.dat\" u 1:($24) w l ls 100 ti \"T_gas hydro Bondi\"\n"
#endif
	  "set lmargin at screen 0.73\n"
	  "set rmargin at screen 0.99\n"
	  "set bmargin at screen .5\n"
	  "set tmargin at screen .9\n"
	  "set format x \"\"\n"
	  "set format y \"%%.1e\"\n" 
	   "set key bottom left\n"
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	   //"unset log y\n"
#ifdef RADIATION
	  "plot \"%s\" u 1:(-$16) w lp ls 4 pt 7 ti \"gas -vr\",  "
	   " \"%s\" u 1:(-($21)) w lp ls 2   ti \"radiation -vr\", "
" \"%s\" u 1:(($21)) w lp ls 2 lc 3  ti \"radiation vr\", "
" \"dumps/out0000.dat\" u 1:(-$16) w l ls 100 ti \"gas -vr hydro Bondi\" "
#else
	  "plot \"%s\" u 1:(-$16) w lp ls 4 pt 7 ti \"gas vr\", "
" \"%s\" u 1:($16) w lp ls 2 pt 7 ps .5  ti \"gas vr\", "
" \"%s\" u 1:($16) w lp ls 2 pt 7 ps .5  ti \"gas vr\", "
" \"dumps/out0000.dat\" u 1:(-$16) w l ls 100 ti \"gas -vr hydro Bondi\" "
#endif

,fname2,t,t/CCC,0.,MDOT,RMIN,RMAX,0.,exp(get_xb(-NG,0)),exp(get_xb(NX+NG,0)),
fname,fname,fname,
	   rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)/calc_mdotEdd()/MDOT,
fname,fname,fname,fname,fname,fname,fname,fname,fname);





/****************************************/
/****************************************/
/****************************************/

//  fprintf(fgnu,"\n");
//  fclose(fgnu);   
//  
//  int i=system("gnuplot plot.gp ");
//  return 0;
//}
