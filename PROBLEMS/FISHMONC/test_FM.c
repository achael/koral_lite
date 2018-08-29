#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define ldouble double
#define FTYPE ldouble
#define BHSPIN 0.9375
#define GAMMA (4./3.)
#define FM_rin 6.
#define FM_rmax 12.
#define FM_rho0 1.


#include "tools.c"


/////////////////////////////////////////////////////////////////

int
main(int argc, char **argv)
{
  int i;
  ldouble r, th = 3.14159265*0.5, a = BHSPIN, rho, uu, ell;
  
  for(i = 0; i < 101; i++)
  {
    r = pow(10., 0.1 + 0.02*i);
    init_dsandvels_fishbone_moncrief(r, th, a, &rho, &uu, &ell);
    printf("%d %e %e %e %e\n", i, r, rho, uu, ell);
  }
  
  return 0;
}


