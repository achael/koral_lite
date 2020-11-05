int if_indisk(int ix, int iy, int iz);

/**********************/
int ret_val = 0;
int gix,giy,giz;
gix = ix + TOI;
giy = iy + TOJ;
giz = iz + TOK;

ret_val = if_indisk(gix,giy,giz);

if(ret_val == 1) 
{
  BCtype = XBCHI;
}
else if(ret_val == 2) //no longer used
{ 
  BCtype = YBCHI;
}
else if(ret_val == 3) //no longer used
{
  BCtype = YBCLO;
}


//if(ret_val != 0 && gix == 144 && giy == 111) printf("BC CHECK - ix iy iz gix giy giz BC : %i %i %i %i %i %i %i \n",ix,iy,iz,gix,giy,giz,BCtype);
