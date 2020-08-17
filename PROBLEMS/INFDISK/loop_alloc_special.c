int if_indisk(int ix, int iy, int iz);

/**********************/
int ret_val = 0;
int gix,giy,giz;
gix = ix + TOI;
giy = iy + TOJ;
giz = iz + TOK;

ret_val = if_indisk(gix,giy,giz);

//if(ret_val != 0 && ret_val != 4) printf("LOOP ALLOC - ix iy iz gix giy giz rv : %i %i %i %i %i %i %i \n\n",ix,iy,iz,gix,giy,giz,ret_val);

//if(ret_val == 4) printf("LOOP ALLOC (Corner) - gix giy giz rv : %i %i %i %i \n\n",gix,giy,giz,ret_val);
