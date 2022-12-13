//returns problem specific BC
//ix, iy, iz hold the ghost cell indices, e.g., (-1,0,0)
//BCtype gives the type of the boundary

int iix,iiy,iiz,iv;

/***********************************************/
//structure of geometry
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);


struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,BLCOORDS);
