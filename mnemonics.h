/*! \file ko.h
 \brief mnemonics for various indices
*/

//rad or hydro
#define RAD 1
#define MHD 2 

//electrons or ions
#define GAS 0
#define IONS 1
#define ELECTRONS 2

//velocities
#define VEL4 1 // lab four-velocity u^i
#define VEL3 2 // lab three-velocity u^i/u^t
#define VELR 3 // relative velocity \tilde u^i

//primitive/conserved ordering
#define RHO 0
#define UU 1
#define VX 2
#define VY 3
#define VZ 4
#define ENTR 5

#define ENTRE 6
#define ENTRI 7

#define NEREL(i) (8+i)

#define B1 (NVHD+0)
#define B2 (NVHD+1)
#define B3 (NVHD+2)

#define EE (NVMHD)
#define FX (NVMHD+1)
#define FY (NVMHD+2)
#define FZ (NVMHD+3)

#if defined(RADIATION) && defined(EVOLVEPHOTONNUMBER)
#define NF (NVMHD+4)
#else
#define NF 0
#endif

//single fluid macros
#define EE0 EE
#define FX0 FX
#define FY0 FY
#define FZ0 FZ
#define NF0 NF

//radiative closures
#define M1CLOSURE 0
#define EDDCLOSURE 1 // no longer used
#define VETCLOSURE 2
//#define M1ORTOCLOSURE 3
//#define MINERBOCLOSURE 4

//u2p inversion types
#define U2P_HOT 0
#define U2P_ENTROPY 1
//#define U2P_HOTMAX 2 // no longer used
//#define U2P_COLD 3
//#define U2P_SLOW 4

#define U2P_EQS_NOBLE 0
#define U2P_EQS_JON 1

#define U2P_SOLVER_WP 0
#define U2P_SOLVER_W 1
//#define U2P_SOLVER_WPPLUS5D 2 // no longer used

//coordinates/metric
#define BLCOORDS 1
#define SCHWCOORDS 1
#define KERRCOORDS 1

#define KSCOORDS 2
#define MKSCOORDS 3
#define MINKCOORDS 4
#define CYLCOORDS 5
#define SPHCOORDS 6
#define MKS1COORDS 7
#define MCYL1COORDS 8
#define MKER1COORDS 9
#define MKS2COORDS 10
#define MSPH1COORDS 11
#define MKS3COORDS 12
#define TKS3COORDS 13
#define TFLATCOORDS 14
#define JETCOORDS 15

//type of boundary
#define XBCLO 1
#define XBCHI 2
#define YBCLO 3
#define YBCHI 4
#define ZBCLO 5
#define ZBCHI 6

//cell flags
#define NFLAGS 7

//boolean
#define ENTROPYFLAG 0
#define RADIMPFIXUPFLAG 1
#define HDFIXUPFLAG 2
#define RADFIXUPFLAG 3
#define ENTROPYFLAG2 4
#define ENTROPYFLAG3 5
#define RADSOURCETYPEFLAG 6

//values for RADFIXUP
#define RADSOURCETYPEEXPLICIT 1
#define RADSOURCETYPEEXPLICITSUBSTEP 2
#define RADSOURCETYPEIMPLICITLAB 3
#define RADSOURCETYPEIMPLICITFF 10

//fluxes
#define LAXF_FLUX 0
#define HLL_FLUX 1
#define HLLC_FLUX 2

//timestepping
#define RK2IMEX 1
#define RK2 0
#define RK1 -1
#define RK2HEUN 2

//types of hd/rad viscosity
#define NOVISCOSITY 0
#define SIMPLEVISCOSITY 1
#define SHEARVISCOSITY 2

//rad.implicit solver parameters
#define RADIMPLICIT_ENERGYEQ 0
#define RADIMPLICIT_ENTROPYEQ 1
#define RADIMPLICIT_LTEEQ 2

#define RADIMPLICIT_LAB 0
#define RADIMPLICIT_FF 1

//global integer slots
#define NGLOBALINTSLOT 18
#define GLOBALINTSLOT_NIMPENERRADFF 13
#define GLOBALINTSLOT_NIMPENERMHDFF 14
#define GLOBALINTSLOT_NIMPENERRAD 0
#define GLOBALINTSLOT_ITERIMPENERRAD 1
#define GLOBALINTSLOT_NIMPENERMHD 2
#define GLOBALINTSLOT_ITERIMPENERMHD 3
#define GLOBALINTSLOT_NIMPENTRRAD 4
#define GLOBALINTSLOT_ITERIMPENTRRAD 5
#define GLOBALINTSLOT_NIMPENTRMHD 6
#define GLOBALINTSLOT_ITERIMPENTRMHD 7
#define GLOBALINTSLOT_NIMPLTE 8
#define GLOBALINTSLOT_ITERIMPLTE 9
#define GLOBALINTSLOT_NRADFIXUPS 10
#define GLOBALINTSLOT_NCRITFAILURES 11
#define GLOBALINTSLOT_NTOTALRADIMPFAILURES 12
#define GLOBALINTSLOT_NTOTALMHDFIXUPS 15
#define GLOBALINTSLOT_NTOTALRADFIXUPS 16
#define GLOBALINTSLOT_NTOTALRADIMPFIXUPS 17

//frames
#define ZAMOFRAME 0
#define FFFRAME 1

//avg quantities
#define AVGBSQ (NV+0)
#define AVGUCON(i) (NV+1+i)
#define AVGUCOV(i) (NV+5+i)
#define AVGBCON(i) (NV+9+i)
#define AVGBCOV(i) (NV+13+i)
#define AVGRHOUCON(i) (NV+17+i)
#define AVGRHOUCOV(i) (NV+21+i)
#define AVGUUUCON(i) (NV+25+i)
#define AVGUUCOV(i) (NV+29+i)
#define AVGBSQUCON(i) (NV+33+i)
#define AVGBSQUCOV(i) (NV+37+i)
#define AVGRHOUCONUCOV(i,j) (NV+41+i*4+j)
#define AVGUUUCONUCOV(i,j) (NV+57+i*4+j)
#define AVGBSQUCONUCOV(i,j) (NV+73+i*4+j)
#define AVGBCONBCOV(i,j) (NV+89+i*4+j)
#define AVGWUCON(i) (NV+105+i)
#define AVGFLUXXL(i) (NV+109+i)
#define AVGFLUXYL(i) (NV+109+NV+i)
#define AVGFLUXZL(i) (NV+109+2*NV+i)
#define AVGRHOURDIFF (NV+109+3*NV)
#define AVGTGAS (NV+110+3*NV)
#define AVGEHAT (NV+3*NV+111)
#define AVGRIJ(i,j) (4*NV+112+i*4+j)
#define AVGEHATUCON(i) (4*NV+128+i)
#define AVGEHATUCOV(i) (4*NV+132+i)
#define AVGURFCON(i) (4*NV+136+i)
#define AVGURFCOV(i) (4*NV+140+i)
#define AVGGHAT(i) (4*NV+144+i)
#define AVGGHATCOMPT(i) (4*NV+148+i)
#define AVGNFHAT (4*NV+152)
#define AVGPGAS (4*NV+153)
#define AVGPE (4*NV+154)
#define AVGPI (4*NV+155)
#define AVGTIJ(i,j) (4*NV+156+i*4+j)
#define AVGVISCHEATING (4*NV+156+16)
#define AVGVISCHEATINGNEGE (4*NV+157+16)
#define AVGVISCHEATINGNEGI (4*NV+158+16)
#define AVGVISCHEATINGTIMESDELTAE (4*NV+159+16)
#define AVGNRELEL (4*NV+176)
#define AVGPRELEL (4*NV+177)
#define AVGURELEL (4*NV+178)
#define AVGNETH (4*NV+179)

//MPI mnemonics
#define MPI_MSG_XLO 100
#define MPI_MSG_XHI 101
#define MPI_MSG_YLO 102
#define MPI_MSG_YHI 103
#define MPI_MSG_ZLO 104
#define MPI_MSG_ZHI 105
#define MPI_MSG_XLOYLO 200
#define MPI_MSG_XLOYHI 201
#define MPI_MSG_XHIYLO 202
#define MPI_MSG_XHIYHI 203
#define MPI_MSG_XLOZLO 204
#define MPI_MSG_XLOZHI 205
#define MPI_MSG_XHIZLO 206
#define MPI_MSG_XHIZHI 207
#define MPI_MSG_YLOZLO 208
#define MPI_MSG_YLOZHI 209
#define MPI_MSG_YHIZLO 210
#define MPI_MSG_YHIZHI 211
#define MPI_MSG_XLOYLOZLO 301
#define MPI_MSG_XLOYLOZHI 302
#define MPI_MSG_XLOYHIZLO 303
#define MPI_MSG_XHIYLOZLO 304
#define MPI_MSG_XLOYHIZHI 305
#define MPI_MSG_XHIYLOZHI 306
#define MPI_MSG_XHIYHIZLO 307
#define MPI_MSG_XHIYHIZHI 308

//fixup types
#define FIXUP_U2PRAD 100
#define FIXUP_U2PMHD 101
#define FIXUP_RADIMP 102

//electron heating types
#define ELECTRONIONHEATTYPE_THROUGHUINT 1
#define ELECTRONIONHEATTYPE_THROUGHENTROPY 2
