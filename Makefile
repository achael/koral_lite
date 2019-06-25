#by default parallel, use 'make SERIAL=1' for serial

ifneq ($(SERIAL),1)
CC=mpicc 
CFLAGS=-O3 -DMPI
else
//CC=gcc
//CFLAGS=-O2 -Wno-unused-result -fopenmp
//-fsanitize=address -g -fno-omit-frame-pointer -Wunused-function 
CC=clang
CFLAGS = -O2 -Wno-unused-result -I/usr/lib/gcc/x86_64-linux-gnu/5.4.0/include -I/usr/include/hdf5/serial -Wunused-function -fopenmp=libiomp5 -g
//-fsanitize=address -fno-omit-frame-pointer
endif

LIBS=-lm -lgsl -lgslcblas -lsiloh5 -lfftw3 -lrt -lhdf5_serial
RM=/bin/rm

OBJS = mpi.o u2prad.o magn.o silo.o postproc.o fileop.o misc.o physics.o finite.o problem.o metric.o relele.o rad.o opacities.o u2p.o frames.o p2u.o nonthermal.o 

all: ko ana avg outavg phisli thsli phiavg regrid dumps2hdf5

ko: ko.o $(OBJS) Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ko ko.o $(OBJS) $(LIBS)

ana: ana.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ana ana.o $(OBJS) $(LIBS)

avg: avg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o avg avg.o $(OBJS) $(LIBS)

outavg: outavg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o outavg outavg.o $(OBJS) $(LIBS)

phisli: phisli.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o phisli phisli.o $(OBJS) $(LIBS)

thsli: thsli.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o thsli thsli.o $(OBJS) $(LIBS)

phiavg: phiavg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o phiavg phiavg.o $(OBJS) $(LIBS)

regrid: regrid.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o regrid regrid.o $(OBJS) $(LIBS)

dumps2hdf5: dumps2hdf5.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o dumps2hdf5 dumps2hdf5.o $(OBJS) $(LIBS)

clean:
	$(RM) -f ko ana avg phiavg phisli thsli outavg regrid *~ *.o *.oo
