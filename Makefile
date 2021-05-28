#by default parallel, use 'make SERIAL=1' for serial

ifneq ($(SERIAL),1)
CC=mpicc 
CFLAGS=-O3 -DMPI

else
CC=gcc
//CC=/usr/bin/h5cc
CFLAGS = -fopenmp -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0/include -I/usr/include/hdf5/serial  
//CFLAGS = -fopenmp -O2 -L/home/achael/software/lib -I/home/achael/software/include

endif

LIBS=-lm -lgsl -lgslcblas -lsiloh5 -lfftw3 -lrt -lhdf5_serial
//LIBS=-lm -lgsl -lgslcblas -lsilo -lfftw3 -lrt 

RM=/bin/rm
//OBJS = mpi.o u2prad.o magn.o silo.o postproc.o fileop.o misc.o physics.o finite.o problem.o metric.o relele.o rad.o opacities.o u2p.o frames.o p2u.o nonthermal.o 
OBJS = mpi.o problem.o finite.o metric.o frames.o relele.o u2p.o p2u.o magn.o physics.o opacities.o misc.o postproc.o fileop.o silo.o


all: ko ana avg phisli thsli phiavg 

ko: ko.o $(OBJS) Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ko ko.o $(OBJS) $(LIBS)

ana: ana.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ana ana.o $(OBJS) $(LIBS)

avg: avg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o avg avg.o $(OBJS) $(LIBS)

phisli: phisli.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o phisli phisli.o $(OBJS) $(LIBS)

thsli: thsli.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o thsli thsli.o $(OBJS) $(LIBS)

phiavg: phiavg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o phiavg phiavg.o $(OBJS) $(LIBS)

clean:
	$(RM) -f ko ana avg phiavg phisli thsli *~ *.o *.oo
