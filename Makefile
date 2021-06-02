#by default parallel, use 'make SERIAL=1' for serial
SERIAL=1
RM=/bin/rm

ifneq ($(SERIAL),1)
CC=mpicc
CFLAGS=-O3 -DMPI

else
CC=gcc
//CC=/usr/bin/h5cc
//CFLAGS = -fopenmp -O2 -I/usr/include/hdf5/serial
//CFLAGS = -O3 -fPIC -I/home/achael/software/include
CFLAGS = -O3 -fPIC 
OMPFLAGS = -fopenmp
endif

//LIBS=-lm -lgsl -lgslcblas -lsiloh5 -lfftw3 -lrt -lhdf5_serial -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0/include
//LIBS=-lm -lgsl -lgslcblas -lfftw3 -lrt -L/home/achael/software/lib 
LIBS=-lm -lgsl -lgslcblas -lfftw3 -lrt 

//OBJS = mpi.o u2prad.o magn.o silo.o postproc.o fileop.o misc.o physics.o finite.o problem.o metric.o relele.o rad.o opacities.o u2p.o frames.o p2u.o nonthermal.o
//OBJS = mpi.o problem.o finite.o metric.o frames.o relele.o u2p.o p2u.o magn.o physics.o opacities.o misc.o postproc.o fileop.o silo.o
OBJS = mpi.o problem.o finite.o metric.o frames.o relele.o u2p.o p2u.o magn.o physics.o opacities.o misc.o postproc.o fileop.o 

all: ko_gpu ko ana avg phisli thsli phiavg



ko_gpu: ko.o $(OBJS) finitegpu.cu Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -c $(OBJS)
	nvcc -arch=sm_80 --compiler-options '$(CFLAGS)' -x cu -c finitegpu.cu
	nvcc -arch=sm_80 $(LIBS) $(OBJS) finitegpu.o -o ko_gpu ko.o

ko: ko.o $(OBJS) Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o ko ko.o $(OBJS) $(LIBS)

ana: ana.o $(OBJS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o ana ana.o $(OBJS) $(LIBS)
avg: avg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o avg avg.o $(OBJS) $(LIBS)
phisli: phisli.o $(OBJS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o phisli phisli.o $(OBJS) $(LIBS)
thsli: thsli.o $(OBJS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o thsli thsli.o $(OBJS) $(LIBS)
phiavg: phiavg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o phiavg phiavg.o $(OBJS) $(LIBS)
clean:
	$(RM) -f ko ko_gpu ana avg phiavg phisli thsli *~ *.o *.oo
