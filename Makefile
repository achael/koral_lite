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
CFLAGSGPU = -O3 -fPIC -DGPUKO

OMPFLAGS = -fopenmp
endif

//LIBS=-lm -lgsl -lgslcblas -lsiloh5 -lfftw3 -lrt -lhdf5_serial -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0/include
//LIBS=-lm -lgsl -lgslcblas -lfftw3 -lrt -L/home/achael/software/lib 
LIBS=-lm -lgsl -lgslcblas -lfftw3 -lrt 

//OBJS = mpi.o u2prad.o magn.o silo.o postproc.o fileop.o misc.o physics.o finite.o problem.o metric.o relele.o rad.o opacities.o u2p.o frames.o p2u.o nonthermal.o
//OBJS = mpi.o problem.o finite.o metric.o frames.o relele.o u2p.o p2u.o magn.o physics.o opacities.o misc.o postproc.o fileop.o silo.o
OBJS = mpi.o problem.o finite.o metric.o frames.o relele.o u2p.o p2u.o magn.o physics.o opacities.o misc.o postproc.o fileop.o 
OBJSGPU = relelegpu.o finitegpu.o
SRCS = mpi.c problem.c finite.c metric.c frames.c relele.c u2p.c p2u.c magn.c physics.c opacities.c misc.c postproc.c fileop.c 

all: ko_gpu ko ana avg phisli thsli phiavg



ko_gpu: ko.o $(SRCS) relelegpu.cu finitegpu.cu Makefile ko.h kogpu.h problem.h mnemonics.h 
	$(CC) $(CFLAGSGPU) -c $(SRCS)
	nvcc -rdc=true -gencode arch=compute_80,code=sm_80 --compiler-options '$(CFLAGSGPU)' -x cu -dc relelegpu.cu finitegpu.cu
	nvcc -rdc=true -gencode arch=compute_80,code=sm_80 -lcudart $(LIBS) $(OBJS) relelegpu.cu finitegpu.cu -o ko_gpu ko.o

ko: ko.o $(SRCS) Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o ko ko.o $(SRCS) $(LIBS)

ana: ana.o $(SRCS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o ana ana.o $(SRCS) $(LIBS)

avg: avg.o $(SRCS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o avg avg.o $(SRCS) $(LIBS)

phisli: phisli.o $(SRCS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o phisli phisli.o $(SRCS) $(LIBS)

thsli: thsli.o $(SRCS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o thsli thsli.o $(SRCS) $(LIBS)

phiavg: phiavg.o $(SRCS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) $(OMPFLAGS) -o phiavg phiavg.o $(SRCS) $(LIBS)

clean:
	$(RM) -f ko ko_gpu ana avg phiavg phisli thsli *~ *.o *.oo
