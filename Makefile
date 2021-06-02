#by default parallel, use 'make SERIAL=1' for serial
SERIAL=1
ifneq ($(SERIAL),1)
CC=mpicc
CFLAGS=-O3 -DMPI
else
CC=gcc
//CC=nvcc
//CC=/usr/bin/h5cc
//CFLAGS = -fopenmp -O2 -L/usr/lib/gcc/x86_64-linux-gnu/5.4.0/include -I/usr/include/hdf5/serial
CFLAGS = -O2 -L/home/achael/software/lib -I/home/achael/software/include
//CFLAGS =  -O2 -L/home/ogreen/koral_lite/software/lib -I/home/ogreen/koral_lite/software/include
endif
//LIBS=-lm -lgsl -lgslcblas -lsiloh5 -lfftw3 -lrt -lhdf5_serial
LIBS=-lm -lgsl -lgslcblas -lsilo -lfftw3 -lrt
RM=/bin/rm
//OBJS = mpi.o u2prad.o magn.o silo.o postproc.o fileop.o misc.o physics.o finite.o problem.o metric.o relele.o rad.o opacities.o u2p.o frames.o p2u.o nonthermal.o
//OBJS = mpi.o problem.o finite.o metric.o frames.o relele.o u2p.o p2u.o magn.o physics.o opacities.o misc.o postproc.o fileop.o silo.o
OBJS = mpi.o problem.o metric.o frames.o relele.o u2p.o p2u.o magn.o physics.o opacities.o misc.o postproc.o fileop.o silo.o
SOURCES		= mpi.c problem.c metric.c frames.c relele.c u2p.c p2u.c magn.c physics.c opacities.c misc.c postproc.c fileop.c silo.c

#nvcc –arch=sm_80 –dlink $(OBJS) –o gpuCode.o
#nvcc $(CFLAGS) -x cu -ccbin=gcc $(SOURCES) finite.cu $(LIBS)

gpu:
	$(CC) $(CFLAGS) -c ko.c $(SOURCES) $(LIBS)
	nvcc --compiler-options '-fPIC' -o finiteCU.o  -c finite.cu
	nvcc -arch=sm_80 -dlink finiteCU.o -o gpucode.o $(OBJS)
	$(CC) gpucode.o $(OBJS) -o kogpu ko.o 
        #$(CC) -shared -o libKO.so $(OBJS)
        #$(CC) $(CFLAGS) gpucode.o -I/usr/local/cuda/include -Wall -g -L. -I. $(LIBS) -L/usr/local/cuda/lib64 ko.c -o main $(OBJS) finiteCU.o -lcudart

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
