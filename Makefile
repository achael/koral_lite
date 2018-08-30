#by default parallel, use 'make SERIAL=1' for serial

ifneq ($(SERIAL),1)
CC=mpicc 
CFLAGS=-O3 -DMPI
else
//CC=gcc
//CFLAGS=-O2 -Wno-unused-result -fopenmp
//-fsanitize=address -g -fno-omit-frame-pointer -Wunused-function 
CC=clang
CFLAGS = -O2 -Wno-unused-result -I/usr/lib/gcc/x86_64-linux-gnu/5.4.0/include -fopenmp=libiomp5 -Wunused-function
//-fsanitize=address -g -fno-omit-frame-pointer
endif

LIBS=-lm -lgsl -lgslcblas -lsiloh5 -lfftw3 -lrt
RM=/bin/rm

OBJS = mpi.o u2prad.o magn.o silo.o postproc.o fileop.o misc.o physics.o finite.o problem.o metric.o relele.o rad.o opacities.o u2p.o frames.o p2u.o nonthermal.o

all: ko ana avg 

ko: ko.o $(OBJS) Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ko ko.o $(OBJS) $(LIBS)

ana: ana.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ana ana.o $(OBJS) $(LIBS)

avg: avg.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o avg avg.o $(OBJS) $(LIBS)

clean:
	$(RM) -f ko ana avg *~ *.o *.oo
