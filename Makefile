
CC=gcc
#CC=pgcc
CFLAGS=-O3 -fopenmp -w -lrt -L/home/maxadmin/silo-4.10.2/lib -I/home/maxadmin/silo-4.10.2/include
#CFLAGS=CFLAGS=-O3 -mp -w -lrt -L/home/maxadmin/silo-4.10.2/lib -I/home/maxadmin/silo-4.10.2/include
LIBS=-lm -lgsl -lgslcblas -lsilo
RM=/bin/rm

OBJS = mpi.o u2prad.o magn.o silo.o postproc.o fileop.o misc.o physics.o finite.o problem.o metric.o relele.o rad.o opacities.o u2p.o frames.o p2u.o nonthermal.o


all: ko 

ko: ko.o $(OBJS) Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ko ko.o $(OBJS) $(LIBS)

ana: ana.o $(OBJS)  Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o ana ana.o $(OBJS) $(LIBS)

avg: avg.o $(OBJS) Makefile ko.h problem.h mnemonics.h 
	$(CC) $(CFLAGS) -o avg avg.o $(OBJS) $(LIBS)


phiavg: phiavg.o
	$(CC) $(CFLAGS) -o phiavg phiavg.o $(LIBS)

phisli: phisli.o
	$(CC) $(CFLAGS) -o phisli phisli.o $(LIBS)

thsli: thsli.o
	$(CC) $(CFLAGS) -o thsli thsli.o $(LIBS)

regrid: regrid.o $(OBJS)  Makefile ko.h problem.h mnemonics.h
	$(CC) $(CFLAGS) -o regrid regrid.o $(OBJS) $(LIBS)

clean:
	$(RM) -f ko ana $(OBJS) *~ *.o *.oo
