#by default parallel, use 'make SERIAL=1' for serial

LANG		= -std=cc14
CXX		= gcc
OPTS		= -lm -O3 -fopenmp
CFLAGS          = -fPIC -O2 -L/home/achael/software/lib -I/home/achael/software/include
SOURCES		= mpi.c problem.c metric.c frames.c relele.c u2p.c p2u.c magn.c physics.c opacities.c misc.c postproc.c fileop.c silo.c
OBJECTS		= $(SOURCES:.c=.o)


LIBS=-lm -lgsl -lgslcblas -lsilo -lfftw3 -lrt

RM=/bin/rm

//OBJS = mpi.o problem.o metric.o frames.o relele.o u2p.o p2u.o magn.o physics.o opacities.o misc.o postproc.o fileop.o silo.o

all: cu_exec

cu_exec: $(SOURCES) finite.cu nvcc $(LIBS) -O3 -std=cc11  -fopenmp $^ -o $@

clean:
	$(RM) -f *.o  cu_exec

