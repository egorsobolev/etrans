CC=gcc
MPICC=mpicc

# icc
#CFLAGS= -O3

# gcc
# in GCC 4.9.2 -ftree-loop-vectorize leads to segmentation fault in sol_update (sol.c)
# -freee-loop-vectorize is enabled by -03, so use -O2
# this has no matter on caclulation speed
CFLAGS= -O2 -funroll-all-loops -ffast-math -fno-trapping-math
LFLAGS= -O3

MPICFLAGS= -O2 -funroll-all-loops -ffast-math -fno-trapping-math
MPILFLAGS= -O3

SRC = $(wildcard *.c)
OBJMPI = $(addprefix .o.mpi/,$(SRC:.c=.o))
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:	serial mpi

serial:	argtable librng $(OBJSRL)
	$(CC) -o ../bin/etrans $(OBJSRL) $(LFLAGS) -lm \
	    ../argtable2-13/lib/libargtable2.a ../rng1-0/librng.a

mpi:	argtable librng $(OBJMPI)
	$(MPICC) -o ../bin/etrans.MPI $(OBJMPI) $(MPILFLAGS) -lm \
            ../argtable2-13/lib/libargtable2.a ../rng1-0/librng.a

.o.mpi/%.o:	%.c
	$(MPICC) -I/usr/include/x86_64-linux-gnu $(MPICFLAGS) \
	    -I../argtable2-13/include -I../rng1-0 \
	    -c -o $@ $<

.o.serial/%.o:	%.c
	$(CC) -I/usr/include/x86_64-linux-gnu ${CFLAGS} \
	    -I../argtable2-13/include -I../rng1-0 \
	    -c -o $@ $<

$(OBJMPI):	| .o.mpi

$(OBJSRL):	| .o.serial

.o.mpi:
	mkdir .o.mpi

.o.serial:
	mkdir .o.serial

argtable:	| ../argtable2-13

../argtable2-13: ../argtable2-13.tar.gz
	( \
	    cd ..; \
	    tar -xzf argtable2-13.tar.gz; \
	    cd argtable2-13; \
	    ./configure --prefix=`pwd`; \
	    make install; \
	    cd ../src; \
	)

../argtable2-13.tar.gz:
	( \
	    cd ..; \
	    wget http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz; \
	    cd src; \
	)

librng:		| ../rng1-0

../rng1-0:
	( \
	    cd ..; \
	    tar -xzf rng1-0.tar.gz; \
	    cd rng1-0; \
            make; \
	    cd ../src; \
	)

.PHONY:	clean
clean:
	-rm -rf ../argtable2-13
	-rm -rf ../rng1-0
	-rm -f ../bin/etrans ../bin/etrans.MPI
	-rm -rf .o.mpi .o.serial

