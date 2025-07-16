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

MPICFLAGS= -O2 -funroll-all-loops -ffast-math -fno-trapping-math -DMPI
MPILFLAGS= -O3

SRC = $(wildcard *.c)
OBJMPI = $(addprefix .o.mpi/,$(SRC:.c=.o))
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:	serial mpi

serial:	argtable librng $(OBJSRL) | bin
	$(CC) -o ../bin/etrans $(OBJSRL) $(LFLAGS) -lm \
	    ../argtable2-13/lib/libargtable2.a ../rng/librng.a

mpi:	argtable librng $(OBJMPI) | bin
	$(MPICC) -o ../bin/etrans.MPI $(OBJMPI) $(MPILFLAGS) -lm \
            ../argtable2-13/lib/libargtable2.a ../rng/librng.a

.o.mpi/%.o:	%.c
	$(MPICC) -I/usr/include/x86_64-linux-gnu $(MPICFLAGS) \
	    -I../argtable2-13/include -I../rng \
	    -c -o $@ $<

.o.serial/%.o:	%.c
	$(CC) -I/usr/include/x86_64-linux-gnu ${CFLAGS} \
	    -I../argtable2-13/include -I../rng \
	    -c -o $@ $<

$(OBJMPI):	| .o.mpi

$(OBJSRL):	| .o.serial

.o.mpi:
	mkdir .o.mpi

.o.serial:
	mkdir .o.serial

bin:
	mkdir -p ../bin

argtable:	| ../argtable2-13

../argtable2-13: ../argtable2-13.tar.gz
	( \
	    cd ..; \
	    tar -xzf argtable2-13.tar.gz; \
	    cd argtable2-13; \
	    CFLAGS='--std=c90' ./configure --prefix=`pwd`; \
	    make install; \
	    cd ../src; \
	)

../argtable2-13.tar.gz:
	( \
	    cd ..; \
	    wget http://prdownloads.sourceforge.net/argtable/argtable2-13.tar.gz; \
	    cd src; \
	)

librng:		| ../rng

../rng:
	( \
	    cd ..; \
	    git clone https://github.com/egorsobolev/rng.git; \
	    cd rng; \
            make; \
	    cd ../src; \
	)

.PHONY:	clean
clean:
	-rm -rf ../argtable2-13
	-rm -rf ../rng
	-rm -rf ../bin
	-rm -rf .o.mpi .o.serial

