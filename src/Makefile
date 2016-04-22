CC=gcc

# icc
#CFLAGS= -O3

# gcc
# in GCC 4.9.2 -ftree-loop-vectorize leads to segmentation fault in sol_update (sol.c)
# -freee-loop-vectorize is enabled by -03, so use -O2
# this has no matter on caclulation speed
CFLAGS= -O2 -funroll-all-loops -ffast-math -fno-trapping-math

# IMPB openmpi
MPICFLAG = -I/usr/lib/openmpi/include -DMPI
MPILFLAG = -L/usr/lib/openmpi/lib -lmpi -lopen-rte -lopen-pal -ldl \
    -Wl,--export-dynamic -lnsl -lutil -ldl

# IMPB mpich2
#MPICFLAG = -I/usr/include/mpich2 -DMPI
#MPILFLAG = -L/usr/lib -lmpich -lopa -lmpl -lrt -lcr -lpthread

# IMPB mvapich2
#MPICFLAG = -I/opt/mvapich2/include -DMPI 
#MPILFLAG = -L/opt/mvapich2/lib -Wl,-rpath,/opt/mvapich2/lib \
#     -lmpich -lopa -lmpl -lpthread -lhwloc -lrdmacm -libverbs -libumad -ldl -lrt


# MVS100k
#MPICFLAG = -I/opt/mvapich/include -DMPI
#MPILFLAG = -Wl,-rpath-link -Wl,/opt/mvapich/lib/shared \
#    -L/opt/mvapich/lib/shared -L/opt/mvapich/lib \
#    -lmpich -L/usr/lib64 -Wl,-rpath=/usr/lib64 -libverbs -libumad

SRC = $(wildcard *.c)
OBJMPI = $(addprefix .o.mpi/,$(SRC:.c=.o))
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:	serial mpi

serial:	argtable librng $(OBJSRL)
	$(CC) -o ../bin/etrans $(OBJSRL) -lm \
	    ../argtable2-13/lib/libargtable2.a ../rng1-0/librng.a

mpi:	argtable librng $(OBJMPI)
	$(CC) -o ../bin/etrans.MPI $(OBJMPI) $(MPILFLAG) -lm \
            ../argtable2-13/lib/libargtable2.a ../rng1-0/librng.a

.o.mpi/%.o:	%.c
	$(CC) -I/usr/include/x86_64-linux-gnu $(CFLAGS) $(MPICFLAG) \
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
	-rm -r ../argtable2-13
	-rm -r ../rng1-0
	-rm -r $(OBJMPI) $(OBJSRL) ../bin/etrans ../bin/etrans.MPI
	-rmdir .o.mpi .o.serial

