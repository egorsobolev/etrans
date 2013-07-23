CC=gcc

# icc
#CFLAGS= -O3

# gcc
CFLAGS= -O3 -funroll-all-loops -ffast-math -fno-trapping-math

# MATHLIB
# GSL
MATHCFLAGS=
MATHLFLAGS= -llapack -lgsl -lgslcblas 

# MKL
#MKLROOT=/opt/intel/mkl
#MATHCFLAGS= -I$(MKLROOT)/include -DUSE_MKL
#MATHLFLAGS= -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
#    $(MKLROOT)/lib/intel64/libmkl_sequential.a \
#    $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -limf

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

serial:	argtable $(OBJSRL)
	$(CC) -o ../bin/etrans $(OBJSRL) ../argtable2-13/lib/libargtable2.a \
            $(MATHLFLAGS) -lm

mpi:	argtable $(OBJMPI)
	$(CC) -o ../bin/etrans.MPI $(OBJMPI) ../argtable2-13/lib/libargtable2.a $(MPILFLAG) \
	    $(MATHLFLAGS) -lm

.o.mpi/%.o:	%.c
	$(CC) -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 \
	    -DHAVE_STDARG_H=1 -DUSE_STDARG=1 -DMALLOC_RET_VOID=1 -c \
	    -I/usr/include/x86_64-linux-gnu $(CFLAGS) $(MPICFLAG) \
	    -I../argtable2-13/include $(MATHCFLAGS) -o $@ $<

.o.serial/%.o:	%.c
	$(CC) -DUSE_STDARG -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 \
	    -DHAVE_STDARG_H=1 -DUSE_STDARG=1 -DMALLOC_RET_VOID=1 -c \
	    -I/usr/include/x86_64-linux-gnu ${CFLAGS} \
	    -I../argtable2-13/include $(MATHCFLAGS) -o $@ $<

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

.PHONY:	clean
clean:
	-rm -r ../argtable2-13
	-rm -r $(OBJMPI) $(OBJSRL) ../bin/etrans ../bin/etrans.MPI
	-rmdir .o.mpi .o.serial

