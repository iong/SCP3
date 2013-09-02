SRCDIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))
OS := $(shell uname -s)

COMPILER ?= GNU

ifdef PROFILE
	include $(SRCDIR)/config/$(PROFILE).mk
endif

include $(SRCDIR)/config/$(COMPILER).mk
include $(SRCDIR)/config/blas.mk

ifdef MPI
CC:=mpicc
CXX:=mpic++
FC:=mpif90
endif

vpath %.cpp $(SRCDIR)
vpath %.h $(SRCDIR)
vpath %.f90 $(SRCDIR)

CPPFLAGS += -I$(SRCDIR)

ifdef DEBUG
$(foreach x,CFLAGS CXXFLAGS,$(eval $(x) += $(CDBG) ) )
	FFLAGS    += $(FDBG)
else
$(foreach x,CFLAGS CXXFLAGS,$(eval $(x) += $(COPT) ) )
	FFLAGS    += $(FOPT)
endif
$(foreach x,CFLAGS CXXFLAGS FFLAGS,$(eval $(x) += $(TARGET_FLAGS) ) )


LIBS :=$(BLAS_LIBRARIES) -lhdf5_cpp -lhdf5 -lz
ifeq ($(OS),Darwin)
LIBS +=-larmadillo
endif

X2O_OBJ =  ps.o qtip4pf.o \
              ttm3f-bits.o ttm3f-bits-smear.o ttm3f.o ttm4-es.o \
              ttm4-smear.o gammq.o
ifdef WHBB
	X2O_OBJ += bowman-bits.o bowman.o bowman-fortran.o ttm4-hbb2-x3b.o x3b.o
	LDFLAGS += -L$(SRCDIR)/bowman
	LIBS += -lpes3bifc -lpes2bifc -lifcore -limf -lsvml -Bstatic -lnetcdf -Bdynamic -ldl
	CPPFLAGS += -DHAVE_BOWMAN
endif

all: SCP3 eigensolver

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $(cxxflags-$*) -o $@ $<

%.o: %.c
	$(CXX) -c $(CPPFLAGS) $(CCFLAGS) $(ccflags-$*) -o $@ $<

sobol_stdnormal.o : sobol.o

cxxflags-SCP3 := $(OPENMP_FLAGS)
cxxflags-ScaledMatrixElements := $(OPENMP_FLAGS)

SCP3: SCP3.o Hessian.o ScaledMatrixElements.o \
	HarmonicOscillatorBasis.o \
	DiskIO.o sobol.o Constants.o beasley_springer_moro.o \
	$(X2O_OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $(OPENMP_FLAGS) $^ $(LDFLAGS) $(LIBS)

SCP1: SCP1_main.o SCP1.o Hessian.o \
	DiskIO.o sobol.o Constants.o beasley_springer_moro.o \
	$(X2O_OBJ)
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBS)


TestMatrixElements: TestMatrixElements.o TestScaledMatrixElements.o \
		ScaledMatrixElements.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

TestSimplexIterator: TestSimplexIterator.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

eigensolver: eigensolver.o DiskIO.o Constants.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

Upot: Upot.o Constants.o DiskIO.o  $(X2O_OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: clean
