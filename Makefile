SRCDIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))
OS := $(shell uname -s)

COMPILER ?= GNU

ifdef PROFILE
	include $(SRCDIR)/config/$(PROFILE).mk
endif

include $(SRCDIR)/config/$(COMPILER).mk
include $(SRCDIR)/config/blas.mk

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
$(foreach x,CFLAGS CXXFLAGS FFLAGS,$(eval $(x) += $(OPENMP_FLAGS) ) )


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
	LIBS += -lpes3bifc -lpes2bifc -lifcore -limf -lsvml -Bstatic -lnetcdf -Bdynamic
	CPPFLAGS += -DHAVE_BOWMAN
endif

all: SCP3

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

sobol_stdnormal.o : sobol.o

SCP3: SCP3.o ScaledMatrixElements.o \
	DiskIO.o sobol.o Constants.o beasley_springer_moro.o \
	$(X2O_OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

TestMatrixElements: TestMatrixElements.o TestScaledMatrixElements.o \
		ScaledMatrixElements.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

TestSimplexIterator: TestSimplexIterator.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

SelectEW: SelectEW.o DiskIO.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

Upot: Upot.o Constants.o DiskIO.o  $(X2O_OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: clean
