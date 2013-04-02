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

FOBJS := water.o sobol.o sobol_stdnormal.o

all: SCP3

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

sobol_stdnormal.o : sobol.o

fobjs: $(FOBJS)

SCP3: SCP3.o ScaledMatrixElements.o DiskIO.o $(FOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(FORTRAN_LIBS)

TestMatrixElements: TestMatrixElements.o TestScaledMatrixElements.o \
		ScaledMatrixElements.o $(FOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(FORTRAN_LIBS)

TestSimplexIterator: TestSimplexIterator.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)
clean:
	$(RM) *.o *.mod

.PHONY: clean
