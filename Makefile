SRCDIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))
OS := $(shell uname -s)

BLAS ?= acml_mp
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
$(foreach x,CFLAGS CXXFLAGS LDFLAGS,$(eval $(x) += $(CDBG) ) )
	FFLAGS    += $(FDBG)
else
$(foreach x,CFLAGS CXXFLAGS LDFLAGS,$(eval $(x) += $(COPT) ) )
	FFLAGS    += $(FOPT)
endif
$(foreach x,CFLAGS CXXFLAGS FFLAGS LDFLAGS,$(eval $(x) += $(TARGET_FLAGS) ) )
$(foreach x,CFLAGS CXXFLAGS FFLAGS LDFLAGS,$(eval $(x) += $(OPENMP_FLAGS) ) )


LIBS :=$(BLAS_LIBRARIES)
ifeq ($(OS),Darwin)
LIBS :=-larmadillo -lhdf5 $(BLAS_LIBRARIES)
endif

FOBJS := water.o sobol.o sobol_stdnormal.o

all: SCP_scaled_nm

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

sobol_stdnormal.o : sobol.o

fobjs: $(FOBJS)

SCP_scaled_nm: SCP_scaled_nm.o ScaledMatrixElements.o $(FOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(FORTRAN_LIBS)

TestMatrixElements: TestMatrixElements.o TestScaledMatrixElements.o \
		ScaledMatrixElements.o $(FOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(FORTRAN_LIBS)

TestSimplexIterator: TestSimplexIterator.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)
clean:
	$(RM) *.o *.mod

.PHONY: clean
