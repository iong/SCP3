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

WHBB_OBJ := math.o smear.o ttm3f_mod2.o mnasa_mod.o mnasa.o ttm3f_mb.o pot_monomer_mod.o pot_monomer.o pes3b.o pes_shell.o

X2O_OBJ =  ps.o qtip4pf.o \
              ttm3f-bits.o ttm3f-bits-smear.o ttm3f.o ttm4-es.o \
              ttm4-smear.o gammq.o
ifdef WHBB
	OMP_SUFFIX:=$(if $(OPENMP_FLAGS),_omp,)

	FFLAGS += -I$(SRCDIR)/water_WHBB/include/mod_3bifc$(OMP_SUFFIX)
	LDFLAGS += -L$(SRCDIR)/water_WHBB/libs
	LIBS += -lpes3bifc$(OMP_SUFFIX) -ldms2bifc$(OMP_SUFFIX) -lpes2bifc$(OMP_SUFFIX) 
	FOBJS += $(WHBB_OBJ)
endif

all: SCP3

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

sobol_stdnormal.o : sobol.o

fobjs: $(FOBJS)

SCP3: SCP3.o ScaledMatrixElements.o DiskIO.o sobol.o Constants.o \
	beasley_springer_moro.o $(X2O_OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

TestMatrixElements: TestMatrixElements.o TestScaledMatrixElements.o \
		ScaledMatrixElements.o $(FOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(FORTRAN_LIBS)

TestSimplexIterator: TestSimplexIterator.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

SelectEW: SelectEW.o DiskIO.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

Upot: $(FOBJS) Upot.o DiskIO.o
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(FORTRAN_LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: clean
