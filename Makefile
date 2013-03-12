#!/usr/bin/make

OS := $(shell uname -s)

CC:=gcc
CXX:=g++
FC:=gfortran

GFORTRAN_LIBDIR := $(subst --libdir=,,$(filter --libdir=%,$(shell $(FC) -v 2>&1)))

#CPPFLAGS += -I/opt/local/include# -DNDEBUG
OPTFLAGS := -O3 -ftree-vectorize -march=native
FFLAGS += $(OPTFLAGS)
#FFLAGS:=-O2 -g
CFLAGS += $(OPTFLAGS) -fopenmp
#CFLAGS:=-O2 -g -fopenmp
CXXFLAGS := $(CFLAGS)

#LDFLAGS += -L/opt/local/lib #-L$(GFORTRAN_LIBDIR)

LIBS :=-lacml_mp -lgfortran
ifeq ($(OS),Darwin)
LIBS :=-larmadillo -lhdf5 -framework Accelerate -lgfortran
endif

FOBJS := water.o sobol.o sobol_stdnormal.o

all: SCP_scaled_nm

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

sobol_stdnormal.o : sobol.o

fobjs: $(FOBJS)

SCP_scaled_nm: SCP_scaled_nm.o ScaledMatrixElements.o $(FOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: clean
