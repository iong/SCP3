#!/usr/bin/make

CC:=gcc
CXX:=g++
FC:=gfortran

GFORTRAN_LIBDIR := $(subst --libdir=,,$(filter --libdir=%,$(shell $(FC) -v 2>&1)))

CPPFLAGS:=-I/opt/local/include# -DNDEBUG
FFLAGS:=-O3 -ftree-vectorize -march=native
#FFLAGS:=-O2 -g
CFLAGS:=-O3 -ftree-vectorize -march=native -fopenmp
#CFLAGS:=-O2 -g -fopenmp
CXXFLAGS:=$(CFLAGS)

LDFLAGS:=-L/opt/local/lib #-L$(GFORTRAN_LIBDIR)
LIBS :=-larmadillo -lhdf5 -framework Accelerate -lgfortran

FOBJS := water.o sobol.o sobol_stdnormal.o

all: SCP_scaled_nm

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

sobol_stdnormal.o : sobol.o

fobjs: $(FOBJS)

SCP_scaled_nm: SCP_scaled_nm.o $(FOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

clean:
	$(RM) *.o *.mod

.PHONY: clean
