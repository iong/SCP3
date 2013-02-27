#!/usr/bin/make

CC:=/usr/bin/gcc
CXX:=/usr/bin/g++
FC:=gfortran

CPPFLAGS:=-I/opt/local/include# -DNDEBUG
FFLAGS:=-O3 -ftree-vectorize
FFLAGS:=-O2 -g
CFLAGS:=-O3 -ftree-vectorize -fopenmp
CFLAGS:=-O2 -g -fopenmp
CXXFLAGS:=$(CFLAGS)

LDFLAGS:=-L/opt/local/lib -larmadillo -lhdf5 -lboost_program_options-mt -framework Accelerate -Wl,-rpath -Wl,$(PWD)

all: SCP_scaled_nm# fortran.dylib

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $<

sobol_stdnormal.o : sobol.o

fortran.dylib: water.o sobol.o sobol_stdnormal.o
	$(FC) $^ -shared -o $@ -install_name @rpath/$@ -static-libgcc

SCP_scaled_nm: SCP_scaled_nm.o fortran.dylib
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

clean:
	$(RM) *.o *.mod

.PHONY: clean
