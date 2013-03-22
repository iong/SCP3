CC:=gcc
FC:=gfortran

FFLAGS+=-ffree-line-length-0 -ffixed-line-length-0 -fimplicit-none

TARGET_FLAGS := -march=native
ifdef TARGET
    TARGET_FLAGS := -march=$(TARGET)
endif

COPT ?= -O3 -ffast-math
FOPT ?= $(COPT)

CDBG ?= -O0 -g -Wall
FDBG ?= -fbounds-check -ffpe-trap=invalid,zero,overflow,denormal \
	-finit-real=SNAN -finit-integer=-1

OPENMP_FLAGS ?= -fopenmp


FORTRAN_LIBS := -lgfortran

GFORTRAN_LIBDIR := $(subst --libdir=,,$(filter --libdir=%,$(shell $(FC) -v 2>&1)))
ifneq ($(GFORTRAN_LIBDIR),)
	FORTRAN_LIBS := -L$(GFORTRAN_LIBDIR) $(FORTRAN_LIBS)
endif
