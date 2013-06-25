CC:=icc
CXX:=icpc
FC:=ifort

FFLAGS += -heap-arrays

TARGET_FLAGS := 
ifdef TARGET
    TARGET_FLAGS := -m$(TARGET)
endif


COPT ?= -O3 -no-prec-div -xHost -DNDEBUG
FOPT ?= -O3 -no-prec-div -xHost

CDBG ?= -O0 -g 
FDBG ?= -fpe0 -traceback -check all -ftrapuv -warn unused

FORTRAN_LIBS := -lifport -lifcore

OPENMP_FLAGS ?= -openmp
