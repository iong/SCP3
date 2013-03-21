CC:=icc
CXX:=icpc
FC:=ifort

FFLAGS += -heap-arrays

TARGET_FLAGS := 
ifdef TARGET
    TARGET_FLAGS := -m$(TARGET)
endif


COPT ?= -fast
FOPT ?= -fast

CDBG ?= -O0 -g 
FDBG ?= -fpe0 -traceback -check all -ftrapuv -warn unused

FORTRAN_LIBS := -lifport -lifcore

OPENMP_FLAGS ?= -openmp
