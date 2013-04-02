#Generic makefile for fortran90

LIBDIR := ../libs
MODDIR := ../include

FC := ifort -r8 
# serial
LIBS := -static -L$(LIBDIR) -lpes3bifc -ldms2bifc -lpes2bifc 
FFLAGS := -O -I$(MODDIR)/mod_3bifc # -heap-arrays 396900 
# parallel  
#LIBS := -static -L$(LIBDIR) -lpes3bifc_omp -ldms2bifc_omp -lpes2bifc_omp
#FFLAGS := -O -openmp -parallel -u -I$(MODDIR)/mod_3bifc_omp 

vpath %.f90 src

OBJ := \
math.o\
smear.o\
ttm3f_mod2.o\
mnasa_mod.o\
mnasa.o\
ttm3f_mb.o\
pot_monomer_mod.o\
pot_monomer.o\
pes3b.o\
pes_shell.o\
getpot.o

# Targets:
.PHONY : clean 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< 
%.o : %.f
	$(FC) -c $(FFLAGS) $< 

all : getpot

getpot :  $(OBJ) 
	$(FC) -o $@  $^ $(FFLAGS) $(LIBS)

pes_shell.o:   pot_monomer_mod.o nbody.f90 
getpot.o:      pes_shell.o
pot_monomer.o: pot_monomer_mod.o
ttm3f_mb.o:    ttm3f_mod2.o
mnasa.o:       ttm3f_mod2.o mnasa_mod.o

clean :
	rm *.o *.mod 
