# The compiler
FC = ifort #gfortran
FCFLAGS = -fopenmp -O3

OBJECTS = general.o abmlib.o linearinterp.o  splineinterp.o  matrices.o propdyn.o Prop.o 

PROGRAMS = Prop


.SUFFIXES: .f90 .for .f .o .a

.f.o :
	$(FC) -c  $(*).f $(FCFLAGS)

.f90.o :
	$(FC) -c  $(*).f90 $(FCFLAGS)


# "make" builds all

all: $(PROGRAMS)

Prop:   $(OBJECTS) 
	$(FC) $(LDFLAGS) -o Prop $(OBJECTS) $(FCFLAGS)


# Utility targets

tar:
	tar cvf Prop.tar *.f* Makefile

clean:
	rm -f *.o *.mod *.MOD Adiab Auger

