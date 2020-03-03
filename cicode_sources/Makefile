#FC = gfortran
FC = ifort
#DEBUGFLAGS = -g -ftrace=full -fbounds-check
#FFLAGS = -fzero -w
#DEBUGFLAGS =  -C -traceback
FFLAGS = -O3 -qopenmp  -mkl  -static-intel
#LIBFLAGS =  -framework vecLib
LIBFLAGS = -O3 -qopenmp  -mkl -static-intel 

#SRCS = CreateCSF.f90
 SRCS = myCI.f90
# SRCS = matColl.f90
#SRCS = Fano.f90
OBJS = $(SRCS:.f90=.o)
EXEC = $(SRCS:.f90=)

all : $(EXEC)

%.o : %.f90
	$(FC) $(DEBUGFLAGS) $(PARFLAGS) $(FFLAGS) -c  $<

%.o : %.f
	$(FC) $(DEBUGFLAGS) $(PARFLAGS) $(FFLAGS)  -c  $<

%: %.o
	$(FC) $(DEBUGFLAGS) $(PARFLAGS) $(LIBFLAGS) $(LDFLAGS) -o $@ $^   

CreateCSF : intindex.o\
       phis.o\
       determinant.o\
       math.o\
       csf.o\
       cidet.o

myCI : intindex.o\
       phis.o\
       determinant.o\
       math.o\
       pythag.o\
       recipies.o\
       lanczos.o\
       diago.o\
       diagolib.o

matColl : intindex.o\
       phis.o\
       determinant.o

clean : 
	rm -f *.o
	rm -f *.mod

tar :
	tar cvf myCI.tar *.f* Makefile 

