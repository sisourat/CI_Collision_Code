#
# Makefile for DOM examples
#
default: all
 all:   getsta # coll# testOv
# all:   getsta_tri  colltri # testOv
#
#---------------------------
MK=$(FLIB_ROOT)/fortran.mk
include $(MK)
#---------------------------
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)
#FFLAGS=$(FFLAGS_CHECK)
#
#LIBS=$(LIB_PREFIX)$(LIB_STD)     -lflib  -I/opt/intel/mklinclude  -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIBS=$(LIB_PREFIX)$(LIB_STD)     -lflib -mkl 
#LIBS=$(LIB_PREFIX)$(LIB_STD) -lflib -L/usr/lib64/atlas/ -llapack -lblas -Bstatic_pgi -Bstatic 
##LIBS=$(LIB_PREFIX)$(LIB_STD) -lflib -L/opt/intel/mkl/lib/intel64 -I/opt/intel/mkl/lib/intel64  -llapack

#LDFLAGS = -p -g -pg
#LDFLAGS = -O3 -parallel
LDFLAGS =  -qopenmp -O3 -pg#-xHost  -no-prec-div 
#FC = fortran-intel
#testov: general.o overlap.o testOv.o       
#	$(FC) $(LDFLAGS) -o testOv general.o overlap.o testOv.o $(LIBS)

#testpot: general.o f11.o nuclear.o special_functions.o  testNuclear.o
#	$(FC) $(LDFLAGS) -o testPot general.o f11.o nuclear.o special_functions.o  testNuclear.o $(LIBS)
#	$(FC) $(LDFLAGS) -o testPot general.o f11.o nuclear.o special_functions.o  testNuclear.o $(LIBS)



getsta: general.o tools.o newtypes.o matrices.o centerlib.o fdn.o setup.o inputlib.o cgto.o overlap.o special_functions.o f11.o nuclear.o repulsion.o  computemat.o matrep.o  computemat_tt.o   diagolib.o  Get2eStaMoints.o
	$(FC) $(LDFLAGS) -o Get2eStaMoints general.o tools.o newtypes.o matrices.o centerlib.o fdn.o setup.o inputlib.o cgto.o overlap.o special_functions.o f11.o nuclear.o  repulsion.o  computemat.o matrep.o  computemat_tt.o  diagolib.o Get2eStaMoints.o $(LIBS)


clean : 
	rm -f *.o
	rm -f *.mod




