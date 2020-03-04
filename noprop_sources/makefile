#
# Makefile for DOM examples
#
default: all
all: testov testpot getsta getpec coll sigma
#
#---------------------------
MK=$(FLIB_ROOT)/fortran.mk
include $(MK)
#---------------------------
#
# Uncomment the following line for debugging support
#
#FFLAGS=$(FFLAGS_DEBUG)  -ffree-line-length-0 -fopenmp -O3
FFLAGS=$(FFLAGS_DEBUG)  -fopenmp -O3
#FFLAGS=$(FFLAGS_CHECK)
#
LIBS=$(LIB_PREFIX)$(LIB_STD) -lflib -mkl

#LDFLAGS = -p -g -pg
LDFLAGS = -O3 -fopenmp #-parallel
#LDFLAGS = -O3 


testov: general.o overlap.o testOv.o
	$(FC) $(LDFLAGS) -o testOv general.o overlap.o testOv.o $(LIBS)

testpot: general.o f11.o nuclear.o special_functions.o  testNuclear.o
	$(FC) $(LDFLAGS) -o testPot general.o f11.o nuclear.o special_functions.o  testNuclear.o $(LIBS)

getsta: general.o tools.o newtypes.o matrices.o centerlib.o fdn.o setup.o inputlib.o cgto.o overlap.o special_functions.o f11.o nuclear.o computemat.o diagolib.o GetSta.o
	$(FC) $(LDFLAGS) -o GetSta general.o tools.o newtypes.o matrices.o centerlib.o fdn.o setup.o inputlib.o cgto.o overlap.o special_functions.o f11.o nuclear.o computemat.o diagolib.o GetSta.o $(LIBS)

getpec: general.o tools.o newtypes.o matrices.o centerlib.o fdn.o setup.o inputlib.o cgto.o overlap.o special_functions.o f11.o nuclear.o computemat.o diagolib.o GetPec.o
	$(FC) $(LDFLAGS) -o GetPec general.o tools.o newtypes.o matrices.o centerlib.o fdn.o setup.o inputlib.o cgto.o overlap.o special_functions.o f11.o nuclear.o computemat.o diagolib.o GetPec.o $(LIBS)

coll: general.o abmlib.o linearinterp.o  splineinterp.o tools.o setup.o newtypes.o matrices.o inputlib.o centerlib.o fdn.o colllib.o cgto.o colldyn.o  overlap.o special_functions.o f11.o nuclear.o  collint.o  Coll.o
	$(FC) $(LDFLAGS) -o Coll general.o abmlib.o linearinterp.o  splineinterp.o  tools.o setup.o newtypes.o matrices.o inputlib.o centerlib.o fdn.o colllib.o cgto.o colldyn.o overlap.o special_functions.o f11.o nuclear.o collint.o  Coll.o $(LIBS)

sigma: general.o newtypes.o gnufor.o sigma.o
	$(FC) $(LDFLAGS) -o Sigma general.o newtypes.o gnufor.o sigma.o $(LIBS)
#
clean: 
	rm -f getsta testov testpot  *.o *.$(MOD_EXT)
#

tar:
	tar cvf coll_jarvis170419.tar *.f* makefile xmlf90-1.2g.tgz ../basis










