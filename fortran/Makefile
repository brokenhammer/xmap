# Use fortran compiler and f2py to generates modules and .so

FCMP := gfortran
EXT_SUFFIX := $(shell python3-config --extension-suffix)
F90OPT := --f90flags=-ffree-line-length-0
FTNFLG := -ffree-line-length-0

.PHONY: all

all: spdata_interface_f$(EXT_SUFFIX) spline_function_f$(EXT_SUFFIX)

spdata_interface_f$(EXT_SUFFIX): smooth.o spdata_interface.F90
	f2py $(F90OPT) -c -m spdata_interface_f spdata_interface.F90

spline_function_f$(EXT_SUFFIX): smooth.o function.F90
	f2py $(F90OPT) -c smooth.F90 function.F90 -m spline_function_f

smooth.o: smooth.F90
	$(FCMP) $(FTNFLG) -c smooth.F90

.PHONY: clean

clean:
	rm -f *.so *.o *.mod 
