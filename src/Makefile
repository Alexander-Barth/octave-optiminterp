
sinclude Makeconf

SOURCES=optiminterp.cc optimal_interpolation.F90 optiminterp_wrapper.F90

OBJECTS=optiminterp.o optimal_interpolation.o optiminterp_wrapper.o

TARGETS=optiminterp.oct

ifndef LAPACK_LIBS
LAPACK_LIBS := $(shell $(MKOCTFILE) -p LAPACK_LIBS)
endif
ifndef BLAS_LIBS
BLAS_LIBS := $(shell $(MKOCTFILE) -p BLAS_LIBS)
endif
ifndef FLIBS
FLIBS := $(shell $(MKOCTFILE) -p FLIBS)
endif
LFLAGS := $(shell $(MKOCTFILE) -p LFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

# Rule for compiling Fortran 90 programs

%.o: %.F90
	FFLAGS="$(FFLAGS)" $(MKOCTFILE) -c $<

ifeq ($(MKOCTFILE_FORTRAN_90),yes)
all: $(TARGETS)
else
all:
	@echo "*** mkoctfile lacks Fortran 90 support"; false
endif


optiminterp.oct: $(OBJECTS)
	LFLAGS="$(LFLAGS)" \
               $(MKOCTFILE) -o $@ $(OBJECTS) $(LIBS) $(LAPACK_LIBS) $(FLIBS)

check:
	octave --silent --norc --eval test_optiminterp

clean:
	rm -f *.o *.oct *~ *.mod
