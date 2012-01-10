
sinclude Makeconf

SOURCES=optiminterp.cc optimal_interpolation.F90 optiminterp_wrapper.F90

OBJECTS=optiminterp.o optimal_interpolation.o optiminterp_wrapper.o

TARGETS=optiminterp.oct

LAPACK_LIBS := $(shell mkoctfile -p LAPACK_LIBS)
FLIBS := $(shell mkoctfile -p FLIBS) 

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
	$(MKOCTFILE) -o $@ $(OBJECTS) $(LIBS) $(LAPACK_LIBS) $(FLIBS)

check:
	octave --silent --norc --eval test_optiminterp

clean:
	rm -f *.o *.oct *~ *.mod
