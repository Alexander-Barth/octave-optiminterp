
sinclude Makeconf

SOURCES=optiminterp.cc optimal_interpolation.F90 optiminterp_wrapper.F90

OBJECTS=optiminterp.o optimal_interpolation.o optiminterp_wrapper.o

TARGETS=optiminterp.oct

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
	$(MKOCTFILE) -o $@ $(OBJECTS) $(LIBS)

check:
	octave --silent --norc --eval test_optiminterp

clean:
	rm -f *.o *.oct *~ *.mod
