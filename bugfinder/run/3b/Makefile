FC=gfortran
MPIFILE=nompi
EXEC=drive3b

MATRIX=matrixstupid

FFLAGS=-O2
F77FLAGS=$(FFLAGS) -fdefault-real-8 -fdefault-double-8
LDFLAGS=-O2
LIBS=

.SUFFIXES:
.SUFFIXES: .o .f .f90

.f90.o:
	$(FC) $(FFLAGS) -c $<

OBJECTS=\
   drive3b.o\
   pot.o\
   ran.o\
   stack.o\
   v6pot.o\
   $(MATRIX).o\
   $(MPIFILE).o\
   cheft.o\
   v3bpot.o\
   correlator.o

$(EXEC): $(OBJECTS)
	$(FC) $(LDFLAGS) -o ./$(EXEC) $(OBJECTS) $(LIBS)

clean:
	rm -f *\.mod *\.o *~

drive3b.o: stack.o ran.o $(MPIFILE).o v6pot.o v3bpot.o $(MATRIX).o\
   correlator.o drive3b.f90
	$(FC) $(FFLAGS) -c drive3b.f90

stack.o: stack.f90
	$(FC) $(FFLAGS) -c stack.f90

ran.o: ran.f90
	$(FC) $(FFLAGS) -c ran.f90

$(MPIFILE).o: $(MPIFILE).f90
	$(FC) $(FFLAGS) -c $(MPIFILE).f90

$(MATRIX).o: $(MATRIX).f90
	$(FC) $(FFLAGS) -c $(MATRIX).f90

v6pot.o: cheft.o v6pot.f90 $(MPIFILE).o
	$(FC) $(FFLAGS) -c v6pot.f90

pot.o: pot.f
	$(FC) $(F77FLAGS) -c pot.f

cheft.o: cheft.f90
	$(FC) $(FFLAGS) -c cheft.f90

correlator.o: $(MATRIX).o correlator.f90
	$(FC) $(FFLAGS) -c correlator.f90
