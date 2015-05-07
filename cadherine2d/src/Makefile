HOME = .
FC = gfortran
MACPORTS = /usr
FCFLAGS = -I$(MACPORTS)/include/fgsl -lfgsl
LDFLAGS = -L$(MACPORTS)/lib -lgsl -lgslcblas
#LDFLAGS = -L$(MACPORTS)/lib
#LOADLIBES = -lgsl -lgslcblas

## ======================================================================
## OPTIONS
## ======================================================================
# -O2 est le niveau par defaut
# utile pour le developpement :
#      -f        Fortran Language Options
#      -g 
#      -pedantic Error and Warning Options
#      -W        Error and Warning Options
# -O0 supprime toute optimisation
#OPT = -O2
OPT = -fPIC -fmax-errors=1 -fcheck=all -fbacktrace -g -pedantic -Wall -Wconversion -Wextra -Wimplicit-interface 
#OPT = -O0

PROGRAMS = exemple

exemple:
	$(FC) $(OPT) exemple.f90 $(LDFLAGS) $(FCFLAGS)
#	$(FC) $(OPT) exemple.f90 $(LDFLAGS) $(LDLIBS) $(FCFLAGS)

cadhe2d: $(INCLUDE) cadhe2d.f90 rock4.o
	$(FC) $(OPT) cadhe2d.f90 rock4.o -o cadhe2d
#	$(FC) $(OPT) cadhe2d.f90 rock4.o -o ../exe/cadhe2d # peut-on faire autrement ?

#all: $(PROGRAMS)                                          # necessaire ?

## ======================================================================
## And now the general rules, these should not require modification
## ======================================================================

## General rule for building prog from prog.o; $^ (GNU extension) is
## used in order to list additional object files on which the
## executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

## General rules for building prog.o from prog.f90 or prog.F90; $< is
## used in order to list only the first prerequisite (the source file)
## and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
#	$(FC) $(OPT) -c $< 

## ======================================================================
## CLEAN
## ======================================================================
clean:
	rm -f *.mod *.o
