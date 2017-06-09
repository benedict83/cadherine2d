#HOME = .

FC = gfortran

FCFLAGS = -I/usr/include/fgsl -lfgsl
#LDFLAGS = -L/usr/lib -lgsl -lgslcblas bug apres ubuntu 15.10 update
LDFLAGS = -L/usr/lib

## ======================================================================
## OPTIONS
## ======================================================================
## -O3 the most optimization
## -O2
## -O1
## -O0 no optimization
## -Og optimize debugging experience
## utiles pour le developpement :
## -f        Fortran Language Options
## -g 	     actually makes debugging possible by requesting the compiler to produce debugging information intended to be used by interactive debugger (GDB)
## -pedantic Error and Warning Options
## -W        Error and Warning Options
OPT = -Og -fPIC -fbacktrace -fcheck=all -fdiagnostics-color=auto -fmax-errors=1 \
      -g -pedantic -Wall -Wextra -Wimplicit-interface

## ======================================================================
## OBJECTS
## ======================================================================
OBJECTS = rock4.o

## ======================================================================
## And now the general rules, these should not require modification
## ======================================================================
exemple: 
	$(FC) -o $@ $@.f90 $(FCFLAGS) $(LDFLAGS) $(OPT)
## General rule for building prog from prog.o
## $^ (GNU extension) is used in order to list
## additional object files on which the executable depends
## I added $(OBJECTS) for cadhe2d
%: $(OBJECTS) %.o
	$(FC) $(FCFLAGS) $(LDFLAGS) $(OPT) $^ -o ../exe/$@
## General rule for building prog.o from prog.f90
## $< is used in order to list only the first prerequisite (the source file)
## and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) $(LDFLAGS) $(OPT) -c $<

## ======================================================================
## CLEAN
## ======================================================================
clean:
	rm -f *.mod *.o
