#gfortran exemple.f90 -L/usr/lib -I/usr/include/fgsl -lfgsl -lgsl -lgslcblas
#
HOME = .
FC = gfortran
MACPORTS = /usr
FCFLAGS = -I$(MACPORTS)/include/fgsl -lfgsl
LDFLAGS = -L$(MACPORTS)/lib -lgsl -lgslcblas
#LDFLAGS = -L$(MACPORTS)/lib
#LOADLIBES = -lgsl -lgslcblas

#----------OPTIONS-----------------------------
# -O2 est le niveau par defaut
# -O1 -g -fbounds-check est utile pour débugger
# -O0 supprime toute optimisation
#----------------------------------------------
#OPT = -O2
OPT = -O1 -g -fbounds-check -Wall -Wextra -Wconversion -pedantic
#OPT = -O0

PROGRAMS = exemple

#----------------------------------------------
# cadhe2d solves the system of reaction-diffusion
# equations for the cadherin model in two dimensions
#----------------------------------------------
exemple:
	$(FC) $(OPT) exemple.f90 $(LDFLAGS) $(FCFLAGS)
#	$(FC) $(OPT) exemple.f90 $(LDFLAGS) $(LDLIBS) $(FCFLAGS)

cadhe2d: $(INCLUDE) cadhe2d.f90 rock4.o
	$(FC) $(OPT) cadhe2d.f90 rock4.o -o ../exe/cadhe2d
#cadhe2d: $(INCLUDE) cadhe2d.f90 rock4.o
#	$(FC) -Wl,--allow-multiple-definition $(OPT)
#	cadhe2d.f90 rock4.o -o ../exe/cadhe2d
#	cp ../exe/cadhe2d ../data/run_input/cadhe2d

all: $(PROGRAMS)

#----------------------------------------------
# integrators
#----------------------------------------------
#exemple.o: fgsl.o integral.o
rock2.o : rock2.f
	$(FC) $(OPT) -c rock2.f
rock4.o : rock4.f90
	$(FC) $(OPT) -c rock4.f90

## ======================================================================
## And now the general rules, these should not require modification
## ======================================================================
#
## General rule for building prog from prog.o; $^ (GNU extension) is
## used in order to list additional object files on which the
## executable depends
#%: %.o
#	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)
#
## General rules for building prog.o from prog.f90 or prog.F90; $< is
## used in order to list only the first prerequisite (the source file)
## and not the additional prerequisites such as module or include files
#%.o: %.f90
#	$(FC) $(FCFLAGS) -c $<
#
#%.o: %.F90
#	$(FC) $(FCFLAGS) -c $<
#

#----------------------------------------------
# clean
#----------------------------------------------
clean:
	rm -f *.mod *.o cadhe2d
