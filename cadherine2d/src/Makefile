HOME = .
INCLUDE =
FC = gfortran
#----------OPTIONS-----------------------------
# -O2 est le niveau par defaut
# -O1 -g -fbounds-check est utile pour débugger
# -O0 supprime toute optimisation
#----------------------------------------------
#OPT = -O2
OPT = -O1 -g -fbounds-check
#OPT = -O0

#----------------------------------------------
# cadhe2d solves the system of reaction-diffusion
# equations for the cadherin model in two dimensions
#----------------------------------------------
cadhe2d: $(INCLUDE) cadhe2d.f90 rock4.o
	$(FC) $(OPT) cadhe2d.f90 rock4.o -o ../exe/cadhe2d
	cp ../exe/cadhe2d ../data/run_input/cadhe2d
#cadhe2d: $(INCLUDE) cadhe2d.f90 rock4.o
#	$(FC) -Wl,--allow-multiple-definition $(OPT)
#	cadhe2d.f90 rock4.o -o ../exe/cadhe2d
#	cp ../exe/cadhe2d ../data/run_input/cadhe2d

#----------------------------------------------
# integrators
#----------------------------------------------
rock2.o : rock2.f
	$(FC) $(OPT) -c rock2.f
rock4.o : rock4.f90
	$(FC) $(OPT) -c rock4.f90

#----------------------------------------------
# clean
#----------------------------------------------
clean:
	rm -f *.mod *.o cadhe2d
