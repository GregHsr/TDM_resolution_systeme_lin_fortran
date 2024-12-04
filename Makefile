# Makefile 
# Nom du compilateur
FC = gfortran

# Options de compilation: optimisation, debug etc...
OPT =  -g  -O3 
# nom de l'executable
EXE = poisson2d
# Options de l'edition de lien..
LINKOPT =  

# Defining the objects (OBJS) variables

OBJS =  \
       poisson2d.o \
       matgen.o \
       Pcg_Icc2.o \
	   ensight_laplacien.o

# Linking object files
exe :   $(OBJS)
	$(FC) $(LINKOPT) $(MODS) $(OBJS)  -o $(EXE) 

poisson2d.o : poisson2d.f90
	$(FC) -c $(OPT)  poisson2d.f90

matgen.o : matgen.f90
	$(FC) -c $(OPT) matgen.f90

Pcg_Icc2.o : Pcg_Icc2.f90
	$(FC) -c $(OPT) Pcg_Icc2.f90

ensight_laplacien.o : ensight_laplacien.f90
	$(FC) -c $(OPT) ensight_laplacien.f90
	

# Removing object files
clean :
	/bin/rm -f $(OBJS) $(EXE)  *.mod

config :
	if [ ! -d obj ] ; then mkdir obj ; fi
	if [ ! -d run ] ; then mkdir run ; fi
