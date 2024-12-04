F90   = gfortran
F90OPTFLAGS= -g -O3  -ffree-form -ffree-line-length-none 
F90FLAGS =$(F90OPTFLAGS)

.SUFFIXES:
.SUFFIXES: .o .f .f90 .c 
#
.f90.o:; $(F90) -c $(F90FLAGS)  $(F90OPTFLAG) $<
.f.o:; $(F90) -c -loglist $(F90FLAGS)  $(F90OPTFLAG) $<
#
OBJS = \
Exact_Solver.o Multigrid_2D_Dirichlet_BC.o  Main.o \

TARGET = poisson2D_mg

all: $(TARGET)

$(TARGET): $(OBJS) 
	$(F90) $(F90FLAGS) -o $(TARGET) \
	$(OBJS) \
	$(F90FLAGS)

clean:
	rm -f *.o $(TARGET) *.mod *~ PI* *.log *.lst *.txt *.plt
