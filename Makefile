# Nom du compilateur
FC= gfortran

# Directives de compilation
FCFLAGS= -g -Og -fcheck=all -Wall

DEBUG = -cpp -std=f2008 -pedantic -Wconversion -Wall \
	-Wcharacter-truncation -Wunderflow -g -fbounds-check -fbacktrace \
	-fimplicit-none -fdump-core -ffpe-trap=invalid,zero,underflow,denormal \
DEBUG_OPTIM = $(DEBUG) -Wextra -Warray-temporaries -ffree-line-length-0 \
	-fcheck=all -finit-real=nan

# Sources (l’ordre est important !)
SRCS= mod1.f90 ter.f90

# Nom de l’exécutable
PROGRAM=TER

# Ne pas éditer en dessous de cette ligne
OBJ=$(SRCS:.f90=.o)


all: $(PROGRAM)

$(PROGRAM): $(OBJ)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

clean:
	rm -f *.o *.mod *.vtk *.dat fort.*
