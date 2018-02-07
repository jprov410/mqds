export MPI_COMPILER=GNU
### suffix rule
.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
		$(F90) $(COMPFLAGS) $*.f90

### compiler
F90 = gfortran
COMMONFLAGS = 
COMPFLAGS   =  -c $(COMMONFLAGS) -llapack
LINKFLAGS   =  $(COMMONFLAGS)
#LDFLAGS     = mpilib

a=kinds.o
b=parameters.o
#d=input_variables.o
c=input_output.o
e=random_numbers.o
f=unit_conversions.o
g=harmonic_bath.o
h=mapping_variables.o
i=hamiltonians.o
j=pldm_master.o
k=calculate_pldm_redmat.o
l=main.o



### objects
OBJ = $(a) $(b) $(c) $(e) $(f) $(g) $(h) $(i) $(j) $(k) $(l) # $(d)
### compile and link

test/mqds: $(OBJ)
	  $(F90) -o $@ $(LINKFLAGS) $(OBJ)

clean:
	rm -f *.o *.mod