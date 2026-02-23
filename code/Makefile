FC = gfortran
FFLAGS = -g -Wall -Wextra -Wimplicit-interface -fmax-errors=5

OBJ = utility.o problemsetup.o discreteFourierTransform.o boundaryValueProblem.o

solveBVP.ex: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

.PHONY: clean

clean:
	rm -f solveBVP.ex *.o *.mod
