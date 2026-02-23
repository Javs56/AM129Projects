#I need to change the object files from hw2 to project


FC = gfortran
FFLAGS = -Wall -Wextra -ffree-line-length-none -fPIC -g -fcheck=all -fbacktrace
LFLAGS = 
 #child.o -> child.o -> parent.o
OBJ = utility.o read_initFile.o problemSetup.o leapfrog.o fput.o  # add changes to this line as needed
											
fput.ex: $(OBJ)  #fput is the name of my program, can use 'main.ex' or whatever
	$(FC) -o $@ $^ $(LFLAGS) $(FFLAGS)

%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

.PHONY: clean

clean:
	rm *.ex *.o *.mod
          