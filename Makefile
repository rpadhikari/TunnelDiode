FC=gfortran -c
LD=gfortran
SRC=params.f90 alloc.f90 dealloc.f90 main.f90
OBJ=params.o alloc.o dealloc.o main.o
diode:
        $(FC) $(SRC)
        $(LD) $(OBJ) -o diode.x -llapack -lblas
        rm -f *.mod *.o
clean:
        rm -f *.o *.mod *.x
