
FC = gfortran
FLAGS = -O3 -ftree-vectorize -march=native 
OMP = -fopenmp
TARGET = fractal 

SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)
LIBS = #-lpthread

%.o: %.f90 
	$(FC) $(OMP) $(FLAGS) -c  $<

all: $(OBJS)
	$(FC) -o $(TARGET) $(OMP) $(OBJS) $(LIBS)


clean:
	rm *.o *.mod $(TARGET)

include make.dep
