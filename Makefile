
FC = gfortran  
FLAGS = -O3 -Wall -fcheck=all
#FLAGS = -g -CB -check all 
#FLAGS = -g -Warray-bounds -fbounds-check
 
TARGET1 = ns

SOURCES1 = common.f90 \
	   gmres.f90\
	   solvers.f90 \
	   main.f90

OBJS1 = $(SOURCES1:.f90=.o)
LIBS = -llapack

%.o: %.f90 
	$(FC) $(FLAGS) -c  $<

all: $(TARGET1) 

$(TARGET1): $(OBJS1)
	$(FC) -o $(TARGET1) $(OBJS1) $(LIBS)



clean:
	rm *.o *.mod $(TARGET1) 

cleandat:
	rm *.dat *.txt *.raw 

solvers.o : common.o
gmres.o : common.o
main.o : common.o solvers.o gmres.o
