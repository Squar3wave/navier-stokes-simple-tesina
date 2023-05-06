
FC = gfortran  
#FLAGS = -O3
FLAGS = -O3 -Wall -fcheck=all
#FLAGS = -g -CB -check all 
#FLAGS = -g -Warray-bounds -fbounds-check
 
TARGET1 = ns

SOURCES1 = common.f90 \
	   gmres.f90\
	   solvers.f90 \
	   main.f90\
	   adj_map.f90 \
	   matdef.f90 \
	   sparsealg.f90

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

sparsealg.o: matdef.o
common.o:    matdef.o sparsealg.o
gmres.o :    matdef.o sparsealg.o common.o
solvers.o :  matdef.o sparsealg.o common.o gmres.o
main.o :     matdef.o sparsealg.o common.o gmres.o solvers.o
