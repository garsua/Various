# Makefile for GSBulk
# Victor M. Garcia-Suarez, September 2015.
FC=gfortran
FLAGS=-g -march=native -funroll-loops -O3 -ffast-math -fgcse-lm -fgcse-sm -ffast-math -ftree-vectorize -fexternal-blas
FLAGS_DEBUG=-g

%.o: %.f90
	$(FC) -c $(FLAGS) $(FLAGS_DEBUG) $<

OBJS=gsbulk.o tools.o ioabc_m.o

all:	$(OBJS)
	$(FC) -o gsbulk $(OBJS) $(LIBS)

clean:
	rm -f *.o *.mod gsbulk

gsbulk.o: ioabc_m.o
