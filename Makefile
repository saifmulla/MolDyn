# Makefile for the molecular dynamics code

#
# C compiler and options
# -mp
CC=     gcc-mp-4.8 -O3
LIB=	-lm -arch x86_64 -fopenmp -march=native

#
# Object files
#
OBJ=    main.o \
	dfill.o \
	domove.o \
        dscal.o \
	fcc.o \
	forces.o \
	mkekin.o \
	mxwell.o \
	prnout.o \
	velavg.o

#
# Compile
#
md:	$(OBJ)
	$(CC) -o $@ $(OBJ) $(LIB)

.c.o:
	$(CC) -c $<

#
# Clean out object files and the executable.
#
clean:
	rm *.o md
