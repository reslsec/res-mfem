default: oil_gas

clean:
	/bin/rm -f *.o  *.m.dat *.vtk  *.m singlepahse  oil_water oil_gas SS IMPES
#include /home/liuming/phg-matmult/Makefile.inc

include	/opt/local-MVAPICH2/phg-0.8.3/share/phg/Makefile.inc

paremater.o: parameter.c
quadfunc.o: quadfunc.c
well.o: well.c 
uzawa.o:uzawa.c

SS.o : SS.c parameter.c quadfunc.c well.c uzawa.c
SS: SS.o  parameter.o well.o quadfunc.o uzawa.o

IMPES.o : IMPES.c parameter.c quadfunc.c well.c uzawa.c
IMPES: IMPES.o  parameter.o well.o quadfunc.o uzawa.o

oil_gas.o : oil_gas.c parameter.c quadfunc.c well.c uzawa.c
oil_gas: oil_gas.o  parameter.o well.o quadfunc.o uzawa.o
#include ${PHG_MAKEFILE_INC}
