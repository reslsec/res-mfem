default: single_phase

clean:
	/bin/rm -f *.o *.m *.m.dat *.vtk single_phase

#include /home/liuming/phg-matmult/Makefile.inc
include	/opt/local-MVAPICH2/phg-0.8.0/share/phg/Makefile.inc

well.o: well.c 
uzawa.o: uzawa.c 
single_phase.o: single_phase.c well.c uzawa.c
single_phase: single_phase.o well.o uzawa.o

new.o: new.c well.c
new: new.o well.o
