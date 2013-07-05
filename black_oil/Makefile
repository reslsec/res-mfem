default: black_oil

clean:
	/bin/rm -f *.o  *.m.dat *.vtk  *.m singlepahse  oil_water_phase black_oil SS 
mixed-poisson: mixed-poisson.o
linear-heat: linear-heat.o
nonlinear-heat: nonlinear-heat.o
heat-first: heat-first.o
heat-sec: heat-sec.o
singlephase-first: singlephase-first.o
singlephase-sec: singlephase-sec.o

#include /home/liuming/phg-0.8.0/Makefile.inc

include /opt/local-MVAPICH2/phg-0.8.3/share/phg/Makefile.inc
paremater.o: parameter.c
quadfunc.o: quadfunc.c
well.o: well.c 
uzawa.o:uzawa.c

black_oil.o : black_oil.c parameter.c well.c quadfunc.c write_matlab.c
black_oil: black_oil.o parameter.o well.o quadfunc.o write_matlab.o

SS.o : SS.c parameter.c well.c quadfunc.c uzawa.c
SS: SS.o parameter.o well.o quadfunc.o uzawa.o

#IMPES.o : IMPES.c parameter.c well.c quadfunc.c
#IMPES: IMPES.o parameter.o well.o quadfunc.o

#include ${PHG_MAKEFILE_INC}
