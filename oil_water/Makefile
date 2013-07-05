default: IMPES

clean:
	/bin/rm -f *.o *.m *.m.dat *.vtk  mixed-poisson linear-heat heat-first heat-sec singlephase-first singlephase-sec nonlinear-heat SS IMPES implicit
mixed-poisson: mixed-poisson.o
linear-heat: linear-heat.o
nonlinear-heat: nonlinear-heat.o
heat-first: heat-first.o
heat-sec: heat-sec.o
singlephase-first: singlephase-first.o
singlephase-sec: singlephase-sec.o
#wophase.o: functions.c wophase.c


include /home/liuming/phg-matmult/Makefile.inc


#wophase: wophase.o functions.o
#	${LINKER} ${LDFLAGS} -o $@ $^ ${LIBS}
paremater.o: parameter.c
utils.o: utils.c
well.o: well.c 
uzawa.o:uzawa.c


##SS.o : SS.c parameter.c well.c utils.c uzawa.c
##SS: SS.o well.o parameter.o utils.o uzawa.o

IMPES.o : IMPES.c parameter.c well.c utils.c uzawa.c 
IMPES: IMPES.o well.o parameter.o utils.o uzawa.o 

#include ${PHG_MAKEFILE_INC}

