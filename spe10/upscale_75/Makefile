default:oil_water up_75 
clean:
	-/bin/rm -f *.o *.tmp oil_water *.vtk  up_75  simplest *.out *.dat 
distclean: clean
	-/bin/rm -f oil_water *.o *~
lib:
	@(cd ../src; $(MAKE))
OBJECTS = rt.o well.o parameter.o quadfunc.o struct.o uzawa.o oil_water.o export.o forlm.o

oil_water: $(OBJECTS) 
	$(LINKER) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

up_75_OBJECTS = rt.o well.o parameter.o quadfunc.o struct.o uzawa.o up_75.o export.o forlm.o
up_75: $(up_75_OBJECTS) 
	$(LINKER) $(LDFLAGS) -o $@ $(up_75_OBJECTS) $(LIBS)
#include /home/liuming/phg-0.8.4/Makefile.inc
#include ${PHG_MAKEFILE_INC}

#include /home/liuming/phg-0.8.4/Makefile.inc
.PHONY: default all clean distclean lib
include /opt/local-MVAPICH2/phg-0.8.3/share/phg/Makefile.inc
