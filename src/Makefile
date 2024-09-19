include ../Makefile.in
include ../A2D_SHELLS.mk

CXX_OBJS = TACSObject.o \
		TACSAssembler.o \
		TACSCreator.o \
		TacsUtilities.o \
		TACSAssembler_thread.o

default: ${CXX_OBJS}

clean:
	rm -f *.o
	rm -f **/*.o
	rm -f **/**/*.o