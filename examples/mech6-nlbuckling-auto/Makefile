include ../../Makefile.in
include ../../A2D_SHELLS.mk

LOCAL_OBJS = mechNLBuckling.o
EXECS = mech-nlbuckling.x

default: ${LOCAL_OBJS} ${EXECS}

mech-nlbuckling.x:
	${CXX} -o mech-nlbuckling.x mechNLBuckling.o ${TACS_LD_FLAGS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

complex: TACS_DEF="-DTACS_USE_COMPLEX"
complex: default

complex_debug: TACS_DEF="-DTACS_USE_COMPLEX"
complex_debug: debug


new:
	$(MAKE) clean
	$(MAKE) debug

clean:
	rm -f *.o
	rm -f *.x

cleanf5:
	rm *.csv
	rm -f **/*.f5
	rm -f **/*.vtk
