include ../../Makefile.in
include ../../A2D_SHELLS.mk

default: thermBuckling.o therm-buckling.x

therm-buckling.x:
	${CXX} -o therm-buckling.x thermBuckling.o ${TACS_LD_FLAGS}

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
