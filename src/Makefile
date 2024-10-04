# ============================================
#
# Make file for TACS_DIR/src/
#
# ============================================

include ../Makefile.in
include ../A2D_SHELLS.mk

CXX_OBJS = TACSAssembler.o \
	TACSCreator.o \
	TACSAssembler_thread.o \
	TACSBuckling.o

DIR=${TACS_DIR}/src

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o
