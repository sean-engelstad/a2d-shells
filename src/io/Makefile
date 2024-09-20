# ============================================
#
# Make file for TACS_DIR/src/io/
#
# ============================================

include ../../Makefile.in
include ../../A2D_SHELLS.mk

CXX_OBJS = TACSFH5.o \
	TACSToFH5.o \
	TACSMeshLoader.o

DIR=${TACS_DIR}/src/io

CXX_OBJS := $(CXX_OBJS:%=$(DIR)/%)

default: ${CXX_OBJS}
	${AR} ${AR_FLAGS} ${TACS_LIB} ${CXX_OBJS}

debug: TACS_CC_FLAGS=${TACS_DEBUG_CC_FLAGS}
debug: default

clean:
	${RM} *.o

