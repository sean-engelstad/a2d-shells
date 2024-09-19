A2D_SHELLS_LIB = ${A2D_SHELLS_DIR}/lib/liba2dshells.${SO_EXT}
MY_SUBDIRS = src \
          src/io \
		  src/constitutive \
		  src/elements \
		  src/elements/shell \
		  src/bpmat
		  
# MY_OBJS := $(addsuffix /*.o, ${MY_SUBDIRS})
MY_OBJS := $(shell find $(A2D_SHELLS_DIR) -name '*.o')
MY_INCLUDE = -I${A2D_SHELLS_DIR}/src \
	-I${A2D_SHELLS_DIR}/src/io \
	-I${A2D_SHELLS_DIR}/src/constitutive \
	-I${A2D_SHELLS_DIR}/src/elements \
	-I${A2D_SHELLS_DIR}/src/elements/basis \
	-I${A2D_SHELLS_DIR}/src/elements/shell \
	-I${A2D_SHELLS_DIR}/src/bpmat

# TACS_EXTERN_LIBS = ${AMD_LIBS} ${METIS_LIB} ${LAPACK_LIBS} ${TECIO_LIBS}
EXTERN_LIBS = ${AMD_LIBS} ${METIS_LIB} ${LAPACK_LIBS}

ALL_FLAGS = ${CFLAGS} ${MY_INCLUDE} ${METIS_INCLUDE}

%.o: %.cpp
	${CC} ${ALL_FLAGS} -std=c++11 -c $< -o $*.o