MY_SUBDIRS = src \
          src/io \
		  src/constitutive \
		  src/elements \
		  src/elements/shell
		  
# MY_OBJS := $(addsuffix /*.o, ${MY_SUBDIRS})
MY_OBJS := $(shell find $(A2D_SHELLS_DIR) -name '*.o')
MY_INCLUDE = -I${A2D_SHELLS_DIR}/src \
	-I${A2D_SHELLS_DIR}/src/io \
	-I${A2D_SHELLS_DIR}/src/constitutive \
	-I${A2D_SHELLS_DIR}/src/elements \
	-I${A2D_SHELLS_DIR}/src/elements/basis \
	-I${A2D_SHELLS_DIR}/src/elements/shell

%.o: %.cpp
	${CC} ${CFLAGS} ${MY_INCLUDE} -std=c++11 -c $< -o $*.o