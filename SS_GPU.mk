MY_SUBDIRS = src \
          src/io
# MY_OBJS := $(addsuffix /*.o, ${MY_SUBDIRS})
MY_OBJS := $(shell find $(SS_GPU_DIR) -name '*.o')
MY_INCLUDE = -I${SS_GPU_DIR}/src \
	-I${SS_GPU_DIR}/src/io

%.o: %.cpp
	${CC} ${CFLAGS} ${MY_INCLUDE} -std=c++11 -c $< -o $*.o