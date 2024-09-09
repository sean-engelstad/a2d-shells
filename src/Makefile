OFLAGS = -Wall -g
# OFLAGS = -Wall "-Wl,-allow-multiple-definition" -g
# OFLAGS = -Wall "-Wl,-allow-multiple-definition" -g
CFLAGS = -c
CC = g++
exec = test_sym

# default: ${OBJECTS}
# 	echo "${OBJECTS}"
# 	$(CC) $(OFLAGS) ${OBJECTS} -o $(exec)

# default: LoadMesh.o TACSMeshLoader.o TACSObject.o
default: TACSObject.o TACSMeshLoader.o LoadMesh.o
	$(CC) $(OFLAGS) TACSObject.o TACSMeshLoader.o LoadMesh.o -o loadmesh

%.o: %.cpp
	${CC} ${CFLAGS} -std=c++11 -c $< -o $*.o

# %.o: %.cpp
# 	${CC} ${CFLAGS} -std=c++11 -c $<

clean:
	rm -f *.o
