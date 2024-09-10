include ../../Makefile.in
include ../../A2D_SHELLS.mk

default: LoadMesh.o RunCRM.o
	$(CC) $(OFLAGS) ${MY_OBJS} LoadMesh.o -o loadmesh
	# $(CC) $(OFLAGS) ${MY_OBJS} RunCRM.o -o runcrm

clean:
	rm -f *.o
