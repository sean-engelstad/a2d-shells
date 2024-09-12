include ../../Makefile.in
include ../../A2D_SHELLS.mk

default: LoadMesh.o loadmesh RunCRM.o runcrm

loadmesh:
	$(CC) $(OFLAGS) ${MY_OBJS} LoadMesh.o -o loadmesh

runcrm:
	$(CC) $(OFLAGS) ${MY_OBJS} RunCRM.o -o runcrm

new:
	$(MAKE) clean
	$(MAKE) default

clean:
	rm -f *.o
	rm loadmesh
	rm runcrm
