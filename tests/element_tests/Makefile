include ../../Makefile.in
include ../../A2D_SHELLS.mk

default: make_shell_element.o make_shell

make_shell:
	$(CC) $(OFLAGS) ${MY_OBJS} make_shell_element.o -o make_shell.x

new:
	$(MAKE) clean
	$(MAKE) default

clean:
	rm -f *.o
	rm -f *.x
