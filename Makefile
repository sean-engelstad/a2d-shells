include Makefile.in.info
include SS_GPU.mk


default: ${MY_OBJS}
	for subdir in $(MY_SUBDIRS) ; do \
		echo "making $@ in $$subdir"; \
		echo; (cd $$subdir && $(MAKE)) || exit 1; \
		done \

clean:
	rm -f *.o
	rm -f **/*.o
	rm -f **/**/*.o