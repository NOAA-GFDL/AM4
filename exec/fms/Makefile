# Makefile created by mkmf 19.3.1

CPPDEFS = -DINTERNAL_FILE_NML -Duse_libMPI -Duse_netCDF -DMAXFIELDMETHODS_=500 -DMAXFIELDS_=500


include $(MK_TEMPLATE)

.DEFAULT:
	-echo $@ does not exist.
all: libFMS.a

libFMS.a: folder configure build 
	echo path to libFMS is ${PWD}/build/libFMS/.libs/libFMS.a
folder:
	mkdir -p build && cd build && autoreconf -i ${SRCROOT}/FMS 
configure: folder
	cd build && ${SRCROOT}/FMS/configure FC="$(FC)" CC="$(CC)" FCFLAGS="$(FFLAGS)" CFLAGS="$(CFLAGS)" CPPFLAGS="$(CPPFLAGS) $(CPPDEFS) $(FPPFLAGS)" FPPFLAGS="$(FPPFLAGS)" LIBS="$(LIBS)" 
build: configure
	cd build && make
clean:
	rm -rf build *.a* *.log
