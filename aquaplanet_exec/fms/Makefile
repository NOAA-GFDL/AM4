# Makefile created by mkmf 19.3.1

CPPDEFS = -Duse_libMPI -Duse_netCDF -DINTERNAL_FILE_NML -DSPMD -DLAND_GRID_FROM_ATMOS -DLAND_BND_TRACERS -DGFDL_COMPATIBLE_MICROP 


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
