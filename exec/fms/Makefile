# Makefile created by mkmf 19.3.1

CPPDEFS = -DINTERNAL_FILE_NML -Duse_libMPI -Duse_netCDF


include $(MK_TEMPLATE)
include $(MK_TEMPLATE)

.DEFAULT:
	-echo $@ does not exist.
all: libfms.a

libfms.a:
	mkdir -p build && cd build && autoreconf -i ${SRCROOT}/FMS 
	cd build && ${SRCROOT}/FMS/configure FC="$(FC)" CC="$(CC)" CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)" CPPFLAGS="$(CPPFLAGS)" FPPFLAGS="$(FPPFLAGS)" 
	cd build && make -j 32

clean:
	rm -rf build
