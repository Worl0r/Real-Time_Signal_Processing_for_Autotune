### Do not edit -- Generated by 'configure --with-whatever' from Makefile.in
### RtAudio library Makefile

RM = /bin/rm
LN = /bin/ln

OBJECTS	= RtAudio.o 

LIBNAME = librtaudio
STATIC = $(LIBNAME).a
SHARED = librtaudio.dylib
RELEASE = 4.1.1
MAJOR = 4
LIBRARIES = $(STATIC) $(SHARED)

CC       = g++
AR       = /usr/bin/ar
RANLIB   = ranlib

DEFS     =   -DHAVE_GETTIMEOFDAY -D__MACOSX_CORE__
CFLAGS   = -O2 -Wall -Wextra -Iinclude -fPIC

PREFIX   = /usr/local

all : $(LIBRARIES)

tests:
	cd tests && $(MAKE) all

$(LIBRARIES): $(OBJECTS)
	$(AR) ruv $(STATIC) $(OBJECTS)
	ranlib $(STATIC)
	$(CC) -fPIC -dynamiclib -o librtaudio.$(RELEASE).dylib $(OBJECTS) -lpthread -framework CoreAudio -framework CoreFoundation
	$(LN) -sf librtaudio.$(RELEASE).dylib $(SHARED)
	$(LN) -sf librtaudio.$(RELEASE).dylib $(SHARED).$(MAJOR)

%.o : %.cpp
	$(CC) $(CFLAGS) $(DEFS) -c $(<) -o $@

%.o : include/%.cpp
	$(CC) $(CFLAGS) $(DEFS) -c $(<) -o $@

install:
	install --mode=755 $(STATIC) $(PREFIX)/lib/
	install --mode=755 librtaudio.$(RELEASE).dylib $(PREFIX)/lib/
	$(LN) -sf librtaudio.$(RELEASE).dylib $(PREFIX)/lib/$(SHARED)
	$(LN) -sf librtaudio.$(RELEASE).dylib $(PREFIX)/lib/$(SHARED).$(MAJOR)
	install --mode=644 $(LIBNAME).pc $(PREFIX)/lib/pkgconfig
	install --mode=644 RtAudio.h $(PREFIX)/include/
	install --mode=755 rtaudio-config $(PREFIX)/bin/

uninstall:
	-@rm -vf $(patsubst %,$(PREFIX)/lib/%, $(LIBRARIES) $(SHARED).$(MAJOR) $(SHARED).$(RELEASE))
	-@rm -vf $(PREFIX)/lib/pkgconfig/$(LIBNAME).pc
	-@rm -vf $(PREFIX)/bin/rtaudio-config

clean : 
	$(RM) -f $(LIBRARIES) librtaudio.$(RELEASE).dylib $(SHARED)*
	$(RM) -f $(OBJECTS)
	$(RM) -f *~
	cd tests && $(MAKE) clean

distclean:
	$(RM) -f $(LIBRARIES) librtaudio.$(RELEASE).dylib $(SHARED)*
	$(RM) -f $(OBJECTS)
	$(RM) -f *~
	$(RM) -rf config.log config.status autom4te.cache Makefile rtaudio-config $(LIBNAME).pc
	cd tests && $(MAKE) distclean

strip : 
	strip $(LIBRARIES)
	ranlib $(LIBRARIES)
	cd tests && $(MAKE) strip

.PHONY: clean distclean strip install uninstall
