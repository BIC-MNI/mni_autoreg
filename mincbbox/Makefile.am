INCLUDES = -I$(top_srcdir)/Proglib @INCLUDEDIRS@

LDADD = ../Proglib/libProglib.a @LIBDIRS@ \
	-lvolume_io -lminc -lnetcdf -lm

bin_PROGRAMS = mincbbox
mincbbox_SOURCES = mincbbox.c mincbbox.h

