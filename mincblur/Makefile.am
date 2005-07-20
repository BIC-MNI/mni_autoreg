INCLUDES = -I$(top_srcdir)/Proglib

LDADD = ../Proglib/libProglib.a -lm

bin_PROGRAMS = mincblur
man_MANS = mincblur.1

EXTRA_DIST = $(man_MANS)

mincblur_SOURCES = \
	apodize_data.c \
	blur_support.c blur_support.h \
	blur_volume.c blur_volume.h \
	fft.c \
	gradient_volume.c \
	gradmag_volume.c gradmag_volume.h \
	kernel.h \
	mincblur.c mincblur.h

