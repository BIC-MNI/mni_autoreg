INCLUDES = -I$(srcdir)/../Include -I$(top_srcdir)/Proglib @MNI_INCLUDES@

LDADD = ../Files/libminctracc_files.a \
	../Optimize/libminctracc_optimize.a \
	../Numerical/libminctracc_numerical.a \
	../Volume/libminctracc_volume.a \
	../../Proglib/libProglib.a \
	@MNI_LIBS@ -lm

bin_PROGRAMS = minctracc
minctracc_SOURCES = \
	make_matlab_data_file.c \
	minctracc.c 

EXTRA_DIST = measure_code.c

