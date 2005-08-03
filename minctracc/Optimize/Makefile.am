INCLUDES = -I$(srcdir)/../Include -I$(top_srcdir)/Proglib

noinst_LIBRARIES = libminctracc_optimize.a
libminctracc_optimize_a_SOURCES = \
	def_obj_functions.c \
	extras.c \
	sub_lattice.c\
	amoeba.c \
	vox_space.c \
	objectives.c \
	optimize.c \
	segment_table.c \
	deform_support.c \
	super_sample_def.c \
	my_grid_support.c \
	obj_fn_mutual_info.c \
	do_nonlinear.c

EXTRA_DIST = switch_obj_func.c \
	louis_splines.h
