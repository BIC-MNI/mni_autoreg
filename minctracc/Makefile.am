SUBDIRS = Files Numerical Optimize Volume Main Extra_progs

# Automake does not handle sources in subdirs (yet).  We kludge our
# way around this by making a bunch of (noninstalled) libraries ---
# one per subdir --- and building the minctracc program in
# subdirectory Main.  The extra programs are likewise built in their
# subdirectory, Extra_progs.

EXTRA_DIST = $(man_MANS) $(docs) $(includes)

man_MANS = minctracc.1

docs = Doc/doc

includes = \
	Include/amoeba.h \
	Include/arg_data.h \
	Include/constants.h \
	Include/cov_to_praxes.h \
	Include/deform_support.h \
	Include/extras.h \
	Include/globals.h \
	Include/init_lattice.h \
	Include/interpolation.h \
	Include/local_macros.h \
	Include/make_rots.h \
	Include/matrix_basics.h \
	Include/minctracc.h \
	Include/objectives.h \
	Include/point_vector.h \
	Include/quad_max_fit.h \
	Include/quaternion.h \
	Include/rotmat_to_ang.h \
	Include/segment_table.h \
	Include/stats.h \
	Include/sub_lattice.h \
	Include/super_sample_def.h \
	Include/vox_space.h

