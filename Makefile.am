AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

bin_SCRIPTS = make_model

SUBDIRS = \
	Proglib \
	make_phantom \
	mincbbox \
	mincblur \
	mincchamfer \
	minctracc \
	perl \
	Testing

EXTRA_DIST = $(m4_files) \
	BUGS FAQ INSTALL.quick README.TESTING \
	Doc/howto_add_obj_fn Doc/howto_add_opt_fn \
	make_model.in \
	epm-header.in

m4_files = \
	m4/mni_REQUIRE_LIB.m4 \
	m4/mni_REQUIRE_MNILIBS.m4 \
	m4/mni_REQUIRE_OPENINVENTOR.m4 \
	m4/mni_cxx_have_koenig_lookup.m4 \
	m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 \
	m4/smr_REQUIRED_LIB.m4 \
	m4/smr_WITH_BUILD_PATH.m4
