INCLUDES = -I$(srcdir)/../Include -I$(top_srcdir)/Proglib @MNI_INCLUDES@
LDADD = ../Optimize/libminctracc_optimize.a \
	../Numerical/libminctracc_numerical.a \
	../Volume/libminctracc_volume.a \
	../../Proglib/libProglib.a \
	@MNI_LIBS@ -lm

bin_PROGRAMS = \
	check_scale \
	param2xfm \
	volume_cog \
	rand_param \
	xcorr_vol \
	xfm2param \
	zscore_vol

check_PROGRAMS = cmpxfm

EXTRA_DIST = $(TESTS)
CLEANFILES = test1.xfm test2.xfm

# The tests are all shell scripts.
TESTS_ENVIRONMENT = $(SHELL)
TESTS = param2xfm.tests

oddities = \
	avg_voxel_vol.c \
	change_range.c \
	correct_sd_vol.c \
	crispify.c \
	def_to_lines.c \
	def_to_mag.c \
	def_to_scale.c \
	def_to_vols.c \
	deform_object.c \
	find_tubes.c \
	inbetween_def.c \
	inbetween_tps.c \
	load_vol.c \
	make_lvv_vol.c \
	make_zscore_vol.c \
	mst_to_lines.c \
	overlap_corr.c \
	random_def.c \
	random_vol.c \
	randomize_tags.c \
	randomize_volume.c \
	reversedef.c \
	rms_dist_lines.c \
	rms_dist_objects.c \
	scan_objects_to_volume.c \
	test.c \
	test_cov.c \
	test_cov_pr.c \
	test_deriv.c \
	test_eigen.c \
	test_kernels.c \
	test_lines.c \
	test_mripet.csh \
	test_prindir.c \
	test_quad.c \
	volume_stat.c \
	xfminvert.c
