INCLUDES = -I$(srcdir)/../Include -I$(top_srcdir)/Proglib
LDADD = ../Numerical/libminctracc_numerical.a \
	../Volume/libminctracc_volume.a \
	../Optimize/libminctracc_optimize.a \
	../../Proglib/libProglib.a

bin_PROGRAMS = \
	check_scale \
	crispify \
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

