INCLUDES = -I$(srcdir)/../Include -I$(top_srcdir)/Proglib @INCLUDEDIRS@
LDADD = \
	../Optimize/libminctracc_optimize.a \
	../Numerical/libminctracc_numerical.a \
	../Volume/libminctracc_volume.a \
	../../Proglib/libProglib.a \
	@LIBDIRS@ \
	-lvolume_io -lminc -lnetcdf -lm

bin_PROGRAMS = \
	check_scale \
	param2xfm \
	volume_cog \
	rand_param \
	xcorr_vol \
	xfm2param \
	zscore_vol

