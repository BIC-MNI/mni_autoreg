perllibdir = $(pkglibdir)
configdir = $(sysconfdir)/@PACKAGE@
modeldir = $(pkgdatadir)

man_MANS = mritotal.1 autocrop.1 mritoself.1
script_in = autocrop.in mritoself.in mritotal.in xfmtool.in
configdata_in = mritotal.cfg.in mritotal.default.cfg.in mritotal.icbm.cfg.in

EXTRA_DIST = $(script_in) $(configdata_in) $(man_MANS)

if PERLGOOD

CLEANFILES = $(bin_SCRIPTS) $(config_DATA)

bin_SCRIPTS = $(script_in:.in=)
config_DATA  = $(configdata_in:.in=)

install-exec-local:

else

install-exec-local:
	@echo 'Perl (version 5.002 or greater) was not found on your system'
	@echo '**** NO PERL SCRIPTS INSTALLED ****'
	@exit 1

endif
