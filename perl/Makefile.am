perllibdir = $(pkglibdir)
configdir = $(sysconfdir)/@PACKAGE@
modeldir = $(pkgdatadir)

# Helper program used by build_perl_stuff.
noinst_PROGRAMS = version

# In future:
#PROGRAMS = stripscalp autocrop assemble_pets\
#           mritopet mritoself mritotal mritotal\
#           do_mritopet do_pettopet do_mritotal

EXTRA_DIST = $(script_in) $(configdata_in) $(man_MANS) $(perllib_DATA) \
	build_perl_stuff 

script_in = autocrop.in mritoself.in mritotal.in xfmtool.in
configdata_in = mritotal.cfg.in mritotal.default.cfg.in mritotal.icbm.cfg.in

man_MANS = mritotal.1 autocrop.1 mritoself.1

perllib_DATA = Startup.pm JobControl.pm \
	misc_utilities.pl \
	file_utilities.pl \
	path_utilities.pl \
	numeric_utilities.pl \
	minc_utilities.pl \
	volume_heuristics.pl


if PERLGOOD

CLEANFILES = $(bin_SCRIPTS) $(config_DATA)

bin_SCRIPTS = autocrop mritotal mritoself xfmtool
config_DATA  = mritotal.cfg mritotal.default.cfg mritotal.icbm.cfg

# Wacky perl script to do substitution
#
$(bin_SCRIPTS) $(config_DATA): $(script_in) $(configdata_in) version Makefile
	$(PERL) $(srcdir)/build_perl_stuff --srcdir $(srcdir) \
	    PERL=$(PERL) LIBDIR=$(perllibdir) \
	    MODELDIR=$(modeldir) CONFIGDIR=$(configdir) \
	    $(bin_SCRIPTS) $(config_DATA)


else

install-exec-local:
	@echo 'Perl (version 5.002 or greater) was not found on your system'
	@echo '**** NO PERL SCRIPTS INSTALLED ****'
	@exit 1

endif
