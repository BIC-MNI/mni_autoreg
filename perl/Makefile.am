perllibdir = $(pkglibdir)
configdir = $(sysconfdir)/@PACKAGE@
modeldir = $(libdir)/model

# Helper program used by build_perl_stuff.
noinst_PROGRAMS = version

# In future:
#PROGRAMS = stripscalp autocrop assemble_pets\
#           mritopet mritoself mritotal mritotal\
#           do_mritopet do_pettopet do_mritotal

man_MANS = mritotal.1 autocrop.1 mritoself.1

EXTRA_DIST = $(man_MANS) $(perllib_DATA) \
	build_perl_stuff \
	autocrop.in mritoself.in mritotal.in \
	mritotal.cfg.in mritotal.default.cfg.in mritotal.icbm.cfg.in 

perllib_DATA = Startup.pm JobControl.pm \
	misc_utilities.pl \
	file_utilities.pl \
	path_utilities.pl \
	numeric_utilities.pl \
	minc_utilities.pl \
	volume_heuristics.pl


if PERLGOOD

CLEANFILES = $(bin_SCRIPTS) $(config_DATA)

bin_SCRIPTS = autocrop mritotal mritoself
config_DATA  = mritotal.cfg mritotal.default.cfg mritotal.icbm.cfg

# Wacky perl script to do substitution
# which fails (sigh) if builddir != srcdir
#
$(bin_SCRIPTS) $(config_DATA): version Makefile
	$(PERL) $(srcdir)/build_perl_stuff \
	    PERL=$(PERL) LIBDIR=$(perllibdir) \
	    MODELDIR=$(modeldir) CONFIGDIR=$(configdir) \
	    $(bin_SCRIPTS) $(config_DATA)

else

install-exec-local:
	@echo 'Perl (version 5.002 or greater) was not found on your system'
	@echo '**** NO PERL SCRIPTS INSTALLED ****'
	@exit 1

endif
