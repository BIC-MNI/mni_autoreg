#
# Makefile for the Perl part of the MNI AutoReg package.  This file assumes
# that it lives in the perl/ subdirectory of the mni_autoreg package,
# i.e. it depends on having file ../Makefile.include to get directories for
# the installation.  (SCRIPTDIR for the scripts, PERLLIBDIR for the Perl
# library, and CONFIGDIR for the config files.)
#

include ../Makefile.include

SHELL = /bin/sh

# The first list of programs here is all the Perl that will eventually be
# in the package; for now, we'll just stick with what's needed for mritotal

#PROGRAMS = stripscalp autocrop assemble_pets\
#           mritopet mritoself mritotal mritotal\
#           do_mritopet do_pettopet do_mritotal
PROGRAMS =  autocrop\
            mritotal\
            mritoself

MANPAGES = mritotal.1\
           autocrop.1\
           mritoself.1

CONFIGS  = mritotal.cfg\
           mritotal.default.cfg\
           mritotal.icbm.cfg

# All this stuff should go in the perl5/site_perl directory.

PERL_LIB = Startup.pm JobControl.pm\
           misc_utilities.pl\
           file_utilities.pl\
           path_utilities.pl\
           numeric_utilities.pl\
           minc_utilities.pl\
           volume_heuristics.pl

# Now the targets, which are specific to this directory (ie.
# we can't just include Louis' Makefile.bottom, because his
# install target doesn't do enough -- here, we need to install
# the Perl library stuff too)

build:	$(MANPAGES) version
	@if test -z "$(PERL)"; then \
	  echo "Perl 5 wasn't found on your system -- skipping this directory"; \
	  echo "(N.B. \"make install\" will fail in this directory -- just ignore it)" ; \
	else \
	  $(PERL) build_perl_stuff \
	    PERL=$(PERL) \
	    MODELDIR=$(MODELDIR) \
	    CONFIGDIR=$(CONFIGDIR) \
	    $(PROGRAMS) $(CONFIGS) ; \
	  chmod 755 $(PROGRAMS) ; \
	  chmod 644 $(CONFIGS) $(PERL_LIB) $(MANPAGES) ; \
	fi

version: version.c ../version.h
	$(CC) -o version version.c

install:
	@if [ ! -d $(SCRIPTDIR) ] ; then mkdir -p $(SCRIPTDIR) ; fi
	@if [ ! -d $(CONFIGDIR) ] ; then mkdir -p $(CONFIGDIR) ; fi
	@if [ ! -d $(PERLLIBDIR) ] ; then mkdir -p $(PERLLIBDIR) ; fi

	-cp $(PROGRAMS) $(SCRIPTDIR)
	-cp $(CONFIGS) $(CONFIGDIR)
	-cp $(PERL_LIB) $(PERLLIBDIR)
	-@$(PERL) checkdir $(PERLLIBDIR)

installman: $(MANPAGES)
	cp -p $(MANPAGES) $(MAN1DIR)

clean:
	\rm -f $(PROGRAMS) $(CONFIGS) $(MANPAGES) version

