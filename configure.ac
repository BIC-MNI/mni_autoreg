# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

# The arguments are package name, and package version.
AC_INIT([mni_autoreg],[0.99.6],[Louis Collins <louis@bic.mni.mcgill.ca>])
AC_CONFIG_SRCDIR([minctracc/Include/minctracc.h])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB

# check for various bits and pieces
AC_C_INLINE
AC_C_CONST
AC_TYPE_SIZE_T
AC_CHECK_HEADERS(float.h limits.h malloc.h math.h stdlib.h)

# Checks for libraries.  See m4/README.
mni_REQUIRE_VOLUMEIO

# for clean MINC2.0 volume_io
AC_DEFINE_UNQUOTED(VIO_PREFIX_NAMES, 1, [Play nice with the other kids volume_io])

# for louis' confangulated history string
AC_CANONICAL_HOST

compile_datetime=`date`
compile_user=`logname`
compile_host=`hostname`
long_version="Package $PACKAGE_STRING, compiled by $compile_user@$compile_host ($host) on $compile_datetime"

AC_DEFINE_UNQUOTED(MNI_AUTOREG_COMPILE_DATETIME,"$compile_datetime",[Configure DateTime])
AC_DEFINE_UNQUOTED(MNI_AUTOREG_COMPILE_USER,"$compile_user",[Configure User])
AC_DEFINE_UNQUOTED(MNI_AUTOREG_COMPILE_HOST,"$compile_host",[Configure Host])
AC_DEFINE_UNQUOTED(MNI_AUTOREG_COMPILE_SYSTEM,"$host",[Configure System])
AC_DEFINE_UNQUOTED(MNI_AUTOREG_LONG_VERSION,"$long_version",[Long Version])

AC_SUBST(MNI_AUTOREG_COMPILE_DATETIME, [$compile_datetime])
AC_SUBST(MNI_AUTOREG_COMPILE_USER, [$compile_user])
AC_SUBST(MNI_AUTOREG_COMPILE_HOST, [$compile_host])
AC_SUBST(MNI_AUTOREG_COMPILE_SYSTEM, [$host])
AC_SUBST(MNI_AUTOREG_LONG_VERSION, [$long_version])

#------------------------------------------------------------------------
# Look for perl (and make sure it's the right version)
#------------------------------------------------------------------------
AC_PATH_PROGS(PERL, [perl], no)

if test "$PERL" = no; then
    AC_MSG_WARN(The mritotal script will not run since Perl 5 is not found.
                You must install Perl 5 before installing MNI AutoReg.)
    PERL=
else
    AC_MSG_CHECKING(Perl version)

    # The perl command will exit with status zero if
    #  (i) $PERL is found, is executable, and
    # (ii) the perl version is NOT < 5.002
    #
    if $PERL -e 'exit ($] < 5.002)'; then
	AC_MSG_RESULT(5.002 or later--good)
    else
	AC_MSG_RESULT(unacceptable)
	AC_MSG_WARN(Perl 5.002 (or greater) is required--you will not be able to run mritotal)
	PERL=
    fi
fi

AM_CONDITIONAL(PERLGOOD, test x"$PERL" != x)


# Finish up by writing output files.
# AC_CONFIG_FILES([Makefile])
AC_OUTPUT(epm-header
make_model
perl/autocrop
perl/mritoself
perl/mritotal
perl/xfmtool
perl/mritotal.cfg
perl/mritotal.default.cfg
perl/mritotal.icbm.cfg
Makefile
Proglib/Makefile
make_phantom/Makefile
mincbbox/Makefile
mincblur/Makefile
mincchamfer/Makefile
minctracc/Makefile
minctracc/Files/Makefile
minctracc/Main/Makefile
minctracc/Numerical/Makefile
minctracc/Optimize/Makefile
minctracc/Volume/Makefile
minctracc/Extra_progs/Makefile
perl/Makefile
Testing/Makefile
)
