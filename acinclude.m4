dnl ####################### -*- Mode: M4 -*- ###########################
dnl smr.m4 -- 
dnl 
dnl Copyright (C) 1999 Matthew D. Langston <langston@SLAC.Stanford.EDU>
dnl Copyright (C) 1998 Steve Robbins <stever@cs.mcgill.ca>
dnl
dnl This file is free software; you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl
dnl This file is distributed in the hope that it will be useful, but
dnl WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this file; if not, write to:
dnl
dnl   Free Software Foundation, Inc.
dnl   Suite 330
dnl   59 Temple Place
dnl   Boston, MA 02111-1307, USA.
dnl ####################################################################


dnl NOTE: The macros in this file are extensively documented in the
dnl       accompanying `smr_macros.texi' Texinfo file.  Please see the
dnl       Texinfo documentation for the definitive specification of how
dnl       these macros are supposed to work.  If the macros work
dnl       differently than the Texinfo documentation says they should,
dnl       then the macros (and not the Texinfo documentation) has the
dnl       bug(s).

dnl This is a convenient macro which translates illegal characters for
dnl bourne shell variables into legal characters.  It has the same
dnl functionality as sed 'y%./+-:%__p__%'.
AC_DEFUN([smr_safe_translation], [patsubst(patsubst([$1], [+], [p]), [./-:], [_])])

AC_DEFUN(smr_SWITCH,
[
  dnl Define convenient aliases for the arguments since there are so
  dnl many of them and I keep confusing myself whenever I have to edit
  dnl this macro.
  pushdef([smr_name],        $1)
  pushdef([smr_help_string], $2)
  pushdef([smr_default],     $3)
  pushdef([smr_yes_define],  $4)
  pushdef([smr_no_define],   $5)

  dnl Do some sanity checking of the arguments.
  ifelse([regexp(smr_default, [^\(yes\|no\)$])], -1, [AC_MSG_ERROR($0: third arg must be either yes or no)])

  dnl Create the help string
  pushdef([smr_lhs], [--ifelse(smr_default, yes, disable, enable)-smr_name])dnl
  pushdef([smr_rhs], [ifelse(smr_default, yes, disable, enable) smr_help_string (default is smr_default)])dnl

  AC_PRETTY_HELP_STRING([smr_lhs], [smr_rhs], smr_name[]_switch_help_string)

  dnl Add the option to `configure --help'.  We don't need to supply the
  dnl 4th argument to AC_ARG_ENABLE (i.e. the code to set the default
  dnl value) because that is done below by AC_CACHE_CHECK.
  AC_ARG_ENABLE([smr_name],
                [$]smr_name[]_switch_help_string,
                smr_cv_enable_[]smr_name=$enableval)

  dnl We cache the result so that the user doesn't have to remember
  dnl which flags they passed to `configure'.
  AC_CACHE_CHECK([whether to enable smr_help_string],
                 smr_cv_enable_[]smr_name,
                 smr_cv_enable_[]smr_name=smr_default)

    ifelse(smr_yes_define, , , test x"[$]smr_cv_enable_[]smr_name" = xyes && AC_DEFINE(smr_yes_define))
    ifelse(smr_no_define, , ,  test x"[$]smr_cv_enable_[]smr_name" = xno  && AC_DEFINE(smr_no_define))

  dnl Sanity check the value assigned to smr_cv_enable_$1 to force it to
  dnl be either `yes' or `no'.
  if test ! x"[$]smr_cv_enable_[]smr_name" = xyes; then
    if test ! x"[$]smr_cv_enable_[]smr_name" = xno; then
      AC_MSG_ERROR([smr_lhs must be either yes or no])
    fi
  fi

  popdef([smr_name])
  popdef([smr_help_string])
  popdef([smr_default])
  popdef([smr_yes_define])
  popdef([smr_no_define])
  popdef([smr_lhs])
  popdef([smr_rhs])
])


AC_DEFUN(smr_ARG_WITHLIB,
[
  dnl Define convenient aliases for the arguments since there are so
  dnl many of them and I keep confusing myself whenever I have to edit
  dnl this macro.
  pushdef([smr_name],        $1)
  pushdef([smr_libname],     ifelse($2, , smr_name, $2))
  pushdef([smr_help_string], $3)
  pushdef([smr_safe_name],   smr_safe_translation(smr_name))

  dnl Create the help string
  AC_PRETTY_HELP_STRING([--with-smr_safe_name-library[[=PATH]]],
                        [use smr_name library ifelse(smr_help_string, , , (smr_help_string))],
                        smr_safe_name[]_library_help_string)

  dnl Add the option to `configure --help'.  We don't need to supply the
  dnl 4th argument to AC_ARG_WITH (i.e. the code to set the default
  dnl value) because that is done below by AC_CACHE_CHECK.
  AC_ARG_WITH(smr_safe_name-library,
              [$]smr_safe_name[]_library_help_string,
              smr_cv_with_[]smr_safe_name[]_library=$withval)

  dnl We cache the result so that the user doesn't have to remember
  dnl which flags they passed to `configure'.
  AC_CACHE_CHECK([whether to use smr_name library],
                 smr_cv_with_[]smr_safe_name[]_library,
                 smr_cv_with_[]smr_safe_name[]_library=maybe)


  case x"[$]smr_cv_with_[]smr_safe_name[]_library" in
      xyes | xmaybe)
          smr_safe_name[]_LIBS="-l[]smr_libname"
          with_[]smr_safe_name=[$]smr_cv_with_[]smr_safe_name[]_library
          ;;
      xno)
          smr_safe_name[]_LIBS=
          with_[]smr_safe_name=no
          ;;
      *)
          if test -f "[$]smr_cv_with_[]smr_safe_name[]_library"; then
            smr_safe_name[]_LIBS=[$]smr_cv_with_[]smr_safe_name[]_library
          elif test -d "[$]smr_cv_with_[]smr_safe_name[]_library"; then
            smr_safe_name[]_LIBS="-L[$]smr_cv_with_[]smr_safe_name[]_library -l[]smr_libname"
          else
            AC_MSG_ERROR([argument must be boolean, file, or directory])
          fi
          with_[]smr_safe_name=yes
          ;;
  esac

  popdef([smr_name])
  popdef([smr_libname])
  popdef([smr_help_string])
  popdef([smr_safe_name])
])


AC_DEFUN(smr_ARG_WITHINCLUDES,
[
  dnl Define convenient aliases for the arguments since there are so
  dnl many of them and I keep confusing myself whenever I have to edit
  dnl this macro.
  pushdef([smr_name],        $1)
  pushdef([smr_header],      $2)
  pushdef([smr_extra_flags], $3)
  pushdef([smr_safe_name],   smr_safe_translation(smr_name))

  dnl Create the help string
  AC_PRETTY_HELP_STRING([--with-smr_safe_name-includes[[=DIR]]],
                        [set directory for smr_name headers],
                        smr_safe_name[]_includes_help_string)

  dnl Add the option to `configure --help'.  We don't need to supply the
  dnl 4th argument to AC_ARG_WITH (i.e. the code to set the default
  dnl value) because that is done below by AC_CACHE_CHECK.
  AC_ARG_WITH(smr_safe_name-includes,
              [$]smr_safe_name[]_includes_help_string,
              smr_cv_with_[]smr_safe_name[]_includes=$withval)

  dnl We cache the result so that the user doesn't have to remember
  dnl which flags they passed to `configure'.
  AC_CACHE_CHECK([where to find the smr_name header files],
                 smr_cv_with_[]smr_safe_name[]_includes,
                 smr_cv_with_[]smr_safe_name[]_includes=)

  if test ! x"[$]smr_cv_with_[]smr_safe_name[]_includes" = x; then
    if test -d "[$]smr_cv_with_[]smr_safe_name[]_includes"; then
      smr_safe_name[]_CFLAGS="-I[$]smr_cv_with_[]smr_safe_name[]_includes"
    else
      AC_MSG_ERROR([argument must be a directory])
    fi
  else
    smr_safe_name[]_CFLAGS=
  fi

  dnl This bit of logic comes from the autoconf AC_PROG_CC macro.  We
  dnl need to put the given include directory into CPPFLAGS temporarily,
  dnl but then restore CPPFLAGS to its old value.
  smr_test_CPPFLAGS="${CPPFLAGS+set}"
  smr_save_CPPFLAGS="$CPPFLAGS"
  CPPFLAGS="$CPPFLAGS [$]smr_safe_name[]_CFLAGS smr_extra_flags"

  dnl If the header file smr_header exists, then define
  dnl HAVE_[]smr_header (in all capitals).
  AC_CHECK_HEADERS([smr_header],
                   smr_have_[]smr_safe_name[]_header=yes,
                   smr_have_[]smr_safe_name[]_header=no)

  if test x"$smr_test_CPPFLAGS" = xset; then
    CPPFLAGS=$smr_save_CPPFLAGS
  else
    unset CPPFLAGS
  fi

  popdef([smr_name])
  popdef([smr_header])
  popdef([smr_extra_flags])
  popdef([smr_safe_name])
])


AC_DEFUN(smr_CHECK_LIB,
[
  dnl Define convenient aliases for the arguments since there are so
  dnl many of them and I keep confusing myself whenever I have to edit
  dnl this macro.
  pushdef([smr_name],        $1)
  pushdef([smr_libname],     ifelse($2, , smr_name, $2))
  pushdef([smr_help_string], $3)
  pushdef([smr_function],    $4)
  pushdef([smr_header],      $5)
  pushdef([smr_extra_libs],  $6)
  pushdef([smr_extra_flags], $7)
  pushdef([smr_prototype],   $8)
  pushdef([smr_safe_name],   smr_safe_translation(smr_name))

  dnl Give the user (via "configure --help") an interface to specify
  dnl whether we should use the library or not, and possibly where we
  dnl should find it.
  smr_ARG_WITHLIB([smr_name], [smr_libname], [smr_help_string])

  if test ! x"$with_[]smr_safe_name" = xno; then

    # If we got this far, then the user didn't explicitly ask not to use
    # the library.

    dnl If the caller of smr_CHECK_LIB specified a header file for this
    dnl library, then give the user (via "configure --help") an
    dnl interface to specify where this header file can be found (if it
    dnl isn't found by the compiler by default).
    ifelse(smr_header, , , [smr_ARG_WITHINCLUDES(smr_name, smr_header, smr_extra_flags)])

    # We need only look for the library if the header has been found
    # (or no header is needed).
    if test [$]smr_have_[]smr_safe_name[]_header != no; then

    	AC_CHECK_LIB(smr_libname,
                     smr_function,
                     smr_have_[]smr_safe_name[]_library=yes,
                     smr_have_[]smr_safe_name[]_library=no,
                     [$]smr_safe_name[]_CFLAGS [smr_extra_flags] [$]smr_safe_name[]_LIBS [smr_extra_libs],
                     [ifelse(smr_prototype, , , [[#]include <smr_header>])],
                     smr_prototype)
    fi

    if test x"[$]smr_have_[]smr_safe_name[]_library" = xyes; then
      AC_MSG_RESULT([using smr_name library])
    else
      smr_safe_name[]_LIBS=
      smr_safe_name[]_CFLAGS=

      if test x"$with_[]smr_safe_name" = xmaybe; then
        AC_MSG_RESULT([not using smr_name library])
      else
        AC_MSG_WARN([requested smr_name library not found!])
      fi
    fi
  fi

  popdef([smr_name])
  popdef([smr_libname])
  popdef([smr_help_string])
  popdef([smr_function])
  popdef([smr_header])
  popdef([smr_extra_libs])
  popdef([smr_extra_flags])
  popdef([smr_prototype])
  popdef([smr_safe_name])
])


dnl NOTE: The AC_PRETTY_HELP_STRING macro (below) has been submitted to
dnl       the Autoconf maintainer, but it hasn't shown up yet in a new
dnl       Autoconf release.  Therefore, I am including it with the
dnl       smr_macros for the time being since they depend on
dnl       AC_PRETTY_HELP_STRING being available.


dnl AC_PRETTY_HELP_STRING
dnl ---------------------
dnl
dnl usage: AC_PRETTY_HELP_STRING(LHS, RHS, HELP-STRING)
dnl
dnl Format an Autoconf macro's help string so that it looks pretty when
dnl the user executes "configure --help".  This macro take three
dnl arguments, a "left hand side" (LHS), a "right hand side" (RHS), and
dnl a variable (HELP-STRING) to set to the pretty-printed concatenation
dnl of LHS and RHS (the new, pretty-printed "help string").
dnl
dnl The resulting string in HELP-STRING is suitable for use in other
dnl macros that require a help string (e.g. AC_ARG_WITH).
dnl 
AC_DEFUN(AC_PRETTY_HELP_STRING,
[
dnl 
dnl Here is the sample string from the Autoconf manual (Node: External
dnl Software) which shows the proper spacing for help strings.
dnl 
dnl    --with-readline         support fancy command line editing
dnl  ^ ^                       ^ 
dnl  | |                       |
dnl  | column 2                column 26     
dnl  |
dnl  column 0
dnl 
dnl A help string is made up of a "left hand side" (LHS) and a "right
dnl hand side" (RHS).  In the example above, the LHS is
dnl "--with-readline", while the RHS is "support fancy command line
dnl editing".
dnl 
dnl If the LHS extends past column 24, then the LHS is terminated with a
dnl newline so that the RHS is on a line of its own beginning in column
dnl 26.
dnl 
dnl Therefore, if the LHS were instead "--with-readline-blah-blah-blah",
dnl then the AC_PRETTY_HELP_STRING macro would expand into:
dnl
dnl
dnl    --with-readline-blah-blah-blah
dnl  ^ ^                       support fancy command line editing
dnl  | |                       ^ 
dnl  | column 2                |
dnl  column 0                  column 26     

dnl We divert everything to AC_DIVERSION_NOTICE (which gets output very
dnl early in the configure script) because we want the user's help
dnl string to be set before it is used.

AC_DIVERT_PUSH(AC_DIVERSION_NOTICE)dnl
# This is from AC_PRETTY_HELP_STRING
lhs="$1"
rhs="$2"

lhs_column=25
rhs_column=`expr $lhs_column + 1`
lhs_pad='                         '      # 25 spaces
rhs_pad='                          '     # 26 spaces


# Insure that the LHS begins with exactly two spaces.
changequote(, )dnl
lhs=`echo "$lhs" | sed -n -e "s/[ ]*\(.*\)/  \1/p"`
changequote([, ])dnl

# Is the length of the LHS greater or equal to the $lhs_column?
if `echo "$lhs" | grep ".\{$lhs_column\}" > /dev/null 2>&1`; then
  # Strip all leading spaces from the RHS.
changequote(, )dnl
  rhs=`echo "$rhs" | sed -n -e "s/[ ]*\(.*\)/\1/p"`
changequote([, ])dnl

$3="$lhs
${rhs_pad}${rhs}"

else
  lhs="${lhs}${lhs_pad}"
changequote(, )dnl
$3=`echo "$lhs" | sed -n -e "/.\{$lhs_column\}[ ][ ]*$/ s/\(.\{$rhs_column\}\).*/\1$rhs/p"`
changequote([, ])dnl

fi 
AC_DIVERT_POP()dnl
])


dnl This is a small fix to the standard `AC_CHECK_LIB' macro which
dnl allows checking for C++ libraries.
dnl
dnl AC_CHECK_LIB(LIBRARY, FUNCTION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND
dnl              [, OTHER-LIBRARIES [, INCLUDES [, FUNC-PTR-PROTOTYPE]]]]])
pushdef([AC_CHECK_LIB],
[AC_MSG_CHECKING([for $2 in -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
ac_lib_var=`echo $1['_']$2 | sed 'y%./+-:%__p__%'`
AC_CACHE_VAL(ac_cv_lib_$ac_lib_var,
[ac_save_LIBS="$LIBS"
LIBS="-l$1 $5 $LIBS"
AC_TRY_LINK(dnl
ifelse([$6], , dnl If there is no INCLUDES then concoct our own declaration.
ifelse(AC_LANG, [FORTRAN77], ,
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
])), [$6]),
dnl If the user passed in a prototype for a pointer to FUNCTION
dnl (i.e. $7) then instead of just calling FUNCTION, construct a
dnl statement that takes the address of FUNCTION.  This is particularly
dnl helpful in checking for C++ class libraries.
ifelse([$7], ,[$2()], [$7 = &$2]),
	    eval "ac_cv_lib_$ac_lib_var=yes",
	    eval "ac_cv_lib_$ac_lib_var=no")
LIBS="$ac_save_LIBS"
])dnl
if eval "test \"`echo '$ac_cv_lib_'$ac_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$3], ,
[changequote(, )dnl
  ac_tr_lib=HAVE_LIB`echo $1 | sed -e 's/[^a-zA-Z0-9_]/_/g' \
    -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($ac_tr_lib)
  LIBS="-l$1 $LIBS"
], [$3])
else
  AC_MSG_RESULT(no)
ifelse([$4], , , [$4
])dnl
fi
])
