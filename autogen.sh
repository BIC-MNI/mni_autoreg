#! /bin/sh

cat <<EOF
Ignore the following three warnings:

	required file \`./config.h.in' not found
	warning: AC_TRY_RUN called without default to allow cross compiling
	warning: AC_TRY_RUN called without default to allow cross compiling



EOF

aclocal
automake --add-missing
autoconf
autoheader

