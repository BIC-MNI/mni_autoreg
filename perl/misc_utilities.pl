# ------------------------------ MNI Header ----------------------------------
#@NAME       : misc_utilities.pl
#@DESCRIPTION: Miscellaneous and unclassifiable (but otherwise useful!)
#              utility routines:
#                 &timestamp
#                 &userstamp
#                 &comp_num_lists
#                 &print_banner
#@GLOBALS    : 
#@CREATED    : 95/08/07, Greg Ward (from gpw_utilities.pl)
#@MODIFIED   : 
#@VERSION    : $Id: misc_utilities.pl,v 1.1.1.1 2000-01-19 14:10:32 louis Exp $
#-----------------------------------------------------------------------------

# --------------------------------------------------------------------
# Copyright (c) 1995 Greg Ward, McConnell Brain Imaging Centre,
# Montreal Neurological Institute, McGill University.  Permission to
# use, copy, modify, and distribute this software and its
# documentation for any purpose and without fee is hereby granted,
# provided that the above copyright notice appear in all copies.  The
# author and McGill University make no representations about the
# suitability of this software for any purpose.  It is provided "as
# is" without express or implied warranty.
# --------------------------------------------------------------------

use POSIX qw/strftime/;
use Cwd;

# used by &userstamp:

my $_host;
BEGIN { chop ($_host = `hostname`) }

# ------------------------------ MNI Header ----------------------------------
#@NAME       : timestamp
#@INPUT      : $tm - [optional] time to use, as seconds since 
#                    1970-01-01 00:00:00 UTC (eg from `time'); 
#                    defaults to the current time
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Generates and returns a timestamp of the form 
#              [95-05-16 22:30:14].
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/05/16, GW (from &doit)
#@MODIFIED   : 1996/05/22, GW: added seconds to time
#              1996/06/17, GW: changed to use strftime from POSIX
#-----------------------------------------------------------------------------
sub timestamp
{
   my ($tm) = @_;

   $tm = time unless defined $tm;
   strftime ('[%Y-%m-%d %H:%M:%S]', localtime ($tm));
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : userstamp
#@INPUT      : $user - [optional] username; defaults to looking up 
#                      login id of $< in password file using getpwuid
#              $host - [optional]; defaults to `hostname`
#              $dir  - [optional]; defaults to current directory, as 
#                      obtained by calling Cwd::getcwd
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Generates and returns a "userstamp" of the form 
#              [greg@bottom:/data/scratch1/greg].
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/05/16, GW
#@MODIFIED   : 1996/05/29, GW: added directory
#-----------------------------------------------------------------------------
sub userstamp
{
   my ($user, $host, $dir) = @_;

   $user = getpwuid ($<) unless defined $user;
#  chop ($host = `hostname`) unless defined $host;
   $host = $_host unless defined $host;
   $dir = getcwd unless defined $dir;
   sprintf ("[%s@%s:%s]", $user, $host, $dir);
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : comp_num_lists
#@INPUT      : $a1 - reference to first list
#              $a2 - reference to second list
#@OUTPUT     : 
#@RETURNS    : 0 if any mismatch between the two lists (i.e. they are of
#                different length, or any corresponding elements aren't
#                equal in the numeric sense)
#              1 if the two lists are numerically identical
#@DESCRIPTION: 
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 96/02/12 GPW (direct copyof comp_arrays in autocrop)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub comp_num_lists
{
   die "comp_num_lists: wrong number of arguments" unless (@_ == 2);
   local ($a1, $a2) = @_;
   
   return 0 unless (@$a1 == @$a2);
   for $i (0 .. $#$a1)
   {
      return 0 unless ($a1->[$i] == $a2->[$i]);
   }
   return 1;
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : make_banner
#@INPUT      : $msg    - the string to print
#              $char   - the character to use when making the "banner"
#                        (optional; defaults to "-")
#              $width  - the width of field to pad to (optional; defaults
#                        to 80, but should default to width of terminal)
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Creates and returns a string of the form 
#              "-- Hello! ----------" (assuming $msg="Hello!", $char="-", 
#              and $width=20)
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1996/05/22, Greg Ward - adapted from do_mritopet
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub make_banner
{
   my ($msg, $char, $width) = @_;

   $width = 80 unless $width;   # should use Term::Cap!!!
   $char = "-" unless $char;

   my $banner = $char x 2 . " " . $msg . " ";
   $banner .= $char x ($width - length ($banner)) . "\n"
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &shellquote
#@INPUT      : @words - list of words to possibly quote or escape
#@OUTPUT     : 
#@RETURNS    : concatenation of @words with necessary quotes and backslashes
#@DESCRIPTION: The inverse of shellwords -- takes a list of arguments 
#              (like @ARGV, or a list passed to system or exec) and 
#              escapes meta-characters or encases in quotes as appropriate
#              to allow later processing by the shell.  (/bin/sh, in 
#              particular -- the list of metacharacters was taken from
#              the Perl source that does an exec().)
#@METHOD     : If a word contains no metacharacters, it is untouched.  
#              If it contains both single and double quotes, all meta-
#              characters are escaped with a backslash, and no quotes 
#              are added.  If it contains just single quotes, it is encased
#              in double quotes.  Otherwise, it is encased in single quotes.
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1996/11/13, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub shellquote
{
   my (@words) = @_;
   
   local $_;
   for (@words)
   {
      # This list of shell metacharacters was taken from the Perl source
      # (do_exec(), in doio.c).  It is, in slightly more readable form:
      # 
      #    $ & * ( ) { } [ ] ' " ; \ | ? < > ~ ` \n
      #
      # (plus whitespace).  This totally screws up cperl-mode's idea of
      # the syntax, unfortunately, so don't expect indenting to work
      # at all in the rest of this function.

      if ($_ eq "" || /[\s\$\&\*\(\)\{\}\[\]\'\";\\\|\?<>~`\n]/)
      {
         # If the word has both " and ' in it, then just backslash all 
         #   metacharacters;
         # if it has just ' then encase it in "";
         # otherwise encase it in ''

         SUBST:
         {
            (s/([\s\$\&\*\(\)\{\}\[\]\'\";\\\|\?<>~`\n])/\\$1/g, last SUBST)
               if (/\"/) && (/\'/);
            ($_ = qq/"$_"/, last SUBST) if (/\'/);
            $_ = qq/'$_'/;
         }
      }
   }

   join (" ", @words);
}


1;
