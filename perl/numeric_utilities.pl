# ------------------------------ MNI Header ----------------------------------
#@NAME       : numeric_utilities.pl
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Short routines for doing common numeric tasks, including
#              &in_range, &labs (list absolute value), &round, ...
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 95/08/07, Greg Ward - broke off from gpw_utilties.pl
#@MODIFIED   : 
#@VERSION    : $Id: numeric_utilities.pl,v 1.1 2000-01-19 14:10:32 louis Exp $
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

use POSIX qw(floor ceil);


# ------------------------------ MNI Header ----------------------------------
#@NAME       : in_range
#@INPUT      : $val - value to test
#              $lo  - lower bound
#              $hi  - upper bound
#@OUTPUT     : 
#@RETURNS    : 0 if $val is within the closed interval [$lo, $hi]
#              -1 if $val is less than $lo
#              +1 if $val is greater than $hi
#@DESCRIPTION: Tests whether a number is within a specified range.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 95/3/8, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub in_range
{
   my ($val, $lo, $hi) = @_;

   return (-1) if ($val < $lo);
   return (+1) if ($val > $hi);
   return (0);
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &labs
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Returns the absolute value of its argument(s).  In an array
#              context, returns a list of absolute values, otherwise just
#              a scalar.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 95/04/12, GW
#@MODIFIED   : 95/09/01, GW: fixed so it doesn't modify its arguments (dohh!)
#              97/04/11, GW: renamed to labs (list absolute value), and
#                            made it call Perl's builtin abs function
#@COMMENTS   : You might think this was obsolete with Perl 5, but
#              the Perl 5 builtin only takes a scalar -- this version
#              will take and return a list if desired.
#-----------------------------------------------------------------------------
sub labs
{
   my (@nums) = @_;
   foreach $elt (@nums)
   {
      $elt = abs ($elt);
   }
   wantarray ? @nums : $nums[0];
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &round
#@INPUT      : $value - number to round
#              $factor - what to round to; eg. 1 to round to an integer,
#                        .5 to round to a half-integer, or 10 to a
#                        factor of 10.
#              $direction - 0 to round to nearest $factor
#                          -1 to round down to next $factor
#                          +1 to round up to next $factor
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Round off a number to a factor of
#              $factor, where $factor is the first argument.  Round either
#              to the nearest factor, or down, or up, depending on the value
#              of $direction.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 95/05/03, Greg Ward (hacked from Alex Zijdenbos' code)
#@MODIFIED   : 95/08/07, GW: added $direction
#              96/10/21, GW: changed "0" to "0.0" in comparison (!?!??!!!?)
#-----------------------------------------------------------------------------
sub round
{
   my($value, $factor, $direction) = @_;
   $factor = 1 unless defined $factor;
   $direction = 0 unless defined $direction;

   $factor = abs ($factor);
   $value /= $factor;
   if ($direction == 0)
   {
      $value += ($value < 0.0) ? (-0.5) : (+0.5);
      $value = int($value) * $factor;
   }
   elsif ($direction == -1)
   {
      $value = floor ($value) * $factor;
   }
   elsif ($direction == +1)
   {
      $value = ceil ($value) * $factor;
   }
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &max
#@INPUT      : List of scalars (numeric).
#@OUTPUT     : 
#@RETURNS    : Maximum value from list.
#@DESCRIPTION: Finds the max value in a list.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : August 95, GW (from Perl man page)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub max
{
   my($max) = pop(@_);
   foreach $foo (@_) {
      $max = $foo if $max < $foo;
   }
   $max;
}

# ------------------------------ MNI Header ----------------------------------
#@NAME       : &min
#@INPUT      : List of scalars (numeric).
#@OUTPUT     : 
#@RETURNS    : Minimum value from list.
#@DESCRIPTION: Finds the min value in a list.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : August 95, GW (from Perl man page)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub min
{
   my($min) = pop(@_);
   foreach $foo (@_) {
      $min = $foo if $min > $foo;
   }
   $min;
}


1;
