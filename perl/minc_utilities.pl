# ------------------------------ MNI Header ----------------------------------
#@NAME       : minc_utilities.pl
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Interfaces to a various MINC utilities for common
#              image- and volume-processing tasks.  Uses &Spawn to run
#              the commands, so be sure you have set the error action
#              to something reasonable (using &SetSpawnOptions).
#              Also, requires the following global variables to give
#              the full location and standard options for the
#              programs:
#                 $MincInfo
#                 $VolumeStats
#                 $VolumeCOG
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : Greg Ward, July 95
#@MODIFIED   : considerably
#@VERSION    : $Id: minc_utilities.pl,v 1.1.1.1 2000-01-19 14:10:31 louis Exp $
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

require 5.001;
use JobControl;
require "numeric_utilities.pl";

sub volume_min
{
   my ($volume) = @_;
   my ($status, $volmin);
   
   ($status, $volmin) = 
      &Spawn ("mincinfo -varvalue image-min $volume | sort -n | head -1");
   if ($Execute) { chop $volmin; return $volmin; }
   0;
}


sub volume_max
{
   my ($volume) = @_;
   my ($status, $volmax);
   
   ($status, $volmax) = 
      &Spawn ("mincinfo -varvalue image-max $volume | sort -n | tail -1");

   if ($Execute) { chop $volmax; return $volmax; }
   0;
}


sub volume_minmax
{
   my ($volume) = @_;
   my ($status, $volmax, $volmin);
#   die "\$MincInfo undefined" unless defined $MincInfo;

   (&volume_min ($volume), &volume_max ($volume));
}


sub volume_threshold
{
   my ($volmin, $volmax, $percent) = @_;
   ($volmax - $volmin) * $percent + $volmin;
}


sub auto_threshold
{
   my ($volume) = @_;
#   die "\$VolumeStats undefined" unless defined $VolumeStats;

   ($status, $threshold) = 
      &Spawn ("volume_stats -biModalT -quiet $volume");

   if ($Execute) { chop $threshold; return $threshold; }
   0;
}


sub volume_cog
{
   my ($volume) = @_;
   my ($cog, @cog);
#   die "\$VolumeCOG undefined" unless defined $VolumeCOG;

   ($status, $cog) = &Spawn ("volume_cog $volume | tail -1");
   if ($Execute) { chop $cog; return $cog; }
   "0 0 0";
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : volume_params
#@INPUT      : $volume - file to query
#              $start  - reference to list to fill with starts
#              $step   - reference to list to fill with step sizes
#              $length - reference to list to fill with dimension lengths
#              $dircos - reference to list to fill with direction cosines
#              $dims   - reference to list to fill with dimension names
#              (If any of these "references" are just false values, then
#              they won't be followed -- that way you can specify "undef"
#              for values you're not interested in.)
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Gets the x, y, and z starts, steps, and lengths for a
#              MINC volume.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 95/04/21, Greg Ward
#@MODIFIED   : 95/08/04, GW: simplified to handle negative steps and extents
#              95/08/16, GW: changed to GetVolumeParams and made so it returns
#                            all three interesting spatial parameters
#              95/08/30, GW: renamed to volume_params and moved from autocrop
#                            to minc_utilities.pl
#              96/02/12, GW: converted to use Perl 5 references
#                            added direction cosines stuff
#              97/04/07, GW: added $dims to get dimension names
#-----------------------------------------------------------------------------
sub volume_params
{
   my ($volume, $start, $step, $length, $dircos, $dims) = @_;
   my ($cmd, $output, @output);

   $cmd = "mincinfo $volume -error_string 0 " .
          "-attval xspace:start -attval yspace:start -attval zspace:start " .
          "-attval xspace:step -attval yspace:step -attval zspace:step " .
	  "-dimlength xspace -dimlength yspace -dimlength zspace " .
	  "-attval xspace:direction_cosines " .
	  "-attval yspace:direction_cosines " .
	  "-attval zspace:direction_cosines " .
          "-dimnames";
   $output = `$cmd`;
   die "Error executing mincinfo on $volume\n" if $?;
   @output = split ("\n", $output);
   grep (s/(^\s+)|(\s+$)//g, @output); # strip leading and trailing whitespace
         
   @$start = @output[0..2] if $start;
   @$step = @output[3..5] if $step;
   @$length = @output[6..8] if $length;
   if ($dircos)
   {
      @$dircos = (1,0,0, 0,1,0, 0,0,1);
      for $i (0,1,2)
      {
	 my ($base) = $i * 3;

	 @$dircos[$base .. $base+2] = split (/\s+/, $output[9+$i])
	    unless ($output[9+$i] eq "0");
      }
   }
   @$dims = split ("", $output[12]) if $dims;         
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &get_history
#@INPUT      : $volume
#@OUTPUT     : 
#@RETURNS    : list containing all the elements from the MINC file
#@DESCRIPTION: Fetches the global "history" attribute from a MINC file,
#              splits it on newlines, and returns the list resulting
#              from that split.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : Nov 1995, Greg Ward (originally in get_flipped_volume)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub get_history
{
   my ($volume) = @_;
   my ($history, @history);

   $history = `mincinfo -error_string "" -attvalue :history $volume`;
   die "Error executing mincinfo on $volume\n" if $?;
   @history = split ("\n", $history);
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &put_history
#@INPUT      : $volume
#              @history
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Joins the @history array with newlines, and puts the 
#              resulting string into a MINC file as the global history
#              attribute.  (This completely replaces the existing
#              history; if you wish to update it, see &update_history).
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : Nov 1995, Greg Ward (originally in get_flipped_volume)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub put_history
{
   my ($volume, @history) = @_;
   my ($history);

   $history = join ("\n", @history) . "\n";
   system ("minc_modify_header", $volume, "-sinsert", ":history=$history");
   die "Error modifying header of $volume\n" if $?;
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &update_history
#@INPUT      : $volume  - name of MINC file to update
#              $replace - number of items to delete from the history list
#                 before adding the new one.  If not given or zero,
#                 no items are removed.
#              $history - entry to add to history; can take one of three
#                 possible forms:
#                   - non-empty string: $history is copied into the MINC
#                     file's history attribute with no changes or additions
#                   - array ref (assumed to be, eg. [$0 @ARGV]): 
#                     a "userstamp" (eg. "[user@host:/foo/bar/baz]") and 
#                     timestamp (eg. "[96-05-29 12:03:36]) are generated,
#                     and printed along with @$history
#                   - empty or undefined (ie. a "false" value): 
#                     same as if you pass an array ref [$0 @ARGV] (so
#                     passing an array ref is really if you want something
#                     *other* than [$0 @ARGV], ie. if -- for whatever
#                     nefarious reason -- you wish to lie about your 
#                     program's name and arguments
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Adds an item to the global history attribute in a MINC file.
#              Optionally deletes other items from the end of the list
#              first, allowing you to (eg., if $replace == 1) replace
#              the last item with your desired item.  This is useful for
#              Perl scripts that are front-ends for a known number of
#              operations by standard MINC tools, each of which results
#              in a single history item.  
#
#              Note: if a timestamp is generated, it will be the time
#              at which minc_utilities.pl was compiled -- presumably,
#              when your script starts -- rather than the current time.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : &get_history, &put_history
#              &userstamp, &timestamp
#@CREATED    : 1996/05/29, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub update_history
{
   require "misc_utilities.pl";               # for &userstamp and &timestamp
   my ($volume, $replace, $history) = @_;

   # First figure out the real history line: either the caller supplied the
   # whole thing ($history is a non-empty string), or they supplied a list
   # with the program name and arguments ($history is an array ref), or
   # they want us to cook it up entirely ($history is a false value)

   if (defined $history && $history ne "")
   {
      # do nothing -- put $history right into the MINC file
   }
   elsif (defined $history && ref $history eq "ARRAY")
   {
      $history = sprintf ("%s %s ", &userstamp(), &timestamp($_start_time)) . 
                 join (" ", @$history);
   }
   elsif (! $history)
   {
      $history = sprintf ("%s %s %s ", 
                          &userstamp(), &timestamp($_start_time), $0) . 
                 join (" ", @ARGV);
   }
   else
   {
      die "don't know what to do with the supplied \$history value";
   }

   # Now fetch the existing history attribute

   my @orig_history = &get_history ($volume);

   # Remove the last $replace items (if $replace was supplied)

   splice (@orig_history, -$replace) if (defined $replace && $replace > 0);

   # Put $history onto the list, and write it to the MINC file

   push (@orig_history, $history);
   &put_history ($volume, @orig_history);
}



# ------------------------------ MNI Header ----------------------------------
#@NAME       : &get_dimension_order
#@INPUT      : $volume - name of MINC file to get dimension names from; 
#                   OR - reference to an array containing the dim names
#@OUTPUT     : 
#@RETURNS    : $order  - ref to dimension order list
#              $perm   - ref to dimension permutation list
#@DESCRIPTION: Computes the dimension order and permutation for a MINC
#              file.  These are two vectors that are very useful when
#              you need to go back and forth between the canonical
#              dimension ordering (x,y,z) and whatever order the dimensions
#              happen to be in in a particular MINC file.
#
#              The dimension order vector is the easy one: order[i] tells
#              you which dimension is the i'th dimension of your volume.
#              For instance, a coronal volume has dimensions (y,z,x); it's
#              order vector is (1,2,0), a simple transcription of (y,z,x)
#              to numerical form.  (Put another way, order[0]==1 means that
#              dimension 0 of the file is canonical dimension 1, or
#              yspace.)
#
#              The permutation vector is a little trickier to wrap your
#              head around, even though it's just the inverse of the order
#              vector.  In short, perm[i] is where to find the i'th
#              canonical dimension in your file's dimension list.  Going
#              with the coronal example again, the permutation vector is
#              (2,0,1): looking up canonical dimension 2 (zspace) in perm[]
#              gives 1, and whaddya know? zspace is at slot 1 in the list
#              of dimensions.
#
#              The main reason that these two are so confusing is that
#              they're usually the same -- the reason I've used the coronal
#              ordering as an example here is that it's the only standard
#              ordering where the order and permutation vectors are
#              different!  (Of the 6 possible orders for three dimensions,
#              only coronal (y,z,x) and the non-standard order (z,x,y) have
#              different order and permutation vectors.)  However, to be
#              truly general, you have to know when to use which one.
#
#              In short: use the order vector when you have something in
#              (x,y,z) order and want it in volume order; use the
#              permutation vector to go from volume to (x,y,z) order.  This
#              is particular easy in Perl using array slices.  Say you have
#              a list of parameters in (x,y,z) order:
#
#                   @step = (x_step, y_step, z_step)
#
#              that you want in volume order.  Again assuming a coronal 
#              volume, the order vector is (1,2,0), and so
#
#                   @step_v = @step[@order]
#                           = @step[1,2,0] 
#                           = (y_step, z_step, x_step)
#
#              which of course is in coronal order.  Neat-o!
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : mincinfo
#@CREATED    : 1996/10/22, GW (from code formerly in autocrop)
#@MODIFIED   : 
#@COMMENTS   : The "order" and "permutation" jargon is entirely my
#                 own invention; I don't know if anybody else uses
#                 the same terms.  Helps me get a grip on this damn 
#                 stuff, at any rate.
#              Shouldn't actually bomb on volumes with < 3 spatial 
#                 dimensions (or with non-spatial dimensions; they will
#                 just be ignored).  However, I really don't know if
#                 it produces useful results in those cases.
#-----------------------------------------------------------------------------
sub get_dimension_order
{
   my ($volume) = @_;
   my (@dimlist, %dim_num, @order, @perm);

   %dim_num = ('xspace', 0, 'yspace', 1, 'zspace', 2);

   if (! ref $volume)           # it's a scalar -- name of MINC file
   {
      @dimlist = split (/\s+/, `mincinfo -dimnames $volume`);
   }
   elsif (ref $volume eq 'ARRAY')
   {
      @dimlist = @$volume;
   }
   else
   {
      die "get_dimension_order: \$volume must be a scalar or array ref";
   }

   @dimlist = grep (/^[xyz]space$/, @dimlist);

   for $i (0 .. $#dimlist)
   {
      $dim_num = $dim_num{$dimlist[$i]};
      $order[$i] = $dim_num;
      $perm[$dim_num] = $i;
   }

   (\@order, \@perm);
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : compute_resample_params
#@INPUT      : @$start  - start (world) (x,y,z) coordinates for resampling
#              @$extent - extent (length in mm) of each dimension (x,y,z)
#              @$step   - step (in mm) of each dimension (x,y,z)
#@OUTPUT     : 
#@RETURNS    : $params - the mincresample command line (except for 
#                        filenames) needed to carry out the bounding
#                        specified by @start and @extent, with steps 
#                        specified by @step
#@DESCRIPTION: Computes the needed dimension lengths to accomodate the
#              desired step and extent, and generates appropriate
#              parameters for mincresample (-start, -step, -nelements)
#              to carry out the volume bounding.
#
#              Does nothing about direction cosines, -like, or
#              -use_input_sampling -- that's up to the caller to figure out.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 95/04/21, Greg Ward
#@MODIFIED   : 96/10/21, GW: changed to use Perl5 ref's and my
#              96/10/22, GW: changed name (from StudlyCaps to lower_case)
#                            and moved from autocrop into minc_utilities.pl
#@COMMENTS   : What should we do about direction cosines?!?!?!?
#-----------------------------------------------------------------------------
sub compute_resample_params
{
   my ($start, $extent, $step) = @_;
   my (@length, $params);

   foreach $i (0,1,2)
   {
      $length[$i] = round ($extent->[$i] / $step->[$i], 1, +1);
   }

   $params = sprintf ("-start %g %g %g -step %g %g %g -nelements %d %d %d",
		      @$start, @$step, @length);
}  


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &compute_reshape_params
#@INPUT      : @$order    - the dimension order vector, as returned by
#                           &get_dimension_order
#              @$oldstart - the current `start' parameters of the volume
#                           to be reshaped
#              @$oldstep  - the current `step' parameters of the volume
#                           to be reshaped
#              @$start    - the desired new `start' parameters (after
#                           reshaping)
#              @$extent   - the desired new extent of each dimension in mm
#              @$step     - the desired new step (may differ from @$oldstep
#                           in sign only!)
#@OUTPUT     : 
#@RETURNS    : $params    - parameters for mincreshape, in 
#                           "-start ... -count ..." form
#@DESCRIPTION: Computes the parameters necessary for mincreshape to
#              give a volume a new spatial extent as described by @$start
#              and @$extent.  If @$oldstart and @$start differ by anything
#              other than integral multiples of @$step (or @$oldstep),
#              then it will only approximate the desired bounds.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : &round (from numeric_utilities.pl)
#@CREATED    : 95/06/25, Greg Ward
#@MODIFIED   : 95/08/18, GW: fixed so it would work
#              95/10/01, GW: fixed finding of vstart 
#              96/10/21, GW: changed to use Perl5 ref's and my
#              96/10/22, GW: changed name (from StudlyCaps to lower_case)
#                            and moved from autocrop into minc_utilities.pl;
#                            changed to take the dimension order as a 
#                            parameter rather than computing it here
#-----------------------------------------------------------------------------
sub compute_reshape_params
{
   my ($order, $oldstart, $oldstep, $start, $extent, $step) = @_;
   my (@perm, @vstart, @count, $params);
   my $i;

   foreach $i (0 .. 2)
   {
      $vstart[$i] = &round (($start->[$i] - $oldstart->[$i]) / $oldstep->[$i]);
      $count[$i] = &round ($extent->[$i] / $step->[$i], 1, +1);
      $count[$i] = -$count[$i] if ($step->[$i] == -$oldstep->[$i]);
   }
   @vstart = @vstart[@$order];
   @count = @count[@$order];

   $params = sprintf ("-start %d,%d,%d -count %g,%g,%g", @vstart, @count);
}


$_start_time = time;            # when script starts running, so 
                                # &update_history will know


1;
