# ------------------------------ MNI Header ----------------------------------
#@NAME       : volume_heuristics.pl
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Routines that make some guesses about the nature of volumetric
#              data sets of the human brain, generally in order to reduce the
#              amount of data in a consistent way.  You should 
#              `require "numeric_utilities.pl"' before calling any of 
#              the routines in this file.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 95/10/27, Greg Ward (from code in mritotal)
#@MODIFIED   : 
#-----------------------------------------------------------------------------

require 5.001;

# ------------------------------ MNI Header ----------------------------------
#@NAME       : &GuessSubsample
#@INPUT      : @steps - the step sizes of the three spatial dimensions
#                       for the input volume
#@OUTPUT     : 
#@RETURNS    : new step sizes to use
#@DESCRIPTION: 
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 95/08/30, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub GuessSubsample
{
   my ($steps) = @_;
   my (@newsteps);

   # Default subsample step values are just 2*step for any dimension
   # with steps < 1.5mm.  This isn't exactly rocket science, and it
   # doesn't have to be since subsampling is done solely to save time
   # and memory.

   @newsteps = @$steps;
   for $i (0,1,2)
   {
      $newsteps[$i] = 2*$steps->[$i] if (abs ($steps->[$i]) < 1.5);
   }
   @newsteps
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &GuessCrop
#@INPUT      : $step, $length - references to arrays containing
#                 the various spatial parameters for the input volume
#@OUTPUT     : 
#@RETURNS    : array of crop specifications to use on the volume
#@DESCRIPTION: Determines how much of a volume we need to remove in order
#              to preserve the top 190mm of the volume.  This is useful
#              to reduce the amount of data in a large sagittal (or coronal)
#              acquisition, as I've seen anecdotally that this seems to
#              be plenty to get in the whole brain, as long as the top of the
#              head is near the top of the volume.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 95/08/30, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub GuessCrop
{
   my ($step, $length) = @_;
   my (@extent, @crop);

   for $i (0,1,2)
   {
      $extent[$i] = $step->[$i] * $length->[$i];
   }
   @extent = &labs (@extent);
   
   # For now, we only crop in the z dimension by default.  The biggest
   # assumption here is that the top of the head is pretty near the
   # top of the scanning volume.  The heuristic used is this: if the
   # z-extent is > 190 mm, we chop off anything outside of the top 190
   # mm.  The reasoning for this is quite arbitrary -- ie. it looks
   # like that's a pretty good figure to use for ICBM data acquired at
   # the MNI, so presumably it applies to the entire world.

   @crop = ("0,0", "0,0", "0,0");
   $crop[2] = sprintf ("%g,0", 190-$extent[2])
      if (abs ($extent[2]) > 190);
   @crop;
}
