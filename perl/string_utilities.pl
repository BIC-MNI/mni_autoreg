# ------------------------------ MNI Header ----------------------------------
#@NAME       : string_utilities.pl
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: String utilties (duh!).
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1996/06/11, Greg Ward
#@MODIFIED   : 
#@VERSION    : $Id: string_utilities.pl,v 1.1 2000-02-28 16:38:25 stever Exp $
#-----------------------------------------------------------------------------


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &parse_num_list
#@INPUT      : $list - a string containing a "list" of numbers, something
#                      like "1-3,4,6,9-11"
#@OUTPUT     : 
#@RETURNS    : a real Perl list with the string expanded as you'd expect
#              (e.g. (1,2,3,4,6,9,10,11))
#              OR empty list if your string is bogus (and prints a warning)
#@DESCRIPTION: 
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : January 1996, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub parse_num_list 
{
   my ($list) = @_;
   my (@chunks, @list);

   $list =~ s/\s//;
   @chunks = split (/,/, $list);
   for $chunk (@chunks)
   {
      if ($chunk =~ /^\d+$/)
      {
	 push (@list, $chunk);
      }
      elsif ($chunk =~ /^(\d+)\-(\d+)$/)
      {
	 push (@list, $1 .. $2);
      }
      else
      {
	 warn "Numeric list $list is illegal";
	 return ();
      }
   }
   @list;
}


1;
