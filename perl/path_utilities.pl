# ------------------------------ MNI Header ----------------------------------
#@NAME       : path_utilities.pl
#@DESCRIPTION: Functions for recognizing, parsing, and tweaking Unix
#              filenames and paths:
#                 &replace_dir
#                 &replace_ext
#                 &split_path
#@GLOBALS    : 
#@CREATED    : 1995/08/07, Greg Ward (from gpw_utilities.pl)
#@MODIFIED   : 
#@VERSION    : $Id: path_utilities.pl,v 1.1.1.1 2000-01-19 14:10:32 louis Exp $
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


sub ReplaceDir { &replace_dir; }
sub ReplaceExt { &replace_ext; }
sub SplitPath { &split_path; }
sub MergePaths { &merge_paths; }


# ------------------------------ MNI Header ----------------------------------
#@NAME       : replace_dir
#@INPUT      : $newpath - directory to replace existing directories with
#              @files   - list of files to have directories replaced
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Replaces the directory component of a list of pathnames.
#              Returns the list of files with substitutions performed.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/05/04, Greg Ward
#@MODIFIED   : 1995/05/23, GW: renamed to &ReplaceDir
#-----------------------------------------------------------------------------
sub replace_dir
{
   local ($newpath, @files) = @_;

   $newpath = "." if $newpath eq "";
   $newpath =~ s|([^/])$|$1/|;	# add a trailing slash if necessary
   foreach $file (@files)
   {
      $file =~ s|.*/||;		# strip off any existing directory
      $file = $newpath . $file;
   }
   return @files if wantarray;
   return $files[0];
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : replace_ext
#@INPUT      : $newext  - extension to replace existing extensions with
#              @files   - list of files to have extensions replaced
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Replaces the final extension (whatever follows the final dot)
#              of a list of pathnames.  Returns the list of files with
#              substitutions performed in array context, or the first filename
#              from the list in a scalar context.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/05/23, Greg Ward (from &ReplaceDir)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub replace_ext
{
   local ($newext, @files) = @_;

   foreach $file (@files)
   {
      $file =~ s/\.[^\.]*$/\.$newext/; # replace existing extension
   }
   return @files if wantarray;
   return $files[0];
}
      


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &split_path
#@INPUT      : $path    - a Unix path specifiction (optional directory + file)
#              $ext_opt - specifies how to deal with file extensions:
#                         if "none", extension is ignored and returned as
#                           part of the base filename
#                         if "first", the *first* dot in a filename denotes
#                           the extension, eg ".mnc.gz" would be an extension
#                         if "last", the *last* dot denotes the extension,
#                           eg. just ".gz" would be the extension
#                         the default is "first"
#@OUTPUT     : 
#@RETURNS    : array: ($dir, $base, $ext)
#@DESCRIPTION: Splits a Unix path specification into directory, base file 
#              name, and extension.  (The extension is chosen based on
#              either the first or last dot in the filename, depending
#              on the $ext_opt argument; by default, it splits on the 
#              first dot in the filename.)
#
#              If there is no directory (ie. the `path' refers to a file
#              in the current directory), then $dir will be the
#              empty string.  Otherwise, $dir will be that portion
#              of $path up to and including the last slash.
#
#              If $ext_opt is 'first' (the default), then $ext will be
#              that portion of $path from the first period after the 
#              last slash, inclusive.  If $ext_opt is 'last', then
#              $ext will go from the last period after the last slash.
#              If $ext_opt is 'none', then $ext will be undefined and
#              whatever extensions are in $path will be rolled into $base.
#              (If there are no extensions in $path, then $ext will be
#              undefined.)
#
#              $base is just whatever portion of $path is leftover after
#              pulling off $dir and $ext -- eg. from the last slash to
#              the first period (if $ext_opt is 'first'), or from the
#              last slash to the last period (if $ext_opt is 'last').
#
#              if $ext_opt is 'first':
#                  "/foo/bar/zap.mnc.gz" --> ("/foo/bar/", "zap", ".mnc.gz")
#              if $ext_opt is 'last':
#                  "/foo/bar/zap.mnc.gz" --> ("/foo/bar/", "zap.mnc", ".gz")
#              if $ext_opt is 'none':
#                  "/foo/bar/zap.mnc.gz" --> ("/foo/bar/", "zap.mnc.gz")
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/05/10, Greg Ward - taken from mritotal and modified
#@MODIFIED   : 1995/08/10, GW: added $ext_opt option to handle splitting off
#                              the extension in different ways
#              1997/02/26, GW: changed to preserve trailing slash and 
#                              empty directory string
#-----------------------------------------------------------------------------
sub split_path
{
    my ($path, $ext_opt) = @_;
    my ($dir, $base, $ext);

    $ext_opt = "first" unless defined $ext_opt;

    # If filename has no extension, don't try to act as though it does
    # (both "last" and "first" options assume there is an extension

#    $ext_opt = "none" if $path !~ m+/?[^/]*\.+;

    if ($ext_opt eq "none")
    {
       ($dir, $base) = $path =~ m+^(.*/)?([^/]*)$+;
    } 
    elsif ($ext_opt eq "first")
    {
       ($dir, $base, $ext) = $path =~ m+^(.*/)?([^/\.]*)(\..*)?$+;
    }
    elsif ($ext_opt eq "last")
    {
       ($dir, $base, $ext) = $path =~ m+^(.*/)?([^/]*)(\.[^\/]*)$+;
    }
    else
    {
       die "split_path: unknown extension option \"$ext_opt\"\n";
    }

    $dir = "" unless ($dir);

    ($dir, $base, $ext);
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &merge_paths
#@INPUT      : $path      - ref to a list of directories
#              $additions - ref to a list of directories to add to $path
#@OUTPUT     : @$path     - updated so elements from @$additions are added
#                           only if they weren't already in @$path
#@RETURNS    : 
#@DESCRIPTION: 
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/12/04 GPW
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub merge_paths
{
   my (%seen, $dir, @path);

   foreach $dir (@_)
   {
      $dir =~ s|/$||;                   # strip trailing slash
      push (@path, $dir) unless $seen{$dir};
      $seen{$dir} = 1;
   }
   @path;
}

1;
