# ------------------------------ MNI Header ----------------------------------
#@NAME       : file_utilities.pl
#@DESCRIPTION: Functions to manipulate/check/validate/search files and
#              directories.
#@GLOBALS    : 
#@CREATED    : 1995/08/07, Greg Ward (from gpw_utilities.pl)
#@MODIFIED   : 
#@VERSION    : $Id: file_utilities.pl,v 1.1.1.1 2000-01-19 14:10:31 louis Exp $
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

use POSIX qw(strftime);
use Config;                     # as soon as we upgrade to 5.003, can lose
                                # this and use $^O...
require 5.001;

# Aliases for the old way of naming these subs:

sub CheckOutputDirs { &check_output_dirs; }
sub CheckOutputPath { &check_output_path; }
sub CheckInputDirs { &check_input_dirs; }
sub CheckDirectories { &check_directories; }
sub SearchDirectories { &search_directories; }
sub CheckFiles { &check_files; }
sub FindProgram { &find_program; }
sub FindPrograms { &find_programs; }
sub Uncompress { &uncompress; }



# ------------------------------ MNI Header ----------------------------------
#@NAME       : check_output_dirs
#@INPUT      : @dirs - list of directories to check
#@OUTPUT     : 
#@RETURNS    : 1 if all directories exist and are writeable, *or* were
#                 successfully created
#              0 otherwise
#@DESCRIPTION: Checks that each directory in a list of directories either
#              exists and is writeable; *or* that it can be created.  
#              Prints a meaningful error message and returns 0 if any fail.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/05/02, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub check_output_dirs
{
   local (@dirs) = @_;
   local ($num_err) = 0;

   foreach $dir (@dirs)
   {
      next unless $dir;		# skip blank strings
      if (-e $dir || -l $dir)   # file exists (or is dangling link)
      {
	 if (! -d $dir)         # but *not* a directory
	 {
	    warn "$dir already exists but is not a directory\n";
	    $num_err++;
	 }
	 elsif (! -w $dir)      # is a directory, but not writeable
	 {
	    warn "$dir is a directory but is not writeable\n";
	    $num_err++;
	 }
      }
      else                      # no file, no dangling link
      {
	 if (! mkdir ($dir, 0755))
	 {
	    warn "Couldn't create \"$dir\": $!\n";
	    $num_err++;
	 }
      }
   }
   return ($num_err == 0);
}



# ------------------------------ MNI Header ----------------------------------
#@NAME       : &check_output_path
#@INPUT      : $path   - path to a file or directory
#@OUTPUT     : 
#@RETURNS    : 0 if any component of $path is not a directory (or
#                cannot be created as one), or if last directory component
#                isn't a writeable directory (or can't be created).
#              1 otherwise
#@DESCRIPTION: If $path is a file (ie. doesn't end with a slash), ensures
#              that conditions are optimal for creating it.  If $path
#              is a directory, ensures that it exists and is writeable.
#
#              In detail, splits a file path up and makes sure that every
#              directory component except the last is indeed a directory,
#              or can be created as one.  Makes sure that the last
#              component before a slash (last directory) is both a
#              directory and writeable; if not, attempts to create it.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1996/03/14, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub check_output_path
{
   my ($path) = @_;
   my (@dirs, $dir, $partial, $last_dir);

   # If $path starts with "/" (absolute path), then the first element of
   # @dirs will be blank, which we just strip off and ignore.  We also
   # strip off the last element of @dirs, unless $path ends in a "/", which
   # means it's a really a directory.  (Ie. if $path is "/foo/bar/baz", we
   # only check if "/foo" and "/foo/bar" are directories; but if it's
   # "/foo/bar/baz/", then we also check "/foo/bar/baz".)

   @dirs = split ("/", $path, -1);
   pop @dirs;                           # strip non-directory component

   if (@dirs)                           # still directories left in the list?
   {
      $last_dir = pop @dirs;            # get the last directory (it's
      if (@dirs && $dirs[0] eq "")      #   checked separately)
      {
         $partial = "/";
         shift @dirs;                   # remove first element (empty string)
      }
      else
      {
         $partial = "";
      }
   }
   else                                 # list empty after popping?
   {                                    # (must have been no slashes, hence
      $last_dir = ".";                  # relative to current dir)
      $partial = "";
   }

   # Check all but the last directory component: each must be a directory.

   for $dir (@dirs)
   {
#      my ($tmp) = "$partial$dir";

      $partial .= $dir;

      if (! (-e $partial || -l $partial))
      {
	 unless (mkdir ($partial, 0755))
	 {
	    warn "Couldn't create \"$partial\": $!\n";
	    return 0;
	 }
      }
      elsif (! -d $partial)
      {
	 warn "\"$path\" is not a writeable path because \"$partial\" is not a directory\n";
	 return 0;
      }
      $partial .= "/";
   }

   # Now check the last directory; must be a directory and writeable.

   $partial .= $last_dir;

   if (! (-e $partial || -l $partial))
   {
      unless (mkdir ($partial, 0755))
      {
	 warn "Couldn't create \"$partial\": $!\n";
	 return 0;
      }
   }
   elsif (! (-d $partial && -w $partial))
   {
      warn "\"$partial\" not a writeable directory\n";
      return 0;
   }

   return 1;
}
      



# ------------------------------ MNI Header ----------------------------------
#@NAME       : check_input_dirs
#@INPUT      : @dirs - list of directories to check
#@OUTPUT     : 
#@RETURNS    : 1 if all directories exist and are readable
#              0 otherwise
#@DESCRIPTION: Checks to see if all desired input directories exist and
#              are readable.  Prints meaningful error messages and returns
#              zero if anything is wrong.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/05/02, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub check_input_dirs
{
   my (@dirs) = @_;
   my ($num_err) = 0;

   foreach $dir (@dirs)
   {
      next unless $dir;		# skip blank strings
      if (-e $dir || -l $dir)
      {
	 if (! -d $dir)
	 {
	    warn "$dir exists but is not a directory\n";
	    $num_err++;
	 }
	 elsif (! -r $dir)
	 {
	    warn "$dir is a directory, but is not readable\n";
	    $num_err++;
	 }
      }
      else
      {
	 warn "directory $dir does not exist\n";
	 $num_err++;
      }
   }
   return ($num_err == 0);   
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : check_directories
#@INPUT      : $source - reference to list of source directories
#              $dest   - reference to list of destination directories
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Goes through the list of source directories, and makes sure
#              that each one exists; goes through list of destination
#              directories, and attempts to create those that don't exist.
#              Dies if any source directories are missing.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1994/08/22, Greg Ward
#@MODIFIED   : 1995/05/02, GW: changed to just call &check_input_dirs and
#                              &check_output_dirs
#-----------------------------------------------------------------------------
sub check_directories
{
   my ($source, $dest) = @_;
   my ($input_ok, $output_ok);

   $input_ok = &check_input_dirs (@$source);
   $output_ok = &check_output_dirs (@$dest);
   return ($input_ok && $output_ok);
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : search_directories
#@INPUT      : $file - name of file to search for
#              @dirs - list of directories to search in.  Note that
#                 the empty string should NOT be used; use "." for 
#                 the current directory.
#@OUTPUT     : 
#@RETURNS    : The directory where $file was found, or undef if it wasn't.
#@DESCRIPTION: Searches a list of directories for a single file.  Note
#              that only the existence test -e is used; doesn't check
#              for readability, executability, or whatever.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1994/09/16, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub search_directories
{
   my ($file, @dirs) = @_;
   my ($dir, $found);

   $found = 0;
   while ($dir = shift @dirs)
   {
      $dir .= "/" unless $dir =~ m|/$|;
      $found = (-e $dir . $file);
      last if $found;
   }

   return ($dir) if $found;
   return undef;
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : check_files
#@INPUT      : @files - a list of filenames to check for readability
#@OUTPUT     : 
#@RETURNS    : 1 if all files exist and are readable
#              0 otherwise (and prints a warning for every bad file)
#@DESCRIPTION: Makes sure that each of a list of files exists and is readable.
#              Existence test done with &test_file, so it will also check
#              for .gz, .pgp, etc. versions of the file (but won't tell
#              which one it found!)
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : Aug 1994, Greg Ward
#@MODIFIED   : 1995/04/25, GW: gave it two possible error messages
#              1995/10/30, GW: changed so it warns and returns 0 instead 
#                              of dying
#-----------------------------------------------------------------------------
sub check_files
{
   my (@files) = @_;
   my ($num_err, $new);

   $num_err = 0;
   foreach (@files)
   {
      $new = &test_file ('-e', $_);
      if (! $new)
      {	
         if (&test_file ('-l', $_))
         {
            warn "$_ is a dangling link (file does not exist)\n";
         }
         else
         {
            warn "$_ does not exist\n";
         }
	 $num_err++;
      }
      elsif (! -r $new)
      {
	 warn "$new not readable\n";
	 $num_err++;
      }
   }
   return ($num_err == 0);
}



# ------------------------------ MNI Header ----------------------------------
#@NAME       : &find_program
#@INPUT      : $program - name of program to search for
#              $path    - colon-delimited list of directories to look in
#                         (if not supplied, $ENV{"PATH"} will be used);
#                         OR a reference to a list of directories
#@OUTPUT     : 
#@RETURNS    : full path of program, or 0 if not found (also prints a 
#              warning if not found)
#@DESCRIPTION: Searches a list of directories for an executable program.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/07/20, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub find_program
{
   my ($program, $path) = @_;
   my (@path, $found, $fullpath);

   $path = $ENV{'PATH'} unless defined $path;

   if (! ref $path)
   {
      @path = split (":", $path);
   }
   elsif (ref $path eq 'ARRAY')
   {
      @path = @$path;
   }
   else 
   {
      warn "\&find_program: \$path must be either a scalar or ref to array\n";
      return undef;
   }

   $found = 0;
   foreach $dir (@path)
   {
      $fullpath = "${dir}/${program}";
      if (-x $fullpath && -f $fullpath)
      {
	 $found = $fullpath;
	 last;
      }
   }

   return $found if $found;
   warn "Couldn't find program \"$program\"\n";
   return 0;
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &find_programs
#@INPUT      : @$programs  - list of programs to search for
#              $path       - reference to list of directories to search
#                            (OR a colon-delimited string)
#@OUTPUT     : 
#@RETURNS    : list of found programs (as full pathnames)
#              OR empty list (if there were any missing programs)
#@DESCRIPTION: Searches for each program specified in the specified search
#              path.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : &find_program
#@CREATED    : 1995/12/08, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub find_programs
{
   my ($programs, $path) = @_;
   my ($missing, $found, @found);

   $missing = 0;
   foreach $prog (@$programs)
   {
      $found = &find_program ($prog, $path);
      $missing++ if ! $found;
      push (@found, $found);
   }
   return @found if $missing == 0;
   return ();
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &uncompress
#@INPUT      : $tmp_dir - directory to uncompress files to
#              @originals - list of files to consider for uncompression
#@OUTPUT     : 
#@RETURNS    : array context: @originals, changed so that any the names
#                of any files that were compressed are now decompressed
#                and in $tmp_dir
#              scalar context: first element of @originals, possibly changed
#                to decompressed version
#@DESCRIPTION: Uncompresses (if applicable) compressed files to $tmp_dir.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/07/31, Greg Ward
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub uncompress
{
   my ($tmp_dir, @originals) = @_;
   my ($uncompressed, @uncompressed, $orig);

   require "path_utilities.pl";	         # for &replace_dir

   foreach $orig (@originals)
   {
      if ($orig =~ /\.(Z|gz|z)$/)
      {
	 ($uncompressed = &replace_dir ($tmp_dir, $orig)) =~ s/\.(Z|gz|z)$//;
         unless (-e $uncompressed)
         {
            system ("gunzip -c $orig > $uncompressed");
            die "gunzip $orig failed\n" if ($?);
         }
	 push (@uncompressed, $uncompressed);
      }
      else
      {
	 push (@uncompressed, $orig);
      }
   }
   return @uncompressed if wantarray;
   return $uncompressed[0];
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &compress
#@INPUT      : $files - ref to list of files to compress
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Compresses (with gzip) each of a list of files.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1997/01/15, GPW
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub compress
{
   my ($files) = @_;
   my $file;

   foreach $file (@$files)
   {
      unless ($file =~ /\.(Z|gz|z)$/)
      {
         system "gzip", $file;
         die "gzip $file failed\n" if ($?);
      }
   }
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &test_file
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : name of the file that ultimately passes the test
#@DESCRIPTION: Performs a single file test (eg. "-e" or "-r") with
#              allowances made for possible compressed or encrypted files.
#              In particular, if the the filename doesn't appear to be
#              compressed already, then we test $file, then $file.gz,
#              $file.pgp, $file.z, and $file.Z.  Only if all of these tests
#              fail does the test as a whole fail.  Returns the name of the
#              file that ultimately passes the test, or 0 if none of
#              them do.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : Greg Ward, March 1996
#@MODIFIED   : 1997/03/06, added "pgp" extension
#-----------------------------------------------------------------------------
sub test_file
{
   my ($test, $file) = @_;

   # If the file as supplied passes the test, immediately return success
   
   return $file if (eval "$test '$file'");
   &Fatal ($@) if $@;

   # If the filename as supplied *is* compressed, then try doing the
   # test on the basename ($` is the portion of $file preceding the
   # match), and if successful return the matched filename.
   
   if ($file =~ /\.(Z|g?z|pgp)$/)
   {
      return $` if (eval "$test '$`'");
   }

   # If the filename as given is not compressed, then test the three
   # possible compressed formats -- return the matching filename as soon as
   # one of them matches.

   else
   {
      for $ending (qw(gz pgp z Z))
      {
	 return "$file.$ending" if (eval "$test '$file.$ending'");
      }
   }
   
   # If we make it to here, then all tests failed -- return 0

   return 0;
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &generate_numbered_filename
#@INPUT      : $base - first part of filename
#@OUTPUT     : $ext  - last part of filename (including "." if wanted)
#@RETURNS    : Empty string on error (conflicting filenames or error 
#              renaming); otherwise, next filename in numbered sequence
#              starting with $base.
#@DESCRIPTION: Generates a numbered filename by incrementing a counter
#              $i until ${base}_${i}${ext} is found not to exist.  If
#              $i is 1 -- i.e. there weren't any files named with $base
#              and $ext -- then "_$i" is left out of the filename entirely.
#
#              For example, if called with $base="foo" and $ext=".log", and
#              neither "foo.log" nor "foo_1.log" exist, returns "foo.log".
#              On the next call, "foo.log" will be renamed to "foo_1.log",
#              and "foo_2.log" is returned.  Subsequent calls return
#              "foo_3.log", "foo_4.log", etc.  If both "foo.log" and
#              "foo_1.log" exist, then we print a warning and return the
#              empty string -- you'll have to deal with this degenerate
#              situation yourself, because it should never arise if you
#              only use &generate_numbered_filename to generate filenames.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1996/08/01, GPW, from code in ICBM.pm
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub generate_numbered_filename
{
   my ($base, $ext, $add_date) = @_;
   my $i = 1;

   $base .= strftime ("_%Y-%m-%d", localtime (time))
      if $add_date;

   if (-e "${base}${ext}")
   {
      if (-e "${base}_1${ext}")
      {
	 warn "conflicting filenames: ${base}${ext} and ${base}_1${ext}";
         return '';
      }
      else
      {
	 rename ("${base}${ext}", "${base}_1${ext}") ||
            (warn ("unable to rename ${base}${ext}: $!\n"), return '');
      }
   }

   $i++ while (-e "${base}_${i}${ext}");
   ($i == 1) ? "${base}${ext}" : "${base}_${i}${ext}";
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &statfs
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : list of useful-sounding values that are common to IRIX 
#              and Linux:
#                 $type    - filesystem (see your header files for definitions)
#                 $bsize   - block size (in bytes)
#                 $blocks  - total number of blocks in filesystem
#                 $bfree   - free blocks in filesystem (under Linux, this
#                            will be the "available block" count -- ie.
#                            number of blocks available to non-superuser)
#                 $files   - number of file nodes
#                 $ffree   - number of available file nodes
#@DESCRIPTION: Attempts to call the statfs(2) system call.  Will only 
#              work on Linux/i86 or IRIX.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1997/03/06, GPW
#@MODIFIED   : 
#@COMMENTS   : there is a File::statfs (or maybe Filesystem::stat?) module
#              mentioned in the module list, but no sign of it on CPAN.  Hmmm.
#-----------------------------------------------------------------------------
sub statfs                      # this ought to be made more standard...
{
   my ($path) = @_;
   my ($buf, $r);
   my ($type, $pad, $bsize, $frsize, $blocks, $bfree, $bavail);
   my ($files, $ffree, $fsid, $namelen, $spare, $fname, $fpack);

   if ($Config{'archname'} =~ /^i[0-9]?86-linux$/)
   {
      # structure size taken from i86 Linux 2.0 man pages; only tested
      # on Linux 2.0/i86 

      require "syscall.ph";
      $buf = " " x ( (7 + 2 + 1 + 6) * 4);
      $r = syscall (&SYS_statfs, $path, $buf);
      ($type,$bsize,$blocks,$bfree,$bavail,$files,$ffree,$fsid,$namelen,$spare)
         = unpack ("lllllll2ll6", $buf);
      $bfree = $bavail;         # ignore the free/available distinction (RTFM)
   }
   elsif ($Config{'osname'} eq 'irix')
   {
      $len = (2 +2 + (6 * 4) + 6 + 6);
      $buf = " " x $len;
      $r = syscall (1035, $path, $buf, $len, 0);
      ($type,$pad,$bsize,$frsize,$blocks,$bfree,$files,$ffree,$fname,$fpack)
         = unpack ("ssllllllc6c6", $buf);
   }
   else
   {
      die "Sorry, I don't know how to do `statfs' on the $Config{'archname'} architecture\n";
   }

   if ($r == 0)                         # success?
   {
      return ($type, $bsize, $blocks, $bfree, $files, $ffree);
   }
   else                                 # failure
   {
      warn "statfs failed on \"$path\": $!\n";
      return;
   }
}


1;
