# ------------------------------ MNI Header ----------------------------------
#@NAME       : JobControl.pm
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: A module for doing process and job stuff -- spawning 
#              children and waiting for them (or detaching them to the
#              background), detaching ourself to the background, running
#              stuff remotely -- that sort of thing.
#@METHOD     : 
#@GLOBALS    : %Options    - package options:
#                  Verbose     print commands as we execute them?
#		   Execute     actually execute commands?
#                  MergeErrors merge stderr with stdout?
#                  ErrorAction what to do when a command fails
#                  Batch       submit commands to batch system?
#                  Clobber     overwrite output files (instead of append)
#                  LogHandle   filehandle to write commands to
#                  Notify      should &Obituary actually send mail?
#              %Programs   - known programs and their full paths
#              %PreOptions - options to come at the start of the command line
#                            for a given program
#              %PostOptions- options for the end of the command line
#@CALLS      : 
#@CREATED    : 1995/03/03, Greg Ward
#@MODIFIED   : 1995/11/14, GW: made into a package
#              1996/04/11, GW: made into a Perl 5 module
#              1996/12/10, GW: thorough overhaul to allow commands as lists,
#                              default options for commands, and a priori
#                              searching for required programs; also major
#                              code reorg and cleanup
#@VERSION    : $Id: JobControl.pm,v 1.1 2000-01-19 14:10:31 louis Exp $
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
require "misc_utilities.pl";    # for &userstamp and &timestamp
require "path_utilities.pl";    # for &split_path
require "file_utilities.pl";    # for &find_programs

package JobControl;

use strict;
use vars qw(@ISA @EXPORT @EXPORT_OK);
use Carp;
use FileHandle;			# for autoflush method
use Cwd;

require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(SetSpawnOptions Backgroundify Spawn Execute Obituary);
@EXPORT_OK = qw(SetProgramOptions AddProgramOptions);

#
# Some terminology that might help you (or me) in deciphering this code:
#
#  a `command' is either a string or a list containing both a program
#    to execute and some or all of the arguments it expects
#  a `word' is either one element of a command (if it's a list) or
#    one whitespace-delimited chunk of characters in the shell sense
#    (if it's a string)
#  a `program' is just the first word of a command
#
# It is important to keep in mind that commands can be lists (well, really
# list refs, for easy passing around).  This has the important advantage of
# guaranteeing that no pesky shells will be involved in executing your
# command (unless, err, batch gets in the way... but that's another
# story...)
#


# ----------------------------------------------------------------------
# Package "globals"

# Here's something deeply weird: the way these variables are declared
# now (ie. as globals, accessible as eg. %JobControl::Options or,
# thanks to the `use vars' here, just as %Options within this package),
# everything works fine.
#
# BUT -- if I try to lexically scope them so they are visible only
# within this file (using the commented-out `my'), then weird things
# start to happen.  In particular, the variables disappear between 
# subroutine calls!  For instance, if I add something to %Programs in
# &FindPrograms, then %Programs maintains that value *only* within
# the loop (ie. scope) where it is set; after that, it disappears.
# In fact, the value of \%Programs changes -- ie. the variable is 
# picking up and moving on elsewhere!  And, of course, the contents
# of the hash are being lost.
#
# Anyone who can explain this to me wins a free muffin at Cafe Neuro.

use vars qw(%Programs %PreOptions %PostOptions %Options);
#my (%Programs, %PreOptions, %PostOptions, %Options);

# $Notify = $ENV{'USER'};
# $ErrorAction = "";
# $MergeErrors = 1;
# $Batch = 0;
# $Clobber = 0;
# $LogHandle = "STDOUT";

# Package options (not to be confused with program options!)
%Options =
   ("Verbose"     => undef,             # print commands as we execute them?
    "Execute"     => undef,             # actually execute commands?
    "Strict"      => 0,                 # complain about unknown programs?
    "SearchPath"  => undef,             # list of directories to search
    "MergeErrors" => 1,                 # merge stderr with stdout?
    "ErrorAction" => "",                # what to do when a command fails
    "Batch"       => 0,                 # submit commands to batch system?
    "Clobber"     => 0,                 # overwrite output files (not append)
    "LogHandle"   => "STDOUT",          # filehandle to write commands to
    "Notify"      => $ENV{'USER'});     # should &Obituary actually send mail?


# Here's where we keep track of what programs we know about and the
# standard options to run 'em with

%Programs = ();                 # hash mapping program name -> full path
%PreOptions = ();               # program name -> list of default options
%PostOptions = ();              # program name -> list of default options
                                #   that come *after* all other arguments



# ----------------------------------------------------------------------
# Option-setting routines:
#   &SetOptions
#   &SetSpawnOptions [obsolete!]
#   &FindPrograms
#   &SetProgramOptions
#   &AddProgramOptions

# ------------------------------ MNI Header ----------------------------------
#@NAME       : SetOptions
#@INPUT      : (option,value) pairs -- as many as you like
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Set various package options.  Carps if you call it in an
#              obsolete way, croaks if you call it with a bad option or 
#              an odd number of arguments.
#@METHOD     : 
#@GLOBALS    : %Options
#@CALLS      : 
#@CREATED    : GPW, 1996/03/27 (from SetSpawnOptions)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub SetOptions
{
   my (@optval) = @_;
   my ($opt, $val);

#   ($filename,$line) = (caller)[1..2]; # for error messages

   if (@optval == 2 && $optval[0] =~ /\d+/)
   {
      carp ("warning: obsolete version of \&SetSpawnOptions called");
      @Options{'MergeErrors','ErrorAction'} = @optval;
      return;
   }

   while (@optval)
   {
      $opt = shift @optval;
      $val = shift @optval;
      if ($opt eq "ClobberLogs")
      {
         $opt = "Clobber";
         carp ("ClobberLogs option is obsolete -- use Clobber");
      }

      croak ("Unmatched option $opt")
	 unless defined $val;

      croak ("Unknown option $opt") unless
	 exists $Options{$opt};

      $Options{$opt} = $val;
   }
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &SetSpawnOptions
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Obsolete version of SetOptions.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : GPW, Summer 1995
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub SetSpawnOptions
{
   carp "Use of \&SetSpawnOptions is deprecated -- " .
        "please use \&JobControl::SetOptions";
   &SetOptions;
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &FindPrograms
#@INPUT      : $programs - [array ref] list of programs to search for;
#              $path     - [optional; array ref or colon-separated string] 
#                          list of directories to search
#@OUTPUT     : none
#@RETURNS    : true if all programs found; false otherwise
#@DESCRIPTION: Searches for programs that you will wish to run at some 
#              point in your program (eg. via &Spawn or &Execute).
#              The full path to each program is stored in a private 
#              structure for future reference by this package; if you
#              need access to these paths, perhaps you should be using
#              &find_programs from file_utilities.pl -- that's all
#              this routine does!
#
#              Do not confuse this with &find_program or &find_programs 
#              (both of which are in file_utilities.pl, and should wind
#              up in the `main' package if you `use JobControl').
#              They are all distinct routines; this &FindPrograms
#              calls &find_programs once, which calls &find_program
#              repeatedly.
#@METHOD     : 
#@GLOBALS    : %Programs
#@CALLS      : 
#@CREATED    : GPW, November 1996
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub FindPrograms
{
   my ($programs, $path) = @_;

   croak "programs must be an array reference"
      unless ref $programs eq 'ARRAY';
   croak "path must be a scalar or an array reference"
      unless (ref $path eq 'ARRAY' || ! ref $path);
   my @fullpaths = &::find_programs ($programs, $path || $Options{'SearchPath'});
   if (@fullpaths)              # all found successfully?
   {
      confess "Wrong number of full paths to go with program list"
         unless (@fullpaths == @$programs);

      map { $Programs{$_} = shift @fullpaths } @$programs;
      return 1;
   }
   else
   {
      return 0;
   }
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &SetProgramOptions
#@INPUT      : $program - the name of the program to set options for
#                         (doesn't have to be registered with &FindPrograms,
#                         but probably a good idea)
#              $options - [array ref] list of options
#              $where   - [optional] where to include these options on
#                         $program's command line -- 'pre' means at the 
#                         beginning, 'post' at the end (defaults to 'pre';
#                         if you need to set both pre and post options for
#                         a program, you need to call this routine twice)
#@OUTPUT     : none
#@RETURNS    : nothing you should care about
#@DESCRIPTION: Set standard options for a single program.
#@METHOD     : 
#@GLOBALS    : %PreOptions, %PostOptions
#@CALLS      : 
#@CREATED    : GPW, November 1996
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub SetProgramOptions
{
   my ($program, $options, $where) = @_;

   my $options_hash = (defined $where && $where eq 'post')
      ? \%PostOptions
      : \%PreOptions;
   my @programs = (ref $program eq 'ARRAY') ? (@$program) : ($program);

   foreach $program (@programs)
   {
      carp ("warning: setting options for unregistered program \"$program\"")
         if ($Options{'Strict'} && ! exists $Programs{$program});
      $options_hash->{$program} = $options;
   }
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &AddProgramOptions
#@INPUT      : (see &SetProgramOptions)
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Adds options to an existing list of program options, such
#              as is set with &SetProgramOptions.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : %PreOptions, %PostOptions
#@CREATED    : GPW, November 1996
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub AddProgramOptions
{
   my ($program, $options, $where) = @_;

   my $options_hash = (defined $where && $where eq 'post')
      ? \%PostOptions
      : \%PreOptions;
   my @programs = (ref $program eq 'ARRAY') ? (@$program) : ($program);

   foreach $program (@programs)
   {
      carp ("warning: adding options for unregistered program \"$program\"")
         if ($Options{'Strict'} && ! exists $Programs{$program});
      push (@{$options_hash->{$program}}, @$options);
   }
}



# ----------------------------------------------------------------------
# Private utility subroutines called by &Backgroundify, &Spawn, 
# and &Execute (and also by each other):
#    &set_undefined_options
#    &check_program
#    &complete_command
#    &exec_command
#    &gather_error
#    &check_return_status
#    &spawn_capture
#    &spawn_redirect
# ----------------------------------------------------------------------


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &set_undefined_options
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Used to inherit values of package options (eg.  Verbose,
#              Execute) from like-named variables in package `main' if
#              they haven't been explicitly set in this package's %Options
#              hash.
#@METHOD     : 
#@GLOBALS    : %Options, and possibly anything from main package 
#@CALLS      : 
#@CREATED    : GPW, March 1996 (copied from Batch module)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub set_undefined_options
{
   package main;
   no strict 'refs';

   my (@options) = @_;
   my $option;

   foreach $option (@options)
   {
      $JobControl::Options{$option} = $$option # neat trick -- but check it!!
         unless defined $JobControl::Options{$option};
   }
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &check_program
#@INPUT      : $program - the first word from a command, ie. the program
#                         to execute
#@OUTPUT     : 
#@RETURNS    : $fullpath    - full path (possibly relative) to the program;
#                             if the input $program was a path, then this
#                             will be identical to it
#              $options     - [array ref] list of options associated with
#                             this program
#              $post_options- [array ref] list of options associated with
#                             this program that are meant to come at the
#                             *end* of the command line
#@DESCRIPTION: Makes sure that a program we wish to execute actually exists.
#              How this is done depends on how the program name is specified;
#              if it is just a filename (no slashes), then we first look
#              in the %Programs hash to see if it was registered via
#              FindProgram; if not, we (optionally) complain and try to
#              find out ourselves via an explicit search.  If this fails,
#              we fail via &check_return_status (meaning we might just
#              return 0, or completely bomb out, or whatever the caller
#              asked for).
#
#              If, however, the program is supplied as a path of some sort
#              (ie. there're slashes in it), then we simply check that the
#              specified file exists and is executable.  If not, we again
#              bomb via &check_return_status.
#@METHOD     : 
#@GLOBALS    : %Programs, %PreOptions, %PostOptions
#@CALLS      : 
#@CREATED    : 1996/11/19, GW (from prototype code in &Execute)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub check_program
{
   my ($command, $program) = @_;
   my ($fullpath, @pre_options, @post_options);

   # If program was supplied as a bare filename, then we poke around a bit
   # to find it.  First try the %Programs hash; if it's not found there,
   # then we issue a warning (as long as the Strict flag is set) and try to
   # find it ourselves with &find_program.  If this also fails, then we
   # call &check_return_status with an artificial failure status (the same
   # status that Perl sets if system() fails because the program doesn't
   # exist) to provoke the appropriate failure action.  It's ok to do this
   # without printing "$program not found" because &find_program does this
   # for us.

   if ($program !~ m|/|)      # just a name, no directories at all
   {
      $fullpath = $Programs{$program};
      unless ($fullpath)
      {
         carp ("warning: program $program not registered with FindPrograms")
            if ($Options{'Strict'});

         unless ($fullpath = &::find_program ($program, $Options{'SearchPath'}))
         {
            $command = &::shellquote (@$command) if ref $command;
            return &check_return_status (255, $program, $command);
         }
      }
   }
   else                      # caller supplied a path (possibly relative,
   {                         # but that doesn't matter)
      unless (-e $program && -x $program)
      {
         warn "Couldn't find program \"$program\" (check permissions)\n";
         $command = &::shellquote (@$command) if ref $command;
         return &check_return_status (255, $program, $command);
      }

      $fullpath = $program;
      $program = (&::split_path ($fullpath, 'none'))[1];
   }

   # Now, $program is base filename, and $fullpath is the full path
   # to it (possibly relative, if that's what we were supplied with).
   # Look up program options in the two options hashes for use
   # by the caller.

   @pre_options = (exists $PreOptions{$program})
      ? @{$PreOptions{$program}}
      : ();
   @post_options = (exists $PostOptions{$program})
      ? @{$PostOptions{$program}}
      : ();

   ($fullpath, \@pre_options, \@post_options);
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &complete_command
#@INPUT      : $command  - [list ref or string]
#              $verbose  - true to print out full command
#@OUTPUT     : 
#@RETURNS    : $command  - [list ref or string -- whichever the input is]
#                          input $command with full path to program and
#                          any standard options included
#              $program  - the program name extracted from $command
#@DESCRIPTION: Extracts the program name from a command (list or string),
#              calls &check_program to get the full path and options
#              for the program, and constructs a new command (list or
#              string) with these goodies included.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : GPW, 1996/12/10, from &Execute
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub complete_command
{
   my ($command, $verbose) = @_;
   my ($program, $pre_options, $post_options);

   if (ref $command)
   {
      croak ("\$command must be either a list ref or a simple string")
         unless (ref $command eq 'ARRAY');

      # Extract the program to run, and then get more details (full path,
      # default options) from &check_program
         
      my @command = @$command;
      $program = shift @command;

      # Skip the fancy stuff if the program name is "escaped"

      if ($program =~ s|^\\||)
      {
         unshift (@command, $program);
      }
      else
      {
         my ($fullpath, $pre_options, $post_options) =
            &check_program ($command, $program);
         
         # Build a new command list using the full path and default options;
         # use this to print and execute
         
         unshift (@command, $fullpath, @$pre_options);
         push (@command, @$post_options);
      }

      # Finally, print the command if the Verbose flag is true

      if ($verbose)
      {
         no strict 'refs';              # so $lh (a string) can act as a FH
         my $lh = $Options{'LogHandle'};
         printf $lh ("[%s] %s %s %s\n",
                     $::ProgramName, &::userstamp(), &::timestamp(),
                     &::shellquote (@command));
      }

      # And return the augmented command list

      return (\@command, $program);
   }
   else                                 # $command is just a string
   {
      # Skip the fancy stuff if the program name is "escaped"

      if ($command =~ s|^\\||)
      {
         ($program) = split (/\s+/, $command, 1);
      }
      else
      {
         # Get the first word (i.e. the program to run) from the command
         # string; we do this via a regexp rather than split so we can
         # remove the program name and it back in later, with default
         # options added

         $command =~ s/^(\S+)\s*//;
         $program = $1;
         my ($fullpath, $pre_options, $post_options) = 
         &check_program ($command, $program);

         # Build a new command string using the full path and default options;
         # print the command and execute it

         $command = join (" ", $fullpath, @$pre_options, $command, @$post_options);
      }

      if ($verbose)
      {
         no strict 'refs';              # so $lh (a string) can act as a FH
         my $lh = $Options{'LogHandle'};
         printf $lh ("[%s] %s %s %s\n",
                     $::ProgramName, &::userstamp(), &::timestamp(), $command);
      }

      return ($command, $program);
   }
}



# ------------------------------ MNI Header ----------------------------------
#@NAME       : &exec_command
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Run by the child process after we fork -- this takes care
#              of redirecting stdout (if appropriate) and stderr (likewise),
#              exec's the command, and bombs if the exec fails.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : GPW, 1996/12/10
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub exec_command
{
   my ($command, $program, $stdout, $stderr) = @_;

   # Redirect stdout and/or stderr if the caller specified file
   # destination(s) (also possibly "capture" stderr via temp file)

   if ($stdout)
   {
      open (STDOUT, $stdout)
         || croak ("Unable to redirect stdout to \"$stdout\": $!\n");
   }

   if ($stderr)
   {
      $stderr = ">/tmp/error$$.log" if $stderr eq "-";
      open (STDERR, $stderr)
         || croak ("Unable to redirect stderr to \"$stderr\": $!");
   }

   # Exec that sucker!

   (ref $command)
      ? exec @$command
      : exec $command;

   # If we get here, the exec failed -- should not happen because of
   # the care we take to find the program in &check_program

   confess ("exec of $program failed: $!");
}



# ------------------------------ MNI Header ----------------------------------
#@NAME       : &gather_error
#@INPUT      : $pid   - process id of the now deceased child, presumed
#                       to have written its stderr to /tmp/error$pid.log#
#                       (in fact, we crash 'n burn if this file is not found)
#@OUTPUT     : 
#@RETURNS    : $error - the contents of that temporary file -- all lines
#                       are concatented together, but the newlines are
#                       preserved
#@DESCRIPTION: Reads in a temporary file that was created to hold a
#              child process' stderr.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : GPW, 1996/12/10, from &spawn_redirect and &spawn_capture
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub gather_error
{
   my ($pid) = @_;
   my ($filename, $error);

   $filename = "/tmp/error${pid}.log";

   open (ERROR, "<$filename") || confess ("Unable to open $filename: $!");
   $error = join ("", <ERROR>);
   close (ERROR);
   unlink ($filename) || warn ("Unable to delete $filename: $!");

   $error;
}
   
   
   
# ------------------------------ MNI Header ----------------------------------
#@NAME       : &check_return_status
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Checks the return status of a program.  If non-zero,
#              prints a warning and carries out the ErrorAction (if
#              defined; if not, just returns 0).
#
#              This is the only place where ErrorAction is interpreted,
#              so I suppose this is as good a place as any to document it.
#              Current two simple values are supported:
#                 `fatal'   - crash the whole program immediately
#                 `notify'  - send email to the address specified in the
#                             Notify option, and *then* crash the whole
#                             program
#
#              Alternately, it can be an arbitrary chunk of Perl code,
#              which is evaluated in the calling package; if there are
#              errors in this code, we again crash by calling croak.
#              Beware of making this code a single string of lowercase
#              letters, as I reserve the right to add more options like
#              `fatal' or `notify' that might cause conflicts.
#
#              Finally, if no ErrorAction is supplied, we print a warning
#              and return 0.  (That's return 0 to the caller, which is
#              most likely &Spawn or &Execute -- what they return is
#              up to them!)
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 96/01/17, GPW (from code in &Spawn)
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub check_return_status
{
   my ($status, $program, $command, $stdout, $stderr, $output, $error) = @_;
   my $ea = $Options{'ErrorAction'};

   # Note to myself: $status is now the full 16-bit value as returned
   # by wait(2).  I should do some bit-twiddling (preferably with the
   # help of "sys/wait.ph", but it's not cooperating right now) 
   # to figure out what kind of crash, exit code, signal number, etc.

   if ($status)
   {
      if ($ea eq "fatal")
      {
         &::Fatal ("crashed while running $program (exit status=$status)");
      }
      elsif ($ea eq "notify")
      {
         &_Obituary ($status, $program, $command,
                     $stdout, $stderr, $output, $error);
      }
      elsif ($ea)                       # some chunk of code to be eval'd
      {
	 my ($i, $this_pkg, $up_pkg, @caller);

         carp "You should be using `notify' rather than hard-coding a call to \&Obituary" 
            if $ea =~ /Obituary/;
         carp "You should be using `fatal' rather than hard-coding a `die' or call to \&Fatal" 
            if $ea =~ /Fatal|die/;

	 $i = 0;
	 $this_pkg = (caller(0))[0];
	 while (@caller = caller($i++))
	 {
	    $up_pkg = $caller[0];
	    last if $up_pkg ne $this_pkg;
	 }

	 eval "package $up_pkg; $ea";
	 croak ("Error in eval'd code: $@\n") if $@;
         return 0;
      }
      else
      {
         warn "$program crashed (exit status=$status)\n";
         return 0;
      }
   }
   else { return 1; }
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &spawn_capture
#@INPUT      : $command - [list ref or string] the full command
#              $program - the program name, as extracted by &complete_command
#              $stderr  - what to do with stderr.  Possible values:
#                  '-'  : capture to a separate variable from stdout
#                  ''   : leave intact
#                  any other string
#                       : redirect to the named file (which could be ">&STDOUT"
#                         to merge stderr with stdout)
#@OUTPUT     : 
#@RETURNS    : [array context] ($result,$output,$error)
#              [scalar context] $result
#                $result - exit status of $command, ie. non-zero means failure
#                $output - contents of child process' stdout
#                $error  - contents of child process' stderr (but only
#                          if stderr was captured! ie. this will be undefined
#                          unless $stderr was '-')
#              Note that calling this routine in a scalar context means
#              that your child process' stdout will be lost.
#@DESCRIPTION: Run a command, capturing stdout to a variable.  stderr might
#              also be captured, but then again it might be redirected
#              (including possibly merged with stdout) or left intact.
#@METHOD     : forks and execs desired command with a pipe between parent
#              and child; captures stderr via a temporary file if it has
#              to
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : GPW, 1996/12/10, from &Spawn
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub spawn_capture
{
   my ($command, $program, $stderr) = @_;
   my ($pid, $output, $status, $error);

   # Run the command through a pipe.  In the child, stdout automagically
   # goes to the pipe; in the parent, we will slurp the child's output from
   # the pipe.  stderr is dealt with according to $stderr -- either it's
   # capture (via a temporary file), or it's redirected to some other file,
   # or it's untouched.

   $pid = open (PIPE, "-|");
   croak ("&Spawn: failed to start child process: $!\n")
      unless defined $pid;

   if ($pid == 0)                       # in the child?
   {
      &exec_command ($command, $program, "", $stderr);
   }

   # $pid is not 0, so we're in the parent.  We read in the child's stdout
   # (and stderr too if it was captured).

   $output = join ("", <PIPE>);         # but leave newlines intact!
   close (PIPE);
   $status = $?;                        # get child's exit code

   if ($stderr eq "-")                  # did we "capture" stderr?
   {                                    # then read it in from temp file
      $error = &gather_error ($pid);
   }

   &check_return_status ($status, $program, $command,
                         undef, $stderr, $output, $error);

   wantarray 
      ? ($status, $output, $error)
      : $status;
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &spawn_redirect
#@INPUT      : $command - [list ref or string] the full command
#              $program - the program name, as extracted by &complete_command
#              $stdout  - what to do with stdout.  Can be empty string
#                         (leave intact) or name of file to redirect to.
#              $stderr  - what to do with stderr.  Possible values:
#                  '-'  : capture to a separate variable from stdout
#                  ''   : leave intact
#                  any other string
#                       : redirect to the named file (which could be ">&STDOUT"
#                         to merge stderr with stdout)
#@OUTPUT     : 
#@RETURNS    : $result  - exit status of $command, ie. non-zero means failure
#@DESCRIPTION: Run a command, redirecting stdout to a file or leaving it
#               intact.  stderr might likewise be redirected to file
#               (including possibly ">&STDOUT", which will merge it with
#               stdout), captured to a variable (via a temporary file), or
#               left intact.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : GPW, 1996/12/10, from &Spawn
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub spawn_redirect
{
   my ($command, $program, $stdout, $stderr) = @_;
   my ($pid, $status, $error);

   $pid = fork;
   croak ("&Spawn: failed to start child process: $!\n") unless defined $pid;

   if ($pid == 0)		# we're in the child process
   {
      &exec_command ($command, $program, $stdout, $stderr);
   }


   # $pid is not 0, so we're in the parent; block and wait for child.

   if (waitpid ($pid, 0) == -1)
   {
      confess ("no child processes even though I just forked!");
   }

   $status = $?;


   # More weirdness here: if I have the code like this, then
   # it actually calls &check_return_status here (or so
   # the debugger claims).

#   ($error = &gather_error ($pid))      # if stderr was put in a temp file,
#      if ($stderr eq "-")               # gather it up to a variable

   if ($stderr eq "-")
   {
      $error = &gather_error ($pid);
   }

   &check_return_status ($status, $program, $command,
                         $stdout, $stderr, undef, $error);

   wantarray 
      ? ($status, undef, $error)
      : $status;
}



# ----------------------------------------------------------------------
# The main workhorses:
#    &Backgroundify
#    &Spawn
#    &Execute
# ----------------------------------------------------------------------


# ------------------------------ MNI Header ----------------------------------
#@NAME       : Backgroundify
#@INPUT      : $logfile - name of file to which the stdout and stderr should
#                         be redirected (if empty, will try to use LogHandle
#                         package option)
#              @args    - arguments program was run with (used to generate
#                         announcement message) (usually just a copy of @ARGV)
#                         (if not supplied, won't print out banner describing
#                         how program was run)
#@OUTPUT     : none
#@RETURNS    : 0 if redirection of stdout or stderr failed (does NOT fork)
#              0 if caller didn't specify anywhere to redirect to, either
#                 with $logfile parameter or LogHandle package option
#              if successful: returns 1 to newly-detached child process
#                             never returns to parent (it terminates!)
#@DESCRIPTION: Redirects stdout and stderr to a file (named by $logfile),
#              and detaches to background execution.  If @args is non-
#              empty, prints two announcements: to stdout, that we're
#              detaching to the background; and to the logfile, the
#              user, host, time, and full argument list (@args).
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1995/02/21, Greg Ward: copied from nearly-identical function
#                 in do_mritopet
#@MODIFIED   : 1996/03/15-26, GW: changed to pay attention to Clobber
#                 and LogHandle options
#              1996/12/10, GW: changed to use &shellquote to print @args
#@COMMENTS   : 
#-----------------------------------------------------------------------------
sub Backgroundify
{
   my ($logfile, @args) = @_;
   my ($time, $user, $log_existed, $stdout);

   &set_undefined_options ("Verbose");

   # If caller didn't supply a log filename, check the package option
   # LogHandle -- if it's anything other than STDOUT, then we'll
   # redirect to it.

   if (! $logfile)
   {
      my $log_handle = $Options{'LogHandle'};
      if ($log_handle ne "STDOUT")
      {
	 print "Redirecting output and detaching to background\n"
	    if $Options{'Verbose'};
	 $stdout = ">&$log_handle";
      }
      else
      {
	 carp ("Backgroundify: must specify either \$logfile or LogHandle");
	 return 0;
      }
   }

   # Caller specified a filename, so use it and ignore LogHandle option.

   else
   {
      print "Redirecting output to $logfile and detaching to background\n"
	 if $Options{'Verbose'};
      $stdout = ($Options{'Clobber'} ? ">" : ">>") . "$logfile";
   }

   unless (open (STDOUT, $stdout))
   {
      warn "Detachment to background failed: couldn't redirect stdout ($!)\n";
      return 0;
   }
   unless (open (STDERR, ">&STDOUT"))
   {
      warn "Detachment to background failed: couldn't redirect stderr ($!)\n";
      return 0;
   }

   autoflush STDERR 1;
   autoflush STDOUT 1;

   if (@args)
   {
      print "\n" . "=" x 80 . "\n"
	 if ($log_existed && !$Options{'Clobber'});
      printf "%s %s opening log file for:\n", &::userstamp(), &::timestamp();
      print &::shellquote ($0, @args) . "\n\n";
   }

   fork && exit;
   return 1;			# to the *child* -- parent terminates!
}



# ------------------------------ MNI Header ----------------------------------
#@NAME       : &Spawn
#@INPUT      : $command - the command (and all arguments) to run
#              $stdout  - tells what to do with the standard output.
#                         Can be either a filename (redirect to that
#                         file), "-" (capture to a variable, which is
#                         returned by &Spawn), or empty (leave it alone).
#                         If $stdout is a filename, it will be either
#                         overwritten or appended depending on the
#                         Clobber option.
#              $stderr  - identical to $stdout.  Note that you
#                         can merge the standard error into the standard
#                         output by using ">&STDOUT" for $stderr (this is
#                         a standard Perl trick).
#              $stdin   - doesn't do anything (yet) (and probably never will)
#@OUTPUT     : 
#@RETURNS    : $status - exit status of the child process
#              $output - standard output of the child process (as long as
#                        $stdout was empty or "-")
#              $error  - standard error of the child process (as long as
#                        $stderr was empty or "-")
#@DESCRIPTION: Runs a child process, optionally capturing its standard
#              output and/or standard error to a file or a Perl variable.
#              To capture either stdout or stderr to a filename, just pass
#              filenames (eg. "foo.log") as $stdout and/or $stderr.  You
#              can merge stderr with stdout either by setting $stderr to
#              ">&STDOUT", or by leaving $stderr empty (or undefined) and
#              setting the MergeErrors option to a true value.
#
#              To capture to a variable, simply evaluate &Spawn in an array
#              context, eg. ($status, $output, $error) = &Spawn (...).  By
#              default, this will make stdout and stderr go to separate
#              variables; you can counteract this (as above) by setting
#              $stderr to ">&STDOUT" or by setting MergeErrors to 2.
#@METHOD     : 
#@GLOBALS    : %Options
#@CALLS      : 
#@CREATED    : 1995/03/03, Greg Ward (didn't really get anywhere)
#@MODIFIED   : 1995/05/16, GW - finished it (well, mostly)
#              1995/06/20, GW - added code for capturing stdout and stderr
#@COMMENTS   : there is a much better explanation of this lurking in
#              my "Tips for Perl Hackers" document: see ~greg/docs/perltips.dvi
#              on the BIC SGI systems
#-----------------------------------------------------------------------------
sub Spawn
{
   my ($command, $stdout, $stderr, $stdin) = @_;
   my ($program);
   my ($pid, $status, $output, $error); # must be dynamically scoped!!!
                                        # must they now...?

   &set_undefined_options ("Verbose", "Execute");

   $stdout = "" if (! defined $stdout); # for guilt-free string comparisons
   $stderr = "" if (! defined $stderr);
   $stdin = "" if (! defined $stdin);

   # Check for unimplemented features

   croak ("&Spawn: sorry, can't redirect stdin at all (yet)\n")
      if ($stdin);


   # Turn the basic command into a somewhat more ornamented beast --
   # ie. with the program turned into a full path, and with standard
   # options inserted fore and aft.

   ($command, $program) = &complete_command
      ($command, ($Options{'Verbose'} && !$Options{'Batch'}));


   # If Batch mode is on, submit this command and be done with it

   if ($Options{'Batch'} && $Options{'Execute'})
   {
      if (! defined &Batch::QueueCommand)
      {
         croak "Batch package not loaded -- you must do \"use Batch\"";
      }

      $stdout = "" if $stdout eq "-";   # can't capture
      $stderr = "" if $stderr eq "-";
      &Batch::QueueCommand ($command, "", $stdout, $stderr, $Options{'MergeErrors'});
      return 0;
   }

   # Capture both stdout and stderr to variables if we're evaluated
   # in an array context and a redirection destination (ie. file)
   # hasn't already been specified for them.

   if (wantarray)
   {
#      print "Spawn called in array context, possibly overriding stdout/stderr\n";
      $stdout = "-" unless $stdout;
      $stderr = "-" unless $stderr;
   }

   # If nothing was explicitly said about what to do with stderr,
   # then we'll merge it with stdout as long as MergeErrors is true.
   # If stderr is to be captured (either through caller explicitly 
   # passing in a "-", or through wantarray trickery above), then
   # we might still merge it with stdout -- but this time, 
   # MergeErrors has to be a little more true (ie. at least 2).

   $stderr = ">&STDOUT" if $stderr eq ""  && $Options{'MergeErrors'} >= 1;
   $stderr = ">&STDOUT" if $stderr eq "-" && $Options{'MergeErrors'} >= 2;

   
   # Figure out how to open files (> to overwrite or >> to append),
   # and hard-code that into $stdout and $stderr if they don't
   # already have such a code in them.

   my $open_prefix = ($Options{'Clobber'}) ? ">" : ">>";
   $stdout = $open_prefix . $stdout if ($stdout && $stdout !~ /^[>-]/);
   $stderr = $open_prefix . $stderr if ($stderr && $stderr !~ /^[>-]/);


   # Here's where we pay attention to the Execute flag: return now unless
   # it is true.  (Note that &Spawn's return value is bass-ackwards: 0
   # implies success, just like a program's exit code [this is not
   # coincidental!].)

   return 0 unless $Options{'Execute'};

   # Do things rather differently depending on whether we are to capture
   # stdout.  In the capture-to-variable case, we will fork-and-exec the
   # desired command, and slurp its output via a pipe; stderr will be
   # captured via a temporary file (if stderr is supposed to be captured --
   # it might be redirected or left alone).  (Yeah, I know stdout could
   # just as well have been captured the same way.  However, I was 1)
   # curious about doing it via a pipe, and 2) concerned about efficiency
   # -- error output is usually pretty small so going via disk won't kill
   # anyone, but stdout could be big.)  I didn't do stderr via a pipe
   # because I don't think you can have *two* pipes (one for stdout and one
   # for stderr) connecting the two processes.  Perhaps it could be done
   # with a named pipe... ah well, another day...

   if ($stdout eq "-")		# capturing to a variable
   {
      &spawn_capture ($command, $program, $stderr);
   }

   # If redirecting to file(s), do it right -- fork, do redirection
   # in the child process, and exec the command; in the parent, just
   # wait for the child to finish

   else
   {
      &spawn_redirect ($command, $program, $stdout, $stderr);
   }
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &Execute
#@INPUT      : 
#@OUTPUT     : 
#@RETURNS    : 1 if the program is found and the command executes succesfully
#              0 otherwise (error message printed to STDERR or mailed
#                to user, depending on ErrorAction)
#@DESCRIPTION: Spawn's simplified little cousin.  The name "spawn" should
#              bring forth images of great shaggy beasts rising from
#              the hellish unseen depths of the universe.  I don't know
#              what "execute" makes you think, but it's a much simpler
#              beast.
#
#              One potentially important advantage of Spawn, though, is
#              that if your child process dies due to being hit by a
#              signal (eg. Ctrl-C from the terminal, or kill command,
#              or resource overrun, or what have you), then that
#              will trigger an error in your program -- usually a Good
#              Thing, at least if you want to be very sensitive to
#              errors and crashes and that sort of thing.
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 1996/11/13, Greg Ward
#@MODIFIED   : 
#@COMMENTS   : this was created as an experimental playground to try out
#              the %Programs and %Options hashes, with eventual incorporation
#              of those features into the great shaggy beast a possibility.
#-----------------------------------------------------------------------------
sub Execute
{
   my ($command) = @_;
   my ($program, @rest, $status);

   &set_undefined_options ("Verbose", "Execute");

   $command = &complete_command ($command, $Options{'Verbose'});
   return 1 unless $Options{'Execute'};
   
   (ref $command)
      ? system @$command                # it's an array ref
      : system $command;                # it's a string

   # Check the return status (check_return_status will fail if
   # it's non-zero)
   
   $status = $?;
   return &check_return_status
      ($status, $program, &::shellquote (@$command));
}


sub Obituary                            # hack for old code that hard-codes
{                                       # "&Obituary(...)" in ErrorAction
   my ($command, $status, $error) = @_;
   my $program;

   carp "Using obsolete form of \&Obituary -- please change your ErrorAction" .
      " to `notify'";
   carp "Bad mix of obsolete \&Obituary and new command-as-list style"
      if (ref $command);

   (ref $command)
      ? ($program = $command->[0])
      : (($program) = split (/\s+/, $command, 1));
   &_Obituary ($command, $program, $status, $error);
}


sub describe_fate
{
   my ($program, $dest, $contents, $what, $was) = @_;
   my ($text);

   $program = "the child program" unless $program;

   # empty means we didn't do anything with it
   if (!$dest)
   {
      $text = <<TEXT;
$program\'s $what $was left intact; check 
$::ProgramName\'s output to see what went wrong.
TEXT
   }

   # "-" means it was captured to a variable, so it should be
   # in $output 
   elsif ($dest eq "-")
   {
      if ($contents)
      {
         $text  = "Here is ${program}'s $what:\n";
         $text .= $contents . "\n";
      }
      else
      {
         $text = "$program\'s $what $was lost";
      }
   }

   # otherwise, it must have gone to a file
   else
   {
#      $dest =~ s/^>>?//;
      $text = "$program\'s $what $was saved in the file $dest";
   }

   $text . "\n";
}


# ------------------------------ MNI Header ----------------------------------
#@NAME       : &_Obituary
#@INPUT      : $status
#              $program
#              $command
#              $stdout
#              $stderr
#              $error
#@OUTPUT     : 
#@RETURNS    : 
#@DESCRIPTION: Mails a reasonably detailed death notice to the address
#              specified in the Notify option, and then crashes (via
#              Fatal).
#@METHOD     : 
#@GLOBALS    : 
#@CALLS      : 
#@CREATED    : 
#@MODIFIED   : 
#-----------------------------------------------------------------------------
sub _Obituary                   # use this for long-running programs that
{				# send mail when done or crashed
   my ($status, $program, $command, $stdout, $stderr, $output, $error) = @_;
   my ($text, $cwd, $stdout_desc, $stderr_desc);

   &set_undefined_options ("Notify", "Execute");

   $text = "";
   $cwd = getcwd();

   $stdout = "" if (! defined $stdout); # for guilt-free string comparisons
   $stderr = "" if (! defined $stderr);
   $stdout =~ s/^>>?//;                 # don't care about clobber state!
   $stderr =~ s/^>>?//;

   if ($Options{'Notify'} && $Options{'Execute'})
   {
      $text .= "$::ProgramName crashed" . 
         ($program ? " while running $program" : "") . "\n";
      $text .= "current dir: $cwd\n";

      if ($command)
      {
         $command = &::shellquote (@$command) if ref $command;
         $text .= "command    : $command\n";
      }

      # Now the fun bit: try to give the poor user some idea of what
      # happened to his stdout and stderr.  First, tell them about
      # stdout -- but if stderr was merged into it, we're really 
      # talking about both of them, so let's make it look that way.

      my ($what, $was);
      if ($stderr eq '&STDOUT')
      {
         $what = "standard output and standard error";
         $was = "were";
         $stderr_desc = "";
      }
      else
      {
         $what = "standard output";
         $was = "was";
         $stderr_desc = &describe_fate
            ($program, $stderr, $error, "standard error", "was");
      }

      $stdout_desc = &describe_fate ($program, $stdout, $output, $what, $was);

      open (MAIL, "|/usr/lib/sendmail $Options{'Notify'}");

      print MAIL <<EOM;
From: ($::ProgramName)
Subject: $::ProgramName crashed!

Dear $Options{'Notify'},

$text
$stdout_desc
$stderr_desc

EOM

      close (MAIL);
   }

   if ($program)
      { &::Fatal ("crashed while running $program"); }
   else
      { &::Fatal ("crashed"); }
}

1;
