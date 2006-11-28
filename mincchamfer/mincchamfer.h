/* ----------------------------- MNI Header -----------------------------------
@NAME       : mincchamfer.h
@DESCRIPTION: this file contains the list of include files and global variables
              needed for mincchamfer.h
@CREATED    : Nov 2, 1998 louis
@MODIFIED   :
---------------------------------------------------------------------------- */

                                /* include list */
#include <volume_io.h>
#include <Proglib.h>
#include <ParseArgv.h>

                                /* globals */
char 
  *prog_name;
int
  first,
  debug,
  verbose;
Real
  max_dist;

static ArgvInfo argTable[] = {
  {"-max_dist", ARGV_FLOAT, (char*)0, (char*)&max_dist,
     "Maximum distance value in transform."},
  {"-first", ARGV_INT, (char *) 0, (char *) &first,
     "Number of initial background structure dilations"},
  {NULL, ARGV_HELP, NULL, NULL,
     "Options for logging progress. Default = -verbose."},
  {"-verbose", ARGV_CONSTANT, (char *) TRUE, (char *) &verbose,
     "Write messages indicating progress"},
  {"-quiet", ARGV_CONSTANT, (char *) FALSE , (char *) &verbose,
     "Do not write log messages"},
  {"-debug", ARGV_CONSTANT, (char *) TRUE, (char *) &debug,
     "Print out debug info."},
  {NULL, ARGV_END, NULL, NULL, NULL}
};


