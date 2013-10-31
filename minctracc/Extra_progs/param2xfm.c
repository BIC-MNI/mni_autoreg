/* ----------------------------- MNI Header -----------------------------------
@NAME       : param2xfm.c
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to check the z-scale of an  MNI transform file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon Nov 29 11:01:47 EST 1993 Louis
@MODIFIED   : $Log: param2xfm.c,v $
@MODIFIED   : Revision 96.6  2006-11-29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.5  2005/07/20 20:45:46  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.4  2002/03/26 14:15:31  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.3  2000/03/15 08:42:37  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 96.2  2000/02/20 04:01:02  stever
@MODIFIED   : * use new history_string() function to generate history strings
@MODIFIED   :   when outputting MNI files (.mnc, .xfm)
@MODIFIED   : * removed unused vax routines from Proglib
@MODIFIED   : * tuned configure script; CPPFLAGS and LDFLAGS are now left alone,
@MODIFIED   :   for the installer to use
@MODIFIED   :
@MODIFIED   : Revision 96.1  1999/10/25 19:52:07  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:36  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:30  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:18  louis
 * Never released version 0.95
 *
 * Revision 1.5  1996/08/12  14:15:09  louis
 * Pre-release
 *
 * Revision 1.4  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.3  94/04/26  12:54:34  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.2  94/04/06  11:48:46  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.1  94/02/21  16:36:03  louis
 * Initial revision
 * 
@COPYRIGHT  :
              Copyright 1993 Louis Collins, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <volume_io.h>
#include <config.h>
#include <Proglib.h>
#include <ParseArgv.h>
#include "make_rots.h"


char *prog_name;


/* Main program */

int main(int argc, char *argv[])
{
   VIO_General_transform new_transform;
   VIO_Transform
     lt;
   char* history = history_string( argc, argv );
   static VIO_Real 
     scales[3], trans[3], rots[3], skews[3], center[3];
   static int clobber = FALSE;
   int i;

   static ArgvInfo argTable[] = {
     {"-center",      ARGV_FLOAT, (char *) 3, (char *)center,
        "Force center of rotation and scale."},
     {"-translation", ARGV_FLOAT, (char *) 3, (char *)trans,
        "Translation x,y,z."},
     {"-rotations",   ARGV_FLOAT, (char *) 3, (char *)rots,
        "Rotation angle (in degrees)."},
     {"-scales",      ARGV_FLOAT, (char *) 3, (char *)scales,
        "Scaling factors."},
     {"-shears",      ARGV_FLOAT, (char *) 3, (char *)skews,
        "Scaling factors."},
     {"-clobber",     ARGV_CONSTANT, (char *) TRUE, (char *) &clobber,
        "Overwrite existing file (default = no clobber)."},
     {"-version", ARGV_FUNC, (char *) print_version_info, (char *)MNI_AUTOREG_LONG_VERSION,
          "Print out version info and exit."},
     {NULL, ARGV_END, NULL, NULL, NULL}
   };
   
   
   prog_name = argv[0];

   for(i=0; i<3; i++) {
     trans[i] = 0.0;
     center[i] = 0.0;
     rots[i] = 0.0;
     scales[i] = 1.0;
     skews[i] = 0.0;
   }


   if (ParseArgv(&argc, argv, argTable, 0) || (argc!=2)) {
      (void) fprintf(stderr, "Usage: %s [options] <result.xfm>\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }


   if (!clobber && file_exists(argv[1])) {
      (void) fprintf(stderr, "%s: file exists already - %s\n",
                     argv[0], argv[1]);
      exit(EXIT_FAILURE);
   }

   for(i=0; i<3; i++) {
     
     rots[i] = rots[i] * 3.1415927/ 180;
   }

   build_transformation_matrix(&lt,
                               center,
                               trans,
                               scales,
                               skews,
                               rots);

   create_linear_transform(&new_transform, &lt);
   

   /* Write out the transform */
   if (output_transform_file(argv[1], history, &new_transform) != VIO_OK) {
      (void) fprintf(stderr, "%s: Error writing transform file %s\n",
                     argv[0], argv[1]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}
