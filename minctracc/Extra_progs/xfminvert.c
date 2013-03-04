/* ----------------------------- MNI Header -----------------------------------
@NAME       : xfminvert
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to invert an MNI transform file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : August 13, 1993 (Peter Neelin)
@MODIFIED   : $Log: xfminvert.c,v $
@MODIFIED   : Revision 96.5  2006-11-29 09:09:32  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.4  2005/07/20 20:45:47  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.3  2004/02/12 05:54:06  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 96.2  2002/03/26 14:15:35  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.1  1999/10/25 19:52:12  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:36  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:32  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:18  louis
 * Never released version 0.95
 *
 * Revision 1.4  1996/08/12  14:15:12  louis
 * Pre-release
 *
 * Revision 1.3  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.2  94/04/06  11:49:01  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.1  94/02/21  16:37:49  louis
 * Initial revision
 * 
 * Revision 1.4  93/10/12  12:52:37  neelin
 * Replaced def_mni.h with volume_io.h
 * 
 * Revision 1.3  93/09/16  09:40:24  neelin
 * Use dave's open_file_with_default_suffix and input_transform_file and
 * output_transform_file to add suffixes to file names.
 * 
 * Revision 1.2  93/09/01  15:59:05  neelin
 * Cast return of fclose to (void).
 * 
 * Revision 1.1  93/08/13  15:26:53  neelin
 * Initial revision
 * 
@COPYRIGHT  :
              Copyright 1993 Peter Neelin, McConnell Brain Imaging Centre, 
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/xfminvert.c,v 96.5 2006-11-29 09:09:32 rotor Exp $";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <volume_io.h>
#include <minc_def.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

/* Main program */

int main(int argc, char *argv[])
{
   VIO_General_transform transform, inverse;

   /* Check arguments */
   if (argc != 3) {
      (void) fprintf(stderr, "Usage: %s <input.xfm> <result.xfm>\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

   /* Read in file to invert */
   if (input_transform_file(argv[1], &transform) != OK) {
      (void) fprintf(stderr, "%s: Error reading transform file %s\n",
                     argv[0], argv[1]);
      exit(EXIT_FAILURE);
   }

   /* Invert the transform */
   create_inverse_general_transform(&transform, &inverse);

   /* Write out the transform */
   if (output_transform_file(argv[2], NULL, &inverse) != OK) {
      (void) fprintf(stderr, "%s: Error writing transform file %s\n",
                     argv[0], argv[2]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}
