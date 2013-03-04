/* ----------------------------- MNI Header -----------------------------------
@NAME       : inbetween_def
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: program designed to calculate an intermediate deformation
              between two input deformations.

@CREATED    : May 5, 1996
@MODIFIED   : $Log: inbetween_def.c,v $
@MODIFIED   : Revision 1.4  2006-11-29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.3  2005/07/20 20:45:46  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.2  2002/03/26 14:15:29  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:06  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
@COPYRIGHT  :
              Copyright 1996 Louis Collins, McConnell Brain Imaging Centre, 
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
static char rcsid[]="$Header:";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <volume_io.h>
#include <sys/types.h>
#include <time.h>



/* Main program */
char *prog_name;

int main(int argc, char *argv[])
{
   VIO_General_transform 
     *grid_transform_ptr,
     transform;
   VIO_Volume 
     volume;
   VIO_Real
     fraction,
     upper_limit, lower_limit;
   int trans_count;
   prog_name = argv[0];



   /* Check arguments */
   if (argc != 4) {
      (void) fprintf(stderr, "Usage: %s <input.xfm> <result.xfm> <fraction>\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

   fraction = atof( argv[3] );

   if (fraction<=0.0) {
      (void) fprintf(stderr, "%s: Fraction cannot be negative <%f>.\n",
                     argv[0], fraction);
      exit(EXIT_FAILURE);
   }

   /* Read deformation field to be used as template */
   if (input_transform_file(argv[1], &transform) != OK) {
      (void) fprintf(stderr, "%s: Error reading transform file %s\n",
                     argv[0], argv[1]);
      exit(EXIT_FAILURE);
   }

   trans_count = 0;
   do {
     grid_transform_ptr = get_nth_general_transform(&transform, trans_count );
     trans_count++;
   } while ((grid_transform_ptr->type != GRID_TRANSFORM) &&
            (trans_count < get_n_concated_transforms(&transform)));
   
   if (grid_transform_ptr->type != GRID_TRANSFORM) {
      (void) fprintf(stderr, "Error: no deformation field to use.\n");
      exit(EXIT_FAILURE);
   }

   volume = grid_transform_ptr->displacement_volume;

   get_volume_real_range(volume, &lower_limit, &upper_limit);
   upper_limit *= fraction;
   lower_limit *= fraction;
   set_volume_real_range(volume, lower_limit, upper_limit);

   /* Write out the random transform */
   if (output_transform_file(argv[2], NULL, &transform) != OK) {
      (void) fprintf(stderr, "%s: Error writing transform file %s\n",
                     argv[0], argv[2]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}
