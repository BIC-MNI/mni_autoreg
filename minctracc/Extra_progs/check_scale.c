/* ----------------------------- MNI Header -----------------------------------
@NAME       : check_scale.c
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program to check the z-scale of an  MNI transform file
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon Nov 29 11:01:47 EST 1993 Louis
@MODIFIED   : $Log: check_scale.c,v $
@MODIFIED   : Revision 1.3  1994-04-26 12:52:49  louis
@MODIFIED   : updated with new versions of make_rots, extract2_parameters_from_matrix 
@MODIFIED   : that include proper interpretation of skew.
@MODIFIED   :
 * Revision 1.2  94/04/06  11:46:45  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.1  94/02/21  16:31:52  louis
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

char *prog_name;

/* Main program */

int main(int argc, char *argv[])
{
   General_transform transform, new_transform;
   Transform
     *lt;
   Real scales[3], trans[3], rots[3], skews[3], center[3];
   int i;

   prog_name = argv[0];

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


   if (get_transform_type(&transform) == CONCATENATED_TRANSFORM) {
     (void) fprintf(stderr, "Error: Cannot deal with concatenated transforms\n");
     exit(EXIT_FAILURE);
   }
   if (get_transform_type(&transform) == THIN_PLATE_SPLINE) {
     (void) fprintf(stderr, "Error: Cannot deal with non-linear transforms\n");
     exit(EXIT_FAILURE);
   }
   if (get_transform_type(&transform) == USER_TRANSFORM) {
     (void) fprintf(stderr, "Error: Cannot deal with user-defined transforms\n");
     exit(EXIT_FAILURE);
   }

   
   /* Extract parameters from transform */


   lt = get_linear_transform_ptr(&transform);
   for_less(i,0,3)
     center[i] = 0.0;
   if (!extract2_parameters_from_matrix(lt,
				       center,
				       trans,
				       scales,
				       skews,
				       rots)) {
     (void) fprintf(stderr, "Error: Cannot extract parameters from matrix\n");
     exit(EXIT_FAILURE);
   }
	
   /* check scaling parameters */


   if (scales[2] > 1.15*(scales[0]+scales[1])/2.0) {
     scales[2] = (scales[0]+scales[1])/2.0;
   }
   
   /* rebuild transformation from parameters */

   build_transformation_matrix(lt,
			       center,
			       trans,
			       scales,
			       skews,
			       rots);

   create_linear_transform(&new_transform, lt);
   

   /* Write out the transform */
   if (output_transform_file(argv[2], NULL, &new_transform) != OK) {
      (void) fprintf(stderr, "%s: Error writing transform file %s\n",
                     argv[0], argv[2]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}
