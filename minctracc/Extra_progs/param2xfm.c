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
@MODIFIED   : Revision 1.1  1994-02-21 16:36:03  louis
@MODIFIED   : Initial revision
@MODIFIED   :
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
#include <ParseArgv.h>

char *prog_name;


/* Main program */

int main(int argc, char *argv[])
{
   General_transform transform, new_transform;
   Transform
     lt;
   static Real 
     scales[3], trans[3], rots[3], skews[3], center[3];
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
     {NULL, ARGV_END, NULL, NULL, NULL}
   };
   
   
   prog_name = argv[0];

   for_less(i,0,3) {
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

   for_less(i,0,3) {
     
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
   if (output_transform_file(argv[1], NULL, &new_transform) != OK) {
      (void) fprintf(stderr, "%s: Error writing transform file %s\n",
                     argv[0], argv[1]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}