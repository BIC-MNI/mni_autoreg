/* ----------------------------- MNI Header -----------------------------------
@NAME       : rand_param.c
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: program to create a list of random parameters that can be fed
              to param2xfm to create a random transformation matrix.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Aug 23 15:32:30 EST 1994 - Louis
@MODIFIED   : $Log: rand_param.c,v $
@MODIFIED   : Revision 1.1  1995-02-22 08:56:06  louis
@MODIFIED   : Initial revision
@MODIFIED   :

---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <volume_io.h>
#include <ParseArgv.h>

#include <sys/types.h>
#include <time.h>


time_t time(time_t *tloc);

char *prog_name;


/* Main program */

int main(int argc, char *argv[])
{
   static Real 
     mag_scales, mag_trans, mag_rots, mag_skews,
     scales[3], trans[3], rots[3], skews[3];
   int i;
   long
     seedval;
   time_t t;
   char tmp, *b1,*b2;

   static ArgvInfo argTable[] = {
     {"-translation", ARGV_FLOAT, (char *) 0, (char *) &mag_trans,
	"Translation x,y,z."},
     {"-rotations",   ARGV_FLOAT, (char *) 0, (char *) &mag_rots,
	"Rotation angle (in degrees)."},
     {"-scales",      ARGV_FLOAT, (char *) 0, (char *) &mag_scales,
	"Scaling factors."},
     {"-shears",      ARGV_FLOAT, (char *) 0, (char *) &mag_skews,
	"Scaling factors."},
     {NULL, ARGV_END, NULL, NULL, NULL}
   };
   
   
   prog_name = argv[0];

   mag_scales = mag_trans =  mag_rots =  mag_skews = 0.0;

   for_less(i,0,3) {
     trans[i]  = 0.0;
     rots[i]   = 0.0;
     scales[i] = 1.0;
     skews[i]  = 0.0;
   }

   t = 2*time(NULL);
   seedval = t; 

   b1 = &seedval;
   b2 = b1+3;
   tmp = *b1;  *b1 = *b2; *b2 = tmp;
   b1++; b2--;
   tmp = *b1;  *b1 = *b2; *b2 = tmp;

   srand48(seedval);
   
   if (ParseArgv(&argc, argv, argTable, 0) || (argc!=1)) {
      (void) fprintf(stderr, "Usage: %s [options]\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

   for_less(i,0,3) {
     rots[i]   = 2.0*(-0.5 + drand48()) * mag_rots;
     scales[i] = 1.0 + 2.0*(-0.5 + drand48()) * mag_scales;
     trans[i]  = 2.0*(-0.5 + drand48()) * mag_trans;
     skews[i]  = 2.0*(-0.5 + drand48()) * mag_skews;
   }

   print ("-rotation %f %f %f -translation %f %f %f -scale %f %f %f -shear %f %f %f\n",
	  rots[0],   rots[1],   rots[2],
	  trans[0],  trans[1],  trans[2],
	  scales[0], scales[1], scales[2],
	  skews[0],  skews[1],  skews[2]);

   exit(EXIT_SUCCESS);
}
