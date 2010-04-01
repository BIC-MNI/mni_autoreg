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
@MODIFIED   : Revision 96.6  2010-04-01 04:49:16  rotor
@MODIFIED   :  * fixed time.h include
@MODIFIED   :
@MODIFIED   : Revision 96.5  2008/11/07 23:54:54  sjschen
@MODIFIED   :
@MODIFIED   :
@MODIFIED   : Change the way the seed value is generated for srand(). Now gives different
@MODIFIED   : values when program in run multiple times under a second. Uses the product
@MODIFIED   : of pid and microsecond time as seed.
@MODIFIED   :
@MODIFIED   : Revision 96.4  2006/11/29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.3  2005/07/20 20:45:46  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.2  2002/03/26 14:15:31  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.1  1999/10/25 19:52:09  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:21:36  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:21:31  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:15:18  louis
 * Never released version 0.95
 *
 * Revision 1.2  1996/08/12  14:15:10  louis
 * Pre-release
 *
 * Revision 1.1  1995/02/22  08:56:06  collins
 * Initial revision
 *

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

#include <sys/types.h>
#include <sys/time.h>


time_t time(time_t *tloc);

char *prog_name;


/* Main program */

int main(int argc, char *argv[])
{
   static VIO_Real 
     mag_scales, mag_trans, mag_rots, mag_skews,
     scales[3], trans[3], rots[3], skews[3];
   int i;

   union
     {
       long   l;
       char   c[4];
     } seedval;
   
   time_t t;
   char tmp;

   static ArgvInfo argTable[] = {
     {"-translation", ARGV_FLOAT, (char *) 0, (char *) &mag_trans,
        "Translation x,y,z."},
     {"-rotations",   ARGV_FLOAT, (char *) 0, (char *) &mag_rots,
        "Rotation angle (in degrees)."},
     {"-scales",      ARGV_FLOAT, (char *) 0, (char *) &mag_scales,
        "Scaling factors."},
     {"-shears",      ARGV_FLOAT, (char *) 0, (char *) &mag_skews,
        "Scaling factors."},
     {"-version", ARGV_FUNC, (char *) print_version_info, (char *)MNI_AUTOREG_LONG_VERSION,
          "Print out version info and exit."},
    {NULL, ARGV_END, NULL, NULL, NULL}
   };
   
   
   prog_name = argv[0];

   mag_scales = mag_trans =  mag_rots =  mag_skews = 0.0;

   for(i=0; i<3; i++) {
     trans[i]  = 0.0;
     rots[i]   = 0.0;
     scales[i] = 1.0;
     skews[i]  = 0.0;
   }

/*    t = 2*time(NULL); */

   struct timeval tv;
   gettimeofday(&tv, NULL);
   t = tv.tv_usec * getpid(); /* use usec time and program id as seed*/

   seedval.l = t; 

   tmp = seedval.c[0]; seedval.c[0] = seedval.c[3]; seedval.c[3] = tmp; 
   tmp = seedval.c[1]; seedval.c[1] = seedval.c[2]; seedval.c[2] = tmp;
   
   srand48(seedval.l);
   
   if (ParseArgv(&argc, argv, argTable, 0) || (argc!=1)) {
      (void) fprintf(stderr, "Usage: %s [options]\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

   for(i=0; i<3; i++)
   {
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
