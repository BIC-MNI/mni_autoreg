/* ----------------------------- MNI Header -----------------------------------
@NAME       : random_vol
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program designed to build a random volume with random intensities 
              that have a gaussian distribution.
@CREATED    : Dec 5 1995 Collins
@MODIFIED   : $Log: random_vol.c,v $
@MODIFIED   : Revision 1.1  1999-10-25 19:52:09  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
@COPYRIGHT  :
              Copyright 1995 Louis Collins, McConnell Brain Imaging Centre, 
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
#include <internal_volume_io.h>
#include <sys/types.h>
#include <time.h>


time_t time(time_t *tloc);

/* Main program */
char *prog_name;


#include <stdlib.h>
#include <math.h>


Real gaussian_random_0_1()
{
  static int iset=0;
  static Real gset;
  Real fac,r,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*drand48()-1.0;
      v2=2.0*drand48()-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0 || r == 0.0);
    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}



int main(int argc, char *argv[])
{
   Volume 
     volume;
   Real
     variability,
     rand_val,
     steps[MAX_DIMENSIONS],
     voxel[MAX_DIMENSIONS];
   int 
     progress_count,
     sizes[MAX_DIMENSIONS],
     index[MAX_DIMENSIONS],
     i,j,k;
   progress_struct
     progress;

   union
     {
       long   l;
       char   c[4];
     } seedval;
   
   time_t t;
   char tmp;

   prog_name = argv[0];

   /* Check arguments */
   if (argc != 5) {
      (void) fprintf(stderr, "Usage: %s output.mnc vol_size vol_seps std_dev_mm\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

   sizes[0] = sizes[1] = sizes[2] = atoi( argv[2] );
   steps[0] = steps[1] = steps[2] = atof( argv[3] );
   variability = atof( argv[4] );

				/* initialize drand function */
   t = 2*time(NULL);
   seedval.l = t; 
   tmp = seedval.c[0]; seedval.c[0] = seedval.c[3]; seedval.c[3] = tmp; 
   tmp = seedval.c[1]; seedval.c[1] = seedval.c[2]; seedval.c[2] = tmp;
   srand48(seedval.l);

				/* create data volume  */
   volume = create_volume(3, (char **) NULL, NC_SHORT, TRUE, 0.0, 0.0);
   set_volume_sizes(volume, sizes);
   set_volume_separations(volume, steps);
   set_volume_real_range(volume, (Real)(-5.0*variability), (Real)(5.0*variability));
   alloc_volume_data(volume);

   initialize_progress_report(&progress, FALSE, 
			      sizes[X]*sizes[Y]*sizes[Z]+1,
			      "Randomizing volume");
   progress_count = 0;


   for_less(i,0,MAX_DIMENSIONS) index[i] = 0;

   /* loop over all voxels */
   for_less( index[ X ], 0, sizes[ X ])
     for_less( index[ Y ], 0, sizes[ Y ])
       for_less( index[ Z ], 0, sizes[ Z ]) {
	     
	 /* get a random value from a gaussian distribution  */

	 rand_val = variability * gaussian_random_0_1();

	 if (rand_val >  5.0*variability) rand_val =  5.0*variability;
	 if (rand_val < -5.0*variability) rand_val = -5.0*variability;
	 
	 set_volume_real_value(volume,
			       index[0],index[1],index[2],index[3],index[4],
			       rand_val);
	 
	 progress_count++;
	 update_progress_report(&progress, progress_count);
       }
   
   terminate_progress_report(&progress);

   

   /* Write out the random volume */
   if (output_volume(argv[1], NC_UNSPECIFIED, FALSE, 0.0, 0.0, volume, 
		     (char *)NULL, (minc_output_options *)NULL) != OK) {
     (void) fprintf(stderr, "%s: Error writing volume file %s\n",
                     argv[0], argv[1]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}
