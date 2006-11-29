/* ----------------------------- MNI Header -----------------------------------
@NAME       : random_def
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program designed to build a random deformation field, given
        an existing field as a template.

@CREATED    : Tue Oct 24 14:22:51 MET 1995 Collins
@MODIFIED   : $Log: random_def.c,v $
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
@MODIFIED   : Revision 1.2  2002/03/26 14:15:31  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:09  louis
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
#include <volume_io.h>
#include <sys/types.h>
#include <time.h>



void get_volume_XYZV_indices(VIO_Volume data, int xyzv[])
{
  
  int 
    axis, i, vol_dims;
  char 
    **data_dim_names;
 
  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);
 
  for(i=0; i<N_DIMENSIONS+1; i++) xyzv[i] = -1;
  for(i=0; i<vol_dims; i++) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i; 
    } 
    else {     /* not a spatial axis */
      xyzv[Z+1] = i;
    }
  }
   delete_dimension_names(data, data_dim_names);

}

time_t time(time_t *tloc);

VIO_Real uniform_random_in_range(VIO_Real min, VIO_Real max);
VIO_Real gaussian_random_w_std(VIO_Real sigma);

/* Main program */
char *prog_name;

int main(int argc, char *argv[])
{
   VIO_General_transform 
     transform, 
     random_trans, 
     *grid_transform_ptr, 
     forward_transform;
   VIO_Volume 
     volume;
   VIO_Real
     variability,
     upper_limit, lower_limit,
     voxel[VIO_MAX_DIMENSIONS],
     def_values[VIO_MAX_DIMENSIONS];
   int 
     prog_count,
     sizes[VIO_MAX_DIMENSIONS],
     xyzv[VIO_MAX_DIMENSIONS],
     index[VIO_MAX_DIMENSIONS],
     i,j,k,
     trans_count;
   VIO_progress_struct
     progress;

   prog_name = argv[0];



   /* Check arguments */
   if (argc != 4) {
      (void) fprintf(stderr, "Usage: %s <input.xfm> <result.xfm> <std dev mm>\n",
                     argv[0]);
      exit(EXIT_FAILURE);
   }

   variability = atof( argv[3] );
   if (variability<=0.0) {
      (void) fprintf(stderr, "%s: Variability cannot be negative <%f>.\n",
                     argv[0], variability);
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
      (void) fprintf(stderr, "Error: no deformation field to use as template in\n");
      (void) fprintf(stderr, "       input file '%s'\n", argv[1]);
      exit(EXIT_FAILURE);
   }

   volume = grid_transform_ptr->displacement_volume;

   get_volume_sizes(volume, sizes);
   get_volume_XYZV_indices(volume,xyzv);
   
   upper_limit =  3.0*variability;
   lower_limit = -3.0*variability;
   set_volume_real_range(volume, lower_limit, upper_limit);


   for(i=0; i<MAX_DIMENSIONS; i++) index[i] = 0;

                                /* initialize drand function */
   init_random();

   for(index[xyzv[VIO_X]]=0; index[xyzv[VIO_X]]<sizes[xyzv[VIO_X]]; index[xyzv[VIO_X]]++)
     for(index[xyzv[VIO_Y]]=0; index[xyzv[VIO_Y]]<sizes[xyzv[VIO_Y]]; index[xyzv[VIO_Y]]++)
       for(index[xyzv[VIO_Z]]=0; index[xyzv[VIO_Z]]<sizes[xyzv[VIO_Z]]; index[xyzv[VIO_Z]]++) 
         for(index[xyzv[Z+1]]=0; index[xyzv[Z+1]]<3; index[xyzv[Z+1]]++)
           set_volume_real_value(volume,
                                 index[0],index[1],index[2],index[3],index[4],
                                 0.0);
   

   initialize_progress_report(&progress, FALSE, 
                              sizes[xyzv[VIO_X]]*sizes[xyzv[VIO_Y]]*sizes[xyzv[VIO_Z]]+1,
                              "Randomizing def field");
   prog_count = 0;

   for(index[xyzv[VIO_X]]=1; index[xyzv[VIO_X]]<sizes[xyzv[VIO_X]]-1; index[xyzv[VIO_X]]++)
     for(index[xyzv[VIO_Y]]=1; index[xyzv[VIO_Y]]<sizes[xyzv[VIO_Y]]-1; index[xyzv[VIO_Y]]++)
       for(index[xyzv[VIO_Z]]=1; index[xyzv[VIO_Z]]<sizes[xyzv[VIO_Z]]-1; index[xyzv[VIO_Z]]++) {
             
         for(i=X; i<=Z; i++) {
           def_values[i] = gaussian_random_w_std(variability);
           if (def_values[i] < lower_limit)
             def_values[i] = lower_limit;
           if (def_values[i] > upper_limit)
             def_values[i] = upper_limit;
         }
           /*or: variablity * uniform_random_in_range(-1.0, 1.0);*/
         
         for(index[xyzv[Z+1]]=0; index[xyzv[Z+1]]<3; index[xyzv[Z+1]]++)
           set_volume_real_value(volume,
                                 index[0],index[1],index[2],index[3],index[4],
                                 def_values[ index[ xyzv[Z+1] ]]);
         
         prog_count++;
         update_progress_report(&progress, prog_count);
       }
   
   terminate_progress_report(&progress);

   

   /* Write out the random transform */
   if (output_transform_file(argv[2], NULL, &transform) != OK) {
      (void) fprintf(stderr, "%s: Error writing transform file %s\n",
                     argv[0], argv[2]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}
