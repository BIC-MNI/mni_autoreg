/* ----------------------------- MNI Header -----------------------------------
@NAME       : flipdef
@INPUT      : argc, argv - command line arguments
@OUTPUT     : (none)
@RETURNS    : status
@DESCRIPTION: Program designed to invert the representation of 
        a deformation field to speed up resampling into the target
        space. 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mon May 29 09:07:14 MET DST 1995 Collins
@MODIFIED   : $Log: reversedef.c,v $
@MODIFIED   : Revision 1.6  2006-11-29 09:09:31  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.5  2005/07/20 20:45:47  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 1.4  2004/02/12 05:54:06  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 1.3  2002/03/26 14:15:32  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 1.2  2000/02/07 19:33:03  stever
@MODIFIED   : replaced HAVE_RECENT_VOLUME_IO with more specific feature tests.
@MODIFIED   :
@MODIFIED   : Revision 1.1  1999/10/25 19:52:10  louis
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/reversedef.c,v 1.6 2006-11-29 09:09:31 rotor Exp $";
#endif

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <volume_io.h>
#include <ParseArgv.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif

void print_usage_and_exit(char *pname);

void get_volume_XYZV_indices(VIO_Volume data, int xyzv[]){
  
  int 
    axis, i, vol_dims;
  char 
    **data_dim_names;
 
  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);
 
  for(i=0; i<VIO_N_DIMENSIONS+1; i++) xyzv[i] = -1;
  for(i=0; i<vol_dims; i++) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i; 
    } 
    else {     /* not a spatial axis */
      xyzv[VIO_Z+1] = i;
    }
  }

  delete_dimension_names(data, data_dim_names);
 
}

 void  general_transform_point_in_trans_plane(
    VIO_General_transform   *transform,
    VIO_Real                x,
    VIO_Real                y,
    VIO_Real                z,
    VIO_Real                *x_transformed,
    VIO_Real                *y_transformed,
    VIO_Real                *z_transformed );

 void  general_inverse_transform_point_in_trans_plane(
    VIO_General_transform   *transform,
    VIO_Real                x,
    VIO_Real                y,
    VIO_Real                z,
    VIO_Real                *x_transformed,
    VIO_Real                *y_transformed,
    VIO_Real                *z_transformed );



/* Main program */
char *prog_name;

int main(int argc, char *argv[])
{
   VIO_General_transform 
     transform, 
     *grid_transform_ptr, 
     forward_transform;
   VIO_Volume 
     target_vol,
     volume;
   volume_input_struct 
     input_info;
   VIO_Real
     voxel[VIO_MAX_DIMENSIONS],
     steps[VIO_MAX_DIMENSIONS],
     start[VIO_N_DIMENSIONS],
     target_steps[VIO_MAX_DIMENSIONS],
     wx,wy,wz, inv_x, inv_y, inv_z,
     def_values[VIO_MAX_DIMENSIONS];

   static int 
     clobber_flag = FALSE,
     verbose      = TRUE,
     debug        = FALSE;
   static char  
     *target_file;

   int 
     parse_flag,
     prog_count,
     sizes[VIO_MAX_DIMENSIONS],
     target_sizes[VIO_MAX_DIMENSIONS],
     xyzv[VIO_MAX_DIMENSIONS],
     target_xyzv[VIO_MAX_DIMENSIONS],
     index[VIO_MAX_DIMENSIONS],
     i,
     trans_count;
   VIO_progress_struct
     progress;


   static ArgvInfo argTable[] = {
     {"-like",       ARGV_STRING,   (char *) 0,     (char *) &target_file,
        "Specify target volume sampling information."},
     {"-no_clobber", ARGV_CONSTANT, (char *) FALSE, (char *) &clobber_flag,
        "Do not overwrite output file (default)."},
     {"-clobber",    ARGV_CONSTANT, (char *) TRUE,  (char *) &clobber_flag,
        "Overwrite output file."},
     {"-verbose",    ARGV_CONSTANT, (char *) TRUE,     (char *) &verbose,
        "Write messages indicating progress (default)"},
     {"-quiet",      ARGV_CONSTANT, (char *) FALSE,    (char *) &verbose,
        "Do not write log messages"},
     {"-debug",      ARGV_CONSTANT, (char *) TRUE,  (char *) &debug,
        "Print out debug info."},
     {NULL, ARGV_END, NULL, NULL, NULL}
   };


   prog_name = argv[0];
   target_file = malloc(1024);
   strcpy(target_file,"");

   /* Call ParseArgv to interpret all command line args (returns TRUE if error) */
   parse_flag = ParseArgv(&argc, argv, argTable, 0);

   /* Check remaining arguments */
   if (parse_flag || argc != 3) print_usage_and_exit(prog_name);

   /* Read in file that has a def field to invert */
   if (input_transform_file(argv[1], &transform) != OK) {
      (void) fprintf(stderr, "%s: Error reading transform file %s\n",
                     argv[0], argv[1]);
      exit(EXIT_FAILURE);
   }


   for(trans_count=0; trans_count<get_n_concated_transforms(&transform); trans_count++ ) {

     grid_transform_ptr = get_nth_general_transform(&transform, trans_count );
     
     if (grid_transform_ptr->type == GRID_TRANSFORM) {

       copy_general_transform(grid_transform_ptr,
                              &forward_transform);

       /* 
          this is the call that should be made
          with the latest version of internal_libvolume_io
        
          invert_general_transform(&forward_transform); */

       forward_transform.inverse_flag = !(forward_transform.inverse_flag);

       volume = grid_transform_ptr->displacement_volume;

       if (strlen(target_file)!=0) {
         if (debug) print ("Def field will be resampled like %s\n",target_file);
         
         if (!file_exists( target_file ) ) {
           (void) fprintf(stderr, "%s: Target file '%s' does not exist\n",
                          prog_name,target_file);
           exit(EXIT_FAILURE);
         }

         start_volume_input(target_file, 3, (char **)NULL, 
                            NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                            TRUE, &target_vol, 
                            (minc_input_options *)NULL,
                            &input_info);
         get_volume_XYZV_indices(volume, xyzv);
         get_volume_separations (volume, steps);
         get_volume_sizes       (volume, sizes);

         get_volume_XYZV_indices(target_vol, target_xyzv);
         get_volume_separations (target_vol, target_steps);
         get_volume_sizes       (target_vol, target_sizes);

         for(i=0; i<VIO_MAX_DIMENSIONS; i++) {
           index[i] = 0;
           voxel[i] = 0.0;
         }
         convert_voxel_to_world(target_vol, voxel, &start[VIO_X], &start[VIO_Y], &start[VIO_Z]);

         if( volume != (void *) NULL ){
           free_volume_data( volume );
         }

         for(i=VIO_X; i<=VIO_Z; i++) {
           steps[ xyzv[i] ] = target_steps[ target_xyzv[i] ] ;
           sizes[ xyzv[i] ] = target_sizes[ target_xyzv[i] ] ;
         }
         set_volume_separations(volume, steps);
         set_volume_sizes(      volume, sizes);
         set_volume_starts(volume, start);
         alloc_volume_data( volume );
       }

       get_volume_sizes(volume, sizes);
       get_volume_XYZV_indices(volume,xyzv);

       for(i=0; i<VIO_MAX_DIMENSIONS; i++){
         index[i] = 0;
       }

       if (verbose){
        initialize_progress_report(&progress, FALSE, 
                                   sizes[xyzv[VIO_X]]*sizes[xyzv[VIO_Y]]*sizes[xyzv[VIO_Z]]+1,
                                   "Inverting def field");
       }
       prog_count = 0;

       for(index[xyzv[VIO_X]]=0; index[xyzv[VIO_X]]<sizes[xyzv[VIO_X]]; index[xyzv[VIO_X]]++)
         for(index[xyzv[VIO_Y]]=0; index[xyzv[VIO_Y]]<sizes[xyzv[VIO_Y]]; index[xyzv[VIO_Y]]++)
           for(index[xyzv[VIO_Z]]=0; index[xyzv[VIO_Z]]<sizes[xyzv[VIO_Z]]; index[xyzv[VIO_Z]]++) {
             
             index[ xyzv[VIO_Z+1] ] = 0;
             for(i=0; i<VIO_MAX_DIMENSIONS; i++) voxel[i] = (VIO_Real)index[i];
       
             convert_voxel_to_world(volume, voxel, &wx, &wy, &wz);
             
             if (sizes[ xyzv[VIO_Z] ] ==1)
                general_inverse_transform_point_in_trans_plane(&forward_transform,
                                            wx, wy, wz,
                                            &inv_x, &inv_y, &inv_z);
             else
               grid_inverse_transform_point(&forward_transform,
                                            wx, wy, wz,
                                            &inv_x, &inv_y, &inv_z);
             def_values[VIO_X] = inv_x - wx;
             def_values[VIO_Y] = inv_y - wy;
             def_values[VIO_Z] = inv_z - wz;

             for(index[xyzv[VIO_Z+1]]=0; index[xyzv[VIO_Z+1]]<3; index[xyzv[VIO_Z+1]]++)
               set_volume_real_value(volume,
                                     index[0],index[1],index[2],index[3],index[4],
                                     def_values[ index[ xyzv[VIO_Z+1] ]]);

             prog_count++;
             if (verbose)
               update_progress_report(&progress, prog_count);
           }
       
       if (verbose)
         terminate_progress_report(&progress);

       delete_general_transform(&forward_transform);

       grid_transform_ptr->inverse_flag = !(grid_transform_ptr->inverse_flag);
       
     }

   }
   

   /* Write out the transform */
   if (output_transform_file(argv[2], NULL, &transform) != OK) {
      (void) fprintf(stderr, "%s: Error writing transform file %s\n",
                     argv[0], argv[2]);
      exit(EXIT_FAILURE);
   }

   exit(EXIT_SUCCESS);
}


void print_usage_and_exit(char *pname) {

  (void) fprintf(stderr, "This program is used to invert the internal representation\n");
  (void) fprintf(stderr, "of a GRID_TRANSFORM in order to speed up resampling.\n\n");
  (void) fprintf(stderr, "Usage: %s [options] <input.xfm> <result.xfm>\n",
                 pname);
  exit(EXIT_FAILURE);

}

