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
@MODIFIED   : Revision 1.1  1999-10-25 19:52:10  louis
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Extra_progs/reversedef.c,v 1.1 1999-10-25 19:52:10 louis Exp $";
#endif

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <internal_volume_io.h>
#include <minc_def.h>
#include <ParseArgv.h>

/* Constants */
#ifndef TRUE
#  define TRUE 1
#  define FALSE 0
#endif
#ifndef public
#  define public
#  define private static
#endif


void print_usage_and_exit(char *pname);

void get_volume_XYZV_indices(Volume data, int xyzv[])
{
  
  int 
    axis, i, vol_dims;
  char 
    **data_dim_names;
 
  vol_dims       = get_volume_n_dimensions(data);
  data_dim_names = get_volume_dimension_names(data);
 
  for_less(i,0,N_DIMENSIONS+1) xyzv[i] = -1;
  for_less(i,0,vol_dims) {
    if (convert_dim_name_to_spatial_axis(data_dim_names[i], &axis )) {
      xyzv[axis] = i; 
    } 
    else {     /* not a spatial axis */
      xyzv[Z+1] = i;
    }
  }

#ifdef HAVE_RECENT_VOLUME_IO
  delete_dimension_names(data, data_dim_names);
#else  
  delete_dimension_names(data_dim_names);
#endif
 
}

public  void  general_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed );

public  void  general_inverse_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed );



/* Main program */
char *prog_name;

int main(int argc, char *argv[])
{
   General_transform 
     transform, 
     *grid_transform_ptr, 
     forward_transform;
   Volume 
     target_vol,
     volume;
   volume_input_struct 
     input_info;
   Real
     voxel[MAX_DIMENSIONS],
     steps[MAX_DIMENSIONS],
     start[N_DIMENSIONS],
     target_steps[MAX_DIMENSIONS],
     wx,wy,wz, inv_x, inv_y, inv_z,
     def_values[MAX_DIMENSIONS];

   static int 
     clobber_flag = FALSE,
     verbose      = TRUE,
     debug        = FALSE;
   static char  
     *target_file;

   int 
     parse_flag,
     prog_count,
     sizes[MAX_DIMENSIONS],
     target_sizes[MAX_DIMENSIONS],
     xyzv[MAX_DIMENSIONS],
     target_xyzv[MAX_DIMENSIONS],
     index[MAX_DIMENSIONS],
     i,
     trans_count;
   progress_struct
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
   ALLOC(target_file,1024);
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


   for_less(trans_count,0, get_n_concated_transforms(&transform) ) {

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

	 for_less(i,0,MAX_DIMENSIONS) {
	   index[i] = 0;
	   voxel[i] = 0.0;
	 }
	 convert_voxel_to_world(target_vol, voxel, &start[X], &start[Y], &start[Z]);

	 if( VOXEL_DATA(volume) != (void *) NULL )
	   free_volume_data( volume );

	 for_inclusive(i,X,Z) {
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

       for_less(i,0,MAX_DIMENSIONS) index[i] = 0;

       if (verbose)
	initialize_progress_report(&progress, FALSE, 
				   sizes[xyzv[X]]*sizes[xyzv[Y]]*sizes[xyzv[Z]]+1,
				   "Inverting def field");
       prog_count = 0;

       for_less( index[ xyzv[X] ], 0, sizes[ xyzv[X] ])
	 for_less( index[ xyzv[Y] ], 0, sizes[ xyzv[Y] ])
	   for_less( index[ xyzv[Z] ], 0, sizes[ xyzv[Z] ]) {
	     
	     index[ xyzv[Z+1] ] = 0;
	     for_less(i,0,MAX_DIMENSIONS) voxel[i] = (Real)index[i];
       
	     convert_voxel_to_world(volume, voxel, &wx, &wy, &wz);
	     
	     if (sizes[ xyzv[Z] ] ==1)
	        general_inverse_transform_point_in_trans_plane(&forward_transform,
					    wx, wy, wz,
					    &inv_x, &inv_y, &inv_z);
	     else
	       grid_inverse_transform_point(&forward_transform,
					    wx, wy, wz,
					    &inv_x, &inv_y, &inv_z);
	     def_values[X] = inv_x - wx;
	     def_values[Y] = inv_y - wy;
	     def_values[Z] = inv_z - wz;

	     for_less(index[ xyzv[Z+1] ],0,3)
	       set_volume_real_value(volume,
				     index[0],index[1],index[2],index[3],index[4],
				     def_values[ index[ xyzv[Z+1] ]]);

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

