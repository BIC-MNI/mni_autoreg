/* ----------------------------- MNI Header -----------------------------------
@NAME       : init_lattice
                set up the START position and element COUNT for the
                sampling lattice used to do the volume comparisons.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
              m1,m2:
                two mask volumes for data (already in memory).
              globals:
                a global data structure containing info from the command line,
                including the input parameters to be optimized, the input matrix,
                and a plethora of flags!
@OUTPUT     : globals->start and globals->count
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: 
@GLOBALS    : 
@CALLS      : 
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

@CREATED    : Wed Jun  9 12:56:08 EST 1993 LC
@MODIFIED   :  $Log: init_lattice.c,v $
@MODIFIED   :  Revision 96.12  2009-03-13 19:51:31  claude
@MODIFIED   :  fixed bug in offsets for minctracc and free memory upon exit
@MODIFIED   :
@MODIFIED   :  Revision 96.11  2006/11/30 09:07:33  rotor
@MODIFIED   :   * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   :  Revision 96.9  2005/07/20 20:45:51  rotor
@MODIFIED   :      * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :      * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :      * Removed old VOLUME_IO cruft #defines
@MODIFIED   :      * Fixed up all Makefile.am's in subdirs
@MODIFIED   :      * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :      * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   :  Revision 96.8  2004/02/12 06:04:52  rotor
@MODIFIED   :   * removed /static defs
@MODIFIED   :
@MODIFIED   :  Revision 96.7  2004/02/04 20:44:42  lenezet
@MODIFIED   :  *** empty log message ***
@MODIFIED   :
@MODIFIED   :  Revision 96.6  2003/02/26 00:50:25  lenezet
@MODIFIED   :  for 2D : now computes all 3 coordinates for the "start" (to take into account the slice position)
@MODIFIED   :
@MODIFIED   :  Revision 96.5  2002/12/13 21:18:45  lenezet
@MODIFIED   :  nonlinear in 2D has changed. The option -2D-non-lin is no more necessary. The grid transform has been adapted to feet on the target volume whatever is size. The Optimization is done on the dimensions for which "count" is greater than 1.
@MODIFIED   :
@MODIFIED   :  Revision 96.4  2002/11/20 21:39:18  lenezet
@MODIFIED   :
@MODIFIED   :  Fix the code to take in consideration the direction cosines especially in the grid transform.
@MODIFIED   :  Add an option to choose the maximum expected deformation magnitude.
@MODIFIED   :
@MODIFIED   :  Revision 96.3  2002/03/26 14:15:47  stever
@MODIFIED   :  Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   :  Revision 96.2  2000/02/07 19:33:07  stever
@MODIFIED   :  replaced HAVE_RECENT_VOLUME_IO with more specific feature tests.
@MODIFIED   :
@MODIFIED   :  Revision 96.1  1999/10/25 19:59:16  louis
@MODIFIED   :  final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.0  1996/08/21  18:22:15  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.5  1996/08/12  14:16:15  louis
 * Release of MNI_AutoReg version 1.0
 *
 * Revision 1.10  1996/08/12  14:16:12  louis
 * Pre-release
 *
 * Revision 1.9  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.8  94/06/02  20:16:00  louis
 * made modifications to allow deformations to be calulated in 2D on slices. 
 * changes had to be made in set_up_lattice, init_lattice when defining
 * the special case of a single slice....
 * Build_default_deformation_field also had to reflect these changes.
 * do_non-linear-optimization also had to check if one of dimensions had
 * a single element.
 * All these changes were made, and slightly tested.  Another type of
 * deformation strategy will be necessary (to replace the deformation 
 * perpendicular to the surface, since it does not work well).
 * 


made change to init lattice to not change start when there is only 1 slice.

 * Revision 1.7  94/04/06  11:48:35  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.6  94/02/21  16:35:29  louis
 * version before feb 22 changes
 * 
 * Revision 1.5  93/11/15  16:25:43  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/VIO_Volume/init_lattice.c,v 96.12 2009-03-13 19:51:31 claude Exp $";
#endif

#include <config.h>
#include <volume_io.h>

#include "matrix_basics.h"
#include "make_rots.h"

#include "constants.h"
#include "minctracc_arg_data.h"

#include "local_macros.h"

#include <Proglib.h>

extern Arg_Data *main_args;

        /* prototype from interpolation.c */
int point_not_masked(VIO_Volume volume, 
                            VIO_Real wx, VIO_Real wy, VIO_Real wz);




void get_volume_XYZV_indices(VIO_Volume data, int xyzv[])
{
  
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


/* 
   find the largest lattice array, aligned with the existing data volume,
   that will fit inside the volumetric space define by the data volume,
   given the user spacing.
   
   return the lattice info in start, count, step and directions. 

   NOTE: the START is in the 'volume world coordinate system', ie, the
   minc coordinate system, where the start is the offset along the
   direction cosines for each axis, and NOT THE ABSOLUTE WORLD COORD!

   count is constrained to be positive,
   step is constrained to have the same sign as the data volume.

   if one of the dimensions of the input volume is of length 1, then the
   output lattice will be of the same length.

*/

void set_up_lattice(VIO_Volume data,       /* in: volume  */
                           double *user_step, /* in: user requested spacing for lattice */
                           double *start,     /* out: starting position of lattice in volume-dircos coords*/
                           double *wstart,     /* out:world starting position of lattice */
                           int    *count,     /* out:number of steps in each direction */
                           double *step,      /* out:step size in each direction */
                           VectorR directions[])/* out: vector directions for each index*/

                                /* note that user_step, start, count, step and directions
                                   are in x,y,z order*/

{
  int 
    xyzv[VIO_MAX_DIMENSIONS],        /* the indices for the x,y,z and vect directions in
                                   the data volume. */
    sizes[VIO_MAX_DIMENSIONS], 
    i,j,verbose;
  VIO_Real 
    sign,
    direction[VIO_N_DIMENSIONS],
    starts[VIO_MAX_DIMENSIONS],
    start_voxel[VIO_MAX_DIMENSIONS],
    start_world[VIO_MAX_DIMENSIONS],
    vect_voxel[VIO_MAX_DIMENSIONS],
    tmp_point[VIO_N_DIMENSIONS],
    num_steps,
    offset[VIO_MAX_DIMENSIONS],
    separations[VIO_MAX_DIMENSIONS];
  VIO_BOOL 
    debug; 

  debug  = main_args->flags.debug;
  verbose= main_args->flags.verbose;
  
  for(i=0; i<VIO_MAX_DIMENSIONS; i++)
    {
      sizes[i]=0;
      separations[i]=0.0;
    }

                                /* get the volume sizes and voxel spacing */
  get_volume_sizes(data, sizes );
  get_volume_separations(data, separations );
                                /* set up the x_ y_ and z_ind indices to be
                                   able to access the volume data below.     */
  get_volume_XYZV_indices(data, xyzv);
  if (debug) {
    print ("In set_up_lattice, xyzv[axes] = %d, %d, %d, %d\n",
           xyzv[VIO_X],xyzv[VIO_Y],xyzv[VIO_Z],xyzv[VIO_Z+1]);
  }
    
  for(i=0; i<3; i++)                /* copy the requested step values for X, Y, Z */
    step[i] = user_step[i]; 
  
  
  if (debug && verbose>1) {
    print ("In set_up_lattice, data volume is (in x y z order):\n");
    print ("sizes: %7d %7d %7d\n",
           xyzv[VIO_X]>-1 ? sizes[xyzv[VIO_X]]: 0,
           xyzv[VIO_Y]>-1 ? sizes[xyzv[VIO_Y]]: 0,
           xyzv[VIO_Z]>-1 ? sizes[xyzv[VIO_Z]]: 0);
    print ("steps: %7.2f %7.2f %7.2f\n",
           xyzv[VIO_X]>-1 ? separations[xyzv[VIO_X]] : 0.0,
           xyzv[VIO_Y]>-1 ? separations[xyzv[VIO_Y]] : 0.0,
           xyzv[VIO_Z]>-1 ? separations[xyzv[VIO_Z]] : 0.0);
    for(i=0; i<VIO_MAX_DIMENSIONS; i++) start_voxel[i] = 0.0;
    convert_voxel_to_world(data, start_voxel,
                           &starts[VIO_X], &starts[VIO_Y], &starts[VIO_Z]);
    print ("start (v->w), 0,0,0 ->(x,y,z)    : %7.2f %7.2f %7.2f\n",starts[VIO_X],starts[VIO_Y],starts[VIO_Z]);
 
    convert_world_to_voxel(data,starts[VIO_X],starts[VIO_Y],starts[VIO_Z],start_voxel);
    print ("start voxel (w->v): %7.2f %7.2f %7.2f\n",start_voxel[xyzv[VIO_X]],start_voxel[xyzv[VIO_Y]],start_voxel[xyzv[VIO_Z]]);

    get_volume_starts(data, starts);
    print ("start (from data struct, !world) : %7.2f %7.2f %7.2f\n",starts[xyzv[VIO_X]],starts[xyzv[VIO_Y]],starts[xyzv[VIO_Z]]);

    
    print ("spatial axes %d %d %d\n",data->spatial_axes[0],data->spatial_axes[1],data->spatial_axes[2]);
    print ("separations %7.2f  %7.2f  %7.2f \n",data->separations[0],data->separations[1],data->separations[2]);
    print ("directions cosines \n %7.2f  %7.2f  %7.2f \n %7.2f  %7.2f  %7.2f \n %7.2f  %7.2f  %7.2f \n",data->direction_cosines[0][VIO_X],data->direction_cosines[0][VIO_Y],data->direction_cosines[0][VIO_Z],data->direction_cosines[1][VIO_X],data->direction_cosines[1][VIO_Y],data->direction_cosines[1][VIO_Z],data->direction_cosines[2][VIO_X],data->direction_cosines[2][VIO_Y],data->direction_cosines[2][VIO_Z]);
         

  }
                                /* given the input step size,
                                   figure out the count and starting offset */



  for(i=0; i<VIO_N_DIMENSIONS; i++) {        
    
    count[i] = 1;
    offset[xyzv[i]] = 0.0;
 
    if (xyzv[i] >= 0 && sizes[ xyzv[i] ] > 1) {
                                /* force step to have the same SIGN as the 
                                   volume voxel spacing. */
      step[i] = fabs( step[i] );
      if (separations[xyzv[i]] < 0.0)  step[i] *= (-1.0);
      
      num_steps = separations[xyzv[i]] * sizes[xyzv[i]] / step[i];
      num_steps = fabs(num_steps);
      count[i] = (int)floor(num_steps + 0.5);
      if (count[i] == 0) count[i] = 1;
    
                                /* this is the offset for the start of the
                                   lattice from the corner of the volume 
                                   in world distance mm                   */
      offset[xyzv[i]] = 0.5 * (separations[xyzv[i]]*sizes[xyzv[i]] - 
                              step[i]*count[i]);

    }
  }

                                /* get the voxel position of the lattice start,
                                   in voxel coordinates of the data volume  */
    
  for(i=0; i<VIO_MAX_DIMENSIONS; i++) {
    start_voxel[i] = 0.0;
    start[i] = 0.0;
    starts[i]=0.0;
  }
  
  

  
  get_volume_starts(data, starts); /* get the volume-world start position (in dir_cos coords */

  for(i=0; i<VIO_N_DIMENSIONS; i++) {
    
    if (separations[xyzv[i]] > 0)
      sign = 1.0;
    else
      sign = -1.0;
    
    if (xyzv[i]>=0 && sizes[xyzv[i]]>1 ) {
      start_voxel[xyzv[i]] = sign*((-0.5)  /* to get to the edge of the voxel,
                                            since the voxel's  coordinates is 
                                            at its center */
                             
                                 + (offset[xyzv[i]]/separations[xyzv[i]]) 
                                         /* the offset to edge of lattice */
                                 + (step[i]/2)/separations[xyzv[i]]);
                                            /* the offset to the center of the 
                                            lattice voxel */

      start_world[xyzv[i]] =  starts[xyzv[i]] - separations[xyzv[i]]/2.0 - sign*offset[xyzv[i]] 
        + step[i]/2.0;
    }
    
    if(xyzv[i]>=0 && sizes[xyzv[i]]==1 ){
      start_world[xyzv[i]]=starts[xyzv[i]];
    }
  }
    
  /* get the absolute world starting coordinates of the origin of the lattice */
  convert_voxel_to_world(data, start_voxel, &wstart[VIO_X], &wstart[VIO_Y], &wstart[VIO_Z]);

  for(i=0; i<VIO_N_DIMENSIONS; i++)
    start[xyzv[i]] = start_world[xyzv[i]];
  
                                /* get direction vectors*/
  for(i=0; i<VIO_N_DIMENSIONS; i++) {

    get_volume_direction_cosine(data, xyzv[i], direction);
    
    fill_Vector(directions[i], 
                direction[VIO_X], direction[VIO_Y], direction[VIO_Z]);
  }
 

  if (debug && verbose>1) {
    
    print ("       for lattice volume:\n");
    print ("sizes: %7d %7d %7d\n",count[VIO_X],count[VIO_Y],count[VIO_Z]);
    print ("steps: %7.2f %7.2f %7.2f\n",step[VIO_X],step[VIO_Y],step[VIO_Z]);
    print ("starts:%7.2f %7.2f %7.2f %7.2f %7.2f <- in volume order\n",starts[0],starts[1],starts[2],starts[3],starts[4]);
    print ("wstart:%7.2f %7.2f %7.2f\n",wstart[VIO_X],wstart[VIO_Y],wstart[VIO_Z]);
    print ("dir_x: %7.2f %7.2f %7.2f\n",directions[VIO_X].coords[VIO_X],
                                        directions[VIO_X].coords[VIO_Y],
                                        directions[VIO_X].coords[VIO_Z]);
    print ("dir_y: %7.2f %7.2f %7.2f\n",directions[VIO_Y].coords[VIO_X],
                                        directions[VIO_Y].coords[VIO_Y],
                                        directions[VIO_Y].coords[VIO_Z]);
    print ("dir_z: %7.2f %7.2f %7.2f\n",directions[VIO_Z].coords[VIO_X],
                                        directions[VIO_Z].coords[VIO_Y],
                                        directions[VIO_Z].coords[VIO_Z]);
    print ("leaving set_up_lattice()\n\n");
  }
 

}

/* 
   in this procedure, the smallest (in number of samples) 3D lattice
   is defined that complete covers either d1 or d2 (taking into account
   the mask volumes m1 and m2).  The lattice spacing is defined by the 
   globals->step variable, where the steps are stored in X,Y,Z order.

   The procedure returns the globals->start, globals->count and 
   globals->step that define the lattice (each in X Y Z order).  

   globals->smallest_vol ==1 or ==2, indicating on which volume space
   the lattice is defined.  (This volume always ==2 when NONLIN 
   transformations are optimized.

*/


void init_lattice(VIO_Volume d1,
                         VIO_Volume d2,
                         VIO_Volume m1,
                         VIO_Volume m2, 
                         Arg_Data *globals)
{

  VectorR
    vector_step;

  PointR 
    starting_position1,
    starting_position2,
    slice,
    row,
    col,
    voxel;

  double
    tx,ty,tz;
  int
    i,r,c,s;

  int 
    vol1,vol2,
    min1_row,   max1_row,
    min1_col,   max1_col,
    min1_slice, max1_slice,
    min2_row,   max2_row,
    min2_col,   max2_col,
    min2_slice, max2_slice;

  VIO_Real
    abs_step,
    factor,
    true_value;
  double
    tmp_threshold;

  int
    count1[VIO_MAX_DIMENSIONS], count2[VIO_MAX_DIMENSIONS];
  double 
    wstart1[VIO_MAX_DIMENSIONS], wstart2[VIO_MAX_DIMENSIONS],
    start1[VIO_MAX_DIMENSIONS], start2[VIO_MAX_DIMENSIONS],
    step1[VIO_MAX_DIMENSIONS],  step2[VIO_MAX_DIMENSIONS];
  VectorR
    scaled_directions1[VIO_MAX_DIMENSIONS],
    scaled_directions2[VIO_MAX_DIMENSIONS],
    directions1[VIO_MAX_DIMENSIONS],
    directions2[VIO_MAX_DIMENSIONS];  

  int voxels_found;

  

  if (globals->flags.debug && globals->flags.verbose>1)
    print ("\n***** entering init_lattice\n");
                                /* build default sampling lattice info
                                   on first data set (d1)               */
  
  /* get start in volume order and wstart,count,step and directions in XYZ order */
  set_up_lattice(d1, globals->step,
                 start1, wstart1, count1, step1, directions1);

  for(i=0; i<3; i++) {
    SCALE_VECTOR( scaled_directions1[i], directions1[i], step1[i]);
  }

  if (globals->flags.debug && globals->flags.verbose>1) {
    print ("in init_lattice: for the source data set, the new lattice is:\n");
    print ("start = %8.2f %8.2f %8.2f %8.2f in vol order\n",start1[0],start1[1],start1[2],start1[3]);
    print ("wstart= %8.2f %8.2f %8.2f \n",wstart1[0],wstart1[1],wstart1[2]);
    print ("count = %8d %8d %8d \n",count1[0],count1[1],count1[2]);
    print ("step  = %8.2f %8.2f %8.2f \n",step1[0],step1[1],step1[2]);
    print ("thres  = %f %f\n",globals->threshold[0],globals->threshold[1]);

    for(i=0; i<3; i++)
      print ("direct= %8.2f %8.2f %8.2f \n",
             Point_x(directions1[i]),
             Point_y(directions1[i]),
             Point_z(directions1[i]));
    
  }


  fill_Point( starting_position1, wstart1[0], wstart1[1], wstart1[2]);
  
                                /* init min1 max1 values */
  min1_col = count1[VIO_X]; max1_col = 0;
  min1_row = count1[VIO_Y]; max1_row = 0;
  min1_slice=count1[VIO_Z]; max1_slice = 0;
  voxels_found = FALSE;

  for(s=0; s<count1[VIO_Z]; s++) {

    SCALE_VECTOR( vector_step, scaled_directions1[VIO_Z], s);
    ADD_POINT_VECTOR( slice, starting_position1, vector_step );

    for(r=0; r<count1[VIO_Y]; r++) {

      SCALE_VECTOR( vector_step, scaled_directions1[VIO_Y], r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for(c=0; c<count1[VIO_X]; c++) {

        convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
        
        fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
        
        if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {        

          if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &true_value )) {
            
            if (true_value > globals->threshold[0]) {

              voxels_found = TRUE;

              if (r > max1_row) max1_row = r; 
              if (r < min1_row) min1_row = r;
              if (c > max1_col) max1_col = c; 
              if (c < min1_col) min1_col = c;
              if (s > max1_slice) max1_slice = s; 
              if (s < min1_slice) min1_slice = s;
            }
          }    
          else {
            if (globals->flags.debug && globals->flags.verbose>3) {
              print ("%3d %3d %3d : %12.5f %12.5f %12.5f -> %12.5f\n",c,r,s,
                      Point_x(col), Point_y(col), Point_z(col),true_value);
            }
          }
                

        } 


        ADD_POINT_VECTOR( col, col, scaled_directions1[VIO_X] );
        
      }
    }
  }

  vol1 = (max1_row - min1_row + 1) *
         (max1_col - min1_col + 1) *
         (max1_slice - min1_slice + 1);

  if (globals->flags.debug && globals->flags.verbose>1) 
    print ("volume =  %d\n",vol1);

  if (globals->flags.debug || (max1_row < min1_row) || (max1_col < min1_col) || 
      (max1_slice < min1_slice)) {

    print ("slice lim %d %d\n",min1_slice, max1_slice);
    print ("row lim   %d %d\n",min1_row, max1_row);
    print ("col lim   %d %d\n",min1_col, max1_col);
    print ("thresh = %10.5f %10.5f\n", globals->threshold[0],globals->threshold[1]);

    if (!voxels_found) {
       print ("No voxels were found in volume 1 with value above threshold (%f).\n",
              globals->threshold[0]);
       print_error_and_line_num("%s", __FILE__, __LINE__,"Cannot calculate size of volume 1\n.");
    }

  }



                                /* build default sampling lattice info
                                   on second data set (d2)               */

  set_up_lattice(d2, globals->step,
                 start2, wstart2, count2, step2, directions2);

  for(i=0; i<3; i++) {
    SCALE_VECTOR( scaled_directions2[i], directions2[i], step2[i]);
  }

  if (globals->flags.debug && globals->flags.verbose>1) {
    print ("in init_lattice: for the target data set\n");
    print ("start = %8.2f %8.2f %8.2f in vol order\n",start2[0],start2[1],start2[2]);
    print ("wstart= %8.2f %8.2f %8.2f \n",wstart2[0],wstart2[1],wstart2[2]);
    print ("count = %8d %8d %8d \n",count2[0],count2[1],count2[2]);
    print ("step  = %8.2f %8.2f %8.2f \n",step2[0],step2[1],step2[2]);
    print ("thresh = %10.5f %10.5f\n", globals->threshold[0],globals->threshold[1]);
    for(i=0; i<3; i++)
      print ("direct= %8.2f %8.2f %8.2f \n",
             Point_x(directions2[i]),
             Point_y(directions2[i]),
             Point_z(directions2[i]));

  }

  fill_Point( starting_position2, wstart2[0], wstart2[1], wstart2[2]);
  
                                /* init min2 max2 values */
  min2_col = count2[VIO_X]; max2_row = 0;
  min2_row = count2[VIO_Y]; max2_col = 0;
  min2_slice=count2[VIO_Z]; max2_slice = 0;

  voxels_found = FALSE;

  for(s=0; s<count2[VIO_Z]; s++) {

    SCALE_VECTOR( vector_step, scaled_directions2[VIO_Z], s);
    ADD_POINT_VECTOR( slice, starting_position2, vector_step );

    for(r=0; r<count2[VIO_Y]; r++) {

      SCALE_VECTOR( vector_step, scaled_directions2[VIO_Y], r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for(c=0; c<count2[VIO_X]; c++) {

        convert_3D_world_to_voxel(d2, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
        
        fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
        
        if (point_not_masked(m2, Point_x(col), Point_y(col), Point_z(col))) {                        /* should be fill_value  */
          
          if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &true_value )) {
            
            if (true_value > globals->threshold[1]) {

              voxels_found = TRUE;

              if (r > max2_row) max2_row = r; 
              if (r < min2_row) min2_row = r;
              if (c > max2_col) max2_col = c; 
              if (c < min2_col) min2_col = c;
              if (s > max2_slice) max2_slice = s; 
              if (s < min2_slice) min2_slice = s;
            }
          }
        } 

        ADD_POINT_VECTOR( col, col, scaled_directions2[VIO_X] );
        
      }
    }
  }

  if ((max2_row < min2_row) || (max2_col < min2_col) || 
      (max2_slice < min2_slice)) {
    print ("slice lim %d %d\n",min2_slice, max2_slice);
    print ("row lim   %d %d\n",min2_row, max2_row);
    print ("col lim   %d %d\n",min2_col, max2_col);
    print ("thresh = %10.5f %10.5f\n", globals->threshold[0],globals->threshold[1]);
    if (!voxels_found) {
      print ("No voxels were found in volume 2 with value above threshold (%f).\n",
             globals->threshold[1]);
      print_error_and_line_num("%s", __FILE__, __LINE__,"Cannot calculate size of volume 2\n." );
    }
  }

  vol2 = (max2_row - min2_row + 1) *
         (max2_col - min2_col + 1) *
         (max2_slice - min2_slice + 1);

  if (globals->flags.debug && globals->flags.verbose>1) 
    print ("volume =  %d\n",vol2);

                                /* set up the lattice on the source only when
                                    vol1<=vol2  or
                                    source is forced, and we are not running NONLIN 
                                   otherwise,
                                    the lattice should be set up on the target*/
  if ( !(globals->trans_info.transform_type==TRANS_NONLIN) && 
      ((vol1<=vol2)  || (globals->force_lattice==1)) && 
      !(globals->force_lattice==2) ) {
    globals->smallest_vol = 1;

    globals->count[VIO_X] = max1_col - min1_col + 1;
    globals->count[VIO_Y] = max1_row - min1_row + 1;
    globals->count[VIO_Z] = max1_slice - min1_slice + 1;

    fill_Point( starting_position1, wstart1[0], wstart1[1], wstart1[2]);
    SCALE_VECTOR( vector_step, scaled_directions1[VIO_Z], min1_slice);
    ADD_POINT_VECTOR( slice, starting_position1, vector_step );
    SCALE_VECTOR( vector_step, scaled_directions1[VIO_Y], min1_row);
    ADD_POINT_VECTOR( row, slice, vector_step );
    SCALE_VECTOR( vector_step, scaled_directions1[VIO_X], min1_col);
    ADD_POINT_VECTOR( col, row, vector_step );

    globals->start[VIO_X] = Point_x(col);
    globals->start[VIO_Y] = Point_y(col);
    globals->start[VIO_Z] = Point_z(col);

    

    for(i=0; i<3; i++) {

      globals->step[i] = step1[i]; /* ensure that global step has same sign as data voxel */

      abs_step = fabs(globals->step[i]);

      Point_x(globals->directions[i]) = Point_x(directions1[i]) * abs_step;
      Point_y(globals->directions[i]) = Point_y(directions1[i]) * abs_step;
      Point_z(globals->directions[i]) = Point_z(directions1[i]) * abs_step;
    }
  }
  else {
                                /* The lattice is always defined on the target volume
                                   when computing non-linear transformations.          */

    globals->smallest_vol = 2;
    globals->count[VIO_X] = max2_col - min2_col + 1;
    globals->count[VIO_Y] = max2_row - min2_row + 1;
    globals->count[VIO_Z] = max2_slice - min2_slice + 1;

    fill_Point( starting_position2, wstart2[0], wstart2[1], wstart2[2]);
    SCALE_VECTOR( vector_step, scaled_directions2[VIO_Z], min2_slice);
    ADD_POINT_VECTOR( slice, starting_position2, vector_step );
    SCALE_VECTOR( vector_step, scaled_directions2[VIO_Y], min2_row);
    ADD_POINT_VECTOR( row, slice, vector_step );
    SCALE_VECTOR( vector_step, scaled_directions2[VIO_X], min2_col);
    ADD_POINT_VECTOR( col, row, vector_step );

    globals->start[VIO_X] = Point_x(col);
    globals->start[VIO_Y] = Point_y(col);
    globals->start[VIO_Z] = Point_z(col);

    for(i=0; i<3; i++) {

      globals->step[i] = step2[i]; /* ensure that global step has same sign as data voxel */

      abs_step = fabs(globals->step[i]);
      
      Point_x(globals->directions[i]) = Point_x(directions2[i]) * abs_step;
      Point_y(globals->directions[i]) = Point_y(directions2[i]) * abs_step;
      Point_z(globals->directions[i]) = Point_z(directions2[i]) * abs_step;
    }
    
    tmp_threshold = globals->threshold[0];
    globals->threshold[0] = globals->threshold[1];
    globals->threshold[1] = tmp_threshold;
  }

  if (globals->flags.debug && globals->flags.verbose>1)
    print ("***** leaving init_lattice\n");

}
