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
@MODIFIED   :  Revision 1.10  1996-08-12 14:16:12  louis
@MODIFIED   :  Pre-release
@MODIFIED   :
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Volume/init_lattice.c,v 1.10 1996-08-12 14:16:12 louis Exp $";
#endif

#include <config.h>
#include <volume_io.h>

#include "matrix_basics.h"
#include "make_rots.h"

#include "constants.h"
#include "arg_data.h"

#include "local_macros.h"

#include <print_error.h>

extern Arg_Data main_args;

	/* prototype from interpolation.c */
public int point_not_masked(Volume volume, 
			    Real wx, Real wy, Real wz);




public void get_volume_XYZV_indices(Volume data, int xyzv[])
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


/* 
   find the largest lattice array, aligned with the existing data volume,
   that will fit inside the volumetric space define by the data volume,
   given the user spacing.
   
   return the lattice info in start, count, step and directions. 

   count is constrained to be positive,
   step is constrained to have the same sign as the data volume.

   if one of the dimensions of the input volume is of length 1, then the
   output lattice will be of the same length.

*/

public void set_up_lattice(Volume data,       /* in: volume  */
			   double *user_step, /* in: user requested spacing for lattice */
			   double *start,     /* out:world starting position of lattice */
			   int    *count,     /* out:number of steps in each direction */
			   double *step,      /* out:step size in each direction */
			   VectorR directions[])/* out: vector directions for each index*/

				/* note that user_step, start, count, step and directions
				   are in x,y,z order*/

{
  int 
    xyzv[MAX_DIMENSIONS],	/* the indices for the x,y,z and vect directions in
				   the data volume. */
    sizes[MAX_DIMENSIONS], 
    i,j,verbose;
  Real 
    sign,
    starts[N_DIMENSIONS],
    start_voxel[MAX_DIMENSIONS],
    vect_voxel[MAX_DIMENSIONS],
    tmp_point[N_DIMENSIONS],
    num_steps,
    offset[MAX_DIMENSIONS],
    separations[MAX_DIMENSIONS];
  BOOLEAN 
    debug; 

  debug  = main_args.flags.debug;
  verbose= main_args.flags.verbose;
  
				/* get the volume sizes and voxel spacing */
  get_volume_sizes(data, sizes );
  get_volume_separations(data, separations );
				/* set up the x_ y_ and z_ind indices to be
				   able to access the volume data below.     */
  get_volume_XYZV_indices(data, xyzv);
  if (debug) {
    print ("In set_up_lattice, xyzv[axes] = %d, %d, %d, %d\n",
	   xyzv[X],xyzv[Y],xyzv[Z],xyzv[Z+1]);
  }
    
  for_less(i,0,3)		/* copy the requested step values for X, Y, Z */
    step[i] = user_step[i]; 
  
  
  if (debug && verbose>1) {
    print ("In set_up_lattice, data volume is (in x y z order):\n");
    print ("sizes: %7d %7d %7d\n",
	   xyzv[X]>-1 ? sizes[xyzv[X]]: 0,
	   xyzv[Y]>-1 ? sizes[xyzv[Y]]: 0,
	   xyzv[Z]>-1 ? sizes[xyzv[Z]]: 0);
    print ("steps: %7.2f %7.2f %7.2f\n",
	   xyzv[X]>-1 ? separations[xyzv[X]] : 0.0,
	   xyzv[Y]>-1 ? separations[xyzv[Y]] : 0.0,
	   xyzv[Z]>-1 ? separations[xyzv[Z]] : 0.0);
    for_less(i,0,MAX_DIMENSIONS) start_voxel[i] = 0.0;
    convert_voxel_to_world(data, start_voxel,
			   &starts[X], &starts[Y], &starts[Z]);
    print ("start: %7.2f %7.2f %7.2f\n",starts[X],starts[Y],starts[Z]); 
  }
				/* given the input step size,
				   figure out the count and starting offset */
  for_less( i, 0, N_DIMENSIONS) {	
    
    count[i] = 1;
    offset[xyzv[i]] = 0.0;
    
    if (xyzv[i] >= 0 && sizes[ xyzv[i] ] > 1) {
				/* force step to have the same SIGN as the 
				   volume voxel spacing. */
      step[i] = ABS( step[i] );
      if (separations[xyzv[i]] < 0.0)  step[i] *= (-1.0);
      
      num_steps = separations[xyzv[i]] * sizes[xyzv[i]] / step[i];
      num_steps = ABS(num_steps);
      
      count[i] = ROUND(num_steps);
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
    
  for_less(i,0,MAX_DIMENSIONS) start_voxel[i] = 0.0;
  
  for_less(i,0,N_DIMENSIONS) {
    
    if (separations[i] > 0)
      sign = 1.0;
    else
      sign = -1.0;
    
    if (xyzv[i]>=0 && sizes[xyzv[i]]>1) {
      start_voxel[xyzv[i]]= sign*((-0.5)  /* to get to the edge of the voxel,
                                            since the voxel's  coordinates is 
					    at its center */
			     
				 + (offset[xyzv[i]]/separations[xyzv[i]]) 
				         /* the offset to edge of lattice */
				 + step[i]/(2*separations[xyzv[i]]));
			   	         /* the offset to the center of the 
					    lattice voxel */
    }
    
  }
    
				/* get the world start.  */
  convert_voxel_to_world(data, start_voxel,
			 &tmp_point[0], &tmp_point[1], &tmp_point[2]);
  
  for_less(i,0,N_DIMENSIONS)
    start[i] = tmp_point[i];
  
				/* build direction vectors*/
  for_less(i,0,N_DIMENSIONS) {
    
    for_less(j,0,N_DIMENSIONS)
      vect_voxel[xyzv[j]] = start_voxel[xyzv[j]];

    vect_voxel[xyzv[i]] += step[i]/separations[xyzv[i]];
    
    convert_voxel_to_world(data, vect_voxel,
			   &tmp_point[0], &tmp_point[1], &tmp_point[2]);
    
    fill_Vector(directions[i], 
		tmp_point[0]-start[0], tmp_point[1]-start[1], tmp_point[2]-start[2]);
  }
  
  if (debug && verbose>1) {
    
    print ("       for lattice volume:\n");
    print ("sizes: %7d %7d %7d\n",count[X],count[Y],count[Z]);
    print ("steps: %7.2f %7.2f %7.2f\n",step[X],step[Y],step[Z]);
    print ("start: %7.2f %7.2f %7.2f\n",start[X],start[Y],start[Z]);
    print ("dir_x: %7.2f %7.2f %7.2f\n",directions[X].coords[X],
	                                directions[X].coords[Y],
	                                directions[X].coords[Z]);
    print ("dir_y: %7.2f %7.2f %7.2f\n",directions[Y].coords[X],
	                                directions[Y].coords[Y],
	                                directions[Y].coords[Z]);
    print ("dir_z: %7.2f %7.2f %7.2f\n",directions[Z].coords[X],
	                                directions[Z].coords[Y],
	                                directions[Z].coords[Z]);
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


public void init_lattice(Volume d1,
			 Volume d2,
			 Volume m1,
			 Volume m2, 
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

  Real
    true_value;
  double
    tmp_threshold;

  int
    count1[MAX_DIMENSIONS], count2[MAX_DIMENSIONS];
  double 
    start1[MAX_DIMENSIONS], start2[MAX_DIMENSIONS],
    step1[MAX_DIMENSIONS],  step2[MAX_DIMENSIONS];
  VectorR
    directions1[MAX_DIMENSIONS],
    directions2[MAX_DIMENSIONS];  

  if (globals->flags.debug && globals->flags.verbose>1)
    print ("\n***** entering init_lattice\n");
				/* build default sampling lattice info
				   on first data set (d1)               */
  set_up_lattice(d1, globals->step,
		 start1, count1, step1, directions1);

  if (globals->flags.debug && globals->flags.verbose>1) {
    print ("in init_lattice: for the source data set\n");
    print ("start = %8.2f %8.2f %8.2f \n",start1[0],start1[1],start1[2]);
    print ("count = %8d %8d %8d \n",count1[0],count1[1],count1[2]);
    print ("step  = %8.2f %8.2f %8.2f \n",step1[0],step1[1],step1[2]);
    print ("thres  = %f %f\n",globals->threshold[0],globals->threshold[1]);

    for_less(i,0,3)
      print ("direct= %8.2f %8.2f %8.2f \n",
	     Point_x(directions1[i]),
	     Point_y(directions1[i]),
	     Point_z(directions1[i]));
    
  }


  fill_Point( starting_position1, start1[0], start1[1], start1[2]);
  
				/* init min1 max1 values */
  min1_col = count1[COL_IND]; max1_row = 0;
  min1_row = count1[ROW_IND]; max1_col = 0;
  min1_slice=count1[SLICE_IND]; max1_slice = 0;

  for_inclusive(s,0,count1[SLICE_IND]) {

    SCALE_VECTOR( vector_step, directions1[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position1, vector_step );

    for_inclusive(r,0,count1[ROW_IND]) {

      SCALE_VECTOR( vector_step, directions1[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,count1[COL_IND]) {

	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {	

	  if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &true_value )) {
	    
	    if (true_value > globals->threshold[0]) {
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


	ADD_POINT_VECTOR( col, col, directions1[COL_IND] );
	
      }
    }
  }

  vol1 = (max1_row - min1_row + 1) *
         (max1_col - min1_col + 1) *
         (max1_slice - min1_slice + 1);

  if (globals->flags.debug && globals->flags.verbose>1) 
    print ("volume =  %d\n",vol1);

  if ((max1_row < min1_row) || (max1_col < min1_col) || 
      (max1_slice < min1_slice)) {

    print ("slice lim %d %d\n",min1_slice, max1_slice);
    print ("row lim   %d %d\n",min1_row, max1_row);
    print ("col lim   %d %d\n",min1_col, max1_col);
    print_error_and_line_num("%s", __FILE__, __LINE__,"Cannot calculate size of volume 1\n.");
  }



				/* build default sampling lattice info
				   on second data set (d2)               */


  set_up_lattice(d2, globals->step,
		 start2, count2, step2, directions2);

  if (globals->flags.debug && globals->flags.verbose>1) {
    print ("in init_lattice: for the target data set\n");
    print ("start = %8.2f %8.2f %8.2f \n",start2[0],start2[1],start2[2]);
    print ("count = %8d %8d %8d \n",count2[0],count2[1],count2[2]);
    print ("step  = %8.2f %8.2f %8.2f \n",step2[0],step2[1],step2[2]);
    
    for_less(i,0,3)
      print ("direct= %8.2f %8.2f %8.2f \n",
	     Point_x(directions2[i]),
	     Point_y(directions2[i]),
	     Point_z(directions2[i]));

  }


  fill_Point( starting_position2, start2[0], start2[1], start2[2]);
  
				/* init min2 max2 values */
  min2_col = count2[COL_IND]; max2_row = 0;
  min2_row = count2[ROW_IND]; max2_col = 0;
  min2_slice=count2[SLICE_IND]; max2_slice = 0;

  for_inclusive(s,0,count2[SLICE_IND]) {

    SCALE_VECTOR( vector_step, directions2[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position2, vector_step );

    for_inclusive(r,0,count2[ROW_IND]) {

      SCALE_VECTOR( vector_step, directions2[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,count2[COL_IND]) {

	convert_3D_world_to_voxel(d2, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m2, Point_x(col), Point_y(col), Point_z(col))) {	                /* should be fill_value  */
	  
	  if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &true_value )) {
	    
	    if (true_value > globals->threshold[1]) {
	      if (r > max2_row) max2_row = r; 
	      if (r < min2_row) min2_row = r;
	      if (c > max2_col) max2_col = c; 
	      if (c < min2_col) min2_col = c;
	      if (s > max2_slice) max2_slice = s; 
	      if (s < min2_slice) min2_slice = s;
	    }
	  }
	} 

	ADD_POINT_VECTOR( col, col, directions2[COL_IND] );
	
      }
    }
  }

  if ((max2_row < min2_row) || (max2_col < min2_col) || 
      (max2_slice < min2_slice)) {
    print ("slice lim %d %d\n",min2_slice, max2_slice);
    print ("row lim   %d %d\n",min2_row, max2_row);
    print ("col lim   %d %d\n",min2_col, max2_col);
    print_error_and_line_num("%s", __FILE__, __LINE__,"Cannot calculate size of volume 2\n." );
  }

  vol2 = (max2_row - min2_row + 1) *
         (max2_col - min2_col + 1) *
         (max2_slice - min2_slice + 1);


  if (globals->flags.debug && globals->flags.verbose>1) 
    print ("volume =  %d\n",vol2);

  if ( !(globals->trans_info.transform_type==TRANS_NONLIN) && vol1<=vol2) {
    globals->smallest_vol = 1;
    globals->count[COL_IND] = max1_col - min1_col + 1;
    globals->count[ROW_IND] = max1_row - min1_row + 1;
    globals->count[SLICE_IND] = max1_slice - min1_slice + 1;

    fill_Point( starting_position1, start1[0], start1[1], start1[2]);
    SCALE_VECTOR( vector_step, directions1[SLICE_IND], min1_slice);
    ADD_POINT_VECTOR( slice, starting_position1, vector_step );
    SCALE_VECTOR( vector_step, directions1[ROW_IND], min1_row);
    ADD_POINT_VECTOR( row, slice, vector_step );
    SCALE_VECTOR( vector_step, directions1[COL_IND], min1_col);
    ADD_POINT_VECTOR( col, row, vector_step );

    globals->start[X] = Point_x(col);
    globals->start[Y] = Point_y(col);
    globals->start[Z] = Point_z(col);

    for_less(i,0,3) {
      Point_x(globals->directions[i]) = Point_x(directions1[i]);
      Point_y(globals->directions[i]) = Point_y(directions1[i]);
      Point_z(globals->directions[i]) = Point_z(directions1[i]);
    }
  }
  else {
				/* The lattice is always defined on the target volume
				   when computing non-linear transformations.          */

    globals->smallest_vol = 2;
    globals->count[COL_IND] = max2_col - min2_col + 1;
    globals->count[ROW_IND] = max2_row - min2_row + 1;
    globals->count[SLICE_IND] = max2_slice - min2_slice + 1;

    fill_Point( starting_position2, start2[0], start2[1], start2[2]);
    SCALE_VECTOR( vector_step, directions2[SLICE_IND], min2_slice);
    ADD_POINT_VECTOR( slice, starting_position2, vector_step );
    SCALE_VECTOR( vector_step, directions2[ROW_IND], min2_row);
    ADD_POINT_VECTOR( row, slice, vector_step );
    SCALE_VECTOR( vector_step, directions2[COL_IND], min2_col);
    ADD_POINT_VECTOR( col, row, vector_step );

    globals->start[X] = Point_x(col);
    globals->start[Y] = Point_y(col);
    globals->start[Z] = Point_z(col);

    for_less(i,0,3) {
      Point_x(globals->directions[i]) = Point_x(directions2[i]);
      Point_y(globals->directions[i]) = Point_y(directions2[i]);
      Point_z(globals->directions[i]) = Point_z(directions2[i]);
    }
    
    tmp_threshold = globals->threshold[0];
    globals->threshold[0] = globals->threshold[1];
    globals->threshold[1] = tmp_threshold;
  }

  if (globals->flags.debug && globals->flags.verbose>1)
    print ("***** leaving init_lattice\n");

}
