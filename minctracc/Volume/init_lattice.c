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
@MODIFIED   :  Revision 1.7  1994-04-06 11:48:35  louis
@MODIFIED   :  working linted version of linear + non-linear registration based on Lvv
@MODIFIED   :  operator working in 3D
@MODIFIED   :
 * Revision 1.6  94/02/21  16:35:29  louis
 * version before feb 22 changes
 * 
 * Revision 1.5  93/11/15  16:25:43  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Volume/init_lattice.c,v 1.7 1994-04-06 11:48:35 louis Exp $";
#endif


#include <volume_io.h>
#include <recipes.h>

#include "matrix_basics.h"
#include "cov_to_praxes.h"
#include "make_rots.h"

#include "constants.h"
#include "arg_data.h"

#include "local_macros.h"

#include <print_error.h>

extern Arg_Data main_args;



/* find the largest lattice array that will fit in the data volume,
   given the user spacing.

   return the lattice info in start, count and step. 
*/

public void set_up_lattice(Volume data, 
			    double *user_step, /* user requested spacing for lattice */
			    double *start,     /* world starting position of lattice */
			    int    *count,     /* number of steps in each direction */
			    double *step,      /* step size in each direction */
			    VectorR directions[])/* array of vector directions for each index*/
{
  int 
    sizes[3], 
    i,j;

  Real 
    sign,
    starts[3],
    start_voxel[3],
    vect_voxel[3],
    tmp_point[3],
    num_steps,
    offset[3],
    separations[3];

  BOOLEAN debug; int verbose;

  debug  = main_args.flags.debug;
  verbose= main_args.flags.verbose;

  for_less(i,0,3)
    step[i] = user_step[i]; 

				/* get the volume sizes and voxel spacing */
  get_volume_sizes(data, sizes );
  get_volume_separations(data, separations );

  

  if (debug && verbose>1) {

    convert_3D_voxel_to_world(data, 0.0,0.0,0.0,
			      &starts[0], &starts[1], &starts[2]);
    
    print ("In set_up_lattice, for orig data volume:\n");
    print ("sizes: %7d %7d %7d\n",sizes[0],sizes[1],sizes[2]);
    print ("steps: %7.2f %7.2f %7.2f\n",separations[0],separations[1],separations[2]);
    print ("start: %7.2f %7.2f %7.2f\n",starts[0],starts[1],starts[2]); 
  }
  
				/* for each direction */
  for_less( i, 0, 3) {	

				/* force step to have the same SIGN as the 
				   volume voxel spacing. */
    step[i] = ABS( step[i] );
    if (separations[i] < 0.0)  step[i] *= (-1.0);

    num_steps = separations[i] * sizes[i] / step[i];
    num_steps = ABS(num_steps);

    count[i] = ROUND(num_steps);
    
				/* this is the offset for the start of the
				   lattice from the corner of the volume 
				   in world distance mm                   */
    offset[i] = (separations[i]*sizes[i] - step[i]*count[i]) / 2.0;

  }

				/* get the voxel position of the lattice start,
				   in voxel coordinates of the data volume       */
  for_less(i,0,3) {

    if (separations[i] > 0)
      sign = 1.0;
    else
      sign = -1.0;

    start_voxel[i] = sign*((-0.5)  	/* to get to the edge of the voxel,
                                           since the voxel's  coordinates is at its center */

                    + (offset[i]/separations[i]) /* the offset to edge of lattice */

		    + step[i]/(2*separations[i]));/* the offset to the center of the 
						    lattice voxel */
    
  }
				/* get the world start.  */

  convert_3D_voxel_to_world(data, start_voxel[0],start_voxel[1], start_voxel[2],
			    &tmp_point[0], &tmp_point[1], &tmp_point[2]);
  
  for_less(i,0,3)
    start[i] = tmp_point[i];

				/* build direction vectors*/
  for_less(i,0,3) {

    for_less(j,0,3)
      vect_voxel[j] = start_voxel[j];

    vect_voxel[i] += step[i]/separations[i];

    convert_3D_voxel_to_world(data, vect_voxel[0], vect_voxel[1], vect_voxel[2],
			      &tmp_point[0], &tmp_point[1], &tmp_point[2]);

    fill_Vector( directions[i], tmp_point[0]-start[0], tmp_point[1]-start[1], tmp_point[2]-start[2]);
  }
  
  if (debug && verbose>1) {
    
    print ("                   for lattice volume:\n");
    print ("sizes: %7d %7d %7d\n",count[0],count[1],count[2]);
    print ("steps: %7.2f %7.2f %7.2f\n",step[0],step[1],step[2]);
    print ("start: %7.2f %7.2f %7.2f\n",start[0],start[1],start[2]);
    print ("leaving set_up_lattice()\n\n");
  }
 

}

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
    count1[3], count2[3];
  double 
    start1[3], start2[3],
    step1[3],  step2[3];
  VectorR
    directions1[3],
    directions2[3];  

				/* build default sampling lattice info
				   on first data set (d1)               */
  set_up_lattice(d1, globals->step,
		 start1, count1, step1, directions1);

  if (globals->flags.debug && globals->flags.verbose>1) {
    print ("start = %8.2f %8.2f %8.2f \n",start1[0],start1[1],start1[2]);
    print ("count = %8d %8d %8d \n",count1[0],count1[1],count1[2]);
    print ("step  = %8.2f %8.2f %8.2f \n",step1[0],step1[1],step1[2]);
    
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
	    if (globals->flags.debug && globals->flags.verbose>2) {
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

  if (globals->flags.debug && globals->flags.verbose>1) {

    print ("slice lim %d %d\n",min1_slice, max1_slice);
    print ("row lim   %d %d\n",min1_row, max1_row);
    print ("col lim   %d %d\n",min1_col, max1_col);
    print ("volume =  %d\n",vol1);
  }

  if ((max1_row < min1_row) || (max1_col < min1_col) || 
      (max1_slice < min1_slice)) {
    print_error("%s", __FILE__, __LINE__,"Cannot calculate size of volume 1\n.");
  }



				/* build default sampling lattice info
				   on second data set (d2)               */


  set_up_lattice(d2, globals->step,
		 start2, count2, step2, directions2);

  if (globals->flags.debug && globals->flags.verbose>1) {

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
    print_error("%s", __FILE__, __LINE__,"Cannot calculate size of volume 2\n." );
  }

  vol2 = (max2_row - min2_row + 1) *
         (max2_col - min2_col + 1) *
         (max2_slice - min2_slice + 1);


  if (globals->flags.debug && globals->flags.verbose>1) {

    print ("slice lim %d %d\n",min2_slice, max2_slice);
    print ("row lim   %d %d\n",min2_row, max2_row);
    print ("col lim   %d %d\n",min2_col, max2_col);
    print ("volume =  %d\n",vol2);
  }

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
  
}
