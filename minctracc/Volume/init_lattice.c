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
@CREATED    : Wed Jun  9 12:56:08 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
#include <def_mni.h>
#include <recipes.h>

#include "matrix_basics.h"
#include "cov_to_praxes.h"
#include "make_rots.h"

#include "def_constants.h"
#include "def_arg_data.h"

#include "local_macros.h"

#include <print_error.h>

extern Arg_Data main_args;


public void init_lattice(Volume d1,
			 Volume d2,
			 Volume m1,
			 Volume m2, 
			 Arg_Data *globals)
{

  VectorR
    vector_step,
    slice_step,
    row_step,
    col_step;

  PointR 
    starting_offset,
    starting_origin,
    starting_position1,
    starting_position2,
    slice,
    row,
    col,
    voxel;

  double
    tx,ty,tz,
    voxel_value;
  int
    i,r,c,s,
    limits[VOL_NDIMS];

  float
    t;

  int 
    vol1,vol2,
    min1_row,   max1_row,
    min1_col,   max1_col,
    min1_slice, max1_slice,
    min2_row,   max2_row,
    min2_col,   max2_col,
    min2_slice, max2_slice;

  Real 
    thickness1[3],thickness2[3];
  int 
    sizes1[3],sizes2[3];
  Real
    rx,ry,rz,
    true_value,
    mask_value;


  get_volume_sizes(d1, sizes1 );
  get_volume_separations(d1, thickness1 );

  get_volume_sizes(d2, sizes2 );
  get_volume_separations(d2, thickness2 );

  
				/* build default sampling lattice info
				   on first data set (d1)               */
  for_less( i, 0, 3) {	
    globals->step[X] *= thickness1[i] / ABS( thickness1[i]);
  }


  fill_Vector( col_step,   globals->step[COL_IND], 0.0,     0.0 );
  fill_Vector( row_step,   0.0,     globals->step[ROW_IND], 0.0 );
  fill_Vector( slice_step, 0.0,     0.0,     globals->step[SLICE_IND] );

  convert_voxel_to_world(d1, 0.0, 0.0, 0.0, &tx, &ty, &tz);

  fill_Point( starting_origin, tx, ty, tz);

  for_less( i, 0, 3) {		/* for each dim, get # of steps in that direction,
				   and set starting offset */
    t = sizes1[i] * thickness1[i] / globals->step[i];
    limits[i] = (int)( ABS( t ) );
    
    Point_coord( starting_offset, 2-i ) = 
      ( (sizes1[i]-1)*thickness1[i] - (limits[i] * globals->step[i] ) ) / 2.0;
  }
  
  ADD_POINTS( starting_position1, starting_origin, starting_offset ); /*  */
  
				/* init min1 max1 values */
  min1_col = limits[COL_IND]; max1_row = 0;
  min1_row = limits[ROW_IND]; max1_col = 0;
  min1_slice=limits[SLICE_IND]; max1_slice = 0;

  for_inclusive(s,0,limits[SLICE_IND]) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position1, vector_step );

    for_inclusive(r,0,limits[ROW_IND]) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,limits[COL_IND]) {

	convert_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (m1 != NULL) {
	  nearest_neighbour_interpolant( m1, &voxel, &mask_value); 
	}
	else
	  mask_value = 1.0;
	
	if (mask_value > 0.0) {	                /* should be fill_value  */
	  
	  if (INTERPOLATE_VOXEL_VALUE( d1, &voxel, &voxel_value )) {
	    
	    true_value = CONVERT_VOXEL_TO_VALUE(d1,voxel_value );
	    
	    if (true_value > globals->threshold) {
	      if (r > max1_row) max1_row = r; 
	      if (r < min1_row) min1_row = r;
	      if (c > max1_col) max1_col = c; 
	      if (c < min1_col) min1_col = c;
	      if (s > max1_slice) max1_slice = s; 
	      if (s < min1_slice) min1_slice = s;
	    }
	  }    
	} /* if mask_value > 0.0 */


	ADD_POINT_VECTOR( col, col, col_step );
	
      }
    }
  }

  vol1 = (max1_row - min1_row + 1) *
         (max1_col - min1_col + 1) *
         (max1_slice - min1_slice + 1);

  if ((max1_row < min1_row) || (max1_col < min1_col) || 
      (max1_slice < min1_slice)) {
    print_error("Cannot calculate size of volume 1\n.", __FILE__, __LINE__);
  }



				/* build default sampling lattice info
				   on second data set (d2)               */
  for_less( i, 0, 3) {	
    globals->step[i] *= thickness2[i] / ABS( thickness2[i]);
  }

  fill_Vector( col_step,   globals->step[COL_IND], 0.0,     0.0 );
  fill_Vector( row_step,   0.0,     globals->step[ROW_IND], 0.0 );
  fill_Vector( slice_step, 0.0,     0.0,     globals->step[SLICE_IND] );

  convert_voxel_to_world(d2, 0.0, 0.0, 0.0, &tx, &ty, &tz);

  fill_Point( starting_origin, tx, ty, tz);

  for_less( i, 0, 3) {		/* for each dim, get # of steps in that direction,
				   and set starting offset */
    t = sizes2[i] * thickness2[i] / globals->step[i];
    limits[i] = (int)( ABS( t ) );
    
    Point_coord( starting_offset, (2-i) ) = 
      ( (sizes2[i]-1)*thickness2[i] - (limits[i] * globals->step[i] ) ) / 2.0;
  }
  
  ADD_POINTS( starting_position2, starting_origin, starting_offset ); /*  */
  
				/* init min2 max2 values */
  min2_row = limits[ROW_IND]; max2_row = 0;
  min2_col = limits[COL_IND]; max2_col = 0;
  min2_slice=limits[SLICE_IND]; max2_slice = 0;

  for_inclusive(s,0,limits[SLICE_IND]) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position2, vector_step );

    for_inclusive(r,0,limits[ROW_IND]) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,limits[COL_IND]) {
	
	convert_world_to_voxel(d2, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (m2 != NULL) {
	  nearest_neighbour_interpolant( m2, &voxel, &mask_value); 
	}
	else
	  mask_value = 1.0;
	
	if (mask_value > 0.0) {	                /* should be fill_value  */
	  
	  if (INTERPOLATE_VOXEL_VALUE( d2, &voxel, &voxel_value )) {
	    
	    true_value = CONVERT_VOXEL_TO_VALUE(d2,voxel_value );
	    
	    if (true_value > globals->threshold) {
	      if (r > max2_row) max2_row = r; 
	      if (r < min2_row) min2_row = r;
	      if (c > max2_col) max2_col = c; 
	      if (c < min2_col) min2_col = c;
	      if (s > max2_slice) max2_slice = s; 
	      if (s < min2_slice) min2_slice = s;
	    }
	  }
	} /* if mask_value > 0.0 */

	ADD_POINT_VECTOR( col, col, col_step );
	
      }
    }
  }

  if ((max2_row < min2_row) || (max2_col < min2_col) || 
      (max2_slice < min2_slice)) {
    print_error("Cannot calculate size of volume 2\n.", __FILE__, __LINE__ );
  }

  vol2 = (max2_row - min2_row + 1) *
         (max2_col - min2_col + 1) *
         (max2_slice - min2_slice + 1);

  if (vol1<vol2) {
    globals->smallest_vol = 1;
    globals->count[COL_IND] = max1_col - min1_col + 1;
    globals->count[ROW_IND] = max1_row - min1_row + 1;
    globals->count[SLICE_IND] = max1_slice - min1_slice + 1;

    convert_voxel_to_world(d1, (Real)min1_col, (Real)min1_row, (Real)min1_slice, &rx, &ry, &rz);
    globals->start[X] = rx;
    globals->start[Y] = ry;
    globals->start[Z] = rz;
  }
  else {
    globals->smallest_vol = 2;
    globals->count[COL_IND] = max2_col - min2_col + 1;
    globals->count[ROW_IND] = max2_row - min2_row + 1;
    globals->count[SLICE_IND] = max2_slice - min2_slice + 1;
    convert_voxel_to_world(d2, (Real)min2_col, (Real)min2_row, (Real)min2_slice, &rx, &ry, &rz);
    globals->start[X] = rx;
    globals->start[Y] = ry;
    globals->start[Z] = rz;

  }
  
}
