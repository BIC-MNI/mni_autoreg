/* ----------------------------- MNI Header -----------------------------------
@NAME       : objectives.c
@DESCRIPTION: File containing objective function routines.
@METHOD     : 
@GLOBALS    : 
@CREATED    : Wed Jun  9 12:56:08 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */


#include <def_mni.h>

#define VOL_NDIMS 3

#include "minctracc.h"


/* ----------------------------- MNI Header -----------------------------------
@NAME       : xcorr_objective
@INPUT      : volumetric data, for use in correlation.  

	      Correlation is accomplished under the currently defined
	      affine transformation in `globals->trans_info.transformation.trans_data'

	      The cross correlation is calculated using sub-sampling
	      in each of the volumes.

@OUTPUT     : none

@RETURNS    : the normalized cross correlation value (min=0.0 max = 1.0)

@DESCRIPTION: this routine does volumetric subsampling in world space.
              These world coordinates are mapped back into each 
	      data volume to get the actual value at the given location.
	      The globals lattice info is used to calculate
              positions of the sub-samples in each vol.

	      The correlation is equal to:

	      r =    f1 / ( sqrt(f2) * sqrt(f3) )
              
	      where:          
	             f1 = sum( d1 .* d2)    (point to point multiply)
		     f2 = sum( d1 .^ 2 )    (sum of all points squared)
		     f3 = sum( d2 .^ 2 )    (sum of all points squared)

		    note: -point refers to the value at a sub-sample
		          -the xmat is applied to d2, before calculating r
		              
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 4, 1992 lc original in fit_vol
@MODIFIED   : Wed Jun 16 10:08:46 EST 1993 LC
        new code for minc tracc, copied from routines in fit_vol.
---------------------------------------------------------------------------- */
public float xcorr_objective(Volume d1,
			     Volume d2,
			     Volume m1,
			     Volume m2, 
			     Arg_Data *globals)
{

  Vector
    vector_step,
    slice_step,
    row_step,
    col_step;

  Point 
    starting_position,
    slice,
    row,
    col,
    col2,
    voxel;

  double
    tx,ty,tz;
  int
    indices[3],
    xi,yi,zi,
    flip_flag,r,c,s;

  Real
    position[3],
    value1, value2,
    voxel_value,
    mask_value;
  
  Real
    s1,s2,s3;                   /* to store the sums for f1,f2,f3 */
  float 
    result;				/* the result */
  int 
    count1,count2;


  fill_Vector( row_step,   globals->step[X], 0.0,              0.0 );
  fill_Vector( col_step,   0.0,              globals->step[Y], 0.0 );
  fill_Vector( slice_step, 0.0,              0.0,              globals->step[Z] );

  fill_Point( starting_position, globals->start[X], globals->start[Y], globals->start[Z]);

  flip_flag = FALSE;

  s1 = s2 = s3 = 0.0;
  count1 = count2 = 0;

  for_inclusive(s,0,globals->count[Z]) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,globals->count[Y]) {
      
      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,globals->count[X]) {
	
	convert_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	position[X] = tx;
	position[Y] = ty;
	position[Z] = tz;
	
	if (voxel_is_within_volume( d1, position ) ) {
	  count1++;
	  fill_Point( voxel, tx, ty, tz );

	  if (m1 != NULL) {
	    xi = ROUND( tx );
	    yi = ROUND( ty );
	    zi = ROUND( tz );
	    
	    GET_VOXEL_3D( mask_value, m1 , xi, yi, zi ); 
	  }
	  else
	    mask_value = 1.0;
	  
	  
	  if (mask_value > 0.0) {	                /* should be fill_value  */
	    INTERPOLATE_VOXEL_VALUE( d1, &voxel, &voxel_value ); 
	    value1 = CONVERT_VOXEL_TO_VALUE( d1, voxel_value);

	    
	    do_linear_transformation_point(&col2, globals->trans_info.transformation.trans_data, 
					   &col);
	    
	    convert_world_to_voxel(d2, Point_x(col2), Point_y(col2), Point_z(col2), &tx, &ty, &tz);
	    
	    position[X] = tx;
	    position[Y] = ty;
	    position[Z] = tz;
	    
	    if (voxel_is_within_volume( d2, position ) ) {
	      count2++;
	      fill_Point( voxel, tx, ty, tz );
	      
	      if (m2 != NULL) {
		xi = ROUND( tx );
		yi = ROUND( ty );
		zi = ROUND( tz );
		
		GET_VOXEL_3D( mask_value, m2 , xi, yi, zi ); 
	      }
	      else
		mask_value = 1.0;
	      
	      if (mask_value > 0.0) {	                /* should be fill_value  */
		
		INTERPOLATE_VOXEL_VALUE( d2, &voxel, &voxel_value ); 
		value2 = CONVERT_VOXEL_TO_VALUE( d1, voxel_value);
		
		if (value1 > globals->threshold || value2 > globals->threshold ) {
		  
		  s1 += value1*value2;
		  s2 += value1*value1;
		  s3 += value2*value2;
		  
		} 
		
	      } /* if mask_value, on volume 2 */
	    } /* if voxel in volume two */
	  } /* if mask_value, on volume 1 */
	} /* if voxel in volume one */
	
	ADD_POINT_VECTOR( col, col, col_step );
	
      } /* for c */
    } /* for r */
  } /* for s */

  result = s1 / (fsqrt(s2)*fsqrt(s3));

  if (globals->flags.debug) printf ("%7d %7d -> %10.8f\n",count1,count2,result);
  
  return (result);
  
}


public float zscore_objective(Volume d1,
			      Volume d2,
			      Volume m1,
			      Volume m2, 
			      Arg_Data *globals)
{
}

public float vr_objective(Volume d1,
			  Volume d2,
			  Volume m1,
			  Volume m2, 
			  Arg_Data *globals)
{
}

public float ssc_objective(Volume d1,
			   Volume d2,
			   Volume m1,
			   Volume m2, 
			   Arg_Data *globals)
{
}

