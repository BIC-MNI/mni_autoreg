#include <def_mni.h>

/* ----------------------------- MNI Header -----------------------------------
@NAME       : volume_functions.c
                collection of routines used to manipulate volume data.
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Jun 15 08:57:23 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */

public void make_zscore_volume(Volume d1, Volume m1, 
			       float threshold); 

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double *start, int *count, double  *size);




public void make_zscore_volume(Volume d1, Volume m1, 
			       float threshold)
{
  int 
    sizes[MAX_DIMENSIONS],
    count,
    s,r,c;
  double
    scale, offset,
    sum, sum2, mean, var, std,
    mask_vox, data_vox;

  count  = 0;
  sum  = 0.0;
  sum2 = 0.0;

  scale  = d1->value_scale;
  offset = d1->value_translation;

  get_volume_sizes( d1, sizes );

				/* do first pass, to get mean and std */
  for_less( s, 0,  sizes[Z])
    for_less( r, 0, sizes[Y])
      for_less( c, 0, sizes[X]) {

	if ( m1 != NULL ) {
	  GET_VOXEL_3D( mask_vox, m1 , c, r, s ); 
	}
	else
	  mask_vox = 1.0;

	if (mask_vox != 0) { /* should be m1->fill_value */
	  
	  GET_VOXEL_3D( data_vox,  d1 , c, r, s );
	  data_vox = data_vox * scale + offset;

	  if (data_vox > threshold) {
	    sum  += data_vox;
	    sum2 += data_vox*data_vox;
	    
	    count++;
	  }
	}
	 
      }

				/* calc mean and std */
  mean = sum / count;
  var  = (count*sum2 - sum*sum) / count*(count-1);
  std  = sqrt(var);
		
  switch( d1->data_type) {
  case UNSIGNED_BYTE:  
    d1->value_scale = (2^8  - 1)/5.0 ; d1->value_translation = 2^7;
    break;  
  case SIGNED_BYTE:  
    d1->value_scale = (2^7  - 1)/5.0 ; d1->value_translation = 0;
    break;  
  case UNSIGNED_SHORT:  
    d1->value_scale = (2^16 - 1)/5.0 ; d1->value_translation = 2^15;
    break;  
  case SIGNED_SHORT:  
    d1->value_scale = (2^15 - 1)/5.0 ; d1->value_translation = 0;
    break;  
  case UNSIGNED_LONG:  
    d1->value_scale = (2^32 - 1)/5.0 ; d1->value_translation = 2^31;
    break;  
  case SIGNED_LONG:  
    d1->value_scale = (2^31 - 1)/5.0 ; d1->value_translation = 0;
    break;  
  case FLOAT:  
    d1->value_scale = 1.0; d1->value_translation = 0.0;
    break;  
  case DOUBLE:  
    d1->value_scale = 1.0; d1->value_translation = 0.0;
    break;  
  }
				/* replace the voxel values */
  for_less( s, 0,  sizes[Z])
    for_less( r, 0, sizes[Y])
      for_less( c, 0, sizes[X]) {

	if ( m1 != NULL ) {
	  GET_VOXEL_3D( mask_vox, m1 , c, r, s ); 
	}
	else
	  mask_vox = 1.0;

	if (mask_vox != 0) { /* should be m1->fill_value */
	  
	  GET_VOXEL_3D( data_vox,  d1, c, r, s );
	  data_vox = data_vox * scale + offset;

	  if (data_vox > threshold) 
	    data_vox = CONVERT_VALUE_TO_VOXEL( d1, data_vox);
	  else
	    data_vox = 0;   /* should be fill_value! */
	
	  SET_VOXEL_3D( d1 , c, r, s, data_vox );
	}
	 
      }

}

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double  *start, int *count, double  *step)
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
    voxel;

  double
    tx,ty,tz,
    voxel_value;
  int
    xi,yi,zi,
    flip_flag,r,c,s;

  fill_Vector( row_step,   step[X], 0.0,     0.0 );
  fill_Vector( col_step,   0.0,     step[Y], 0.0 );
  fill_Vector( slice_step, 0.0,     0.0,     step[Z] );

  fill_Point( starting_position, start[X], start[Y], start[Z]);

  flip_flag = FALSE;

  for_inclusive(s,0,count[Z]) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,count[Y]) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,count[X]) {

	convert_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

	fill_Point( voxel, tx, ty, tz );

	xi = ROUND( tx );
	yi = ROUND( ty );
	zi = ROUND( tz );

	GET_VOXEL_3D( voxel_value, d1 , xi, yi, zi ); 

/*	INTERPOLATE_VOXEL_VALUE( d1, &voxel, &voxel_value ); */
	
	if (flip_flag)
	  voxel_value = voxel_value * (1 + 0.01*speckle);
	else
	  voxel_value = voxel_value * (1 - 0.01*speckle);

	flip_flag = !flip_flag;

	SET_VOXEL_3D( d1 , xi, yi, zi, voxel_value );

	ADD_POINT_VECTOR( col, col, col_step );
	
      }
    }
  }


}

