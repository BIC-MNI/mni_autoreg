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
#include <def_mni.h>
#include <limits.h>

public void make_zscore_volume(Volume d1, Volume m1, 
			       float threshold); 

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double *start, int *count, double  *size);

public void make_scale_and_translation_with_offset(Volume volume, 
						   long offset_bottom, 
						   long offset_top, 
						   Real real_min, 
						   Real real_max);

public void make_scale_and_translation(Volume volume, Real real_min, Real real_max);

public void save_volume(Volume d, char *filename);



public void make_zscore_volume(Volume d1, Volume m1, 
			       float threshold)
{
  unsigned long
    count;
  int 
    sizes[MAX_DIMENSIONS],
    s,r,c;
  Real
    valid_min_mvoxel, valid_max_mvoxel,
    valid_min_dvoxel, valid_max_dvoxel,
    min,max,
    scale, offset,
    sum, sum2, mean, var, std,
    mask_vox, data_vox,data_val;

  count  = 0;
  sum  = 0.0;
  sum2 = 0.0;
  min = FLT_MAX;
  max = -FLT_MAX;

  scale  = d1->value_scale;
  offset = d1->value_translation;

  get_volume_sizes( d1, sizes );

  get_volume_voxel_range(d1, &valid_min_dvoxel, &valid_max_dvoxel);
  get_volume_voxel_range(m1, &valid_min_mvoxel, &valid_max_mvoxel);

				/* do first pass, to get mean and std */
  for_less( s, 0,  sizes[Z]) {
    for_less( r, 0, sizes[Y]) {
      for_less( c, 0, sizes[X]) {

	if ( m1 != NULL ) {
	  GET_VOXEL_3D( mask_vox, m1 , c, r, s ); 
	}
	else
	  mask_vox = -DBL_MAX;

	if (mask_vox >= valid_min_mvoxel && mask_vox <= valid_max_mvoxel) { 
	  
	  GET_VOXEL_3D( data_vox,  d1 , c, r, s );

				/* instead of   
				   data_val = CONVERT_VOXEL_TO_VALUE(d1, data_vox);

				   I use the scale and offset values saved at the beginning
				   since these values are changed with respect to the
				   new z-score volume*/

	  if (data_vox >= valid_min_dvoxel && data_vox <= valid_max_dvoxel) { 
	    data_val = data_vox*scale + offset;
	    
	    if (data_val > threshold) {
	      sum  += data_val;
	      sum2 += data_val*data_val;
	      
	      count++;
	      
	      if (data_val < min)
		min = data_val;
	      else
		if (data_val > max)
		  max = data_val;
	    }
	  }
	}
      }
    }
  }


				/* calc mean and std */
  mean = sum / (float)count;
  var  = ((float)count*sum2 - sum*sum) / ((float)count*((float)count-1));
  std  = sqrt(var);

  make_scale_and_translation(d1,-5.0, 5.0);

  min = FLT_MAX;
  max = -FLT_MAX;

				/* replace the voxel values */
  for_less( s, 0,  sizes[Z]) {
    for_less( r, 0, sizes[Y]) {
      for_less( c, 0, sizes[X]) {
	
	GET_VOXEL_3D( data_vox,  d1, c, r, s );
	
	if (data_vox >= valid_min_dvoxel && data_vox <= valid_max_dvoxel) { 
	  
	  data_val = data_vox*scale + offset;
	  
	  if (data_val > threshold) {
	    data_val = (data_val - mean) / std;
	    data_vox = CONVERT_VALUE_TO_VOXEL( d1, data_val);
	    
	    if (data_val < min) {
	      min = data_val;
	    }
	    else {
	      if (data_val > max)
		max = data_val;
	    }
	  }
	  else
	    data_vox = -DBL_MAX;   /* should be fill_value! */
	  
	  SET_VOXEL_3D( d1 , c, r, s, data_vox );
	}
	
      }
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
    col;
  Real valid_min_voxel, valid_max_voxel;

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

  get_volume_voxel_range(d1,&valid_min_voxel, &valid_max_voxel);

  for_inclusive(s,0,count[Z]) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,count[Y]) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,count[X]) {

	convert_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

	xi = ROUND( tx );
	yi = ROUND( ty );
	zi = ROUND( tz );

	GET_VOXEL_3D( voxel_value, d1 , xi, yi, zi ); 


	if (voxel_value >= valid_min_voxel && voxel_value <= valid_max_voxel) {
	  if (flip_flag)
	    voxel_value = voxel_value * (1 + 0.01*speckle);
	  else
	    voxel_value = voxel_value * (1 - 0.01*speckle);

	  SET_VOXEL_3D( d1 , xi, yi, zi, voxel_value );
	}

	flip_flag = !flip_flag;

	ADD_POINT_VECTOR( col, col, col_step );
	
      }
    }
  }


}

public void make_scale_and_translation_with_offset(Volume volume, 
						   long offset_bottom, 
						   long offset_top, 
						   Real real_min, 
						   Real real_max)
{
  Real
    min_value,
    max_value;


  switch( volume->data_type) {
  case UNSIGNED_BYTE:  
    min_value = 0+offset_bottom ; max_value = UCHAR_MAX-offset_top;
    break;  
  case SIGNED_BYTE:  
    min_value =  SCHAR_MIN+offset_bottom; max_value = SCHAR_MAX-offset_top;
    break;  
  case UNSIGNED_SHORT:  
    min_value = 0+offset_bottom ; max_value = USHRT_MAX-offset_top;
    break;  
  case SIGNED_SHORT:  
    min_value = SHRT_MIN+offset_bottom ; max_value = SHRT_MAX-offset_top;
    break;  
  case UNSIGNED_LONG:  
    min_value = 0+offset_bottom ; max_value = ULONG_MAX-offset_top;
    break;  
  case SIGNED_LONG:  
    min_value = LONG_MIN+offset_bottom ; max_value = LONG_MAX-offset_top;
    break;  
  case FLOAT:  
    volume->value_scale = 1.0; volume->value_translation = 0.0;
    return;
    break;  
  case DOUBLE:  
    volume->value_scale = 1.0; volume->value_translation = 0.0;
    return;
    break;  
  }

  volume->min_voxel = min_value;
  volume->max_voxel = max_value;

  if( real_min == real_max )
    volume->value_scale = 1.0;
  else
    volume->value_scale = (real_max - real_min) / (max_value - min_value);

  volume->value_translation = real_min - min_value * volume->value_scale;

}


public void make_scale_and_translation(Volume volume, Real real_min, Real real_max)
{

  make_scale_and_translation_with_offset(volume, (long)0, (long)0, real_min, real_max);

}


public void save_volume(Volume d, char *filename)
{
  Minc_file minc_fp;
  Real vox,val,min_val,max_val;
  int i,j,k,sizes[3];
  Status status;

  get_volume_sizes(d, sizes);
  min_val = FLT_MAX;
  max_val = -FLT_MAX;

  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]) {
	GET_VOXEL_3D(vox, d, i,j,k);
	val = CONVERT_VOXEL_TO_VALUE(d, vox);

	if (val < min_val)
	  min_val = val;
	else
	  if (val > max_val)
	    max_val = val;

      }


  minc_fp = initialize_minc_output(filename, 3, d->dimension_names, d->sizes,
                                   d->nc_data_type, FALSE, min_val, max_val,
                                   &(d->voxel_to_world_transform));

  status = output_minc_volume(minc_fp, d);

  if (status == OK)
    close_minc_output(minc_fp);
  else
    print_error("problems writing  volume `%s'.",__FILE__, __LINE__, filename,0,0,0,0);

}
