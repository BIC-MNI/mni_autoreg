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
#include "def_point_vector.h"
#include "def_constants.h"

public void make_zscore_volume(Volume d1, Volume m1, 
			       float threshold); 

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double *start, int *count, double  *size);

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
    sum, sum2, mean, var, std,
    mask_vox, data_vox,data_val,
    thick[MAX_DIMENSIONS];

  Volume vol;
  General_transform *transform;

  /* get default information from data and mask */

  get_volume_sizes(d1, sizes);
  get_volume_separations(d1, thick);
  transform = get_voxel_to_world_transform(d1);
  get_volume_voxel_range(d1,&valid_min_dvoxel, &valid_max_dvoxel);

  get_volume_voxel_range(m1, &valid_min_mvoxel, &valid_max_mvoxel);

  /* build temporary working volume */
 
  vol = create_volume(3, d1->dimension_names,  NC_FLOAT, TRUE);
  set_volume_size(vol, NC_FLOAT, TRUE, sizes);
  set_volume_separations(vol, thick);
  set_volume_voxel_range(vol, valid_min_dvoxel, valid_max_dvoxel);
  set_volume_real_range(vol, -5.0, 5.0);
  set_voxel_to_world_transform(vol, transform);


  /* initialize counters and sums */

  count  = 0;
  sum  = 0.0;
  sum2 = 0.0;
  min = FLT_MAX;
  max = -FLT_MAX;

				/* do first pass, to get mean and std */
  for_less( s, 0,  sizes[0]) {
    for_less( r, 0, sizes[1]) {
      for_less( c, 0, sizes[2]) {

	if ( m1 != NULL ) {
	  GET_VOXEL_3D( mask_vox, m1 , s, r, c ); 
	}
	else
	  mask_vox = -DBL_MAX;
	
	if (mask_vox >= valid_min_mvoxel && mask_vox <= valid_max_mvoxel) { 
	  
	  GET_VOXEL_3D( data_vox,  d1 , s, r, c );

	  if (data_vox >= valid_min_dvoxel && data_vox <= valid_max_dvoxel) { 

	    data_val = CONVERT_VOXEL_TO_VALUE(d1, data_vox);
	    
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

  min = FLT_MAX;
  max = -FLT_MAX;

				/* replace the voxel values */
  for_less( s, 0,  sizes[0]) {
    for_less( r, 0, sizes[1]) {
      for_less( c, 0, sizes[2]) {
	
	GET_VOXEL_3D( data_vox,  d1, s, r, c );
	
	if (data_vox >= valid_min_dvoxel && data_vox <= valid_max_dvoxel) { 
	  
	  data_val = CONVERT_VOXEL_TO_VALUE(d1, data_vox);
	  
	  if (data_val > threshold) {

				/* instead of   
				   data_val = CONVERT_VALUE_TO_VOXEL(d1, data_vox);
				   i will use
				   data_val = CONVERT_VALUE_TO_VOXEL(d1, vol);

				   since the values in vol are changed with respect to the
				   new z-score volume */

	    data_val = (data_val - mean) / std;
	    data_vox = CONVERT_VALUE_TO_VOXEL( vol, data_val);
	    
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
	  
	  SET_VOXEL_3D( d1 , s, r, c, data_vox );
	}
	
      }
    }
  }

  delete_volume(vol);
  
}

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double  *start, int *count, double  *step)
{
  VectorR
    vector_step,
    slice_step,
    row_step,
    col_step;

  PointR 
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

  fill_Vector( col_step,   step[COL_IND], 0.0,     0.0 );
  fill_Vector( row_step,   0.0,     step[ROW_IND], 0.0 );
  fill_Vector( slice_step, 0.0,     0.0,     step[SLICE_IND] );

  fill_Point( starting_position, start[X], start[Y], start[Z]);

  flip_flag = FALSE;

  get_volume_voxel_range(d1,&valid_min_voxel, &valid_max_voxel);

  for_inclusive(s,0,count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, slice_step, s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,count[ROW_IND]) {

      SCALE_VECTOR( vector_step, row_step, r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,count[COL_IND]) {

	convert_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

	xi = ROUND( tx );
	yi = ROUND( ty );
	zi = ROUND( tz );

	GET_VOXEL_3D( voxel_value, d1 , zi, yi, xi ); 


	if (voxel_value >= valid_min_voxel && voxel_value <= valid_max_voxel) {
	  if (flip_flag)
	    voxel_value = voxel_value * (1 + 0.01*speckle);
	  else
	    voxel_value = voxel_value * (1 - 0.01*speckle);

	  SET_VOXEL_3D( d1 , zi, yi, xi, voxel_value );
	}

	flip_flag = !flip_flag;

	ADD_POINT_VECTOR( col, col, col_step );
	
      }
    }
  }


}


public void save_volume(Volume d, char *filename)
{
  Status status;

  status = output_volume(filename,FALSE, d, NULL);

  if (status != OK)
    print_error("problems writing  volume `%s'.",__FILE__, __LINE__, filename);

}
