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
@MODIFIED   :  $Log: volume_functions.c,v $
@MODIFIED   :  Revision 1.6  1993-11-15 16:27:13  louis
@MODIFIED   :  working version, with new library, with RCS revision stuff,
@MODIFIED   :  before deformations included
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Volume/volume_functions.c,v 1.6 1993-11-15 16:27:13 louis Exp $";
#endif

#include <mni.h>
#include <limits.h>
#include "point_vector.h"
#include "constants.h"

public void make_zscore_volume(Volume d1, Volume m1, 
			       float threshold); 

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double  *start, int *count, VectorR directions[]);

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
    wx,wy,wz,
    valid_min_mvoxel, valid_max_mvoxel,
    valid_min_dvoxel, valid_max_dvoxel,
    min,max,
    sum, sum2, mean, var, std,
    data_vox,data_val,
    thick[MAX_DIMENSIONS];

  PointR 
    voxel;

  Volume 
    vol;

  /* get default information from data and mask */

  /* build temporary working volume */
 
  vol = copy_volume_definition(d1, NC_UNSPECIFIED, FALSE, 0.0, 0.0);
  set_volume_real_range(vol, -5.0, 5.0);
  get_volume_sizes(d1, sizes);
  get_volume_separations(d1, thick);
  get_volume_voxel_range(d1, &valid_min_dvoxel, &valid_max_dvoxel);
  get_volume_voxel_range(m1, &valid_min_mvoxel, &valid_max_mvoxel);

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

	convert_3D_voxel_to_world(d1, (Real)s, (Real)r, (Real)c, &wx, &wy, &wz);
	convert_3D_world_to_voxel(m1, wx, wy, wz, &Point_x(voxel), &Point_y(voxel), &Point_z(voxel));

	if (point_not_masked(m1, wx,wy,wz)) {
	  
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
				  double  *start, int *count, VectorR directions[])
{
  VectorR
    vector_step;

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


  flip_flag = FALSE;

  get_volume_voxel_range(d1, &valid_min_voxel, &valid_max_voxel);
  fill_Point( starting_position, start[0], start[1], start[2]);
  
  for_inclusive(s,0,count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,count[ROW_IND]) {

      SCALE_VECTOR( vector_step, directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );

      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,count[COL_IND]) {

	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

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


	ADD_POINT_VECTOR( col, col, directions[COL_IND] );

      }
	


    }
  }


}


public void save_volume(Volume d, char *filename)
{
  Status status;

  status = output_volume(filename,NC_UNSPECIFIED, FALSE, 0.0, 0.0, d, (char *)NULL,
			 (minc_output_options *)NULL);

  if (status != OK)
    print_error("problems writing  volume `%s'.",__FILE__, __LINE__, filename);

}
