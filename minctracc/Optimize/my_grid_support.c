#include <volume_io.h>
#include <louis_splines.h>
#include <print_error.h>

public void get_volume_XYZV_indices(Volume data, int xyzv[]);

public void interpolate_deformation_slice(Volume volume, 
					  Real wx,Real wy,Real wz,
					  Real def[]);

public  void  grid_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed );

public  void  grid_inverse_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed );

private void transform_or_inverse_point_in_trans_plane(General_transform *transform,
					  BOOLEAN           inverse_flag,
					  Real              x, 
					  Real              y, 
					  Real              z,
					  Real              *x_transformed,  
					  Real              *y_transformed,  
					  Real              *z_transformed);

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







/*********************************************************************** 
  this routine will use interpolation on the 2D deformation
  field to calculate the inverse position */

private void transform_or_inverse_point_in_trans_plane(General_transform *transform,
					  BOOLEAN           inverse_flag,
					  Real              x, 
					  Real              y, 
					  Real              z,
					  Real              *x_transformed,  
					  Real              *y_transformed,  
					  Real              *z_transformed) {
  
  int      trans;
  
  switch( transform->type )
    {
    case LINEAR:
      if( inverse_flag )
	transform_point( transform->inverse_linear_transform,
			x, y, z,
			x_transformed, y_transformed, z_transformed );
      else
	transform_point( transform->linear_transform,
			x, y, z,
			x_transformed, y_transformed, z_transformed );
      break;
      
    case THIN_PLATE_SPLINE:
      if( inverse_flag )
        {
	  thin_plate_spline_inverse_transform( transform->n_dimensions,
					      transform->n_points,
					      transform->points,
					      transform->displacements,
					      x, y, z,
					      x_transformed, y_transformed,
					      z_transformed );
        }
      else
        {
	  thin_plate_spline_transform( transform->n_dimensions,
				      transform->n_points,
				      transform->points,
				      transform->displacements,
				      x, y, z,
				      x_transformed, y_transformed,
				      z_transformed );
        }
      break;
      
    case GRID_TRANSFORM:
      if( inverse_flag )
        {
	  grid_inverse_transform_point_in_trans_plane( transform,
					  x, y, z,
					  x_transformed, y_transformed,
					  z_transformed );
        }
      else
        {
	  grid_transform_point_in_trans_plane( transform,
				  x, y, z,
				  x_transformed, y_transformed,
				  z_transformed );
        }
      break;
      
    case USER_TRANSFORM:
      if( inverse_flag )
        {
	  transform->user_inverse_transform_function(
		transform->user_data, x, y, z,
		x_transformed, y_transformed, z_transformed );
        }
      else
        {
	  transform->user_transform_function(
		transform->user_data, x, y, z,
		x_transformed, y_transformed, z_transformed );
        }
      break;
      
    case CONCATENATED_TRANSFORM:
      *x_transformed = x;
      *y_transformed = y;
      *z_transformed = z;
      
      if( transform->inverse_flag )
	inverse_flag = !inverse_flag;
      
      if( inverse_flag )
        {
	  for( trans = transform->n_transforms-1;  trans >= 0;  --trans )
            {
	      general_inverse_transform_point_in_trans_plane(&transform->transforms[trans],
		  *x_transformed, *y_transformed, *z_transformed,
		  x_transformed, y_transformed, z_transformed );
            }
        }
      else
        {
	  for_less( trans, 0, transform->n_transforms )
            {
	      general_transform_point_in_trans_plane( &transform->transforms[trans],
		  *x_transformed, *y_transformed, *z_transformed,
		  x_transformed, y_transformed, z_transformed );
            }
        }
      break;
      
    default:
      handle_internal_error( "transform_or_invert_point" );
      break;
    }
}

public  void  general_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed )
{

    transform_or_inverse_point_in_trans_plane( transform, transform->inverse_flag, 
				       x, y, z,
				       x_transformed, 
				       y_transformed, 
				       z_transformed );
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : general_inverse_transform_point
@INPUT      : transform
              x
              y
              z
@OUTPUT     : x_transformed
              y_transformed
              z_transformed
@RETURNS    : 
@DESCRIPTION: Transforms a point by the inverse of the general transform.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 1993            David MacDonald
@MODIFIED   : 1995 louis
---------------------------------------------------------------------------- */

public  void  general_inverse_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed )
{

    transform_or_inverse_point_in_trans_plane( transform, !transform->inverse_flag, 
				       x, y, z,
				       x_transformed, 
				       y_transformed, 
				       z_transformed );
}



public  void  grid_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed )
{
    Real    displacements[MAX_DIMENSIONS];
    Volume  volume;

    /* --- the volume that defines the transform is an offset vector,
           so evaluate the volume at the given position and add the
           resulting offset to the given position */

    volume = (Volume) transform->displacement_volume;

    interpolate_deformation_slice(volume, x, y, z, displacements);

    *x_transformed = x + displacements[X];
    *y_transformed = y + displacements[Y];
    *z_transformed = z + displacements[Z];
}


#define  NUMBER_TRIES  10

public  void  grid_inverse_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed )
{
    int    tries;
    Real   best_x, best_y, best_z;
    Real   tx, ty, tz;
    Real   gx, gy, gz;
    Real   error_x, error_y, error_z, error, smallest_e;
    Real   ftol;

    ftol = 0.05;

    grid_transform_point_in_trans_plane( transform, x, y, z, &tx, &ty, &tz );

    tx = x - (tx - x);
    ty = y - (ty - y);
    tz = z - (tz - z);

    grid_transform_point_in_trans_plane( transform, tx, ty, tz, &gx, &gy, &gz );

    error_x = x - gx;
    error_y = y - gy;
    error_z = z - gz;

    tries = 0;

    error = smallest_e = ABS(error_x) + ABS(error_y) + ABS(error_z);

    best_x = tx;
    best_y = ty;
    best_z = tz;

    while( ++tries < NUMBER_TRIES && smallest_e > ftol )
    {
        tx += 0.67 * error_x;
        ty += 0.67 * error_y;
        tz += 0.67 * error_z;

        grid_transform_point_in_trans_plane( transform, tx, ty, tz, &gx, &gy, &gz );

        error_x = x - gx;
        error_y = y - gy;
        error_z = z - gz;
    
        error = ABS(error_x) + ABS(error_y) + ABS(error_z);

        if( error < smallest_e )
        {
            smallest_e = error;
            best_x = tx;
            best_y = ty;
            best_z = tz;
        }
    }

    *x_transformed = best_x;
    *y_transformed = best_y;
    *z_transformed = best_z;
}

public void interpolate_deformation_slice(Volume volume, 
					  Real wx,Real wy,Real wz,
					  Real def[])
{
  Real
    world[N_DIMENSIONS],
    v00,v01,v02,v03, 
    v10,v11,v12,v13, 
    v20,v21,v22,v23, 
    v30,v31,v32,v33;
  long 
    ind0, ind1, ind2;
  Real
    voxel[MAX_DIMENSIONS],
    frac[MAX_DIMENSIONS];
  int 
    i,
    xyzv[MAX_DIMENSIONS],
    sizes[MAX_DIMENSIONS];
  double temp_result;
  
  
  def[X] = def[Y] = def[Z] = 0.0;
  world[X] = wx; world[Y] = wy; world[Z] = wz;
  convert_world_to_voxel(volume, wx, wy, wz, voxel);
  
  /* Check that the coordinate is inside the volume */
  
  get_volume_sizes(volume, sizes);
  get_volume_XYZV_indices(volume, xyzv);

  if ((voxel[ xyzv[X] ] < 0) || (voxel[ xyzv[X] ] >=sizes[ xyzv[X]]) ||
      (voxel[ xyzv[Y] ] < 0) || (voxel[ xyzv[Y] ] >=sizes[ xyzv[Y]]) ||
      (voxel[ xyzv[Z] ] < 0) || (voxel[ xyzv[Z] ] >=sizes[ xyzv[Z]]))  {
    return ;
  }
  
  if (/*(sizes[ xyzv[Z] ] == 1) &&*/ xyzv[Z]==1) {
    
    /* Get the whole and fractional part of the coordinate */
    ind0 = FLOOR( voxel[ xyzv[Z] ] );
    ind1 = FLOOR( voxel[ xyzv[Y] ] );
    ind2 = FLOOR( voxel[ xyzv[X] ] );
    frac[Y] = voxel[ xyzv[Y] ] - ind1;
    frac[X] = voxel[ xyzv[X] ] - ind2;
    
    if (sizes[xyzv[X]] < 4 || sizes[xyzv[Y]] < 4) {
      def[X] = get_volume_real_value(volume, 0, ind0, ind1, ind2, 0);
      def[Y] = get_volume_real_value(volume, 1, ind0, ind1, ind2, 0);
      def[Z] = get_volume_real_value(volume, 2, ind0, ind1, ind2, 0);
      return;
    }

    if (ind1==0) 
      frac[Y] -= 1.0;
    else {
      ind1--;
      while ( ind1+3 >= sizes[xyzv[Y]] ) {
	ind1--;
	frac[Y] += 1.0;
      }
    }
    if (ind2==0) 
      frac[X] -= 1.0;
    else {
      ind2--;
      while ( ind2+3 >= sizes[xyzv[X]] ) {
	ind2--;
	frac[X] += 1.0;
      }
    }
    for_less(i,0,3) {
      GET_VOXEL_4D(v00, volume, i, ind0, ind1, ind2);
      GET_VOXEL_4D(v01, volume, i, ind0, ind1, ind2+1);
      GET_VOXEL_4D(v02, volume, i, ind0, ind1, ind2+2);
      GET_VOXEL_4D(v03, volume, i, ind0, ind1, ind2+3);
      
      GET_VOXEL_4D(v10, volume, i, ind0, ind1+1, ind2);
      GET_VOXEL_4D(v11, volume, i, ind0, ind1+1, ind2+1);
      GET_VOXEL_4D(v12, volume, i, ind0, ind1+1, ind2+2);
      GET_VOXEL_4D(v13, volume, i, ind0, ind1+1, ind2+3);
      
      GET_VOXEL_4D(v20, volume, i, ind0, ind1+2, ind2);
      GET_VOXEL_4D(v21, volume, i, ind0, ind1+2, ind2+1);
      GET_VOXEL_4D(v22, volume, i, ind0, ind1+2, ind2+2);
      GET_VOXEL_4D(v23, volume, i, ind0, ind1+2, ind2+3);
      
      GET_VOXEL_4D(v30, volume, i, ind0, ind1+3, ind2);
      GET_VOXEL_4D(v31, volume, i, ind0, ind1+3, ind2+1);
      GET_VOXEL_4D(v32, volume, i, ind0, ind1+3, ind2+2);
      GET_VOXEL_4D(v33, volume, i, ind0, ind1+3, ind2+3);
      
      CUBIC_BIVAR(v00,v01,v02,v03, \
		  v10,v11,v12,v13, \
		  v20,v21,v22,v23, \
		  v30,v31,v32,v33, \
		  frac[Y],frac[X], temp_result);
      def[i] = CONVERT_VOXEL_TO_VALUE(volume,temp_result);
    }
    
    

  }
  else {
    print ("\n\nxyzv[] = %d %d %d %d\n",xyzv[0],xyzv[1],xyzv[2],xyzv[3]);
    print ("\n\nsizes[]= %d %d %d %d\n",sizes[0],sizes[1],sizes[2],sizes[3]);
    print_error_and_line_num("error in slice_interpolate", __FILE__, __LINE__);
  }

}

