#include <internal_volume_io.h>
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

/* ----------------------------- MNI Header -----------------------------------
@NAME       : interpolate_volume
@INPUT      : n_dims        - number of dimensions to interpolate
              parameters[]  - 0 to 1 parameters for each dim
              n_values      - number of values to interpolate
              degree        - degree of interpolation, 4 == cubic
              coefs         - [degree*degree*degree... *n_values] coeficients
@OUTPUT     : values        - pass back values
              first_deriv   - pass first derivs [n_values][n_dims]
              second_deriv  - pass back values  [n_values][n_dims][n_dims]
@RETURNS    : 
@DESCRIPTION: Computes the interpolation of the box specified by coefs and
              its derivatives, if necessary.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Mar. 20, 1995    David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

private  void   my_interpolate_volume(
    int      n_dims,
    Real     parameters[],
    int      n_values,
    int      degree,
    Real     coefs[],
    Real     values[])
{
    int       v    ;
    Real      *derivs;

    /*--- evaluate the interpolating spline */
    
    
    ALLOC( derivs, n_values );

    evaluate_interpolating_spline( n_dims, parameters, degree, n_values, coefs,
                                   1 , derivs );

    for_less( v, 0, n_values )
            values[v] = derivs[v];

    FREE(derivs);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : evaluate_volume
@INPUT      : volume
              voxel
              interpolating_dimensions - whether each dimension is interpolated
              degrees_continuity
              use_linear_at_edge
              outside_value
@OUTPUT     : values
              first_deriv
              second_deriv
@RETURNS    : 
@DESCRIPTION: Takes a voxel space position and evaluates the value within
              the volume by nearest_neighbour, linear, quadratic, or
              cubic interpolation. degrees_continuity == 2 corresponds to
              cubic, 1 for quadratic, etc.
              If first_deriv is not a null pointer, then the first derivatives
              are passed back.  Similarly for the second_deriv.
              If use_linear_at_edge is TRUE, then near the boundaries, either
              linear or nearest neighbour interpolation is used, even if cubic
              is specified by the degrees_continuity.
              If use_linear_at_edge is FALSE, then the 'outside_value' is used
              to provide coefficients for outside the volume, and the degree
              specified by degrees_continuity is used.

              Each dimension may or may not be interpolated, specified by the
              interpolating_dimensions parameter.  For instance, a 4D volume
              of x,y,z,RGB may be interpolated in 3D (x,y,z) for each of the
              3 RGB components, with one call to evaluate_volume.
@CREATED    : Mar   1993           David MacDonald
@MODIFIED   : 
---------------------------------------------------------------------------- */

public  int   my_evaluate_volume(
    Volume         volume,
    Real           voxel[],
    BOOLEAN        interpolating_dimensions[],
    int            degrees_continuity,
    BOOLEAN        use_linear_at_edge,
    Real           outside_value,
    Real           values[],
    int            n_dims,
    int            sizes[] )
{
  int      inc0, inc1, inc2, inc3, inc4, inc[MAX_DIMENSIONS];
  int      ind0, spline_degree;
  int      start0, start1, start2, start3, start4;
  int      end0, end1, end2, end3, end4;
  int      v0, v1, v2, v3, next_d;
  int      n, v, d, n_values;
  int      start[MAX_DIMENSIONS], n_interp_dims;
  int      end[MAX_DIMENSIONS];
  int      interp_dims[MAX_DIMENSIONS];
  int      n_coefs;
  Real     fraction[MAX_DIMENSIONS], bound, *coefs, pos;
  BOOLEAN  fully_inside, fully_outside;
  
  /*--- check if the degrees continuity is between nearest neighbour
    and cubic */
  
  if( degrees_continuity < -1 || degrees_continuity > 2 )
    {
      print_error( "Warning: evaluate_volume(), degrees invalid: %d\n",
		  degrees_continuity );
      degrees_continuity = 0;
    }
    
  bound = (Real) degrees_continuity / 2.0;
  
  /*--- if we must use linear interpolation near the boundaries, then
        check if we are near the boundaries, and adjust the 
	degrees_continuity accordingly */

  if( use_linear_at_edge ) {
    for_less( d, 0, n_dims )  {
      if( interpolating_dimensions == NULL || interpolating_dimensions[d]) {
	while( degrees_continuity >= 0 &&
	      (voxel[d] < bound  ||
	       voxel[d] > (Real) sizes[d] - 1.0 - bound) ) {
	  --degrees_continuity;
	  if( degrees_continuity == 1 )
	    degrees_continuity = 0;
	  bound = (Real) degrees_continuity / 2.0;
	}
      }
    }
  }

  /*--- now check which dimensions are being interpolated.  Also, compute
        how many values must be interpolated, which are all the values not
        in the interpolated dimensions */

  n_interp_dims = 0;
  n_values = 1;
  n_coefs = 1;
  spline_degree = degrees_continuity + 2;
  
  fully_inside = TRUE;
  fully_outside = FALSE;
  
  for_less( d, 0, n_dims ) {
    if( interpolating_dimensions == NULL || interpolating_dimensions[d]) {
      interp_dims[n_interp_dims] = d;
      pos = voxel[d] - bound;
      start[d] =       FLOOR( pos );
      fraction[n_interp_dims] = pos - start[d];
      
      if( voxel[d] == (Real) sizes[d] - 1.0 - bound ) {
	--start[d];
	fraction[n_interp_dims] = 1.0;
      }

      end[d] = start[d] + spline_degree;
      n_coefs *= spline_degree;

      if( start[d] < 0 || end[d] > sizes[d] ) {

	fully_inside = FALSE;
	
	if( end[d] <= 0 || start[d] >= sizes[d] ) {
	  
	  fully_outside = TRUE;
	  break;
	}
      }

      ++n_interp_dims;
    }
    else
      n_values *= sizes[d];
    }

    /*--- check if the values are uncomputable, i.e., outside volume */

  if( fully_outside ) {

    if( values != NULL ) {
      
      for_less( v, 0, n_values )
	values[v] = outside_value;
    }

    return( n_values );
  }

  /*--- add the non-interpolated dimensions to the list of dimensions, in
        order, after the interpolated dimensions */

  n = 0;
  for_less( d, 0, n_dims ) {
    
    if( interpolating_dimensions != NULL && !interpolating_dimensions[d] ) {
      
      interp_dims[n_interp_dims+n] = d;
      start[d] = 0;
      end[d] = sizes[d];
      ++n;
    }
  }

  /*--- make room for the coeficients */

  ALLOC( coefs, n_values * n_coefs );

  /*--- compute the increments in the coefs[] array for each dimension,
        in order to simulate a multidimensional array with a single dim
        array, coefs */

  inc[interp_dims[n_dims-1]] = 1;
  for_down( d, n_dims-2, 0 ) {
    next_d = interp_dims[d+1];
    inc[interp_dims[d]] = inc[next_d] * (end[next_d] - start[next_d]);
  }

  /*--- figure out the offset within coefs.  If we are inside, the offset
        is zero, since all coefs must be filled in.  If we are partially
        inside, set the offset to the first coef within the volume. */

  ind0 = 0;

  if( !fully_inside ) {

    for_less( d, 0, n_dims ) {
      
      if( start[d] < 0 ) {
	ind0 += -start[d] * inc[d];
	start[d] = 0;
      }

      if( end[d] > sizes[d] )
	end[d] = sizes[d];
    }

    for_less( v, 0, n_values * n_coefs )
      coefs[v] = outside_value;

    /*--- adjust the inc stride for stepping through coefs to account
          for the additions of the inner loops */
    
    for_less( d, 0, n_dims-1 )
      inc[d] -= inc[d+1] * (end[d+1] - start[d+1]);
  }
  else {

        /*--- adjust the inc stride for stepping through coefs to account
              for the additions of the inner loops */

    for_less( d, 0, n_dims-1 )
      inc[d] -= inc[d+1] * spline_degree;
  }

  /*--- for speed, use non-array variables for the loops */

  start0 = start[0];
  start1 = start[1];
  start2 = start[2];
  start3 = start[3];
  start4 = start[4];
  
  end0 = end[0];
  end1 = end[1];
  end2 = end[2];
  end3 = end[3];
  end4 = end[4];
  
  inc0 = inc[0];
  inc1 = inc[1];
  inc2 = inc[2];
  inc3 = inc[3];
  inc4 = inc[4];

  /*--- get the coefs[] from the volume.  For speed, do each dimension
        separately */

  for_less( v0, start0, end0 ) {
    for_less( v1, start1, end1 ) {
       for_less( v2, start2, end2 ) {
	 for_less( v3, start3, end3 ) {

	   GET_VALUE_4D( coefs[ind0], volume, v0, v1, v2, v3 );
	   ind0 += inc3;
	 }
	 ind0 += inc2;
       }
       ind0 += inc1;
     }
    ind0 += inc0;
  }

  /*--- now that we have the coeficients, do the interpolation */

  switch( degrees_continuity ) {

  case -1:                        /*--- nearest neighbour interpolation */
    for_less( v, 0, n_values )
      values[v] = coefs[v];
    break;

  case 0:
  case 1:
  case 2:

/*
  if (n_values != 3) 
      print ("n_values= %d\n",n_values);
    evaluate_interpolating_spline(n_interp_dims, fraction, spline_degree, 
				  n_values, coefs,
				  1 , values );
*/
    my_interpolate_volume(n_interp_dims, fraction, n_values,
			  spline_degree, coefs,
			  values );
    
  break;
  }

  FREE( coefs );

  return( n_values );
}

public  void   my_evaluate_volume_in_world(
    Volume         volume,
    Real           x,
    Real           y,
    Real           z,
    int            degrees_continuity,
    BOOLEAN        use_linear_at_edge,
    Real           outside_value,
    Real           values[],
    int            n_dims,
    int            sizes[])
{
  Real      voxel[MAX_DIMENSIONS];
  int       d, n_values;
  BOOLEAN   interpolating_dimensions[MAX_DIMENSIONS];
  
  /*--- convert the world space to a voxel coordinate */
  
  convert_world_to_voxel( volume, x, y, z, voxel );

  /*--- initialize all dimensions to not being interpolated */
  
  for_less( d, 0, n_dims )
    interpolating_dimensions[d] = FALSE;
  
  /*--- set each spatial dimension to being interpolated */
  
  for_less( d, 0, N_DIMENSIONS )
    interpolating_dimensions[volume->spatial_axes[d]] = TRUE;
  
  /*--- evaluate the volume in voxel space */
  
  n_values = my_evaluate_volume( volume, voxel, interpolating_dimensions,
				degrees_continuity, use_linear_at_edge, 
				outside_value,
				values, n_dims, sizes);
  
}





