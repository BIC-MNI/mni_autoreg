#define BRUTE_FORCE     1
#define SECANT_METHOD   0
#define TEST_LENGTH     1.0

#define EQUAL(a,b) ( ABS( (a) - (b)) < 0.000001)

#include "line_data.h"

private float interp_data_along_gradient(Real dist_from, int inter_type,
					 Real x, Real y, Real z, 
					 Real dx, Real dy, Real dz,
					 Volume data);

private Real get_normalized_world_gradient_direction(Real xworld, Real yworld, Real zworld, 
						     Volume data_x,Volume data_y,Volume data_z,
						     Real *DX_dir, Real *DY_dir, Real *DZ_dir,
						     Real thresh);

private Real optimize_1D_deformation_for_single_node(Real spacing, 
						     Real threshold1, Real threshold2,
						     Real wx, Real wy, Real wz,
						     Real mx, Real my, Real mz,
						     Real *def_x, Real *def_y, Real *def_z,
						     int match_type);
private BOOLEAN build_second_derivative_data(Line_data *line, int num_samples, Real spacing,
					     Real wx, Real wy, Real wz, 
					     Real dx_dir, Real dy_dir, Real dz_dir, 
					     Volume intensity);

private BOOLEAN build_first_derivative_data(Line_data *line, int num_samples, Real spacing,
					     Real wx, Real wy, Real wz, 
					     Real dx_dir, Real dy_dir, Real dz_dir, 
					     Volume intensity);



/*******************************************************************************/
/* this routine will find the best deformation that has to be applied 
   to a single node to increase to overall objective function value.

   note:
   fwhm = 2.35*sigma
   fwtm = 4.3 *sigma

   spacing should be equal to half the fwhm value used to blur the data.


   wx,wy,wz is the position of the point in the source volume
   mx,my,mz is the position of the average neighbourhood warp, in the target volume.

*/




private Real optimize_1D_deformation_for_single_node(Real spacing, Real threshold1, Real threshold2,
						  Real wx, Real wy, Real wz,
						  Real mx, Real my, Real mz,
						  Real *def_x, Real *def_y, Real *def_z,
						  int match_type)
{
  Line_data data1, data2;
  Real
    mag_normal2,mag,
    offset, 
    result,
    tx,ty,tz, tempx, tempy, tempz, temp1x, temp1y, temp1z,
    px_dir, py_dir, pz_dir,
    d1x_dir, d1y_dir, d1z_dir,
    d2x_dir, d2y_dir, d2z_dir;
  int 
    flag1, flag2;

  result = 0.0;			/* assume no additional deformation */
  *def_x = 0.0;
  *def_y = 0.0;
  *def_z = 0.0;


  if (!get_best_start_from_neighbours(threshold1,  
				      wx,  wy,  wz,
				      mx,  my,  mz,
				      &tx, &ty,  &tz,
				      &d1x_dir, &d1y_dir, &d1z_dir,
				      def_x,  def_y,  def_z)) {
    return(result);
  }

  mag_normal2 =
    get_normalized_world_gradient_direction(tx,ty,tz,
					    Gd2_dx, Gd2_dy, Gd2_dz,
					    &d2x_dir, &d2y_dir, &d2z_dir,
					    threshold2);

                                /* project normal 1 in space of normal 2 */
    tempx = wx + d1x_dir;
    tempy = wy + d1y_dir;
    tempz = wz + d1z_dir;
    general_transform_point(Gglobals->trans_info.transformation,
                            tempx,tempy,tempz,
                            &temp1x,&temp1y,&temp1z);
    tempx = temp1x - tx;
    tempy = temp1y - ty;
    tempz = temp1z - tz;

    mag = sqrt(tempx*tempx + tempy*tempy + tempz*tempz);

    if (mag>0) {
      px_dir = tempx / mag;
      py_dir = tempy / mag;
      pz_dir = tempz / mag;
    }
    else {
      px_dir = tempx;
      py_dir = tempy;
      pz_dir = tempz;
    }

                                /* if mag_normal2 is too small, than use projection of d1?_dir
                                   in the target space for the direction vector*/
    if (mag_normal2 < threshold2) {
      d2x_dir = px_dir;
      d2y_dir = py_dir;
      d2z_dir = pz_dir;
    }
    else {                      /* make sure that the both normal directions are  aligned,
                                   and are not pointing in opposite directions!*/
      tempx = d2x_dir + px_dir;
      tempy = d2y_dir + py_dir;
      tempz = d2z_dir + pz_dir;
      mag = tempx*tempx + tempy*tempy + tempz*tempz;

      tempx = d2x_dir - px_dir;
      tempy = d2y_dir - py_dir;
      tempz = d2z_dir - pz_dir;

                                /* if the d2?_dir + p?_dir < d2?_dir - p?_dir,
                                   then d2?_dir is pointing the wrong way! x*/

      if (mag < (tempx*tempx + tempy*tempy + tempz*tempz)) {
        d2x_dir = -(d2x_dir);
        d2y_dir = -(d2y_dir);
        d2z_dir = -(d2z_dir);
      }

    }



  /* we now have direction vectors and positions for both volumes.
     Now, get samples of 2nd derivative along gradient direction 
     from both volumes. */


				/* from data1, get samples from -spacing to spacing */
  if (Gglobals->trans_info.use_magnitude==TRUE) {

    flag1 = build_first_derivative_data(&data1, 5, spacing, 
			    wx,wy,wz, d1x_dir, d1y_dir, d1z_dir, Gd1_dxyz);

				/* from data2, get samples from -2*spacing to 2*spacing */
    flag2 = build_first_derivative_data(&data2, 10, spacing, 
			    tx,ty,tz, d2x_dir, d2y_dir, d2z_dir, Gd2_dxyz);

  }
  else {
    flag1 = build_second_derivative_data(&data1, 5, spacing, 
			    wx,wy,wz, d1x_dir, d1y_dir, d1z_dir, Gd1);

				/* from data2, get samples from -2*spacing to 2*spacing */
    flag2 = build_second_derivative_data(&data2, 10, spacing, 
			    tx,ty,tz, d2x_dir, d2y_dir, d2z_dir, Gd2);
  }

  if (flag1 && flag2 && data1.count>=5 && data2.count>=10) {
    offset = find_offset_to_match(&data1,&data2,spacing,match_type);

    if (ABS(offset)>0.95*spacing) {
      offset = 0.0; 
    }
    
    
    *def_x += offset*px_dir;	/* use the projection of dir 1 */
    *def_y += offset*py_dir;
    *def_z += offset*pz_dir;

    result = offset;
	
  }
  else {
    if (!flag1) print("can't build_second_derivative_data at (v1) %f %f %f\n", wx,wy,wz);
    if (!flag2) print("can't build_second_derivative_data at (v2) %f %f %f\n", tx,ty,tz);
    if (data1.count<5) print("can't data1.count<5 at (v1) %f %f %f\n", wx,wy,wz);
    if (data2.count<9) print("can't data2.count<9 at (v2) %f %f %f\n", tx,ty,tz);

  }

  return(result);
}


/*******************************************************************************/
/*  procedure: build_second_derivative_data

    build an array of data corresponding to samples of the second derivative along the 
    direction specified by dx_dir, dy_dir and dz_dir (the normalized direction).

    This is done by extracting samples from the intensity volume and fitting a cubic 
    spline to the data, and using the cubic interpolant to estimate the second derivative.

    line = (1+num_samples*2) data samples are returned in line that span a space defined 
           from from -spacing*(num_samples-1)/4  to spacing*(num_samples-1)/4;

	   so when num_samples = 5, the data is defined from -spacing to spacing, centered
           on the coordinate wx,wy, wz


    a cubic spline is used to interpolate the second derivative from four
    function values, f1,f2,f3 & f4, where:
    f1 = f(x=-1), f2 = f(x=0), f3 = f(x=1), f4 = f(x=2), then
       fpp(x=0)   =  f1 - 2*f2 + f3      
       fpp(x=0.5) = (f1 - f2   - f3 + f4)/2  
      



*/


private BOOLEAN build_second_derivative_data(Line_data *line, int num_samples, Real spacing,
					  Real wx, Real wy, Real wz, 
					  Real dx_dir, Real dy_dir, Real dz_dir, 
					  Volume intensity)
{

  int near_edge,i,start_i, end_i;
  Real frac;
  Real data[30];
  
				/* set up Line_data structure */
  line->step = spacing/4.0;
  line->start = -num_samples * line->step;
  line->count = num_samples*2 +1;

 				/* make space for intensity data */

				/* get intensity data */
  near_edge = FALSE;
  for_inclusive(i,-num_samples-2,num_samples+2) {
    frac = ((Real)i - 0.5)*line->step ;
    data[i+num_samples+2] = interp_data_along_gradient(frac, 3,  
						       wx,wy,wz, 
						       dx_dir, dy_dir, dz_dir, intensity);
    if (data[i+num_samples+2] == -DBL_MAX) {
      near_edge = TRUE;
    }
  }

  if (near_edge) {

    for_inclusive(i,-num_samples-2,num_samples+2) {
      start_i = i;
      if (data[i+num_samples+2] != -DBL_MAX) break;
    }
    for( i = num_samples+2; i >= -num_samples-2; --i ) {
      end_i = i;
      if (data[i+num_samples+2] != -DBL_MAX) break;
    }

    if (end_i - start_i < 4) {
       return(FALSE);
    }


    if (start_i > -num_samples) {
      line->data[start_i-1+num_samples] = CUBIC_PP(data[start_i-1+num_samples+1],
						   data[start_i-1+num_samples+2],
						   data[start_i-1+num_samples+3],
						   data[start_i-1+num_samples+4], -0.5);
    }

    for_inclusive(i,start_i,end_i-3) {
      line->data[i+num_samples] = (data[i+num_samples+1] - 
				   data[i+num_samples+2] - 
				   data[i+num_samples+3] + 
				   data[i+num_samples+4])/2.0;
    }
    if (end_i < num_samples)
      line->data[end_i-3+1+num_samples] = CUBIC_PP(data[end_i-1+num_samples+1],
						   data[end_i-1+num_samples+2],
						   data[end_i-1+num_samples+3],
						   data[end_i-1+num_samples+4], 1.5);
    
    line->count = end_i - start_i + 1;
    line->start = start_i * line->step;
    
    if (start_i > -num_samples) {
      line->count++;
      line->start = line->start - line->step;
    }

    if (end_i<num_samples) 
      line->count++;

    for_less(i,0,line->count) {
      line->data[i] = line->data[start_i+num_samples];
    }

    return(TRUE);
  }
  else {
    for_inclusive(i,-num_samples,num_samples) {
      
      /* fpp(xpos=0.5) = (f1 - f2   - f3 + f4)/2   */
      
      line->data[i+num_samples] = (data[i+num_samples+1] - 
				   data[i+num_samples+2] - 
				   data[i+num_samples+3] + 
				   data[i+num_samples+4])/2.0;
    }
    return(TRUE);
  }

}

private BOOLEAN build_first_derivative_data(Line_data *line, int num_samples, Real spacing,
					    Real wx, Real wy, Real wz, 
					    Real dx_dir, Real dy_dir, Real dz_dir, 
					    Volume intensity)
{
  Real frac;
  int i;
				/* set up Line_data structure */
  line->step = spacing/4.0;
  line->start = -num_samples * line->step;
  line->count = num_samples*2 +1;

 				/* make space for intensity data */

				/* get intensity data */


  for_inclusive(i,-num_samples,num_samples) {
    frac = ((Real)i - 0.5)*line->step ;
    line->data[i+num_samples] = interp_data_along_gradient(frac,  1, 
						     wx,wy,wz, 
						     dx_dir, dy_dir, dz_dir, intensity);
  }
  return(TRUE);

}

/*******************************************************************************/
/*
  return a vector in 'DX_dir, DY_dir, DZ_dir' pointing uphill in the intensity
  field (point towards brighter areas). This vector is returned in the World 
  Coordinate system

  note: xworld, yworld, zworld  <- world coordinates
        data_x, data_y, data_z  <- derivatives in voxel coordinates! 
	                           along xspace, yspace and zspace, respectively.

        DX_dir, DY_dir, DZ_dir  <- direction of gradient magnitude in WORLD COORDS!

*/
private Real get_normalized_world_gradient_direction(Real xworld, Real yworld, Real zworld, 
						     Volume data_x,Volume data_y,Volume data_z,
						     Real *DX_dir, Real *DY_dir, Real *DZ_dir,
						     Real thresh)
     
{

  Real
    mag,mag1,
    dx_dir, dy_dir, dz_dir,
    xvox, yvox, zvox;
  
  PointR 
    voxel;

  convert_3D_world_to_voxel(data_x, xworld, yworld, zworld,
			 &zvox, &yvox, &xvox);
  fill_Point( voxel, zvox, yvox, xvox );

				/* interpolate the real values from */
				/* the volume data structures */
  if (!trilinear_interpolant(data_x, &voxel, &dx_dir)) return(0.0); 
  if (!trilinear_interpolant(data_y, &voxel, &dy_dir)) return(0.0); 
  if (!trilinear_interpolant(data_z, &voxel, &dz_dir)) return(0.0); 

  mag = sqrt(dx_dir*dx_dir + dy_dir*dy_dir + dz_dir*dz_dir);
  if (mag < thresh) {
    mag1 = (float)mag;
    return(mag1);
  }
				/* check to see if we are on a ridge line: If the first 
				   derivative changes sign on the point we are interested in,
				   then we are on an intensity ridge (in that direction) */

/*  if ((ABS(dx_dir) < thresh/20) || (ABS(dy_dir) < thresh/20) || (ABS(dz_dir) < thresh/20) )
    return(0.0);


  if (ABS(dx_dir) < thresh/10.0) {
    print ("x");
    fill_Point( voxel, zvox, yvox, xvox+1 );
    if (!trilinear_interpolant(data_x, &voxel, &d1))  return(0.0); 
    fill_Point( voxel, zvox, yvox, xvox-1 );
    if (!trilinear_interpolant(data_x, &voxel, &d2))  return(0.0); 
    if (d1 * d2 < 0.0)	{ / * if derivative changes sign, we're on a ridge * /
      if ( ABS(d1) > ABS(d2 ) )
	dx_dir = d1;	/ * so, take a value nearby in order to point at ridge *  /
      else
	dx_dir = d2;
      print ("X");
    }
  }


  if (ABS(dy_dir) < thresh/10.0) {
    print ("y");
    fill_Point( voxel, zvox, yvox+1, xvox );
    if (!trilinear_interpolant(data_y, &voxel, &d1))  return(0.0); 
    fill_Point( voxel, zvox, yvox-1, xvox );
    if (!trilinear_interpolant(data_y, &voxel, &d2))  return(0.0); 
    if (d1 * d2 < 0.0)	{ / * if derivative changes sign, we're on a ridge * / 
      if ( ABS(d1) > ABS(d2) )
	dy_dir = d1;	/ * so, take a value nearby in order to point at ridge * /
      else
	dy_dir = d2;
      print ("Y");
    }
  }
  if (ABS(dz_dir) < thresh/10.0) {
    print ("z");
    fill_Point( voxel, zvox+1, yvox, xvox );
    if (!trilinear_interpolant(data_z, &voxel, &d1))  return(0.0); 
    fill_Point( voxel, zvox-1, yvox, xvox );
    if (!trilinear_interpolant(data_z, &voxel, &d2))  return(0.0); 
    if (d1 * d2 < 0.0) {	/ * if derivative changes sign, we're on a ridge * /
      if ( ABS(d1) > ABS(d2) )
	dz_dir = d1;	/ * so, take a value nearby in order to point at ridge * /
      else
	dz_dir = d2;
      print ("Z");
    }
  }
*/				/* get normal vector in voxel coordinates: */

  mag = sqrt(dx_dir*dx_dir + dy_dir*dy_dir + dz_dir*dz_dir);
  if (mag > 0) {
    dx_dir /= mag;
    dy_dir /= mag;
    dz_dir /= mag;
  }
  else {
    dx_dir = 0.0;
    dy_dir = 0.0;
    dz_dir = 0.0;
  }

  mag1 = (float) mag;

  *DX_dir = dx_dir;
  *DY_dir = dy_dir;
  *DZ_dir = dz_dir;

  return(mag1);
}

/*******************************************************************************/
/* return the value of the gradient magnitude of data volume a distance 'x'
   from the global position Gxpf, Gypf, Gzpf along the gradient direction
   stored in DX_dir, DY_dir, DZ_dir */

private float interp_data_along_gradient(Real dist_from, int inter_type,
					 Real x, Real y, Real z, 
					 Real dx, Real dy, Real dz,
					 Volume data)  
{

  float
    value_flt;
  Real 
    value,
    xvox, yvox, zvox;
  PointR
    voxel;
  int f;

  convert_3D_world_to_voxel(data,
			    x+dist_from*dx, y+dist_from*dy, z+dist_from*dz,
			    &xvox, &yvox, &zvox);
  
  fill_Point( voxel, xvox, yvox, zvox );

  if (inter_type==3)
    f = tricubic_interpolant(data, &voxel, &value);
  else
    f = trilinear_interpolant(data, &voxel, &value);

  
  if ( f  )
    value_flt = (float)value;
  else {
    value_flt = -DBL_MAX;
  }

  return(value_flt);
}




typedef Real (*Difference_Function) 
     (Line_data *m, Line_data *d, Real offset);



/*
% return the value of the 1st or second derivative of the data, based on cubic
% interpolation through the four equally spaced data function value
% points.
% xpos is the fractional distance from the x position of f2 for which the
% value is to be returned.
function v = cubic_p(f1,f2,f3,f4,xpos)

%  for f1 = f(-1),
%      f2 = f(0),
%      f3 = f(1),
%      f4 = f(2), 
%
%      a1 =                f2;
%      a2 =  -1.0*f1/3.0 - f2/2.0 + f3     - f4/6.0;
%      a3 =       f1/2.0 - f2     + f3/2.0 ;
%      a4 =  -1.0*f1/6.0 + f2/2.0 - f3/2.0 + f4/6.0;
%
%  to estimate function f() at xpos, using cubic interpolation:
%
%      f   = a1 + a2*xpos + a3*xpos^2 + a4*xpos^3
%      fp  = a2 + 2*a3*xpos + 3*a4*xpos*xpos;      first derivitive
%      fpp = 2*a3 + 6*a4*xpos;	                   second derivitive 
*/

public Real CUBIC_PP(Real f1, Real f2, Real f3, Real f4, Real xpos) {
  Real v,a3,a4;
  
  a3 =      f1/2.0 - f2     + f3/2.0 ;
  a4 = -1.0*f1/6.0 + f2/2.0 - f3/2.0 + f4/6.0;
  
  v =  2*a3 + 6*a4*xpos;
  
  return(v);
}

public Real CUBIC_P(Real f1, Real f2, Real f3, Real f4, Real xpos) {

  Real v,a2,a3,a4;
  
  a2 = -1.0*f1/3.0 - f2/2.0 + f3     - f4/6.0;
  a3 =      f1/2.0 - f2     + f3/2.0 ;
  a4 = -1.0*f1/6.0 + f2/2.0 - f3/2.0 + f4/6.0;
  
  v =  a2 + 2*a3*xpos + 3*a4*xpos*xpos;
  
  return(v);
}



/*
  real space 2nd derivative value interpolation
*/
public Real inter_2p(Line_data *line, Real real_p) {

  Real pos, u, v, v0,v1,v2,v3;
  int p,m;
				/* establish voxel position */
  
  pos = (real_p - line->start) / line->step;
  m = line->count;
  
  p = FLOOR(pos);
  u = pos-p;
  
  v=0;
  
  if (pos>=1 &&  pos<(m-2)) {
    v0 = line->data[p-1];
    v1 = line->data[p];
    v2 = line->data[p+1];
    v3 = line->data[p+2];
    
    v = CUBIC_PP(v0,v1,v2,v3,u);
    
  } else if (pos>=0 &&  pos<1) {
    v0 = line->data[p];
    v1 = line->data[p+1];
    v2 = line->data[p+2];
    v3 = line->data[p+3];
    u=u+1;

    v = CUBIC_PP(v0,v1,v2,v3,u);

  } else if (pos>=(m-2) &&  pos<(m-1)) {
    v0 = line->data[p-2];
    v1 = line->data[p-1];
    v2 = line->data[p];
    v3 = line->data[p+1];
    u=u-1;
    
    v = CUBIC_PP(v0,v1,v2,v3,u);
  }

  return(v);
}




/* 
  real space 1st derivative value interpolation
*/ 
public Real inter_1p(Line_data *line, Real real_p) {

  Real pos, u, v, v0,v1,v2,v3;
  int m,p;
				/* establish voxel position */
  
  pos = (real_p - line->start) / line->step;
  m = line->count;
  
  p = FLOOR(pos);
  u = pos-p;
  
  v=0;
  
  if (pos>=1 &&  pos<(m-2)) {
    v0 = line->data[p-1];
    v1 = line->data[p];
    v2 = line->data[p+1];
    v3 = line->data[p+2];
    
    v = CUBIC_P(v0,v1,v2,v3,u);
    
  } else if (pos>=0 &&  pos<1) {
    v0 = line->data[p];
    v1 = line->data[p+1];
    v2 = line->data[p+2];
    v3 = line->data[p+3];
    u=u+1;

    v = CUBIC_P(v0,v1,v2,v3,u);

  } else if (pos>=(m-2) &&  pos<(m-1)) {
    v0 = line->data[p-2];
    v1 = line->data[p-1];
    v2 = line->data[p];
    v3 = line->data[p+1];
    u=u-1;
    
    v = CUBIC_P(v0,v1,v2,v3,u);
  }

  return(v);
}







/*
  voxel space value interpolation
*/
public Real inter(Line_data *line, Real pos) {

  Real v,u,v0,v1,v2,v3,f;
  int m,p;

  m = line->count;
  v = 0;

  if (pos>=1 &&  pos<(m-2)) {
    p = FLOOR(pos);
    u = pos-p;
    
    v0 = line->data[p-1];
    v1 = line->data[p];
    v2 = line->data[p+1];
    v3 = line->data[p+2];
    
    v = ( (v1) + (u) * ( 0.5 * ((v2)-(v0)) + (u) * ( (v0) - 2.5 * (v1) + 2.0 * (v2) - 0.5 * (v3) + (u) * ( -0.5 * (v0) + 1.5 * (v1) - 1.5 * (v2) + 0.5 * (v3)  ) ) ) );
  }
  else {
    if (pos<-0.5 || pos>=(m-0.5))
      v = 0;
    else {
      if (pos<0 || pos>(m-1)) {
	if (pos<0)
	  v = 0.5 * line->data[0];
	else
	  v = 0.5 * line->data[m-1];
      }
      else {
	if (pos<1) {
	  f = pos;
	  v = (1.0-f)*line->data[0] + f*line->data[1];
	} else
	  if (pos>=m-2) {
	    f = pos-(m-2);
	    v = (1.0-f)*line->data[m-2] + f*line->data[m-1];
	  }
      }
    }
  }
  
  return(v);
}

/*
  real space value interpolation
*/
public Real inter_p(Line_data *data, Real real_p) {

  Real pos,v;

  pos = (real_p - data->start) / data->step;

  v = inter(data, pos);

  return(v);
}



/*  return the difference between the overlapping area
    of the two functions of model+offset wrt data.
*/
public Real difference(Line_data *model,
		      Line_data *data,
		      Real offset) {

  Real t, d, val, minx, maxx;
  int i;
  
  val = 0;
  minx = data->start;
  maxx = data->start + data->count*data->step;

  for_less(i,0,model->count) {

    t = model->start + i*model->step + offset; /* new position */

    if (t>=minx && t <=maxx) {
      d = model->data[i] - inter_p(data, t);
      val = val + d*d;
    }
  }
  return(val);
}

/* procedure: brute
     use brute force to search for the offset that minimizes the difference
     between the data in m and d, where m is a spatial subset of d.
   method:
     search with offsets from -limit to limit, stepping by the value in step
     and return the best offset position.
*/
public Real brute(Difference_Function differ,
		  Line_data *m,
		  Line_data *d,
		  Real step,
		  Real limit)
{

  Real 
    v_min, v, pos,
    dist;

  v_min = FLT_MAX;

  dist = 0;

  pos = -limit;

  while (pos<=limit) {


    v = (*(differ))(m,d,pos);

/*
print ("%6.2f %10.8f\n",pos,v*1000);
*/
    if (v<v_min) {
      v_min = v;
      dist = pos;
    }
    pos += step;
  }

  pos = dist;

  if (dist<limit) {
    v = (*(differ))(m,d,pos+step/2);

    if (v<v_min) {
      v_min = v;
      dist = pos+step/2;
    }
  }

  if (dist> -limit) {
    v = (*(differ))(m,d,pos-step/2);
    if (v<v_min) {
      v_min = v;
      dist = pos-step/2;
    }
  }

  return(dist);
}

#define R  0.61803399
#define C  1.0-R

public Real gold(Difference_Function differ,
		 Line_data *m,
		 Line_data *d,
		 Real limit)
{

  Real /*f0,*/f1,f2,/*f3,*/d1,d2,x0,x3,x1,x2,ax,bx,cx,dist;
  
  int i;

  i = 0;
  dist = 0;

  ax = -0.6*limit;
  cx = 0.6*limit;


  
  x0 = ax;
  bx = 0;
  x3 = cx;

  d1=cx-bx;
  d2=bx=ax;

  if ( ABS(d1) > ABS(d2)) {
    x1=bx;
    x2=bx+C*(cx-bx);
  }
  else {
    x2=bx;
    x1=bx-C*(bx-ax);
  }

  f1= (*(differ))(m,d,x1);
  f2= (*(differ))(m,d,x2);

  while (abs(x3-x0) > limit/16 && i<20) {
    i++;
    if (f2 < f1) {
      x0=x1; x1=x2; x2=R*x1+C*x3;
      /* f0=f1;  */
      f1=f2; f2= (*(differ))(m,d,x2);
    }
    else {
      x3=x2; x2=x1; x1=R*x2+C*x0;
      /* f3=f2; */
      f2=f1; f1= (*(differ))(m,d,x1);
    }
  }
  
  if (f1 < f2)
    dist=x1;
  else
    dist=x2;
  
  return(dist);
  
  
}


/*
  find the  offset to best match model (m) with data (d)

  offset returned is constrained to: 
       | offset | <= limit

  if op_type = 1, do brute force, otherwise use secant
*/
public Real find_offset_to_match(Line_data *m,
				 Line_data *d,
				 Real      limit,
				 int       op_type) 
{

  Real dist,step;

  if (op_type==1) {
    step = limit/10;
    dist = brute(difference,m,d,step,limit);
  }
  else
    dist = gold(difference,m,d,limit);

  return(dist);
}

/*******************************************************************************/
/*  procedure: smooth_the_warp

    desc: this procedure will smooth the current warp stored in warp_d? and 
          return the smoothed warp in smooth_d?

    meth: smoothing is accomplished by averaging the 6 neighbour of each node
          with the value at that node.

	  new_val = FRAC1*old_val + (1-FRAC1)*sum_6_neighbours(val);
    
*/
public void smooth2_the_warp(Volume smooth_dx, Volume smooth_dy, Volume smooth_dz, 
			     Volume warp_dx, Volume warp_dy, Volume warp_dz) {

  int i,j,k,sizes[3];
  Real voxel, value, value2;
  progress_struct
    progress;

  get_volume_sizes(warp_dx, sizes);
  
  initialize_progress_report( &progress, FALSE, sizes[0]*sizes[1] + 1,
			     "Smoothing deformations" );

  for_less(i,1,sizes[0]-1)
    for_less(j,1,sizes[1]-1) {
      for_less(k,1,sizes[2]-1){
	
	GET_VOXEL_3D(voxel, warp_dx, i,j,k); 
	value = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value *= FRAC1;
	GET_VOXEL_3D(voxel, warp_dx, i+1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i-1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i,j+1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i,j-1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i,j,k+1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dx, i,j,k-1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dx, voxel ); value += value2*FRAC2; 
	voxel = CONVERT_VALUE_TO_VOXEL(smooth_dx, value);
	SET_VOXEL_3D(smooth_dx, i,j,k, voxel); 
	
	GET_VOXEL_3D(voxel, warp_dy, i,j,k); 
	value = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value *= FRAC1;
	GET_VOXEL_3D(voxel, warp_dy, i+1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i-1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i,j+1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i,j-1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i,j,k+1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dy, i,j,k-1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dy, voxel ); value += value2*FRAC2; 
	voxel = CONVERT_VALUE_TO_VOXEL(smooth_dy, value);
	SET_VOXEL_3D(smooth_dy, i,j,k, voxel); 
	
	GET_VOXEL_3D(voxel, warp_dz, i,j,k); 
	value = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value *= FRAC1;
	GET_VOXEL_3D(voxel, warp_dz, i+1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i-1,j,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i,j+1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i,j-1,k); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i,j,k+1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	GET_VOXEL_3D(voxel, warp_dz, i,j,k-1); 
	value2 = CONVERT_VOXEL_TO_VALUE( warp_dz, voxel ); value += value2*FRAC2; 
	voxel = CONVERT_VALUE_TO_VOXEL(smooth_dz, value);
	SET_VOXEL_3D(smooth_dz, i,j,k, voxel); 
	
      }
      update_progress_report( &progress, sizes[1]*i+j+1 );

    }
    terminate_progress_report( &progress );


  
}
