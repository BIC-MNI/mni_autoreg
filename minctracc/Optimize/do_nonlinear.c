/* ----------------------------- MNI Header -----------------------------------
@NAME       : do_nonlinear.c
@DESCRIPTION: routines to do non-linear registration by local linear
              correlation.
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

@CREATED    : Thu Nov 18 11:22:26 EST 1993 LC
@MODIFIED   : $Log: do_nonlinear.c,v $
@MODIFIED   : Revision 1.1  1994-02-21 16:33:41  louis
@MODIFIED   : Initial revision
@MODIFIED   :
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/do_nonlinear.c,v 1.1 1994-02-21 16:33:41 louis Exp $";
#endif

#include <volume_io.h>
#include "arg_data.h"
#include "deform_field.h"


static Volume   Gd1;
static Volume   Gd1_dx; 
static Volume   Gd1_dy; 
static Volume   Gd1_dz; 
static Volume   Gd1_dxyz;
static Volume   Gd2;
static Volume   Gd2_dx; 
static Volume   Gd2_dy; 
static Volume   Gd2_dz; 
static Volume   Gd2_dxyz;
static Volume   Gm1;
static Volume   Gm2; 
static Arg_Data *Gglobals;

#define BAD_VALUE -10000000

/* prototypes */

private Real find_offset_to_match(Real data1[],Real data2[]);

private Real optimize_deformation_for_single_node(Real sigma, Real threshold,
						  Real wx, Real wy, Real wz,
						  Real *def_x, Real *def_y, Real *def_z);

private Real get_normalized_world_gradient_direction(Real xworld, Real yworld, Real zworld, 
						     Volume data_x,Volume data_y,Volume data_z,
						     Real *DX_dir, Real *DY_dir, Real *DZ_dir);

private float interp_data_along_gradient(Real dist_from, Real x, Real y, Real z, 
					 Real dx, Real dy, Real dz,
					 Volume data);

public Status do_non_linear_optimization(Volume d1,
					 Volume d1_dx, 
					 Volume d1_dy, 
					 Volume d1_dz, 
					 Volume d1_dxyz,
					 Volume d2,
					 Volume d2_dx, 
					 Volume d2_dy, 
					 Volume d2_dz, 
					 Volume d2_dxyz,
					 Volume m1,
					 Volume m2, 
					 Arg_Data *globals)
{
  Volume            dx,dy,dz;
  Deform_field *def_data;
  Real steps[3];
  int i,j,k,sizes[3];
  Real voxel, def_x, def_y, def_z, wx,wy,wz;
  Real 
    threshold,
    result;
  progress_struct
    progress;

				/* set up globals for comminication with other routines */
  Gd1     = d1;
  Gd1_dx  = d1_dx; 
  Gd1_dy  = d1_dy; 
  Gd1_dz  = d1_dz; 
  Gd1_dxyz= d1_dxyz;
  Gd2     = d2;
  Gd2_dx  = d2_dx; 
  Gd2_dy  = d2_dy; 
  Gd2_dz  = d2_dz; 
  Gd2_dxyz= d2_dxyz;
  Gm1     = m1;
  Gm2     = m2; 
  Gglobals= globals;
  
				/* copy current deformation for this iteration */

  def_data = (Deform_field *)NULL;
  for_less(i,0,get_n_concated_transforms(globals->trans_info.transformation))
    if (get_transform_type( get_nth_general_transform(globals->trans_info.transformation,i) ) 
	       == USER_TRANSFORM)
      def_data = (Deform_field *)get_nth_general_transform(globals->trans_info.transformation,i)->user_data;

  if (def_data == (Deform_field *)NULL) {
    print_error("Cannot find the deformation field transform to optimize",
		__FILE__, __LINE__);
  }

  dx = copy_volume_definition(def_data->dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  dy = copy_volume_definition(def_data->dy, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  dz = copy_volume_definition(def_data->dz, NC_UNSPECIFIED, FALSE, 0.0,0.0);

  alloc_volume_data(dx);
  alloc_volume_data(dy);
  alloc_volume_data(dz);

  get_volume_sizes(dx, sizes);
  get_volume_separations(dx, steps);

  threshold = 0.15 * get_volume_real_max(dx);

  initialize_progress_report( &progress, FALSE, sizes[0]*sizes[1] + 1,
			     "Estimating deformations" );
 
  for_less(i,0,sizes[0]) {
    for_less(j,0,sizes[1]) {
      for_less(k,0,sizes[2]){

	convert_3D_voxel_to_world(def_data->dx, (Real)i, (Real)j, (Real)k,
				  &wx, &wy, &wz);

	/* check masking here! */

	GET_VALUE_3D( def_x, def_data->dx, i,j,k );
	GET_VALUE_3D( def_y, def_data->dy, i,j,k );
	GET_VALUE_3D( def_z, def_data->dz, i,j,k );

	result = optimize_deformation_for_single_node(steps[0], threshold,
						      wx,wy,wz,
						      &def_x, &def_y, &def_z);
	def_z = result;

	def_x = CONVERT_VALUE_TO_VOXEL(dx, def_x); 
	SET_VOXEL_3D(dx, i,j,k, def_x);
	def_y = CONVERT_VALUE_TO_VOXEL(dy, def_y); 
	SET_VOXEL_3D(dy, i,j,k, def_y);
	def_z = CONVERT_VALUE_TO_VOXEL(dz, def_z); 
	SET_VOXEL_3D(dz, i,j,k, def_z);
	
      }
      update_progress_report( &progress, sizes[1]*i+j+1 );
    }
  }
  terminate_progress_report( &progress );

				/* copy data from temporary volume, back into
				   complete general transformation */
  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]){
	GET_VOXEL_3D(voxel, dx, i,j,k); 
	SET_VOXEL_3D(def_data->dx, i,j,k, voxel); 
	GET_VOXEL_3D(voxel, dy, i,j,k); 
	SET_VOXEL_3D(def_data->dy, i,j,k, voxel); 
	GET_VOXEL_3D(voxel, dz, i,j,k); 
	SET_VOXEL_3D(def_data->dz, i,j,k, voxel); 
      }
  

  return (OK);
}

/* this routine will find the best deformation that has to be applied 
   to a single node to increase to overall objective function value.

   note:
   fwhm = 2.35*sigma
   fwtm = 4.3 *sigma
*/

private Real optimize_deformation_for_single_node(Real spacing, Real threshold,
						  Real wx, Real wy, Real wz,
						  Real *def_x, Real *def_y, Real *def_z)
{
  Real 
    offset,
    data1[9], data2[9],
    frac,
    result,mag_normal1,mag_normal2,
    mag,
    tx,ty,tz, tempx,tempy,tempz,temp1x,temp1y,temp1z,
    d1x_dir, d1y_dir, d1z_dir,
    d2x_dir, d2y_dir, d2z_dir;

  int i;

  result = 0.0;			/* assume no additional deformation */


  mag_normal1 = 
    get_normalized_world_gradient_direction(wx,wy,wz, 
					    Gd1_dx, Gd1_dy, Gd1_dz, 
					    &d1x_dir, &d1y_dir, &d1z_dir);


  if (mag_normal1 < threshold) {
				/* if the mag is too small, then the direction is
				   unreliable, so try to get a direction from the
				   target volume */

    general_transform_point(Gglobals->trans_info.transformation, wx,wy,wz, &tx,&ty,&tz);
    mag_normal2 = 
      get_normalized_world_gradient_direction(tx,ty,tz, 
					      Gd2_dx, Gd2_dy, Gd2_dz, 
					      &d2x_dir, &d2y_dir, &d2z_dir);
    if (mag_normal2 < threshold)
      return(result);		/* RETURN since no directions can be found! */

				/* since the first mag_normal is too small,  project d2?_dir
				   into the source space, and use that for a direction */
    tempx = tx+d2x_dir;
    tempy = ty+d2y_dir;
    tempz = tz+d2z_dir;
    general_inverse_transform_point(Gglobals->trans_info.transformation,
				    tempx,tempy,tempz,
				    &temp1x,&temp1y,&temp1z);
    tempx = temp1x - wx;
    tempy = temp1y - wy;
    tempz = temp1z - wz;
    
    mag = sqrt(tempx*tempx + tempy*tempy + tempz*tempz);
    
    if (mag>0) {
      d1x_dir = tempx / mag;
      d1y_dir = tempy / mag;
      d1z_dir = tempz / mag;
    } 
    else {
      d1x_dir = tempx;
      d1y_dir = tempy;
      d1z_dir = tempz;
    }
      
  }
  else {
				/* first direction is fine, get its equivalent in the 
				   other volume */

    general_transform_point(Gglobals->trans_info.transformation, wx,wy,wz, &tx,&ty,&tz);

    mag_normal2 = 
      get_normalized_world_gradient_direction(tx,ty,tz, 
					      Gd2_dx, Gd2_dy, Gd2_dz, 
					      &d2x_dir, &d2y_dir, &d2z_dir);

				/* if mag_normal is too small, than project d1?_dir
				   into the target space, and use that for a direction */
    if (mag_normal2 < threshold) {
      tempx = wx+d1x_dir;
      tempy = wy+d1y_dir;
      tempz = wz+d1z_dir;
      general_transform_point(Gglobals->trans_info.transformation,
			      tempx,tempy,tempz,
			      &temp1x,&temp1y,&temp1z);
      tempx = temp1x - tx;
      tempy = temp1y - ty;
      tempz = temp1z - tz;

      mag = sqrt(tempx*tempx + tempy*tempy + tempz*tempz);

      if (mag>0) {
	d2x_dir = tempx / mag;
	d2y_dir = tempy / mag;
	d2z_dir = tempz / mag;
      } 
      else {
	d2x_dir = tempx;
	d2y_dir = tempy;
	d2z_dir = tempz;
      }
      
    }
  }

  /* we now have direction vectors and positions for both volumes.
     Now, get samples of 2nd derivative along gradient direction 
     from both volumes. */

  for_less(i,-4,5) {
    frac = (Real)i*spacing/4.0;
    data1[i+3] = interp_data_along_gradient(frac, wx,wy,wz, d1x_dir, d1y_dir, d1z_dir, Gd1);
  }
				/* do the same for the target data set */
  for_less(i,-4,5) {
    frac = (Real)i*spacing/4.0;
    data2[i+3] = interp_data_along_gradient(frac, tx,ty,tz, d2x_dir, d2y_dir, d2z_dir, Gd2);
  }


  offset = find_offset_to_match(data1,data2);
  result = offset*spacing/4.0;

  

  *def_x += offset*d2x_dir;
  *def_y += offset*d2y_dir;
  *def_z += offset*d2z_dir;

  result = 1.0;
  return(result);
}

private Real find_offset_to_match(Real data1[],Real data2[])
{
  Real blur1[9],blur2[9];
  int i, num_zeros;
  
  for_less(i,1,8) {
    blur1[i] = (data1[i-1] + data1[i] + data1[i+1]) / 3.0;
    blur2[i] = (data2[i-1] + data2[i] + data2[i+1]) / 3.0;
  }
  blur1[0] = (data1[0]+data1[1])/2.0;
  blur1[8] = (data1[7]+data1[8])/2.0;
  blur2[0] = (data2[0]+data2[1])/2.0;
  blur2[8] = (data2[7]+data2[8])/2.0;

/*
  sign = blur1[0] > 0;

  for_less(i,1,8) {
*/
}


/*
  return a vector in 'DX_dir, DY_dir, DZ_dir' pointing uphill in the intensity
  field (point towards brighter areas)
*/
private Real get_normalized_world_gradient_direction(Real xworld, Real yworld, Real zworld, 
						     Volume data_x,Volume data_y,Volume data_z,
						     Real *DX_dir, Real *DY_dir, Real *DZ_dir)
     
{

  Real
    dist,mag1,
    dx_dir, dy_dir, dz_dir,
    xvox, yvox, zvox;
  
  PointR 
    voxel,
    voxel2;

  convert_3D_world_to_voxel(data_x, xworld, yworld, zworld,
			 &xvox, &yvox, &zvox);
  fill_Point( voxel, xvox, yvox, zvox );
  
  trilinear_interpolant(data_x, &voxel, &dx_dir); /* interpolate the real values from */
  trilinear_interpolant(data_y, &voxel, &dy_dir); /* the volume data structures */
  trilinear_interpolant(data_z, &voxel, &dz_dir); 


				/* get normal vector in voxel coordinates: */

  dist = sqrt(dx_dir*dx_dir + dy_dir*dy_dir + dz_dir*dz_dir);
  if (dist > 0) {
    dx_dir /= dist;
    dy_dir /= dist;
    dz_dir /= dist;
  }
  else {
    dx_dir = 0.0;
    dy_dir = 0.0;
    dz_dir = 0.0;
  }

  mag1 = (float) dist;

  if (dist > 0.0) {

			 /*----- now, get normal vector in world coordinates: */

				/*                             --->           */
				/* get a voxel point, exactly  grad vec  away */
    fill_Point( voxel2, xvox+dx_dir, yvox+dy_dir, zvox+dz_dir );
    
				/* get this point in world coords             */
    
    convert_3D_voxel_to_world(data_x, Point_x(voxel2), Point_y(voxel2), Point_z(voxel2), 
			   &dx_dir,&dy_dir,&dz_dir);
    
    dx_dir = dx_dir-xworld;	/* now calculate the grad vec in world coords */
    dy_dir = dy_dir-yworld;
    dz_dir = dz_dir-zworld;

				/* normalize it in the world coords           */
    dist = sqrt(dx_dir*dx_dir + dy_dir*dy_dir + dz_dir*dz_dir);
    dx_dir /= dist;
    dy_dir /= dist;
    dz_dir /= dist;

  }

  *DX_dir = -dx_dir;
  *DY_dir = -dy_dir;
  *DZ_dir = -dz_dir;

  return(mag1);
}

/* return the value of the gradient magnitude of data volume a distance 'x'
   from the global position Gxpf, Gypf, Gzpf along the gradient direction
   stored in DX_dir, DY_dir, DZ_dir */

private float interp_data_along_gradient(Real dist_from, Real x, Real y, Real z, 
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

  convert_3D_world_to_voxel(data,
			    x+dist_from*dx, y+dist_from*dy, z+dist_from*dz,
			    &xvox, &yvox, &zvox);
  
  fill_Point( voxel, xvox, yvox, zvox );
  
  if (trilinear_interpolant(data, &voxel, &value) )
    value_flt = (float)value;
  else {
    value_flt = BAD_VALUE;
  }

  return(value_flt);
}

