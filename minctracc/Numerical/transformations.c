/* ----------------------------- MNI Header -----------------------------------
@NAME       : transformations.c
@DESCRIPTION: routines to apply the forward and inverse transformations
              of the non-linear deformation field.
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

@CREATED    : Tue Nov 16 14:51:04 EST 1993 lc
                    based on transformations.c from fit_vol
@MODIFIED   : $Log: transformations.c,v $
@MODIFIED   : Revision 1.6  1994-05-28 16:18:29  louis
@MODIFIED   : working version before modification of non-linear optimiation
@MODIFIED   :
 * Revision 1.5  94/04/06  11:48:52  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.4  94/02/21  16:37:02  louis
 * version before feb 22 changes
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Numerical/transformations.c,v 1.6 1994-05-28 16:18:29 louis Exp $";
#endif


#include <volume_io.h>

#include "arg_data.h"
#include "local_macros.h"

#define NUMBER_TRIES 10


#include "deform_field.h"

#include "constants.h"

extern Arg_Data main_args;


static char *default_dim_names[N_DIMENSIONS] =
   { MIzspace, MIyspace, MIxspace };


public void
  apply_deformation_field_to_point(void *trans_data, 
				   Real x_in, Real y_in, Real z_in, 
				   Real *x_out, Real *y_out, Real *z_out)
{
   
  Real
    xv,yv,zv,
    delta_x,
    delta_y,
    delta_z;
  PointR
    voxel;
  Deform_field
    *def;

  def = (Deform_field *)trans_data;

  if (def->dx==NULL) {
    delta_x = 0.0;
    delta_y = 0.0;
    delta_z = 0.0;
  }
  else {

    convert_3D_world_to_voxel(def->dx, x_in,y_in,z_in, &xv, &yv, &zv);
    
    fill_Point( voxel, xv, yv, zv ); /* build the voxel POINT */

/*    
    if ( !INTERPOLATE_TRUE_VALUE( def->dx, &voxel, &delta_x ))
      delta_x = 0.0;
    if ( !INTERPOLATE_TRUE_VALUE( def->dy, &voxel, &delta_y ))
      delta_y = 0.0;
    if ( !INTERPOLATE_TRUE_VALUE( def->dz, &voxel, &delta_z ))
      delta_z = 0.0;
*/

    if ( !tricubic_interpolant( def->dx, &voxel, &delta_x ))
      delta_x = 0.0;
    if ( !tricubic_interpolant( def->dy, &voxel, &delta_y ))
      delta_y = 0.0;
    if ( !tricubic_interpolant( def->dz, &voxel, &delta_z ))
      delta_z = 0.0;
    
  }

  *x_out = x_in + delta_x;
  *y_out = y_in + delta_y;
  *z_out = z_in + delta_z;
  
}

/* apply the non-linear transformation inversely, mapping points 
   of volume2 into  the space of volume 1   using numerical inversion.
*/

public void
  apply_inv_deformation_field_to_point(void *trans_data, 
				       Real x_in, Real y_in, Real z_in, 
				       Real *x_out, Real *y_out, Real *z_out)

{

  Real
    xv,yv,zv,
    delta_x,
    delta_y,
    delta_z;
  PointR
    voxel;
  int
    tries;
  Real
    ftol,
    x,y,z,                                     /* xyz for first volume                   */
    forward_x, forward_y, forward_z,           /* xyz mapped forward into second volume. */
    best_x,best_y,best_z, smallest_e, 
    error, error_x,  error_y,  error_z;
  Deform_field
    *def;

  def = (Deform_field *)trans_data;
  ftol = 0.05;
  
  x = x_in;
  y = y_in;
  z = z_in;

  if (def->dx!=NULL)  {
    
    convert_3D_world_to_voxel(def->dx, x,y,z, &xv, &yv, &zv);
    
    fill_Point( voxel, xv, yv, zv ); /* build the voxel POINT */
    
    /* if there is a deformation for this point,
       subtract it to form the first guess... */
/*
    if ( !INTERPOLATE_TRUE_VALUE( def->dx, &voxel, &delta_x ))
      delta_x = 0.0;
    if ( !INTERPOLATE_TRUE_VALUE( def->dy, &voxel, &delta_y ))
      delta_y = 0.0;
    if ( !INTERPOLATE_TRUE_VALUE( def->dz, &voxel, &delta_z ))
      delta_z = 0.0;
*/
    if ( !tricubic_interpolant( def->dx, &voxel, &delta_x ))
      delta_x = 0.0;
    if ( !tricubic_interpolant( def->dy, &voxel, &delta_y ))
      delta_y = 0.0;
    if ( !tricubic_interpolant( def->dz, &voxel, &delta_z ))
      delta_z = 0.0;
       
    x -= delta_x;
    y -= delta_y;
    z -= delta_z;
  }

   				/* map the guess forward into the second volume            */

  apply_deformation_field_to_point(trans_data, x,y,z,&forward_x,&forward_y,&forward_z);

				/* test it to see how close we are in the second volume    */
  error_x = x_in - forward_x;
  error_y = y_in - forward_y;
  error_z = z_in - forward_z;
  
  /* x,y,z - stores the position in volume one of the first guess */

  tries = 0;
  
  error = smallest_e = ABS(error_x) + ABS(error_y) + ABS(error_z);
  best_x = x;
  best_y = y;
  best_z = z;
  
  while ((++tries<NUMBER_TRIES) && 
	 (smallest_e > ftol) ) {
    
    x += 0.67*error_x;
    y += 0.67*error_y;
    z += 0.67*error_z;
    
    apply_deformation_field_to_point(trans_data, x,y,z,
				     &forward_x,&forward_y,&forward_z);
    
    error_x = x_in - forward_x;
    error_y = y_in - forward_y;
    error_z = z_in - forward_z;
    
    error = ABS(error_x) + ABS(error_y) + ABS(error_z);
    
    if (error<smallest_e) {
      smallest_e = error;
      best_x = x;
      best_y = y;
      best_z = z;
    }
    
  } 
  
  *x_out = best_x;
  *y_out = best_y;
  *z_out = best_z;


/*
   printf ("t=%2d %6.1f %6.1f %6.1f (%6.1f %6.1f %6.1f) %6.1f %6.1f %6.1f %8.6f)\n",
	   tries,x1,y1,z1,
	   forward_x,forward_y,forward_z, 
	   *x_out, *y_out, *z_out
	   smallest_e);
*/

}


private General_transform *get_linear_part_of_transformation(General_transform *trans)
{
  General_transform *result,*concated,*current_lin;
  int i;

  ALLOC(result, 1);
  ALLOC(concated,1 );
  ALLOC(current_lin,1);

  create_linear_transform(result, NULL); /* start with identity */

  for_less(i,0,get_n_concated_transforms(trans)) {
    if (get_transform_type( get_nth_general_transform(trans, i-0) ) == LINEAR){

      copy_general_transform( get_nth_general_transform(trans, i-0), current_lin);
      concat_general_transforms(result, current_lin, concated);

      delete_general_transform(result);
      delete_general_transform(current_lin);
      copy_general_transform(concated, result);
      delete_general_transform(concated);

   }
  }

  return(result);
}

/*
   concatenate a zero-valued non-linear deformation field transformation 
   to the end of the globals->trans_info.transformation
   
   use the lattice information (in *globals) to build the deformation volumes.
*/
public void build_default_deformation_field(Arg_Data *globals)
{
  Volume *dx,*dy,*dz;
  Real zero, voxel[3], point[3], start2[3];
  int i,j,k;
  int count_extended[3];
  Deform_field *def;
  General_transform *trans;
/*
  int count[3];
  Real step[3];
*/


  /* build the three deformation volumes: */
  ALLOC(dx,1);
  ALLOC(dy,1);
  ALLOC(dz,1);

  *dx = create_volume(3, default_dim_names, NC_SHORT, FALSE, 0.0, 0.0);
  set_volume_voxel_range(*dx, 1.0, 32767.0);  
  set_volume_real_range(*dx, -2.0*globals->step[0], 2.0*globals->step[0] );

/*   see below for count_extended 

  for_less(i,0,3) {
    count[i] = 2*globals->count[i];
    step[i] = globals->step[i]/2.0;
  }

  set_volume_sizes(*dx, count);
  set_volume_separations(*dx,step);
*/


  set_volume_sizes(*dx, globals->count);
  set_volume_separations(*dx,globals->step);

  voxel[0] = 0.0;
  voxel[1] = 0.0;
  voxel[2] = 0.0;
  if (globals->smallest_vol==1 || globals->trans_info.transform_type==TRANS_NONLIN) { /* lattice was defined on volume 1 */
    set_volume_translation(*dx, voxel, globals->start);
  } 
  else {			/* lattice was defined on volume 2, but it must be 
				   mapped back into volume 1, using only the linear
				   part of the transformation */
    
    trans = get_linear_part_of_transformation(globals->trans_info.transformation);
    general_inverse_transform_point(trans,
				    globals->start[0], globals->start[1], globals->start[2],
				    &start2[0], &start2[1], &start2[2]);
    set_volume_translation(*dx, voxel, start2);

    general_inverse_transform_point(trans,
				    globals->start[0]+1, globals->start[1], globals->start[2],
				    &point[0], &point[1], &point[2]);
    for_less(i,0,3) point[i] -= start2[i];
    set_volume_direction_cosine(*dx, 2, point);
    
    general_inverse_transform_point(trans,
				    globals->start[0], globals->start[1]+1, globals->start[2],
				    &point[0], &point[1], &point[2]);
    for_less(i,0,3) point[i] -= start2[i];
    set_volume_direction_cosine(*dx, 1, point);
    
    general_inverse_transform_point(trans,
				    globals->start[0], globals->start[1], globals->start[2]+1,
				    &point[0], &point[1], &point[2]);
    for_less(i,0,3) point[i] -= start2[i];
    set_volume_direction_cosine(*dx, 0, point);

  }

				/* now make the volume four voxels larger, all around */

  convert_3D_voxel_to_world(*dx, -4.0, -4.0, -4.0, &point[X], &point[Y], &point[Z]);

  for_less(i,0,3)
    count_extended[i] = globals->count[i]+8;
  set_volume_sizes(*dx, count_extended);

  voxel[0] = 0.0;
  voxel[1] = 0.0;
  voxel[2] = 0.0;
				/* reset the first voxel position */
  set_volume_translation(*dx, voxel, point);

  *dy = copy_volume_definition(*dx, NC_UNSPECIFIED, FALSE, 0.0, 0.0);
  *dz = copy_volume_definition(*dx, NC_UNSPECIFIED, FALSE, 0.0, 0.0);
    
  alloc_volume_data(*dx);
  alloc_volume_data(*dy);
  alloc_volume_data(*dz);

				/* set them to zero */

  zero = CONVERT_VALUE_TO_VOXEL(*dx, 0.0);

  for_less(i,0,count_extended[0])
    for_less(j,0,count_extended[1])
      for_less(k,0,count_extended[2]){
	SET_VOXEL_3D(*dx, i,j,k, zero);
	SET_VOXEL_3D(*dy, i,j,k, zero);
	SET_VOXEL_3D(*dz, i,j,k, zero);
      }
  
				/* put the volumes into the Deform_field */

  ALLOC(def, 1);

  (void)strcpy(def->basename,globals->filenames.output_trans);
  (void)strcpy(def->xfmname,"\0");
  def->dx = *dx;
  def->dy = *dy;
  def->dz = *dz;
  
  ALLOC(trans,1);
  
  create_user_transform(trans, def, sizeof(*def),
			apply_deformation_field_to_point,
			apply_inv_deformation_field_to_point);
  
  concat_general_transforms(globals->trans_info.transformation, trans, 
			    globals->trans_info.transformation);
}
