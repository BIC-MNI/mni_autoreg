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
@MODIFIED   : Revision 1.9  1995-02-22 08:56:06  louis
@MODIFIED   : Montreal Neurological Institute version.
@MODIFIED   : compiled and working on SGI.  this is before any changes for SPARC/
@MODIFIED   : Solaris.
@MODIFIED   :
 * Revision 1.8  94/06/06  09:37:56  louis
 * modifed the voxel and real range calls, in build_default_deformation_field
 * where the new voxel range: 0.0, 32767.0, new real range: -50.0, 50.0.
 * 
 * These ranges are copied for all other deformation fields created from
 * the 1st scale deformation field.
 * 
 * Revision 1.7  94/06/02  20:15:59  louis
 * made modifications to allow deformations to be calulated in 2D on slices. 
 * changes had to be made in set_up_lattice, init_lattice when defining
 * the special case of a single slice....
 * Build_default_deformation_field also had to reflect these changes.
 * do_non-linear-optimization also had to check if one of dimensions had
 * a single element.
 * All these changes were made, and slightly tested.  Another type of
 * deformation strategy will be necessary (to replace the deformation 
 * perpendicular to the surface, since it does not work well).
 * 
 * Revision 1.6  94/05/28  16:18:29  louis
 * working version before modification of non-linear optimiation
 * 
 * Revision 1.5  94/04/06  11:48:52  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.4  94/02/21  16:37:02  louis
 * version before feb 22 changes
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Numerical/transformations.c,v 1.9 1995-02-22 08:56:06 louis Exp $";
#endif


#include <volume_io.h>

#include "arg_data.h"
#include "local_macros.h"

#include "deform_field.h"

#include "constants.h"

#define NUMBER_TRIES 10
#define MY_MAX_VOX 32766.0
#define MY_MAX_REAL 50.0

public void set_up_lattice(Volume data, 
			    double *user_step, /* user requested spacing for lattice */
			    double *start,     /* world starting position of lattice */
			    int    *count,     /* number of steps in each direction */
			    double *step,      /* step size in each direction */
			    VectorR directions[]);/* array of vector directions for each index*/

public int tricubic_interpolant(Volume volume, 
                                PointR *coord, double *result);

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
  This routine will do one of two tasks.  If there is NO non-linear deformation
  (as defined by a USER_TRANSFORM) in the globals->trans_info.transformation,
  then the routine will:
     concatenate a zero-valued non-linear deformation field transformation 
     to the end of the globals->trans_info.transformation, using the lattice 
     information (in *globals) to build the deformation volumes.
  If there is already a deformation stored, then 
     this deformation will be replaced by its equivalent, by resampling 
     each of the deformation field volumes using the information in *globals 
     to get the STEPping information.
     The new deformation volume will then replace those that were used in the
     general_transform_point routines.

*/
public void build_default_deformation_field(Arg_Data *globals)
{
  Volume *dx,*dy,*dz;
  Real zero, voxel_val, voxel[3], point[3], start2[3];
  int i,j,k;
  int count_extended[3];
  Deform_field *def;
  General_transform *trans, *non_lin_part;
  int count[3];
  Real start[3],step[3], step2[3];
  VectorR directions[3];
  PointR voxel_point; 
  Real del_x, del_y, del_z, nx,ny,nz, wx, wy,wz;
  progress_struct
    progress;
  int n_transforms;

  ALLOC(dx,1);
  ALLOC(dy,1);
  ALLOC(dz,1);
  
  n_transforms = get_n_concated_transforms(globals->trans_info.transformation);
  if (n_transforms==1) {


    /* build the three deformation volumes: */
    
    *dx = create_volume(3, default_dim_names, NC_SHORT, TRUE, 0.0, 0.0);
    set_volume_voxel_range(*dx, -MY_MAX_VOX, MY_MAX_VOX);
    set_volume_real_range(*dx, -MY_MAX_REAL, MY_MAX_REAL);

    set_volume_sizes(*dx, globals->count);
    set_volume_separations(*dx,globals->step);
    
    voxel[0] = 0.0;
    voxel[1] = 0.0;
    voxel[2] = 0.0;
    if (globals->smallest_vol==1 || 
	globals->trans_info.transform_type==TRANS_NONLIN) { /* lattice was defined on volume 1 */
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
    
    for_less(i,0,3) {
      if (globals->count[i]>1) {
	voxel[i] = -2.5;
	count_extended[i] = globals->count[i]+5;
      }
      else {
	voxel[i] = 0.0;
	count_extended[i] = 1;
      }
	
    }

    if (main_args.flags.debug) {
      print ("in build def transformation count = %d %d %d\n",count_extended[0],count_extended[1],count_extended[2]);
    }


    convert_3D_voxel_to_world(*dx, voxel[0], voxel[1], voxel[2], &point[X], &point[Y], &point[Z]);

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
  else {			/* we are starting with concatenated transforms,
				   so go get a point to the linear part....*/

      trans = get_linear_part_of_transformation(globals->trans_info.transformation);

				/* go get the pointer to the non-linear part */

      non_lin_part = get_nth_general_transform(
			  globals->trans_info.transformation,
                          n_transforms-1); 

      if (get_transform_type( non_lin_part ) == USER_TRANSFORM)
	def = (Deform_field *)non_lin_part->user_data;
      else {
	for_less(i,0,n_transforms)
	  print ("Transform %d is of type %d\n",i, 
		 get_transform_type(
				    get_nth_general_transform(
							      globals->trans_info.transformation,
							      i) ));

	print_error("Cannot find the deformation field transform to resample",
		    __FILE__, __LINE__);
      }



      *dx = copy_volume_definition(def->dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);
      set_volume_voxel_range(*dx, -MY_MAX_VOX, MY_MAX_VOX);  
      set_volume_real_range(*dx, -MY_MAX_REAL, MY_MAX_REAL);

				/* now redefine the new volume, using the new STEP
				   values stored in *globals->step along with the 
				   existing start and count */

      for_less(i,0,3) {
	if(globals->count[i]==1) 
	  step2[i] = 1000;
	else
	  step2[i] = globals->step[i];
      }
      
      set_up_lattice(*dx, step2, start, count, step, directions);

      for_less(i,0,3)		/* use the sign of the step returned to set the true step size */
	if (step[i]<0) step[i] = -1.0 * ABS(globals->step[i]); else step[i] = ABS(globals->step[i]);
      
      set_volume_sizes(*dx, count);
      set_volume_separations(*dx,step);
      voxel[0] = 0.0;  voxel[1] = 0.0;  voxel[2] = 0.0;
      set_volume_translation(*dx, voxel, start);

      *dy = copy_volume_definition(*dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);
      *dz = copy_volume_definition(*dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);

      alloc_volume_data(*dx);
      alloc_volume_data(*dy);
      alloc_volume_data(*dz);


      if (globals->flags.verbose>0)
	initialize_progress_report( &progress, FALSE, count[0],
				   "Interpolating new field" );
    
				/* now resample the values from the input deformation */
      for_less(i,0,count[0]) {
	for_less(j,0,count[1]) {
	  for_less(k,0,count[2]){
	    convert_3D_voxel_to_world(*dx, (Real)i, (Real)j, (Real)k,
				      &wx,&wy,&wz);
	    convert_3D_world_to_voxel(def->dx,wx,wy,wz,&nx,&ny,&nz);

	    fill_Point( voxel_point, nx, ny, nz  ); /* build the voxel POINT */
	    
	    if ( !tricubic_interpolant( def->dx, &voxel_point, &del_x ))
	      del_x = 0.0;
	    if ( !tricubic_interpolant( def->dy, &voxel_point, &del_y ))
	      del_y = 0.0;
	    if ( !tricubic_interpolant( def->dz, &voxel_point, &del_z ))
	      del_z = 0.0;
	    voxel_val = CONVERT_VALUE_TO_VOXEL(*dx, del_x); SET_VOXEL_3D(*dx, i,j,k, voxel_val);
	    voxel_val = CONVERT_VALUE_TO_VOXEL(*dy, del_y); SET_VOXEL_3D(*dy, i,j,k, voxel_val);
	    voxel_val = CONVERT_VALUE_TO_VOXEL(*dz, del_z); SET_VOXEL_3D(*dz, i,j,k, voxel_val);


	  }
	}
	if (globals->flags.verbose>0) 
	  update_progress_report( &progress, i+1 );
      }
       if (globals->flags.verbose>0) 
	terminate_progress_report( &progress );

				/* delete and free up old data */
      delete_volume(def->dx);
      delete_volume(def->dy);
      delete_volume(def->dz);
				/* set new volumes into transform` */
      def->dx = *dx;
      def->dy = *dy;
      def->dz = *dz;

      strcpy(def->basename, globals->filenames.output_trans);
      print ("globals output file name is %s\n", globals->filenames.output_trans);
      print ("def_basename file name is   %s\n", def->basename);
      
      
    }
}
