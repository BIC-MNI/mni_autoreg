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
@MODIFIED   : Revision 1.2  1994-04-06 11:47:47  louis
@MODIFIED   : working linted version of linear + non-linear registration based on Lvv
@MODIFIED   : operator working in 3D
@MODIFIED   :
 * Revision 1.1  94/02/21  16:33:41  louis
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/do_nonlinear.c,v 1.2 1994-04-06 11:47:47 louis Exp $";
#endif

#include <volume_io.h>
#include "arg_data.h"
#include "deform_field.h"
#include "line_data.h"
#include <print_error.h>

/*
 # define ITERATION_LIMIT 30
 # define WEIGHTING       0.4
 # define FRAC1           0.5
 # define FRAC2           0.0833333
 # define BRUTE_FORCE     1
 # define SECANT_METHOD   0

 hunt for !!!
*/
#define ITERATION_LIMIT 30
#define TEST_LENGTH     1.0
#define WEIGHTING       0.3
#define FRAC1           0.5
#define FRAC2           0.0833333
#define BRUTE_FORCE     1
#define SECANT_METHOD   0

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

public  Status  output_deformation_file(
    char                filename[],
    char                comments[],
    General_transform   *transform );

public Status save_deform_data(Volume dx,
			       Volume dy,
			       Volume dz,
			       char *name,
			       char *history);

public Real find_offset_to_match(Line_data *m,
				 Line_data *d,
				 Real      limit,
				 int       op_type);

private Real optimize_deformation_for_single_node(Real spacing, Real threshold,
						  Real wx, Real wy, Real wz,
						  Real *def_x, Real *def_y, Real *def_z,
						  int match_type);

private Real get_normalized_world_gradient_direction(Real xworld, Real yworld, Real zworld, 
						     Volume data_x,Volume data_y,Volume data_z,
						     Real *DX_dir, Real *DY_dir, Real *DZ_dir);

private float interp_data_along_gradient(Real dist_from, Real x, Real y, Real z, 
					 Real dx, Real dy, Real dz,
					 Volume data);

private void smooth_the_warp(Volume smooth_dx, Volume smooth_dy, Volume smooth_dz, 
			     Volume warp_dx, Volume warp_dy, Volume warp_dz);

private void build_second_derivative_data(Line_data *line, int num_samples, Real spacing,
					  Real wx, Real wy, Real wz, 
					  Real dx_dir, Real dy_dir, Real dz_dir, 
					  Volume intensity);





public void save_data(char *basename, int i, int j,
		      Volume dx,Volume dy,Volume dz);



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
  Volume            
    warps_dx,warps_dy,warps_dz,
    warps1_dx,warps1_dy,warps1_dz;
  Deform_field 
    *def_data;
  Real steps[3];
  int iters,i,j,k,sizes[3];
  Real voxel, value,value2,  def_x, def_y, def_z, wx,wy,wz,tx,ty,tz;
  Real 
    displace,
    zero,
    threshold,
    result;
  progress_struct
    progress;
  int 
    nodes_done, nodes_seen, matching_type;
  General_transform *all_until_last, *tmp_trans;

				/* set up globals for communication with other routines */
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

  /*******************************************************************************/
  /*  copy all but the last general transformation (which should be non-lin) to  */
  /*  form a general transformation 'all_until_last'                             */
  
  ALLOC(all_until_last,1)
				/* copy the first transform from the global data struct  */
  copy_general_transform(get_nth_general_transform(globals->trans_info.transformation, 0),
			 all_until_last);

				/* copy and concat the rest of tem, stopping before the end */
  for_less(i,1,get_n_concated_transforms(globals->trans_info.transformation)-1){
    ALLOC(tmp_trans,1);
    copy_general_transform(get_nth_general_transform(globals->trans_info.transformation, i),
			   tmp_trans);
    concat_general_transforms(all_until_last, 
			     tmp_trans,
			     all_until_last);
  }
  if (globals->flags.debug) {	/* print some debugging info    */
    print("orig transform is %d long\n",
	  get_n_concated_transforms(globals->trans_info.transformation));
    print("all_until_last is %d long\n",
	  get_n_concated_transforms(all_until_last));
  }

  /*******************************************************************************/
  /*  get a pointer to the last non-linear transform in the globals struct.      */
  /*  so that it can be used to update the global transformation at each         */
  /*  iterative step.                                                            */

  def_data = (Deform_field *)NULL;
  for_less(i,0,get_n_concated_transforms(globals->trans_info.transformation)) {
    if (get_transform_type( get_nth_general_transform(globals->trans_info.transformation,i) ) 
	       == USER_TRANSFORM)
      def_data = (Deform_field *)get_nth_general_transform(globals->trans_info.transformation,i)->user_data;
  }
  if (def_data == (Deform_field *)NULL) { /* exit if no deformation */
    print_error("Cannot find the deformation field transform to optimize",
		__FILE__, __LINE__);
  }

  /*******************************************************************************/
  /*   build and allocate the temporary deformation volume needed                */

  warps1_dx = copy_volume_definition(def_data->dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  warps1_dy = copy_volume_definition(def_data->dy, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  warps1_dz = copy_volume_definition(def_data->dz, NC_UNSPECIFIED, FALSE, 0.0,0.0);

  alloc_volume_data(warps1_dx);
  alloc_volume_data(warps1_dy);
  alloc_volume_data(warps1_dz);

  warps_dx = copy_volume_definition(def_data->dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  warps_dy = copy_volume_definition(def_data->dy, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  warps_dz = copy_volume_definition(def_data->dz, NC_UNSPECIFIED, FALSE, 0.0,0.0);

  alloc_volume_data(warps_dx);
  alloc_volume_data(warps_dy);
  alloc_volume_data(warps_dz);

  get_volume_sizes(warps1_dx, sizes);
  get_volume_separations(warps1_dx, steps);

  /*******************************************************************************/
  /*   initialize this iterations warp                                           */

  zero = CONVERT_VALUE_TO_VOXEL(warps_dx, 0.0);	
  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]){
	SET_VOXEL_3D(warps_dx, i,j,k, zero);
	SET_VOXEL_3D(warps_dy, i,j,k, zero);
	SET_VOXEL_3D(warps_dz, i,j,k, zero);
      }


  /*******************************************************************************/
  /*    set the threshold to be 15% of the maximum gradient magnitude            */

  threshold = 0.15 * get_volume_real_max(d1_dx); 


  /*******************************************************************************/
  /*    start the iterations to estimate the warp                     

	for each iteration {
	   for each voxel in the deformation field {
	      get the original coordinate node
	      find the best deformation for the node, taking into consideration
	        the previous deformation
	      store this best deformation
	   }

	   for each voxel in the deformation field {
 	      add WEIGHTING*best_deformation to the current deformation
	   }

	   smooth the warp;
	}
  */

  for_less(iters,0,ITERATION_LIMIT) {

    print("Iteration %2d of %2d\n",iters+1, ITERATION_LIMIT);
  
    initialize_progress_report( &progress, FALSE, sizes[0]*sizes[1] + 1,
			       "Estimating deformations" );
    
    if (iters<3)
      matching_type = BRUTE_FORCE;
    else
      matching_type = SECANT_METHOD;

    matching_type = BRUTE_FORCE; /* force matching type to be BRUTE_FORCE for now!!! */


    nodes_done = 0; nodes_seen=0; displace = 0.0;

    for_less(i,0,sizes[0]) {
      for_less(j,0,sizes[1]) {
	for_less(k,0,sizes[2]){
	  
	  nodes_seen++;
				/* remember that lattice stored in def_data was defined on 
				   the target volume! and now I have to get the corresponding
				   point back in the original space */

	  convert_3D_voxel_to_world(def_data->dx, (Real)i, (Real)j, (Real)k,
				    &wx, &wy, &wz);
	  general_inverse_transform_point(all_until_last,
					  wx,wy,wz,
					  &tx,&ty,&tz);

				/* check masking here !!!< */
	  
	  result = optimize_deformation_for_single_node(steps[0], threshold,
							tx,ty,tz,
							&def_x, &def_y, &def_z, matching_type);
	  
	  if (result != 0.0) {
	    nodes_done++;
	    displace += ABS(result);
	  }
	  else {
	    def_x = def_y = def_z = 0.0;
	  }
	  def_x = CONVERT_VALUE_TO_VOXEL(warps1_dx, def_x); 
	  SET_VOXEL_3D(warps1_dx, i,j,k, def_x);
	  def_y = CONVERT_VALUE_TO_VOXEL(warps1_dy, def_y); 
	  SET_VOXEL_3D(warps1_dy, i,j,k, def_y);
	  def_z = CONVERT_VALUE_TO_VOXEL(warps1_dz, def_z); 
	  SET_VOXEL_3D(warps1_dz, i,j,k, def_z);
	  
	}
	update_progress_report( &progress, sizes[1]*i+j+1 );
      }
    }
    terminate_progress_report( &progress );

    if (globals->flags.debug) {
      if (nodes_done>0)
	displace = displace/nodes_done;
      else
	displace=0.0;
      print ("Nodes seen = %d  , done = %d, avg displacement = %f\n",
	     nodes_seen, nodes_done, displace);
    }

				/* add this iteration to the current warp at this scale */
    for_less(i,0,sizes[0])
      for_less(j,0,sizes[1])
	for_less(k,0,sizes[2]){
	  GET_VOXEL_3D(voxel, warps1_dx, i,j,k); 
	  value = CONVERT_VOXEL_TO_VALUE( warps1_dx, voxel );
	  GET_VOXEL_3D(voxel, warps_dx, i,j,k); 
	  value2 = CONVERT_VOXEL_TO_VALUE( warps_dx, voxel );
	  value = WEIGHTING*value + value2; 
	  voxel = CONVERT_VALUE_TO_VOXEL(warps_dx, value);
	  SET_VOXEL_3D(warps_dx, i,j,k, voxel); 

	  GET_VOXEL_3D(voxel, warps1_dy, i,j,k); 
	  value = CONVERT_VOXEL_TO_VALUE( warps1_dy, voxel );
	  GET_VOXEL_3D(voxel, warps_dy, i,j,k); 
	  value2 = CONVERT_VOXEL_TO_VALUE( warps_dy, voxel );
	  value = WEIGHTING*value + value2; 
	  voxel = CONVERT_VALUE_TO_VOXEL(warps_dy, value);
	  SET_VOXEL_3D(warps_dy, i,j,k, voxel); 

	  GET_VOXEL_3D(voxel, warps1_dz, i,j,k); 
	  value = CONVERT_VOXEL_TO_VALUE( warps1_dz, voxel );
	  GET_VOXEL_3D(voxel, warps_dz, i,j,k); 
	  value2 = CONVERT_VOXEL_TO_VALUE( warps_dz, voxel );
	  value = WEIGHTING*value + value2; 
	  voxel = CONVERT_VALUE_TO_VOXEL(warps_dz, value);
	  SET_VOXEL_3D(warps_dz, i,j,k, voxel); 
	}
    

				/* smooth the current warp  
				   (smooth warps_ , returning warps1_) */

    smooth_the_warp(warps1_dx, warps1_dy, warps1_dz, 
		    warps_dx,  warps_dy,  warps_dz);
    
    
				/* copy data from temporary volume, back into
				   complete general transformation, and into
				   temp current warps... */
    for_less(i,0,sizes[0])
      for_less(j,0,sizes[1])
	for_less(k,0,sizes[2]){
	  GET_VOXEL_3D(voxel, warps1_dx, i,j,k); 
	  SET_VOXEL_3D(def_data->dx, i,j,k, voxel); 
	  SET_VOXEL_3D(warps_dx, i,j,k, voxel); 

	  GET_VOXEL_3D(voxel, warps1_dy, i,j,k); 
	  SET_VOXEL_3D(def_data->dy, i,j,k, voxel); 
	  SET_VOXEL_3D(warps_dy, i,j,k, voxel); 

	  GET_VOXEL_3D(voxel, warps1_dz, i,j,k); 
	  SET_VOXEL_3D(def_data->dz, i,j,k, voxel); 
	  SET_VOXEL_3D(warps_dz, i,j,k, voxel); 
	}

    
				/* reset the next iteration's warp. */
    zero = CONVERT_VALUE_TO_VOXEL(warps1_dx, 0.0);
    for_less(i,0,sizes[0])
      for_less(j,0,sizes[1])
	for_less(k,0,sizes[2]){
	  SET_VOXEL_3D(warps1_dx, i,j,k, zero);
	  SET_VOXEL_3D(warps1_dy, i,j,k, zero);
	  SET_VOXEL_3D(warps1_dz, i,j,k, zero);

	}
  

    if (globals->flags.debug && 
	globals->flags.verbose == 3)
      save_data(globals->filenames.output_trans, iters+1, ITERATION_LIMIT, 
		def_data->dx, def_data->dy, def_data->dz);

  }

				/* free up allocated temporary deformation volumes */

  (void)free_volume_data(warps1_dx); (void)delete_volume(warps1_dx);
  (void)free_volume_data(warps1_dy); (void)delete_volume(warps1_dy);
  (void)free_volume_data(warps1_dz); (void)delete_volume(warps1_dz);

  (void)free_volume_data(warps_dx); (void)delete_volume(warps_dx);
  (void)free_volume_data(warps_dy); (void)delete_volume(warps_dy);
  (void)free_volume_data(warps_dz); (void)delete_volume(warps_dz);

  return (OK);

}



/*******************************************************************************/
/*  procedure: smooth_the_warp

    desc: this procedure will smooth the current warp stored in warp_d? and 
          return the smoothed warp in smooth_d?

    meth: smoothing is accomplished by averaging the 6 neighbour of each node
          with the value at that node.

	  new_val = FRAC1*old_val + frac2*sum_6_neighbours(val);
    
*/
private void smooth_the_warp(Volume smooth_dx, Volume smooth_dy, Volume smooth_dz, 
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
	voxel = CONVERT_VALUE_TO_VOXEL(warp_dx, value);
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



/*******************************************************************************/
/* this routine will find the best deformation that has to be applied 
   to a single node to increase to overall objective function value.

   note:
   fwhm = 2.35*sigma
   fwtm = 4.3 *sigma

   spacing should be equal to half the fwhm value used to blur the data.

*/

private Real optimize_deformation_for_single_node(Real spacing, Real threshold,
						  Real wx, Real wy, Real wz,
						  Real *def_x, Real *def_y, Real *def_z,
						  int match_type)
{
  Line_data 
    data1, data2;
  Real
    offset,
    result,mag_normal1,mag_normal2,
    mag,
    tx,ty,tz, tempx,tempy,tempz,temp1x,temp1y,temp1z,
    px_dir, py_dir, pz_dir,
    d1x_dir, d1y_dir, d1z_dir,
    d2x_dir, d2y_dir, d2z_dir;

  result = 0.0;			/* assume no additional deformation */
  *def_x = 0.0;
  *def_y = 0.0;
  *def_z = 0.0;


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

    px_dir = d2x_dir;
    py_dir = d2y_dir;
    pz_dir = d2z_dir;


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

				/* project normal 1 in space of normal 2 */
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
      px_dir = tempx / mag;
      py_dir = tempy / mag;
      pz_dir = tempz / mag;
    } 
    else {
      px_dir = tempx;
      py_dir = tempy;
      pz_dir = tempz;
    }
      
				/* if mag_normal is too small, than use projection of d1?_dir
				   in the target space for the direction vector */
    if (mag_normal2 < threshold) {
      d2x_dir = px_dir;
      d2y_dir = py_dir;
      d2z_dir = pz_dir;
    }
    else {			/* make sure that the both normal directions are aligned,
				   and are not pointing in opposite directions! */
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
	d2x_dir = -d2x_dir;
	d2y_dir = -d2y_dir;
	d2z_dir = -d2z_dir;
      }

    }
  }

  /* check to see that the magnitude of both normals is 
     above threshold.  */


  if (mag_normal1 < threshold || mag_normal2<threshold)
    return(result);


  /* we now have direction vectors and positions for both volumes.
     Now, get samples of 2nd derivative along gradient direction 
     from both volumes. */


				/* from data1, get samples from -spacing to spacing */
  build_second_derivative_data(&data1, 5, spacing, wx,wy,wz, d1x_dir, d1y_dir, d1z_dir, Gd1);

				/* from data2, get samples from -2*spacing to 2*spacing */
  build_second_derivative_data(&data2, 9, spacing, tx,ty,tz, d2x_dir, d2y_dir, d2z_dir, Gd2);

  offset = find_offset_to_match(&data1,&data2,spacing,match_type);


/*
  *def_x = offset*d2x_dir;
  *def_y = offset*d2y_dir;
  *def_z = offset*d2z_dir;
*/
  *def_x = offset*px_dir;	/* use the projection of dir 1 */
  *def_y = offset*py_dir;
  *def_z = offset*pz_dir;

  result = offset;

  return(offset);
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
*/
private void build_second_derivative_data(Line_data *line, int num_samples, Real spacing,
					  Real wx, Real wy, Real wz, 
					  Real dx_dir, Real dy_dir, Real dz_dir, 
					  Volume intensity)
{

  int i;
  Real frac;
  Real data[30];
  
				/* set up Line_data structure */
  line->step = spacing/4.0;
  line->start = -num_samples * line->step;
  line->count = num_samples*2 +1;

				/* make space for intensity data */

				/* get intensity data */
  for_inclusive(i,-num_samples-2,num_samples+2) {
    frac = ((Real)i - 0.5)*line->step ;
    data[i+num_samples+2] = interp_data_along_gradient(frac, 
						       wx,wy,wz, 
						       dx_dir, dy_dir, dz_dir, intensity);
  }

  for_inclusive(i,-num_samples,num_samples) {

				/* fpp(xpos=0)   =  f1 - 2*f2 + f3       */
				/* fpp(xpos=0.5) = (f1 - f2   - f3 + f4)/2   */
    
    line->data[i+num_samples] = (data[i+num_samples+1] - 
				 data[i+num_samples+2] - 
				 data[i+num_samples+3] + 
				 data[i+num_samples+4])/2.0;
  }


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
						     Real *DX_dir, Real *DY_dir, Real *DZ_dir)
     
{

  Real
    dist,mag1,
    dx_dir, dy_dir, dz_dir,
    xvox, yvox, zvox;
  
  PointR 
    voxel;

  convert_3D_world_to_voxel(data_x, xworld, yworld, zworld,
			 &xvox, &yvox, &zvox);
  fill_Point( voxel, xvox, yvox, zvox );

				/* interpolate the real values from */
				/* the volume data structures */
  if (!trilinear_interpolant(data_x, &voxel, &dx_dir))  return(0.0); 
  if (!trilinear_interpolant(data_y, &voxel, &dy_dir))  return(0.0); 
  if (!trilinear_interpolant(data_z, &voxel, &dz_dir))  return(0.0); 


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

  *DX_dir = dx_dir;
  *DY_dir = dy_dir;
  *DZ_dir = dz_dir;

  return(mag1);
}

/*******************************************************************************/
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
  
  if (tricubic_interpolant(data, &voxel, &value) )
    value_flt = (float)value;
  else {
    value_flt = 0.0;
  }

  return(value_flt);
}

/*******************************************************************************/
/* debugging procedure called to save the deformation at each itereative step  */

public void save_data(char *basename, int i, int j,
		      Volume dx,Volume dy,Volume dz)
{

  Status status;
  STRING comments,name;

  (void)sprintf(comments,"step %d of %d of the non-linear estimation",i,j);
  (void)sprintf(name,"%s%d",basename,i);

  status = save_deform_data(dx,dy,dz,name,comments);
  
  if (status==OK)
    status = output_deformation_file(name,
				     comments,
				     Gglobals->trans_info.transformation);

  if (status!=OK)
    print ("Error saving %s\n",name);
}

