/* ----------------------------- MNI Header -----------------------------------
@NAME       : deform_support.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
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

@CREATED    : Tue Feb 22 08:37:49 EST 1994
@MODIFIED   : $Log: deform_support.c,v $
@MODIFIED   : Revision 96.4  2000-02-07 19:33:05  stever
@MODIFIED   : replaced HAVE_RECENT_VOLUME_IO with more specific feature tests.
@MODIFIED   :
@MODIFIED   : Revision 96.3  1999/10/25 19:59:06  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.2  1997/11/03  20:05:41  louis
 * reorganized deform_support...
 * putting:
 *   init_the_volume_to_zero()
 *   get_volume_maximum_real_value()
 *   save_data()
 *  into extras.c
 *
 * and putting:
 *   general_transform_point_in_trans_plane()
 *   build_source_lattice()
 *   go_get_samples_in_source()
 *   go_get_samples_with_offset()
 *   build_target_lattice()
 *   build_target_lattice_using_super_sampled_def()
 *  into sub_lattice.c
 *
 * Revision 96.2  1997/11/03  20:05:41  louis
 * reorganized deform_support...
 * putting:
 *   init_the_volume_to_zero()
 *   get_volume_maximum_real_value()
 *   save_data()
 *  into extras.c
 *
 * and putting:
 *   general_transform_point_in_trans_plane()
 *   build_source_lattice()
 *   go_get_samples_in_source()
 *   go_get_samples_with_offset()
 *   build_target_lattice()
 *   build_target_lattice_using_super_sampled_def()
 *  into sub_lattice.c
 *
 * Revision 96.1  1997/11/03  15:06:29  louis
 * working version, before creation of mni_animal package, and before inserting
 * distance transforms
 *
 * Revision 96.0  1996/08/21  18:22:10  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:01  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:03  louis
 * Never released version 0.95
 *
 * Revision 1.16  1996/08/12  14:15:53  louis
 * Pre-release
 *
 * Revision 1.14  1996/04/01  09:16:29  collins
 * removed the code to super sample the deformation field,
 * and placed it in super_sample_def.c.
 *
 * Revision 1.13  1996/04/01  09:02:10  collins
 * added optimized code to super-sample the deformation field for the special
 * case when super_sample = 2.  interpolate_super_sampled_data_by2() is
 * approximately 7 times faster than interpolate_super_sampled_data().
 *
 * Revision 1.12  1996/03/25  10:33:15  collins
 * used inter_type to specify degress_continuity in the call to
 * evaluate_volume_in_world() in procedure go_get_samples_in_source().
 *
 * ----
 * added support for tri-linear interpolation in
 * go_get_samples_with_offset() since the local objective function is not
 * smooth (but is actually step-wise) with NN interpolation.  This is
 * specified by an extra boolean parameter to go_get_samples_with_offset
 * that should be TRUE to use NN interpolation.  This is controled from
 * the command line by the -nearest_neighbour / -trilinear options.
 *
 * NOTE: NN is faster, but the local obj function will have steps (ie it
 *       will have plateaus and discontinuous changes).  Also, there
 *       seems to be a `drapery' effect, such that the function is _not_
 *       monotonically increasing with mis-registration -> thus yielding
 *       local minima.
 *       -> trilin is the default, and should be used unless you really
 *          know what you are doing!
 *
 * ----
 *
 * removed dx += 0.5; dy += 0.5; dz += 0.5.; in
 * go_get_samples_with_offset().  this was originally supposed to give a
 * `rounding' effect for NN interpolation, but is no longer needed with
 * tri-linear interpolation.  IN FACT: keeping it caused a certain bias
 * so that a deformation was returned _even_ when an object was matched
 * to itself; after removal- there are only small sporadic deviations
 * from a null deformation field.
 *
 * Revision 1.11  1996/03/07  13:25:19  collins
 * small reorganisation of procedures and working version of non-isotropic
 * smoothing.
 *
 * Revision 1.10  1995/10/06  09:25:02  collins
 * removed references to line_data.h since it hos not been used in a while.
 *
 * included "constants.h" to have access to NONLIN_* similarity func ids.
 *
 * modified go_get_samples_with_offset to account for different similarity
 * functions.
 *
 * Revision 1.9  1995/09/07  10:05:11  collins
 * All references to numerical recipes routines are being removed.  At this
 * stage, any num rec routine should be local in the file.  All memory
 * allocation calls to vector(), matrix(), free_vector() etc... have been
 * replaced with ALLOC and FREE from the volume_io library of routines.
 *
 * Revision 1.8  1995/06/12  14:29:46  collins
 * working version - 2d,3d w/ simplex and -direct.
 *
 * Revision 1.7  1995/05/04  14:25:18  collins
 * compilable version, seems to run a bit with GRID_TRANSFORM, still
 * needs work for the super sampled volumes... and lots of testing.
 *
 * Revision 1.6  1995/05/02  11:30:27  collins
 * started clean up of code, separation of non used procedures into
 * old_methods.c.  This version was working, but I am now going to
 * rewrite everything to use GRID_TRANSFORM.
 *
 * Revision 1.5  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.4  94/06/21  10:59:28  louis
 * working optimized version of program.  when compiled with -O, this
 * code is approximately 60% faster than the previous version.
 * 
 * 
 * Revision 1.3  94/06/06  18:46:53  louis
 * working version: clamp and blur of deformation lattice now ensures
 * a smooth recovered deformation.  Unfortunately, the r = cost-similarity
 * function used in the optimization is too heavy on the cost_fn.  This has
 * to get fixed...
 * 
 * 

 * Revision 1.2  94/06/02  20:15:56  louis
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
 * Revision 1.1  94/04/06  11:47:27  louis
 * Initial revision
 * 

---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/deform_support.c,v 96.4 2000-02-07 19:33:05 stever Exp $";
#endif

#include <config.h>
#include <internal_volume_io.h>
#include "arg_data.h"
#include <louis_splines.h>
#include <print_error.h>
#include "local_macros.h"
#include "constants.h"

extern Arg_Data main_args;

#define DERIV_FRAC      0.6
#define FRAC1           0.5
#define FRAC2           0.0833333
#define ABSOLUTE_MAX_DEFORMATION       50.0

extern double smoothing_weight;
extern char *my_XYZ_dim_names;

public void get_volume_XYZV_indices(Volume data, int xyzv[]);

public int trilinear_interpolant(Volume volume, 
                                 PointR *coord, double *result);
 
public int nearest_neighbour_interpolant(Volume volume, 
                                         PointR *coord, double *result);

public void init_the_volume_to_zero(Volume volume);

public Real get_volume_maximum_real_value(Volume volume);

 

				/* define the limits for three nested
				   for loops, so that we loop through
				   each spatial dimension.
				   'start' and 'end' are returned in
				   XYZ order.  */

public void  get_voxel_spatial_loop_limits(Volume volume,
					   int start[],	
					   int end[])
{
  int 
    i,
    sizes[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS];

  get_volume_sizes(volume, sizes);
  get_volume_XYZV_indices(volume, xyzv);

  for_less(i,0,N_DIMENSIONS) {
    if (sizes[ xyzv[i] ]>3) {
      start[ i ] = 1;
      end[ i ] = sizes[ xyzv[i] ]-1;
    }
    else {
      start[ i ] = 0;
      end[ i ] = sizes[ xyzv[i] ];
    }
  }

  if (get_volume_n_dimensions(volume)>3) {
    start[Z+1] = 0;
    end[Z+1] = sizes[ xyzv[Z+1] ] ;
  }
  else {
    start[Z+1] = 0;
    end[Z+1] = 0;
  }
    
}

public BOOLEAN get_average_warp_vector_from_neighbours(General_transform *trans,
						       int voxel[],
						       int avg_type,
						       Real *mx, Real *my, Real *mz)
{
  int
    start[MAX_DIMENSIONS],
    end[MAX_DIMENSIONS],
    count, i, 
    voxel2[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    sizes[MAX_DIMENSIONS];
  Real 
    def_vector[N_DIMENSIONS];
  Volume volume;
  
  if (trans->type != GRID_TRANSFORM) {
    print_error_and_line_num("get_average_warp_vector_from_neighbours not called with GRID_TRANSFORM",
			     __FILE__, __LINE__);
    return (FALSE);
  }

  volume = trans->displacement_volume;
  
  get_volume_sizes(volume, sizes);
  get_volume_XYZV_indices(volume, xyzv);
  count = 0;
  *mx = 0.0; *my = 0.0; *mz = 0.0; /* assume no warp vector in volume */

				/* make sure voxel is in volume */
  for_less(i,0,3 ) {
    if (voxel[ xyzv[i]]<0 || voxel[ xyzv[i] ]>=sizes[ xyzv[i]] ) {
       return (FALSE);
    }
  }
  
  for_less(i,0,MAX_DIMENSIONS ) { /* copy the voxel position */
    voxel2[i] = voxel[i];
  }
  
  switch (avg_type) {
  case 1:  {			/* get the 6 neighbours, 2 along
				   each of the spatial axes,
				   if they exist 
				   ie the 4-connected immediate neighbours only */

    for_less(i, 0 , N_DIMENSIONS) {
      
      if ((voxel[ xyzv[i] ]+1) < sizes[ xyzv[i] ]) {
	
	voxel2[ xyzv[i] ] = voxel[ xyzv[i] ] + 1;
	
	for_less(voxel2[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
	  def_vector[voxel2[ xyzv[Z+1] ]] = 
	    get_volume_real_value(volume,
				  voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
	}
	
	voxel2[ xyzv[i] ] = voxel[ xyzv[i] ];
	
	*mx += def_vector[X]; *my += def_vector[Y]; *mz += def_vector[Z];
	++count;
      }
      if ((voxel[ xyzv[i] ]-1) >= 0) {
	voxel2[ xyzv[i] ] = voxel[ xyzv[i] ] - 1;
	
	for_less(voxel2[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
	  def_vector[voxel2[ xyzv[Z+1] ]] = 
	    get_volume_real_value(volume,
				  voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
	}
	
	voxel2[ xyzv[i] ] = voxel[ xyzv[i] ];
	
	*mx += def_vector[X]; *my += def_vector[Y]; *mz += def_vector[Z];
	++count;
      }
      
    }
    break;
  }
  case 2: {			/* 3x3x3 ie the 26 8-connected
				   immediate neighbours only */
    
    for_less(i, 0 , N_DIMENSIONS) {
      start[i] = voxel[ xyzv[i] ] - 1;
      if (start[i]<0) start[i]=0;
      end[i] = voxel[ xyzv[i] ] + 1;
      if (end[i]>sizes[ xyzv[i] ]-1) end[i] = sizes[ xyzv[i] ]-1;
    }
    
    
    for_inclusive(voxel2[ xyzv[X] ],start[X], end[X])
      for_inclusive(voxel2[ xyzv[Y] ],start[Y], end[Y])
	for_inclusive(voxel2[ xyzv[Z] ],start[Z], end[Z]) {

	  if ((voxel2[ xyzv[X]] != voxel[ xyzv[X] ]) ||
	      (voxel2[ xyzv[Y]] != voxel[ xyzv[Y] ]) ||
	      (voxel2[ xyzv[Z]] != voxel[ xyzv[Z] ])) {
	    for_less(voxel2[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
	      def_vector[voxel2[ xyzv[Z+1] ]] = 
		get_volume_real_value(volume,
				      voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
	    }
	    *mx += def_vector[X]; *my += def_vector[Y]; *mz += def_vector[Z];
	    ++count;
	  }
	}
    break;
  }
  case 3: {			/* 5x5x5 */
    for_less(i, 0 , N_DIMENSIONS) {
      start[i] = voxel[ xyzv[i] ] - 2;
      if (start[i]<0) start[i]=0;
      end[i] = voxel[ xyzv[i] ] + 2;
      if (end[i]>sizes[ xyzv[i] ]-1) end[i] = sizes[ xyzv[i] ]-1;
    }
    
    
    for_inclusive(voxel2[ xyzv[X] ],start[X], end[X])
      for_inclusive(voxel2[ xyzv[Y] ],start[Y], end[Y])
	for_inclusive(voxel2[ xyzv[Z] ],start[Z], end[Z]) {
	  
	  if ((voxel2[ xyzv[X]] != voxel[ xyzv[X] ]) ||
	      (voxel2[ xyzv[Y]] != voxel[ xyzv[Y] ]) ||
	      (voxel2[ xyzv[Z]] != voxel[ xyzv[Z] ])) {

	    for_less(voxel2[ xyzv[Z+1] ], 0, sizes[ xyzv[Z+1] ]) {
	      def_vector[voxel2[ xyzv[Z+1] ]] = 
		get_volume_real_value(volume,
				      voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
	    }
	    *mx += def_vector[X]; *my += def_vector[Y]; *mz += def_vector[Z];
	    ++count;

	  }
	}
    break;
    
  }
  }    

  if (count>0) {		/* average deformation vector */
    *mx /= count; 
    *my /= count; 
    *mz /= count; 
    return(TRUE);
  }
  else {
    return(FALSE);
  }
}


public BOOLEAN get_average_warp_of_neighbours(General_transform *trans,
					      int voxel[],
					      Real mean_pos[])
{
  int       i;
  Real      voxel_real[MAX_DIMENSIONS],
            dx, dy, dz;
  Volume     volume;
  
  if (trans->type != GRID_TRANSFORM) {
    print_error_and_line_num("get_average_warp_of_neighbours not called with GRID_TRANSFORM",
			     __FILE__, __LINE__);
    return (FALSE);
  }

  volume = trans->displacement_volume;
  
  for_less(i,0,get_volume_n_dimensions(volume) ) {
    voxel_real[i] = (Real)voxel[i];
  }
  convert_voxel_to_world(volume, voxel_real, 
			 &(mean_pos[X]),&(mean_pos[Y]),&(mean_pos[Z]) );

  if ( ! get_average_warp_vector_from_neighbours(trans, voxel, 1, &dx, &dy, &dz)) {
    return(FALSE);
  }
  else {
    mean_pos[X] += dx; mean_pos[Y] += dy; mean_pos[Z] += dz;
    return(TRUE);
  }

}

/* add additional to current, return
   answer in additional 

   I assume for this procedure that the displacement_volume
   for both the additional and current transform have the 
   same volumetric definition ie, both have the same dimension
   order and same length along each dimension.
*/

public void add_additional_warp_to_current(General_transform *additional,
					   General_transform *current,
					   Real weight)
{
  int
    count[MAX_DIMENSIONS],
    count_current[MAX_DIMENSIONS],
    xyzv_additional[MAX_DIMENSIONS],
    xyzv_current[MAX_DIMENSIONS],
    index[MAX_DIMENSIONS],
    i;
  Real 
    additional_value, current_value;


  if (get_volume_n_dimensions(additional->displacement_volume) != 
      get_volume_n_dimensions(current->displacement_volume)) {
    print_error_and_line_num("add_additional_warp_to_current: warp dim error",
			     __FILE__, __LINE__);
  }

  get_volume_sizes(additional->displacement_volume, count);
  get_volume_sizes(current->displacement_volume, count_current);
  for_less(i,0,get_volume_n_dimensions(current->displacement_volume)) {
    if (count_current[i] != count[i]) {
      print_error_and_line_num("add_additional_warp_to_current: dim count error",
			       __FILE__, __LINE__);
    }
  }

  get_volume_XYZV_indices(additional->displacement_volume, xyzv_additional);
  get_volume_XYZV_indices(current->displacement_volume, xyzv_current);
  for_less(i,0,get_volume_n_dimensions(current->displacement_volume)) {
    if (xyzv_current[i] != xyzv_additional[i]) {
      print_error_and_line_num("add_additional_warp_to_current: dim match error",
			       __FILE__, __LINE__);
    }
  }

  for_less(i,0,MAX_DIMENSIONS) index[i]=0;

  for_less(index[ xyzv_additional[X] ], 0, count[ xyzv_additional[X] ])
    for_less(index[ xyzv_additional[Y] ], 0, count[ xyzv_additional[Y] ])
      for_less(index[ xyzv_additional[Z] ], 0, count[ xyzv_additional[Z] ])
	for_less(index[ xyzv_additional[Z+1] ], 0, count[ xyzv_additional[Z+1] ]) {

	  additional_value = get_volume_real_value(
				 additional->displacement_volume,
				 index[0],index[1],index[2],index[3],index[4]);
	  current_value = get_volume_real_value(
                                 current->displacement_volume,
				 index[0],index[1],index[2],index[3],index[4]);

	  additional_value = current_value + additional_value*weight;

	  if (additional_value >  40.0) additional_value =  40.0;
	  if (additional_value < -40.0) additional_value = -40.0;

	  set_volume_real_value(additional->displacement_volume,
				index[0],index[1],index[2],index[3],index[4],
				additional_value);

	  
	}

}
					    


/*******************************************************************
  procedure: smooth_the_warp

    desc: this procedure will smooth the current warp stored in
          current and return the smoothed warp in smoothed

    meth: smoothing is accomplished by averaging the value of the 
          node's deformation vector with the mean deformation vector
	  of it's neighbours.

	  def'  = sw*mean  + (1-sw)*def

	  where: sw   = smoothing_weight
                 mean = neighbourhood mean deformation
		 def  = estimate def for current node
*/

public void smooth_the_warp(General_transform *smoothed,
			    General_transform *current,
			    Volume warp_mag, Real thres) 
{
  int
    count_smoothed[MAX_DIMENSIONS],
    count_current[MAX_DIMENSIONS],
    count_mag[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    xyzv_current[MAX_DIMENSIONS],
    xyzv_mag[MAX_DIMENSIONS],
    mag_index[MAX_DIMENSIONS],
    index[MAX_DIMENSIONS],
    start[MAX_DIMENSIONS], 
    end[MAX_DIMENSIONS],
    i;
  Real 
    value[3], 
    wx, wy, wz, 
    mx, my, mz;
  progress_struct
    progress;
  
  
  if (get_volume_n_dimensions(smoothed->displacement_volume) != 
      get_volume_n_dimensions(current->displacement_volume)) {
    print_error_and_line_num("smooth_the_warp: warp dim error",
			     __FILE__, __LINE__);
  }
  
  get_volume_sizes(smoothed->displacement_volume, count_smoothed);
  get_volume_sizes(current->displacement_volume, count_current);
  for_less(i,0,get_volume_n_dimensions(current->displacement_volume)) {
    if (count_current[i] != count_smoothed[i]) {
      print_error_and_line_num("smooth_the_warp: dim count error",
			       __FILE__, __LINE__);
    }
  }

  get_volume_XYZV_indices(smoothed->displacement_volume, xyzv);
  get_volume_XYZV_indices(current->displacement_volume, xyzv_current);
  for_less(i,0,get_volume_n_dimensions(current->displacement_volume)) {
    if (xyzv_current[i] != xyzv[i]) {
      print_error_and_line_num("smooth_the_warp: dim match error",
			       __FILE__, __LINE__);
    }
  }
  
  get_volume_XYZV_indices(warp_mag, xyzv_mag);
  get_volume_sizes(warp_mag, count_mag);

  for_less(i,0,get_volume_n_dimensions(warp_mag)) {
    if (count_current[xyzv_current[i]] != count_mag[i]) {
      print_error_and_line_num("smooth_the_warp: dim count error w/mag (%d: %d != %d)\n",
			       __FILE__, __LINE__, i, count_current[xyzv_current[i]], count_mag[i] );
    }
  }
  
  for_less(i,0,MAX_DIMENSIONS) {
    index[i]=0;
    start[i] = 0;
    end[i] = 0;
  }
  
  get_voxel_spatial_loop_limits(smoothed->displacement_volume, start, end);
  start[Z+1] = 0;
  end[Z+1] = 3;
 

  initialize_progress_report( &progress, FALSE, 
			     (end[X]-start[X])*
			     (end[Y]-start[Y]) + 1,
			     "Smoothing deformations" );


  for_less(index[ xyzv[X] ], start[ X ], end[ X ]) {
    for_less(index[ xyzv[Y] ], start[ Y ], end[ Y ]) {
      for_less(index[ xyzv[Z] ], start[ Z ], end[ Z ]) {

	for_less(i,0,get_volume_n_dimensions(warp_mag))
	  mag_index[ xyzv_mag[i] ] = index[ xyzv[i] ];
	
	if (  get_volume_real_value(warp_mag, 
				    mag_index[ X ],
				    mag_index[ Y ],
				    mag_index[ Z ],
				    0, 0) >= thres) {

				/* go get the current warp vector for
                                   this node. */
	  
	  for_less(index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) {

	    value[index[ xyzv[Z+1] ]] = 
	      get_volume_real_value(current->displacement_volume,
				    index[0],index[1],index[2],
				    index[3],index[4]);

	  }
				/* store the current warp in wx, wy,
                                   wz */

	  wx = value[X]; wy = value[Y]; wz = value[Z]; 

	  index[ xyzv[Z+1]] = 0;

				/* if we can get a neighbourhood mean
				   warp vector, then we average it
				   with the current warp vector */

	  if ( get_average_warp_vector_from_neighbours(current,
						      index, 2 ,
						      &mx, &my, &mz) ) {
	    
	    wx = (1.0 - smoothing_weight) * value[X] + smoothing_weight * mx;
	    wy = (1.0 - smoothing_weight) * value[Y] + smoothing_weight * my;
	    wz = (1.0 - smoothing_weight) * value[Z] + smoothing_weight * mz;
	    value[X] = wx; 
	    value[Y] = wy; 
	    value[Z] = wz; 

	  } 
	  
				/* now put the averaged vector into
				   the smoothed volume */
 
	  for_less(index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ])  
	    set_volume_real_value(smoothed->displacement_volume,
				  index[0],index[1],index[2],
				  index[3],index[4],
				  value[index[ xyzv[Z+1] ]] );  
	}
	else {

	  for_less(index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ])  {
	    
	    value[index[ xyzv[Z+1] ]] = 
	      get_volume_real_value(current->displacement_volume,
				    index[0],index[1],index[2],
				    index[3],index[4]);

	    set_volume_real_value(smoothed->displacement_volume,
				  index[0],index[1],index[2],
				  index[3],index[4],
				  value[index[ xyzv[Z+1] ]] );
	  }

	} /* mag > thresh */
	  
      }
      update_progress_report( &progress,
			     ((end[ Y ]-start[ Y ])*
			      (index[ xyzv[X]  ]-start[ X ])) +
			     (index[ xyzv[Y] ]-start[ X ])  +    1  );
    }
  }

  terminate_progress_report( &progress );
}


/*

   We want to extrapolate (and smooth) the estimated deformations to
   nodes where no estimation was possible (and thus no local smoothing
   completed).  This procedure should only be called when
   Gglobals->trans_info.use_local_smoothing is true.  (When false,
   extrapolation is not needed, since it is addressed in the global
   smoothing process.)

   additional  = sw*mean    + (1-sw)current
             i          i-1                i-1 
               
   Here we must store the mean_vector location in additional_warp for
   all nodes where there was no estimatation possible (and thus no
   local smoothing completed).  This will achieve a homogeneous
   smoothing throughout the entire field.  (Here, smoothing_weight is
   effectively equal to 1.0 - since the local node has no information
   available to it.

   The idea is to set the node value in additional_vol, so that when
   added to the corresponding node value in current_vol, a smoothed
   vector field results.
   ie 
      additional = required_absolute_mean_vector - absolute_current_vector


   note estimated_flag_vol is created to be accessed in [X][Y][Z] order.

      */

public void extrapolate_to_unestimated_nodes(General_transform *current,
					     General_transform *additional,
					     Volume estimated_flag_vol) 
{

  int 
    many,
    total,
    extrapolated,
    count,
    count_additional[MAX_DIMENSIONS],
    count_current[MAX_DIMENSIONS],
    count_flag[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    xyzv_current[MAX_DIMENSIONS],
    xyzv_flag[MAX_DIMENSIONS],
    flag_index[MAX_DIMENSIONS],
    index[MAX_DIMENSIONS],
    start[MAX_DIMENSIONS], 
    end[MAX_DIMENSIONS],
    start2[MAX_DIMENSIONS], 
    end2[MAX_DIMENSIONS],
    voxel2[MAX_DIMENSIONS],
    i;
  Real 
    current_deform[N_DIMENSIONS], 
    additional_deform[N_DIMENSIONS], 
    mx, my, mz;
  progress_struct
    progress;

  extrapolated = many = total = 0;

                                /* verify that the deformation volumes for
                                   current and additional warps are
                                   compatible in dimension, size and
                                   order. */

  if (get_volume_n_dimensions(additional->displacement_volume) != 
      get_volume_n_dimensions(current->displacement_volume)) {
    print_error_and_line_num("extrapolate_the_warp: warp dim error",
			     __FILE__, __LINE__);
  }
  
  get_volume_sizes(additional->displacement_volume, count_additional);
  get_volume_sizes(current->displacement_volume, count_current);
  for_less(i,0,get_volume_n_dimensions(current->displacement_volume)) {
    if (count_current[i] != count_additional[i]) {
      print_error_and_line_num("extrapolate_the_warp: dim count error",
			       __FILE__, __LINE__);
    }
  }

  get_volume_XYZV_indices(additional->displacement_volume, xyzv);
  get_volume_XYZV_indices(current->displacement_volume, xyzv_current);
  for_less(i,0,get_volume_n_dimensions(current->displacement_volume)) {
    if (xyzv_current[i] != xyzv[i]) {
      print_error_and_line_num("extrapolate_the_warp: dim match error",
			       __FILE__, __LINE__);
    }
  }
  
  get_volume_XYZV_indices(estimated_flag_vol, xyzv_flag);
  get_volume_sizes(estimated_flag_vol, count_flag);

  for_less(i,0,get_volume_n_dimensions(estimated_flag_vol)) {
    if (count_current[xyzv_current[i]] != count_flag[i]) {
      print_error_and_line_num("extrapolate_the_warp: dim count error w/flag (%d: %d != %d)\n",
			       __FILE__, __LINE__, i, count_current[xyzv_current[i]], count_flag[i] );
    }
  }

                                /* verification completed, now get loop
                                   limits to go through the deformation
                                   volume and extrapolate the estimated
                                   vectors from the additional volume */
  for_less(i,0,MAX_DIMENSIONS) {
    index[i] = 0;
    start[i] = 0;
    end[i]   = 0;
  }
  
  get_voxel_spatial_loop_limits(additional->displacement_volume, start, end);
  start[Z+1] = 0;
  end[Z+1]   = 3;
 
  initialize_progress_report( &progress, FALSE, 
			     (end[X]-start[X])*
			     (end[Y]-start[Y]) + 1,
			     "Extrapolating estimations" );


  for_less(index[ xyzv[X] ], start[ X ], end[ X ]) {
    for_less(index[ xyzv[Y] ], start[ Y ], end[ Y ]) {
      for_less(index[ xyzv[Z] ], start[ Z ], end[ Z ]) {

	for_less(i,0,get_volume_n_dimensions(estimated_flag_vol))
	  flag_index[ xyzv_flag[i] ] = index[ xyzv[i] ];
	
	total++;


	if (  get_volume_real_value(estimated_flag_vol, 
				    index[ xyzv[X] ],
				    index[ xyzv[Y] ],
				    index[ xyzv[Z] ],
				    0, 0) < 1.0) {

	  /* 
             then, this node has not been estimated at this iteration,
             so we have to interpolate a deformation value from
             neighbouring nodes, and then apply the proper local
	     smoothing. 

             local smoothing is defined to be the average deformation of
             the 26 8-connected neighbours and is computed in the loop
             below
          */

	  many++;
				/* go get the current warp vector for
                                   this node. */
	  
	  for_less(index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) {

	    current_deform[index[ xyzv[Z+1] ]] = 
	      get_volume_real_value(current->displacement_volume,
				    index[0],index[1],index[2],
				    index[3],index[4]);

	  }
				/* get an average of the current additional 
				   deformation */

	  for_less(i, 0 , N_DIMENSIONS) {
	    start2[i] = index[ xyzv[i] ] - 1;
	    if (start2[i]<0) start2[i]=0;
	    end2[i] = index[ xyzv[i] ] + 1;
	    if (end2[i]> count_flag[ xyzv_flag[i] ]-1) end2[i] =  count_flag[ xyzv_flag[i] ]-1;
	  }
    
	  for_less(i,0,N_DIMENSIONS)   additional_deform[i]= 0.0;
	  for_less(i,0,MAX_DIMENSIONS) voxel2[i]=0;
	  count = 0;

	  for_inclusive(voxel2[ xyzv[X] ],start2[X], end2[X])
	    for_inclusive(voxel2[ xyzv[Y] ],start2[Y], end2[Y])
	      for_inclusive(voxel2[ xyzv[Z] ],start2[Z], end2[Z]) {

		if (get_volume_real_value(estimated_flag_vol, 
				    voxel2[ xyzv[X] ],
				    voxel2[ xyzv[Y] ],
				    voxel2[ xyzv[Z] ],
				    0, 0) >= 0.5) {

		  if (!((voxel2[ xyzv[X]] == index[ xyzv[X] ]) &&
			(voxel2[ xyzv[Y]] == index[ xyzv[Y] ]) &&
			(voxel2[ xyzv[Z]] == index[ xyzv[Z] ]))   ) {
		    
		    for_less(voxel2[ xyzv[Z+1] ], 0, N_DIMENSIONS) {
		      additional_deform[ voxel2[ xyzv[Z+1] ] ] += 
			get_volume_real_value(additional->displacement_volume,
					      voxel2[0],voxel2[1],voxel2[2],
					      voxel2[3],voxel2[4]);
		    }
		    ++count;
		  }
		}

		
	      }

	  if (count>0) {
	    extrapolated++;
	    for_less(i,0,N_DIMENSIONS)
	      additional_deform[i] /= 26.0;
	  }

				/* if we can get a neighbourhood mean
				   warp vector from the previous iterations, 
				   then we average it with the previous warp 
				   vector */
	  index[ xyzv[Z+1]] = 0;
	  if ( get_average_warp_vector_from_neighbours(current,
						       index, 2 ,
						       &mx, &my, &mz) ) {

	    /* additional_deform += sw*mean + (1-sw)*current - current

	       with sw = 0.5 gives: */
	    
	    additional_deform[X] += (mx - current_deform[X])/2.0; 
	    additional_deform[Y] += (my - current_deform[Y])/2.0; 
	    additional_deform[Z] += (mz - current_deform[Z])/2.0; 
	  } 

	  
				/* now put the averaged vector into
				   the additional volume */

	  for_less(index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ])  
	    set_volume_real_value(additional->displacement_volume,
				  index[0],index[1],index[2],
				  index[3],index[4],
				  additional_deform[index[ xyzv[Z+1] ]] );  
	}
	  
      }
      update_progress_report( &progress,
			     ((end[ Y ]-start[ Y ])*
			      (index[ xyzv[X]  ]-start[ X ])) +
			     (index[ xyzv[Y] ]-start[ X ])  +    1  );
    }
  }

  terminate_progress_report( &progress );

  print ("There were %d out of %d extrapolated (%d left) (%d extrapolated)\n",many,total,total-many, extrapolated);

}


/* copied from tricubic_interpolant in interpolation.c */

private int tricubic_slice_interpolant(Volume volume, 
				      PointR *coord, double *result)
{
  Real
    v00,v01,v02,v03, 
    v10,v11,v12,v13, 
    v20,v21,v22,v23, 
    v30,v31,v32,v33;
   long 
     ind0, ind1, ind2, max[3];
   double 
     frac[3];
   int 
     sizes[3];
   int 
     flag;
   double temp_result;

   /* Check that the coordinate is inside the volume */

   get_volume_sizes(volume, sizes);
   max[0] = sizes[0];
   max[1] = sizes[1];
   max[2] = sizes[2];
   
   if ((Point_y( *coord ) < 0) || (Point_y( *coord ) >= max[1]-1) ||
       (Point_z( *coord ) < 0) || (Point_z( *coord ) >= max[2]-1)) {

     flag = nearest_neighbour_interpolant(volume, coord, &temp_result) ;
     *result = temp_result;
     return(flag);
   }


   /* Get the whole and fractional part of the coordinate */
   ind0 = (long) Point_x( *coord );
   ind1 = (long) Point_y( *coord );
   ind2 = (long) Point_z( *coord );
   frac[1] = Point_y( *coord ) - ind1;
   frac[2] = Point_z( *coord ) - ind2;
   ind1--;
   ind2--;

   /* Check for edges - do linear interpolation at edges */
   if ((ind1 >= max[1]-3) || (ind1 < 0) ||
       (ind2 >= max[2]-3) || (ind2 < 0)) {
      return trilinear_interpolant(volume, coord, result);
   }
   
   GET_VOXEL_3D(v00, volume, ind0, ind1, ind2);
   GET_VOXEL_3D(v01, volume, ind0, ind1, ind2+1);
   GET_VOXEL_3D(v02, volume, ind0, ind1, ind2+2);
   GET_VOXEL_3D(v03, volume, ind0, ind1, ind2+3);

   GET_VOXEL_3D(v10, volume, ind0, ind1+1, ind2);
   GET_VOXEL_3D(v11, volume, ind0, ind1+1, ind2+1);
   GET_VOXEL_3D(v12, volume, ind0, ind1+1, ind2+2);
   GET_VOXEL_3D(v13, volume, ind0, ind1+1, ind2+3);

   GET_VOXEL_3D(v20, volume, ind0, ind1+2, ind2);
   GET_VOXEL_3D(v21, volume, ind0, ind1+2, ind2+1);
   GET_VOXEL_3D(v22, volume, ind0, ind1+2, ind2+2);
   GET_VOXEL_3D(v23, volume, ind0, ind1+2, ind2+3);

   GET_VOXEL_3D(v30, volume, ind0, ind1+3, ind2);
   GET_VOXEL_3D(v31, volume, ind0, ind1+3, ind2+1);
   GET_VOXEL_3D(v32, volume, ind0, ind1+3, ind2+2);
   GET_VOXEL_3D(v33, volume, ind0, ind1+3, ind2+3);

   CUBIC_BIVAR( v00,v01,v02,v03, v10,v11,v12,v13, v20,v21,v22,v23, v30,v31,v32,v33, frac[1],frac[2], temp_result);


  *result = CONVERT_VOXEL_TO_VALUE(volume,temp_result);

   return(TRUE);

}




public Real get_value_of_point_in_volume(Real xw, Real yw, Real zw, 
					  Volume data)
     
{

  Real
    mag,
    xvox, yvox, zvox;
  
  PointR 
    voxel;

  convert_3D_world_to_voxel(data, xw, yw, zw,
			    &zvox, &yvox, &xvox);
  fill_Point( voxel, zvox, yvox, xvox );

  if (!trilinear_interpolant(data, &voxel, &mag)) 
    return(-DBL_MAX); 
  else 
    return(mag);
}



/*********************************************************************** 
   split the General_transform stored in 'total' to extract two parts,
      1 - the last warp (that will be optimized)
      2 - everything else,
   so that:

   total = all_until_last + last_warp 
*/

public void split_up_the_transformation(General_transform *total,
					General_transform **all_until_last,
					General_transform **last_warp) 
{
  General_transform
    *tmp_trans;
  int i;

  ALLOC(*all_until_last,1);
				/* copy the first transform from the global data struct  */
  copy_general_transform(get_nth_general_transform(total, 0),
			 *all_until_last);

				/* copy and concat the rest of tem, stopping before the end */
  for_less(i,1,get_n_concated_transforms(total)-1){
    ALLOC(tmp_trans,1);
    copy_general_transform(get_nth_general_transform(total, i), tmp_trans);
    concat_general_transforms(*all_until_last, tmp_trans, *all_until_last);
  }

  *last_warp = (General_transform *)NULL;
  for_less(i,0,get_n_concated_transforms(total)) {
    if (get_transform_type( get_nth_general_transform(total,i) ) 
	       == GRID_TRANSFORM)
      *last_warp = get_nth_general_transform(total,i);
  }

}



/*   return the maximum value stored in the data volume */

private Real get_maximum_magnitude(Volume dxyz)
{

  Real max;

  max = get_volume_maximum_real_value(dxyz);

  return(max);
  
}

/***************************************************************************/
/*    set the threshold to be 10% of the maximum gradient magnitude        */
/*    for each source and target volumes                                   */

public void  set_feature_value_threshold(Volume d1, 
					 Volume d2,
					 Real *global_thres1, 
					 Real *global_thres2, 
					 Real *threshold1, 
					 Real *threshold2) 

{
  if (*global_thres1==0.0)
    *threshold1 = 0.10 * get_maximum_magnitude(d1);
  else
    *threshold1 = *global_thres1 * get_maximum_magnitude(d1);

  if (*global_thres2==0.0)
    *threshold2 = 0.10 * get_maximum_magnitude(d2);
  else
    *threshold2 = *global_thres2 * get_maximum_magnitude(d2);
}



public void build_two_perpendicular_vectors(Real orig[], 
					     Real p1[], 
					     Real p2[])
{
  Vector
    v,v1,v2;
  Real
    len;
  fill_Vector(v, orig[X], orig[Y], orig[Z]);
  create_two_orthogonal_vectors( &v, &v1, &v2);

  len = MAGNITUDE( v1 );
  if (len>0) {
    p1[X] = Vector_x(v1) / len;
    p1[Y] = Vector_y(v1) / len;
    p1[Z] = Vector_z(v1) / len;
  }
  else
    print_error_and_line_num("Null length for vector normalization\n", 
		__FILE__, __LINE__);

  len = MAGNITUDE( v2 );
  if (len>0) {
    p2[X] = Vector_x(v2) / len;
    p2[Y] = Vector_y(v2) / len;
    p2[Z] = Vector_z(v2) / len;
  }
  else
    print_error_and_line_num("Null length for vector normalization\n", 
		__FILE__, __LINE__);
}


public float xcorr_objective_with_def(Volume d1,
                                      Volume d2,
                                      Volume m1,
                                      Volume m2, 
                                      Arg_Data *globals)
{

  VectorR
    vector_step;

  PointR
    starting_position,
    slice,
    row,
    col,
    pos2,
    voxel;

  double
    tx,ty,tz;
  int
    r,c,s;

  Real
    value1, value2;
  
  Real
    s1,s2,s3;                   /* to store the sums for f1,f2,f3 */
  float 
    result;				/* the result */
  int 
    count1,count2;


  fill_Point( starting_position, globals->start[X], globals->start[Y], globals->start[Z]);

  s1 = s2 = s3 = 0.0;
  count1 = count2 = 0;

  for_inclusive(s,0,globals->count[SLICE_IND]) {

    SCALE_VECTOR( vector_step, globals->directions[SLICE_IND], s);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );

    for_inclusive(r,0,globals->count[ROW_IND]) {
      
      SCALE_VECTOR( vector_step, globals->directions[ROW_IND], r);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */
      for_inclusive(c,0,globals->count[COL_IND]) {
	
	convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);
	
	fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
	  
	  if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &value1 )) {

	    count1++;

	    DO_TRANSFORM(pos2, globals->trans_info.transformation, col);
	    
	    convert_3D_world_to_voxel(d2, Point_x(pos2), Point_y(pos2), Point_z(pos2), &tx, &ty, &tz);
	    
	    fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
	
	    if (point_not_masked(m2, Point_x(pos2), Point_y(pos2), Point_z(pos2))) {
	      
	      if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &value2 )) {


		if (value1 > globals->threshold[0] && value2 > globals->threshold[1] ) {
		  
		  count2++;

		  s1 += value1*value2;
		  s2 += value1*value1;
		  s3 += value2*value2;
		  
		} 
		
	      } /* if voxel in d2 */
	    } /* if point in mask volume two */
	  } /* if voxel in d1 */
	} /* if point in mask volume one */
	
	ADD_POINT_VECTOR( col, col, globals->directions[COL_IND] );
	
      } /* for c */
    } /* for r */
  } /* for s */
  
  result = 1.0 - s1 / (sqrt((double)s2)*sqrt((double)s3));
  
  if (globals->flags.debug) (void)print ("%7d %7d -> %10.8f\n",count1,count2,result);
  
  return (result);
  
}

