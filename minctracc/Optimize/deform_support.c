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
@MODIFIED   : Revision 1.15  1996-06-24 13:54:58  louis
@MODIFIED   : working version, although there may be a bug in the smoohting code.
@MODIFIED   :
 * Revision 1.14  1996/04/01  09:16:29  louis
 * removed the code to super sample the deformation field,
 * and placed it in super_sample_def.c.
 *
 * Revision 1.13  1996/04/01  09:02:10  louis
 * added optimized code to super-sample the deformation field for the special
 * case when super_sample = 2.  interpolate_super_sampled_data_by2() is
 * approximately 7 times faster than interpolate_super_sampled_data().
 *
 * Revision 1.12  1996/03/25  10:33:15  louis
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
 * Revision 1.11  1996/03/07  13:25:19  louis
 * small reorganisation of procedures and working version of non-isotropic
 * smoothing.
 *
 * Revision 1.10  1995/10/06  09:25:02  louis
 * removed references to line_data.h since it hos not been used in a while.
 *
 * included "constants.h" to have access to NONLIN_* similarity func ids.
 *
 * modified go_get_samples_with_offset to account for different similarity
 * functions.
 *
 * Revision 1.9  1995/09/07  10:05:11  louis
 * All references to numerical recipes routines are being removed.  At this
 * stage, any num rec routine should be local in the file.  All memory
 * allocation calls to vector(), matrix(), free_vector() etc... have been
 * replaced with ALLOC and FREE from the volume_io library of routines.
 *
 * Revision 1.8  1995/06/12  14:29:46  louis
 * working version - 2d,3d w/ simplex and -direct.
 *
 * Revision 1.7  1995/05/04  14:25:18  louis
 * compilable version, seems to run a bit with GRID_TRANSFORM, still
 * needs work for the super sampled volumes... and lots of testing.
 *
 * Revision 1.6  1995/05/02  11:30:27  louis
 * started clean up of code, separation of non used procedures into
 * old_methods.c.  This version was working, but I am now going to
 * rewrite everything to use GRID_TRANSFORM.
 *
 * Revision 1.5  1995/02/22  08:56:06  louis
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/deform_support.c,v 1.15 1996-06-24 13:54:58 louis Exp $";
#endif

#include <config.h>             /* for VOXEL_DATA macro */
#include <limits.h>
#include <volume_io.h>
#include <louis_splines.h>
#include <print_error.h>

#include "point_vector.h"

#include "constants.h"


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

public void init_the_volume_to_zero(Volume volume)
{
    int             v0, v1, v2, v3, v4;
    Real            zero;
  

    zero = CONVERT_VALUE_TO_VOXEL(volume,0.0);
    
    BEGIN_ALL_VOXELS( volume, v0, v1, v2, v3, v4 )
      
      set_volume_voxel_value( volume, v0, v1, v2, v3, v4, zero );
    
    END_ALL_VOXELS

}
 



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
				   if they exist */

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
  case 2: {			/* 3x3x3 */
    
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
  
  for_less(i,0,MAX_DIMENSIONS) {
    index[i]=0;
    start[i] = 0;
    end[i] = 0;
  }
  
  get_voxel_spatial_loop_limits(additional->displacement_volume, start, end);
  start[Z+1] = 0;
  end[Z+1] = 3;
 
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

	  /* then, this node has not been estimated at this iteration,
             so we have to interpolate a deformation value from
             neighbouring nodes, and then apply the proper local
	     smoothing. */

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


public void clamp_warp_deriv(Volume dx, Volume dy, Volume dz)
{
  int clamp_needed, clamped_once, i,j,k,sizes[3];
  Real steps[3], diff, voxel, value1, value2, old2 ;
  progress_struct
    progress;

  get_volume_sizes(dx, sizes);
  get_volume_separations(dx, steps);
  
  initialize_progress_report( &progress, FALSE, sizes[0]+sizes[1]+sizes[2] + 1,
			     "Deriv check" );

  clamped_once = FALSE;
  voxel = 0.0;

  if (sizes[0]>2)
  for_less(j,0,sizes[1]) {
    for_less(k,0,sizes[2]){      

      
      do {
	
	clamp_needed = FALSE;
	
	GET_VOXEL_3D(voxel, dz, 0 ,j,k); value1 = CONVERT_VOXEL_TO_VALUE( dz, voxel ); 
	
	for_less(i,1,sizes[0]) {
	  GET_VOXEL_3D(voxel, dz, i,j,k);   value2 = CONVERT_VOXEL_TO_VALUE( dz, voxel ); 
	  old2 = value2;
	  
	  diff = value1 - (value2+steps[0]);
	  
	  if (diff>0.0) {
	    clamp_needed = TRUE; clamped_once = TRUE;
	    
	    value1 -= diff * DERIV_FRAC;
	    value2 += diff * DERIV_FRAC;
	    
	    voxel = CONVERT_VALUE_TO_VOXEL( dz, value1); SET_VOXEL_3D(dz, i-1,j,k, voxel);
	    voxel = CONVERT_VALUE_TO_VOXEL( dz, value2); SET_VOXEL_3D(dz, i,j,k, voxel);
	  }
	  
	  value1 = old2;
	}

      } while (clamp_needed); 
    }

    update_progress_report( &progress, j+1 );
  }


  if (sizes[1]>2)
  for_less(k,0,sizes[2]){      
    for_less(i,0,sizes[0]) {

      
      do { 
	
      clamp_needed = FALSE;

	GET_VOXEL_3D(voxel, dy, i ,0,k); value1 = CONVERT_VOXEL_TO_VALUE( dy, voxel ); 
	
	for_less(j,1,sizes[1]) {
	  GET_VOXEL_3D(voxel, dy, i,j,k);   value2 = CONVERT_VOXEL_TO_VALUE( dy, voxel ); 
	  old2 = value2;
	  
	  diff = value1 - (value2+steps[1]);
	  
	  if (diff>0.0) {
	    clamp_needed = TRUE; clamped_once = TRUE;

	    value1 -= diff * DERIV_FRAC;
	    value2 += diff * DERIV_FRAC;
	    
	    voxel = CONVERT_VALUE_TO_VOXEL( dy, value1); SET_VOXEL_3D(dy, i,j-1,k, voxel);
	    voxel = CONVERT_VALUE_TO_VOXEL( dy, value2); SET_VOXEL_3D(dy, i,j,k, voxel);
	  }
	  
	  value1 = old2;
	}

      } while (clamp_needed); 
    }
    update_progress_report( &progress, sizes[1]+k+1 );
  }


  if (sizes[2]>2)
  for_less(i,0,sizes[0]) {
    for_less(j,0,sizes[1]) {

      
      do {
	
      clamp_needed = FALSE;

	GET_VOXEL_3D(voxel, dx, i ,j,0); value1 = CONVERT_VOXEL_TO_VALUE( dx, voxel ); 
	
	for_less(k,0,sizes[2]){      
	  GET_VOXEL_3D(voxel, dx, i,j,k);   value2 = CONVERT_VOXEL_TO_VALUE( dx, voxel ); 
	  old2 = value2;
	  
	  diff = value1 - (value2+steps[2]);
	  
	  if (diff>0.0) {
	    clamp_needed = TRUE; clamped_once = TRUE;

	    value1 -= diff * DERIV_FRAC;
	    value2 += diff * DERIV_FRAC;
	    
	    voxel = CONVERT_VALUE_TO_VOXEL( dx, value1); SET_VOXEL_3D(dx, i,j,k-1, voxel);
	    voxel = CONVERT_VALUE_TO_VOXEL( dx, value2); SET_VOXEL_3D(dx, i,j,k, voxel);
	  }
	  
	  value1 = old2;
	}

      } while (clamp_needed);
    }
    update_progress_report( &progress, sizes[2]+sizes[1]+i+1 );

  }


  if (clamped_once)
    print ("Clamped once!\n");

}



/* copied from tricubic_interpolant in interpolation.c */

public int tricubic_slice_interpolant(Volume volume, 
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

/*************************************************************************/
/* debug procedure called to save the deformation at each iterative step */

public void save_data(char *basename, int i, int j,
		      General_transform *transform)
{

  Status 
    status;
  STRING 
    comments,name;
  FILE *file;

  (void)sprintf(comments,"step %d of %d of the non-linear estimation",i,j);
  
  (void)sprintf(name,"%s%d",basename,i);
  status = open_file_with_default_suffix(name,
					 get_default_transform_file_suffix(),
					 WRITE_FILE, ASCII_FORMAT, &file );
  
  if( status == OK )
    status = output_transform(file,
			      basename,
			      &i,
			      comments,
			      transform);
  
  if( status == OK )
    status = close_file( file );
  
  if (status!=OK)
    print ("Error saving %s%d\n",basename,i);
}

/*********************************************************************** 
   build a regular (2D) 3D lattice of coordinates to represent the
   local (circular) spherical neighbourhood surrounding the point
   x,y,z.

   - *length coordinates are returned in PX[], PY[], PZ[]. 
   - The radius of the lattice is defined by width_{x,y,z}.
   - The equivalent retangular lattice has (nx)(ny)(nz) samples,
     but the round lattice returned has *length samples.  */

public void    build_source_lattice(Real x, Real y, Real z,
				    float PX[], float PY[], float PZ[],
				    Real width_x, Real width_y, Real width_z, 
				    int nx, int ny, int nz,
				    int ndim, int *length)
{
  int 
    c, 
    i,j,k;
  float 
    radius_squared,
    tx,ty,tz;

  *length = 0; 
  c = 1;
  radius_squared = 0.55 * 0.55;	/* a bit bigger than .5^2 */
  
  if (ndim==2) {
    for_less(i,0,nx)
      for_less(j,0,ny) {
	tx = -0.5 + (float)(i)/(float)(nx-1);
	ty = -0.5 + (float)(j)/(float)(ny-1);
	if ((tx*tx + ty*ty) <= radius_squared) {
	  PX[c] = (float)(x + width_x * tx);
	  PY[c] = (float)(y + width_y * ty);
	  PZ[c] = (float)z;
	  c++;
	  (*length)++;
	}
      }
  }
  else {
    for_less(i,0,nx)
      for_less(j,0,ny)
	for_less(k,0,nz) {
	  tx = -0.5 + (float)(i)/(float)(nx-1);
	  ty = -0.5 + (float)(j)/(float)(ny-1);
	  tz = -0.5 + (float)(k)/(float)(nz-1);
	  if ((tx*tx + ty*ty + tz*tz) <= radius_squared) {
	    PX[c] = (float)(x + width_x * tx);
	    PY[c] = (float)(y + width_y * ty);
	    PZ[c] = (float)(z + width_z * tz);
	    c++;
	    (*length)++;
	  }
	}
  }
}

/*********************************************************************** 
   use the world coordinates stored in x[],y[],z[] to interpolate len
   samples from the volume 'data' */

public void go_get_samples_in_source(Volume data,
				     float x[], float y[], float z[],
				     float samples[],
				     int len,
				     int inter_type) 
{
  int 
    c;
  Real 
    val[MAX_DIMENSIONS];
  
  
  for_inclusive(c,1,len) {
    evaluate_volume_in_world(data,
			     (Real)x[c], (Real)y[c], (Real)z[c], 
			     inter_type,
			     TRUE,
			     0.0, 
			     val,
			     NULL, NULL, NULL, 
			     NULL, NULL, NULL, 
			     NULL, NULL, NULL );
  
    samples[c] = val[0];
  }
  
}

/*********************************************************************** 
   use the list of voxel coordinates stored in x[], y[], z[] and the
   voxel offset stored in dx, dy, dz to interpolate len samples from
   the volume 'data' 

   return the value of the normalized correlation between the list of
   values in *a1 and the values interpolated from data.

   note: the volume is assumed to be in x,y,z order.

   actually, the order does not matter, except that dx corresponds to
   the displacement along the 1st dimension x[], dy corresponds to
   second dimension y[] and dz to z[], the last.  When doing 2D
   processing, the first dimension (x[]) is assumed to be the slowest
   varying, and the deformation in dx=0.

   CAVEAT 1: only nearest neighbour and tri-linear interpolation
             are supported.
   CAVEAT 2: only Volume data types of UNSIGNED_BYTE, SIGNED_SHORT, and
             UNSIGNED_SHORT are supported.

*/

public float go_get_samples_with_offset(Volume data,
					float *x, float *y, float *z,
					Real  dx, Real  dy, Real dz,
					int obj_func,
					int len, 
					float sqrt_s1, float *a1,
					BOOLEAN use_nearest_neighbour)
{
  float
    tmp,
    sample, r,
    s1,s3;			/* to store the sums for f1,f2,f3 */
  int 
    sizes[3],
    ind0, ind1, ind2, 
    c;  
  int xs,ys,zs;
  float
    f_trans, f_scale;

  static double v0, v1, v2;
  static double f0, f1, f2, r0, r1, r2, r1r2, r1f2, f1r2, f1f2;
  static double v000, v001, v010, v011, v100, v101, v110, v111;

  unsigned char  ***byte_ptr;
  unsigned short ***ushort_ptr;
           short ***sshort_ptr;

  Data_types 
    data_type;
  

  data_type = get_volume_data_type(data);
  get_volume_sizes(data, sizes);  
  xs = sizes[0];  
  ys = sizes[1];  
  zs = sizes[2];

  f_trans = data->real_value_translation;
  f_scale = data->real_value_scale;


  s1 = 0.0;
  s3 = 0.0;

  ++a1; 


  if (use_nearest_neighbour) {
				/* then do fast NN interpolation */

    dx += 0.5;			/* to achieve `rounding' for ind0, ind1 and */
    dy += 0.5;			/* ind2 below */
    dz += 0.5;
    
    switch( data_type ) {
    case UNSIGNED_BYTE:  
      
      byte_ptr = VOXEL_DATA (data);
      
      ++x; ++y; ++z; 
      
      for_inclusive(c,1,len) {
	
	/*
	   this is code to test timing of David's evaluate_volume_in_world()
	   interpolation code.  While it is _VERY_ general, the eval_vol_in_world()'s
	   NN interpolation is approximately 8-9 times slower than the bit of 
	   fast NN code below.
	   
	   Real
	   sampleR,index[5], wx, wy,wz;
	   
	   index[0] = (Real) ( *x++ + dx );
	   index[1] = (Real) ( *y++ + dy );
	   index[2] = (Real) ( *z++ + dz );
	   index[3] = 0.0;
	   index[4] = 0.0;
	   
	   convert_voxel_to_world(data, index, &wx, &wy, &wz );
	   
	   evaluate_volume_in_world(data, wx, wy, wz,
				    0, TRUE, 0.0, &sampleR,
				    NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
	   sample = (float)sampleR;
	*/
	
	ind0 = (int) ( *x++ + dx );
	ind1 = (int) ( *y++ + dy );
	ind2 = (int) ( *z++ + dz );
	
	if (ind0>=0 && ind0<xs &&
	    ind1>=0 && ind1<ys &&
	    ind2>=0 && ind2<zs) {
	  sample = (float)(byte_ptr[ind0][ind1][ind2]) * f_scale + f_trans;
	}
	else
	  sample = 0.0;
		
	
#include "switch_obj_func.c"

      }
      break;
    case SIGNED_SHORT:  
      
      sshort_ptr = VOXEL_DATA (data);
      
      ++x; ++y; ++z; 
      
      for_inclusive(c,1,len) {
	
	ind0 = (int) ( *x++ + dx );
	ind1 = (int) ( *y++ + dy );
	ind2 = (int) ( *z++ + dz );
	
	if (ind0>=0 && ind0<xs &&
	    ind1>=0 && ind1<ys &&
	    ind2>=0 && ind2<zs) {
	  sample = sshort_ptr[ind0][ind1][ind2] * f_scale + f_trans;
	}
	else
	  sample = 0.0;
	
#include "switch_obj_func.c"
	
      }
      break;
    case UNSIGNED_SHORT:  
      
      ushort_ptr = VOXEL_DATA (data);
      
      ++x; ++y; ++z; 
      
      for_inclusive(c,1,len) {
	
	ind0 = (int) ( *x++ + dx );
	ind1 = (int) ( *y++ + dy );
	ind2 = (int) ( *z++ + dz );
	
	if (ind0>=0 && ind0<xs &&
	    ind1>=0 && ind1<ys &&
	    ind2>=0 && ind2<zs) {
	  sample = (float)(ushort_ptr[ind0][ind1][ind2]) * f_scale + f_trans;
	}
	else
	  sample = 0.0;
	
#include "switch_obj_func.c"

      }
      break;
    default:
      print_error_and_line_num("Data type not supported in go_get_samples_with_offset (only signed_byte, signed_short, unsigned_short allowed)",__FILE__, __LINE__);
    }
    
  }
  else {			/* then do fast trilinear interpolation */
    
    switch( data_type ) {
    case UNSIGNED_BYTE:  
      
      byte_ptr = VOXEL_DATA (data);
      
      ++x; ++y; ++z; 
      
      for_inclusive(c,1,len) {
	
	/*  fast tri-linear interplation */
	
	v0 = (Real) ( *x++ + dx );
	v1 = (Real) ( *y++ + dy );
	v2 = (Real) ( *z++ + dz );
	
	ind0 = (int)v0;
	ind1 = (int)v1;
	ind2 = (int)v2;
	
	if (ind0>=0 && ind0<(xs-1) &&
	    ind1>=0 && ind1<(ys-1) &&
	    ind2>=0 && ind2<(zs-1)) {
	  
	  /* get the data */
	  v000 = (Real)(byte_ptr[ind0  ][ind1  ][ind2  ]);
	  v001 = (Real)(byte_ptr[ind0  ][ind1  ][ind2+1]);
	  v010 = (Real)(byte_ptr[ind0  ][ind1+1][ind2  ]);
	  v011 = (Real)(byte_ptr[ind0  ][ind1+1][ind2+1]);
	  v100 = (Real)(byte_ptr[ind0+1][ind1  ][ind2  ]);
	  v101 = (Real)(byte_ptr[ind0+1][ind1  ][ind2+1]);
	  v110 = (Real)(byte_ptr[ind0+1][ind1+1][ind2  ]);
	  v111 = (Real)(byte_ptr[ind0+1][ind1+1][ind2+1]);
	  
	  /* Get the fraction parts */
	  f0 = v0 - ind0;
	  f1 = v1 - ind1;
	  f2 = v2 - ind2;
	  r0 = 1.0 - f0;
	  r1 = 1.0 - f1;
	  r2 = 1.0 - f2;
	  
	  /* Do the interpolation */
	  r1r2 = r1 * r2;
	  r1f2 = r1 * f2;
	  f1r2 = f1 * r2;
	  f1f2 = f1 * f2;
	  
	  sample   = 
	    r0 *  (r1r2 * v000 +
		   r1f2 * v001 +
		   f1r2 * v010 +
		   f1f2 * v011);
	  sample  +=
	    f0 *  (r1r2 * v100 +
		   r1f2 * v101 +
		   f1r2 * v110 +
		   f1f2 * v111);
	  
	  sample = sample * f_scale + f_trans;
	  
	}
	else
	  sample = 0.0;
	
	
#include "switch_obj_func.c"

      }
      break;
    case SIGNED_SHORT:  
      
      sshort_ptr = VOXEL_DATA (data);
      
      ++x; ++y; ++z; 
      
      for_inclusive(c,1,len) {
	
	/*  fast tri-linear interpolation */
	
	v0 = (Real) ( *x++ + dx );
	v1 = (Real) ( *y++ + dy );
	v2 = (Real) ( *z++ + dz );
	
	ind0 = (int)v0;
	ind1 = (int)v1;
	ind2 = (int)v2;
	
	if (ind0>=0 && ind0<(xs-1) &&
	    ind1>=0 && ind1<(ys-1) &&
	    ind2>=0 && ind2<(zs-1)) {
	  
	  /* get the data */
	  v000 = (Real)(sshort_ptr[ind0  ][ind1  ][ind2  ]);
	  v001 = (Real)(sshort_ptr[ind0  ][ind1  ][ind2+1]);
	  v010 = (Real)(sshort_ptr[ind0  ][ind1+1][ind2  ]);
	  v011 = (Real)(sshort_ptr[ind0  ][ind1+1][ind2+1]);
	  v100 = (Real)(sshort_ptr[ind0+1][ind1  ][ind2  ]);
	  v101 = (Real)(sshort_ptr[ind0+1][ind1  ][ind2+1]);
	  v110 = (Real)(sshort_ptr[ind0+1][ind1+1][ind2  ]);
	  v111 = (Real)(sshort_ptr[ind0+1][ind1+1][ind2+1]);
	  
	  /* Get the fraction parts */
	  f0 = v0 - ind0;
	  f1 = v1 - ind1;
	  f2 = v2 - ind2;
	  r0 = 1.0 - f0;
	  r1 = 1.0 - f1;
	  r2 = 1.0 - f2;
	  
	  /* Do the interpolation */
	  r1r2 = r1 * r2;
	  r1f2 = r1 * f2;
	  f1r2 = f1 * r2;
	  f1f2 = f1 * f2;
	  
	  sample   = 
	    r0 *  (r1r2 * v000 +
		   r1f2 * v001 +
		   f1r2 * v010 +
		   f1f2 * v011);
	  sample  +=
	    f0 *  (r1r2 * v100 +
		   r1f2 * v101 +
		   f1r2 * v110 +
		   f1f2 * v111);
	  
	  sample = sample * f_scale + f_trans;
	  
	}
	else
	  sample = 0.0;
	
#include "switch_obj_func.c"
	
      }
      break;
    case UNSIGNED_SHORT:  
      
      ushort_ptr = VOXEL_DATA (data);
      
      ++x; ++y; ++z; 
      
      for_inclusive(c,1,len) {
	
	/*  fast tri-linear interplation */
	
	v0 = (Real) ( *x++ + dx );
	v1 = (Real) ( *y++ + dy );
	v2 = (Real) ( *z++ + dz );
	
	ind0 = (int)v0;
	ind1 = (int)v1;
	ind2 = (int)v2;
	
	if (ind0>=0 && ind0<(xs-1) &&
	    ind1>=0 && ind1<(ys-1) &&
	    ind2>=0 && ind2<(zs-1)) {
	  
	  /* get the data */
	  v000 = (Real)(ushort_ptr[ind0  ][ind1  ][ind2  ]);
	  v001 = (Real)(ushort_ptr[ind0  ][ind1  ][ind2+1]);
	  v010 = (Real)(ushort_ptr[ind0  ][ind1+1][ind2  ]);
	  v011 = (Real)(ushort_ptr[ind0  ][ind1+1][ind2+1]);
	  v100 = (Real)(ushort_ptr[ind0+1][ind1  ][ind2  ]);
	  v101 = (Real)(ushort_ptr[ind0+1][ind1  ][ind2+1]);
	  v110 = (Real)(ushort_ptr[ind0+1][ind1+1][ind2  ]);
	  v111 = (Real)(ushort_ptr[ind0+1][ind1+1][ind2+1]);
	  
	  /* Get the fraction parts */
	  f0 = v0 - ind0;
	  f1 = v1 - ind1;
	  f2 = v2 - ind2;
	  r0 = 1.0 - f0;
	  r1 = 1.0 - f1;
	  r2 = 1.0 - f2;
	  
	  /* Do the interpolation */
	  r1r2 = r1 * r2;
	  r1f2 = r1 * f2;
	  f1r2 = f1 * r2;
	  f1f2 = f1 * f2;
	  
	  sample   = 
	    r0 *  (r1r2 * v000 +
		   r1f2 * v001 +
		   f1r2 * v010 +
		   f1f2 * v011);
	  sample  +=
	    f0 *  (r1r2 * v100 +
		   r1f2 * v101 +
		   f1r2 * v110 +
		   f1f2 * v111);
	  
	  sample = sample * f_scale + f_trans;
	  
	}
	else
	  sample = 0.0;
	
#include "switch_obj_func.c"

      }
      break;
    default:
      print_error_and_line_num("Data type not supported in go_get_samples_with_offset (only signed_byte, signed_short, unsigned_short allowed)",__FILE__, __LINE__);
    }
    
  }
    
  
				/* do the last bits of the similarity function
				   calculation here: */
  r = 0.0;
  switch (obj_func) {
  case NONLIN_XCORR:
    if ( sqrt_s1 < 0.001 && s3 < 0.00001) {
      r = 1.0;
    }
    else {
      if ( sqrt_s1 < 0.001 || s3 < 0.00001) {
	r = 0.0;
      }
      else {
	r = s1 / (sqrt_s1*sqrt((double)s3));
      }
    }
    
    break;
  case NONLIN_DIFF:
  case NONLIN_LABEL:
    r = s1 / sqrt_s1;		/* sqrt_s1 can't be 0, since = Glen */
    break;
  default:
    print_error_and_line_num("Objective function %d not supported in go_get_samples_with_offset",__FILE__, __LINE__,obj_func);
  }
  
  
  return(r);
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
  int 
    i,j,k,sizes[3];
  Real voxel, val, max;

  max = -DBL_MAX;
  get_volume_sizes(dxyz, sizes);
  
  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]){
	
	GET_VOXEL_3D(voxel, dxyz, i,j,k); 
	val = CONVERT_VOXEL_TO_VALUE(dxyz, voxel);

	if (val > max) max = val;

      }
  
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

