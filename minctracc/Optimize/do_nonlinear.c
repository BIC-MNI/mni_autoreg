/* ----------------------------- MNI Header -----------------------------------
@NAME       : nonlin_fit.c
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
@MODIFIED   : Revision 1.16  1995-05-04 14:25:18  louis
@MODIFIED   : compilable version, seems to run a bit with GRID_TRANSFORM, still
@MODIFIED   : needs work for the super sampled volumes... and lots of testing.
@MODIFIED   :
 * Revision 1.14  1995/02/22  08:56:06  louis
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.13  94/06/23  09:18:04  louis
 * working version before testing of increased cpu time when calculating
 * deformations on a particular slice.
 * 
 * Revision 1.12  94/06/21  10:58:32  louis
 * working optimized version of program.  when compiled with -O, this
 * code is approximately 60% faster than the previous version.
 * 
 * 
 * Revision 1.11  94/06/19  15:43:50  louis
 * clean working version of 3D local deformation with simplex optimization
 * (by default) on magnitude data (default).  No more FFT stuff.
 * This working version before change of deformation field in do_nonlinear.c
 * 
 * 
 * Revision 1.10  94/06/19  15:00:20  louis
 * working version with 3D simplex optimization to find local deformation
 * vector.  time for cleanup.  will remove all fft stuff from code in 
 * next version.
 * 
 * 
 * Revision 1.9  94/06/18  12:29:12  louis
 * temporary version.  trying to find previously working simplex
 * method.
 * 
 * Revision 1.8  94/06/17  10:29:20  louis
 * this version has both a simplex method and convolution method 
 * working to find the best local offset for a given lattice node.
 * 
 * The simplex method works well, is stable and yields good results.
 * The convolution method also works, but the discretization of the 
 * possible deformed vector yields results that are not always stable.
 * 
 * The conv technique is only 2 times faster than the simplex method, 
 * and therefore does not warrent further testing....
 * 
 * I will now modify the simplex method:
 *    - use nearest neighbour interpolation only for the moving lattice
 *    - use tri-cubic for the non-moving lattice
 *    - swap the moving and stable grids.  now the source grid moves,
 *      the target (deformed) grid should move to increase the effective
 *      resolution of the deformation vector (over the linear interpolant
 *      limit).
 * 
 * 
 * Revision 1.7  94/06/15  09:46:47  louis
 * non-working version using FFT for local def.
 *
 * Revision 1.6  94/06/06  18:45:24  louis
 * working version: clamp and blur of deformation lattice now ensures
 * a smooth recovered deformation.  Unfortunately, the r = cost-similarity
 * function used in the optimization is too heavy on the cost_fn.  This has
 * to get fixed...
 * 
 * 
 * Revision 1.5  94/06/06  09:32:41  louis
 * working version: 2d deformations based on local neighbourhood correlation.
 * numerous small bugs over 1.4.  Major bug fix: the routine will now 
 * properly calculate the additional warp (instead of the absolute amount)
 * that needs to be added.  This was due to a mis-type in a call to
 * go_get_values_with_offset.  It was being called with points from the
 * target volume, when they should have been from the source vol.
 * 
 * Some fixes still need to be done: smoothing, clamp 1st deriv, balence
 * of r = cost +  similarity.
 * 
 * Revision 1.4  94/06/02  20:12:07  louis
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
 * Revision 1.3  94/05/28  15:57:33  louis
 * working version, before removing smoothing and adding springs!
 * 
 * Revision 1.2  94/04/06  11:47:47  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.1  94/02/21  16:33:41  louis
 * Initial revision
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/do_nonlinear.c,v 1.16 1995-05-04 14:25:18 louis Exp $";
#endif



#include <volume_io.h>		/* tools to deal with volumes */
#include "arg_data.h"		/* definition of the global data struct */
#include "deform_field.h"	/* my definition of the deformation
				   field.  This should be removed after
				   inclusion of grad_transforms */
#include "local_macros.h"	
#include <print_error.h>	/* def of print_error_and_..  */
#include <limits.h>
#include <recipes.h>
#include "deform_support.h"	/* prototypes for routines called
				   from deformation procedures.   */

#include <sys/types.h>		/* for timing the deformaitons */
#include <time.h>
time_t time(time_t *tloc);


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
static Real     Gsimplex_size;
static Real     Gcost_radius;
static General_transform 
                *Glinear_transform;
static  Volume  Gsuper_dx,Gsuper_dy,Gsuper_dz;
static float    Gsqrt_s1;

static Arg_Data *Gglobals;

static char *my_XYZ_dim_names[] = { MIxspace, MIyspace, MIzspace };


static float		/* these are used to communicate to the correlation */
  *Ga1xyz,		/* functions over top the SIMPLEX optimization */
  *Ga2xyz,		/* routine */
  *TX, *TY, *TZ,
  *SX, *SY, *SZ;
static int
  Glen;


extern double iteration_weight;
extern double similarity_cost_ratio;
extern int    iteration_limit;
extern int    number_dimensions;
extern double ftol;

#define ABSOLUTE_MAX_DEFORMATION 50.0

/* prototypes */


int amoeba2(float **p, 
	    float y[], 
	    int ndim, 
	    float ftol, 
	    float (*funk)(), 
	    int *nfunk);

public void    build_target_lattice1(float px[], float py[], float pz[],
				    float tx[], float ty[], float tz[],
				    int len);

public void    build_target_lattice2(float px[], float py[], float pz[],
				    float tx[], float ty[], float tz[],
				    int len);


private Real cost_fn(float x, float y, float z, Real max_length);

private float xcorr_fitting_function(float *x);

public Status save_deform_data(Volume dx,
			       Volume dy,
			       Volume dz,
			       char *name,
			       char *history);

private Real optimize_3D_deformation_for_single_node(Real spacing, Real threshold1, 
						     Real src_x, Real src_y, Real src_z,
						     Real mx, Real my, Real mz,
						     Real *def_x, Real *def_y, Real *def_z,
						     int iteration, int total_iters,
						     int *nfunks,
						     int ndim);

private BOOLEAN get_best_start_from_neighbours(Real threshold1, 
					       Real wx, Real wy, Real wz,
					       Real mx, Real my, Real mz,
					       Real *tx, Real *ty, Real *tz,
					       Real *d1x_dir, Real *d1y_dir, Real *d1z_dir,
					       Real *def_x, Real *def_y, Real *def_z);
  



/**************************************************************************/

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
  Volume			/* these are to be deleted !!! */
    additional_dx,additional_dy,additional_dz;
  Deform_field 
    *current;

  General_transform
    *all_until_last, 
    *additional_warp,
    *current_warp;

  Volume
    additional_vol,
    current_vol,
    additional_mag;
  
  long
    timer1,timer2,
    nfunk_total;

  int 
    additional_count[MAX_DIMENSIONS],
    current_count[MAX_DIMENSIONS],
    mag_count[MAX_DIMENSIONS],
    sizes[MAX_DIMENSIONS],    
    index[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    start[MAX_DIMENSIONS], 
    end[MAX_DIMENSIONS],
    iters,
    i,j,k,m,
    nodes_done, nodes_tried, nodes_seen, matching_type, over,
    nfunks,
    nfunk1, nodes1;

  Real 
    voxel[MAX_DIMENSIONS],
    steps[MAX_DIMENSIONS],
    current_steps[MAX_DIMENSIONS],
    steps_data[MAX_DIMENSIONS],
    min, max, sum, sum2, mag, mean_disp_mag, std, var, 
    val_x, val_y, val_z,
    def_vector[3],
    wx,wy,wz,
    mx,my,mz,
    tx,ty,tz, 
    displace,
    zero,
    threshold1,
    threshold2,
    result;
  progress_struct
    progress;

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


  Ga1xyz = vector(1,512);	/* allocate space for the global data for */
  Ga2xyz = vector(1,512);	/* the local neighborhood lattice values */
  
  SX = vector(1,512);		/* and coordinates in source volume  */
  SY = vector(1,512);
  SZ = vector(1,512);
  TX = vector(1,512);		/* and coordinates in target volume  */
  TY = vector(1,512);
  TZ = vector(1,512);



				/* split the total transformation into
				   the first linear part and the last
				   non-linear def.                    */
  split_up_the_transformation(globals->trans_info.transformation,
			      &all_until_last,
			      &current_warp);

				/* exit if no deformation */
  if (current_warp == (General_transform *)NULL) { 
    print_error_and_line_num("Cannot find the deformation field transform to optimize",
		__FILE__, __LINE__);

  }
				/* set up linear part of transformation   */
  Glinear_transform = all_until_last; 

				/* print some debugging info    */
  if (globals->flags.debug) {	
    print("orig transform is %d long\n",
	  get_n_concated_transforms(globals->trans_info.transformation));
    print("all_until_last is %d long\n",
	  get_n_concated_transforms(all_until_last));
  }

				/* make a copy of the current warp 
				   that will be used to store the
				   additional deformation needed to
				   optimize the current warp.      */
  ALLOC(additional_warp,1);
  copy_general_transform(current_warp, additional_warp);
  current_vol = current_warp->displacement_volume;
  additional_vol = additional_warp->displacement_volume;

  get_volume_sizes(additional_vol, additional_count);
  get_volume_XYZV_indices(additional_vol, xyzv);

				/* reset the additional warp to zero */
  zero = CONVERT_VALUE_TO_VOXEL(additional_vol, 0.0);	
  for_less (i, 0, additional_count[0])
    for_less (j, 0, additional_count[1])
      for_less (k, 0, additional_count[2])
	for_less (m, 0, additional_count[3])
	  set_volume_real_value(additional_vol, i,j,k,m,0,zero);

				/* build a temporary volume that will
				   be used to store the magnitude of
				   the deformation at each iteration */

  additional_mag = create_volume(3, my_XYZ_dim_names, NC_SHORT, TRUE, 0.0, 0.0);
  for_less(i,0,3)
    mag_count[i] = additional_count[ xyzv[i] ];
  set_volume_sizes(additional_mag, mag_count);
  alloc_volume_data(additional_mag);
  set_volume_real_range(additional_mag, 0.0, ABSOLUTE_MAX_DEFORMATION);

				/* reset additional mag to zero */
  zero = CONVERT_VALUE_TO_VOXEL(additional_mag, 0.0);	
  for_less (i, 0, mag_count[0])
    for_less (j, 0, mag_count[1])
      for_less (k, 0, mag_count[2])
	  set_volume_real_value(additional_mag, i,j,k,0,0,zero);
  mean_disp_mag = 0.0;

				/* test simplex size against size
				   of data voxels */

  get_volume_separations(additional_vol, steps);
  get_volume_separations(d2_dxyz, steps_data);

  if (steps_data[0]!=0.0) {
    Gsimplex_size= ABS(steps[xyzv[X]]) / 2.0; /* /steps_data[0]); */
    if (ABS(Gsimplex_size) < ABS(steps_data[0])) {
      print ("*** WARNING ***\n");
      print ("Simplex size will be smaller than data voxel size (%f < %f)\n",
	     Gsimplex_size,steps_data[0]);
    }
  }
  else
    print_error_and_line_num("Zero step size for gradient data2: %f %f %f\n", 
		__FILE__, __LINE__,steps_data[0],steps_data[1],steps_data[2]);


				/* set up the loop limits to
				   step through the deformation field,
				   node by node.                        */

  get_voxel_spatial_loop_limits(additional_vol, 
				start, end);
  if (globals->flags.debug) {
    print("loop: (%d %d) (%d %d) (%d %d)\n",
	  start[0],end[0],start[1],end[1],start[2],end[2]);
  }

				/* build a super-sampled version of the
				   current transformation, if needed     */
  if (globals->trans_info.use_super) {
    make_super_sampled_data(current->dx, &Gsuper_dx,  &Gsuper_dy,  &Gsuper_dz);
  }
  else {
    /*
       Gsuper_dx = current->dx;
       Gsuper_dy = current->dy;
       Gsuper_dz = current->dz;
    */
  }
				/* set up other parameters needed
				   for non linear fitting */

  Gcost_radius = 8*Gsimplex_size*Gsimplex_size*Gsimplex_size;

  set_feature_value_threshold(d1_dxyz, d2_dxyz,
			      &(globals->threshold[0]), 
			      &(globals->threshold[0]),
			      &threshold1,
			      &threshold2);			      
  if (threshold1<0.0 || threshold2<0.0) {
    print_error_and_line_num("Gradient magnitude threshold error: %f %f\n",
			     __FILE__, __LINE__, 
			     threshold1, threshold2);
  }


  if (globals->flags.debug) {	
    print("Source vol threshold = %f\n", threshold1);
    print("Target vol threshold = %f\n", threshold2);
    print("Iteration limit      = %d\n", iteration_limit);
    print("Iteration weight     = %f\n", iteration_weight);
    print("xyzv                 = %3d %3d %3d %3d \n",
	  xyzv[X], xyzv[Y], xyzv[Z], xyzv[Z+1]);
  }


  /*************************************************************************/
  /*    start the iterations to estimate the warp                     

	for each iteration {
	   for each voxel in the deformation field {
	      get the original coordinate node
	      find the best deformation for the node, taking into 
	        consideration the previous deformation
	      store this best deformation
	   }

	   for each voxel in the deformation field {
 	      add WEIGHTING*best_deformation to the current deformation
	   }

	}
  */

  for_less(iters,0,iteration_limit) {

    if (globals->trans_info.use_super)
      interpolate_super_sampled_data(current->dx, current->dy, current->dz,
				     Gsuper_dx,  Gsuper_dy,  Gsuper_dz,
				     number_dimensions);
  
    
    print("Iteration %2d of %2d\n",iters+1, iteration_limit);

    nodes_done = 0; nodes_tried = 0; nodes_seen=0; 
    displace = 0.0; over = 0;        nfunk_total = 0;

    sum = sum2 = 0.0;
    min = DBL_MAX;
    max = -DBL_MAX;

    initialize_progress_report( &progress, FALSE, 
			       (end[X]-start[X])*(end[Y]-start[Y]) + 1,
			       "Estimating deformations" );

    for_less(i,0,MAX_DIMENSIONS) index[i]=0;

    for_less( index[ xyzv[X] ] , start[ X ], end[ X ]) {

      timer1 = time(NULL);
      nfunk1 = 0; nodes1 = 0;

      for_less( index[ xyzv[Y] ] , start[ Y ], end[ Y ]) {

	for_less( index[ xyzv[Z] ] , start[ Z ], end[ Z ]) {
	 
print ("%3d %3d %3d %3d %3d \n",index[0],index[1],index[2],index[3],index[4]);

	  nodes_seen++;	  
				/* get the lattice coordinate 
				   of the current index node  */
	  for_less(i,0,MAX_DIMENSIONS) voxel[i]=index[i];
	  convert_voxel_to_world(current_vol, 
				 voxel,
				 &wx, &wy, &wz);

				/* get the warp that needs to be added 
				   to get the target world coordinate
				   position */

	  for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
	    def_vector[ index[ xyzv[Z+1] ] ] = 
	      get_volume_real_value(current_vol,
				    index[0],index[1],index[2],index[3],index[4]);

				/* add the warp to get the target 
				   lattice position in world coords */

	  wx += def_vector[X]; wy += def_vector[Y]; wz += def_vector[Z];

	  if (point_not_masked(m2, wx, wy, wz) &&
	      get_value_of_point_in_volume(wx,wy,wz, Gd2_dxyz)>threshold2) {

				/* now get the mean warped position of 
				   the target's neighbours */

	    get_average_warp_of_neighbours(current_warp,
					   index,
					   &mx, &my, &mz);
	    
				/* get the targets homolog in the
				   world coord system of the source
				   data volume                      */

	    general_inverse_transform_point(globals->trans_info.transformation,
					    wx,wy,wz,
					    &tx,&ty,&tz);


				/* find the best deformation for
				   this node                        */

	    def_vector[X] = def_vector[Y] = def_vector[Z] = 0.0;
	    result = optimize_3D_deformation_for_single_node(steps[0], 
							     threshold1,
							     tx,ty,tz,
							     mx, my, mz,
							     &def_vector[X],
							     &def_vector[Y],
							     &def_vector[Z],
							     iters, 
							     iteration_limit,
							     &nfunks,
							     number_dimensions);
	    

	    if (result == -40.0) {
	      nodes_tried++;
	      result = 0.0;
	    } else {
	      if (ABS(result) > 0.95*steps[0])
		over++;

	      nodes_done++;

	      displace += ABS(result);
	      sum += ABS(result);
	      sum2 += ABS(result) * ABS(result);

	      if (ABS(result)>max) max = ABS(result);
	      if (ABS(result)<min) min = ABS(result);

	      nfunk_total += nfunks;

	      nfunk1 += nfunks; nodes1++;

	    }
	    
	  }
	  
	  mag = sqrt(def_vector[X]*def_vector[X] +
		     def_vector[Y]*def_vector[Y] +
		     def_vector[Z]*def_vector[Z]);
	  
	  set_volume_real_value(additional_mag,
				index[xyzv[X]],index[xyzv[Y]],index[xyzv[Z]],0,0,
				mag);

	  for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
	    set_volume_real_value(additional_vol,
				  index[0],index[1],index[2],index[3],index[4],
				  def_vector[ index[ xyzv[Z+1] ] ]);
	}

	update_progress_report( &progress, 
			       (end[Y]-start[Y])*(i-start[X])+
			       (j-start[Y])+1 );
      }
      timer2 = time(NULL);

      if (globals->flags.debug) 
	print ("xslice: (%d : %d) = %d sec -- nodes=%d av funks %f\n",
	       i+1-start[X], 
	       end[X]-start[X]+1, 
	       timer2-timer1, 
	       nodes1,
	       nodes1==0? 0.0:(float)nfunk1/(float)nodes1);
      
    }
    terminate_progress_report( &progress );

    if (globals->flags.debug) {
      if (nodes_done>0) {
	mean_disp_mag = displace/nodes_done;
	var = ((sum2 * nodes_done) - sum*sum) / ((float)nodes_done*(float)(nodes_done-1));
	std = sqrt(var);
	nfunks = nfunk_total / nodes_done;
      }
      else {
	mean_disp_mag=0.0; std = 0.0;
      }
      print ("Nodes seen = %d, tried = %d, done = %d, avg disp = %f +/- %f\n",
	     nodes_seen, nodes_tried, nodes_done, mean_disp_mag, std);
      print ("av nfunks = %d , over = %d, max disp = %f, min disp = %f\n", nfunks, over, max, min);

      nodes_tried = 0;
      for_less(i,0,mag_count[0])
	for_less(j,0,mag_count[1])
	  for_less(k,0,mag_count[2]){
	    mag = get_volume_real_value(additional_mag,i,j,k,0,0);

	    if (mag >= (mean_disp_mag+std))
	      nodes_tried++;
	  }
      print ("there are %d of %d over (mean+1std) disp.\n", nodes_tried, nodes_done);
  

    }
				/* update the current warp, so that the
				   next iteration will use all the data
				   calculated thus far.
				   (the result goes into additional)  */

    add_additional_warp_to_current(additional_warp,
				   current_warp,
				   iteration_weight);

				/* smooth the warp (result->current)  */

    smooth_the_warp(current_warp,
		    additional_warp,
		    additional_mag, 0.0);

    
/*                                 clamp the data so that the 1st derivative of
				   the deformation field does not exceed 1.0*step
				   in magnitude 

     clamp_warp_deriv(current->dx, current->dy, current->dz); 
*/

    
    if (iters<iteration_limit-1) {

      for_less(i,0,MAX_DIMENSIONS) index[i]=0;

      for_less( index[ xyzv[X] ] , start[ X ], end[ X ]) {

	for_less( index[ xyzv[Y] ] , start[ Y ], end[ Y ]) {

	  for_less( index[ xyzv[Z] ] , start[ Z ], end[ Z ]) {

	    
	    mag = get_volume_real_value(additional_mag,
					index[xyzv[X]],index[xyzv[Y]],index[xyzv[Z]],0,0);

	    if (mag >= (mean_disp_mag+std)) {
	      
				/* get the lattice coordinate 
				   of the current index node  */
	      for_less(i,0,MAX_DIMENSIONS) voxel[i]=index[i];
	      convert_voxel_to_world(current_vol, 
				     voxel,
				     &wx, &wy, &wz);

				/* get the warp that needs to be added 
				   to get the target world coordinate
				   position */

	      for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
		def_vector[ index[ xyzv[Z+1] ] ] = 
		  get_volume_real_value(current_vol,
					index[0],index[1],index[2],index[3],index[4]);
	
				/* add the warp to get the target 
				   lattice position in world coords */

	      wx += def_vector[X]; wy += def_vector[Y]; wz += def_vector[Z];
      
	      if ( point_not_masked(m2, wx, wy, wz)) {
		
		/* now get the mean warped position of the target's neighbours */
		
				/* now get the mean warped position of 
				   the target's neighbours */

		get_average_warp_of_neighbours(current_warp,
					       index,
					       &mx, &my, &mz);

				/* get the targets homolog in the
				   world coord system of the source
				   data volume                      */

		general_inverse_transform_point(globals->trans_info.transformation,
						wx,wy,wz,
						&tx,&ty,&tz);
		
				/* find the best deformation for
				   this node                        */

		def_vector[X] = def_vector[Y] = def_vector[Z] = 0.0;

		result = optimize_3D_deformation_for_single_node(steps[0], 
								 threshold1,
								 tx,ty,tz,
								 mx, my, mz,
								 &def_vector[X],
								 &def_vector[Y],
								 &def_vector[Z],
								 iters, iteration_limit,
								 &nfunks,
								 number_dimensions);
		
	      }
	    }
	    for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
	      set_volume_real_value(additional_vol,
				    index[0],index[1],index[2],index[3],index[4],
				    def_vector[ index[ xyzv[Z+1] ] ]);
	    
	  }
	}
      }
      
      
      add_additional_warp_to_current(additional_warp,
				     current_warp,
				     iteration_weight);
      
      /* smooth the warp (result into current->d*)  */
      
      smooth_the_warp(current_warp,
		      additional_warp,
		      additional_mag, (Real)(mean_disp_mag+std));

    }

				/* reset the next iteration's warp. */

    zero = CONVERT_VALUE_TO_VOXEL(additional_vol, 0.0);	
    for_less (i, 0, additional_count[0])
      for_less (j, 0, additional_count[1])
	for_less (k, 0, additional_count[2])
	  for_less (m, 0, additional_count[3])
	    set_volume_real_value(additional_vol, i,j,k,m,0,zero);
    
    if (globals->flags.debug && 
	globals->flags.verbose == 3)
      save_data(globals->filenames.output_trans, 
		iters+1, iteration_limit, 
		current_warp);
    
  }
    
    /* free up allocated temporary deformation volumes */

  if (globals->trans_info.use_super) {
    (void)free_volume_data(Gsuper_dx); (void)delete_volume(Gsuper_dx);
    (void)free_volume_data(Gsuper_dy); (void)delete_volume(Gsuper_dy);
    (void)free_volume_data(Gsuper_dz); (void)delete_volume(Gsuper_dz);
  }

  delete_general_transform(additional_warp);

  (void)free_volume_data(additional_mag); (void)delete_volume(additional_mag);

  free_vector(Ga1xyz ,1,512);
  free_vector(Ga2xyz ,1,512);
  free_vector(TX ,1,512);
  free_vector(TY ,1,512);
  free_vector(TZ ,1,512);
  free_vector(SX ,1,512);
  free_vector(SY ,1,512);
  free_vector(SZ ,1,512);
    

  return (OK);

}





private BOOLEAN get_best_start_from_neighbours(Real threshold1, 
					       Real wx, Real wy, Real wz,
					       Real mx, Real my, Real mz,
					       Real *tx, Real *ty, Real *tz,
					       Real *d1x_dir, Real *d1y_dir, Real *d1z_dir,
					       Real *def_x, Real *def_y, Real *def_z)
     
{
  Real
    mag_normal1,
    nx, ny, nz;


				/* map point from source, forward into target space */

  general_transform_point(Gglobals->trans_info.transformation, wx,wy,wz, tx,ty,tz);

				/* average out target point with the mean position 
				   of its neightbours */

  nx = (*tx+mx)/2.0; 
  ny = (*ty+my)/2.0; 
  nz = (*tz+mz)/2.0; 
				/* what is the deformation needed to achieve this 
				   displacement */
  *def_x = nx - *tx;
  *def_y = ny - *ty;
  *def_z = nz - *tz;
  
  *tx = nx; *ty = ny; *tz = nz;
  
  mag_normal1 = get_value_of_point_in_volume(wx,wy,wz, Gd1_dxyz);

  if (mag_normal1 < threshold1)
    return(FALSE);	
  else
    return(TRUE);
  
}




/***********************************************************/
/* note that the value of the spacing coming in is FWHM/2 for the data
   used to do the correlation. */

private Real optimize_3D_deformation_for_single_node(Real spacing, 
						     Real threshold1, 
						     Real src_x, Real src_y, Real src_z,
						     Real mx, Real my, Real mz,
						     Real *def_x, Real *def_y, Real *def_z,
						     int iteration, int total_iters,
						     int *num_functions,
						     int ndim)
{
  Real
    voxel[3],
    pos[3],
    steps[3],
    simplex_size,
    result,
    tx,ty,tz, 
    xp,yp,zp, 
    d1x_dir, d1y_dir, d1z_dir;
  float 
     *y, **p;
    
    
  int 
    nfunk,len,
    numsteps,
    i,j;

  result = 0.0;			/* assume no additional deformation */
  *num_functions = 0;		/* assume no optimization done      */

  if (!get_best_start_from_neighbours(threshold1, 
				      src_x,  src_y,  src_z,
				      mx,  my,  mz,
				      &tx,  &ty,  &tz,
				      &d1x_dir, &d1y_dir, &d1z_dir,
				      def_x,  def_y,  def_z)) {

				/* then there is no gradient magnitude strong enough
				   to grab onto...*/

    result = -40.0;		/* set result to a flag value used above */

  }
  else {

    /* we now have the info needed to continue...

       src_x, src_y, src_z - point in source volume
       tx, ty, tz          - best point in target volume, so far.
       def_x,def_y,def_z   - currently contains the additional def needed to 
                             take src_x,src_y,src_z mid-way to the neighbour's mean-point
			     (which is now stored in tx,ty,tz).
                    
       */   
    
    
    numsteps = 7;

    len = numsteps+1; 
    for_less(i,1,ndim) {
      len *= (numsteps+1);
    }
				/* get the node point in the source volume taking
				/* into consideration the warp from neighbours  */

    general_inverse_transform_point(Gglobals->trans_info.transformation, 
				    tx,  ty,  tz,
				    &xp, &yp, &zp);


				/* build a lattice of points, in source volume  */
    build_source_lattice(xp, yp, zp, 
			 SX,    SY,    SZ,
			 spacing*3, spacing*3, spacing*3, /* sub-lattice 
							     diameter= 1.5*fwhm */
			 numsteps+1,  numsteps+1,  numsteps+1,
			 ndim, &Glen);

    
				/* map this lattice forward into the target space,
				   using the current transformation, in order to 
				   build a deformed lattice */
    if (Gglobals->trans_info.use_super) {
      build_target_lattice2(SX,SY,SZ,
			    TX,TY,TZ,
			    Glen);
    }
    else {
      build_target_lattice1(SX,SY,SZ,
			    TX,TY,TZ,
			    Glen);
    }

    for_inclusive(i,1,Glen) {
      convert_3D_world_to_voxel(Gd2_dxyz, 
				(Real)TX[i],(Real)TY[i],(Real)TZ[i], 
				&pos[0], &pos[1], &pos[2]);

      if (ndim>2)		/* jiggle the points off center. */
				/* so that nearest-neighbour interpolation can be used */
	TX[i] = pos[0] -1.0 + 2.0*drand48();
      else
	TX[i] = pos[0];
      TY[i] = pos[1] -1.0 + 2.0*drand48();
      TZ[i] = pos[2] -1.0 + 2.0*drand48(); 
    }

				/* build the source lattice (without local neighbour warp),
				   that will be used in the optimization below */
    for_inclusive(i,1,Glen) {
      SX[i] += src_x - xp;
      SY[i] += src_y - yp;
      SZ[i] += src_z - zp;
    }
    
    go_get_samples_in_source(Gd1_dxyz, SX,SY,SZ, Ga1xyz, Glen, 1);

    Gsqrt_s1 = 0.0;
    for_inclusive(i,1,Glen)
      Gsqrt_s1 += Ga1xyz[i]*Ga1xyz[i];

    Gsqrt_s1 = sqrt((double)Gsqrt_s1);

				/* set up SIMPLEX OPTIMIZATION */

    nfunk = 0;
    p = matrix(1,ndim+1,1,ndim);	/* simplex */
    y = vector(1,ndim+1);        /* value of correlation at simplex vertices */
    
    
    p[1][1] = 0.0;
    p[1][2] = 0.0;
    if (ndim > 2)
      p[1][3] = 0.0;
    
    
    for (i=2; i<=(ndim+1); ++i)	/* copy initial guess to all points of simplex */
      for (j=1; j<=ndim; ++j)
	p[i][j] = p[1][j];

    get_volume_separations(Gd2_dxyz,steps);

    simplex_size = Gsimplex_size * (0.5 + 
				    0.5*((Real)(total_iters-iteration)/(Real)total_iters));

    
    p[2][1] += simplex_size;	/* set up all vertices of simplex */
    p[3][2] += simplex_size;
    if (ndim > 2) {
      p[4][3] += simplex_size;
    }
    
    for (i=1; i<=(ndim+1); ++i)	{ /* set up value of correlation at all points of simplex */
      y[i] = xcorr_fitting_function(p[i]);
    }


    if (amoeba2(p,y,ndim,ftol,xcorr_fitting_function,&nfunk)) {    /* do optimization */

      if ( y[1] < y[2] + 0.00001 )
	i=1;
      else
	i=2;
      
      if ( y[i] > y[3] + 0.00001)
	i=3;
      
      if ((ndim > 2) && ( y[i] > y[4] + 0.00001))
	i=4;

      *num_functions = nfunk;
      

      convert_3D_world_to_voxel(Gd2_dxyz, tx,ty,tz, &voxel[0], &voxel[1], &voxel[2]);

      if (ndim>2)
	convert_3D_voxel_to_world(Gd2_dxyz, 
				  (Real)(voxel[0]+p[i][3]), 
				  (Real)(voxel[1]+p[i][2]), 
				  (Real)(voxel[2]+p[i][1]),
				  &pos[0], &pos[1], &pos[2]);
      else
	convert_3D_voxel_to_world(Gd2_dxyz, 
				  (Real)(voxel[0]), 
				  (Real)(voxel[1]+p[i][2]), 
				  (Real)(voxel[2]+p[i][1]),
				  &pos[0], &pos[1], &pos[2]);


      *def_x += pos[0]-tx;
      *def_y += pos[1]-ty;
      if (ndim > 2)
	*def_z += pos[2]-tz;
      else
	*def_z = 0.0;
      
      result = sqrt((*def_x * *def_x) + (*def_y * *def_y) + (*def_z * *def_z)) ;      

/*******begin debug stuff*************

      if (result>2.0) {

	p[1][1] = p[i][1];
	p[1][2] = p[i][2];
	if (ndim > 2)
	  p[1][3] = p[i][3];
	
	print ("dx=[");
	for_inclusive(i,-4,4) {
	  p[2][1] = p[1][1];
	  p[2][2] = p[1][2];
	  if (ndim > 2)
	    p[2][3] = p[1][3];
	  
	  p[2][1] += (float)i * Gsimplex_size/4.0;
	  
	  y[1] = xcorr_fitting_function(p[2]);
	  print ("%9.7f ",y[1]);
	}
	print ("];\n");
	
	print ("dy=[");
	for_inclusive(i,-4,4) {
	  p[2][1] = p[1][1];
	  p[2][2] = p[1][2];
	  if (ndim > 2)
	    p[2][3] = p[1][3];
	  
	  p[2][2] += (float)i * Gsimplex_size/4.0;
	  
	  y[1] = xcorr_fitting_function(p[2]);
	  print ("%9.7f ",y[1]);
	}
	print ("];\n");
	
	if (ndim>2) {
	  print ("dz=[");
	  for_inclusive(i,-4,4) {
	    p[2][1] = p[1][1];
	    p[2][2] = p[1][2];
	    if (ndim > 2)
	      p[2][3] = p[1][3];
	    
	    p[2][3] += (float)i * Gsimplex_size/4.0;
	    
	    y[1] = xcorr_fitting_function(p[2]);
	    print ("%9.7f ",y[1]);
	  }
	print ("];\n");
	}
	p[1][1] = 0.0;
	p[1][2] = 0.0;
	if (ndim > 2)
	  p[1][3] = 0.0;
	
	print ("dx=[");
	for_inclusive(i,-4,4) {
	  p[2][1] = p[1][1];
	  p[2][2] = p[1][2];
	  if (ndim > 2)
	    p[2][3] = p[1][3];
	  
	  p[2][1] += (float)i * Gsimplex_size/4.0;
	  
	  y[1] = xcorr_fitting_function(p[2]);
	  print ("%9.7f ",y[1]);
	}
	print ("];\n");
	
	print ("dy=[");
	for_inclusive(i,-4,4) {
	  p[2][1] = p[1][1];
	  p[2][2] = p[1][2];
	  if (ndim > 2)
	    p[2][3] = p[1][3];
	  
	  p[2][2] += (float)i * Gsimplex_size/4.0;
	  
	  y[1] = xcorr_fitting_function(p[2]);
	  print ("%9.7f ",y[1]);
	}
	print ("];\n");
	
	if (ndim>2) {
	  print ("dz=[");
	  for_inclusive(i,-4,4) {
	    p[2][1] = p[1][1];
	    p[2][2] = p[1][2];
	    if (ndim > 2)
	      p[2][3] = p[1][3];
	    
	    p[2][3] += (float)i * Gsimplex_size/4.0;
	    
	    y[1] = xcorr_fitting_function(p[2]);
	    print ("%9.7f ",y[1]);
	  }
	print ("];\n");
	}

	print ("\n");


      }

*******end debug stuff  *************/



    }
    else {
      *num_functions = 0;
      *def_x += 0.0;
      *def_y += 0.0;
      *def_z += 0.0;     
      result = 0.0;
    }

    free_matrix(p, 1,ndim+1,1,ndim);    
    free_vector(y, 1,ndim+1);  
    

  }
  
  return(result);
}


/* Build the target lattice by transforming the source points through the
   current non-linear transformation stored in:

        Gglobals->trans_info.transformation                        

*/
public void    build_target_lattice1(float px[], float py[], float pz[],
				     float tx[], float ty[], float tz[],
				     int len)
{
  int i;
  Real x,y,z;


  for_inclusive(i,1,len) {
    general_transform_point(Gglobals->trans_info.transformation, 
			    (Real)px[i],(Real) py[i], (Real)pz[i], 
			    &x, &y, &z);
    tx[i] = (float)x;
    ty[i] = (float)y;
    tz[i] = (float)z;
  }

}

/* Build the target lattice by transforming the source points through the
   current non-linear transformation stored in:

        Glinear_transform and Gsuper_d{x,y,z} deformation volumes  

*/
public void    build_target_lattice2(float px[], float py[], float pz[],
				     float tx[], float ty[], float tz[],
				     int len)
{
  int i,sizes[3];
  Real x,y,z;
  Real vx,vy,vz;
  Real dx,dy,dz;
  long 
     ind0, ind1, ind2;

  get_volume_sizes(Gsuper_dx,sizes);

  for_inclusive(i,1,len) {
    general_transform_point(Glinear_transform,
			    (Real)px[i],(Real) py[i], (Real)pz[i], 
			    &x, &y, &z);
    convert_3D_world_to_voxel(Gsuper_dx, x,y,z, &vx, &vy, &vz);

    if ((vx >= -0.5) && (vx < sizes[0]-0.5) &&
	(vy >= -0.5) && (vy < sizes[1]-0.5) &&
	(vz >= -0.5) && (vz < sizes[2]-0.5)    ) {

      ind0 = (long) (vx+0.5);
      ind1 = (long) (vy+0.5);
      ind2 = (long) (vz+0.5);

      GET_VALUE_3D( dx ,  Gsuper_dx, ind0  , ind1  , ind2  );
      GET_VALUE_3D( dy ,  Gsuper_dy, ind0  , ind1  , ind2  );
      GET_VALUE_3D( dz ,  Gsuper_dz, ind0  , ind1  , ind2  );

      x += dx;
      y += dy;
      z += dz;
    }


    tx[i] = (float)x;
    ty[i] = (float)y;
    tz[i] = (float)z;

  }

}




private Real cost_fn(float x, float y, float z, Real max_length)
{
  Real v2,v,d;

  v2 = x*x + y*y + z*z;
  v = sqrt(v2);

  v *= v2;

  if (v<max_length)
    d = 0.2 * v / (max_length - v);
  else
    d = FLT_MAX;

  return(d);
}

private float xcorr_fitting_function(float *d)

{
   float
      similarity,cost, r;

   similarity = go_get_samples_with_offset(Gd2_dxyz, TX,TY,TZ,
					   d[1], d[2], d[3],
					   Glen, Gsqrt_s1, Ga1xyz);

   cost = (float)cost_fn(d[1], d[2], d[3], Gcost_radius);

   r = 1.0 - similarity*similarity_cost_ratio + cost*(1.0-similarity_cost_ratio);
   
   return(r);
}





