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
@MODIFIED   : Revision 1.11  1994-06-19 15:43:50  louis
@MODIFIED   : clean working version of 3D local deformation with simplex optimization
@MODIFIED   : (by default) on magnitude data (default).  No more FFT stuff.
@MODIFIED   : This working version before change of deformation field in do_nonlinear.c
@MODIFIED   :
@MODIFIED   :
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/do_nonlinear.c,v 1.11 1994-06-19 15:43:50 louis Exp $";
#endif


#include <volume_io.h>
#include "arg_data.h"
#include "deform_field.h"
#include "line_data.h"
#include <print_error.h>
#include <limits.h>
#include <recipes.h>


#include <sys/types.h>
#include <time.h>

time_t time(time_t *tloc);


#define BRUTE_FORCE     1
#define SECANT_METHOD   0
#define TEST_LENGTH     1.0

#define EQUAL(a,b) ( ABS( (a) - (b)) < 0.000001)

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

static Arg_Data *Gglobals;

static float			/* these are used to communicate to the correlation */
  *Ga1xyz,			/* functions over top the SIMPLEX optimization */
  *Ga2xyz,			/* routine */
  *TX, *TY, *TZ,
  *SX, *SY, *SZ;
static int
  Glen;


extern double iteration_weight;
extern double similarity_cost_ratio;
extern int    iteration_limit;
extern int    number_dimensions;

/* prototypes */


int amoeba2(float **p, 
	    float y[], 
	    int ndim, 
	    float ftol, 
	    float (*funk)(), 
	    int *nfunk);

public void    build_source_lattice(Real x, Real y, Real z,
				    float PX[], float PY[], float PZ[],
				    Real fwhm_x, Real fwhm_y, Real fwhm_z, 
				    int nx, int ny, int nz,
				    int ndim);

public void    build_target_lattice(float px[], float py[], float pz[],
				    float tx[], float ty[], float tz[],
				    int len);

public void go_get_samples_in_source(Volume data,
			   float x[], float y[], float z[],
			   float samples[],
			   int len,
			   int inter_type);

public void go_get_samples_with_offset(Volume data,
				       float x[], float y[], float z[],
				       float samples[],
				       Real dx, Real dy, Real dz,
				       int len,
				       int inter_type) ;

private float gauss_3d(float c,
		       float fwhm_x,float fwhm_y,float fwhm_z,
		       float mux,float muy,float muz,
		       float x,float y,float z,
		       int ndim);

public float calc_scalar_correlation(float *a1,float *a2, int len);

private Real cost_fn(float x, float y, float z, Real max_length);

private float xcorr_fitting_function(float *x);

public  Status  output_deformation_file(
    char                filename[],
    char                comments[],
    General_transform   *transform );

public Status save_deform_data(Volume dx,
			       Volume dy,
			       Volume dz,
			       char *name,
			       char *history);

public  Real CUBIC_PP(Real f1, Real f2, Real f3, Real f4, Real xpos);
public  Real CUBIC_P(Real f1, Real f2, Real f3, Real f4, Real xpos);

public Real find_offset_to_match(Line_data *m,
				 Line_data *d,
				 Real      limit,
				 int       op_type);

private Real optimize_1D_deformation_for_single_node(Real spacing, 
						     Real threshold1, Real threshold2,
						     Real wx, Real wy, Real wz,
						     Real mx, Real my, Real mz,
						     Real *def_x, Real *def_y, Real *def_z,
						     int match_type);

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
  

private Real get_normalized_world_gradient_direction(Real xworld, Real yworld, Real zworld, 
						     Volume data_x,Volume data_y,Volume data_z,
						     Real *DX_dir, Real *DY_dir, Real *DZ_dir,
						     Real thresh);

private float interp_data_along_gradient(Real dist_from, int inter_type,
					 Real x, Real y, Real z, 
					 Real dx, Real dy, Real dz,
					 Volume data);


public void clamp_warp_deriv(Volume dx, Volume dy, Volume dz);


public void get_average_warp_of_neighbours(Volume dx, Volume dy, Volume dz, 
					   int i,int j,int k,
					   Real *mx, Real *my, Real *mz);

public void smooth_the_warp(Volume smooth_dx, Volume smooth_dy, Volume smooth_dz, 
			    Volume warp_dx, Volume warp_dy, Volume warp_dz,
			    Volume warp_mag, Real thres);

public void add_additional_warp_to_current(Volume dx, Volume dy, Volume dz,
					    Volume adx, Volume ady, Volume adz,
					    Real weight);

public Real get_maximum_magnitude(Volume dx, Volume dy, Volume dz);




private BOOLEAN build_second_derivative_data(Line_data *line, int num_samples, Real spacing,
					     Real wx, Real wy, Real wz, 
					     Real dx_dir, Real dy_dir, Real dz_dir, 
					     Volume intensity);

private BOOLEAN build_first_derivative_data(Line_data *line, int num_samples, Real spacing,
					     Real wx, Real wy, Real wz, 
					     Real dx_dir, Real dy_dir, Real dz_dir, 
					     Volume intensity);

private BOOLEAN build_line_data(Line_data *line, int num_samples, Real spacing,
				Real wx, Real wy, Real wz, 
				Real dx_dir, Real dy_dir, Real dz_dir, 
				Volume intensity);

public void save_data(char *basename, int i, int j,
		      Volume dx,Volume dy,Volume dz);




/******************************************************************************************/

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
    additional_dx,additional_dy,additional_dz,additional_mag;
  Deform_field 
    *current;
  long
    timer1,timer2,
    nfunk_total;
  int 
    loop_start[3], loop_end[3],
    iters,
    i,j,k,
    nodes_done, nodes_tried, nodes_seen, matching_type, over,
    nfunks,
    sizes[3];
  Real 
    min, max, sum, sum2, mag, mean_disp_mag, std, var, 
    voxel,val_x, val_y, val_z,
    steps[3],steps_data[3],
    def_x, def_y, def_z, 
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
  General_transform 
    *all_until_last, *tmp_trans;

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

  current = (Deform_field *)NULL;
  for_less(i,0,get_n_concated_transforms(globals->trans_info.transformation)) {
    if (get_transform_type( get_nth_general_transform(globals->trans_info.transformation,i) ) 
	       == USER_TRANSFORM)
      current = (Deform_field *)get_nth_general_transform(globals->trans_info.transformation,i)->user_data;
  }
  if (current == (Deform_field *)NULL) { /* exit if no deformation */
    print_error("Cannot find the deformation field transform to optimize",
		__FILE__, __LINE__);
  }


  /*******************************************************************************/
  /*   build and allocate the temporary deformation volume needed                */

  additional_dx = copy_volume_definition(current->dx, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  additional_dy = copy_volume_definition(current->dy, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  additional_dz = copy_volume_definition(current->dz, NC_UNSPECIFIED, FALSE, 0.0,0.0);
  additional_mag = copy_volume_definition(current->dz, NC_UNSPECIFIED, FALSE, 0.0,0.0);

  alloc_volume_data(additional_dx);
  alloc_volume_data(additional_dy);
  alloc_volume_data(additional_dz);
  alloc_volume_data(additional_mag);

  get_volume_sizes(additional_dx, sizes);
  get_volume_separations(additional_dx, steps);

  get_volume_separations(d2_dxyz, steps_data);

  if (steps_data[0]!=0.0) {
    Gsimplex_size= ABS(steps[0]/steps_data[0]);
    if (ABS(Gsimplex_size) < ABS(steps_data[0])) {
      print ("*** WARNING ***\n");
      print ("Simplex size will be smaller than data voxel size (%f < %f)\n",
	     Gsimplex_size,steps_data[0]);
    }
  }
  else
    print_error("Zero step size for gradient data2: %f %f %f\n", 
		__FILE__, __LINE__,steps_data[0],steps_data[1],steps_data[2]);


  Gcost_radius = 8*Gsimplex_size*Gsimplex_size*Gsimplex_size;


  /*******************************************************************************/
  /*   initialize this iterations warp                                           */

  zero = CONVERT_VALUE_TO_VOXEL(additional_dx, 0.0);	
  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]){
	SET_VOXEL_3D(additional_dx, i,j,k, zero);
	SET_VOXEL_3D(additional_dy, i,j,k, zero);
	SET_VOXEL_3D(additional_dz, i,j,k, zero);
	SET_VOXEL_3D(additional_mag, i,j,k, zero);
      }
  mean_disp_mag = 0.0;


  /*******************************************************************************/
  /*    set the threshold to be 10% of the maximum gradient magnitude            */
  /*    for each source and target volumes                                       */

  threshold1 = 0.10 * get_maximum_magnitude(d1_dx, d1_dy, d1_dz); 
  threshold2 = 0.10 * get_maximum_magnitude(d2_dx, d2_dy, d2_dz); 

  if (threshold1<0.0 || threshold2<0.0) {
    print_error("Gradient magnitude threshold error: %f %f\n", __FILE__, __LINE__, threshold1, threshold2);
  }
  if (globals->flags.debug) {	
    print("Source vol threshold = %f\n", threshold1);
    print("Target vol threshold = %f\n", threshold2);
    print("Iteration limit      = %d\n", iteration_limit);
    print("Iteration weight     = %f\n", iteration_weight);
  }

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

	}
  */

  for_less(i,0,3) {
    if (sizes[i]>3) {
      loop_start[i] = 1;
      loop_end[i] = sizes[i]-1;
    }
    else {
      loop_start[i]=0;
      loop_end[i] = sizes[i];
    }
  }


/*  loop_start[0] += 7;
  loop_end[0] -= 7;
*/

  if (globals->flags.debug) {
    print("loop: (%d %d) (%d %d) (%d %d)\n",
	  loop_start[0],loop_end[0],loop_start[1],loop_end[1],loop_start[2],loop_end[2]);
  }

  for_less(iters,0,iteration_limit) {

    print("Iteration %2d of %2d\n",iters+1, iteration_limit);

    if (Gglobals->trans_info.use_simplex==TRUE) 
      print ("will use simplex\n");
    else
      print ("will use FFT\n");

    
    /* if (iters<3) matching_type = BRUTE_FORCE; else matching_type = SECANT_METHOD; */

    matching_type = BRUTE_FORCE; /* force matching type to be BRUTE_FORCE for now!!! */

    nodes_done = 0; nodes_tried = 0; nodes_seen=0; displace = 0.0; over = 0;
    nfunk_total = 0;

    sum = sum2 = 0.0;
    min = 1000.0;
    max = -1000.0;

    initialize_progress_report( &progress, FALSE, (loop_end[0]-loop_start[0])*(loop_end[1]-loop_start[1]) + 1,
			       "Estimating deformations" );

    for_less(i,loop_start[0],loop_end[0]) {

      timer1 = time(NULL);

      for_less(j,loop_start[1],loop_end[1]) {
	for_less(k,loop_start[2],loop_end[2]){
	  
	  nodes_seen++;
	  def_x = def_y = def_z = 0.0;
	  
				/* get the lattice coordinate */
	  convert_3D_voxel_to_world(current->dx, 
				    (Real)i, (Real)j, (Real)k,
				    &wx, &wy, &wz);

				/* get the warp to be added to the target point */
	  GET_VOXEL_3D(voxel, current->dx, i,j,k); 
	  val_x = CONVERT_VOXEL_TO_VALUE(current->dx , voxel ); 
	  GET_VOXEL_3D(voxel, current->dy, i,j,k); 
	  val_y = CONVERT_VOXEL_TO_VALUE(current->dy , voxel ); 
	  GET_VOXEL_3D(voxel, current->dz, i,j,k); 
	  val_z = CONVERT_VOXEL_TO_VALUE(current->dz , voxel ); 

				/* add the warp to the target lattice point */
	  wx += val_x; wy += val_y; wz += val_z;

	  if ( point_not_masked(m2, wx, wy, wz)) {

				/* now get the mean warped position of the target's neighbours */

	    get_average_warp_of_neighbours(current->dx, current->dy, current->dz, 
					   i,j,k,
					   &mx, &my, &mz);
	    
	    general_inverse_transform_point(globals->trans_info.transformation,
					    wx,wy,wz,
					    &tx,&ty,&tz);

	    if (Gglobals->trans_info.use_simplex==TRUE)
	      result = optimize_3D_deformation_for_single_node(steps[0], 
							       threshold1,
							       tx,ty,tz,
							       mx, my, mz,
							       &def_x, &def_y, &def_z, 
							       iters, iteration_limit,
							       &nfunks,
							       number_dimensions);
	    else
	      result = optimize_1D_deformation_for_single_node(steps[0], 
							       threshold1,threshold2,
							       tx,ty,tz,
							       mx, my, mz,
							       &def_x, &def_y, &def_z, 
							       matching_type);
	    
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
	    }
	    
	  }
	  
	  mag = sqrt(def_x*def_x + def_y*def_y + def_z*def_z);
	  mag = CONVERT_VALUE_TO_VOXEL(additional_dx, mag); 
	  SET_VOXEL_3D(additional_mag, i,j,k, mag);

	  def_x = CONVERT_VALUE_TO_VOXEL(additional_dx, def_x); 
	  SET_VOXEL_3D(additional_dx, i,j,k, def_x);
	  def_y = CONVERT_VALUE_TO_VOXEL(additional_dy, def_y); 
	  SET_VOXEL_3D(additional_dy, i,j,k, def_y);
	  def_z = CONVERT_VALUE_TO_VOXEL(additional_dz, def_z); 
	  SET_VOXEL_3D(additional_dz, i,j,k, def_z);
	  
	}

	update_progress_report( &progress, 
			       (loop_end[1]-loop_start[1])*(i-loop_start[0])+(j-loop_start[1])+1 );
      }
      timer2 = time(NULL);

print ("time: %d - %d = %d seconds\n",timer2, timer1, timer2-timer1);

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
      for_less(i,0,sizes[0])
	for_less(j,0,sizes[1])
	  for_less(k,0,sizes[2]){
	    GET_VOXEL_3D(voxel, additional_mag, i,j,k); 
	    mag = CONVERT_VOXEL_TO_VALUE(additional_mag , voxel ); 
	    if (mag >= (mean_disp_mag+std))
	      nodes_tried++;
	  }
      print ("there are %d of %d over (mean+1std) disp.\n", nodes_tried, nodes_done);
  

    }


				/* update the current warp, so that the
				   next iteration will use all the data
				   calculated thus far.
				   (the result goes into additional_d* )   */

    add_additional_warp_to_current(additional_dx, additional_dy, additional_dz,
				   current->dx, current->dy, current->dz,
				   iteration_weight);

				/* smooth the warp (result into current->d*)  */

    smooth_the_warp(current->dx, current->dy, current->dz,
		    additional_dx, additional_dy, additional_dz,
		    additional_mag, 0.0);

    
/*                                 clamp the data so that the 1st derivative of
				   the deformation field does not exceed 1.0*step
				   in magnitude 

     clamp_warp_deriv(current->dx, current->dy, current->dz); 
*/

    
    if (iters<iteration_limit-1 && Gglobals->trans_info.use_simplex==TRUE) {
      for_less(i,loop_start[0],loop_end[0]) {
	for_less(j,loop_start[1],loop_end[1]) {
	  for_less(k,loop_start[2],loop_end[2]){
	    
	    def_x = def_y = def_z = 0.0;
	    
	    GET_VOXEL_3D(voxel, additional_mag, i,j,k); 
	    mag = CONVERT_VOXEL_TO_VALUE(additional_mag , voxel ); 
	    if (mag >= (mean_disp_mag+std)) {
	      
	      /* get the lattice coordinate */
	      convert_3D_voxel_to_world(current->dx, 
					(Real)i, (Real)j, (Real)k,
					&wx, &wy, &wz);
	      
	      /* get the warp to be added to the target point */
	      GET_VOXEL_3D(voxel, current->dx, i,j,k); 
	      val_x = CONVERT_VOXEL_TO_VALUE(current->dx , voxel ); 
	      GET_VOXEL_3D(voxel, current->dy, i,j,k); 
	      val_y = CONVERT_VOXEL_TO_VALUE(current->dy , voxel ); 
	      GET_VOXEL_3D(voxel, current->dz, i,j,k); 
	      val_z = CONVERT_VOXEL_TO_VALUE(current->dz , voxel ); 
	      
	      /* add the warp to the target lattice point */
	      wx += val_x; wy += val_y; wz += val_z;
	      
	      if ( point_not_masked(m2, wx, wy, wz)) {
		
		/* now get the mean warped position of the target's neighbours */
		
		get_average_warp_of_neighbours(current->dx, current->dy, current->dz, 
					       i,j,k,
					       &mx, &my, &mz);
		
		general_inverse_transform_point(globals->trans_info.transformation,
						wx,wy,wz,
						&tx,&ty,&tz);
		
		if (Gglobals->trans_info.use_simplex==TRUE)
		  result = optimize_3D_deformation_for_single_node(steps[0], 
								   threshold1,
								   tx,ty,tz,
								   mx, my, mz,
								   &def_x, &def_y, &def_z, 
								   iters, iteration_limit,
								   &nfunks,
								   number_dimensions);
		else
		  result = optimize_1D_deformation_for_single_node(steps[0], 
								   threshold1,threshold2,
								   tx,ty,tz,
								   mx, my, mz,
								   &def_x, &def_y, &def_z, 
								   matching_type);
		
	      }
	    }
	    
	    def_x = CONVERT_VALUE_TO_VOXEL(additional_dx, def_x); 
	    SET_VOXEL_3D(additional_dx, i,j,k, def_x);
	    def_y = CONVERT_VALUE_TO_VOXEL(additional_dy, def_y); 
	    SET_VOXEL_3D(additional_dy, i,j,k, def_y);
	    def_z = CONVERT_VALUE_TO_VOXEL(additional_dz, def_z); 
	    SET_VOXEL_3D(additional_dz, i,j,k, def_z);
	    
	  }
	}
      }
      
      
      add_additional_warp_to_current(additional_dx, additional_dy, additional_dz,
				     current->dx, current->dy, current->dz,
				     iteration_weight);
      
      /* smooth the warp (result into current->d*)  */
      
      smooth_the_warp(current->dx, current->dy, current->dz,
		      additional_dx, additional_dy, additional_dz,
		      additional_mag, (Real)(mean_disp_mag+std));


    }


				/* reset the next iteration's warp. */


    zero = CONVERT_VALUE_TO_VOXEL(additional_dx, 0.0);
    for_less(i,0,sizes[0])
      for_less(j,0,sizes[1])
	for_less(k,0,sizes[2]){
	  SET_VOXEL_3D(additional_dx, i,j,k, zero);
	  SET_VOXEL_3D(additional_dy, i,j,k, zero);
	  SET_VOXEL_3D(additional_dz, i,j,k, zero);

	}
  

    if (globals->flags.debug && 
	globals->flags.verbose == 3)
      save_data(globals->filenames.output_trans, iters+1, iteration_limit, 
		current->dx, current->dy, current->dz);

  }

				/* free up allocated temporary deformation volumes */

  (void)free_volume_data(additional_dx); (void)delete_volume(additional_dx);
  (void)free_volume_data(additional_dy); (void)delete_volume(additional_dy);
  (void)free_volume_data(additional_dz); (void)delete_volume(additional_dz);
  (void)free_volume_data(additional_mag); (void)delete_volume(additional_mag);


  return (OK);



}





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
  
  mag_normal1 = 
    get_normalized_world_gradient_direction(wx,wy,wz, 
					    Gd1_dx, Gd1_dy, Gd1_dz, 
					    d1x_dir, d1y_dir, d1z_dir,
					    threshold1);

  if (mag_normal1 < threshold1)
    return(FALSE);	
  else
    return(TRUE);
  
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

/*******************************************************************************/
/* debugging procedure called to save the deformation at each iterative step  */

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
    ftol, *y, **p;
    
    
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

    Glen = len;
    
    Ga1xyz = vector(1,len);
    Ga2xyz = vector(1,len);
    
    SX = vector(1,len);
    SY = vector(1,len);
    SZ = vector(1,len);
    TX = vector(1,len);
    TY = vector(1,len);
    TZ = vector(1,len);

				/* get the node point in the source volume taking   */
				/* into consideration the warp from neighbours      */

    general_inverse_transform_point(Gglobals->trans_info.transformation, 
				    tx,  ty,  tz,
				    &xp, &yp, &zp);

				/* build a lattice of points, in the source volume  */
    build_source_lattice(xp, yp, zp, 
			 SX,    SY,    SZ,
			 spacing*3, spacing*3, spacing*3, /* lattice size= 1.5*fwhm */
			 numsteps+1,  numsteps+1,  numsteps+1,
			 ndim);
    
				/* map this lattice forward into the target space,
				   using the current transformation, in order to 
				   build a deformed lattice */
    build_target_lattice(SX,SY,SZ,
			 TX,TY,TZ,
			 len);

    for_inclusive(i,1,len) {
      convert_3D_world_to_voxel(Gd2_dxyz, 
				(Real)TX[i],(Real)TY[i],(Real)TZ[i], 
				&pos[0], &pos[1], &pos[2]);
      TX[i] = pos[0] -1.0 + 2.0*drand48(); /* jiggle the point off center. */
      TY[i] = pos[1] -1.0 + 2.0*drand48();
      if (ndim>2)
	TZ[i] = pos[2] -1.0 + 2.0*drand48();
      else
	TZ[i] = pos[2];
    }

				/* build the source lattice (without local neighbour warp),
				   that will be used in the optimization below */
    build_source_lattice(src_x, src_y, src_z,
			 SX,    SY,    SZ,
			 spacing*3, spacing*3, spacing*3, /* lattice size= 1.5*fwhm */
			 numsteps+1,  numsteps+1,  numsteps+1,
			 ndim);
    go_get_samples_in_source(Gd1_dxyz, SX,SY,SZ, Ga1xyz, len, 1);


				/* set up SIMPLEX OPTIMIZATION */

    ftol  = MAX(0.004,spacing/3000);

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

/* print ("corr: %12.9f %12.9f %12.9f  \n",y[1], y[2], y[3]); */
    
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
    
    free_vector(Ga1xyz ,1,len);
    free_vector(Ga2xyz ,1,len);
    free_vector(TX ,1,len);
    free_vector(TY ,1,len);
    free_vector(TZ ,1,len);
    free_vector(SX ,1,len);
    free_vector(SY ,1,len);
    free_vector(SZ ,1,len);
    

  }
  
  return(result);
}


public void    build_source_lattice(Real x, Real y, Real z,
				    float PX[], float PY[], float PZ[],
				    Real width_x, Real width_y, Real width_z, 
				    int nx, int ny, int nz,
				    int ndim)
{
  int c, i,j,k;

  c = 1;
  if (ndim==2) {
    for_less(i,0,nx)
      for_less(j,0,ny) {
	PX[c] = (float)(x + width_x * (-0.5 + (float)(i)/(float)(nx-1)));
	PY[c] = (float)(y + width_y * (-0.5 + (float)(j)/(float)(ny-1)));
        PZ[c] = (float)z;
	c++;
      }
  }
  else {
    for_less(i,0,nx)
      for_less(j,0,ny)
	for_less(k,0,nz) {
	  PX[c] = (float)(x + width_x * (-0.5 + (float)(i)/(float)(nx-1)));
	  PY[c] = (float)(y + width_y * (-0.5 + (float)(j)/(float)(ny-1)));
	  PZ[c] = (float)(z + width_z * (-0.5 + (float)(k)/(float)(nz-1)));
	  c++;
	}
  }
}

public void    build_target_lattice(float px[], float py[], float pz[],
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


public void go_get_samples_in_source(Volume data,
				     float x[], float y[], float z[],
				     float samples[],
				     int len,
				     int inter_type) 
{
  int 
    flag, c;
  Real 
    val;
  PointR
    point;

  for_inclusive(c,1,len) {
    convert_3D_world_to_voxel(data, (Real)x[c], (Real)y[c], (Real)z[c], 
			      &Point_x(point), &Point_y(point), &Point_z(point));

    if (inter_type==0)
	flag = nearest_neighbour_interpolant(data,  &point , &val);
    else {
      if (inter_type==1)
	flag = trilinear_interpolant(data,  &point , &val);
      else
	flag = tricubic_interpolant(data,  &point , &val);
    }

    if (flag)
      samples[c] = (float)val;
    else
      samples[c] = 0.0;
  }

}


public void go_get_samples_with_offset(Volume data,
				       float x[], float y[], float z[],
				       float samples[],
				       Real dx, Real dy, Real dz,
				       int len,
				       int inter_type) 
{
  int 
    flag, c;
  Real 
    val;
  PointR
    point;

  for_inclusive(c,1,len) {

    Point_x(point) = x[c]+dz;
    Point_y(point) = y[c]+dy;
    Point_z(point) = z[c]+dx;
    
    if (inter_type==0)
      flag = nearest_neighbour_interpolant(data,  &point , &val);
    else {
      if (inter_type==1)
	flag = trilinear_interpolant(data,  &point , &val);
      else
	flag = tricubic_interpolant(data,  &point , &val);
    }
    if (flag)
      samples[c] = (float)val;
    else
      samples[c] = (float)0.0;
  }

}


/************************************************************/
/* return the value of the normal dist at x, given c,sigma  */
/* and mu ----   all in mm                                  */
/************************************************************/
private float gauss_3d(float c,
		       float fwhm_x,float fwhm_y,float fwhm_z,
		       float mux,float muy,float muz,
		       float x,float y,float z,
		       int ndim)
{
  float 
    dx,dy,dz,sigma_x,sigma_y,sigma_z,t1,t2,t3,f;
  
  
  if (ndim==2) {
    
    if (EQUAL(fwhm_x,0.0) || EQUAL(fwhm_y,0.0) ) {
      if ( EQUAL(x,mux) && EQUAL(y,muy)  )
	f = c;
      else
	f = 0;
    }
    else {
      
      sigma_x = fwhm_x/2.0;
      sigma_y = fwhm_y/2.0;
      
      dx = mux-x;
      dy = muy-y;
      
      t1 = c; /* /(sigma_x*sigma_y*fsqrt(two_pi3));*/
      
      sigma_x = sigma_x*sigma_x;  /* square the differences and the sigmas */
      sigma_y = sigma_y*sigma_y;
      
      dx = dx*dx;
      dy = dy*dy;
      
      t2 = (-1.0/4.0) * 
	(sigma_y*dx + sigma_x*dy ) / (sigma_x*sigma_y);
      
      t3 = fexp(t2);
      
      f = t1*t3;
      
    }
    
    
  }
  else {
    
    if (EQUAL(fwhm_x,0.0) || EQUAL(fwhm_y,0.0) || EQUAL(fwhm_z,0.0) ) {
      if ( EQUAL(x,mux) && EQUAL(y,muy) && EQUAL(z,muz)  )
	f = c;
      else
	f = 0;
    }
    else {
      
      sigma_x = fwhm_x/2.0;
      sigma_y = fwhm_y/2.0;
      sigma_z = fwhm_z/2.0;
      
      dx = mux-x;
      dy = muy-y;
      dz = muz-z;
      
      
      t1 = c; /*/(sigma_x*sigma_y*sigma_z*fsqrt(two_pi3)); */
      
      sigma_x = sigma_x*sigma_x;  /* square the differences and the sigmas */
      sigma_y = sigma_y*sigma_y;
      sigma_z = sigma_z*sigma_z;
      
      dx = dx*dx;
      dy = dy*dy;
      dz = dz*dz;
      
      t2 = (-1.0/8.0) * 
	(sigma_y*sigma_z*dx + sigma_x*sigma_z*dy + sigma_x*sigma_y*dz) / (sigma_x*sigma_y*sigma_z);
      
      t3 = fexp(t2);
      
      f = t1*t3;
    }
    
  }
  return(f);
}


/* ------------------------------------------------------------
  calculate the scalar correlation value from arrays a1 and a2
  return the value as a float
*/
public float calc_scalar_correlation(float *a1,float *a2, int len)
{
   
   float
      r,
      s1,s2,s3;  /* to store the sums for f1,f2,f3 */

   int
      i;

   s1 = 0.0;
   s2 = 0.0;
   s3 = 0.0;

   for (i=1; i<=len; ++i) {
     s1 += a1[i]*a2[i];
     s2 += a1[i]*a1[i];
     s3 += a2[i]*a2[i];
   }

   if ( s2 < 0.0001 && s3 < 0.0001) {
      r = 1.0;
   }
   else {
      if ( s2 < 0.0001 || s3 < 0.0001) {
	 r = 0.0;
      }
      else {
	 r = s1 / (fsqrt(s2)*fsqrt(s3));
      }
   }

   return(r);
}


private Real cost_fn(float x, float y, float z, Real max_length)
{
  Real v2,v,d;

  v2 = x*x + y*y + z*z;
  v = sqrt(v2);

  v *= v2;

  if (v<max_length)
    d = 0.01 * v / (max_length - v);
  else
    d = FLT_MAX;

  return(d);
}

private float xcorr_fitting_function(float *d)

{
   float
      similarity,cost, r;

/*
   int size[3];
   Real 
     pos[3],voxel;

   get_volume_sizes(Gd2_dxyz, size);


   for_inclusive(c,1,Glen) {
     pos[0] = TX[c]+d[1];
     pos[1] = TY[c]+d[2];
     pos[2] = TZ[c]+d[3];

     convert_3D_world_to_voxel(Gd2_dxyz, 
			       (Real)(TX[c]+d[1]),(Real)(TY[c]+d[2]),(Real)(TZ[c]+d[3]), 
			       &pos[0], &pos[1], &pos[2]);


     if (pos[0]<-0.5 || pos[0]>=size[0]-0.5 ||
	 pos[1]<-0.5 || pos[1]>=size[1]-0.5 ||
	 pos[2]<-0.5 || pos[2]>=size[2]-0.5) {



       Ga2xyz[c] = 0.0;
     }
     else {
       GET_VOXEL_3D(voxel, Gd2_dxyz, (int)(pos[0]+0.5), (int)(pos[1]+0.5), (int)(pos[2]+0.5));
       Ga2xyz[c] = Gd2_dxyz->real_value_scale * voxel + Gd2_dxyz->real_value_translation;
     }
   }
*/


   go_get_samples_with_offset(Gd2_dxyz, TX,TY,TZ,
			      Ga2xyz, d[1], d[2], d[3],
			      Glen, 0);



   cost = (float)cost_fn(d[1], d[2], d[3], Gcost_radius);
   similarity = calc_scalar_correlation(Ga1xyz, Ga2xyz, Glen);

   r = 1.0 - similarity*similarity_cost_ratio + cost*(1.0-similarity_cost_ratio);
   
/*
print ("s c g: %12.8f %12.8f (%6.1f)-- %6.2f %6.2f \n", similarity, cost, Gsimplex_size, x[1], x[2]);
*/

   return(r);
}



