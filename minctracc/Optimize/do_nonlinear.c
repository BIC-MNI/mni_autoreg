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
@MODIFIED   : Revision 1.19  1995-08-21 11:36:43  louis
@MODIFIED   : This version of get_deformation_vector_for_node with the 3x3x3 weights
@MODIFIED   : of 1.0, 1.0, 1.0, 1.0 for all nodes used in the quadratic fit and works
@MODIFIED   : well with return_3D_disp_from_quad_fit() in ver 1.2 of quad_max_fit.c
@MODIFIED   :
@MODIFIED   : This version of the quadratic fit seems to work reasonably well in 3D
@MODIFIED   : and compares favorably to the simplex optimization for hoge/jacob fitting
@MODIFIED   : at 16mm.
@MODIFIED   :
@MODIFIED   : The use of quadratic fitting is twice as fast as simplex (6.4 vs 15.2 min)
@MODIFIED   : for the 16mm fit (again with {hoge,jacob}_16_dxyz.mnc and -step 8 8 8.
@MODIFIED   :
 * Revision 1.18  1995/08/17  12:17:59  louis
 * bug fixed for warping problems at the scalp region.  This was
 * caused by adding a deoformation vector to the field, even when
 * it was not estimated.  In this case, the deformation vector was
 * an average of the neighbourhood warp and caused an amplification
 * of the deformation field where the field should be simply smoothed.
 *
 * Revision 1.17  1995/06/12  14:28:57  louis
 * working version - 2d,3d w/ simplex and -direct.
 *
 * Revision 1.16  1995/05/04  14:25:18  louis
 * compilable version, seems to run a bit with GRID_TRANSFORM, still
 * needs work for the super sampled volumes... and lots of testing.
 *
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/do_nonlinear.c,v 1.19 1995-08-21 11:36:43 louis Exp $";
#endif

#include <volume_io.h>		/* structs & tools to deal with volumes data */
#include "arg_data.h"		/* definition of the global data struct      */
#include <print_error.h>	/* def of print_error_and_..                 */
#include <limits.h>		/* MAXtype and MIN defs                      */
#include <recipes.h>		/* numerical recipes defs                    */
#include "deform_support.h"	/* prototypes for routines called
				   from deformation procedures.              */

#include <sys/types.h>		/* for timing the deformations               */
#include <time.h>
time_t time(time_t *tloc);


         /* these Globals are used to communicate to the correlation */
         /* functions over top the SIMPLEX optimization  routine     */ 
static float	
  Gsqrt_s1,			/* normalization const for correlation       */
  *Ga1xyz,			/* samples in source sub-lattice             */
  *Ga2xyz,			/* samples in target sub-lattice             */
  *TX, *TY, *TZ,		/* sample sub-lattice positions in target    */
  *SX, *SY, *SZ;		/* sample sub-lattice positions in source    */
static int 
  Glen;				/* # of samples in sub-lattice               */
 

         /* these Globals are used to communicate the projection */
         /* values over top the SIMPLEX optimization  routine    */ 
static Real
  Gtarget_vox_x, Gtarget_vox_y, Gtarget_vox_z,
  Gproj_d1,  Gproj_d1x,  Gproj_d1y,  Gproj_d1z, 
  Gproj_d2,  Gproj_d2x,  Gproj_d2y,  Gproj_d2z;

	/* Global pointers to the different data volumes */
static Volume   Gd1;    /* source: blurred volume                            */
static Volume   Gd1_dx;	        /* deriv along x, not always loaded          */
static Volume   Gd1_dy;	        /* deriv along y, not always loaded          */
static Volume   Gd1_dz;	        /* deriv along z, not always loaded          */
static Volume   Gd1_dxyz;       /* Gradient magnitude vol                    */
static Volume   Gm1;		/* mask vol for source                       */
static Volume   Gd2;    /* target: blurred volume                            */
static Volume   Gd2_dx;         /* deriv along x, not always loaded          */
static Volume   Gd2_dy;         /* deriv along y, not always loaded          */
static Volume   Gd2_dz;         /* deriv along z, not always loaded          */
static Volume   Gd2_dxyz;       /* Gradient magnitude vol                    */
static Volume   Gm2;		/* mask vol for target                       */

	/* Globals used for lcoal simplex Optimization  */
static Real     Gsimplex_size;	/* the radius of the local simplex           */
static Real     Gcost_radius;	/* constant used in the cost function        */

        /* Globals used to split the input transformation into a
	   linear part and a super-sampled non-linear part */
static General_transform 
                *Gsuper_sampled_warp,
                *Glinear_transform;
static  Volume  Gsuper_sampled_vol;

	/* Volume order definition for super sampled data */
static char *my_XYZ_dim_names[] = { MIxspace, MIyspace, MIzspace };


        /* program Global data used to store all info regarding data
           and transformations  */
static Arg_Data *Gglobals;


       /* constants defined on command line to control optimization */
extern double     smoothing_weight;      /* weight given to neigbours        */
extern double     iteration_weight;      /* wght given to a singer iteration */
extern double     similarity_cost_ratio; /* obj fn = sim * s+c+r -
					             cost * (1-s_c_r)        */
extern int        iteration_limit;       /* total number of iterations       */
extern int        number_dimensions;     /* ==2 or ==3                       */
extern double     ftol;		         /* stopping tolerence for simplex   */
extern Real       initial_corr;	         /* value of correlation before
					    optimization                     */

				/* absolute maximum range for deformation
				   allowed */
#define ABSOLUTE_MAX_DEFORMATION       50.0

				/* diameter of the local neighbourhood
				   sub-lattice, in number of elements-1 */
#define DIAMETER_OF_LOCAL_LATTICE   7
#define MAX_G_LEN (DIAMETER_OF_LOCAL_LATTICE+1)*\
                  (DIAMETER_OF_LOCAL_LATTICE+1)*\
                  (DIAMETER_OF_LOCAL_LATTICE+1)

				/* weighting for 3D quadratic fitting */
/*
#define CENTER_WEIGHT        1.00
#define CENTER_FACE_WEIGHT   0.99
#define CENTER_EDGE_WEIGHT   0.975
#define CORNER_WEIGHT        0.90

#define CENTER_WEIGHT        1.00
#define CENTER_FACE_WEIGHT   0.99
#define CENTER_EDGE_WEIGHT   0.9875
#define CORNER_WEIGHT        0.95
*/

#define CENTER_WEIGHT        1.00
#define CENTER_FACE_WEIGHT   1.0
#define CENTER_EDGE_WEIGHT   1.0
#define CORNER_WEIGHT        1.0

				/* weighting for 2D quadratic fitting */
#define CENTER_WEIGHT_2D     1.00
#define EDGE_WEIGHT_2D       0.97
#define CORNER_WEIGHT_2D     0.89



        /* prototypes function definitions */


int amoeba2(float **p, 
	    float y[], 
	    int ndim, 
	    float ftol, 
	    float (*funk)(), 
	    int *nfunk);

public float xcorr_objective(Volume d1,
                             Volume d2,
                             Volume m1,
                             Volume m2, 
                             Arg_Data *globals);

public void    build_target_lattice1(float px[], float py[], float pz[],
				    float tx[], float ty[], float tz[],
				    int len);

public void    build_target_lattice2(float px[], float py[], float pz[],
				    float tx[], float ty[], float tz[],
				    int len);


private Real cost_fn(float x, float y, float z, Real max_length);

private Real similarity_fn(float *d);

private float xcorr_fitting_function(float *x);

public Status save_deform_data(Volume dx,
			       Volume dy,
			       Volume dz,
			       char *name,
			       char *history);

private Real get_deformation_vector_for_node(Real spacing, Real threshold1, 
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
					       Real *def_x, Real *def_y, Real *def_z);
  
public BOOLEAN return_3D_disp_from_quad_fit(Real r[3][3][3], /* the values used in the quad fit */
					    Real *dispu, /* the displacements returned */
					    Real *dispv, 
					    Real *dispw);	

public BOOLEAN return_2D_disp_from_quad_fit(Real r[3][3], /* the values used in the quad fit */
					    Real *dispu, /* the displacements returned */
					    Real *dispv);	

public  void  louis_general_transform_point(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed );

public  void  louis_general_inverse_transform_point(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed );




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
  General_transform
    *all_until_last, 
    *additional_warp,
    *current_warp;

  Volume
    estimated_flag_vol,
    additional_vol,
    current_vol,
    additional_mag;
  
  long
    iteration_start_time,
    iteration_end_time,
    timer1,timer2,
    nfunk_total;
  STRING 
    time_total_string;
  Real 
    time_total;

  int 
    additional_count[MAX_DIMENSIONS],
    current_count[MAX_DIMENSIONS],
    mag_count[MAX_DIMENSIONS],
    sizes[MAX_DIMENSIONS],    
    index[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS],
    data_xyzv[MAX_DIMENSIONS],
    start[MAX_DIMENSIONS], 
    end[MAX_DIMENSIONS],
    debug_sizes[MAX_DIMENSIONS],
    iters,
    i,j,k,m,
    nodes_done, nodes_tried, nodes_seen, matching_type, over,
    nfunks,
    nfunk1, nodes1;

  Real 
    debug_steps[MAX_DIMENSIONS],
    voxel[MAX_DIMENSIONS],
    steps[MAX_DIMENSIONS],
    current_steps[MAX_DIMENSIONS],
    steps_data[MAX_DIMENSIONS],
    min, max, sum, sum2, mag, mean_disp_mag, std, var, 
    val_x, val_y, val_z,
    def_vector[3],
    wx,wy,wz,
    node_x, node_y, node_z,
    test_x, test_y, test_z,
    mx,my,mz,
    tx,ty,tz, 
    displace,
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


  Ga1xyz = vector(1,MAX_G_LEN);	/* allocate space for the global data for */
  Ga2xyz = vector(1,MAX_G_LEN);	/* the local neighborhood lattice values */
  
  SX = vector(1,MAX_G_LEN);		/* and coordinates in source volume  */
  SY = vector(1,MAX_G_LEN);
  SZ = vector(1,MAX_G_LEN);
  TX = vector(1,MAX_G_LEN);		/* and coordinates in target volume  */
  TY = vector(1,MAX_G_LEN);
  TZ = vector(1,MAX_G_LEN);

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
  init_the_volume_to_zero(additional_vol);

				/* build a temporary volume that will
				   be used to store the magnitude of
				   the deformation at each iteration */

  additional_mag = create_volume(3, my_XYZ_dim_names, NC_SHORT, FALSE, 0.0, 0.0);
  for_less(i,0,3)
    mag_count[i] = additional_count[ xyzv[i] ];
  set_volume_sizes(additional_mag, mag_count);
  alloc_volume_data(additional_mag);
  set_volume_real_range(additional_mag, 0.0, ABSOLUTE_MAX_DEFORMATION);

				/* reset additional mag to zero */
  init_the_volume_to_zero(additional_mag);
  mean_disp_mag = 0.0;

				/* build a volume to flag all voxels
				   where the deformation has been estimated */
  estimated_flag_vol = create_volume(3, my_XYZ_dim_names, NC_BYTE, FALSE, 0.0, 0.0);
  for_less(i,0,3)
    mag_count[i] = additional_count[ xyzv[i] ];
  set_volume_sizes(estimated_flag_vol, mag_count);
  alloc_volume_data(estimated_flag_vol);
  set_volume_real_range(estimated_flag_vol, 0.0, 255.0);
  init_the_volume_to_zero(estimated_flag_vol);

				/* test simplex size against size
				   of data voxels */

  get_volume_separations(additional_vol, steps);
  get_volume_separations(d2_dxyz, steps_data);

  if (steps_data[0]!=0.0) {
				/* Gsimplex_size is in voxel units
				   in the data volume.             */
    Gsimplex_size= ABS(steps[xyzv[X]]) / MAX3(steps_data[0],steps_data[1],steps_data[2]);  

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
  if (globals->trans_info.use_super>0) {
    ALLOC(Gsuper_sampled_warp,1);
    create_super_sampled_data_volumes(current_warp, 
				      Gsuper_sampled_warp,
				      globals->trans_info.use_super);
    Gsuper_sampled_vol = Gsuper_sampled_warp->displacement_volume;

    if (globals->flags.debug) {
      for_less(i,0,MAX_DIMENSIONS) voxel[i]=0.0;
      convert_voxel_to_world(current_warp->displacement_volume, 
			     voxel,
			     &wx, &wy, &wz);
      get_volume_sizes(current_warp->displacement_volume, 
		       debug_sizes);
      get_volume_separations(current_warp->displacement_volume, 
			     debug_steps);
      print ("After super sampling:\n");
      print ("curnt sizes: %7d  %7d  %7d  %7d  %7d\n",
	     debug_sizes[0],debug_sizes[1],debug_sizes[2],debug_sizes[3],debug_sizes[4]);
      print ("curnt steps: %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n",
	     debug_steps[0],debug_steps[1],debug_steps[2],debug_steps[3],debug_steps[4]);
      print ("curnt start: %7.2f  %7.2f  %7.2f  \n", wx, wy, wz);

      convert_voxel_to_world(Gsuper_sampled_warp->displacement_volume, 
			     voxel,
			     &wx, &wy, &wz);
      get_volume_sizes(Gsuper_sampled_warp->displacement_volume, 
		       debug_sizes);
      get_volume_separations(Gsuper_sampled_warp->displacement_volume, 
			     debug_steps);
      print ("After super sampling:\n");
      print ("super sizes: %7d  %7d  %7d  %7d  %7d\n",
	     debug_sizes[0],debug_sizes[1],debug_sizes[2],debug_sizes[3],debug_sizes[4]);
      print ("super steps: %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n",
	     debug_steps[0],debug_steps[1],debug_steps[2],debug_steps[3],debug_steps[4]);
      print ("super start: %7.2f  %7.2f  %7.2f  \n", wx, wy, wz);


    }
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

  initial_corr = xcorr_objective(d1, d2, m1, m2, globals );

  if (globals->flags.debug) {	
    print("Initial corr         = %f\n",initial_corr);
    print("Source vol threshold = %f\n", threshold1);
    print("Target vol threshold = %f\n", threshold2);
    print("Iteration limit      = %d\n", iteration_limit);
    print("Iteration weight     = %f\n", iteration_weight);
    print("xyzv                 = %3d %3d %3d %3d \n",
	  xyzv[X], xyzv[Y], xyzv[Z], xyzv[Z+1]);
    print("number_dimensions    = %d\n",number_dimensions);
    print("smoothing_weight     = %f\n",smoothing_weight);
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

    if (globals->trans_info.use_super>0)
      interpolate_super_sampled_data(current_warp,
				     Gsuper_sampled_warp,
				     number_dimensions);
  
    
    print("Iteration %2d of %2d\n",iters+1, iteration_limit);

    nodes_done = 0; nodes_tried = 0; nodes_seen=0; 
    displace = 0.0; over = 0;        nfunk_total = 0;

    sum = sum2 = 0.0;
    min = DBL_MAX;
    max = -DBL_MAX;

    iteration_start_time = time(NULL);

    initialize_progress_report( &progress, FALSE, 
			       (end[X]-start[X])*(end[Y]-start[Y]) + 1,
			       "Estimating deformations" );

    for_less(i,0,MAX_DIMENSIONS) index[i]=0;

    for_less( index[ xyzv[X] ] , start[ X ], end[ X ]) {

      timer1 = time(NULL);
      nfunk1 = 0; nodes1 = 0;

      for_less( index[ xyzv[Y] ] , start[ Y ], end[ Y ]) {

	for_less( index[ xyzv[Z] ] , start[ Z ], end[ Z ]) {

	  nodes_seen++;	  
				/* get the lattice coordinate 
				   of the current index node  */
	  for_less(i,0,MAX_DIMENSIONS) voxel[i]=index[i];
	  convert_voxel_to_world(current_vol, 
				 voxel,
				 &node_x, &node_y, &node_z);

				/* get the warp that needs to be added 
				   to get the target world coordinate
				   position */

	  for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
	    def_vector[ index[ xyzv[Z+1] ] ] = 
	      get_volume_real_value(current_vol,
				    index[0],index[1],index[2],index[3],index[4]);

				/* add the warp to get the target 
				   lattice position in world coords */

	  wx = node_x + def_vector[X]; 
	  wy = node_y + def_vector[Y]; 
	  wz = node_z + def_vector[Z];

	  if (point_not_masked(m2, wx, wy, wz) &&
	      get_value_of_point_in_volume(wx,wy,wz, Gd2_dxyz)>threshold2) {

				/* now get the mean warped position of 
				   the target's neighbours */

	    index[ xyzv[Z+1] ] = 0;
	    get_average_warp_of_neighbours(current_warp,
					   index,
					   &mx, &my, &mz);
	    
				/* get the targets homolog in the
				   world coord system of the source
				   data volume                      */

	    general_inverse_transform_point(Glinear_transform,
					    node_x, node_y, node_z,
					    &tx,&ty,&tz);

				/* find the best deformation for
				   this node                        */

	    def_vector[X] = def_vector[Y] = def_vector[Z] = 0.0;


	    result = get_deformation_vector_for_node(steps[xyzv[X]], 
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
	    
	    if (result < 0.0) {
	      nodes_tried++;
	      result = 0.0;
	    } else {
				/* store the deformation vector */
	      for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
		set_volume_real_value(additional_vol,
				      index[0],index[1],index[2],
				      index[3],index[4],
				      def_vector[ index[ xyzv[Z+1] ] ]);
				/* store the def magnitude */
	      set_volume_real_value(additional_mag,
		       index[xyzv[X]],index[xyzv[Y]],index[xyzv[Z]],0,0,
		       result);
				/* set the 'node estimated' flag */
	      set_volume_real_value(estimated_flag_vol,
		       index[xyzv[X]],index[xyzv[Y]],index[xyzv[Z]],0,0,
		       1.0);

				/* calculate statistics  */
	      if (ABS(result) > 0.95*steps[0]) over++;
	      nodes_done++;
	      displace += ABS(result);
	      sum      += ABS(result);
	      sum2     += ABS(result) * ABS(result);
	      
	      if (ABS(result)>max) max = ABS(result);
	      if (ABS(result)<min) min = ABS(result);
	      
	      nfunk_total += nfunks;
	      nfunk1      += nfunks; 
	      nodes1++;	      

	    } /* of else (result<0) */
	    
	    
	  } /* val in volume2 > /* threshold 2 */
	  

	}  /* forless on X index */

	update_progress_report( &progress, 
			       (end[Y]-start[Y])*(index[ xyzv[X]]-start[X])+
			       (index[ xyzv[Y]]-start[Y])+1 );
      } /* forless on Y index */
      timer2 = time(NULL);

      if (globals->flags.debug) 
	print ("xslice: (%3d:%3d) = %d sec -- nodes=%d av funks %f\n",
	       index[ xyzv[X] ] +1-start[X], 
	       end[X]-start[X], 
	       timer2-timer1, 
	       nodes1,
	       nodes1==0? 0.0:(float)nfunk1/(float)nodes1);
      
    } /* forless on X index */

    terminate_progress_report( &progress );

    iteration_end_time = time(NULL);

    if (globals->flags.debug) {

      if (nodes_done>0) {
	mean_disp_mag = displace/nodes_done;
	var = ((sum2 * nodes_done) - sum*sum) / 
	      ((float)nodes_done*(float)(nodes_done-1));
	std = sqrt(var);
	nfunks = nfunk_total / nodes_done;
      }
      else {
	mean_disp_mag=0.0; std = 0.0;
      }
      print ("Nodes seen = %d, tried = %d, done = %d, avg disp = %f +/- %f\n",
	     nodes_seen, nodes_tried, nodes_done, mean_disp_mag, std);
      print ("av nfunks = %d , over = %d, max disp = %f, min disp = %f\n", 
	     nfunks, over, max, min);

      nodes_tried = 0; nodes_seen = 0;
      for_less(i,0,mag_count[0])
	for_less(j,0,mag_count[1])
	  for_less(k,0,mag_count[2]){
	    mag = get_volume_real_value(additional_mag,i,j,k,0,0);
	    if (mag >= 0.000001)
	      nodes_seen++;
	    if (mag >= (mean_disp_mag+std))
	      nodes_tried++;
	  }
      print ("there are %d of %d over (mean+1std) out of %d.\n", 
	     nodes_tried, nodes_done, nodes_seen);

      time_total = (Real)(iteration_end_time-iteration_start_time);
      format_time( time_total_string,"%g %s", time_total);
      print ("This iteration took %s (%d seconds)\n", 
	     time_total_string, 
	     iteration_end_time-iteration_start_time);
    }
				/* update the current warp, so that the
				   next iteration will use all the data
				   calculated thus far.

				   additional = additional + current    */

    add_additional_warp_to_current(additional_warp,
				   current_warp,
				   iteration_weight);


				/* smooth the warp in additional,
				   leaving the result in current 

				   current = smooth(additional) */

    smooth_the_warp(current_warp,
		    additional_warp,
		    additional_mag, -1.0);


    /* clamp the data so that the 1st derivative of the deformation
       field does not exceed 1.0*step in magnitude

    clamp_warp_deriv(current->dx, current->dy, current->dz); */

    
    if (iters<iteration_limit-1) {

				/* reset additional warp to zero */
      init_the_volume_to_zero(additional_vol);

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
				     &node_x, &node_y, &node_z);

				/* get the warp that needs to be added 
				   to get the target world coordinate
				   position */

	      for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
		def_vector[ index[ xyzv[Z+1] ] ] = 
		  get_volume_real_value(current_vol,
					index[0],index[1],index[2],index[3],index[4]);
	
				/* add the warp to get the target 
				   lattice position in world coords */

	      wx = node_x + def_vector[X]; 
	      wy = node_y + def_vector[Y]; 
	      wz = node_z + def_vector[Z];
      
	      if ( point_not_masked(m2, wx, wy, wz)) {
		
		/* now get the mean warped position of the target's neighbours */
		
				/* now get the mean warped position of 
				   the target's neighbours */

		index[ xyzv[Z+1] ] = 0;
		get_average_warp_of_neighbours(current_warp,
					       index,
					       &mx, &my, &mz);

				/* get the targets homolog in the
				   world coord system of the source
				   data volume                      */


		general_inverse_transform_point(Glinear_transform,
						node_x, node_y, node_z,
						&tx,&ty,&tz);

		
				/* find the best deformation for
				   this node                        */

		def_vector[X] = def_vector[Y] = def_vector[Z] = 0.0;

		result = get_deformation_vector_for_node(steps[xyzv[X]],
							 threshold1,
							 tx,ty,tz,
							 mx, my, mz,
							 &def_vector[X],
							 &def_vector[Y],
							 &def_vector[Z],
							 iters, iteration_limit,
							 &nfunks,
							 number_dimensions);
		
    
		if (result>=0) 
		  for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
		    set_volume_real_value(additional_vol,
					  index[0],index[1],index[2],
					  index[3],index[4],
					  def_vector[ index[ xyzv[Z+1] ] ]);
	      
	      } /* point not masked */
	      
	    } /* if mag > mead + std */

	  } /* forless on z index */
	} /* forless on y index */
      } /* forless on x index */
      

				/* additional = additional + current */
      add_additional_warp_to_current(additional_warp,
				     current_warp,
				     iteration_weight);

				/* current = smooth(additional) */
      smooth_the_warp(current_warp,
		      additional_warp,
		      additional_mag, (Real)(mean_disp_mag+std));

    } /* if (iters< iteration_limit-1) */

				/* reset the next iteration's warp. */
    init_the_volume_to_zero(additional_vol);
    init_the_volume_to_zero(additional_mag);

    if (globals->flags.debug && 
	globals->flags.verbose == 3)
      save_data(globals->filenames.output_trans, 
		iters+1, iteration_limit, 
		globals->trans_info.transformation);
    
    if (globals->flags.debug) { 
      print("initial corr %f ->  this step %f\n",initial_corr,xcorr_objective(d1, d2, m1, m2, globals) );
    }
  }
    
    /* free up allocated temporary deformation volumes */

  if (globals->trans_info.use_super>0) {
    delete_general_transform(Gsuper_sampled_warp);
  }

  (void)delete_general_transform(additional_warp);
  (void)delete_volume(additional_mag);

  free_vector(Ga1xyz ,1,MAX_G_LEN);
  free_vector(Ga2xyz ,1,MAX_G_LEN);
  free_vector(TX ,1,MAX_G_LEN);
  free_vector(TY ,1,MAX_G_LEN);
  free_vector(TZ ,1,MAX_G_LEN);
  free_vector(SX ,1,MAX_G_LEN);
  free_vector(SY ,1,MAX_G_LEN);
  free_vector(SZ ,1,MAX_G_LEN);
    

  return (OK);

}




private BOOLEAN get_best_start_from_neighbours(Real threshold1, 
					       Real wx, Real wy, Real wz,
					       Real mx, Real my, Real mz,
					       Real *tx, Real *ty, Real *tz,
					       Real *def_x, Real *def_y, Real *def_z)
     
{
  Real
    mag_normal1,
    nx, ny, nz;



  mag_normal1 = get_value_of_point_in_volume(wx,wy,wz, Gd1_dxyz);

  if (mag_normal1 < threshold1)
    return(FALSE);	
  else {

				/* map point from source, forward into
                                   target space */

    general_transform_point(Gglobals->trans_info.transformation, 
			    wx,wy,wz, tx,ty,tz);

				/* average out target point with the
				   mean position of its neightbours */

    nx = (*tx+mx)/2.0; 
    ny = (*ty+my)/2.0; 
    nz = (*tz+mz)/2.0; 
				/* what is the deformation needed to
				   achieve this displacement */
    *def_x = nx - *tx;
    *def_y = ny - *ty;
    *def_z = nz - *tz;
  
    *tx = nx; *ty = ny; *tz = nz;
  
    return(TRUE);
  }  
}




/***********************************************************/
/* note that the value of the spacing coming in is FWHM/2 for the data
   used to do the correlation. */

private Real get_deformation_vector_for_node(Real spacing, 
					     Real threshold1, 
					     Real src_x, Real src_y, Real src_z,
					     Real mx, Real my, Real mz,
					     Real *def_x, Real *def_y, Real *def_z,
					     int iteration, int total_iters,
					     int *num_functions,
					     int ndim)
{
  Real
    du,dv,dw,
    xt, yt, zt,
    local_corr3D[3][3][3],
    local_corr2D[3][3],
    voxel_displacement[3],
    voxel[3],
    val[MAX_DIMENSIONS],
    pos[3],
    simplex_size,
    result,
    tx,ty,tz, 
    xp,yp,zp;
  float 
    pos_vector[4],
    tmp_val,
     *y, **p;
  int 
    flag,
    nfunk,
    numsteps,
    i,j,k;


  result = 0.0;			/* assume no additional deformation */
  *num_functions = 0;		/* assume no optimization done      */

  
  /* if there is no gradient magnitude strong enough to grab onto,
     then set a negative magnitude deformation and skip the optimization */

  if (!get_best_start_from_neighbours(threshold1, 
				      src_x,  src_y,  src_z,
				      mx,  my,  mz,
				      &tx,  &ty,  &tz,
				      def_x,  def_y,  def_z)) {

    result = -DBL_MAX;		/* set result to a flag value used above */

  }
  else {   
    /* we now have the info needed to continue...

       src_x, src_y, src_z - point in source volume
       tx, ty, tz          - best point in target volume, so far.
       def_x,def_y,def_z   - currently contains the additional def needed to 
                             take src_x,src_y,src_z mid-way to the neighbour's
			     mean-point (which is now stored in tx,ty,tz).   */

    /* get the world coord position of the node in the source volume
       taking into consideration the warp from neighbours */

    if (ndim==3)
      general_inverse_transform_point(Gglobals->trans_info.transformation, 
				      tx,  ty,  tz,
				      &xp, &yp, &zp);
    else
      louis_general_inverse_transform_point(Gglobals->trans_info.transformation, 
					 tx,  ty,  tz,
					 &xp, &yp, &zp);

    if (Gglobals->trans_info.use_magnitude) {

      /* build a spherical sub-lattice of Glen points in the source
         volume, note: sub-lattice diameter= 1.5*fwhm */

      numsteps = DIAMETER_OF_LOCAL_LATTICE;	
      build_source_lattice(xp, yp, zp, 
			   SX, SY, SZ,
			   spacing*3, spacing*3, spacing*3,
			   numsteps+1,  numsteps+1,  numsteps+1,
			   ndim, &Glen);
    }
    else {

      /* store the coordinate of the center of the neighbour only,
         since we are using projections. */

      SX[1] = xp; SY[1] = yp; SZ[1] = zp;
      Glen = 1;
      numsteps = 0;
    }

		/* BUILD THE TARGET INFO */

    if ( Gglobals->trans_info.use_magnitude ) {

      /* map this lattice forward into the target space, using the
	 current transformation, in order to build a deformed lattice
	 (in the WORLD COORDS of the target volume) */

      if (Gglobals->trans_info.use_super>0) 
	build_target_lattice2(SX,SY,SZ, TX,TY,TZ, Glen);
      else 
	build_target_lattice1(SX,SY,SZ, TX,TY,TZ, Glen);

    }
    else {

      /* get the target voxel position */

      if (number_dimensions==3)
	general_transform_point(Gglobals->trans_info.transformation, 
				xp,yp,zp,    &xt, &yt, &zt);
      else
	louis_general_transform_point(Gglobals->trans_info.transformation, 
					xp,yp,zp,  &xt, &yt, &zt);
					
      convert_3D_world_to_voxel(Gd2_dx, xt, yt, zt,
				&Gtarget_vox_x, &Gtarget_vox_y, &Gtarget_vox_z);
    }


    /* for the objective function used in the actual optimization, I
       need the voxel coordinates of the target lattice.  

       note: I assume that the volume is stored in ZYX order!  */

    for_inclusive(i,1,Glen) {
      convert_3D_world_to_voxel(Gd2_dxyz, 
				(Real)TX[i],(Real)TY[i],(Real)TZ[i], 
				&pos[0], &pos[1], &pos[2]);

      if (Gglobals->trans_info.use_magnitude) {

	/* jiggle the voxel coordinates of the sub-lattice off center
	   so that nearest-neighbour interpolation can be used */

	if (ndim>2)		
	  TX[i] = pos[0] - 0.5 + drand48(); 
	else
	  TX[i] = pos[0];
	TY[i] = pos[1] - 0.5 + drand48();
	TZ[i] = pos[2] - 0.5 + drand48();
      }
      else {

	/* otherwise, store only the voxel position of the center of
           target lattice, used below to get the projections */

	TX[i] = pos[0];
	TY[i] = pos[1];
	TZ[i] = pos[2];
      }
    }

    /* re-build the source lattice (without local neighbour warp),
       that will be used in the optimization below */

    if (Gglobals->trans_info.use_magnitude) {
      for_inclusive(i,1,Glen) {
	SX[i] += src_x - xp;
	SY[i] += src_y - yp;
	SZ[i] += src_z - zp;
      }
    }
    
    if (Gglobals->trans_info.use_magnitude) {
      go_get_samples_in_source(Gd1_dxyz, SX,SY,SZ, Ga1xyz, Glen, 1);
    }
    else {			

      /* build projection data, ie go get f, df/dx, df/dy, df/dz to be
         able to reconstruct the equivalent of a sub-lattice, based
         only on these first four terms of the Taylor expansion. */

      evaluate_volume_in_world(Gd1,
			       xp, yp, zp, 
			       0, TRUE, 0.0, val,
			       NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d1 = val[0];
      evaluate_volume_in_world(Gd1_dx,
			       xp, yp, zp, 
			       0, TRUE, 0.0, val,
			       NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d1x = 2.0*val[0];
      evaluate_volume_in_world(Gd1_dy,
			       xp, yp, zp, 
			       0, TRUE, 0.0, val,
			       NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d1y = 2.0*val[0];
      evaluate_volume_in_world(Gd1_dz,
			       xp, yp, zp, 
			       0, TRUE, 0.0, val,
			       NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d1z = 2.0*val[0];

    }

    /* calc one of the normalization coefficients for correlation,
       when using the magnitude data.  Note that this variable is not
       used when using projections. */

    Gsqrt_s1 = 0.0;
    for_inclusive(i,1,Glen)
      Gsqrt_s1 += Ga1xyz[i]*Ga1xyz[i];

    Gsqrt_s1 = sqrt((double)Gsqrt_s1);

    /* now do the optimization to find the best local deformation that
       maximises the local neighbourhood correlation between the source
       values stored in Ga1xyz and the values at positions TX,TY,TZ 
       in the target volume */

    if ( !Gglobals->trans_info.use_simplex) {
				/* do quad fit optimization */

      if (ndim==3) {
	/* build up the 3x3x3 matrix of local correlation values */

	for_inclusive(i,-1,1)
	  for_inclusive(j,-1,1)
	    for_inclusive(k,-1,1) {
	      pos_vector[1] = (float) i * Gsimplex_size/2.0;
	      pos_vector[2] = (float) j * Gsimplex_size/2.0;
	      pos_vector[3] = (float) k * Gsimplex_size/2.0;
	      local_corr3D[i+1][j+1][k+1] = similarity_fn(pos_vector); 
	    }
	*num_functions = 27;


/*
      print ("with proj data:\n");
      print ("%8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f\n",
	     local_corr3D[0][0][0],local_corr3D[1][0][0],local_corr3D[2][0][0], 
	     local_corr3D[0][0][1],local_corr3D[1][0][1],local_corr3D[2][0][1],
	     local_corr3D[0][0][2],local_corr3D[1][0][2],local_corr3D[2][0][2]);
      print ("%8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f\n",
	     local_corr3D[0][1][0],local_corr3D[1][1][0],local_corr3D[2][1][0], 
	     local_corr3D[0][1][1],local_corr3D[1][1][1],local_corr3D[2][1][1],
	     local_corr3D[0][1][2],local_corr3D[1][1][2],local_corr3D[2][1][2]);
      print ("%8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f\n\n\n",
	     local_corr3D[0][2][0],local_corr3D[1][2][0],local_corr3D[2][2][0], 
	     local_corr3D[0][2][1],local_corr3D[1][2][1],local_corr3D[2][2][1],
	     local_corr3D[0][2][2],local_corr3D[1][2][2],local_corr3D[2][2][2]);
*/


				/*  the central node  */
	local_corr3D[1][1][1] *= CENTER_WEIGHT;

				/*  the six face centers  */

	local_corr3D[1][1][0] *= CENTER_FACE_WEIGHT;
	local_corr3D[1][1][2] *= CENTER_FACE_WEIGHT;
	local_corr3D[1][0][1] *= CENTER_FACE_WEIGHT;
	local_corr3D[1][2][1] *= CENTER_FACE_WEIGHT;
	local_corr3D[0][1][1] *= CENTER_FACE_WEIGHT;
	local_corr3D[2][1][1] *= CENTER_FACE_WEIGHT;

				/*  the twelve edge centers  */
	local_corr3D[1][0][0] *= CENTER_EDGE_WEIGHT;
	local_corr3D[1][0][2] *= CENTER_EDGE_WEIGHT;
	local_corr3D[1][2][0] *= CENTER_EDGE_WEIGHT;
	local_corr3D[1][2][2] *= CENTER_EDGE_WEIGHT;

	local_corr3D[0][1][0] *= CENTER_EDGE_WEIGHT;
	local_corr3D[0][1][2] *= CENTER_EDGE_WEIGHT;
	local_corr3D[2][1][0] *= CENTER_EDGE_WEIGHT;
	local_corr3D[2][1][2] *= CENTER_EDGE_WEIGHT;

	local_corr3D[0][0][1] *= CENTER_EDGE_WEIGHT;
	local_corr3D[0][2][1] *= CENTER_EDGE_WEIGHT;
	local_corr3D[2][0][1] *= CENTER_EDGE_WEIGHT;
	local_corr3D[2][2][1] *= CENTER_EDGE_WEIGHT;

				/*  the eight corners  */
	local_corr3D[0][0][0] *= CORNER_WEIGHT;
	local_corr3D[0][0][2] *= CORNER_WEIGHT;
	local_corr3D[0][2][0] *= CORNER_WEIGHT;
	local_corr3D[0][2][2] *= CORNER_WEIGHT;
	local_corr3D[2][0][0] *= CORNER_WEIGHT;
	local_corr3D[2][0][2] *= CORNER_WEIGHT;
	local_corr3D[2][2][0] *= CORNER_WEIGHT;
	local_corr3D[2][2][2] *= CORNER_WEIGHT;


	/* use a quadratic fit to estimate the local maxima */
	
	flag = return_3D_disp_from_quad_fit(local_corr3D, &du, &dv, &dw);
	
      }
      else {
	/* build up the 3x3 matrix of local correlation values */
	
	for_inclusive(i,-1,1)
	  for_inclusive(j,-1,1) {
	    pos_vector[1] = (float) i * Gsimplex_size/2.0;
	    pos_vector[2] = (float) j * Gsimplex_size/2.0;
	    pos_vector[3] = 0.0;
	    local_corr2D[i+1][j+1] = similarity_fn(pos_vector); 
	  }
	*num_functions = 9;
	
	local_corr2D[1][1] *= CENTER_WEIGHT_2D;

	local_corr2D[0][1] *= EDGE_WEIGHT_2D;
	local_corr2D[2][1] *= EDGE_WEIGHT_2D;
	local_corr2D[1][0] *= EDGE_WEIGHT_2D;
	local_corr2D[1][2] *= EDGE_WEIGHT_2D;

	local_corr2D[0][0] *= CORNER_WEIGHT_2D;
	local_corr2D[0][2] *= CORNER_WEIGHT_2D;
	local_corr2D[2][0] *= CORNER_WEIGHT_2D;
	local_corr2D[2][2] *= CORNER_WEIGHT_2D;
	
	/* use a quadratic fit to estimate the local maxima */


	flag = return_2D_disp_from_quad_fit(local_corr2D,  &du, &dv);
	dw = 0.0;
	
      }

      if ( flag ) {
	voxel_displacement[0] = dw * Gsimplex_size/2.0;
	voxel_displacement[1] = dv * Gsimplex_size/2.0;
	voxel_displacement[2] = du * Gsimplex_size/2.0;
      }
      else {
	result = -DBL_MAX;
	voxel_displacement[0] = 0.0;
	voxel_displacement[1] = 0.0;
	voxel_displacement[2] = 0.0;
      }


    }
    else {	       /* set up SIMPLEX OPTIMIZATION */

      nfunk = 0;
      p = matrix(1,ndim+1,1,3);	/* simplex structure */
      y = vector(1,ndim+1);	/* value of correlation at simplex vertices */
    
				/* init simplex for no deformation */
      p[1][1] = p[1][2] = p[1][3] = 0.0;
      for (i=2; i<=(ndim+1); ++i)	
	for (j=1; j<=3; ++j)
	  p[i][j] = p[1][j];
				/* set the simplex diameter so as to 
				   reduce the size of the simplex, and
				   hence reduce the search space with
				   each iteration.                      */
      simplex_size = Gsimplex_size * 
	(0.5 + 
	 0.5*((Real)(total_iters-iteration)/(Real)total_iters));
      
      p[2][1] += simplex_size;	/* set up all vertices of simplex */
      p[3][2] += simplex_size;
      if (ndim > 2) {
	p[4][3] += simplex_size;
      }
				/* set up value of correlation at 
				   each vertex of simplex */
      for (i=1; i<=(ndim+1); ++i)	{ 
	y[i] = xcorr_fitting_function(p[i]);
      }

				/* do the actual SIMPLEX optimization */
      if (amoeba2(p,y,ndim,ftol,xcorr_fitting_function,&nfunk)) {    
	
				/* find the best result */
	if ( y[1] < y[2] + 0.00001 )
	  i=1;
	else
	  i=2;
	
	if ( y[i] > y[3] + 0.00001)
	  i=3;
	
	if ((ndim > 2) && ( y[i] > y[4] + 0.00001))
	  i=4;
	
	*num_functions = nfunk;
    
	if (ndim>2) {
	  voxel_displacement[0] = p[i][3];
	  voxel_displacement[1] = p[i][2];
	  voxel_displacement[2] = p[i][1];
	}
	else {
	  voxel_displacement[0] = 0.0;
	  voxel_displacement[1] = p[i][2];
	  voxel_displacement[2] = p[i][1];
	}

      }
      else {
	      /* simplez optimization found nothing, so set the additional
		 displacement to 0 */

	voxel_displacement[0] = 0.0;
	voxel_displacement[1] = 0.0;
	voxel_displacement[2] = 0.0;
	result                = 0.0;
	*num_functions        = 0;

      } /*  if ameoba2 */

      free_matrix(p, 1,ndim+1,1,3);    
      free_vector(y, 1,ndim+1);  

    } /* else use_magnitude */


    convert_3D_world_to_voxel(Gd2_dxyz, 
			      tx,ty,tz, 
			      &voxel[0], &voxel[1], &voxel[2]);
    
    convert_3D_voxel_to_world(Gd2_dxyz, 
			      (Real)(voxel[0]+voxel_displacement[0]), 
			      (Real)(voxel[1]+voxel_displacement[1]), 
			      (Real)(voxel[2]+voxel_displacement[2]),
			      &pos[X], &pos[Y], &pos[Z]);
    
    *def_x += pos[X]-tx;
    *def_y += pos[Y]-ty;
    *def_z += pos[Z]-tz;
    
    result = sqrt((*def_x * *def_x) + (*def_y * *def_y) + (*def_z * *def_z)) ;      
    if (result > 50.0) {
      print ("??? displacement too big: %f\n",result);
      print ("vox_disp %f %f %f\n",
	     voxel_displacement[0],
	     voxel_displacement[1],
	     voxel_displacement[2]);
      print ("deform   %f %f %f\n", def_x, def_y, def_z);
      print ("Gsimplex_size = %f\n",Gsimplex_size);
    }


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


    if (number_dimensions==3)
      general_transform_point(Gglobals->trans_info.transformation, 
			      (Real)px[i],(Real) py[i], (Real)pz[i], 
			      &x, &y, &z);
    else
      louis_general_transform_point(Gglobals->trans_info.transformation, 
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
  int 
    i,j,
    n_dim,
    sizes[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS];
  Real 
    def_vector[N_DIMENSIONS],
    voxel[MAX_DIMENSIONS],
    x,y,z;
  Real vx,vy,vz;
  Real dx,dy,dz;
  long 
    index[MAX_DIMENSIONS],
    ind0, ind1, ind2;


  get_volume_sizes(Gsuper_sampled_vol,sizes);
  get_volume_XYZV_indices(Gsuper_sampled_vol,xyzv);
  n_dim = get_volume_n_dimensions(Gsuper_sampled_vol);

  for_inclusive(i,1,len) {

    general_transform_point(Glinear_transform,
			    (Real)px[i], (Real)py[i], (Real)pz[i], 
			    &x, &y, &z);
    convert_world_to_voxel(Gsuper_sampled_vol, 
			   x,y,z, voxel);


    if ((voxel[ xyzv[X] ] >= -0.5) && (voxel[ xyzv[X] ] < sizes[xyzv[X]]-0.5) &&
	(voxel[ xyzv[Y] ] >= -0.5) && (voxel[ xyzv[Y] ] < sizes[xyzv[Y]]-0.5) &&
	(voxel[ xyzv[Z] ] >= -0.5) && (voxel[ xyzv[Z] ] < sizes[xyzv[Z]]-0.5) ) {

      for_less(j,0,3) index[ xyzv[j] ] = (long) (voxel[ xyzv[j] ]+0.5);
      
      for_less(index[ xyzv[Z+1] ],0, sizes[ xyzv[Z+1] ]) 
	GET_VALUE_4D(def_vector[ index[ xyzv[Z+1] ]  ], \
		     Gsuper_sampled_vol, \
		     index[0], index[1], index[2], index[3]);


      x += def_vector[X];
      y += def_vector[Y];
      z += def_vector[Z];
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

private Real similarity_fn(float *d)
{
  Real 
    val[MAX_DIMENSIONS],
    voxel[MAX_DIMENSIONS],
    xw,yw,zw,
    s, s1, s2;
				/* note: here the displacement order
				   is 3,2,1 since the same function is
				   used for 2D and 3D optimization.
				   In 2D, simplex will only modify
				   d[1] and d[2].
				   (d[3] stays const=0) */

  if (Gglobals->trans_info.use_magnitude) {
    s = (Real)go_get_samples_with_offset(Gd2_dxyz, TX,TY,TZ,
					 d[3], d[2], d[1],
					 Glen, Gsqrt_s1, Ga1xyz);
  } else {
	/* calc correlation based on projection data */
    voxel[0] = Gtarget_vox_x + d[3];
    voxel[1] = Gtarget_vox_y + d[2]; 
    voxel[2] = Gtarget_vox_z + d[1];
    voxel[3] = 0.0;
    voxel[4] = 0.0;

    convert_voxel_to_world(Gd2, voxel, &xw, &yw, &zw);

    evaluate_volume_in_world(Gd2, xw, yw, zw, 
			     0, TRUE, 0.0, val,
			     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    Gproj_d2 = val[0];
    evaluate_volume_in_world(Gd2_dx, xw, yw, zw, 
			     0, TRUE, 0.0, val,
			     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    Gproj_d2x = 2.0*val[0];
    evaluate_volume_in_world(Gd2_dy, xw, yw, zw, 
			     0, TRUE, 0.0, val,
			     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    Gproj_d2y = 2.0*val[0];
    evaluate_volume_in_world(Gd2_dz, xw, yw, zw, 
			     0, TRUE, 0.0, val,
			     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    Gproj_d2z = 2.0*val[0];
    
    s  = (Gproj_d1*Gproj_d2   + Gproj_d1x*Gproj_d2x + 
	  Gproj_d1y*Gproj_d2y + Gproj_d1z*Gproj_d2z);
    s1 = (Gproj_d1*Gproj_d1   +  Gproj_d1x*Gproj_d1x + 
	  Gproj_d1y*Gproj_d1y + Gproj_d1z*Gproj_d1z);
    s2 = (Gproj_d2*Gproj_d2   + Gproj_d2x*Gproj_d2x + 
	  Gproj_d2y*Gproj_d2y + Gproj_d2z*Gproj_d2z);

    if ((s1!=0.0) && (s2 !=0.0)) {
      s = s / ( sqrt(s1) * sqrt(s2) );
    }
    else
      s = 0.0;
    
    
  }
  return( s );
}



private float xcorr_fitting_function(float *d)

{
  float
    similarity,
    cost, 
    r;

   similarity = (float)similarity_fn( d );
   cost =       (float)cost_fn( d[1], d[2], d[3], Gcost_radius );

   r = 1.0 - 
       similarity * similarity_cost_ratio + 
       cost       * (1.0-similarity_cost_ratio);
   
   return(r);
}

