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
@MODIFIED   : Revision 9.6  1996-08-21 18:22:02  louis
@MODIFIED   : Pre-release
@MODIFIED   :
 * Revision 9.5  1996/08/12  14:16:03  louis
 * Never released version 0.95
 *
 * Revision 1.26  1996/08/12  14:15:54  louis
 * Pre-release
 *
 * Revision 1.25  1996/03/25  10:33:15  collins
 * removed jiggle in do_nonlinear.c - The jiggle was causing multiple
 * local minima in the local objective function.  This in turn was the
 * MAJOR  CAUSE of errors in the estimation of the deformation field when
 * using simplex optimization (and in fact was caused error for quad fit
 * as well).
 *
 * ----
 * go_get_samples_in_source called with -1 for inter_type (before
 * the value for inter_type was ignored).
 *
 * ----
 * completed non-isotropic smoothing in return_locally_smoothed_def()
 * that is based on the eigen vectors of the hessian matrix (that is
 * derived from the differential geometry of the shape of the local
 * objective function).  The estimated deformation vectors are modulated
 * by a confidence function (see next para) such that the vectors are
 * smoothed less in the direction where the obj function has greater
 * curvature.
 *
 * ----
 * confidence_function has been implemented to give a confidence of
 *      -> 1.0 to eigen values >= the mean largest eigen_val (eval'd over
 *         the previous iteration)
 *      -> 0.0 to eigen values <= the mean smallest eigen_val
 *      -> [0.0 .. 0.5] to eigen_vals between smallest and middle mean
 *         eigen val
 *      -> [0.5 .. 1.0] to eigen_vals between middle and largest mean
 *         eigen val
 * This confidence_function is used to weight the amount of the additional
 * deformation vector added to the current def in return_locally_smoothed_def().
 *
 * Revision 1.24  1996/03/07  13:25:19  collins
 * small reorganisation of procedures and working version of non-isotropic
 * smoothing.
 *
 * Revision 1.23  1995/10/13  11:16:55  collins
 * added comments for local fitting process.
 *
 * changed the objective function used in the quadratic fitting routines to
 * be the same as that used for simplex, both should (and now do use)
 * xcorr_fitting_function() so that both the similarity of the data and the
 * cost of the deformation are taken into account.
 *
 * started to add code for non-isotropic smoothing, but it is not completely
 * working in this version.
 *
 * Revision 1.22  1995/10/06  09:25:02  collins
 * working version.  Added multiple feature support.  The extra features
 * can included in the global similarity function with user-selectable
 * individual similarity functions (straight diff, xcorr, label diff).
 *
 * I removed all references to the following global variables (since they
 * are no longer needed with the feature structure: (Gsqrt_s1,*Ga1xyz,
 * *Ga2xyz,Gd1,Gd1_dx,Gd1_dy,Gd1_dz,Gd1_dxyz,Gm1,Gd2,Gd2_dx,Gd2_dy,
 * Gd2_dz,Gd2_dxyz,Gm2)
 *
 * I removed all the volume parameters from do_non_linear_optimization()
 * since they are now all represented in the globals->features struc.
 *
 * changed the code wherever Gsqrt_s1, Ga1xyz, Ga2xyz, Gd1_dxyz, etc were
 * used, replacing the code with the proper reference to globals->features.
 *
 * the similarity function is replaced by the weighted sum of the individual
 * feature similarity functions.
 *
 * this version has been tested and seems to work with multiple features in
 * 2D and in 3D.
 *
 * Revision 1.21  1995/09/14  09:53:49  collins
 * working version, using David's simplex optimization replacing the
 * numerical recipes amoeba() routine.
 *
 * Revision 1.20  1995/09/07  10:05:11  collins
 * All references to numerical recipes routines are being removed.  At this
 * stage, any num rec routine should be local in the file.  All memory
 * allocation calls to vector(), matrix(), free_vector() etc... have been
 * replaced with ALLOC and FREE from the volume_io library of routines.
 *
 * Revision 1.19  1995/08/21  11:36:43  collins
 * This version of get_deformation_vector_for_node with the 3x3x3 weights
 * of 1.0, 1.0, 1.0, 1.0 for all nodes used in the quadratic fit and works
 * well with return_3D_disp_from_quad_fit() in ver 1.2 of quad_max_fit.c
 *
 * This version of the quadratic fit seems to work reasonably well in 3D
 * and compares favorably to the simplex optimization for hoge/jacob fitting
 * at 16mm.
 *
 * The use of quadratic fitting is twice as fast as simplex (6.4 vs 15.2 min)
 * for the 16mm fit (again with {hoge,jacob}_16_dxyz.mnc and -step 8 8 8.
 *
 * Revision 1.18  1995/08/17  12:17:59  collins
 * bug fixed for warping problems at the scalp region.  This was
 * caused by adding a deoformation vector to the field, even when
 * it was not estimated.  In this case, the deformation vector was
 * an average of the neighbourhood warp and caused an amplification
 * of the deformation field where the field should be simply smoothed.
 *
 * Revision 1.17  1995/06/12  14:28:57  collins
 * working version - 2d,3d w/ simplex and -direct.
 *
 * Revision 1.16  1995/05/04  14:25:18  collins
 * compilable version, seems to run a bit with GRID_TRANSFORM, still
 * needs work for the super sampled volumes... and lots of testing.
 *
 * Revision 1.14  1995/02/22  08:56:06  collins
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/do_nonlinear.c,v 9.6 1996-08-21 18:22:02 louis Exp $";
#endif

#include <config.h>		/* MAXtype and MIN defs                      */
#include <volume_io.h>		/* structs & tools to deal with volumes data */
#include <amoeba.h>		/* simplex optimization struct               */

#include <stdlib.h>		/* to get header info for drand48()          */
#include <arg_data.h>		/* definition of the global data struct      */
#include <print_error.h>	/* def of print_error_and_..                 */
#include <deform_support.h>	/* prototypes for routines called
				   from deformation procedures.              */

#include "constants.h"		/* internal constant definitions             */

#include <sys/types.h>		/* for timing the deformations               */
#include <time.h>
time_t time(time_t *tloc);

#include <stats.h>


int FLAG_HAHA = FALSE;

int stat_quad_total;
int stat_quad_zero;
int stat_quad_two;
int stat_quad_plus;
int stat_quad_minus;
int stat_quad_semi;

         /* these Globals are used to communicate to the correlation */
         /* functions over top the SIMPLEX optimization  routine     */ 
static float	
  *Gsqrt_features,		/* normalization const for correlation       */
  **Ga1_features,		/* samples in source sub-lattice             */
  *TX, *TY, *TZ,		/* sample sub-lattice positions in target    */
  *SX, *SY, *SZ;		/* sample sub-lattice positions in source    */
static int 
  Glen;				/* # of samples in sub-lattice               */
 
static stats_struct
   stat_def_mag,
   stat_num_funks,
   stat_conf0,
   stat_conf1,
   stat_conf2,
   stat_eigval0,
   stat_eigval1,
   stat_eigval2;


         /* these Globals are used to communicate the projection */
         /* values over top the SIMPLEX optimization  routine    */ 
static Real
  Gtarget_vox_x, Gtarget_vox_y, Gtarget_vox_z,
  Gproj_d1,  Gproj_d1x,  Gproj_d1y,  Gproj_d1z, 
  Gproj_d2,  Gproj_d2x,  Gproj_d2y,  Gproj_d2z;

	/* Globals used for local simplex Optimization  */
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
extern double     smoothing_weight;      /* weight given to neighbours       */
extern double     iteration_weight;      /* wght given to a singer iteration */
extern double     similarity_cost_ratio; /* obj fn = sim * s+c+r -
					             cost * (1-s_c_r)        */
extern int        iteration_limit;       /* total number of iterations       */
extern int        number_dimensions;     /* ==2 or ==3                       */
extern double     ftol;		         /* stopping tolerence for simplex   */
extern Real       initial_corr, final_corr;
				         /* value of correlation before/after
					    optimization                     */

				/* diameter of the local neighbourhood
				   sub-lattice, in number of elements-1 */

extern int        Diameter_of_local_lattice;
#define MAX_G_LEN (Diameter_of_local_lattice)*\
                  (Diameter_of_local_lattice)*\
                  (Diameter_of_local_lattice)


				/* absolute maximum range for deformation
				   allowed */
#define ABSOLUTE_MAX_DEFORMATION       50.0

				/* constants for non-isotropic smoothing, used
				   for first iteration only.  They are updated
				   for each following iteration.              */

/*
#define DEFAULT_MEAN_E0  0.048457
#define DEFAULT_MEAN_E1  0.008863
#define DEFAULT_MEAN_E2  0.004401
#define DEFAULT_STD_E0   0.024583
#define DEFAULT_STD_E1   0.006246
#define DEFAULT_STD_E2   0.002482
*/

#define DEFAULT_MEAN_E0  0.023836
#define DEFAULT_MEAN_E1  0.005587
#define DEFAULT_MEAN_E2  0.000871
#define DEFAULT_STD_E0   0.018712
#define DEFAULT_STD_E1   0.010252
#define DEFAULT_STD_E2   0.002711

static Real
   previous_mean_eig_val[3] = {DEFAULT_MEAN_E0,DEFAULT_MEAN_E1,DEFAULT_MEAN_E2},
   previous_std_eig_val[3]  = {DEFAULT_STD_E0,DEFAULT_STD_E1,DEFAULT_STD_E2};

        /* prototypes function definitions */


public  BOOLEAN  perform_amoeba(amoeba_struct  *amoeba, int *num_funks );
public  void  initialize_amoeba(amoeba_struct     *amoeba,
				int               n_parameters,
				Real              initial_parameters[],
				Real              parameter_delta,
				amoeba_function   function,
				void              *function_data,
				Real              tolerance );

public  Real  get_amoeba_parameters(amoeba_struct  *amoeba,
				    Real           parameters[] );

public  void  terminate_amoeba( amoeba_struct  *amoeba );

private Real amoeba_obj_function(void * dummy, float d[]);

public float xcorr_objective(Volume d1,
                             Volume d2,
                             Volume m1,
                             Volume m2, 
                             Arg_Data *globals);

public void    build_target_lattice(float px[], float py[], float pz[],
				    float tx[], float ty[], float tz[],
				    int len);

public void    build_target_lattice_using_super_sampled_def(
                                    float px[], float py[], float pz[],
				    float tx[], float ty[], float tz[],
				    int len);

private Real cost_fn(float x, float y, float z, Real max_length);

private Real similarity_fn(float *d);

private Real xcorr_fitting_function(float *x);

private Real get_deformation_vector_for_node(Real spacing, Real threshold1, 
					     Real source_coord[],
					     Real mean_target[],
					     Real def_vector[],
					     Real voxel_displacement[],
					     int iteration, int total_iters,
					     int *nfunks,
					     int ndim);

private double return_locally_smoothed_def(int  isotropic_smoothing,
					 int  ndim,
					 Real smoothing_wght,
					 Real iteration_wght,
					 Real smoothed_result[],
					 Real previous_def[],
					 Real neighbour_mean[],
					 Real additional_def[],
					 Real another_vector[],
					 Real voxel_displacement[]);

public BOOLEAN return_local_eigen(Real r[3][3][3],
				  Real dir_1[3],
				  Real dir_2[3],
				  Real dir_3[3],
				  Real val[3]);

public BOOLEAN return_local_eigen_from_hessian(Real r[3][3][3],
				  Real dir_1[3],
				  Real dir_2[3],
				  Real dir_3[3],
				  Real val[3]);

private BOOLEAN get_best_start_from_neighbours(
					       Real threshold1, 
					       Real source[],
					       Real mean_target[],
					       Real target[],
					       Real def[]);

  
public BOOLEAN return_3D_disp_from_quad_fit(Real r[3][3][3], /* the values used in the quad fit */
					    Real *dispu, /* the displacements returned */
					    Real *dispv, 
					    Real *dispw);	

public BOOLEAN return_3D_disp_from_min_quad_fit(Real r[3][3][3], /* the values used in the quad fit */
					    Real *dispu, /* the displacements returned */
					    Real *dispv, 
					    Real *dispw);	

public BOOLEAN return_2D_disp_from_quad_fit(Real r[3][3], /* the values used in the quad fit */
					    Real *dispu, /* the displacements returned */
					    Real *dispv);	

public  void  general_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed );

public  void  general_inverse_transform_point_in_trans_plane(
    General_transform   *transform,
    Real                x,
    Real                y,
    Real                z,
    Real                *x_transformed,
    Real                *y_transformed,
    Real                *z_transformed );

public void get_volume_XYZV_indices(Volume data, int xyzv[]);

public int point_not_masked(Volume volume, 
			    Real wx, Real wy, Real wz);



public int nearest_neighbour_interpolant(Volume volume, 
                                 PointR *coord, double *result);



/**************************************************************************/
public Status do_non_linear_optimization(Arg_Data *globals)
{
  General_transform
    *all_until_last,		/* will contain the first (linear) part of the xform  */
    *additional_warp,		/* storage of estimates of the needed additional warp */
    *another_warp,		/* storage of estimates of the needed additional warp */
    *current_warp;		/* storage of the current (best) warp so far          */

  Volume
    estimated_flag_vol,		/* flags indicating node estimated or not             */
    additional_vol,		/* volume pointer to additional_warp transform        */
    another_vol,		/* volume pointer to additional_warp transform        */
    current_vol,		/* volume pointer to current_warp transform           */
    additional_mag;		/* volume storing mag of additional_warp vectors      */
  
  long
    iteration_start_time,	/* variables to time each iteration                   */
    iteration_end_time,
    temp_start_time,
    temp_end_time,
    timer1,timer2,
    nfunk_total;
  STRING 
    time_total_string;
  Real 
    time_total;

  int 
    additional_count[MAX_DIMENSIONS], /* size (in voxels) of  additional_vol          */
    mag_count[MAX_DIMENSIONS],	      /* size (in voxels) of  additional_mag          */
    xyzv[MAX_DIMENSIONS],	      /* order of voxel indices                       */
    index[MAX_DIMENSIONS],	      /* used to step through all nodes of def field  */
    start[MAX_DIMENSIONS],	      /* starting limit of index[]                    */
    end[MAX_DIMENSIONS],	      /* ending limit of index[]                      */
    debug_sizes[MAX_DIMENSIONS],
    iters,			      /* iteration counter */
    i,j,k,
    nodes_done, nodes_tried,	 /* variables to calc stats on deformation estim  */
    nodes_seen, over,
    nfunks, nfunk1, nodes1;

  Real 
    debug_steps[MAX_DIMENSIONS], 
    voxel[MAX_DIMENSIONS],	/* voxel position used for interpolation         */
    mag_steps[MAX_DIMENSIONS],
    steps[MAX_DIMENSIONS],	/* voxel size (width,height,length) of def field */
    steps_data[MAX_DIMENSIONS],	/* voxel size of data                            */

				/* variables to calc stats on deformation estim  */
    mag, mean_disp_mag, std, 

    another_vector[3],	
    current_def_vector[3],	/* the current deformation vector for a  node    */
    result_def_vector[3],	/* the smoothed deformation vector for a node    */
    def_vector[3],		/* the additional deformation estimated for node,
				   for these three vectors in real world coords
				   index [0] stores the x disp, and [2] the z disp */
    voxel_displacement[3],      /* the additional displacement, in voxel coords 
				   with [0] storing the displacement of the
				   fastest varying index, and [2] the slowest in
				   the data voluming = xdim and zdim respectively*/
    wx,wy,wz,			/* temporary storage for a world coordinate      */
    source_node[3],		/* world coordinate of source node               */
    target_node[3],		/* world coordinate of corresponding target node */
    mean_target[3],		/* mean deformed pos, determined by neighbors    */
    mean_vector[3],		/* mean deformed vector, determined by neighbors */
    threshold1,			/* intensity thresh for source vol               */
    threshold2,			/* intensity thresh for target vol               */
    result,			/* magnitude of additional warp vector           */
    eig1;

  progress_struct		/* to print out program progress report */
    progress;
  STRING filenamestring;

  /*******************************************************************************/

		    /* set up globals for communication with other routines */


  Gglobals= globals;


  if (Gglobals->features.number_of_features > 0) {
    ALLOC2D (Ga1_features, Gglobals->features.number_of_features, MAX_G_LEN+1);
    ALLOC( Gsqrt_features, Gglobals->features.number_of_features);
  }
  else {
    print_error_and_line_num("There are no features to match in non_lin optimization\n", 
		__FILE__, __LINE__);
  }
  ALLOC(SX,MAX_G_LEN+1);		/* and coordinates in source volume  */
  ALLOC(SY,MAX_G_LEN+1);
  ALLOC(SZ,MAX_G_LEN+1);
  ALLOC(TX,MAX_G_LEN+1);		/* and coordinates in target volume  */
  ALLOC(TY,MAX_G_LEN+1);
  ALLOC(TZ,MAX_G_LEN+1);

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


  ALLOC(another_warp,1);	/* build a debugging volume/xform */
  copy_general_transform(current_warp, another_warp);
  another_vol = another_warp->displacement_volume;
  set_volume_real_range(another_vol, -20.0 , 20.0);
  init_the_volume_to_zero(another_vol);


				/* reset the additional warp to zero */
  init_the_volume_to_zero(additional_vol);


				/* build a temporary volume that will
				   be used to store the magnitude of
				   the deformation at each iteration */


  additional_mag = create_volume(3, my_XYZ_dim_names, NC_SHORT, FALSE, 0.0, 0.0);
  for_less(i,0,3)
    mag_count[i] = additional_count[ xyzv[i] ];
  set_volume_sizes(additional_mag, mag_count);
  get_volume_separations(additional_vol, steps);
  for_less(i,0,3)
    mag_steps[i] = steps[ xyzv[i] ];
  set_volume_separations(additional_mag, mag_steps);
  for_less(i,0,MAX_DIMENSIONS)
    voxel[i] = 0.0;
  convert_voxel_to_world(additional_vol, 
			 voxel,
			 &(target_node[X]), &(target_node[Y]), &(target_node[Z]));
  set_volume_translation(additional_mag, voxel, target_node);
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
  set_volume_separations(estimated_flag_vol, mag_steps);
  set_volume_translation(estimated_flag_vol, voxel, target_node);
  alloc_volume_data(estimated_flag_vol);
  set_volume_real_range(estimated_flag_vol, 0.0, 1.0);
 
	 			/* test simplex size against size
		 		   of data voxels */

  get_volume_separations(additional_vol, steps);
  get_volume_separations(Gglobals->features.model[0], steps_data);

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

				/* build a super-sampled version of the
				   current transformation, if needed     */
  if (globals->trans_info.use_super>0) {
    ALLOC(Gsuper_sampled_warp,1);

    if (globals->trans_info.use_super==2) {
      create_super_sampled_data_volumes_by2(current_warp, 
					    Gsuper_sampled_warp);
    }
    else
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
      print ("Before super sampling, orig def field:\n");
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

  set_feature_value_threshold(Gglobals->features.data[0], 
			      Gglobals->features.model[0],
			      &(globals->threshold[0]), 
			      &(globals->threshold[1]),
			      &threshold1,
			      &threshold2);			      

  if (threshold1<0.0 || threshold2<0.0) {
    print_error_and_line_num("Gradient magnitude threshold error: %f %f\n",
			     __FILE__, __LINE__, 
			     threshold1, threshold2);
  }

  initial_corr = xcorr_objective(Gglobals->features.data[0], Gglobals->features.model[0],
				 Gglobals->features.data_mask[0], Gglobals->features.model_mask[0],
				 globals );

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


    if ( Gglobals->trans_info.use_simplex) 
      print ("This fit will use local simplex optimization and\n");
    else
      print ("This fit will use local quadratic fitting and\n");

    if (Gglobals->trans_info.use_local_smoothing) {
      if ( Gglobals->trans_info.use_local_isotropic) 
	print ("  local isotroptic smoothing.\n");
      else
	print ("  local non-isotroptic smoothing.\n");
    }
    else {
      print ("  global smoothing.\n");
    }
    if (Gglobals->interpolant==nearest_neighbour_interpolant) 
      print ("The similarity function will be evaluated using NN interpolation\n");
    else
      print ("The similarity function will be evaluated using tri-linear interpolation\n");
    print ("  on a spherical sub-lattice with a diameter of %d nodes.\n",
	   Diameter_of_local_lattice);

    print("loop: (%d %d) (%d %d) (%d %d)\n",
	  start[0],end[0],start[1],end[1],start[2],end[2]);

  }


  /*************************************************************************/
  /*    
        start the iterations to estimate the warp                     

	for each iteration {

	   supersample the deformation if needed

	   for each voxel in the deformation field defined on the target volume {
	      get the original coordinate node in the source volume
	      estimate the best deformation for the node, taking into 
	        consideration the previous deformation
	      store this best deformation
	   }

	   for each voxel in the deformation field that had a large deformation {
	      get the original coordinate node in the source volume
	      esitmate the best deformation for this node, taking into 
	        consideration the previous deformation field
	      update the deformation for this node
	   }

	   for each voxel in the deformation field {
 	      add WEIGHTING*best_deformation to the current deformation
	   }

	   smooth the deformation field

	}
  */

  for_less(iters,0,iteration_limit) {

    if (globals->trans_info.use_super>0) {

      temp_start_time = time(NULL);

      if (globals->trans_info.use_super==2) {
	interpolate_super_sampled_data_by2(current_warp,
					       Gsuper_sampled_warp,
					       number_dimensions);
      }
      else
	interpolate_super_sampled_data(current_warp,
				     Gsuper_sampled_warp,
				     number_dimensions);
      temp_end_time = time(NULL);

      if (globals->flags.debug) {	
	time_total = (Real)(temp_end_time-temp_start_time);

#ifdef HAVE_RECENT_VOLUME_IO
	time_total_string = format_time("%g %s", time_total);
#else
	format_time( time_total_string,"%g %s", time_total);
#endif

	print ("\nInterpolating super-sampled data %s (%d seconds)\n", 
	       time_total_string, 
	       temp_end_time-temp_start_time);
      }
    }  

    init_the_volume_to_zero(estimated_flag_vol);
    
    print("Iteration %2d of %2d\n",iters+1, iteration_limit);


				/* for various stats on this iteration*/
    stat_quad_total = 0;
    stat_quad_zero = 0;
    stat_quad_two  = 0;
    stat_quad_plus = 0;
    stat_quad_minus= 0;
    stat_quad_semi = 0;

    nodes_done  = 0; 
    nodes_tried = 0; 
    nodes_seen  = 0; 
    over        = 0;        
    nfunk_total = 0;
    std = 0.0;

    init_stats(&stat_def_mag,  "def_mag");
    init_stats(&stat_num_funks,"num_funks");
    init_stats(&stat_eigval0,  "eigval[0]");
    init_stats(&stat_eigval1,  "eigval[1]");
    init_stats(&stat_eigval2,  "eigval[2]");
    init_stats(&stat_conf0,    "conf[0]");
    init_stats(&stat_conf1,    "conf[1]");
    init_stats(&stat_conf2,    "conf[2]");


    iteration_start_time = time(NULL);

    initialize_progress_report( &progress, FALSE, 
			       (end[X]-start[X])*(end[Y]-start[Y]) + 1,
			       "Estimating deformations" );

    for_less(i,0,MAX_DIMENSIONS) index[i]=0;

    /* FLAG_HAHA stuff:   

       start[X] = 7;  end[X] = start[X]+1;
       start[Y] = 13; end[Y] = start[Y]+1;
       start[Z] = 15; end[Z] = start[Z]+1;
       */

    /* step index[] through all the nodes in the deformation field. */

    for_less( index[ xyzv[X] ] , start[ X ], end[ X ]) {

      timer1 = time(NULL);	/* for stats on this slice */
      nfunk1 = 0; nodes1 = 0;

      for_less( index[ xyzv[Y] ] , start[ Y ], end[ Y ]) {

	for_less( index[ xyzv[Z] ] , start[ Z ], end[ Z ]) {

	  nodes_seen++;	  
				/* get the lattice coordinate 
				   of the current index node  */
	  for_less(i,0,MAX_DIMENSIONS) voxel[i]=index[i];

	  convert_voxel_to_world(current_vol, 
				 voxel,
				 &(target_node[X]), &(target_node[Y]), &(target_node[Z]));

				/* get the warp that needs to be added 
				   to get the target world coordinate
				   position */

	  for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
	    current_def_vector[ index[ xyzv[Z+1] ] ] = 
	      get_volume_real_value(current_vol,
				    index[0],index[1],index[2],index[3],index[4]);

				/* add the warp to get the target 
				   lattice position in world coords */

	  wx = target_node[X] + current_def_vector[X]; 
	  wy = target_node[Y] + current_def_vector[Y]; 
	  wz = target_node[Z] + current_def_vector[Z];

	  if (point_not_masked(Gglobals->features.model_mask[0], wx, wy, wz) &&
	      get_value_of_point_in_volume(wx,wy,wz, Gglobals->features.model[0])>threshold2) {


				/* now get the mean warped position of 
				   the target's neighbours */

	    index[ xyzv[Z+1] ] = 0;
	    if (get_average_warp_of_neighbours(current_warp,
					   index, mean_target)) {
	    
				/* what is the offset to the mean_target? */

	      for_inclusive(i,X,Z)
		mean_vector[i] = mean_target[i] - target_node[i];

				/* get the targets homolog in the
				   world coord system of the source
				   data volume                      */

	      general_inverse_transform_point(Glinear_transform,
			target_node[X], target_node[Y], target_node[Z],
			&(source_node[X]),&(source_node[Y]),&(source_node[Z]));

				/* find the best deformation for
				   this node                        */


/*	      FLAG_HAHA = (index[xyzv[X]] == 7 && index[xyzv[Y]] == 13 && index[xyzv[Z]] == 15);
*/

	      result = get_deformation_vector_for_node(steps[xyzv[X]], 
						       threshold1,
						       source_node,
						       mean_target,
						       def_vector,
						       voxel_displacement,
						       iters, iteration_limit, 
						       &nfunks,
						       number_dimensions);
	    
	      if (result < 0.0) {
		nodes_tried++;
		result = 0.0;
	      } else {
				/* store the deformation vector */

		eig1 = 0.0;
		if (Gglobals->trans_info.use_local_smoothing) {

		  eig1 = return_locally_smoothed_def(
				Gglobals->trans_info.use_local_isotropic,
				number_dimensions,
				smoothing_weight,
				iteration_weight,
				result_def_vector,
				current_def_vector,
				mean_vector,
				def_vector,
				another_vector,
				voxel_displacement);

		  /* Remember that I can't modify current_vol just
		     yet, so I have to set additional_vol to a value,
		     that when added to current_vol (below) I will
		     have the correct result!  */

		  for_inclusive(i,X,Z)
		    result_def_vector[ i ] -= current_def_vector[ i ];

		  for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ])  {
		    set_volume_real_value(additional_vol,
					  index[0],index[1],index[2],
					  index[3],index[4],
					  result_def_vector[ index[ xyzv[Z+1] ] ]);
		    set_volume_real_value(another_vol,
					  index[0],index[1],index[2],
					  index[3],index[4],
					  another_vector[ index[ xyzv[Z+1] ] ]);
		  }

		}
		else {		/* then prepare for global smoothing, (this will
				   actually be done after all nodes 
				   have been estimated  */

		  for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
		    set_volume_real_value(additional_vol,
					  index[0],index[1],index[2],
					  index[3],index[4],
					  def_vector[ index[ xyzv[Z+1] ] ]);

		}
				/* store the def magnitude */

		set_volume_real_value(additional_mag,
		       index[xyzv[X]],index[xyzv[Y]],index[xyzv[Z]],0,0,
		       result);
				/* set the 'node estimated' flag */
		set_volume_real_value(estimated_flag_vol,
		       index[xyzv[X]],index[xyzv[Y]],index[xyzv[Z]],0,0,
		       1.0);

				/* tally up some statistics for this iteration */
		if (ABS(result) > 0.95*steps[xyzv[X]]) over++;
		
		nfunk_total += nfunks;
		nfunk1      += nfunks; 
		nodes1++;	
		nodes_done++;

		tally_stats(&stat_def_mag,   result);
		tally_stats(&stat_num_funks, nfunks);

	      } /* of else (result<0) */

	    } /* if get_average_warp_of_neighbours */
	    
	  } /* point not masked and val in volume2 > threshold 2 */
	  
	    
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


    if (globals->flags.debug) {

/*
      sprintf (filenamestring,"pr_axes2_%d.xfm",iters);
      (void)output_transform_file(filenamestring,NULL,another_warp);
*/

      stat_title();
      report_stats(&stat_num_funks);
      report_stats(&stat_def_mag);
      if (Gglobals->trans_info.use_local_smoothing && 
	  !Gglobals->trans_info.use_local_isotropic) {
	report_stats(&stat_eigval0);
	report_stats(&stat_eigval1);
	report_stats(&stat_eigval2);
	report_stats(&stat_conf0);
	report_stats(&stat_conf1);
	report_stats(&stat_conf2);

      }
      if (!Gglobals->trans_info.use_simplex) {
	print ("quad fit stats: tot + ~ 0 2 -: %5d %5d %5d %5d %5d %5d\n",
	       stat_quad_total,
	       stat_quad_plus,
	       stat_quad_semi,
	       stat_quad_zero,
	       stat_quad_two,
	       stat_quad_minus);
      }

      print ("Nodes seen = %d, tried = %d, done = %d, over = %d\n",
	     nodes_seen, nodes_tried, nodes_done, over);

      mean_disp_mag = stat_get_mean(&stat_def_mag);
      std           = stat_get_standard_deviation(&stat_def_mag);
      
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

      print ("there are %d of %d over (mean+1std) out of %d estimated.\n", 
	     nodes_tried, nodes_done, nodes_seen);
      
      
    }


    if (Gglobals->trans_info.use_local_smoothing && 
	!Gglobals->trans_info.use_local_isotropic) {
      previous_mean_eig_val[0] = stat_get_mean(&stat_eigval0);
      previous_mean_eig_val[1] = stat_get_mean(&stat_eigval1);
      previous_mean_eig_val[2] = stat_get_mean(&stat_eigval2);
      previous_std_eig_val[0]  = stat_get_standard_deviation(&stat_eigval0);
      previous_std_eig_val[1]  = stat_get_standard_deviation(&stat_eigval1);
      previous_std_eig_val[2]  = stat_get_standard_deviation(&stat_eigval2);
    }


				/* update the current warp, so that the
				   next iteration will use all the data
				   calculated thus far.   	*/

    if (Gglobals->trans_info.use_local_smoothing) {

				/* extrapolate (and smooth) the warp to 
				   un-estimated nodes in additional_vol */
      
      extrapolate_to_unestimated_nodes(current_warp,
				       additional_warp,
				       estimated_flag_vol);

				/* current = current + additional */

      add_additional_warp_to_current(current_warp,
				     additional_warp,
				     1.0);



    }
    else {
				/* additional = additional + current    */

      temp_start_time = time(NULL);
      add_additional_warp_to_current(additional_warp,
				     current_warp,
				     iteration_weight);
      temp_end_time = time(NULL);

      if (globals->flags.debug) {	
	time_total = (Real)(temp_end_time-temp_start_time);
#ifdef HAVE_RECENT_VOLUME_IO
	time_total_string = format_time("%g %s", time_total);
#else
	format_time( time_total_string,"%g %s", time_total);
#endif
	print ("\nAdding additional to current took  %s (%d seconds)\n", 
	       time_total_string, 
	       temp_end_time-temp_start_time);
      }
       
      
				/* smooth the warp in additional,
				   leaving the result in current 

				   current = smooth(additional) */

      temp_start_time = time(NULL);

      smooth_the_warp(current_warp,
		      additional_warp,
		      additional_mag, -1.0);
      
      temp_end_time = time(NULL);

      if (globals->flags.debug) {	
	time_total = (Real)(temp_end_time-temp_start_time);
#ifdef HAVE_RECENT_VOLUME_IO
	time_total_string = format_time("%g %s", time_total);
#else
	format_time( time_total_string,"%g %s", time_total);
#endif
	print ("\nSmoothing the current warp took  %s (%d seconds)\n", 
	       time_total_string, 
	       temp_end_time-temp_start_time);
      }
    }

    /* clamp the data so that the 1st derivative of the deformation
       field does not exceed 1.0*step in magnitude

    clamp_warp_deriv(current->dx, current->dy, current->dz); */

    
    if (/*FALSE && */iters<iteration_limit-1) {

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
				     &(target_node[X]), &(target_node[Y]), &(target_node[Z]));

				/* get the warp that needs to be added 
				   to get the target world coordinate
				   position */

	      for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
		current_def_vector[ index[ xyzv[Z+1] ] ] = 
		  get_volume_real_value(current_vol,
					index[0],index[1],index[2],index[3],index[4]);
	
				/* add the warp to get the target 
				   lattice position in world coords */

	      wx = target_node[X] + current_def_vector[X]; 
	      wy = target_node[Y] + current_def_vector[Y]; 
	      wz = target_node[Z] + current_def_vector[Z];
      
	      if ( point_not_masked(Gglobals->features.model_mask[0], wx, wy, wz)) {
		
				/* now get the mean warped position of 
				   the target's neighbours */

		index[ xyzv[Z+1] ] = 0;

		if (get_average_warp_of_neighbours(current_warp,
					       index, mean_target)) {

				/* what is the offset to the mean_target? */

		  for_inclusive(i,X,Z)
		    mean_vector[i] = mean_target[i] - target_node[i];

				/* get the targets homolog in the
				   world coord system of the source
				   data volume                      */

		  general_inverse_transform_point(Glinear_transform,
			     target_node[X], target_node[Y], target_node[Z],
			     &(source_node[X]),&(source_node[Y]),&(source_node[Z]));
		
				/* find the best deformation for
				   this node                        */


		  result = get_deformation_vector_for_node(steps[xyzv[X]],
							   threshold1,
							   source_node,
							   mean_target,
							   def_vector,
							   voxel_displacement,
							   iters, iteration_limit, 
							   &nfunks,
							   number_dimensions);

		  if (result>=0) {

		    if (Gglobals->trans_info.use_local_smoothing) {

		      eig1 = return_locally_smoothed_def(
				Gglobals->trans_info.use_local_isotropic,
				number_dimensions,
				smoothing_weight,
				iteration_weight,
				result_def_vector,
				current_def_vector,
				mean_vector,
				def_vector,
				another_vector,
				voxel_displacement);

		      /* 
			 Remember that I can't modify current_vol just yet, so I
			 have to set additional_vol to a value, that when added to
			 current_vol (below) I will have the correct result!
  		      */

		      for_inclusive(i,X,Z) {
			result_def_vector[ i ] -= current_def_vector[ i ];

/*			current_def_vector[i] =((1.0 - smoothing_weight)*
		                                (current_def_vector[i] + def_vector[i])) + 
					       (smoothing_weight * mean_vector[i])       -
					       current_def_vector[i];
*/
		      }
		      for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
			set_volume_real_value(additional_vol,
					      index[0],index[1],index[2],
					      index[3],index[4],
					      result_def_vector[ index[ xyzv[Z+1] ] ]);

		    }
		    else {
		      for_less( index[ xyzv[Z+1] ], start[ Z+1 ], end[ Z+1 ]) 
			set_volume_real_value(additional_vol,
					      index[0],index[1],index[2],
					      index[3],index[4],
					      def_vector[ index[ xyzv[Z+1] ] ]);
		    }
		  } /* if (result>=0) */
		} /* if  get_average_warp_of_neighbours */
	      } /* if point not masked */
	      
	    } /* if mag > mead + std */

	  } /* forless on z index */
	} /* forless on y index */
      } /* forless on x index */
      

      if (Gglobals->trans_info.use_local_smoothing) {
				/* current = current + additional */

	add_additional_warp_to_current(current_warp,
				       additional_warp,
				       1.0);

      }
      else {
				/* additional = additional + current */

	add_additional_warp_to_current(additional_warp,
				       current_warp,
				       iteration_weight);
	
				/* current = smooth(additional) */
	smooth_the_warp(current_warp,
			additional_warp,
			additional_mag, (Real)(mean_disp_mag+std));
      }
    } /* if (iters< iteration_limit-1) */

				/* reset the next iteration's warp. */


    init_the_volume_to_zero(additional_vol);
    init_the_volume_to_zero(additional_mag);

    if (globals->flags.debug && 
	globals->flags.verbose == 3 && FALSE)
      save_data(globals->filenames.output_trans, 
		iters+1, iteration_limit, 
		globals->trans_info.transformation);
    
    if (globals->flags.debug) {

      iteration_end_time = time(NULL);
      time_total = (Real)(iteration_end_time-iteration_start_time);
#ifdef HAVE_RECENT_VOLUME_IO
	time_total_string = format_time("%g %s", time_total);
#else
	format_time( time_total_string,"%g %s", time_total);
#endif
      print ("\nThis iteration took %s (%d seconds)\n", 
	     time_total_string, 
	     iteration_end_time-iteration_start_time);

      final_corr = xcorr_objective(Gglobals->features.data[0], Gglobals->features.model[0],
				   Gglobals->features.data_mask[0], Gglobals->features.model_mask[0],
				   globals );
      print("initial corr %f ->  this step %f\n",
	    initial_corr,final_corr);

    }
  }
    
  if (!globals->flags.debug)
    final_corr = xcorr_objective(Gglobals->features.data[0], Gglobals->features.model[0],
				 Gglobals->features.data_mask[0], Gglobals->features.model_mask[0],
				 globals );


    /* free up allocated temporary deformation volumes */

  if (globals->trans_info.use_super>0) {
    delete_general_transform(Gsuper_sampled_warp);
  }

  (void)delete_general_transform(additional_warp);
  (void)delete_volume(additional_mag);

  if (Gglobals->features.number_of_features>0) {
    FREE2D(Ga1_features);
    FREE(Gsqrt_features);
  }
  FREE(TX );
  FREE(TY );
  FREE(TZ );
  FREE(SX );
  FREE(SY );
  FREE(SZ );
    

  return (OK);

}

/*
   if local isotropic smoothing:
          n+1    ___n                                 ___n
       def    =  def  + (1-sw)(def + iw*additional -  def )
          where sw is smoothing_weight,
                iw iteration_wieght
   otherwise:
      use directionally sensitive smoothing, by fitting a quadratic
      function to the local neighbourhood of the objective function
      to determine the maximum and minimum principal curvature directions
      in order to smooth only in these directions, and not normal to the
      obj func isosurf.

          n+1    ___n                              ___n    ^      ^
       def    =  def  + (c_max/(1+c_max))( (def -  def ) . e_max) e_max
                                                   ___n    ^      ^
                      + (c_min/(1+c_min))( (def -  def ) . e_min) e_min

     where c_max, c_min are 'confidence measures':
                      ^
                    | k1 |
      c_max =  ------------------; and c_min likewise
               a + b*Smin +c*|k1|

     with a,b and c determined experimentally to have mean( c_max )
     around 1.0.

     note that this is only possible when the principal directions are
     defined.  When ????
*/
#define G_MEAN  0.05
#define K_MEAN  5.0
#define A_const K_MEAN
#define B_const 0.0
#define C_const 0.0

#define P1_MEAN 0.33
#define P2_MEAN 0.33
#define P3_MEAN 0.33

#define CONST_P 0.5
#define CONST_P_MIN .30
#define CONST_P_WIDTH 0.07

/*  #define CONF_FN( x ) ( ((x) < 2.0*CONST_P) ? \
      ( (1.0 + cos( (0.5*PI * (x)/CONST_P)))/2.0 ) : 0.0 )

      c_max = ABS(k1) / (A_const + B_const*Smin  +C_const*ABS(k1));
      c_min = ABS(k2) / (A_const + B_const*Smin  +C_const*ABS(k2));
      c_norm= gradient_magnitude / G_MEAN;
      
      c_norm= ABS(gradient_magnitude) / (P1_MEAN + B_const*Smin  +C_const*ABS(gradient_magnitude));
      c_max = ABS(k1) / (P2_MEAN + B_const*Smin  +C_const*ABS(k1));
      c_min = ABS(k2) / (P3_MEAN + B_const*Smin  +C_const*ABS(k2));
      
	 */
      


/* return a value representing confidence, with the value between 0 and 1 */

private double confidence_function(double x) {

  
  double t;


  t = 0.5;
				/* double linear */
  if 
    (x > previous_mean_eig_val[0]) t = 1.0;
  else { 
    if  

      (x < previous_mean_eig_val[2]) t = 0.0;

    else {

      if (x > previous_mean_eig_val[1]) /* first linear part */

	t = 0.5 + 0.5 * (x -  previous_mean_eig_val[1]) / 
	  (  previous_mean_eig_val[0] - previous_mean_eig_val[1]);

      else			      /* second linear part */

	t = 0.5 * (x -  previous_mean_eig_val[2]) / 
	  (  previous_mean_eig_val[1] - previous_mean_eig_val[2]);

    }
  }

  return(t);

}

private double return_locally_smoothed_def(int isotropic_smoothing,
					   int  ndim,
					   Real smoothing_wght,
					   Real iteration_wght,
					   Real smoothed_result[],
					   Real previous_def[],
					   Real neighbour_mean[],
					   Real additional_def[],
					   Real another_vector[],
					   Real voxel_displacement[])
{

  Real
    mean_vector[3],
    local_corr3D[3][3][3],
    mag_eig_vals,
    conf[3],			/* confidence values                    */
    eig_vals[3],		/* eigen values                         */
    eig_vecs[3][3],		/* eigen vectors (STORED IN ROWS HERE!) */
    diff[3],			/* to represent def - mean_def          */
    len[3],			/* length of projection onto eig_vecs   */
    Smin,			/* best local correlation value         */
    eps;			/* epsilon value                        */
  int
    flag,i,j,k;

  float 
    pos_vector[4];

  eps = 0.0001; /* SMALL_EPSILON_VALUE*/
  eig_vals[0] = 0.0;

				/* let the mean be the average of the 
				   previous_def and the previous neighbour_mean */
  for_inclusive(i,X,Z) 
    mean_vector[i] =  (previous_def[i] + neighbour_mean[i]) / 2.0;

				/* calc the diff vec = estimated - mean */
  for_inclusive(i,X,Z) {
    diff[i] = previous_def[i] + iteration_wght*additional_def[i] - 
              mean_vector[i];
  }

  mag_eig_vals = 0.0;

  if (isotropic_smoothing) {

    /* use the user-defined (command-line option -stiff)
       smoothing-weight to smooth the estimated deformation vector.
       Higher stiffness means more smoothing: 
         smoothing_weight = 1.0
             -> no estimates ever returned, since they are eliminated in the
	        (1.0 -smoothing_weight) term.  
         smoothing_weight = 0.5
	     -> the vector returned is the average of the estimated
                deformation vector and the neighbourhood mean def
         smoothing_weight = 0.0
	     -> no smoothing at all.  the vector returned is simply the
	        estimated deformation vector. */

    for_inclusive(i,X,Z) 
      smoothed_result[i] = mean_vector[i] + (1.0 - smoothing_wght)*diff[i]; 

  }
  else {

    /* use a quadratic approximation to the objective function around
       the local neighbourhood to determine the smoothing directions.
       This will use the global information already stored for
       evaluation of the similarity function, with the exception of
       one change: TX, TY, TZ store the voxel coordinates of the
       target lattice centered on the previous target node (remember
       that TX holds the slowest varying data index).  These values
       have to be updated with the value in voxel_displacement in
       order to approximate the obj function of the source node with
       the _new_ target position */

		      /* take into account the weighting for this iteration */
    for_inclusive(i,X,Z) {
      voxel_displacement[i] *= iteration_wght;
    }
		      /* update target lattice position */
    for_less(i,1,Glen) {
      TX[i] += voxel_displacement[2]; /* slowest varying index for data */
      TY[i] += voxel_displacement[1];
      TZ[i] += voxel_displacement[0]; /* fastest index */
    }

    flag = FALSE;
    if (ndim==2) {
      /* build up the 3x3 matrix of local correlation values,
	 and get the directions */
      
      print_error_and_line_num("2D non-isotropic smoothing not yet supported\n", 
			       __FILE__, __LINE__);
    } else {

      /* build up the 3x3x3 matrix of local correlation values,
         and get the principal directions */

      Smin = DBL_MAX;
      for_inclusive(i,-1,1)
	for_inclusive(j,-1,1)
	  for_inclusive(k,-1,1) {
	    pos_vector[1] = (float) (i * Gsimplex_size)/4.0;
	    pos_vector[2] = (float) (j * Gsimplex_size)/4.0;
	    pos_vector[3] = (float) (k * Gsimplex_size)/4.0;

	    local_corr3D[i+1][j+1][k+1] = xcorr_fitting_function(pos_vector); 

	    if ( local_corr3D[i+1][j+1][k+1] < Smin)
	      Smin = local_corr3D[i+1][j+1][k+1];

	  }

      flag = return_local_eigen_from_hessian(local_corr3D, 
					     eig_vecs[0], eig_vecs[1], eig_vecs[2], eig_vals);

    }
    
    if ( flag ) {		/* if eigen vectors found, then use them */
				/* get associated confidence values */


      for_less(i,0,3)
	conf[i] = confidence_function( eig_vals[i] );

      conf[0] = 0.2;
      conf[1] = 0.5;
      conf[2] = 0.5;


      tally_stats(&stat_eigval0, eig_vals[0]);
      tally_stats(&stat_eigval1, eig_vals[1]);
      tally_stats(&stat_eigval2, eig_vals[2]);
      tally_stats(&stat_conf0, conf[0]);
      tally_stats(&stat_conf1, conf[1]);
      tally_stats(&stat_conf2, conf[2]);
		
				/* project the diff onto each of the 
				   eigen vecs [i] */
      for_less(i,0,3) {
	len[i] = (diff[X]*eig_vecs[i][X] + 
		  diff[Y]*eig_vecs[i][Y] + 
		  diff[Z]*eig_vecs[i][Z]);
      }

				/* calculate the non-iso smoothed
                                   result, using the confidence values
                                   as a weighting factor: 
				   see eq. above for def(n+1) */
      for_inclusive(i,X,Z) {
	smoothed_result[i] = mean_vector[i] + 
	                     conf[2] * len[2] * eig_vecs[2][i] +
			     conf[1] * len[1] * eig_vecs[1][i] +
			     conf[0] * len[0] * eig_vecs[0][i];
      }

				/* for debugging volume... */

      another_vector[X] = 200.0*eig_vals[0]*eig_vecs[0][X];
      another_vector[Y] = 200.0*eig_vals[0]*eig_vecs[0][Y];
      another_vector[Z] = 200.0*eig_vals[0]*eig_vecs[0][Z];
/*
      another_vector[X] = 100.0*eig_vals[0];
      another_vector[Y] = 100.0*eig_vals[1];
      another_vector[Z] = 100.0*eig_vals[2];
*/

    }
    else {			/* default to local isotropic smoothing */
      for_inclusive(i,X,Z) 
	smoothed_result[i] = mean_vector[i]  +(1.0 - smoothing_wght)*diff[i]; 

      for_inclusive(i,X,Z)	/* debugging volume */
	another_vector[i] = 0.0;
    }

    
  } /* else nonisotropic smoothing */

  return(mag_eig_vals);
}


/*
   get_best_start_from_neighbours: find the 'best' point in the target
   volume to start searching for the true 'best' position.

   returns 
      FALSE if the source coordinate does not fall within the
            threshold range, otherwise
      TRUE  and
            the starting target = (current_target + mean_target) / 2
	        def = deformation needed to bring current_target to 
		      the best starting target.
*/
private BOOLEAN get_best_start_from_neighbours(
			   Real threshold1, 
			   Real source[],
			   Real mean_target[],
			   Real target[],
			   Real def[])
     
{
  Real
    mag_normal1,
    nx, ny, nz;


  mag_normal1 = get_value_of_point_in_volume(source[X],source[Y],source[Z], 
					     Gglobals->features.data[0]);

  if (mag_normal1 < threshold1)
    return(FALSE);	
  else {
				/* map point from source, forward into
                                   target_ space */

    general_transform_point(Gglobals->trans_info.transformation, 
			    source[X],source[Y],source[Z], 
			    &(target[X]),&(target[Y]),&(target[Z]));

				/* average out target_ point with the
				   mean position of its neightbours */
    nx = (target[X] + mean_target[X])/2.0; 
    ny = (target[Y] + mean_target[Y])/2.0; 
    nz = (target[Z] + mean_target[Z])/2.0; 
				/* what is the deformation needed to
				   achieve this displacement */
    def[X] = nx - target[X];
    def[Y] = ny - target[Y];
    def[Z] = nz - target[Z];


				/* reset the target location to be
				   halfway to the mean position of its
				   neighbours */

    target[X] = nx; target[Y] = ny; target[Z] = nz;

    return(TRUE);
  }  
}




/**********************************************************

  get_deformation_vector_for_node will return the magnitude of the
  deformation, if one can be estimated by either simplex or quadratic
  fitting.  Otherwise, the procedure will return -DBL_MAX.

  The procedure begins by calling calling build_source_lattice() to
  define the local neighbourhood with a sub-lattice defined on the
  source volume.  The coordinates of the source lattice are mapped
  into the target space for eventual use when the objective function
  must be evaluated for different possible offsets.

  Either "Nelder Mead Simplex" or "Quadratic Obj Func Fitting" will be
  used to determine the best deformation vector (additional offset)
  that maximises local neighboughood correlation between source and
  target volumes.

  If Quadratic: then the objective function is evaluated on a 3x3x3
  neighbourhood and a quadratic function is fit to the data. The
  minimum of the 3D obj function is found directly in
  return_3D_disp_from_min_quad_fit().

  If Simplex: the objective function is evaluated on the 4 vertices of
  a 3D simplex. perform_amoeba() is called repreatedly until the
  volume of the ameoba has been reduced below a pre-selected
  tolerence.

  The necessary additional offset is returned in def_vector[].

  note that the value of the spacing coming in is FWHM/2 for the data
  used to do the correlation.


*/

private Real get_deformation_vector_for_node(Real spacing, 
					     Real threshold1, 
					     Real source_coord[],
					     Real mean_target[],
					     Real def_vector[],
					     Real voxel_displacement[],
					     int iteration, int total_iters,
					     int *num_functions,
					     int ndim)
{
  Real
    du,dv,dw,
    xt, yt, zt,
    local_corr3D[3][3][3],
    local_corr2D[3][3],
    voxel[3],
    val[MAX_DIMENSIONS],
    pos[3],
    simplex_size,
    result,
    target_coord[3],
    xp,yp,zp;
  float 
    pos_vector[4];
  int 
    flag,
    nfunk,
    i,j,k;
  amoeba_struct
    the_amoeba;
  Real
    *parameters;
  FILE *tmp_fp;



  result = 0.0;			/* assume no additional deformation */
  def_vector[X] = def_vector[Y] = def_vector[Z] = 0.0;
  *num_functions = 0;		/* assume no optimization done      */

  
  /* if there is no gradient magnitude strong enough to grab onto,
     then set a negative magnitude deformation and skip the optimization 

     otherwise

     bias the initial starting search position for the local
     deformation vector by calculating a target position that is equal
     to the vector average of 
     (current_transformation(source) + average_position(source in target)

     */


  if (!get_best_start_from_neighbours(Gglobals->features.thresh_data[0], 
				      source_coord, mean_target, target_coord,
				      def_vector)) {

    result = -DBL_MAX;		/* set result to a flag value used
				   above when looping through all nodes */

  }
  else {   
    /* we now have the info needed to continue...

       source_coord[] - point in source volume
       target_coord[] - best point in target volume, so far.
       def_vector[]   - currently contains the additional def needed to 
                        take source_coord mid-way to the neighbour's
		        mean-point (which is now stored in target_coord[]).   */


    /* -------------------------------------------------------------- */
    /* get the world coord position of the node in the source volume,
       corresponding to the target node in question, taking into
       consideration the current warp  

       xp,yp,zp is the position in the source volume, corrsponding to
       the 'best target' position */


    if (ndim==3)
      general_inverse_transform_point(Gglobals->trans_info.transformation, 
				      target_coord[X],  target_coord[Y],  target_coord[Z],
				      &xp, &yp, &zp);
    else
      general_inverse_transform_point_in_trans_plane(Gglobals->trans_info.transformation, 
					 target_coord[X],  target_coord[Y],  target_coord[Z],
					 &xp, &yp, &zp);

    /* -------------------------------------------------------------- */
    /* BUILD THE SOURCE VOLUME LOCAL NEIGHBOURHOOD INFO:
          build the list of world coordinates representing nodes in the
          sub-lattice within the source volume                        */

    if (Gglobals->trans_info.use_magnitude) {

      /* build a spherical sub-lattice of Glen points in the source
         volume, note: sub-lattice diameter= 1.5*fwhm 
	 (here specified as 3*spacing = 2*(fwhm/2) to specify radius 
	  in build_source_lattice) 
      */

      build_source_lattice(xp, yp, zp, 
			   SX, SY, SZ,
			   spacing*3, spacing*3, spacing*3,
			   Diameter_of_local_lattice,  
			   Diameter_of_local_lattice,  
			   Diameter_of_local_lattice,
			   ndim, &Glen);
    }
    else {			/* use projections */

      /* store the coordinate of the center of the neighbourhood only,
         since we are using projections. */

      SX[1] = xp; SY[1] = yp; SZ[1] = zp;
      Glen = 1;
    }

    /* -------------------------------------------------------------- */
    /* BUILD THE TARGET VOLUME LOCAL NEIGHBOURHOOD INFO */

    if ( Gglobals->trans_info.use_magnitude ) {

      /* map this lattice forward into the target space, using the
	 current transformation, in order to build a deformed lattice
	 (in the WORLD COORDS of the target volume) */

      if (Gglobals->trans_info.use_super>0) 
	build_target_lattice_using_super_sampled_def(
                             SX,SY,SZ, TX,TY,TZ, Glen);
      else 
	build_target_lattice(SX,SY,SZ, TX,TY,TZ, Glen);

    }
    else {

      /* get only the target voxel position */

      if (number_dimensions==3)
	general_transform_point(Gglobals->trans_info.transformation, 
				xp,yp,zp,    &xt, &yt, &zt);
      else
	general_transform_point_in_trans_plane(Gglobals->trans_info.transformation, 
					xp,yp,zp,  &xt, &yt, &zt);
					
      convert_3D_world_to_voxel(Gglobals->features.model[0], xt, yt, zt,
				&Gtarget_vox_x, &Gtarget_vox_y, &Gtarget_vox_z);
    }


    /* -------------------------------------------------------------- */
    /* GET THE VOXEL COORDINATE LIST:
          for the objective function used in the actual optimization,
	  need the voxel coordinates of the target lattice.

	  note: I assume that the volume is stored in ZYX order !
	        (since I load the features in ZYX order in main() and
		in get_feature_volume()                               

		so, TX[] will store the voxel zdim position, TY with
		ydim, and TZ the voxel xdim coordinate.  BIZARRE I know,
		but it works... */

    for_inclusive(i,1,Glen) {
      convert_3D_world_to_voxel(Gglobals->features.model[0], 
				(Real)TX[i],(Real)TY[i],(Real)TZ[i], 
				&pos[0], &pos[1], &pos[2]);
      TX[i] = pos[0];
      TY[i] = pos[1];
      TZ[i] = pos[2];
    }

    /* -------------------------------------------------------------- */
    /* re-build the source lattice (without local neighbour warp),
       that will be used in the optimization below */

    if (Gglobals->trans_info.use_magnitude) {
      for_inclusive(i,1,Glen) {
	SX[i] += source_coord[X] - xp;
	SY[i] += source_coord[Y] - yp;
	SZ[i] += source_coord[Z] - zp;
      }
    }
    
    /* -------------------------------------------------------------- */
    /* GO GET FEATURES IN SOURCE VOLUME
          actually get the feature data from the source volume
	  local neighbourhood.  The target volume features are 
	  retrieved in go_get_samples_with_offset() called
	  from 1.0 - xcorr_fitting_function()                                       */

    if (Gglobals->trans_info.use_magnitude) {
      for_less(i,0, Gglobals->features.number_of_features) {
	go_get_samples_in_source(Gglobals->features.data[i], 
				 SX,SY,SZ, Ga1_features[i], Glen, 
				 (Gglobals->interpolant==nearest_neighbour_interpolant ? -1 : 0)
				 );
      }
    }
    else {			

      /* build projection data, ie go get f, df/dx, df/dy, df/dz to be
         able to reconstruct the equivalent of a sub-lattice, based
         only on these first four terms of the Taylor expansion. */

      evaluate_volume_in_world(Gglobals->features.data[0],
			       xp, yp, zp, 
			       0, TRUE, 0.0, val,
			       NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d1 = val[0];
      evaluate_volume_in_world(Gglobals->features.data[1],
			       xp, yp, zp, 
			       0, TRUE, 0.0, val,
			       NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d1x = 2.0*val[0];
      evaluate_volume_in_world(Gglobals->features.data[2],
			       xp, yp, zp, 
			       0, TRUE, 0.0, val,
			       NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d1y = 2.0*val[0];
      evaluate_volume_in_world(Gglobals->features.data[3],
			       xp, yp, zp, 
			       0, TRUE, 0.0, val,
			       NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
      Gproj_d1z = 2.0*val[0];

    }

    /* -------------------------------------------------------------- */
    /* calc one of the normalization coefficients for correlation,
       when using the magnitude data.  This saves a few CPU cycles in
       go_get_samples_with_offset(), since the constants only have to be
       eval'd once for the source volume. Note that this variable is not
       used when using projections. */

    for_less(i,0,Gglobals->features.number_of_features) {

      switch (Gglobals->features.obj_func[i]) {
      case NONLIN_XCORR:
	Gsqrt_features[i] = 0.0;
	for_inclusive(j,1,Glen)
	  Gsqrt_features[i] += Ga1_features[i][j]*Ga1_features[i][j];
	Gsqrt_features[i] = sqrt((double)Gsqrt_features[i]);
	break;
      case NONLIN_DIFF:
	Gsqrt_features[i] = (Real)Glen;
	break;
      case NONLIN_LABEL:
	Gsqrt_features[i] = (Real)Glen;
	break;
      default:
	print_error_and_line_num("Objective function %d not supported in go_get_samples_with_offset",__FILE__, __LINE__,Gglobals->features.obj_func[i]);
      }
    }

    /* -------------------------------------------------------------- */
    /*  FIND BEST DEFORMATION VECTOR
          now find the best local deformation that maximises the local
	  neighbourhood correlation between the source values stored in
	  **Ga1_features at positions Sx, SY, SZ with the homologous 
	  values at positions TX,TY,TZ in the target volume */

    if ( !Gglobals->trans_info.use_simplex) {

      /* ----------------------------------------------------------- */
      /*  USE QUADRATIC FITTING to find best deformation vector      */

      if (ndim==3) {
	/* build up the 3x3x3 matrix of local correlation values */

	for_inclusive(i,-1,1)
	  for_inclusive(j,-1,1)
	    for_inclusive(k,-1,1) {
	      pos_vector[1] = (float) i * Gsimplex_size/2.0;
	      pos_vector[2] = (float) j * Gsimplex_size/2.0;
	      pos_vector[3] = (float) k * Gsimplex_size/2.0;
	      local_corr3D[i+1][j+1][k+1] = xcorr_fitting_function(pos_vector); 
	    }
	*num_functions = 27;
	flag = return_3D_disp_from_min_quad_fit(local_corr3D, &du, &dv, &dw);
	
      }
      else {
	/* build up the 3x3 matrix of local correlation values */
	
	for_inclusive(i,-1,1)
	  for_inclusive(j,-1,1) {
	    pos_vector[1] = (float) i * Gsimplex_size/2.0;
	    pos_vector[2] = (float) j * Gsimplex_size/2.0;
	    pos_vector[3] = 0.0;
	    local_corr2D[i+1][j+1] = 1.0 - xcorr_fitting_function(pos_vector); 
	  }
	*num_functions = 9;
	
	flag = return_2D_disp_from_quad_fit(local_corr2D,  &du, &dv);
	dw = 0.0;
	
      }

      if ( flag ) {
	voxel_displacement[0] = dw * Gsimplex_size/2.0;	/* fastest (X) data index */
	voxel_displacement[1] = dv * Gsimplex_size/2.0;	/* Y */
	voxel_displacement[2] = du * Gsimplex_size/2.0;	/* slowest, Z */
      }
      else {
	result = -DBL_MAX;
	voxel_displacement[0] = 0.0;
	voxel_displacement[1] = 0.0;
	voxel_displacement[2] = 0.0;
      }
    }
    else {

      /* ----------------------------------------------------------- */
      /*  USE SIMPLEX OPTIMIZATION to find best deformation vector   */

				/* set up SIMPLEX OPTIMIZATION */
      nfunk = 0;

      ALLOC(parameters, ndim);	
      for_less(i,0,ndim)	/* init parameters for _NO_ deformation  */
	parameters[i] = 0.0;
				/* set the simplex diameter so as to 
				   reduce the size of the simplex, and
				   hence reduce the search space with
				   each iteration.                      */
      simplex_size = Gsimplex_size * 
	(0.5 + 
	 0.5*((Real)(total_iters-iteration)/(Real)total_iters));

      initialize_amoeba(&the_amoeba, ndim, parameters, 
			simplex_size, amoeba_obj_function, 
			NULL, (Real)ftol);
      
      /*                do the actual SIMPLEX optimization        */

      nfunk = 4;		/* since 4 eval's needed to init
				   the amoeba                     */

				/*  nfunk is incremented inside 
				    perform_amoeba  */
      while (nfunk < 400 && 
	     perform_amoeba(&the_amoeba, &nfunk) );
    
      if (nfunk<400) {

	get_amoeba_parameters(&the_amoeba,parameters);
	
	if (FLAG_HAHA) {
	  for_less( i, 0, the_amoeba.n_parameters+1 )  {
	    print ("%2d: %10.8f ",i,the_amoeba.values[i]);
	    for_less( j, 0, the_amoeba.n_parameters )
	      print ("%8.6f ",the_amoeba.parameters[i][j]);
	    print ("\n");
	  }
	}

	*num_functions = nfunk;
    
	if (ndim>2) {
	  voxel_displacement[0] = parameters[2]; /* fastest data index */
	  voxel_displacement[1] = parameters[1];
	  voxel_displacement[2] = parameters[0];
	}
	else {
	  voxel_displacement[0] = 0.0; /* corresponds to z-dim in data */
	  voxel_displacement[1] = parameters[1];
	  voxel_displacement[2] = parameters[0];
	}

      }
      else {
	      /* simplex optimization found nothing, so set the additional
		 displacement to 0 */

	voxel_displacement[0] = 0.0;
	voxel_displacement[1] = 0.0;
	voxel_displacement[2] = 0.0;
	result                = -DBL_MAX;
	*num_functions        = 0;

      } /*  if perform_amoeba */

      terminate_amoeba(&the_amoeba);      
      FREE(parameters);

    } /* else use_magnitude */


    /* -------------------------------------------------------------- */
    /* RETURN DEFORMATION FOUND 
          deformation calculated above is in voxel coordinates, it
	  has to be transformed into real world coordinates so that it
	  can be saved in the GRID_TRANSFORM                          */

    if ((voxel_displacement[0] == 0.0 &&
	 voxel_displacement[1] == 0.0 &&
	 voxel_displacement[2] == 0.0)) {
      
      def_vector[X] += 0.0;
      def_vector[Y] += 0.0;
      def_vector[Z] += 0.0;
    }
    else {

      convert_3D_world_to_voxel(Gglobals->features.model[0], 
				target_coord[X],target_coord[Y],target_coord[Z], 
				&voxel[0], &voxel[1], &voxel[2]);
      
      convert_3D_voxel_to_world(Gglobals->features.model[0], 
				(Real)(voxel[0]+voxel_displacement[0]), 
				(Real)(voxel[1]+voxel_displacement[1]), 
				(Real)(voxel[2]+voxel_displacement[2]),
				&pos[X], &pos[Y], &pos[Z]);
      
      def_vector[X] += pos[X]-target_coord[X];
      def_vector[Y] += pos[Y]-target_coord[Y];
      def_vector[Z] += pos[Z]-target_coord[Z];
      
    }

    result = sqrt((def_vector[X] * def_vector[X]) + 
		  (def_vector[Y] * def_vector[Y]) + 
		  (def_vector[Z] * def_vector[Z])) ;      
    
    if (result > 50.0 && 
	voxel_displacement[0]!=-20.0 && 
	voxel_displacement[1]!=-20.0 && 
	voxel_displacement[2]!=-20.0   ) {	/* this is a check, and should never
				   occur when using normal brain data,
				   unless we are scanning trauma
				   patients!  */

      print ("??? displacement way too big: %f\n",result);
      print ("vox_disp %f %f %f\n",
	     voxel_displacement[0],
	     voxel_displacement[1],
	     voxel_displacement[2]);
      print ("deform   %f %f %f\n", def_vector[X], def_vector[Y], def_vector[Z]);
      print ("Gsimplex_size = %f\n",Gsimplex_size);
    }


    /* build output file x,y & z to look at neighbourhood correlation data */

    if (FLAG_HAHA==TRUE) {
      print ("Flag_haha in get def for node\n");
      print ("def found was: %f %f %f -> %f %f %f\n", 
	     voxel_displacement[0], voxel_displacement[1], voxel_displacement[2],
	     def_vector[X], def_vector[Y], def_vector[Z]);
      
      pos_vector[3] = 0.0;
      pos_vector[2] = 0.0;
      pos_vector[1] = 0.0;
      print ("At 0,0,0 obj_fn        = %f\n",xcorr_fitting_function(pos_vector)); 
      pos_vector[3] = voxel_displacement[0];
      pos_vector[2] = voxel_displacement[1];
      pos_vector[1] = voxel_displacement[2];
      print ("At returned def obj_fn = %f\n",xcorr_fitting_function(pos_vector)); 
      if (open_file("x", WRITE_FILE, ASCII_FORMAT, &tmp_fp) != OK) {
	print_error_and_line_num("Can't open 'x' file\n",
				 __FILE__, __LINE__);
      }
      else {
	i = j = k = 0;
	for_inclusive(i,-200,200) {
	  pos_vector[1] = ((float)i/100.0) * Gsimplex_size/4.0;
	  pos_vector[2] = ((float)j/100.0) * Gsimplex_size/4.0;
	  pos_vector[3] = ((float)k/100.0) * Gsimplex_size/4.0;
	  fprintf (tmp_fp, "%f %f\n",(float)i/400.0,xcorr_fitting_function(pos_vector)); 
	}
	close_file(tmp_fp);
      }
      
      if (open_file("y", WRITE_FILE, ASCII_FORMAT, &tmp_fp) != OK) {
	print_error_and_line_num("Can't open 'y' file\n",
				 __FILE__, __LINE__);
      }
      else {
	i = j = k = 0;
	for_inclusive(j,-200,200) {
	  pos_vector[1] = ((float)i/100.0) * Gsimplex_size/4.0;
	  pos_vector[2] = ((float)j/100.0) * Gsimplex_size/4.0;
	  pos_vector[3] = ((float)k/100.0) * Gsimplex_size/4.0;
	  fprintf (tmp_fp, "%f %f\n",(float)j/400.0,xcorr_fitting_function(pos_vector)); 
	}
	close_file(tmp_fp);
      }
      
      if (open_file("z", WRITE_FILE, ASCII_FORMAT, &tmp_fp) != OK) {
	print_error_and_line_num("Can't open 'z' file\n",
				 __FILE__, __LINE__);
      }
      else {
	i = j = k = 0;
	for_inclusive(k,-200,200) {
	  pos_vector[1] = ((float)i/100.0) * Gsimplex_size/4.0;
	  pos_vector[2] = ((float)j/100.0) * Gsimplex_size/4.0;
	  pos_vector[3] = ((float)k/100.0) * Gsimplex_size/4.0;
	  fprintf (tmp_fp, "%f %f\n",(float)k/400.0,xcorr_fitting_function(pos_vector)); 
	}
	close_file(tmp_fp);
      }
    }

    

  } /*  if (!get_best_start_from_neighbours()) */
  
  return(result);
}


/* Build the target lattice by transforming the source points through the
   current non-linear transformation stored in:

        Gglobals->trans_info.transformation                        

*/
public void    build_target_lattice(float px[], float py[], float pz[],
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
      general_transform_point_in_trans_plane(Gglobals->trans_info.transformation, 
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
public void    build_target_lattice_using_super_sampled_def(
                                     float px[], float py[], float pz[],
				     float tx[], float ty[], float tz[],
				     int len)
{
  int 
    i,j,
    sizes[MAX_DIMENSIONS],
    xyzv[MAX_DIMENSIONS];
  Real 
    def_vector[N_DIMENSIONS],
    voxel[MAX_DIMENSIONS],
    x,y,z;
  long 
    index[MAX_DIMENSIONS];

  get_volume_sizes(Gsuper_sampled_vol,sizes);
  get_volume_XYZV_indices(Gsuper_sampled_vol,xyzv);

  for_inclusive(i,1,len) {

				/* apply linear part of the transformation */

    general_transform_point(Glinear_transform,
			    (Real)px[i], (Real)py[i], (Real)pz[i], 
			    &x, &y, &z);

				/* now get the non-linear part, using
                                   nearest neighbour interpolation in
                                   the super-sampled deformation
                                   volume. */

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

/* This is the COST FUNCTION TO BE MINIMIZED.
   so that very large displacements are impossible */

private Real cost_fn(float x, float y, float z, Real max_length)
{
  Real v2,v,d;

  v2 = x*x + y*y + z*z;
  v = sqrt(v2);

  v *= v2;

  if (v<max_length)
    d = 0.2 * v / (max_length - v);
  else
    d = 1e38;

  return(d);
}

/* This is the SIMILARITY FUNCTION TO BE MAXIMIZED.

      it is maximum when source and target data are most similar.
      The value is normalized by the weights associated with each feature,
      so that 0 <= similarity_fn() <= 1.0 

      the input parameter D stores the displacement offsets for the
      target lattice:
         D[1] stores the xdisp in voxels
         D[2] stores the ydisp
         D[3] stores the zdisp
*/

private Real similarity_fn(float *d)
{
  int i;
  Real
    norm,
    val[MAX_DIMENSIONS],
    voxel[MAX_DIMENSIONS],
    xw,yw,zw,
    s, s1, s2;

  /* note: here the displacement order for go_get_samples_with_offset
     is 3,2,1 (=Z,Y,X) since the source and target volumes are stored in
     Z,Y,X order.

     This backward ordering is necessary since the same function is
     used for 2D and 3D optimization, and simplex will only modify d[1]
     and d[2] in 2D (d[3] stays const=0). */

  if (Gglobals->trans_info.use_magnitude) {

    s = norm = 0.0;
    
    for_less(i,0,Gglobals->features.number_of_features)  {
      norm += ABS(Gglobals->features.weight[i]);
      s += Gglobals->features.weight[i] * 
	(Real)go_get_samples_with_offset(Gglobals->features.model[i],
					 TX,TY,TZ,
					 d[3], d[2], d[1],
					 Gglobals->features.obj_func[i],
					 Glen, 
					 Gsqrt_features[i], Ga1_features[i],
					 Gglobals->interpolant==nearest_neighbour_interpolant);

    }

    if (norm != 0.0) 
      s = s / norm;
    else
      print_error_and_line_num("The feature weights are null.", 
		__FILE__, __LINE__);


  } else {
	/* calc correlation based on projection data */
    voxel[0] = Gtarget_vox_x + d[3];
    voxel[1] = Gtarget_vox_y + d[2]; 
    voxel[2] = Gtarget_vox_z + d[1];
    voxel[3] = 0.0;
    voxel[4] = 0.0;

    convert_voxel_to_world(Gglobals->features.model[0], voxel, &xw, &yw, &zw);

    evaluate_volume_in_world(Gglobals->features.model[0], xw, yw, zw, 
			     0, TRUE, 0.0, val,
			     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    Gproj_d2 = val[0];
    evaluate_volume_in_world(Gglobals->features.model[1], xw, yw, zw, 
			     0, TRUE, 0.0, val,
			     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    Gproj_d2x = 2.0*val[0];
    evaluate_volume_in_world(Gglobals->features.model[2], xw, yw, zw, 
			     0, TRUE, 0.0, val,
			     NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    Gproj_d2y = 2.0*val[0];
    evaluate_volume_in_world(Gglobals->features.model[3], xw, yw, zw, 
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


private Real amoeba_obj_function(void * dummy, float d[])
{
  int i;
  float p[4];
  Real obj_func_val;


  for_less(i,0,number_dimensions)
    p[i+1] = d[i];
  obj_func_val =  xcorr_fitting_function(p);

  if (FLAG_HAHA) {
    print ("\n%15.12f: %12.8f  %12.8f  %12.8f: ", obj_func_val, d[0], d[1], d[2]);
  }

  return ( obj_func_val );
  
}

private Real xcorr_fitting_function(float *d)
     
{
  Real
    similarity,
    cost, 
    r;
  
  similarity = (Real)similarity_fn( d );
  cost       = (Real)cost_fn( d[1], d[2], d[3], Gcost_radius );
  
  r = 1.0 - 
      similarity * similarity_cost_ratio + 
      cost       * (1.0-similarity_cost_ratio);
  
  return(r);
}

