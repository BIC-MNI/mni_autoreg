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
@MODIFIED   : Revision 96.8  2000-03-15 08:42:46  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 96.7  2000/01/18 18:53:37  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.6  1999/10/25  19:59:07  louis
 * final checkin before switch to CVS
 *
 * Revision 96.5  1999/06/10  12:51:23  louis
 * update with optical flow working in addition to xcorr, label, and diff
 * sub-lattice computed only if needed
 *
 * Revision 96.4  1999/06/09  13:11:08  louis
 * working version with optical flow (working by itself).
 *
 * Revision 96.3  1997/11/12  21:07:43  louis
 * - added support for chamfer distance as a local obj func
 * - moved all procedures used to compute the local obj function
 *   into def_obj_funcitons.c
 *   cost_fn(), similarity_fn(), go_get_samples_with_offset(),
 *   local_objective_function(), amoeba_NL_obj_function()
 *
 * Revision 96.2  1997/11/03  20:05:41  louis
 * reorganized do_nonlinear.c
 *  - added prototypes in header files
 *    sub_lattice.h
 *    extras.h
 *    quad_max_fit.h
 *  - moved functions
 *    build_target_lattice()
 *    build_target_lattice_using_super_sampled_def()
 *    build_source_lattice()
 *    into sub_lattice.c
 *
 * Revision 96.1  1997/11/03  15:06:29  louis
 * working version, before creation of mni_animal package, and before inserting
 * distance transforms
 *
 * Revision 96.0  1996/08/21  18:22:10  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:02  louis
 * Pre-release
 *
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
 * local_objective_function() so that both the similarity of the data and the
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/do_nonlinear.c,v 96.8 2000-03-15 08:42:46 stever Exp $";
#endif

#include <config.h>		/* MAXtype and MIN defs                      */
#include <internal_volume_io.h>	/* structs & tools to deal with volumes data */
#include <amoeba.h>		/* simplex optimization struct               */

#include <stdlib.h>		/* to get header info for drand48()          */
#include <arg_data.h>		/* definition of the global data struct      */
#include <print_error.h>	/* def of print_error_and_..                 */
#include <deform_support.h>	/* prototypes for routines called
				   from deformation procedures.              */

#include "constants.h"		/* internal constant definitions             */
#include "interpolation.h"
#include "super_sample_def.h"
#include <sys/types.h>		/* for timing the deformations               */
#include <time.h>
time_t time(time_t *tloc);

#include <stats.h>              /* for local stats computations              */
#include <sub_lattice.h>        /* prototypes for sub_lattice manipulation   */
#include <extras.h>             /* prototypes for extra convienience routines*/
#include <quad_max_fit.h>       /* prototypes for quadratic fitting routines */

int stat_quad_total;            /* these are used as globals to tally stats  */
int stat_quad_zero;             /* in Numerical/quad_max_stats.c             */
int stat_quad_two;              /* (mostly for debugging)                    */
int stat_quad_plus;
int stat_quad_minus;
int stat_quad_semi;

                                /* these globals are used to tally stats over
                                   do_non_linear_optimization() and 
                                   return_locally_smoothed_def               */
static stats_struct
   stat_def_mag,
   stat_num_funks,
   stat_conf0,
   stat_conf1,
   stat_conf2,
   stat_eigval0,
   stat_eigval1,
   stat_eigval2;

                                /* these Globals are used to communicate to
                                   the correlation functions over top the
                                   SIMPLEX optimization routine */
float	
  *Gsqrt_features,		/* normalization const for correlation       */
  **Ga1_features,		/* samples in source sub-lattice             */
  *TX, *TY, *TZ,		/* sample sub-lattice positions in target    */
  *SX, *SY, *SZ;		/* sample sub-lattice positions in source    */

int 
  Glen;				/* # of samples in sub-lattice               */


         /* these Globals are used to communicate the projection */
         /* values over top the SIMPLEX optimization  routine    */ 
Real
  Gtarget_vox_x, Gtarget_vox_y, Gtarget_vox_z,
  Gproj_d1,  Gproj_d1x,  Gproj_d1y,  Gproj_d1z, 
  Gproj_d2,  Gproj_d2x,  Gproj_d2y,  Gproj_d2z;

	/* Globals used for local simplex Optimization  */
static Real     Gsimplex_size;	/* the radius of the local simplex           */
Real     Gcost_radius;	/* constant used in the cost function        */

        /* Globals used to split the input transformation into a
	   linear part and a super-sampled non-linear part */
General_transform 
                *Gsuper_sampled_warp,
                *Glinear_transform;
Volume  Gsuper_sampled_vol;

	/* Volume order definition for super sampled data */
static char *my_XYZ_dim_names[] = { MIxspace, MIyspace, MIzspace };


        /* program Global data used to store all info regarding data
           and transformations  */
Arg_Data *Gglobals;


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

public  Real amoeba_NL_obj_function(void * dummy, float d[]);

#define AMOEBA_ITERATION_LIMIT  400 /* max number of iterations for amoeba */

public  Real local_objective_function(float *x);

private Real get_deformation_vector_for_node(Real spacing, Real threshold1, 
					     Real source_coord[],
					     Real mean_target[],
					     Real def_vector[],
					     Real voxel_displacement[],
					     int iteration, int total_iters,
					     int *nfunks,
					     int ndim,
					     BOOLEAN sub_lattice_needed);

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



private BOOLEAN get_best_start_from_neighbours(
					       Real threshold1, 
					       Real source[],
					       Real mean_target[],
					       Real target[],
					       Real def[]);

  

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

private BOOLEAN is_a_sub_lattice_needed (char obj_func[],
					 int  number_of_features);

private BOOLEAN build_lattices(Real spacing, 
			       Real threshold, 
			       Real source_coord[],
			       Real mean_target[],
			       Real target_coord[],
			       Real def_vector[],
			       int ndim);

private Real get_optical_flow_vector(Real threshold1, 
				     Real source_coord[],
				     Real mean_target[],
				     Real def_vector[],
				     Real voxel_displacement[],
				     Volume data,
				     Volume model,
				     int ndim);

private Real get_chamfer_vector(Real threshold1, 
				Real source_coord[],
				Real mean_target[],
				Real def_vector[],
				Real voxel_displacement[],
				Volume data,
				Volume chamfer,
				int ndim);

/**************************************************************************/
public Status do_non_linear_optimization(Arg_Data *globals)
{
   General_transform
      *all_until_last,		/* will contain the first (linear) part of the xform  */
      *additional_warp,		/* storage of estimates of the needed additional warp */
      *another_warp,		/* storage of estimates of the needed additional warp */
      *current_warp;		/* storage of the current (best) warp so far          */
   
   Volume
      estimated_flag_vol,	/* flags indicating node estimated or not             */
      additional_vol,		/* volume pointer to additional_warp transform        */
      another_vol,		/* volume pointer to additional_warp transform        */
      current_vol,		/* volume pointer to current_warp transform           */
      additional_mag;		/* volume storing mag of additional_warp vectors      */
  
   long
      iteration_start_time,	/* variables to time each iteration                   */
      temp_start_time,
      timer1,timer2,
      nfunk_total;

   int 
      additional_count[MAX_DIMENSIONS], /* size (in voxels) of  additional_vol  */
      mag_count[MAX_DIMENSIONS],/* size (in voxels) of  additional_mag          */
      xyzv[MAX_DIMENSIONS],	/* order of voxel indices                       */
      index[MAX_DIMENSIONS],	/* used to step through all nodes of def field  */
      start[MAX_DIMENSIONS],	/* starting limit of index[]                    */
      end[MAX_DIMENSIONS],	/* ending limit of index[]                      */
      debug_sizes[MAX_DIMENSIONS],
      iters,			/* iteration counter */
      i,j,k,
      nodes_done, nodes_tried,	/* variables to calc stats on deformation estim  */
      nodes_seen, over,
      nfunks, nfunk1, nodes1,
      sub_lattice_needed;

   Real 
      debug_steps[MAX_DIMENSIONS], 
      voxel[MAX_DIMENSIONS],	/* voxel position used for interpolation         */
      mag_steps[MAX_DIMENSIONS],
      steps[MAX_DIMENSIONS],	/* voxel size (width,height,length) of def field */
      steps_data[MAX_DIMENSIONS], /* voxel size of data                          */
      
				/* variables to calc stats on deformation estim  */
      mag, mean_disp_mag, std, 

      another_vector[3],	
      current_def_vector[3],	/* the current deformation vector for a  node    */
      result_def_vector[3],	/* the smoothed deformation vector for a node    */
      def_vector[3],		/* the additional deformation estimated for node,
				   for these three vectors in real world coords
				   index [0] stores the x disp, and [2] the z disp */
      voxel_displacement[3],    /* the additional displacement, in voxel coords 
				   with [0] storing the displacement of the
				   fastest varying index, and [2] the slowest in
				   the data voluming = xdim and zdim respectively*/
      wx,wy,wz,			/* temporary storage for a world coordinate      */
      source_node[3],		/* world coordinate of source node               */
      target_node[3],		/* world coordinate of corresponding target node */
      mean_target[3],		/* mean deformed pos, determined by neighbors    */
      mean_vector[3],		/* mean deformed vector, determined by neighbors */
      threshold1,		/* intensity thresh for source vol               */
      threshold2,		/* intensity thresh for target vol               */
      result,			/* magnitude of additional warp vector           */
      eig1;

   progress_struct		/* to print out program progress report */
      progress;

  /*******************************************************************************/

	   /* set up globals for communication with other routines */

   Gglobals= globals;

   /* allocate space required for some globals */

   if (Gglobals->features.number_of_features > 0) {
      ALLOC2D (Ga1_features, Gglobals->features.number_of_features, MAX_G_LEN+1);
      ALLOC( Gsqrt_features, Gglobals->features.number_of_features);

      sub_lattice_needed = is_a_sub_lattice_needed (Gglobals->features.obj_func,
						    Gglobals->features.number_of_features);

      if (globals->flags.debug) {
        print ("There are %d feature pairs\n",Gglobals->features.number_of_features);
        for_less(i,0,Gglobals->features.number_of_features) {
	  print ("%d: [%d] [%7.5f] %s <-> %s\n", i, 
		 Gglobals->features.obj_func[i],
		 Gglobals->features.weight[i],
		 Gglobals->features.data_name[i],
		 Gglobals->features.model_name[i]);
        }
	if (sub_lattice_needed) 
	  print ("A sub-lattice is needed for at least one feature\n");
	else
	  print ("No sub-lattice needed.  Should be a fast run!\n");

      } 
   }
   else {			/* we should never get here */
     print_error_and_line_num("There are no features to match in non_lin optimization\n", 
			      __FILE__, __LINE__);
   }

   ALLOC(SX,MAX_G_LEN+1);        /* and coordinates in source volume  */
   ALLOC(SY,MAX_G_LEN+1);
   ALLOC(SZ,MAX_G_LEN+1);
   ALLOC(TX,MAX_G_LEN+1);        /* and coordinates in target volume  */
   ALLOC(TY,MAX_G_LEN+1);
   ALLOC(TZ,MAX_G_LEN+1);

   /* split the total transformation into the first linear part and the
      last non-linear def.  */
   split_up_the_transformation(globals->trans_info.transformation,
                               &all_until_last,
                               &current_warp);

				/* exit if no deformation */
   if (current_warp == (General_transform *)NULL) { 
      print_error_and_line_num("Cannot find the deformation field to optimize at end of input transform",
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

   /* make a copy of the current warp that will be used to store the
      additional deformation needed to optimize the current warp.  */
   ALLOC(additional_warp,1);
   copy_general_transform(current_warp, additional_warp);
   current_vol = current_warp->displacement_volume;
   additional_vol = additional_warp->displacement_volume;
   get_volume_sizes(additional_vol, additional_count);
   get_volume_XYZV_indices(additional_vol, xyzv);
   init_the_volume_to_zero(additional_vol);    /* reset it to zero */


   /* build a debugging volume/xform */
   ALLOC(another_warp,1);	
   copy_general_transform(current_warp, another_warp);
   another_vol = another_warp->displacement_volume;
   set_volume_real_range(another_vol, -20.0 , 20.0);
   init_the_volume_to_zero(another_vol);

   /* build a temporary volume that will be used to store the magnitude of
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

   /* build a volume to flag all voxels where the deformation has been
      estimated */
   estimated_flag_vol = create_volume(3, my_XYZ_dim_names, NC_BYTE, FALSE, 0.0, 0.0);
   for_less(i,0,3)
      mag_count[i] = additional_count[ xyzv[i] ];
   set_volume_sizes(estimated_flag_vol, mag_count);
   set_volume_separations(estimated_flag_vol, mag_steps);
   set_volume_translation(estimated_flag_vol, voxel, target_node);
   alloc_volume_data(estimated_flag_vol);
   set_volume_real_range(estimated_flag_vol, 0.0, 1.0);
 
   /* test simplex size against size of data voxels */

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
      print_error_and_line_num("Zero voxel size for model data: %f %f %f\n", 
                               __FILE__, __LINE__,
                               steps_data[0],steps_data[1],steps_data[2]);


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

  initial_corr = xcorr_objective_with_def(Gglobals->features.data[0], 
                                          Gglobals->features.model[0],
                                          Gglobals->features.data_mask[0], 
                                          Gglobals->features.model_mask[0],
                                          globals );

  if (globals->flags.debug) {	
    print("\n\nDebug info from do_nonlinear_optimization---------------\n");
    print("Initial corr         = %f\n",initial_corr);
    print("Source vol threshold = %f\n", threshold1);
    print("Target vol threshold = %f\n", threshold2);
    print("Iteration limit      = %d\n", iteration_limit);
    print("Iteration weight     = %f\n", iteration_weight);
    print("xyzv                 = %3d %3d %3d %3d \n",
	  xyzv[X], xyzv[Y], xyzv[Z], xyzv[Z+1]);
    print("number_dimensions    = %d\n",number_dimensions);
    print("smoothing_weight     = %f\n",smoothing_weight);
    print("loop                 = (%d %d) (%d %d) (%d %d)\n",
	  start[0],end[0],start[1],end[1],start[2],end[2]);


    print ("\nFitting STRATEGY ----------\n");
    
    if ( Gglobals->trans_info.use_magnitude) {

      if ( Gglobals->trans_info.use_simplex) {
	print ("  This fit will use local simplex optimization and\n");
	print ("  Simplex radius = %7.2f (voxels) or %7.2f(mm)\n",
	       Gsimplex_size, 
	       Gsimplex_size * MAX3(steps_data[0],steps_data[1],steps_data[2]));      }
      else {
	print ("  This fit will use local quadratic fitting and\n");
	print ("  Search/quad fit radius= %7.2f (data voxels) or %7.2f(mm)\n",
	       Gsimplex_size /2.0, 
	       Gsimplex_size * MAX3(steps_data[0],steps_data[1],steps_data[2])/2.0);
      }
    }
    else {
	print ("  This fit will use optical flow for direct fitting\n");        
    }

    if (Gglobals->trans_info.use_local_smoothing) {
      if ( Gglobals->trans_info.use_local_isotropic) 
	print ("    local isotroptic smoothing.\n");
      else
	print ("    local non-isotroptic smoothing.\n");
    }
    else {
      print ("    global smoothing.\n");
    }
    

    if (Gglobals->interpolant==nearest_neighbour_interpolant) 
      print ("  The similarity function will be evaluated using NN interpolation\n");
    else
      print ("  The similarity function will be evaluated using tri-linear interpolation\n");
    
    if ( Gglobals->trans_info.use_magnitude) {
      print ("    on a spherical sub-lattice with a diameter of %7.2f (data voxels),\n",
	     3.0*steps[xyzv[X]]/MAX3(steps_data[0],steps_data[1],steps_data[2]) );          
      print ("    %7.2f(mm) or %d nodes (%7.2f vox/node or%7.2f mm/node) \n", 
	     3.0*steps[xyzv[X]], Diameter_of_local_lattice,
	     3.0*steps[xyzv[X]]/((Diameter_of_local_lattice-1)*MAX3(steps_data[0],steps_data[1],steps_data[2])),
	     3.0*steps[xyzv[X]]/(Diameter_of_local_lattice-1));
    }
    
    print ("-----------------------\n");
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

   mean_disp_mag = 0.0;

  for_less(iters,0,iteration_limit) {

    iteration_start_time = time(NULL);

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
      if (globals->flags.debug)
         report_time(temp_start_time, "TIME:Interpolating super-sampled data");

    }  

    init_the_volume_to_zero(estimated_flag_vol);
    
    print("Iteration %2d of %2d\n",iters+1, iteration_limit);


				/* for various stats on this iteration*/
    stat_quad_total = 0;
    stat_quad_zero  = 0;
    stat_quad_two   = 0;
    stat_quad_plus  = 0;
    stat_quad_minus = 0;
    stat_quad_semi  = 0;

    nodes_done      = 0; 
    nodes_tried     = 0; 
    nodes_seen      = 0; 
    over            = 0;        
    nfunk_total     = 0;
    std             = 0.0;

    init_stats(&stat_def_mag,  "def_mag");
    init_stats(&stat_num_funks,"num_funks");
    init_stats(&stat_eigval0,  "eigval[0]");
    init_stats(&stat_eigval1,  "eigval[1]");
    init_stats(&stat_eigval2,  "eigval[2]");
    init_stats(&stat_conf0,    "conf[0]");
    init_stats(&stat_conf1,    "conf[1]");
    init_stats(&stat_conf2,    "conf[2]");



    initialize_progress_report( &progress, FALSE, 
			       (end[X]-start[X])*(end[Y]-start[Y]) + 1,
			       "Estimating deformations" );

    temp_start_time = time(NULL);

    for_less(i,0,MAX_DIMENSIONS) index[i]=0;

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


	      result = get_deformation_vector_for_node(steps[xyzv[X]], 
						       threshold1,
						       source_node,
						       mean_target,
						       def_vector,
						       voxel_displacement,
						       iters, iteration_limit, 
						       &nfunks,
						       number_dimensions,
						       sub_lattice_needed);
	    
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
          ALLOC(filenamestring,512);
          sprintf (filenamestring,"pr_axes2_%d.xfm",iters);
          (void)output_transform_file(filenamestring,NULL,another_warp);
          FREE(filenamestring);
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
       
      print ("Nodes seen = %d: [no def = %d], [w/def = %d (over = %d)]\n",
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
      
       report_time(temp_start_time, "TIME:Estimation defs");
      
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

       /* extrapolate (and smooth) the newly estimated deformation vectors
          (stored in additional vol) to un-estimated nodes, leaving the
          extrapolated result in additional_vol (note that it still has to
          be added to current */


       temp_start_time = time(NULL);

       extrapolate_to_unestimated_nodes(current_warp,
                                        additional_warp,
                                        estimated_flag_vol);
       if (globals->flags.debug) 
         report_time(temp_start_time, "TIME:Extrapolating the current warp");


       /* current = current + additional */

       temp_start_time = time(NULL);

       add_additional_warp_to_current(current_warp,
                                      additional_warp,
                                      1.0);
       if (globals->flags.debug) 
         report_time(temp_start_time, "TIME:Adding additional to current");


    }
    else {

       /* additional = additional + current, 
          and then apply global  smoothing */
       
       temp_start_time = time(NULL);
       add_additional_warp_to_current(additional_warp,
                                      current_warp,
                                      iteration_weight);
       if (globals->flags.debug) 
         report_time(temp_start_time, "TIME:Adding additional to current");
       

				/* smooth the warp in additional,
				   leaving the result in current 

				   current = smooth(additional) */

       temp_start_time = time(NULL);
       
       smooth_the_warp(current_warp,
                       additional_warp,
                       additional_mag, -1.0);
      
       if (globals->flags.debug) 
         report_time(temp_start_time, "TIME:Smoothing the current warp");

    }
    

    /* 
       clamp the data so that the 1st derivative of the deformation
       field does not exceed 1.0*step in magnitude

       clamp_warp_deriv(current->dx, current->dy, current->dz); 
    */

    
/*    if (FALSE && iters<iteration_limit-1) {  !!! */

    if (FALSE) { /* !!! */
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
							   number_dimensions,
							   sub_lattice_needed);

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


       final_corr = xcorr_objective_with_def(Gglobals->features.data[0], Gglobals->features.model[0],
                                    Gglobals->features.data_mask[0], Gglobals->features.model_mask[0],
                                    globals );
       print("initial corr %f ->  this step %f\n",
             initial_corr,final_corr);

       report_time(iteration_start_time, "TIME:This iteration");

    }
  }
  

                                /* now we are done with the iterations,
                                   compute the final correlation if we
                                   haven't already done so just above in
                                   the debug statement */
  if (!globals->flags.debug)
    final_corr = xcorr_objective_with_def(Gglobals->features.data[0], 
                                 Gglobals->features.model[0],
				 Gglobals->features.data_mask[0], 
                                 Gglobals->features.model_mask[0],
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



/*   look though the list of object functions requested,
     and set is_a_sub_lattice_needed=TRUE if any obj function
     is used other than Optical Flow
*/
private BOOLEAN is_a_sub_lattice_needed (char obj_func[],
					 int  number_of_features) {
  BOOLEAN needed;
  int i;

  needed = FALSE;
  for_less (i,0,number_of_features) {
    if (obj_func[i] != NONLIN_OPTICALFLOW && obj_func[i] != NONLIN_CHAMFER) 
      needed = TRUE;
  }
  return(needed);
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
	    pos_vector[1] = (float) (i * Gsimplex_size)/2.0;
	    pos_vector[2] = (float) (j * Gsimplex_size)/2.0;
	    pos_vector[3] = (float) (k * Gsimplex_size)/2.0;

	    local_corr3D[i+1][j+1][k+1] = local_objective_function(pos_vector); 

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
      another_vector[X] = conf[0];
      another_vector[Y] = conf[1];
      another_vector[Z] = conf[2]; 
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

/* use optical flow to compute the deformation (def_vector) required
to warp the source_coord onto the mean_target using the intensities of
the source_coord in the source volume and the mean_target coord in the
target volume, as well as the 1st intensity derivatives in the source
volume.  

inputs:
    threshold1,
    source_coord
    mean_target
    ndim

outputs:
    def_vector         (in world coords)
    voxel_displacement (in voxel coords)
    num_functions      (for number of function evaluations)

Based on Horn and Schunck Artificial Intell 17 (1981) 185-203
*/

#define Min_deriv  0.02

private Real get_optical_flow_vector(Real threshold1, 
				     Real source_coord[],
				     Real mean_target[],
				     Real def_vector[],
				     Real voxel_displacement[],
				     Volume data,
				     Volume model,
				     int ndim)
{ 
  Real
    val[MAX_DIMENSIONS],	/* the interpolated intensity value     */
    result,			/* the magnitude of the estimated def   */
    min,max,thresh,             /* volume real min, max and estimate on smallest
				   derivative */
    xp, yp, zp,			/* temp storage for coordinate position */
    dx[MAX_DIMENSIONS],	        /* derivative in X (world-coord)        */
    dy[MAX_DIMENSIONS],	        /*     "         Y                      */
    dz[MAX_DIMENSIONS],	        /*     "         Z                      */
    steps[MAX_DIMENSIONS];
  int
    i;				/* a counter                            */

    get_volume_separations(data, steps);
  
    xp = mean_target[0];	/* get intensity and derivatives        */
    yp = mean_target[1];        /* in target volume                     */
    zp = mean_target[2];
    evaluate_volume_in_world(model,
			     xp, yp, zp,
			     0, TRUE, 0.0, val,
			     dx,dy,dz,
			     NULL,NULL,NULL,NULL,NULL,NULL);
    Gproj_d2 = val[0];
    
    xp = source_coord[0];	/* get intensity only                   */
    yp = source_coord[1];	/* in source volume                     */
    zp = source_coord[2];
    evaluate_volume_in_world(data,
			     xp, yp, zp,
			     0, TRUE, 0.0, val,
                             NULL,NULL,NULL,
			     NULL,NULL,NULL,
			     NULL,NULL,NULL);
    Gproj_d1 = val[0];
    
				/* compute deformations directly!       */

    get_volume_real_range(model, &min, &max);
    thresh = Min_deriv * (max - min);

    if (ABS(dx[0]) > thresh)   /* fastest (X) */
      def_vector[0] = ((Gproj_d1 - Gproj_d2) /  dx[0]);
    else
      def_vector[0] = 0.0;
    
    if (ABS(dy[0]) > thresh)
      def_vector[1] = ((Gproj_d1 - Gproj_d2) /  dy[0]);
    else
      def_vector[1] = 0.0;
    
    if (ABS(dz[0]) > thresh && ndim==3)  /* slowest  (Z) */
      def_vector[2] = ((Gproj_d1 - Gproj_d2) /  dz[0]);
    else
      def_vector[2] = 0.0;
    
    
    for_less(i,0,3)            /* build the real-world displacement */
      voxel_displacement[i] = def_vector[i]  * steps[i];
    
    result = sqrt (def_vector[0]*def_vector[0] +
		   def_vector[1]*def_vector[1] +
		   def_vector[2]*def_vector[2]);
    
    return(result);

}

/* compute the deformation only if the source_coord is on a surface voxel

   use the chamfer volume (an approximation for distance) to determine how
   far (and in which direction) the mean_target is from the sync in the
   target


   modified 9/21/1999 LC:

      compute the deformation for this node from the nearest surface.
      ie - if there is a surface within the threshold distance (= .5
      spacing between nodes) then figure out what displacement would
      be applied to that coordinate, and return that value for this
      node.

*/

#define MAX_CAPTURE 3.8		

private Real get_chamfer_vector(Real capture_limit, 
				Real source_coord[],
				Real mean_target[],
				Real def_vector[],
				Real voxel_displacement[],
				Volume data,
				Volume chamfer,
				int ndim)
{ 
  Real
    
    sx,sy,sz,tx,ty,tz,          /* source and target coords             */
    dist, thresh, min, max,
    dist_weight,
    zero,
    result,			/* the magnitude of the estimated def   */
    val[MAX_DIMENSIONS],	/* the interpolated intensity value     */
    mag, gx,gy,gz,
    dx[MAX_DIMENSIONS],	        /* derivative in X (world-coord)        */
    dy[MAX_DIMENSIONS],	        /*     "         Y                      */
    dz[MAX_DIMENSIONS],	        /*     "         Z                      */
    steps[MAX_DIMENSIONS];      /* voxel size of data volume            */
  int
    i;				/* a counter                            */

  
                                /* get intensity in source volume       */
  evaluate_volume_in_world(data,
                           source_coord[0], source_coord[1], source_coord[2],
                           0, TRUE, 0.0, val,
                           NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  
  dist = val[0];

  if (dist > MAX_CAPTURE*capture_limit) { /* we are not near a surface */
    result = 0.0;
    for_less(i,0,3) {
        voxel_displacement[i] = 0.0;
	def_vector[i] = 0.0;
    }
  }
  else {     

    if (dist < capture_limit) {
      dist_weight = 1.0;
    }
    else {
      dist_weight = 1.0 - ((dist - capture_limit) / (capture_limit* (MAX_CAPTURE-1.0)));
    }

    sx = source_coord[0];
    sy = source_coord[1];
    sz = source_coord[2];
    tx = mean_target[0];
    ty = mean_target[1];
    tz = mean_target[2];
    
    zero = CONVERT_VOXEL_TO_VALUE(data, 0.0);

    if (val[0]!=zero) {		/* we are not already on  the surface */

				/* get derivative of the data volume */
      evaluate_volume_in_world(data,
			       sx, sy, sz,
			       0, TRUE, 0.0, val,
			       dx,dy,dz,
			       NULL,NULL,NULL,
			       NULL,NULL,NULL);
      
      mag = sqrt (dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0]);
				/* if mag > 0, then we can compute
				   a normalized vector in the gradient
				   direction */
      if (mag > 0.0) {
	
	gx = dx[0] / mag;
	gy = dy[0] / mag;
	gz = dz[0] / mag;

	/* use derivative info to find nearest surface point.   */


	/*
print ("%5.3f: %7.2f %7.2f %7.2f -> %7.2f %7.2f %7.2f [%7.2f %7.2f %7.2f ]",
       dist,
       sx,sy,sz, tx,ty,tz, gx,gy,gz);
	*/

	sx += -1.0 * dist * gx;
	sy += -1.0 * dist * gy;
	sz += -1.0 * dist * gz;
      
	/* sx,sy,sz is now on the closest surface in the data volume,
	   we now need the equivalent target coord */

	general_transform_point(Gglobals->trans_info.transformation, 
			      sx,sy,sz,  &tx,&ty,&tz);


  evaluate_volume_in_world(data,
                           sx,sy,sz,
                           0, TRUE, 0.0, val,
                           NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
  
  /*print ("%10.5f: %7.2f %7.2f %7.2f -> %7.2f %7.2f %7.2f: %5.3f:  ",
       dist,sx,sy,sz, tx,ty,tz, val[0]);
  */
        /* Now, both sx,sy,sz AND tx,ty,tz are offset by the same amount, so
	   that the following optimization will work. 

           Essentially, the deformation vector associated with the offset
           point will be returned for the original point sx,sy,sz   */

      }
    }
     
     /* get intensity in target volume       */
     evaluate_volume_in_world(chamfer,
                              tx, ty, tz,
                              0, TRUE, 0.0, val,
                              dx,dy,dz,
                              NULL,NULL,NULL,
                              NULL,NULL,NULL);
     

     get_volume_real_range(chamfer, &min, &max);
     thresh = Min_deriv * (max - min);

     thresh = 0.1;		/* gradient threshold... testing */

     /* compute deformations directly!       */

     /* print ("%7.2f: %7.2f %7.2f\n",dist,dist_weight,capture_limit); */
     
     if (ABS(dx[0]) > thresh)              /* fastest (X) */
        def_vector[0] = -1.0 * dist_weight * val[0] /  dx[0];
     else
        def_vector[0] = 0.0;
     
     if (ABS(dy[0]) > thresh)
        def_vector[1] = -1.0 * dist_weight * val[0] /  dy[0];
     else
        def_vector[1] = 0.0;
     
     if (ABS(dz[0]) > thresh && ndim==3)  /* slowest  (Z) */
        def_vector[2] = -1.0 * dist_weight * val[0] /  dz[0];
     else
        def_vector[2] = 0.0;


     get_volume_separations(data, steps);
     for_less(i,0,3)            /* build the real-world displacement */
        voxel_displacement[i] = def_vector[i]  * steps[i];

     /*print ("der:%5.2f %5.2f %5.2f def:%5.2f %5.2f %5.2f (%7.3f %7.3f  %7.3f)\n", 
       dx[0],dy[0],dz[0],
       def_vector[0],def_vector[1], def_vector[2],
       dist_weight,thresh,val[0]);
     */
     result = sqrt (def_vector[0]*def_vector[0] +
                    def_vector[1]*def_vector[1] +
                    def_vector[2]*def_vector[2]);
  }
  return(result);
}

/* procedure build_lattices will construct the sub-lattice on the
   source and target volumes for xcorr, diff and label objective
   functions that need to be optimized with the similarity_function 

   inputs:
   outputs:

*/

private BOOLEAN build_lattices(Real spacing, 
			       Real threshold, 
			       Real source_coord[],
			       Real mean_target[],
			       Real target_coord[],
			       Real def_vector[],
			       int ndim)
{

  BOOLEAN
    result;
  Real
    pos[3],
    xp,yp,zp;
  int 
    i,j;

  /* if there is no gradient magnitude strong enough to grab onto,
     then set a negative magnitude deformation and skip the optimization 
     
     otherwise
     
     Bias the initial starting search position for the local
     deformation vector by calculating a target position that is equal
     to the vector average of 
     (current_transformation(source) + average_position(source in target)
     
  */
  

  result = TRUE;

  if (!get_best_start_from_neighbours(threshold,
				      source_coord, mean_target, target_coord,
				      def_vector)) {
    
    for_less (i,0,3)
      target_coord[i] = mean_target[i];

    result = FALSE;
    
  }
  else {   

    /* we now have the info needed to continue...
       
       source_coord[] - point in source volume
       target_coord[] - best point in target volume, so far.
       def_vector[]   - currently contains the additional def needed to 
                        take source_coord mid-way to the neighbour's
                        mean-point (which is now stored in target_coord[]).  

       now,
       get the world coord position of the node in the source volume,
       corresponding to the target node in question, taking into
       consideration the current warp  
       
       xp,yp,zp will be the position in the source volume, corresponding to
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
          sub-lattice within the source volume                       

       The spherical sub-lattice will have Glen points in the source
       volume, note: sub-lattice diameter= 1.5*fwhm 
       (here specified as 3*spacing = 2*(fwhm/2) to specify radius 
       in build_source_lattice) 

       note that SX, SY, SZ, TX,TY,TZ, and Glen are all globals
    */

    build_source_lattice(xp, yp, zp, 
			 SX, SY, SZ,
			 spacing*3*2, spacing*3*2, spacing*3*2,
			 Diameter_of_local_lattice,  
			 Diameter_of_local_lattice,  
			 Diameter_of_local_lattice,
			 ndim, &Glen);

    /* -------------------------------------------------------------- */
    /* BUILD THE TARGET VOLUME LOCAL NEIGHBOURHOOD INFO */

    /* map this lattice forward into the target space, using the
       current transformation, in order to build a deformed lattice
       (in the WORLD COORDS of the target volume) */

    if (Gglobals->trans_info.use_super>0) 
      build_target_lattice_using_super_sampled_def(
		  SX,SY,SZ, TX,TY,TZ, Glen, ndim);
    else 
      build_target_lattice(SX,SY,SZ, TX,TY,TZ, Glen, ndim);
      

    /* -------------------------------------------------------------- */
    /* GET THE VOXEL COORDINATE LIST:
	  need the voxel coordinates of the target lattice
          for the objective function used in the actual 
	  optimization,

	  (world coords are use for the source, since it is only 
	   interpolated once)

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
       that will be used in the optimization below                    */

    if (Gglobals->trans_info.use_magnitude) {
      for_inclusive(i,1,Glen) {
	SX[i] += source_coord[X] - xp;
	SY[i] += source_coord[Y] - yp;
	SZ[i] += source_coord[Z] - zp;
      }
    }

    /* -------------------------------------------------------------- */
    /* GO GET FEATURES IN SOURCE VOLUME actually get the feature data from
       the source volume local neighbourhood and compute the required
       similarity function.  The target volume features are retrieved and
       the similarity function computed in go_get_samples_with_offset() 

       note that I only have to get features from the volumes that
       will use the sublattice in the optimization 
    */

    for_less(i,0, Gglobals->features.number_of_features) {

      if (Gglobals->features.obj_func[i] != NONLIN_OPTICALFLOW && Gglobals->features.obj_func[i] != NONLIN_CHAMFER)

	go_get_samples_in_source(Gglobals->features.data[i], 
				 SX,SY,SZ, Ga1_features[i], Glen, 
				 (Gglobals->interpolant==nearest_neighbour_interpolant ? -1 : 0)
				 );
    }

    /* -------------------------------------------------------------- */
    /* calc one of the normalization coefficients for the similarity
       measure, when using the magnitude data.  This saves a few CPU cycles
       in go_get_samples_with_offset(), since the constants only have to be
       eval'd once for the source volume. Note that this variable is not
       used when doing OPTICAL FLOW. */

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
      case NONLIN_CHAMFER:
	Gsqrt_features[i] = 0;
	break;
      case NONLIN_OPTICALFLOW:
	Gsqrt_features[i] = 0;
	break;
      default:
	print_error_and_line_num("Objective function %d not supported in build_lattices",
                                 __FILE__, __LINE__,Gglobals->features.obj_func[i]);
      }
    }
    
  }
  return(result );
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
					     int ndim,
					     BOOLEAN sub_lattice_needed)
{

  Real
    real_def[3], vox_def[3],
    temp_total_weight,
    optical_partial_weight,
    other_partial_weight,
    total_weight,
    du,dv,dw,
    local_corr3D[3][3][3],
    local_corr2D[3][3],
    optical_def_vector[3],
    optical_voxel_displacement[3],
    voxel[3],
    pos[3],
    simplex_size,
    result,
    target_coord[3];
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

				/* initialize for no deformation */
  result = 0.0;			
  def_vector[X] = def_vector[Y] = def_vector[Z] = 0.0;
  *num_functions = 0;		

				/* build sub-lattice if necessary */
  if (sub_lattice_needed) {

    if ( ! build_lattices(spacing, threshold1, 
			  source_coord, mean_target, target_coord, def_vector,
			  ndim) ){
      result = -DBL_MAX;
      return(result);		/* return if we don't make the threshold */
    }

  }
				/* compute weighting factors */

  optical_partial_weight = other_partial_weight = total_weight = 0.0;

  for_less(i,0,Gglobals->features.number_of_features) {

    if ((Gglobals->features.obj_func[i] == NONLIN_OPTICALFLOW) || 
        (Gglobals->features.obj_func[i] == NONLIN_CHAMFER) )
      optical_partial_weight += Gglobals->features.weight[i];
    else
      other_partial_weight += Gglobals->features.weight[i];

    total_weight += Gglobals->features.weight[i];
  }

  if (total_weight == 0.0) {
    print_error_and_line_num("Objective functions have no total weight in get_deformation_vector_for_node",
			     __FILE__, __LINE__);
  }

  
  /* estimate deformations using optimization over the sub-lattice if
     required */

  if (other_partial_weight > 0.0) {
    /* -------------------------------------------------------------- */
    /*  FIND BEST DEFORMATION VECTOR
	now find the best local deformation that maximises the local
	neighbourhood correlation between the source values stored in
	**Ga1_features at positions Sx, SY, SZ with the homologous 
	values at positions TX,TY,TZ in the target volume */
    
    if ( !Gglobals->trans_info.use_simplex) {
      
      /* ----------------------------------------------------------- */
      /*  USE QUADRATIC FITTING to find best deformation vector      */
      
      if (ndim==3) { /* build up the 3x3x3 matrix of local correlation values */
	
	for_inclusive(i,-1,1) {
	  
	  pos_vector[1] = (float) i * Gsimplex_size/2.0;
	  for_inclusive(j,-1,1) {
	    
	    pos_vector[2] = (float) j * Gsimplex_size/2.0;
	    for_inclusive(k,-1,1) {
	      pos_vector[3] = (float) k * Gsimplex_size/2.0;
	      local_corr3D[i+1][j+1][k+1] = local_objective_function(pos_vector); 
	    }
	  }
	}
	*num_functions += 27;
	flag = return_3D_disp_from_min_quad_fit(local_corr3D, &du, &dv, &dw);
	
      }
      else {
	/* build up the 3x3 matrix of local correlation values */
	
	pos_vector[3] = 0.0;	/* since 2D */
	
	for_inclusive(i,-1,1) {
	  pos_vector[1] = (float) i * Gsimplex_size/2.0;
	  for_inclusive(j,-1,1) {
	    pos_vector[2] = (float) j * Gsimplex_size/2.0;
	    local_corr2D[i+1][j+1] = 1.0 - local_objective_function(pos_vector); 
	  }
	}
	*num_functions += 9;
	
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
			simplex_size, amoeba_NL_obj_function, 
			NULL, (Real)ftol);
      
      
      nfunk = 4;		/* since 4 eval's needed to init the amoeba */
      
      /*   do the actual SIMPLEX optimization,
	   note that nfunk is incremented inside perform_amoeba  */

      while (nfunk < AMOEBA_ITERATION_LIMIT && 
	     perform_amoeba(&the_amoeba, &nfunk) );
      
      *num_functions += nfunk;

      if (nfunk < AMOEBA_ITERATION_LIMIT) {
	
	get_amoeba_parameters(&the_amoeba,parameters);
	
	
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
	
      } /*  if perform_amoeba */
      
      terminate_amoeba(&the_amoeba);      
      FREE(parameters);
      
    } /* else use_simplex */
    
  
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
      
      for_less (j,0,3) {	/* weight these displacements properly */
	def_vector[j]         *= other_partial_weight / total_weight;
	voxel_displacement[j] *= other_partial_weight / total_weight;
      }

    }

  }
  
  /* at this point, we have the deformations coming from optimization
     of correlation, label or difference objective functions (if any)
     stored in def_vector and voxel_displacement, now time to find
     displacements with OPTICAL FLOW if needed */

  if (optical_partial_weight > 0.0) {

    for_less(i,0,3) {		/* init optical/chamfer to zero */
      optical_def_vector[i] = 0.0;
      optical_voxel_displacement[i] = 0.0;
    }

    temp_total_weight = 0;

    for_less(i,0,Gglobals->features.number_of_features) {
      
      if (Gglobals->features.obj_func[i] == NONLIN_OPTICALFLOW ||  
          Gglobals->features.obj_func[i] == NONLIN_CHAMFER)  {
	
	if (Gglobals->features.obj_func[i] == NONLIN_OPTICALFLOW)
	  result =  get_optical_flow_vector(threshold1, 
					    source_coord, mean_target,
					    real_def, vox_def,
					    Gglobals->features.data[i],
					    Gglobals->features.model[i],
					    ndim);
	else                   /* must be CHAMFER */
	  result =  get_chamfer_vector(spacing,   
				       source_coord, mean_target,
				       real_def, vox_def,
				       Gglobals->features.data[i],
				       Gglobals->features.model[i],
				       ndim);
	if (result > 0.0) {
	  *num_functions += 1;
                                /* add in the weighted deformations */

	  temp_total_weight += Gglobals->features.weight[i];
	  
	  for_less (j,0,3) {
	    optical_def_vector[j]         += real_def[j] * Gglobals->features.weight[i];
	    optical_voxel_displacement[j] += vox_def[j]  * Gglobals->features.weight[i];
	  }
	} 

      }

    }
    				/* add in the weighted defs from optical/chamfer */
    if (temp_total_weight > 0.0) {
      for_less (j,0,3) {
	def_vector[j]         += optical_def_vector[j] / temp_total_weight;
	voxel_displacement[j] += optical_voxel_displacement[j] / temp_total_weight;
      }
    }

  }

  result = sqrt((def_vector[X] * def_vector[X]) + 
		(def_vector[Y] * def_vector[Y]) + 
		(def_vector[Z] * def_vector[Z])) ;      

  return(result);
}


