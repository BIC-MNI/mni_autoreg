/* ----------------------------- MNI Header -----------------------------------
@NAME       : optimize.c
@DESCRIPTION: collection of routines for user-specified optimization.
@METHOD     : now, only simplex method is used.
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

@MODIFIED   : $Log: optimize.c,v $
@MODIFIED   : Revision 96.5  2000-05-05 17:57:04  louis
@MODIFIED   : addes volume intensity normalization code
@MODIFIED   :
@MODIFIED   : Revision 96.4  2000/02/16 22:09:36  stever
@MODIFIED   : fine tuning basic nonlinear tests
@MODIFIED   :
@MODIFIED   : Revision 96.3  2000/02/15 19:02:08  stever
@MODIFIED   : Add tests for param2xfm, minctracc -linear.
@MODIFIED   :
@MODIFIED   : Revision 96.2  1997/11/12 21:07:43  louis
@MODIFIED   : no changes, other than rcsid...
@MODIFIED   :
 * Revision 96.1  1997/11/03  15:06:29  louis
 * working version, before creation of mni_animal package, and before inserting
 * distance transforms
 *
 * Revision 96.1  1997/11/03  15:06:29  louis
 * working version, before creation of mni_animal package, and before inserting
 * distance transforms
 *
 * Revision 96.0  1996/08/21  18:22:10  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:05  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:03  louis
 * Never released version 0.95
 *
 * Revision 1.17  1996/08/12  14:15:57  louis
 * Pre-release
 *
 * Revision 1.16  1996/05/02  19:39:12  collins
 * fixed the matrix inversion bug.  When source and target volumes were
 * swapped, the call to fit_function, to get a value for final_corr,
 * rebuilt the inverse registration matrix.  And this, just before the
 * matrix was to be written out to disk.
 *
 * Now, the matrix is correctly computed, after the last call to
 * fit_function.
 *
 * Revision 1.15  1996/03/25  10:33:15  collins
 * modifications apply to the implementation of mutual information
 * and the constraints of using byte-only data and 256 groups.
 *
 * also, I now call  get_volume_data_type() to get the volume data type.
 *
 * also, calls to set initial_corr and final_corr are done in
 * optimize_linear_transformation() instead of calling xcorr_fitting_fn
 * in main(), so that the user-selected obj-fn is used to measure
 * the before and after fit.
 *
 * Revision 1.14  1996/03/07  13:25:19  collins
 * small reorganisation of procedures and working version of non-isotropic
 * smoothing.
 *
 * Revision 1.13  1995/10/06  09:25:02  collins
 * removed all the volume parameters from do_non_linear_optimization() and
 * optimize_non_linear_transformation() since they are all represented in
 * globals->features struc.
 *
 * made the necesary changes to access globals->features instead of previous
 * ly used volumes in input parameters.
 *
 * Revision 1.12  1995/09/07  13:02:13  collins
 * removed numerical recipes call to amoeba, the simplex optimization
 * procedure.  It has been replaced with Davids version of the simplex
 * optimization code.
 *
 * Revision 1.11  1995/09/07  10:05:11  collins
 * All references to numerical recipes routines are being removed.  At this
 * stage, any num rec routine should be local in the file.  All memory
 * allocation calls to vector(), matrix(), free_vector() etc... have been
 * replaced with ALLOC and FREE from the volume_io library of routines.
 *
 * Revision 1.10  1995/03/17  11:15:35  collins
 * added prototype for amoeba2 in optimize.c - I think that Greg may
 * have changed the call from amoeba2 to amoeba so that he could take
 * advantage of the amoeba routine in librecipes.  I don't think he
 * knew that there are small differences between the original NR amoeba
 * and the version I use.  So I put it back to the way it was, or in
 * any case, to the way that works.
 *
 * Revision 1.9  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.8  94/04/26  12:54:33  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.7  94/04/06  11:48:44  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.6  94/02/21  16:35:59  louis
 * version before feb 22 changes
 * 
 * Revision 1.5  93/11/15  16:27:08  louis
 * working version, with new library, with RCS revision stuff,
 * before deformations included
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/optimize.c,v 96.5 2000-05-05 17:57:04 louis Exp $";
#endif

#include <config.h>
#include <internal_volume_io.h>
#include <print_error.h>
#include <amoeba.h>

#include "constants.h"
#include "arg_data.h"
#include "objectives.h"
#include "make_rots.h"
#include "segment_table.h"

extern Arg_Data main_args;

Volume   Gdata1, Gdata2, Gmask1, Gmask2;
int      Ginverse_mapping_flag, Gndim;

extern   double   ftol ;        
extern   double   simplex_size ;
extern   Real     initial_corr, final_corr;

         Segment_Table  *segment_table;	/* for variance of ratios */

         Real            **prob_hash_table;   /* for mutual information */
         Real            *prob_fn1;      /*     for vol 1 */
         Real            *prob_fn2;      /*     for vol 2 */


/* external calls: */

public  BOOLEAN  perform_amoeba(amoeba_struct  *amoeba, int *num_funks );
public  void  initialize_amoeba(
    amoeba_struct     *amoeba,
    int               n_parameters,
    Real              initial_parameters[],
    Real              parameter_delta,
    amoeba_function   function,
    void              *function_data,
    Real              tolerance );

public  Real  get_amoeba_parameters(
    amoeba_struct  *amoeba,
    Real           parameters[] );

public  void  terminate_amoeba(
    amoeba_struct  *amoeba );


public void make_zscore_volume(Volume d1, Volume m1, 
			       Real *threshold); 

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double  *start, int *count, 
				  VectorR directions[]);

public Status do_non_linear_optimization(Arg_Data *globals);


public void parameters_to_vector(double *trans, 
				 double *rots, 
				 double *scales,
				 double *shears,
				 float  *op_vector,
				 double *weights) 
{
  int i;

  i = 1;

  if (weights[0]  != 0.0) { op_vector[i] = trans[0] /weights[0]; ++i; }
  if (weights[1]  != 0.0) { op_vector[i] = trans[1] /weights[1]; ++i; }
  if (weights[2]  != 0.0) { op_vector[i] = trans[2] /weights[2]; ++i; }

  if (weights[3]  != 0.0) { op_vector[i] = rots[0] / weights[3]; ++i; }
  if (weights[4]  != 0.0) { op_vector[i] = rots[1] / weights[4]; ++i; }
  if (weights[5]  != 0.0) { op_vector[i] = rots[2] / weights[5]; ++i; }

  if (weights[6]  != 0.0) { op_vector[i] = scales[0]/weights[6]; ++i; }
  if (weights[7]  != 0.0) { op_vector[i] = scales[1]/weights[7]; ++i; }
  if (weights[8]  != 0.0) { op_vector[i] = scales[2]/weights[8]; ++i; }

  if (weights[9]  != 0.0) { op_vector[i] = shears[0]/weights[9]; ++i; }
  if (weights[10] != 0.0) { op_vector[i] = shears[1]/weights[10];++i; }
  if (weights[11] != 0.0) { op_vector[i] = shears[2]/weights[11];     }


}

private void vector_to_parameters(double *trans, 
				  double *rots, 
				  double *scales,
				  double *shears,
				  float  *op_vector,
				  double *weights) 
{
  int i;

  i = 1;

  if (weights[0]  != 0.0) { trans[0] = op_vector[i] * weights[0]; ++i; }
  if (weights[1]  != 0.0) { trans[1] = op_vector[i] * weights[1]; ++i; }
  if (weights[2]  != 0.0) { trans[2] = op_vector[i] * weights[2]; ++i; }

  if (weights[3]  != 0.0) { rots[0]  = op_vector[i] * weights[3]; ++i; }
  if (weights[4]  != 0.0) { rots[1]  = op_vector[i] * weights[4]; ++i; }
  if (weights[5]  != 0.0) { rots[2]  = op_vector[i] * weights[5]; ++i; }

  if (weights[6]  != 0.0) { scales[0]= op_vector[i] * weights[6]; ++i; }
  if (weights[7]  != 0.0) { scales[1]= op_vector[i] * weights[7]; ++i; }
  if (weights[8]  != 0.0) { scales[2]= op_vector[i] * weights[8]; ++i; }

  if (weights[9]  != 0.0) { shears[0]= op_vector[i] * weights[9]; ++i; }
  if (weights[10] != 0.0) { shears[1]= op_vector[i] * weights[10];++i; }
  if (weights[11] != 0.0) { shears[2]= op_vector[i] * weights[11];     }
}

inline private BOOLEAN in_limits(double x,double lower,double upper)
{
    return lower <= x && x <= upper;
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : fit_function
@INPUT      : params - a variable length array of floats
@OUTPUT     :               
@RETURNS    : a float value of the user requested objective function,
              measuring the similarity between two data sets.
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : 
@MODIFIED   : 
---------------------------------------------------------------------------- */

public float fit_function(float *params) 
{

  Transform *mat;
  int i;
  float r;


  double trans[3];
  double cent[3];
  double rots[3];
  double scale[3];
  double shear[6];


  for_less( i, 0, 3 ) {		/* set default values from GLOBAL MAIN_ARGS */
    shear[i] = main_args.trans_info.shears[i];
    scale[i] = main_args.trans_info.scales[i];
    trans[i] = main_args.trans_info.translations[i];
    rots[i]  = main_args.trans_info.rotations[i];
    cent[i]  = main_args.trans_info.center[i];
  }


				/* modify the parameters to be optimized */
  vector_to_parameters(trans, rots, scale, shear, params, main_args.trans_info.weights);
  
  if (main_args.trans_info.transform_type==TRANS_LSQ7) { /* adjust scaley and scalez only */
                                                         /* if 7 parameter fit.  */
    scale[1] = scale[0];
    scale[2] = scale[0];
  }

  if (!in_limits(rots[0], (double)-3.1415927/2.0, (double)3.1415927/2.0) ||
      !in_limits(rots[1], (double)-3.1415927/2.0, (double)3.1415927/2.0) ||
      !in_limits(rots[2], (double)-3.1415927/2.0, (double)3.1415927/2.0) ||
      !in_limits(scale[0],(double)0.0, (double)3.0) ||
      !in_limits(scale[1],(double)0.0, (double)3.0) ||
      !in_limits(scale[2],(double)0.0, (double)3.0) ||
      !in_limits(shear[0],(double)-2.0, (double)2.0) ||
      !in_limits(shear[1],(double)-2.0, (double)2.0) ||
      !in_limits(shear[2],(double)-2.0, (double)2.0))

    {

 (void)printf("out : %7.4f=%c %7.4f=%c %7.4f=%c   %7.4f=%c %7.4f=%c %7.4f=%c   %7.4f=%c %7.4f=%c %7.4f=%c \n",
       rots[0], in_limits(rots[0], (double)-3.1415927/2.0, (double)3.1415927/2.0) ? 'T': 'F' , 
       rots[1], in_limits(rots[1], (double)-3.1415927/2.0, (double)3.1415927/2.0) ? 'T': 'F' , 
       rots[2], in_limits(rots[2], (double)-3.1415927/2.0, (double)3.1415927/2.0) ? 'T': 'F' , 
       scale[0],in_limits(scale[0], (double)0.0, (double)3.0)? 'T': 'F' , 
       scale[1],in_limits(scale[1], (double)0.0, (double)3.0)? 'T': 'F' , 
       scale[2],in_limits(scale[2], (double)0.0, (double)3.0)? 'T': 'F' , 
       shear[0],in_limits(shear[0], (double)-2.0, (double)2.0)? 'T': 'F' , 
       shear[1],in_limits(shear[1], (double)-2.0, (double)2.0)? 'T': 'F' , 
       shear[2],in_limits(shear[2], (double)-2.0, (double)2.0)? 'T': 'F' );

    r = 1e10;
  }
  else {
				/* get the linear transformation ptr */

    if (get_transform_type(main_args.trans_info.transformation) == CONCATENATED_TRANSFORM) {
      mat = get_linear_transform_ptr(
             get_nth_general_transform(main_args.trans_info.transformation,0));
    }
    else
      mat = get_linear_transform_ptr(main_args.trans_info.transformation);
    
    if (Ginverse_mapping_flag)
      build_inverse_transformation_matrix(mat, cent, trans, scale, shear, rots);
    else
      build_transformation_matrix(mat, cent, trans, scale, shear, rots);
    
    /* call the needed objective function */
    
    r = (main_args.obj_function)(Gdata1,Gdata2,Gmask1,Gmask2,&main_args);
  }

  return(r);
}


public Real amoeba_obj_function(void *dummy, float d[])
{
  int i;
  float p[13];

  for_less(i,0,Gndim)
    p[i+1] = d[i];
  
  return ( (Real)fit_function(p) );
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : optimize_simplex
                get the parameters necessary to map volume 1 to volume 2
		using the simplex optimizaqtion algorithm and a user specified
		objective function.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
	      m1,m2:
                two mask volumes for data (already in memory).
	      globals:
	        a global data structure containing info from the command line,
		including the input parameters to be optimized, the input matrix,
		and a plethora of flags!
@OUTPUT     : 
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: 
@METHOD     : uses the simplex algorithm extracted from BICPL
@GLOBALS    : 
@CALLS      : 
@CREATED    : Fri Jun 11 11:16:25 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
public BOOLEAN optimize_simplex(Volume d1,
				Volume d2,
				Volume m1,
				Volume m2, 
				Arg_Data *globals)
{
  BOOLEAN 
    stat;
  float 
    local_ftol,
    *p;
  amoeba_struct 
    the_amoeba;
  int 
    iteration_number,
    max_iters,
    i,j, 
    ndim;

  Transform
    *mat;

  Real
    *parameters;

  double trans[3];
  double cent[3];
  double rots[3];
  double scale[3];
  double shear[6];

  
  stat = TRUE;
  local_ftol = ftol;
				/* find number of dimensions for optimization */
  ndim = 0;
  for_less(i,0,12)
    if (globals->trans_info.weights[i] != 0.0) ndim++;

				/* set GLOBALS to communicate with the
				   function to be fitted!              */
  if (stat && ndim>0) {
    Gndim = ndim;

    ALLOC(p,ndim+1+1);		/* my parameters for the simplex 
				   [1..ndim+1]*/

    ALLOC(parameters, ndim+1);	/* David's parmaters for the simplex
				   [0..ndim] */
    
				/* build the parameter vector from the 
				   initial transformation parameters   */
    parameters_to_vector(globals->trans_info.translations,
			 globals->trans_info.rotations,
			 globals->trans_info.scales,
			 globals->trans_info.shears,
			 p,
			 globals->trans_info.weights);

    for_less(i,0,ndim+1)		/* copy initial guess into parameter list */
      parameters[i] = (Real)p[i+1];

    initialize_amoeba(&the_amoeba, ndim, parameters, 
		      simplex_size, amoeba_obj_function, 
		      NULL, (Real)local_ftol);

    max_iters = 400;
    iteration_number = 0;
				/* do the ameoba optimization */
    while ( iteration_number<max_iters && perform_amoeba(&the_amoeba, &iteration_number) ) 
      /* empty */ ;

    
    if (globals->flags.debug) {
      
      (void)print("done with simplex after %d iterations\n",iteration_number);
      for_less(i,0,the_amoeba.n_parameters+1) {
	
	(void)print ("%d %7.5f:",i,the_amoeba.values[i]);
	for_less(j,0,the_amoeba.n_parameters) {
	  (void)print ("%8.5f ", the_amoeba.parameters[i][j]);
	}
	(void)print ("\n");
	
      }
    }

				/* copy result into main data structure */
    get_amoeba_parameters(&the_amoeba,parameters);
    for_less(i,0,ndim+1)		
      p[i+1] = (float)parameters[i];
    
    vector_to_parameters(globals->trans_info.translations,
			 globals->trans_info.rotations,
			 globals->trans_info.scales,
			 globals->trans_info.shears,
			 p,
			 globals->trans_info.weights);
    terminate_amoeba(&the_amoeba);

    if (globals->trans_info.transform_type==TRANS_LSQ7) { /* adjust scaley and scalez only */
      /* if 7 parameter fit.  */
      globals->trans_info.scales[1] = globals->trans_info.scales[0];
      globals->trans_info.scales[2] = globals->trans_info.scales[0];
    }
    
    for_less( i, 0, 3 ) {		/* set translations */
      trans[i] = globals->trans_info.translations[i]; 
      rots[i]  = globals->trans_info.rotations[i];
      scale[i] = globals->trans_info.scales[i];
      cent[i]  = globals->trans_info.center[i];
      shear[i] = globals->trans_info.shears[i];
    }


    if (globals->flags.debug) {
      print("after parameter optimization\n");
      print("-center      %10.5f %10.5f %10.5f\n", cent[0], cent[1], cent[2]);
      print("-translation %10.5f %10.5f %10.5f\n", trans[0], trans[1], trans[2]);
      print("-rotation    %10.5f %10.5f %10.5f\n", 
	    rots[0]*180.0/3.1415927, rots[1]*180.0/3.1415927, rots[2]*180.0/3.1415927);
      print("-scale       %10.5f %10.5f %10.5f\n", scale[0], scale[1], scale[2]);
      print("-shear       %10.5f %10.5f %10.5f\n", shear[0], shear[1], shear[2]);
    }
  
    if (get_transform_type(globals->trans_info.transformation) == CONCATENATED_TRANSFORM) {
      mat = get_linear_transform_ptr(
              get_nth_general_transform(globals->trans_info.transformation,0));
    }
    else
      mat = get_linear_transform_ptr(globals->trans_info.transformation);
    
    build_transformation_matrix(mat, cent, trans, scale, shear, rots);

    FREE(p);
    FREE(parameters);
  }

  return( stat );
}


public BOOLEAN replace_volume_data_with_ubyte(Volume data)
{
  Volume tmp_vol;
  int sizes[MAX_DIMENSIONS];
  int i,j,k,count,n_dim;
  progress_struct		
    progress;

  n_dim = get_volume_n_dimensions(data);
  get_volume_sizes(data, sizes);

  if (n_dim != 3) {
    print ("Volume must have 3 dimensions for byte copy\n");
    return(FALSE);
  }
				/* build a matching temporary ubyte 
				   volume */

  tmp_vol = copy_volume_definition(data, NC_BYTE, FALSE, 0.0, 0.0);

				/* copy the original voxel data into
				   the byte voxels */
  count = 0;
  initialize_progress_report(&progress, FALSE, sizes[0]*sizes[1]*sizes[2] + 1,
			     "Converting" );
  for_less(i,0,sizes[0])
    for_less(j,0,sizes[1])
      for_less(k,0,sizes[2]) {
	set_volume_real_value(tmp_vol, i,j,k,0,0,
			      get_volume_real_value(data,i,j,k,0,0));
	count++;
	update_progress_report( &progress, count);
      }
  terminate_progress_report( &progress );

  
  free_volume_data( data );	/* get rid of original data */

   /* 
   * BLEAGHH!! This is an evil and nasty hack, made worse by the fact that
   * we have to support two versions of Volume_io for it to work!  (At
   * least this is done by the VOXEL_DATA hack^H^H^H^Hmacro, in <config.h>.)
   * -GPW 96/04/34
   */

  VOXEL_DATA (data) = VOXEL_DATA (tmp_vol);
  
  set_volume_type(data, NC_BYTE, FALSE, 0.0, 0.0);

  VOXEL_DATA (tmp_vol) = NULL;  /* is this really necessary?!?!? */


  return(TRUE);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : optimize_linear_transformation
                get the parameters necessary to map volume 1 to volume 2
		using a user specified optimization strategy and objective
		function.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
	      m1,m2:
                two mask volumes for data (already in memory).
	      globals:
	        a global data structure containing info from the command line,
		including the input parameters to be optimized, the input matrix,
		and a plethora of flags!
@OUTPUT     : 
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: 
@METHOD     :
                1- this routine begins by initializing the volume data structures
		to be used by the objective functions.

		2- optimization function is called
		
		3- the optimized parameters are returned in globals...
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun  9 12:56:08 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
public BOOLEAN optimize_linear_transformation(Volume d1,
					      Volume d2,
					      Volume m1,
					      Volume m2, 
					      Arg_Data *globals)
{
  BOOLEAN 
    stat;
  int i;
  Data_types
    data_type;
  float *p;
  Transform
    *mat;

  double trans[3];
  double cent[3];
  double rots[3];
  double scale[3];
  double shear[6];


  stat = TRUE;

          /* --------------------------------------------------------------*/
          /*----------------- prepare data for optimization -------------- */
  
  if (globals->obj_function == zscore_objective) 
				/* normalize volumes before correlation */
    { 
      /* replace volume d1 and d2 by zscore volume  */

      make_zscore_volume(d1,m1,&globals->threshold[0]);
      make_zscore_volume(d2,m2,&globals->threshold[1]);
    } else
  if (globals->obj_function == ssc_objective)
				/* Stocastic sign change (or zero-crossings) */
    {
      /* add speckle to the data set, after making both data sets
         comparable in mean and sd...                             */

      make_zscore_volume(d1,m1,&globals->threshold[0]); 
      make_zscore_volume(d2,m2,&globals->threshold[1]); 

      if (globals->smallest_vol == 1)
	add_speckle_to_volume(d1, 
			      globals->speckle,
			      globals->start, globals->count, 
			      globals->directions);
      else
	add_speckle_to_volume(d2, 
			      globals->speckle,
			      globals->start, globals->count, 
			      globals->directions);    
    } else
  if (globals->obj_function == vr_objective)
				/* Woods' variance of ratios */
    {

      if (globals->smallest_vol == 1) {
	if (!build_segment_table(&segment_table, d1, globals->groups))
	  print_error_and_line_num(
	    "Could not build segment table for SOURCE volume\n",
	    __FILE__, __LINE__);
      }
      else {
	if (!build_segment_table(&segment_table, d2, globals->groups))
	  print_error_and_line_num(
	    "Could not build segment table for TARGET volume\n",
	    __FILE__, __LINE__);	
      }
      
      if (globals->flags.debug && globals->flags.verbose>1) {
	print ("groups = %d\n",segment_table->groups);
	for_less(i, segment_table->min, segment_table->max+1) {
	  print ("%5d: table = %5d, function = %5d\n",i,segment_table->table[i],
		 (segment_table->segment)(i,segment_table) );
	}
      }
    } else
  if (globals->obj_function == mutual_information_objective)
				/* Collignon's mutual information */
    {

      if ( globals->groups != 256 ) {
	print ("WARNING: -groups was %d, but will be forced to 256 in this run\n",globals->groups);
	globals->groups = 256;
      }

      data_type = get_volume_data_type (d1);
      if (data_type != UNSIGNED_BYTE) {
	print ("WARNING: source volume not UNSIGNED_BYTE, will do conversion now.\n");
	if (!replace_volume_data_with_ubyte(d1)) {
	  print_error_and_line_num("Can't replace volume data with unsigned bytes\n",
			     __FILE__, __LINE__);
	}
      }

      data_type = get_volume_data_type (d2);
      if (data_type != UNSIGNED_BYTE) {
	print ("WARNING: target volume not UNSIGNED_BYTE, will do conversion now.\n");
	if (!replace_volume_data_with_ubyte(d2)) {
	  print_error_and_line_num("Can't replace volume data with unsigned bytes\n",
			     __FILE__, __LINE__);
	}
      }

      ALLOC(   prob_fn1,   globals->groups);
      ALLOC(   prob_fn2,   globals->groups);
      ALLOC2D( prob_hash_table, globals->groups, globals->groups);

    } else
  if (globals->obj_function == xcorr_objective) {
				/*EMPTY*/
  }
  else  {
    print_error_and_line_num("Unknown objective function value\n",
			     __FILE__, __LINE__);
  }    



  /* --------------------------------------------------------------*/
  /* -------- prepare the weighting array for optimization --------*/
  
  switch (globals->trans_info.transform_type) {
  case TRANS_LSQ3: 
    for_less(i,3,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
      globals->trans_info.scales[i] = 1.0;
      globals->trans_info.rotations[i] = 0.0;
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ6: 
    for_less(i,6,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
      globals->trans_info.scales[i] = 1.0;
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ7: 
    for_less(i,7,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ9: 
    for_less(i,9,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ10: 
    for_less(i,10,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,1,3) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ: 
				/* nothing to be zeroed */
    break;
  case TRANS_LSQ12: 
				/* nothing to be zeroed */
    break;
  default:
    (void)fprintf(stderr, "Unknown type of transformation requested (%d)\n",
		   globals->trans_info.transform_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    stat = FALSE;
  }

	   /* ---------------- swap the volumes, so that the smallest
	                       is first (to save on CPU)             ---------*/

  if (globals->smallest_vol == 1) {
    Gdata1 = d1;      Gdata2 = d2;
    Gmask1 = m1;      Gmask2 = m2;
    Ginverse_mapping_flag = FALSE;
  }
  else {
    Gdata1 = d2;      Gdata2 = d1;
    Gmask1 = m2;      Gmask2 = m1;
    Ginverse_mapping_flag = TRUE;
  }


	   /* ---------------- call the requested obj_function to 
	                       establish the initial fitting value  ---------*/

  ALLOC(p,13);
  parameters_to_vector(globals->trans_info.translations,
		       globals->trans_info.rotations,
		       globals->trans_info.scales,
		       globals->trans_info.shears,
		       p,
		       globals->trans_info.weights);

  initial_corr = fit_function(p);

	   /* ---------------- call requested optimization strategy ---------*/

  switch (globals->optimize_type) {
  case OPT_SIMPLEX:
    stat = optimize_simplex(d1, d2, m1, m2, globals);
    break;
  default:
    (void)fprintf(stderr, "Unknown type of optimization requested (%d)\n",
		  globals->optimize_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    stat = FALSE;
  }
  
  parameters_to_vector(globals->trans_info.translations,
		       globals->trans_info.rotations,
		       globals->trans_info.scales,
		       globals->trans_info.shears,
		       p,
		       globals->trans_info.weights);

  final_corr = fit_function(p);

  FREE(p);

  /*--------- set up final transformation matrix ------------------*/

  if (get_transform_type(globals->trans_info.transformation) == CONCATENATED_TRANSFORM) {
    mat = get_linear_transform_ptr(
	   get_nth_general_transform(globals->trans_info.transformation,0));
  }
  else
    mat = get_linear_transform_ptr(globals->trans_info.transformation);
  
  for_less( i, 0, 3 ) {		/* set translations */
    trans[i] = globals->trans_info.translations[i]; 
    rots[i]  = globals->trans_info.rotations[i];
    scale[i] = globals->trans_info.scales[i];
    cent[i]  = globals->trans_info.center[i];
    shear[i] = globals->trans_info.shears[i];
  }

  build_transformation_matrix(mat, cent, trans, scale, shear, rots);


          /* ----------------finish up parameter/matrix manipulations ------*/

  if (globals->obj_function == vr_objective)
    {
      stat = stat && free_segment_table(segment_table);
    } else
  if (globals->obj_function == mutual_information_objective)
				/* Collignon's mutual information */
    {
      FREE(   prob_fn1 );
      FREE(   prob_fn2 );
      FREE2D( prob_hash_table);
    }


  return(stat);
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : measure_fit
                evaluate the objective function comparing d1 to d2 a single
		time.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
	      m1,m2:
                two mask volumes for data (already in memory).
	      globals:
	        a global data structure containing info from the command line,
		including the input parameters, the input matrix,
		and a plethora of flags!
@OUTPUT     : 
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: 
@METHOD     :
                1- this routine begins by initializing the volume data structures
		to be used by the objective functions.

		2- the objective function is evaluated
		
		3- its value is returned
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Jul 13 11:15:57 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
public float measure_fit(Volume d1,
			 Volume d2,
			 Volume m1,
			 Volume m2, 
			 Arg_Data *globals)
{
  BOOLEAN 
    stat;
  float 
    **p, y;
  int 
    i, 
    ndim;
  Data_types
    data_type;

/*  Transform
    *mat;
    double trans[3];  double cent[3];  double rots[3];  double scale[3];  double shear[6]; */

  
  stat = TRUE;

	     /*----------------- prepare data for objective function evaluation ------------ */

  
  if (globals->obj_function == zscore_objective) { /* replace volume d1 and d2 by zscore volume  */
    make_zscore_volume(d1,m1,&globals->threshold[0]);
    make_zscore_volume(d2,m2,&globals->threshold[1]);
  } 
  else  if (globals->obj_function == ssc_objective) {	/* add speckle to the data set */

    make_zscore_volume(d1,m1,&globals->threshold[0]); /* need to make data sets comparable */
    make_zscore_volume(d2,m2,&globals->threshold[1]); /* in mean and sd...                 */

    if (globals->smallest_vol == 1)
      add_speckle_to_volume(d1, 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);
    else
      add_speckle_to_volume(d2, 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);    
  } else if (globals->obj_function == vr_objective) {

    if (globals->smallest_vol == 1) {
      if (!build_segment_table(&segment_table, d1, globals->groups))
	print_error_and_line_num("Could not build segment table for source volume\n",__FILE__, __LINE__);
    }
    else {
      if (!build_segment_table(&segment_table, d2, globals->groups))
	print_error_and_line_num("Could not build segment table for target volume\n",__FILE__, __LINE__);
    }

  } else if (globals->obj_function == mutual_information_objective)
				/* Collignon's mutual information */
    {

      if ( globals->groups != 256 ) {
	print ("WARNING: -groups was %d, but will be forced to 256 in this run\n",globals->groups);
	globals->groups = 256;
      }

      data_type = get_volume_data_type (d1);
      if (data_type != UNSIGNED_BYTE) {
	print ("WARNING: source volume not UNSIGNED_BYTE, will do conversion now.\n");
	if (!replace_volume_data_with_ubyte(d1)) {
	  print_error_and_line_num("Can't replace volume data with unsigned bytes\n",
			     __FILE__, __LINE__);
	}
      }

      data_type = get_volume_data_type (d2);
      if (data_type != UNSIGNED_BYTE) {
	print ("WARNING: target volume not UNSIGNED_BYTE, will do conversion now.\n");
	if (!replace_volume_data_with_ubyte(d2)) {
	  print_error_and_line_num("Can't replace volume data with unsigned bytes\n",
			     __FILE__, __LINE__);
	}
      }

      ALLOC(   prob_fn1,   globals->groups);
      ALLOC(   prob_fn2,   globals->groups);
      ALLOC2D( prob_hash_table, globals->groups, globals->groups);

    } 
          /* ---------------- prepare the weighting array for obj func evaluation  ---------*/

          /* ---------------- prepare the weighting array for optimization ---------*/

  switch (globals->trans_info.transform_type) {
  case TRANS_LSQ3: 
    for_less(i,3,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
      globals->trans_info.scales[i] = 1.0;
      globals->trans_info.rotations[i] = 0.0;
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ6: 
    for_less(i,6,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
      globals->trans_info.scales[i] = 1.0;
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ7: 
    for_less(i,7,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ9: 
    for_less(i,9,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ10: 
    for_less(i,10,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,1,3) {
      globals->trans_info.shears[i] = 0.0;
    }
    break;
  case TRANS_LSQ: 
				/* nothing to be zeroed */
    break;
  case TRANS_LSQ12: 
				/* nothing to be zeroed */
    break;
  default:
    (void)fprintf(stderr, "Unknown type of transformation requested (%d)\n",
		   globals->trans_info.transform_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    stat = FALSE;
  }

	

  ndim = 0;
  for_less(i,0,12)
    if (globals->trans_info.weights[i] != 0.0) ndim++;

				/* set GLOBALS to communicate with the
				   function to be fitted!              */
  y = -1e10; 

  if (stat) {
    Gndim = ndim;
    if (globals->smallest_vol == 1) {
      Gdata1 = d1;
      Gdata2 = d2;
      Gmask1 = m1;
      Gmask2 = m2;
      Ginverse_mapping_flag = FALSE;
    }
    else {
      Gdata1 = d2;
      Gdata2 = d1;
      Gmask1 = m2;
      Gmask2 = m1;
      Ginverse_mapping_flag = TRUE;
    }

    ALLOC2D(p,ndim+1+1,ndim+1); /* simplex */
    
    parameters_to_vector(globals->trans_info.translations,
			 globals->trans_info.rotations,
			 globals->trans_info.scales,
			 globals->trans_info.shears,
			 p[1],
			 globals->trans_info.weights);

    y = fit_function(p[1]);	/* evaluate the objective  function */

    FREE2D(p); /* simplex */

  }
  
          /* ----------------finish up parameter/matrix manipulations ------*/

  if (globals->obj_function == vr_objective) {
    if (!free_segment_table(segment_table)) {
      (void)fprintf(stderr, "Can't free segment table.\n");
      (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    }
  } else
  if (globals->obj_function == mutual_information_objective)
				/* Collignon's mutual information */
    {
      FREE(   prob_fn1 );
      FREE(   prob_fn2 );
      FREE2D( prob_hash_table);
    }


  return(y);
}






/* ----------------------------- MNI Header -----------------------------------
@NAME       : optimize_non_linear_transformation
                get the parameters necessary to map volume 1 to volume 2
		using non-linear deformation fields.
@INPUT      : d1,d2:
                two volumes of data (already in memory).
	      m1,m2:
                two mask volumes for data (already in memory).
	      d1_dx, d1_dy, d1_dz, d1_dxyz,
	      d2_dx, d2_dy, d2_dz, d2_dxyz:
	        data sets corresponding to the intensity derivative in the x,y, 
		and z dirs.
	      globals:
	        a global data structure containing info from the command line,
		including the input parameters to be optimized, the input matrix,
		and a plethora of flags!
@OUTPUT     : 
@RETURNS    : TRUE if ok, FALSE if error.
@DESCRIPTION: 
@METHOD     :
                1- this routine begins by initializing the volume data structures
		to be used by the objective functions.

		2- optimization function is called
		
		3- the optimized parameters are returned in globals...
@GLOBALS    : 
@CALLS      : 
@CREATED    : Tue Nov 16 14:27:10 EST 1993 LC
@MODIFIED   : 
---------------------------------------------------------------------------- */
public BOOLEAN optimize_non_linear_transformation(Arg_Data *globals)
{
  BOOLEAN 
    stat;
  int i;

  stat = TRUE;
  
	     /*----------------- prepare data for optimization ------------ */
  
  if (globals->obj_function == zscore_objective) { /* replace volumes 
						      globals->features.data[0] and 
						      globals->features.model[0] 
						      by zscore volume  */

    make_zscore_volume(globals->features.data[0],
		       globals->features.data_mask[0],
		       &globals->threshold[0]);
    make_zscore_volume(globals->features.model[0],
		       globals->features.model_mask[0],
		       &globals->threshold[1]);

  } 
  else  if (globals->obj_function == ssc_objective) {	/* add speckle to the data set */

    make_zscore_volume(globals->features.data[0],        /* need to make data sets comparable */
		       globals->features.data_mask[0],   /* in mean and sd...                 */
		       &globals->threshold[0]); 
    make_zscore_volume(globals->features.model[0],
		       globals->features.model_mask[0],
		       &globals->threshold[1]); 

    if (globals->smallest_vol == 1)
      add_speckle_to_volume(globals->features.data[0], 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);
    else
      add_speckle_to_volume(globals->features.model[0], 
			    globals->speckle,
			    globals->start, globals->count, globals->directions);    
  } else if (globals->obj_function == vr_objective) {
    if (globals->smallest_vol == 1) {
      if (!build_segment_table(&segment_table, globals->features.data[0], globals->groups))
	print_error_and_line_num("Could not build segment table for source volume\n",__FILE__, __LINE__);
    }
    else {
      if (!build_segment_table(&segment_table, globals->features.model[0], globals->groups))
	print_error_and_line_num("Could not build segment table for target volume\n",__FILE__, __LINE__);

    }

    if (globals->flags.debug && globals->flags.verbose>1) {
      print ("groups = %d\n",segment_table->groups);
      for_less(i, segment_table->min, segment_table->max+1) {
	print ("%5d: table = %5d, function = %5d\n",i,segment_table->table[i],
	       (segment_table->segment)(i,segment_table) );
      }
    }

  }

  if (!globals->trans_info.use_magnitude) { /* then we are doing optical flow */
    normalize_data_to_match_target(globals->features.data[0],
				   globals->features.data_mask[0],
				   &globals->threshold[0],
				   globals->features.model[0],
				   globals->features.model_mask[0],
				   &globals->threshold[1],
				   globals);
				   
  }

	   /* ---------------- call requested optimization strategy ---------*/

  stat = ( do_non_linear_optimization(globals)==OK );


  
          /* ----------------finish up parameter/matrix manipulations ------*/

  if (stat && globals->obj_function == vr_objective) {
    stat = free_segment_table(segment_table);
  }


  return(stat);
}

