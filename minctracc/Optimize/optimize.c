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
@MODIFIED   : Revision 1.9  1995-02-22 08:56:06  louis
@MODIFIED   : Montreal Neurological Institute version.
@MODIFIED   : compiled and working on SGI.  this is before any changes for SPARC/
@MODIFIED   : Solaris.
@MODIFIED   :
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
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Optimize/optimize.c,v 1.9 1995-02-22 08:56:06 louis Exp $";
#endif

#include <volume_io.h>
#include <recipes.h>
#include <limits.h>
#include <print_error.h>

#include "constants.h"
#include "arg_data.h"
#include "objectives.h"
#include "make_rots.h"
#include "segment_table.h"

/* external calls: */

public void make_zscore_volume(Volume d1, Volume m1, 
			       float threshold); 

public void add_speckle_to_volume(Volume d1, 
				  float speckle,
				  double  *start, int *count, VectorR directions[]);

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
					 Arg_Data *globals);

extern Arg_Data main_args;

Volume   Gdata1, Gdata2, Gmask1, Gmask2;
int      Ginverse_mapping_flag, Gndim;

extern   double   ftol ;        
extern   double   simplex_size ;

         Segment_Table  *segment_table;


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

private BOOLEAN in_limits(double x,double lower,double upper)
{

  if ( x>=lower && x <= upper ) 
    return(TRUE);
  else
    return(FALSE);
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

printf("out : %7.4f=%c %7.4f=%c %7.4f=%c   %7.4f=%c %7.4f=%c %7.4f=%c   %7.4f=%c %7.4f=%c %7.4f=%c \n",
       rots[0], in_limits(rots[0], (double)-3.1415927/2.0, (double)3.1415927/2.0) ? 'T': 'F' , 
       rots[1], in_limits(rots[1], (double)-3.1415927/2.0, (double)3.1415927/2.0) ? 'T': 'F' , 
       rots[2], in_limits(rots[2], (double)-3.1415927/2.0, (double)3.1415927/2.0) ? 'T': 'F' , 
       scale[0],in_limits(scale[0], (double)0.0, (double)3.0)? 'T': 'F' , 
       scale[1],in_limits(scale[1], (double)0.0, (double)3.0)? 'T': 'F' , 
       scale[2],in_limits(scale[2], (double)0.0, (double)3.0)? 'T': 'F' , 
       shear[0],in_limits(shear[0], (double)-2.0, (double)2.0)? 'T': 'F' , 
       shear[1],in_limits(shear[1], (double)-2.0, (double)2.0)? 'T': 'F' , 
       shear[2],in_limits(shear[2], (double)-2.0, (double)2.0)? 'T': 'F' );

    r = FLT_MAX;
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
@METHOD     : uses the ameoba simplex algorithm from numerical recipes. 
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
    **p, *y;
  int 
    i,j, 
    ndim, nfunk;

  Transform
    *mat;

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

  }


  if (stat && ndim>0) {

    p = matrix(1,ndim+1,1,ndim); /* simplex */
    y = vector(1,ndim+1);        /* value of correlation at simplex vertices */
    
    parameters_to_vector(globals->trans_info.translations,
			 globals->trans_info.rotations,
			 globals->trans_info.scales,
			 globals->trans_info.shears,
			 p[1],
			 globals->trans_info.weights);

    for_inclusive(i,2,ndim+1)	/* copy initial guess to all points of simplex */
      for_inclusive(j,1,ndim)
	p[i][j] = p[1][j];
    
				/* set up all vertices of simplex */
    for_inclusive(j,1,ndim) {
	p[j+1][j] = p[1][j] + simplex_size;
    }

    for (i=1; i<=(ndim+1); ++i)	{   /* set up value of correlation at all points of simplex */
      
      y[i] = fit_function(p[i]);

      (void)print ("corr = %6.4f",y[i]);
      for (j=1; j<=ndim; ++j)  {
	(void)print (" %9.4f",p[i][j]);
      }
      (void)print ("\n");
      
    }
	
    amoeba(p,y,ndim,local_ftol,fit_function,&nfunk);
    
    (void)print("after %d iterations\n",nfunk);
    
    for (i=1; i<=(ndim+1); ++i)	{   /* print out value of correlation at all points of simplex */
      if (i==1) {
	(void)print ("end corr = %f",y[i]);
	for (j=1; j<=ndim; ++j)  {
	  (void)print (" %9.4f",p[i][j]);
	}
	(void)print ("\n");
       } 
    }   

				/* copy result into main data structure */

    vector_to_parameters(globals->trans_info.translations,
			 globals->trans_info.rotations,
			 globals->trans_info.scales,
			 globals->trans_info.shears,
			 p[1],
			 globals->trans_info.weights);

print ("%d == %d\n",globals->trans_info.transform_type,TRANS_LSQ7);

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

    free_matrix(p,1,ndim+1,1,ndim); /* simplex */
    free_vector(y,1,ndim+1);        /* value of correlation at simplex vertices */

  }

  return( stat );
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

  stat = TRUE;
  
	     /*----------------- prepare data for optimization ------------ */
  
  if (globals->obj_function == zscore_objective) { /* replace volume d1 and d2 by zscore volume  */
    make_zscore_volume(d1,m1,(float)globals->threshold[0]);
    make_zscore_volume(d2,m2,(float)globals->threshold[1]);
  } 
  else  if (globals->obj_function == ssc_objective) {	/* add speckle to the data set */

    make_zscore_volume(d1,m1,(float)globals->threshold[0]); /* need to make data sets comparable */
    make_zscore_volume(d2,m2,(float)globals->threshold[1]); /* in mean and sd...                 */

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
	print_error("Could not build segment table for source volume\n",__FILE__, __LINE__);
    }
    else {
      if (!build_segment_table(&segment_table, d2, globals->groups))
	print_error("Could not build segment table for target volume\n",__FILE__, __LINE__);

    }

    if (globals->flags.debug && globals->flags.verbose>1) {
      print ("groups = %d\n",segment_table->groups);
      for_less(i, segment_table->min, segment_table->max+1) {
	print ("%5d: table = %5d, function = %5d\n",i,segment_table->table[i],
	       (segment_table->segment)(i,segment_table) );
      }
    }

  }

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
  case TRANS_PROCRUSTES: 
    for_less(i,7,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
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

	   /* ---------------- call requested optimization strategy ---------*/

  switch (globals->optimize_type) {
  case OPT_SIMPLEX:
    stat = optimize_simplex(d1, d2, m1, m2, globals);
    break;
  default:
    (void)fprintf(stderr, "Unknown type of optimization requested (%d)\n",globals->optimize_type);
    (void)fprintf(stderr, "Error in line %d, file %s\n",__LINE__, __FILE__);
    stat = FALSE;
  }
  
          /* ----------------finish up parameter/matrix manipulations ------*/

  if (globals->obj_function == vr_objective) {
    stat = free_segment_table(segment_table);
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

/*  Transform
    *mat;
    double trans[3];  double cent[3];  double rots[3];  double scale[3];  double shear[6]; */

  
  stat = TRUE;

	     /*----------------- prepare data for objective function evaluation ------------ */

  
  if (globals->obj_function == zscore_objective) { /* replace volume d1 and d2 by zscore volume  */
    make_zscore_volume(d1,m1,(float)globals->threshold[0]);
    make_zscore_volume(d2,m2,(float)globals->threshold[1]);
  } 
  else  if (globals->obj_function == ssc_objective) {	/* add speckle to the data set */

    make_zscore_volume(d1,m1,(float)globals->threshold[0]); /* need to make data sets comparable */
    make_zscore_volume(d2,m2,(float)globals->threshold[1]); /* in mean and sd...                 */

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
	print_error("Could not build segment table for source volume\n",__FILE__, __LINE__);
    }
    else {
      if (!build_segment_table(&segment_table, d2, globals->groups))
	print_error("Could not build segment table for target volume\n",__FILE__, __LINE__);
    }

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
  case TRANS_PROCRUSTES: 
    for_less(i,7,12) globals->trans_info.weights[i] = 0.0;
    for_less(i,0,3) {
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
  }



  if (stat) {

    p = matrix(1,ndim+1,1,ndim); /* simplex */
    
    parameters_to_vector(globals->trans_info.translations,
			 globals->trans_info.rotations,
			 globals->trans_info.scales,
			 globals->trans_info.shears,
			 p[1],
			 globals->trans_info.weights);

    y = fit_function(p[1]);	/* evaluate the objective  function */

    free_matrix(p,1,ndim+1,1,ndim); /* simplex */

  }
  
          /* ----------------finish up parameter/matrix manipulations ------*/

  if (globals->obj_function == vr_objective) {
    stat = free_segment_table(segment_table);
  }

  if (stat)
    return(y);
  else
    return(-FLT_MAX);
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
public BOOLEAN optimize_non_linear_transformation(Volume d1,
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
  BOOLEAN 
    stat;
  int i;

  stat = TRUE;
  
	     /*----------------- prepare data for optimization ------------ */
  
  if (globals->obj_function == zscore_objective) { /* replace volume d1 and d2 by zscore volume  */
    make_zscore_volume(d1,m1,(float)globals->threshold[0]);
    make_zscore_volume(d2,m2,(float)globals->threshold[1]);
  } 
  else  if (globals->obj_function == ssc_objective) {	/* add speckle to the data set */

    make_zscore_volume(d1,m1,(float)globals->threshold[0]); /* need to make data sets comparable */
    make_zscore_volume(d2,m2,(float)globals->threshold[1]); /* in mean and sd...                 */

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
	print_error("Could not build segment table for source volume\n",__FILE__, __LINE__);
    }
    else {
      if (!build_segment_table(&segment_table, d2, globals->groups))
	print_error("Could not build segment table for target volume\n",__FILE__, __LINE__);

    }

    if (globals->flags.debug && globals->flags.verbose>1) {
      print ("groups = %d\n",segment_table->groups);
      for_less(i, segment_table->min, segment_table->max+1) {
	print ("%5d: table = %5d, function = %5d\n",i,segment_table->table[i],
	       (segment_table->segment)(i,segment_table) );
      }
    }

  }
	   /* ---------------- call requested optimization strategy ---------*/

  stat = (do_non_linear_optimization(d1,d1_dx, d1_dy, d1_dz, d1_dxyz,
				    d2,d2_dx, d2_dy, d2_dz, d2_dxyz,
				    m1,m2, 
				    globals)==OK);


  
          /* ----------------finish up parameter/matrix manipulations ------*/

  if (stat && globals->obj_function == vr_objective) {
    stat = free_segment_table(segment_table);
  }


  return(stat);
}

