/* ----------------------------- MNI Header -----------------------------------
@NAME       : minctracc.h
@DESCRIPTION: Header file for minctracc.c
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

@CREATED    : Thu May 20 14:20:21 EST 1993 Louis Collins
@MODIFIED   : $Log: minctracc.h,v $
@MODIFIED   : Revision 9.5  1996-08-12 14:15:34  louis
@MODIFIED   : Release of MNI_AutoReg version 1.0
@MODIFIED   :
 * Revision 1.15  1996/08/12  14:15:26  louis
 * Pre-release
 *
 * Revision 1.14  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.13  94/06/19  15:42:15  louis
 * clean working version of 3D local deformation with simplex optimization
 * (by default) on magnitude data (default).  No more FFT stuff.
 * This working version before change of deformation field in do_nonlinear.c
 * 
 * 
 * Revision 1.12  94/06/06  09:46:52  louis
 * modified the initialization of main_args to reflect the use_magnitude
 * field in the trans_info struct.
 * 
 * Revision 1.11  94/05/28  16:18:55  louis
 * working version before modification of non-linear optimiation
 * 
 * Revision 1.10  94/04/26  12:55:22  louis
 * updated with new versions of make_rots, extract2_parameters_from_matrix 
 * that include proper interpretation of skew.
 * 
 * Revision 1.9  94/04/06  11:49:49  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.8  94/02/21  16:38:49  louis
 * version before feb 22 changes
 * 
 * Revision 1.7  93/11/15  13:12:53  louis
 * working version, deform deform installation
 * 
---------------------------------------------------------------------------- */

#include <minc.h>
#include <Proglib.h>


/* ------------------------  Constants used in program  -------------------- */

#include "constants.h"

/* ------------------------  Types used in program  ------------------------ */

#include "arg_data.h"

/*  ------------------------ Function prototypes  ------------------------ */

public int trilinear_interpolant(Volume volume, 
                                 PointR *coord, double *result);

public int tricubic_interpolant(Volume volume, 
                                PointR *coord, double *result);

public void do_Ncubic_interpolation(Volume volume, 
                                    long index[], int cur_dim, 
                                    double frac[], double *result);

public int nearest_neighbour_interpolant(Volume volume, 
                                         PointR *coord, double *result);
public int point_not_masked(Volume volume,
                            Real wx, Real wy, Real wz);

public int get_transformation(char *dst, char *key, char *nextArg);

public int get_mask_file(char *dst, char *key, char *nextArg);

public int get_feature_volumes(char *dst, char *key, int argc, char **argv);

public void procrustes(int npoints, int ndim, 
                       float **Apoints, float **Bpoints,
                       float *translation, float *centre_of_rotation,
                       float **rotation, float *scale);

public void transformations_to_homogeneous(int ndim, 
                  float *translation, float *centre_of_rotation,
                  float **rotation, float scale,
                  float **transformation);

public void translation_to_homogeneous(int ndim, float *translation,
                                       float **transformation);

public void rotation_to_homogeneous(int ndim, float **rotation,
                                       float **transformation);

public float fit_function(float *x);        /* apply cross correlation to the data sets    */

public float zscore_function(float *x);     /* calculate rms z-score difference.           */

public float check_function(float *x);      /* calculate the squared error between points2 */

public void invertmatrix(int n, float **mat, float **mat_invert);

public BOOLEAN init_params(Volume d1,
			   Volume d2,
			   Volume m1,
			   Volume m2, 
			   Arg_Data *globals);

public void init_lattice(Volume d1,
			 Volume d2,
			 Volume m1,
			 Volume m2, 
			 Arg_Data *globals);

public BOOLEAN optimize_linear_transformation(Volume d1,
					      Volume d2,
					      Volume m1,
					      Volume m2, 
					      Arg_Data *globals);

public BOOLEAN optimize_non_linear_transformation(Arg_Data *globals);

#include "objectives.h"

public float measure_fit(Volume d1,
			 Volume d2,
			 Volume m1,
			 Volume m2, 
			 Arg_Data *globals);

public void make_matlab_data_file(Volume d1,
				  Volume d2,
				  Volume m1,
				  Volume m2, 
				  char *comments,
				  Arg_Data *globals);

Status read_all_data(Volume *dblur,
		     Volume *dx,
		     Volume *dy,
		     Volume *dz,
		     Volume *dxyz, 
		     char *name);

public void build_default_deformation_field(Arg_Data *globals);


int allocate_a_new_feature(Feature_volumes *features);

void add_a_feature_for_matching(Feature_volumes *features,
				Volume data,
				Volume model,
				Volume data_mask,
				Volume model_mask,
				char *data_name,
				char *model_name,
				char *mask_data_name,
				char *mask_model_name,
				char obj_func,
				Real weight,
				Real thresh_data,
				Real thresh_model);

/*  ------------------------ Macros used in program  ------------------------ */

#include "local_macros.h"

/*  ------------------------ Global data structure for program  ------------------------ */

Arg_Data main_args = {
  {"","","","","","",""},	/* filenames           */
  {1,FALSE},			/* verbose, debug      */
  {				/* transformation info */
    FALSE,			/*   use identity tranformation to start */
    TRUE,			/*   do default tranformation (PAT) to start */
    TRUE,			/*   use_mag=TRUE; do not use projections by default */
    TRUE,			/*   use_simplex=TRUE ie use 3d simplex by default */
    2,				/*   use super sampling of deformation field  */
    FALSE,			/* use local smoothing       */
    TRUE,			/* use isotropic smoothing */
    "",			/*   filename */
    NULL,			/*   file_contents */
    0,                          /* buffer_length   */
    (General_transform *)NULL,	/*   General transform */
    (General_transform *)NULL,	/*   General transform copy of input */
    TRANS_PROCRUSTES,		/*   default type      */
    {-DBL_MAX, -DBL_MAX, -DBL_MAX},		/*   center            */
    {1.0, 1.0, 1.0},		/*   scale             */
    {0.0, 0.0, 0.0},		/*   shears            */
    {0.0, 0.0, 0.0},		/*   rotations         */
    {0.0, 0.0, 0.0},		/*   translations      */
    {1.0, 1.0, 1.0,  3.1415927/180.0, 3.1415927/180.0, 3.1415927/180.0,   0.02, 0.02, 0.02,  0.02, 0.02, 0.02}, /* optimization weights*/
    FALSE},			/*   invert_mapping_flag                  */
  {0,NULL, NULL, NULL, NULL, NULL, NULL},	/* FEATURE VOL */
  trilinear_interpolant,	/* use trilinear interpolation by default */
  xcorr_objective,              /* use cross-correlation by default       */
  OPT_SIMPLEX,                  /* use simplex optimization strategy      */
  {4.0,4.0,4.0},		/* default step sizes for lattice         */
  {0.0,0.0,0.0},		/* default start for lattice, reset in init_lattice */
  {0,0,0},                      /* default number of element in lattice, also reset */

  {{{1.0,0.0,0.0}},		/* default sampling lattice axes directions */
   {{0.0,1.0,0.0}},
   {{0.0,0.0,1.0}}},

  1,                            /* use first volume as default smallest volume      */
  {FALSE, FALSE, FALSE, FALSE},	/* Transform flags: est_cent, _scale, _rots, _trans */
  {0.0,0.0},			/* lower limit of voxels considered                 */
  5.0,				/* percent noise speckle                            */
  64				/* number of groups to use for ratio of variance    */
};


 
