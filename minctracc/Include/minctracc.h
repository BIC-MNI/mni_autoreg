#ifndef MINCTRACC_H
#define MINCTRACC_H

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
@MODIFIED   : Revision 96.11  2006-11-29 09:09:32  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.10  2006/06/04 07:02:35  rotor
@MODIFIED   :  * Fixed 64 bit function pointer and ParseArgv problem with an enum for
@MODIFIED   :       objective function type and interpolation type (thanks jason)
@MODIFIED   :
@MODIFIED   : Revision 96.9  2004/02/12 05:54:16  rotor
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 96.8  2003/02/04 06:08:44  stever
@MODIFIED   : Add support for correlation coefficient and sum-of-squared difference.
@MODIFIED   :
@MODIFIED   : Revision 96.7  2002/11/20 21:38:02  lenezet
@MODIFIED   :
@MODIFIED   : Fix the code to take in consideration the direction cosines especially in the grid transform.
@MODIFIED   : Add an option to choose the maximum expected deformation magnitude.
@MODIFIED   :
@MODIFIED   : Revision 96.6  2002/08/14 19:54:49  lenezet
@MODIFIED   :  quaternion option added for the rotation
@MODIFIED   :
@MODIFIED   : Revision 96.5  2002/03/07 19:07:51  louis
@MODIFIED   : Added -lattice_diameter as an optionto minctracc to account for a
@MODIFIED   : problem with the automated calculation of the sub-lattice diameter.
@MODIFIED   : It used to be step*3*2 - which was pretty big, when step = 8mm.
@MODIFIED   :
@MODIFIED   : Now, the sub lattice diameter can be input on the command line, and I
@MODIFIED   : suggest a lattice size 3 times greater than the step size.
@MODIFIED   :
@MODIFIED   : If not on the command line, the default is = 24mm.
@MODIFIED   :
@MODIFIED   : Revision 96.4  2000/03/15 08:42:39  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 96.3  2000/02/15 19:02:06  stever
@MODIFIED   : Add tests for param2xfm, minctracc -linear.
@MODIFIED   :
@MODIFIED   : Revision 96.2  1999/10/25 19:52:16  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.1  1997/11/03  19:52:55  louis
 * changed the number of default groups from 64 to 256 for -mi.
 * added the default pdf blurring size for -mi
 *
 * Revision 96.0  1996/08/21  18:21:43  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.5  1996/08/12  14:15:34  louis
 * Release of MNI_AutoReg version 1.0
 *
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

#include "minctracc_arg_data.h"

/*  ------------------------ Function prototypes  ------------------------ */

#include "interpolation.h"

int get_transformation(char *dst, char *key, char *nextArg);

int get_mask_file(char *dst, char *key, char *nextArg);

int get_nonlinear_objective(char *dst, char *key, char *nextArg);

int get_feature_volumes(char *dst, char *key, int argc, char **argv);

void procrustes(int npoints, int ndim, 
                       float **Apoints, float **Bpoints,
                       float *translation, float *centre_of_rotation,
                       float **rotation, float *scale);

void transformations_to_homogeneous(int ndim, 
                  float *translation, float *centre_of_rotation,
                  float **rotation, float scale,
                  float **transformation);

void translation_to_homogeneous(int ndim, float *translation,
                                       float **transformation);

void rotation_to_homogeneous(int ndim, float **rotation,
                                       float **transformation);

float fit_function(float *x);        /* apply cross correlation to the data sets    */

float zscore_function(float *x);     /* calculate rms z-score difference.           */

float check_function(float *x);      /* calculate the squared error between points2 */

void invertmatrix(int n, float **mat, float **mat_invert);

VIO_BOOL vol_to_cov(VIO_Volume d1, VIO_Volume m1, float *centroid, float **covar, double *step);


VIO_BOOL init_params(VIO_Volume d1,
                           VIO_Volume d2,
                           VIO_Volume m1,
                           VIO_Volume m2, 
                           Arg_Data *globals);

VIO_BOOL init_params_quater(VIO_Volume d1,
                                  VIO_Volume d2,
                                  VIO_Volume m1,
                                  VIO_Volume m2, 
                                  Arg_Data *globals);


void init_lattice(VIO_Volume d1,
                         VIO_Volume d2,
                         VIO_Volume m1,
                         VIO_Volume m2, 
                         Arg_Data *globals);

VIO_BOOL optimize_linear_transformation(VIO_Volume d1,
                                              VIO_Volume d2,
                                              VIO_Volume m1,
                                              VIO_Volume m2, 
                                              Arg_Data *globals);

VIO_BOOL optimize_linear_transformation_quater(VIO_Volume d1,
                                                     VIO_Volume d2,
                                                     VIO_Volume m1,
                                                     VIO_Volume m2, 
                                                     Arg_Data *globals);

VIO_BOOL optimize_non_linear_transformation(Arg_Data *globals);

#include "objectives.h"

float measure_fit(VIO_Volume d1,
                         VIO_Volume d2,
                         VIO_Volume m1,
                         VIO_Volume m2, 
                         Arg_Data *globals);

void make_matlab_data_file(VIO_Volume d1,
                                  VIO_Volume d2,
                                  VIO_Volume m1,
                                  VIO_Volume m2, 
                                  char *comments,
                                  Arg_Data *globals);

VIO_Status read_all_data(VIO_Volume *dblur,
                     VIO_Volume *dx,
                     VIO_Volume *dy,
                     VIO_Volume *dz,
                     VIO_Volume *dxyz, 
                     char *name);

void build_default_deformation_field(Arg_Data *globals);


int allocate_a_new_feature(Feature_volumes *features);

void add_a_feature_for_matching(Feature_volumes *features,
                                VIO_Volume data,
                                VIO_Volume model,
                                VIO_Volume data_mask,
                                VIO_Volume model_mask,
                                char *data_name,
                                char *model_name,
                                char *mask_data_name,
                                char *mask_model_name,
                                char obj_func,
                                VIO_Real weight,
                                VIO_Real thresh_data,
                                VIO_Real thresh_model);


int minctraccOldFashioned ( int argc, char* argv[] );


/*---------------------- functions relatives to quaternions-------------------------------------------*/
#include "quaternion.h" 



/*  ------------------------ Macros used in program  ------------------------ */

#include "local_macros.h"

/*  ------------------------ Global data structure for program  ------------------------ */

extern Arg_Data *main_args;

 
#endif