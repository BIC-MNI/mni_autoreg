/* ----------------------------- MNI Header -----------------------------------
@NAME       : transformations.c
@DESCRIPTION: routines to apply the forward and inverse transformations
              of the non-linear deformation field.
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

@CREATED    : Tue Nov 16 14:51:04 EST 1993 lc
                    based on transformations.c from fit_vol
@MODIFIED   : $Log: transformations.c,v $
@MODIFIED   : Revision 1.13  2006-11-29 09:09:33  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 1.12  2004/02/12 16:22:43  louis
@MODIFIED   :  * removed public/private defs
@MODIFIED   :
@MODIFIED   : Revision 1.11  1995/09/11 12:37:16  louis
@MODIFIED   : All refs to numerical recipes routines have been replaced.
@MODIFIED   : this is an updated working version - corresponds to mni_reg-0.1g
@MODIFIED   :
 * Revision 1.10  1995/07/04  11:42:45  louis
 * removed apply_deformation_field and apply_inverse_deform... since they
 * are now included in volume_io
 *
 * Revision 1.9  1995/02/22  08:56:06  louis
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.8  94/06/06  09:37:56  louis
 * modifed the voxel and real range calls, in build_default_deformation_field
 * where the new voxel range: 0.0, 32767.0, new real range: -50.0, 50.0.
 * 
 * These ranges are copied for all other deformation fields created from
 * the 1st scale deformation field.
 * 
 * Revision 1.7  94/06/02  20:15:59  louis
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
 * Revision 1.6  94/05/28  16:18:29  louis
 * working version before modification of non-linear optimiation
 * 
 * Revision 1.5  94/04/06  11:48:52  louis
 * working linted version of linear + non-linear registration based on Lvv
 * operator working in 3D
 * 
 * Revision 1.4  94/02/21  16:37:02  louis
 * version before feb 22 changes
 * 
---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /private-cvsroot/registration/mni_autoreg/minctracc/Numerical/transformations.c,v 1.13 2006-11-29 09:09:33 rotor Exp $";
#endif


#include <volume_io.h>

#include "minctracc_arg_data.h"
#include "local_macros.h"

#include "constants.h"

#define NUMBER_TRIES 10

void build_default_deformation_field(Arg_Data *globals);


void set_up_lattice(VIO_Volume data, 
                            double *user_step, /* user requested spacing for lattice */
                            double *start,     /* world starting position of lattice */
                            int    *count,     /* number of steps in each direction */
                            double *step,      /* step size in each direction */
                            VectorR directions[]);/* array of vector directions for each index*/

 int tricubic_interpolant(VIO_Volume volume, 
                                PointR *coord, double *result);

extern Arg_Data *main_args;


 VIO_General_transform *get_linear_part_of_transformation(VIO_General_transform *trans)
{
  VIO_General_transform *result,*concated,*current_lin;
  int i;

  ALLOC(result, 1);
  ALLOC(concated,1 );
  ALLOC(current_lin,1);

  create_linear_transform(result, NULL); /* start with identity */

  for(i=0; i<get_n_concated_transforms(trans; i++)) {
    if (get_transform_type( get_nth_general_transform(trans, i-0) ) == LINEAR){

      copy_general_transform( get_nth_general_transform(trans, i-0), current_lin);
      concat_general_transforms(result, current_lin, concated);

      delete_general_transform(result);
      delete_general_transform(current_lin);
      copy_general_transform(concated, result);
      delete_general_transform(concated);

   }
  }

  return(result);
}

