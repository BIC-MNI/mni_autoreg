/* ----------------------------- MNI Header -----------------------------------
@NAME       : deform_support.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
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

@CREATED    : Tue Feb 22 08:37:49 EST 1994
@MODIFIED   : $Log: deform_support.c,v $
@MODIFIED   : Revision 96.15  2009-04-03 18:36:59  louis
@MODIFIED   : made changes to use only DOUBLES for input source and model volumes, and for all estimation of deformation fields
@MODIFIED   :
@MODIFIED   : Revision 96.14  2006/11/30 17:23:43  rotor
@MODIFIED   :  * fixed a small bug in init_volume_to_zero
@MODIFIED   :
@MODIFIED   : Revision 96.13  2006/11/30 09:07:32  rotor
@MODIFIED   :  * many more changes for clean minc 2.0 build
@MODIFIED   :
@MODIFIED   : Revision 96.12  2006/11/29 09:09:34  rotor
@MODIFIED   :  * first bunch of changes for minc 2.0 compliance
@MODIFIED   :
@MODIFIED   : Revision 96.11  2005/07/20 20:45:50  rotor
@MODIFIED   :     * Complete rewrite of the autoconf stuff (configure.in -> configure.am)
@MODIFIED   :     * Many changes to includes of files (float.h, limits.h, etc)
@MODIFIED   :     * Removed old VOLUME_IO cruft #defines
@MODIFIED   :     * Fixed up all Makefile.am's in subdirs
@MODIFIED   :     * Removed all things in Proglib that are now part of MINC proper
@MODIFIED   :     * Still working on fixing up perl subdirectory - removing mni_perllib
@MODIFIED   :
@MODIFIED   : Revision 96.10  2004/02/12 06:08:19  rotor
@MODIFIED   :  * removed /static defs
@MODIFIED   :
@MODIFIED   : Revision 96.9  2004/02/04 20:44:11  lenezet
@MODIFIED   : *** empty log message ***
@MODIFIED   :
@MODIFIED   : Revision 96.8  2002/12/13 21:18:20  lenezet
@MODIFIED   :
@MODIFIED   : A memory leak has been repaired
@MODIFIED   :
@MODIFIED   : Revision 96.7  2002/11/20 21:39:14  lenezet
@MODIFIED   :
@MODIFIED   : Fix the code to take in consideration the direction cosines especially in the grid transform.
@MODIFIED   : Add an option to choose the maximum expected deformation magnitude.
@MODIFIED   :
@MODIFIED   : Revision 96.6  2002/03/26 14:15:43  stever
@MODIFIED   : Update includes to <volume_io/foo.h> style.
@MODIFIED   :
@MODIFIED   : Revision 96.5  2000/03/15 08:42:45  stever
@MODIFIED   : Code cleanup: all functions prototyped (except ParseArgs.c), no useless declarations, etc
@MODIFIED   :
@MODIFIED   : Revision 96.4  2000/02/07 19:33:05  stever
@MODIFIED   : replaced HAVE_RECENT_VOLUME_IO with more specific feature tests.
@MODIFIED   :
@MODIFIED   : Revision 96.3  1999/10/25 19:59:06  louis
@MODIFIED   : final checkin before switch to CVS
@MODIFIED   :
 * Revision 96.2  1997/11/03  20:05:41  louis
 * reorganized deform_support...
 * putting:
 *   init_the_volume_to_zero()
 *   get_volume_maximum_real_value()
 *   save_data()
 *  into extras.c
 *
 * and putting:
 *   general_transform_point_in_trans_plane()
 *   build_source_lattice()
 *   go_get_samples_in_source()
 *   go_get_samples_with_offset()
 *   build_target_lattice()
 *   build_target_lattice_using_super_sampled_def()
 *  into sub_lattice.c
 *
 * Revision 96.2  1997/11/03  20:05:41  louis
 * reorganized deform_support...
 * putting:
 *   init_the_volume_to_zero()
 *   get_volume_maximum_real_value()
 *   save_data()
 *  into extras.c
 *
 * and putting:
 *   general_transform_point_in_trans_plane()
 *   build_source_lattice()
 *   go_get_samples_in_source()
 *   go_get_samples_with_offset()
 *   build_target_lattice()
 *   build_target_lattice_using_super_sampled_def()
 *  into sub_lattice.c
 *
 * Revision 96.1  1997/11/03  15:06:29  louis
 * working version, before creation of mni_animal package, and before inserting
 * distance transforms
 *
 * Revision 96.0  1996/08/21  18:22:10  louis
 * Release of MNI_AutoReg version 0.96
 *
 * Revision 9.6  1996/08/21  18:22:01  louis
 * Pre-release
 *
 * Revision 9.5  1996/08/12  14:16:03  louis
 * Never released version 0.95
 *
 * Revision 1.16  1996/08/12  14:15:53  louis
 * Pre-release
 *
 * Revision 1.14  1996/04/01  09:16:29  collins
 * removed the code to super sample the deformation field,
 * and placed it in super_sample_def.c.
 *
 * Revision 1.13  1996/04/01  09:02:10  collins
 * added optimized code to super-sample the deformation field for the special
 * case when super_sample = 2.  interpolate_super_sampled_data_by2() is
 * approximately 7 times faster than interpolate_super_sampled_data().
 *
 * Revision 1.12  1996/03/25  10:33:15  collins
 * used inter_type to specify degress_continuity in the call to
 * evaluate_volume_in_world() in procedure go_get_samples_in_source().
 *
 * ----
 * added support for tri-linear interpolation in
 * go_get_samples_with_offset() since the local objective function is not
 * smooth (but is actually step-wise) with NN interpolation.  This is
 * specified by an extra boolean parameter to go_get_samples_with_offset
 * that should be TRUE to use NN interpolation.  This is controled from
 * the command line by the -nearest_neighbour / -trilinear options.
 *
 * NOTE: NN is faster, but the local obj function will have steps (ie it
 *       will have plateaus and discontinuous changes).  Also, there
 *       seems to be a `drapery' effect, such that the function is _not_
 *       monotonically increasing with mis-registration -> thus yielding
 *       local minima.
 *       -> trilin is the default, and should be used unless you really
 *          know what you are doing!
 *
 * ----
 *
 * removed dx += 0.5; dy += 0.5; dz += 0.5.; in
 * go_get_samples_with_offset().  this was originally supposed to give a
 * `rounding' effect for NN interpolation, but is no longer needed with
 * tri-linear interpolation.  IN FACT: keeping it caused a certain bias
 * so that a deformation was returned _even_ when an object was matched
 * to itself; after removal- there are only small sporadic deviations
 * from a null deformation field.
 *
 * Revision 1.11  1996/03/07  13:25:19  collins
 * small reorganisation of procedures and working version of non-isotropic
 * smoothing.
 *
 * Revision 1.10  1995/10/06  09:25:02  collins
 * removed references to line_data.h since it hos not been used in a while.
 *
 * included "constants.h" to have access to NONLIN_* similarity func ids.
 *
 * modified go_get_samples_with_offset to account for different similarity
 * functions.
 *
 * Revision 1.9  1995/09/07  10:05:11  collins
 * All references to numerical recipes routines are being removed.  At this
 * stage, any num rec routine should be local in the file.  All memory
 * allocation calls to vector(), matrix(), free_vector() etc... have been
 * replaced with ALLOC and FREE from the volume_io library of routines.
 *
 * Revision 1.8  1995/06/12  14:29:46  collins
 * working version - 2d,3d w/ simplex and -direct.
 *
 * Revision 1.7  1995/05/04  14:25:18  collins
 * compilable version, seems to run a bit with GRID_TRANSFORM, still
 * needs work for the super sampled volumes... and lots of testing.
 *
 * Revision 1.6  1995/05/02  11:30:27  collins
 * started clean up of code, separation of non used procedures into
 * old_methods.c.  This version was working, but I am now going to
 * rewrite everything to use GRID_TRANSFORM.
 *
 * Revision 1.5  1995/02/22  08:56:06  collins
 * Montreal Neurological Institute version.
 * compiled and working on SGI.  this is before any changes for SPARC/
 * Solaris.
 *
 * Revision 1.4  94/06/21  10:59:28  louis
 * working optimized version of program.  when compiled with -O, this
 * code is approximately 60% faster than the previous version.
 * 
 * 
 * Revision 1.3  94/06/06  18:46:53  louis
 * working version: clamp and blur of deformation lattice now ensures
 * a smooth recovered deformation.  Unfortunately, the r = cost-similarity
 * function used in the optimization is too heavy on the cost_fn.  This has
 * to get fixed...
 * 
 * 

 * Revision 1.2  94/06/02  20:15:56  louis
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
 * Revision 1.1  94/04/06  11:47:27  louis
 * Initial revision
 * 

---------------------------------------------------------------------------- */

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/registration/mni_autoreg/minctracc/Optimize/deform_support.c,v 96.15 2009-04-03 18:36:59 louis Exp $";
#endif

#include <config.h>
#include <float.h>
#include <volume_io.h>
#include "minctracc_arg_data.h"
#include "louis_splines.h"
#include <Proglib.h>
#include "local_macros.h"
#include "constants.h"
#include "interpolation.h"

extern Arg_Data *main_args;

#define DERIV_FRAC      0.6
#define FRAC1           0.5
#define FRAC2           0.0833333
#define ABSOLUTE_MAX_DEFORMATION       50.0

extern double smoothing_weight;
extern char *my_XYZ_dim_names;

void get_volume_XYZV_indices(VIO_Volume data, int xyzv[]);

int trilinear_interpolant(VIO_Volume volume, 
                                 PointR *coord, double *result);
 
int nearest_neighbour_interpolant(VIO_Volume volume, 
                                         PointR *coord, double *result);

VIO_Real get_volume_maximum_real_value(VIO_Volume volume);

VIO_Real get_coeff_from_neighbours(VIO_General_transform *trans,
                                      int voxel[],
                                      int avg_type); 

                                /* define the limits for three nested
                                   for loops, so that we loop through
                                   each spatial dimension.
                                   'start' and 'end' are returned in
                                   XYZ order.  */

void  get_voxel_spatial_loop_limits(VIO_Volume volume,
                                           int start[],        
                                           int end[])
{
  int 
    i,
    sizes[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS];

  get_volume_sizes(volume, sizes);
  get_volume_XYZV_indices(volume, xyzv);

  for(i=0; i<VIO_N_DIMENSIONS; i++) {
    if (sizes[ xyzv[i] ]>3) {
      start[ i ] = 1;
      end[ i ] = sizes[ xyzv[i] ]-1;
    }
    else {
      start[ i ] = 0;
      end[ i ] = sizes[ xyzv[i] ];
    }
  }

  if (get_volume_n_dimensions(volume)>3) {
    start[VIO_Z+1] = 0;
    end[VIO_Z+1] = sizes[ xyzv[VIO_Z+1] ] ;
  }
  else {
    start[VIO_Z+1] = 0;
    end[VIO_Z+1] = 0;
  }
    
}

VIO_BOOL get_average_warp_vector_from_neighbours(VIO_General_transform *trans,
                                                       int voxel[],
                                                       int avg_type,
                                                       VIO_Real *mx, VIO_Real *my, VIO_Real *mz)
{
  int
    start[VIO_MAX_DIMENSIONS],
    end[VIO_MAX_DIMENSIONS],
    count, i, 
    voxel2[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    sizes[VIO_MAX_DIMENSIONS];
  VIO_Real 
    def_vector[VIO_N_DIMENSIONS];
  VIO_Volume volume;
  
  if (trans->type != GRID_TRANSFORM) {
    print_error_and_line_num("get_average_warp_vector_from_neighbours not called with GRID_TRANSFORM",
                             __FILE__, __LINE__);
    return (FALSE);
  }

  volume = trans->displacement_volume;
  
  get_volume_sizes(volume, sizes);
  get_volume_XYZV_indices(volume, xyzv);
  count = 0;
  *mx = 0.0; *my = 0.0; *mz = 0.0; /* assume no warp vector in volume */

                                /* make sure voxel is in volume */
  for(i=0; i<3; i++) {
    if (voxel[ xyzv[i]]<0 || voxel[ xyzv[i] ]>=sizes[ xyzv[i]] ) {
       return (FALSE);
    }
  }
  
  for(i=0; i<VIO_MAX_DIMENSIONS; i++) { /* copy the voxel position */
    voxel2[i] = voxel[i];
  }
  
  switch (avg_type) {
  case 1:  {                        /* get the 6 neighbours, 2 along
                                   each of the spatial axes,
                                   if they exist 
                                   ie the 4-connected immediate neighbours only */

    for(i=0; i<VIO_N_DIMENSIONS; i++) {
      
      if ((voxel[ xyzv[i] ]+1) < sizes[ xyzv[i] ]) {
        
        voxel2[ xyzv[i] ] = voxel[ xyzv[i] ] + 1;
        
        for(voxel2[xyzv[VIO_Z+1]]=0; voxel2[xyzv[VIO_Z+1]]<sizes[xyzv[VIO_Z+1]]; voxel2[xyzv[VIO_Z+1]]++) {
          def_vector[voxel2[ xyzv[VIO_Z+1] ]] = 
            get_volume_real_value(volume,
                                  voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
        }
        
        voxel2[ xyzv[i] ] = voxel[ xyzv[i] ];
        
        *mx += def_vector[VIO_X]; *my += def_vector[VIO_Y]; *mz += def_vector[VIO_Z];
        ++count;
      }
      if ((voxel[ xyzv[i] ]-1) >= 0) {
        voxel2[ xyzv[i] ] = voxel[ xyzv[i] ] - 1;
        
        for(voxel2[xyzv[VIO_Z+1]]=0; voxel2[xyzv[VIO_Z+1]]<sizes[xyzv[VIO_Z+1]]; voxel2[xyzv[VIO_Z+1]]++) {
          def_vector[voxel2[ xyzv[VIO_Z+1] ]] = 
            get_volume_real_value(volume,
                                  voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
        }
        
        voxel2[ xyzv[i] ] = voxel[ xyzv[i] ];
        
        *mx += def_vector[VIO_X]; *my += def_vector[VIO_Y]; *mz += def_vector[VIO_Z];
        ++count;
      }
      
    }
    break;
  }
  case 2: {                        /* 3x3x3 ie the 26 8-connected
                                   immediate neighbours only */
    
    for(i=0; i<VIO_N_DIMENSIONS; i++) {
      start[i] = voxel[ xyzv[i] ] - 1;
      if (start[i]<0) start[i]=0;
      end[i] = voxel[ xyzv[i] ] + 1;
      if (end[i]>sizes[ xyzv[i] ]-1) end[i] = sizes[ xyzv[i] ]-1;
    }
    
    
    for(voxel2[xyzv[VIO_X]]=start[VIO_X]; voxel2[xyzv[VIO_X]]<=end[VIO_X]; voxel2[xyzv[VIO_X]]++)
      for(voxel2[xyzv[VIO_Y]]=start[VIO_Y]; voxel2[xyzv[VIO_Y]]<=end[VIO_Y]; voxel2[xyzv[VIO_Y]]++)
        for(voxel2[xyzv[VIO_Z]]=start[VIO_Z]; voxel2[xyzv[VIO_Z]]<=end[VIO_Z]; voxel2[xyzv[VIO_Z]]++) {

          if ((voxel2[ xyzv[VIO_X]] != voxel[ xyzv[VIO_X] ]) ||
              (voxel2[ xyzv[VIO_Y]] != voxel[ xyzv[VIO_Y] ]) ||
              (voxel2[ xyzv[VIO_Z]] != voxel[ xyzv[VIO_Z] ])) {
            for(voxel2[xyzv[VIO_Z+1]]=0; voxel2[xyzv[VIO_Z+1]]<sizes[xyzv[VIO_Z+1]]; voxel2[xyzv[VIO_Z+1]]++) {
              def_vector[voxel2[ xyzv[VIO_Z+1] ]] = 
                get_volume_real_value(volume,
                                      voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
            }
            *mx += def_vector[VIO_X]; *my += def_vector[VIO_Y]; *mz += def_vector[VIO_Z];
            ++count;
          }
        }
    break;
  }
  case 3: {                        /* 5x5x5 */
    for(i=0; i<VIO_N_DIMENSIONS; i++) {
      start[i] = voxel[ xyzv[i] ] - 2;
      if (start[i]<0) start[i]=0;
      end[i] = voxel[ xyzv[i] ] + 2;
      if (end[i]>sizes[ xyzv[i] ]-1) end[i] = sizes[ xyzv[i] ]-1;
    }
    
    
    for(voxel2[xyzv[VIO_X]]=start[VIO_X]; voxel2[xyzv[VIO_X]]<=end[VIO_X]; voxel2[xyzv[VIO_X]]++)
      for(voxel2[xyzv[VIO_Y]]=start[VIO_Y]; voxel2[xyzv[VIO_Y]]<=end[VIO_Y]; voxel2[xyzv[VIO_Y]]++)
        for(voxel2[xyzv[VIO_Z]]=start[VIO_Z]; voxel2[xyzv[VIO_Z]]<=end[VIO_Z]; voxel2[xyzv[VIO_Z]]++) {
          
          if ((voxel2[ xyzv[VIO_X]] != voxel[ xyzv[VIO_X] ]) ||
              (voxel2[ xyzv[VIO_Y]] != voxel[ xyzv[VIO_Y] ]) ||
              (voxel2[ xyzv[VIO_Z]] != voxel[ xyzv[VIO_Z] ])) {

            for(voxel2[xyzv[VIO_Z+1]]=0; voxel2[xyzv[VIO_Z+1]]<sizes[xyzv[VIO_Z+1]]; voxel2[xyzv[VIO_Z+1]]++) {
              def_vector[voxel2[ xyzv[VIO_Z+1] ]] = 
                get_volume_real_value(volume,
                                      voxel2[0],voxel2[1],voxel2[2],voxel2[3],voxel2[4]);
            }
            *mx += def_vector[VIO_X]; *my += def_vector[VIO_Y]; *mz += def_vector[VIO_Z];
            ++count;

          }
        }
    break;
    
  }
  }    

  if (count>0) {                /* average deformation vector */
    *mx /= count; 
    *my /= count; 
    *mz /= count; 
    return(TRUE);
  }
  else {
    return(FALSE);
  }
}


VIO_BOOL get_average_warp_of_neighbours(VIO_General_transform *trans,
                                              int voxel[],
                                              VIO_Real mean_pos[])
{
  int       i;
  VIO_Real      voxel_real[VIO_MAX_DIMENSIONS],
            dx, dy, dz;
  VIO_Volume     volume;
  
  if (trans->type != GRID_TRANSFORM) {
    print_error_and_line_num("get_average_warp_of_neighbours not called with GRID_TRANSFORM",
                             __FILE__, __LINE__);
    return (FALSE);
  }

  volume = trans->displacement_volume;
  
  for(i=0; i<get_volume_n_dimensions(volume); i++ ) {
    voxel_real[i] = (VIO_Real)voxel[i];
  }
  convert_voxel_to_world(volume, voxel_real, 
                         &(mean_pos[VIO_X]),&(mean_pos[VIO_Y]),&(mean_pos[VIO_Z]) );

  if ( ! get_average_warp_vector_from_neighbours(trans, voxel, 1, &dx, &dy, &dz)) {
    return(FALSE);
  }
  else {
    mean_pos[VIO_X] += dx; mean_pos[VIO_Y] += dy; mean_pos[VIO_Z] += dz;
    return(TRUE);
  }

}

/* add additional to current, return
   answer in additional 

   I assume for this procedure that the displacement_volume
   for both the additional and current transform have the 
   same volumetric definition ie, both have the same dimension
   order and same length along each dimension.
*/

void add_additional_warp_to_current(VIO_General_transform *additional,
                                           VIO_General_transform *current,
                                           VIO_Real weight)
{
  int
    count[VIO_MAX_DIMENSIONS],
    count_current[VIO_MAX_DIMENSIONS],
    xyzv_additional[VIO_MAX_DIMENSIONS],
    xyzv_current[VIO_MAX_DIMENSIONS],
    index[VIO_MAX_DIMENSIONS],
    i;
  VIO_Real 
    additional_value, current_value;


  if (get_volume_n_dimensions(additional->displacement_volume) != 
      get_volume_n_dimensions(current->displacement_volume)) {
    print_error_and_line_num("add_additional_warp_to_current: warp dim error",
                             __FILE__, __LINE__);
  }

  get_volume_sizes(additional->displacement_volume, count);
  get_volume_sizes(current->displacement_volume, count_current);
  for(i=0; i<get_volume_n_dimensions(current->displacement_volume); i++) {
    if (count_current[i] != count[i]) {
      print_error_and_line_num("add_additional_warp_to_current: dim count error",
                               __FILE__, __LINE__);
    }
  }

  get_volume_XYZV_indices(additional->displacement_volume, xyzv_additional);
  get_volume_XYZV_indices(current->displacement_volume, xyzv_current);
  for(i=0; i<get_volume_n_dimensions(current->displacement_volume); i++) {
    if (xyzv_current[i] != xyzv_additional[i]) {
      print_error_and_line_num("add_additional_warp_to_current: dim match error",
                               __FILE__, __LINE__);
    }
  }

  for(i=0; i<VIO_MAX_DIMENSIONS; i++) index[i]=0;

  for(index[xyzv_additional[VIO_X]]=0; index[xyzv_additional[VIO_X]]<count[xyzv_additional[VIO_X]]; index[xyzv_additional[VIO_X]]++)
    for(index[xyzv_additional[VIO_Y]]=0; index[xyzv_additional[VIO_Y]]<count[xyzv_additional[VIO_Y]]; index[xyzv_additional[VIO_Y]]++)
      for(index[xyzv_additional[VIO_Z]]=0; index[xyzv_additional[VIO_Z]]<count[xyzv_additional[VIO_Z]]; index[xyzv_additional[VIO_Z]]++)
        for(index[xyzv_additional[VIO_Z+1]]=0; index[xyzv_additional[VIO_Z+1]]<count[xyzv_additional[VIO_Z+1]]; index[xyzv_additional[VIO_Z+1]]++) {

          additional_value = get_volume_real_value(
                                 additional->displacement_volume,
                                 index[0],index[1],index[2],index[3],index[4]);
          current_value = get_volume_real_value(
                                 current->displacement_volume,
                                 index[0],index[1],index[2],index[3],index[4]);

          additional_value = current_value + additional_value*weight;

          set_volume_real_value(additional->displacement_volume,
                                index[0],index[1],index[2],index[3],index[4],
                                additional_value);

          
        }

}
                                            


/*******************************************************************
  procedure: smooth_the_warp

    desc: this procedure will smooth the current warp stored in
          current and return the smoothed warp in smoothed

    meth: smoothing is accomplished by averaging the value of the 
          node's deformation vector with the mean deformation vector
          of it's neighbours.

          def'  = sw*mean  + (1-sw)*def

          where: sw   = smoothing_weight
                 mean = neighbourhood mean deformation
                 def  = estimate def for current node
*/

void smooth_the_warp(VIO_General_transform *smoothed,
                            VIO_General_transform *current,
                            VIO_Volume warp_mag, VIO_Real thres) 
{
  int
    count_smoothed[VIO_MAX_DIMENSIONS],
    count_current[VIO_MAX_DIMENSIONS],
    count_mag[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    xyzv_current[VIO_MAX_DIMENSIONS],
    xyzv_mag[VIO_MAX_DIMENSIONS],
    mag_index[VIO_MAX_DIMENSIONS],
    index[VIO_MAX_DIMENSIONS],
    start[VIO_MAX_DIMENSIONS], 
    end[VIO_MAX_DIMENSIONS],
    i;
  VIO_Real 
    value[3], 
    wx, wy, wz, 
    mx, my, mz,smoothing;
  VIO_progress_struct
    progress;
  
  
  if (get_volume_n_dimensions(smoothed->displacement_volume) != 
      get_volume_n_dimensions(current->displacement_volume)) {
    print_error_and_line_num("smooth_the_warp: warp dim error",
                             __FILE__, __LINE__);
  }
  
  get_volume_sizes(smoothed->displacement_volume, count_smoothed);
  get_volume_sizes(current->displacement_volume, count_current);
  for(i=0; i<get_volume_n_dimensions(current->displacement_volume); i++) {
    if (count_current[i] != count_smoothed[i]) {
      print_error_and_line_num("smooth_the_warp: dim count error",
                               __FILE__, __LINE__);
    }
  }
  
  get_volume_XYZV_indices(smoothed->displacement_volume, xyzv);
  get_volume_XYZV_indices(current->displacement_volume, xyzv_current);
  for(i=0; i<get_volume_n_dimensions(current->displacement_volume); i++) {
    if (xyzv_current[i] != xyzv[i]) {
      print_error_and_line_num("smooth_the_warp: dim match error",
                               __FILE__, __LINE__);
    }
  }
  
  get_volume_XYZV_indices(warp_mag, xyzv_mag);
  get_volume_sizes(warp_mag, count_mag);
  
  for(i=0; i<get_volume_n_dimensions(warp_mag); i++) {
    if (count_current[xyzv_current[i]] != count_mag[i]) {
      print_error_and_line_num("smooth_the_warp: dim count error w/mag (%d: %d != %d)\n",
                               __FILE__, __LINE__, i, count_current[xyzv_current[i]], count_mag[i] );
    }
  }
  
  for(i=0; i<VIO_MAX_DIMENSIONS; i++) {
    index[i]=0;
    start[i] = 0;
    end[i] = 0;
  }
  
  get_voxel_spatial_loop_limits(smoothed->displacement_volume, start, end);
  start[VIO_Z+1] = 0;
  end[VIO_Z+1] = 3;
  
  
  initialize_progress_report( &progress, FALSE, 
			      (end[VIO_X]-start[VIO_X])*
			      (end[VIO_Y]-start[VIO_Y]) + 1,
			      "Smoothing deformations" );
  
  
  for(index[xyzv[VIO_X]]=start[VIO_X]; index[xyzv[VIO_X]]<end[VIO_X]; index[xyzv[VIO_X]]++) {
    for(index[xyzv[VIO_Y]]=start[VIO_Y]; index[xyzv[VIO_Y]]<end[VIO_Y]; index[xyzv[VIO_Y]]++) {
      for(index[xyzv[VIO_Z]]=start[VIO_Z]; index[xyzv[VIO_Z]]<end[VIO_Z]; index[xyzv[VIO_Z]]++) {
	
        for(i=0; i<get_volume_n_dimensions(warp_mag); i++){
          mag_index[ xyzv_mag[i] ] = index[ xyzv[i] ];
	}
	
	
	/* go get the current warp vector for
	   this node. */
	
	for(index[xyzv[VIO_Z+1]]=start[VIO_Z+1]; index[xyzv[VIO_Z+1]]<end[VIO_Z+1]; index[xyzv[VIO_Z+1]]++) {
	  
	  value[index[ xyzv[VIO_Z+1] ]] = 
	    get_volume_real_value(current->displacement_volume,
				  index[0],index[1],index[2],
				  index[3],index[4]);
	  
	}
	/* store the current warp in wx, wy,wz */
	
	wx = value[VIO_X]; wy = value[VIO_Y]; wz = value[VIO_Z]; 
	
	index[ xyzv[VIO_Z+1]] = 0;
	
	/* if we can get a neighbourhood mean
	   warp vector, then we average it
	   with the current warp vector */
	
	if ( get_average_warp_vector_from_neighbours(current,
						     index, 2 ,
						     &mx, &my, &mz) ) {
	  
	  wx = (1.0 - smoothing_weight) * value[VIO_X] + smoothing_weight * mx;
	  wy = (1.0 - smoothing_weight) * value[VIO_Y] + smoothing_weight * my;
	  wz = (1.0 - smoothing_weight) * value[VIO_Z] + smoothing_weight * mz;
	  value[VIO_X] = wx; 
	  value[VIO_Y] = wy; 
	  value[VIO_Z] = wz; 
	  
	} 
          
                                /* now put the averaged vector into
                                   the smoothed volume */
 
	for(index[xyzv[VIO_Z+1]]=start[VIO_Z+1]; index[xyzv[VIO_Z+1]]<end[VIO_Z+1]; index[xyzv[VIO_Z+1]]++)  
	  set_volume_real_value(smoothed->displacement_volume,
				index[0],index[1],index[2],
				index[3],index[4],
				value[index[ xyzv[VIO_Z+1] ]] );  
      }

          
    }
    update_progress_report( &progress,
			    ((end[ VIO_Y ]-start[ VIO_Y ])*
			     (index[ xyzv[VIO_X]  ]-start[ VIO_X ])) +
			    (index[ xyzv[VIO_Y] ]-start[ VIO_X ])  +    1  );
  }


  terminate_progress_report( &progress );
}


/*

   We want to extrapolate (and smooth) the estimated deformations to
   nodes where no estimation was possible (and thus no local smoothing
   completed).  This procedure should only be called when
   Gglobals->trans_info.use_local_smoothing is true.  (When false,
   extrapolation is not needed, since it is addressed in the global
   smoothing process.)

   additional  = sw*mean    + (1-sw)current
             i          i-1                i-1 
               
   Here we must store the mean_vector location in additional_warp for
   all nodes where there was no estimatation possible (and thus no
   local smoothing completed).  This will achieve a homogeneous
   smoothing throughout the entire field.  (Here, smoothing_weight is
   effectively equal to 1.0 - since the local node has no information
   available to it.

   The idea is to set the node value in additional_vol, so that when
   added to the corresponding node value in current_vol, a smoothed
   vector field results.
   ie 
      additional = required_absolute_mean_vector - absolute_current_vector


   note estimated_flag_vol is created to be accessed in [VIO_X][VIO_Y][VIO_Z] order.

      */

void extrapolate_to_unestimated_nodes(VIO_General_transform *current,
                                             VIO_General_transform *additional,
                                             VIO_Volume estimated_flag_vol) 
{

  int 
    many,
    total,
    extrapolated,
    count,
    count_additional[VIO_MAX_DIMENSIONS],
    count_current[VIO_MAX_DIMENSIONS],
    count_flag[VIO_MAX_DIMENSIONS],
    xyzv[VIO_MAX_DIMENSIONS],
    xyzv_current[VIO_MAX_DIMENSIONS],
    xyzv_flag[VIO_MAX_DIMENSIONS],
    flag_index[VIO_MAX_DIMENSIONS],
    index[VIO_MAX_DIMENSIONS],
    start[VIO_MAX_DIMENSIONS], 
    end[VIO_MAX_DIMENSIONS],
    start2[VIO_MAX_DIMENSIONS], 
    end2[VIO_MAX_DIMENSIONS],
    voxel2[VIO_MAX_DIMENSIONS],
    i;
  VIO_Real 
    current_deform[VIO_N_DIMENSIONS], 
    additional_deform[VIO_N_DIMENSIONS], 
    mx, my, mz;
  VIO_progress_struct
    progress;

  extrapolated = many = total = 0;

                                /* verify that the deformation volumes for
                                   current and additional warps are
                                   compatible in dimension, size and
                                   order. */

  if (get_volume_n_dimensions(additional->displacement_volume) != 
      get_volume_n_dimensions(current->displacement_volume)) {
    print_error_and_line_num("extrapolate_the_warp: warp dim error",
                             __FILE__, __LINE__);
  }
  
  get_volume_sizes(additional->displacement_volume, count_additional);
  get_volume_sizes(current->displacement_volume, count_current);
  for(i=0; i<get_volume_n_dimensions(current->displacement_volume); i++) {
    if (count_current[i] != count_additional[i]) {
      print_error_and_line_num("extrapolate_the_warp: dim count error",
                               __FILE__, __LINE__);
    }
  }

  get_volume_XYZV_indices(additional->displacement_volume, xyzv);
  get_volume_XYZV_indices(current->displacement_volume, xyzv_current);
  for(i=0; i<get_volume_n_dimensions(current->displacement_volume); i++) {
    if (xyzv_current[i] != xyzv[i]) {
      print_error_and_line_num("extrapolate_the_warp: dim match error",
                               __FILE__, __LINE__);
    }
  }
  
  get_volume_XYZV_indices(estimated_flag_vol, xyzv_flag);
  get_volume_sizes(estimated_flag_vol, count_flag);

  for(i=0; i<get_volume_n_dimensions(estimated_flag_vol); i++) {
    if (count_current[xyzv_current[i]] != count_flag[i]) {
      print_error_and_line_num("extrapolate_the_warp: dim count error w/flag (%d: %d != %d)\n",
                               __FILE__, __LINE__, i, count_current[xyzv_current[i]], count_flag[i] );
    }
  }

                                /* verification completed, now get loop
                                   limits to go through the deformation
                                   volume and extrapolate the estimated
                                   vectors from the additional volume */
  for(i=0; i<VIO_MAX_DIMENSIONS; i++) {
    index[i] = 0;
    start[i] = 0;
    end[i]   = 0;
  }
  
  get_voxel_spatial_loop_limits(additional->displacement_volume, start, end);
  start[VIO_Z+1] = 0;
  end[VIO_Z+1]   = 3;
 
  initialize_progress_report( &progress, FALSE, 
                             (end[VIO_X]-start[VIO_X])*
                             (end[VIO_Y]-start[VIO_Y]) + 1,
                             "Extrapolating estimations" );


  for(index[xyzv[VIO_X]]=start[VIO_X]; index[xyzv[VIO_X]]<end[VIO_X]; index[xyzv[VIO_X]]++) {
    for(index[xyzv[VIO_Y]]=start[VIO_Y]; index[xyzv[VIO_Y]]<end[VIO_Y]; index[xyzv[VIO_Y]]++) {
      for(index[xyzv[VIO_Z]]=start[VIO_Z]; index[xyzv[VIO_Z]]<end[VIO_Z]; index[xyzv[VIO_Z]]++) {

        for(i=0; i<get_volume_n_dimensions(estimated_flag_vol); i++)
          flag_index[ xyzv_flag[i] ] = index[ xyzv[i] ];
        
        total++;


        if (  get_volume_real_value(estimated_flag_vol, 
                                    index[ xyzv[VIO_X] ],
                                    index[ xyzv[VIO_Y] ],
                                    index[ xyzv[VIO_Z] ],
                                    0, 0) < 1.0) {

          /* 
             then, this node has not been estimated at this iteration,
             so we have to interpolate a deformation value from
             neighbouring nodes, and then apply the proper local
             smoothing. 

             local smoothing is defined to be the average deformation of
             the 26 8-connected neighbours and is computed in the loop
             below
          */

          many++;
                                /* go get the current warp vector for
                                   this node. */
          
          for(index[xyzv[VIO_Z+1]]=start[VIO_Z+1]; index[xyzv[VIO_Z+1]]<end[VIO_Z+1]; index[xyzv[VIO_Z+1]]++) {

            current_deform[index[ xyzv[VIO_Z+1] ]] = 
              get_volume_real_value(current->displacement_volume,
                                    index[0],index[1],index[2],
                                    index[3],index[4]);

          }
                                /* get an average of the current additional 
                                   deformation */

          for(i=0; i<VIO_N_DIMENSIONS; i++) {
            start2[i] = index[ xyzv[i] ] - 1;
            if (start2[i]<0) start2[i]=0;
            end2[i] = index[ xyzv[i] ] + 1;
            if (end2[i]> count_flag[ xyzv_flag[i] ]-1) end2[i] =  count_flag[ xyzv_flag[i] ]-1;
          }
    
          for(i=0; i<VIO_N_DIMENSIONS; i++)   additional_deform[i]= 0.0;
          for(i=0; i<VIO_MAX_DIMENSIONS; i++) voxel2[i]=0;
          count = 0;

          for(voxel2[xyzv[VIO_X]]=start2[VIO_X]; voxel2[xyzv[VIO_X]]<=end2[VIO_X]; voxel2[xyzv[VIO_X]]++)
            for(voxel2[xyzv[VIO_Y]]=start2[VIO_Y]; voxel2[xyzv[VIO_Y]]<=end2[VIO_Y]; voxel2[xyzv[VIO_Y]]++)
              for(voxel2[xyzv[VIO_Z]]=start2[VIO_Z]; voxel2[xyzv[VIO_Z]]<=end2[VIO_Z]; voxel2[xyzv[VIO_Z]]++) {

                if (get_volume_real_value(estimated_flag_vol, 
                                    voxel2[ xyzv[VIO_X] ],
                                    voxel2[ xyzv[VIO_Y] ],
                                    voxel2[ xyzv[VIO_Z] ],
                                    0, 0) >= 0.5) {

                  if (!((voxel2[ xyzv[VIO_X]] == index[ xyzv[VIO_X] ]) &&
                        (voxel2[ xyzv[VIO_Y]] == index[ xyzv[VIO_Y] ]) &&
                        (voxel2[ xyzv[VIO_Z]] == index[ xyzv[VIO_Z] ]))   ) {
                    
                    for(voxel2[xyzv[VIO_Z+1]]=0; voxel2[xyzv[VIO_Z+1]]<VIO_N_DIMENSIONS; voxel2[xyzv[VIO_Z+1]]++) {
                      additional_deform[ voxel2[ xyzv[VIO_Z+1] ] ] += 
                        get_volume_real_value(additional->displacement_volume,
                                              voxel2[0],voxel2[1],voxel2[2],
                                              voxel2[3],voxel2[4]);
                    }
                    ++count;
                  }
                }

                
              }

          if (count>0) {
            extrapolated++;
            for(i=0; i<VIO_N_DIMENSIONS; i++)
              additional_deform[i] /= 26.0;
          }

                                /* if we can get a neighbourhood mean
                                   warp vector from the previous iterations, 
                                   then we average it with the previous warp 
                                   vector */
          index[ xyzv[VIO_Z+1]] = 0;
          if ( get_average_warp_vector_from_neighbours(current,
                                                       index, 2 ,
                                                       &mx, &my, &mz) ) {

            /* additional_deform += sw*mean + (1-sw)*current - current

               with sw = 0.5 gives: */
            
            additional_deform[VIO_X] += (mx - current_deform[VIO_X])/2.0; 
            additional_deform[VIO_Y] += (my - current_deform[VIO_Y])/2.0; 
            additional_deform[VIO_Z] += (mz - current_deform[VIO_Z])/2.0; 
          } 

          
                                /* now put the averaged vector into
                                   the additional volume */

          for(index[xyzv[VIO_Z+1]]=start[VIO_Z+1]; index[xyzv[VIO_Z+1]]<end[VIO_Z+1]; index[xyzv[VIO_Z+1]]++)  
            set_volume_real_value(additional->displacement_volume,
                                  index[0],index[1],index[2],
                                  index[3],index[4],
                                  additional_deform[index[ xyzv[VIO_Z+1] ]] );  
        }
          
      }
      update_progress_report( &progress,
                             ((end[ VIO_Y ]-start[ VIO_Y ])*
                              (index[ xyzv[VIO_X]  ]-start[ VIO_X ])) +
                             (index[ xyzv[VIO_Y] ]-start[ VIO_X ])  +    1  );
    }
  }

  terminate_progress_report( &progress );

  print ("There were %d out of %d extrapolated (%d left) (%d extrapolated)\n",many,total,total-many, extrapolated);

}


VIO_Real get_value_of_point_in_volume(VIO_Real xw, VIO_Real yw, VIO_Real zw, 
                                          VIO_Volume data)
     
{

  VIO_Real
    mag,
    xvox, yvox, zvox;
  
  PointR 
    voxel;

  convert_3D_world_to_voxel(data, xw, yw, zw,
                            &zvox, &yvox, &xvox);
  fill_Point( voxel, zvox, yvox, xvox );

  if (!trilinear_interpolant(data, &voxel, &mag)) 
    return(-DBL_MAX); 
  else 
    return(mag);
}



/*********************************************************************** 
   split the VIO_General_transform stored in 'total' to extract two parts,
      1 - the last warp (that will be optimized)
      2 - everything else,
   so that:

   total = all_until_last + last_warp 
*/

void split_up_the_transformation(VIO_General_transform *total,
                                        VIO_General_transform **all_until_last,
                                        VIO_General_transform **last_warp) 
{
  VIO_General_transform
    *tmp_trans;
  int i;

  ALLOC(*all_until_last,1);
                                /* copy the first transform from the global data struct  */
  copy_general_transform(get_nth_general_transform(total, 0),
                         *all_until_last);

                                /* copy and concat the rest of tem, stopping before the end */
  for(i=1; i<get_n_concated_transforms(total)-1; i++){
    ALLOC(tmp_trans,1);
    copy_general_transform(get_nth_general_transform(total, i), tmp_trans);
    concat_general_transforms(*all_until_last, tmp_trans, *all_until_last);
  }

  *last_warp = (VIO_General_transform *)NULL;
  for(i=0; i<get_n_concated_transforms(total); i++) {
    if (get_transform_type( get_nth_general_transform(total,i) ) 
               == GRID_TRANSFORM)
      *last_warp = get_nth_general_transform(total,i);
  }

}



/*   return the maximum value stored in the data volume */

static VIO_Real get_maximum_magnitude(VIO_Volume dxyz)
{

  VIO_Real max;

  max = get_volume_maximum_real_value(dxyz);

  return(max);
  
}

/***************************************************************************/
/*    set the threshold to be 10% of the maximum gradient magnitude        */
/*    for each source and target volumes                                   */

void  set_feature_value_threshold(VIO_Volume d1, 
                                         VIO_Volume d2,
                                         VIO_Real *global_thres1, 
                                         VIO_Real *global_thres2, 
                                         VIO_Real *threshold1, 
                                         VIO_Real *threshold2) 

{
  if (*global_thres1==0.0)
    *threshold1 = 0.10 * get_maximum_magnitude(d1);
  else
    *threshold1 = *global_thres1 * get_maximum_magnitude(d1);

  if (*global_thres2==0.0)
    *threshold2 = 0.10 * get_maximum_magnitude(d2);
  else
    *threshold2 = *global_thres2 * get_maximum_magnitude(d2);
}



void build_two_perpendicular_vectors(VIO_Real orig[], 
                                             VIO_Real p1[], 
                                             VIO_Real p2[])
{
  VIO_Vector
    v,v1,v2;
  VIO_Real
    len;
  fill_Vector(v, orig[VIO_X], orig[VIO_Y], orig[VIO_Z]);
  create_two_orthogonal_vectors( &v, &v1, &v2);

  len = MAGNITUDE( v1 );
  if (len>0) {
    p1[VIO_X] = Vector_x(v1) / len;
    p1[VIO_Y] = Vector_y(v1) / len;
    p1[VIO_Z] = Vector_z(v1) / len;
  }
  else
    print_error_and_line_num("Null length for vector normalization\n", 
                __FILE__, __LINE__);

  len = MAGNITUDE( v2 );
  if (len>0) {
    p2[VIO_X] = Vector_x(v2) / len;
    p2[VIO_Y] = Vector_y(v2) / len;
    p2[VIO_Z] = Vector_z(v2) / len;
  }
  else
    print_error_and_line_num("Null length for vector normalization\n", 
                __FILE__, __LINE__);
}

float xcorr_objective_with_def(VIO_Volume d1,
                                      VIO_Volume d2,
                                      VIO_Volume m1,
                                      VIO_Volume m2, 
                                      Arg_Data *globals)
{

  VectorR
    slice_dir,
    row_dir,
    col_dir,
    vector_step;

  PointR
    starting_position,
    slice,
    row,
    col,
    pos2,
    voxel;

  double
    tx,ty,tz;
  int
    r,c,s;

  VIO_Real
    sign_x,sign_y,sign_z,
    value1, value2;
  
  VIO_Real
    s1,s2,s3;                   /* to store the sums for f1,f2,f3 */
  float 
    result;                                /* the result */
  int 
    count1,count2;





  fill_Point( starting_position, globals->start[VIO_X], globals->start[VIO_Y], globals->start[VIO_Z]);


  s1 = s2 = s3 = 0.0;
  count1 = count2 = 0;

  if (globals->step[VIO_X] > 0 ) sign_x = 1.0; else sign_x = -1.0;
  if (globals->step[VIO_Y] > 0 ) sign_y = 1.0; else sign_y = -1.0;
  if (globals->step[VIO_Z] > 0 ) sign_z = 1.0; else sign_z = -1.0;

  for(s=0; s<=globals->count[VIO_Z]; s++) {

    SCALE_VECTOR( vector_step, globals->directions[VIO_Z], s*sign_z);
    ADD_POINT_VECTOR( slice, starting_position, vector_step );


    for(r=0; r<=globals->count[VIO_Y]; r++) {
      
      SCALE_VECTOR( vector_step, globals->directions[VIO_Y], r*sign_y);
      ADD_POINT_VECTOR( row, slice, vector_step );
      
      SCALE_POINT( col, row, 1.0); /* init first col position */
      for(c=0; c<=globals->count[VIO_X]; c++) {
        

        convert_3D_world_to_voxel(d1, Point_x(col), Point_y(col), Point_z(col), &tx, &ty, &tz);

  
        
        fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
        
        if (point_not_masked(m1, Point_x(col), Point_y(col), Point_z(col))) {
          
          if (INTERPOLATE_TRUE_VALUE( d1, &voxel, &value1 )) {

            count1++;

            DO_TRANSFORM(pos2, globals->trans_info.transformation, col);
            
            convert_3D_world_to_voxel(d2, Point_x(pos2), Point_y(pos2), Point_z(pos2), &tx, &ty, &tz);
            
            fill_Point( voxel, tx, ty, tz ); /* build the voxel POINT */
        
            if (point_not_masked(m2, Point_x(pos2), Point_y(pos2), Point_z(pos2))) {
              
              if (INTERPOLATE_TRUE_VALUE( d2, &voxel, &value2 )) {


                if (value1 > globals->threshold[0] && value2 > globals->threshold[1] ) {
                  
                  count2++;

                  s1 += value1*value2;
                  s2 += value1*value1;
                  s3 += value2*value2;
                  
                } 
                
              } /* if voxel in d2 */
            } /* if point in mask volume two */
          } /* if voxel in d1 */
        } /* if point in mask volume one */
        
        if (sign_x > 0) {
          ADD_POINT_VECTOR( col, col, globals->directions[VIO_X] );
        }
        else {
          SUB_POINT_VECTOR( col, col, globals->directions[VIO_X] );
        }
        
      } /* for c */
    } /* for r */
  } /* for s */
  
  result = 1.0 - s1 / (sqrt((double)s2)*sqrt((double)s3));
  
  if (globals->flags.debug) (void)print ("%7d %7d -> %10.8f\n",count1,count2,result);
  
  return (result);
  
}

